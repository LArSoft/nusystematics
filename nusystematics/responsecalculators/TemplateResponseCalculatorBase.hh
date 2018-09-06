#ifndef nusystematics_RESPONSE_CALCULATORS_TEMPLATE_RESPONSE_BASE_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_TEMPLATE_RESPONSE_BASE_HH_SEEN

#include "systematicstools/interface/types.hh"

#include "systematicstools/utility/ROOTUtility.hh"
#include "systematicstools/utility/exceptions.hh"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

#ifndef NO_ART
#include "cetlib/filepath_maker.h"
#endif

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSpline.h"

// #define TemplateResponseCalculatorBase_DEBUG

namespace nusyst {
template <size_t N> struct THType {};
template <> struct THType<1> {
  typedef TH1D type;
  static size_t GetNbins(type *H) { return H->GetNbinsX(); }
  static size_t GetNbins(std::unique_ptr<type> &H) { return H->GetNbinsX(); }
};
template <> struct THType<2> {
  typedef TH2D type;
  static size_t GetNbins(type *H) { return H->GetNbinsX() * H->GetNbinsY(); }
  static size_t GetNbins(std::unique_ptr<type> &H) {
    return H->GetNbinsX() * H->GetNbinsY();
  }
};
template <> struct THType<3> {
  typedef TH3D type;
  static size_t GetNbins(type *H) {
    return H->GetNbinsX() * H->GetNbinsY() * H->GetNbinsZ();
  }
  static size_t GetNbins(std::unique_ptr<type> &H) {
    return H->GetNbinsX() * H->GetNbinsY() * H->GetNbinsZ();
  }
};

NEW_SYSTTOOLS_EXCEPT(no_responses_loaded);
NEW_SYSTTOOLS_EXCEPT(incompatible_number_of_bins);

template <size_t NDims, bool Continuous = true>
class TemplateResponseCalculatorBase {

protected:
  std::map<systtools::paramId_t, std::vector<TSpline3>> SplinedBinResponses;
  std::map<systtools::paramId_t,
           std::map<double, std::unique_ptr<typename THType<NDims>::type>>>
      BinnedResponses;

  void ValidateInputHistograms();
  /// Reads and loads input fhicl
  ///
  /// Expected fhicl like:
  /// > # Optional default root file name
  /// > InputTemplatesFile: "file.root"
  /// > InputTemplates: [
  /// >   { parameter_name: a
  /// >     # Optional default root file name for this parameters inputs
  /// >     input_file: "file.root"
  /// >     inputs: [
  /// >       { value: 0
  /// >         # Optional if either of the less-specific
  /// >         input_file: "file.root"
  /// >         hist_name: "histo_name"
  /// >       }
  /// >     ]
  /// >   }
  /// > ]
  void LoadInputHistograms(
      fhicl::ParameterSet const &ps,
      std::map<std::string, systtools::paramId_t> const &ParamNames);

public:
  TemplateResponseCalculatorBase();

  typedef size_t bin_it_t;

  bin_it_t kBinOutsideRange = std::numeric_limits<bin_it_t>::max();

  virtual bin_it_t GetBin(systtools::paramId_t,
                          std::array<double, NDims> const &) = 0;

  virtual std::string GetCalculatorName() = 0;

  /// Extra template parameters required for SINFAE
  template <bool Enable = Continuous>
  typename std::enable_if<Enable, double>::type
  GetVariation(systtools::paramId_t, double,
               typename std::enable_if<Enable, bin_it_t>::type);
  /// Extra template parameters required for SINFAE
  template <bool Enable = Continuous>
  typename std::enable_if<!Enable, double>::type
  GetVariation(systtools::paramId_t, double,
               typename std::enable_if<!Enable, bin_it_t>::type);

  double GetVariation(systtools::paramId_t param, double val,
                      std::array<double, NDims> const &kinematics);

  /// Extra template parameters required for SINFAE
  template <bool Enable = Continuous>
  typename std::enable_if<Enable, std::pair<double, double>>::type
      GetValidVariations(
          typename std::enable_if<Enable, systtools::paramId_t>::type);
  /// Extra template parameters required for SINFAE
  template <bool Enable = Continuous>
  typename std::enable_if<!Enable, std::vector<double>>::type
      GetValidVariations(
          typename std::enable_if<!Enable, systtools::paramId_t>::type);
};

//*********************** Implementations

template <size_t NDims, bool Continuous>
void TemplateResponseCalculatorBase<NDims,
                                    Continuous>::ValidateInputHistograms() {
  for (auto &param_resps : BinnedResponses) {
    if (!param_resps.second.size()) {
      throw no_responses_loaded()
          << "[ERROR]: Expected to find some responses to parameter "
          << param_resps.first << " but found none.";
    }
    size_t NBins = THType<NDims>::GetNbins(param_resps.second.begin()->second);
    for (auto &val_resp : param_resps.second) {
      if (THType<NDims>::GetNbins(val_resp.second) != NBins) {
        throw incompatible_number_of_bins()
            << "[ERROR]: The first histogram at parameter value "
            << param_resps.second.begin()->first << " has a response in "
            << NBins << " bins, at parameter value " << val_resp.first
            << " found " << THType<NDims>::GetNbins(val_resp.second)
            << " bins.";
      }
    }
  }
}

template <size_t NDims, bool Continuous>
void TemplateResponseCalculatorBase<NDims, Continuous>::LoadInputHistograms(
    fhicl::ParameterSet const &ps,
    std::map<std::string, systtools::paramId_t> const &ParamNames) {

  std::string provider_default_root_file =
      ps.get<std::string>("InputTemplatesFile", "");

  std::vector<fhicl::ParameterSet> const &provider_config =
      ps.get<std::vector<fhicl::ParameterSet>>("InputTemplates");

  for (auto const &param_config : provider_config) {
    std::string pname = param_config.get<std::string>("parameter_name");

    if (ParamNames.find(pname) == ParamNames.end()) {
      std::stringstream ss("");
      ss << "[";
      for (auto const &p : ParamNames) {
        ss << p.second << ", ";
      }
      std::string valid_pnames = ss.str();
      valid_pnames = valid_pnames.substr(0, valid_pnames.size() - 2) + " ]";
      throw systtools::invalid_parameter_name()
          << "[ERROR]: Template reweight requested for parameter named "
          << std::quoted(pname) << ", but " << GetCalculatorName()
          << " knows about parameters: " << valid_pnames;
    }

    systtools::paramId_t pId = ParamNames.at(pname);

    std::string parameter_default_root_file =
        param_config.get<std::string>("input_file", provider_default_root_file);

    std::vector<fhicl::ParameterSet> const &inputs_config =
        param_config.get<std::vector<fhicl::ParameterSet>>("inputs");

    for (auto const &val_config : inputs_config) {
      double pval = val_config.get<double>("value");
      std::string input_file = val_config.get<std::string>(
          "input_file", parameter_default_root_file);
      std::string hist_name = val_config.get<std::string>("hist_name");

      BinnedResponses[pId][pval] =
          std::unique_ptr<typename THType<NDims>::type>(
              GetHistogram<typename THType<NDims>::type>(input_file,
                                                         hist_name));
    }
  }
}
template <size_t NDims, bool Continuous>
TemplateResponseCalculatorBase<NDims,
                               Continuous>::TemplateResponseCalculatorBase() {}

template <size_t NDims, bool Continuous>
template <bool Enable>
typename std::enable_if<Enable, double>::type
TemplateResponseCalculatorBase<NDims, Continuous>::GetVariation(
    systtools::paramId_t param, double val,
    typename std::enable_if<Enable, bin_it_t>::type bin) {
  return SplinedBinResponses[param][bin]->Eval(val);
}

template <size_t NDims, bool Continuous>
template <bool Enable>
typename std::enable_if<!Enable, double>::type
TemplateResponseCalculatorBase<NDims, Continuous>::GetVariation(
    systtools::paramId_t param, double val,
    typename std::enable_if<!Enable, bin_it_t>::type bin) {

  if (bin == kBinOutsideRange) {
    return 1;
  }

  for (auto &resp : BinnedResponses[param]) {
    if (fabs(val - resp.first) <
        (std::numeric_limits<double>::epsilon() * 1E4)) {
#ifdef TemplateResponseCalculatorBase_DEBUG
      std::cout << "[INFO]: Getting bin content for bin: " << bin
                << " for parameter: " << param << " at value: " << val << " = "
                << resp.second->GetBinContent(bin)
                << " from hist: " << resp.second->GetName() << std::endl;
#endif
      return resp.second->GetBinContent(bin);
    }
  }

  std::stringstream ss("");
  ss << "[";
  for (auto const &v : GetValidVariations(param)) {
    ss << v << ", ";
  }
  std::string valid_vals = ss.str();
  valid_vals = valid_vals.substr(0, valid_vals.size() - 2) + " ]";
  throw systtools::invalid_parameter_value()
      << "[ERROR]: Invalid parameter value, " << val
      << " used for template response " << GetCalculatorName()
      << ", configured values: " << valid_vals;
}

template <size_t NDims, bool Continuous>
double TemplateResponseCalculatorBase<NDims, Continuous>::GetVariation(
    systtools::paramId_t param, double val,
    std::array<double, NDims> const &kinematics) {
  return GetVariation(param, val, GetBin(param, kinematics));
}

template <size_t NDims, bool Continuous>
template <bool Enable>
typename std::enable_if<Enable, std::pair<double, double>>::type
TemplateResponseCalculatorBase<NDims, Continuous>::GetValidVariations(
    typename std::enable_if<Enable, systtools::paramId_t>::type param) {
  std::pair<double, double> minmax{std::numeric_limits<double>::max(),
                                   -std::numeric_limits<double>::max()};
  for (auto const &var : BinnedResponses[param]) {
    minmax.first = std::min(minmax.first, var.first);
    minmax.second = std::max(minmax.second, var.first);
  }

  return minmax;
}

template <size_t NDims, bool Continuous>
template <bool Enable>
typename std::enable_if<!Enable, std::vector<double>>::type
TemplateResponseCalculatorBase<NDims, Continuous>::GetValidVariations(
    typename std::enable_if<!Enable, systtools::paramId_t>::type param) {
  std::vector<double> vals;
  for (auto const &var : BinnedResponses[param]) {
    vals.push_back(var.first);
  }
  return vals;
}
} // namespace nusyst

#endif
