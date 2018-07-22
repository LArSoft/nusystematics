#include "larsyst/interface/ISystProvider_tool.hh"

using namespace larsyst;
using namespace fhicl;

class BeRPAWeightProvider : public ISystProvider_tool {
public:
  explicit BeRPAWeightProvider(ParameterSet const &);

  SystMetaData ConfigureFromFHICL(ParameterSet const &, paramId_t);

  bool Configure();
  std::unique_ptr<EventResponse> GetEventResponse(art::Event &);
  std::string AsString();

private:
  enum Coeffs { kA = 0, kB, kD, kE, kU, kNCoeffs };
  std::map<Coeffs, std::pair<double, double>> DefaultValues;
  std::map<Coeffs, paramId_t> VariedParameters;
  size_t responseParam_idx;
};

BeRPAWeightProvider::BeRPAWeightProvider(ParameterSet const &params)
    : ISystProvider_tool(params) {
  DefaultValues[kA] = {0.59, 0.59 * 0.20};
  DefaultValues[kB] = {1.05, 1.05 * 0.20};
  DefaultValues[kD] = {1.13, 1.13 * 0.15};
  DefaultValues[kE] = {0.88, 0.88 * 0.40};
  DefaultValues[kE] = {1.2, 0};

  responseParam_idx = kParamUnhandled<size_t>;
}

SystMetaData BeRPAWeightProvider::ConfigureFromFHICL(ParameterSet const &ps,
                                                     paramId_t firstId) {

  SystMetaData smd;

  SystParamHeader resp;
  resp.prettyName = "BeRPA_Response";
  // Don't increment this, if we use it, increment all the others.
  resp.systParamId = firstId;

  if (ps.has_key("ATweakDefinition")) {
    double ACV = kDefaultDouble;
    if (ps.has_key("ACentralValue")) {
      ACV = ps.get<double>("ACentralValue");
    }
    SystParamHeader A = init_header_from_tweak_definition_string(
        ACV, ps.get<std::string>("ATweakDefinition"));
    A.prettyName = "BeRPA_A";
    A.systParamId = firstId++;
    smd.push_back(std::move(A));
  }
  if (ps.has_key("BTweakDefinition")) {
    double BCV = kDefaultDouble;
    if (ps.has_key("BCentralValue")) {
      BCV = ps.get<double>("BCentralValue");
    }
    SystParamHeader B = init_header_from_tweak_definition_string(
        BCV, ps.get<std::string>("BTweakDefinition"));
    B.prettyName = "BeRPA_B";
    B.systParamId = firstId++;
    smd.push_back(std::move(B));
  }
  if (ps.has_key("DTweakDefinition")) {
    double DCV = kDefaultDouble;
    if (ps.has_key("DCentralValue")) {
      DCV = ps.get<double>("DCentralValue");
    }
    SystParamHeader D = init_header_from_tweak_definition_string(
        DCV, ps.get<std::string>("DTweakDefinition"));
    D.prettyName = "BeRPA_D";
    D.systParamId = firstId++;
    smd.push_back(std::move(D));
  }
  if (ps.has_key("ETweakDefinition")) {
    double ECV = kDefaultDouble;
    if (ps.has_key("ECentralValue")) {
      ECV = ps.get<double>("ECentralValue");
    }
    SystParamHeader E = init_header_from_tweak_definition_string(
        ECV, ps.get<std::string>("ETweakDefinition"));
    E.prettyName = "BeRPA_E";
    E.systParamId = firstId++;
    smd.push_back(std::move(E));
  }

  // If multiple parameters were used then add the response parameter and
  // increment the id of all of the other parameters.
  if (smd.size() > 1) {
    for (auto &h : smd) {
      h.systParamId++;
      h.isResponselessParam = true;
      h.responseParamId = resp.systParamId;
    }

    smd.insert(0, std::move(resp));
  }

  return smd;
}

bool BeRPAWeightProvider::Configure() {

  paramId_t AId = GetParamId(fMetaData, "BeRPA_A");
  paramId_t BId = GetParamId(fMetaData, "BeRPA_B");
  paramId_t DId = GetParamId(fMetaData, "BeRPA_D");
  paramId_t EId = GetParamId(fMetaData, "BeRPA_E");
  paramId_t RespId = GetParamId(fMetaData, "BeRPA_Response");

  if (AId != kParamUnhandled<paramId_t>) {
    VariedParameters[kA] = GetParamIndex(fMetaData, "BeRPA_A");
    responseParam_idx = VariedParameters[kA];
  }
  if (BId != kParamUnhandled<paramId_t>) {
    VariedParameters[kB] = GetParamIndex(fMetaData, "BeRPA_B");
    responseParam_idx = VariedParameters[kB];
  }
  if (DId != kParamUnhandled<paramId_t>) {
    VariedParameters[kD] = GetParamIndex(fMetaData, "BeRPA_D");
    responseParam_idx = VariedParameters[kD];
  }
  if (EId != kParamUnhandled<paramId_t>) {
    VariedParameters[kE] = GetParamIndex(fMetaData, "BeRPA_E");
    responseParam_idx = VariedParameters[kE];
  }
  if (RespId != kParamUnhandled<paramId_t>) {
    VariedParameters[kNCoeffs] = GetParamIndex(fMetaData, "BeRPA_Response");
    responseParam_idx = VariedParameters[kNCoeffs];
  }
}

std::unique_ptr<EventResponse>
BeRPAWeightProvider::GetEventResponse(art::Event &e) {

  art::Handle<std::vector<simb::GTruth>> gTruthHandle;
  e.getByLabel(fGENIEModuleLabel, gTruthHandle);

  std::unique_ptr<EventResponse> er = std::make_unique<EventResponse>();

  size_t NEventUnits = gTruthHandle->size();
  er->responses.resize(NEventUnits);
  size_t NParamVars =
      fMetaData[responseParam_idx].paramVariations.size();

  std::map<Coeffs, double> pVals;
  for (size_t ev_it = 0; ev_it < NEventUnits; ++ev_it) {
    std::vector<double> RPAWeights;
    for (size_t v_it = 0; v_it < NParamVars; ++v_it) {
      if (VariedParameters.find(kA) != VariedParameters.end()) {
        pVals[kA] =
            fMetaData[VariedParameters[kA]].paramVariations[v_it];
      } else {
        pVals[kA] = DefaultValues[kA].first;
      }
      if (VariedParameters.find(kB) != VariedParameters.end()) {
        pVals[kB] =
            fMetaData[VariedParameters[kB]].paramVariations[v_it];
      } else {
        pVals[kB] = DefaultValues[kB].first;
      }
      if (VariedParameters.find(kD) != VariedParameters.end()) {
        pVals[kD] =
            fMetaData[VariedParameters[kD]].paramVariations[v_it];
      } else {
        pVals[kD] = DefaultValues[kD].first;
      }
      if (VariedParameters.find(kE) != VariedParameters.end()) {
        pVals[kE] =
            fMetaData[VariedParameters[kE]].paramVariations[v_it];
      } else {
        pVals[kE] = DefaultValues[kE].first;
      }
      RPAWeights.push_back(nusyst::BeRPA(gTruthHandle->at(ev_it).fgQ2,
                                         pVals[kA], pVals[kB], pVals[kD],
                                         pVals[kE]));
    }
    er->responses[ev_it].push_back(
        {fMetaData[responseParam_idx].systParamId,
         std::move(RPAWeights)});
  }
  return er;
}
std::string BeRPAWeightProvider::AsString() { return ""; }
