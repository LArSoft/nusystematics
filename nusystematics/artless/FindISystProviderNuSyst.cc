#include "nusystematics/artless/response_helper.hh"

#include "fhiclcpp/ParameterSet.h"

#include <iomanip>
#include <iostream>

namespace cliopts {
std::string provider_name = "";
bool quiet = false;
bool dump_example_config = false;
} // namespace cliopts

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n" << std::endl;
  std::cout << "\t-?|--help          : Show this message.\n"
               "\t-p <provider_name> : Full name of provider to search for.\n"
               "\t-q                 : Quieten all output. Use returncode to\n"
               "\t                     determine if the provider exists.\n"
               "\t-d                 : Dump provider example documentation.\n"
            << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-p") {
      cliopts::provider_name = argv[++opt];
    } else if (std::string(argv[opt]) == "-q") {
      cliopts::quiet = true;
    } else if (std::string(argv[opt]) == "-d") {
      cliopts::dump_example_config = true;
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

int main(int argc, char const *argv[]) {
  HandleOpts(argc, argv);
  if (!cliopts::provider_name.size()) {
    std::cout << "[ERROR]: Expected to be passed a -p option." << std::endl;
    SayUsage(argv);
    exit(1);
  }

  fhicl::ParameterSet dummy_ps;
  dummy_ps.put("tool_type", cliopts::provider_name);

  std::unique_ptr<nusyst::IGENIESystProvider_tool> tool;
  try {
    tool = nusyst::make_instance(dummy_ps);
  } catch (nusyst::unknown_nusyst_systprovider const &e) {
    if (!cliopts::quiet) {
      std::cout << "[INFO]: Failed to find ISystProvider named "
                << std::quoted(cliopts::provider_name) << std::endl;
    }
    return 0;
  }

  if (!cliopts::quiet) {
    std::cout << "[INFO]: Found ISystProvider named "
              << std::quoted(cliopts::provider_name) << std::endl;
    if (cliopts::dump_example_config) {
      fhicl::ParameterSet ex_cfg;
      ex_cfg.put<fhicl::ParameterSet>("example_tool_config",
                                      tool->GetExampleToolConfiguration());
      std::cout << ex_cfg.to_indented_string() << std::endl;
    }
  }

  return 1;
}
