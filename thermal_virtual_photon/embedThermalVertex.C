#include "NA6PMCEventHeader.h"
#include <TFile.h>
#include <TTree.h>
#include <cstdlib>
#include <stdexcept>
#include <string>

void embedThermalVertex(float& x, float& y, float& z)
{
  const auto requiredEnv = [](const char* name) {
    const auto* value = std::getenv(name);
    if (!value || !*value) throw std::runtime_error(std::string("Missing embedding variable ") + name);
    return value;
  };
  static TFile background(requiredEnv("NA6P_EMBED_BACKGROUND"));
  static auto* tree = background.Get<TTree>("mckine");
  static auto first = std::stoll(requiredEnv("NA6P_EMBED_FIRST"));
  static const auto replicas = std::stoll(requiredEnv("NA6P_EMBED_REPLICAS"));
  static NA6PMCEventHeader* header = nullptr;
  if (!tree || tree->GetEntries() == 0 || first < 0 || replicas <= 0) {
    throw std::runtime_error("Invalid embedding background or replica count");
  }
  if (!header && tree->SetBranchAddress("header", &header) < 0) {
    throw std::runtime_error("Missing background event header");
  }
  if (tree->GetEntry(first++ / replicas % tree->GetEntries()) <= 0 || !header) {
    throw std::runtime_error("Cannot read background vertex");
  }
  x = header->getVX();
  y = header->getVY();
  z = header->getVZ();
}
