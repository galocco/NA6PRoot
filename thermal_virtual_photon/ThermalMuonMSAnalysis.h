#ifndef THERMAL_MUON_MS_ANALYSIS_H
#define THERMAL_MUON_MS_ANALYSIS_H

#include "NA6PLayoutParam.h"
#include "NA6PMCEventHeader.h"
#include "NA6PRecoParam.h"
#include "Propagator.h"
#include <TParticle.h>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace thermalMuonMS
{
inline Propagator::PropOpt loadPropagation(const std::string& run)
{
  na6p::conf::ConfigurableParam::updateFromFile(run + "/reco/na6pLayout.ini", "layout");
  na6p::conf::ConfigurableParam::updateFromFile(run + "/reco/na6pRecoParam.ini", "reco");
  if (!Propagator::loadField() || !Propagator::loadGeometry(run + "/sim/geometry.root")) {
    throw std::runtime_error("Cannot initialize field or geometry for propagation");
  }
  Propagator::PropOpt propagation;
  propagation.maxStep = NA6PRecoParam::Instance().maxPropagationStep;
  if (!std::isfinite(propagation.maxStep) || propagation.maxStep <= 0.f) {
    throw std::runtime_error("Invalid maximum propagation step");
  }
  return propagation;
}

inline std::vector<bool> findThermalMuons(const NA6PMCEventHeader& header, const std::vector<TParticle>& particles)
{
  std::vector<bool> thermal(particles.size());
  bool hasThermalSource = false;
  for (const auto& generator : header.getGenHeaders()) {
    const auto& name = generator.getInfo();
    if (name.find("VirtualPhoton_pdg23") != 0 &&
        name.find("thermalMuPlus_pdg-13") != 0 &&
        name.find("thermalMuMinus_pdg13") != 0) {
      continue;
    }
    hasThermalSource = true;
    const size_t begin = generator.getPrimariesOffset();
    const size_t count = generator.getNPrimaries();
    if (begin > thermal.size() || count > thermal.size() - begin) {
      throw std::runtime_error("Thermal generator range exceeds MC tracks");
    }
    for (size_t i = begin; i < begin + count; ++i) {
      thermal[i] = std::abs(particles[i].GetPdgCode()) == 13;
    }
  }
  if (!hasThermalSource) {
    throw std::runtime_error("Missing thermal generator header");
  }
  return thermal;
}
}

#endif
