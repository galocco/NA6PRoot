// Compare the simulated W step length with the MS-only backward extrapolation.
#include "NA6PLayoutParam.h"
#include "NA6PBeamParam.h"
#include "NA6PMCComposedLabel.h"
#include "NA6PMCEventHeader.h"
#include "NA6PRecoParam.h"
#include "NA6PTrack.h"
#include "NA6PMuonSpecModularHit.h"
#include "Propagator.h"
#include <TFile.h>
#include <TParticle.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

R__LOAD_LIBRARY(librecLib)

void compareMSMuonPlugTruth(const char* runDir, const char* outputDir = "", Long64_t maxEvents = -1)
{
  const std::string run = runDir;
  const std::string out = *outputDir ? outputDir : run + "/plug_truth";
  TFile reco((run + "/reco/TracksMuonSpec.root").c_str());
  TFile kine((run + "/sim/MCKine.root").c_str());
  auto* msTree = reco.Get<TTree>("tracksMuonSpec");
  auto* mcTree = kine.Get<TTree>("mckine");
  if (!msTree || !mcTree || msTree->GetEntries() != mcTree->GetEntries()) {
    throw std::runtime_error("Missing or inconsistent MS and MC trees");
  }

  na6p::conf::ConfigurableParam::updateFromFile(run + "/reco/na6pLayout.ini", "layout");
  na6p::conf::ConfigurableParam::updateFromFile(run + "/reco/na6pLayout.ini", "beam");
  na6p::conf::ConfigurableParam::updateFromFile(run + "/reco/na6pRecoParam.ini", "reco");
  if (!Propagator::loadField() || !Propagator::loadGeometry(run + "/sim/geometry.root")) {
    throw std::runtime_error("Cannot initialize field or geometry");
  }
  Propagator::PropOpt options;
  options.maxStep = NA6PRecoParam::Instance().maxPropagationStep;
  if (!std::isfinite(options.maxStep) || options.maxStep <= 0.f) {
    throw std::runtime_error("Invalid maximum propagation step");
  }
  const auto* propagator = Propagator::Instance();
  TTreeReader ms(msTree), mc(mcTree);
  TTreeReaderValue<std::vector<NA6PTrack>> tracks(ms, "MuonSpec");
  TTreeReaderValue<std::vector<NA6PMCComposedLabel>> labels(ms, "MuonSpecMCTruth");
  TTreeReaderValue<NA6PMCEventHeader> header(mc, "header");
  TTreeReaderValue<std::vector<TParticle>> particles(mc, "tracks");
  struct Slice {
    double lo, hi, rIn, rOut;
  };
  const auto& layout = NA6PLayoutParam::Instance();
  std::vector<Slice> slices;
  double z = layout.posZStartAbsorber;
  for (int i = 0; i < layout.nAbsorberSlices - layout.nAbsorberSlicesWall; ++i) {
    const double end = z + layout.thicknessAbsorber[i];
    if (layout.medAbsorberPlug[i] == "Tungsten" && layout.radPlug[i] > 0.) {
      slices.push_back({z, end, layout.radPlugInt[i], layout.radPlug[i]});
    }
    z = end;
  }
  const auto lengthInPlug = [&](const std::array<double, 3>& a, const std::array<double, 3>& b) {
    const double dx = b[0] - a[0], dy = b[1] - a[1], dz = b[2] - a[2], transverse2 = dx * dx + dy * dy;
    if (dz == 0.) {
      return 0.;
    }
    double length = 0.;
    for (const auto& sl : slices) {
      if (std::max(a[2], b[2]) <= sl.lo || std::min(a[2], b[2]) >= sl.hi) {
        continue;
      }
      const double t1 = (sl.lo - a[2]) / dz, t2 = (sl.hi - a[2]) / dz;
      const double tl = std::max(0., std::min(t1, t2)), th = std::min(1., std::max(t1, t2));
      const auto disk = [&](double r) {
        if (r == 0.) {
          return 0.;
        }
        if (transverse2 == 0.) {
          return a[0] * a[0] + a[1] * a[1] <= r * r ? th - tl : 0.;
        }
        const double center = -(dx * a[0] + dy * a[1]) / transverse2;
        const double x = a[0] + center * dx, y = a[1] + center * dy;
        const double radial = r * r - x * x - y * y;
        if (radial <= 0.) {
          return 0.;
        }
        const double dt = std::sqrt(radial / transverse2);
        return std::max(0., std::min(th, center + dt) - std::max(tl, center - dt));
      };
      length += std::hypot(dx, dy, dz) * (disk(sl.rOut) - disk(sl.rIn));
    }
    return length;
  };
  if (slices.empty() || !mcTree->GetBranch("muonPlugLengthCm")) {
    throw std::runtime_error("Need a tungsten plug and simulation with --record-muon-plug-length");
  }
  TTreeReaderValue<std::vector<float>> truthLength(mc, "muonPlugLengthCm");
  TFile hitFile((run + "/sim/HitsMuonSpecModular.root").c_str());
  auto* hitTree = hitFile.Get<TTree>("hitsMuonSpecModular");
  if (!hitTree || hitTree->GetEntries() != mcTree->GetEntries()) {
    throw std::runtime_error("Missing or inconsistent MS hits");
  }
  TTreeReader hitReader(hitTree);
  TTreeReaderValue<std::vector<NA6PMuonSpecModularHit>> hits(hitReader, "MuonSpecModular");
  if (gSystem->mkdir(out.c_str(), true) != 0 && gSystem->AccessPathName(out.c_str(), kWritePermission)) {
    throw std::runtime_error("Cannot create the plug truth output directory");
  }
  std::ofstream csv(out + "/tracks.csv");
  csv << std::setprecision(9)
      << "event,track_id,y_gen,e_gen_GeV,e_ms_mc_GeV,e_ms_fit_GeV,e_pv_fit_GeV,w_mc_cm,w_fit_cm,r_pv_fit_to_beam_cm\n";
  const double mass = PID::getMass(PID::Muon);
  const auto& beam = NA6PBeamParam::Instance();
  size_t used = 0, failed = 0;
  const Long64_t nEvents = maxEvents < 0 ? mcTree->GetEntries() : std::min(maxEvents, mcTree->GetEntries());
  for (Long64_t event = 0; event < nEvents; ++event) {
    if (!ms.Next() || !mc.Next() || !hitReader.Next() || tracks->size() != labels->size() ||
        truthLength->size() != particles->size()) {
      throw std::runtime_error("Inconsistent MS and MC data at entry " + std::to_string(event));
    }
    const double xPV = header->getVX(), yPV = header->getVY(), zPV = header->getVZ();
    for (size_t i = 0; i < tracks->size(); ++i) {
      const auto& label = labels->at(i);
      if (!label.isValid() || label.isFake() || label.getSourceID() != 0) {
        continue;
      }
      const int id = label.getTrackID();
      // Event IDs in composed labels are stored modulo 2^nbitsEvID.
      if (label.getEventID() != (header->getEventID() & NA6PMCComposedLabel::maskEvID) ||
          id < 0 || static_cast<size_t>(id) >= particles->size()) {
        throw std::runtime_error("Invalid MC track label");
      }
      const auto& particle = particles->at(id);
      const auto& track = tracks->at(i);
      if (std::abs(particle.GetPdgCode()) != 13 ||
          !(std::hypot(particle.Vx() - xPV, particle.Vy() - yPV, particle.Vz() - zPV) <= .001) ||
          !std::isfinite(particle.Energy()) || particle.Energy() < mass || !track.isValid() ||
          !(track.getZ() > zPV) || !std::isfinite(track.getP()) || track.getP() <= 0.) {
        continue;
      }
      const auto length = truthLength->at(id);
      if (!std::isfinite(length) || length < 0.f) {
        throw std::runtime_error("Invalid simulated plug length");
      }

      // Compare the MC and fitted momentum at a common plane.
      const NA6PMuonSpecModularHit* closest = nullptr;
      double distance = 1.;
      for (const auto& hit : *hits) {
        const double dz = std::abs(track.getZ() - hit.getZIn());
        if (hit.getTrackID() == id && hit.isEnteringStart() && hit.getPZIn() > 0. && dz < distance) {
          distance = dz;
          closest = &hit;
        }
      }
      if (!closest) {
        ++failed;
        continue;
      }
      NA6PTrackPar atHit = track;
      atHit.setPID(PID::Muon);
      if (!propagator->propagateToZ(atHit, closest->getZIn(), options)) {
        ++failed;
        continue;
      }
      const double eMC = std::hypot(std::hypot(closest->getPXIn(), closest->getPYIn(), closest->getPZIn()), mass);
      const double eFit = std::hypot(atHit.getP(), mass);

      // Integrate the chord through the tungsten slices in each propagation step.
      NA6PTrackPar state = track;
      state.setPID(PID::Muon);
      double wFit = 0.;
      bool ok = true;
      while (state.getZ() > zPV) {
        const auto start = state.getXYZ<double>();
        const float target = std::max(static_cast<float>(zPV), state.getZ() - options.maxStep);
        if (!propagator->propagateToZ(state, target, options)) {
          ok = false;
          break;
        }
        if (start[2] > slices.front().lo && target < slices.back().hi) {
          wFit += lengthInPlug(start, state.getXYZ<double>());
        }
      }
      if (!ok || !(state.getP() > 0.) || !std::isfinite(state.getP()) || !std::isfinite(eFit) || !std::isfinite(eMC)) {
        ++failed;
        continue;
      }
      const double y = std::asinh(particle.Pz() / std::hypot(particle.Pt(), mass));
      const double dzBeam = zPV - beam.zFocus;
      const double rBeam = std::hypot(state.getX() - (beam.meanX + dzBeam * beam.meanSlopeX),
                                      state.getY() - (beam.meanY + dzBeam * beam.meanSlopeY));
      csv << event << ',' << id << ',' << y << ',' << particle.Energy() << ',' << eMC << ',' << eFit << ','
          << std::hypot(state.getP(), mass) << ',' << length << ',' << std::max(0., wFit) << ',' << rBeam << '\n';
      ++used;
    }
  }
  csv.close();
  if (!csv || !used) {
    throw std::runtime_error("No usable comparisons or failed writing the CSV");
  }
  std::cout << "Compared " << used << " tracks; excluded " << failed << " for missing hits or failed propagation\n";
}
