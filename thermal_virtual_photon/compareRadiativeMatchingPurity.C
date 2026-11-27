#include "ThermalMuonMSAnalysis.h"
#include "NA6PMCComposedLabel.h"
#include "NA6PMatch.h"
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

R__LOAD_LIBRARY(librecLib)

namespace
{
struct MatchingInput {
  explicit MatchingInput(const std::string& path) : file(path.c_str())
  {
    tree = file.Get<TTree>("tracksMatching");
    if (!tree) {
      throw std::runtime_error("Missing matching tree in " + path);
    }
    reader.SetTree(tree);
  }
  TFile file;
  TTree* tree = nullptr;
  TTreeReader reader;
  TTreeReaderValue<std::vector<NA6PMatch>> matches{reader, "Matching"};
  TTreeReaderValue<std::vector<NA6PMCComposedLabel>> labels{reader, "MatchingMCTruth"};
};

struct PurityCount {
  Long64_t matched = 0, correct = 0;

  void add(bool isCorrect)
  {
    ++matched;
    correct += isCorrect;
  }

  double percent() const { return matched ? 100. * correct / matched : 0.; }
};

struct MatchingScore {
  PurityCount all, thermal;
  std::map<int, PurityCount> thermalByY;
};

struct PairChanges {
  Long64_t same = 0, changed = 0, lost = 0, newMatch = 0;
};

struct PairedPurity {
  // Background-level cross products of (correct off, matched off, correct on, matched on).
  std::array<std::array<double, 4>, 4> products{};

  void add(const PurityCount& off, const PurityCount& on)
  {
    const std::array<Long64_t, 4> counts{off.correct, off.matched, on.correct, on.matched};
    for (size_t i = 0; i < counts.size(); ++i) {
      for (size_t j = 0; j < counts.size(); ++j) {
        products[i][j] += static_cast<double>(counts[i]) * counts[j];
      }
    }
  }

  double standardErrorPercent(const PurityCount& off, const PurityCount& on, size_t groups) const
  {
    if (!off.matched || !on.matched || groups < 2) {
      return 0.;
    }
    const std::array<double, 4> gradients{
      -1. / off.matched, off.percent() / (100. * off.matched),
      1. / on.matched, -on.percent() / (100. * on.matched)};
    double variance = 0.;
    for (size_t i = 0; i < gradients.size(); ++i) {
      for (size_t j = 0; j < gradients.size(); ++j) {
        variance += gradients[i] * products[i][j] * gradients[j];
      }
    }
    return 100. * std::sqrt(std::max(0., variance) * groups / (groups - 1.));
  }
};

struct BackgroundCounts {
  std::array<PurityCount, 2> all, thermal;
};
}

void compareRadiativeMatchingPurity(const char* runDir)
{
  const std::string run = runDir;
  const auto out = run + "/radiative_matching_purity";
  TFile vtFile((run + "/reco/TracksVerTel.root").c_str());
  TFile msFile((run + "/reco/TracksMuonSpec.root").c_str());
  TFile mcFile((run + "/sim/MCKine.root").c_str());
  auto* vtTree = vtFile.Get<TTree>("tracksVerTel");
  auto* msTree = msFile.Get<TTree>("tracksMuonSpec");
  auto* mcTree = mcFile.Get<TTree>("mckine");
  if (!vtTree || !msTree || !mcTree || vtTree->GetEntries() != mcTree->GetEntries() ||
      msTree->GetEntries() != mcTree->GetEntries()) {
    throw std::runtime_error("Missing or inconsistent matching inputs");
  }
  const auto nEvents = mcTree->GetEntries();
  TTreeReader vt(vtTree), ms(msTree), mc(mcTree);
  TTreeReaderValue<std::vector<NA6PMCComposedLabel>> vtLabels(vt, "VerTelMCTruth"), msLabels(ms, "MuonSpecMCTruth");
  TTreeReaderValue<NA6PMCEventHeader> header(mc, "header");
  TTreeReaderValue<std::vector<TParticle>> particles(mc, "tracks");
  constexpr std::array<int, 3> thresholds{3, 4, 5};
  constexpr std::array<const char*, 2> models{"ionization", "radiative"};
  constexpr double dy = .2;
  const auto mass = PID::getMass(PID::Muon);
  std::array<std::array<std::unique_ptr<MatchingInput>, models.size()>, thresholds.size()> inputs;
  std::array<std::array<MatchingScore, models.size()>, thresholds.size()> scores;
  std::array<std::array<PairChanges, 2>, thresholds.size()> changes;
  std::map<int, std::array<BackgroundCounts, thresholds.size()>> byBackground;
  for (size_t cut = 0; cut < thresholds.size(); ++cut) {
    for (size_t model = 0; model < models.size(); ++model) {
      const std::string dir = out + "/" + models[model] + "_vt" + std::to_string(thresholds[cut]);
      inputs[cut][model] = std::make_unique<MatchingInput>(dir + "/TracksMatching.root");
      if (inputs[cut][model]->tree->GetEntries() != nEvents) {
        throw std::runtime_error("Unexpected matching event count in " + dir);
      }
    }
  }
  const auto genuine = [](const NA6PMCComposedLabel& label) { return label.isValid() && !label.isFake(); };
  for (Long64_t event = 0; event < nEvents; ++event) {
    if (!vt.Next() || !ms.Next() || !mc.Next()) {
      throw std::runtime_error("Cannot read original inputs");
    }
    const auto& generators = header->getGenHeaders();
    const auto provenance = std::find_if(generators.begin(), generators.end(),
                                         [](const auto& g) { return g.getInfo() == "EmbeddedBackground"; });
    if (provenance == generators.end() || !provenance->hasKey("entry")) {
      throw std::runtime_error("Missing provenance for the embedded background");
    }
    auto& background = byBackground[static_cast<int>(provenance->getValue("entry"))];
    std::optional<std::vector<bool>> thermal;
    std::vector<std::optional<int>> rapidityBin(msLabels->size());
    const auto thermalBin = [&](size_t msIndex) -> std::optional<int> {
      const auto& label = msLabels->at(msIndex);
      if (!genuine(label) || label.getSourceID() != 0) {
        return std::nullopt;
      }
      const int id = label.getTrackID();
      if (label.getEventID() != header->getEventID() || id < 0 || static_cast<size_t>(id) >= particles->size()) {
        throw std::runtime_error("Invalid MC label for MS track");
      }
      if (!thermal) {
        thermal = thermalMuonMS::findThermalMuons(*header, *particles);
      }
      if (!thermal->at(id)) {
        return std::nullopt;
      }
      if (!rapidityBin[msIndex]) {
        const auto& particle = particles->at(id);
        const double y = std::asinh(particle.Pz() / std::hypot(particle.Pt(), mass));
        if (!std::isfinite(y)) {
          throw std::runtime_error("Invalid generated muon rapidity");
        }
        rapidityBin[msIndex] = static_cast<int>(std::floor(y / dy));
      }
      return rapidityBin[msIndex];
    };
    for (size_t cut = 0; cut < thresholds.size(); ++cut) {
      std::array<std::vector<int>, models.size()> matchedVT;
      for (size_t model = 0; model < models.size(); ++model) {
        auto& input = *inputs[cut][model];
        if (!input.reader.Next() || input.matches->size() != input.labels->size()) {
          throw std::runtime_error("Cannot read matching event");
        }
        auto& score = scores[cut][model];
        auto& partner = matchedVT[model];
        partner.assign(msLabels->size(), -1);
        for (size_t matchID = 0; matchID < input.matches->size(); ++matchID) {
          const auto& match = input.matches->at(matchID);
          const int vi = match.getIndexVT(), mi = match.getIndexMS();
          if (vi < 0 || mi < 0 || static_cast<size_t>(vi) >= vtLabels->size() ||
              static_cast<size_t>(mi) >= msLabels->size() || partner[mi] >= 0) {
            throw std::runtime_error("Invalid or duplicate matching pair");
          }
          partner[mi] = vi;
          const auto& vl = vtLabels->at(vi);
          const auto& ml = msLabels->at(mi);
          if (input.labels->at(matchID) != ml) {
            throw std::runtime_error("Matching output is out of sync with MS input");
          }
          const bool correct = genuine(vl) && genuine(ml) && vl == ml;
          score.all.add(correct);
          background[cut].all[model].add(correct);
          const auto bin = thermalBin(mi);
          if (bin) {
            score.thermal.add(correct);
            background[cut].thermal[model].add(correct);
            score.thermalByY[*bin].add(correct);
          }
        }
      }
      for (size_t mi = 0; mi < msLabels->size(); ++mi) {
        const int off = matchedVT[0][mi], on = matchedVT[1][mi];
        if (off < 0 && on < 0) {
          continue;
        }
        for (size_t sample = 0; sample <= static_cast<size_t>(rapidityBin[mi].has_value()); ++sample) {
          auto& c = changes[cut][sample];
          if (off < 0) {
            ++c.newMatch;
          } else if (on < 0) {
            ++c.lost;
          } else if (off == on) {
            ++c.same;
          } else {
            ++c.changed;
          }
        }
      }
    }
    if ((event + 1) % 10000 == 0) {
      std::cout << "Events " << event + 1 << '/' << nEvents << std::endl;
    }
  }
  std::ofstream summary(out + "/purity.csv"), differential(out + "/thermal_purity_vs_y.csv"),
    transitions(out + "/pair_changes.csv"), delta(out + "/purity_delta.csv");
  summary << std::setprecision(10) << "min_vt_hits,model,sample,matched,correct,purity_percent\n";
  differential << std::setprecision(10) << "min_vt_hits,model,y_low,y_high,matched,correct,purity_percent\n";
  transitions << "min_vt_hits,sample,same_pair,changed_vt_partner,lost_match,new_match\n";
  delta << std::setprecision(10) << "min_vt_hits,sample,delta_pp,background_groups,paired_se_pp,ci95_low_pp,ci95_high_pp\n";
  for (size_t cut = 0; cut < thresholds.size(); ++cut) {
    for (size_t model = 0; model < models.size(); ++model) {
      const auto& score = scores[cut][model];
      for (const auto& [name, counts] : {std::pair{"all", score.all}, std::pair{"thermal", score.thermal}}) {
        summary << thresholds[cut] << ',' << models[model] << ',' << name << ',' << counts.matched << ','
                << counts.correct << ',' << counts.percent() << '\n';
      }
      for (const auto& [bin, counts] : score.thermalByY) {
        differential << thresholds[cut] << ',' << models[model] << ',' << bin * dy << ',' << (bin + 1) * dy << ','
                     << counts.matched << ',' << counts.correct << ',' << counts.percent() << '\n';
      }
    }
    for (size_t sample = 0; sample < changes[cut].size(); ++sample) {
      const auto name = sample == 0 ? "all" : "thermal";
      const auto& c = changes[cut][sample];
      transitions << thresholds[cut] << ',' << name << ',' << c.same << ','
                  << c.changed << ',' << c.lost << ',' << c.newMatch << '\n';
      const auto& off = sample == 0 ? scores[cut][0].all : scores[cut][0].thermal;
      const auto& on = sample == 0 ? scores[cut][1].all : scores[cut][1].thermal;
      const double diff = on.percent() - off.percent();
      PairedPurity paired;
      for (const auto& [id, counts] : byBackground) {
        const auto& pair = sample == 0 ? counts[cut].all : counts[cut].thermal;
        paired.add(pair[0], pair[1]);
      }
      const double error = paired.standardErrorPercent(off, on, byBackground.size());
      delta << thresholds[cut] << ',' << name << ',' << diff << ',' << byBackground.size() << ',' << error << ','
            << diff - 1.96 * error << ',' << diff + 1.96 * error << '\n';
    }
  }
  summary.close();
  differential.close();
  transitions.close();
  delta.close();
  if (!summary || !differential || !transitions || !delta) {
    throw std::runtime_error("Cannot write matching statistics");
  }
}
