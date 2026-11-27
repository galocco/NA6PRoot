#include <TFile.h>
#include <TEfficiency.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

// Run matchingQuality.C first to produce matchingQuality.root from the three passes.
void fiveHitComparison(const char* runDir = ".")
{
  const std::string run = runDir;
  TFile input((run + "/matchingQuality.root").c_str(), "READ");
  if (input.IsZombie()) throw std::runtime_error("Cannot open matchingQuality.root");
  std::ofstream csv(run + "/fiveHitComparison.csv");
  if (!csv) throw std::runtime_error("Cannot create fiveHitComparison.csv");
  csv << "threshold,thermal_ms,matched_thermal_ms,correct_thermal_matches,all_matches,fake_matches,efficiency_percent,purity_percent,fake_rate_percent\n";
  std::cout << "threshold | matched thermal / thermal MS | efficiency | correct / thermal matched | purity | fake / all matches | fake rate\n";
  for (int threshold : {3, 4, 5}) {
    const std::string suffix = "_minVT" + std::to_string(threshold);
    auto* efficiency = input.Get<TEfficiency>(("efficiency" + suffix).c_str());
    auto* purity = input.Get<TEfficiency>(("purity" + suffix).c_str());
    auto* fake = input.Get<TEfficiency>(("fakeRate" + suffix).c_str());
    if (!efficiency || !purity || !fake) throw std::runtime_error("Missing results for threshold " + std::to_string(threshold));
    constexpr int bin = 3; // VT track has five hits.
    const auto thermal = efficiency->GetTotalHistogram()->GetBinContent(bin);
    const auto matched = efficiency->GetPassedHistogram()->GetBinContent(bin);
    const auto correct = purity->GetPassedHistogram()->GetBinContent(bin);
    const auto thermalMatches = purity->GetTotalHistogram()->GetBinContent(bin);
    const auto all = fake->GetTotalHistogram()->GetBinContent(bin);
    const auto fakes = fake->GetPassedHistogram()->GetBinContent(bin);
    if (matched != thermalMatches || all - fakes != correct)
      throw std::runtime_error("Inconsistent five-hit counts");
    const auto percent = [](double numerator, double denominator) {
      return denominator > 0 ? std::to_string(100. * numerator / denominator) : std::string{};
    };
    const auto efficiencyPercent = percent(matched, thermal);
    const auto purityPercent = percent(correct, thermalMatches);
    const auto fakePercent = percent(fakes, all);
    csv << threshold << ',' << thermal << ',' << matched << ',' << correct << ',' << all << ',' << fakes << ','
        << efficiencyPercent << ',' << purityPercent << ',' << fakePercent << '\n';
    std::cout << threshold << " | " << matched << '/' << thermal << " | " << (efficiencyPercent.empty() ? "n/a" : efficiencyPercent + "%")
              << " | " << correct << '/' << thermalMatches << " | " << (purityPercent.empty() ? "n/a" : purityPercent + "%")
              << " | " << fakes << '/' << all << " | " << (fakePercent.empty() ? "n/a" : fakePercent + "%") << '\n';
  }
}
