#include "NA6PRecoTiming.h"
#include <cstdio>
#include <cstdlib>

namespace NA6PRecoTiming
{
namespace
{
struct Totals {
  std::chrono::steady_clock::duration transport{}, tgeo{}, correction{}, candidateFit{};
  std::uint64_t steps = 0, fits = 0, fitFailures = 0;
  ~Totals()
  {
    if (!enabled()) {
      return;
    }
    auto seconds = [](auto d) { return std::chrono::duration<double>(d).count(); };
    std::fprintf(stderr, "NA6P_RECO_TIMING steps=%llu candidateFits=%llu fitFailures=%llu transport=%.6f s TGeo=%.6f s materialCorrection=%.6f s candidateFit=%.6f s\n",
                 static_cast<unsigned long long>(steps), static_cast<unsigned long long>(fits),
                 static_cast<unsigned long long>(fitFailures), seconds(transport), seconds(tgeo),
                 seconds(correction), seconds(candidateFit));
  }
};
Totals& totals() { static Totals value; return value; }
} // namespace

bool enabled()
{
  static const bool value = [] { const char* env = std::getenv("NA6P_RECO_TIMING"); return env && env[0] && env[0] != '0'; }();
  return value;
}
void add(Part part, std::chrono::steady_clock::duration elapsed)
{
  auto& t = totals();
  switch (part) {
    case Part::Transport: t.transport += elapsed; break;
    case Part::TGeo: t.tgeo += elapsed; break;
    case Part::MaterialCorrection: t.correction += elapsed; break;
    case Part::CandidateFit: t.candidateFit += elapsed; break;
  }
}
void countPropagationStep() { if (enabled()) { ++totals().steps; } }
void countCandidateFit(bool success)
{
  if (enabled()) {
    ++totals().fits;
    totals().fitFailures += !success;
  }
}
} // namespace NA6PRecoTiming
