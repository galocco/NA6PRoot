// Optional reconstruction CPU profile. Set NA6P_RECO_TIMING=1 to print it at exit.
#ifndef NA6P_RECO_TIMING_H
#define NA6P_RECO_TIMING_H

#include <chrono>
#include <cstdint>

namespace NA6PRecoTiming
{
enum class Part { Transport, TGeo, MaterialCorrection, CandidateFit };

bool enabled();
void add(Part part, std::chrono::steady_clock::duration elapsed);
void countPropagationStep();
void countCandidateFit(bool success);

class Timer
{
 public:
  explicit Timer(Part part) : mPart(part), mActive(enabled())
  {
    if (mActive) {
      mStart = std::chrono::steady_clock::now();
    }
  }
  ~Timer()
  {
    if (mActive) {
      add(mPart, std::chrono::steady_clock::now() - mStart);
    }
  }
  Timer(const Timer&) = delete;
  Timer& operator=(const Timer&) = delete;

 private:
  Part mPart;
  bool mActive;
  std::chrono::steady_clock::time_point mStart{};
};
} // namespace NA6PRecoTiming

#endif
