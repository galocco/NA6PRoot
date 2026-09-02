// NA6PCCopyright

#ifndef NA6P_TOF_MATCH_H
#define NA6P_TOF_MATCH_H

#include "NA6PTrackParCov.h"

#include <Rtypes.h>
#include <string>

class NA6PTOFMatch : public NA6PTrackParCov
{
 public:
  using NA6PTrackParCov::NA6PTrackParCov;

  NA6PTOFMatch() = default;
  explicit NA6PTOFMatch(const NA6PTrackParCov& track) : NA6PTrackParCov(track) {}

  float getTOF() const { return mTOF; }
  void setTOF(float time) { mTOF = time; }

  float getPathLength() const { return mPathLength; }
  void setPathLength(float length) { mPathLength = length; }
  float getBeta() const
  {
    constexpr float SpeedOfLightCmPerSecond = 2.99792458e10f;
    return mTOF > 0.f ? mPathLength / (SpeedOfLightCmPerSecond * mTOF) : 0.f;
  }

  float getMatchChi2() const { return mMatchChi2; }
  void setMatchChi2(float chi2) { mMatchChi2 = chi2; }

  int getIndexTrack() const { return mIndexTrack; }
  void setIndexTrack(int index) { mIndexTrack = index; }

  int getIndexTOF() const { return mIndexTOF; }
  void setIndexTOF(int index) { mIndexTOF = index; }

  void print() const;
  std::string asString() const;

 private:
  float mTOF = 0.f;           // measured arrival time [s]
  float mPathLength = 0.f;    // helical path length from the primary-vertex Z to the TOF [cm]
  float mMatchChi2 = 0.f;     // position matching chi2 (X,Y)
  int mIndexTrack = -1;       // index of the source track
  int mIndexTOF = -1;         // index of the matched TOF cluster

  ClassDefNV(NA6PTOFMatch, 2);
};

#endif
