// NA6PCCopyright

#ifndef NA6P_TOF_MATCHING_H
#define NA6P_TOF_MATCHING_H

#include "NA6PReconstruction.h"
#include "NA6PMatch.h"
#include "NA6PTOFCluster.h"
#include "NA6PTOFMatch.h"
#include "NA6PTrack.h"

#include <Rtypes.h>
#include <vector>

class NA6PTOFHit;
class TFile;
class TTree;

class NA6PTOFMatching : public NA6PReconstruction
{
 public:
  explicit NA6PTOFMatching(bool init = true);
  ~NA6PTOFMatching() override = default;

  bool initTOFMatching();
  void configureFromRecoParam();

  void setClusterSpaceResolutionX(float res) { mCluResX = res; }
  void setClusterSpaceResolutionY(float res) { mCluResY = res; }
  void setClusterTimeResolution(float res) { mTimeRes = res; }
  void hitsToRecPoints(const std::vector<NA6PTOFHit>& hits);

  void createClustersOutput() override;
  void clearClusters() override { mClusters.clear(); }
  void writeClusters() override;
  void closeClustersOutput() override;

  void createTracksOutput() override;
  void clearTracks() override { mMatchedTracks.clear(); }
  void writeTracks() override;
  void closeTracksOutput() override;

  void setTracks(const std::vector<NA6PTrackParCov>& tracks);
  void setMatchedTracks(const std::vector<NA6PMatch>& tracks);
  void setVerTelTracks(const std::vector<NA6PTrack>& tracks);
  void setClusters(std::vector<NA6PTOFCluster>& clusters);
  void setMaxMatchChi2(float chi2) { mMaxMatchChi2 = chi2; }
  void runTOFMatching();

  const std::vector<NA6PTOFCluster>& getClusters() const { return mClusters; }
  const std::vector<NA6PTOFMatch>& getTracks() const { return mMatchedTracks; }

 private:
  std::vector<NA6PTOFCluster> mClusters;
  std::vector<NA6PTOFCluster>* hClusPtr = &mClusters;
  TFile* mClusFile = nullptr;
  TTree* mClusTree = nullptr;
  float mCluResX = 100.e-4f / 3.464101615f;
  float mCluResY = 300.e-4f / 3.464101615f;
  float mTimeRes = 20.e-12f; // 20 ps from ALICE3 TOF

  std::vector<NA6PTOFMatch> mMatchedTracks;
  std::vector<NA6PTOFMatch>* hMatchedTrackPtr = &mMatchedTracks;
  TFile* mTrackFile = nullptr;
  TTree* mTrackTree = nullptr;

  std::vector<NA6PTrackParCov> mInputTracks;
  float mMaxMatchChi2 = 25.f;

  ClassDefNV(NA6PTOFMatching, 1);
};

#endif
