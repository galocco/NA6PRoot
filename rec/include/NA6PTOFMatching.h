// NA6PCCopyright

#ifndef NA6P_TOF_MATCHING_H
#define NA6P_TOF_MATCHING_H

#include <string>
#include <vector>
#include <Rtypes.h>

#include "NA6PReconstruction.h"
#include "NA6PTOFCluster.h"
#include "NA6PTrack.h"

class TFile;
class TTree;
class NA6PTOFHit;

class NA6PTOFMatching : public NA6PReconstruction
{
 public:
  NA6PTOFMatching();
  NA6PTOFMatching(const char* geofile, const char* geoname = "NA6P");
  ~NA6PTOFMatching() override = default;

  bool initTOFMatching();

  // ── Hit → cluster conversion ────────────────────────────────────
  void setClusterSpaceResolutionX(double res) { mCluResX = res; }
  void setClusterSpaceResolutionY(double res) { mCluResY = res; }
  void hitsToRecPoints(const std::vector<NA6PTOFHit>& hits);

  // ── Cluster I/O ─────────────────────────────────────────────────
  void createClustersOutput() override;
  void clearClusters() override { mClusters.clear(); }
  void writeClusters() override;
  void closeClustersOutput() override;

  // ── Track I/O ───────────────────────────────────────────────────
  void createTracksOutput() override;
  void clearTracks() override { mTracks.clear(); }
  void writeTracks() override;
  void closeTracksOutput() override;

  // ── Matching ────────────────────────────────────────────────────
  void setVerTelTracks(std::vector<NA6PTrack>& tracks) { hVerTelTrackPtr = &tracks; }
  void setClusters(std::vector<NA6PTOFCluster>& clusters);
  void setMaxMatchDist(double d) { mMaxMatchDist = d; }
  void runTOFMatching();

  const std::vector<NA6PTOFCluster>& getClusters() const { return mClusters; }
  const std::vector<NA6PTrack>& getTracks() const { return mTracks; }

 private:
  // clusters
  std::vector<NA6PTOFCluster> mClusters, *hClusPtr = &mClusters;
  TFile* mClusFile = nullptr;
  TTree* mClusTree = nullptr;
  double mCluResX = 100.e-4/3.464101615; // cluster spatial resolution [cm]
  double mCluResY = 300.e-4/3.464101615; // cluster spatial resolution [cm]
  double mTimeRes = 20.e-12; // cluster time resolution: 20 ps [s]

  // output tracks (VT tracks with TOF assigned)
  std::vector<NA6PTrack> mTracks, *hTrackPtr = &mTracks;
  TFile* mTrackFile = nullptr;
  TTree* mTrackTree = nullptr;

  // input VT tracks (external pointer)
  std::vector<NA6PTrack>* hVerTelTrackPtr = nullptr;

  double mMaxMatchDist = 0.32; // max XY distance for TOF matching [cm]

  ClassDefNV(NA6PTOFMatching, 1);
};

#endif
