// NA6PCCopyright

// Minimal matching reconstruction driver

#ifndef NA6P_MATCHING_H
#define NA6P_MATCHING_H

#include <string>
#include <Rtypes.h>
#include <TVector3.h>

#include "NA6PReconstruction.h"
#include "NA6PVerTelCluster.h"
#include "NA6PMuonSpecCluster.h"
#include "NA6PFastTrackFitter.h"
#include "NA6PRecoParam.h"
#include "NA6PMCComposedLabel.h"
#include "NA6PTreeStreamRedirector.h"

#include "NA6PMatch.h"
#include "NA6PTrack.h"
#include "NA6PVertex.h"

class TFile;
class TTree;
class NA6PFastTrackFitter;

// #define _CHI2_TUNING_MODE_

class NA6PMatching : public NA6PReconstruction
{
 public:
  struct MatchRecord {
    int vtID = -1; // reference of the VT track
    int msID = -1; // reference of the MS track
    float chi2Match = 1e9;
  };

  NA6PMatching(bool init = true);
  ~NA6PMatching() override;

  void setMCMatching(bool mc) { mMCMatching = mc; }

  void setVerTelTracks(std::vector<NA6PTrack>& tracks);
  void setMuonSpecTracks(std::vector<NA6PTrack>& tracks);
  void setVerTelTrackMCLabels(std::vector<NA6PMCComposedLabel>& trLab);
  void setMuonSpecTrackMCLabels(std::vector<NA6PMCComposedLabel>& trLab);

  void setVerTelClusters(std::vector<NA6PVerTelCluster>& clusters);
  void setMuonSpecClusters(std::vector<NA6PMuonSpecCluster>& clusters);

  bool initMatching();
  void createTracksOutput() override;
  void clearTracks() override
  {
    mMatchedTracks.clear();
    mMatchedTrkMCLabels.clear();
  }
  void writeTracks() override;
  void closeTracksOutput() override;
  void propToZMatching(std::vector<NA6PTrack>& tracks, float z, bool outer = false);
  void runMatching();

  bool fitAndStoreMatchedTrack(int vtIdx, int msIdx, float matchChi2);

 private:
  void addClustersToFitter(const NA6PTrack& trk, const auto* clusPtr);
  void runMCMatching();
  void runDataMatching();
  void prefilterTracks();
  void buildMatchingCandidates(int msID);

  std::unordered_map<NA6PMCComposedLabel, int> buildMCMatchingIndex();

  const NA6PRecoParam* mRecoParam = nullptr;

  std::vector<NA6PVerTelCluster>* hVerTelClusPtr = nullptr;
  std::vector<NA6PMuonSpecCluster>* hMuonSpecClusPtr = nullptr;

  int mNDOF = 5;
  bool mMCMatching = false;

  std::unique_ptr<NA6PFastTrackFitter> mTrackFitter;
  std::vector<NA6PTrack>* hVerTelTrackPtr = nullptr;                   // Vertex telescope tracks
  std::vector<NA6PTrack>* hMuonSpecTrackPtr = nullptr;                 // Muon spectrometer tracks
  std::vector<NA6PMCComposedLabel>* hVerTelTrkMCLabelsPtr = nullptr;   // vertel track MC labels
  std::vector<NA6PMCComposedLabel>* hMuonSpecTrkMCLabelsPtr = nullptr; // vertel track MC labels

  std::vector<int> mSelIDVT;
  std::vector<int> mSelIDMS;
  std::vector<NA6PMatch> mMatchedTracks, *hMatchedTrackPtr = &mMatchedTracks;                           // Matched tracks
  std::vector<NA6PMCComposedLabel> mMatchedTrkMCLabels, *hMatchedTrkMCLabelsPtr = &mMatchedTrkMCLabels; // track MC labels
  std::vector<MatchRecord> mMatchRecords;

  TFile* mMatchedTrackFile = nullptr; // file with Matched tracks
  TTree* mMatchedTrackTree = nullptr; // tree of Matched tracks

  std::unique_ptr<NA6PTreeStreamRedirector> dbgStream;

  ClassDefNV(NA6PMatching, 1);
};

#endif
