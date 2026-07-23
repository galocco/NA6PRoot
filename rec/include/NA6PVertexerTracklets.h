// NA6PCCopyright

// Based on:
// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef NA6P_VERTEXER_TRACKLETS_H
#define NA6P_VERTEXER_TRACKLETS_H

#include <Rtypes.h>
#include <vector>
#include <array>
#include "NA6PLayoutParam.h"
#include "NA6PLine.h"
#include "NA6PRecoParam.h"
#include "NA6PMCTruthContainer.h"

class NA6PVertex;
class NA6PVerTelCluster;

// Vertex reconstruction from tracklets in the telescope

// structures for temporary objects used in track finding

struct TrackletForVertex {
  int startingLayer;
  int firstClusterIndex;
  int secondClusterIndex;
  float tanL;
  float phi;
  float pxpz;
  float pypz;
};

struct TracklIntersection {
  float zeta;
  float sigmazeta;
  float tanl;
  int firstClusterIndex;
  int secondClusterIndex;
};

struct ClusterLines {
  ClusterLines(int firstLabel, const NA6PLine& firstLine, int secondLabel, const NA6PLine& secondLine);
  void add(int lineLabel, const NA6PLine& line);
  void computeClusterCentroid();
  const std::vector<int>& getLabels() { return lineLabels; }
  int getSize() const { return lineLabels.size(); }
  const std::array<float, 3>& getVertex() const { return lineCluVertex; }
  const std::array<float, 6>& getRMS2() const { return lineCluRMS2; }
  inline float getAvgDistance2() const { return lineCluAvgDistance2; }

  bool operator==(const ClusterLines&) const;

  std::array<double, 6> lineCluAMatrix = {0.0};  // AX=B
  std::array<double, 3> lineCluBMatrix = {0.0};  // AX=B
  std::vector<int> lineLabels;                   // labels
  std::array<float, 9> lineWeightMatrix = {0.f}; // weight matrix
  std::array<float, 3> lineCluVertex = {0.f};    // cluster centroid position
  std::array<float, 6> lineCluRMS2 = {0.f};      // symmetric matrix: diagonal is RMS2
  float lineCluAvgDistance2 = 0.f;               // substitute for chi2
};

class NA6PVertexerTracklets
{
 public:
  enum { kYZ,
         kRZ,
         k3D,
         kXZ };
  enum { kNoWeight,
         kTanL,
         kSigma };
  enum { kMaxPileupVertices = 5,
         kMaxBinsForPeakFind = 1000 };
  enum { kKDE,
         kHistoPeak,
         kPairs };
  enum { kStandardKDE,
         kAdaptiveKDE };
  enum { kMultiVertOff,
         kMultiVertIterative,
         kAllVerticesInOneGo };

  NA6PVertexerTracklets();
  ~NA6PVertexerTracklets() = default;

  void configurePeakFinding()
  {
    mZBinWidth = (mRecoParam->vertexerZMax - mRecoParam->vertexerZMin) / mRecoParam->vertexerNBinsForPeakFind;
    mHistIntersec.assign(mRecoParam->vertexerNBinsForPeakFind, 0);
  }
  // Setters
  void setNLayersVT(int n) { mNLayersVT = n; }
  void setLayerToStart(int lay) { mLayerToStart = lay; }
  void setBeamX(float x) { mBeamX = x; }
  void setBeamY(float y) { mBeamY = y; }
  void setRecoInYZPlane() { mRecoType = kYZ; }
  void setRecoInXZPlane() { mRecoType = kXZ; }
  void setRecoInRZPlane() { mRecoType = kRZ; }
  void setRecoIn3D() { mRecoType = k3D; }
  void setUseHistoForPeakFinding() { mMethod = kHistoPeak; }
  void setUseKDEForPeakFinding() { mMethod = kKDE; }
  void setUseTrackletPairs() { mMethod = kPairs; }
  void setUseNoWeight() { mWeightedMeanOption = kNoWeight; }
  void setUseTanLWeight() { mWeightedMeanOption = kTanL; }
  void setUseSigmaWeight() { mWeightedMeanOption = kSigma; }

  void setUseStandardKDE() { mKDEOption = kStandardKDE; }
  void setUseAdaptiveKDE() { mKDEOption = kAdaptiveKDE; }
  void setMultiVertexOff() { mMultiVertexMode = kMultiVertOff; }
  void setMultiVertexInOneGo() { mMultiVertexMode = kAllVerticesInOneGo; }
  void setMultiVertexIterative() { mMultiVertexMode = kMultiVertIterative; }
  void setVerbosity(bool opt = true) { mVerbose = opt; }
  void setClusterMCTruth(NA6PMCTruthContainer* cont) { mCluMCLabels = cont; }

  void configureFromRecoParam();
  void printConfiguration() const;

  // methods for vertex selection
  void initTargets();
  bool isVertexInTarget(float zRecoVert, float tolerance = 0.1);

  // methods for vertex calculation
  void findVertices(std::vector<NA6PVerTelCluster>& cluArr,
                    std::vector<NA6PVertex>& vertices);

  void sortClustersByLayerAndEta(std::vector<NA6PVerTelCluster>& cluArr,
                                 std::vector<int>& firstIndex,
                                 std::vector<int>& lastIndex);
  void sortTrackletsByLayerAndIndex(std::vector<TrackletForVertex>& tracklets,
                                    std::vector<int>& firstIndex,
                                    std::vector<int>& lastIndex);
  void computeLayerTracklets(const std::vector<NA6PVerTelCluster>& cluArr,
                             const std::vector<int>& firstIndex,
                             const std::vector<int>& lastIndex,
                             std::vector<TrackletForVertex>& tracklets);
  void selectTracklets(const std::vector<TrackletForVertex>& tracklets,
                       const std::vector<int>& firstIndex,
                       const std::vector<int>& lastIndex,
                       std::vector<TrackletForVertex>& selTracklets);
  void filterOutUsedTracklets(const std::vector<TrackletForVertex>& tracklets,
                              std::vector<TrackletForVertex>& usableTracklets);
  void computeIntersections(const std::vector<TrackletForVertex>& selTracklets,
                            const std::vector<NA6PVerTelCluster>& cluArr,
                            std::vector<TracklIntersection>& zIntersec);
  bool findVertexHistoPeak(std::vector<TracklIntersection>& zIntersec,
                           std::vector<NA6PVertex>& vertices);

  float gaussKernel(float u)
  {
    return std::exp(-0.5 * u * u) / std::sqrt(2.0 * M_PI);
  }
  void resetClusters(int nClus)
  {
    mIsClusterUsed.resize(nClus);
    for (int jClu = 0; jClu < nClus; jClu++)
      mIsClusterUsed[jClu] = false;
  }
  bool findVertexKDE(const std::vector<TracklIntersection>& zIntersec,
                     std::vector<NA6PVertex>& vertices);
  bool compute3DVertices(const std::vector<TrackletForVertex>& selTracklets,
                         const std::vector<NA6PVerTelCluster>& cluArr,
                         std::vector<NA6PVertex>& vertices);
  void printStats(const std::vector<TrackletForVertex>& candidates,
                  const std::vector<NA6PVerTelCluster>& cluArr,
                  const std::string& label);
  short getMethodForPeakFinding() const { return mMethod; }
  short getMultiVertexMode() const { return mMultiVertexMode; }

 private:
  int mNLayersVT = 5;                                  // number of layers in the VT
  int mLayerToStart = 0;                               // innermost layer used in tracklets
  std::vector<bool> mIsClusterUsed = {};               // flag for used clusters
  float mBeamX = 0.;                                   // beam transverse coordindates
  float mBeamY = 0.;                                   // beam transverse coordindates
  short mRecoType = kYZ;                               // method to compute tracklet intersections with beam axis
  short mMethod = kKDE;                                // method for peak finding (KDE vs histo)
  int mWeightedMeanOption = kNoWeight;                 // option for weighted mean
  float mZBinWidth = 0.1;                              // bin width (recalculated from n bins)
  std::vector<int> mHistIntersec;                      // histogram for the peak finding method
  int mKDEOption = kStandardKDE;                       // option for using uncertainties in KDE sigma
  short mMultiVertexMode = kMultiVertIterative;        // method for multiple vertices
  bool mVerbose = false;                               // flag for verbosity
  int mNTargets = 0;                                   // number of targets
  float mZPosTarg[NA6PLayoutParam::MaxTargets] = {};   // target positions
  float mZThickTarg[NA6PLayoutParam::MaxTargets] = {}; // target thickesses
  const NA6PRecoParam* mRecoParam = nullptr;           // reconstruction parameters
  NA6PMCTruthContainer* mCluMCLabels = nullptr;        // MC truth for clusters

  ClassDefNV(NA6PVertexerTracklets, 1);
};

#endif
