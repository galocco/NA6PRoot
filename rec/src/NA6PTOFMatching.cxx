// NA6PCCopyright

#include "NA6PTOFMatching.h"

#include "NA6PTOFHit.h"
#include "NA6PRecoParam.h"
#include "NA6PVertex.h"
#include "Propagator.h"

#include <TFile.h>
#include <TRandom.h>
#include <TTree.h>
#include <fairlogger/Logger.h>
#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
bool getHelixLengthToZ(const NA6PTrackPar& track, float z, float by, float& length)
{
  const float dz = z - track.getZ();
  if (std::abs(dz) < 1.e-6f) {
    length = 0.f;
    return true;
  }

  const float s0 = track.getTx();
  const float c0 = track.getCosPsi();
  const float kappa = track.getCurvature(by);
  const float pitch = track.getP2Pxz();
  if (c0 < NA6PTrackPar::kTinyF) {
    return false;
  }
  if (std::abs(kappa) < NA6PTrackPar::kSmallKappa) {
    length = std::abs(dz) * pitch / c0;
    return true;
  }

  const float s1 = s0 + kappa * dz;
  if (std::abs(s1) >= NA6PTrackPar::kAlmost1F) {
    return false;
  }
  const float deltaPsi = std::asin(s1) - std::asin(s0);
  length = std::abs(deltaPsi / kappa) * pitch;
  return std::isfinite(length);
}

bool propagateToZWithHelixLength(NA6PTrackParCov& track, float z, float& pathLength)
{
  constexpr float MaxStep = 2.f;
  constexpr float ZTolerance = 1.e-5f;
  auto* propagator = Propagator::Instance();
  pathLength = 0.f;

  while (std::abs(z - track.getZ()) > ZTolerance) {
    const float dz = z - track.getZ();
    const float zNext = track.getZ() + std::copysign(std::min(std::abs(dz), MaxStep), dz);
    const float by = propagator->getBy(track.getXYZ());
    float stepLength = 0.f;
    if (!getHelixLengthToZ(track, zNext, by, stepLength) ||
        !propagator->propagateToZ(track, zNext)) {
      return false;
    }
    pathLength += stepLength;
  }
  return true;
}
} // namespace

NA6PTOFMatching::NA6PTOFMatching(bool init) : NA6PReconstruction("TOFMatching")
{
  if (init) {
    initTOFMatching();
  }
}

bool NA6PTOFMatching::initTOFMatching()
{
  configureFromRecoParam();
  createTracksOutput();
  return true;
}

void NA6PTOFMatching::configureFromRecoParam()
{
  const auto& recoParam = NA6PRecoParam::Instance();
  mCluResX = recoParam.tofClusterResolutionX;
  mCluResY = recoParam.tofClusterResolutionY;
  mTimeRes = recoParam.tofTimeResolution;
  mMaxMatchChi2 = recoParam.tofMaxMatchChi2;
  LOGP(info, "TOF configuration: sigmaX={} cm, sigmaY={} cm, sigmaT={} s, maxChi2={}",
       mCluResX, mCluResY, mTimeRes, mMaxMatchChi2);
}

void NA6PTOFMatching::hitsToRecPoints(const std::vector<NA6PTOFHit>& hits)
{
  for (size_t hitID = 0; hitID < hits.size(); ++hitID) {
    const auto& hit = hits[hitID];
    float x = hit.getX();
    float y = hit.getY();
    const float z = hit.getZ();
    float errX2 = 5.e-4f;
    float errY2 = 5.e-4f;
    if (mCluResX > 0.f && mCluResY > 0.f) {
      x = gRandom->Gaus(x, mCluResX);
      y = gRandom->Gaus(y, mCluResY);
      errX2 = mCluResX * mCluResX;
      errY2 = mCluResY * mCluResY;
    }

    int clusterSize = 1;
    
    float time = hit.getTime();
    if (mTimeRes > 0.f) {
      time = gRandom->Gaus(time, mTimeRes);
    }

    const int clusterID = static_cast<int>(mClusters.size());
    mClusters.emplace_back(x, y, z, clusterSize, 0, time);
    auto& cluster = mClusters.back();
    cluster.setErr(errX2, 0.f, errY2);
    cluster.setDetectorID(hit.getDetectorID());
    cluster.setHitID(static_cast<int>(hitID));
    cluster.setClusterIndex(clusterID);
  }
}

void NA6PTOFMatching::createClustersOutput()
{
  const auto name = fmt::format("ClustersTOF.root");
  mClusFile = TFile::Open(name.c_str(), "recreate");
  mClusTree = new TTree("clustersTOF", "TOF Clusters");
  mClusTree->Branch("TOF", &hClusPtr);
  LOGP(info, "Will store TOF clusters in {}", name);
}

void NA6PTOFMatching::writeClusters()
{
  if (mClusTree) {
    mClusTree->Fill();
    LOGP(info, "Saved {} clusters in tree with {} entries", mClusters.size(), mClusTree->GetEntries());
  }
}

void NA6PTOFMatching::closeClustersOutput()
{
  if (mClusTree && mClusFile) {
    mClusFile->cd();
    mClusTree->Write();
    delete mClusTree;
    mClusTree = nullptr;
    mClusFile->Close();
    delete mClusFile;
    mClusFile = nullptr;
  }
}

void NA6PTOFMatching::createTracksOutput()
{
  const auto name = fmt::format("TracksTOFMatching.root");
  mTrackFile = TFile::Open(name.c_str(), "recreate");
  mTrackTree = new TTree("tracksTOFMatching", "TOF-matched VT Tracks");
  mTrackTree->Branch("TOFMatching", &hMatchedTrackPtr);
  LOGP(info, "Will store TOF-matched tracks in {}", name);
}

void NA6PTOFMatching::writeTracks()
{
  if (mTrackTree) {
    mTrackTree->Fill();
    LOGP(info, "Saved {} TOF-matched tracks in tree with {} entries",
         mMatchedTracks.size(), mTrackTree->GetEntries());
  }
}

void NA6PTOFMatching::closeTracksOutput()
{
  if (mTrackTree && mTrackFile) {
    mTrackFile->cd();
    mTrackTree->Write();
    delete mTrackTree;
    mTrackTree = nullptr;
    mTrackFile->Close();
    delete mTrackFile;
    mTrackFile = nullptr;
  }
}

void NA6PTOFMatching::setClusters(std::vector<NA6PTOFCluster>& clusters)
{
  hClusPtr = &clusters;
  for (size_t i = 0; i < clusters.size(); ++i) {
    clusters[i].setClusterIndex(static_cast<int>(i));
  }
}

void NA6PTOFMatching::setTracks(const std::vector<NA6PTrackParCov>& tracks)
{
  mInputTracks = tracks;
}

void NA6PTOFMatching::setMatchedTracks(const std::vector<NA6PMatch>& tracks)
{
  mInputTracks.clear();
  mInputTracks.reserve(tracks.size());
  for (const auto& track : tracks) {
    mInputTracks.emplace_back(static_cast<const NA6PTrackParCov&>(track));
  }
}

void NA6PTOFMatching::setVerTelTracks(const std::vector<NA6PTrack>& tracks)
{
  mInputTracks.clear();
  mInputTracks.reserve(tracks.size());
  for (const auto& track : tracks) {
    mInputTracks.emplace_back(track.getOuterParam());
  }
}

void NA6PTOFMatching::runTOFMatching()
{
  clearTracks();
  if (mInputTracks.empty()) {
    LOGP(warn, "No input tracks for TOF matching");
    writeTracks();
    return;
  }
  if (!hClusPtr || hClusPtr->empty()) {
    LOGP(warn, "No TOF clusters for this event");
    writeTracks();
    return;
  }
  if (!mPrimaryVertex) {
    LOGP(warn, "No primary vertex for TOF matching");
    writeTracks();
    return;
  }

  LOGP(info, "Process event with nVTTracks {} nTOFClusters {}, primary vertex in z = {} cm",
       mInputTracks.size(), hClusPtr->size(), mPrimaryVertex->getZ());

  const float zTOF = hClusPtr->front().getZ();
  int nPropagationFailures = 0;
  int nVertexPropagationFailures = 0;

  for (size_t trackID = 0; trackID < mInputTracks.size(); ++trackID) {
    NA6PTrackParCov trackAtTOF = mInputTracks[trackID];
    if (!trackAtTOF.isValid() ||
        !Propagator::Instance()->propagateToZ(trackAtTOF, zTOF)) {
      ++nPropagationFailures;
      continue;
    }

    float bestChi2 = std::numeric_limits<float>::max();
    int bestClusterID = -1;
    for (size_t clusterID = 0; clusterID < hClusPtr->size(); ++clusterID) {
      const auto& cluster = (*hClusPtr)[clusterID];
      const float chi2 = trackAtTOF.getPredictedChi2(cluster);
      if (chi2 < bestChi2) {
        bestChi2 = chi2;
        bestClusterID = static_cast<int>(clusterID);
      }
    }

    if (bestClusterID >= 0 && bestChi2 < mMaxMatchChi2) {
      const auto& cluster = (*hClusPtr)[bestClusterID];
      NA6PTrackParCov trackAtVertex = trackAtTOF;
      float pathLength = 0.f;
      if (!propagateToZWithHelixLength(trackAtVertex, mPrimaryVertex->getZ(), pathLength)) {
        ++nVertexPropagationFailures;
        continue;
      }
      auto& match = mMatchedTracks.emplace_back(trackAtVertex);
      match.setTOF(cluster.getTime());
      match.setPathLength(pathLength);
      match.setMatchChi2(bestChi2);
      match.setIndexTrack(static_cast<int>(trackID));
      match.setIndexTOF(bestClusterID);
    }
  }

  LOGP(info,
       "TOF matching: {} / {} tracks matched (maxChi2={:.2f}, {} TOF and {} vertex propagation failures)",
       mMatchedTracks.size(), mInputTracks.size(), mMaxMatchChi2,
       nPropagationFailures, nVertexPropagationFailures);
  writeTracks();
}
