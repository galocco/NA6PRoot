// NA6PCCopyright

#include "NA6PTOFMatching.h"
#include "NA6PTOFHit.h"

#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <fairlogger/Logger.h>
#include <cmath>
#include <limits>

ClassImp(NA6PTOFMatching)

// ============================================================
//  Constructors
// ============================================================

NA6PTOFMatching::NA6PTOFMatching() : NA6PReconstruction("TOFMatching")
{
}

NA6PTOFMatching::NA6PTOFMatching(const char* geofile,
                                 const char* geoname)
  : NA6PReconstruction("TOFMatching")
{
  mGeoFilName = geofile;
  mGeoObjName = geoname;
  initTOFMatching();
}

bool NA6PTOFMatching::initTOFMatching()
{
  NA6PReconstruction::init(mGeoFilName.c_str(), mGeoObjName.c_str());
  createTracksOutput();
  return true;
}

// ============================================================
//  Hit → Cluster conversion
// ============================================================

void NA6PTOFMatching::hitsToRecPoints(const std::vector<NA6PTOFHit>& hits)
{
  int nHits = hits.size();
  for (int jHit = 0; jHit < nHits; ++jHit) {
    const auto& hit = hits[jHit];
    double x = hit.getX();
    double y = hit.getY();
    double z = hit.getZ();
    double ex2 = 5.e-4;
    double ey2 = 5.e-4;
    if (mCluResX > 0 && mCluResY > 0) {
      x = gRandom->Gaus(hit.getX(), mCluResX);
      y = gRandom->Gaus(hit.getY(), mCluResY);
      ex2 = mCluResX * mCluResX;
      ey2 = mCluResY * mCluResY;
    }
    double eloss = hit.getHitValue();
    int clusiz = 2;
    if (eloss > 2.e-5 && eloss < 5.e-5)
      clusiz = 3;
    else if (eloss > 5.e-5)
      clusiz = 4;

    int nDet = hit.getDetectorID();
    int idPart = hit.getTrackID();
    int layer = nDet / 4; // single plane → layer 0
    float time = hit.getTime();
    if (mTimeRes > 0.f)
      time = gRandom->Gaus(time, mTimeRes);
    

    mClusters.emplace_back(x, y, z, clusiz, layer, time);
    auto& clu = mClusters.back();
    clu.setErr(ex2, 0., ey2);
    clu.setDetectorID(nDet);
    clu.setParticleID(idPart);
    clu.setHitID(jHit);
  }
}

// ============================================================
//  Cluster I/O
// ============================================================

void NA6PTOFMatching::createClustersOutput()
{
  auto nm = fmt::format("ClustersTOF.root");
  mClusFile = TFile::Open(nm.c_str(), "recreate");
  mClusTree = new TTree("clustersTOF", "TOF Clusters");
  mClusTree->Branch("TOF", &hClusPtr);
  LOGP(info, "Will store TOF clusters in {}", nm);
}

void NA6PTOFMatching::writeClusters()
{
  if (mClusTree) {
    mClusTree->Fill();
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

// ============================================================
//  Track I/O
// ============================================================

void NA6PTOFMatching::createTracksOutput()
{
  auto nm = fmt::format("TracksTOFMatching.root");
  mTrackFile = TFile::Open(nm.c_str(), "recreate");
  mTrackTree = new TTree("tracksTOFMatching", "VT Tracks with TOF");
  mTrackTree->Branch("TOFMatching", &hTrackPtr);
  LOGP(info, "Will store TOF-matched tracks in {}", nm);
}

void NA6PTOFMatching::writeTracks()
{
  if (mTrackTree) {
    mTrackTree->Fill();
    LOGP(info, "Saved {} TOF-matched tracks in tree with {} entries",
         mTracks.size(), mTrackTree->GetEntries());
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

// ============================================================
//  setClusters
// ============================================================

void NA6PTOFMatching::setClusters(std::vector<NA6PTOFCluster>& clusters)
{
  hClusPtr = &clusters;
  for (size_t i = 0; i < clusters.size(); ++i) {
    clusters[i].setClusterIndex(static_cast<int>(i));
  }
}

// ============================================================
//  TOF matching
//
//  For each VT track:
//    1. Copy the outer param and propagate to the TOF plane z.
//    2. Find the closest TOF cluster in (x,y).
//    3. If distance < mMaxMatchDist, assign the TOF time to the track.
//    4. Store the track in the output.
// ============================================================

void NA6PTOFMatching::runTOFMatching()
{
  if (!hVerTelTrackPtr) {
    LOGP(error, "No VT tracks set for TOF matching");
    return;
  }
  if (!hClusPtr || hClusPtr->empty()) {
    LOGP(warn, "No TOF clusters for this event");
    // still save all tracks without TOF
    for (auto& trk : *hVerTelTrackPtr) {
      mTracks.push_back(trk);
    }
    writeTracks();
    return;
  }

  // Determine the TOF plane z from the first cluster
  const double zTOF = (*hClusPtr)[0].getZLab();

  int nMatched = 0;
  for (auto& vtTrack : *hVerTelTrackPtr) {
    NA6PTrack trk = vtTrack; // work on a copy

    // Propagate the outer param to the TOF z
    if (!trk.propagateToZBxByBz(zTOF, 1.0, 0., 0., true)) {
      LOGP(debug, "Failed to propagate VT track to TOF z={:.2f}", zTOF);
      mTracks.push_back(vtTrack); // keep original, no TOF
      continue;
    }

    // Track position at TOF in lab frame (outer param)
    const double trkX = trk.getOuterParam().getY(); // Y_TF = X_lab
    const double trkY = trk.getOuterParam().getZ(); // Z_TF = Y_lab

    // Find closest TOF cluster
    double bestDist2 = std::numeric_limits<double>::max();
    int bestIdx = -1;
    for (size_t ic = 0; ic < hClusPtr->size(); ++ic) {
      const auto& clu = (*hClusPtr)[ic];
      double dx = clu.getXLab() - trkX;
      double dy = clu.getYLab() - trkY;
      double d2 = dx * dx + dy * dy;
      if (d2 < bestDist2) {
        bestDist2 = d2;
        bestIdx = static_cast<int>(ic);
      }
    }

    double bestDist = std::sqrt(bestDist2);
    if (bestIdx >= 0 && bestDist < mMaxMatchDist) {
      vtTrack.setTOF((*hClusPtr)[bestIdx].getTime());
      ++nMatched;
    }
    mTracks.push_back(vtTrack);
  }

  LOGP(info, "TOF matching: {} / {} VT tracks matched (maxDist={:.2f} cm)",
       nMatched, hVerTelTrackPtr->size(), mMaxMatchDist);
  writeTracks();
}
