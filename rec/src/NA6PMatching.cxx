// NA6PCCopyright

#include "NA6PMatching.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <unordered_map>
#include <fairlogger/Logger.h>

#include <TMatrixD.h>
#include <TVectorD.h>

NA6PMatching::NA6PMatching(bool init) : NA6PReconstruction("Matching")
{
#ifdef _CHI2_TUNING_MODE_
  dbgStream = std::make_unique<NA6PTreeStreamRedirector>("matching_dbg.root");
#endif
  if (init) {
    initMatching();
  }
}

NA6PMatching::~NA6PMatching()
{
#ifdef _CHI2_TUNING_MODE_
  if (dbgStream) {
    dbgStream->Close();
  }
#endif
}

bool NA6PMatching::initMatching()
{
  mRecoParam = &NA6PRecoParam::Instance();
  mTrackFitter = std::make_unique<NA6PFastTrackFitter>();
  mTrackFitter->enableMaterialCorrections();
  mTrackFitter->setPropagateToPrimaryVertex(false);
  mTrackFitter->setMaxPropagationStep(mRecoParam->maxPropagationStep);

  createTracksOutput();
  mNDOF = Propagator::Instance()->getNDOFTrack();
  return true;
}

void NA6PMatching::createTracksOutput()
{
  auto nm = fmt::format("Tracks{}.root", getName());
  mMatchedTrackFile = TFile::Open(nm.c_str(), "recreate");
  mMatchedTrackTree = new TTree(fmt::format("tracks{}", getName()).c_str(), fmt::format("{} Tracks", getName()).c_str());
  mMatchedTrackTree->Branch(getName().c_str(), &hMatchedTrackPtr);
  if (mReadMCTruth)
    mMatchedTrackTree->Branch(fmt::format("{}MCTruth", getName()).c_str(), &hMatchedTrkMCLabelsPtr);
  LOGP(info, "Will store {} tracks in {}", getName(), nm);
}

void NA6PMatching::writeTracks()
{
  if (mMatchedTrackTree) {
    mMatchedTrackTree->Fill();
    LOGP(info, "Saved {} tracks in tree with {} entries", mMatchedTracks.size(), mMatchedTrackTree->GetEntries());
  }
}

void NA6PMatching::closeTracksOutput()
{
  if (mMatchedTrackTree && mMatchedTrackFile) {
    mMatchedTrackFile->cd();
    mMatchedTrackTree->Write();
    delete mMatchedTrackTree;
    mMatchedTrackTree = nullptr;
    mMatchedTrackFile->Close();
    delete mMatchedTrackFile;
    mMatchedTrackFile = nullptr;
  }
}

void NA6PMatching::setVerTelTracks(std::vector<NA6PTrack>& tracks)
{
  hVerTelTrackPtr = &tracks;
}

void NA6PMatching::setMuonSpecTracks(std::vector<NA6PTrack>& tracks)
{
  hMuonSpecTrackPtr = &tracks;
}
void NA6PMatching::setVerTelTrackMCLabels(std::vector<NA6PMCComposedLabel>& trLab)
{
  hVerTelTrkMCLabelsPtr = &trLab;
}

void NA6PMatching::setMuonSpecTrackMCLabels(std::vector<NA6PMCComposedLabel>& trLab)
{
  hMuonSpecTrkMCLabelsPtr = &trLab;
}

void NA6PMatching::setMuonSpecClusters(std::vector<NA6PMuonSpecCluster>& clusters)
{
  hMuonSpecClusPtr = &clusters;
  for (size_t i = 0; i < clusters.size(); ++i) {
    clusters[i].setClusterIndex(static_cast<int>(i));
  }
}

void NA6PMatching::setVerTelClusters(std::vector<NA6PVerTelCluster>& clusters)
{
  hVerTelClusPtr = &clusters;
  for (size_t i = 0; i < clusters.size(); ++i) {
    clusters[i].setClusterIndex(static_cast<int>(i));
  }
}

void NA6PMatching::propToZMatching(std::vector<NA6PTrack>& tracks, float z, bool outer)
{
  // RSTODO consider propagation track copies. Or better, define outer ITS and inward MS tracks at the reference Z
  for (auto& track : tracks) {
    NA6PTrackParCov& trc = outer ? track.getOuterParam() : track;
    if (trc.isValid() && std::abs(trc.getZ() - z) > 1e-4) {
      if (!Propagator::Instance()->propagateToZ(trc, z, mTrackFitter->getPropOpt())) {
        trc.invalidate();
      }
    }
  }
};

void NA6PMatching::addClustersToFitter(const NA6PTrack& trk, const auto* clusPtr)
{
  uint32_t cmap = trk.getClusterMap();
  for (int lr = 0; lr < NA6PTrack::kMaxLr; ++lr) {
    if (cmap & (1u << lr)) {
      mTrackFitter->addCluster((*clusPtr)[trk.getClusterIndex(lr)]);
    }
  }
}

void NA6PMatching::buildMatchingCandidates(int iMS)
{
  // find at most mMaxMatchesPerTrack best matches for given MS track
  const auto& msTrack = (*hMuonSpecTrackPtr)[mSelIDMS[iMS]];
  const float curv = msTrack.getQ2Pxz();
  const float err = std::sqrt(msTrack.getSigmaQ2Pxz2() + mRecoParam->mtPrelDiffMargin2[4]) * mRecoParam->mtPrelCurveNSigma;
  auto lower = std::partition_point(mSelIDVT.begin(), mSelIDVT.end(), [&](size_t id) { return (*hVerTelTrackPtr)[id].getOuterParam().getQ2Pxz() < curv - err; });
  auto upper = std::partition_point(lower, mSelIDVT.end(), [&](size_t id) { return (*hVerTelTrackPtr)[id].getOuterParam().getQ2Pxz() < curv + err; });

  float chi2Worst = -999.f;
  int idWorst = -1, nCand = 0;
  for (auto it = lower; it != upper; ++it) {
    const auto& vtTrackOut = (*hVerTelTrackPtr)[*it].getOuterParam();
    // crude check in diagonals only
    float chi2Crude = 0.f, chi2Match = 0.f;
    std::array<float, 5> nsig2 = {};
    for (int i = 0; i < mNDOF; i++) {
      auto dif = msTrack.getParam(i) - vtTrackOut.getParam(i);
      nsig2[i] = dif * dif / (msTrack.getCovMatElem(i, i) + vtTrackOut.getCovMatElem(i, i) + mRecoParam->mtPrelDiffMargin2[i]);
      if (nsig2[i] > mRecoParam->mtNSigma2Crude[i]) {
        chi2Crude = mRecoParam->mtMaxChi2Crude + 1.f;
        break;
      }
      chi2Crude += nsig2[i];
    }
    if (chi2Crude > mRecoParam->mtMaxChi2Crude) { // full chi2 calculation
      continue;
    }
    chi2Match = msTrack.getPredictedChi2(vtTrackOut);
    if (chi2Match < mRecoParam->mtMaxChi2Match) {
      if (chi2Match > chi2Worst) {
        chi2Worst = chi2Match;
        idWorst = mMatchRecords.size();
      }
      if (nCand < mRecoParam->mtMaxMatchesPerTrack) {
        mMatchRecords.emplace_back(MatchRecord{.vtID = int(std::distance(mSelIDVT.begin(), it)), .msID = iMS, .chi2Match = chi2Match});
        nCand++;
      } else if (idWorst >= (int)mMatchRecords.size()) { // substitute the old worst match
        mMatchRecords[idWorst] = MatchRecord{.vtID = int(std::distance(mSelIDVT.begin(), it)), .msID = iMS, .chi2Match = chi2Match};
      }
    }
#ifdef _CHI2_TUNING_MODE_
    std::vector<float> nsigcrude2(5);
    NA6PMCComposedLabel lblMS = (*hMuonSpecTrkMCLabelsPtr)[mSelIDMS[iMS]];
    for (int i = 0; i < 5; i++) {
      nsigcrude2[i] = nsig2[i];
    }
    (*dbgStream) << "match"
                 << "vtTrack=" << vtTrackOut << "msTrack=" << ((NA6PTrackParCov&)msTrack)
                 << "nsig2crude=" << nsigcrude2 << "chi2Match=" << chi2Match << "chi2Crude=" << chi2Crude
                 << "vtPartID=" << (*hVerTelTrkMCLabelsPtr)[*it] << "msPartID=" << lblMS
                 << "\n";
#endif
  }
}

std::unordered_map<NA6PMCComposedLabel, int> NA6PMatching::buildMCMatchingIndex()
{
  std::unordered_map<NA6PMCComposedLabel, int> vtByPid;
  vtByPid.reserve(hVerTelTrackPtr->size());

  for (size_t i = 0; i < hVerTelTrackPtr->size(); ++i) {
    const auto& vt = (*hVerTelTrackPtr)[i];
    NA6PMCComposedLabel lab = (*hVerTelTrkMCLabelsPtr)[i];
    if (lab.isFake() || vt.getNHits() < mRecoParam->mtMinVTHits) {
      continue;
    }
    // Prefer track with more hits if there are duplicates with same particle ID
    auto it = vtByPid.find(lab);
    if (it == vtByPid.end() || vt.getNHits() > (*hVerTelTrackPtr)[it->second].getNHits()) {
      vtByPid[lab] = static_cast<int>(i);
    }
  }
  return vtByPid;
}

void NA6PMatching::prefilterTracks()
{
  // the optimal is to work with measured parameters whose errors are normal
  // apply selections
  for (int i = 0; i < (int)hVerTelTrackPtr->size(); ++i) {
    const auto& track = (*hVerTelTrackPtr)[i];
    if (track.getOuterParam().isValid() && track.getNHits() >= mRecoParam->mtMinVTHits) {
      mSelIDVT.push_back(i);
    }
  }
  for (int i = 0; i < (int)hMuonSpecTrackPtr->size(); ++i) {
    const auto& track = (*hMuonSpecTrackPtr)[i];
    if (track.getNHits() >= mRecoParam->mtMinMSHits) {
      mSelIDMS.push_back(i);
    }
  }
  // sort in increasing curvature
  std::sort(mSelIDVT.begin(), mSelIDVT.end(), [&](const int i, const int j) { return (*hVerTelTrackPtr)[i].getOuterParam().getQ2Pxz() < (*hVerTelTrackPtr)[j].getOuterParam().getQ2Pxz(); });
  std::sort(mSelIDMS.begin(), mSelIDMS.end(), [&](const int i, const int j) { return (*hMuonSpecTrackPtr)[i].getQ2Pxz() < (*hMuonSpecTrackPtr)[j].getQ2Pxz(); });
}

bool NA6PMatching::fitAndStoreMatchedTrack(int vtIdx, int msIdx, float chi2Match)
{
  const NA6PTrack& vtTrk = (*hVerTelTrackPtr)[vtIdx];
  const NA6PTrack& msTrk = (*hMuonSpecTrackPtr)[msIdx];
  mTrackFitter->cleanupAndStartFit();
  addClustersToFitter(vtTrk, hVerTelClusPtr);
  addClustersToFitter(msTrk, hMuonSpecClusPtr);

  auto& matchedTrack = mMatchedTracks.emplace_back();
  mTrackFitter->setMaxChi2Cl(mRecoParam->mtMaxChi2TrCl); // RSTODO make this general init
  int ndf = std::max(1, mTrackFitter->getNumberOfClusters() * 2 - mNDOF);
  float chi2Out = -1.f, chi2Inw = -1.f;
  if (mRecoParam->mtDoOutwardInwardFit) {
    ((NA6PTrackParCov&)matchedTrack) = vtTrk;
    matchedTrack.setPID(PID::Muon);
    chi2Out = mTrackFitter->fitSeedOutward(matchedTrack, true, mRecoParam->mtUseLinRefOut);
    if (chi2Out < 0 || chi2Out >= mRecoParam->mtMaxChi2NormRefit * ndf) {
      LOGP(warn, "Forward fit of matched track failed, skipping (MCTruth: VTID:{} MSID:{})", (*hVerTelTrkMCLabelsPtr)[vtIdx].asString(), (*hMuonSpecTrkMCLabelsPtr)[msIdx].asString());
      mMatchedTracks.pop_back();
      return false;
    }
  } else {
    ((NA6PTrackParCov&)matchedTrack) = msTrk.getOuterParam();
    matchedTrack.setPID(PID::Muon);
  }
  chi2Inw = mTrackFitter->fitSeedInward(matchedTrack, true, mRecoParam->mtUseLinRefInw);
  if (chi2Inw < 0 || chi2Inw >= mRecoParam->mtMaxChi2NormRefit * ndf) {
    LOGP(warn, "Inward fit of matched track failed, skipping (MCTruth: VTID:{} MSID:{})", (*hVerTelTrkMCLabelsPtr)[vtIdx].asString(), (*hMuonSpecTrkMCLabelsPtr)[msIdx].asString());
    mMatchedTracks.pop_back();
    return false;
  }
  if (mRecoParam->mtPropagateMatchedTracksToPV && mPrimaryVertex) {
    Propagator::Instance()->propagateToZ(matchedTrack, mPrimaryVertex->getZ(), mTrackFitter->getPropOpt());
  }
  matchedTrack.setChi2Match(chi2Match);
  matchedTrack.setChi2Refit(chi2Inw);
  matchedTrack.setNClusters(mTrackFitter->getNumberOfClusters());
  matchedTrack.setIndexVT(vtIdx);
  matchedTrack.setIndexMS(msIdx);
  if (mReadMCTruth) {
    NA6PMCComposedLabel lblVT = (*hVerTelTrkMCLabelsPtr)[vtIdx];
    NA6PMCComposedLabel lblMS = (*hMuonSpecTrkMCLabelsPtr)[msIdx];
    NA6PMCComposedLabel lblGlo = lblMS;
    if (lblVT != lblMS)
      lblGlo.setFakeFlag();
    mMatchedTrkMCLabels.push_back(lblGlo);
  }

#ifdef _CHI2_TUNING_MODE_
  (*dbgStream) << "refit"
               << "vtTrack=" << ((NA6PTrackParCov&)vtTrk) << "msTrack=" << ((NA6PTrackParCov&)msTrk)
               << "chi2vec=" << mTrackFitter->getChi2Buffer() << "chi2Match=" << chi2Match
               << "chi2Out=" << chi2Out << "chi2Inw=" << chi2Inw
               << "vtPartID=" << (*hVerTelTrkMCLabelsPtr)[vtIdx]() << "msPartID=" << (*hMuonSpecTrkMCLabelsPtr)[msIdx]()
               << "\n";
#endif
  return true;
}

void NA6PMatching::runMCMatching()
{
  auto vtByPid = buildMCMatchingIndex();
  for (size_t msIdx = 0; msIdx < hMuonSpecTrackPtr->size(); ++msIdx) {
    const auto& msTrack = (*hMuonSpecTrackPtr)[msIdx];
    if (msTrack.getNHits() < mRecoParam->mtMinMSHits)
      continue;
    NA6PMCComposedLabel lblMS = (*hMuonSpecTrkMCLabelsPtr)[msIdx];
    if (lblMS.isFake())
      continue;

    auto it = vtByPid.find(lblMS);
    if (it == vtByPid.end())
      continue;

    int vtIdx = it->second;
    fitAndStoreMatchedTrack(vtIdx, msIdx, 0.0);
  }
}

void NA6PMatching::runDataMatching()
{
  mSelIDVT.clear();
  mSelIDMS.clear();
  mMatchRecords.clear();

  propToZMatching(*hVerTelTrackPtr, mRecoParam->mtZMatching, true);
  propToZMatching(*hMuonSpecTrackPtr, mRecoParam->mtZMatching, false);
  prefilterTracks();

  for (int iMS = 0; iMS < (int)mSelIDMS.size(); iMS++) { // build matching candidates for each muon
    buildMatchingCandidates(iMS);
  }
  std::sort(mMatchRecords.begin(), mMatchRecords.end(), [](const auto& a, const auto& b) { return a.chi2Match < b.chi2Match; });
  // now we have select the best mutual matches, discarding already validated tracks
  for (const auto& mtc : mMatchRecords) {
    if (mSelIDVT[mtc.vtID] < 0 || mSelIDMS[mtc.msID] < 0) { // one of these tracke was already used for a match
      continue;
    }
    if (fitAndStoreMatchedTrack(mSelIDVT[mtc.vtID], mSelIDMS[mtc.msID], mtc.chi2Match)) { // register if good fit and flag tracks as used
      mSelIDVT[mtc.vtID] = -1;
      mSelIDMS[mtc.msID] = -1;
    }
  }
}

void NA6PMatching::runMatching()
{
  clearTracks();
  mMatchedTracks.reserve(hMuonSpecTrackPtr ? hMuonSpecTrackPtr->size() : 0);
  LOGP(info, "Process event with nVTTracks {} nMSTracks {}, primary vertex in z = {} cm",
       hVerTelTrackPtr ? hVerTelTrackPtr->size() : 0,
       hMuonSpecTrackPtr ? hMuonSpecTrackPtr->size() : 0,
       mPrimaryVertex ? mPrimaryVertex->getZ() : -999.);

  if (mMCMatching) {
    runMCMatching();
  } else {
    runDataMatching();
  }

  writeTracks();
}
