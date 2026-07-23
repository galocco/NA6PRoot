#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TVector3.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "NA6PMuonSpecCluster.h"
#include "NA6PTrack.h"
#include "NA6PLine.h"
#include "NA6PVertex.h"
#include "NA6PTrackerCA.h"
#include "MagneticField.h"
#include "NA6PMCTruthContainer.h"
#include "NA6PReconstruction.h"
#endif

void runMSTrackMIDTrackletMatching(int firstEv = 0,
                                   int lastEv = 9999,
                                   const char* dirSimu = ".")
{

  auto magField = new MagneticField();
  magField->loadField();
  magField->setAsGlobalField();

  int nMomBins = 40;
  float pMax = 30.0;

  NA6PTrackerCA* tracker = new NA6PTrackerCA();
  tracker->setNLayers(4);
  tracker->setStartLayer(5);
  tracker->setDoOutwardPropagation(true);
  tracker->setZForOutwardPropagation(810.);
  // tracker->setVerbosity(true);
  if (!tracker->loadGeometry(Form("%s/geometry.root", dirSimu)))
    return;
  // pass here the configuration of the tracker via an ini file
  // tracker->configureFromRecoParam(/* filename.ini */);
  // alternatively the configuration can be set calling setters for the iterations
  tracker->setNumberOfIterations(2);
  tracker->setIterationParams(0, 0.06, 0.1, 6., 0.6, 0.05, 0.05, 5., 5., 5., 4);
  tracker->setIterationParams(1, 0.1, 0.6, 9., 0.8, 0.08, 0.08, 10., 10., 10., 4);
  tracker->printConfiguration();

  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fc = new TFile(Form("%s/ClustersMuonSpec.root", dirSimu));
  printf("Open cluster file: %s\n", fc->GetName());
  TTree* tc = (TTree*)fc->Get("clustersMuonSpec");
  std::vector<NA6PMuonSpecCluster> msClus, *msClusPtr = &msClus;
  NA6PMCTruthContainer vtCluMCLabels, *vtCluMCLabelsPtr = &vtCluMCLabels;
  tc->SetBranchAddress("MuonSpec", &msClusPtr);
  tc->SetBranchAddress("MuonSpecMCTruth", &vtCluMCLabelsPtr);

  int nEv = tc->GetEntries();
  if (lastEv > nEv || lastEv < 0)
    lastEv = nEv;
  if (firstEv < 0)
    firstEv = 0;
  NA6PVertex primVert;

  int nIterationsCA = tracker->getNIterations();
  TH1F* hDistXYGood = new TH1F("hDistXYGood", "", 100, 0., 100.);
  TH1F* hCosThGood = new TH1F("hCosThXYGood", "", 100, 0., 1.);
  TH1F* h1mCosThGood = new TH1F("h1mCosThXYGood", "", 100, 0., 0.01);
  TH1F* hDistXYComb = new TH1F("hDistXYComb", "", 100, 0., 100.);
  TH1F* hCosThComb = new TH1F("hCosThXYComb", "", 100, 0., 1.);
  TH1F* hNGoodFound = new TH1F("hNGoodFound", "", 10, -0.5, 9.5);
  TH1F* hNMatches = new TH1F("hNMatches", "", 10, -0.5, 9.5);
  TH1F* hNGoodMatches = new TH1F("hNGoodMatches", "", 10, -0.5, 9.5);
  TH1F* hMatchType = new TH1F("hMatchType", "", 5, -0.5, 4.5);
  hMatchType->GetXaxis()->SetBinLabel(1, "Correct match");
  hMatchType->GetXaxis()->SetBinLabel(2, "Fake track, correct match");
  hMatchType->GetXaxis()->SetBinLabel(3, "Good track, fake match");
  hMatchType->GetXaxis()->SetBinLabel(4, "Fake track, fake match");
  hMatchType->GetXaxis()->SetBinLabel(5, "No match");

  TH1F* hEtaOrigTracks = new TH1F("hEtaOrigTracks", ";#eta;counts", 20, 1., 5.);
  TH1F* hEtaRefitTracks = new TH1F("hEtaRefitTracks", ";#eta;counts", 20, 1., 5.);
  TH1F* hNClusOrigTracks = new TH1F("hNClusOrigTracks", ";n_{clusters};counts", 7, -0.5, 6.5);
  TH1F* hNClusRefitTracks = new TH1F("hNClusRefitTracks", ";n_{clusters};counts", 7, -0.5, 6.5);

  std::vector<NA6PTrack> ms42Tracks, *hTrackPtr = &ms42Tracks;
  std::vector<NA6PMCComposedLabel> ms42TrackLabs, *hTrackLabPtr = &ms42TrackLabs;
  TFile* fouttr = TFile::Open("TracksFourPlusTwo.root", "recreate");
  TTree* trackTree = new TTree("tracksMuonSpec", "MuonSpec Tracks");
  trackTree->Branch("MuonSpec", &hTrackPtr);
  trackTree->Branch("MuonSpecMCTruth", &hTrackLabPtr);

  NA6PReconstruction rec("dummy");
  std::vector<NA6PMCComposedLabel> trkMCLabs;

  int nFailed = 0;

  for (int jEv = firstEv; jEv < lastEv; jEv++) {
    ms42Tracks.clear();
    ms42TrackLabs.clear();
    mcTree->GetEvent(jEv);
    int nPart = mcArr->size();
    double xvert = 0.;
    double yvert = 0.;
    double zvert = 0;
    // get primary vertex position from the Kine Tree
    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      if (curPart.IsPrimary()) {
        xvert = curPart.Vx();
        yvert = curPart.Vy();
        zvert = curPart.Vz();
        break;
      }
    }
    primVert.setXYZ(xvert, yvert, zvert);
    tc->GetEvent(jEv);
    for (size_t i = 0; i < msClus.size(); ++i) {
      msClus[i].setClusterIndex(static_cast<int>(i));
    }
    tracker->findTracks(msClus, &primVert);
    std::vector<NA6PTrack> trks = tracker->getTracks();
    rec.assignMCLabels(trks, trkMCLabs, vtCluMCLabels);
    int nTrks = trks.size();
    std::vector<std::pair<NA6PMuonSpecCluster, NA6PMuonSpecCluster>> trkltsMID = tracker->findTracklets(9, 10, msClus, &primVert);
    int nTrklets = trkltsMID.size();
    NA6PFastTrackFitter* fitter = tracker->getTrackFitter();
    auto clusterMatch = [](const auto& clu, int trackID, const auto& cluMCLabels) {
      int cluInd = clu.getClusterIndex();
      std::span labels = cluMCLabels.getLabels(cluInd);
      int nLabels = labels.size();
      for (int jLab = 0; jLab < nLabels; jLab++) {
        NA6PMCComposedLabel lbl = labels[jLab];
        if (lbl.getTrackID() == trackID)
          return true;
      }
      return false;
    };

    // first loop for counting and histos
    for (int jT = 0; jT < nTrks; jT++) {
      NA6PTrack tr = trks[jT];
      NA6PMCComposedLabel& lblTr = trkMCLabs[jT];
      int nMatches = 0;
      int nGoodFound = 0;
      int nGoodMatches = 0;
      for (int jS = 0; jS < nTrklets; jS++) {
        const auto& tracklet = trkltsMID[jS];
        NA6PMuonSpecCluster clu1 = tracklet.first;
        NA6PMuonSpecCluster clu2 = tracklet.second;
        float xyzClu1[3] = {clu1.getX(), clu1.getY(), clu1.getZ()};
        float xyzClu2[3] = {clu2.getX(), clu2.getY(), clu2.getZ()};
        if (xyzClu1[2] > xyzClu2[2]) {
          float tmp[3] = {xyzClu1[0], xyzClu1[1], xyzClu1[2]};
          for (int j = 0; j < 3; ++j) {
            xyzClu1[j] = xyzClu2[j];
            xyzClu2[j] = tmp[j];
          }
        }
        float zToProp = xyzClu1[2];
        float dirSegm[3];
        float norm = 0.f;
        for (int j = 0; j < 3; ++j) {
          dirSegm[j] = xyzClu2[j] - xyzClu1[j];
          norm += dirSegm[j] * dirSegm[j];
        }
        norm = std::sqrt(norm);
        for (int j = 0; j < 3; ++j)
          dirSegm[j] /= norm;

        auto trOuter = tr.getOuterParam();
        if (!Propagator::Instance()->propagateToZ(trOuter, zToProp, fitter->getPropOpt())) {
          continue;
        }
        auto xyzTr = trOuter.getXYZ<float>();
        auto pxyzTr = trOuter.getPXYZ<float>();
        norm = 0.f;
        for (int j = 0; j < 3; ++j)
          norm += pxyzTr[j] * pxyzTr[j];
        norm = std::sqrt(norm);
        for (int j = 0; j < 3; ++j)
          pxyzTr[j] /= norm;
        float distXY = std::sqrt((xyzClu1[0] - xyzTr[0]) * (xyzClu1[0] - xyzTr[0]) + (xyzClu1[1] - xyzTr[1]) * (xyzClu1[1] - xyzTr[1]));
        float costh = dirSegm[0] * pxyzTr[0] + dirSegm[1] * pxyzTr[1] + dirSegm[2] * pxyzTr[2];
        int idTrack = lblTr.getTrackID();
        bool ok1 = clusterMatch(clu1, idTrack, vtCluMCLabels);
        bool ok2 = clusterMatch(clu2, idTrack, vtCluMCLabels);
        if (!lblTr.isFake() && ok1 && ok2) {
          nGoodFound++;
          hDistXYGood->Fill(distXY);
          hCosThGood->Fill(costh);
          h1mCosThGood->Fill(1 - costh);
        } else {
          hDistXYComb->Fill(distXY);
          hCosThComb->Fill(costh);
        }
        if (distXY < 7.5 && costh > 0.999) {
          nMatches++;
          if (!lblTr.isFake() && ok1 && ok2) {
            nGoodMatches++;
          }
        }
      }
      hNMatches->Fill(nMatches);
      hNGoodMatches->Fill(nGoodMatches);
      hNGoodFound->Fill(nGoodFound);
    }

    // matching loop
    struct MatchCandidate {
      int track;
      int segment;
      float score;
    };
    std::vector<MatchCandidate> candidates;
    for (int jT = 0; jT < nTrks; jT++) {
      NA6PTrack tr = trks[jT];
      NA6PMCComposedLabel& lblTr = trkMCLabs[jT];
      for (int jS = 0; jS < nTrklets; jS++) {
        const auto& tracklet = trkltsMID[jS];
        NA6PMuonSpecCluster clu1 = tracklet.first;
        NA6PMuonSpecCluster clu2 = tracklet.second;
        float xyzClu1[3] = {clu1.getX(), clu1.getY(), clu1.getZ()};
        float xyzClu2[3] = {clu2.getX(), clu2.getY(), clu2.getZ()};
        if (xyzClu1[2] > xyzClu2[2]) {
          float tmp[3] = {xyzClu1[0], xyzClu1[1], xyzClu1[2]};
          for (int j = 0; j < 3; ++j) {
            xyzClu1[j] = xyzClu2[j];
            xyzClu2[j] = tmp[j];
          }
        }
        float zToProp = xyzClu1[2];
        float dirSegm[3];
        float norm = 0.f;
        for (int j = 0; j < 3; ++j) {
          dirSegm[j] = xyzClu2[j] - xyzClu1[j];
          norm += dirSegm[j] * dirSegm[j];
        }
        norm = std::sqrt(norm);
        for (int j = 0; j < 3; ++j)
          dirSegm[j] /= norm;

        auto trOuter = tr.getOuterParam();
        if (!Propagator::Instance()->propagateToZ(trOuter, zToProp, fitter->getPropOpt())) {
          continue;
        }
        auto xyzTr = trOuter.getXYZ<float>();
        auto pxyzTr = trOuter.getPXYZ<float>();
        norm = 0.f;
        for (int j = 0; j < 3; ++j)
          norm += pxyzTr[j] * pxyzTr[j];
        norm = std::sqrt(norm);
        for (int j = 0; j < 3; ++j)
          pxyzTr[j] /= norm;
        float distXY = std::sqrt((xyzClu1[0] - xyzTr[0]) * (xyzClu1[0] - xyzTr[0]) + (xyzClu1[1] - xyzTr[1]) * (xyzClu1[1] - xyzTr[1]));
        float costh = dirSegm[0] * pxyzTr[0] + dirSegm[1] * pxyzTr[1] + dirSegm[2] * pxyzTr[2];
        if (distXY < 7.5 && costh > 0.999) {
          costh = std::clamp(costh, -1.f, 1.f);
          float score = distXY / 2. + (1 - costh) / 1.e-4;
          candidates.push_back({jT, jS, score});
        }
      }
    }
    std::sort(candidates.begin(), candidates.end(),
              [](auto& a, auto& b) { return a.score < b.score; });

    std::vector<bool> trackUsed(nTrks, false);
    std::vector<bool> segmentUsed(nTrklets, false);
    std::vector<int> trackMatch(nTrks, -1);

    for (auto& c : candidates) {
      if (trackUsed[c.track])
        continue;
      if (segmentUsed[c.segment])
        continue;
      trackMatch[c.track] = c.segment;
      trackUsed[c.track] = true;
      segmentUsed[c.segment] = true;
    }

    for (int jT = 0; jT < nTrks; jT++) {
      NA6PTrack tr = trks[jT];
      NA6PMCComposedLabel& lblTr = trkMCLabs[jT];
      hNClusOrigTracks->Fill(tr.getNHits());
      auto pxyz = tr.getPXYZ<float>();
      float momtr = tr.getP();
      float thetatr = std::acos(pxyz[2] / momtr);
      float etatr = -std::log(std::tan(thetatr / 2.));
      hEtaOrigTracks->Fill(etatr);

      int jS = trackMatch[jT];
      if (jS < 0) {
        hMatchType->Fill(4);
        continue;
      }
      const auto& tracklet = trkltsMID[jS];
      NA6PMuonSpecCluster clu1 = tracklet.first;
      NA6PMuonSpecCluster clu2 = tracklet.second;
      int idTrack = lblTr.getTrackID();
      bool ok1 = clusterMatch(clu1, idTrack, vtCluMCLabels);
      bool ok2 = clusterMatch(clu2, idTrack, vtCluMCLabels);
      if (!lblTr.isFake()) {
        if (ok1 && ok2)
          hMatchType->Fill(0);
        else
          hMatchType->Fill(2);
      } else {
        if (ok1 && ok2)
          hMatchType->Fill(1);
        else
          hMatchType->Fill(3);
      }

      fitter->cleanupAndStartFit();
      std::vector<int> clusterLookup(msClus.size(), -1);
      for (size_t i = 0; i < msClus.size(); ++i) {
        int originalID = msClus[i].getClusterIndex();
        if (originalID >= 0 && (size_t)originalID < msClus.size())
          clusterLookup[originalID] = static_cast<int>(i);
      }

      for (int jl = 0; jl < 4; ++jl) {
        int originalID = tr.getClusterIndex(jl + 5);
        if (originalID >= 0 && (size_t)originalID < clusterLookup.size()) {
          int jNewPos = clusterLookup[originalID];
          const auto& cl = msClus[jNewPos];
          fitter->addCluster(cl);
        }
      }

      NA6PTrack outTr = tr;
      static_cast<NA6PTrackParCov&>(outTr) = tr.getOuterParam();
      outTr.setChi2(tr.getChi2Outer());
      float zToProp = 810.;
      if (!Propagator::Instance()->propagateToZ(outTr, zToProp, fitter->getPropOpt())) {
        continue;
      }

      float zFirst = -999.;
      float zSecond = -999.;
      if (clu1.getZ() < clu2.getZ()) {
        fitter->addCluster(clu1);
        fitter->addCluster(clu2);
        zFirst = clu1.getZ();
        zSecond = clu2.getZ();
      } else {
        fitter->addCluster(clu2);
        fitter->addCluster(clu1);
        zFirst = clu2.getZ();
        zSecond = clu1.getZ();
      }
      //      fitter->printClusters();
      if (!Propagator::Instance()->propagateToZ(outTr, zFirst, fitter->getPropOpt())) {
        continue;
      }
      bool success = false;
      if (outTr.update(clu1)) {
        if (!Propagator::Instance()->propagateToZ(outTr, zSecond, fitter->getPropOpt())) {
          continue;
        }
        if (outTr.update(clu2)) {
          success = true;
        }
      }
      if (!success)
        continue;

      NA6PTrack refitInw = outTr;
      float chi2Refit = fitter->fitSeedInward(refitInw, true);
      if (chi2Refit < 0.f)
        continue;
      refitInw.setChi2(chi2Refit);
      refitInw.setOuterParam(outTr);
      fitter->addClustersToTrack(refitInw);

      hNClusRefitTracks->Fill(refitInw.getNHits());
      pxyz = refitInw.getPXYZ<float>();
      momtr = refitInw.getP();
      thetatr = std::acos(pxyz[2] / momtr);
      etatr = -std::log(std::tan(thetatr / 2.));
      hEtaRefitTracks->Fill(etatr);
      fitter->constrainTrackToVertex(refitInw, primVert);
      ms42Tracks.push_back(refitInw);
    }
    rec.assignMCLabels(ms42Tracks, ms42TrackLabs, vtCluMCLabels);
    trackTree->Fill();
  }
  fouttr->cd();
  trackTree->Write();
  fouttr->Close();

  hDistXYComb->Scale(1. / hDistXYComb->GetEntries());
  hDistXYGood->Scale(1. / hDistXYGood->GetEntries());
  hCosThComb->Scale(1. / hCosThComb->GetEntries());
  hCosThGood->Scale(1. / hCosThGood->GetEntries());

  TCanvas* cmatch = new TCanvas("cmatch", "", 1400, 800);
  cmatch->Divide(3, 2);
  cmatch->cd(1);
  hDistXYComb->SetLineColor(2);
  hDistXYComb->Draw("histo");
  hDistXYGood->Draw("histo,sames");
  cmatch->cd(2);
  hCosThComb->SetLineColor(2);
  hCosThComb->Draw("histo");
  hCosThGood->Draw("histo,sames");
  cmatch->cd(3);
  hMatchType->Draw("histo");
  hMatchType->Draw("text,same");
  cmatch->cd(4);
  hNMatches->Draw();
  hNMatches->Draw("text,same");
  cmatch->cd(5);
  hNGoodFound->Draw();
  hNGoodFound->Draw("text,same");
  cmatch->cd(6);
  hNGoodMatches->Draw();
  hNGoodMatches->Draw("text,same");

  TCanvas* ctr = new TCanvas("ctr", "", 1400, 800);
  ctr->Divide(3, 2);
  ctr->cd(1);
  hNClusOrigTracks->Draw();
  ctr->cd(2);
  hEtaOrigTracks->Draw();
  ctr->cd(4);
  hNClusRefitTracks->Draw();
  ctr->cd(5);
  hEtaRefitTracks->Draw();
  ctr->cd(6);
  TH1F* hRatioEta = (TH1F*)hEtaRefitTracks->Clone("hRatioEta");
  hRatioEta->Divide(hEtaRefitTracks, hEtaOrigTracks, 1., 1., "B");
  hRatioEta->Draw();
}
