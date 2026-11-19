#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <map>
#include <cmath>
#include "NA6PTrack.h"
#include "MagneticField.h"
#include "NA6PMuonSpecModularHit.h"
#include "NA6PLayoutParam.h"
#endif

void checkMSPurity(const char* dirSimu = ".")
{
  auto magField = new MagneticField();
  magField->loadField();
  magField->setAsGlobalField();

  TFile* ft = new TFile("TracksMuonSpec.root");
  if (!ft)
    return;
  TTree* trTree = (TTree*)ft->Get("tracksMuonSpec");
  std::vector<NA6PTrack>* trArr = nullptr;
  trTree->SetBranchAddress("MuonSpec", &trArr);

  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fh = new TFile(Form("%s/HitsMuonSpecModular.root", dirSimu));
  TTree* th = (TTree*)fh->Get("hitsMuonSpecModular");
  std::vector<NA6PMuonSpecModularHit> msHits, *msHitsPtr = &msHits;
  th->SetBranchAddress("MuonSpecModular", &msHitsPtr);

  TH1F* hNBadClusters = new TH1F("hNBadClusters", "Number of bad clusters per track;N bad clusters;Entries", 10, -0.5, 9.5);
  TH1F* hPDGBadCluster = new TH1F("hPDGBadCluster", "PDG of tracks with bad clusters;PDG code;Entries", 51, -25.5, 25.5);
  TH1F* hPDGTrack = new TH1F("hPDGTrack", "Track PDG;PDG code;Entries", 51, -25.5, 25.5);
  TH1F* hLayerBadCluster = new TH1F("hLayerBadCluster", "Layer of bad clusters;MS layer;Wrong cluster frequency (%)", 6, -0.5, 5.5);
  TH1F* hPDGBadClusterMother = new TH1F("hPDGBadClusterMother", "PDG of tracks with bad clusters;PDG code;Entries", 51, -25.5, 25.5);
  TH1F* hBadClusterXOrigin = new TH1F("hBadClusterXOri", "X origin of bad clusters;X (cm);Entries", 200, -200, 200);
  TH1F* hBadClusterYOrigin = new TH1F("hBadClusterYOri", "Y origin of bad clusters;Y (cm);Entries", 200, -200, 200);
  TH1F* hBadClusterZOrigin = new TH1F("hBadClusterZOri", "Z origin of bad clusters;Z (cm);Entries", 100000, 290, 900);
  float pMax = 30;
  TH1F* hPurityVsP = new TH1F("hPurityVsP", ";#it{p} (GeV/#it{c});Purity", 10, 0, pMax);
  TH1F* hRecTrackVsP = new TH1F("hRecTrackVsP", "Reconstructed tracks vs P;P (GeV/c);Entries", 10, 0, pMax);
  // do one histogram per layer
  std::vector<TH1F*> hHitsDeltaPerLayer;
  hHitsDeltaPerLayer.reserve(NA6PTrack::kMaxLr);
  std::vector<TH1F*> hHitsDeltaXPerLayer;
  hHitsDeltaXPerLayer.reserve(NA6PTrack::kMaxLr);
  std::vector<TH1F*> hHitsDeltaYPerLayer;
  hHitsDeltaYPerLayer.reserve(NA6PTrack::kMaxLr);
  for (int lr = 0; lr < NA6PTrack::kMaxLr; ++lr) {
    hHitsDeltaPerLayer.push_back(new TH1F(Form("hHitsDelta_L%d", lr),
                                          Form("Delta between associated hit and true hit (L%d);#DeltaR (cm);Entries", lr - 5),
                                          100, 0, 0.5));
    hHitsDeltaXPerLayer.push_back(new TH1F(Form("hHitsDeltaX_L%d", lr),
                                           Form("Delta X between associated hit and true hit (L%d);#DeltaX (cm);Entries", lr - 5),
                                           100, -0.5, 0.5));
    hHitsDeltaYPerLayer.push_back(new TH1F(Form("hHitsDeltaY_L%d", lr),
                                           Form("Delta Y between associated hit and true hit (L%d);#DeltaY (cm);Entries", lr - 5),
                                           100, -0.5, 0.5));
  }
  int nEv = trTree->GetEntries();
  printf("Number of events = %d\n", nEv);
  for (int jEv = 0; jEv < nEv; jEv++) {
    mcTree->GetEvent(jEv);
    trTree->GetEvent(jEv);
    th->GetEvent(jEv);
    int nPart = mcArr->size();
    int nTracks = trArr->size();
    int nHits = msHits.size();

    for (int jp = 0; jp < nTracks; jp++) {
      auto curTrack = trArr->at(jp);
      int partID = curTrack.getParticleID();
      float pTrack = curTrack.getP();
      hRecTrackVsP->Fill(pTrack);
      if (partID < 0) {
        partID = abs(partID);
        // Check the PDG code of the particles associated to the clusters
        uint32_t clusterMap = curTrack.getClusterMap();
        std::map<int, int> pdgCounts; // PDG code -> count
        int badClusterCount = 0;
        // printf("Track %d (partID=%d, p=%.2f GeV/c) cluster PDG codes:\n", jp, partID, pTrack);
        for (int lr = 0; lr < NA6PTrack::kMaxLr; ++lr) {
          if (!(clusterMap & (1u << lr)))
            continue;

          int particleLabel = curTrack.getParticleLabel(lr);
          int clusterIndex = curTrack.getClusterIndex(lr);
          // check if partID left an hit in the layer
          bool hasHitFromPartID = false;

          double corrHitX = 0;
          double corrHitY = 0;
          const auto& layout = NA6PLayoutParam::Instance();
          // find a hit from the same MC track in this layer
          for (const auto& msHit : msHits) {
            if (msHit.getTrackID() != partID) {
              continue;
            }
            const float z = msHit.getZ();
            float minDist = std::abs(z - layout.posMSPlaneZ[0]);
            int hitLayerGlobal = layout.nVerTelPlanes;
            for (int i = 1; i < layout.nMSPlanes; ++i) {
              float dist = std::abs(z - layout.posMSPlaneZ[i]);
              if (dist < minDist) {
                minDist = dist;
                hitLayerGlobal = i + layout.nVerTelPlanes;
              }
            }
            const int hitLayerMS = hitLayerGlobal - layout.nVerTelPlanes;
            if (lr != hitLayerGlobal && lr != hitLayerMS) {
              continue;
            }
            hasHitFromPartID = true;
            corrHitX = msHit.getX();
            corrHitY = msHit.getY();
            break;
          }

          if (particleLabel >= 0 && particleLabel < mcArr->size()) {

            int pdgCode = mcArr->at(particleLabel).GetPdgCode();
            pdgCounts[pdgCode]++;
            // printf("  Layer %d: cluster=%d, particle=%d, PDG=%d\n",
            //      lr, clusterIndex, particleLabel, pdgCode);
            // get the cluster position and compute deltaR wrt the cluster of the tracks associated particle in the layer

            if (particleLabel != partID) {
              const auto& hit = msHits[clusterIndex];
              if (hasHitFromPartID) {
                double hitX = hit.getX();
                double hitY = hit.getY();

                float deltaR = std::sqrt(std::pow(hitX - corrHitX, 2) +
                                         std::pow(hitY - corrHitY, 2));
                if (lr >= 0 && lr < static_cast<int>(hHitsDeltaPerLayer.size())) {
                  hHitsDeltaPerLayer[lr]->Fill(deltaR);
                  hHitsDeltaXPerLayer[lr]->Fill(hitX - corrHitX);
                  hHitsDeltaYPerLayer[lr]->Fill(hitY - corrHitY);
                }
              }
              badClusterCount++;
              hPDGBadCluster->Fill(pdgCode);
              hLayerBadCluster->Fill(lr - 5);
              // retrieve mother particle PDG
              if (particleLabel >= 0 && particleLabel < mcArr->size()) {
                int motherLabel = mcArr->at(particleLabel).GetFirstMother();

                float xorigin = mcArr->at(particleLabel).Vx();
                float yorigin = mcArr->at(particleLabel).Vy();
                float zorigin = mcArr->at(particleLabel).Vz();
                hBadClusterXOrigin->Fill(xorigin);
                hBadClusterYOrigin->Fill(yorigin);
                hBadClusterZOrigin->Fill(zorigin);
                if (motherLabel >= 0 && motherLabel < mcArr->size()) {
                  int motherPDG = mcArr->at(motherLabel).GetPdgCode();
                  hPDGBadClusterMother->Fill(motherPDG);
                  // printf("    -> Bad cluster from particle %d (PDG=%d), mother %d (PDG=%d)\n",
                  //   particleLabel, pdgCode, motherLabel, motherPDG);
                }
              }
            }
          }
        }
        hNBadClusters->Fill(badClusterCount);
        // get pdg code from the PID
        hPDGTrack->Fill(mcArr->at(partID).GetPdgCode());

      } else
        hPurityVsP->Fill(pTrack);
    }
  }

  hPurityVsP->Divide(hRecTrackVsP);

  gStyle->SetOptStat(0);
  TFile* fout = new TFile(Form("%s/PurityMuonSpec.root", dirSimu), "RECREATE");

  auto setAxisStyle = [](TH1* hist) {
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetLabelFont(42);
  };


  setAxisStyle(hPDGBadCluster);
  setAxisStyle(hPDGBadClusterMother);
  setAxisStyle(hPurityVsP);
  setAxisStyle(hLayerBadCluster);
  setAxisStyle(hNBadClusters);
  setAxisStyle(hPDGTrack);
  setAxisStyle(hBadClusterXOrigin);
  setAxisStyle(hBadClusterYOrigin);
  setAxisStyle(hBadClusterZOrigin);


  hPDGBadCluster->Write();
  hPDGBadClusterMother->Write();
  hPurityVsP->Write();
  hLayerBadCluster->Scale(100. / hRecTrackVsP->GetEntries());
  hLayerBadCluster->Write();
  hNBadClusters->Write();
  hPDGTrack->Write();
  hBadClusterXOrigin->Write();
  hBadClusterYOrigin->Write();
  hBadClusterZOrigin->Write();
  for (auto* h : hHitsDeltaPerLayer) {
    setAxisStyle(h);
    h->Write();
  }
  for (auto* h : hHitsDeltaXPerLayer) {
    setAxisStyle(h);
    h->Write();
  }
  for (auto* h : hHitsDeltaYPerLayer) {
    setAxisStyle(h);
    h->Write();
  }

  TCanvas* c1 = new TCanvas("c1", "Purity vs P", 1200, 800);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.15);
  gPad->SetTopMargin(0.08);
  hPurityVsP->SetMarkerStyle(20);
  hPurityVsP->SetMarkerSize(1.5);
  hPurityVsP->SetMarkerColor(kBlue);
  hPurityVsP->Draw();
  //set the margin e axis title size
  c1->SaveAs(Form("%s/PurityVsP.png", dirSimu));

  hPDGBadCluster->Draw();
  c1->SaveAs(Form("%s/PDGBadCluster.png", dirSimu));
  hPDGBadClusterMother->Draw();
  c1->SaveAs(Form("%s/PDGBadClusterMother.png", dirSimu));
  hLayerBadCluster->Draw();
  c1->SaveAs(Form("%s/LayerBadCluster.png", dirSimu));
  hNBadClusters->Draw();
  c1->SaveAs(Form("%s/NBadClusters.png", dirSimu));
  hPDGTrack->Draw();
  c1->SaveAs(Form("%s/PDGTrack.png", dirSimu));
  for(int lr = 0; lr < NA6PTrack::kMaxLr; ++lr) {
    if (hHitsDeltaPerLayer[lr]->GetEntries() == 0)
      continue;
    hHitsDeltaPerLayer[lr]->Draw();
    c1->SaveAs(Form("%s/HitsDelta_L%d.png", dirSimu, lr));
    hHitsDeltaXPerLayer[lr]->Draw();
    c1->SaveAs(Form("%s/HitsDeltaX_L%d.png", dirSimu, lr));
    hHitsDeltaYPerLayer[lr]->Draw();
    c1->SaveAs(Form("%s/HitsDeltaY_L%d.png", dirSimu, lr));
  }
  hBadClusterZOrigin->Draw();
  c1->SaveAs(Form("%s/BadClusterZOrigin.png", dirSimu));
  hBadClusterXOrigin->Draw();
  c1->SaveAs(Form("%s/BadClusterXOrigin.png", dirSimu));
  hBadClusterYOrigin->Draw();
  c1->SaveAs(Form("%s/BadClusterYOrigin.png", dirSimu));
  fout->Close();
}
