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
#include <TGeoManager.h>
#include <TStyle.h>
#include <TSystem.h>
#include "NA6PTrack.h"
#include "MagneticField.h"
#include "NA6PMuonSpecModularHit.h"
#include "NA6PLayoutParam.h"
#include "NA6PFastTrackFitter.h"
#endif

void plotMSDimuons2(const char* dirSimu = ".", bool Matching = kFALSE, bool FullMSTracking = kFALSE, int NRecoClusters = 4, bool VertexConstrained = kTRUE)
{

  TGeoManager::Import(Form("%s/geometry.root", dirSimu));
  if (!gGeoManager) {
    std::cerr << "Error: Could not load geometry" << std::endl;
    return;
  }

  // Load magnetic field from .ini file
  TString fileName = Form("%s/na6pLayout.ini", dirSimu);

  if (gSystem->AccessPathName(fileName)) {
    std::cerr << "ERROR: File not found: " << fileName << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto& param = NA6PLayoutParam::Instance();
  param.updateFromFile(std::string(fileName.Data()), "", true);

  auto magField = new MagneticField();
  magField->loadField();
  magField->setAsGlobalField();

  TFile* ft = nullptr;

  if (FullMSTracking && !Matching)
    ft = TFile::Open(Form("%s/TracksFourPlusTwo.root", dirSimu));
  else if (!FullMSTracking && !Matching)
    ft = TFile::Open(Form("%s/TracksMuonSpec.root", dirSimu));
  else if (Matching)
    ft = TFile::Open(Form("%s/TracksMatching.root", dirSimu));
  else if (!ft || ft->IsZombie()) {
    std::cerr << "ERROR: Could not open tracks file" << std::endl;
    return;
  }

  TTree* trTree;
  if (Matching)
    trTree = dynamic_cast<TTree*>(ft->Get("tracksMatching"));
  else if (!Matching)
    trTree = dynamic_cast<TTree*>(ft->Get("tracksMuonSpec"));
  if (!trTree) {
    std::cerr << "ERROR: tree not found in file" << std::endl;
    ft->Close();
    return;
  }
  std::vector<NA6PTrack>* trArr = nullptr;
  if (Matching)
    trTree->SetBranchAddress("Matching", &trArr);
  else if (!Matching)
    trTree->SetBranchAddress("MuonSpec", &trArr);

  TFile* fk = TFile::Open(Form("%s/MCKine.root", dirSimu));
  if (!fk || fk->IsZombie()) {
    std::cerr << "ERROR: Could not open MCKine.root" << std::endl;
    ft->Close();
    return;
  }
  TTree* mcTree = dynamic_cast<TTree*>(fk->Get("mckine"));
  if (!mcTree) {
    std::cerr << "ERROR: mckine tree not found in MCKine.root" << std::endl;
    fk->Close();
    ft->Close();
    return;
  }
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fh = TFile::Open(Form("%s/HitsMuonSpecModular.root", dirSimu));
  if (!fh || fh->IsZombie()) {
    std::cerr << "ERROR: Could not open HitsMuonSpecModular.root" << std::endl;
    fk->Close();
    ft->Close();
    return;
  }
  TTree* th = dynamic_cast<TTree*>(fh->Get("hitsMuonSpecModular"));
  if (!th) {
    std::cerr << "ERROR: hitsMuonSpecModular tree not found" << std::endl;
    fh->Close();
    fk->Close();
    ft->Close();
    return;
  }
  std::vector<NA6PMuonSpecModularHit> msHits, *msHitsPtr = &msHits;
  th->SetBranchAddress("MuonSpecModular", &msHitsPtr);

  TH1D* hdimuy = new TH1D("hdimuy", "dimuon y", 100, 0., 5.);
  TH1D* hdimupt = new TH1D("hdimupt", "dimuon pT", 100, 0., 5.);
  TH1D* hdimumass = new TH1D("hdimumass", "dimuon m", 250, 0., 5.);
  TH1D* hdimumass_MCtruth = new TH1D("hdimumass_MCtruth", "dimuon m MC truth", 200, 0., 4.);
  TH1D* hdimuygen = new TH1D("hdimuygen", "dimuon y gen", 100, 0., 5.);
  TH1D* hdimuptgen = new TH1D("hdimuptgen", "dimuon pT gen", 100, 0., 5.);
  TH1D* hdimumassgen = new TH1D("hdimumassgen", "dimuon m gen", 100, 0., 5.);
  TH2D* hdimuyrecvsgen = new TH2D("hdimuyrecvsgen", "dimuon y rec vs gen", 100, 1., 5., 100, 1., 5.);
  TH2D* hdimuptrecvsgen = new TH2D("hdimuptrecvsgen", "dimuon pt rec vs gen", 100, 0., 5., 100, 0., 5.);
  TH1D* hdimudeltayrecvsgen = new TH1D("hdimudeltayrecvsgen", "dimuon y rec - gen", 100, -1., 1.);
  TH1D* hdimudeltaptrecvsgen = new TH1D("hdimudeltaptrecvsgen", "dimuon pt rec - gen", 100, -6., 6.);
  TH2D* hdimurecmassvspt = new TH2D("hdimurecmassvspt", "dimuon rec mass vs pt", 200, 0., 4., 100, 0., 5.);
  TH2D* hdimurecmassvseta = new TH2D("hdimurecmassvseta", "dimuon rec mass vs eta", 200, 0., 4., 100, 1., 5.);
  TH2D* hdimurecptvseta = new TH2D("hdimurecptvseta", "dimuon rec pt vs eta", 100, 0., 5., 100, 1., 5.);
  TH2D* hdimurecpvseta = new TH2D("hdimurecpvseta", "dimuon rec p vs eta", 100, 0., 30., 100, 1., 5.);

  int nEv = trTree->GetEntries();
  printf("Number of events = %d\n", nEv);

  for (int jEv = 0; jEv < nEv; jEv++) {
    if(jEv % 10000 == 0) printf("Event %d\n", jEv);
    mcTree->GetEvent(jEv);
    trTree->GetEvent(jEv);
    th->GetEvent(jEv);
    size_t nPart = mcArr->size();
    size_t nTracks = trArr->size();
    size_t nHits = msHits.size();
    // printf("\nEvent %d particles = %d tracks = %d nHits = %d\n", jEv, nPart, nTracks, nHits);
    double zOrig = 0.;
    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      if (curPart.IsPrimary()) {
        zOrig = curPart.Vz();
        break;
      }
    }
    int nTrackEventNCl = 0;

    for (int jTr = 0; jTr < nTracks; ++jTr) {
      NA6PTrack tr = trArr->at(jTr);
      int nClusters = tr.getNHits();
      uint32_t clustermap = tr.getClusterMap();

      if (!Matching) {
        if (nClusters < NRecoClusters)
          continue;
        uint32_t maskClusters;
        if (NRecoClusters == 4)
          maskClusters = ((1u << 4) - 1) << 5;
        else if (NRecoClusters == 6)
          maskClusters = ((1u << 6) - 1) << 5;
        if ((clustermap & maskClusters) != maskClusters) // last six layers all have ≥1 hit, independentely from the first 5
          continue;
        nTrackEventNCl++;
      }
    }

    TLorentzVector mothgen;
    TLorentzVector moth;
    TLorentzVector decprod0;
    TLorentzVector decprod1;
    TLorentzVector decprodgen0;
    TLorentzVector decprodgen1;
    int ndec = 0;
    double RecXYZ[3];
    double RecP;
    double RecMom[3];

    if (!Matching && nTrackEventNCl < 2)
      continue;
    if (Matching && nTracks < 2)
      continue;

    int mcLabel0 = -9999;
    int nskip = 0;
    for (int jTr1 = 0; jTr1 < nTracks; ++jTr1) {
      NA6PTrack tr = trArr->at(jTr1);
      NA6PFastTrackFitter fitter;
      if (!VertexConstrained && !Matching)
        fitter.propagateToZ(&tr, zOrig); // to correct for material budget
      if (VertexConstrained && !Matching) {
        tr.getPXYZVertexConstrained(RecMom);
        tr.getXYZVertexConstrained(RecXYZ);
        RecP = tr.getPVertexConstrained();
        if (!tr.getStatusConstrained())
          continue;
      } else {
        tr.getPXYZ(RecMom);
        RecXYZ[0] = tr.getXLab();
        RecXYZ[1] = tr.getYLab();
        RecXYZ[2] = tr.getZLab();
        RecP = tr.getP();
      }

      mcLabel0 = tr.getParticleID();

      if (TMath::Abs(mcLabel0) >= nPart) {
        printf("mcLabel0 = %d out of range (nPart = %d), %d skipping track\n", mcLabel0, nPart, nskip);
        nskip++;
        continue; // CHECK SE HA SENSO
      }
        decprod0.SetPxPyPzE(RecMom[0], RecMom[1], RecMom[2], TMath::Sqrt(RecP * RecP + 0.10566 * 0.10566));

        if (mcLabel0 > 0) { // check if MC truth exists
          TParticle part = mcArr->at(std::abs(mcLabel0));
          decprodgen0.SetPxPyPzE(part.Px(), part.Py(), part.Pz(), part.Energy());
        }
        double RecMom2[3];
        double RecXYZ2[3];
        double RecP2;

        for (int jTr2 = jTr1 + 1; jTr2 < nTracks; ++jTr2) {
          NA6PTrack tr = trArr->at(jTr2);
          NA6PFastTrackFitter fitter;
          if (!VertexConstrained && !Matching)
            fitter.propagateToZ(&tr, zOrig); // to correct for material budget
          if (VertexConstrained && !Matching) {
            tr.getPXYZVertexConstrained(RecMom2);
            tr.getXYZVertexConstrained(RecXYZ2);
            RecP2 = tr.getPVertexConstrained();
            if (!tr.getStatusConstrained())
              continue;
          } else {
            tr.getPXYZ(RecMom2);
            RecXYZ2[0] = tr.getXLab();
            RecXYZ2[1] = tr.getYLab();
            RecXYZ2[2] = tr.getZLab();
            RecP2 = tr.getP();
          }

          decprod1.SetPxPyPzE(RecMom2[0], RecMom2[1], RecMom2[2], TMath::Sqrt(RecP2 * RecP2 + 0.10566 * 0.10566));

          int mcLabel1 = tr.getParticleID();

          if (TMath::Abs(mcLabel1) >= nPart) {
            printf("mcLabel1 = %d out of range (nPart = %zu), %d skipping track\n", mcLabel1, nPart, nskip);
            nskip++;
            continue; // CHECK SE HA SENSO
          }
          if (mcLabel1 > 0) {
            TParticle part = mcArr->at(std::abs(mcLabel1));
            decprodgen1.SetPxPyPzE(part.Px(), part.Py(), part.Pz(), part.Energy());
          }

          // Build dimuon and fill histos
          moth = decprod0 + decprod1;
          if (mcLabel0 > 0 && mcLabel1 > 0)
            mothgen = decprodgen0 + decprodgen1;

          hdimuy->Fill(moth.Rapidity());
          hdimupt->Fill(moth.Pt());
          hdimumass->Fill(moth.M());
          if (mcLabel0 > 0 && mcLabel1 > 0) {
            hdimumass_MCtruth->Fill(moth.M());
            hdimumassgen->Fill(mothgen.M());
            hdimuyrecvsgen->Fill(moth.Rapidity(), mothgen.Rapidity());
            hdimuptrecvsgen->Fill(moth.Pt(), mothgen.Pt());
            hdimudeltayrecvsgen->Fill(moth.Rapidity() - mothgen.Rapidity());
            hdimudeltaptrecvsgen->Fill(moth.Pt() - mothgen.Pt());
          }
          hdimurecmassvspt->Fill(moth.M(), moth.Pt());
          hdimurecmassvseta->Fill(moth.M(), moth.Rapidity());
          hdimurecptvseta->Fill(moth.Pt(), moth.Rapidity());
          hdimurecpvseta->Fill(moth.P(), moth.Rapidity());
        }
      }
      //}
    }

    TCanvas* cDimu = new TCanvas("cDimu", "Dimum mass", 1200, 500);
    cDimu->Divide(2, 1);
    cDimu->cd(1);
    hdimumass->SetMarkerColor(kAzure + 8);
    hdimumass->SetLineColor(kBlue + 1);
    hdimumass->SetMarkerStyle(20);
    hdimumass->SetTitle("Dimuon Mass Spectrum");
    hdimumass->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    hdimumass->GetYaxis()->SetTitle("Counts / 0.1 GeV/c^{2}");
    hdimumass->Fit("gaus", "", " ", 2.85, 3.35);
    hdimumass->Draw("EP");
    gStyle->SetOptFit(1111);

    cDimu->cd(2);
    hdimumass_MCtruth->SetMarkerColor(kAzure + 8);
    hdimumass_MCtruth->SetLineColor(kBlue + 1);
    hdimumass_MCtruth->SetMarkerStyle(20);
    hdimumass_MCtruth->SetTitle("Dimuon Mass Spectrum");
    hdimumass_MCtruth->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    hdimumass_MCtruth->GetYaxis()->SetTitle("Counts / 0.1 GeV/c^{2}");
    hdimumass_MCtruth->Fit("gaus", "", " ", 2.85, 3.35);
    hdimumass_MCtruth->Draw("EP");

    TCanvas* cDimu2 = new TCanvas("cDimu2", "Dimuon rec vs gen", 1400, 900);
    cDimu2->Divide(2, 2);
    cDimu2->cd(1);
    hdimuyrecvsgen->Draw("COLZ");
    cDimu2->cd(2);
    hdimuptrecvsgen->Draw("COLZ");
    cDimu2->cd(3);
    hdimudeltayrecvsgen->Draw();
    cDimu2->cd(4);
    hdimudeltaptrecvsgen->Draw();

    TCanvas* cDimu3 = new TCanvas("cDimu3", "Dimuon Rec", 1400, 900);
    cDimu3->Divide(2, 2);
    cDimu3->cd(1);
    hdimurecmassvspt->Draw("COLZ");
    cDimu3->cd(2);
    hdimurecmassvseta->Draw("COLZ");
    cDimu3->cd(3);
    hdimurecptvseta->Draw("COLZ");
    cDimu3->cd(4);
    hdimurecpvseta->Draw("COLZ");

    TFile* fout = new TFile(Form("%s/DimuonHistos.root", dirSimu), "recreate");
    hdimumass->Write();
    fout->Close();
  }