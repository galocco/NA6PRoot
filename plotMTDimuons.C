#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TGeoManager.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include "NA6PTrack.h"
#include "MagneticField.h"
#include "NA6PLayoutParam.h"
#include "NA6PFastTrackFitter.h"
#endif

bool isMuonFrom(const NA6PTrack& trk, const std::vector<TParticle>& mcArr, int motherPdg) {
  int pid = trk.getParticleID();
  if (pid < 0 || pid >= static_cast<int>(mcArr.size())) return false;
  const TParticle& part = mcArr.at(pid);
  if (TMath::Abs(part.GetPdgCode()) != 13) return false;
  int motherIndex = part.GetFirstMother();
  if (motherIndex < 0 || motherIndex >= static_cast<int>(mcArr.size())) return false;
  const TParticle& mother = mcArr[motherIndex];
  return TMath::Abs(mother.GetPdgCode()) == motherPdg;
}

TH1D* getMassSigmaVsRapidity(TH2D* h2d,
                             const char* outName = "hMassSigmaVsRapidity",
                             int minEntries = 20,
                             std::vector<TH1D*>* outProjections = nullptr)
{
  if (!h2d) return nullptr;

  int nYBins = h2d->GetYaxis()->GetNbins();
  TH1D* hSigma = new TH1D(outName,
    ";Rapidity;#sigma_{M} (GeV/c^{2})",
    nYBins, h2d->GetYaxis()->GetXmin(), h2d->GetYaxis()->GetXmax());
  hSigma->SetStats(0);
  for (int iy = 1; iy <= nYBins; iy++) {
    TH1D* hproj = h2d->ProjectionX(Form("%s_projYbin%d", outName, iy), iy, iy);
    if (!hproj) continue;
    if (hproj->GetEntries() < minEntries || hproj->GetMaximum() <= 0.) {
      if (outProjections) {
        TH1D* hProjClone = static_cast<TH1D*>(hproj->Clone(Form("%s_saved_projYbin%d", outName, iy)));
        hProjClone->SetDirectory(nullptr);
        outProjections->push_back(hProjClone);
      }
      delete hproj;
      continue;
    }

    const double peak = hproj->GetMean();
    double rms = hproj->GetRMS();

    double halfWindow = (rms > 0.) ? 1.5 * rms : 0.03;
    if (halfWindow < 0.01) halfWindow = 0.01;
    double fitMin = TMath::Max(peak - halfWindow, hproj->GetXaxis()->GetXmin());
    double fitMax = TMath::Min(peak + halfWindow, hproj->GetXaxis()->GetXmax());

    TF1 fitFunc(Form("%s_fitYbin%d", outName, iy), "gaus", fitMin, fitMax);
    fitFunc.SetParameters(hproj->GetMaximum(), peak, (rms > 0.) ? rms : 0.02);
    if (hproj->Fit(&fitFunc, "QR") == 0) {
      hSigma->SetBinContent(iy, TMath::Abs(fitFunc.GetParameter(2)));
      hSigma->SetBinError(iy, fitFunc.GetParError(2));
    }
    if (outProjections) {
      TH1D* hProjClone = static_cast<TH1D*>(hproj->Clone(Form("%s_saved_projYbin%d", outName, iy)));
      hProjClone->SetDirectory(nullptr);
      outProjections->push_back(hProjClone);
    }
    delete hproj;
  }
  return hSigma;
}

void plotMTDimuons(int pdg = 443, const char* dirSimu = ".", int NRecoClusters = 9)
{
  // --- Initialization & Field Loading ---
  auto& param = NA6PLayoutParam::Instance();
  param.updateFromFile(Form("%s/na6pLayout.ini", dirSimu), "", true);

  auto magField = new MagneticField();
  magField->loadField();
  magField->setAsGlobalField();

  NA6PFastTrackFitter fitter;
  fitter.loadGeometry();
  fitter.enableMaterialCorrections();

  // --- File Loading ---
  TFile* ft = TFile::Open(Form("%s/TracksMatching.root", dirSimu));
  if (!ft || ft->IsZombie()) return;
  TTree* trTree = (TTree*)ft->Get("tracksMatching");
  std::vector<NA6PTrack>* trArr = nullptr;
  trTree->SetBranchAddress("Matching", &trArr);

  TFile* fk = TFile::Open(Form("%s/MCKine.root", dirSimu));
  if (!fk || fk->IsZombie()) return;
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);
  float mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  float rangeMMin = mass - 0.15*3;
  float rangeMMax = mass + 0.15*3;
  float rangeYMin = 2.1;
  float rangeYMax = 4.2;
  if (pdg == 443) { // J/psi
    rangeMMin = mass - 1;
    rangeMMax = mass + 1;
    rangeYMin = 1.9;
    rangeYMax = 3.0;
  }
  // --- Histogram Definitions (Restored) ---
  TH1D* hdimuy = new TH1D("hdimuy", "dimuon y", 100, 0., 5.);
  TH1D* hdimupt = new TH1D("hdimupt", "dimuon pT", 100, 0., 5.);
  TH1D* hdimumass = new TH1D("hdimumass", "dimuon m", 200, rangeMMin, rangeMMax);
  TH1D* hdimuygen = new TH1D("hdimuygen", "dimuon y gen", 100, 0., 5.);
  TH1D* hdimuptgen = new TH1D("hdimuptgen", "dimuon pT gen", 100, 0., 5.);
  TH1D* hdimumassgen = new TH1D("hdimumassgen", "dimuon m gen", 100, 0., 5.);
  TH2D* hdimuyrecvsgen = new TH2D("hdimuyrecvsgen", "dimuon y rec vs gen", 100, 1., 5., 100, 1., 5.);
  TH2D* hdimuptrecvsgen = new TH2D("hdimuptrecvsgen", "dimuon pt rec vs gen", 100, 0., 5., 100, 0., 5.);
  TH1D* hdimudeltayrecvsgen = new TH1D("hdimudeltayrecvsgen", "dimuon y rec - gen", 100, -1., 1.);
  TH1D* hdimudeltaptrecvsgen = new TH1D("hdimudeltaptrecvsgen", "dimuon pt rec - gen", 100, -6., 6.);
  TH2D* hdimurecmassvspt = new TH2D("hdimurecmassvspt", "dimuon rec mass vs pt", 200, 0., 4., 100, 0., 5.);
  TH2D* hdimurecmassvsxovx0 = new TH2D("hdimurecmassvsxovx0", "dimuon rec mass vs x/X0", 200, 0., 8., 100, 0., 100.);


  int nbinsM = 100;
  int nbinsY = 20;
  TH2D* hdimurecmassvseta = new TH2D("hdimurecmassvseta", "dimuon rec mass vs eta", nbinsM, rangeMMin, rangeMMax, nbinsY, rangeYMin, rangeYMax);
  TH2D* hdimurecptvseta = new TH2D("hdimurecptvseta", "dimuon rec pt vs eta", 100, 0., 5., 100, 1., 5.);
  TH2D* hdimurecpvseta = new TH2D("hdimurecpvseta", "dimuon rec p vs eta", 100, 0., 30., 100, 1., 5.);
  TH1D* hxovx0 = new TH1D("hxovx0", "x/X0 of dimuon tracks", 200, 0., 2.);
  TH1D* hzstep = new TH1D("hzstep", "z of dimuon tracks", 950, -50., 900.);
  TH2D* hxovx0VsZ = new TH2D("hxovx0VsZ", "x/X0 per step (1 cm) vs z; x/X0; z (cm)", 200, 0., 2., 900, -50., 850.);
  TH2D* hAVsZ = new TH2D("hAVsZ", "A mean per step (1 cm) vs z; A; z (cm)", 200, 0., 200., 900, -50., 850.);
  TH2D* hDensityVsZ = new TH2D("hDensityVsZ", "Density mean per step (1 cm) vs z; density (g/cm3); z (cm)", 200, 0., 10., 900, -50., 850.);
  TH1D* hRapidity = new TH1D("hRapidity", "dimuon rapidity", 100, -2, 7.);
  TH1D* hRapidityNeg = new TH1D("hRapidityNeg", "dimuon rapidity", 100, -2, 7.);
  TH1D* hRapidityPos = new TH1D("hRapidityPos", "dimuon rapidity", 100, -2, 7.);
  TH1D* hRapidityMC = new TH1D("hRapidityMC", "dimuon rapidity", 100, -2, 7.);
  int nEv = trTree->GetEntries();
  printf("Number of events = %d\n", nEv);

  for (int jEv = 0; jEv < nEv; jEv++) {
    mcTree->GetEvent(jEv);
    trTree->GetEvent(jEv);

    int nPart = (mcArr) ? mcArr->size() : 0;
    int nTracks = (trArr) ? trArr->size() : 0;
    printf("\nEvent %d particles = %d tracks = %d\n", jEv, nPart, nTracks);

    double zOrig = 0.;
    for (int jp = 0; jp < nPart; jp++) {
      auto& curPart = mcArr->at(jp);
      if (curPart.IsPrimary()) {
        zOrig = curPart.Vz();
        printf("zOrigin from primary vertex: %.2f\n", zOrig);
        break;
      }
    }

    //loop over the MCkine to fille the rapidity of muons
    for (int jp = 0; jp < nPart; jp++) {
      auto& curPart = mcArr->at(jp);
      if (TMath::Abs(curPart.GetPdgCode()) == 13) {
        double px = curPart.Px();
        double py = curPart.Py();
        double pz = curPart.Pz();
        double p = TMath::Sqrt(px*px + py*py + pz*pz);
        double rapidity = 0.5 * TMath::Log((p + pz) / (p - pz));
        hRapidityMC->Fill(rapidity);
      }
    }

    // Pre-filter tracks based on ClusterMap
    std::vector<int> filteredIndices;
    for (int jTr = 0; jTr < nTracks; ++jTr) {
      NA6PTrack& tr = trArr->at(jTr);
      double pxyz[3];
      tr.getPXYZ(pxyz);
      double px = pxyz[0];
      double py = pxyz[1];
      double pz = pxyz[2];
      double p = TMath::Sqrt(px*px + py*py + pz*pz);
      double rapidity = 0.5 * TMath::Log((p + pz) / (p - pz));
      int nClusters = tr.getNHits();

      if (nClusters < NRecoClusters) continue;
      if (TMath::Abs(tr.getParticleID()) > nPart) continue;
      if (!isMuonFrom(tr, *mcArr, pdg)) continue;
      // check if the track is a muon 
      hRapidity->Fill(rapidity);
      if (tr.getCharge() < 0) hRapidityNeg->Fill(rapidity);
      if (tr.getCharge() > 0) hRapidityPos->Fill(rapidity);
      uint32_t clustermap = tr.getClusterMap();
      filteredIndices.push_back(jTr);
    }

    // --- Dimuon Pairing ---
    if (filteredIndices.size() >= 2) {
      for (size_t i = 0; i < filteredIndices.size(); ++i) {
        NA6PTrack tr1 = trArr->at(filteredIndices[i]);
        double RecMom1[3], RecXYZ1[3], RecP1;

        fitter.propagateToZ(&tr1, zOrig);
          tr1.getPXYZ(RecMom1);
          RecP1 = tr1.getP();
        

        TLorentzVector decprod0;
        decprod0.SetPxPyPzE(RecMom1[0], RecMom1[1], RecMom1[2], TMath::Sqrt(RecP1 * RecP1 + 0.10566 * 0.10566));

        // MC Truth safety for Track 1
        int mcLabel1 = TMath::Abs(tr1.getParticleID());
        bool hasMC1 = (mcLabel1 >= 0 && mcLabel1 < nPart);
        TLorentzVector decprodgen0;
        if (hasMC1) {
          TParticle& p1 = mcArr->at(mcLabel1);
          decprodgen0.SetPxPyPzE(p1.Px(), p1.Py(), p1.Pz(), p1.Energy());
        }

        for (size_t j = i + 1; j < filteredIndices.size(); ++j) {
          NA6PTrack tr2 = trArr->at(filteredIndices[j]);
          double RecMom2[3], RecXYZ2[3], RecP2;

          fitter.propagateToZ(&tr2, zOrig);
          tr2.getPXYZ(RecMom2);
          RecP2 = tr2.getP();
          RecP2 = tr2.getP();
          

          TLorentzVector decprod1;
          decprod1.SetPxPyPzE(RecMom2[0], RecMom2[1], RecMom2[2], TMath::Sqrt(RecP2 * RecP2 + 0.10566 * 0.10566));

          TLorentzVector moth = decprod0 + decprod1;
          hdimuy->Fill(moth.Rapidity());
          hdimupt->Fill(moth.Pt());
          hdimumass->Fill(moth.M());
          hdimurecmassvspt->Fill(moth.M(), moth.Pt());
          hdimurecmassvseta->Fill(moth.M(), moth.Rapidity());
          hdimurecptvseta->Fill(moth.Pt(), moth.Rapidity());
          hdimurecpvseta->Fill(moth.P(), moth.Rapidity());

          // Fill Gen info if both have MC truth
          int mcLabel2 = TMath::Abs(tr2.getParticleID());
          if (hasMC1 && (mcLabel2 >= 0 && mcLabel2 < nPart)) {
            TParticle& p2 = mcArr->at(mcLabel2);
            TLorentzVector decprodgen1;
            decprodgen1.SetPxPyPzE(p2.Px(), p2.Py(), p2.Pz(), p2.Energy());
            TLorentzVector mothgen = decprodgen0 + decprodgen1;
            
            hdimumassgen->Fill(mothgen.M());
            hdimuyrecvsgen->Fill(moth.Rapidity(), mothgen.Rapidity());
            hdimuptrecvsgen->Fill(moth.Pt(), mothgen.Pt());
            hdimudeltayrecvsgen->Fill(moth.Rapidity() - mothgen.Rapidity());
            hdimudeltaptrecvsgen->Fill(moth.Pt() - mothgen.Pt());
          }
        }
      }
    }
  }

  // --- Plotting & Saving (All Canvases Restored) ---
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  TCanvas* cDimu = new TCanvas("cDimu", "Dimuon mass", 1200, 1200);
  hdimumass->SetMarkerColor(kAzure + 8);
  hdimumass->SetLineColor(kBlue + 1);
  hdimumass->SetMarkerStyle(20);
  hdimumass->SetTitle(Form(";m_{#mu#mu} (GeV/c^{2});Counts / %g MeV/c^{2}", 1000*hdimumass->GetXaxis()->GetBinWidth(2)));
  // get mass from pdgcode 
  hdimumass->Fit("gaus", "", "", mass - 0.3, mass + 0.3);
  hdimumass->Draw("EP");
  cDimu->SaveAs(Form("%s/dimuon_mass_spectrumMT.png", dirSimu));
  
  TFile fOut(Form("%s/dimuon_analysisVT.root", dirSimu), "RECREATE");
  hdimumass->Write();
  hRapidity->Write();
  hRapidityNeg->Write();
  hRapidityPos->Write();
  hRapidityMC->Write();
  TH1F * hEffVsRapidity = (TH1F*)hRapidity->Clone("hEffVsRapidity");
  hEffVsRapidity->Divide(hRapidityMC);
  hEffVsRapidity->SetTitle("Dimuon Reconstruction Efficiency vs Rapidity;Rapidity;Efficiency");
  hEffVsRapidity->Write();
  TH1F * hEffVsRapidityNeg = (TH1F*)hRapidityNeg->Clone("hEffVsRapidityNeg");
  hEffVsRapidityNeg->Divide(hRapidityMC);
  hEffVsRapidityNeg->SetTitle("Dimuon Reconstruction Efficiency vs Rapidity (Negative);Rapidity;Efficiency");
  hEffVsRapidityNeg->Write();
  TH1F * hEffVsRapidityPos = (TH1F*)hRapidityPos->Clone("hEffVsRapidityPos");
  hEffVsRapidityPos->Divide(hRapidityMC);
  hEffVsRapidityPos->SetTitle("Dimuon Reconstruction Efficiency vs Rapidity (Positive);Rapidity;Efficiency");
  hEffVsRapidityPos->Write();

  std::vector<TH1D*> massProjectionsVsY;
  TH1D* hMassSigmaVsRapidity = getMassSigmaVsRapidity(hdimurecmassvseta, "hMassSigmaVsRapidity", 20, &massProjectionsVsY);
  hMassSigmaVsRapidity->SetMarkerStyle(20);
  hMassSigmaVsRapidity->SetMarkerColor(kRed + 1);
  hMassSigmaVsRapidity->SetLineColor(kRed + 1);
  TCanvas* cSigmaRapidity = new TCanvas("cSigmaRapidity", "Mass Sigma vs Rapidity", 1200, 1200);
  cSigmaRapidity->SetLeftMargin(0.14);
  hMassSigmaVsRapidity->Draw("EP");
  hMassSigmaVsRapidity->GetYaxis()->SetRangeUser(0., 0.1);
  cSigmaRapidity->SaveAs(Form("%s/dimuon_mass_sigma_vs_rapidity.png", dirSimu));
  hMassSigmaVsRapidity->Write();
  if (!massProjectionsVsY.empty()) {
    TDirectory* dProj = fOut.mkdir("MassProjectionsVsRapidity");
    if (dProj) {
      dProj->cd();
      for (auto* hProj : massProjectionsVsY) {
        if (hProj) hProj->Write();
      }
      fOut.cd();
    }
  }

  fOut.Close();
  for (auto* hProj : massProjectionsVsY) {
    delete hProj;
  }
}