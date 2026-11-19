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
#include "NA6PVerTelHit.h"
#include "NA6PLayoutParam.h"
#include "NA6PFastTrackFitter.h"
#endif

bool isFromD0(const TParticle part, const std::vector<TParticle>& mcArr) {
  if (TMath::Abs(part.GetPdgCode()) != 211) return false; // Not a pion
  int motherIndex = part.GetFirstMother();
  if (motherIndex < 0 || motherIndex >= static_cast<int>(mcArr.size())) return false; // Invalid mother index
  const TParticle& mother = mcArr[motherIndex];
  return TMath::Abs(mother.GetPdgCode()) == 421; // Check if mother is D0
}

bool isPiFromD0(const NA6PTrack& trk, const std::vector<TParticle>& mcArr) {
  int pid = trk.getParticleID();
  auto part = mcArr.at(pid);
  if (TMath::Abs(part.GetPdgCode()) != 211) return false; // Not a pion
  int motherIndex = part.GetFirstMother();
  if (motherIndex < 0 || motherIndex >= static_cast<int>(mcArr.size())) return false; // Invalid mother index
  const TParticle& mother = mcArr[motherIndex];
  return TMath::Abs(mother.GetPdgCode()) == 421; // Check if mother is D0
}
bool isKaonFromD0(const NA6PTrack& trk, const std::vector<TParticle>& mcArr) {
  int pid = trk.getParticleID();
  auto part = mcArr.at(pid);
  if (TMath::Abs(part.GetPdgCode()) != 321) return false; // Not a kaon
  int motherIndex = part.GetFirstMother();
  if (motherIndex < 0 || motherIndex >= static_cast<int>(mcArr.size())) return false; // Invalid mother index
  const TParticle& mother = mcArr[motherIndex];
  return TMath::Abs(mother.GetPdgCode()) == 421; // Check if mother is D0
}

TH1D* getMassSigmaVsRapidity(TH2D* hMassVsRapidity,
                             const char* outName = "hMassSigmaVsRapidity",
                             int minEntriesPerBin = 20,
                             std::vector<TH1D*>* outProjections = nullptr)
{
  if (!hMassVsRapidity) {
    return nullptr;
  }

  const int nRapidityBins = hMassVsRapidity->GetYaxis()->GetNbins();
  const double rapidityMin = hMassVsRapidity->GetYaxis()->GetXmin();
  const double rapidityMax = hMassVsRapidity->GetYaxis()->GetXmax();

  TH1D* hSigmaVsRapidity = new TH1D(outName,
                                    "Mass #sigma from gaussian fit vs rapidity;Rapidity;#sigma_{M} (GeV/c^{2})",
                                    nRapidityBins, rapidityMin, rapidityMax);

  for (int iRap = 1; iRap <= nRapidityBins; ++iRap) {
    TH1D* hMassProj = hMassVsRapidity->ProjectionX(Form("%s_projYbin%d", outName, iRap), iRap, iRap);
    if (!hMassProj) {
      continue;
    }

    if (hMassProj->GetEntries() < minEntriesPerBin || hMassProj->GetMaximum() <= 0.) {
      delete hMassProj;
      continue;
    }

    const int maxBin = hMassProj->GetMaximumBin();
    const double peak = hMassProj->GetMean();
    const double rms = hMassProj->GetRMS();

    double halfWindow = (rms > 0.) ? 1.5 * rms : 0.03;
    if (halfWindow < 0.01) {
      halfWindow = 0.01;
    }

    double fitMin = peak - halfWindow;
    double fitMax = peak + halfWindow;
    const double xMin = hMassProj->GetXaxis()->GetXmin();
    const double xMax = hMassProj->GetXaxis()->GetXmax();
    if (fitMin < xMin) {
      fitMin = xMin;
    }
    if (fitMax > xMax) {
      fitMax = xMax;
    }

    TF1 fitFunc(Form("%s_fitYbin%d", outName, iRap), "gaus", fitMin, fitMax);
    fitFunc.SetParameters(hMassProj->GetMaximum(), peak, (rms > 0.) ? rms : 0.02);

    int fitStatus = hMassProj->Fit(&fitFunc, "QR");
    if (fitStatus == 0) {
      hSigmaVsRapidity->SetBinContent(iRap, TMath::Abs(fitFunc.GetParameter(2)));
      hSigmaVsRapidity->SetBinError(iRap, fitFunc.GetParError(2));
      if (outProjections) {
        TH1D* hProjClone = static_cast<TH1D*>(hMassProj->Clone(Form("%s_saved_projYbin%d", outName, iRap)));
        hProjClone->SetDirectory(nullptr);
        outProjections->push_back(hProjClone);
      }
    }

    delete hMassProj;
  }

  return hSigmaVsRapidity;
}

void plotVTD0(const char* dirSimu = ".", int NRecoClusters = 4)
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
  TFile* ft = TFile::Open(Form("%s/TracksVerTel.root", dirSimu));
  if (!ft || ft->IsZombie()) return;
  TTree* trTree = (TTree*)ft->Get("tracksVerTel");
  std::vector<NA6PTrack>* trArr = nullptr;
  trTree->SetBranchAddress("VerTel", &trArr);

  TFile* fk = TFile::Open(Form("%s/MCKine.root", dirSimu));
  if (!fk || fk->IsZombie()) return;
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  // Note: Hits tree is opened but not used in pairing logic in original script
  TFile* fh = TFile::Open(Form("%s/HitsVerTel.root", dirSimu));
  TTree* th = (fh) ? (TTree*)fh->Get("hitsVerTel") : nullptr;

  float d0mass = 1.86483; // D0 mass in GeV/c^2
  int nbinsM = 100;
  int nbinsY = 20;
  // --- Histogram Definitions (Restored) ---
  TH1D* hdimuy = new TH1D("hdimuy", "dimuon y", 100, 0., 5.);
  TH1D* hdimupt = new TH1D("hdimupt", "dimuon pT", 100, 0., 5.);
  TH1D* hdimumass = new TH1D("hdimumass", "dimuon m", nbinsM, d0mass - 0.15, d0mass + 0.15);
  TH1D* hdimuygen = new TH1D("hdimuygen", "dimuon y gen", 100, 0., 5.);
  TH1D* hdimuptgen = new TH1D("hdimuptgen", "dimuon pT gen", 100, 0., 5.);
  TH1D* hdimumassgen = new TH1D("hdimumassgen", "dimuon m gen", 100, 0., 5.);
  TH2D* hdimuyrecvsgen = new TH2D("hdimuyrecvsgen", "dimuon y rec vs gen", 100, 1., 5., 100, 1., 5.);
  TH1D* hdimuyrec = new TH1D("hdimuyrec", "dimuon y rec", 100, 0., 5.);
  TH2D* hdimuptrecvsgen = new TH2D("hdimuptrecvsgen", "dimuon pt rec vs gen", 100, 0., 5., 100, 0., 5.);
  TH1D* hdimudeltayrecvsgen = new TH1D("hdimudeltayrecvsgen", "dimuon y rec - gen", 100, -1., 1.);
  TH1D* hdimudeltaptrecvsgen = new TH1D("hdimudeltaptrecvsgen", "dimuon pt rec - gen", 100, -6., 6.);
  TH2D* hdimurecmassvspt = new TH2D("hdimurecmassvspt", "dimuon rec mass vs pt", nbinsM, d0mass - 0.15, d0mass + 0.15, 100, 0., 5.);
  TH2D* hdimurecmassvsxovx0 = new TH2D("hdimurecmassvsxovx0", "dimuon rec mass vs x/X0", 200, 0., 8., 100, 0., 100.);
  TH2D* hdimurecmassvsy = new TH2D("hdimurecmassvsy", "dimuon rec mass vs y", nbinsM, d0mass - 0.15, d0mass + 0.15, nbinsY, 0, 5.);
  TH2D* hdimurecptvsy = new TH2D("hdimurecptvsy", "dimuon rec pt vs eta", 100, 0., 5., 100, 1., 5.);
  TH2D* hdimurecpvsy = new TH2D("hdimurecpvsy", "dimuon rec p vs eta", 100, 0., 30., 100, 1., 5.);
  TH1D* hRapidity = new TH1D("hRapidity", "dimuon rapidity", 100, -2, 7.);
  TH1D* hRapidityNeg = new TH1D("hRapidityNeg", "dimuon rapidity", 100, -2, 7.);
  TH1D* hRapidityPos = new TH1D("hRapidityPos", "dimuon rapidity", 100, -2, 7.);
  TH1D* hRapidityMC = new TH1D("hRapidityMC", "dimuon rapidity MC", 100, -2, 7.);
  // Efficiency: reco D0 tracks / generated D0 tracks vs momentum
  TH1D* hPTrue = new TH1D("hPTrue", "D0 reco momentum (true)", 50, 0., 30.);
  TH1D* hPMC   = new TH1D("hPMC",  "D0 generated momentum",  50, 0., 30.);
  int nEv = trTree->GetEntries();
  printf("Number of events = %d\n", nEv);

  for (int jEv = 0; jEv < nEv; jEv++) {
    mcTree->GetEvent(jEv);
    trTree->GetEvent(jEv);
    if(th) th->GetEvent(jEv);

    int nPart = (mcArr) ? mcArr->size() : 0;
    int nTracks = (trArr) ? trArr->size() : 0;
    printf("\nEvent %d particles = %d tracks = %d\n", jEv, nPart, nTracks);

    double zOrig = 0.;
    for (int jp = 0; jp < nPart; jp++) {
      auto& curPart = mcArr->at(jp);
      if (isFromD0(curPart, *mcArr)) {
        zOrig = curPart.Vz();
        printf("zOrigin from D0 vertex: %.2f\n", zOrig);
        break;
      }
    }

    //loop over the MCkine to fille the rapidity of muons
    for (int jp = 0; jp < nPart; jp++) {
      auto& curPart = mcArr->at(jp);
      if (TMath::Abs(curPart.GetPdgCode()) == 421) {
        hRapidityMC->Fill(curPart.Y());
        hPMC->Fill(curPart.P());
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
      if (tr.getParticleID() < 0 || TMath::Abs(tr.getParticleID()) > nPart) continue; // check if the track is a muon
      // check if the track is a muon 
      hRapidity->Fill(rapidity);
      if (tr.getCharge() < 0) hRapidityNeg->Fill(rapidity);
      if (tr.getCharge() > 0) hRapidityPos->Fill(rapidity);
      // Fill true D0 track momentum for efficiency numerator
      int pid = tr.getParticleID();
      if (pid >= 0 && pid < nPart) {
        auto& mcPart = mcArr->at(pid);
        int motherIdx = mcPart.GetFirstMother();
        if (motherIdx >= 0 && motherIdx < nPart &&
            TMath::Abs(mcArr->at(motherIdx).GetPdgCode()) == 421) {
          hPTrue->Fill(p);
        }
      }
      uint32_t clustermap = tr.getClusterMap();
      filteredIndices.push_back(jTr);
    }

    // --- Dimuon Pairing ---
    if (filteredIndices.size() >= 2) {
      for (size_t i = 0; i < filteredIndices.size(); ++i) {
        NA6PTrack tr1 = trArr->at(filteredIndices[i]);
        if(!isPiFromD0(tr1, *mcArr)) continue;
        double RecMom1[3], RecXYZ1[3], RecP1;

        fitter.propagateToZ(&tr1, zOrig);
          tr1.getPXYZ(RecMom1);
          RecP1 = tr1.getP();
        
        float mpion = 0.10566; // GeV/c^2
        float mkaon = 0.49367; // GeV/c^2

        // MC Truth safety for Track 1
        int mcLabel1 = tr1.getParticleID();
        bool hasMC1 = (mcLabel1 >= 0 && mcLabel1 < nPart);
        TLorentzVector decprodgen0;
        if (hasMC1) {
          TParticle& p1 = mcArr->at(mcLabel1);
          decprodgen0.SetPxPyPzE(p1.Px(), p1.Py(), p1.Pz(), p1.Energy());
        }

        for (size_t j = i + 1; j < filteredIndices.size(); ++j) {
          NA6PTrack tr2 = trArr->at(filteredIndices[j]);
          if(!isKaonFromD0(tr2, *mcArr)) continue;
          double RecMom2[3], RecXYZ2[3], RecP2;

          fitter.propagateToZ(&tr2, zOrig);
          tr2.getPXYZ(RecMom2);
          RecP2 = tr2.getP();
          RecP2 = tr2.getP();
          

          TLorentzVector decprod0;
          decprod0.SetPxPyPzE(RecMom1[0], RecMom1[1], RecMom1[2], TMath::Sqrt(RecP1 * RecP1 + mpion * mpion));
          TLorentzVector decprod1;
          decprod1.SetPxPyPzE(RecMom2[0], RecMom2[1], RecMom2[2], TMath::Sqrt(RecP2 * RecP2 + mkaon * mkaon));

          TLorentzVector moth = decprod0 + decprod1;

          hdimuy->Fill(moth.Rapidity());
          hdimupt->Fill(moth.Pt());
          hdimumass->Fill(moth.M());
          hdimurecmassvspt->Fill(moth.M(), moth.Pt());
          hdimurecmassvsy->Fill(moth.M(), moth.Rapidity());
            /*
          decprod0.SetPxPyPzE(RecMom1[0], RecMom1[1], RecMom1[2], TMath::Sqrt(RecP1 * RecP1 + mkaon * mkaon));
          decprod1.SetPxPyPzE(RecMom2[0], RecMom2[1], RecMom2[2], TMath::Sqrt(RecP2 * RecP2 + mpion * mpion));
          moth = decprod0 + decprod1;
          

          hdimuy->Fill(moth.Rapidity());
          hdimupt->Fill(moth.Pt());
          hdimumass->Fill(moth.M());
          hdimurecmassvspt->Fill(moth.M(), moth.Pt());
          hdimurecmassvseta->Fill(moth.M(), moth.Rapidity());
          hdimurecptvseta->Fill(moth.Pt(), moth.Rapidity());
          hdimurecpvseta->Fill(moth.P(), moth.Rapidity());
*/
          // Fill Gen info if both have MC truth
          int mcLabel2 = tr2.getParticleID();
          if (hasMC1 && (mcLabel2 >= 0 && mcLabel2 < nPart)) {
            TParticle& p2 = mcArr->at(mcLabel2);
            TLorentzVector decprodgen1;
            decprodgen1.SetPxPyPzE(p2.Px(), p2.Py(), p2.Pz(), p2.Energy());
            TLorentzVector mothgen = decprodgen0 + decprodgen1;
            
            hdimumassgen->Fill(mothgen.M());
            hdimuyrecvsgen->Fill(moth.Rapidity(), mothgen.Rapidity());
            hdimuyrec->Fill(moth.Rapidity());
            hdimuptrecvsgen->Fill(moth.Pt(), mothgen.Pt());
            hdimudeltayrecvsgen->Fill(moth.Rapidity() - mothgen.Rapidity());
            hdimudeltaptrecvsgen->Fill(moth.Pt() - mothgen.Pt());
          }
        }
      }
    }
  }

  // --- Plotting & Saving (All Canvases Restored) ---
  gStyle->SetOptFit(1111);
  TCanvas* cDimu = new TCanvas("cDimu", "Dimuon mass", 1200, 1200);
  hdimumass->SetMarkerColor(kAzure + 8);
  hdimumass->SetLineColor(kBlue + 1);
  hdimumass->SetMarkerStyle(20);
  hdimumass->SetTitle("Dimuon Mass Spectrum;Mass (GeV/c^{2});Counts / 0.1 GeV/c^{2}");
  hdimumass->Fit("gaus", "", "", 2.85, 3.35);
  hdimumass->Draw("EP");
  cDimu->SaveAs(Form("%s/D0_mass.png", dirSimu));

  std::vector<TH1D*> massProjectionsVsY;
  TH1D* hMassSigmaVsRapidity = getMassSigmaVsRapidity(hdimurecmassvsy, "hMassSigmaVsRapidity", 20, &massProjectionsVsY);
  hMassSigmaVsRapidity->SetMarkerColor(kRed + 1);
  hMassSigmaVsRapidity->SetLineColor(kRed + 1);
  TCanvas* cSigmaRapidity = new TCanvas("cSigmaRapidity", "Mass Sigma vs Rapidity", 1200, 1200);
  hMassSigmaVsRapidity->Draw("EP");
  hMassSigmaVsRapidity->GetYaxis()->SetRangeUser(0., 0.1);
  cSigmaRapidity->SaveAs(Form("%s/D0_mass_sigma_vs_rapidity.png", dirSimu));
  // --- Efficiency vs momentum canvas ---
  TH1D* hEfficiency = (TH1D*)hPTrue->Clone("hEfficiency");
  hEfficiency->Divide(hPMC);
  hEfficiency->SetTitle("D^{0} reconstruction efficiency vs p;p (GeV/c);Efficiency");
  hEfficiency->SetMarkerColor(kBlue+1);
  hEfficiency->SetLineColor(kBlue+1);
  hEfficiency->SetMarkerStyle(kFullCircle);
  hEfficiency->SetMarkerSize(1.2);
  hEfficiency->GetYaxis()->SetRangeUser(0., 1.1);
  hEfficiency->SetStats(0);
  TCanvas* cEff = new TCanvas("cEff", "Efficiency", 1200, 1200);
  cEff->SetGrid();
  cEff->SetLeftMargin(0.14);
  cEff->SetRightMargin(0.05);
  cEff->SetBottomMargin(0.13);
  cEff->SetTopMargin(0.07);
  hEfficiency->Draw("E P");
  cEff->SaveAs(Form("%s/D0_efficiency_vs_p.png", dirSimu));

  TFile fOut(Form("%s/D0_analysis.root", dirSimu), "RECREATE");
  hdimumass->Write();
  hdimurecmassvsy->Write();
  if (hMassSigmaVsRapidity) hMassSigmaVsRapidity->Write();
  hRapidity->Write();
  hRapidityMC->Write();
  hEfficiency->Write();
  hdimuyrec->Write();
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