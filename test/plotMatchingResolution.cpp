#include <vector>
#include <string>
#include <iostream>

void plotMatchingResolution()
{
  // === Lista dei valori N ===
  std::vector<double> Nvalues = {38.1175, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300};

  // Numero di istogrammi residuals
  const int nResiduals = 5; // hResiduals0 ... hResiduals4

  // list of the track parameter names
  //   external param0:   local Y-coordinate of a track (cm)
  //   external param1:   local Z-coordinate of a track (cm)
  //   external param2:   local sine of the track momentum azimuthal angle
  //   external param3:   tangent of the track momentum dip angle
  //   external param4:   1/pt (1/(GeV/c))
  const char* paramNames[nResiduals] = {
    "local Y-coordinate (cm)",
    "local Z-coordinate (cm)",
    "sine of azimuthal angle",
    "tangent of dip angle",
    "1/p_{T} (1/(GeV/c))"};

  for (int ires = 0; ires < nResiduals; ++ires) {

    std::string histName = Form("hResidualsTrue%d", ires);

    std::vector<double> N;
    std::vector<double> sigma;
    std::vector<double> sigmaErr;
    std::vector<TH1D*> histograms;

    // ======================================================
    // Overlay istogrammi
    // ======================================================
    TCanvas* c1 = new TCanvas(
      Form("c_res_%d", ires),
      Form("Residuals %d overlay", ires),
      1400, 900);

    bool first = true;
    int color = 1;

    for (double n : Nvalues) {

      std::string fname = Form("zprop_%g/QAMatchingMatching.root", n);
      TFile* f = TFile::Open(fname.c_str());

      if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << fname << std::endl;
        continue;
      }

      TH1D* h = (TH1D*)f->Get(histName.c_str());
      if (!h) {
        std::cerr << "Histogram " << histName
                  << " not found in " << fname << std::endl;
        f->Close();
        continue;
      }

      TH1D* hClone = (TH1D*)h->Clone(
        Form("h_%s_%d", histName.c_str(), (int)n));
      hClone->SetDirectory(nullptr); // <<< FIX CHIAVE
      histograms.push_back(hClone);

      hClone->SetLineColor(color);
      hClone->SetLineWidth(2);
      hClone->SetStats(false);

      c1->cd();

      if (first) {
        hClone->Draw("hist");
        hClone->GetXaxis()->SetTitle("Residual");
        hClone->GetYaxis()->SetTitle("Entries");
        first = false;
      } else {
        hClone->Draw("hist same");
      }

      // === Fit gaussiano ===
      TF1* g = new TF1(
        Form("g_%s_%d", histName.c_str(), (int)n),
        "gaus",
        hClone->GetMean() - 2 * hClone->GetRMS(),
        hClone->GetMean() + 2 * hClone->GetRMS());

      hClone->Fit(g, "Q0");

      N.push_back(n);
      sigma.push_back(g->GetParameter(2));
      sigmaErr.push_back(g->GetParError(2));

      color++;
      f->Close();
    }

    c1->SaveAs(Form("ResidualsTrueOverlay_%s.png", histName.c_str()));
    delete c1;

    // ======================================================
    // Risoluzione vs N
    // ======================================================
    TCanvas* c2 = new TCanvas(
      Form("c_res_vs_Z_%d", ires),
      Form("Resolution vs z_{match} (cm) (%s)", histName.c_str()),
      1600, 1400);
    c2->SetGrid();

    TGraphErrors* gr = new TGraphErrors(
      N.size(), N.data(), sigma.data(), nullptr, sigmaErr.data());

    gr->SetTitle(Form("%s vs z_{match} (cm)", paramNames[ires]));
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.2);
    gr->SetLineWidth(2);

    gr->GetXaxis()->SetTitle("z_{prop} (cm)");
    gr->GetYaxis()->SetTitle("#sigma_{residual}");
    // set left margin
    c2->SetLeftMargin(0.15);
    c2->SetBottomMargin(0.15);
    gr->GetYaxis()->SetTitleOffset(1.5);
    // set axis title size
    gr->GetXaxis()->SetTitleSize(0.05);
    gr->GetYaxis()->SetTitleSize(0.05);

    gr->Draw("APL");
    c2->SaveAs(Form("ResolutionVsN_%s.png", histName.c_str()));
    delete c2;
  }

  const std::string chi2HistNameList[3] = {"hChi2MatchTrue", "hChi2MatchFake", "hChi2MatchSemiFake"};
  for (const auto& chi2HistName : chi2HistNameList) {
  TCanvas* cChi2 = new TCanvas(
    "c_chi2",
    "Chi2 match overlay",
    1400, 900);
  TLegend* leg = new TLegend(0.35, 0.45, 0.8, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);
  bool first = true;
  int color = 1;

  std::vector<TH1D*> chi2Hists;

  for (double n : Nvalues) {

    std::string fname = Form("zprop_%g/QAMatchingMatching.root", n);
    TFile* f = TFile::Open(fname.c_str());

    if (!f || f->IsZombie()) {
      std::cerr << "Cannot open " << fname << std::endl;
      continue;
    }

    TH1D* h = (TH1D*)f->Get(chi2HistName.c_str());
    if (!h) {
      std::cerr << "Histogram " << chi2HistName
                << " not found in " << fname << std::endl;
      f->Close();
      continue;
    }

    TH1D* hClone = (TH1D*)h->Clone(
      Form("h_%s_%d", chi2HistName.c_str(), (int)n));
    hClone->SetDirectory(nullptr); // FONDAMENTALE
    chi2Hists.push_back(hClone);

    hClone->SetLineColor(color);
    hClone->SetLineWidth(2);
    hClone->SetStats(false);

    cChi2->cd();

    if (first) {
      hClone->Draw("hist");
      hClone->GetXaxis()->SetTitle("#chi^{2}_{match}/nDoF");
      hClone->GetYaxis()->SetTitle("Entries");
      first = false;
    } else {
      hClone->Draw("hist same");
    }
    leg->AddEntry(hClone, Form("Z_{match} = %.0f, mean = %.2f, std dev = %.2f", n, hClone->GetMean(), hClone->GetStdDev()), "l");

    color++;
    f->Close();
  }

  cChi2->Modified();
  cChi2->Update();
  leg->Draw();
  cChi2->SaveAs(Form("Chi2MatchTrueOverlay_%s.png", chi2HistName.c_str()));
  delete cChi2;
}

  // Loop to plot the counts of the TH2D vs zprop
  const std::string countsHistList[3] = {"hTrue", "hFake", "hSemiFake"};
  for (const auto& countsHistName : countsHistList) {
    
    std::vector<double> zprop_vals;
    std::vector<double> counts;
    std::vector<double> countsErr;
    
    for (double n : Nvalues) {
      std::string fname = Form("zprop_%g/QAMatchingMatching.root", n);
      TFile* f = TFile::Open(fname.c_str());
      
      if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << fname << std::endl;
        continue;
      }
      
      TH2D* h = (TH2D*)f->Get(countsHistName.c_str());
      if (!h) {
        std::cerr << "Histogram " << countsHistName
                  << " not found in " << fname << std::endl;
        f->Close();
        continue;
      }
      
      double totalCounts = h->GetEntries();
      double totalCountsErr = 0; //std::sqrt(totalCounts);
      
      zprop_vals.push_back(n);
      counts.push_back(totalCounts);
      countsErr.push_back(totalCountsErr);
      
      f->Close();
    }
    
    // Plot counts vs zprop
    TCanvas* cCounts = new TCanvas(
      Form("c_counts_%s", countsHistName.c_str()),
      Form("Counts vs z_{prop} (%s)", countsHistName.c_str()),
      1600, 1400);
    cCounts->SetGrid();
    
    TGraphErrors* grCounts = new TGraphErrors(
      zprop_vals.size(), zprop_vals.data(), counts.data(), 
      nullptr, countsErr.data());
    
    grCounts->SetTitle(Form("%s counts vs z_{prop}", countsHistName.c_str()));
    grCounts->SetMarkerStyle(20);
    grCounts->SetMarkerSize(1.2);
    grCounts->SetLineWidth(2);
    
    grCounts->GetXaxis()->SetTitle("z_{prop} (cm)");
    grCounts->GetYaxis()->SetTitle("Number of tracks");
    cCounts->SetLeftMargin(0.15);
    cCounts->SetBottomMargin(0.15);
    grCounts->GetYaxis()->SetTitleOffset(1.5);
    grCounts->GetXaxis()->SetTitleSize(0.05);
    grCounts->GetYaxis()->SetTitleSize(0.05);
    
    grCounts->Draw("APL");
    cCounts->SaveAs(Form("CountsVsZprop_%s.png", countsHistName.c_str()));
    delete cCounts;
  }
}
