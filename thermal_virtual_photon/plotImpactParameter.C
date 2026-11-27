#include "NA6PMCEventHeader.h"
#include "NA6PMCComposedLabel.h"
#include "NA6PTrack.h"
#include "NA6PMatch.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
R__LOAD_LIBRARY(librecLib)

void plotImpactParameter(const char* runDir = ".", double maxImpactUm = 2000.)
{
  const std::string run = runDir;
  if (maxImpactUm <= 0.) throw std::runtime_error("Invalid impact parameter range");
  TFile vf((run + "/reco/TracksVerTel.root").c_str()), mf((run + "/reco/TracksMuonSpec.root").c_str());
  TFile mc((run + "/sim/MCKine.root").c_str()), matched((run + "/reco_5hits/TracksMatching.root").c_str());
  if (vf.IsZombie() || mf.IsZombie() || mc.IsZombie() || matched.IsZombie()) throw std::runtime_error("Cannot open input ROOT files");
  TTreeReader vr("tracksVerTel", &vf), mr("tracksMuonSpec", &mf), cr("mckine", &mc), tr("tracksMatching", &matched);
  TTreeReaderValue<std::vector<NA6PTrack>> vt(vr, "VerTel");
  TTreeReaderValue<std::vector<NA6PMCComposedLabel>> vl(vr, "VerTelMCTruth"), ml(mr, "MuonSpecMCTruth");
  TTreeReaderValue<NA6PMCEventHeader> vertex(cr, "header");
  TTreeReaderValue<std::vector<NA6PMatch>> pairs(tr, "Matching");
  TH1D correct("correct", ";Transverse distance at primary-vertex z (#mum);Matches / bin", 100, 0., maxImpactUm);
  TH1D wrong("wrong", ";Transverse distance at primary-vertex z (#mum);Matches / bin", 100, 0., maxImpactUm);
  correct.SetDirectory(nullptr);
  wrong.SetDirectory(nullptr);
  int events=0, offPlane=0;
  while (cr.Next()) {
    if (!vr.Next() || !mr.Next() || !tr.Next()) throw std::runtime_error("Mismatched tree entries");
    for (const auto& pair : *pairs) {
      const int vi=pair.getIndexVT(), mi=pair.getIndexMS();
      if (vi<0 || mi<0 || size_t(vi)>=vt->size() || size_t(vi)>=vl->size() || size_t(mi)>=ml->size())
        throw std::runtime_error("Invalid match index");
      if (vt->at(vi).getNHits()<5) throw std::runtime_error("Match below VT threshold");
      if (std::abs(pair.getZ()-vertex->getVZ())>1.e-3) {
        ++offPlane;
        continue;
      }
      const auto& v=vl->at(vi);
      const auto& m=ml->at(mi);
      const bool good=v.isValid() && !v.isFake() && m.isValid() && !m.isFake() && v==m;
      const double displacement=1.e4*std::hypot(pair.getX()-vertex->getVX(),pair.getY()-vertex->getVY());
      if (!std::isfinite(displacement)) throw std::runtime_error("Non-finite impact parameter");
      (good ? correct : wrong).Fill(displacement);
    }
    ++events;
  }
  if (vr.Next() || mr.Next() || tr.Next()) throw std::runtime_error("Mismatched tree entries");
  gSystem->mkdir((run+"/impact_parameter").c_str(),true);
  TFile output((run+"/impact_parameter/impact_parameter.root").c_str(),"RECREATE");
  if (output.IsZombie()) throw std::runtime_error("Cannot create output ROOT file");
  correct.Write();
  wrong.Write();
  gStyle->SetOptStat(0);
  TCanvas canvas("impactParameter", "Impact parameter, VT threshold 5", 900, 700);
  canvas.SetMargin(.13,.05,.13,.16);
  canvas.SetTicks(1,1);
  canvas.SetLogy();
  correct.SetLineColor(kAzure+2);
  correct.SetLineWidth(3);
  wrong.SetLineColor(kOrange+7);
  wrong.SetLineWidth(3);
  correct.SetMinimum(.5);
  correct.SetMaximum(2.5*std::max(correct.GetMaximum(),wrong.GetMaximum()));
  correct.Draw("HIST");
  wrong.Draw("HIST SAME");
  TLegend legend(.58,.69,.92,.82);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(.028);
  legend.AddEntry(&correct,Form("Correct match (%lld)",static_cast<Long64_t>(correct.GetEntries())),"l");
  legend.AddEntry(&wrong,Form("False match (%lld)",static_cast<Long64_t>(wrong.GetEntries())),"l");
  legend.Draw();
  TLatex label;
  label.SetNDC();
  label.SetTextFont(62);
  label.SetTextSize(.038);
  label.DrawLatex(.13,.95,"Matched-track impact parameter");
  label.SetTextFont(42);
  label.SetTextSize(.028);
  label.DrawLatex(.13,.90,"VT threshold #geq5  |  MC primary vertex  |  all matched tracks");
  canvas.SaveAs((run+"/impact_parameter/impactParameter.png").c_str());
  canvas.Write();
  std::cout<<"Events "<<events<<", correct "<<correct.GetEntries()<<", false "<<wrong.GetEntries()
           <<", overflow correct "<<correct.GetBinContent(correct.GetNbinsX()+1)
           <<", overflow false "<<wrong.GetBinContent(wrong.GetNbinsX()+1)
           <<", rejected off-plane "<<offPlane<<"\n";
}
