#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include "NA6PTrack.h"
#include "NA6PLayoutParam.h"
#endif

void checkMuonTracks(const char* dirSimu = ".", int minNClusters = 6)
{
  auto& param = NA6PLayoutParam::Instance();
  param.updateFromFile(Form("%s/na6pLayout.ini", dirSimu), "", true);

  // --- Open files ---
  TFile* ft = TFile::Open(Form("%s/TracksMuonSpec.root", dirSimu));
  if (!ft || ft->IsZombie()) { printf("Cannot open TracksMuonSpec.root\n"); return; }
  TTree* trTree = (TTree*)ft->Get("tracksMuonSpec");
  if (!trTree) { printf("Cannot find tracksMuonSpec tree\n"); return; }
  std::vector<NA6PTrack>* trArr = nullptr;
  trTree->SetBranchAddress("MuonSpec", &trArr);

  TFile* fk = TFile::Open(Form("%s/MCKine.root", dirSimu));
  if (!fk || fk->IsZombie()) { printf("Cannot open MCKine.root\n"); return; }
  TTree* mcTree = (TTree*)fk->Get("mckine");
  if (!mcTree) { printf("Cannot find mckine tree\n"); return; }
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  // --- PDG label helper: ordered list of expected PDG codes ---
  // Each entry: {pdg code, display name}
  const std::vector<std::pair<int,std::string>> kPDGList = {
    {-13,   "#mu^{+}"},
    { 13,   "#mu^{-}"},
    {-211,  "#pi^{-}"},
    { 211,  "#pi^{+}"},
    {-321,  "K^{-}"},
    { 321,  "K^{+}"},
    {-2212, "#bar{p}"},
    { 2212, "p"},
    {-11,   "e^{+}"},
    { 11,   "e^{-}"},
    { 22,   "#gamma"},
    { 2112, "n"},
    { 111,  "#pi^{0}"},
    { 311,  "K^{0}"},
    { 130,  "K^{0}_{L}"},
    { 310,  "K^{0}_{S}"},
    { 443,  "J/#psi"},
  };
  int nPDGBins = (int)kPDGList.size();

  auto makePDGHist = [&](const char* name, const char* title) -> TH1D* {
    TH1D* h = new TH1D(name, title, nPDGBins, 0, nPDGBins);
    for (int i = 0; i < nPDGBins; i++)
      h->GetXaxis()->SetBinLabel(i + 1, kPDGList[i].second.c_str());
    return h;
  };

  auto fillPDGHist = [&](TH1D* h, int pdg) {
    for (int i = 0; i < nPDGBins; i++) {
      if (kPDGList[i].first == pdg) { h->Fill(i + 0.5); return; }
    }
    // unknown: fill last+1 overflow or just ignore
  };

  TH1D* hPDG       = makePDGHist("hPDG",      "PDG of track origin;;Counts");
  TH1D* hMotherPDG = makePDGHist("hMotherPDG","PDG of mother of track origin;;Counts");
  TH1D* hNClusters  = new TH1D("hNClusters", "N clusters per track;N clusters;Counts", 11, -0.5, 10.5);
  TH1D* hTrackPurity= new TH1D("hTrackPurity","Track purity (most common particle label);Purity;Counts", 110, 0., 1.1);

  int nEv = trTree->GetEntries();
  printf("Number of events = %d\n", nEv);

  int nTracksTotal = 0, nTracksWithMC = 0, nMuons = 0;

  for (int jEv = 0; jEv < nEv; jEv++) {
    mcTree->GetEvent(jEv);
    trTree->GetEvent(jEv);

    int nPart   = mcArr  ? (int)mcArr->size()  : 0;
    int nTracks = trArr  ? (int)trArr->size()  : 0;

    for (int jTr = 0; jTr < nTracks; jTr++) {
      NA6PTrack& tr = trArr->at(jTr);
      if(TMath::Abs(tr.getParticleID()) >= nPart) {
        printf("Pid above MC array size: %d >= %d\n", TMath::Abs(tr.getParticleID()), nPart);
      }
      hNClusters->Fill(tr.getNHits());
      if (tr.getNHits() < minNClusters) continue;
      nTracksTotal++;

      // --- Find majority particle label across cluster layers ---
      std::map<int,int> labelCount;
      for (int lr = 0; lr < NA6PLayoutParam::Instance().nVerTelPlanes + NA6PLayoutParam::Instance().nMSPlanes; lr++) {
        int lbl = tr.getParticleLabel(lr);
        if (lbl >= 0) labelCount[lbl]++;
      }
      int majorityLabel = tr.getParticleID(); // default: assigned MC label
      int majorityCount = 0;
      int totalLabeled  = 0;
      for (auto& [lbl, cnt] : labelCount) {
        totalLabeled += cnt;
        if (cnt > majorityCount) { majorityCount = cnt; majorityLabel = lbl; }
      }
      if (totalLabeled > 0) hTrackPurity->Fill((double)majorityCount / totalLabeled);

      // --- Look up particle in MCKine ---
      int mcLabel = majorityLabel;
      if (mcLabel < 0 || mcLabel >= nPart) {
        printf("  Ev%d Tr%d: no MC label (%d)\n", jEv, jTr, mcLabel);
        continue;
      }
      nTracksWithMC++;

      TParticle& part = mcArr->at(mcLabel);
      int pdg         = part.GetPdgCode();
      int motherIdx   = part.GetFirstMother();

      // particle name
      TParticlePDG* pdgInfo = TDatabasePDG::Instance()->GetParticle(pdg);
      const char* partName  = pdgInfo ? pdgInfo->GetName() : "unknown";

      fillPDGHist(hPDG, pdg);
      if (TMath::Abs(pdg) == 13) nMuons++;

      // --- Mother ---
      int motherPDG = 0;
      const char* motherName = "none";
      if (motherIdx >= 0 && motherIdx < nPart) {
        TParticle& mother = mcArr->at(motherIdx);
        motherPDG = mother.GetPdgCode();
        TParticlePDG* mInfo = TDatabasePDG::Instance()->GetParticle(motherPDG);
        motherName = mInfo ? mInfo->GetName() : "unknown";
        fillPDGHist(hMotherPDG, motherPDG);
      }

      printf("  Ev%3d Tr%3d  NClu=%d  particle: [%6d] %-12s  mother idx=%d  mother: [%6d] %s\n",
             jEv, jTr, tr.getNHits(), pdg, partName, motherIdx, motherPDG, motherName);
    }
  }

  printf("\n=== Summary ===\n");
  printf("Total tracks (>= %d clusters): %d\n", minNClusters, nTracksTotal);
  printf("Tracks with MC label:          %d\n", nTracksWithMC);
  printf("Of which muons (|pdg|==13):    %d\n", nMuons);

  // --- Save output ---
  TFile fOut(Form("%s/checkMuonTracks.root", dirSimu), "RECREATE");
  hPDG->Write();
  hMotherPDG->Write();
  hNClusters->Write();
  hTrackPurity->Write();
  fOut.Close();
  printf("Output saved to %s/checkMuonTracks.root\n", dirSimu);



  gStyle->SetTitleSize(0.06, "XYT");
  gStyle->SetLabelSize(0.05, "XY");
  gStyle->SetOptStat(0);
  hPDG->Scale(1. / nEv);
  hMotherPDG->Scale(1. / nEv);

  TCanvas* cPDG = new TCanvas("cPDG", "PDG of track origin", 1800, 1200);
  cPDG->SetLogy();
  cPDG->SetLeftMargin(0.2);
  cPDG->SetBottomMargin(0.18);
  hPDG->SetTitle("PDG of track origin;;Counts/N_{ev}");
  hPDG->GetXaxis()->SetTitleSize(0.06);
  hPDG->GetYaxis()->SetTitleSize(0.06);
  hPDG->GetXaxis()->SetLabelSize(0.06);
  hPDG->GetYaxis()->SetLabelSize(0.05);
  hPDG->GetXaxis()->LabelsOption("h");
  hPDG->SetFillColor(kAzure + 7);
  hPDG->SetLineColor(kBlue + 1);
  hPDG->Draw("bar");
  cPDG->SaveAs(Form("%s/track_origin_pdg.png", dirSimu));

  TCanvas* cMotherPDG = new TCanvas("cMotherPDG", "PDG of mother of track origin", 1800, 1200);
  cMotherPDG->SetLogy();
  cMotherPDG->SetLeftMargin(0.2);
  cMotherPDG->SetBottomMargin(0.18);
  hMotherPDG->SetTitle("PDG of mother of track origin;;Counts/N_{ev}");
  hMotherPDG->GetXaxis()->SetTitleSize(0.06);
  hMotherPDG->GetYaxis()->SetTitleSize(0.06);
  hMotherPDG->GetXaxis()->SetLabelSize(0.06);
  hMotherPDG->GetYaxis()->SetLabelSize(0.05);
  hMotherPDG->GetXaxis()->LabelsOption("h");
  hMotherPDG->SetFillColor(kOrange - 3);
  hMotherPDG->SetLineColor(kOrange + 2);
  hMotherPDG->Draw("bar");
  cMotherPDG->SaveAs(Form("%s/track_origin_mother_pdg.png", dirSimu));

}
