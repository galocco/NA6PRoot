// resoWithMC.C - Resolution analysis with MC truth matching
#include <unordered_map>
#include <algorithm>

void computeEff(
  bool useMC = false,
  const char* dir = ".",
  const char* dirMC = "../jpsi_bkg_40gev")
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);

  // ---- Open files with error checking
  auto openFile = [](const char* name) -> TFile* {
    TFile* f = TFile::Open(name);
    if (!f || f->IsZombie()) {
      printf("Error: Cannot open file %s\n", name);
      return nullptr;
    }
    return f;
  };

  TFile* fVT = openFile(Form("%s/TracksVerTel.root", dir));
  TFile* fMS = openFile(Form("%s/TracksMuonSpec.root", dir));
  TFile* fM = openFile(Form("%s/TracksMatching%s.root", dir, useMC ? "_MC" : ""));
  TFile* fMC = openFile(Form("%s/MCKine.root", dirMC));

  if (!fVT || !fMS || !fM || !fMC) {
    printf("Failed to open required files. Exiting.\n");
    return;
  }

  // ---- Get trees with error checking
  auto getTree = [](TFile* f, const char* name) -> TTree* {
    TTree* t = (TTree*)f->Get(name);
    if (!t) {
      printf("Error: Cannot find tree %s\n", name);
    }
    return t;
  };

  TTree* tVT = getTree(fVT, "tracksVerTel");
  TTree* tMS = getTree(fMS, "tracksMuonSpec");
  TTree* tM = getTree(fM, "tracksMatching");
  TTree* tMC = getTree(fMC, "mckine");

  if (!tVT || !tMS || !tM || !tMC) {
    printf("Failed to retrieve required trees. Exiting.\n");
    return;
  }

  std::vector<NA6PTrack>* vVT = nullptr;
  std::vector<NA6PTrack>* vMS = nullptr;
  std::vector<NA6PTrack>* vM = nullptr;
  std::vector<TParticle>* vMC = nullptr;

  tVT->SetBranchAddress("VerTel", &vVT);
  tMS->SetBranchAddress("MuonSpec", &vMS);
  tM->SetBranchAddress("Matching", &vM);
  tMC->SetBranchAddress("tracks", &vMC);
  float maxP = 25.0; // max momentum for histograms
  float minP = 0.0;  // min momentum for histograms
  int nPBins = 25;   // number of momentum bins

  TH1F* hEff = new TH1F("hEff", ";#it{p} (GeV/#it{c});# correct match / # matchable", nPBins, minP, maxP);
  TH1F* hEffFromJpsi = new TH1F("hEffFromJpsi", ";#it{p} (GeV/#it{c});# correct match / # matchable", nPBins, minP, maxP);
  
  TH1F* hEffTot = new TH1F("hEffTot", ";#it{p} (GeV/#it{c});# correct match / # gen", nPBins, minP, maxP);
  TH1F* hMatchable = new TH1F("hMatchable", ";#it{p} (GeV/#it{c});# matchable", nPBins, minP, maxP);
  TH1F* hMatchableFromJpsi = new TH1F("hMatchableFromJpsi", ";#it{p} (GeV/#it{c});# matchable from J/psi", nPBins, minP, maxP);
  TH1F* hMC = new TH1F("hMC", ";#it{p} (GeV/#it{c});# gen", nPBins, minP, maxP);
  TH1F* hFake = new TH1F("hFake", ";#it{p} (GeV/#it{c});# fake matches", nPBins, minP, maxP);
  TH1F* hFakeFromJpsi = new TH1F("hFakeFromJpsi", ";#it{p} (GeV/#it{c});# fake matches from J/psi", nPBins, minP, maxP);
  TH1F* hRefitted = new TH1F("hRefitted", ";#it{p} (GeV/#it{c});# refitted matches", nPBins, minP, maxP);
  TH1F* hPID = new TH1F("hPID", ";Particle ID;Entries", 10000, -0.5, 9999.5);  
  TH1F* hPIDmuJpsi = new TH1F("hPIDmuJpsi", ";Particle ID;Entries", 10000, -0.5, 9999.5);  
  TH1F* hMSFromJpsi = new TH1F("hMSFromJpsi", ";#it{p} (GeV/#it{c});# MS tracks from J/psi", nPBins, minP, maxP);
  TH1F* hMotherPDG = new TH1F("hMotherPDG", ";Mother PDG;Entries", 1000, -500.5, 499.5);
  TH1F* hMuonPDG = new TH1F("hMuonPDG", ";Muon PDG;Entries", 100, -50.5, 49.5);
  TH1F* hMCEta = new TH1F("hMCEta", ";#eta;# gen muons from J/psi", 100, -5, 5);
  TH1F* hMSEta = new TH1F("hMSEta", ";#eta;# MS muons from J/psi", 100, -5, 5);

  // Loop over events
  int nEvents = tVT->GetEntries();
  int debugCount = 0;
  int debugCountMS = 0;
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if ((iEvent + 1) % 100 == 0)
      printf("Processing event %d / %d\r", iEvent + 1, nEvents);
    tVT->GetEntry(iEvent);
    tMS->GetEntry(iEvent);
    tM->GetEntry(iEvent);
    tMC->GetEntry(iEvent);
    // Build MC truth map once per event
    std::unordered_map<int, const TParticle*> truth;
    truth.reserve(vMC->size());
    for (size_t idx = 0; idx < vMC->size(); ++idx) {
      auto& p = (*vMC)[idx];
      int trackId = static_cast<int>(idx);  // Use vector index as track ID
      
      hPID->Fill(trackId);
      truth[trackId] = &p;
      
      // Check all muons and their mothers
      if (std::abs(p.GetPdgCode()) == 13) {
        hMuonPDG->Fill(p.GetPdgCode());
        int motherPdg = -1;
        if (p.GetMother(0) >= 0 && p.GetMother(0) < vMC->size()) {
          motherPdg = (*vMC)[p.GetMother(0)].GetPdgCode();
          hMotherPDG->Fill(motherPdg);
          
          if (motherPdg == 443) {
            hMC->Fill(p.P());
            hPIDmuJpsi->Fill(trackId);
            hMCEta->Fill(p.Eta());
            if (debugCount < 5) {
              printf("Event %d: Found mu from J/psi, trackId=%d, P=%.2f, eta=%.2f, mother_idx=%d\n", 
                     iEvent, trackId, p.P(), p.Eta(), p.GetMother(0));
              debugCount++;
            }
          }
        }
      }
    }
    // Create map of particle IDs in VT tracks
    std::unordered_map<int, int> vtTracks;
    for (auto& track : *vVT) {
      vtTracks[track.getParticleID()] = 0;
    }
    // Create map of particle IDs in VT tracks
    std::unordered_map<int, int> matchedTracks;
    for (auto& track : *vM) {
      int pid = track.getParticleID();
      matchedTracks[pid] = 0;

      hRefitted->Fill(track.getP());
      if (pid > 0) continue;
      hFake->Fill(track.getP());
      // check if the fake match comes from muon from jpsi
      int absPid = std::abs(pid);
      auto itTruth = truth.find(absPid);
      if (itTruth != truth.end()) {
        const TParticle* part = itTruth->second;
        int motherPdg = -1;
        if (part->GetMother(0) >= 0 && part->GetMother(0) < vMC->size()) {
          motherPdg = (*vMC)[part->GetMother(0)].GetPdgCode();
        }
        
        if (std::abs(part->GetPdgCode()) == 13 && motherPdg == 443){
          hFakeFromJpsi->Fill(track.getP());
        }
      }
    }

    // Check for matching particle IDs in MS tracks
    for (auto& msTrack : *vMS) {
      int msPid = msTrack.getParticleID();
      
      if (debugCountMS < 10 && msPid > 0) {
        printf("Event %d: MS track with PID=%d, P=%.2f\n", iEvent, msPid, msTrack.getP());
        debugCountMS++;
      }
      
      if (msPid <= 0)
        continue;
      
      // Check if MS track is from J/psi muon
      auto itTruthMS = truth.find(msPid);
      if (itTruthMS != truth.end()) {
        const TParticle* partMS = itTruthMS->second;
        int motherPdgMS = -1;
        if (partMS->GetMother(0) >= 0 && partMS->GetMother(0) < vMC->size()) {
          motherPdgMS = (*vMC)[partMS->GetMother(0)].GetPdgCode();
        }
        if (std::abs(partMS->GetPdgCode()) == 13 && motherPdgMS == 443){
          hMSFromJpsi->Fill(msTrack.getP());
          hMSEta->Fill(partMS->Eta());
        }
      }
      
      if (vtTracks.find(msPid) != vtTracks.end()) {
        hMatchable->Fill(msTrack.getP());
        
        // Check if the matchable comes from muon from J/psi
        auto itTruth = truth.find(msPid);
        bool fromJpsi = false;
        if (itTruth != truth.end()) {
          const TParticle* part = itTruth->second;
          int motherPdg = -1;
          if (part->GetMother(0) >= 0 && part->GetMother(0) < vMC->size()) {
            motherPdg = (*vMC)[part->GetMother(0)].GetPdgCode();
          }
          
          if (std::abs(part->GetPdgCode()) == 13 && motherPdg == 443){
            hMatchableFromJpsi->Fill(msTrack.getP());
            fromJpsi = true;
          }
        }

        // Check if this track was matched
        if (matchedTracks.find(msPid) != matchedTracks.end()) {
          hEff->Fill(msTrack.getP());
          hEffTot->Fill(msTrack.getP());
          // Only fill hEffFromJpsi if this matched track comes from J/psi
          if (fromJpsi) {
            hEffFromJpsi->Fill(msTrack.getP());
          }
        }
      }
    }
  }
  
  printf("\n");
  TFile out(Form("efficiencyMatching%s.root", useMC ? "_MC" : ""), "RECREATE");
  printf("# matched / # matchable = %f\n", hEff->GetEntries() / hMatchable->GetEntries());
  hEff->Divide(hMatchable);
  hEff->Write();
  hEffTot->Divide(hMC);
  hEffTot->Write();
  hEffFromJpsi->Divide(hMatchableFromJpsi);
  hEffFromJpsi->Write();
  hMC->Write();
  hMatchable->Write();
  hMatchableFromJpsi->Write();
  hFake->Divide(hRefitted);
  hFake->Write();
  hFakeFromJpsi->Divide(hRefitted);
  hRefitted->Write();
  hFakeFromJpsi->Write();
  hPID->Write();
  hPIDmuJpsi->Write();
  hMSFromJpsi->Write();
  hMotherPDG->Write();
  hMuonPDG->Write();
  hMCEta->Write();
  hMSEta->Write();
  printf("hPIDmuJpsi entries: %f\n", hPIDmuJpsi->GetEntries());
  printf("hMSFromJpsi entries: %f\n", hMSFromJpsi->GetEntries());
  printf("hMatchableFromJpsi entries: %f\n", hMatchableFromJpsi->GetEntries());
  printf("Acceptance MS: %.4f%%\n", 100.0 * hMSFromJpsi->GetEntries() / hPIDmuJpsi->GetEntries());
  out.Close();

  fVT->Close();
  fMS->Close();
  fM->Close();
  fMC->Close();

  delete fVT;
  delete fMS;
  delete fM;
  delete fMC;
  printf("Done! Saved efficiencyMatching.root\n");
}