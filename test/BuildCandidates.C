#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include "NA6PTrack.h"
#include "NA6PMCEventHeader.h"
#include "NA6PDCAFitterN.h"
#include "NA6PLayoutParam.h"
#include "MagneticField.h"
#include "ConfigurableParam.h"
#endif

// ─── helpers ──────────────────────────────────────────────────────────────────

void FillMeanAndRms(TH2F* h2D, TH1F* hMean, TH1F* hRms, TH1F* hSig)
{
  if (!h2D || !hMean || !hRms || !hSig) {
    printf("One of pointers is not set: h2D=%p hMean=%p hRms=%p hSig=%p\n", h2D, hMean, hRms, hSig);
    return;
  }
  int nBins = h2D->GetXaxis()->GetNbins();
  hMean->SetStats(0);
  hRms->SetStats(0);
  hSig->SetStats(0);
  for (int jp = 1; jp <= nBins; jp++) {
    TH1D* htmp1 = h2D->ProjectionY("htmp1", jp, jp);
    double mean = htmp1->GetMean();
    double emean = htmp1->GetMeanError();
    double rms = htmp1->GetRMS();
    double erms = htmp1->GetRMSError();
    hMean->SetBinContent(jp, mean);
    hMean->SetBinError(jp, emean);
    hRms->SetBinContent(jp, rms);
    hRms->SetBinError(jp, erms);
    double gw = rms, egw = erms;
    if (htmp1->GetEntries() > 20.) {
      int fitStatus = htmp1->Fit("gaus", "N0Q", "", mean - 3. * rms, mean + 3. * rms);
      if (fitStatus == 0) {
        TF1* fg = (TF1*)htmp1->GetListOfFunctions()->FindObject("gaus");
        if (fg) {
          gw = fg->GetParameter(2);
          egw = fg->GetParError(2);
        }
      }
    }
    hSig->SetBinContent(jp, gw);
    hSig->SetBinError(jp, egw);
    delete htmp1;
  }
}

// ─────────────────────────────────────────────────────────────────────────────

bool CheckCharmDecay(TParticle& mothPart,
                     std::vector<TParticle>* mcArr,
                     int pdgProngs[5],
                     int idDau[5],
                     float decvert[3])
{
  bool okDau[5] = {false, false, false, false, false};
  int nDau = mothPart.GetNDaughters();
  int finalProngs = nDau;
  int nExpectedProngs = 0;
  for (int i = 0; i < 5; ++i)
    if (pdgProngs[i] != 0)
      nExpectedProngs++;

  bool isDecVertSet = false;
  for (int jd = 0; jd < nDau; ++jd) {
    int indexDau = mothPart.GetFirstDaughter() + jd;
    auto& dauPart = mcArr->at(indexDau);
    int pdgDau = dauPart.GetPdgCode();
    int absPdg = std::abs(pdgDau);
    bool isExp = false;
    for (int ld = 0; ld < 5; ++ld)
      if (pdgDau == pdgProngs[ld]) {
        isExp = true;
        break;
      }
    if (absPdg == 313 || absPdg == 333)
      isExp = true;
    if (isExp && !isDecVertSet) {
      decvert[0] = dauPart.Vx();
      decvert[1] = dauPart.Vy();
      decvert[2] = dauPart.Vz();
      isDecVertSet = true;
    }
  }
  if (!isDecVertSet)
    return false;

  for (int jd = 0; jd < nDau; ++jd) {
    int indexDau = mothPart.GetFirstDaughter() + jd;
    auto& dauPart = mcArr->at(indexDau);
    int pdgDau = dauPart.GetPdgCode();
    int absPdg = std::abs(pdgDau);
    float dist2 = (dauPart.Vx() - decvert[0]) * (dauPart.Vx() - decvert[0]) + (dauPart.Vy() - decvert[1]) * (dauPart.Vy() - decvert[1]) + (dauPart.Vz() - decvert[2]) * (dauPart.Vz() - decvert[2]);
    if (pdgDau == 11 && dist2 > 1e-8) {
      finalProngs--;
      continue;
    }
    bool found = false;
    for (int ld = 0; ld < 5; ++ld) {
      if (pdgDau == pdgProngs[ld] && !okDau[ld]) {
        okDau[ld] = true;
        idDau[ld] = indexDau;
        found = true;
        break;
      }
    }
    if (!found && (absPdg == 313 || absPdg == 333)) {
      int nDauRes = dauPart.GetNDaughters();
      finalProngs = finalProngs - 1 + nDauRes;
      for (int kd = 0; kd < nDauRes; ++kd) {
        int indexDau2 = dauPart.GetFirstDaughter() + kd;
        auto dauPart2 = mcArr->at(indexDau2);
        int pdgDau2 = dauPart2.GetPdgCode();
        for (int ld = 0; ld < 5; ++ld)
          if (pdgDau2 == pdgProngs[ld] && !okDau[ld]) {
            okDau[ld] = true;
            idDau[ld] = indexDau2;
            break;
          }
      }
    }
  }
  if (finalProngs != nExpectedProngs)
    return false;
  for (int ld = 0; ld < 5; ++ld)
    if (pdgProngs[ld] != 0 && !okDau[ld])
      return false;
  return true;
}

// ─────────────────────────────────────────────────────────────────────────────
// Armenteros-Podolanski variables for a 2-prong decay
// alpha = (pL+ - pL-) / (pL+ + pL-)  where pL is longitudinal momentum
//         of each prong w.r.t. the mother direction
// qT    = transverse momentum of one prong w.r.t. the mother direction
// ─────────────────────────────────────────────────────────────────────────────
void ArmenterosPodolanski(float px0, float py0, float pz0, // prong 0
                          float px1, float py1, float pz1, // prong 1
                          float& alpha, float& qt)
{
  float pxM = px0 + px1, pyM = py0 + py1, pzM = pz0 + pz1;
  float pM = std::sqrt(pxM * pxM + pyM * pyM + pzM * pzM);
  if (pM < 1e-9) {
    alpha = 0.;
    qt = 0.;
    return;
  }
  // unit vector along mother
  float nx = pxM / pM, ny = pyM / pM, nz = pzM / pM;
  float pL0 = px0 * nx + py0 * ny + pz0 * nz;
  float pL1 = px1 * nx + py1 * ny + pz1 * nz;
  float sum = pL0 + pL1;
  alpha = (std::abs(sum) > 1e-9) ? (pL0 - pL1) / sum : 0.;
  // qT from prong 0
  float p02 = px0 * px0 + py0 * py0 + pz0 * pz0;
  qt = std::sqrt(std::max(0.f, p02 - pL0 * pL0));
}

// ─────────────────────────────────────────────────────────────────────────────

void BuildCandidates(const char* part = "D0",
                     const char* dirSimu = ".",
                     int minNClusters = 4,
                     int maxEvents = 1000,
                     bool makeCompCanvas = true)
{
  na6p::conf::ConfigurableParam::updateFromFile(Form("%s/na6pLayout.ini", dirSimu), "", true);
  if (TGeoGlobalMagField::Instance()->GetField() == nullptr) {
    auto magField = new MagneticField();
    magField->loadField();
    magField->setAsGlobalField();
  }

  // ── particle / decay setup ─────────────────────────────────────────────────
  int pdgPart = 421;
  int nProngs = 2;
  int pdgProngs[5] = {0, 0, 0, 0, 0};
  int nPosProngs = 0;
  int nNegProngs = 0;
  float maxDL = 5000.;
  float minDL = 0.;
  float maxDLXY = 5000.;
  float minDLXY = 0.;
  float maxCosPA = 0.999999f;
  float maxCosPAXY = 0.9999f;
  float invMassMin = 1.7;
  float invMassMax = 2.5;
  bool isLongLivedV0 = false;
  std::string decay;

  if (strcmp(part, "D0") == 0) {
    pdgPart = 421;
    nProngs = 2;
    pdgProngs[0] = 211;
    pdgProngs[1] = -321;
    maxDL = 1000.;
    decay = "D^{0} #rightarrow K^{-} #pi^{+}";
  } else if (strcmp(part, "D0bar") == 0) {
    pdgPart = -421;
    nProngs = 2;
    pdgProngs[0] = -211;
    pdgProngs[1] = 321;
    maxDL = 1000.;
    decay = "#bar{D^{0}} #rightarrow K^{+} #pi^{-}";
  } else if (strcmp(part, "Dplus") == 0) {
    pdgPart = 411;
    nProngs = 3;
    pdgProngs[0] = 211;
    pdgProngs[1] = 211;
    pdgProngs[2] = -321;
    maxDL = 3000.;
    decay = "D^{+} #rightarrow K^{-} #pi^{+} #pi^{+}";
  } else if (strcmp(part, "Dminus") == 0) {
    pdgPart = -411;
    nProngs = 3;
    pdgProngs[0] = -211;
    pdgProngs[1] = -211;
    pdgProngs[2] = 321;
    maxDL = 3000.;
    decay = "D^{-} #rightarrow K^{+} #pi^{-} #pi^{-}";
  } else if (strcmp(part, "Ds") == 0) {
    pdgPart = 431;
    nProngs = 3;
    pdgProngs[0] = 211;
    pdgProngs[1] = 321;
    pdgProngs[2] = -321;
    maxDL = 1000.;
    decay = "D_{s}^{+} #rightarrow K^{-} K^{+} #pi^{+}";
  } else if (strcmp(part, "DsMinus") == 0) {
    pdgPart = -431;
    nProngs = 3;
    pdgProngs[0] = -211;
    pdgProngs[1] = 321;
    pdgProngs[2] = -321;
    maxDL = 1000.;
    decay = "D_{s}^{-} #rightarrow K^{-} K^{+} #pi^{-}";
  } else if (strcmp(part, "LambdaC") == 0 || strcmp(part, "LcPlus") == 0) {
    pdgPart = 4122;
    nProngs = 3;
    pdgProngs[0] = 2212;
    pdgProngs[1] = 211;
    pdgProngs[2] = -321;
    maxDL = 1000.;
    decay = "#Lambda_{c}^{+} #rightarrow p K^{-} #pi^{+}";
  } else if (strcmp(part, "antiLambdaC") == 0 || strcmp(part, "LcMinus") == 0) {
    pdgPart = -4122;
    nProngs = 3;
    pdgProngs[0] = -2212;
    pdgProngs[1] = -211;
    pdgProngs[2] = 321;
    maxDL = 1000.;
    decay = "#bar{#Lambda}_{c}^{-} #rightarrow #bar{p} K^{+} #pi^{-}";
  } else if (strcmp(part, "Lambda0") == 0 || strcmp(part, "Lambda") == 0) {
    pdgPart = 3122;
    nProngs = 2;
    pdgProngs[0] = 2212;
    pdgProngs[1] = -211;
    isLongLivedV0 = true;
    maxDL = 200000.;
    minDL = 2.;
    maxDLXY = 50000.;
    minDLXY = 0.1;
    invMassMin = TDatabasePDG::Instance()->GetParticle(pdgPart)->Mass() - 0.03;
    invMassMax = TDatabasePDG::Instance()->GetParticle(pdgPart)->Mass() + 0.03;
    decay = "#Lambda^{0} #rightarrow p #pi^{-}";
  } else if (strcmp(part, "antiLambda0") == 0 || strcmp(part, "antiLambda") == 0) {
    pdgPart = -3122;
    nProngs = 2;
    pdgProngs[0] = -2212;
    pdgProngs[1] = 211;
    isLongLivedV0 = true;
    maxDL = 200000.;
    minDL = 2.;
    maxDLXY = 50000.;
    minDLXY = 0.1;
    invMassMin = TDatabasePDG::Instance()->GetParticle(pdgPart)->Mass() - 0.03;
    invMassMax = TDatabasePDG::Instance()->GetParticle(pdgPart)->Mass() + 0.03;
    decay = "#bar{#Lambda}^{0} #rightarrow #bar{p} #pi^{+}";
  } else if (strcmp(part, "K0s") == 0 || strcmp(part, "Kshort") == 0 || strcmp(part, "Ks") == 0) {
    pdgPart = 310;
    nProngs = 2;
    pdgProngs[0] = 211;
    pdgProngs[1] = -211;
    isLongLivedV0 = true;
    maxDL = 200000.;
    invMassMin = 0.44;
    invMassMax = 0.56;
    decay = "K_{S}^{0} #rightarrow #pi^{+} #pi^{-}";
  }

  for (int ip = 0; ip < nProngs; ++ip) {
    if (pdgProngs[ip] > 0)
      ++nPosProngs;
    else if (pdgProngs[ip] < 0)
      ++nNegProngs;
  }

  // ── open input files / trees ───────────────────────────────────────────────
  TFile* ft = new TFile(Form("%s/TracksVerTel.root", dirSimu));
  if (!ft)
    return;
  TTree* trTree = (TTree*)ft->Get("tracksVerTel");
  std::vector<NA6PTrack>* trArr = nullptr;
  trTree->SetBranchAddress("VerTel", &trArr);

  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  NA6PMCEventHeader* mcHead = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);
  mcTree->SetBranchAddress("header", &mcHead);

  // ── fitters ────────────────────────────────────────────────────────────────
  NA6PDCAFitterN<2> df2;
  df2.setPropagateToPCA(true);
  df2.setUseAbsDCA(true);
  NA6PDCAFitterN<3> df3;
  df3.setPropagateToPCA(true);
  df3.setUseAbsDCA(true);
  if (isLongLivedV0) {
    df2.setMaxDZIni(500.);
    df2.setMinVertZ(-500.);
    df2.setMaxVertZ(200.);
  }

  // ── layout ────────────────────────────────────────────────────────────────
  auto& layoutPar = NA6PLayoutParam::Instance();
  int nTargets = std::max(1, std::min((int)layoutPar.nTargets, (int)NA6PLayoutParam::MaxTargets));
  std::vector<float> zTarget(nTargets, 0.f);
  for (int it = 0; it < nTargets; ++it)
    zTarget[it] = layoutPar.posTargetZ[it] + layoutPar.shiftTargets[2];

  auto getTargetIndex = [&](float zVtx) -> int {
    int best = 0;
    float bestDz = 1.e9;
    for (int it = 0; it < nTargets; ++it) {
      float dz = std::abs(zVtx - zTarget[it]);
      float halfThick = 0.5f * layoutPar.thicknessTarget[it];
      if (dz <= halfThick + 1e-5f)
        return it;
      if (dz < bestDz) {
        bestDz = dz;
        best = it;
      }
    }
    return best;
  };

  // ── output tree ───────────────────────────────────────────────────────────
  TFile* fout = new TFile(Form("%s-Candidates.root", part), "recreate");
  TTree* candTree = new TTree("candidates", Form("%s candidates", part));
  TTree* genTree = new TTree("generatedMothers", Form("%s generated mothers", part));

  // branch variables
  float t_mass = 0.f;
  float t_pt = 0.f;
  float t_rap = 0.f;
  // DCAs between all prong pairs (up to 3 prongs -> 3 pairs: 01,02,12)
  float t_dca01 = 0.f, t_dca02 = 0.f, t_dca12 = 0.f;
  // cosine of pointing angle (2D and 3D) w.r.t. primary vertex
  float t_cpa = 0.f;   // 3D
  float t_cpaXY = 0.f; // transverse
  // secondary vertex position
  float t_svx = 0.f, t_svy = 0.f, t_svz = 0.f;
  // decay length
  float t_decLen = 0.f;
  float t_decLenXY = 0.f;
  // Armenteros (only filled for 2-prong)
  float t_armAlpha = 0.f;
  float t_armQt = 0.f;
  // signal flag: 1 = signal, 0 = background
  int t_isSignal = 0;
  // primary vertex (for pointing angle)
  float t_pvx = 0.f, t_pvy = 0.f, t_pvz = 0.f;
  float g_pt = 0.f;
  float g_y = 0.f;
  float g_pvx = 0.f, g_pvy = 0.f, g_pvz = 0.f;
  float g_dvx = 0.f, g_dvy = 0.f, g_dvz = 0.f;

  candTree->Branch("mass", &t_mass, "mass/F");
  candTree->Branch("pt", &t_pt, "pt/F");
  candTree->Branch("rap", &t_rap, "rap/F");
  candTree->Branch("dca01", &t_dca01, "dca01/F");
  candTree->Branch("dca02", &t_dca02, "dca02/F");
  candTree->Branch("dca12", &t_dca12, "dca12/F");
  candTree->Branch("cpa", &t_cpa, "cpa/F");
  candTree->Branch("cpaXY", &t_cpaXY, "cpaXY/F");
  candTree->Branch("svx", &t_svx, "svx/F");
  candTree->Branch("svy", &t_svy, "svy/F");
  candTree->Branch("svz", &t_svz, "svz/F");
  candTree->Branch("decLen", &t_decLen, "decLen/F");
  candTree->Branch("decLenXY", &t_decLenXY, "decLenXY/F");
  candTree->Branch("armAlpha", &t_armAlpha, "armAlpha/F");
  candTree->Branch("armQt", &t_armQt, "armQt/F");
  candTree->Branch("isSignal", &t_isSignal, "isSignal/I");
  candTree->Branch("pvx", &t_pvx, "pvx/F");
  candTree->Branch("pvy", &t_pvy, "pvy/F");
  candTree->Branch("pvz", &t_pvz, "pvz/F");

  genTree->Branch("pt", &g_pt, "pt/F");
  genTree->Branch("y", &g_y, "y/F");
  genTree->Branch("pvx", &g_pvx, "pvx/F");
  genTree->Branch("pvy", &g_pvy, "pvy/F");
  genTree->Branch("pvz", &g_pvz, "pvz/F");
  genTree->Branch("dvx", &g_dvx, "dvx/F");
  genTree->Branch("dvy", &g_dvy, "dvy/F");
  genTree->Branch("dvz", &g_dvz, "dvz/F");

  // ── monitoring histos (kept for quick QA plots) ───────────────────────────
  TH1F* hInvMassSig = new TH1F("hInvMassSig", ";M (GeV/c^{2});entries", 200, invMassMin, invMassMax);
  TH1F* hInvMassBkg = new TH1F("hInvMassBkg", ";M (GeV/c^{2});entries", 200, invMassMin, invMassMax);

  // ── event loop ────────────────────────────────────────────────────────────
  int nEv = trTree->GetEntries();
  if (maxEvents > 0 && maxEvents < nEv)
    nEv = maxEvents;
  printf("Number of events = %d\n", nEv);

  for (int jEv = 0; jEv < nEv; jEv++) {
    mcTree->GetEvent(jEv);
    trTree->GetEvent(jEv);
    printf("--- Event %d ---\n", jEv);

    int nPart = mcArr->size();
    int nTracks = trArr->size();

    double xPV = mcHead->getVX();
    double yPV = mcHead->getVY();
    double zPV = mcHead->getVZ();
    t_pvx = (float)xPV;
    t_pvy = (float)yPV;
    t_pvz = (float)zPV;
    g_pvx = static_cast<float>(xPV);
    g_pvy = static_cast<float>(yPV);
    g_pvz = static_cast<float>(zPV);

    // field at primary vertex
    double bxyz[3], xyz[3] = {xPV, yPV, zPV};
    TGeoGlobalMagField::Instance()->Field(xyz, bxyz);
    df2.setBy(bxyz[1]);
    df3.setBy(bxyz[1]);

    // ── build MC truth signal map for this event ──────────────────────────
    // For each (sorted) set of track MC labels that form a true signal decay,
    // we store the set so we can flag combinatorial candidates.
    // Key = sorted vector of idDau indices for the true decay prongs.
    std::vector<std::set<int>> signalSets;

    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      if (curPart.GetPdgCode() != pdgPart)
        continue;
      int idDau[5] = {-1, -1, -1, -1, -1};
      float decvert[3] = {0., 0., -9999.};
      if (!CheckCharmDecay(curPart, mcArr, pdgProngs, idDau, decvert))
        continue;
      g_pt = static_cast<float>(curPart.Pt());
      g_y = static_cast<float>(curPart.Y());
      g_dvx = decvert[0];
      g_dvy = decvert[1];
      g_dvz = decvert[2];
      genTree->Fill();
      std::set<int> sigSet;
      for (int kd = 0; kd < nProngs; ++kd)
        if (idDau[kd] >= 0)
          sigSet.insert(idDau[kd]);
      if ((int)sigSet.size() == nProngs)
        signalSets.push_back(sigSet);
    }

    // ── select tracks passing quality cut ─────────────────────────────────
    std::vector<NA6PTrack> goodTracks;
    goodTracks.reserve(nTracks);
    for (int jTr = 0; jTr < nTracks; ++jTr) {
      NA6PTrack tr = trArr->at(jTr);
      if (tr.getNVTHits() >= minNClusters)
        goodTracks.push_back(tr);
    }
    int nGood = goodTracks.size();

    // invariant mass (use PDG masses of expected prongs)
    float m0 = TDatabasePDG::Instance()->GetParticle(pdgProngs[0])->Mass();
    float m1 = TDatabasePDG::Instance()->GetParticle(pdgProngs[1])->Mass();
    // ── combinatorial loop ─────────────────────────────────────────────────
    // 2-prong
    if (nProngs == 2) {
      for (int i0 = 0; i0 < nGood - 1; ++i0) {
        float ch0 = goodTracks[i0].getCharge();
        for (int i1 = i0 + 1; i1 < nGood; ++i1) {
          float ch1 = goodTracks[i1].getCharge();
          if (ch0 * ch1 > 0)
            continue; // only opposite charge pairs
          // compute the mass of the mother for a preselection (before vertex fit) to speed up processing of long-lived V0s
          
          double pxyz0[3], pxyz1[3];
          goodTracks[i0].getPXYZ(pxyz0);
          goodTracks[i1].getPXYZ(pxyz1);
          double px0 = pxyz0[0], py0 = pxyz0[1], pz0 = pxyz0[2];
          double px1 = pxyz1[0], py1 = pxyz1[1], pz1 = pxyz1[2];
          double eSum = std::sqrt(m0 * m0 + px0 * px0 + py0 * py0 + pz0 * pz0) + std::sqrt(m1 * m1 + px1 * px1 + py1 * py1 + pz1 * pz1);
          double p2Sum = (px0 + px1) * (px0 + px1) + (py0 + py1) * (py0 + py1) + (pz0 + pz1) * (pz0 + pz1);
          t_mass = std::sqrt(std::max(0., eSum * eSum - p2Sum));
          if (t_mass < invMassMin || t_mass > invMassMax+0.02)
            continue;
          
          // vertex fit
          int n = df2.process(goodTracks[i0], goodTracks[i1]);
          if (n == 0)
            continue;

          // pick best candidate (closest to beam axis as proxy)
          int bestCand = 0;
          double bestR2 = 1.e30;
          for (int ic = 0; ic < df2.getNCandidates(); ++ic) {
            const auto& v = df2.getPCACandidate(ic);
            double r2 = v[0] * v[0] + v[1] * v[1];
            if (r2 < bestR2) {
              bestR2 = r2;
              bestCand = ic;
            }
          }
          const auto& sv = df2.getPCACandidate(bestCand);
          t_svx = (float)sv[0];
          t_svy = (float)sv[1];
          t_svz = (float)sv[2];

          NA6PTrack p0 = goodTracks[i0], p1 = goodTracks[i1];
          p0.propagateToZBxByBz(sv[2]);
          p1.propagateToZBxByBz(sv[2]);
          double xyz0[3], xyz1[3];
          p0.getXYZ(xyz0);
          p1.getXYZ(xyz1);
          p0.getPXYZ(pxyz0);
          p1.getPXYZ(pxyz1);
          double ddx = xyz0[0] - xyz1[0], ddy = xyz0[1] - xyz1[1], ddz = xyz0[2] - xyz1[2];
          t_dca01 = (float)std::sqrt(ddx * ddx + ddy * ddy + ddz * ddz) * 1e4; // µm

          t_dca02 = 0.f;
          t_dca12 = 0.f; // not applicable for 2-prong

          float px = (float)(pxyz0[0] + pxyz1[0]);
          float py = (float)(pxyz0[1] + pxyz1[1]);
          float pz = (float)(pxyz0[2] + pxyz1[2]);
          float p2 = px * px + py * py + pz * pz;

          float mom0 = p0.getP(), mom1 = p1.getP();
          float e0 = std::sqrt(m0 * m0 + mom0 * mom0);
          float e1 = std::sqrt(m1 * m1 + mom1 * mom1);
          eSum = e0 + e1;
          t_mass = std::sqrt(eSum * eSum - p2);

          if (t_mass < invMassMin || t_mass > invMassMax)
            continue;

          t_pt = std::sqrt(px * px + py * py);

          // rapidity (use particle mass from PDG)
          float mPart = TDatabasePDG::Instance()->GetParticle(pdgPart)->Mass();
          float ePart = std::sqrt(mPart * mPart + p2);
          t_rap = (std::abs(ePart + pz) > 1e-9 && std::abs(ePart - pz) > 1e-9) ? 0.5f * std::log((ePart + pz) / (ePart - pz)) : 0.f;

          // pointing angle
          float dx = t_svx - (float)xPV;
          float dy = t_svy - (float)yPV;
          float dz = t_svz - (float)zPV;
          float lxyz = std::sqrt(dx * dx + dy * dy + dz * dz);
          if (lxyz < minDL)
            continue; // decay length cut
          float lxy = std::sqrt(dx * dx + dy * dy);
          if (lxy < minDLXY)
            continue; // decay length cut
          float ptot = std::sqrt(p2);
          t_cpa = (lxyz > 1e-9 && ptot > 1e-9) ? (dx * px + dy * py + dz * pz) / (lxyz * ptot) : 0.f;
          if (t_cpa < maxCosPA)
            continue; // loose pointing angle cut to reduce combinatorics for 3-prong

          t_cpaXY = (lxy > 1e-9 && t_pt > 1e-9) ? (dx * px + dy * py) / (lxy * t_pt) : 0.f;
          if (t_cpaXY < maxCosPAXY)
            continue;            // loose transverse pointing angle cut to reduce combinatorics for 3-prong
          t_decLen = lxyz * 1e4; // µm
          t_decLenXY = lxy * 1e4;

          // Armenteros
          ArmenterosPodolanski((float)pxyz0[0], (float)pxyz0[1], (float)pxyz0[2],
                               (float)pxyz1[0], (float)pxyz1[1], (float)pxyz1[2],
                               t_armAlpha, t_armQt);

          // signal flag
          int lab0 = goodTracks[i0].getParticleID();
          int lab1 = goodTracks[i1].getParticleID();
          t_isSignal = 0;
          for (auto& sigSet : signalSets) {
            if (sigSet.count(lab0) && sigSet.count(lab1)) {
              t_isSignal = 1;
              break;
            }
          }

          candTree->Fill();
          if (t_isSignal)
            hInvMassSig->Fill(t_mass);
          else
            hInvMassBkg->Fill(t_mass);

        } // i1
      } // i0
    } // nProngs==2

    // ── 3-prong combinatorial loop ─────────────────────────────────────────
    else if (nProngs == 3) {
      for (int i0 = 0; i0 < nGood - 2; ++i0) {
        float ch0 = goodTracks[i0].getCharge();
        for (int i1 = i0 + 1; i1 < nGood - 1; ++i1) {
          float ch1 = goodTracks[i1].getCharge();
          for (int i2 = i1 + 1; i2 < nGood; ++i2) {
            float ch2 = goodTracks[i2].getCharge();
            int nPosCand = (ch0 > 0.f) + (ch1 > 0.f) + (ch2 > 0.f);
            int nNegCand = (ch0 < 0.f) + (ch1 < 0.f) + (ch2 < 0.f);
            if (nPosCand != nPosProngs || nNegCand != nNegProngs)
              continue;

            int n = df3.process(goodTracks[i0], goodTracks[i1], goodTracks[i2]);
            if (n == 0)
              continue;

            int bestCand = 0;
            double bestR2 = 1.e30;
            for (int ic = 0; ic < df3.getNCandidates(); ++ic) {
              const auto& v = df3.getPCACandidate(ic);
              double r2 = v[0] * v[0] + v[1] * v[1];
              if (r2 < bestR2) {
                bestR2 = r2;
                bestCand = ic;
              }
            }
            const auto& sv = df3.getPCACandidate(bestCand);
            t_svx = (float)sv[0];
            t_svy = (float)sv[1];
            t_svz = (float)sv[2];

            // propagate prongs to secondary vertex z
            NA6PTrack p0 = goodTracks[i0], p1 = goodTracks[i1], p2t = goodTracks[i2];
            p0.propagateToZBxByBz(sv[2]);
            p1.propagateToZBxByBz(sv[2]);
            p2t.propagateToZBxByBz(sv[2]);
            double xyz0[3], xyz1[3], xyz2[3];
            p0.getXYZ(xyz0);
            p1.getXYZ(xyz1);
            p2t.getXYZ(xyz2);

            // pairwise DCAs at vertex (µm)
            auto pairDCA = [&](double* a, double* b) -> float {
              double ddx = a[0] - b[0], ddy = a[1] - b[1], ddz = a[2] - b[2];
              return (float)std::sqrt(ddx * ddx + ddy * ddy + ddz * ddz) * 1e4;
            };
            t_dca01 = pairDCA(xyz0, xyz1);
            t_dca02 = pairDCA(xyz0, xyz2);
            t_dca12 = pairDCA(xyz1, xyz2);

            double pxyz0[3], pxyz1[3], pxyz2[3];
            p0.getPXYZ(pxyz0);
            p1.getPXYZ(pxyz1);
            p2t.getPXYZ(pxyz2);

            float px = (float)(pxyz0[0] + pxyz1[0] + pxyz2[0]);
            float py = (float)(pxyz0[1] + pxyz1[1] + pxyz2[1]);
            float pz = (float)(pxyz0[2] + pxyz1[2] + pxyz2[2]);
            float p2 = px * px + py * py + pz * pz;

            // use PDG masses of expected prongs
            int pids[3] = {pdgProngs[0], pdgProngs[1], pdgProngs[2]};
            float masses[3];
            for (int k = 0; k < 3; ++k)
              masses[k] = (pids[k] != 0) ? TDatabasePDG::Instance()->GetParticle(pids[k])->Mass() : 0.14f;
            float moms[3] = {static_cast<float>(p0.getP()), static_cast<float>(p1.getP()), static_cast<float>(p2t.getP())};
            float eSum = 0.f;
            for (int k = 0; k < 3; ++k)
              eSum += std::sqrt(masses[k] * masses[k] + moms[k] * moms[k]);
            t_mass = std::sqrt(std::max(0.f, eSum * eSum - p2));
            t_pt = std::sqrt(px * px + py * py);

            float mPart = TDatabasePDG::Instance()->GetParticle(pdgPart) ? TDatabasePDG::Instance()->GetParticle(pdgPart)->Mass() : t_mass;
            float ePart = std::sqrt(mPart * mPart + p2);
            t_rap = (std::abs(ePart + pz) > 1e-9 && std::abs(ePart - pz) > 1e-9) ? 0.5f * std::log((ePart + pz) / (ePart - pz)) : 0.f;

            float dx = t_svx - (float)xPV, dy = t_svy - (float)yPV, dz = t_svz - (float)zPV;
            float lxyz = std::sqrt(dx * dx + dy * dy + dz * dz);
            float lxy = std::sqrt(dx * dx + dy * dy);
            float ptot = std::sqrt(p2), ptxy = std::sqrt(px * px + py * py);
            t_cpa = (lxyz > 1e-9 && ptot > 1e-9) ? (dx * px + dy * py + dz * pz) / (lxyz * ptot) : 0.f;
            if (t_cpa < 0.9f)
              continue; // loose pointing angle cut to reduce combinatorics for 3-prong
            t_cpaXY = (lxy > 1e-9 && ptxy > 1e-9) ? (dx * px + dy * py) / (lxy * ptxy) : 0.f;
            t_decLen = lxyz * 1e4;
            t_decLenXY = lxy * 1e4;

            // Armenteros: use first two prongs (convention)
            ArmenterosPodolanski((float)pxyz0[0], (float)pxyz0[1], (float)pxyz0[2],
                                 (float)pxyz1[0], (float)pxyz1[1], (float)pxyz1[2],
                                 t_armAlpha, t_armQt);

            int lab0 = goodTracks[i0].getParticleID();
            int lab1 = goodTracks[i1].getParticleID();
            int lab2 = goodTracks[i2].getParticleID();
            t_isSignal = 0;
            for (auto& sigSet : signalSets) {
              if (sigSet.count(lab0) && sigSet.count(lab1) && sigSet.count(lab2)) {
                t_isSignal = 1;
                break;
              }
            }

            if (t_mass >= invMassMin && t_mass <= invMassMax) {
              candTree->Fill();
              if (t_isSignal)
                hInvMassSig->Fill(t_mass);
              else
                hInvMassBkg->Fill(t_mass);
            }
          } // i2
        } // i1
      } // i0
    } // nProngs==3
  } // event loop

  // ── quick QA canvas ───────────────────────────────────────────────────────
  TCanvas* cqm = new TCanvas("cqm", "Inv Mass sig vs bkg", 700, 500);
  hInvMassBkg->SetLineColor(kBlue + 1);
  hInvMassBkg->SetFillColorAlpha(kBlue + 1, 0.25);
  hInvMassSig->SetLineColor(kRed + 1);
  hInvMassSig->SetFillColorAlpha(kRed + 1, 0.35);
  hInvMassBkg->SetTitle(Form(";M (GeV/c^{2});entries  [%s]", part));
  hInvMassBkg->Draw("hist");
  hInvMassSig->Draw("hist same");
  TLegend* leg = new TLegend(0.60, 0.65, 0.88, 0.82);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hInvMassSig, "Signal", "f");
  leg->AddEntry(hInvMassBkg, "Background", "f");
  leg->Draw();
  TLatex tdec;
  tdec.SetNDC();
  tdec.DrawLatex(0.13, 0.92, decay.data());
  cqm->SaveAs(Form("%s-InvMassSigBkg.png", part));

  if (makeCompCanvas) {
    // ── normalized signal vs background comparison for all tree quantities ───
    struct VarDef {
      const char* name;
      const char* title;
      int nbins;
      float vmin;
      float vmax;
    };

    float dCosPA = (1 - maxCosPA) / 10.f;
    float dCosPAXY = (1 - maxCosPAXY) / 10.f;
    std::vector<VarDef> vars = {
      {"mass", ";M (GeV/c^{2});a.u.", 120, invMassMin, invMassMax},
      {"pt", ";p_{T} (GeV/c);a.u.", 120, 0.f, 8.f},
      {"rap", ";y;a.u.", 120, 0.f, 6.f},
      {"dca01", ";DCA_{01} (#mum);a.u.", 120, 0.f, 300.f},
      {"dca02", ";DCA_{02} (#mum);a.u.", 120, 0.f, 300.f},
      {"dca12", ";DCA_{12} (#mum);a.u.", 120, 0.f, 300.f},
      {"cpa", ";CPA;a.u.", 120, maxCosPA - dCosPA, 1 + dCosPA},
      {"cpaXY", ";CPA_{XY};a.u.", 120, maxCosPAXY - dCosPAXY, 1 + dCosPAXY},
      {"svx", ";x_{SV} (cm);a.u.", 120, -5.f, 5.f},
      {"svy", ";y_{SV} (cm);a.u.", 120, -5.f, 5.f},
      {"svz", ";z_{SV} (cm);a.u.", 120, -80.f, 80.f},
      {"decLen", ";L_{xyz} (#mum);a.u.", 120, minDL, maxDL},
      {"decLenXY", ";L_{xy} (#mum);a.u.", 120, minDLXY, maxDLXY},
      {"armAlpha", ";#alpha;a.u.", 120, -1.2f, 1.2f},
      {"armQt", ";q_{T} (GeV/c);a.u.", 120, 0.f, 1.2f},
      {"pvx", ";x_{PV} (cm);a.u.", 120, -1.f, 1.f},
      {"pvy", ";y_{PV} (cm);a.u.", 120, -1.f, 1.f},
      {"pvz", ";z_{PV} (cm);a.u.", 120, -10.f, 10.f},
    };

    int nVars = vars.size();
    int nCols = 6;
    int nRows = (nVars + nCols - 1) / nCols;
    TCanvas* ccmp = new TCanvas("ccmp", "Signal vs Background (normalized)", 450 * nCols, 330 * nRows);
    ccmp->Divide(nCols, nRows);

    TLegend* legCmp = new TLegend(0.48, 0.72, 0.88, 0.88);
    legCmp->SetBorderSize(0);
    legCmp->SetFillStyle(0);

    for (int iv = 0; iv < nVars; ++iv) {
      ccmp->cd(iv + 1);
      gPad->SetTickx();
      gPad->SetTicky();

      auto& v = vars[iv];
      TH1F* hSig = new TH1F(Form("hSig_%s", v.name), v.title, v.nbins, v.vmin, v.vmax);
      TH1F* hBkg = new TH1F(Form("hBkg_%s", v.name), v.title, v.nbins, v.vmin, v.vmax);

      candTree->Draw(Form("%s>>hSig_%s", v.name, v.name), "isSignal==1", "goff");
      candTree->Draw(Form("%s>>hBkg_%s", v.name, v.name), "isSignal==0", "goff");

      double intSig = hSig->Integral();
      double intBkg = hBkg->Integral();
      if (intSig > 0.)
        hSig->Scale(1. / intSig);
      if (intBkg > 0.)
        hBkg->Scale(1. / intBkg);

      hSig->SetLineColor(kRed + 1);
      hSig->SetLineWidth(2);
      hBkg->SetLineColor(kBlue + 1);
      hBkg->SetLineWidth(2);
      hBkg->SetLineStyle(2);
      hSig->SetStats(0);
      hBkg->SetStats(0);

      double maxY = std::max(hSig->GetMaximum(), hBkg->GetMaximum());
      hSig->SetMaximum(maxY > 0. ? 1.2 * maxY : 1.);
      hSig->Draw("hist");
      hBkg->Draw("hist same");

      if (iv == 0) {
        legCmp->AddEntry(hSig, "Signal (norm)", "l");
        legCmp->AddEntry(hBkg, "Background (norm)", "l");
        legCmp->Draw();
      }
    }
    ccmp->SaveAs(Form("%s-CandidatesSigBkgComparison.png", part));
  }

  // ── write everything ──────────────────────────────────────────────────────
  fout->cd();
  candTree->Write();
  genTree->Write();
  hInvMassSig->Write();
  hInvMassBkg->Write();
  fout->Close();

  printf("\nDone. Output: %s-Candidates.root\n", part);
  printf("Additional tree: generatedMothers(pt, y, pvx, pvy, pvz, dvx, dvy, dvz)\n");
  printf("Tree branches: mass, pt, rap, dca01, dca02, dca12,\n");
  printf("               cpa, cpaXY, svx, svy, svz, decLen, decLenXY,\n");
  printf("               armAlpha, armQt, isSignal, pvx, pvy, pvz\n");
}