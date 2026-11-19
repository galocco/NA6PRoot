// plot_rapidity.C
// Reads MCKine.root (TTree "mckine", branch "tracks")
// and plots rapidity, pseudorapidity, multiplicity and PDG composition.
//
// The tracks branch can be either std::vector<TParticle> or TClonesArray.
// TParticle provides: Px(), Py(), Pz(), Energy(), GetPdgCode(),
//                     GetStatusCode(), GetMother(0), IsPrimary()
//
// Usage:
//   root -l -b -q 'plot_rapidity.C("MCKine.root")'
//   root -l -b -q 'plot_rapidity.C("MCKine.root", 211)'           // only pions
//   root -l -b -q 'plot_rapidity.C("MCKine.root", 0, true)'       // split by charge
//   root -l -b -q 'plot_rapidity.C("MCKine.root", 0, false, -4., 4.)' // custom window

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLine.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
double rapidity(double px, double py, double pz, double e) {
    double denom = e - pz;
    if (std::abs(denom) < 1e-10 || denom <= 0.0) return -999.0;
    return 0.5 * std::log((e + pz) / denom);
}

double pseudorapidity(double px, double py, double pz) {
    double pmag = std::sqrt(px*px + py*py + pz*pz);
    if (pmag < 1e-10 || std::abs(pmag - std::abs(pz)) < 1e-10) return -999.0;
    return 0.5 * std::log((pmag + pz) / (pmag - pz));
}

bool isNeutralPDG(int pidabs) {
    return (pidabs == 22  || pidabs == 111 || pidabs == 130 ||
            pidabs == 310 || pidabs == 311 || pidabs == 2112 ||
            pidabs == 12  || pidabs == 14  || pidabs == 16);
}

std::string pdgLabel(int pid) {
    static std::map<int,std::string> names = {
        { 211,  "#pi^{+}"},   {-211,  "#pi^{-}"},  { 111,  "#pi^{0}"},
        { 321,  "K^{+}"},     {-321,  "K^{-}"},
        { 130,  "K^{0}_{L}"}, { 310,  "K^{0}_{S}"},{ 311,  "K^{0}"},
        { 2212, "p"},         {-2212, "#bar{p}"},
        { 2112, "n"},         {-2112, "#bar{n}"},
        {  11,  "e^{-}"},     {  -11, "e^{+}"},
        {  13,  "#mu^{-}"},   {  -13, "#mu^{+}"},
        {  22,  "#gamma"},
        { 3122, "#Lambda"},   {-3122, "#bar{#Lambda}"},
        { 3112, "#Sigma^{-}"},{ 3222, "#Sigma^{+}"},  { 3212, "#Sigma^{0}"},
        { 3312, "#Xi^{-}"},   { 3322, "#Xi^{0}"},
        { 3334, "#Omega^{-}"},
        { 4122, "#Lambda_{c}^{+}"},{-4122,"#bar{#Lambda}_{c}^{-}"},
        {  421, "D^{0}"},     {  411, "D^{+}"},    { -411, "D^{-}"},
        {  431, "D_{s}^{+}"}, { -431, "D_{s}^{-}"},
        {  443, "J/#psi"},    {  553, "#Upsilon"},
        {  223, "#omega"},    {  333, "#phi"},      {  113, "#rho^{0}"},
        {  213, "#rho^{+}"}, { -213, "#rho^{-}"},
        {   12, "#nu_{e}"},   {   14, "#nu_{#mu}"},{   16, "#nu_{#tau}"},
        {  -12, "#bar{#nu}_{e}"},{-14,"#bar{#nu}_{#mu}"}
    };
    if (names.count(pid)) return names.at(pid);
    TDatabasePDG* db = TDatabasePDG::Instance();
    TParticlePDG* part = db->GetParticle(pid);
    if (part) return std::string(part->GetName());
    return std::to_string(pid);
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void plot_rapidity(const std::string& inputFile   = "MCKine.root",
                   int    selectPDG   = 0,         // 0=all charged, >0=specific |PDG|
                   bool   splitCharge = false,
                   double yMin  = -6.0,
                   double yMax  =  6.0,
                   int    nBins = 120)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadRightMargin(0.06);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetHistLineWidth(2);
    TH1::AddDirectory(kFALSE);

    // --- Open file and tree ---
    TFile* f = TFile::Open(inputFile.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open file: " << inputFile << "\n";
        return;
    }

    TTree* tree = nullptr;
    for (const auto& nm : {"mckine", "MCKine", "kine", "Kine"}) {
        tree = (TTree*)f->Get(nm);
        if (tree) { std::cout << "Found tree: " << nm << "\n"; break; }
    }
    if (!tree) { std::cerr << "ERROR: TTree not found.\n"; f->Close(); return; }

    // --- Connect tracks branch (supports both vector<TParticle> and TClonesArray) ---
    TBranch* bTracks = tree->GetBranch("tracks");
    if (!bTracks) {
        std::cerr << "ERROR: branch 'tracks' not found.\n";
        f->Close();
        return;
    }

    std::string tracksClass = bTracks->GetClassName();
    std::vector<TParticle>* tracksVec = nullptr;
    TClonesArray* tracksCA = nullptr;
    bool useVectorTracks = false;

    if (tracksClass.find("vector<") != std::string::npos) {
        tree->SetBranchAddress("tracks", &tracksVec);
        useVectorTracks = true;
        std::cout << "Using tracks branch as " << tracksClass << "\n";
    } else {
        tracksCA = new TClonesArray("TParticle");
        tree->SetBranchAddress("tracks", &tracksCA);
        std::cout << "Using tracks branch as " << tracksClass << "\n";
    }

    // Also read mNTracks from header for safety cross-check
    UInt_t mNTracks = 0;
    if (tree->GetBranch("mNTracks")) {
        tree->SetBranchAddress("mNTracks", &mNTracks);
    }

    // --- Histograms ---
    TH1D* hRap    = new TH1D("hRap",    "", nBins, yMin, yMax);
    TH1D* hRapPos = new TH1D("hRapPos", "", nBins, yMin, yMax);
    TH1D* hRapNeg = new TH1D("hRapNeg", "", nBins, yMin, yMax);
    TH1D* hEta    = new TH1D("hEta",    "", nBins, yMin, yMax);
    TH1I* hMult   = new TH1I("hMult",   "", 300, 0, 3000);

    hRap->SetLineColor(kBlack);    hRap->SetFillColorAlpha(kGray, 0.5);
    hRapPos->SetLineColor(kRed+1);
    hRapNeg->SetLineColor(kBlue+1);
    hEta->SetLineColor(kGreen+2);
    hMult->SetLineColor(kBlack);   hMult->SetFillColorAlpha(kAzure+1, 0.3);

    std::map<int,long> pdgCount;

    std::map<int,std::string> axisNames = {
        {0,"all charged"}, {211,"#pi^{#pm}"}, {111,"#pi^{0}"},
        {321,"K^{#pm}"}, {2212,"p/#bar{p}"}, {2112,"n/#bar{n}"},
        {11,"e^{#pm}"}, {13,"#mu^{#pm}"}
    };
    std::string partLabel = axisNames.count(selectPDG)
        ? axisNames[selectPDG]
        : "PDG |" + std::to_string(selectPDG) + "|";

    // -----------------------------------------------------------------------
    // Event loop
    // -----------------------------------------------------------------------
    long nev   = 0;
    long npart = 0;
    long nEntries = tree->GetEntries();
    std::cout << "Processing " << nEntries << " events...\n";

    for (long iev = 0; iev < 1000; ++iev) {
        tree->GetEntry(iev);

        int nTrk = 0;
        if (useVectorTracks) {
            nTrk = tracksVec ? static_cast<int>(tracksVec->size()) : 0;
        } else {
            nTrk = tracksCA ? tracksCA->GetEntriesFast() : 0;
        }

        if (mNTracks > 0 && static_cast<UInt_t>(nTrk) != mNTracks && iev < 5) {
            std::cout << "Warning: event " << iev << " has tracks=" << nTrk
                      << " but mNTracks=" << mNTracks << "\n";
        }
        int multEv = 0;

        for (int i = 0; i < nTrk; ++i) {
            TParticle* p = nullptr;
            if (useVectorTracks) {
                if (!tracksVec) continue;
                p = &tracksVec->at(i);
            } else {
                if (!tracksCA) continue;
                p = (TParticle*)tracksCA->At(i);
            }
            if (!p) continue;

            int pid    = p->GetPdgCode();
            int pidabs = std::abs(pid);
            // Primary: mother index < 0
            bool primary = (p->GetMother(0) < 0);
            if (!primary) continue;
            pdgCount[pid]++;

            // --- filter for rapidity plots ---
            if (selectPDG != 0 && pidabs != selectPDG) continue;
            if (selectPDG == 0 && isNeutralPDG(pidabs)) continue;

            double px = p->Px(), py = p->Py(), pz = p->Pz(), e = p->Energy();
            double y   = rapidity(px, py, pz, e);
            double eta = pseudorapidity(px, py, pz);
            if (y < -990.0 || eta < -990.0) continue;

            hRap->Fill(y);
            hEta->Fill(eta);
            if (pid > 0) hRapPos->Fill(y);
            else         hRapNeg->Fill(y);

            ++multEv;
            ++npart;
        }

        hMult->Fill(multEv);
        ++nev;
        if (nev % 500 == 0) std::cout << "  " << nev << " events...\n";
    }

    f->Close();

    if (nev == 0) { std::cerr << "ERROR: No events processed!\n"; return; }

    std::cout << "\n=== Summary ===\n"
              << "  Events          : " << nev               << "\n"
              << "  Particles (sel) : " << npart             << "\n"
              << "  <N> per event   : " << (double)npart/nev << "\n"
              << "  Unique PDG codes: " << pdgCount.size()   << "\n\n";

    std::cout << "[dbg] after summary\n";

    // Normalize to dN/dy per event
    double binW = (yMax - yMin) / nBins;
    hRap->Scale(1.0 / (nev * binW));
    hRapPos->Scale(1.0 / (nev * binW));
    hRapNeg->Scale(1.0 / (nev * binW));
    hEta->Scale(1.0 / (nev * binW));

    std::string evLabel = std::to_string(nev) + " events";

    // -----------------------------------------------------------------------
    // Build sorted PDG histogram
    // -----------------------------------------------------------------------
    std::vector<std::pair<long,int>> sorted;
    for (const auto& kv : pdgCount) sorted.push_back({kv.second, kv.first});
    std::sort(sorted.begin(), sorted.end(), [](const auto& a, const auto& b){
        return a.first > b.first;
    });

    std::cout << "[dbg] after sorting PDG map, sorted size=" << sorted.size() << "\n";

    int maxSpecies = std::min((int)sorted.size(), 30);
    int canvasW    = std::max(900, std::max(1, maxSpecies) * 45);

    TH1F* hPDG = nullptr;
    if (maxSpecies > 0) {
        std::cout << "[dbg] building hPDG maxSpecies=" << maxSpecies << "\n";
        hPDG = new TH1F("hPDG", "", maxSpecies, 0, maxSpecies);
        hPDG->SetFillColorAlpha(kOrange+1, 0.75);
        hPDG->SetLineColor(kOrange+3);
        hPDG->GetYaxis()->SetTitle("<N_{primary}> per event");
        hPDG->GetXaxis()->SetLabelSize(0.048);
        hPDG->GetYaxis()->SetTitleOffset(1.4);

        std::cout << "Primary PDG composition (top " << maxSpecies << "):\n";
        for (int i = 0; i < maxSpecies; ++i) {
            long  cnt = sorted[i].first;
            int   pid = sorted[i].second;
            std::string lbl = pdgLabel(pid);
            hPDG->SetBinContent(i+1, (double)cnt/nev);
            hPDG->GetXaxis()->SetBinLabel(i+1, lbl.c_str());
            std::cout << "  " << lbl << " (pid=" << pid << "): "
                      << cnt << "  (<N>=" << (double)cnt/nev << ")\n";
        }
    } else {
        std::cout << "Primary PDG composition: no entries after selection.\n";
    }

    std::cout << "[dbg] before canvas c1\n";

    TLatex tex;
    tex.SetNDC();

    // -----------------------------------------------------------------------
    // Canvas 1 — rapidity
    // -----------------------------------------------------------------------
    TCanvas* c1 = new TCanvas("c1", "Rapidity", 800, 650);
    double yMaxPlot = std::max(1e-6, hRap->GetMaximum()) * 1.28;
    hRap->GetXaxis()->SetTitle("y");
    hRap->GetYaxis()->SetTitle("dN/dy");
    hRap->GetYaxis()->SetRangeUser(0, yMaxPlot);
    hRap->Draw("HIST");
    if (splitCharge) {
        hRapPos->Draw("HIST SAME");
        hRapNeg->Draw("HIST SAME");
    }
    TLine l0(0, 0, 0, yMaxPlot);
    l0.SetLineStyle(2); l0.SetLineColor(kGray+2); l0.Draw();

    TLegend* leg1 = new TLegend(0.65, 0.72, 0.93, 0.90);
    leg1->SetBorderSize(0); leg1->SetFillStyle(0); leg1->SetTextSize(0.034);
    leg1->AddEntry(hRap, (partLabel + " (all)").c_str(), "f");
    if (splitCharge) {
        leg1->AddEntry(hRapPos, "Positive", "l");
        leg1->AddEntry(hRapNeg, "Negative", "l");
    }
    leg1->Draw();

    tex.SetTextSize(0.038); tex.SetTextColor(kBlack);
    tex.DrawLatex(0.15, 0.93, ("dN/dy  #minus  " + partLabel).c_str());
    tex.SetTextSize(0.030); tex.SetTextColor(kGray+2);
    tex.DrawLatex(0.15, 0.875, evLabel.c_str());

    c1->SaveAs("rapidity.pdf"); c1->SaveAs("rapidity.png");
    std::cout << "\nSaved: rapidity.pdf / rapidity.png\n";

    std::cout << "[dbg] before canvas c2\n";

    // -----------------------------------------------------------------------
    // Canvas 2 — y vs eta
    // -----------------------------------------------------------------------
    TCanvas* c2 = new TCanvas("c2", "y vs eta", 800, 650);
    double yMaxPlot2 = std::max(1e-6, std::max(hRap->GetMaximum(), hEta->GetMaximum())) * 1.28;
    hRap->GetYaxis()->SetTitle("dN/dy  or  dN/d#eta");
    hRap->GetYaxis()->SetRangeUser(0, yMaxPlot2);
    hRap->Draw("HIST");
    hEta->Draw("HIST SAME");
    TLine l0b(0, 0, 0, yMaxPlot2);
    l0b.SetLineStyle(2); l0b.SetLineColor(kGray+2); l0b.Draw();

    TLegend* leg2 = new TLegend(0.65, 0.75, 0.93, 0.90);
    leg2->SetBorderSize(0); leg2->SetFillStyle(0); leg2->SetTextSize(0.034);
    leg2->AddEntry(hRap, "dN/dy",    "f");
    leg2->AddEntry(hEta, "dN/d#eta", "l");
    leg2->Draw();

    tex.SetTextColor(kBlack); tex.SetTextSize(0.038);
    tex.DrawLatex(0.15, 0.93, ("y vs #eta  #minus  " + partLabel).c_str());
    tex.SetTextSize(0.030); tex.SetTextColor(kGray+2);
    tex.DrawLatex(0.15, 0.875, evLabel.c_str());

    c2->SaveAs("rapidity_vs_eta.pdf"); c2->SaveAs("rapidity_vs_eta.png");
    std::cout << "Saved: rapidity_vs_eta.pdf / rapidity_vs_eta.png\n";

    std::cout << "[dbg] before canvas c3\n";

    // -----------------------------------------------------------------------
    // Canvas 3 — multiplicity
    // -----------------------------------------------------------------------
    TCanvas* c3 = new TCanvas("c3", "Multiplicity", 800, 650);
    c3->SetLogy();
    hMult->GetXaxis()->SetTitle("N_{ch}");
    hMult->GetYaxis()->SetTitle("Events");
    hMult->Draw("HIST");

    tex.SetTextColor(kBlack); tex.SetTextSize(0.038);
    tex.DrawLatex(0.15, 0.93, ("Charged multiplicity  #minus  " + partLabel).c_str());
    tex.SetTextSize(0.030); tex.SetTextColor(kGray+2);
    std::string meanStr = evLabel + ",  <N> = " +
        std::to_string((int)std::round((double)npart/nev));
    tex.DrawLatex(0.15, 0.875, meanStr.c_str());

    c3->SaveAs("multiplicity.pdf"); c3->SaveAs("multiplicity.png");
    std::cout << "Saved: multiplicity.pdf / multiplicity.png\n";

    std::cout << "[dbg] before canvas c4\n";

    // -----------------------------------------------------------------------
    // Canvas 4 — PDG composition
    // -----------------------------------------------------------------------
    if (hPDG) {
        TCanvas* c4 = new TCanvas("c4", "PDG composition", canvasW, 650);
        c4->SetLogy();
        c4->SetBottomMargin(0.22);

        double pdgMax = std::max(1e-6, hPDG->GetMaximum()) * 8.0;
        double pdgMin = hPDG->GetMinimum(0) * 0.5;
        if (pdgMin <= 0) pdgMin = 1e-3;
        hPDG->GetYaxis()->SetRangeUser(pdgMin, pdgMax);
        hPDG->GetXaxis()->SetLabelSize(0.048);
        hPDG->GetXaxis()->LabelsOption("v");
        hPDG->Draw("HIST text");

        tex.SetTextColor(kBlack); tex.SetTextSize(0.038);
        tex.DrawLatex(0.15, 0.93, "Primary particle composition");
        tex.SetTextSize(0.030); tex.SetTextColor(kGray+2);
        tex.DrawLatex(0.15, 0.875, (evLabel + ",  top " +
            std::to_string(maxSpecies) + " species").c_str());

        c4->SaveAs("pdg_composition.pdf"); c4->SaveAs("pdg_composition.png");
        std::cout << "Saved: pdg_composition.pdf / pdg_composition.png\n";
    }

    std::cout << "\nDone.\n";
}