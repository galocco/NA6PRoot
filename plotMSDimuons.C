#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TGeoManager.h>
#include <TStyle.h>
#include <TString.h>
#include <TLatex.h>
#include <TLorentzVector.h>
#include <cstdio>
#include <exception>
#include <vector>
#include "NA6PMCComposedLabel.h"
#include "NA6PMCTruthContainer.h"
#include "NA6PEventReader.h"
#include "NA6PTrack.h"
#include "MagneticField.h"
#include "NA6PMuonSpecModularHit.h"
#include "NA6PBeamParam.h"
#include "NA6PLayoutParam.h"
#include "NA6PFastTrackFitter.h"
#include "Propagator.h"
#endif

void fillMassMeanVsYBins(TH2D* source, TH1D* fittedMean, TH1D* arithmeticMean, const char* projPrefix, const char* fitPrefix)
{
  if (!source || !fittedMean || !arithmeticMean) {
    return;
  }

  for (int i = 1; i <= source->GetNbinsY(); i++) {
    double yCenter = source->GetYaxis()->GetBinCenter(i);
    TH1D* projection = source->ProjectionX(Form("%s_bin%d", projPrefix, i), i, i);
    if (projection->GetEntries() > 10) {
      TF1* fit = new TF1(Form("%s_bin%d", fitPrefix, i), "gaus", 2.5, 2.75);
      projection->Fit(fit, "Q", "", 2.85, 2.735);
      double meanMass = fit->GetParameter(1);
      double meanMassErr = fit->GetParError(1);
      int bin = fittedMean->FindBin(yCenter);
      fittedMean->SetBinContent(bin, meanMass);
      fittedMean->SetBinError(bin, meanMassErr);
      arithmeticMean->SetBinContent(bin, projection->GetMean());
      arithmeticMean->SetBinError(bin, projection->GetMeanError());

      delete fit;
    }
    delete projection;
  }
}

void fillMassResolutionVsBins(TH2D* source, TH1D* resolution, const char* projPrefix, const char* fitPrefix)
{
  if (!source || !resolution) return;

  for (int i = 1; i <= source->GetNbinsY(); ++i) {
    TH1D* projection = source->ProjectionX(Form("%s_bin%d", projPrefix, i), i, i);
    if (projection->GetEntries() > 10 && projection->GetRMS() > 0.) {
      const double fitMin = TMath::Max(projection->GetXaxis()->GetXmin(), projection->GetMean() - 1.5 * projection->GetRMS());
      const double fitMax = TMath::Min(projection->GetXaxis()->GetXmax(), projection->GetMean() + 1.5 * projection->GetRMS());
      TF1* fit = new TF1(Form("%s_bin%d", fitPrefix, i), "gaus", fitMin, fitMax);
      projection->Fit(fit, "QN", "", fitMin, fitMax);
      const double sigma = TMath::Abs(fit->GetParameter(2));
      const double sigmaError = fit->GetParError(2);
      if (sigma > 0. && sigmaError >= 0.) {
        resolution->SetBinContent(i, 1000. * sigma);
        resolution->SetBinError(i, 1000. * sigmaError);
      }
      delete fit;
    }
    delete projection;
  }
}

void plotMSDimuons(int pdg = 443, int NRecoClusters = 6, const char* dirSimu = ".",
                   bool VertexConstrained = kTRUE, const char* outputFileName = "dimuon_analysisMS.root")
{
  // Keep histograms out of the input files' directories
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // --- Initialization & Field Loading ---
  auto& param = NA6PLayoutParam::Instance();
  try {
    param.updateFromFile(Form("%s/na6pLayout.ini", dirSimu), "", true);
  } catch (const std::exception& error) {
    std::fprintf(stderr, "[plotMSDimuons] ERROR: Failed to load configuration: %s\n", error.what());
    return;
  } catch (...) {
    std::fprintf(stderr, "[plotMSDimuons] ERROR: Failed to load configuration: unknown error\n");
    return;
  }

  auto magField = new MagneticField();
  magField->loadField();
  magField->setAsGlobalField();

  NA6PFastTrackFitter fitter;
  const bool hasGeometry = fitter.loadGeometry(Form("%s/geometry.root", dirSimu), "NA6P");
  if (hasGeometry) {
    fitter.enableMaterialCorrections();
  } else {
    fitter.disableMaterialCorrections();
    printf("[plotMSDimuons] WARNING: geometry not loaded from %s/geometry.root, disabling material corrections.\n", dirSimu);
  }

  // --- Event loading ---
  NA6PEventReader reader(Form("%s/VerticesVerTel.root", dirSimu),
                         Form("%s/TracksVerTel.root", dirSimu),
                         Form("%s/TracksMuonSpec.root", dirSimu),
                         Form("%s/TracksMatching.root", dirSimu),
                         Form("%s/MCKine.root", dirSimu),
                         false, false, true, false, true);
  if (!reader.hasTracksMuonSpec() || !reader.hasMSTrackMCLabels() || !reader.hasMC()) {
    printf("[plotMSDimuons] ERROR: required event data are unavailable in %s\n", dirSimu);
    return;
  }

  // Note: Hits tree is opened but not used in pairing logic
  TFile* fh = TFile::Open(Form("%s/HitsMuonSpecModular.root", dirSimu));
  TTree* th = (fh) ? (TTree*)fh->Get("hitsMuonSpecModular") : nullptr;
  std::vector<NA6PMuonSpecModularHit>* hitArr = nullptr;
  if (!th || !th->GetBranch("MuonSpecModular")) {
    printf("[plotMSDimuons] ERROR: required hit tree/branch not found in %s/HitsMuonSpecModular.root\n", dirSimu);
    return;
  }
  th->SetBranchAddress("MuonSpecModular", &hitArr);

  // The track-level MC label only says that at least one associated cluster
  // has a different truth label.  Read the cluster truth container to locate
  // the offending MS plane(s).
  TFile* fClusters = TFile::Open(Form("%s/ClustersMuonSpec.root", dirSimu));
  TTree* tClusters = (fClusters) ? (TTree*)fClusters->Get("clustersMuonSpec") : nullptr;
  NA6PMCTruthContainer* msClusterMCTruth = nullptr;
  if (!tClusters || !tClusters->GetBranch("MuonSpecMCTruth")) {
    printf("[plotMSDimuons] ERROR: required cluster-truth tree/branch not found in %s/ClustersMuonSpec.root\n", dirSimu);
    return;
  }
  tClusters->SetBranchAddress("MuonSpecMCTruth", &msClusterMCTruth);

  constexpr double muonMass = 0.1056583755; // GeV/c^2
  constexpr double maxMuonP = 30.; // GeV/c
  const double yCMShift = NA6PBeamParam::Instance().getYCM();
  float mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  float rangeMMin = mass - 0.15 * 3;
  float rangeMMax = mass + 0.15 * 3;
  if (pdg == 443) { // J/psi
    rangeMMin = mass - 1;
    rangeMMax = mass + 1;
  }

  // --- Histogram Definitions ---
  TH1D* hdimuy = new TH1D("hdimuy", "dimuon y", 100, 0., 5.);
  TH1D* hdimupt = new TH1D("hdimupt", "dimuon pT", 100, 0., 5.);
  TH1D* hdimumassUnconstrained = new TH1D("hdimumassUnconstrained", "Dimuon mass; m_{#mu#mu} (GeV/c^{2});Normalized counts", 200, rangeMMin, rangeMMax);
  TH1D* hdimumassConstrained = new TH1D("hdimumassConstrained", "Dimuon mass; m_{#mu#mu} (GeV/c^{2});Normalized counts", 200, rangeMMin, rangeMMax);
  TH1D* hdimumassUnconstrainedGood = new TH1D("hdimumassUnconstrainedGood", "Good-label pairs; m_{#mu#mu} (GeV/c^{2});Counts", 200, rangeMMin, rangeMMax);
  TH1D* hdimumassConstrainedGood = new TH1D("hdimumassConstrainedGood", "Good-label pairs; m_{#mu#mu} (GeV/c^{2});Counts", 200, rangeMMin, rangeMMax);
  TH1D* hdimumassUnconstrainedFake = new TH1D("hdimumassUnconstrainedFake", "Pairs containing fake labels; m_{#mu#mu} (GeV/c^{2});Counts", 200, rangeMMin, rangeMMax);
  TH1D* hdimumassConstrainedFake = new TH1D("hdimumassConstrainedFake", "Pairs containing fake labels; m_{#mu#mu} (GeV/c^{2});Counts", 200, rangeMMin, rangeMMax);
  TH1D* hMSTrackFakeHitCategory = new TH1D("hMSTrackFakeHitCategory",
                                            "MS-track MC-truth category;Category;Tracks", 9, 0.5, 9.5);
  hMSTrackFakeHitCategory->GetXaxis()->SetBinLabel(1, "true");
  for (int iMS = 0; iMS < 6; ++iMS) {
    hMSTrackFakeHitCategory->GetXaxis()->SetBinLabel(iMS + 2, Form("fake: MS%d", iMS));
  }
  hMSTrackFakeHitCategory->GetXaxis()->SetBinLabel(8, "fake: multiple");
  hMSTrackFakeHitCategory->GetXaxis()->SetBinLabel(9, "fake: unavailable");
  TH1D* hdimuygen = new TH1D("hdimuygen", "dimuon y gen", 100, 0., 5.);
  TH1D* hdimuptgen = new TH1D("hdimuptgen", "dimuon pT gen", 100, 0., 5.);
  TH1D* hdimumassgen = new TH1D("hdimumassgen", "dimuon m gen", 100, 0., 5.);

  TH2D* hdimuyrecvsgen = new TH2D("hdimuyrecvsgen", "dimuon y rec vs gen", 100, 1., 5., 100, 1., 5.);
  TH2D* hdimuptrecvsgen = new TH2D("hdimuptrecvsgen", "dimuon pt rec vs gen", 100, 0., 5., 100, 0., 5.);
  TH1D* hdimudeltayrecvsgen = new TH1D("hdimudeltayrecvsgen", "dimuon y rec - gen", 100, -1., 1.);
  TH1D* hdimudeltaptrecvsgen = new TH1D("hdimudeltaptrecvsgen", "dimuon pt rec - gen", 100, -6., 6.);
  TH2D* hdimurecmassvspt = new TH2D("hdimurecmassvspt", "dimuon rec mass vs pt", 200, 0., 2.7, 100, 0., 5.);
  TH2D* hdimurecmassvseta = new TH2D("hdimurecmassvseta", "dimuon rec mass vs eta", 200, 0., 2.7, 100, 1., 5.);
  TH1D* hmassvsetafit = new TH1D("hmassvsetafit", "dimuon rec mass", 200, 0., 2.7);
  TH1D* hmassvsetamean = new TH1D("hmassvsetamean", "dimuon rec mass mean", 200, 0., 2.7);

  TH2D* hdimurecptvseta = new TH2D("hdimurecptvseta", "dimuon rec pt vs eta", 100, 0., 5., 100, 1., 5.);
  TH2D* hdimurecpvseta = new TH2D("hdimurecpvseta", "dimuon rec p vs eta", 100, 0., 30., 100, 1., 5.);

  constexpr int nMassResolutionPBins = 15;
  constexpr int nMassResolutionYBins = 12;
  TH2D* hMassVsPUnconstrained = new TH2D("hMassVsPUnconstrained", "Unconstrained;m_{#mu#mu}^{rec} (GeV/c^{2});p_{#mu#mu}^{MC} (GeV/c)",
                                         200, rangeMMin, rangeMMax, nMassResolutionPBins, 10., 30.);
  TH2D* hMassVsPConstrained = new TH2D("hMassVsPConstrained", "Constrained;m_{#mu#mu}^{rec} (GeV/c^{2});p_{#mu#mu}^{MC} (GeV/c)",
                                       200, rangeMMin, rangeMMax, nMassResolutionPBins, 10., 30.);
  TH2D* hMassVsYUnconstrained = new TH2D("hMassVsYUnconstrained", "Unconstrained;m_{#mu#mu}^{rec} (GeV/c^{2});y_{#mu#mu}^{MC}",
                                         200, rangeMMin, rangeMMax, nMassResolutionYBins, 2., 2.7);
  TH2D* hMassVsYConstrained = new TH2D("hMassVsYConstrained", "Constrained;m_{#mu#mu}^{rec} (GeV/c^{2});y_{#mu#mu}^{MC}",
                                       200, rangeMMin, rangeMMax, nMassResolutionYBins, 2., 2.7);

  // Single-muon rapidity (reco numerator, MC denominator), split by charge
  TH1D* hRapidity = new TH1D("hRapidity", "muon rapidity", 100, 1., 7.);
  TH1D* hRapidityNeg = new TH1D("hRapidityNeg", "muon rapidity", 100, 1., 7.);
  TH1D* hRapidityPos = new TH1D("hRapidityPos", "muon rapidity", 100, 1., 7.);
  TH1D* hRapidityMC = new TH1D("hRapidityMC", "muon rapidity", 100, 1., 7.);
  TH1D* hRapidityMCNeg = new TH1D("hRapidityMCNeg", "muon rapidity", 100, 1., 7.);
  TH1D* hRapidityMCPos = new TH1D("hRapidityMCPos", "muon rapidity", 100, 1., 7.);
  TH1D* hRapidityReconstructable = new TH1D("hRapidityReconstructable", "Reconstructable muons;MC rapidity;Counts", 100, 1., 7.);
  TH1D* hRapidityReconstructableNeg = new TH1D("hRapidityReconstructableNeg", "Reconstructable #mu^{-};MC rapidity;Counts", 100, 1., 7.);
  TH1D* hRapidityReconstructablePos = new TH1D("hRapidityReconstructablePos", "Reconstructable #mu^{+};MC rapidity;Counts", 100, 1., 7.);
  TH1D* hRapidityReconstructedMC = new TH1D("hRapidityReconstructedMC", "Reconstructed muons;MC rapidity;Counts", 100, 1., 7.);
  TH1D* hRapidityReconstructedMCNeg = new TH1D("hRapidityReconstructedMCNeg", "Reconstructed #mu^{-};MC rapidity;Counts", 100, 1., 7.);
  TH1D* hRapidityReconstructedMCPos = new TH1D("hRapidityReconstructedMCPos", "Reconstructed #mu^{+};MC rapidity;Counts", 100, 1., 7.);

  TH1D* hMomentumMC = new TH1D("hMomentumMC", "Generated muons;MC p (GeV/c);Counts", 100, 0., maxMuonP);
  TH1D* hMomentumMCNeg = new TH1D("hMomentumMCNeg", "Generated #mu^{-};MC p (GeV/c);Counts", 100, 0., maxMuonP);
  TH1D* hMomentumMCPos = new TH1D("hMomentumMCPos", "Generated #mu^{+};MC p (GeV/c);Counts", 100, 0., maxMuonP);
  TH1D* hMomentumReconstructable = new TH1D("hMomentumReconstructable", "Reconstructable muons;MC p (GeV/c);Counts", 100, 0., maxMuonP);
  TH1D* hMomentumReconstructableNeg = new TH1D("hMomentumReconstructableNeg", "Reconstructable #mu^{-};MC p (GeV/c);Counts", 100, 0., maxMuonP);
  TH1D* hMomentumReconstructablePos = new TH1D("hMomentumReconstructablePos", "Reconstructable #mu^{+};MC p (GeV/c);Counts", 100, 0., maxMuonP);
  TH1D* hMomentumReconstructedMC = new TH1D("hMomentumReconstructedMC", "Reconstructed muons;MC p (GeV/c);Counts", 100, 0., maxMuonP);
  TH1D* hMomentumReconstructedMCNeg = new TH1D("hMomentumReconstructedMCNeg", "Reconstructed #mu^{-};MC p (GeV/c);Counts", 100, 0., maxMuonP);
  TH1D* hMomentumReconstructedMCPos = new TH1D("hMomentumReconstructedMCPos", "Reconstructed #mu^{+};MC p (GeV/c);Counts", 100, 0., maxMuonP);

  // Single-muon momentum residuals at the production vertex, split by MC charge
  // and reconstructed with or without the vertex constraint.
  TH1D* hDeltaPxNegUnconstrained = new TH1D("hDeltaPxNegUnconstrained", "#mu^{-};p_{x}^{rec}-p_{x}^{MC} (GeV/c);Normalized counts", 200, -1., 1.);
  TH1D* hDeltaPyNegUnconstrained = new TH1D("hDeltaPyNegUnconstrained", "#mu^{-};p_{y}^{rec}-p_{y}^{MC} (GeV/c);Normalized counts", 200, -1., 1.);
  TH1D* hDeltaPzNegUnconstrained = new TH1D("hDeltaPzNegUnconstrained", "#mu^{-};p_{z}^{rec}-p_{z}^{MC} (GeV/c);Normalized counts", 200, -5., 5.);
  TH1D* hDeltaPxPosUnconstrained = new TH1D("hDeltaPxPosUnconstrained", "#mu^{+};p_{x}^{rec}-p_{x}^{MC} (GeV/c);Normalized counts", 200, -1., 1.);
  TH1D* hDeltaPyPosUnconstrained = new TH1D("hDeltaPyPosUnconstrained", "#mu^{+};p_{y}^{rec}-p_{y}^{MC} (GeV/c);Normalized counts", 200, -1., 1.);
  TH1D* hDeltaPzPosUnconstrained = new TH1D("hDeltaPzPosUnconstrained", "#mu^{+};p_{z}^{rec}-p_{z}^{MC} (GeV/c);Normalized counts", 200, -5., 5.);
  TH1D* hDeltaPxNegConstrained = new TH1D("hDeltaPxNegConstrained", "#mu^{-};p_{x}^{rec}-p_{x}^{MC} (GeV/c);Normalized counts", 200, -1., 1.);
  TH1D* hDeltaPyNegConstrained = new TH1D("hDeltaPyNegConstrained", "#mu^{-};p_{y}^{rec}-p_{y}^{MC} (GeV/c);Normalized counts", 200, -1., 1.);
  TH1D* hDeltaPzNegConstrained = new TH1D("hDeltaPzNegConstrained", "#mu^{-};p_{z}^{rec}-p_{z}^{MC} (GeV/c);Normalized counts", 200, -5., 5.);
  TH1D* hDeltaPxPosConstrained = new TH1D("hDeltaPxPosConstrained", "#mu^{+};p_{x}^{rec}-p_{x}^{MC} (GeV/c);Normalized counts", 200, -1., 1.);
  TH1D* hDeltaPyPosConstrained = new TH1D("hDeltaPyPosConstrained", "#mu^{+};p_{y}^{rec}-p_{y}^{MC} (GeV/c);Normalized counts", 200, -1., 1.);
  TH1D* hDeltaPzPosConstrained = new TH1D("hDeltaPzPosConstrained", "#mu^{+};p_{z}^{rec}-p_{z}^{MC} (GeV/c);Normalized counts", 200, -5., 5.);

  // Pulls of the five native track parameters at the production vertex:
  // {x, y, tx, ty, q/pxz}. The covariance is the diagonal element of the
  // propagated fitted-state covariance matrix.
  const char* parNames[5] = {"x", "y", "t_{x}", "t_{y}", "q/p_{XZ}"};
  const char* parAxes[5] = {"x (cm)", "y (cm)", "t_{x}", "t_{y}", "q/p_{XZ} ((GeV/c)^{-1})"};
  TH1D* hTrackPullsUnconstrained[5] = {};
  TH1D* hTrackPullsConstrained[5] = {};
  for (int ip = 0; ip < 5; ++ip) {
    hTrackPullsUnconstrained[ip] = new TH1D(
      Form("hTrackPull_%s_Unconstrained", ip == 0 ? "X" : ip == 1 ? "Y" : ip == 2 ? "Tx" : ip == 3 ? "Ty" : "Q2Pxz"),
      Form("Unconstrained %s pull;(%s^{rec}-%s^{MC})/#sigma_{%s};Tracks", parNames[ip], parAxes[ip], parAxes[ip], parNames[ip]),
      160, -8., 8.);
    hTrackPullsConstrained[ip] = new TH1D(
      Form("hTrackPull_%s_Constrained", ip == 0 ? "X" : ip == 1 ? "Y" : ip == 2 ? "Tx" : ip == 3 ? "Ty" : "Q2Pxz"),
      Form("Constrained %s pull;(%s^{rec}-%s^{MC})/#sigma_{%s};Tracks", parNames[ip], parAxes[ip], parAxes[ip], parNames[ip]),
      160, -8., 8.);
  }

  long long nPairCandidates = 0;
  long long nPairFilled = 0;
  long long nMSTracks = 0;
  long long nPassNHits = 0;
  long long nPassValidMCLabel = 0;
  long long nPassMuonPDG = 0;
  long long nPassClusterMap = 0;
  long long nPassMotherPDG = 0;
  long long nEventsWithTwoFilteredTracks = 0;

  const std::int64_t readerEntries = reader.entries();
  const std::int64_t clusterEntries = tClusters->GetEntries();
  const std::int64_t nEv = readerEntries < clusterEntries ? readerEntries : clusterEntries;
  if (clusterEntries != readerEntries) {
    printf("[plotMSDimuons] WARNING: cluster-truth tree has %lld entries, reader has %lld; using %lld events\n",
           static_cast<long long>(clusterEntries), static_cast<long long>(readerEntries), static_cast<long long>(nEv));
  }
  printf("Number of events = %lld\n", static_cast<long long>(nEv));

  for (int jEv = 0; jEv < nEv; jEv++) {
    printf("Processing event %d/%lld\r", jEv + 1, static_cast<long long>(nEv));
    if (!reader.loadEvent(jEv)) continue;
    tClusters->GetEntry(jEv);
    const auto& trArr = reader.tracksMuonSpec();
    const auto& trMCLabels = reader.trackLabelsMuonSpec();
    const auto& mcArr = reader.mcParticles();
    double primaryVertex[3] = {0., 0., 0.};
    if (const auto* header = reader.mcHeader()) {
      primaryVertex[0] = header->getVX();
      primaryVertex[1] = header->getVY();
      primaryVertex[2] = header->getVZ();
    }

    th->GetEvent(jEv);

    int nPart = mcArr.size();
    int nTracks = trArr.size();
    nMSTracks += nTracks;
    if (trMCLabels.size() != static_cast<size_t>(nTracks)) {
      printf("[plotMSDimuons] ERROR: track/MC-label size mismatch in event %d\n", jEv);
      return;
    }

    // Categorize every reconstructed MS track with a valid MC label.  The
    // reconstruction marks a track fake when its dominant MC label is absent
    // from one or more associated clusters.  Reproduce that criterion plane
    // by plane to identify the single wrong hit, or the multiple-hit case.
    for (int jTr = 0; jTr < nTracks; ++jTr) {
      const NA6PTrack& tr = trArr.at(jTr);
      const auto& trackLabel = trMCLabels.at(jTr);
      if (!trackLabel.isValid()) continue;
      if (!trackLabel.isFake()) {
        hMSTrackFakeHitCategory->Fill(1.);
        continue;
      }

      int nWrongMSHits = 0;
      int wrongMSLayer = -1;
      bool unavailable = msClusterMCTruth == nullptr;
      for (int iMS = 0; iMS < 6; ++iMS) {
        const int trackLayer = param.nVerTelPlanes + iMS;
        if ((tr.getClusterMap() & (1u << trackLayer)) == 0u) {
          unavailable = true;
          continue;
        }
        const int clusterIndex = tr.getClusterIndex(trackLayer);
        if (clusterIndex < 0) {
          unavailable = true;
          continue;
        }
        const auto labels = msClusterMCTruth->getLabels(clusterIndex);
        bool matchesTrackLabel = false;
        for (const auto& clusterLabel : labels) {
          if (clusterLabel == trackLabel) {
            matchesTrackLabel = true;
            break;
          }
        }
        if (!matchesTrackLabel) {
          ++nWrongMSHits;
          wrongMSLayer = iMS;
        }
      }

      if (unavailable || nWrongMSHits == 0) {
        hMSTrackFakeHitCategory->Fill(9.);
      } else if (nWrongMSHits == 1) {
        hMSTrackFakeHitCategory->Fill(wrongMSLayer + 2.);
      } else {
        hMSTrackFakeHitCategory->Fill(8.);
      }
    }

    // Build the MC-track hit mask in physical MS-plane numbering (0 ... nMSPlanes-1).
    std::vector<uint32_t> mcHitMask(nPart, 0u);
    if (hitArr) {
      for (const auto& hit : *hitArr) {
        const int mcTrackID = hit.getTrackID();
        if (mcTrackID < 0 || mcTrackID >= nPart || param.nMSPlanes <= 0) continue;
        const double hitZ = hit.getZ();
        int closestPlane = 0;
        double minDistance = TMath::Abs(hitZ - (param.shiftMS[2] + param.posMSPlaneZ[0]));
        for (int iPlane = 1; iPlane < param.nMSPlanes; ++iPlane) {
          const double distance = TMath::Abs(hitZ - (param.shiftMS[2] + param.posMSPlaneZ[iPlane]));
          if (distance < minDistance) {
            minDistance = distance;
            closestPlane = iPlane;
          }
        }
        if (closestPlane < 32) mcHitMask[mcTrackID] |= (1u << closestPlane);
      }
    }
    const int nRequiredPlanes = TMath::Max(0, TMath::Min(NRecoClusters, param.nMSPlanes));
    const uint32_t requiredHitMask = nRequiredPlanes >= 32 ? 0xffffffffu : ((1u << nRequiredPlanes) - 1u);
    std::vector<bool> isReconstructable(nPart, false);
    std::vector<bool> isReconstructed(nPart, false);

    const double zOrig = primaryVertex[2];

    // Loop over MC kine: muons from the requested mother only (efficiency denominator)
    for (int jp = 0; jp < nPart; jp++) {
      auto& curPart = mcArr.at(jp);
      if (TMath::Abs(curPart.GetPdgCode()) != 13) continue;
      const int iMother = curPart.GetFirstMother();
      if (iMother < 0 || iMother >= nPart) continue;
      if (mcArr.at(iMother).GetPdgCode() != pdg) continue;

      const double y = curPart.Y();
      const double p = curPart.P();
      hRapidityMC->Fill(y);
      hMomentumMC->Fill(p);
      if (curPart.GetPdgCode() == 13) hRapidityMCNeg->Fill(y);  // mu-
      if (curPart.GetPdgCode() == -13) hRapidityMCPos->Fill(y); // mu+
      if (curPart.GetPdgCode() == 13) hMomentumMCNeg->Fill(p);
      if (curPart.GetPdgCode() == -13) hMomentumMCPos->Fill(p);
      isReconstructable[jp] = requiredHitMask != 0u && (mcHitMask[jp] & requiredHitMask) == requiredHitMask;
      if (isReconstructable[jp]) {
        hRapidityReconstructable->Fill(y);
        hMomentumReconstructable->Fill(p);
        if (curPart.GetPdgCode() == 13) hRapidityReconstructableNeg->Fill(y);
        if (curPart.GetPdgCode() == -13) hRapidityReconstructablePos->Fill(y);
        if (curPart.GetPdgCode() == 13) hMomentumReconstructableNeg->Fill(p);
        if (curPart.GetPdgCode() == -13) hMomentumReconstructablePos->Fill(p);
      }
    }

    // Pre-filter tracks based on ClusterMap
    std::vector<int> filteredIndices;
    for (int jTr = 0; jTr < nTracks; ++jTr) {
      const NA6PTrack& tr = trArr.at(jTr);
      if (tr.getNHits() < NRecoClusters) continue;
      ++nPassNHits;
      const auto& mcLabel = trMCLabels.at(jTr);
      const int mcTrackID = mcLabel.getTrackID();
      if (!mcLabel.isValid() || mcTrackID < 0 || mcTrackID >= nPart) continue;
      ++nPassValidMCLabel;
      if (TMath::Abs(mcArr.at(mcTrackID).GetPdgCode()) != 13) continue;
      ++nPassMuonPDG;

      uint32_t clustermap = tr.getClusterMap();
      uint32_t maskClusters = (NRecoClusters == 4) ? (((1u << 4) - 1) << 5) : (((1u << 6) - 1) << 5);
      if ((clustermap & maskClusters) != maskClusters) continue;
      ++nPassClusterMap;

      // Same mother selection as the denominator
      const int iMother = mcArr.at(mcTrackID).GetFirstMother();
      if (iMother < 0 || iMother >= nPart || mcArr.at(iMother).GetPdgCode() != pdg) continue;
      ++nPassMotherPDG;

      // Momentum resolutions at the production vertex for both track states.
      NA6PTrackParCov inward = tr.getInwardParam();
      const bool hasInward = Propagator::Instance()->propagateToZ(inward, zOrig, fitter.getPropOpt());
      NA6PTrackParCov constrained = tr.getVertexConstrainedParam();
      const bool constrainedAtVertex = TMath::Abs(constrained.getZ() - zOrig) < 1.e-5;
      const bool hasConstrained = tr.getStatusConstrained() &&
                                  (constrainedAtVertex || Propagator::Instance()->propagateToZ(constrained, zOrig, fitter.getPropOpt()));
      const TParticle& mcMuon = mcArr.at(mcTrackID);
      const double mcPxz = std::hypot(mcMuon.Px(), mcMuon.Pz());
      const double mcCharge = mcMuon.GetPdgCode() == 13 ? -1. : 1.;
      const double mcPar[5] = {
        primaryVertex[0], primaryVertex[1],
        mcPxz > 0. ? mcMuon.Px() / mcPxz : 0.,
        mcPxz > 0. ? mcMuon.Py() / mcPxz : 0.,
        mcPxz > 0. ? mcCharge / mcPxz : 0.};
      auto fillMomentumResiduals = [&](const NA6PTrackParCov& state, bool isConstrained) {
        const auto recMom = state.getPXYZ();
        TH1D* hPx = nullptr;
        TH1D* hPy = nullptr;
        TH1D* hPz = nullptr;
        if (mcMuon.GetPdgCode() == 13) {
          hPx = isConstrained ? hDeltaPxNegConstrained : hDeltaPxNegUnconstrained;
          hPy = isConstrained ? hDeltaPyNegConstrained : hDeltaPyNegUnconstrained;
          hPz = isConstrained ? hDeltaPzNegConstrained : hDeltaPzNegUnconstrained;
        } else {
          hPx = isConstrained ? hDeltaPxPosConstrained : hDeltaPxPosUnconstrained;
          hPy = isConstrained ? hDeltaPyPosConstrained : hDeltaPyPosUnconstrained;
          hPz = isConstrained ? hDeltaPzPosConstrained : hDeltaPzPosUnconstrained;
        }
        hPx->Fill(recMom[0] - mcMuon.Px());
        hPy->Fill(recMom[1] - mcMuon.Py());
        hPz->Fill(recMom[2] - mcMuon.Pz());
      };
      auto fillTrackPulls = [&](const NA6PTrackParCov& state, bool isConstrained) {
        TH1D** histograms = isConstrained ? hTrackPullsConstrained : hTrackPullsUnconstrained;
        for (int ip = 0; ip < 5; ++ip) {
          const double variance = state.getCovMatElem(ip, ip);
          if (variance > 0. && std::isfinite(variance)) {
            const double pull = (state.getParam(ip) - mcPar[ip]) / std::sqrt(variance);
            if (std::isfinite(pull)) {
              histograms[ip]->Fill(pull);
            }
          }
        }
      };
      // Fill both from the same track sample so their normalized shapes are directly comparable.
      if (hasInward && hasConstrained) {
        fillMomentumResiduals(inward, false);
        fillMomentumResiduals(constrained, true);
        fillTrackPulls(inward, false);
        fillTrackPulls(constrained, true);
      }

      // Rapidity with the muon mass hypothesis
      const auto pxyz = tr.getPXYZ();
      const double p = TMath::Sqrt(pxyz[0] * pxyz[0] + pxyz[1] * pxyz[1] + pxyz[2] * pxyz[2]);
      const double E = TMath::Sqrt(p * p + muonMass * muonMass);
      const double rapidity = 0.5 * TMath::Log((E + pxyz[2]) / (E - pxyz[2]));

      hRapidity->Fill(rapidity);
      if (tr.getCharge() < 0) hRapidityNeg->Fill(rapidity);
      if (tr.getCharge() > 0) hRapidityPos->Fill(rapidity);
      isReconstructed[mcTrackID] = true;
      filteredIndices.push_back(jTr);
    }

    // Fill the reconstructed numerator once per MC muon, using MC rapidity to
    // avoid migration between numerator and denominator bins.
    for (int jp = 0; jp < nPart; ++jp) {
      if (!isReconstructed[jp] || !isReconstructable[jp]) continue;
      const auto& curPart = mcArr.at(jp);
      if (TMath::Abs(curPart.GetPdgCode()) != 13) continue;
      const int iMother = curPart.GetFirstMother();
      if (iMother < 0 || iMother >= nPart || mcArr.at(iMother).GetPdgCode() != pdg) continue;
      hRapidityReconstructedMC->Fill(curPart.Y());
      hMomentumReconstructedMC->Fill(curPart.P());
      if (curPart.GetPdgCode() == 13) hRapidityReconstructedMCNeg->Fill(curPart.Y());
      if (curPart.GetPdgCode() == -13) hRapidityReconstructedMCPos->Fill(curPart.Y());
      if (curPart.GetPdgCode() == 13) hMomentumReconstructedMCNeg->Fill(curPart.P());
      if (curPart.GetPdgCode() == -13) hMomentumReconstructedMCPos->Fill(curPart.P());
    }

    // --- Dimuon Pairing ---
    if (filteredIndices.size() < 2) continue;
    ++nEventsWithTwoFilteredTracks;
    nPairCandidates += static_cast<long long>(filteredIndices.size()) * static_cast<long long>(filteredIndices.size() - 1) / 2;

    for (size_t i = 0; i < filteredIndices.size(); ++i) {
      NA6PTrack tr1 = trArr.at(filteredIndices[i]);
      NA6PTrackParCov inward1 = tr1.getInwardParam();
      const bool hasInward1 = Propagator::Instance()->propagateToZ(inward1, zOrig, fitter.getPropOpt());
      if (!hasInward1) continue;
      NA6PTrackParCov constrained1 = tr1.getVertexConstrainedParam();
      const bool constrained1AtVertex = TMath::Abs(constrained1.getZ() - zOrig) < 1.e-5;
      const bool hasConstrained1 = tr1.getStatusConstrained() &&
                                   (constrained1AtVertex || Propagator::Instance()->propagateToZ(constrained1, zOrig, fitter.getPropOpt()));
      if (!hasConstrained1) continue;

      const auto recMomInward1 = inward1.getPXYZ();
      const auto recMomConstrained1 = constrained1.getPXYZ();
      TLorentzVector decprodInward0;
      decprodInward0.SetPxPyPzE(recMomInward1[0], recMomInward1[1], recMomInward1[2],
                               TMath::Sqrt(inward1.getP() * inward1.getP() + muonMass * muonMass));
      TLorentzVector decprodConstrained0;
      decprodConstrained0.SetPxPyPzE(recMomConstrained1[0], recMomConstrained1[1], recMomConstrained1[2],
                                    TMath::Sqrt(constrained1.getP() * constrained1.getP() + muonMass * muonMass));

      // MC truth for track 1 (validity already ensured by the pre-filter)
      const int mcLabel1 = trMCLabels.at(filteredIndices[i]).getTrackID();
      const TParticle& mu1 = mcArr.at(mcLabel1);
      TLorentzVector decprodgen0;
      decprodgen0.SetPxPyPzE(mu1.Px(), mu1.Py(), mu1.Pz(), mu1.Energy());

      for (size_t j = i + 1; j < filteredIndices.size(); ++j) {
        NA6PTrack tr2 = trArr.at(filteredIndices[j]);
        const int mcLabel2 = trMCLabels.at(filteredIndices[j]).getTrackID();

        if (tr1.getCharge() * tr2.getCharge() >= 0) continue;
        if (mcLabel1 == mcLabel2) continue;

        const TParticle& mu2 = mcArr.at(mcLabel2);
        const int iMother = mu1.GetFirstMother();
        if (iMother < 0 || iMother >= nPart || iMother != mu2.GetFirstMother()) continue;
        if (mcArr.at(iMother).GetPdgCode() != pdg) continue;

        NA6PTrackParCov inward2 = tr2.getInwardParam();
        const bool hasInward2 = Propagator::Instance()->propagateToZ(inward2, zOrig, fitter.getPropOpt());
        if (!hasInward2) continue;
        NA6PTrackParCov constrained2 = tr2.getVertexConstrainedParam();
        const bool constrained2AtVertex = TMath::Abs(constrained2.getZ() - zOrig) < 1.e-5;
        const bool hasConstrained2 = tr2.getStatusConstrained() &&
                                     (constrained2AtVertex || Propagator::Instance()->propagateToZ(constrained2, zOrig, fitter.getPropOpt()));
        if (!hasConstrained2) continue;

        const auto recMomInward2 = inward2.getPXYZ();
        const auto recMomConstrained2 = constrained2.getPXYZ();
        TLorentzVector decprodInward1;
        decprodInward1.SetPxPyPzE(recMomInward2[0], recMomInward2[1], recMomInward2[2],
                                 TMath::Sqrt(inward2.getP() * inward2.getP() + muonMass * muonMass));
        TLorentzVector decprodConstrained1;
        decprodConstrained1.SetPxPyPzE(recMomConstrained2[0], recMomConstrained2[1], recMomConstrained2[2],
                                      TMath::Sqrt(constrained2.getP() * constrained2.getP() + muonMass * muonMass));

        const TLorentzVector mothUnconstrained = decprodInward0 + decprodInward1;
        const TLorentzVector mothConstrained = decprodConstrained0 + decprodConstrained1;
        hdimumassUnconstrained->Fill(mothUnconstrained.M());
        hdimumassConstrained->Fill(mothConstrained.M());
        const bool pairIsFake = trMCLabels.at(filteredIndices[i]).isFake() ||
                                trMCLabels.at(filteredIndices[j]).isFake();
        if (pairIsFake) {
          hdimumassUnconstrainedFake->Fill(mothUnconstrained.M());
          hdimumassConstrainedFake->Fill(mothConstrained.M());
        } else {
          hdimumassUnconstrainedGood->Fill(mothUnconstrained.M());
          hdimumassConstrainedGood->Fill(mothConstrained.M());
        }
        const TLorentzVector& moth = VertexConstrained ? mothConstrained : mothUnconstrained;
        hdimuy->Fill(moth.Rapidity());
        hdimupt->Fill(moth.Pt());
        ++nPairFilled;
        hdimurecmassvspt->Fill(moth.M(), moth.Pt());
        hdimurecmassvseta->Fill(moth.M(), moth.Rapidity());
        hdimurecptvseta->Fill(moth.Pt(), moth.Rapidity());
        hdimurecpvseta->Fill(moth.P(), moth.Rapidity());

        TLorentzVector decprodgen1;
        decprodgen1.SetPxPyPzE(mu2.Px(), mu2.Py(), mu2.Pz(), mu2.Energy());
        TLorentzVector mothgen = decprodgen0 + decprodgen1;

        hMassVsPUnconstrained->Fill(mothUnconstrained.M(), mothgen.P());
        hMassVsPConstrained->Fill(mothConstrained.M(), mothgen.P());
        hMassVsYUnconstrained->Fill(mothUnconstrained.M(), mothgen.Rapidity());
        hMassVsYConstrained->Fill(mothConstrained.M(), mothgen.Rapidity());

        hdimumassgen->Fill(mothgen.M());
        hdimuygen->Fill(mothgen.Rapidity());
        hdimuptgen->Fill(mothgen.Pt());
        hdimuyrecvsgen->Fill(moth.Rapidity(), mothgen.Rapidity());
        hdimuptrecvsgen->Fill(moth.Pt(), mothgen.Pt());
        hdimudeltayrecvsgen->Fill(moth.Rapidity() - mothgen.Rapidity());
        hdimudeltaptrecvsgen->Fill(moth.Pt() - mothgen.Pt());
      }
    }
  }

  // --- Plotting & Saving ---
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  TCanvas* cDimu = new TCanvas("cDimu", "Dimuon mass", 1200, 1200);
  for (auto* hist : {hdimumassUnconstrained, hdimumassConstrained}) {
    if (hist->GetEntries() > 10 && hist->GetRMS() > 0.) {
      hist->Fit("gaus", "Q0", "", mass - 1.5 * hist->GetRMS(), mass + 1.5 * hist->GetRMS());
    }
  }
  printf("[plotMSDimuons] Pair summary: candidates=%lld comparison sample=%lld\n", nPairCandidates, nPairFilled);
  printf("[plotMSDimuons] MS-track prefilter (NRecoClusters=%d): total=%lld, "
         "nHits=%lld, valid MC label=%lld, muons=%lld, full cluster map=%lld, "
         "J/#psi daughters=%lld, events with >=2 selected tracks=%lld\n",
         NRecoClusters, nMSTracks, nPassNHits, nPassValidMCLabel, nPassMuonPDG,
         nPassClusterMap, nPassMotherPDG, nEventsWithTwoFilteredTracks);
  auto* hMassUnconstrainedDraw = static_cast<TH1D*>(hdimumassUnconstrained->Clone("hdimumassUnconstrained_draw"));
  auto* hMassConstrainedDraw = static_cast<TH1D*>(hdimumassConstrained->Clone("hdimumassConstrained_draw"));
  hMassUnconstrainedDraw->SetDirectory(nullptr);
  hMassConstrainedDraw->SetDirectory(nullptr);
  if (hMassUnconstrainedDraw->Integral() > 0.) hMassUnconstrainedDraw->Scale(1. / hMassUnconstrainedDraw->Integral());
  if (hMassConstrainedDraw->Integral() > 0.) hMassConstrainedDraw->Scale(1. / hMassConstrainedDraw->Integral());
  hMassUnconstrainedDraw->SetLineColor(kBlue + 1);
  hMassUnconstrainedDraw->SetLineWidth(2);
  hMassConstrainedDraw->SetLineColor(kRed + 1);
  hMassConstrainedDraw->SetLineWidth(2);
  hMassUnconstrainedDraw->SetMaximum(1.15 * TMath::Max(hMassUnconstrainedDraw->GetMaximum(), hMassConstrainedDraw->GetMaximum()));
  hMassUnconstrainedDraw->Draw("HIST");
  hMassConstrainedDraw->Draw("HIST SAME");
  auto* massLegend = new TLegend(0.52, 0.72, 0.88, 0.88);
  massLegend->SetBorderSize(0);
  const auto* massFitUnconstrained = hdimumassUnconstrained->GetFunction("gaus");
  const auto* massFitConstrained = hdimumassConstrained->GetFunction("gaus");
  massLegend->AddEntry(hMassUnconstrainedDraw,
                       massFitUnconstrained ? Form("Unconstrained, #sigma=%.3g GeV/c^{2}", massFitUnconstrained->GetParameter(2)) : "Unconstrained",
                       "l");
  massLegend->AddEntry(hMassConstrainedDraw,
                       massFitConstrained ? Form("Constrained, #sigma=%.3g GeV/c^{2}", massFitConstrained->GetParameter(2)) : "Constrained",
                       "l");
  massLegend->Draw();
  cDimu->SaveAs(Form("%s/dimuon_mass_spectrumMS_comparison.png", dirSimu));

  TCanvas* cDimuFake = new TCanvas("cDimuFake", "Dimuon mass: fake and good labels", 1200, 1200);
  auto* hGoodDraw = static_cast<TH1D*>(hdimumassConstrainedGood->Clone("hdimumassConstrainedGood_draw"));
  auto* hFakeDraw = static_cast<TH1D*>(hdimumassConstrainedFake->Clone("hdimumassConstrainedFake_draw"));
  hGoodDraw->SetDirectory(nullptr);
  hFakeDraw->SetDirectory(nullptr);
  hGoodDraw->SetStats(kFALSE);
  hFakeDraw->SetStats(kFALSE);
  if (hGoodDraw->Integral() > 0.) hGoodDraw->Scale(1. / hGoodDraw->Integral());
  if (hFakeDraw->Integral() > 0.) hFakeDraw->Scale(1. / hFakeDraw->Integral());
  hGoodDraw->SetLineColor(kGreen + 2);
  hGoodDraw->SetLineWidth(2);
  hFakeDraw->SetLineColor(kMagenta + 1);
  hFakeDraw->SetLineWidth(2);
  if (hdimumassConstrainedGood->GetEntries() > 10 && hdimumassConstrainedGood->GetRMS() > 0.) {
    hdimumassConstrainedGood->Fit("gaus", "Q0", "", mass - 1.5 * hdimumassConstrainedGood->GetRMS(),
                                  mass + 1.5 * hdimumassConstrainedGood->GetRMS());
  }
  if (hdimumassConstrainedFake->GetEntries() > 10 && hdimumassConstrainedFake->GetRMS() > 0.) {
    hdimumassConstrainedFake->Fit("gaus", "Q0", "", mass - 1.5 * hdimumassConstrainedFake->GetRMS(),
                                  mass + 1.5 * hdimumassConstrainedFake->GetRMS());
  }
  hGoodDraw->SetTitle("Vertex-constrained dimuon mass by MC-label quality");
  hGoodDraw->SetMaximum(1.15 * TMath::Max(hGoodDraw->GetMaximum(), hFakeDraw->GetMaximum()));
  hGoodDraw->Draw("HIST");
  hFakeDraw->Draw("HIST SAME");
  auto* fakeLegend = new TLegend(0.52, 0.72, 0.88, 0.88);
  fakeLegend->SetBorderSize(0);
  const auto* goodFit = hdimumassConstrainedGood->GetFunction("gaus");
  const auto* fakeFit = hdimumassConstrainedFake->GetFunction("gaus");
  fakeLegend->AddEntry(hGoodDraw,
                       goodFit ? Form("Good labels only: N=%lld, #sigma_{fit}=%.3g, RMS=%.3g GeV/c^{2}",
                                      static_cast<Long64_t>(hdimumassConstrainedGood->GetEntries()), goodFit->GetParameter(2),
                                      hdimumassConstrainedGood->GetRMS())
                               : Form("Good labels only: N=%lld, RMS=%.3g GeV/c^{2}",
                                      static_cast<Long64_t>(hdimumassConstrainedGood->GetEntries()), hdimumassConstrainedGood->GetRMS()),
                       "l");
  fakeLegend->AddEntry(hFakeDraw,
                       fakeFit ? Form("At least one fake label: N=%lld, #sigma_{fit}=%.3g, RMS=%.3g GeV/c^{2}",
                                      static_cast<Long64_t>(hdimumassConstrainedFake->GetEntries()), fakeFit->GetParameter(2),
                                      hdimumassConstrainedFake->GetRMS())
                               : Form("At least one fake label: N=%lld, RMS=%.3g GeV/c^{2}",
                                      static_cast<Long64_t>(hdimumassConstrainedFake->GetEntries()), hdimumassConstrainedFake->GetRMS()),
                       "l");
  fakeLegend->Draw();
  cDimuFake->SaveAs(Form("%s/dimuon_mass_spectrumMS_fake_vs_good.png", dirSimu));

  TH1D* hMassResolutionVsPUnconstrained = new TH1D("hMassResolutionVsPUnconstrained", "Mass resolution;p_{#mu#mu}^{MC} (GeV/c);#sigma_{m} (MeV/c^{2})",
                                                    nMassResolutionPBins, 10., 30.);
  TH1D* hMassResolutionVsPConstrained = new TH1D("hMassResolutionVsPConstrained", "Mass resolution;p_{#mu#mu}^{MC} (GeV/c);#sigma_{m} (MeV/c^{2})",
                                                  nMassResolutionPBins, 10., 30.);
  TH1D* hMassResolutionVsYUnconstrained = new TH1D("hMassResolutionVsYUnconstrained", "Mass resolution;y_{#mu#mu}^{MC};#sigma_{m} (MeV/c^{2})",
                                                    nMassResolutionYBins, 2., 2.7);
  TH1D* hMassResolutionVsYConstrained = new TH1D("hMassResolutionVsYConstrained", "Mass resolution;y_{#mu#mu}^{MC};#sigma_{m} (MeV/c^{2})",
                                                  nMassResolutionYBins, 2., 2.7);
  fillMassResolutionVsBins(hMassVsPUnconstrained, hMassResolutionVsPUnconstrained, "massProjPUnconstrained", "massFitPUnconstrained");
  fillMassResolutionVsBins(hMassVsPConstrained, hMassResolutionVsPConstrained, "massProjPConstrained", "massFitPConstrained");
  fillMassResolutionVsBins(hMassVsYUnconstrained, hMassResolutionVsYUnconstrained, "massProjYUnconstrained", "massFitYUnconstrained");
  fillMassResolutionVsBins(hMassVsYConstrained, hMassResolutionVsYConstrained, "massProjYConstrained", "massFitYConstrained");

  TCanvas* cMassResolution = new TCanvas("cMassResolution", "Dimuon mass resolution", 1400, 600);
  cMassResolution->Divide(2, 1);
  TH1D* massResolutionUnconstrained[2] = {hMassResolutionVsPUnconstrained, hMassResolutionVsYUnconstrained};
  TH1D* massResolutionConstrained[2] = {hMassResolutionVsPConstrained, hMassResolutionVsYConstrained};
  for (int i = 0; i < 2; ++i) {
    cMassResolution->cd(i + 1);
    auto* unconstrained = massResolutionUnconstrained[i];
    auto* constrained = massResolutionConstrained[i];
    unconstrained->SetMinimum(0.);
    unconstrained->SetMaximum(1.2 * TMath::Max(unconstrained->GetMaximum(), constrained->GetMaximum()));
    unconstrained->SetMarkerStyle(20);
    unconstrained->SetMarkerColor(kBlue + 1);
    unconstrained->SetLineColor(kBlue + 1);
    constrained->SetMarkerStyle(24);
    constrained->SetMarkerColor(kRed + 1);
    constrained->SetLineColor(kRed + 1);
    unconstrained->Draw("EP");
    constrained->Draw("EP SAME");
    auto* legend = new TLegend(0.58, 0.75, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->AddEntry(unconstrained, "Unconstrained", "lp");
    legend->AddEntry(constrained, "Constrained", "lp");
    legend->Draw();
  }
  cMassResolution->SaveAs(Form("%s/dimuon_mass_resolution_vs_p_yMS.png", dirSimu));

  TCanvas* cMomRes = new TCanvas("cMomRes", "Single-muon momentum resolution", 1500, 900);
  cMomRes->Divide(3, 2);
  TH1D* momentumResidualsUnconstrained[6] = {
    hDeltaPxNegUnconstrained, hDeltaPyNegUnconstrained, hDeltaPzNegUnconstrained,
    hDeltaPxPosUnconstrained, hDeltaPyPosUnconstrained, hDeltaPzPosUnconstrained};
  TH1D* momentumResidualsConstrained[6] = {
    hDeltaPxNegConstrained, hDeltaPyNegConstrained, hDeltaPzNegConstrained,
    hDeltaPxPosConstrained, hDeltaPyPosConstrained, hDeltaPzPosConstrained};
  for (int i = 0; i < 6; ++i) {
    cMomRes->cd(i + 1);
    auto* hUnconstrained = momentumResidualsUnconstrained[i];
    auto* hConstrained = momentumResidualsConstrained[i];
    for (auto* hist : {hUnconstrained, hConstrained}) {
      if (hist->GetEntries() > 10 && hist->GetRMS() > 0.) {
        const double fitMin = hist->GetMean() - 2. * hist->GetRMS();
        const double fitMax = hist->GetMean() + 2. * hist->GetRMS();
        hist->Fit("gaus", "Q0", "", fitMin, fitMax);
      }
    }
    auto* hUnconstrainedDraw = static_cast<TH1D*>(hUnconstrained->Clone(Form("%s_draw", hUnconstrained->GetName())));
    auto* hConstrainedDraw = static_cast<TH1D*>(hConstrained->Clone(Form("%s_draw", hConstrained->GetName())));
    hUnconstrainedDraw->SetDirectory(nullptr);
    hConstrainedDraw->SetDirectory(nullptr);
    if (hUnconstrainedDraw->Integral() > 0.) hUnconstrainedDraw->Scale(1. / hUnconstrainedDraw->Integral());
    if (hConstrainedDraw->Integral() > 0.) hConstrainedDraw->Scale(1. / hConstrainedDraw->Integral());
    hUnconstrainedDraw->SetLineColor(kBlue + 1);
    hUnconstrainedDraw->SetLineWidth(2);
    hConstrainedDraw->SetLineColor(kRed + 1);
    hConstrainedDraw->SetLineWidth(2);
    hUnconstrainedDraw->SetMaximum(1.15 * TMath::Max(hUnconstrainedDraw->GetMaximum(), hConstrainedDraw->GetMaximum()));
    hUnconstrainedDraw->Draw("HIST");
    hConstrainedDraw->Draw("HIST SAME");
    auto* legend = new TLegend(0.48, 0.72, 0.88, 0.88);
    legend->SetBorderSize(0);
    const auto* fitUnconstrained = hUnconstrained->GetFunction("gaus");
    const auto* fitConstrained = hConstrained->GetFunction("gaus");
    legend->AddEntry(hUnconstrainedDraw,
                     fitUnconstrained ? Form("Unconstrained, #sigma=%.3g", fitUnconstrained->GetParameter(2)) : "Unconstrained",
                     "l");
    legend->AddEntry(hConstrainedDraw,
                     fitConstrained ? Form("Constrained, #sigma=%.3g", fitConstrained->GetParameter(2)) : "Constrained",
                     "l");
    legend->Draw();
  }
  cMomRes->SaveAs(Form("%s/muon_momentum_resolutionMS.png", dirSimu));

  TCanvas* cTrackPulls = new TCanvas("cTrackPulls", "Single-muon track-parameter pulls", 1500, 900);
  cTrackPulls->Divide(3, 2);
  for (int ip = 0; ip < 5; ++ip) {
    cTrackPulls->cd(ip + 1);
    auto* hUnconstrained = hTrackPullsUnconstrained[ip];
    auto* hConstrained = hTrackPullsConstrained[ip];
    hUnconstrained->SetLineColor(kBlue + 1);
    hUnconstrained->SetLineWidth(2);
    hConstrained->SetLineColor(kRed + 1);
    hConstrained->SetLineWidth(2);
    hUnconstrained->SetTitle(Form("%s pull", parNames[ip]));
    hUnconstrained->SetMaximum(1.15 * TMath::Max(hUnconstrained->GetMaximum(), hConstrained->GetMaximum()));
    hUnconstrained->Draw("HIST");
    hConstrained->Draw("HIST SAME");
    auto* legend = new TLegend(0.48, 0.75, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->AddEntry(hUnconstrained, "Unconstrained", "l");
    legend->AddEntry(hConstrained, "Constrained", "l");
    legend->AddEntry(static_cast<TObject*>(nullptr),
                     Form("Unc. under/over: %.0f / %.0f",
                          hUnconstrained->GetBinContent(0),
                          hUnconstrained->GetBinContent(hUnconstrained->GetNbinsX() + 1)), "");
    legend->AddEntry(static_cast<TObject*>(nullptr),
                     Form("Con. under/over: %.0f / %.0f",
                          hConstrained->GetBinContent(0),
                          hConstrained->GetBinContent(hConstrained->GetNbinsX() + 1)), "");
    legend->Draw();
  }
  cTrackPulls->SaveAs(Form("%s/muon_track_parameter_pullsMS.png", dirSimu));

  TCanvas* cMSTrackFakeHitCategory = new TCanvas("cMSTrackFakeHitCategory", "MS-track fake-hit category", 1400, 700);
  cMSTrackFakeHitCategory->SetBottomMargin(0.20);
  hMSTrackFakeHitCategory->SetLineColor(kBlue + 1);
  hMSTrackFakeHitCategory->SetLineWidth(2);
  hMSTrackFakeHitCategory->SetFillColor(kAzure - 9);
  hMSTrackFakeHitCategory->GetXaxis()->LabelsOption("v");
  const double nCategorizedTracks = hMSTrackFakeHitCategory->Integral(1, hMSTrackFakeHitCategory->GetNbinsX());
  const double maxCategoryCount = hMSTrackFakeHitCategory->GetMaximum();
  if (maxCategoryCount > 0.) {
    hMSTrackFakeHitCategory->SetMaximum(1.15 * maxCategoryCount);
  }
  hMSTrackFakeHitCategory->Draw("HIST");
  if (nCategorizedTracks > 0.) {
    TLatex categoryPercentage;
    categoryPercentage.SetTextAlign(22);
    categoryPercentage.SetTextSize(0.030);
    for (int bin = 1; bin <= hMSTrackFakeHitCategory->GetNbinsX(); ++bin) {
      const double count = hMSTrackFakeHitCategory->GetBinContent(bin);
      categoryPercentage.DrawLatex(hMSTrackFakeHitCategory->GetXaxis()->GetBinCenter(bin),
                                   count + 0.035 * maxCategoryCount,
                                   Form("%.1f%%", 100. * count / nCategorizedTracks));
    }
  }
  cMSTrackFakeHitCategory->SaveAs(Form("%s/ms_track_fake_hit_category.png", dirSimu));

  const TString outputPath = Form("%s/%s", dirSimu, outputFileName);
  TFile fOut(outputPath, "RECREATE");
  if (fOut.IsZombie()) {
    std::fprintf(stderr, "[plotMSDimuons] ERROR: Failed to create histogram file %s\n", outputPath.Data());
    return;
  }
  fOut.cd();

  // Dimuon reconstructed and generated kinematics.
  hdimuy->Write();
  hdimupt->Write();
  hdimumassUnconstrained->Write();
  hdimumassConstrained->Write();
  hdimumassUnconstrainedGood->Write();
  hdimumassConstrainedGood->Write();
  hdimumassUnconstrainedFake->Write();
  hdimumassConstrainedFake->Write();
  hMSTrackFakeHitCategory->Write();
  hdimuygen->Write();
  hdimuptgen->Write();
  hdimumassgen->Write();
  hdimuyrecvsgen->Write();
  hdimuptrecvsgen->Write();
  hdimudeltayrecvsgen->Write();
  hdimudeltaptrecvsgen->Write();
  hdimurecmassvspt->Write();
  hdimurecmassvseta->Write();
  hmassvsetafit->Write();
  hmassvsetamean->Write();
  hdimurecptvseta->Write();
  hdimurecpvseta->Write();

  // Dimuon mass-resolution inputs and fitted summaries.
  hMassVsPUnconstrained->Write();
  hMassVsPConstrained->Write();
  hMassVsYUnconstrained->Write();
  hMassVsYConstrained->Write();
  hMassResolutionVsPUnconstrained->Write();
  hMassResolutionVsPConstrained->Write();
  hMassResolutionVsYUnconstrained->Write();
  hMassResolutionVsYConstrained->Write();
  hRapidity->Write();
  hRapidityNeg->Write();
  hRapidityPos->Write();
  hRapidityMC->Write();
  hRapidityMCNeg->Write();
  hRapidityMCPos->Write();
  hRapidityReconstructable->Write();
  hRapidityReconstructableNeg->Write();
  hRapidityReconstructablePos->Write();
  hRapidityReconstructedMC->Write();
  hRapidityReconstructedMCNeg->Write();
  hRapidityReconstructedMCPos->Write();
  hMomentumMC->Write();
  hMomentumMCNeg->Write();
  hMomentumMCPos->Write();
  hMomentumReconstructable->Write();
  hMomentumReconstructableNeg->Write();
  hMomentumReconstructablePos->Write();
  hMomentumReconstructedMC->Write();
  hMomentumReconstructedMCNeg->Write();
  hMomentumReconstructedMCPos->Write();
  for (auto* hist : momentumResidualsUnconstrained) hist->Write();
  for (auto* hist : momentumResidualsConstrained) hist->Write();
  for (auto* hist : hTrackPullsUnconstrained) hist->Write();
  for (auto* hist : hTrackPullsConstrained) hist->Write();

  TH1D* hRecoEffVsRapidity = static_cast<TH1D*>(hRapidityReconstructedMC->Clone("hRecoEffVsRapidity"));
  TH1D* hRecoEffVsRapidityNeg = static_cast<TH1D*>(hRapidityReconstructedMCNeg->Clone("hRecoEffVsRapidityNeg"));
  TH1D* hRecoEffVsRapidityPos = static_cast<TH1D*>(hRapidityReconstructedMCPos->Clone("hRecoEffVsRapidityPos"));
  hRecoEffVsRapidity->Divide(hRapidityReconstructedMC, hRapidityReconstructable, 1., 1., "B");
  hRecoEffVsRapidityNeg->Divide(hRapidityReconstructedMCNeg, hRapidityReconstructableNeg, 1., 1., "B");
  hRecoEffVsRapidityPos->Divide(hRapidityReconstructedMCPos, hRapidityReconstructablePos, 1., 1., "B");

  TH1D* hEffTimesAccVsRapidity = static_cast<TH1D*>(hRapidityReconstructedMC->Clone("hEffTimesAccVsRapidity"));
  TH1D* hEffTimesAccVsRapidityNeg = static_cast<TH1D*>(hRapidityReconstructedMCNeg->Clone("hEffTimesAccVsRapidityNeg"));
  TH1D* hEffTimesAccVsRapidityPos = static_cast<TH1D*>(hRapidityReconstructedMCPos->Clone("hEffTimesAccVsRapidityPos"));
  hEffTimesAccVsRapidity->Divide(hRapidityReconstructedMC, hRapidityMC, 1., 1., "B");
  hEffTimesAccVsRapidityNeg->Divide(hRapidityReconstructedMCNeg, hRapidityMCNeg, 1., 1., "B");
  hEffTimesAccVsRapidityPos->Divide(hRapidityReconstructedMCPos, hRapidityMCPos, 1., 1., "B");

  TH1D* recoEfficiencies[3] = {hRecoEffVsRapidity, hRecoEffVsRapidityNeg, hRecoEffVsRapidityPos};
  TH1D* efficienciesTimesAcceptance[3] = {hEffTimesAccVsRapidity, hEffTimesAccVsRapidityNeg, hEffTimesAccVsRapidityPos};
  const char* efficiencyTitles[3] = {"Single muons", "#mu^{-}", "#mu^{+}"};
  TCanvas* cEfficiency = new TCanvas("cEfficiency", "Single-muon reconstruction efficiency", 1800, 600);
  cEfficiency->Divide(3, 1);
  fOut.cd();
  for (int i = 0; i < 3; ++i) {
    cEfficiency->cd(i + 1);
    auto* recoEfficiency = recoEfficiencies[i];
    auto* efficiencyTimesAcceptance = efficienciesTimesAcceptance[i];
    recoEfficiency->SetStats(kFALSE);
    efficiencyTimesAcceptance->SetStats(kFALSE);
    recoEfficiency->SetTitle(Form("%s;MC rapidity;Efficiency", efficiencyTitles[i]));
    recoEfficiency->SetMinimum(0.);
    recoEfficiency->SetMaximum(1.05);
    recoEfficiency->SetMarkerStyle(20);
    recoEfficiency->SetMarkerColor(kBlue + 1);
    recoEfficiency->SetLineColor(kBlue + 1);
    efficiencyTimesAcceptance->SetMarkerStyle(24);
    efficiencyTimesAcceptance->SetMarkerColor(kRed + 1);
    efficiencyTimesAcceptance->SetLineColor(kRed + 1);
    recoEfficiency->Draw("EP");
    efficiencyTimesAcceptance->Draw("EP SAME");
    auto* legend = new TLegend(0.52, 0.16, 0.88, 0.30);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.035);
    legend->AddEntry(recoEfficiency, "N_{reco}/N_{reconstructable}", "lp");
    legend->AddEntry(efficiencyTimesAcceptance, "N_{reco}/N_{generated}", "lp");
    legend->Draw();
    recoEfficiency->Write();
    efficiencyTimesAcceptance->Write();
  }
  cEfficiency->SaveAs(Form("%s/single_muon_efficiency_vs_rapidityMS.png", dirSimu));

  TH1D* hRecoEffVsRapidityCM = static_cast<TH1D*>(hRecoEffVsRapidity->Clone("hRecoEffVsRapidityCM"));
  TH1D* hRecoEffVsRapidityCMNeg = static_cast<TH1D*>(hRecoEffVsRapidityNeg->Clone("hRecoEffVsRapidityCMNeg"));
  TH1D* hRecoEffVsRapidityCMPos = static_cast<TH1D*>(hRecoEffVsRapidityPos->Clone("hRecoEffVsRapidityCMPos"));
  TH1D* hEffTimesAccVsRapidityCM = static_cast<TH1D*>(hEffTimesAccVsRapidity->Clone("hEffTimesAccVsRapidityCM"));
  TH1D* hEffTimesAccVsRapidityCMNeg = static_cast<TH1D*>(hEffTimesAccVsRapidityNeg->Clone("hEffTimesAccVsRapidityCMNeg"));
  TH1D* hEffTimesAccVsRapidityCMPos = static_cast<TH1D*>(hEffTimesAccVsRapidityPos->Clone("hEffTimesAccVsRapidityCMPos"));
  TH1D* recoEfficienciesCM[3] = {hRecoEffVsRapidityCM, hRecoEffVsRapidityCMNeg, hRecoEffVsRapidityCMPos};
  TH1D* efficienciesTimesAcceptanceCM[3] = {hEffTimesAccVsRapidityCM, hEffTimesAccVsRapidityCMNeg, hEffTimesAccVsRapidityCMPos};
  TCanvas* cEfficiencyCM = new TCanvas("cEfficiencyCM", "Single-muon efficiency versus center-of-mass rapidity", 1800, 600);
  cEfficiencyCM->Divide(3, 1);
  fOut.cd();
  for (int i = 0; i < 3; ++i) {
    cEfficiencyCM->cd(i + 1);
    auto* recoEfficiency = recoEfficienciesCM[i];
    auto* efficiencyTimesAcceptance = efficienciesTimesAcceptanceCM[i];
    recoEfficiency->SetStats(kFALSE);
    efficiencyTimesAcceptance->SetStats(kFALSE);
    recoEfficiency->GetXaxis()->SetLimits(1. - yCMShift, 7. - yCMShift);
    efficiencyTimesAcceptance->GetXaxis()->SetLimits(1. - yCMShift, 7. - yCMShift);
    recoEfficiency->SetTitle(Form("%s;MC y-y_{CM};Efficiency", efficiencyTitles[i]));
    recoEfficiency->SetMinimum(0.);
    recoEfficiency->SetMaximum(1.05);
    recoEfficiency->SetMarkerStyle(20);
    recoEfficiency->SetMarkerColor(kBlue + 1);
    recoEfficiency->SetLineColor(kBlue + 1);
    efficiencyTimesAcceptance->SetMarkerStyle(24);
    efficiencyTimesAcceptance->SetMarkerColor(kRed + 1);
    efficiencyTimesAcceptance->SetLineColor(kRed + 1);
    recoEfficiency->Draw("EP");
    efficiencyTimesAcceptance->Draw("EP SAME");
    auto* legend = new TLegend(0.52, 0.16, 0.88, 0.30);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.035);
    legend->AddEntry(recoEfficiency, "N_{reco}/N_{reconstructable}", "lp");
    legend->AddEntry(efficiencyTimesAcceptance, "N_{reco}/N_{generated}", "lp");
    legend->Draw();
    recoEfficiency->Write();
    efficiencyTimesAcceptance->Write();
  }
  printf("[plotMSDimuons] Fixed-target rapidity shift at %.1f GeV/nucleon: yCM=%.6f\n",
         NA6PBeamParam::Instance().energyPerNucleon, yCMShift);
  cEfficiencyCM->SaveAs(Form("%s/single_muon_efficiency_vs_rapidityCM_MS.png", dirSimu));

  TH1D* hRecoEffVsMomentum = static_cast<TH1D*>(hMomentumReconstructedMC->Clone("hRecoEffVsMomentum"));
  TH1D* hRecoEffVsMomentumNeg = static_cast<TH1D*>(hMomentumReconstructedMCNeg->Clone("hRecoEffVsMomentumNeg"));
  TH1D* hRecoEffVsMomentumPos = static_cast<TH1D*>(hMomentumReconstructedMCPos->Clone("hRecoEffVsMomentumPos"));
  hRecoEffVsMomentum->Divide(hMomentumReconstructedMC, hMomentumReconstructable, 1., 1., "B");
  hRecoEffVsMomentumNeg->Divide(hMomentumReconstructedMCNeg, hMomentumReconstructableNeg, 1., 1., "B");
  hRecoEffVsMomentumPos->Divide(hMomentumReconstructedMCPos, hMomentumReconstructablePos, 1., 1., "B");

  TH1D* hEffTimesAccVsMomentum = static_cast<TH1D*>(hMomentumReconstructedMC->Clone("hEffTimesAccVsMomentum"));
  TH1D* hEffTimesAccVsMomentumNeg = static_cast<TH1D*>(hMomentumReconstructedMCNeg->Clone("hEffTimesAccVsMomentumNeg"));
  TH1D* hEffTimesAccVsMomentumPos = static_cast<TH1D*>(hMomentumReconstructedMCPos->Clone("hEffTimesAccVsMomentumPos"));
  hEffTimesAccVsMomentum->Divide(hMomentumReconstructedMC, hMomentumMC, 1., 1., "B");
  hEffTimesAccVsMomentumNeg->Divide(hMomentumReconstructedMCNeg, hMomentumMCNeg, 1., 1., "B");
  hEffTimesAccVsMomentumPos->Divide(hMomentumReconstructedMCPos, hMomentumMCPos, 1., 1., "B");

  TH1D* momentumRecoEfficiencies[3] = {hRecoEffVsMomentum, hRecoEffVsMomentumNeg, hRecoEffVsMomentumPos};
  TH1D* momentumEfficienciesTimesAcceptance[3] = {hEffTimesAccVsMomentum, hEffTimesAccVsMomentumNeg, hEffTimesAccVsMomentumPos};
  TCanvas* cMomentumEfficiency = new TCanvas("cMomentumEfficiency", "Single-muon efficiency versus momentum", 1800, 600);
  cMomentumEfficiency->Divide(3, 1);
  fOut.cd();
  for (int i = 0; i < 3; ++i) {
    cMomentumEfficiency->cd(i + 1);
    auto* recoEfficiency = momentumRecoEfficiencies[i];
    auto* efficiencyTimesAcceptance = momentumEfficienciesTimesAcceptance[i];
    recoEfficiency->SetStats(kFALSE);
    efficiencyTimesAcceptance->SetStats(kFALSE);
    recoEfficiency->SetTitle(Form("%s;MC p (GeV/c);Efficiency", efficiencyTitles[i]));
    recoEfficiency->SetMinimum(0.);
    recoEfficiency->SetMaximum(1.05);
    recoEfficiency->SetMarkerStyle(20);
    recoEfficiency->SetMarkerColor(kBlue + 1);
    recoEfficiency->SetLineColor(kBlue + 1);
    efficiencyTimesAcceptance->SetMarkerStyle(24);
    efficiencyTimesAcceptance->SetMarkerColor(kRed + 1);
    efficiencyTimesAcceptance->SetLineColor(kRed + 1);
    recoEfficiency->Draw("EP");
    efficiencyTimesAcceptance->Draw("EP SAME");
    auto* legend = new TLegend(0.52, 0.16, 0.88, 0.30);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.035);
    legend->AddEntry(recoEfficiency, "N_{reco}/N_{reconstructable}", "lp");
    legend->AddEntry(efficiencyTimesAcceptance, "N_{reco}/N_{generated}", "lp");
    legend->Draw();
    recoEfficiency->Write();
    efficiencyTimesAcceptance->Write();
  }
  cMomentumEfficiency->SaveAs(Form("%s/single_muon_efficiency_vs_momentumMS.png", dirSimu));
  fOut.Close();
  printf("[plotMSDimuons] Histograms saved to %s\n", outputPath.Data());
}
