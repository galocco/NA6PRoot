// NA6PCCopyright

#include "ConfigurableParam.h"
#include <fairlogger/Logger.h>
#include <TApplication.h>
#include <TMCManager.h>
#include <TGeoGlobalMagField.h>
#include <boost/program_options.hpp>
#include <filesystem>
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TStopwatch.h>
#include "NA6PVerTelHit.h"
#include "NA6PMuonSpecModularHit.h"
#include "MagneticField.h"
#include "StringUtils.h"
#include "NA6PVerTelReconstruction.h"
#include "NA6PMuonSpecReconstruction.h"
#include "NA6PMatching.h"
#include "NA6PTOFMatching.h"
#include "NA6PTOFHit.h"

double getPrimaryVertexZ(TTree* mcTree, std::vector<TParticle>* mcArr, int eventID)
{
  if (!mcTree || !mcArr) {
    return 0.0;
  }
  mcTree->GetEvent(eventID);
  for (const auto& part : *mcArr) {
    if (part.IsPrimary()) {
      return part.Vz();
    }
  }
  return 0.0;
}

int main(int argc, char** argv)
{
  const std::string paramsIni = "na6pRecoParam.ini";
  fair::Logger::OnFatal([]() { throw std::runtime_error("Fatal error"); });

  namespace bpo = boost::program_options;
  bpo::variables_map vm;
  bpo::options_description opt_general("Usage:\n  " + std::string(argv[0]));
  bpo::options_description opt_hidden("");
  bpo::options_description opt_all;
  bpo::positional_options_description opt_pos;

  try {
    auto add_option = opt_general.add_options();
    add_option("help,h", "Print this help message");
    add_option("verbosity,v", bpo::value<int>()->default_value(0), "verbosity level [0 = no output]");
    add_option("configKeyValues", bpo::value<std::string>()->default_value(""), "comma-separated configKeyValues");
    add_option("load-ini", bpo::value<std::string>()->default_value(""), "load configurables from ini file (if defined), overridden by configKeyValues");
    add_option("load-recoparam", bpo::value<std::string>()->default_value(""), "load reco parameters from ini file (if defined), overridden by configKeyValues");
    add_option("disable-write-ini", bpo::value<bool>()->default_value(false)->implicit_value(true), "do not write reco parameters ini file");
    add_option("firstevent,f", bpo::value<int32_t>()->default_value(0), "first event");
    add_option("lastevent,l", bpo::value<int32_t>()->default_value(-1), "last event");
    add_option("doHitsToRecPoints,cl", bpo::value<bool>()->default_value(true), "run hits->clusters");
    add_option("doTrackletVertex,vert", bpo::value<bool>()->default_value(true), "run tracklet vertexer");
    add_option("doVTTracking,vt", bpo::value<bool>()->default_value(true), "run VT tracker");
    add_option("doMSTracking,ms", bpo::value<bool>()->default_value(true), "run MS tracker");
    add_option("doMatching,mt", bpo::value<bool>()->default_value(true), "run matching between VT and MS");
    add_option("doTOFHitsToRecPoints", bpo::value<bool>()->default_value(true), "run TOF hits->clusters");
    add_option("doTOFMatching", bpo::value<bool>()->default_value(true), "run TOF matching to VT tracks");
    opt_all.add(opt_general).add(opt_hidden);
    bpo::store(bpo::command_line_parser(argc, argv).options(opt_all).positional(opt_pos).run(), vm);

    if (vm.count("help")) {
      std::cout << opt_general << std::endl;
      exit(0);
    }

    bpo::notify(vm);
  } catch (bpo::error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl
              << std::endl;
    std::cerr << opt_general << std::endl;
    exit(1);
  } catch (std::exception& e) {
    std::cerr << e.what() << ", application will now exit" << std::endl;
    exit(2);
  }

  auto flini = vm["load-ini"].as<std::string>();
  na6p::conf::ConfigurableParam::updateFromString(vm["configKeyValues"].as<std::string>());
  if (!flini.empty()) {
    na6p::conf::ConfigurableParam::updateFromFile(flini, "", true);
  }
  auto flrp = vm["load-recoparam"].as<std::string>();
  if (!flrp.empty()) {
    na6p::conf::ConfigurableParam::updateFromFile(flrp, "", true);
  }
  LOGP(info, "Printing all configs");
  na6p::conf::ConfigurableParam::printAllKeyValuePairs();

  const bool doHitsToRecPoints = vm["doHitsToRecPoints"].as<bool>();
  const bool doTrackletVertex = vm["doTrackletVertex"].as<bool>();
  const bool doVTTracking = vm["doVTTracking"].as<bool>();
  const bool doMSTracking = vm["doMSTracking"].as<bool>();
  const bool doMatching = vm["doMatching"].as<bool>();
  const bool doTOFHitsToRecPoints = vm["doTOFHitsToRecPoints"].as<bool>();
  const bool doTOFMatching = vm["doTOFMatching"].as<bool>();

  int firstEv = vm["firstevent"].as<int32_t>();
  int lastEv = vm["lastevent"].as<int32_t>();

  auto magField = new MagneticField();
  magField->loadField();
  magField->setAsGlobalField();

  NA6PVerTelReconstruction* vtrec = new NA6PVerTelReconstruction();
  vtrec->setGeometryFile("geometry.root");
  NA6PMuonSpecReconstruction* msrec = new NA6PMuonSpecReconstruction();
  NA6PMatching* matching = new NA6PMatching();
  NA6PTOFMatching* tofMatching = new NA6PTOFMatching();

  if (doHitsToRecPoints) {
    TFile* fhVT = TFile::Open("HitsVerTel.root");
    if (!fhVT || fhVT->IsZombie()) {
      LOGP(error, "Cannot open HitsVerTel.root");
      if (fhVT)
        delete fhVT;
      return -1;
    }
    TTree* thVT = (TTree*)fhVT->Get("hitsVerTel");
    if (!thVT) {
      LOGP(error, "Cannot find tree 'hitsVerTel' in HitsVerTel.root");
      fhVT->Close();
      delete fhVT;
      return -1;
    }
    std::vector<NA6PVerTelHit> vtHits, *vtHitsPtr = &vtHits;
    thVT->SetBranchAddress("VerTel", &vtHitsPtr);
    int nEvVT = thVT->GetEntriesFast();

    vtrec->createClustersOutput();
    for (int jEv = 0; jEv < nEvVT; jEv++) {
      thVT->GetEvent(jEv);
      int nHits = vtHits.size();
      LOGP(info, "VerTel Event {} nHits= {}", jEv, nHits);
      vtrec->clearClusters();
      vtrec->hitsToRecPoints(vtHits);
      vtrec->writeClusters();
    }
    vtrec->closeClustersOutput();
    fhVT->Close();
    delete fhVT;

    // Muon Spectrometer hits -> clusters
    TFile* fhMS = TFile::Open("HitsMuonSpecModular.root");
    if (!fhMS || fhMS->IsZombie()) {
      LOGP(error, "Cannot open HitsMuonSpecModular.root");
      if (fhMS)
        delete fhMS;
      return -1;
    }
    TTree* thMS = (TTree*)fhMS->Get("hitsMuonSpecModular");
    if (!thMS) {
      LOGP(error, "Cannot find tree 'hitsMuonSpecModular' in HitsMuonSpecModular.root");
      fhMS->Close();
      delete fhMS;
      return -1;
    }
    std::vector<NA6PMuonSpecModularHit> msHits, *msHitsPtr = &msHits;
    thMS->SetBranchAddress("MuonSpecModular", &msHitsPtr);
    int nEvMS = thMS->GetEntriesFast();

    msrec->createClustersOutput();
    for (int jEv = 0; jEv < nEvMS; jEv++) {
      thMS->GetEvent(jEv);
      int nHits = msHits.size();
      LOGP(info, "MuonSpec Event {} nHits= {}", jEv, nHits);
      msrec->clearClusters();
      msrec->hitsToRecPoints(msHits);
      msrec->writeClusters();
    }
    msrec->closeClustersOutput();
    fhMS->Close();
    delete fhMS;

    // ── TOF hits → clusters ──────────────────────────────────────
    if (doTOFHitsToRecPoints) {
      TFile* fhTOF = TFile::Open("HitsTOF.root");
      if (!fhTOF || fhTOF->IsZombie()) {
        LOGP(error, "Cannot open HitsTOF.root");
        if (fhTOF)
          delete fhTOF;
        return -1;
      }
      TTree* thTOF = (TTree*)fhTOF->Get("hitsTOF");
      if (!thTOF) {
        LOGP(error, "Cannot find tree 'hitsTOF' in HitsTOF.root");
        fhTOF->Close();
        delete fhTOF;
        return -1;
      }
      std::vector<NA6PTOFHit> tofHits, *tofHitsPtr = &tofHits;
      thTOF->SetBranchAddress("TOF", &tofHitsPtr);
      int nEvTOF = thTOF->GetEntriesFast();

      tofMatching->createClustersOutput();
      for (int jEv = 0; jEv < nEvTOF; jEv++) {
        thTOF->GetEvent(jEv);
        LOGP(info, "TOF Event {} nHits= {}", jEv, (int)tofHits.size());
        tofMatching->clearClusters();
        tofMatching->hitsToRecPoints(tofHits);
        tofMatching->writeClusters();
      }
      tofMatching->closeClustersOutput();
      fhTOF->Close();
      delete fhTOF;
    }
  } else {
    LOGP(info, "Hits -> Recpoints disabled from input options");
  }

  const bool needsMCKine = doTrackletVertex || doVTTracking || doMSTracking || doMatching;
  TFile* fk = nullptr;
  TTree* mcTree = nullptr;
  std::vector<TParticle>* mcArr = nullptr;
  if (needsMCKine) {
    fk = TFile::Open("MCKine.root");
    if (!fk || fk->IsZombie()) {
      LOGP(error, "Cannot open MCKine.root");
      if (fk)
        delete fk;
      return -1;
    }
    mcTree = (TTree*)fk->Get("mckine");
    if (!mcTree) {
      LOGP(error, "Cannot find tree 'mckine' in MCKine.root");
      fk->Close();
      delete fk;
      return -1;
    }
    mcTree->SetBranchAddress("tracks", &mcArr);
  }

  TFile* fcVT = nullptr;
  TTree* tcVT = nullptr;
  std::vector<NA6PVerTelCluster> vtClus, *vtClusPtr = &vtClus;

  if (doTrackletVertex || doVTTracking) {
    fcVT = TFile::Open("ClustersVerTel.root");
    if (!fcVT || fcVT->IsZombie()) {
      LOGP(error, "Cannot open ClustersVerTel.root");
      if (fcVT)
        delete fcVT;
      return -1;
    }
    tcVT = (TTree*)fcVT->Get("clustersVerTel");
    if (!tcVT) {
      LOGP(error, "Cannot find tree 'clustersVerTel' in ClustersVerTel.root");
      fcVT->Close();
      delete fcVT;
      return -1;
    }
    tcVT->SetBranchAddress("VerTel", &vtClusPtr);
  }

  TFile* fcMS = nullptr;
  TTree* tcMS = nullptr;
  std::vector<NA6PMuonSpecCluster> msClus, *msClusPtr = &msClus;

  if (doMSTracking) {
    fcMS = TFile::Open("ClustersMuonSpec.root");
    if (!fcMS || fcMS->IsZombie()) {
      LOGP(error, "Cannot open ClustersMuonSpec.root");
      if (fcMS)
        delete fcMS;
      return -1;
    }
    tcMS = (TTree*)fcMS->Get("clustersMuonSpec");
    if (!tcMS) {
      LOGP(error, "Cannot find tree 'clustersMuonSpec' in ClustersMuonSpec.root");
      fcMS->Close();
      delete fcMS;
      return -1;
    }
    tcMS->SetBranchAddress("MuonSpec", &msClusPtr);
  }

  if (doTrackletVertex || doVTTracking || doMSTracking) {
    int nEv = -1;
    if (tcVT) {
      nEv = tcVT->GetEntries();
    }
    if (tcMS) {
      int nEvMS = tcMS->GetEntries();
      nEv = (nEv < 0) ? nEvMS : std::min(nEv, nEvMS);
    }
    if (nEv < 0) {
      LOGP(error, "No input trees for tracking stage");
      return -1;
    }
    if (lastEv > nEv || lastEv < 0)
      lastEv = nEv;
    if (firstEv < 0)
      firstEv = 0;

    if (doTrackletVertex)
      vtrec->initVertexer();
    if (doVTTracking)
      vtrec->initTracker();
    if (doMSTracking)
      msrec->initTracker();
    
    TStopwatch timer;
    timer.Start();

    NA6PVertex pvert;
    for (int jEv = firstEv; jEv < lastEv; jEv++) {
      LOGP(info, "Process event {}", jEv);
      const double zvert = getPrimaryVertexZ(mcTree, mcArr, jEv);
      pvert.setXYZ(0.f, 0.f, zvert);

      if (doTrackletVertex || doVTTracking) {
        LOGP(info, "VT reconstruction");
        vtrec->clearEvent();
        tcVT->GetEvent(jEv);
        vtrec->setPrimaryVertex(&pvert);
        vtrec->setClusters(vtClus);
        if (doTrackletVertex)
          vtrec->runVertexerTracklets();
        if (doVTTracking)
          vtrec->runTracking();
      }

      if (doMSTracking) {
        LOGP(info, "MS reconstruction");
        msrec->clearEvent();
        tcMS->GetEvent(jEv);
        msrec->setPrimaryVertex(&pvert);
        msrec->setClusters(msClus);
        msrec->runTracking();
      }
    }

    if (doTrackletVertex || doVTTracking) {
      vtrec->closeVerticesOutput();
      vtrec->closeTracksOutput();
    }
    if (doMSTracking) {
      msrec->closeTracksOutput();
    }

    timer.Stop();
    LOGP(info, "------ Tracking time ------");
    timer.Print();
    LOGP(info, "------------------------------");
  }

  if (doMatching) {
    matching->initMatching();

    TFile* ft = TFile::Open("TracksVerTel.root");
    if (!ft || ft->IsZombie()) {
      LOGP(error, "Cannot open TracksVerTel.root");
      if (ft)
        delete ft;
      return -1;
    }
    TTree* tt = (TTree*)ft->Get("tracksVerTel");
    if (!tt) {
      LOGP(error, "Cannot find tree 'tracksVerTel' in TracksVerTel.root");
      ft->Close();
      delete ft;
      return -1;
    }
    std::vector<NA6PTrack> vtTracks, *vtTracksPtr = &vtTracks;
    tt->SetBranchAddress("VerTel", &vtTracksPtr);

    TFile* fm = TFile::Open("TracksMuonSpec.root");
    if (!fm || fm->IsZombie()) {
      LOGP(error, "Cannot open TracksMuonSpec.root");
      if (fm)
        delete fm;
      return -1;
    }
    TTree* tm = (TTree*)fm->Get("tracksMuonSpec");
    if (!tm) {
      LOGP(error, "Cannot find tree 'tracksMuonSpec' in TracksMuonSpec.root");
      fm->Close();
      delete fm;
      return -1;
    }
    std::vector<NA6PTrack> msTracks, *msTracksPtr = &msTracks;
    tm->SetBranchAddress("MuonSpec", &msTracksPtr);

    TFile* fvc = TFile::Open("ClustersVerTel.root");
    if (!fvc || fvc->IsZombie()) {
      LOGP(error, "Cannot open ClustersVerTel.root");
      if (fvc)
        delete fvc;
      return -1;
    }
    TTree* tvertelc = (TTree*)fvc->Get("clustersVerTel");
    if (!tvertelc) {
      LOGP(error, "Cannot find tree 'clustersVerTel' in ClustersVerTel.root");
      fvc->Close();
      delete fvc;
      return -1;
    }
    std::vector<NA6PVerTelCluster> vtClusMatch, *vtClusMatchPtr = &vtClusMatch;
    tvertelc->SetBranchAddress("VerTel", &vtClusMatchPtr);

    TFile* fmc = TFile::Open("ClustersMuonSpec.root");
    if (!fmc || fmc->IsZombie()) {
      LOGP(error, "Cannot open ClustersMuonSpec.root");
      if (fmc)
        delete fmc;
      return -1;
    }
    TTree* tmuonspec = (TTree*)fmc->Get("clustersMuonSpec");
    if (!tmuonspec) {
      LOGP(error, "Cannot find tree 'clustersMuonSpec' in ClustersMuonSpec.root");
      fmc->Close();
      delete fmc;
      return -1;
    }
    std::vector<NA6PMuonSpecCluster> msClusMatch, *msClusMatchPtr = &msClusMatch;
    tmuonspec->SetBranchAddress("MuonSpec", &msClusMatchPtr);

    int nEvVT = tt->GetEntries();
    int nEvMS = tm->GetEntries();
    int nEv = std::min(nEvVT, nEvMS);
    if (lastEv > nEv || lastEv < 0)
      lastEv = nEv;
    if (firstEv < 0)
      firstEv = 0;

    TStopwatch timer;
    timer.Start();

    NA6PVertex pvert;
    for (int jEv = firstEv; jEv < lastEv; jEv++) {
      LOGP(info, "Process event {}", jEv);
      const double zvert = getPrimaryVertexZ(mcTree, mcArr, jEv);
      pvert.setXYZ(0.f, 0.f, zvert);

      tvertelc->GetEvent(jEv);
      tmuonspec->GetEvent(jEv);
      tt->GetEvent(jEv);
      tm->GetEvent(jEv);

      matching->setPrimaryVertex(&pvert);
      matching->setVerTelClusters(vtClusMatch);
      matching->setMuonSpecClusters(msClusMatch);
      matching->setVerTelTracks(vtTracks);
      matching->setMuonSpecTracks(msTracks);

      matching->runMatching();
    }
    matching->closeTracksOutput();

    timer.Stop();
    LOGP(info, "------ Matching time ------");
    timer.Print();
    LOGP(info, "------------------------------");

    ft->Close();
    delete ft;
    fm->Close();
    delete fm;
    fvc->Close();
    delete fvc;
    fmc->Close();
    delete fmc;
  }

  // ── TOF matching ────────────────────────────────────────────────
  if (doTOFMatching) {
    tofMatching->initTOFMatching();

    TFile* ftVT = TFile::Open("TracksVerTel.root");
    if (!ftVT || ftVT->IsZombie()) {
      LOGP(error, "Cannot open TracksVerTel.root for TOF matching");
      if (ftVT)
        delete ftVT;
      return -1;
    }
    TTree* ttVT = (TTree*)ftVT->Get("tracksVerTel");
    if (!ttVT) {
      LOGP(error, "Cannot find tree 'tracksVerTel' in TracksVerTel.root");
      ftVT->Close();
      delete ftVT;
      return -1;
    }
    std::vector<NA6PTrack> vtTracksTOF, *vtTracksTOFPtr = &vtTracksTOF;
    ttVT->SetBranchAddress("VerTel", &vtTracksTOFPtr);

    TFile* fcTOF = TFile::Open("ClustersTOF.root");
    if (!fcTOF || fcTOF->IsZombie()) {
      LOGP(error, "Cannot open ClustersTOF.root");
      if (fcTOF)
        delete fcTOF;
      return -1;
    }
    TTree* tcTOF = (TTree*)fcTOF->Get("clustersTOF");
    if (!tcTOF) {
      LOGP(error, "Cannot find tree 'clustersTOF' in ClustersTOF.root");
      fcTOF->Close();
      delete fcTOF;
      return -1;
    }
    std::vector<NA6PTOFCluster> tofClus, *tofClusPtr = &tofClus;
    tcTOF->SetBranchAddress("TOF", &tofClusPtr);

    int nEvVTtof = ttVT->GetEntries();
    int nEvTOFc  = tcTOF->GetEntries();
    int nEvTOF   = std::min(nEvVTtof, nEvTOFc);
    int firstEvTOF = (firstEv < 0) ? 0 : firstEv;
    int lastEvTOF  = (lastEv > nEvTOF || lastEv < 0) ? nEvTOF : lastEv;

    TStopwatch timerTOF;
    timerTOF.Start();

    for (int jEv = firstEvTOF; jEv < lastEvTOF; jEv++) {
      LOGP(info, "TOF matching event {}", jEv);
      ttVT->GetEvent(jEv);
      tcTOF->GetEvent(jEv);

      tofMatching->clearTracks();
      tofMatching->setVerTelTracks(vtTracksTOF);
      tofMatching->setClusters(tofClus);
      tofMatching->runTOFMatching();
    }
    tofMatching->closeTracksOutput();

    timerTOF.Stop();
    LOGP(info, "------ TOF Matching time ------");
    timerTOF.Print();
    LOGP(info, "------------------------------");

    ftVT->Close();
    delete ftVT;
    fcTOF->Close();
    delete fcTOF;
  }

  if (fcVT) {
    fcVT->Close();
    delete fcVT;
  }
  if (fcMS) {
    fcMS->Close();
    delete fcMS;
  }
  if (fk) {
    fk->Close();
    delete fk;
  }

  if (!vm["disable-write-ini"].as<bool>()) {
    try {
      std::string pth = na6p::utils::Str::rectifyDirectory(na6p::conf::ConfigurableParam::getOutputDir());
      if (pth.empty()) {
        pth = std::filesystem::current_path();
      }
      na6p::conf::ConfigurableParam::writeINI(paramsIni, "reco");
      LOGP(info, "Stored configurable params to {}/{}", pth, paramsIni);
    } catch (std::exception e) {
      LOGP(error, "Failed to store configurable params, the reason is {}", e.what());
    }
  }

  delete vtrec;
  delete msrec;
  delete matching;
  delete tofMatching;
  return 0;
}
