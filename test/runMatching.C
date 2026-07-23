#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TStopwatch.h>
#include "NA6PVerTelCluster.h"
#include "NA6PMuonSpecCluster.h"
#include "NA6PMatching.h"
#include "ConfigurableParam.h"
#endif

void runMatching(bool useMC = false,
                 double zmatch = 38.1175,
                 int firstEv = 0,
                 int lastEv = 99999,
                 const char* dirSimu = ".",
                 const char* dirRec = ".")
{
  NA6PRecoParam::setValue("reco.mtZMatching", std::to_string(zmatch));
  na6p::conf::ConfigurableParam::updateFromFile(Form("%s/na6pLayout.ini", dirSimu), "", true);
  if (!Propagator::loadField() || !Propagator::loadGeometry(Form("%s/geometry.root", dirSimu))) {
    return;
  }

  std::unique_ptr<NA6PMatching> matching = std::make_unique<NA6PMatching>();
  matching->setMCMatching(useMC);

  // kine file needed to get the primary vertex position (temporary)
  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  // io cluster file for VerTel
  TFile* fvertelc = new TFile(Form("%s/ClustersVerTel.root", dirSimu), "READ");
  printf("Open cluster file: %s\n", fvertelc->GetName());
  TTree* tvertelc = (TTree*)fvertelc->Get("clustersVerTel");
  std::vector<NA6PVerTelCluster> vtClus, *vtClusPtr = &vtClus;
  tvertelc->SetBranchAddress("VerTel", &vtClusPtr);

  // io cluster file for MuonSpec
  TFile* fmuonspec = new TFile(Form("%s/ClustersMuonSpec.root", dirSimu), "READ");
  printf("Open cluster file: %s\n", fmuonspec->GetName());
  TTree* tmuonspec = (TTree*)fmuonspec->Get("clustersMuonSpec");
  std::vector<NA6PMuonSpecCluster> msClus, *msClusPtr = &msClus;
  tmuonspec->SetBranchAddress("MuonSpec", &msClusPtr);

  // io track file for VerTel
  TFile* fvttracks = new TFile(Form("%s/TracksVerTel.root", dirRec), "READ");
  printf("Open track file: %s\n", fvttracks->GetName());
  TTree* tvttracks = (TTree*)fvttracks->Get("tracksVerTel");
  std::vector<NA6PTrack> vtTracks, *vtTrackPtr = &vtTracks;
  std::vector<NA6PMCComposedLabel> vtTrMCLabs, *vtTrMCLabsPtr = &vtTrMCLabs;
  tvttracks->SetBranchAddress("VerTel", &vtTrackPtr);
  tvttracks->SetBranchAddress("VerTelMCTruth", &vtTrMCLabsPtr);
  // io track file for MuonSpec
  TFile* fmstracks = new TFile(Form("%s/TracksMuonSpec.root", dirRec), "READ");
  printf("Open track file: %s\n", fmstracks->GetName());
  TTree* tmstracks = (TTree*)fmstracks->Get("tracksMuonSpec");
  std::vector<NA6PTrack> msTracks, *msTrackPtr = &msTracks;
  std::vector<NA6PMCComposedLabel> msTrMCLabs, *msTrMCLabsPtr = &msTrMCLabs;
  tmstracks->SetBranchAddress("MuonSpec", &msTrackPtr);
  tmstracks->SetBranchAddress("MuonSpecMCTruth", &msTrMCLabsPtr);

  int nEv = tmstracks->GetEntries();
  if (lastEv > nEv || lastEv < 0)
    lastEv = nEv;
  if (firstEv < 0)
    firstEv = 0;

  TStopwatch timer;
  timer.Start();

  matching->createTracksOutput();
  for (int jEv = firstEv; jEv < lastEv; jEv++) {
    mcTree->GetEvent(jEv);
    int nPart = mcArr->size();
    double zvert = 0;
    // get primary vertex position from the Kine Tree
    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      if (curPart.IsPrimary()) {
        zvert = curPart.Vz();
        break;
      }
    }
    tvertelc->GetEvent(jEv);
    tmuonspec->GetEvent(jEv);
    tvttracks->GetEvent(jEv);
    tmstracks->GetEvent(jEv);

    NA6PVertex pvert;
    pvert.setXYZ(0., 0., zvert);
    matching->setPrimaryVertex(&pvert);

    matching->setVerTelClusters(vtClus);
    matching->setMuonSpecClusters(msClus);
    matching->setVerTelTracks(vtTracks);
    matching->setVerTelTrackMCLabels(vtTrMCLabs);
    matching->setMuonSpecTracks(msTracks);
    matching->setMuonSpecTrackMCLabels(msTrMCLabs);

    matching->runMatching();
  }
  matching->closeTracksOutput();
  timer.Stop();
  timer.Print();
}
