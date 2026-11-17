#include <TFile.h>
#include <TTree.h>
#include <iostream>

// A successful exit alone does not prove that the event reached every stage.
int main()
{
  const char* files[] = {"ClustersVerTel.root", "ClustersMuonSpec.root",
                         "VerticesVerTel.root", "TracksVerTel.root",
                         "TracksMuonSpec.root", "TracksMatching.root"};
  const char* trees[] = {"clustersVerTel", "clustersMuonSpec", "verticesVerTel",
                         "tracksVerTel", "tracksMuonSpec", "tracksMatching"};
  for (int i = 0; i < 6; ++i) {
    TFile file(files[i], "READ");
    TTree* tree = nullptr;
    if (!file.IsZombie()) {
      file.GetObject(trees[i], tree);
    }
    if (!tree || tree->GetEntries() != 1) {
      std::cerr << "Expected one event in " << files[i] << ":" << trees[i] << std::endl;
      return 1;
    }
  }
  return 0;
}
