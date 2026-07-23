// NA6PCCopyright

#include <fairlogger/Logger.h>
#include <TFile.h>
#include <TSystem.h>
#include "NA6PTrack.h"
#include "NA6PMCTruthContainer.h"
#include "NA6PReconstruction.h"

void NA6PReconstruction::createClustersOutput()
{
  LOGP(warning, "createClustersOutput called for {}, this should not happen", getName());
}

void NA6PReconstruction::clearClusters()
{
  LOGP(warning, "clearClusters called for {}, this should not happen", getName());
}

void NA6PReconstruction::writeClusters()
{
  LOGP(warning, "writeClusters called for {}, this should not happen", getName());
}

void NA6PReconstruction::closeClustersOutput()
{
  LOGP(warning, "closeClustersOutput called for {}, this should not happen", getName());
}

void NA6PReconstruction::createVerticesOutput()
{
  LOGP(warning, "createVerticesOutput called for {}, this should not happen", getName());
}

void NA6PReconstruction::clearVertices()
{
  LOGP(warning, "clearVertices called for {}, this should not happen", getName());
}

void NA6PReconstruction::writeVertices()
{
  LOGP(warning, "writeVertices called for {}, this should not happen", getName());
}

void NA6PReconstruction::closeVerticesOutput()
{
  LOGP(warning, "closeVerticesOutput called for {}, this should not happen", getName());
}

void NA6PReconstruction::createTracksOutput()
{
  LOGP(warning, "createTracksOutput called for {}, this should not happen", getName());
}

void NA6PReconstruction::clearTracks()
{
  LOGP(warning, "clearTracks called for {}, this should not happen", getName());
}

void NA6PReconstruction::writeTracks()
{
  LOGP(warning, "writeTracks called for {}, this should not happen", getName());
}

void NA6PReconstruction::closeTracksOutput()
{
  LOGP(warning, "closeTracksOutput called for {}, this should not happen", getName());
}

void NA6PReconstruction::assignMCLabels(std::vector<NA6PTrack>& trk,
                                        std::vector<NA6PMCComposedLabel>& mcTrkLabels,
                                        const NA6PMCTruthContainer& mcCluLabels)
{
  mcTrkLabels.clear();
  mcTrkLabels.reserve(trk.size());
  for (auto& tr : trk) {
    // assign label using MCTruthContainer
    std::vector<std::pair<NA6PMCComposedLabel, int>> countLabs;
    for (int jLay = 0; jLay < NA6PTrack::kMaxLr; ++jLay) {
      int cluID = tr.getClusterIndex(jLay);
      if (cluID >= 0) {
        std::span labels = mcCluLabels.getLabels(cluID);
        int nLabels = labels.size();
        for (int jLab = 0; jLab < nLabels; jLab++) {
          NA6PMCComposedLabel lbl = labels[jLab];
          bool found = false;
          for (auto& p : countLabs) {
            if (p.first == lbl) {
              ++p.second;
              found = true;
              break;
            }
          }
          if (!found) {
            countLabs.push_back({lbl, 1});
          }
        }
      }
    }
    NA6PMCComposedLabel lblTrack;
    lblTrack.unset();
    int maxCountLabs = 0;
    for (const auto& p : countLabs) {
      if (p.second > maxCountLabs) {
        maxCountLabs = p.second;
        lblTrack = p.first;
      }
    }
    // assign fake label to tracks with misassociations
    if (lblTrack.isSet()) {
      for (int jLay = 0; jLay < NA6PTrack::kMaxLr; ++jLay) {
        int cluID = tr.getClusterIndex(jLay);
        if (cluID >= 0) {
          std::span labels = mcCluLabels.getLabels(cluID);
          int nLabels = labels.size();
          bool found = false;
          for (int jLab = 0; jLab < nLabels; jLab++) {
            NA6PMCComposedLabel lbl = labels[jLab];
            if (lbl == lblTrack) {
              found = true;
              break;
            }
          }
          if (!found) {
            lblTrack.setFakeFlag();
            break;
          }
        }
      }
    }
    mcTrkLabels.push_back(lblTrack);
  }
}
