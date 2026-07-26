#ifndef NA6P_EVENT_READER_H
#define NA6P_EVENT_READER_H

#include "TFile.h"
#include "TParticle.h"
#include "TTree.h"

#include "NA6PMatch.h"
#include "NA6PMCComposedLabel.h"
#include "NA6PMCEventHeader.h"
#include "NA6PTrack.h"
#include "NA6PVertex.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

class NA6PEventReader
{
 public:
  using TrackIndex = std::size_t;
  using LabelToTrackIndices = std::unordered_map<NA6PMCComposedLabel, std::vector<TrackIndex>>;

  NA6PEventReader(
    const char* fileNameVerTel = "VerticesVerTel.root",
    const char* fileNameTracksVerTel = "TracksVerTel.root",
    const char* fileNameTracksMuonSpec = "TracksMuonSpec.root",
    const char* fileNameTracksMatching = "TracksMatching.root",
    const char* fileNameMC = "MCKine.root",
    bool readVerTel = true,
    bool readTracksVerTel = true,
    bool readTracksMuonSpec = true,
    bool readTracksMatching = true,
    bool readMC = true)
  {
    mReadMC = readMC;

    if (readVerTel) {
      openFile(mFileVerTel, fileNameVerTel);
      setupTree(mFileVerTel.get(), "verticesVerTel", mTreeVerTel);
    }

    if (readTracksVerTel) {
      openFile(mFileTracksVerTel, fileNameTracksVerTel);
      setupTree(mFileTracksVerTel.get(), "tracksVerTel", mTreeTracksVerTel);
    }

    if (readTracksMuonSpec) {
      openFile(mFileTracksMuonSpec, fileNameTracksMuonSpec);
      setupTree(mFileTracksMuonSpec.get(), "tracksMuonSpec", mTreeTracksMuonSpec);
    }

    if (readTracksMatching) {
      openFile(mFileTracksMatching, fileNameTracksMatching);
      setupTree(mFileTracksMatching.get(), "tracksMatching", mTreeTracksMatching);
    }

    if (readMC) {
      openFile(mFileMC, fileNameMC);
      setupTree(mFileMC.get(), "mckine", mTreeMC);
    }

    // Vertex telescope vertices
    if (mTreeVerTel) {
      setupRequiredBranch(mTreeVerTel, "VerTel", mVerticesVerTel);
    }

    // Vertex telescope tracks and corresponding MC labels
    if (mTreeTracksVerTel) {
      setupRequiredBranch(mTreeTracksVerTel, "VerTel", mTracksVerTel);
      if (readMC) {
        setupMCLabelBranch(mTreeTracksVerTel, "VerTelMCTruth", mTrackLabelsVerTel);
      }
    }

    // Muon spectrometer tracks and corresponding MC labels
    if (mTreeTracksMuonSpec) {
      setupRequiredBranch(mTreeTracksMuonSpec, "MuonSpec", mTracksMuonSpec);
      if (readMC) {
        setupMCLabelBranch(mTreeTracksMuonSpec, "MuonSpecMCTruth", mTrackLabelsMuonSpec);
      }
    }

    // Matched tracks and corresponding MC labels
    if (mTreeTracksMatching) {
      setupRequiredBranch(mTreeTracksMatching, "Matching", mMatches);
      if (readMC) {
        setupMCLabelBranch(mTreeTracksMatching, "MatchingMCTruth", mMatchLabels);
      }
    }

    // Generator-level MC information
    if (mTreeMC) {
      setupRequiredBranch(mTreeMC, "header", mMCHeader);
      setupRequiredBranch(mTreeMC, "tracks", mMCParticles);
    }

    determineEntries();
  }

  std::int64_t entries() const { return mEntries; }
  std::int64_t currentEntry() const { return mCurrentEntry; }

  bool loadEvent(std::int64_t i)
  {
    if (i < 0 || i >= mEntries) {
      return false;
    }

    if (mTreeMC) {
      mTreeMC->GetEntry(i);
    }

    if (mTreeVerTel) {
      mTreeVerTel->GetEntry(i);
    }

    if (mTreeTracksVerTel) {
      mTreeTracksVerTel->GetEntry(i);
    }

    if (mTreeTracksMuonSpec) {
      mTreeTracksMuonSpec->GetEntry(i);
    }

    if (mTreeTracksMatching) {
      mTreeTracksMatching->GetEntry(i);
    }

    mCurrentEntry = i;

    mVTTrackIndicesByLabel.clear();
    mMSTrackIndicesByLabel.clear();

    if (mReadMC) {
      buildTrackLabelMaps(); // after GetEntry(i);
    }
    return true;
  }

  // ------------------------------------------------------------------
  // Availability
  // ------------------------------------------------------------------

  bool hasVerticesVerTel() const { return mTreeVerTel != nullptr; }
  bool hasTracksVerTel() const { return mTreeTracksVerTel != nullptr; }
  bool hasTracksMuonSpec() const { return mTreeTracksMuonSpec != nullptr; }
  bool hasMatching() const { return mTreeTracksMatching != nullptr; }
  bool hasMC() const { return mTreeMC != nullptr; }
  bool hasVTTrackMCLabels() const { return mTrackLabelsVerTel != nullptr; }
  bool hasMSTrackMCLabels() const { return mTrackLabelsMuonSpec != nullptr; }
  bool hasMatchedTrackMCLabels() const { return mMatchLabels != nullptr; }

  // ------------------------------------------------------------------
  // Event-level getters
  // ------------------------------------------------------------------

  NA6PMCEventHeader* mcHeader() const { return mMCHeader; }

  const std::vector<NA6PVertex>& verticesVerTel() const
  {
    static const std::vector<NA6PVertex> empty;
    return mVerticesVerTel ? *mVerticesVerTel : empty;
  }

  const std::vector<NA6PTrack>& tracksVerTel() const
  {
    static const std::vector<NA6PTrack> empty;
    return mTracksVerTel ? *mTracksVerTel : empty;
  }

  const std::vector<NA6PTrack>& tracksMuonSpec() const
  {
    static const std::vector<NA6PTrack> empty;
    return mTracksMuonSpec ? *mTracksMuonSpec : empty;
  }

  const std::vector<NA6PMatch>& matches() const
  {
    static const std::vector<NA6PMatch> empty;
    return mMatches ? *mMatches : empty;
  }

  const std::vector<TParticle>& mcParticles() const
  {
    static const std::vector<TParticle> empty;
    return mMCParticles ? *mMCParticles : empty;
  }

  const std::vector<NA6PMCComposedLabel>& trackLabelsVerTel() const
  {
    static const std::vector<NA6PMCComposedLabel> empty;
    return mTrackLabelsVerTel ? *mTrackLabelsVerTel : empty;
  }

  const std::vector<NA6PMCComposedLabel>& trackLabelsMuonSpec() const
  {
    static const std::vector<NA6PMCComposedLabel> empty;
    return mTrackLabelsMuonSpec ? *mTrackLabelsMuonSpec : empty;
  }

  const std::vector<NA6PMCComposedLabel>& matchLabels() const
  {
    static const std::vector<NA6PMCComposedLabel> empty;
    return mMatchLabels ? *mMatchLabels : empty;
  }

  // ------------------------------------------------------------------
  // Individual reconstructed objects
  // ------------------------------------------------------------------

  const NA6PVertex& getVTVertex(std::size_t i) const
  {
    return verticesVerTel().at(i); // if index is out of range, at returns std::out_of_range;
  }

  const NA6PTrack& getVTTrack(std::size_t i) const
  {
    return tracksVerTel().at(i); // if index is out of range, at returns std::out_of_range;
  }

  const NA6PTrack& getMSTrack(std::size_t i) const
  {
    return tracksMuonSpec().at(i); // if index is out of range, at returns std::out_of_range;
  }

  const NA6PMatch& getMatchedTrack(std::size_t i) const
  {
    return matches().at(i); // if index is out of range, at returns std::out_of_range;
  }

  const NA6PTrack* getVTTrack(const NA6PMCComposedLabel& label) const
  {
    const auto* indices = getVTTrackIndices(label);
    if (!indices || indices->empty()) {
      return nullptr;
    }
    return &getVTTrack(indices->front());
  }

  const NA6PTrack* getMSTrack(const NA6PMCComposedLabel& label) const
  {
    const auto* indices = getMSTrackIndices(label);
    if (!indices || indices->empty()) {
      return nullptr;
    }
    return &getMSTrack(indices->front());
  }

  // ------------------------------------------------------------------
  // MC-label access
  // ------------------------------------------------------------------

  const NA6PMCComposedLabel& getVTTrackLabel(std::size_t i) const
  {
    return trackLabelsVerTel().at(i);
  }

  const NA6PMCComposedLabel& getMSTrackLabel(std::size_t i) const
  {
    return trackLabelsMuonSpec().at(i);
  }

  const NA6PMCComposedLabel& getMatchedTrackLabel(std::size_t i) const
  {
    return matchLabels().at(i);
  }

  const LabelToTrackIndices& vtTrackIndicesByLabel() const { return mVTTrackIndicesByLabel; }
  const LabelToTrackIndices& msTrackIndicesByLabel() const { return mMSTrackIndicesByLabel; }

  const std::vector<TrackIndex>* getVTTrackIndices(const NA6PMCComposedLabel& label) const
  {
    const auto it = mVTTrackIndicesByLabel.find(label);
    return it != mVTTrackIndicesByLabel.end() ? &it->second : nullptr;
  }

  const std::vector<TrackIndex>* getMSTrackIndices(const NA6PMCComposedLabel& label) const
  {
    const auto it = mMSTrackIndicesByLabel.find(label);
    return it != mMSTrackIndicesByLabel.end() ? &it->second : nullptr;
  }

  // ------------------------------------------------------------------
  // MC-particle access
  // ------------------------------------------------------------------

  const TParticle& getMCParticle(std::size_t i) const
  {
    return mcParticles().at(i); // if index is out of range, at returns std::out_of_range;
  }

  const TParticle* getMCParticle(const NA6PMCComposedLabel& label) const
  {
    if (!label.isValid()) {
      return nullptr;
    }

    /*
     * Do not use getTrackIDSigned() here.
     *
     * getTrackIDSigned() returns a negative value for fake labels.
     * The index of the TParticle vector must always be obtained using
     * getTrackID().
     */
    return &getMCParticle(label.getTrackID());
  }

 private:
  static void openFile(std::unique_ptr<TFile>& file, const char* fileName)
  {
    file.reset(TFile::Open(fileName, "READ"));

    if (!file || file->IsZombie()) {
      throw std::runtime_error(std::string("Cannot open file: ") + fileName);
    }
  }

  static void setupTree(TFile* file, const char* treeName, TTree*& tree)
  {
    tree = nullptr;

    if (!file) {
      return;
    }

    file->GetObject(treeName, tree);

    if (!tree) {
      throw std::runtime_error(std::string("Cannot find tree: ") + treeName);
    }
  }

  template <typename T>
  static void setupRequiredBranch(TTree* tree, const char* branchName, T*& object)
  {
    object = nullptr;

    if (!tree) {
      throw std::runtime_error(std::string("Null tree while setting branch: ") + branchName);
    }

    if (!tree->GetBranch(branchName)) {
      throw std::runtime_error(std::string("Cannot find branch '") + branchName + "' in tree '" + tree->GetName() + "'");
    }

    const int status = tree->SetBranchAddress(branchName, &object);

    if (status < 0) {
      throw std::runtime_error(std::string("Cannot set branch address for '") + branchName + "' in tree '" + tree->GetName() + "'");
    }
  }

  static bool setupMCLabelBranch(TTree* tree, const char* branchName, std::vector<NA6PMCComposedLabel>*& labels)
  {
    labels = nullptr;

    if (!tree) {
      return false;
    }

    if (!tree->GetBranch(branchName)) {
      std::cerr << "WARNING: tree '" << tree->GetName() << "' does not contain MC-label branch '" << branchName << "'" << std::endl;
      return false;
    }

    const int status = tree->SetBranchAddress(branchName, &labels);

    if (status < 0) {
      throw std::runtime_error(std::string("Cannot set branch address for '") + branchName + "' in tree '" + tree->GetName() + "'");
    }

    return true;
  }

  void buildTrackLabelMaps()
  {
    mVTTrackIndicesByLabel.clear();
    mMSTrackIndicesByLabel.clear();
    buildTrackLabelMap(tracksVerTel().size(), trackLabelsVerTel(), mVTTrackIndicesByLabel, "VT");
    buildTrackLabelMap(tracksMuonSpec().size(), trackLabelsMuonSpec(), mMSTrackIndicesByLabel, "MS");
  }

  static void buildTrackLabelMap(std::size_t numberOfTracks, const std::vector<NA6PMCComposedLabel>& labels, LabelToTrackIndices& output, const char* detectorName)
  {
    output.clear();

    if (numberOfTracks != labels.size()) {
      std::cerr << "WARNING: " << detectorName << " track/label size mismatch: tracks=" << numberOfTracks << ", labels=" << labels.size() << std::endl;
    }

    const std::size_t numberOfEntries = std::min(numberOfTracks, labels.size());

    for (std::size_t trackIndex = 0; trackIndex < numberOfEntries; ++trackIndex) {
      const auto& label = labels[trackIndex];
      if (!label.isValid() || label.isFake()) {
        continue;
      }
      output[label].push_back(trackIndex);
    }
  }

  void determineEntries()
  {
    std::vector<std::int64_t> numbersOfEntries;

    if (mTreeMC) {
      numbersOfEntries.push_back(mTreeMC->GetEntries());
    }

    if (mTreeVerTel) {
      numbersOfEntries.push_back(mTreeVerTel->GetEntries());
    }

    if (mTreeTracksVerTel) {
      numbersOfEntries.push_back(mTreeTracksVerTel->GetEntries());
    }

    if (mTreeTracksMuonSpec) {
      numbersOfEntries.push_back(mTreeTracksMuonSpec->GetEntries());
    }

    if (mTreeTracksMatching) {
      numbersOfEntries.push_back(mTreeTracksMatching->GetEntries());
    }

    if (numbersOfEntries.empty()) {
      throw std::runtime_error("No input trees were loaded");
    }

    const auto [minIt, maxIt] = std::minmax_element(numbersOfEntries.begin(), numbersOfEntries.end());

    if (*minIt != *maxIt) {
      std::cerr << "WARNING: input trees have different numbers of entries. Using minimum entries = " << *minIt << ", maximum entries = " << *maxIt << std::endl;
    }

    mEntries = *minIt;
  }

 private:
  // Files
  std::unique_ptr<TFile> mFileVerTel;
  std::unique_ptr<TFile> mFileTracksVerTel;
  std::unique_ptr<TFile> mFileTracksMuonSpec;
  std::unique_ptr<TFile> mFileTracksMatching;
  std::unique_ptr<TFile> mFileMC;
  bool mReadMC{false};

  // Trees
  TTree* mTreeVerTel = nullptr;
  TTree* mTreeTracksVerTel = nullptr;
  TTree* mTreeTracksMuonSpec = nullptr;
  TTree* mTreeTracksMatching = nullptr;
  TTree* mTreeMC = nullptr;

  // Reconstructed objects
  std::vector<NA6PVertex>* mVerticesVerTel = nullptr;
  std::vector<NA6PTrack>* mTracksVerTel = nullptr;
  std::vector<NA6PTrack>* mTracksMuonSpec = nullptr;
  std::vector<NA6PMatch>* mMatches = nullptr;

  // Reconstructed-object MC labels
  std::vector<NA6PMCComposedLabel>* mTrackLabelsVerTel = nullptr;
  std::vector<NA6PMCComposedLabel>* mTrackLabelsMuonSpec = nullptr;
  std::vector<NA6PMCComposedLabel>* mMatchLabels = nullptr;

  // Generator-level MC objects
  std::vector<TParticle>* mMCParticles = nullptr;
  NA6PMCEventHeader* mMCHeader = nullptr;

  LabelToTrackIndices mVTTrackIndicesByLabel;
  LabelToTrackIndices mMSTrackIndicesByLabel;

  // Entry bookkeeping
  std::int64_t mEntries = 0;
  std::int64_t mCurrentEntry = -1;
};

#endif // NA6P_EVENT_READER_H
