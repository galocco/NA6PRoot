#include "NA6PMCEventHeader.h"
#include "NA6PMuonSpecModularHit.h"
#include "NA6PVerTelHit.h"
#include <TFile.h>
#include <TParticle.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

R__LOAD_LIBRARY(libsimLib)

namespace
{
struct EmbeddingInput {
  TFile kine, vt, ms;
  TTreeReader kr{"mckine", &kine}, vr{"hitsVerTel", &vt}, mr{"hitsMuonSpecModular", &ms};
  TTreeReaderValue<NA6PMCEventHeader> header{kr, "header"};
  TTreeReaderValue<std::vector<TParticle>> tracks{kr, "tracks"};
  TTreeReaderValue<std::vector<NA6PVerTelHit>> vtHits{vr, "VerTel"};
  TTreeReaderValue<std::vector<NA6PMuonSpecModularHit>> msHits{mr, "MuonSpecModular"};
  std::vector<int> owner, remap;
  std::vector<NA6PMCGenHeader> generators;

  explicit EmbeddingInput(const std::string& dir)
    : kine((dir + "/MCKine.root").c_str()), vt((dir + "/HitsVerTel.root").c_str()),
      ms((dir + "/HitsMuonSpecModular.root").c_str())
  {
    if (entries() <= 0 || vr.GetEntries() != entries() || mr.GetEntries() != entries()) {
      throw std::runtime_error("Missing or inconsistent input trees in " + dir);
    }
  }

  Long64_t entries() { return kr.GetEntries(); }

  void load(Long64_t entry, bool removeThermal)
  {
    if (kr.SetEntry(entry) != TTreeReader::kEntryValid || vr.SetEntry(entry) != TTreeReader::kEntryValid ||
        mr.SetEntry(entry) != TTreeReader::kEntryValid || header->getNTracks() != tracks->size() ||
        header->getNPrimaries() > tracks->size()) {
      throw std::runtime_error("Cannot read embedding event");
    }
    owner.assign(tracks->size(), -2);
    remap.resize(tracks->size());
    generators.clear();
    unsigned primaries = 0;
    for (const auto& generator : header->getGenHeaders()) {
      const auto count = generator.getNPrimaries();
      if (!count) continue;
      if (generator.getPrimariesOffset() != primaries || count > header->getNPrimaries() - primaries) {
        throw std::runtime_error("Generator primary offsets do not cover the primary block");
      }
      const auto& name = generator.getInfo();
      const bool thermal = name.find("VirtualPhoton_pdg23") == 0 ||
                           name.find("thermalMuPlus_pdg-13") == 0 || name.find("thermalMuMinus_pdg13") == 0;
      if (removeThermal && thermal) {
        for (unsigned i = 0; i < count; ++i) owner[primaries++] = -1;
        continue;
      }
      const int generatorID = generators.size();
      auto& copy = generators.emplace_back(generator);
      copy.setNPrimaries(0);
      copy.setNSecondaries(0);
      copy.setPrimariesOffset(0);
      copy.setSecondariesOffset(std::numeric_limits<uint32_t>::max());
      for (unsigned i = 0; i < count; ++i) owner[primaries++] = generatorID;
    }
    if (primaries != header->getNPrimaries()) {
      throw std::runtime_error("Missing generator metadata for an embedding primary");
    }
    for (size_t i = primaries; i < owner.size(); ++i) {
      const int mother = tracks->at(i).GetFirstMother();
      if (mother < 0 || static_cast<size_t>(mother) >= i) {
        throw std::runtime_error("Unexpected MC ancestry ordering");
      }
      owner[i] = owner[mother];
    }
  }

  void map(bool primary, unsigned& next)
  {
    const size_t begin = primary ? 0 : header->getNPrimaries();
    const size_t end = primary ? header->getNPrimaries() : owner.size();
    for (size_t i = begin; i < end; ++i) {
      remap[i] = -1;
      if (owner[i] < 0) continue;
      remap[i] = next++;
    }
  }

  void append(bool primary, std::vector<TParticle>& output)
  {
    const size_t begin = primary ? 0 : header->getNPrimaries();
    const size_t end = primary ? header->getNPrimaries() : owner.size();
    for (size_t i = begin; i < end; ++i) {
      if (remap[i] < 0) continue;
      auto& generator = generators[owner[i]];
      if (primary) {
        if (generator.getNPrimaries() == 0) generator.setPrimariesOffset(output.size());
        generator.incNPrimaries();
      } else {
        if (generator.getNSecondaries() == 0) generator.setSecondariesOffset(output.size());
        generator.incNSecondaries();
      }
      auto& p = output.emplace_back(tracks->at(i));
      // SecondMother stores the original stack ID in NA6PMC, not a saved-tree index.
      if (p.GetFirstMother() >= 0) p.SetFirstMother(remap.at(p.GetFirstMother()));
      if (p.GetFirstDaughter() >= 0) p.SetFirstDaughter(remap.at(p.GetFirstDaughter()));
      if (p.GetLastDaughter() >= 0) p.SetLastDaughter(remap.at(p.GetLastDaughter()));
    }
  }
};

template <typename Hit>
void appendHits(const std::vector<Hit>& input, const std::vector<int>& remap, std::vector<Hit>& output)
{
  for (const auto& hit : input) {
    const int id = remap.at(hit.getTrackID());
    if (id >= 0) output.emplace_back(hit).setTrackID(id);
  }
}
}

void embedThermalHits(const char* background, const char* signal, const char* output, int first, int replicas)
{
  if (first < 0 || replicas <= 0) throw std::runtime_error("Invalid embedding event range");
  EmbeddingInput bg(background), sig(signal);
  const std::string dir = output;
  TFile kine((dir + "/MCKine.root").c_str(), "RECREATE");
  TFile vt((dir + "/HitsVerTel.root").c_str(), "RECREATE");
  TFile ms((dir + "/HitsMuonSpecModular.root").c_str(), "RECREATE");
  if (kine.IsZombie() || vt.IsZombie() || ms.IsZombie()) throw std::runtime_error("Cannot create embedding outputs");
  NA6PMCEventHeader header;
  std::vector<TParticle> tracks;
  std::vector<NA6PVerTelHit> vtHits;
  std::vector<NA6PMuonSpecModularHit> msHits;
  kine.cd();
  TTree kr("mckine", "Embedded kinematics");
  kr.Branch("header", &header);
  kr.Branch("tracks", &tracks);
  vt.cd();
  TTree vr("hitsVerTel", "Embedded VerTel hits");
  vr.Branch("VerTel", &vtHits);
  ms.cd();
  TTree mr("hitsMuonSpecModular", "Embedded MuonSpecModular hits");
  mr.Branch("MuonSpecModular", &msHits);

  Long64_t lastBackground = -1;
  for (Long64_t i = 0; i < sig.entries(); ++i) {
    const auto index = (first + i) / replicas % bg.entries();
    if (index != lastBackground) {
      bg.load(index, true);
      lastBackground = index;
    }
    sig.load(i, false);
    if (sig.header->getNPrimaries() != 3 || sig.generators.size() != 1 ||
        sig.generators.front().getInfo().find("VirtualPhoton_pdg23") != 0) {
      throw std::runtime_error("Expected one virtual photon and its two primary muons");
    }
    if (sig.header->getEventID() != first + i ||
        sig.header->getVX() != bg.header->getVX() || sig.header->getVY() != bg.header->getVY() ||
        sig.header->getVZ() != bg.header->getVZ()) {
      throw std::runtime_error("Signal event or vertex does not match the selected background");
    }
    header = *bg.header;
    header.getGenHeaders().clear();
    header.setEventID(sig.header->getEventID());
    tracks.clear();
    vtHits.clear();
    msHits.clear();
    unsigned next = 0;
    bg.map(true, next);
    sig.map(true, next);
    header.setNPrimaries(next);
    bg.map(false, next);
    sig.map(false, next);
    header.setNTracks(next);
    for (bool primary : {true, false}) {
      bg.append(primary, tracks);
      sig.append(primary, tracks);
    }
    for (const auto& generator : bg.generators) header.addGenHeader(generator);
    for (const auto& generator : sig.generators) header.addGenHeader(generator);
    auto& provenance = header.getGenHeaders().emplace_back(0, 0, 0, 0, "EmbeddedBackground");
    provenance.addKeyVal("entry", index);
    provenance.addKeyVal("eventID", bg.header->getEventID());
    appendHits(*bg.vtHits, bg.remap, vtHits);
    appendHits(*sig.vtHits, sig.remap, vtHits);
    appendHits(*bg.msHits, bg.remap, msHits);
    appendHits(*sig.msHits, sig.remap, msHits);
    if (kr.Fill() < 0 || vr.Fill() < 0 || mr.Fill() < 0) throw std::runtime_error("Failed to write embedded event");
    // Background can be reused on the next iteration, including its filtered owner map.
    for (auto& generator : bg.generators) {
      generator.setNPrimaries(0);
      generator.setNSecondaries(0);
      generator.setSecondariesOffset(std::numeric_limits<uint32_t>::max());
    }
  }
  kine.cd();
  if (kr.Write() <= 0) throw std::runtime_error("Cannot write mixed kinematics tree");
  vt.cd();
  if (vr.Write() <= 0) throw std::runtime_error("Cannot write mixed VT hits tree");
  ms.cd();
  if (mr.Write() <= 0) throw std::runtime_error("Cannot write mixed MS hits tree");
  std::cout << "Embedded " << sig.entries() << " events starting at " << first << '\n';
}
