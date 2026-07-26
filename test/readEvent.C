#include <cstdint>
#include "NA6PEventReader.h"

void readEvent(
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
  NA6PEventReader reader(fileNameVerTel, fileNameTracksVerTel, fileNameTracksMuonSpec, fileNameTracksMatching, fileNameMC, readVerTel, readTracksVerTel, readTracksMuonSpec, readTracksMatching, readMC);

  for (std::int64_t iev = 0; iev < reader.entries(); ++iev) {
    reader.loadEvent(iev);

    const auto& vertices = reader.verticesVerTel();
    const auto& vtTracks = reader.tracksVerTel();
    const auto& msTracks = reader.tracksMuonSpec();
    const auto& matches = reader.matches();
    const auto& mcParts = reader.mcParticles();
    const auto& mcHeader = reader.mcHeader();

    std::cout << "event " << iev << " nVtx=" << vertices.size() << " nVT=" << vtTracks.size() << " nMS=" << msTracks.size() << " nMatch=" << matches.size() << " nMC=" << mcParts.size() << " mcVZ=" << mcHeader->getVZ() << std::endl;

    for (std::size_t iTrack = 0; iTrack < vtTracks.size(); ++iTrack) {
      const auto& track = vtTracks[iTrack];
      const auto& label = reader.getVTTrackLabel(iTrack);

      if (!label.isValid() || label.isFake()) {
        continue;
      }

      const auto* mcParticle = reader.getMCParticle(label);
      if (!mcParticle) {
        continue;
      }

      const int pdgCode = mcParticle->GetPdgCode();

      printf("iTrack = %zu, track.getPt() = %f, mcParticle->Pt() = %f, pdgCode = %d\n", iTrack, track.getPt(), mcParticle->Pt(), pdgCode);
    } // end of track loop
  } // end of mc collision loop
}
