#include "NA6PMCEventHeader.h"
#include "NA6PMCTruthContainer.h"
#include "NA6PTrack.h"
#include "NA6PMatch.h"
#include <TFile.h>
#include <TParticle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
R__LOAD_LIBRARY(librecLib)

void auditFakeMSOrigins(const char* directory)
{
  const std::string run=directory;
  TFile vtFile((run+"/reco/TracksVerTel.root").c_str()), msFile((run+"/reco/TracksMuonSpec.root").c_str());
  TFile matchFile((run+"/reco_5hits/TracksMatching.root").c_str()), mcFile((run+"/sim/MCKine.root").c_str());
  TFile cluFile((run+"/reco/ClustersMuonSpec.root").c_str());
  TTreeReader vr("tracksVerTel",&vtFile), mr("tracksMuonSpec",&msFile), rr("tracksMatching",&matchFile);
  TTreeReader cr("mckine",&mcFile), hr("clustersMuonSpec",&cluFile);
  TTreeReaderValue<std::vector<NA6PMCComposedLabel>> vl(vr,"VerTelMCTruth"), ml(mr,"MuonSpecMCTruth");
  TTreeReaderValue<std::vector<NA6PTrack>> ms(mr,"MuonSpec");
  TTreeReaderValue<std::vector<NA6PMatch>> matches(rr,"Matching");
  TTreeReaderValue<NA6PMCEventHeader> header(cr,"header");
  TTreeReaderValue<std::vector<TParticle>> particles(cr,"tracks");
  TTreeReaderValue<NA6PMCTruthContainer> truth(hr,"MuonSpecMCTruth");
  enum Category {Thermal, PiDecay, KDecay, OtherMuon, OtherParticle, Missing, NCategories};
  std::array<int,NCategories> tracksWith{}, clustersWith{}, labelledContributions{};
  int events=0, selected=0, fakeLabel=0, missingLabel=0, clusterCount=0;
  auto genuine=[](const NA6PMCComposedLabel& l){return l.isValid() && !l.isFake();};
  while(cr.Next()) {
    if(!vr.Next()||!mr.Next()||!rr.Next()||!hr.Next()) throw std::runtime_error("Mismatched event trees");
    std::vector<bool> thermal(particles->size(),false);
    for(const auto& gen:header->getGenHeaders()) {
      const auto& name=gen.getInfo();
      if(name.find("VirtualPhoton_pdg23")!=0 && name.find("thermalMuPlus_pdg-13")!=0 && name.find("thermalMuMinus_pdg13")!=0) continue;
      size_t end=gen.getPrimariesOffset()+gen.getNPrimaries();
      if(end>particles->size()) throw std::runtime_error("MC generator range");
      for(size_t i=gen.getPrimariesOffset();i<end;++i) if(std::abs(particles->at(i).GetPdgCode())==13) thermal[i]=true;
    }
    for(const auto& match:*matches) {
      int vi=match.getIndexVT(), mi=match.getIndexMS();
      if(vi<0||mi<0||size_t(vi)>=vl->size()||size_t(mi)>=ml->size()||size_t(mi)>=ms->size()) throw std::runtime_error("Invalid match index");
      const auto& v=vl->at(vi); const auto& label=ml->at(mi);
      if(genuine(v)&&genuine(label)&&v==label) continue;
      if(genuine(label)) continue;
      ++selected;
      if(!label.isValid()) ++missingLabel; else ++fakeLabel;
      std::array<bool,NCategories> trackFlags{};
      const auto& track=ms->at(mi);
      for(int layer=0;layer<NA6PTrack::kMaxLr;++layer) {
        int index=track.getClusterIndex(layer);
        if(index<0) continue;
        if(size_t(index)>=truth->getIndexedSize()) throw std::runtime_error("Invalid cluster truth index");
        ++clusterCount;
        std::array<bool,NCategories> clusterFlags{};
        auto labels=truth->getLabels(index);
        if(labels.empty()) clusterFlags[Missing]=true;
        for(const auto& l:labels) {
          int id=l.getTrackID();
          if(!l.isValid()||l.getSourceID()!=0||l.getEventID()!=header->getEventID()||id<0||size_t(id)>=particles->size()) {
            clusterFlags[Missing]=true; continue;
          }
          int pdg=std::abs(particles->at(id).GetPdgCode());
          Category category=OtherParticle;
          if(pdg==13) {
            if(thermal[id]) category=Thermal;
            else {
              category=OtherMuon;
              int mother=particles->at(id).GetFirstMother();
              int depth=0;
              while(mother>=0 && size_t(mother)<particles->size() && depth++<100) {
                int ancestor=std::abs(particles->at(mother).GetPdgCode());
                if(ancestor==211||ancestor==111) {category=PiDecay; break;}
                if(ancestor==321||ancestor==311||ancestor==130||ancestor==310) {category=KDecay; break;}
                mother=particles->at(mother).GetFirstMother();
              }
            }
          }
          clusterFlags[category]=true;
          ++labelledContributions[category];
        }
        for(int i=0;i<NCategories;++i) if(clusterFlags[i]) {++clustersWith[i]; trackFlags[i]=true;}
      }
      for(int i=0;i<NCategories;++i) if(trackFlags[i]) ++tracksWith[i];
    }
    ++events;
  }
  if(vr.Next()||mr.Next()||rr.Next()||hr.Next()) throw std::runtime_error("Extra tree entries");
  const std::array<const char*,NCategories> names{"thermal_muon","muon_from_pion","muon_from_kaon","other_muon","non_muon","unidentified"};
  std::ofstream out(run+"/impact_parameter/fakeMSOrigins.csv");
  out<<"category,tracks_with_contribution,clusters_with_contribution,cluster_label_contributions\n";
  for(int i=0;i<NCategories;++i) {
    out<<names[i]<<","<<tracksWith[i]<<","<<clustersWith[i]<<","<<labelledContributions[i]<<"\n";
    std::cout<<names[i]<<" tracks="<<tracksWith[i]<<" clusters="<<clustersWith[i]<<" labels="<<labelledContributions[i]<<"\n";
  }
  std::cout<<"events="<<events<<" selected="<<selected<<" fake_labels="<<fakeLabel<<" invalid_labels="<<missingLabel<<" clusters="<<clusterCount<<"\n";
}
