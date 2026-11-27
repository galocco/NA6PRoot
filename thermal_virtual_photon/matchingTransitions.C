#include "NA6PMCEventHeader.h"
#include "NA6PMCComposedLabel.h"
#include "NA6PTrack.h"
#include "NA6PMatch.h"
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TParticle.h>
#include <TSystem.h>
#include <array>
#include <fstream>
#include <stdexcept>
#include <vector>
R__LOAD_LIBRARY(librecLib)

void matchingTransitions(const char* runDir = ".")
{
  std::string path = runDir;
  TFile vf((path + "/reco/TracksVerTel.root").c_str()), mf((path + "/reco/TracksMuonSpec.root").c_str());
  TFile mc((path + "/sim/MCKine.root").c_str()), f3((path + "/reco/TracksMatching.root").c_str());
  TFile f4((path + "/reco_4hits/TracksMatching.root").c_str()), f5((path + "/reco_5hits/TracksMatching.root").c_str());
  TTreeReader vr("tracksVerTel",&vf), mr("tracksMuonSpec",&mf), cr("mckine",&mc);
  TTreeReader r3("tracksMatching",&f3), r4("tracksMatching",&f4), r5("tracksMatching",&f5);
  TTreeReaderValue<std::vector<NA6PTrack>> v(vr,"VerTel"), m(mr,"MuonSpec");
  TTreeReaderValue<std::vector<NA6PMCComposedLabel>> vl(vr,"VerTelMCTruth"), ml(mr,"MuonSpecMCTruth");
  TTreeReaderValue<NA6PMCEventHeader> header(cr,"header");
  TTreeReaderValue<std::vector<TParticle>> particles(cr,"tracks");
  TTreeReaderValue<std::vector<NA6PMatch>> a(r3,"Matching"), b(r4,"Matching"), c(r5,"Matching");
  constexpr int nStates=7;
  std::array<std::array<std::array<int,nStates>,nStates>,3> flow{};
  std::array<std::array<int,nStates>,3> populations{};
  int events=0, selected=0;
  gSystem->mkdir((path + "/transitions").c_str(), true);
  while(cr.Next()) {
    if(!vr.Next() || !mr.Next() || !r3.Next() || !r4.Next() || !r5.Next()) throw std::runtime_error("Mismatched event trees");
    if(v->size()!=vl->size() || m->size()!=ml->size()) throw std::runtime_error("Mismatched track labels");
    std::vector<bool> thermal(particles->size(),false);
    bool found=false;
    for(const auto& gen:header->getGenHeaders()) {
      const auto& name=gen.getInfo();
      if(name.find("VirtualPhoton_pdg23")!=0 && name.find("thermalMuPlus_pdg-13")!=0 && name.find("thermalMuMinus_pdg13")!=0) continue;
      found=true;
      size_t end=gen.getPrimariesOffset()+gen.getNPrimaries();
      if(end>particles->size()) throw std::runtime_error("Invalid MC generator range");
      for(size_t i=gen.getPrimariesOffset();i<end;++i) if(std::abs(particles->at(i).GetPdgCode())==13) thermal[i]=true;
    }
    if(!found) throw std::runtime_error("Missing thermal generator");
    auto genuine=[](const NA6PMCComposedLabel& label){return label.isValid() && !label.isFake();};
    std::array<std::vector<int>,3> state{std::vector<int>(m->size()),std::vector<int>(m->size()),std::vector<int>(m->size())};
    std::array<std::vector<bool>,3> seen{std::vector<bool>(m->size()),std::vector<bool>(m->size()),std::vector<bool>(m->size())};
    const std::array<const std::vector<NA6PMatch>*,3> matches{&*a,&*b,&*c};
    for(size_t pass=0;pass<matches.size();++pass) for(const auto& match:*matches[pass]) {
      int vi=match.getIndexVT(), mi=match.getIndexMS();
      if(vi<0 || mi<0 || size_t(vi)>=v->size() || size_t(mi)>=m->size()) throw std::runtime_error("Invalid match index");
      if(seen[pass][mi]) throw std::runtime_error("Multiple matches for one MS track");
      seen[pass][mi]=true;
      const int hits=v->at(vi).getNHits();
      if(hits<3 || hits>5) throw std::runtime_error("Unexpected VT hit count");
      const bool correct=genuine(vl->at(vi)) && genuine(ml->at(mi)) && vl->at(vi)==ml->at(mi);
      state[pass][mi]=(correct ? 1 : 4)+(hits-3);
    }
    for(size_t mi=0;mi<m->size();++mi) {
      const auto& lbl=ml->at(mi);
      const int id=lbl.getTrackID();
      if(m->at(mi).getNHits()<6 || !genuine(lbl) || lbl.getSourceID()!=0 ||
         lbl.getEventID()!=header->getEventID() || id<0 || size_t(id)>=thermal.size() || !thermal[id]) continue;
      ++selected;
      for(int pass=0;pass<3;++pass) ++populations[pass][state[pass][mi]];
      for(int pass=0;pass<2;++pass) ++flow[pass][state[pass][mi]][state[pass+1][mi]];
      ++flow[2][state[2][mi]][state[0][mi]];
    }
    ++events;
  }
  if(vr.Next()||mr.Next()||r3.Next()||r4.Next()||r5.Next()) throw std::runtime_error("Mismatched event counts");
  std::ofstream out(path+"/transitions/flows.csv");
  if(!out) throw std::runtime_error("Cannot create flows.csv");
  out<<"from_threshold,to_threshold,from_state,to_state,count\n";
  for(int step=0;step<3;++step) for(int i=0;i<nStates;++i) for(int j=0;j<nStates;++j)
    if(flow[step][i][j]) out<<(step==0?3:step==1?4:5)<<","<<(step==0?4:step==1?5:3)<<","<<i<<","<<j<<","<<flow[step][i][j]<<"\n";
  std::ofstream stats(path+"/transitions/populations.csv");
  if(!stats) throw std::runtime_error("Cannot create populations.csv");
  stats<<"threshold,state,count\n";
  for(int pass=0;pass<3;++pass) for(int stateId=0;stateId<nStates;++stateId)
    stats<<(pass+3)<<","<<stateId<<","<<populations[pass][stateId]<<"\n";
  std::cout<<"Events "<<events<<"; genuine thermal MS tracks with at least 6 hits "<<selected<<"\n";
}
