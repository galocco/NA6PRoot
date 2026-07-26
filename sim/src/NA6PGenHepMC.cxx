// NA6PCCopyright

#include "NA6PGenHepMC.h"
#include "NA6PMCStack.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/Attribute.h"
#include "HepMC3/Units.h"
#include <TDatabasePDG.h>
#include <TMath.h>
#include <cmath>

namespace
{
constexpr double SpeedOfLightCmPerS = 2.99792458e10;

double momentumScaleToGeV(HepMC3::Units::MomentumUnit unit)
{
  switch (unit) {
    case HepMC3::Units::GEV:
      return 1.;
    case HepMC3::Units::MEV:
      return 1.e-3;
  }
  LOGP(fatal, "Unsupported HepMC3 momentum unit {}", static_cast<int>(unit));
  return 1.;
}

double lengthScaleToCm(HepMC3::Units::LengthUnit unit)
{
  switch (unit) {
    case HepMC3::Units::CM:
      return 1.;
    case HepMC3::Units::MM:
      return 0.1;
  }
  LOGP(fatal, "Unsupported HepMC3 length unit {}", static_cast<int>(unit));
  return 1.;
}
} // namespace

NA6PGenHepMC::NA6PGenHepMC(const std::string& name, const std::string& filname, bool storeDecayed) : NA6PGenerator(name), mFileName(filname), mStoreDecayedPrimaries(storeDecayed)
{
  LOGP(info, "create NA6PGenHepMC {} {}, decayed primaries will {}be stored", name, filname, storeDecayed ? "" : "not ");
}

long NA6PGenHepMC::canGenerateMaxEvents() const
{
  if (mNEvInTree < 0 || isInitDone()) {
    LOGP(fatal, "Initialization was not done");
  }
  return mNEvInTree;
}

void NA6PGenHepMC::init()
{
  LOGP(info, "NA6PGenHepMC init");
  if (isInitDone()) {
    return;
  }
  try {
    auto fl = TFile::Open(mFileName.c_str());
    assert(fl);
    TTree* tree = (TTree*)fl->Get(HEPTreeName.c_str());
    assert(tree);
    mNEvInTree = tree->GetEntries();
    delete tree;
    fl->Close();
  } catch (...) {
    LOGP(fatal, "Failed to extract the HEPMC tree {} from file {}", HEPTreeName, mFileName);
  }

  mHEPRootFileReader = std::make_unique<HepMC3::ReaderRootTree>(mFileName);
  if (mHEPRootFileReader->failed()) {
    LOGP(fatal, "No HepMC input file found {}", mFileName);
  }

  NA6PGenerator::init();
}

void NA6PGenHepMC::skipEvents(long nEvents)
{
  if (nEvents <= 0) {
    return;
  }
  if (!mHEPRootFileReader) {
    LOGP(fatal, "Cannot skip {} HepMC events before reader initialization", nEvents);
  }
  HepMC3::GenEvent evt;
  for (long i = 0; i < nEvents; ++i) {
    mHEPRootFileReader->read_event(evt);
    if (mHEPRootFileReader->failed()) {
      LOGP(fatal, "Failed to skip HepMC event {} out of {} from {}", i, nEvents, mFileName);
    }
    evt.clear();
  }
  mReadEvents += nEvents;
  LOGP(info, "Skipped {} HepMC events", nEvents);
}

void NA6PGenHepMC::generate()
{
  if (mReadEvents >= mNEvInTree) {
    LOGP(error, "Something is wrong: event {} beyond the available range [0:{}] is requested, wrapping to 0", mReadEvents, mNEvInTree - 1);
    mReadEvents = 0;
  }
  generatePrimaryVertex(); // will generate if not generated yet, otherwise, will set the origin from the MCHeader of the stack.
  auto mcHead = getStack()->getEventHeader();
  double xv = mcHead->getVX();
  double yv = mcHead->getVY();
  double zv = mcHead->getVZ();
  std::unordered_map<int, int> partID2StoreID;

  HepMC3::GenEvent evt;
  mHEPRootFileReader->read_event(evt); // Read one event
  if (mHEPRootFileReader->failed()) {
    LOGP(info, " End of file reached . Exit .\n");
    return;
  }

  const double momentumScale = momentumScaleToGeV(evt.momentum_unit());
  const double lengthScale = lengthScaleToCm(evt.length_unit());
  const double timeScale = lengthScale / SpeedOfLightCmPerS; // HepMC stores t in length/c units.

  // Global variables for heavy-ion running
  double b;
  int ncoll, npart, npart_proj, npart_targ;
  long int id1, id2;
  double e1, e2;
  int zproj, aproj, ztarg, atarg;
  float za_ratio_proj, za_ratio_targ;

// Centrality-related info
  auto attr_b_ptr = evt.attribute<HepMC3::DoubleAttribute>("b");
  if (attr_b_ptr) {
    b = attr_b_ptr->value();
    LOGP(info, "b = {} fm", attr_b_ptr->value());
  }
  auto attr_ncoll_ptr  = evt.attribute<HepMC3::IntAttribute>("nColl");
  if (attr_ncoll_ptr) {
    ncoll = attr_ncoll_ptr->value();
    LOGP(info, "Ncoll = {}", attr_ncoll_ptr->value());
  }
  auto attr_npart_ptr  = evt.attribute<HepMC3::IntAttribute>("nPartTot");
  if (attr_npart_ptr) {
    npart = attr_npart_ptr->value();
    LOGP(info, "Npart = {}", attr_npart_ptr->value());
  }
  auto attr_npart_proj_ptr = evt.attribute<HepMC3::IntAttribute>("nPartProj");
  if (attr_npart_proj_ptr) {
    npart_proj = attr_npart_proj_ptr->value();
    LOGP(info, "Npart_proj = {}", attr_npart_proj_ptr->value());
  }
  auto attr_npart_targ_ptr = evt.attribute<HepMC3::IntAttribute>("nPartTarg");
  if (attr_npart_targ_ptr) {
    npart_targ = attr_npart_targ_ptr->value();
    LOGP(info, "Npart_targ = {}", attr_npart_targ_ptr->value());
  }

// Proj/target species and energy (needs to be stored at each event because storing in GenRunInfo is not supported for HepMC3 v3.03.00 currently linked from O2
  auto attr_id1_ptr = evt.attribute<HepMC3::LongAttribute>("idbeam1");
  auto attr_id2_ptr = evt.attribute<HepMC3::LongAttribute>("idbeam2");
  auto attr_e1_ptr = evt.attribute<HepMC3::DoubleAttribute>("ebeam1");
  auto attr_e2_ptr = evt.attribute<HepMC3::DoubleAttribute>("ebeam2");
  if (attr_id1_ptr && attr_id2_ptr && attr_e1_ptr && attr_e2_ptr) {
    id1 = attr_id1_ptr->value();
    id2 = attr_id2_ptr->value();
    e1 = attr_e1_ptr->value();
    e2 = attr_e2_ptr->value();
    if (id1> 1000000000) {
      zproj = (id1 / 10000) % 1000; // Extract Z (atomic number)
      aproj = (id1 / 10) % 1000; // Extract A (mass number)
      za_ratio_proj = (float) zproj/aproj;
    }
    if (id2> 1000000000) {
      ztarg = (id2 / 10000) % 1000; // Extract Z (atomic number)
      atarg = (id2 / 10) % 1000; // Extract A (mass number)
      za_ratio_targ = (float) ztarg/atarg;
    }
  }

  int nparticlesgood = 0;
  int pdgCode, status;
  float phi, pt, pZ, en, sn, cs;
  float vx, vy, vz, tof; // NA6PMCStack expects positions in cm and tof in s.
  int tobetracked = 0;

  std::vector<int> partons = {1, 2, 3, 4, 5, 6, 21}; // to remove quarks and gluons
  int cnt = 0;
  // Loop on particles
  for (const auto p : evt.particles()) {
    cnt++;
    auto isParton = [](int pabs) { return (pabs >= 1 & pabs <= 6) || pabs == 21; };
    auto isDiquark = [](int pabs) { return pabs > 1000 && pabs < 10000 && ((pabs / 10) % 10) == 0; };
    auto isColorState = [](int pabs) { return pabs > 1e6; };

    // Upstream and downstream vertices
    auto prod_vtx = p->production_vertex();
    auto end_vtx = p->end_vertex();
    auto pidabs = std::abs(p->pid());

    if (prod_vtx) {
      if (isParton(pidabs) ||  // Remove quarks and gluons
          isDiquark(pidabs) || // Remove diquarks
          isColorState(pidabs) // Remove particles like octet states for charmonia and similar
      ) {
        continue;
      }
      // Intermediate particles have both a production vertex and an end vertex
      if (end_vtx) {
        if (!mStoreDecayedPrimaries) {
          continue;
        }
        tobetracked = 0; // decayed particles should not be tracked by GEANT4
      } else {
        tobetracked = 1;
      }
      partID2StoreID[p->id()] = nparticlesgood++;

      int mothertobestored = -1;
      for (const auto& mother : prod_vtx->particles_in()) {
        auto mothEntry = partID2StoreID.find(mother->id());
        if (mothEntry != partID2StoreID.end()) {
          mothertobestored = mothEntry->second;
          break;
        }
      }

      phi = p->momentum().phi();
      pt = p->momentum().perp() * momentumScale;
      pZ = p->momentum().pz() * momentumScale;
      en = p->momentum().e() * momentumScale;
      status = p->status();
      sn = std::sin(phi);
      cs = std::cos(phi);
      pdgCode = p->pid();
      HepMC3::FourVector pos = prod_vtx->position();
      vx = pos.x() * lengthScale + xv; // shift to account for generated position of the primary vertex
      vy = pos.y() * lengthScale + yv;
      vz = pos.z() * lengthScale + zv;
      tof = pos.t() * timeScale;
      int dummy = 0;

      getStack()->PushTrack(tobetracked, mothertobestored, pdgCode, pt * cs, pt * sn, pZ, en, vx, vy, vz, tof, 0., 0., 0., TMCProcess::kPPrimary, dummy, 1., status);
      getStack()->setNColl(ncoll);
      getStack()->setNPart(npart);
      getStack()->setImpPar(b);
    } else {
      LOGP(info, "particle ID# {} with no production vertex\n", p->id());
      continue;
    }
  }

// Tracking spectator neutrons and protons
  if(mTrackSpectators && attr_npart_proj_ptr && id1 > 1000000000 && id2 > 1000000000){
    int nSpecProtons = (aproj - npart_proj) * za_ratio_proj;
    int nSpecNeutrons = (aproj - npart_proj) - nSpecProtons;
    int dummy = 0;
    LOGP(info, "Tracking {} spectator protons",nSpecProtons);
    for (int i=0; i<nSpecProtons; i++){
      getStack()->PushTrack(1, -1, 2212, 0., 0., std::sqrt(e1*e1 - 0.9383*0.9383), e1, xv, yv, zv, 0., 0., 0., 0., TMCProcess::kPPrimary, dummy, 1., status);
      nparticlesgood++;
    }
    LOGP(info, "Tracking {} spectator neutrons",nSpecNeutrons);
    for (int i=0; i<nSpecNeutrons; i++){
      getStack()->PushTrack(1, -1, 2112, 0., 0., std::sqrt(e1*e1 - 0.9396*0.9396), e1, xv, yv, zv, 0., 0., 0., 0., TMCProcess::kPPrimary, dummy, 1., status);
      nparticlesgood++;
    }
  }

  auto info = fmt::format("{}_x{}", getName(), nparticlesgood);
  mcHead->getGenHeaders().emplace_back(nparticlesgood, 0, mcHead->getNPrimaries(), 0, info);
  mcHead->incNPrimaries(nparticlesgood);
  ++mReadEvents;
}
