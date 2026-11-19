#if !defined(__CINT__) || defined(__MAKECINT__)
#include "NA6PGenParam.h"
#include "NA6PBeamParam.h"
#include "NA6PGenCocktail.h"
#include <TMath.h>
#include <fairlogger/Logger.h>
#endif

bool getParams(float& T_piM, float& y0_piM, float& ysig_piM, float& dNdY_piM,
               float& T_piP, float& y0_piP, float& ysig_piP, float& dNdY_piP,
               float& T_KM, float& y0_KM, float& ysig_KM, float& dNdY_KM,
               float& T_KP, float& y0_KP, float& ysig_KP, float& dNdY_KP,
               float& T_Prot, float& y0_Prot, float& ysig_Prot, float& dNdY_Prot,
               float& yboxhalfw_Jpsi, float& ysig_Jpsi, float& pt1_Jpsi, float& pt2_Jpsi, float& pt3_Jpsi,
               float& sigySG, float& TSG_Phi, float& TSG_Omega, float& MotherMass_Phi, float& MotherMass_Omega)
{
  float en = NA6PBeamParam::Instance().energyPerNucleon;

  // macro to generate a fixed number of signals (J/Psi, phi, omega2body) in a kinematic range
  //
  // JPsi y and pT parametrization are based on PYTHIA8 (Enrico, NA60+ meeting 7/10/2025)
  // Phi and omega2body parametrizations are taken from fast sim (available only at 40 and 150 GeV)
  // Background parametrization is available only for 40 and 150 Gev and it can be switched on/off through the argument Bg

  if (std::abs(en - 40.) < 3) {
    LOGP(info, "Setting hadronic event parameters for PbPb at 40 GeV/c");
    // pions and Kaons from  NA49 nucl-ex/0205002
    T_piM = T_piP = 0.169;
    T_KM = 0.232;
    T_KP = 0.226;
    T_Prot = 0.25;

    y0_piM = y0_piP = 0.666; // pions 2-gaussian poles half distance
    y0_KM = 0.694;           // K- 2-gaussian poles half distance
    y0_KP = 0.694;           // K+ 2-gaussian poles half distance
    y0_Prot = 0.907;         // ???

    ysig_piM = ysig_piP = 0.872;
    ysig_KM = 0.635;
    ysig_KP = 0.725;
    ysig_Prot = 0.798;

    dNdY_piM = 106.1;
    dNdY_piP = 96.6;
    dNdY_KM = 7.58;
    dNdY_KP = 20.1;
    dNdY_Prot = 39.9;

  } else if (std::abs(en - 158.) < 3) {
    LOGP(info, "Setting hadronic event parameters for PbPb at 158 GeV/c");

    T_piM = T_piP = 0.18;
    T_KM = 0.226;
    T_KP = 0.232;
    T_Prot = 0.31;

    y0_piM = y0_piP = 0.72; // pions 2-gaussian poles half distance
    y0_KM = 0.727;          // K- 2-gaussian poles half distance
    y0_KP = 0.839;          // K+ 2-gaussian poles half distance
    y0_Prot = 39.8;         // ???

    ysig_piM = ysig_piP = 1.18;
    ysig_KM = 0.81;
    ysig_KP = 0.88;
    ysig_Prot = 8.07;

    dNdY_piM = 175.4;
    dNdY_piP = 170.1;
    dNdY_KM = 16.8;
    dNdY_KP = 29.6;
    dNdY_Prot = 27.;


  } else {
    LOGP(fatal, "No parameters defined for beam energy {}", en);
  }
  return true;
}

NA6PGenerator* addMuonBgEventGenerator(NA6PGenCocktail* genCocktail, int nMuons, float ptMin, float ptMax, float yMin, float yMax, bool Bg)
{
  float T_piM = 0., y0_piM = 0., ysig_piM = 0., dNdY_piM = 0.;
  float T_piP = 0., y0_piP = 0., ysig_piP = 0., dNdY_piP = 0.;
  float T_KM = 0., y0_KM = 0., ysig_KM = 0., dNdY_KM = 0.;
  float T_KP = 0., y0_KP = 0., ysig_KP = 0., dNdY_KP = 0.;
  float T_Prot = 0., y0_Prot = 0., ysig_Prot = 0., dNdY_Prot = 0.;
  float yboxhalfw_Jpsi = 0, ysig_Jpsi = 0., pt1_Jpsi = 0., pt2_Jpsi = 0., pt3_Jpsi = 0.;
  float sigySG = 0., TSG_Phi = 0., TSG_Omega = 0., MotherMass_Phi = 0, MotherMass_Omega = 0;

  const float ycm = NA6PBeamParam::Instance().getYCM();

  std::string dndptBgFun = "x*exp(-sqrt(x*x+[0]*[0])/[1])";
  std::string dndyBgFun = "exp(-0.5*pow((x-[0]-[1])/[2],2))+exp(-0.5*pow((x-[0]+[1])/[2],2))";


  getParams(T_piM, y0_piM, ysig_piM, dNdY_piM,
            T_piP, y0_piP, ysig_piP, dNdY_piP,
            T_KM, y0_KM, ysig_KM, dNdY_KM,
            T_KP, y0_KP, ysig_KP, dNdY_KP,
            T_Prot, y0_Prot, ysig_Prot, dNdY_Prot,
            yboxhalfw_Jpsi, ysig_Jpsi, pt1_Jpsi, pt2_Jpsi, pt3_Jpsi,
            sigySG, TSG_Phi, TSG_Omega, MotherMass_Phi, MotherMass_Omega);

  bool IsPoisson = true; 	

  if (Bg) {
    auto genpiM = new NA6PGenParam("pi-", -211, 0, dndptBgFun, dndyBgFun, ptMin, ptMax, yMin, yMax, true, true, IsPoisson);
    genpiM->setParametersTrans({0.1396f, T_piM});
    genpiM->setParametersLong({ycm, y0_piM, ysig_piM});
    genpiM->setdNdY(dNdY_piM);
    genCocktail->addGenerator(genpiM);

    auto genpiP = new NA6PGenParam("pi+", 211, 0, dndptBgFun, dndyBgFun, ptMin, ptMax, yMin, yMax, true, true, IsPoisson);
    genpiP->setParametersTrans({0.1396f, T_piP});
    genpiP->setParametersLong({ycm, y0_piP, ysig_piP});
    genpiP->setdNdY(dNdY_piP);
    genCocktail->addGenerator(genpiP);

    auto genKM = new NA6PGenParam("K-", -321, 0, dndptBgFun, dndyBgFun, ptMin, ptMax, yMin, yMax, true, true, IsPoisson);
    genKM->setParametersTrans({0.4937f, T_KM});
    genKM->setParametersLong({ycm, y0_KM, ysig_KM});
    genKM->setdNdY(dNdY_KM);
    genCocktail->addGenerator(genKM);

    auto genKP = new NA6PGenParam("K+", 321, 0, dndptBgFun, dndyBgFun, ptMin, ptMax, yMin, yMax, true, true, IsPoisson);
    genKP->setParametersTrans({0.4937f, T_KP});
    genKP->setParametersLong({ycm, y0_KP, ysig_KP});
    genKP->setdNdY(dNdY_KP);
    genCocktail->addGenerator(genKP);

    auto genProt = new NA6PGenParam("Proton", 2212, 0, dndptBgFun, dndyBgFun, ptMin, ptMax, yMin, yMax, true, true, IsPoisson);
    genProt->setParametersTrans({0.9383f, T_Prot});
    genProt->setParametersLong({ycm, y0_Prot, ysig_Prot});
    genProt->setdNdY(dNdY_Prot);
    genCocktail->addGenerator(genProt);
  }

  std::string dndptFun = "1"; 
  std::string dndyFun = "1";        

  const int nPlus = (nMuons + 1) / 2;
  const int nMinus = nMuons / 2;
  if (nPlus > 0) {
    auto genMuP = new NA6PGenParam("mu+", -13, nPlus, dndptFun, dndyFun, ptMin, ptMax, yMin, yMax, true, true, false);
    genCocktail->addGenerator(genMuP);
  }
  if (nMinus > 0) {
    auto genMuM = new NA6PGenParam("mu-", 13, nMinus, dndptFun, dndyFun, ptMin, ptMax, yMin, yMax, true, true, false);
    genCocktail->addGenerator(genMuM);
  }
  

  return genCocktail;
}

NA6PGenerator* genMuonBgEvent(float ptMin = 0., float ptMax = 5., float yMin = 0, float yMax = 6, bool Bg = true, int nMuons = 2)
{
  if (nMuons < 1 || ptMin < 0.f || ptMax <= ptMin || yMax <= yMin) {
    LOGP(error, "Invalid muon generator parameters: nMuons={}, pT=[{}, {}], y=[{}, {}]", nMuons, ptMin, ptMax, yMin, yMax);
    return nullptr;
  }
  NA6PGenCocktail* genCockt = new NA6PGenCocktail("cocktail");
  return addMuonBgEventGenerator(genCockt, nMuons, ptMin, ptMax, yMin, yMax, Bg);
}
