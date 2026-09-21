// NA6PCCopyright

#ifndef NA6P_MWPC_PARAMS_H_
#define NA6P_MWPC_PARAMS_H_

#include "ConfigurableParam.h"
#include "ConfigurableParamHelper.h"

struct NA6PMWPCParam : public na6p::conf::ConfigurableParamHelper<NA6PMWPCParam> {
  static constexpr int MaxStations = 7;

  // Detailed chamber dimensions, cm, in the chamber-local mechanical convention
  // inherited from the prototype drawing/model:
  //   local x = 53.2 cm short body side
  //   local y = 68.56 cm long body side
  // The complete chamber is rotated around z when installed in an MS station so
  // that detector X is horizontal/along the long side (wires) and detector Y is
  // vertical/across the wires. With the default +90 deg rotation the standard
  // sensitive gas footprint in detector coordinates is therefore 66.56 x 51.2 cm.
  float bodyX = 53.2f;
  float bodyY = 68.56f;
  float chamberRotationZDeg = 90.0f;
  float innerFrameWidth = 1.0f;
  float honeycombEdgeWall = 0.1f;
  float readoutEndMargin = 0.78f;
  float electronicsPlaneOverhang = 6.8f;
  float dividerPlaneLengthMargin = 0.2f;
  float dividerPlaneOverhang = 4.0f;
  float dividerPlaneThickness = 0.2f;

  // Layer thicknesses, cm
  float outerSkinFront = 0.015f;
  float honeycombFront = 0.6f;
  float readoutFR4 = 0.08f;
  float readoutCopper = 0.0017f;
  float gasGap = 0.6f;
  float coverFR4 = 0.015f;
  float coverCopper = 0.0017f;
  float honeycombBack = 0.6f;
  float outerSkinBack = 0.015f;

  bool coverCopperOnGasFace = false;
  bool readoutCopperOnElectronicsExtension = true;
  bool includeDividerPlane = true;

  // Approximate required detector working rectangles (global detector X,Y), cm,
  // from the updated MS layout study. They define the present finite station
  // footprints used when comparing regular staggered grids.
  float stationWorkingAreaX[MaxStations] = {220.f, 230.f, 310.f, 320.f, 440.f, 500.f, 0.f};
  float stationWorkingAreaY[MaxStations] = {220.f, 240.f, 310.f, 320.f, 410.f, 440.f, 0.f};

  // Current provisional regular grids in detector coordinates. NX counts
  // chambers horizontally (detector X), NY vertically (detector Y):
  //   MS0 4x8, MS1 4x5, MS2 6x7, MS3 6x7, MS4 8x9, MS5 9x10.
  // Ox=Oy=3 cm and adjacent A/B/C/D z spacing is 4 cm. These grids retain the
  // requested finite working areas while adding one horizontal column at MS2-MS5
  // relative to the minimum-covering baseline for the current acceptance study.
  // They are a working point, not a frozen final design; the arrays remain
  // configurable for further acceptance studies.
  int stationGridNX[MaxStations] = {4, 4, 6, 6, 8, 9, 0};
  int stationGridNY[MaxStations] = {8, 5, 7, 7, 9, 10, 0};

  // Strip readout definition shared by digitization and reconstruction.
  // A strip count of zero lets the sensitive chamber extent define the range.
  float stripPitch[2] = {0.1f, 0.1f};
  // Station-dependent pitches chosen for sigma = pitch/sqrt(12), in cm.
  // MS0-MS3: 0.1 mm (X), 0.5 mm (Y); MS4-MS5: 5 mm in both coordinates.
  float stationStripPitch[MaxStations][2] = {
    {0.034641f, 0.173205f}, {0.034641f, 0.173205f},
    {0.034641f, 0.173205f}, {0.034641f, 0.173205f},
    {1.73205f, 1.73205f}, {1.73205f, 1.73205f}, {0.f, 0.f}};
  float stripAngleDeg[2] = {0.f, 90.f};
  int stripCount[2] = {0, 0};
  float maxClusterChargeAsymmetry = 0.5f;

  int getNModules(int nStations = MaxStations) const
  {
    if (nStations < 0) {
      nStations = 0;
    } else if (nStations > MaxStations) {
      nStations = MaxStations;
    }
    int nModules = 0;
    for (int i = 0; i < nStations; ++i) {
      nModules += stationGridNX[i] * stationGridNY[i];
    }
    return nModules;
  }

  // Active-gas overlap, not mechanical-envelope overlap. Chamber centre pitches
  // are globalGasX-activeOverlapX and globalGasY-activeOverlapY.
  float activeOverlapX = 3.0f;
  float activeOverlapY = 3.0f;

  // Four-level checkerboard longitudinal stagger. For a chamber at (row,col):
  //   q = 2*(row%2) + (col%2) = A,B,C,D -> 0,1,2,3
  //   zGas - zStation = (q - 1.5)*staggerZStep.
  // With the default 4 cm step the gas centres are at -6,-2,+2,+6 cm.
  float staggerZStep = 4.0f;

  // MS0 uses the requested shorter vertical active dimension. Because the
  // detailed chamber builder keeps the prototype-local axes and the installed
  // chamber is rotated by 90 deg, detector Y corresponds to chamber-local x.
  bool useNarrowMS0 = true;
  float ms0GasY = 30.72f;

  // Binary Ar/CO2 gas mixture. The CO2 fraction is 1 - argonMoleFraction.
  float argonMoleFraction = 0.7f;
  float gasTemperatureK = 293.15f;
  float gasPressurePa = 101325.f;

  // Base material names. The parent detector module adds its module suffix.
  std::string medFR4 = "FR4";
  std::string medCopper = "Copper";
  std::string medHoneycomb = "Honeycomb";
  std::string medGas = "ArCO2";

  NA6PParamDef(NA6PMWPCParam, "mwpc");
};

#endif
