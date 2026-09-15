// NA6PCCopyright

#ifndef NA6P_MWPC_PARAMS_H_
#define NA6P_MWPC_PARAMS_H_

#include "ConfigurableParam.h"
#include "ConfigurableParamHelper.h"

struct NA6PMWPCParam : public na6p::conf::ConfigurableParamHelper<NA6PMWPCParam> {
  static constexpr int MaxStations = 7;

  // Chamber dimensions, cm.  These are the standard MS1--MS5 chamber values.
  float bodyX = 53.2f;
  float bodyY = 68.56f;
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

  // Preliminary full Muon-Spectrometer layout used for the stagger study.
  // Rows/columns are chamber counts in y/x respectively.  The defaults are
  // the present balanced 308-chamber working point:
  //   MS0 10x4, MS1 5x4, MS2 8x5, MS3 8x5, MS4 12x7, MS5 12x7.
  int stationGridNX[MaxStations] = {10, 5, 8, 8, 12, 12, 0};
  int stationGridNY[MaxStations] = {4, 4, 5, 5, 7, 7, 0};

  // Active-gas overlap, not mechanical-envelope overlap.  Chamber centre
  // pitches are gasX-activeOverlapX and gasY-activeOverlapY.
  float activeOverlapX = 3.0f;
  float activeOverlapY = 3.0f;

  // Four-level checkerboard longitudinal stagger.  For a chamber at (row,col):
  //   q = 2*(row%2) + (col%2) = A,B,C,D -> 0,1,2,3
  //   zGas - zStation = (q - 1.5)*staggerZStep.
  // With the default 4 cm step the gas centres are at -6,-2,+2,+6 cm.
  float staggerZStep = 4.0f;

  // MS0 uses the current first-pass narrow chamber form factor.  The detailed
  // mechanics are obtained by replacing bodyX while preserving the same frame
  // width and longitudinal stack, so ms0GasX is exactly the sensitive gas width.
  bool useNarrowMS0 = true;
  float ms0GasX = 25.6f;

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
