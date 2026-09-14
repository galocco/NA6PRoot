// NA6PCCopyright

#ifndef NA6P_MWPC_PARAMS_H_
#define NA6P_MWPC_PARAMS_H_

#include "ConfigurableParam.h"
#include "ConfigurableParamHelper.h"

struct NA6PMWPCParam : public na6p::conf::ConfigurableParamHelper<NA6PMWPCParam> {
  // Chamber dimensions, cm
  float bodyX = 53.2f;
  float bodyY = 68.56f;
  float innerFrameWidth = 1.0f;
  float honeycombEdgeWall = 0.1f;
  float readoutEndMargin = 0.78f;
  float blueOverhang = 6.8f;
  float greenLengthMargin = 0.2f;
  float greenOverhang = 4.0f;
  float greenThickness = 0.2f;

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
  bool readoutCopperOnBlueExtension = true;
  bool includeGreenExternalStrip = true;

  // Gas conditions
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
