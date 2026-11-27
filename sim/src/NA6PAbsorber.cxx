// NA6PCCopyright

#include "NA6PAbsorber.h"
#include "NA6PDetector.h"
#include "NA6PTGeoHelper.h"
#include "NA6PLayoutParam.h"

#include <TGeoVolume.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TColor.h>
#include <fairlogger/Logger.h>

void NA6PAbsorber::createMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  std::string nameM;
  nameM = addName("Iron");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 55.845, 26, 7.874);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kGray);
  }
  nameM = addName("BeO");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 2, 3.01);
    mixt->AddElement(new TGeoElement("Be", "Beryllium",  4, 9.012), 1);
    mixt->AddElement(new TGeoElement("O",  "Oxygen",     8, 16.00), 1);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM, "", kYellow + 2);
  }
  nameM = addName("Graphite");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 12.01, 6, 2.267);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kGray + 2);
  }
  nameM = addName("SyntGraphite");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 12.01, 6, 1.9);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kBlue - 5);
  }
  nameM = addName("Tungsten");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 183.84, 74, 19.3);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kBlack);
  }
  nameM = addName("Air");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 2, 0.001);
    mixt->AddElement(new TGeoElement("N", "Nitrogen", 7, 14.01), 0.78);
    mixt->AddElement(new TGeoElement("O", "Oxygen", 8, 16.00), 0.22);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM);
  }
  nameM = addName("Al2O3");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 2, 3.97);
    mixt->AddElement(new TGeoElement("Al", "Aluminium", 13, 26.98), 2);
    mixt->AddElement(new TGeoElement("O",  "Oxygen",     8, 16.00), 3);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM, "", kBlue - 7);
  }
  nameM = addName("Concrete");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 6, 2.3);
    mixt->AddElement(new TGeoElement("H",  "Hydrogen",  1,  1.01),  0.01);
    mixt->AddElement(new TGeoElement("O",  "Oxygen",    8, 16.00),  0.529);
    mixt->AddElement(new TGeoElement("Al", "Aluminium",13, 26.98), 0.034);
    mixt->AddElement(new TGeoElement("Si", "Silicon", 14, 28.09),  0.337);
    mixt->AddElement(new TGeoElement("Ca", "Calcium", 20, 40.08),  0.044);
    mixt->AddElement(new TGeoElement("Fe", "Iron",    26, 55.85),  0.046);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM, "", kBlue - 10);
  }
}

void NA6PAbsorber::createGeometry(TGeoVolume* world)
{
  const auto& param = NA6PLayoutParam::Instance();

  createMaterials();

  float zplace = param.posZStartAbsorber;
  if (param.nAbsorberSlices < param.nAbsorberSlicesWall) {
    LOGP(fatal, "Number of MuonWall slices {} exceeds total number of absorber slices {}", param.nAbsorberSlicesWall, param.nAbsorberSlices);
  }
  for (int isl = 0; isl < param.nAbsorberSlices; isl++) {
    if (isl == param.nAbsorberSlices - param.nAbsorberSlicesWall) { // special setting for last slice (muon wall)
      zplace = param.posMuonWall;
    }
    TGeoShape* slice = nullptr;
    auto slnm = fmt::format("AbsSl{}{}", isl, param.medAbsorber[isl]);
    if (param.dimYAbsorber[isl] > 0) {
      slice = new TGeoBBox((slnm + "SH").c_str(), param.dimXAbsorber[isl] / 2, param.dimYAbsorber[isl] / 2, param.thicknessAbsorber[isl] / 2);
    } else {
      slice = new TGeoTube((slnm + "SH").c_str(), 0.0, param.dimXAbsorber[isl] / 2, param.thicknessAbsorber[isl] / 2);
    }

    double zpos = zplace + param.thicknessAbsorber[isl] / 2;

    // Determine what exists
    bool hasAir = (param.radAir[isl] > 0);
    bool hasPlug = (!param.medAbsorberPlug[isl].empty() && param.radPlug[isl] > 0);

    // Determine hole radius (important!)
    double holeRadius = 0.0;
    if (hasAir) {
      holeRadius = param.radAir[isl];
    } else if (hasPlug) {
      holeRadius = param.radPlug[isl];
    }

    // --- 1. Subtract hole if needed
    if (holeRadius > 0) {
      auto hole = new TGeoTube(
        (slnm + "_H").c_str(),
        0.0,
        holeRadius,
        param.thicknessAbsorber[isl]);

      slice = new TGeoCompositeShape(
        (slnm + "_HS").c_str(),
        new TGeoSubtraction(slice, hole));
    }

    // --- 2. Create AIR (if requested)
    TGeoVolume* airVol = nullptr;

    if (hasAir) {
      auto* airShape = new TGeoTube(
        (slnm + "_AIR_SH").c_str(),
        0.0,
        param.radAir[isl],
        param.thicknessAbsorber[isl] / 2);

      airVol = new TGeoVolume(
        (slnm + "_AIR").c_str(),
        airShape,
        NA6PTGeoHelper::instance().getMedium("Air"));

      airVol->SetLineColor(kBlue - 10);

      world->AddNode(
        airVol,
        composeNonSensorVolID(isl + 100), // (you may want to fix IDs as discussed before)
        new TGeoTranslation(0.0, 0.0, zpos));
    }

    // --- 3. Create PLUG (independent of air)
    if (hasPlug) {

      auto* plugShape = new TGeoTube(
        (slnm + "_PLUG_SH").c_str(),
        param.radPlugInt[isl],
        param.radPlug[isl],
        param.thicknessAbsorber[isl] / 2);

      auto* plugVol = new TGeoVolume(
        (slnm + "_PLUG").c_str(),
        plugShape,
        NA6PTGeoHelper::instance().getMedium(addName(param.medAbsorberPlug[isl])));

      plugVol->SetLineColor(
        NA6PTGeoHelper::instance().getMediumColor(addName(param.medAbsorberPlug[isl])));

      if (hasAir) {
        // place inside air
        airVol->AddNode(
          plugVol,
          composeNonSensorVolID(isl + PlugCopyOffset),
          new TGeoTranslation(0.0, 0.0, 0.0));
      } else {
        // place directly in world
        world->AddNode(
          plugVol,
          composeNonSensorVolID(isl + PlugCopyOffset),
          new TGeoTranslation(0.0, 0.0, zpos));
      }
    }
    auto sliceVol = new TGeoVolume(slnm.c_str(), slice, NA6PTGeoHelper::instance().getMedium(addName(param.medAbsorber[isl])));
    sliceVol->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName(param.medAbsorber[isl])));
    world->AddNode(sliceVol, composeNonSensorVolID(isl), new TGeoTranslation(0.0, 0.0, zplace + param.thicknessAbsorber[isl] / 2));
    zplace += param.thicknessAbsorber[isl];
  }
}
