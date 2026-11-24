// NA6PCCopyright

#include "MagneticField.h"
#include "NA6PLayoutParam.h"
#include <fairlogger/Logger.h>
#include <TVirtualMagField.h>
#include <TSystem.h>
#include "NA6PLayoutParam.h"

MagneticField::~MagneticField()
{
  if (TGeoGlobalMagField::Instance()->GetField() == this) {
    LOGP(error, "This field {} was set as global field, cannot be deleted explicitly", GetName());
  }
}

void MagneticField::loadField()
{
  const auto& param = NA6PLayoutParam::Instance();
  auto loadOne = [](auto& dipole, const std::string& path, const float pos[3]) {
    std::string fname = gSystem->ExpandPathName(path.c_str());
    if (fname.find(".inp") != std::string::npos) {
      // FLUKA format
      dipole.loadFlukaField(fname);
    } else {
      // Alternative format: OPERA3D or Ansys
      dipole.loadOpera3DField(fname);
    }
    dipole.setRefPosition(pos[0], pos[1], pos[2]);
  };
  loadOne(mDipoleVT, param.fieldMapDipVT, param.posDipIP);
  loadOne(mDipoleMS, param.fieldMapDipMS, param.posDipMS);
  mInitDone = true;
}

void MagneticField::setAsGlobalField()
{
  if (!isInitDone()) {
    LOGP(fatal, "The magnetic field was not initialized yet");
  }
  TGeoGlobalMagField::Instance()->SetField(this);
}