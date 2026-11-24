// NA6PCCopyright
#ifndef _MAGFIELD_H
#define _MAGFIELD_H

#include "MagneticFieldRegion.h"
#include <TGeoGlobalMagField.h>
#include <TVirtualMagField.h>

class MagneticField : public TVirtualMagField
{
 public:
  MagneticField() : TVirtualMagField("NA6PMagField") {}
  ~MagneticField() override;
  void Field(const Double_t* x, Double_t* B) override
  {
    getField(x, B);
  }

  auto isInitDone() const { return mInitDone; }
  void loadField();

  void setScale2Unit(float v) { mScale2Unit = v; }
  auto getScale2Unit() const { return mScale2Unit; }

  template <typename T = float>
  bool getField(const T* xyz, T* bxbybz) const
  {
    bxbybz[0] = bxbybz[1] = bxbybz[2] = 0;
    auto res = mDipoleVT.addField(xyz, bxbybz);
    res |= mDipoleMS.addField(xyz, bxbybz);
    if (res && mScale2Unit) {
      for (int i = 0; i < 3; i++) {
        bxbybz[i] *= mScale2Unit;
      }
    }
    return res;
  }

  void setAsGlobalField();

  MagneticFieldRegion& getMagFieldIP() { return mDipoleVT; }
  MagneticFieldRegion& getMagFieldMS() { return mDipoleMS; }

 private:
  MagneticFieldRegion mDipoleVT{};
  MagneticFieldRegion mDipoleMS{};
  float mScale2Unit = 10.f; // multiplier to scale to needed unit (e.g. T -> kGaus)
  bool mInitDone = false;

  ClassDefOverride(MagneticField, 1);
};

#endif
