// NA6PCCopyright

/// \file Propagator
/// \brief Singleton class for track propagation routines
/// \author ruben.shahoyan@cern.ch

#ifndef NA6P_PROPAGATOR_
#define NA6P_PROPAGATOR_

#include "NA6PTrackParCov.h"
#include "MagneticField.h"
#include <array>
#include <string>

struct MatBudget {
  float meanRho = 0.f;  ///< mean density, g/cm^3
  float meanX2X0 = 0.f; ///< fraction of radiaton lenght
  float meanA = 0.f;    ///< mean A: sum(x_i*A_i)/sum(x_i)
  float meanZ = 0.f;    ///< mean Z: sum(x_i*Z_i)/sum(x_i)
  float meanL = 0.f;    ///< length: sum(x_i)
  float meanZ2A = 0.f;  ///< Z/A mean: sum(x_i*Z_i/A_i)/sum(x_i)
  int nBoundaries = 0;  ///< number of boundary crosses
  ClassDefNV(MatBudget, 1);
};

class Propagator
{
 public:
  static Propagator* Instance(bool uninitialized = false)
  {
    static Propagator instance(uninitialized);
    return &instance;
  }
  static bool loadField();
  static bool loadGeometry(const std::string& path = "geometry.root", const std::string geomName = "NA6P");

  enum class MatCorrType : int8_t {
    USEMatCorrNONE, // flag to not use material corrections
    USEMatCorrTGeo, // flag to use TGeo for material queries
    //    USEMatCorrLUT // TODO
  }; // flag to use LUT for material queries (user must provide a pointer

  struct PropOpt {
    NA6PTrackPar* linRef = nullptr; // optional Kalman linearization reference
    float maxStep{2.0f};
    MatCorrType matCorr{MatCorrType::USEMatCorrTGeo};
    bool byOnly{false};
    bool fullFieldJacobian{true};
    bool fieldGradientJacobian{false};
    ClassDefNV(PropOpt, 3);
  };

  bool propagateToZ(NA6PTrackParCov& track, float z) const;
  bool propagateToZ(NA6PTrackPar& track, float z) const;
  bool propagateToZ(NA6PTrackParCov& track, float z, const PropOpt& opt) const;
  bool propagateToZ(NA6PTrackPar& track, float z, const PropOpt& opt) const;

  template <class T>
  bool propagatePCAToLine(T& track, float x, float y, float tolerance, const PropOpt& opt) const;
  template <class T>
  bool propagatePCAToLine(T& track, float x, float y, float tolerance) const;

  Propagator(Propagator const&) = delete;
  Propagator(Propagator&&) = delete;
  Propagator& operator=(Propagator const&) = delete;
  Propagator& operator=(Propagator&&) = delete;

  template <typename T = float>
  MatBudget getMeanMaterial(const T* start, const T* end) const;
  template <typename T = float>
  MatBudget getMeanMaterial(const std::array<T, 3>& start, const std::array<T, 3>& end) const
  {
    return getMeanMaterial(start.data(), end.data());
  }

  void updateField();
  bool hasMagFieldSet() const { return mField != nullptr; }
  bool hasGeometryLoaded() const;

  template <typename T = float>
  void getFieldXYZ(const T* xyz, T* bxyz) const
  {
    mField->getField(xyz, bxyz);
  }

  template <typename T = float>
  void getFieldXYZ(const std::array<T, 3>& xyz, std::array<T, 3>& bxyz) const
  {
    getFieldXYZ(xyz.data(), bxyz.data());
  }

  template <typename T = float>
  std::array<T, 3> getFieldXYZ(const std::array<T, 3>& xyz) const
  {
    std::array<T, 3> bxyz;
    getFieldXYZ(xyz.data(), bxyz.data());
    return bxyz;
  }

  template <typename T = float>
  std::array<T, 3> getFieldXYZ(const T* xyz) const
  {
    std::array<T, 3> bxyz;
    getFieldXYZ(xyz, bxyz.data());
    return bxyz;
  }

  template <typename T = float>
  float getBy(const T* xyz) const
  {
    T b[3];
    mField->getField(xyz, b);
    return b[1];
  }

  template <typename T = float>
  float getBy(const std::array<T, 3>& xyz) const
  {
    return getBy(xyz.data());
  }

  int getNDOFTrack() const { return mNDOF; }

 private:
  static constexpr float Epsilon = 0.00001; // precision of propagation to Z

  Propagator(bool uninitialized = false);
  ~Propagator() = default;

  template <typename T = float>
  MatBudget getMeanMaterialBudgetFromGeom(const T* start, const T* end) const;
  template <typename T = float>
  MatBudget getMeanMaterialBudgetFromGeom(const std::array<T, 3>& start, const std::array<T, 3>& end) const
  {
    return getMeanMaterialBudgetFromGeom(start.data(), end.data());
  }

  MagneticField* mField = nullptr; ///< External or own nominal field map
  int mNDOF = 5;

  ClassDefNV(Propagator, 0);
};

inline bool Propagator::propagateToZ(NA6PTrackParCov& track, float z) const
{
  return propagateToZ(track, z, PropOpt{});
}

inline bool Propagator::propagateToZ(NA6PTrackPar& track, float z) const
{
  return propagateToZ(track, z, PropOpt{});
}

template <typename T>
inline MatBudget Propagator::getMeanMaterial(const T* start, const T* end) const
{
  // at the moment TGeo only, in future LUT should be added
  return getMeanMaterialBudgetFromGeom(start, end);
}

template <class T>
inline bool Propagator::propagatePCAToLine(T& track, float x, float y, float tolerance) const
{
  return propagatePCAToLine(track, x, y, tolerance, PropOpt{});
}

#endif
