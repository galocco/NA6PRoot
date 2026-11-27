// NA6PCCopyright
#include "NA6PTrackParCov.h"
#include <TMath.h>
#include <fmt/format.h>
#include <fairlogger/Logger.h>

// ----------------------- ctor / init -----------------------
NA6PTrackParCov::NA6PTrackParCov(const float* xyz, const float* pxyz, int sign, float errLoose)
  : NA6PTrackPar(xyz, pxyz, sign)
{
  // same “reasonable defaults” as your original
  setCov(CovArray{1e-6, 0., 1e-6, 0., 0., 1e-6, 0., 0., 0., 1e-6, 0., 0., 0., 0., 1e-2});
  if (errLoose >= MaxErrSelRescale) {
    resetCovariance(errLoose);
  }
}
NA6PTrackParCov::NA6PTrackParCov(const std::array<float, 3>& xyz, const std::array<float, 3>& pxyz, int sign, float errLoose)
  : NA6PTrackPar(xyz, pxyz, sign)
{
  // same “reasonable defaults” as your original
  setCov(CovArray{1e-6, 0., 1e-6, 0., 0., 1e-6, 0., 0., 0., 1e-6, 0., 0., 0., 0., 1e-2});
  if (errLoose >= MaxErrSelRescale) {
    resetCovariance(errLoose);
  }
}

void NA6PTrackParCov::init(const float* xyz, const float* pxyz, int sign, float errLoose)
{
  initParam(xyz, pxyz, sign);
  setCov(CovArray{1e-6, 0., 1e-6, 0., 0., 1e-6, 0., 0., 0., 1e-6, 0., 0., 0., 0., 1e-2});
  if (errLoose >= MaxErrSelRescale) {
    resetCovariance(errLoose);
  }
}

std::string NA6PTrackParCov::asString() const
{
  return fmt::format("{} Cov=[{:+.3e}, {:+.3e},{:+.3e}, {:+.3e},{:+.3e},{:+.3e}, {:+.3e},{:+.3e},{:+.3e},{:+.3e}, {:+.3e},{:+.3e},{:+.3e},{:+.3e},{:+.3e}]",
                     NA6PTrackPar::asString(),
                     mC[0], mC[1], mC[2], mC[3], mC[4], mC[5], mC[6], mC[7], mC[8], mC[9], mC[10], mC[11], mC[12], mC[13], mC[14]);
}

// ------------------------------------------------------------------
// propagateToZ: state + Jacobian + covariance transport
// State: {x, y, s=sin(psi), ty, q/pxz}
// kappa = K * q,  where K = kB2C*By
// s1 = s0 + kappa*dz
// dx = dz*(s0+s1)/(c0+c1)
// dy = ty*dz*(dpsi/ds),  dpsi = asin(s1)-asin(s0), ds = s1-s0
// ------------------------------------------------------------------
bool NA6PTrackParCov::propagateToZ(float z, float by)
{
  NA6PTrackPar linRef{*this};
  return propagateToZ(z, by, linRef);
}

// propagate with linearization wrt provided reference rather than the track itself
bool NA6PTrackParCov::propagateToZ(float z, float by, NA6PTrackPar& linRef0)
{
#ifdef _SAVE_TRACK_FOR_DEBUG_
  auto sav = *this;
#endif
  const float dz = z - mZ;
  if (std::abs(dz) < 1e-6f) {
    setZ(z);
    linRef0.setZ(z);
    return true;
  }
  // propagate reference track
  NA6PTrackPar linRef1 = linRef0;
  if (!linRef1.propagateParamToZ(z, by)) {
    return false;
  }
  // The state stores q/pXZ for charged tracks, hence the charge is already
  // included in the fifth parameter. For neutral tracks +1/pXZ is only a
  // momentum placeholder and must not generate magnetic bending.
  const float K = linRef0.getCharge() == 0 ? 0.f : kB2C * by;

  prec_t snpRef0 = linRef0.getTx(), cspRef0 = linRef0.getCosPsi();
  prec_t snpRef1 = linRef1.getTx(), cspRef1 = linRef1.getCosPsi();
  if (cspRef0 < kTinyF || cspRef1 < kTinyF) {
    return false;
  }
  prec_t cc = cspRef0 + cspRef1;
  if (cc < kTinyF) {
    return false;
  }
  prec_t cspRef0Inv = 1. / cspRef0, cspRef1Inv = 1. / cspRef1, ccInv = 1. / cc, dx2dz = (snpRef0 + snpRef1) * ccInv;
  prec_t dzccInv = dz * ccInv, hh = dzccInv * cspRef1Inv * (1 + cspRef0 * cspRef1 + snpRef0 * snpRef1), jj = dz * (dx2dz - snpRef1 * cspRef1Inv);
  // --- Jacobian for covariance transport for lin ref---
  prec_t f02 = hh * cspRef0Inv;
  prec_t f04 = hh * dzccInv * K;
  prec_t f24 = dz * K;
  prec_t f12 = linRef0.getTy() * (f02 * snpRef1 + jj);
  prec_t f13 = dz * (cspRef1 + snpRef1 * dx2dz); // dS
  prec_t f14 = linRef0.getTy() * (f04 * snpRef1 + jj * f24);

  // difference between the current and reference state
  float diff[5];
  for (int i = 0; i < 5; i++) {
    diff[i] = getParam(i) - linRef0.getParam(i);
  }
  float snpUpd = snpRef1 + diff[kTx] + f24 * diff[kQ2Pxz];
  if (std::abs(snpUpd) > kAlmost1F) {
    return false;
  }
  setZ(z);
  setX(linRef1.getX() + diff[kX] + f02 * diff[kTx] + f04 * diff[kQ2Pxz]);
  setY(linRef1.getY() + diff[kY] + f12 * diff[kTx] + f13 * diff[kTy] + f14 * diff[kQ2Pxz]);
  setTx(snpUpd);
  setTy(linRef1.getTy() + diff[kTy]);
  setQ2Pxz(linRef1.getQ2Pxz() + diff[kQ2Pxz]);
  linRef0 = linRef1; // update reference track after transporting this track, as in O2

  transportCovariance(f02, f04, f12, f13, f14, f24);
  return true;
}

bool NA6PTrackParCov::propagateToZ(float z, const float* bxyz)
{
  NA6PTrackPar linRef{*this};
  return propagateToZ(z, bxyz, linRef);
}

bool NA6PTrackParCov::propagateToZ(float z, const float* bxyz, NA6PTrackPar& linRef0)
{
  //----------------------------------------------------------------
  // Extrapolate this track to the plane z in the field bxyz. Cov matrix is trasported with by only. Linearization wrt externally provided linRef
  //----------------------------------------------------------------
#ifdef _SAVE_TRACK_FOR_DEBUG_
  auto sav = *this;
#endif
  const float dz = z - mZ;
  if (std::abs(dz) < 1e-6f) {
    setZ(z);
    linRef0.setZ(z);
    return true;
  }
  const float kappa = (std::abs(bxyz[1]) < kTinyF) ? 0.f : linRef0.getCurvature(bxyz[1]); // kB2C*By*(q/pxz)
  if (std::abs(kappa) < kSmallKappa) {
    return propagateToZ(z, 0.f, linRef0); // for the straight-line propagation use 1D field method
  }
  const float K = kB2C * bxyz[1], bend = kappa * dz, abend = std::abs(bend);
  const float s0 = linRef0.getTx(), s1 = s0 + bend;
  if (std::abs(s0) > kAlmost1F || std::abs(s1) > kAlmost1F) {
    return false;
  }
  const float c0 = linRef0.getCosPsi(), c1 = getCosFromSin(s1), denom = c0 + c1;
  if (denom < kTinyF) {
    return false;
  }
  const float dx2dz = (s0 + s1) / denom, cps_dx2dz = (c1 + s1 * dx2dz);
  const float step = linRef0.getP2Pxz() * (abend < 0.05f ? dz * std::abs(cps_dx2dz) :                                      // chord
                                             2.f * std::asin(0.5f * dz * std::sqrt(1.f + dx2dz * dx2dz) * kappa) / kappa); // arc
  //
  //
  // get the track x,y,z,px/p,py/p,pz/p,p
  std::array<float, 7> vecLab{0.f};
  if (!linRef0.getPosDirGlo(vecLab)) {
    return false;
  }
  // rotate to the system where Bx=By=0.
  float bxy2 = bxyz[0] * bxyz[0] + bxyz[1] * bxyz[1];
  float bt = std::sqrt(bxy2);
  float cosphi = 1.f, sinphi = 0.f;
  if (bt > kTinyF) {
    cosphi = bxyz[0] / bt;
    sinphi = bxyz[1] / bt;
  }
  float bb = std::sqrt(bxy2 + bxyz[2] * bxyz[2]);
  float costet = 1.f, sintet = 0.f;
  if (bb > kTinyF) {
    costet = bxyz[2] / bb;
    sintet = bt / bb;
  }
  std::array<float, 7> vect{costet * cosphi * vecLab[0] + costet * sinphi * vecLab[1] - sintet * vecLab[2],
                            -sinphi * vecLab[0] + cosphi * vecLab[1],
                            sintet * cosphi * vecLab[0] + sintet * sinphi * vecLab[1] + costet * vecLab[2],
                            costet * cosphi * vecLab[3] + costet * sinphi * vecLab[4] - sintet * vecLab[5],
                            -sinphi * vecLab[3] + cosphi * vecLab[4],
                            sintet * cosphi * vecLab[3] + sintet * sinphi * vecLab[4] + costet * vecLab[5],
                            vecLab[6]};

  // Do the helix step
  float q = getCharge();
  g3helx3(q * bb, step, vect);

  // rotate back to the Global System
  vecLab[0] = cosphi * costet * vect[0] - sinphi * vect[1] + cosphi * sintet * vect[2];
  vecLab[1] = sinphi * costet * vect[0] + cosphi * vect[1] + sinphi * sintet * vect[2];
  vecLab[2] = -sintet * vect[0] + costet * vect[2];

  vecLab[3] = cosphi * costet * vect[3] - sinphi * vect[4] + cosphi * sintet * vect[5];
  vecLab[4] = sinphi * costet * vect[3] + cosphi * vect[4] + sinphi * sintet * vect[5];
  vecLab[5] = -sintet * vect[3] + costet * vect[5];

  // Do the final correcting step to the target plane (linear approximation)
  float x = vecLab[0], y = vecLab[1];
  auto dzFin = z - vecLab[2];
  if (std::abs(dzFin) > kTinyF) {
    if (std::abs(vecLab[5]) < kTinyF) {
      return false;
    }
    x += dzFin * vecLab[3] / vecLab[5]; // dz * px/pz
    y += dzFin * vecLab[4] / vecLab[5]; // dz * py/pz
  }
  // Calculate the linRef updated track parameters
  auto linRef1 = linRef0;
  auto t = 1.f / std::sqrt(vecLab[3] * vecLab[3] + vecLab[5] * vecLab[5]); // p / pxz
  linRef1.setZ(z);
  linRef1.setX(x);
  linRef1.setY(y);
  linRef1.setTx(vecLab[3] * t);
  linRef1.setTy(vecLab[4] * t);
  linRef1.setQ2Pxz(q * t / vecLab[6]);
  //
  // transport cov matrix
  // Recompute the Jacobian from the transported linearization reference, matching O2.
  prec_t snpRef0 = linRef0.getTx();
  prec_t snpRef1 = linRef1.getTx();
  prec_t cspRef0 = linRef0.getCosPsi();
  prec_t cspRef1 = linRef1.getCosPsi();
  prec_t cc = cspRef0 + cspRef1;
  if (cspRef0 < kTinyF || cspRef1 < kTinyF || cc < kTinyF) {
    return false;
  }
  prec_t cspRef0Inv = 1. / cspRef0;
  prec_t cspRef1Inv = 1. / cspRef1;
  prec_t ccInv = 1. / cc;
  prec_t dx2dzRef = (snpRef0 + snpRef1) * ccInv;
  prec_t dzccInv = dz * ccInv;
  prec_t hh = dzccInv * cspRef1Inv * (1. + cspRef0 * cspRef1 + snpRef0 * snpRef1);
  prec_t jj = dz * (dx2dzRef - snpRef1 * cspRef1Inv);
  prec_t f02 = hh * cspRef0Inv;
  prec_t f04 = hh * dzccInv * K;
  prec_t f24 = dz * K;
  prec_t f12 = linRef0.getTy() * (f02 * snpRef1 + jj);
  prec_t f13 = dz * (cspRef1 + snpRef1 * dx2dzRef);
  prec_t f14 = linRef0.getTy() * (f04 * snpRef1 + jj * f24);

  float diff[5];
  for (int i = 0; i < 5; i++) {
    diff[i] = getParam(i) - linRef0.getParam(i);
  }
  float snpUpd = linRef1.getTx() + diff[kTx] + f24 * diff[kQ2Pxz];
  if (std::abs(snpUpd) > kAlmost1F) {
    return false;
  }
  setZ(z);
  setX(linRef1.getX() + diff[kX] + f02 * diff[kTx] + f04 * diff[kQ2Pxz]);
  setY(linRef1.getY() + diff[kY] + f12 * diff[kTx] + f13 * diff[kTy] + f14 * diff[kQ2Pxz]);
  setTx(snpUpd);
  setTy(linRef1.getTy() + diff[kTy]);
  setQ2Pxz(linRef1.getQ2Pxz() + diff[kQ2Pxz]);
  linRef0 = linRef1; // update reference track after transporting this track, as in O2

  transportCovariance(f02, f04, f12, f13, f14, f24);
  return true;
}

void NA6PTrackParCov::transportCovariance(prec_t f02, prec_t f04, prec_t f12, prec_t f13, prec_t f14, prec_t f24)
{
  auto &C00 = mC[kXX], &C10 = mC[kYX], &C20 = mC[kTxX], &C30 = mC[kTyX], &C40 = mC[kQ2PxzX], &C11 = mC[kYY], &C21 = mC[kTxY],
       &C31 = mC[kTyY], &C41 = mC[kQ2PxzY], &C22 = mC[kTxTx], &C32 = mC[kTyTx], &C42 = mC[kQ2PxzTx], &C33 = mC[kTyTy], &C43 = mC[kQ2PxzTy], &C44 = mC[kQ2PxzQ2Pxz];

  // b = C*ft
  prec_t b00 = f02 * C20 + f04 * C40, b01 = f12 * C20 + f14 * C40 + f13 * C30;
  prec_t b02 = f24 * C40;
  prec_t b10 = f02 * C21 + f04 * C41, b11 = f12 * C21 + f14 * C41 + f13 * C31;
  prec_t b12 = f24 * C41;
  prec_t b20 = f02 * C22 + f04 * C42, b21 = f12 * C22 + f14 * C42 + f13 * C32;
  prec_t b22 = f24 * C42;
  prec_t b40 = f02 * C42 + f04 * C44, b41 = f12 * C42 + f14 * C44 + f13 * C43;
  prec_t b42 = f24 * C44;
  prec_t b30 = f02 * C32 + f04 * C43, b31 = f12 * C32 + f14 * C43 + f13 * C33;
  prec_t b32 = f24 * C43;

  // a = f*b = f*C*ft
  prec_t a00 = f02 * b20 + f04 * b40, a01 = f02 * b21 + f04 * b41, a02 = f02 * b22 + f04 * b42;
  prec_t a11 = f12 * b21 + f14 * b41 + f13 * b31, a12 = f12 * b22 + f14 * b42 + f13 * b32;
  prec_t a22 = f24 * b42;

  // F*C*Ft = C + (b + bt + a)
  C00 += b00 + b00 + a00;
  C10 += b10 + b01 + a01;
  C20 += b20 + b02 + a02;
  C30 += b30;
  C40 += b40;
  C11 += b11 + b11 + a11;
  C21 += b21 + b12 + a12;
  C31 += b31;
  C41 += b41;
  C22 += b22 + b22 + a22;
  C32 += b32;
  C42 += b42;

  checkCovariance();
}

// ---------------- chi2 and update ----------------
// Copied/compatible with your original (state still x,y observed).

float NA6PTrackParCov::getPredictedChi2(float xm, float ym, float sx2, float syx, float sy2) const
{
  auto exx = static_cast<prec_t>(getSigmaX2()) + static_cast<prec_t>(sx2);
  auto eyx = static_cast<prec_t>(getSigmaYX()) + static_cast<prec_t>(syx);
  auto eyy = static_cast<prec_t>(getSigmaY2()) + static_cast<prec_t>(sy2);
  auto det = exx * eyy - eyx * eyx;

  if (det < 1e-16)
    return 1.e16;

  float dx = getX() - xm;
  float dy = getY() - ym;
  float chi2 = (dx * (eyy * dx - eyx * dy) + dy * (exx * dy - dx * eyx)) / det;
  if (chi2 < 0.) {
    LOGP(warning, "Negative chi2={}, Cluster Sig2: {} {} {} Dx:{} Dy:{} | exx:{} eyx:{} eyy:{} det:{}",
         chi2, sx2, syx, sy2, dx, dy, exx, eyx, eyy, det);
    LOGP(warning, "Track: {}", asString());
  }
  return chi2;
}

bool NA6PTrackParCov::update(float xm, float ym, float sx2, float sxy, float sy2)
{
#ifdef _SAVE_TRACK_FOR_DEBUG_
  auto sav = *this;
#endif
  auto& C00 = mC[kXX];
  auto& C10 = mC[kYX];
  auto& C20 = mC[kTxX];
  auto& C30 = mC[kTyX];
  auto& C40 = mC[kQ2PxzX];
  auto& C11 = mC[kYY];
  auto& C21 = mC[kTxY];
  auto& C31 = mC[kTyY];
  auto& C41 = mC[kQ2PxzY];
  auto& C22 = mC[kTxTx];
  auto& C32 = mC[kTyTx];
  auto& C42 = mC[kQ2PxzTx];
  auto& C33 = mC[kTyTy];
  auto& C43 = mC[kQ2PxzTy];
  auto& C44 = mC[kQ2PxzQ2Pxz];

  prec_t r00 = static_cast<prec_t>(sx2) + static_cast<prec_t>(C00);
  prec_t r01 = static_cast<prec_t>(sxy) + static_cast<prec_t>(C10);
  prec_t r11 = static_cast<prec_t>(sy2) + static_cast<prec_t>(C11);
  const prec_t det = r00 * r11 - r01 * r01;
  if (std::abs(det) < 1e-20) {
    return false;
  }

  const prec_t tmp = r00;
  r00 = r11 / det;
  r11 = tmp / det;
  r01 = -r01 / det;

  const prec_t k00 = C00 * r00 + C10 * r01;
  const prec_t k01 = C00 * r01 + C10 * r11;
  const prec_t k10 = C10 * r00 + C11 * r01;
  const prec_t k11 = C10 * r01 + C11 * r11;
  const prec_t k20 = C20 * r00 + C21 * r01;
  const prec_t k21 = C20 * r01 + C21 * r11;
  const prec_t k30 = C30 * r00 + C31 * r01;
  const prec_t k31 = C30 * r01 + C31 * r11;
  const prec_t k40 = C40 * r00 + C41 * r01;
  const prec_t k41 = C40 * r01 + C41 * r11;

  const prec_t dx = static_cast<prec_t>(xm) - static_cast<prec_t>(getX());
  const prec_t dy = static_cast<prec_t>(ym) - static_cast<prec_t>(getY());
  const prec_t txUpd = static_cast<prec_t>(getTx()) + k20 * dx + k21 * dy;
  if (std::abs(txUpd) > kAlmost1F) {
    return false;
  }

  mP[kX] += k00 * dx + k01 * dy;
  mP[kY] += k10 * dx + k11 * dy;
  mP[kTx] = txUpd;
  mP[kTy] += k30 * dx + k31 * dy;
  mP[kQ2Pxz] += k40 * dx + k41 * dy;

  const prec_t c01 = C10;
  const prec_t c02 = C20;
  const prec_t c03 = C30;
  const prec_t c04 = C40;
  const prec_t c12 = C21;
  const prec_t c13 = C31;
  const prec_t c14 = C41;

  C00 -= k00 * C00 + k01 * C10;
  C10 -= k00 * c01 + k01 * C11;
  C20 -= k00 * c02 + k01 * c12;
  C30 -= k00 * c03 + k01 * c13;
  C40 -= k00 * c04 + k01 * c14;

  C11 -= k10 * c01 + k11 * C11;
  C21 -= k10 * c02 + k11 * c12;
  C31 -= k10 * c03 + k11 * c13;
  C41 -= k10 * c04 + k11 * c14;

  C22 -= k20 * c02 + k21 * c12;
  C32 -= k20 * c03 + k21 * c13;
  C42 -= k20 * c04 + k21 * c14;

  C33 -= k30 * c03 + k31 * c13;
  C43 -= k30 * c04 + k31 * c14;

  C44 -= k40 * c04 + k41 * c14;
  checkCovariance();
  return true;
}

// ---------------- correlation checks / reset ----------------
// For brevity: keep identical semantics; adapt if you had extra project-specific logic.

void NA6PTrackParCov::resetCovariance(float s2)
{
  // Reset the covarince matrix to "something big"
  double d0(kCX2Ini), d1(kCY2Ini), d2(kCTX2Ini), d3(kCTY2Ini), d4(kC1Pxz2Ini);
  if (s2 > kTinyF) { // if negative, leave at ini error
    d0 = getSigmaX2() * s2;
    d1 = getSigmaY2() * s2;
    d2 = getSigmaTx2() * s2;
    d3 = getSigmaTy2() * s2;
    d4 = getSigmaQ2Pxz2() * s2;
    if (d0 > kCX2max) {
      d0 = kCX2max;
    }
    if (d1 > kCY2max) {
      d1 = kCY2max;
    }
    if (d2 > kCTX2max) {
      d2 = kCTX2max;
    }
    if (d3 > kCTY2max) {
      d3 = kCTY2max;
    }
    if (d4 > kC1Pxz2max) {
      d4 = kC1Pxz2max;
    }
  }
  for (int i = 0; i < 15; i++) {
    mC[i] = 0;
  }
  mC[kXX] = d0;
  mC[kYY] = d1;
  mC[kTxTx] = d2;
  mC[kTyTy] = d3;
  mC[kQ2PxzQ2Pxz] = d4;
  if (s2 == MaxErrSelRescale) {
    auto useQ2Pxz = getParam(4);
    useQ2Pxz += TMath::Sign(MinQ2PXZBias, useQ2Pxz);
    mC[kQ2PxzQ2Pxz] *= useQ2Pxz * useQ2Pxz;
  }
}

void NA6PTrackParCov::printCorr() const
{
  // optional: implement like your original if needed
}

bool NA6PTrackParCov::isCovariancePositiveDefinite(const NA6PTrackParCov::CovArray& c)
{
  constexpr int map[5][5] = {
    {0, 1, 3, 6, 10},
    {1, 2, 4, 7, 11},
    {3, 4, 5, 8, 12},
    {6, 7, 8, 9, 13},
    {10, 11, 12, 13, 14}};

  double a[5][5];
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      a[i][j] = a[j][i] = static_cast<double>(c[map[i][j]]);
    }
  }

  // In-place Cholesky. No heap allocation.
  constexpr double epsAbs = 1e-20;
  constexpr double epsRel = 1e-12;

  for (int i = 0; i < 5; ++i) {
    double s = a[i][i];
    for (int k = 0; k < i; ++k) {
      s -= a[i][k] * a[i][k];
    }

    const double minPivot = std::max(epsAbs, epsRel * std::abs(static_cast<double>(c[map[i][i]])));
    if (!(s > minPivot)) {
#ifdef _PRINT_COV_PROBLEMS_
      LOGP(warn, "ANOMALY: s={}  < minPivot {}", s, minPivot);
#endif
      return false;
    }

    const double lii = std::sqrt(s);
    a[i][i] = lii;

    for (int j = i + 1; j < 5; ++j) {
      double v = a[j][i];
      for (int k = 0; k < i; ++k) {
        v -= a[j][k] * a[i][k];
      }
      a[j][i] = v / lii;
    }
  }
  return true;
}

void NA6PTrackParCov::checkCovariancePosDef()
{
  if (!isCovariancePositiveDefinite(mC)) {
#ifdef _PRINT_COV_PROBLEMS_
    LOGP(warning, "Non positive-definite covariance detected: {}", asString());
#endif //_PRINT_COV_PROBLEMS_

#ifdef _FIX_COV_PROBLEMS_
    constexpr int OffDiag[10] = {1, 3, 6, 10, 4, 7, 11, 8, 12, 13};

    // Preserve variances, progressively reduce all correlations.
    float scale = 0.5f;
    for (int iter = 0; iter < 12 && !isCovariancePositiveDefinite(mC); ++iter) {
      for (int k = 0; k < 10; ++k) {
        mC[OffDiag[k]] *= scale;
      }
    }
    // Last-resort repair: diagonal covariance is positive definite after
    // the earlier diagonal abs/bounds pass.
    if (!isCovariancePositiveDefinite(mC)) {
      for (int k = 0; k < 10; ++k) {
        mC[OffDiag[k]] = 0.f;
      }
    }
#endif //_FIX_COV_PROBLEMS_
  }
}

void NA6PTrackParCov::checkCovariance()
{
#ifdef _CHECK_COV_POSDEF_
  checkCovariancePosDef();
#endif //_CHECK_COV_POSDEF_
  // In case the diagonal element is bigger than the maximal allowed value, it is set to
  // the limit and the off-diagonal elements that correspond to it are scaled accordingly
#ifdef _LIMIT_COV_TO_MAX_ERR_
  for (int i = 0; i < 5; i++) {
    if (mC[CovarDiag[i]] < 0.f) {
#ifdef _PRINT_COV_PROBLEMS_
      LOGP(warning, "Negative covariance diag detected: {}", asString());
#endif //_PRINT_COV_PROBLEMS_
#ifdef _FIX_COV_PROBLEMS_
      mC[CovarDiag[i]] = -mC[CovarDiag[i]];
#endif //_PRINT_COV_PROBLEMS_
    }
    const auto maxD = i < 4 ? kCMaxDiag[i] : kCMaxDiag[i] * mP[4] * mP[4];
    if (mC[CovarDiag[i]] > maxD) {
      auto scl = std::sqrt(maxD / mC[CovarDiag[i]]);
      //      LOGP(info, "DBGR covLimit diag={} before={} max={} scale={} state={}", i, mC[CovarDiag[i]], maxD, scl, asString());
      for (int j = 0; j < 5; j++) {
        mC[CovarMap[i][j]] *= scl;
      }
      mC[CovarDiag[i]] = maxD;
    }
  }
#endif // _LIMIT_COV_TO_MAX_ERR_

#ifdef _CHECK_BAD_CORRELATIONS_
  // Keep all 2x2 principal minors positive. Without this, propagation or
  // diagonal limiting can leave an indefinite covariance with positive
  // diagonal elements; the next Kalman update then legitimately produces
  // negative variances.
  constexpr float MaxCorr = 0.99f;
  for (int i = 1; i < 5; ++i) {
    for (int j = 0; j < i; ++j) {
      const auto sig2 = mC[CovarDiag[i]] * mC[CovarDiag[j]];
      auto& cov = mC[CovarMap[i][j]];
      if (sig2 <= 0.f) {
        cov = 0.f;
        continue;
      }
      const auto maxCov = MaxCorr * std::sqrt(sig2);
      if (std::abs(cov) > maxCov) {
#ifdef _PRINT_COV_PROBLEMS_
        LOGP(warning, "Bad covariance correlation detected for {} {}: {} vs diag {} {} -> {}",
             i, j, cov, mC[CovarDiag[i]], mC[CovarDiag[j]], cov / std::sqrt(sig2));
#endif //_PRINT_COV_PROBLEMS_
#ifdef _FIX_COV_PROBLEMS_
        cov = std::copysign(maxCov, cov);
#endif //_FIX_COV_PROBLEMS_
      }
    }
  }
#endif //_CHECK_BAD_CORRELATIONS_
}

bool NA6PTrackParCov::correctForMaterial(float x2x0, float xrho, float density, float atomicZ, float zOverA, bool anglecorr)
{
  //------------------------------------------------------------------
  // This function corrects the track parameters for the crossed material.
  // "x2x0"   - X/X0, the thickness in units of the radiation length.
  // "xrho" - is the product length*density (g/cm^2).
  //     It should be passed as negative when propagating tracks
  //     from the intreaction point to the outside of the central barrel.
  // "density" - mean density (g/cm^3)
  // "atomicZ" - mean Z
  // "zOverA" - mean Z/A
  // "anglecorr" - switch for the angular correction
  //------------------------------------------------------------------
#ifdef _SAVE_TRACK_FOR_DEBUG_
  auto sav = *this;
#endif
  constexpr float kMSConst2 = 0.0136f * 0.0136f;
  constexpr float kMinP = 0.01f; // kill below this momentum

  float csp2 = getCosPsi2();
  float cst2I = (1.f + getTy() * getTy()); // 1/cos(lambda)^2
  if (anglecorr) {                         // Apply angle correction, if requested
    float angle = std::sqrt(cst2I / csp2);
    x2x0 *= angle;
    xrho *= angle;
  }
  auto m = getPID().getMass();
  int charge = getCharge();
  int charge2 = charge * charge; // in case we introduce charge > 1 particle: getAbsCharge() * getAbsCharge();
  float p = getP(), p0 = p, p02 = p * p, e2 = p02 + getPID().getMass2(), massInv = 1.f / m, bg = p * massInv, dETot = 0.f;
  float e = std::sqrt(e2), e0 = e;
  if (m > 0 && xrho != 0.f) {
    float mI = (atomicZ < 13.f) ? (12.f * atomicZ + 7.f) * 1.e-9f : (9.76f * atomicZ + 58.8f * std::pow(atomicZ, -0.19f)) * 1.e-9f;
    float ekin = e - m, dedx = BetheBlochSolid(bg, density, 0.2f, 3.f, mI, zOverA);
#ifdef _BB_NONCONST_CORR_
    float dedxDer = 0.f, dedx1 = dedx;
#endif
    if (charge2 != 1) {
      dedx *= charge2;
    }
    float dE = dedx * xrho;
    int na = 1 + int(std::abs(dE) / ekin * ELoss2EKinThreshInv);
    if (na > MaxELossIter) {
      na = MaxELossIter;
    }
    if (na > 1) {
      dE /= na;
      xrho /= na;
#ifdef _BB_NONCONST_CORR_
      dedxDer = getBetheBlochSolidDerivativeApprox(dedx1, bg, mI); // require correction for non-constantness of dedx vs betagamma
      if (charge2 != 1) {
        dedxDer *= charge2;
      }
#endif
    }
    while (na--) {
#ifdef _BB_NONCONST_CORR_
      if (dedxDer != 0.) { // correction for non-constantness of dedx vs beta*gamma (in linear approximation): for a single step dE -> dE * [(exp(dedxDer) - 1)/dedxDer]
        if (xrho < 0) {
          dedxDer = -dedxDer; // E.loss ( -> positive derivative)
        }
        auto corrC = (std::exp(dedxDer) - 1.) / dedxDer;
        dE *= corrC;
      }
#endif
      e += dE;
      if (e > m) { // stopped
        p = std::sqrt(e * e - getPID().getMass2());
      } else {
        return false;
      }
      if (na) {
        bg = p * massInv;
        dedx = getdEdxBBOpt(bg);
#ifdef _BB_NONCONST_CORR_
        dedxDer = getBetheBlochSolidDerivativeApprox(dedx, bg, mI);
#endif
        if (charge2 != 1) {
          dedx *= charge2;
#ifdef _BB_NONCONST_CORR_
          dedxDer *= charge2;
#endif
        }
        dE = dedx * xrho;
      }
    }

    if (p < kMinP) {
      return false;
    }
    dETot = e - e0;
  } // end of e.loss correction

  // Calculating the multiple scattering corrections******************
  precCov_t& fC22 = mC[kTxTx];
  precCov_t& fC33 = mC[kTyTy];
  precCov_t& fC43 = mC[kQ2PxzTy];
  precCov_t& fC44 = mC[kQ2PxzQ2Pxz];
  //
  precCov_t cC22(0.), cC33(0.), cC43(0.), cC44(0.);
  if (x2x0 != 0.f) {
    float beta2 = p02 / e2, theta2 = kMSConst2 / (beta2 * p02) * std::abs(x2x0);
    float fp34 = getTy() * getCharge2Pxz();
    if (charge2 != 1) {
      theta2 *= charge2;
    }
    if (theta2 > phys_const::PI * phys_const::PI) {
      return false;
    }
    float t2c2I = theta2 * cst2I;
    cC22 = t2c2I * csp2;
    cC33 = t2c2I * cst2I;
    cC43 = t2c2I * fp34;
    cC44 = theta2 * fp34 * fp34;
  }

  // the energy loss correction contribution to cov.matrix: approximate energy loss fluctuation (M.Ivanov)
  float sigmadE = ELoss2SigmaE * std::sqrt(std::abs(dETot)) * e0 / p02 * getCharge2Pxz();
  cC44 += sigmadE * sigmadE;

  // Applying the corrections*****************************
  fC22 += cC22;
  fC33 += cC33;
  fC43 += cC43;
  fC44 += cC44;

  mP[kQ2Pxz] *= p0 / p;

  checkCovariance();

  return true;
}

bool NA6PTrackParCov::correctForMaterial(float x2x0, float xrho, float density, float atomicZ, float zOverA,
                                       NA6PTrackPar& linRef, bool anglecorr, float radiativeDEDX)
{
  //------------------------------------------------------------------
  // This function corrects the track parameters for the crossed material.
  // "x2x0"   - X/X0, the thickness in units of the radiation length.
  // "xrho" - is the product length*density (g/cm^2).
  //     It should be passed as negative when propagating tracks
  //     from the intreaction point to the outside of the central barrel.
  // "radiativeDEDX" - mean radiative stopping power in GeV/(g/cm^2).
  // "anglecorr" - switch for the angular correction
  //------------------------------------------------------------------
#ifdef _SAVE_TRACK_FOR_DEBUG_
  auto sav = *this;
#endif
  constexpr float kMSConst2 = 0.0136f * 0.0136f;
  constexpr float kMinP = 0.01f; // kill below this momentum

  float csp2 = linRef.getCosPsi2();
  float cst2I = (1.f + linRef.getTy() * linRef.getTy()); // 1/cos(lambda)^2
  if (anglecorr) {                                       // Apply angle correction, if requested
    float angle = std::sqrt(cst2I / csp2);
    x2x0 *= angle;
    xrho *= angle;
  }
  auto m = getPID().getMass();
  int charge = getPID().getCharge();
  int charge2 = charge * charge;
  float p = linRef.getP(), p0 = p, p02 = p * p, e2 = p02 + getPID().getMass2(), massInv = 1. / m, bg = p * massInv, dETot = 0.;
  float e = std::sqrt(e2), e0 = e;
  if (m > 0 && xrho != 0.f) {
    float mI = (atomicZ < 13.f) ? (12.f * atomicZ + 7.f) * 1.e-9f : (9.76f * atomicZ + 58.8f * std::pow(atomicZ, -0.19f)) * 1.e-9f;
    float ekin = e - m, dedx = BetheBlochSolid(bg, density, 0.2, 3, mI, zOverA);
#ifdef _BB_NONCONST_CORR_
    float dedxDer = 0., dedx1 = dedx;
#endif
    if (charge2 != 1) {
      dedx *= charge2;
    }
    float dE = dedx * xrho;
    const float dERadTotal = radiativeDEDX * xrho;
    int na = 1 + int(std::abs(dE + dERadTotal) / ekin * ELoss2EKinThreshInv);
    if (na > MaxELossIter) {
      na = MaxELossIter;
    }
    if (na > 1) {
      dE /= na;
      xrho /= na;
#ifdef _BB_NONCONST_CORR_
      dedxDer = getBetheBlochSolidDerivativeApprox(dedx1, bg, mI); // require correction for non-constantness of dedx vs betagamma
      if (charge2 != 1) {
        dedxDer *= charge2;
      }
#endif
    }
    const float dERad = radiativeDEDX * xrho;
    while (na--) {
#ifdef _BB_NONCONST_CORR_
      if (dedxDer != 0.) { // correction for non-constantness of dedx vs beta*gamma (in linear approximation): for a single step dE -> dE * [(exp(dedxDer) - 1)/dedxDer]
        if (xrho < 0) {
          dedxDer = -dedxDer; // E.loss ( -> positive derivative)
        }
        auto corrC = (std::exp(dedxDer) - 1.) / dedxDer;
        dE *= corrC;
      }
#endif
      e += dE + dERad;
      if (e > m) { // stopped
        p = std::sqrt(e * e - getPID().getMass2());
      } else {
        return false;
      }
      if (na) {
        bg = p * massInv;
        dedx = getdEdxBBOpt(bg);
#ifdef _BB_NONCONST_CORR_
        dedxDer = getBetheBlochSolidDerivativeApprox(dedx, bg, mI);
#endif
        if (charge2 != 1) {
          dedx *= charge2;
#ifdef _BB_NONCONST_CORR_
          dedxDer *= charge2;
#endif
        }
        dE = dedx * xrho;
      }
    }

    if (p < kMinP) {
      return false;
    }
    dETot = e - e0 - dERadTotal; // Keep the existing fluctuation estimate tied to ionization.
  } // end of e.loss correction

  // Calculating the multiple scattering corrections******************
  precCov_t& fC22 = mC[kTxTx];
  precCov_t& fC33 = mC[kTyTy];
  precCov_t& fC43 = mC[kQ2PxzTy];
  precCov_t& fC44 = mC[kQ2PxzQ2Pxz];
  //
  precCov_t cC22(0.), cC33(0.), cC43(0.), cC44(0.);
  if (x2x0 != 0.f) {
    float beta2 = p02 / e2, theta2 = kMSConst2 / (beta2 * p02) * std::abs(x2x0);
    float fp34 = linRef.getTy() * linRef.getCharge2Pxz();
    if (charge2 != 1) {
      theta2 *= charge2;
    }
    if (theta2 > phys_const::PI * phys_const::PI) {
      return false;
    }
    float t2c2I = theta2 * cst2I;
    cC22 = t2c2I * csp2;
    cC33 = t2c2I * cst2I;
    cC43 = t2c2I * fp34;
    cC44 = theta2 * fp34 * fp34;
  }

  // the energy loss correction contribution to cov.matrix: approximate energy loss fluctuation (M.Ivanov)
  float sigmadE = ELoss2SigmaE * std::sqrt(std::abs(dETot)) * e0 / p02 * getCharge2Pxz();
  cC44 += sigmadE * sigmadE;

  // Applying the corrections*****************************
  fC22 += cC22;
  fC33 += cC33;
  fC43 += cC43;
  fC44 += cC44;
  auto fact = p0 / p;
  linRef.setQ2Pxz(linRef.getQ2Pxz() * fact);
  if (&linRef != this) {
    mP[kQ2Pxz] *= fact;
  }
  checkCovariance();

  return true;
}

//______________________________________________
float NA6PTrackParCov::getPredictedChi2(const NA6PTrackParCov& rhs) const
{
  MatrixD5Sym cov; // perform matrix operations in double!
  return getPredictedChi2(rhs, cov);
}

//______________________________________________
float NA6PTrackParCov::getPredictedChi2Quiet(const NA6PTrackParCov& rhs) const
{
  if (std::abs(getZ() - rhs.getZ()) > kTinyF) {
    return 2.f / kTinyF;
  }
  MatrixD5Sym cov;
  buildCombinedCovMatrix(rhs, cov);
  if (!cov.Invert()) {
    return 2.f / kTinyF;
  }
  double chi2diag = 0., chi2ndiag = 0., diff[5];
  for (int i = 5; i--;) {
    diff[i] = getParam(i) - rhs.getParam(i);
    chi2diag += diff[i] * diff[i] * cov(i, i);
  }
  for (int i = 5; i--;) {
    for (int j = i; j--;) {
      chi2ndiag += diff[i] * diff[j] * cov(i, j);
    }
  }
  return chi2diag + 2. * chi2ndiag;
}

//______________________________________________
void NA6PTrackParCov::buildCombinedCovMatrix(const NA6PTrackParCov& rhs, MatrixD5Sym& cov) const
{
  // fill combined cov.matrix (NOT inverted)
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j <= i; j++) {
      cov(i, j) = getCovMatElem(i, j) + rhs.getCovMatElem(i, j);
    }
  }
}

//______________________________________________
float NA6PTrackParCov::getPredictedChi2(const NA6PTrackParCov& rhs, MatrixD5Sym& covToSet) const
{
  // get chi2 wrt other track, which must be defined at the same parameters X,alpha
  // Supplied non-initialized covToSet matrix is filled by inverse combined matrix for further use
  if (std::abs(getZ() - rhs.getZ()) > kTinyF) {
    LOGP(error, "The reference Z of the tracks differ: {} : {}", getZ(), rhs.getZ());
    return 2.f / kTinyF;
  }
  buildCombinedCovMatrix(rhs, covToSet);
  if (!covToSet.Invert()) {
    LOG(warning) << "Cov.matrix inversion failed: " << covToSet;
    return 2.f / kTinyF;
  }
  double chi2diag = 0., chi2ndiag = 0., diff[5];
  for (int i = 5; i--;) {
    diff[i] = getParam(i) - rhs.getParam(i);
    chi2diag += diff[i] * diff[i] * covToSet(i, i);
  }
  for (int i = 5; i--;) {
    for (int j = i; j--;) {
      chi2ndiag += diff[i] * diff[j] * covToSet(i, j);
    }
  }
  return chi2diag + 2. * chi2ndiag;
}

//______________________________________________
bool NA6PTrackParCov::update(const NA6PTrackParCov& rhs, const MatrixD5Sym& covInv)
{
  // update track with other track, the inverted combined cov matrix should be supplied

  // consider skipping this check, since it is usually already done upstream
  if (std::abs(this->getZ() - rhs.getZ()) > kTinyF) {
    LOGP(error, "The reference Z of the tracks differ: {} : {}", getZ(), rhs.getZ());
    return false;
  }
  // gain matrix K = Cov0*H*(Cov0+Cov0)^-1 (for measurement matrix H=I)
  MatrixD5Sym matC0;
  double diff[5];
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j <= i; j++) {
      matC0(i, j) = getCovMatElem(i, j);
    }
    diff[i] = rhs.getParam(i) - getParam(i);
  }
  MatrixD5 matK = matC0 * covInv;

  for (int i = 4; i >= 0; --i) {
    for (int j = 4; j >= 0; --j) {
      incParam(i, matK(i, j) * diff[j]);
    }
  }

  // updated covariance: Cov0 = Cov0 - K*Cov0
  matK *= MatrixD5(matC0);
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      incCovMatElem(i, j, -matK(i, j));
    }
  }
  return true;
}

//______________________________________________
bool NA6PTrackParCov::update(const NA6PTrackParCov& rhs)
{
  // update track with other track
  MatrixD5Sym covI; // perform matrix operations in double!
  buildCombinedCovMatrix(rhs, covI);
  if (!covI.Invert()) {
    LOG(warning) << "Cov.matrix inversion failed: " << covI;
    return false;
  }
  return update(rhs, covI);
}
