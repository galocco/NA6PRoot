// NA6PCCopyright

#include "NA6PTrackPar.h"
#include <fairlogger/Logger.h>

void NA6PTrackPar::initParam(const float* xyz, const float* pxyz, int charge)
{
  mZ = xyz[2];
  mP[kX] = xyz[0];
  mP[kY] = xyz[1];

  const float px = pxyz[0];
  const float py = pxyz[1];
  const float pz = pxyz[2]; // assume pz>0 in your convention
  if (pz <= kTinyF) {
    LOGP(fatal, "Track model requires positive pZ > {}, {} was provided", kTinyF, pz);
  }
  const float pxz = std::hypot(px, pz);
  mP[kTx] = px / pxz; // sin(psi)
  mP[kTy] = py / pxz; // py/pxz
  mP[kQ2Pxz] = (charge == 0 ? 1.f : static_cast<float>(charge)) / pxz;
  if (charge == 0) {
    mPID = PID::PI0;
  } else if (std::abs(charge) == 2) {
    mPID = PID::Alpha;
  }
}

bool NA6PTrackPar::propagateParamToZ(float z, float by)
{
  const float dz = z - mZ;
  if (std::abs(dz) < 1e-6f) {
    mZ = z;
    return true;
  }

  const float kappa = getCurvature(by); // kB2C*By*(q/pxz)
  const float bend = kappa * dz, abend = std::abs(bend);

  if (std::abs(kappa) < kSmallKappa || abend < kSmallBend) { // Straight/near-straight transport: tx, ty, q/pxz unchanged to this order
    const float s0 = getTx(), c0Idz = dz / getCosPsi();
    mP[kX] += s0 * c0Idz;      // dx/dz = tan(psi) = sin/cos = tx/cos
    mP[kY] += getTy() * c0Idz; // dy/dz = (py/pz) = (py/pxz)/(pz/pxz) = ty/cos
    mZ = z;
    return true;
  }

  // Exact dipole transport using z-geometry:  // sin(psi1) = sin(psi0) + kappa*dz
  const float s0 = getTx(), s1 = s0 + bend;

  if (std::abs(s1) >= kAlmost1F) { // Physically this would correspond to psi leaving [-pi/2,pi/2] (pz changing sign).
    return false;
  }

  const float c0 = getCosPsi(), c1 = getCosFromSin(s1), denom = c0 + c1;
  if (denom < kTinyF) {
    return false;
  }
  const float dx2dz = (s0 + s1) / denom;

  if (abend < 0.05f) {
    mP[kY] += dz * (c1 + s1 * dx2dz) * getTy();
  } else {
    // for small bends the linear apporximation of the arc by the segment is OK, but at large bends need precise value
    // angle traversed delta = 2*asin(dist_start_end * |kappa| / 2), hence the arc is: deltaPhi/|kappa|
    // The dist_start_end is obtained from sqrt(dx^2+dz^2) = z/(c0+c1)*sqrt(2+s0*s1+c0*c1)
    auto arg = c0 * s1 - c1 * s0;
    if (std::abs(arg) > kAlmost1F) {
      return false; // loop
    }
    float rot = std::asin(arg);
    if (s0 * s0 + s1 * s1 > 1.f && s0 * s1 < 0) {
      if (s1 > 0.f) {
        rot = phys_const::PI - rot;
      } else {
        rot = -phys_const::PI - rot;
      }
    }
    mP[kY] += getTy() / kappa * rot;
  }

  mP[kX] += dx2dz * dz;
  mP[kTx] += bend;
  mZ = z;
  return true; // ty and q/pxz invariant in pure By
}

bool NA6PTrackPar::propagateParamToZ(float z, const float* bxyz)
{
  //----------------------------------------------------------------
  // Extrapolate this track params (w/o cov matrix) to the plane z in the field bxyz.
  //----------------------------------------------------------------
  const float dz = z - mZ;
  if (std::abs(dz) < 1e-6f) {
    mZ = z;
    return true;
  }
  const float kappa = getCurvature(bxyz[1]); // kB2C*By*(q/pxz)
  const float field2 = bxyz[0] * bxyz[0] + bxyz[1] * bxyz[1] + bxyz[2] * bxyz[2];
  if (getCharge() == 0 || field2 < kTinyF * kTinyF) {
    return propagateParamToZ(z, 0.f); // for the straight-line propagation use 1D field method
  }
  const float bend = kappa * dz, abend = std::abs(bend);
  const float s0 = getTx(), s1 = s0 + bend;
  if (std::abs(s0) > kAlmost1F || std::abs(s1) > kAlmost1F) {
    return false;
  }
  const float c0 = getCosPsi(), c1 = getCosFromSin(s1), denom = c0 + c1;
  if (denom < kTinyF) {
    return false;
  }
  const float dx2dz = (s0 + s1) / denom;
  const float cpsDx2dz = c1 + s1 * dx2dz;
  const float step = getP2Pxz() * (abend < 0.05f ? dz * std::abs(cpsDx2dz) :                                       // chord
                                     2.f * std::asin(0.5f * dz * std::sqrt(1.f + dx2dz * dx2dz) * kappa) / kappa); // arc
  //
  // get the track x,y,z,px/p,py/p,pz/p,p
  std::array<float, 7> vecLab{0.f};
  if (!getPosDirGlo(vecLab)) {
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
  mP[kX] = x;
  mP[kY] = y;

  // Calculate the track parameters
  auto t = 1.f / std::sqrt(vecLab[3] * vecLab[3] + vecLab[5] * vecLab[5]); // p / pxz
  mZ = z;
  mP[kTx] = vecLab[3] * t;
  mP[kTy] = vecLab[4] * t;
  mP[kQ2Pxz] = q * t / vecLab[6];
  return true;
}

bool NA6PTrackPar::getPosDirGlo(std::array<float, 7>& posdirp) const
{
  // fill vector with lab x,y,z,px/p,py/p,pz/p,p,sinAlpha,cosAlpha
  posdirp[0] = getX();
  posdirp[1] = getY();
  posdirp[2] = getZ();
  auto p2pxz = getP2Pxz();
  auto pxz2p = 1.f / p2pxz;
  posdirp[3] = getTx() * pxz2p;     // px/p
  posdirp[4] = getTy() * pxz2p;     // py/p
  posdirp[5] = getCosPsi() * pxz2p; // pz/p
  posdirp[6] = getPxz() * p2pxz;    // p
  return true;
}

bool NA6PTrackPar::correctForELoss(float xrho, float density, float atomicZ, float zOverA, bool anglecorr)
{
  //------------------------------------------------------------------
  // This function corrects the track parameters for the energy loss in crossed material.
  // "xrho" - is the product length*density (g/cm^2).
  //     It should be passed as negative when propagating tracks
  //     from the intreaction point to the outside of the central barrel.
  // "density" - mean density (g/cm^3)
  // "atomicZ" - mean Z
  // "zOverA" - mean Z/A
  // "anglecorr" - switch for the angular correction
  //------------------------------------------------------------------
  constexpr float kMinP = 0.01f; // kill below this momentum

  auto m = getPID().getMass();
  if (m > 0 && xrho != 0.f) {
    // Apply angle correction, if requested
    if (anglecorr) {
      float csp2 = getCosPsi2();               // cos(psi)^2
      float cst2I = (1.f + getTy() * getTy()); // 1/cos(lambda)^2
      float angle = std::sqrt(cst2I / (csp2));
      xrho *= angle;
    }
    int charge2 = 1; // in case we introduce charge > 1 particle: getAbsCharge() * getAbsCharge();
    float p = getP(), p0 = p, p2 = p * p, e2 = p2 + getPID().getMass2(), massInv = 1. / m, bg = p * massInv;
    float mI = (atomicZ < 13.f) ? (12.f * atomicZ + 7.f) * 1.e-9f : (9.76f * atomicZ + 58.8f * std::pow(atomicZ, -0.19f)) * 1.e-9f;
    float e = std::sqrt(e2), ekin = e - m, dedx = BetheBlochSolid(bg, density, 0.2f, 3.f, mI, zOverA);
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
      dedxDer = getBetheBlochSolidDerivative(dedx1, bg, mI); // require correction for non-constantness of dedx vs betagamma
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
        dedxDer = getBetheBlochSolidDerivative(dedx, bg, mI);
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
    mP[kQ2Pxz] *= p0 / p;
  }

  return true;
}

std::string NA6PTrackPar::asString() const
{
  return fmt::format("Z:{:+.4e} Par: {:+.4e} {:+.4e} {:+.4e} {:+.4e} {:+.4e} {:s}", getZ(), getX(), getY(), getTx(), getTy(), getQ2Pxz(), getPID().getName());
}

float NA6PTrackPar::BetheBlochSolid(float bg, float rho, float kp1, float kp2, float meanI, float meanZA)
{
  //
  // This is the parameterization of the Bethe-Bloch formula inspired by Geant.
  //
  // bg  - beta*gamma
  // rho - density [g/cm^3]
  // kp1 - density effect first junction point
  // kp2 - density effect second junction point
  // meanI - mean excitation energy [GeV]
  // meanZA - mean Z/A
  //
  // The default values for the kp* parameters are for silicon.
  // The returned value is in [GeV/(g/cm^2)].
  //
  constexpr float mK = 0.307075e-3; // [GeV*cm^2/g]
  constexpr float me = 0.511e-3;    // [GeV/c^2]
  kp1 *= 2.303f;
  kp2 *= 2.303f;
  clampBetaGamma(bg);
  float bg2 = bg * bg, beta2 = bg2 / (1 + bg2);
  float maxT = 2.f * me * bg2; // neglecting the electron mass

  //*** Density effect
  float d2 = 0.;
  const float x = std::log(bg);
  const float lhwI = std::log(28.816f * 1e-9f * std::sqrt(rho * meanZA) / meanI);
  if (x > kp2) {
    d2 = lhwI + x - 0.5f;
  } else if (x > kp1) {
    double r = (kp2 - x) / (kp2 - kp1);
    d2 = lhwI + x - 0.5f + (0.5f - lhwI - kp1) * r * r * r;
  }
  auto dedx = mK * meanZA / beta2 * (0.5f * std::log(2 * me * bg2 * maxT / (meanI * meanI)) - beta2 - d2);
  return dedx > 0. ? dedx : 0.;
}

float NA6PTrackPar::BetheBlochSolidDerivative(float dedx, float bg, float meanZA)
{
  //
  // This is approximate derivative of the BB over betagamma, NO check for the consistency of the provided dedx and bg is done
  // Charge 1 particle is assumed for the provied dedx. For charge > 1 particles dedx/q^2 should be provided and obtained value must be scaled by q^2
  // The call should be usually done as
  // auto dedx = BetheBlochSolid(bg, rho, kp1, kp2, meanI, meanZA);
  // // if derivative needed
  // auto ddedx = BetheBlochSolidDerivative(dedx, bg, bg*bg)
  //
  // dedx - precalculated dedx for bg
  // bg  - beta*gamma
  //
  constexpr float mK = 0.307075e-3; // [GeV*cm^2/g]
  clampBetaGamma(bg);
  auto bg2 = bg * bg;
  auto t1 = 1 + bg2;
  //  auto derH = (mK * meanZA * (t1+bg2) - dedx*bg2)/(bg*t1);
  auto derH = (mK * meanZA * (t1 + 1. / bg2) - dedx) / (bg * t1);
  return derH + derH;
}

void NA6PTrackPar::g3helx3(float qfield, float step, std::array<float, 7>& vect)
{
  /******************************************************************
   *                                                                *
   *       GEANT3 tracking routine in a constant field oriented     *
   *       along axis 3                                             *
   *       Tracking is performed with a conventional                *
   *       helix step method                                        *
   *                                                                *
   *       Authors    R.Brun, M.Hansroul  *********                 *
   *       Rewritten  V.Perevoztchikov                              *
   *                                                                *
   *       Rewritten in C++ by I.Belikov                            *
   *                                                                *
   *  qfield (kG)       - particle charge times magnetic field      *
   *  step   (cm)       - step length along the helix               *
   *  vect[7](cm,GeV/c) - input/output x, y, z, px/p, py/p ,pz/p, p *
   *                                                                *
   ******************************************************************/
  const int ix = 0, iy = 1, iz = 2, ipx = 3, ipy = 4, ipz = 5, ipp = 6;
  constexpr float kOvSqSix = 0.408248f; // std::sqrt(1./6.);

  float cosx = vect[ipx], cosy = vect[ipy], cosz = vect[ipz];

  float rho = qfield * kB2C / vect[ipp];
  float tet = rho * step;

  float tsint, sintt, sint, cos1t;
  if (std::abs(tet) > 0.03f) {
    sint = std::sin(tet);
    sintt = sint / tet;
    tsint = (tet - sint) / tet;
    float t = std::sin(0.5f * tet);
    cos1t = 2 * t * t / tet;
  } else {
    tsint = tet * tet / 6.f;
    sintt = (1.f - tet * kOvSqSix) * (1.f + tet * kOvSqSix); // 1.- tsint;
    sint = tet * sintt;
    cos1t = 0.5f * tet;
  }

  float f1 = step * sintt;
  float f2 = step * cos1t;
  float f3 = step * tsint * cosz;
  float f4 = -tet * cos1t;
  float f5 = sint;

  vect[ix] += f1 * cosx - f2 * cosy;
  vect[iy] += f1 * cosy + f2 * cosx;
  vect[iz] += f1 * cosz + f3;

  vect[ipx] += f4 * cosx - f5 * cosy;
  vect[ipy] += f4 * cosy + f5 * cosx;
}

bool NA6PTrackPar::getZPCAToLine(float XL, float YL, float by, float& zca) const
{
  // Compute the Z of closest approach of the helical track to a vertical line (XL,YL)
  // in a constant By field. The problem is reduced to a 1D root in dz = z - z0,
  // solving F = s*(x-XL) + ty*(y-YL) = 0. The implementation uses numerically
  // stable expressions (no small differences like s-s0, cos differences, or
  // asin differences): x is evaluated via dz*(s0+s)/(c0+c) and the bending angle
  // via atan2. A straight-line approximation provides the initial guess, followed
  // by a few Newton iterations.
  constexpr double epsC = 1e-14, epsStep = 1.e-5, epsF = 1.e-12;
  constexpr int maxIter = 20;

  const double x0 = getX(), y0 = getY(), z0 = getZ(), s0 = getTx(), ty = getTy(), c0 = std::sqrt((1. - s0) * (1. + s0)), kappa = getCurvature(by);

  if (std::abs(kappa) < kSmallKappa) { // Straight-line case
    const double dx0 = x0 - XL, dy0 = y0 - YL, denom = s0 * s0 + ty * ty;
    zca = denom > epsC ? z0 - c0 * (s0 * dx0 + ty * dy0) / denom : z0;
    return true;
  }

  auto eval = [&](double dz, double& F, double& dF) -> bool { // for Newton iteration
    const double s = s0 + kappa * dz;
    if (abs(s) >= 1.0) {
      return false;
    }
    const double c = std::sqrt((1. - s) * (1. + s)), denomX = c0 + c;
    if (c < epsC || denomX < epsC) {
      return false;
    }
    const double x = x0 + dz * (s0 + s) / denomX;

    // Stable delta psi = atan2( sin(psi-psi0), cos(psi-psi0) )
    const double sinDP = s * c0 - s0 * c, cosDP = c * c0 + s * s0, dpsi = std::atan2(sinDP, cosDP), y = y0 + ty * dpsi / kappa;
    const double dx = x - XL, dy = y - YL;
    F = s * dx + ty * dy;
    dF = kappa * dx + (s * s + ty * ty) / c;
    return true;
  };

  double dz = 0.0; // Straight-line seed
  {
    const double dx0 = x0 - XL, dy0 = y0 - YL, denom = s0 * s0 + ty * ty;
    if (denom > epsC) {
      dz = -c0 * (s0 * dx0 + ty * dy0) / denom;
    }
  }
  double F = 0.0, dF = 0.0;
  bool ok = false;
  for (int it = 0; it < maxIter; ++it) {
    if (!eval(dz, F, dF)) {
      return false;
    }
    if (std::abs(F) < epsF) {
      ok = true;
      break;
    }
    if (std::abs(dF) < epsC) {
      break;
    }
    const double step = -F / dF;
    dz += step;
    if (std::abs(step) < epsStep) {
      ok = true;
      break;
    }
  }
  if (!ok) {
    if (!eval(dz, F, dF)) {
      return false;
    }
  }
  zca = z0 + dz;
  return true;
}

float NA6PTrackPar::getRapidity() const
{
  float en = getEnergy(), pz = getPz();
  return en > pz ? 0.5f * std::log((en + pz) / (en - pz)) : 999.f;
}
