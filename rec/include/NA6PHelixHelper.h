// NA6PCCopyright

// Based on:
// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef _NA6P_HELIX_HELPER_
#define _NA6P_HELIX_HELPER_

#include <Rtypes.h>
#include "NA6PTrackPar.h"
#include "NA6PLine.h"

struct NA6PHelixHelper {
  using CircleXZ = std::array<float, 3>;
  enum { kR,
         kX,
         kZ };
  static constexpr float MaxDistXZDef = 10.f;
  float xDCA[2] = {};
  float zDCA[2] = {};
  int nDCA = 0;

  int circlesCrossInfo(const CircleXZ& trax0, const CircleXZ& trax1, float maxDistXZ = MaxDistXZDef, bool isCollinear = false)
  {
    const auto& trcA = trax0[kR] > trax1[kR] ? trax0 : trax1; // designate the largest circle as A
    const auto& trcB = trax0[kR] > trax1[kR] ? trax1 : trax0;
    nDCA = 0;
    float xDist = trcB[kX] - trcA[kX], zDist = trcB[kZ] - trcA[kZ];
    float dist2 = xDist * xDist + zDist * zDist, dist = std::sqrt(dist2), rsum = trcA[kR] + trcB[kR];
    if (dist < 1.e-12f) {
      return nDCA; // circles are concentric?
    }
    if (dist > rsum) { // circles don't touch, chose a point in between
      // the parametric equation of lines connecting the centers is
      // x = x0 + t/dist * (x1-x0), z = z0 + t/dist * (z1-z0)
      if (dist - rsum > maxDistXZ) { // too large distance
        return nDCA;
      }
      notTouchingXZ(dist, xDist, zDist, trcA, trcB[kR], isCollinear);
    } else if (auto dfr = dist + trcB[kR] - trcA[kR]; dfr < 0.) { // the small circle is nestled into large one w/o touching
      if (dfr > -maxDistXZ) {
        // select the point of closest approach of 2 circles
        notTouchingXZ(dist, xDist, zDist, trcA, -trcB[kR], isCollinear);
      } else {
        return nDCA;
      }
    } else { // 2 intersection points
      if (isCollinear) {
        /// collinear tracks, e.g. electrons from photon conversion
        /// if there are 2 crossings of the circle it is better to take
        /// a weighted average of the crossing points as a radius
        float r2r = trcA[kR] + trcB[kR];
        float r1_r = trcA[kR] / r2r;
        float r2_r = trcB[kR] / r2r;
        xDCA[0] = r2_r * trcA[kX] + r1_r * trcB[kX];
        zDCA[0] = r2_r * trcA[kZ] + r1_r * trcB[kZ];
        nDCA = 1;
      } else if (std::abs(zDist) < std::abs(xDist)) {
        // to simplify calculations, we move to new frame x->x+Xc0, z->z+Zc0, so that the 1st one is centered in origin
        float a = (trcA[kR] * trcA[kR] - trcB[kR] * trcB[kR] + dist2) / (2. * xDist), b = -zDist / xDist, ab = a * b, bb = b * b;
        float det = ab * ab - (1. + bb) * (a * a - trcA[kR] * trcA[kR]);
        if (det > 0.) {
          det = std::sqrt(det);
          zDCA[0] = (-ab + det) / (1. + b * b);
          xDCA[0] = a + b * zDCA[0] + trcA[kX];
          zDCA[0] += trcA[kZ];
          zDCA[1] = (-ab - det) / (1. + b * b);
          xDCA[1] = a + b * zDCA[1] + trcA[kX];
          zDCA[1] += trcA[kZ];
          nDCA = 2;
        } else { // due to the finite precision the det<=0, i.e. the circles are barely touching, fall back to this special case
          notTouchingXZ(dist, xDist, zDist, trcA, trcB[kR]);
        }
      } else {
        float a = (trcA[kR] * trcA[kR] - trcB[kR] * trcB[kR] + dist2) / (2. * zDist), b = -xDist / zDist, ab = a * b, bb = b * b;
        float det = ab * ab - (1. + bb) * (a * a - trcA[kR] * trcA[kR]);
        if (det > 0.) {
          det = std::sqrt(det);
          xDCA[0] = (-ab + det) / (1. + bb);
          zDCA[0] = a + b * xDCA[0] + trcA[kZ];
          xDCA[0] += trcA[kX];
          xDCA[1] = (-ab - det) / (1. + bb);
          zDCA[1] = a + b * xDCA[1] + trcA[kZ];
          xDCA[1] += trcA[kX];
          nDCA = 2;
        } else { // due to the finite precision the det<=0, i.e. the circles are barely touching, fall back to this special case
          notTouchingXZ(dist, xDist, zDist, trcA, trcB[kR]);
        }
      }
    }
    return nDCA;
  }

  void notTouchingXZ(float dist, float xDist, float zDist, const CircleXZ& trcA, float rBSign, bool isCollinear = false)
  {
    if (isCollinear) {
      /// for collinear tracks it is better to take
      /// a weighted average of the crossing points as a radius
      float r2r = trcA[kR] + rBSign;
      float r1_r = trcA[kR] / r2r;
      float r2_r = rBSign / r2r;
      xDCA[0] = r2_r * trcA[kX] + r1_r * (xDist + trcA[kX]);
      zDCA[0] = r2_r * trcA[kZ] + r1_r * (zDist + trcA[kZ]);
    } else {
      // fast method to calculate DCA between 2 circles, assuming that they don't touch each outer:
      // the parametric equation of lines connecting the centers is x = xA + t/dist * xDist, z = zA + t/dist * zDist
      // with xA,zA being the center of the circle A ( = trcA[kX], trcA[kZ] ), xDist = trcB[kX] = trcA[kX] ...
      // There are 2 special cases:
      // (a) small circle is inside the large one: provide rBSign as -trcB[kR]
      // (b) circle are side by side: provide rBSign as trcB[kR]
      auto t2d = (dist + trcA[kR] - rBSign) / dist;
      xDCA[0] = trcA[kX] + 0.5 * (xDist * t2d);
      zDCA[0] = trcA[kZ] + 0.5 * (zDist * t2d);
    }
    nDCA = 1;
  }

  int linesCrossInfo(const NA6PTrackPar& tr0, const NA6PTrackPar& tr1, float maxDistXZ = MaxDistXZDef)
  {
    /// closest approach of 2 straight lines
    /// exploit NA6PLine
    nDCA = 0;
    NA6PLine line0(tr0), line1(tr1);
    std::array<float, 3> p1, p2;
    bool isOk = line0.getClosestPoints(line1, p1, p2);
    if (!isOk) {
      return nDCA;
    }
    float dx = p1[0] - p2[0];
    float dz = p1[2] - p2[2];
    if (dx * dx + dz * dz > maxDistXZ * maxDistXZ) {
      return nDCA; // lines too far apart at closest approach
    }
    xDCA[0] = 0.5f * (p1[0] + p2[0]);
    zDCA[0] = 0.5f * (p1[2] + p2[2]);
    nDCA = 1;

    /* Strictly speaking, the code above does not minimize in XZ but rather in 3D, which is probably OK. For the record, the XZ version would be
       float x0 = tr0.getX(), z0 = tr0.getZ(), x1 = tr1.getX(), z1 = tr1.getZ();
       float dx0 = l0.mCosinesDirector[0], dz0 = l0.mCosinesDirector[2], dx1 = l1.mCosinesDirector[0], dz1 = l1.mCosinesDirector[2];
       float a00 = dx0 * dx0 + dz0 * dz0, a01 = -(dx0 * dx1 + dz0 * dz1), a11 = dx1 * dx1 + dz1 * dz1;
       float b0 = dx0 * (x1 - x0) + dz0 * (z1 - z0), b1 = dx1 * (x1 - x0) + dz1 * (z1 - z0), det = a00 * a11 - a01 * a01;
       if (std::abs(det) < 1e-12f) {
         return 0;
       }
       float t0 = ( b0 * a11 + a01 * b1) / det, t1 = ( a00 * b1 + a01 * b0) / det;
       float xp0 = x0 + dx0 * t0, zp0 = z0 + dz0 * t0, xp1 = x1 + dx1 * t1, zp1 = z1 + dz1 * t1;
       float d2 = (xp0 - xp1) * (xp0 - xp1) + (zp0 - zp1) * (zp0 - zp1);
       if (d2 > maxDistXZ * maxDistXZ) {
         return 0;
       }
       xDCA[0] = 0.5f * (xp0 + xp1);
       zDCA[0] = 0.5f * (zp0 + zp1);
       nDCA = 1;
     */
    return nDCA;
  }

  int circleLineCrossInfo(const CircleXZ& circ, const NA6PTrackPar& tr, float maxDistXZ = MaxDistXZDef)
  {
    /// closest approach of line and circle
    /// find intersection of a straight line with a circle in the XZ plane
    /// line: x = x0 + dx*t,  z = z0 + dz*t  (direction cosines dx, dz)
    /// circle: (x - xC)^2 + (z - zC)^2 = rC^2
    ///  (x0 - xC + dx*t)^2 + (z0 -zC + dz*t)^2 - rC^2 = 0
    ///  define: x0C = x0 - xC;  z0C = z0 - zC
    ///  x0C^2 + (dx*t)^2 + 2*x0C*dx*t + z0C^2 + (dz*t)^2 + 2*z0C*dz*t -rC^2 = 0
    ///  t^2 * (dx ^2 + dz^2) + t * (2*x0C*dx + 2*z0C*dz) + (x0C^2 + z0C^2 - rC^2) = 0
    ///  define: aA = (dx ^2 + dz^2); bB = x0C*dx + z0C*dz; cC = x0C^2 + z0C^2 - rC^2; det = bB^2 - aA*cC
    ///  t_1,2 = -bB / aA +- sqrt(det) / aA
    ///  xCross = x0 + t_1,2 * dx
    ///  zCross = z0 + t_1,2 * dz

    nDCA = 0;
    NA6PLine line(tr);
    float dx = line.mCosinesDirector[0], dz = line.mCosinesDirector[2];
    float x0C = line.mOriginPoint[0] - circ[kX], z0C = line.mOriginPoint[2] - circ[kZ];
    float aA = dx * dx + dz * dz, bB = (x0C * dx + z0C * dz), cC = (x0C * x0C + z0C * z0C - circ[kR] * circ[kR]);
    float det = bB * bB - aA * cC;
    if (det > 0.f) {
      float t1 = -bB / aA;
      float t2 = std::sqrt(det) / aA;
      int nCand = (t2 < 1e-6f * std::abs(t1)) ? 1 : 2;
      float t[2] = {t1 - t2, t1 + t2};
      for (int i = 0; i < nCand; i++) {
        xDCA[nDCA] = line.mOriginPoint[0] + t[i] * dx;
        zDCA[nDCA] = line.mOriginPoint[2] + t[i] * dz;
        nDCA++;
      }
    } else {
      // there is no crossing, find the point of the closest approach on the line which is closest to the circle center
      float tClosest = -bB / aA, xL = line.mOriginPoint[0] + tClosest * dx, zL = line.mOriginPoint[2] + tClosest * dz;
      float dxc = xL - circ[kX], dzc = zL - circ[kZ], dist = std::hypot(dxc, dzc);
      if (dist - circ[kR] > maxDistXZ) {
        return nDCA;
      }

      float drcf = circ[kR] / dist; // radius / distance to circle center
      float xH = circ[kX] + dxc * drcf, zH = circ[kZ] + dzc * drcf;
      xDCA[0] = 0.5f * (xL + xH);
      zDCA[0] = 0.5f * (zL + zH);
      nDCA = 1;
    }
    return nDCA;
  }

  int set(const CircleXZ& trax0, const NA6PTrack& tr0, const CircleXZ& trax1, const NA6PTrack& tr1, float maxDistXZ = MaxDistXZDef, bool isCollinear = false)
  {
    // calculate up to 2 crossings between 2 tracks
    nDCA = 0;
    bool curv0 = trax0[kR] > 0, curv1 = trax1[kR] > 0;

    if (curv0 && curv1) { // both are not straight lines
      nDCA = circlesCrossInfo(trax0, trax1, maxDistXZ, isCollinear);
    } else if (!curv0 && !curv1) { // both are straight lines
      nDCA = linesCrossInfo(tr0, tr1, maxDistXZ);
    } else {
      // curved track -> circle
      const auto& circ = curv0 ? trax0 : trax1;
      // not curved track -> line
      const auto& trLin = curv0 ? tr1 : tr0;
      nDCA = circleLineCrossInfo(circ, trLin, maxDistXZ);
    }
    //
    return nDCA;
  }

  NA6PHelixHelper() = default;

  NA6PHelixHelper(const CircleXZ& trax0, const NA6PTrack& tr0, const CircleXZ& trax1, const NA6PTrack& tr1, float maxDistXZ = MaxDistXZDef, bool isCollinear = false)
  {
    set(trax0, tr0, trax1, tr1, maxDistXZ, isCollinear);
  }

  ClassDefNV(NA6PHelixHelper, 1);
};

#endif
