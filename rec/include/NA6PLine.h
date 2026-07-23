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

#ifndef NA6P_LINE_H
#define NA6P_LINE_H

#include <cmath>
#include <array>
#include <numeric>
#include <Rtypes.h>

class NA6PTrackPar; // fwd declaration only

// class for operations with straight lines (tracklets)

struct NA6PLine {
  NA6PLine() = default;
  NA6PLine(const NA6PLine& source) = default;
  NA6PLine& operator=(const NA6PLine& source) = default;
  NA6PLine(const NA6PTrackPar& trc);

  std::array<float, 3> evalAt(float t) const { return propagate(mOriginPoint, mCosinesDirector, t); }

  template <typename T = float>
  NA6PLine(const std::array<T, 3>& start, const std::array<T, 3>& end);

  template <typename T = float>
  static std::array<T, 3> propagate(const std::array<T, 3>& origin, const std::array<T, 3>& dirs, T t)
  {
    return getSum(origin, getScaled(dirs, t));
  }

  template <typename T = float>
  static NA6PLine fromTwoPoints(const std::array<T, 3>& p0, const std::array<T, 3>& p1)
  {
    return NA6PLine(p0, p1);
  }

  template <typename T = float>
  static NA6PLine fromPointAndDirection(const std::array<T, 3>& xyz, const std::array<T, 3>& dir);

  template <typename T1, typename T2>
  static std::array<T1, 3> getDiff(const std::array<T1, 3>& start, const std::array<T2, 3>& end)
  {
    return {end[0] - start[0], end[1] - start[1], end[2] - start[2]};
  }

  template <typename T1, typename T2>
  static std::array<T1, 3> getSum(const std::array<T1, 3>& start, const std::array<T2, 3>& end)
  {
    return {end[0] + start[0], end[1] + start[1], end[2] + start[2]};
  }

  template <typename T1, typename T2>
  static void add(std::array<T1, 3>& p, const std::array<T2, 3>& d)
  {
    for (int i = 0; i < 3; i++)
      p[i] += d[i];
  }

  template <typename T1, typename T2>
  static void subtract(std::array<T1, 3>& p, const std::array<T2, 3>& d)
  {
    for (int i = 0; i < 3; i++)
      p[i] -= d[i];
  }

  template <typename T>
  static T getNorm2(const std::array<T, 3>& v)
  {
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  }

  template <typename T>
  static void scale(std::array<T, 3>& v, T s)
  {
    for (int i = 0; i < 3; i++) {
      v[i] *= s;
    }
  }

  template <typename T>
  static std::array<T, 3> getScaled(const std::array<T, 3>& v, T s)
  {
    return {v[0] * s, v[1] * s, v[2] * s};
  }

  template <typename T>
  static T getNorm(const std::array<T, 3>& v)
  {
    return std::sqrt(getNorm2(v));
  }

  template <typename T>
  static std::array<T, 3> getCross(const std::array<T, 3>& v1, const std::array<T, 3>& v2)
  {
    return {v1[1] * v2[2] - v1[2] * v2[1], -v1[0] * v2[2] + v1[2] * v2[0], v1[0] * v2[1] - v1[1] * v2[0]};
  }

  template <typename T>
  static T getScalar(const std::array<T, 3>& v1, const std::array<T, 3>& v2)
  {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
  }

  template <typename T>
  float getDistanceFromPoint(const std::array<T, 3>& point) const;

  template <typename T>
  std::array<float, 6> getDCAComponents(const std::array<T, 3>& point) const;

  float getDCA(const NA6PLine&, float precision = 1e-7f) const;

  bool getClosestPoints(const NA6PLine& line2, std::array<float, 3>& p1, std::array<float, 3>& p2, float precision = 1e-7f) const;

  bool getCrossingPoint(const NA6PLine& line2, std::array<float, 3>& p, float precision = 1e-7f) const;

  bool areParallel(const NA6PLine&, float precision = 1e-7f) const;

  bool isSameLine(const NA6PLine& other, float precision) const { return !areParallel(other, precision) ? false : getDistanceFromPoint(other.mOriginPoint) <= precision; }

  bool isEmpty() const { return (mOriginPoint[0] == 0.f && mOriginPoint[1] == 0.f && mOriginPoint[2] == 0.f) && (mCosinesDirector[0] == 0.f && mCosinesDirector[1] == 0.f && mCosinesDirector[2] == 0.f); }

  std::array<float, 3> mOriginPoint{};
  std::array<float, 3> mCosinesDirector{};

  ClassDefNV(NA6PLine, 1);
};

template <typename T>
inline NA6PLine::NA6PLine(const std::array<T, 3>& start, const std::array<T, 3>& end)
{
  *this = fromPointAndDirection(start, getDiff(start, end));
}

template <typename T>
inline NA6PLine NA6PLine::fromPointAndDirection(const std::array<T, 3>& xyz, const std::array<T, 3>& dir)
{
  T norm = std::hypot(dir[0], dir[1], dir[2]), inv = (norm > T(1e-12)) ? T(1) / norm : T(0);
  NA6PLine l;
  l.mOriginPoint = xyz;
  l.mCosinesDirector = getScaled(dir, inv);
  return std::move(l);
}

template <typename T>
inline float NA6PLine::getDistanceFromPoint(const std::array<T, 3>& point) const
{
  const auto dif = getDiff(mOriginPoint, point);
  const float t = getScalar(dif, mCosinesDirector);
  return getNorm(getDiff(dif, getScaled(mCosinesDirector, t)));
}

template <typename T>
std::array<float, 6> NA6PLine::getDCAComponents(const std::array<T, 3>& point) const
{
  const auto dif = getDiff(mOriginPoint, point);
  const auto t = getScalar(dif, mCosinesDirector);

  std::array<float, 6> components{dif[0] - t * mCosinesDirector[0], 0., 0.,
                                  dif[1] - t * mCosinesDirector[1], 0.,
                                  dif[2] - t * mCosinesDirector[2]};
  const float cc0 = components[0] * components[0], cc3 = components[3] * components[3], cc5 = components[5] * components[5];
  components[1] = std::sqrt(cc0 + cc3);
  components[2] = std::sqrt(cc0 + cc5);
  components[4] = std::sqrt(cc3 + cc5);
  return components;
}

inline float NA6PLine::getDCA(const NA6PLine& secondLine, const float precision) const
{
  auto nxyz = getCross(mCosinesDirector, secondLine.mCosinesDirector);
  const float norm2 = getNorm2(nxyz);
  if (norm2 <= precision * precision) {
    return getDistanceFromPoint(secondLine.mOriginPoint);
  }
  const auto dif = getDiff(secondLine.mOriginPoint, mOriginPoint);
  return std::abs(getScalar(dif, nxyz)) / std::sqrt(norm2);
}

#endif // NA6P_LINE_H
