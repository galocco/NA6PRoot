#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

namespace
{
TGeoNode* findStationNode(TGeoVolume* top, int station)
{
  const TString wanted = Form("MS%d", station);
  for (int i = 0; i < top->GetNdaughters(); ++i) {
    auto* node = top->GetNode(i);
    if (node && node->GetVolume() && wanted == node->GetVolume()->GetName()) {
      return node;
    }
  }
  return nullptr;
}

TGeoNode* findGasNode(TGeoVolume* chamber)
{
  for (int i = 0; i < chamber->GetNdaughters(); ++i) {
    auto* node = chamber->GetNode(i);
    if (node && node->GetVolume() && std::strstr(node->GetVolume()->GetName(), "_Gas")) {
      return node;
    }
  }
  return nullptr;
}

std::vector<double> uniqueSorted(std::vector<double> values, double tol = 1.e-4)
{
  std::sort(values.begin(), values.end());
  std::vector<double> out;
  for (double v : values) {
    if (out.empty() || std::abs(v - out.back()) > tol) {
      out.push_back(v);
    }
  }
  return out;
}

double regularPitch(const std::vector<double>& centres)
{
  if (centres.size() < 2) {
    return 0.;
  }
  double pitch = 1.e30;
  for (std::size_t i = 1; i < centres.size(); ++i) {
    pitch = std::min(pitch, centres[i] - centres[i - 1]);
  }
  return pitch;
}
}

void checkMWPCStaggeredLayout(const char* geometryFile = "geometry.root")
{
  // Nx and Ny are inferred from the generated geometry so this checker remains
  // valid when the station grids are changed during later efficiency studies.
  const double expectedZ[4] = {-6., -2., 2., 6.};
  constexpr double overlapX = 3.;
  constexpr double overlapY = 3.;
  constexpr double tol = 1.e-4;

  TGeoManager::Import(geometryFile);
  if (!gGeoManager || !gGeoManager->GetTopVolume()) {
    printf("FAIL: could not import %s\n", geometryFile);
    return;
  }

  int errors = 0;
  int total = 0;
  auto* top = gGeoManager->GetTopVolume();
  for (int st = 0; st < 6; ++st) {
    auto* stationNode = findStationNode(top, st);
    if (!stationNode) {
      printf("FAIL MS%d: station node not found\n", st);
      ++errors;
      continue;
    }
    auto* station = stationNode->GetVolume();

    const double expectedGasX = 66.56;
    const double expectedGasY = st == 0 ? 30.72 : 51.2;
    const double expectedPitchX = expectedGasX - overlapX;
    const double expectedPitchY = expectedGasY - overlapY;

    std::vector<double> centresX;
    std::vector<double> centresY;
    double gasGlobalX = -1., gasGlobalY = -1.;

    // First pass: infer the detector-coordinate grid from the actual sensitive
    // gas centres and measure the installed gas size after chamber rotation.
    for (int ich = 0; ich < station->GetNdaughters(); ++ich) {
      auto* chamberNode = station->GetNode(ich);
      auto* gasNode = chamberNode ? findGasNode(chamberNode->GetVolume()) : nullptr;
      if (!chamberNode || !gasNode) {
        printf("FAIL MS%d chamber %d: gas node missing\n", st, ich);
        ++errors;
        continue;
      }

      const double* gasLocal = gasNode->GetMatrix()->GetTranslation();
      double gasStation[3] = {0., 0., 0.};
      chamberNode->GetMatrix()->LocalToMaster(gasLocal, gasStation);
      centresX.push_back(gasStation[0]);
      centresY.push_back(gasStation[1]);

      auto* box = dynamic_cast<TGeoBBox*>(gasNode->GetVolume()->GetShape());
      if (!box) {
        printf("FAIL MS%d chamber %d: gas is not a box\n", st, ich);
        ++errors;
      } else {
        const double vxLocal[3] = {box->GetDX(), 0., 0.};
        const double vyLocal[3] = {0., box->GetDY(), 0.};
        double vxGlobal[3] = {0., 0., 0.};
        double vyGlobal[3] = {0., 0., 0.};
        chamberNode->GetMatrix()->LocalToMasterVect(vxLocal, vxGlobal);
        chamberNode->GetMatrix()->LocalToMasterVect(vyLocal, vyGlobal);
        gasGlobalX = 2. * (std::abs(vxGlobal[0]) + std::abs(vyGlobal[0]));
        gasGlobalY = 2. * (std::abs(vxGlobal[1]) + std::abs(vyGlobal[1]));
      }
    }

    const auto xs = uniqueSorted(centresX);
    const auto ys = uniqueSorted(centresY);
    const int nx = static_cast<int>(xs.size());
    const int ny = static_cast<int>(ys.size());
    const double pitchX = regularPitch(xs);
    const double pitchY = regularPitch(ys);

    if (nx <= 0 || ny <= 0 || station->GetNdaughters() != nx * ny) {
      printf("FAIL MS%d: inferred grid %dx%d is inconsistent with %d chambers\n",
             st, nx, ny, station->GetNdaughters());
      ++errors;
      continue;
    }
    if (nx > 1 && std::abs(pitchX - expectedPitchX) > tol) {
      printf("FAIL MS%d: X pitch %.6f cm, expected %.6f cm\n", st, pitchX, expectedPitchX);
      ++errors;
    }
    if (ny > 1 && std::abs(pitchY - expectedPitchY) > tol) {
      printf("FAIL MS%d: Y pitch %.6f cm, expected %.6f cm\n", st, pitchY, expectedPitchY);
      ++errors;
    }

    // Second pass: node ordering must still correspond to row-major (row,col),
    // and the A/B/C/D parity must have the correct z level.
    for (int ich = 0; ich < station->GetNdaughters(); ++ich) {
      auto* chamberNode = station->GetNode(ich);
      auto* gasNode = chamberNode ? findGasNode(chamberNode->GetVolume()) : nullptr;
      if (!chamberNode || !gasNode) {
        continue;
      }

      const int row = ich / nx;
      const int col = ich % nx;
      const int q = 2 * (row & 1) + (col & 1);
      const double* gasLocal = gasNode->GetMatrix()->GetTranslation();
      double gasStation[3] = {0., 0., 0.};
      chamberNode->GetMatrix()->LocalToMaster(gasLocal, gasStation);

      const double expectedX = (static_cast<double>(col) - 0.5 * (nx - 1)) * expectedPitchX;
      const double expectedY = (static_cast<double>(row) - 0.5 * (ny - 1)) * expectedPitchY;
      if (std::abs(gasStation[0] - expectedX) > tol ||
          std::abs(gasStation[1] - expectedY) > tol ||
          std::abs(gasStation[2] - expectedZ[q]) > tol) {
        printf("FAIL MS%d row %d col %d: gas=(%g,%g,%g), expected=(%g,%g,%g)\n",
               st, row, col,
               gasStation[0], gasStation[1], gasStation[2],
               expectedX, expectedY, expectedZ[q]);
        ++errors;
      }

      if (ich < 4) {
        printf("MS%d chamber %3d row=%d col=%d q=%d gas=(%8.3f,%8.3f,%6.3f) cm\n",
               st, ich, row, col, q,
               gasStation[0], gasStation[1], gasStation[2]);
      }
      ++total;
    }

    if (std::abs(gasGlobalX - expectedGasX) > tol ||
        std::abs(gasGlobalY - expectedGasY) > tol) {
      printf("FAIL MS%d: installed gas size %.6f x %.6f cm, expected %.6f x %.6f cm\n",
             st, gasGlobalX, gasGlobalY, expectedGasX, expectedGasY);
      ++errors;
    }

    const double coverageX = expectedGasX + (nx - 1) * expectedPitchX;
    const double coverageY = expectedGasY + (ny - 1) * expectedPitchY;
    printf("MS%d OK summary: inferred grid=%dx%d N=%d detectorGas=%.2fx%.2f cm coverage=%.2fx%.2f cm\n",
           st, nx, ny, station->GetNdaughters(),
           gasGlobalX, gasGlobalY, coverageX, coverageY);
  }

  printf("\nChecked %d chambers. %s (%d errors)\n",
         total, errors == 0 ? "PASS" : "FAIL", errors);
}
