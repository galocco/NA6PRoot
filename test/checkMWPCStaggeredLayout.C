#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>

#include <cmath>
#include <cstdio>
#include <cstring>

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
}

void checkMWPCStaggeredLayout(const char* geometryFile = "geometry.root")
{
  // Detector-coordinate minimal regular grids for the corrected orientation.
  const int nx[6] = {4, 4, 5, 5, 7, 8};
  const int ny[6] = {8, 5, 7, 7, 9, 10};
  const double expectedZ[4] = {-6., -2., 2., 6.};
  constexpr double overlapX = 3.;
  constexpr double overlapY = 3.;
  constexpr double tol = 1.e-5;

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
    const int expectedN = nx[st] * ny[st];
    if (station->GetNdaughters() != expectedN) {
      printf("FAIL MS%d: expected %d chambers, found %d\n",
             st, expectedN, station->GetNdaughters());
      ++errors;
    }

    const double expectedGasX = 66.56;
    const double expectedGasY = st == 0 ? 30.72 : 51.2;
    const double pitchX = expectedGasX - overlapX;
    const double pitchY = expectedGasY - overlapY;
    double gasGlobalX = -1., gasGlobalY = -1.;

    for (int ich = 0; ich < station->GetNdaughters(); ++ich) {
      auto* chamberNode = station->GetNode(ich);
      auto* gasNode = chamberNode ? findGasNode(chamberNode->GetVolume()) : nullptr;
      if (!chamberNode || !gasNode) {
        printf("FAIL MS%d chamber %d: gas node missing\n", st, ich);
        ++errors;
        continue;
      }

      const int row = ich / nx[st];
      const int col = ich % nx[st];
      const int q = 2 * (row & 1) + (col & 1);

      // Gas-node translation is expressed in chamber-local coordinates. Apply
      // the chamber placement matrix (including the 90-degree z rotation) to
      // obtain detector/station coordinates.
      const double* gasLocal = gasNode->GetMatrix()->GetTranslation();
      double gasStation[3] = {0., 0., 0.};
      chamberNode->GetMatrix()->LocalToMaster(gasLocal, gasStation);

      const double expectedX = (static_cast<double>(col) - 0.5 * (nx[st] - 1)) * pitchX;
      const double expectedY = (static_cast<double>(row) - 0.5 * (ny[st] - 1)) * pitchY;
      if (std::abs(gasStation[0] - expectedX) > tol ||
          std::abs(gasStation[1] - expectedY) > tol ||
          std::abs(gasStation[2] - expectedZ[q]) > tol) {
        printf("FAIL MS%d row %d col %d: gas=(%g,%g,%g), expected=(%g,%g,%g)\n",
               st, row, col,
               gasStation[0], gasStation[1], gasStation[2],
               expectedX, expectedY, expectedZ[q]);
        ++errors;
      }

      auto* box = dynamic_cast<TGeoBBox*>(gasNode->GetVolume()->GetShape());
      if (!box) {
        printf("FAIL MS%d chamber %d: gas is not a box\n", st, ich);
        ++errors;
      } else {
        // Project the chamber-local half-axis vectors into station coordinates
        // to measure the installed sensitive footprint after rotation.
        const double vxLocal[3] = {box->GetDX(), 0., 0.};
        const double vyLocal[3] = {0., box->GetDY(), 0.};
        double vxGlobal[3] = {0., 0., 0.};
        double vyGlobal[3] = {0., 0., 0.};
        chamberNode->GetMatrix()->LocalToMasterVect(vxLocal, vxGlobal);
        chamberNode->GetMatrix()->LocalToMasterVect(vyLocal, vyGlobal);
        gasGlobalX = 2. * (std::abs(vxGlobal[0]) + std::abs(vyGlobal[0]));
        gasGlobalY = 2. * (std::abs(vxGlobal[1]) + std::abs(vyGlobal[1]));
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

    const double coverageX = expectedGasX + (nx[st] - 1) * pitchX;
    const double coverageY = expectedGasY + (ny[st] - 1) * pitchY;
    printf("MS%d OK summary: grid=%dx%d N=%d detectorGas=%.2fx%.2f cm coverage=%.2fx%.2f cm\n",
           st, nx[st], ny[st], station->GetNdaughters(),
           gasGlobalX, gasGlobalY, coverageX, coverageY);
  }

  printf("\nChecked %d chambers. %s (%d errors)\n",
         total, errors == 0 ? "PASS" : "FAIL", errors);
}
