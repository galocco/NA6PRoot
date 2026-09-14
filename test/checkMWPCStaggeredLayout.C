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
  const int nx[6] = {10, 5, 8, 8, 12, 12};
  const int ny[6] = {4, 4, 5, 5, 7, 7};
  const double expectedZ[4] = {-6., -2., 2., 6.};
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

    double gasX = -1., gasY = -1.;
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
      const double* chamberTr = chamberNode->GetMatrix()->GetTranslation();
      const double* gasTr = gasNode->GetMatrix()->GetTranslation();
      const double gx = chamberTr[0] + gasTr[0];
      const double gy = chamberTr[1] + gasTr[1];
      const double gz = chamberTr[2] + gasTr[2];

      if (std::abs(gz - expectedZ[q]) > tol) {
        printf("FAIL MS%d row %d col %d: gas z=%g, expected %g\n",
               st, row, col, gz, expectedZ[q]);
        ++errors;
      }

      auto* box = dynamic_cast<TGeoBBox*>(gasNode->GetVolume()->GetShape());
      if (!box) {
        printf("FAIL MS%d chamber %d: gas is not a box\n", st, ich);
        ++errors;
      } else {
        gasX = 2. * box->GetDX();
        gasY = 2. * box->GetDY();
      }

      if (ich < 4) {
        printf("MS%d chamber %3d row=%d col=%d q=%d gas=(%8.3f,%8.3f,%6.3f) cm\n",
               st, ich, row, col, q, gx, gy, gz);
      }
      ++total;
    }

    const double expectedGasX = st == 0 ? 25.6 : 51.2;
    if (std::abs(gasX - expectedGasX) > tol || std::abs(gasY - 66.56) > tol) {
      printf("FAIL MS%d: gas size %.6f x %.6f cm, expected %.6f x 66.56 cm\n",
             st, gasX, gasY, expectedGasX);
      ++errors;
    }

    printf("MS%d OK summary: grid=%dx%d N=%d gas=%.2fx%.2f cm\n",
           st, nx[st], ny[st], station->GetNdaughters(), gasX, gasY);
  }

  printf("\nChecked %d chambers. %s (%d errors)\n",
         total, errors == 0 ? "PASS" : "FAIL", errors);
}
