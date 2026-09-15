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

void checkMWPCWorkingAreaCoverage(const char* geometryFile = "geometry.root")
{
  // Updated working rectangles from the current MS layout study, detector X/Y.
  const double targetX[6] = {220., 230., 310., 320., 440., 500.};
  const double targetY[6] = {220., 240., 310., 320., 410., 440.};
  constexpr double tol = 1.e-4;

  TGeoManager::Import(geometryFile);
  if (!gGeoManager || !gGeoManager->GetTopVolume()) {
    printf("FAIL: could not import %s\n", geometryFile);
    return;
  }

  int errors = 0;
  int totalChambers = 0;
  auto* top = gGeoManager->GetTopVolume();

  printf("Station  grid      chambers  active coverage [cm]  target [cm]       excess [cm]      minimal\n");
  printf("------------------------------------------------------------------------------------------------\n");

  for (int st = 0; st < 6; ++st) {
    auto* stationNode = findStationNode(top, st);
    if (!stationNode) {
      printf("MS%d FAIL: station node not found\n", st);
      ++errors;
      continue;
    }

    auto* station = stationNode->GetVolume();
    totalChambers += station->GetNdaughters();
    double minX = 1.e30, maxX = -1.e30;
    double minY = 1.e30, maxY = -1.e30;
    std::vector<double> centresX;
    std::vector<double> centresY;

    for (int ich = 0; ich < station->GetNdaughters(); ++ich) {
      auto* chamberNode = station->GetNode(ich);
      auto* gasNode = chamberNode ? findGasNode(chamberNode->GetVolume()) : nullptr;
      if (!chamberNode || !gasNode) {
        printf("MS%d FAIL chamber %d: gas node missing\n", st, ich);
        ++errors;
        continue;
      }
      auto* box = dynamic_cast<TGeoBBox*>(gasNode->GetVolume()->GetShape());
      if (!box) {
        printf("MS%d FAIL chamber %d: gas is not TGeoBBox\n", st, ich);
        ++errors;
        continue;
      }

      const double gasOrigin[3] = {0., 0., 0.};
      double centreChamber[3] = {0., 0., 0.};
      double centreStation[3] = {0., 0., 0.};
      gasNode->GetMatrix()->LocalToMaster(gasOrigin, centreChamber);
      chamberNode->GetMatrix()->LocalToMaster(centreChamber, centreStation);
      centresX.push_back(centreStation[0]);
      centresY.push_back(centreStation[1]);

      for (int sx : {-1, +1}) {
        for (int sy : {-1, +1}) {
          const double cornerGas[3] = {sx * box->GetDX(), sy * box->GetDY(), 0.};
          double cornerChamber[3] = {0., 0., 0.};
          double cornerStation[3] = {0., 0., 0.};
          gasNode->GetMatrix()->LocalToMaster(cornerGas, cornerChamber);
          chamberNode->GetMatrix()->LocalToMaster(cornerChamber, cornerStation);
          minX = std::min(minX, cornerStation[0]);
          maxX = std::max(maxX, cornerStation[0]);
          minY = std::min(minY, cornerStation[1]);
          maxY = std::max(maxY, cornerStation[1]);
        }
      }
    }

    const auto xs = uniqueSorted(centresX);
    const auto ys = uniqueSorted(centresY);
    const int nx = static_cast<int>(xs.size());
    const int ny = static_cast<int>(ys.size());
    const double spanX = maxX - minX;
    const double spanY = maxY - minY;
    const double pitchX = regularPitch(xs);
    const double pitchY = regularPitch(ys);

    const bool coversX = spanX + tol >= targetX[st];
    const bool coversY = spanY + tol >= targetY[st];
    const bool minimalX = nx <= 1 || spanX - pitchX < targetX[st] - tol;
    const bool minimalY = ny <= 1 || spanY - pitchY < targetY[st] - tol;
    const bool pass = coversX && coversY && minimalX && minimalY;
    if (!pass) {
      ++errors;
    }

    printf("MS%d      %dx%-3d   %4d      %7.2f x %-7.2f  %6.1f x %-6.1f   %+6.2f x %+6.2f   %s\n",
           st, nx, ny, station->GetNdaughters(),
           spanX, spanY, targetX[st], targetY[st],
           spanX - targetX[st], spanY - targetY[st],
           pass ? "YES" : "NO");

    if (!coversX || !coversY) {
      printf("  FAIL MS%d: active footprint does not cover the required working rectangle\n", st);
    }
    if (!minimalX || !minimalY) {
      printf("  FAIL MS%d: grid is not minimal; a complete %s could be removed and still cover the target\n",
             st, !minimalX && !minimalY ? "column/row" : (!minimalX ? "column" : "row"));
    }
  }

  printf("------------------------------------------------------------------------------------------------\n");
  printf("Total chambers: %d\n", totalChambers);
  printf("%s: working-area coverage and minimal regular tiling (%d errors)\n",
         errors == 0 ? "PASS" : "FAIL", errors);
}
