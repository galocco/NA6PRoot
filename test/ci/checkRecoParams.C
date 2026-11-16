R__LOAD_LIBRARY(librecLib)

#include "ConfigurableParam.h"
#include "NA6PRecoParam.h"
#include <TSystem.h>
#include <cmath>
#include <cstdlib>
#include <iostream>

void checkRecoParams()
{
  constexpr float expected = 41.25f;
  na6p::conf::ConfigurableParam::updateFromString("reco.mtZMatching=41.25");
  na6p::conf::ConfigurableParam::writeINI("na6pRecoParam.ini", "reco");

  na6p::conf::ConfigurableParam::updateFromString("reco.mtZMatching=17.0");
  na6p::conf::ConfigurableParam::updateFromFile("na6pRecoParam.ini", "reco", false);
  if (std::abs(NA6PRecoParam::Instance().mtZMatching - expected) > 1.e-5f) {
    std::cerr << "Reco parameter INI round-trip failed" << std::endl;
    gSystem->Exit(EXIT_FAILURE);
  }
}
