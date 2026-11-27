#include "NA6PGenCocktail.h"
#include "NA6PGenParam.h"

// One muon per charge at the primary vertex, with independent flat pT and lab y.
NA6PGenerator* genMSPlugMuons(float ptMin = .15f, float ptMax = 1.2f, float yMin = 2.3f, float yMax = 4.6f)
{
  auto* generator = new NA6PGenCocktail("msPlugMuons");
  generator->setOwner(true);
  generator->addGenerator(new NA6PGenParam("plugMuPlus", -13, 1.f, "1", "1", ptMin, ptMax, yMin, yMax, true, true, false));
  generator->addGenerator(new NA6PGenParam("plugMuMinus", 13, 1.f, "1", "1", ptMin, ptMax, yMin, yMax, true, true, false));
  return generator;
}
