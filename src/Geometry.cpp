
#include <Geometry.h>
#include <math.h>


Geometry::Geometry() {};


bool Geometry::calculate(CDTReadout & readout, CDTGeometry & dg) {

  int factor_a = floor(readout.cathode / 16);
  int factor_b = floor(readout.anode / 16);

  if (factor_b > 3) {
    return false;
  }

  if (factor_a > 3  or (readout.sumo == 6 and factor_a > 7 )) {
    return false;
  }

  switch (readout.sumo) {

  case 3:
    dg.nw_layer = s3[factor_b][factor_a];
    break;

  case 4:
    dg.nw_layer = s4[factor_b][factor_a];
    break;

  case 5:
    dg.nw_layer = s5[factor_b][factor_a];
    break;

  case 6:
    dg.nw_layer = s6[factor_b][factor_a];
    break;

  default:  // All other values are invalid
    return false;
    break;
  }

  dg.nsegment = floor(dg.nw_layer / 2) + 1;
  dg.ncounter = dg.nw_layer % 2 + 1;

  ///\todo replace hardcoded value 16 with const variable
  /// \todo seems wrong always gives 1 ??
  dg.nwire = readout.anode - 16 * factor_b + 1;
  dg.nstrip = readout.cathode - 16 * factor_a + 1;

  return true;
}
