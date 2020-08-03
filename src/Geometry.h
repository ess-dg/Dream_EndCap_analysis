

#pragma once
#include <analysis_cdt.h>

class Geometry {
public:
  Geometry();

  bool calculate(CDTReadout & readout, CDTGeometry & dg);

private:

  // used for the mapping of the electronics channels to detector voxels
  // integers > 0 correspond to the wire layer number given in the 2nd column
  // in the CDT file with channel mapping

  int s3[4][4] = {
      {0, 2, -1, -1}, {1, 3, -1, -1}, {-1, -1, 4, 6}, {-1, -1, 5, 7}}; // sumo3
  int s4[4][4] = {
      {-1, 2, 6, 10}, {-1, 3, 7, 11}, {0, 4, 8, -1}, {1, 5, 9, -1}}; // sumo4
  int s5[4][4] = {
      {0, 4, 8, 12}, {1, 5, 9, 13}, {2, 6, 10, 14}, {3, 7, 11, 15}}; // sumo5
  int s6[4][8] = {{0, -1, 4, 8, 12, 16, -1, -1},
                  {1, -1, 5, 9, 13, 17, -1, -1},
                  {-1, -1, 2, 6, 10, 14, -1, 18},
                  {-1, -1, 3, 7, 11, 15, -1, 19}}; // sumo6
};
