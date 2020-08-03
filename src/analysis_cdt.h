
#pragma once

const double pi = 3.14159265;
const double mn = 1.674927471e-27; // in kg
const double hp = 6.62607004e-34;  // in kg * m2/s
const double m2a = 1e+10;          // unit to convert m in A

// Moved constants

struct CDTReadout {
  int anode;
  int cathode;
  int subID;
  int sumo;
  int module;
  unsigned long int time;
  unsigned long int chopperTime;
  unsigned int boardID;
};

struct CDTGeometry {
  int nw_layer;
  int nstrip;
  int ncounter;
  int nsegment;
  int nwire;
  int nindexC;
  int nindexS;
};
