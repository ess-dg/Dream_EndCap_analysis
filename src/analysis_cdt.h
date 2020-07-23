
#pragma once

const double pi = 3.14159265;
const double mn = 1.674927471e-27; // in kg
const double hp = 6.62607004e-34;  // in kg * m2/s
const double m2a = 1e+10;          // unit to convert m in A

// Moved constants
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

// time shifts and delays for the tof correction for events collected in WFM
// mode

signed long int tdelay = 0;
signed long int tshift[6] = {0, 2600, 4400, 6500, 9270, 11700};
signed long int tmin[6] = {9600, 21940, 33000, 42540, 51670, 60300};
signed long int tmax[6] = {20680, 31200, 40420, 48850, 57800, 70000};
