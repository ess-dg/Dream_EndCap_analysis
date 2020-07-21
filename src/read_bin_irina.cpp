/* ********************************************
 * binary file read for CDT SUMO detectors
 * date: 09.07.2019
 * attention: format could change
 * ********************************************/

#include <inttypes.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

typedef unsigned int uint32;
typedef unsigned long long uint64;

const uint64 kNeutronDataType = (0x0ULL) << 62;
const uint64 kMetaDataType = (0x2ULL) << 62;
const uint64 kAdcDataType = (0x3ULL) << 62;

const uint64 kMetaDataIndex0 = (0x0ULL) << 58;
const uint64 kMetaDataIndex1 = (0x1ULL) << 58;
const uint64 kMetaDataIndex2 = (0x2ULL) << 58;
const uint64 kMetaDataIndex3 = (0x3ULL) << 58;
const uint64 kMetaDataIndex4 = (0x4ULL) << 58;
const uint64 kMetaDataIndex5 = (0x5ULL) << 58;
const uint64 kMetaDataIndex6 = (0x6ULL) << 58;
const uint64 kMetaDataIndex7 = (0x7ULL) << 58;
const uint64 kMetaDataIndex8 = (0x8ULL) << 58;

const uint64 kMetaDataSubIndex0 = (0x0ULL) << 54;
const uint64 kMetaDataSubIndex1 = (0x1ULL) << 54;
const uint64 kMetaDataSubIndex2 = (0x2ULL) << 54;
const uint64 kMetaDataSubIndex3 = (0x3ULL) << 54;
const uint64 kMetaDataSubIndex4 = (0x4ULL) << 54;
const uint64 kMetaDataSubIndex5 = (0x5ULL) << 54;
const uint64 kMetaDataSubIndex6 = (0x6ULL) << 54;
const uint64 kMetaDataSubIndex7 = (0x7ULL) << 54;
const uint64 kMetaDataSubIndex8 = (0x8ULL) << 54;
const uint64 kMetaDataSubIndex9 = (0x9ULL) << 54;
const uint64 kMetaDataSubIndex10 = (0xAULL) << 54;

int main(int argc, char **argv) {

  //  remove("myfile_info.txt");

  if (remove("data_for_root.txt") != 0)
    puts("File not found, continue....");
  else
    puts("Output file data_for_root.txt succesfully deleted at the begining of "
         "the run, continue....");

  FILE *fp;
  FILE *ofile;
  int header[8];

  if (argc != 2) {
    printf("usage: read_bin filename\n");
    return -1;
  }
  const char *fn = argv[1];

  fp = fopen(fn, "rb");
  ofile = fopen("data_for_root.txt", "a");

  if (fp == NULL) {
    printf("ERROR: Can't open file %s\n", fn);
    return 0;
  }

  // read 32 Byte file header
  fread(header, 8 * sizeof(int), 1, fp);

  cout << "header...." << header << endl;

  // read rest of file
  while (!feof(fp)) {
    // readbuffer
    uint64 data = 0;

    // read and interpret data
    fread((char *)&data, sizeof(uint64), 1, fp);

    switch (data & ((0x3ULL) << 62)) {
    // NEUTRON
    case kNeutronDataType: {
      uint32 anode, cathode, subId;
      uint64 neutronTime = 0;

      anode = data & 0x3F; // read the last 6 bits

      cathode = (data & 0x1FC0) >>
                6; // read the last 7 but 6 bits and shift by 6 bits

      // skip the first 3 bits and read the next 16; shifts to the end by 13
      // positions

      neutronTime = (data & (0x1FFFFFFFFFFFE000ULL)) >> 13;

      // skip the first 3 bits and read the next 16; shifts to the end by 13
      // positions

      subId = (data & ((1ULL) << 58));

      if (ofile != NULL) {

        fprintf(ofile, "111   %llu    %d    %d     %d\n", neutronTime, cathode,
                anode, subId);
      }
      printf("NEUTRON DATA : time: 0x%016llX  cathode: 0x%0X  anode: 0x%0X "
             "subId; %d\n",
             neutronTime, cathode, anode, subId);
      //				printf("NEUTRON DATA : time: %llu
      // cathode: %d  anode: %d subId; %d\n",neutronTime,cathode,anode,subId);
      printf(" \n");
    } break;
    // METADATA
    case kMetaDataType:
      switch (data & ((0xFULL) << 58)) // this is 0, so this is a neutron event,
                                       // see definition at the beginning
      {
      case kMetaDataIndex0:
        break;
      case kMetaDataIndex1:
        switch (data & (uint64(0xF) << 54)) {
        case kMetaDataSubIndex0: {
          uint64 dataTime = 0;
          dataTime = (data & (uint64(0xFFFFFFFFFFFF)));

          if (ofile != NULL) {

            fprintf(ofile, "11  %llu   -10   -10   -10   \n", dataTime);
          }

          printf("CHOPPER TIMESTAMP : time: 0x%" PRIx64 "08X\n", dataTime);
        } break;
        case kMetaDataSubIndex1: {
          uint32 lastBoard = (data & 0xffffff);

          if (ofile != NULL) {

            fprintf(ofile, "1   %u    -1   -1   -1  \n", lastBoard);
          }

          printf("BOARD-ID: id: %u\n", lastBoard);
        } break;
        default:
          break;
        }
        break;
      case kMetaDataIndex2:
        switch (data & (uint64(0xF) << 54)) {
        case kMetaDataSubIndex0:
        case kMetaDataSubIndex1:
        case kMetaDataSubIndex2:
        case kMetaDataSubIndex3:
        case kMetaDataSubIndex4:
        case kMetaDataSubIndex5:
        case kMetaDataSubIndex6:
        case kMetaDataSubIndex7:
        case kMetaDataSubIndex8:
        case kMetaDataSubIndex9:
        case kMetaDataSubIndex10:
        default:
          break;
        }
        break;
      case kMetaDataIndex3:
        switch (data & (uint64(0xF) << 54)) {
        case kMetaDataSubIndex0:
        case kMetaDataSubIndex1:
        case kMetaDataSubIndex2:
        case kMetaDataSubIndex3:
        case kMetaDataSubIndex4:
        case kMetaDataSubIndex5:
        case kMetaDataSubIndex6:
        case kMetaDataSubIndex7:
        case kMetaDataSubIndex8:
        case kMetaDataSubIndex9:
        case kMetaDataSubIndex10:
        default:
          break;
        }
        break;
      case kMetaDataIndex4:
        switch (data & (uint64(0xF) << 54)) {
        case kMetaDataSubIndex0:
        case kMetaDataSubIndex1:
        case kMetaDataSubIndex2:
        default:
          break;
        }
        break;
      case kMetaDataIndex5:
        switch (data & (uint64(0xF) << 54)) {
        case kMetaDataSubIndex0:
        case kMetaDataSubIndex1:
        case kMetaDataSubIndex2:
        case kMetaDataSubIndex3:
        case kMetaDataSubIndex4:
        case kMetaDataSubIndex5:
        case kMetaDataSubIndex6:
        default:
          break;
        }
        break;
      case kMetaDataIndex6:
        switch (data & (uint64(0xF) << 54)) {
        case kMetaDataSubIndex0:
          break;
        }
        break;
      case kMetaDataIndex7:
        switch (data & (uint64(0xF) << 54)) {
        case kMetaDataSubIndex0:
        case kMetaDataSubIndex1:
        case kMetaDataSubIndex2:
        case kMetaDataSubIndex3:
        default:
          break;
        }
        break;
      case kMetaDataIndex8:
        switch (data & (uint64(0xF) << 54)) {
        case kMetaDataSubIndex0:
          break;
        default:
          break;
        }
        break;
      default:
        break;
      }
      break;
      // ADC DATA
    case kAdcDataType:
      break;
    default:
      break;
    }
  }

  fclose(fp);
  fclose(ofile);

  return 1;
}
