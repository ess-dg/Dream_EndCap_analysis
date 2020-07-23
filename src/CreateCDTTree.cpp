/* Copyright (C) 2020 European Spallation Source, ERIC. See LICENSE file      */
//===----------------------------------------------------------------------===//
///
/// \file CreateCDTTree.cpp
///
/// Usage: .L CreateCDTTree.cpp /CreateCDTRootFile()
///
//===----------------------------------------------------------------------===//

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

void CreateCDTRootFile() {
  TTree *t = new TTree("cdt_ev", "cdt_ev");

  Int_t cathode, anode, subID, module, sumo;
  UInt_t boardID;
  ULong64_t chopperTime, neutronTime;
  int index;

  // see below the description of the variables

  t->Branch("neutronTime", &neutronTime, "neutronTime/l");
  t->Branch("cathode", &cathode, "cathode/I");
  t->Branch("anode", &anode, "anode/I");
  t->Branch("subID", &subID, "subID/I");
  t->Branch("boardID", &boardID, "boardID/i");
  t->Branch("module", &module, "module/I");
  t->Branch("sumo", &sumo, "sumo/I");
  t->Branch("chopperTime", &chopperTime, "chopperTime/l");

  gSystem->Exec("rm cdt.root");
  TFile *rf = new TFile("cdt.root", "new");

  char hname[20];

  Char_t a[5], b[15], c[5], d[5], e[5];

  sprintf(hname, "data_for_root.txt");
  FILE *fp;
  printf("Data file: %s\n", hname);
  fp = fopen(hname, "r");

  if (fp == NULL) {
    printf("Error opening the input file!");
    exit(1);
  }

  while (!feof(fp)) {
    fscanf(fp, "%s %s %s %s %s", a, b, c, d, e);

    //        std::cout<<"a = "<<a<<", b = "<<b<<", c = "<<c<<", d = "<<d<<", e
    //        = "<<e<<std::endl;

    index = strtod(a, NULL);

    switch (index) {
      // this is a neutron event
    case 111:

      neutronTime = strtol(b, NULL, 10);
      cathode = strtod(c, NULL);
      anode = strtod(d, NULL);
      subID = strtod(e, NULL);
      chopperTime = 0;
      boardID = 0;
      module = 1;
      sumo = 0;

      break;
      // this is a chopper event
    case 11:

      chopperTime = strtol(b, NULL, 10);
      neutronTime = 0;
      cathode = 0;
      anode = 0;
      subID = 0;
      boardID = 0;
      module = 1;
      sumo = 0;

      break;
      // this is the board ID
    case 1:

      boardID = strtol(b, NULL, 10);
      chopperTime = 0;
      neutronTime = 0;
      cathode = 0;
      anode = 0;
      subID = 0;
      module = 1;
      if (boardID == 1418045) {
        sumo = 3;
      } else if (boardID == 1416964) {
        sumo = 4;
      } else if (boardID == 1416799) {
        sumo = 5;
      } else if (boardID == 1416697) {
        sumo = 6;
      }

      break;
    }

    t->Fill();
  }

  //	else {printf("Couldn't open data file %s!\n",hname);}

  fclose(fp);
  t->Write();
  std::cout << "Fill data tree...done" << std::endl;
  std::cout << "The data tree has " << t->GetEntries() << " entries."
            << std::endl;
  std::cout << "  " << std::endl;

  rf->Close();
}

int main() {
  CreateCDTRootFile();
  return 0;
}
