#include "TF1.h"
#include "TH1.h"
#include "TH3.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLeaf.h"
#include "TStyle.h"
#include "TRandom.h"
#include "stdio.h"
#include "stdlib.h"
#include </Library/Developer/CommandLineTools/usr/include/c++/v1/iostream>

//Usage: .L analysis_heimdal.C / analysis() 
// to be used with the cdt.root file created by using CreateCDTtree.C
//***********************************
// Written by Irina Stefanescu from ESS DG to analyse the data collected with the 
// DREAM EndCap detector in July 2019 at V20-HZB, may the beamline rest in peace. 
//****************************

const double pi = 3.14159265; 
const double mn = 1.674927471e-27;  // in kg
const double hp = 6.62607004e-34;   // in kg * m2/s
const double m2a = 1e+10;    // unit to convert m in A

void analysis()
{

TFile* file = TFile::Open("cdt.root","update");

TTree* tree = (TTree*)file->Get("cdt_ev"); 

const int nevents = tree->GetEntries();

int* vanode=new int[nevents]();
int* vcathode=new int[nevents]();
int* vsubID=new int[nevents]();
int* vsumo=new int[nevents]();
int* vmodule=new int[nevents]();
unsigned long int* vtime=new unsigned long int[nevents]();
unsigned long int* vchopperTime=new unsigned long int[nevents]();
unsigned int* vboardID=new unsigned int[nevents]();

Int_t seg,scopy,strip,voxel,coat;

Int_t n_wires = 16;
Int_t n_strips = 16;   

Double_t energy,posx,posy,posz;
Int_t cathode,anode,subID;
Int_t module, sumo;
UInt_t boardID;
ULong64_t chopperTime,neutronTime;
Int_t mult,index;

//Tree containing the raw events (created at the previous step by running CreateCDTTree.C)

tree->SetBranchAddress("neutronTime",&neutronTime);
tree->SetBranchAddress("boardID",&boardID);
tree->SetBranchAddress("chopperTime",&chopperTime);
tree->SetBranchAddress("sumo",&sumo);
tree->SetBranchAddress("module",&module);
tree->SetBranchAddress("subID",&subID);
tree->SetBranchAddress("anode",&anode);
tree->SetBranchAddress("cathode",&cathode);

// New event tree for the 'calibrated events'

if ((TTree*) file->Get("cdt_new"))
	{std::cout<<"Found an old tree, I am going to delete it!"<<std::endl; 
	 gDirectory->Delete("cdt_new");
	 std::cout<<"....Done!"<<std::endl; std::cout<<" "<<std::endl; }
	
std::cout<<"now I am going to create the new tree!"<<std::endl; 
TTree* newt=new TTree("cdt_new","cdt_new");
std::cout<<"....Done!"<<std::endl; 
std::cout<<" "<<std::endl; 
     
int ncathode,nanode,nmult,nsubID,nevent;
int nmodule,nsumo,nw_layer,nstrip;
int ncounter,nsegment,factor_a,factor_b,nwire;
int nindexC,nindexS,ccycle;

signed int nboardID;
signed long int nchopperTime,ntime;
signed long int ntof_wfm,tdiff_wfm,ntof,tdiff,ntof_wfm_corr;
     

/*    supported types in ROOT (reminder)

    C : a character string terminated by the 0 character
    B : an 8 bit signed integer (Char_t)
    b : an 8 bit unsigned integer (UChar_t)
    S : a 16 bit signed integer (Short_t)
    s : a 16 bit unsigned integer (UShort_t)
    I : a 32 bit signed integer (Int_t)
    i : a 32 bit unsigned integer (UInt_t)
    F : a 32 bit floating point (Float_t)
    f : a 24 bit floating point with truncated mantissa (Float16_t)
    D : a 64 bit floating point (Double_t)
    d : a 24 bit truncated floating point (Double32_t)
    L : a 64 bit signed integer (Long64_t)
    l : a 64 bit unsigned integer (ULong64_t)
    O : [the letter o, not a zero] a boolean (Bool_t)
    
*/


newt->Branch("ntime",&ntime,"ntime/L");   //neutron time stamp
newt->Branch("nchopperTime",&nchopperTime,"nchopperTime/L");  //chopper timestamp
newt->Branch("ntof",&ntof,"ntof/L");     //neutron time of flight
newt->Branch("ntof_wfm",&ntof_wfm,"ntof_wfm/L");    // corrected tof; divide by 10000 to get it in ms
newt->Branch("ntof_wfm_corr",&ntof_wfm_corr,"ntof_wfm_corr/L");    // divide by 10000 to get it in ms
newt->Branch("tdiff",&tdiff,"tdiff/L");             // divide by 10000 to get it in ms
newt->Branch("tdiff_wfm",&tdiff_wfm,"tdiff_wfm/L");  // divide by 10000 to get it in ms
newt->Branch("ncathode",&ncathode,"ncathode/I");   //cathode plane number 
newt->Branch("nanode",&nanode,"nanode/I");   //anode plane number
newt->Branch("nsubID",&nsubID,"nsubID/I");  // nsubID, not used
newt->Branch("nsumo",&nsumo,"nsumo/I");    // sumo number
newt->Branch("nmodule",&nmodule,"nmodule/I");   // module number (12-deg sector)
newt->Branch("nmult",&nmult,"nmult/I");  //multiplicity event, not used
newt->Branch("nsegment",&nsegment,"nsegment/I");   //segment number
newt->Branch("ncounter",&ncounter,"ncounter/I");   //counter number (1 for counter left of common cathode, 2 for right)
newt->Branch("nwire",&nwire,"nwire/I");   // wire number 
newt->Branch("nstrip",&nstrip,"nstrip/I");   //strip number 
newt->Branch("nevent",&nevent,"nevent/I");    //event counter
newt->Branch("nboardID",&nboardID,"nboardID/i");   // board no => sumo number 
newt->Branch("nw_layer",&nw_layer,"nw_layer/I");  // wire layer#
newt->Branch("nindexS",&nindexS,"nindexS/I");   // index1 for mapping
newt->Branch("nindexC",&nindexC,"nindexC/I");    //index2 for mapping
newt->Branch("ccycle",&ccycle,"ccycle/I");    // chopper cycle

// tree for mapping the real events with the voxel positions from GEANT4
 
TFile* ffile = TFile::Open("cdt_new_cal.root","recreate");

if ((TTree*) ffile->Get("cdt_new_cal"))
	{std::cout<<"Found an old cal tree, I am going to delete it!"<<std::endl; 
	 gDirectory->Delete("cdt_new_cal");
	 std::cout<<"....Done!"<<std::endl; std::cout<<" "<<std::endl; }
	
std::cout<<"now I am going to create the new calibration tree!"<<std::endl; 
TTree* fnewt=new TTree("cdt_new_cal","cdt_new_cal");
std::cout<<"....Done!"<<std::endl; 
std::cout<<" "<<std::endl; 

double rad_c,lambda_c,lambda_c_wfm_corr,twotheta_c,phi_c,ntof_c;
double nvoxel_x,nvoxel_y,nvoxel_z,nangle;
int csumo,nc,nikS,nikC;
float dspacing_c_wfm,dspacing_c;

fnewt->Branch("nikS",&nikS,"nikS/I");  //index1 used for mapping events
fnewt->Branch("nikC",&nikC,"nikC/I");   //index2 used for mapping events
fnewt->Branch("nc",&nc,"nc/I");     //counter number
fnewt->Branch("csumo",&csumo,"csumo/I");   //sumo number
fnewt->Branch("nvoxel_x",&nvoxel_x,"nvoxel_x/D");  //GEANT4 x-pos of voxel centre
fnewt->Branch("nvoxel_y",&nvoxel_y,"nvoxel_y/D");   //GEANT4 y-pos of voxel centre
fnewt->Branch("nvoxel_z",&nvoxel_z,"nvoxel_z/D");   //GEANT4 z-pos of voxel centre
fnewt->Branch("twotheta_c",&twotheta_c,"twotheta_c/D");  //theta-angle of the GEANT4 voxel 
fnewt->Branch("nangle",&nangle,"nangle/D");  // inclination angle of the Boron layer 
fnewt->Branch("phi_c",&phi_c,"phi_c/D");   //phi-angle of the GEANT4 voxel 
fnewt->Branch("rad_c",&rad_c,"rad_c/D");   // GEANT4 distance of the voxel centre to the sample
fnewt->Branch("ntof_c",&ntof_c,"ntof_c/D");   //defined, but not used
fnewt->Branch("lambda_c",&lambda_c,"lambda_c/D");  //lambda calculated by using the G4 position of the voxel (normal operation mode)
fnewt->Branch("dspacing_c",&dspacing_c,"dspacing_c/F");  //d_spacing calculated by using the G4 position of the voxel (normal operation mode)
fnewt->Branch("lambda_c_wfm_corr",&lambda_c_wfm_corr,"lambda_c_wfm_corr/D"); //lambda calculated by using the G4 position of the voxel (WFM mode)
fnewt->Branch("dspacing_c_wfm",&dspacing_c_wfm,"dspacing_c_wfm/F");  //d_spacing calculated by using the G4 position of the voxel (WFM mode)

//used for the mapping of the electronics channels to detector voxels 
//integers > 0 correspond to the wire layer number given in the 2nd column 
// in the CDT file with channel mapping

int s3[4][4] = {{0,2,-1,-1}, {1,3,-1,-1}, {-1,-1,4,6}, {-1,-1,5,7}};  //sumo3
int s4[4][4] = {{-1,2,6,10}, {-1,3,7,11}, {0,4,8,-1}, {1,5,9,-1}};  //sumo4
int s5[4][4] = {{0,4,8,12}, {1,5,9,13}, {2,6,10,14}, {3,7,11,15}};  //sumo5
int s6[4][8] = {{0,-1,4,8,12,16,-1,-1},{1,-1,5,9,13,17,-1,-1},{-1,-1,2,6,10,14,-1,18},{-1,-1,3,7,11,15,-1,19}};  //sumo6

//time shifts and delays for the tof correction for events collected in WFM mode 

signed long int tdelay = 0;
signed long int tshift[6] = {0,2600,4400,6500,9270,11700};
signed long int tmin[6] = {9600,21940,33000,42540,51670,60300};
signed long int tmax[6] = {20680,31200,40420,48850,57800,70000};


for (Long64_t i=0; i<=nevents; i++){

  	tree->GetEntry(i);
		
	    vanode[i] = anode;
     	vcathode[i] = cathode;
        vtime[i] = neutronTime;     
        vsumo[i] = sumo;   
        vmodule[i] = module;
        vboardID[i] = boardID;  
        vchopperTime[i] = chopperTime;     
        vsubID[i] = subID;       

//      std::cout<<"sumo ="<<vsumo[i]<<", module ="<<vmodule[i]<<" , boardID ="<<boardID<<std::endl;
                
      if (i%100000 == 0) printf("reading procent out the original tree.....%2.3f\n",i*100./nevents);  
	
}

//the raw data included in the 'cdt_ev' tree in cdt.root file must be sorted and reorganised. 
// In the raw data the neutron and chopper events belong to different instances. 
// In the following each neutron event will get attached the timestamp of the most recent chopper event
// the new data will be stored an a new tree called 'cdt_new'. 

  for (Long64_t i=1; i<=nevents; i++){
  
 	if (vboardID[i] == 0) {vboardID[i] = vboardID[i-1]; vsumo[i] = vsumo[i-1];}
 	
 	if (vchopperTime[i] == 0) {vchopperTime[i] = vchopperTime[i-1];}
        
 		}
 		
 // time corrections and conversion from ASIC channels to segment, counter, wire and strip number

 file->cd();
 
  for (Long64_t i=0; i<=nevents; i++){
      
       if (vchopperTime[i] !=0 && vtime[i] !=0) {
       
			ntime = vtime[i];
			nanode = vanode[i];
			ncathode = vcathode[i];
			nboardID = vboardID[i];
			nchopperTime = vchopperTime[i];
			nsubID = vsubID[i];
			nsumo = vsumo[i];
			nmodule = vmodule[i];
			nevent = i;
			nmult = 0;
			
			tdiff_wfm = ntime-nchopperTime - 71428;  // 71428 = 1/14 
			ntof_wfm = tdiff_wfm; // tof for calculating lambda in wfm mode
			
//			tdiff = ntime-nchopperTime + 145230; // tcorr from FP's file (Dtt1 in normal mode)
//			tdiff = ntime-nchopperTime + 145230; // tcorr from FP's file (Dtt1 in normal mode)
			tdiff = ntime-nchopperTime + 238100; // tcorr from FP's file (Dtt1 in normal mode)
//			tdiff = ntime-nchopperTime + 238100 + 145230; // tcorr from FP's file (Dtt1 in normal mode)
			ntof = tdiff;
			ccycle = nchopperTime/744468;

			factor_a = floor(ncathode/16);
			factor_b = floor(nanode/16);
			
			switch (nsumo){
			
				case 3:
					nw_layer = s3[factor_b][factor_a];
					break;
					
				case 4:
					nw_layer = s4[factor_b][factor_a];
					break;

				case 5:
					nw_layer = s5[factor_b][factor_a];
					break;

				case 6:
					nw_layer = s6[factor_b][factor_a];
					break;
				
			}
						
			nsegment = floor(nw_layer/2) + 1;
			ncounter = nw_layer % 2 + 1;
			
			nwire = nanode - 16*factor_b + 1;
			nstrip = ncathode - 16*factor_a + 1;
			
			
			nindexS = nsumo*100000+nmodule*100+ncounter;
		    nindexC = nsegment*100000+nstrip*100+nwire;

			
//			this is the neutron tof used when running in normal mode
			
		if (ntof > 744468) { 
		
			ntof = ntof - 744468; 
			
			}

//			this is the neutron tof used when running in WFM mode

		if (ntof_wfm < 0) { 
		
			ntof_wfm = ntof_wfm + 744468; 
			
			}

        for (int ii=0; ii<=5; ii++){
            
             if (ntof_wfm>= 10*tmin[ii] && ntof_wfm<=10*tmax[ii]){
             
             	ntof_wfm_corr = ntof_wfm - 10*tshift[ii];}
             	
             }
	     
	    newt->Fill();
	    
	    }
	     
	   if (i%100000 == 0) printf("    now writing the new tree.....%2.3f\n",i*100./nevents);   
	     
}

// the nmult variable isn't working yet 
	
/*unsigned int b1,b2;
	 
  for (Long64_t i=0; i<=nevents; i++){
  
        newt->GetEntry(i);
        
        b1 = nboardID[i];
        
        newt->GetEntry(i+1);
        
        b2 = nboardID[i+1];
        
		if (b1 == b2) {
			mult++;
			
			} else {
				mult = 1;
				}
		
		nmult[i+1] = mult;	
		
		}*/
		
		

	newt->Write(0,TObject::kOverwrite); //because I don't want to create a new key for the new tree every time I execute the code
	std::cout<<"Writing the new cal tree...done"<<std::endl; 
	
	std::cout<<"number of entries in the new tree is "<<newt->GetEntries()<<std::endl; 
	std::cout<<" "<<std::endl; 
	
	delete[] vanode; delete[] vcathode; delete[] vsumo; delete[] vmodule; delete[] vboardID; delete[] vsubID;
	delete[] vchopperTime; delete[] vtime;
	
// in this part of the code the coordinates x,y,z of the voxel centre become part of the neutron event. 
// The 'endcap_lookup.root' file was created in the GEANT4 simulation of EndCap.
// In the simulation the detector was positioned to match the position of the EndCap detector in the V20 test at HZB. 
// this was the calculated voxels positions match the positions of the detector voxels in the real detector. 


    TFile* lfile = TFile::Open("endcap_lookup.root","read");
    
    if ( lfile->IsOpen() ) printf("endcap_lookup.root file created in GEANT4 opened successfully\n");
    
    printf("...now matching the data....please be patient!\n");

// call the event tree 

	TTree* ltree = (TTree*)lfile->Get("lookup_tree"); 
	
// each entry in the lookup_tree corresponds to a voxel that is uniquely identified 
// through the following information (included as branches (leaves) in the root tree)
// 			sumo            
// 			wire            
// 			seg             
// 			counter         
// 			strip           
// 			module         
// 			indexi          
// 			indexj          
// 			posx            
// 			posy            
// 			posz            

	TBranch *nx = ltree->GetBranch("posx");
	TLeaf *mx = nx->GetLeaf("posx");
	TBranch *ny = ltree->GetBranch("posy");
	TLeaf *my = ny->GetLeaf("posy");
	TBranch *nz = ltree->GetBranch("posz");
	TLeaf *mz = nz->GetLeaf("posz");
	TBranch *zz = ltree->GetBranch("counter");
	TLeaf *zza = zz->GetLeaf("counter");

// next line creates an index using the leaves "indexj" and "indexi" of the 'lookup_tree' tree.
// indexi and indexj were defined ReadLookupTable.C (used for the GEANT4 data) as: 
//	indexi = sumo*100000+module*100+counter;   
//	indexj = seg*100000+strip*100+wire;


	ltree->BuildIndex("indexj","indexi");
	
// The next for-loop goes over all events in the event tree with the new calibrated data ('cdt_new')
// created above and selects the leaves 'nikC' and 'nikS'.
// In case 'nikC' is equal with 'indexi' and 'nikS' is equal to 'indexj', the event in the lookup_tree containing the 
// information from the GEANT4 simulation and the event in 'cdt_new' tree containing the experimental data 
// will be synchronised. 
// Unfortunately this is the most time-consuming part of the code, but I don't know a better way of doing it.

	const int noev = newt->GetEntries();

	for (Long64_t k=0; k<noev; k++){
	
        newt->GetEntry(k);
	    
	    nikS = nindexS;
	    nikC = nindexC;
	    csumo = nsumo;
	    
	    Int_t ev = ltree->GetEntryNumberWithIndex(nikC,nikS);
	    	    
 	    ltree->GetEntry(ev);
		
//	        std::cout<<"event number in lookup rootfile = "<<ev<<std::endl;
//	        std::cout<<"event number in cdt rootfile = "<<k<<std::endl;
	     
	     ffile->cd();

//gets the coordinates of the detector voxel
	                
	     nvoxel_x = mx->GetValue(0);
	     nvoxel_y = my->GetValue(0);
	     nvoxel_z = mz->GetValue(0);	
	     nc = zza->GetValue(0);	
	     
//	     std::cout<<"x ="<<mx->GetValue(0)<<", y ="<< my->GetValue(0)<<", z ="<<mz->GetValue(0)<<std::endl; 

// calculate the distance from the voxel centre to the sample
	
	     rad_c = sqrt(nvoxel_x*nvoxel_x+nvoxel_y*nvoxel_y+nvoxel_z*nvoxel_z);  // in mm
	     
//	     twotheta_c = 90.+atan(nvoxel_z/nvoxel_x)*180/pi; // for GPS 

// calculate the theta angle of the voxel centre wrt the sample
		   
	     twotheta_c = 180 - acos(-nvoxel_z/rad_c)*180/pi; // initially
//	     twotheta_c = 180 - atan(nvoxel_y/rad_c)*180/pi; 
//	     twotheta_c = 180 - acos((rad_c*cos(nvoxel_x/rad_c))/sqrt(rad_c*rad_c+nvoxel_y*nvoxel_y))*180/pi; // try this for mantel
	     
// calculate the phi angle of the voxel centre wrt the sample

	     phi_c = 2*atan(nvoxel_x/rad_c)*180/pi; // for GPS and event file
//		 phi_c = acos(-nvoxel_y/sqrt(nvoxel_y*nvoxel_y+pow(rad_c*sin(nvoxel_x/rad_c),2)))*180./pi;

// calculate wavelength for the normal operation mode and the WFM mode	
     
//	     lambda_c = hp*ntof/mn/rad_c*m2a; //  for GPS
	     lambda_c_wfm_corr = hp*ntof_wfm_corr/mn/(22602+rad_c)*m2a/10e3; // 22.602 m from FP's file, WFM mode 
	     lambda_c = hp*ntof/mn/(29602+rad_c)*m2a/10e3; // 29.602 m from FP's file, normal operation mode
	     
// calculate d-spacing for the normal operation mode and the WFM mode	

	     dspacing_c_wfm = lambda_c_wfm_corr/2/sin(twotheta_c/2*pi/180.);
	     dspacing_c = lambda_c/2/sin(twotheta_c/2*pi/180.);

// save the leaves in the event tree 
	     
	     fnewt->Fill();
	     
}
	
	ffile->cd(); 
	fnewt->Write(0,TObject::kOverwrite); //because I don't want to create a new key for the new tree every time I execute the code
	std::cout<<"Writing the new cal tree...done"<<std::endl; 
	
	std::cout<<"number of entries in the new tree is "<<noev<<std::endl; 
	std::cout<<" "<<std::endl; 
	
// the are ROOT command lines to plot various branches 

//		cdt_new->AddFriend("cdt_new_cal","cdt_new_cal.root")
//      cdt_new->Scan("ntime:nchopperTime:tdiff:ccycle:ntof","")
//		cdt_new->Draw("tdiff/10000>>tt(800,0,80)","","")
//		cdt_new->Draw("ntof/10000>>tt(900,-10,80)","","")
//		cdt_new->Draw("ntof_wfm/10000>>tt(900,-10,80)","","")
//		cdt_new->Draw("lambda_c>>tt(1000,0,10)","","")
//		cdt_new->Draw("lambda_c_wfm_corr>>tt(1000,0,10)","","")
//		cdt_new->Draw("nvoxel_y:nvoxel_z","nsumo<=4","")
//		cdt_new->Draw("dspacing_c>>tt(600,0,6)","nsumo==6","")
//		cdt_new->Draw("twotheta_c>>tt(350,135,170)","","")
//		cdt_new->Draw("nvoxel_y:nvoxel_z>>tt(500,-1600,-1100,1000,200,1200)","","colz")
//		cdt_new->Draw("nvoxel_y:nvoxel_x>>tt(300,-100,200,1000,200,1200)","","colz")

// *************and some automatic plotting 
	
	
	TCanvas *can=new TCanvas("can","can",100,100,700,700);

	gStyle->SetOptTitle(1);
	gStyle->SetOptStat(1);

	can->SetFillColor(0);
	can->SetGrid();
	Float_t small=1e-5;

	can->Divide(1,3,small,small);
	can->ToggleEventStatus();
	
	can->cd(1);
	newt->Draw("ntof_wfm/10000>>tt(800,-1,80)","","");
//	newt->Draw("ncathode:nanode","nsumo==3","colz");
	can->cd(2);
	newt->Draw("ntof_wfm_corr/10000>>tt1(800,-1,80)","","");
	can->cd(3);
	newt->Draw("ntof/10000>>tt2(1000,1,100)","","");
//	newt->Draw("ncathode:nanode","nsumo==4","colz");
//	can->cd(3);
//	newt->Draw("ncathode:nanode","nsumo==5","colz");
//	can->cd(4);
//	newt->Draw("ncathode:nanode","nsumo==6","colz");

/*	TCanvas *canT=new TCanvas("canT","canT",100,100,800,800);

	gStyle->SetOptTitle(1);
	gStyle->SetOptStat(1);

	canT->SetFillColor(0);
	canT->SetGrid();

	canT->Divide(2,2,small,small);
	
	canT->cd(1);
	newt->Draw("(ntime-nchopperTime)/10e6>>tt(2000,-0.1,0.1)","nboardID==1418045");
	canT->cd(2);
	newt->Draw("(ntime-nchopperTime)/10e6","nboardID==1416964");
	canT->cd(3);
	newt->Draw("(ntime-nchopperTime)/10e6","nboardID==1416799");
	canT->cd(4);
	newt->Draw("(ntime-nchopperTime)/10e6","nboardID==1416697");*/

}

  int main() {
  analysis();  
  return 0;
   }
