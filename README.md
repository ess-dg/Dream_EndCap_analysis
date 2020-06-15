# Dream_EndCap_analysis

This repository contains the ROOT scripts for analyzing the DREAM EndCap data collected at V20 in July 2019.

Author: Irina Stefanescu, ESS DG. Email: irina.stefanescu@ess.eu.

First version of the codes pulled in June 2020. 

Description of the repository
-----------------------------

The /data/ directory contains some data that can be used to test the codes. 

The /documentation/ directory contains some documentation relevant to the experiment and data analysis. 

The /src/ directory contains three C++ files that read out the binary data files, outputs the event information in ascii format. this ascii file will be use to create a ROOT raw event tree. The data analysis will be perfomed on the raw event tree and will consists  of corrections of the time-of-flights and mapping of the anode and cathode IDs with the x, y, z positions of the detector voxels calculated in the EndCap simulation project (with GEANT4).   

Additional information on the data format and analysis is given inside the directories. 
