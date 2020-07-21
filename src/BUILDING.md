# Building

## Building with cmake
_cmake builds are under development and may not yet work as intended._

    > mkdir build
    > cd build
    > cmake ..
    > make

After this the binaries will be in _build/bin_

## Running
and execute it by typing

    > ./bin/read_bin_irina /path_data_directory/filename.bin >> test.out
    > ./bin/CreateCDTRootFile
    > ./bin/analysis

## Building with ROOT
For the ROOT macros there are two ways of running the macros

### Using root interpreter
(1) the usual way, just load the macro, execute the function inside and quit ROOT

>root -l

root[1] .L CreateCDTTree.C
root[2] CreateCDTRootFile()
root[3].q

>root -l

root[1] .L analysis_cdt.C
root[2] analysis()
root[3].q

// at the end of these steps one should have two rootfiles available:
// cdt.root -> the main file with the data
// cdt_new_cal.root  -> the file with the positions of the voxels

Open the ROOT file and look at the data:

>root -l cdt.root

### Command line
(This should be deleted a cmake does this???)
**********(2) the advanced way, by compiling the macros


To create a stand-alone program from the macro called analysis_cdt.C, simply type in a terminal window:

> g++ -o analysis analysis_cdt.C `root-config --cflags --libs`

 and execute it by typing:

>./analysis

Do the same for the macro 'CreateCDTTree.C'

>g++ -o CreateCDTRootFile CreateCDTTree.C `root-config --cflags --libs`

> ./CreateCDTRootFile

// same as above, at the end of these steps one should have two rootfiles available:
// cdt.root -> the main file with the data
// cdt_new_cal.root  -> the file with the positions of the voxels

Open the ROOT file and look at the data:

root -l cdt.root
