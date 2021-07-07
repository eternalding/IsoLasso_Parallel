# IsoLasso2.0
The parallel version of IsoLasso, a bioinformatic tool for transciptome assembly and expression level estimation.

## Install
* Clone Isolasso to your workspace
* cd IsoLasso_Parallel
* mkdir build
* cd build
* cmake ../
* make
## Unittest
* ./unittest (if unittest passed, IsoLasso shall work properly)
## Normal Usage
* ./IsoLasso_Parallel -i "YOUR SAM FILE" -o "OUTPUTNAME" </br>
An file named "OUTPUTNAME_RGStats.txt" will be generated, recording assembled read group information.
