# TnP_LowPtElectrons
## Description
This page contains instructions and examples of a way to perform Tag and Probe on LowPt electrons
## Usage instructions

* Set up your area:
```
cmsrel CMSSW_10_2_10
cd CMSSW_10_2_10/src
cmsenv
mkdir TnP_fits
cd TnP_fits
mkdir plots
```
* Get the TnP tool:
```
git clone https://github.com/apetkovi1/TnP_LowPtElectrons.git
```
* Cretate TnP pairs:
```
root -q -b CreateTnPpairs.cpp "()"
```
* Perform fits and get the efficiencies (figures will be stored in plots directory). Bins and CMS part (Barrel or Endcap) can be selected inside TnPAnalyzer.cpp
```
root -q -b TnPAnalyzer.cpp "()"
```
