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
* Cretate TnP pairs from a ntuple (possible to modify the included branches, the pair selection cuts and the ID of interest for efficiency calculation):
```
root -q -b CreateTnPpairsData.cpp "()"
root -q -b CreateTnPpairsMC.cpp "()"

root -q -b CreateTnPpairsData_Run3.cpp "()"
root -q -b CreateTnPpairsMC_Run3.cpp "()"

```
* Perform fits and get the efficiencies (figures will be stored in plots directory). pT and rel_Iso bins can be defined inside TnPAnalyzer_*.cpp files. A corresponding TnPanalyzer file is present for noIso and Iso ID efficiency calculations.
```
root -q -b TnPAnalyzer_noIsoID.cpp "()"
root -q -b TnPAnalyzer_isoID.cpp "()"
```
* If the TnPAnalyzer for isoID is run, then additional extrapolation step is needed to be done to calcualte the SF values. A simplistic code is availabe and can be run here with python 3.
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc8-opt/setup.sh
python isoID_eff_extrapolator.py
```