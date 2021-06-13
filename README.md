# LHE Toy MC

A toy MC kinematics model to emulate proton-proton (pp) collisions, (S) intermediate scalars and final hidden particles. This package includes emulations of pp to ppXZ or pp to ppSX where S decays to H(bbar). The configuration files placed in "test/" folder are able to generate LHE files, which are used within CMSSW framework to produce a fully detector simulated sample, including a fast simulation of the protons detected by PPS roman pots. 

The package to generate a full CMS sample is placed [here](https://github.com/dfigueiredo/PPSMCProduction).

## Installation

Install CMSSW, clone the repo and compile the package.

```bash
cmsrel CMSSW_10_6_17
cd CMSSW_10_6_17/src
cmsenv
git clone https://github.com/dfigueiredo/ppXToyMC
scram b -j 8
```

## Usage

```bash
cd test/
cmsRun HXToyMC_LHE.py 
```
