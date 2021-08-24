# LHE Toy MC

A toy MC kinematics model to emulate proton-proton (pp) collisions, (S) intermediate scalars and final hidden particles. This package includes emulations of pp to ppS(ZX) or pp to ppS(HX) where S decays to [H(bbar) || Z(MuMu||EE)] + X(NuNu). The configuration files placed in "test/" folder are able to generate LHE files, which are used within CMSSW framework to produce a fully detector simulated sample, including a fast simulation of the protons detected by PPS roman pots. 

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
cmsRun XHToyMC_LHE.py (or XZToyMC_LHE.py) 
```

## Validation

There is a script "test/GenParticleShow.py" that can be used to cross check your CMSSW GEN collection. It is possible to printout the protons, particle flow (PF) and GenJets kinematics. Furthermore, a tree file is produced for futher analysis.

```bash
cd test/
cmsRun GenParticleShow.py
```

There are options flags that can be used:

```python
process.genParticle = cms.EDAnalyzer("GenParticleShow",
					    GenPartTag = cms.InputTag('genParticles'||'prunedGenParticles'),
					    GenJetTag = cms.InputTag('ak4GenJets'||'slimmedGenJets'),
					    EBeam = cms.double(6500.),
					    DebugProtons = cms.bool(True),
					    DebugPF = cms.bool(True),
					    DebugJets = cms.bool(True)
					)
```
