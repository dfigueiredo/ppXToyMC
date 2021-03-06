import sys
import FWCore.ParameterSet.Config as cms

miniAOD = False

fileinput = 'file:///afs/cern.ch/user/d/dmf/private/work/private/CMSPhysicsAnalysis/PrivateMCProduction/PPSMCProduction/working/RunIISummer20UL17RECO.root'
#fileinput = 'file:///afs/cern.ch/user/d/dmf/private/work/private/CMSPhysicsAnalysis/PrivateMCProduction/PPSMCProduction/working/RunIISummer20UL17GEN.root'
#fileinput = '/store/user/dmf/ZXToyMC-13TeV-miniAOD_2020-03-23_UTC13-55-30/ZXToyMC-GEN_Muon_PostTS2_950_150/ZXToyMC-13TeV-miniAOD_2020-03-23_UTC13-55-30/200323_125619/0000/output_22.root'


gentag = ''
genjettag = ''

if miniAOD:
  gentag = 'prunedGenParticles'
  genjettag = 'slimmedGenJets'
else:
  gentag = 'genParticles'
  genjettag = 'ak4GenJets'


process = cms.Process('Analysis')

#-----------------
# General options
#-----------------

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True),
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

#-------------
# Input files    
#-------------
process.source = cms.Source("PoolSource",
				fileNames = cms.untracked.vstring(fileinput),
)

#-------------
# Preskimming      
#-------------
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v17') # 94X MC campaing
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

#---------------
# Proton Filter     
#---------------
process.genParticle = cms.EDAnalyzer("GenParticleShow",
					    GenPartTag = cms.InputTag(gentag),
					    GenJetTag = cms.InputTag(genjettag),
					    EBeam = cms.double(6500.),
					    DebugProtons = cms.bool(True),
					    DebugPF = cms.bool(True),
					    DebugJets = cms.bool(False)
					)

# prepare the output file
process.TFileService = cms.Service('TFileService',
    fileName = cms.string("gen_output.root"),
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(process.genParticle)
