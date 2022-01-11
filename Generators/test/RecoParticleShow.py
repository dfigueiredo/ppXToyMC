import sys
import FWCore.ParameterSet.Config as cms

#fileinput = 'file:///afs/cern.ch/user/d/dmf/private/work/private/CMSPhysicsAnalysis/PrivateMCProduction/PPSMCProduction/working/RunIISummer20UL17RECO.root'
fileinput = '/store/user/lpagliai/ZXToyMC-RECO_Electron_PostTS2_950_150_2020-03-11_UTC14-49-56/ZXToyMC-GEN_Electron_PostTS2_950_150/ZXToyMC-RECO_Electron_PostTS2_950_150_2020-03-11_UTC14-49-56/200311_135024/0000/output_24.root'

recotag = 'particleFlow'
recojettag = 'ak4PFJets'

gentag = 'genParticles'
genjettag = 'ak4GenJets'

process = cms.Process('Analysis')

#-----------------
# Recoeral options
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
process.recoParticle = cms.EDAnalyzer("RecoParticleShow",
					    RecoPartTag = cms.InputTag(recotag),
					    RecoJetTag = cms.InputTag(recojettag),
                                            GenPartTag = cms.InputTag(gentag),
                                            GenJetTag = cms.InputTag(genjettag),
					    EBeam = cms.double(6500.),
					    MatchingMC = cms.bool(True),
					    DebugProtons = cms.bool(False),
					    DebugPF = cms.bool(False),
					    DebugJets = cms.bool(False)
					)

# prepare the output file
process.TFileService = cms.Service('TFileService',
    fileName = cms.string("reco_output.root"),
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(process.recoParticle)
