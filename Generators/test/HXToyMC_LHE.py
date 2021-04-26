#flake8: noqa

import sys
from FWCore.ParameterSet.VarParsing import VarParsing

# Setting Input Parameters from Line Command
options = VarParsing ('analysis')
options.register('Mass',850,VarParsing.multiplicity.singleton, VarParsing.varType.int,"X Mass")
options.parseArguments()

print("")
print("Mass: %s"%options.Mass)
print("")

import FWCore.ParameterSet.Config as cms

# random seeds
RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
  sourceSeed = cms.PSet(initialSeed = cms.untracked.uint32(98765)),
  generator = cms.PSet(initialSeed = cms.untracked.uint32(98766)),
  beamDivergenceVtxGenerator = cms.PSet(initialSeed = cms.untracked.uint32(3849)),
  ctppsDirectProtonSimulation = cms.PSet(initialSeed = cms.untracked.uint32(4981))
)

# number of events
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(500)
)

# redefine particle generator
process.load("ppXToyMC.Generators.PPSToyMcLHE_cfi")
process.generator.verbosity = False
process.generator.m_X_mean = options.Mass
process.generator.m_S_mean = options.Mass + 200
process.generator.m_f = (options.Mass - 200)/2
process.generator.decayX = True
process.generator.decayHToBottoms = True

# distribution plotter
process.ctppsTrackDistributionPlotter = cms.EDAnalyzer("CTPPSTrackDistributionPlotter",
  tagTracks = cms.InputTag("ctppsLocalTrackLiteProducer"),
  outputFile = cms.string("output.root")
)

# acceptance plotter
process.ctppsAcceptancePlotter = cms.EDAnalyzer("CTPPSAcceptancePlotter",
  tagHepMC = cms.InputTag("generator", "unsmeared"),
  tagTracks = cms.InputTag("ctppsLocalTrackLiteProducer"),

  rpId_45_F = process.rpIds.rp_45_F,
  rpId_45_N = process.rpIds.rp_45_N,
  rpId_56_N = process.rpIds.rp_56_N,
  rpId_56_F = process.rpIds.rp_56_F,

  outputFile = cms.string("acceptance.root")
)

# generator plots
process.load("ppXToyMC.Generators.PPXHGeneratorValidation_cfi")
process.ppxhGeneratorValidation.tagHepMC = cms.InputTag("generator", "unsmeared")
process.ppxhGeneratorValidation.tagRecoTracks = cms.InputTag("ctppsLocalTrackLiteProducer")
process.ppxhGeneratorValidation.tagRecoProtonsSingleRP = cms.InputTag("ctppsProtons", "singleRP")
process.ppxhGeneratorValidation.tagRecoProtonsMultiRP = cms.InputTag("ctppsProtons", "multiRP")
process.ppxhGeneratorValidation.referenceRPDecId_45 = process.rpIds.rp_45_F
process.ppxhGeneratorValidation.referenceRPDecId_56 = process.rpIds.rp_56_F
process.ppxhGeneratorValidation.outputFile = "ppxhGeneratorValidation.root"

# processing path
process.p = cms.Path(
  process.generator
  #*process.ppxhGeneratorValidation
  #*process.ctppsAcceptancePlotter
)

# output configuration
#process.output = cms.OutputModule("PoolOutputModule",
#  fileName = cms.untracked.string("output.root"),
#  splitLevel = cms.untracked.int32(0),
#  eventAutoFlushCompressedSize=cms.untracked.int32(-900),
#  compressionAlgorithm=cms.untracked.string("LZMA"),
#  compressionLevel=cms.untracked.int32(9),
#  outputCommands = cms.untracked.vstring(
#    'drop *',
#    'keep edmHepMCProduct_*_*_*',
#    'keep LHE*_*_*_*',
#    'keep CTPPSLocalTrackLites_*_*_*',
#    'keep recoForwardProtons_*_*_*'
#  )
#)

process.outpath = cms.EndPath(process.output)

