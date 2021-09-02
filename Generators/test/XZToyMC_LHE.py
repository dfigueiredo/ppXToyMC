#flake8: noqa

import sys
from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

# Setting Input Parameters from Line Command
options = VarParsing ('analysis')
options.register('Mass',950,VarParsing.multiplicity.singleton, VarParsing.varType.int,"X Mass")
options.parseArguments()

process = cms.Process("CTPPSDirectSimulation")

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",

    # Include a PSet for each module label that needs a
    # random engine.  The name is the module label.
    # You must supply a seed or seeds.
    # Optionally an engine type can be specified

    generator = cms.PSet(
        initialSeed = cms.untracked.uint32(123456789),
        engineName = cms.untracked.string('TRandom3') #HepJamesRandom, RanecuEngine 
    )

    # This is optional.  If you want the service to save the state
    # of all engines to a separate text file which is overwritten before
    # modules begin processing on each event. The text file is only
    # needed for one type of replay. The filename can be anything
    # you want but the replay process will need to reference it.
    #,saveFileName = cms.untracked.string('RandomEngineStates.txt')
)

# redefine particle generator
process.load("ppXToyMC.Generators.PPSToyMcLHE_XZ_cfi")
process.generator.verbosity = 0
process.generator.m_X = options.Mass
process.generator.m_XZ_min = options.Mass + 100
process.generator.m_X_pr1 = options.Mass - 100
process.generator.decayX = True
process.generator.m_S_mean = options.Mass + 100
process.generator.m_S_gamma = 0.01 * (options.Mass + 100)
process.generator.useResonantIntermediateState = False

# This is optional.  If you want the service to save the state
# of all engines to each Event and LuminosityBlock, then
# include this producer here and in the path.  This is only
# needed for one type of replay.  The module label can be
# anything, but the replay process will need to reference it.
process.randomEngineStateProducer = cms.EDProducer("RandomEngineStateProducer")

# The RandomNumberGeneratorService should work with
# any kind of source
process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500000)
)

process.p = cms.Path(process.generator+process.randomEngineStateProducer)

