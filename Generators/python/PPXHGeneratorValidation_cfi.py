import FWCore.ParameterSet.Config as cms

ppxhGeneratorValidation = cms.EDAnalyzer("PPXHGeneratorValidation",
  verbosity = cms.untracked.uint32(0),

  tagHepMC = cms.InputTag("generator", "unsmeared"),
  tagRecoTracks = cms.InputTag("ctppsLocalTrackLiteProducer"),
  tagRecoProtonsSingleRP = cms.InputTag("ctppsProtons", "singleRP"),
  tagRecoProtonsMultiRP = cms.InputTag("ctppsProtons", "multiRP"),

  outputFile = cms.string("ppxhGeneratorValidation.root")
)
