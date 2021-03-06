import FWCore.ParameterSet.Config as cms

generator = cms.EDProducer("PPSToyMcLHE_XH",
  verbosity = cms.untracked.uint32(1),
  decayX = cms.bool(True),
  decayHToBottoms = cms.bool(True),
  m_S_mean = cms.double(1200), #Scalar particle produced in pp collision
  m_S_gamma = cms.double(10), 
  m_H_mean = cms.double(125.18),
  m_H_gamma = cms.double(0.013), #Upper Limit @ 95% CL
  m_X_mean = cms.double(1050), #Particle from S->H+X decay
  m_X_gamma = cms.double(10),
  m_f = cms.double(525),
  m_b = cms.double(4.18),
  p_beam = cms.double(6500),
  p_z_LAB_2p_min = cms.double(-1500),
  p_z_LAB_2p_max = cms.double(+1500),
  p_T_H_min = cms.double(0.),
  lheOutputFile = cms.untracked.string("toy_mc_pps_XH.lhe")
)
