import FWCore.ParameterSet.Config as cms

generator = cms.EDProducer("PPXHGenerator",
  verbosity = cms.untracked.uint32(1),
  decayX = cms.bool(True),
  decayHToBottoms = cms.bool(True),
  m_S_mean = cms.double(1200), #Scalar particle produced in pp collision
  m_S_gamma = cms.double(10), 
  m_Z_mean = cms.double(91.1876),
  m_Z_gamma = cms.double(2.4952),
  m_H_mean = cms.double(125.18),
  m_H_gamma = cms.double(0.013), #Upper Limit @ 95% CL
  m_X_mean = cms.double(1050), #Particle from S->H+X decay
  m_X_gamma = cms.double(10),

  m_f = cms.double(525),

  m_e = cms.double(0.5109e-3),
  m_mu = cms.double(105.658e-3),
  m_b = cms.double(4.18),

  p_beam = cms.double(6500),

  p_z_LAB_2p_min = cms.double(-1500),
  p_z_LAB_2p_max = cms.double(+1500),

  p_T_H_min = cms.double(0.)
)
