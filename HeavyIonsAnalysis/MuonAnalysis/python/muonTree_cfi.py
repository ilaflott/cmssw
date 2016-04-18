import FWCore.ParameterSet.Config as cms

muonTree = cms.EDAnalyzer(
    "muonTree",
    doGenParticles     = cms.bool(False),
    runOnParticleGun   = cms.bool(False),
    minPt     = cms.double(1.0),
    minAbsEta = cms.double(0.0),
    maxAbsEta = cms.double(100.0),
    pileupCollection   = cms.InputTag("addPileupInfo"),
    genParticleSrc     = cms.InputTag("genParticles"),
    recoMuonSrc        = cms.InputTag("muons"),
    VtxLabel           = cms.InputTag("offlinePrimaryVertices")
    #beamSpot           = cms.InputTag('offlineBeamSpot'),
    #particleFlowCollection = cms.InputTag("particleFlowTmp"),
)
