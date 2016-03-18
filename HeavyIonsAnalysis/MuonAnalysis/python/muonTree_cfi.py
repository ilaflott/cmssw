import FWCore.ParameterSet.Config as cms

muonTree = cms.EDAnalyzer(
    "muonTree",
    doGenParticles     = cms.bool(False),
    runOnParticleGun   = cms.bool(False),
    pileupCollection   = cms.InputTag("addPileupInfo"),
    genParticleSrc     = cms.InputTag("genParticles"),
    recoMuonSrc        = cms.InputTag("muons"),
    VtxLabel           = cms.InputTag("offlinePrimaryVertices")
    #beamSpot           = cms.InputTag('offlineBeamSpot'),
    #particleFlowCollection = cms.InputTag("particleFlowTmp"),
)
