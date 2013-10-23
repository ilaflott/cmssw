
import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.patHeavyIonSequences_cff import *
from CmsHi.JetAnalysis.inclusiveJetAnalyzer_cff import *


akVs4Calomatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akVs4CaloJets"),
    matched = cms.InputTag("ak4HiGenJets")
    )

akVs4Caloparton = patJetPartonMatch.clone(src = cms.InputTag("akVs4CaloJets"),
                                                        matched = cms.InputTag("hiGenParticles")
                                                        )

akVs4Calocorr = patJetCorrFactors.clone(
    useNPV = False,
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),                                                                
    src = cms.InputTag("akVs4CaloJets"),
    payload = "AK4Calo_hiIterativeTracks"
    )

akVs4CalopatJets = patJets.clone(jetSource = cms.InputTag("akVs4CaloJets"),
                                               jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akVs4Calocorr")),
                                               genJetMatch = cms.InputTag("akVs4Calomatch"),
                                               genPartonMatch = cms.InputTag("akVs4Caloparton"),
                                               jetIDMap = cms.InputTag("akVs4CaloJetID"),
                                               addBTagInfo         = False,
                                               addTagInfos         = False,
                                               addDiscriminators   = False,
                                               addAssociatedTracks = False,
                                               addJetCharge        = False,
                                               addJetID            = True,
                                               getJetMCFlavour     = False,
                                               addGenPartonMatch   = True,
                                               addGenJetMatch      = True,
                                               embedGenJetMatch    = True,
                                               embedGenPartonMatch = True,
                                               embedCaloTowers     = False,
				            )

akVs4CaloAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akVs4CalopatJets"),
                                                             genjetTag = 'ak4HiGenJets',
                                                             rParam = 0.5,
                                                             matchJets = cms.untracked.bool(True),
                                                             matchTag = 'akVs4CalopatJets',
                                                             pfCandidateLabel = cms.untracked.InputTag('particleFlowTmp'),
                                                             trackTag = cms.InputTag("hiGeneralTracks"),
                                                             fillGenJets = True,
                                                             isMC = True,
                                                             genParticles = cms.untracked.InputTag("hiGenParticles")
                                                             )


akVs4CaloJetSequence_mc = cms.Sequence(akVs4Calomatch
                                                  *
                                                  akVs4Caloparton
                                                  *
                                                  akVs4Calocorr
                                                  *
                                                  akVs4CalopatJets
                                                  *
                                                  akVs4CaloAnalyzer
                                                  )

akVs4CaloJetSequence_data = cms.Sequence(akVs4Calocorr
                                                    *
                                                    akVs4CalopatJets
                                                    *
                                                    akVs4CaloAnalyzer
                                                    )

akVs4CaloJetSequence = cms.Sequence(akVs4CaloJetSequence_mc)
