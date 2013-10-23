
import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.patHeavyIonSequences_cff import *
from CmsHi.JetAnalysis.inclusiveJetAnalyzer_cff import *


akVs7Calomatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akVs7CaloJets"),
    matched = cms.InputTag("ak7HiGenJets")
    )

akVs7Caloparton = patJetPartonMatch.clone(src = cms.InputTag("akVs7CaloJets"),
                                                        matched = cms.InputTag("hiGenParticles")
                                                        )

akVs7Calocorr = patJetCorrFactors.clone(
    useNPV = False,
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),                                                                
    src = cms.InputTag("akVs7CaloJets"),
    payload = "AK7Calo_HI"
    )

akVs7CalopatJets = patJets.clone(jetSource = cms.InputTag("akVs7CaloJets"),
                                               jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akVs7Calocorr")),
                                               genJetMatch = cms.InputTag("akVs7Calomatch"),
                                               genPartonMatch = cms.InputTag("akVs7Caloparton"),
                                               jetIDMap = cms.InputTag("akVs7CaloJetID"),
                                               addBTagInfo         = False,
                                               addTagInfos         = False,
                                               addDiscriminators   = False,
                                               addAssociatedTracks = False,
                                               addJetCharge        = False,
                                               addJetID            = False,
                                               getJetMCFlavour     = False,
                                               addGenPartonMatch   = False,
                                               addGenJetMatch      = False,
                                               embedGenJetMatch    = False,
                                               embedGenPartonMatch = False,
                                               embedCaloTowers     = False,
                                               embedPFCandidates = False
				            )

akVs7CaloJetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akVs7CalopatJets"),
                                                             genjetTag = 'ak7HiGenJets',
                                                             rParam = 0.7,
                                                             matchJets = cms.untracked.bool(False),
                                                             matchTag = 'akVs7CalopatJets',
                                                             pfCandidateLabel = cms.untracked.InputTag('particleFlowTmp'),
                                                             trackTag = cms.InputTag("hiGeneralTracks"),
                                                             fillGenJets = False,
                                                             isMC = False,
                                                             genParticles = cms.untracked.InputTag("hiGenParticles")
                                                             )

akVs7CaloJetSequence_mc = cms.Sequence(akVs7Calomatch
                                                  *
                                                  akVs7Caloparton
                                                  *
                                                  akVs7Calocorr
                                                  *
                                                  akVs7CalopatJets
                                                  *
                                                  akVs7CaloJetAnalyzer
                                                  )

akVs7CaloJetSequence_data = cms.Sequence(akVs7Calocorr
                                                    *
                                                    akVs7CalopatJets
                                                    *
                                                    akVs7CaloJetAnalyzer
                                                    )

akVs7CaloJetSequence = cms.Sequence(akVs7CaloJetSequence_data)
