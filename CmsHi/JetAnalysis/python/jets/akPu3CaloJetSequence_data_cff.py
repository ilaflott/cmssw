
import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.patHeavyIonSequences_cff import *
from CmsHi.JetAnalysis.inclusiveJetAnalyzer_cff import *


akPu3Calomatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akPu3CaloJets"),
    matched = cms.InputTag("ak3HiGenJets")
    )

akPu3Caloparton = patJetPartonMatch.clone(src = cms.InputTag("akPu3CaloJets"),
                                                        matched = cms.InputTag("hiGenParticles")
                                                        )

akPu3Calocorr = patJetCorrFactors.clone(
    useNPV = False,
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),                                                                
    src = cms.InputTag("akPu3CaloJets"),
    payload = "AK3Calo_hiIterativeTracks"
    )

akPu3CalopatJets = patJets.clone(jetSource = cms.InputTag("akPu3CaloJets"),
                                               jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu3Calocorr")),
                                               genJetMatch = cms.InputTag("akPu3Calomatch"),
                                               genPartonMatch = cms.InputTag("akPu3Caloparton"),
                                               jetIDMap = cms.InputTag("akPu3CaloJetID"),
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

akPu3CaloAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akPu3CalopatJets"),
                                                             genjetTag = 'ak3HiGenJets',
                                                             rParam = 0.5,
                                                             matchJets = cms.untracked.bool(True),
                                                             matchTag = 'akPu3CalopatJets',
                                                             pfCandidateLabel = cms.untracked.InputTag('particleFlowTmp'),
                                                             trackTag = cms.InputTag("hiGeneralTracks"),
                                                             fillGenJets = True,
                                                             isMC = True,
                                                             genParticles = cms.untracked.InputTag("hiGenParticles")
                                                             )


akPu3CaloJetSequence_mc = cms.Sequence(akPu3Calomatch
                                                  *
                                                  akPu3Caloparton
                                                  *
                                                  akPu3Calocorr
                                                  *
                                                  akPu3CalopatJets
                                                  *
                                                  akPu3CaloAnalyzer
                                                  )

akPu3CaloJetSequence_data = cms.Sequence(akPu3Calocorr
                                                    *
                                                    akPu3CalopatJets
                                                    *
                                                    akPu3CaloAnalyzer
                                                    )

akPu3CaloJetSequence = cms.Sequence(akPu3CaloJetSequence_data)
