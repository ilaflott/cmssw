
import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.patHeavyIonSequences_cff import *
from CmsHi.JetAnalysis.inclusiveJetAnalyzer_cff import *


akPu2PFmatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akPu2PFJets"),
    matched = cms.InputTag("ak2HiGenJets")
    )

akPu2PFparton = patJetPartonMatch.clone(src = cms.InputTag("akPu2PFJets"),
                                                        matched = cms.InputTag("hiGenParticles")
                                                        )

akPu2PFcorr = patJetCorrFactors.clone(
    useNPV = False,
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),                                                                
    src = cms.InputTag("akPu2PFJets"),
    payload = "AK2PF_hiIterativeTracks"
    )

akPu2PFpatJets = patJets.clone(jetSource = cms.InputTag("akPu2PFJets"),
                                               jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu2PFcorr")),
                                               genJetMatch = cms.InputTag("akPu2PFmatch"),
                                               genPartonMatch = cms.InputTag("akPu2PFparton"),
                                               jetIDMap = cms.InputTag("akPu2PFJetID"),
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

akPu2PFAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akPu2PFpatJets"),
                                                             genjetTag = 'ak2HiGenJets',
                                                             rParam = 0.5,
                                                             matchJets = cms.untracked.bool(True),
                                                             matchTag = 'akPu2PFpatJets',
                                                             pfCandidateLabel = cms.untracked.InputTag('particleFlowTmp'),
                                                             trackTag = cms.InputTag("hiGeneralTracks"),
                                                             fillGenJets = True,
                                                             isMC = True,
                                                             genParticles = cms.untracked.InputTag("hiGenParticles")
                                                             )


akPu2PFJetSequence_mc = cms.Sequence(akPu2PFmatch
                                                  *
                                                  akPu2PFparton
                                                  *
                                                  akPu2PFcorr
                                                  *
                                                  akPu2PFpatJets
                                                  *
                                                  akPu2PFAnalyzer
                                                  )

akPu2PFJetSequence_data = cms.Sequence(akPu2PFcorr
                                                    *
                                                    akPu2PFpatJets
                                                    *
                                                    akPu2PFAnalyzer
                                                    )

akPu2PFJetSequence = cms.Sequence(akPu2PFJetSequence_data)
