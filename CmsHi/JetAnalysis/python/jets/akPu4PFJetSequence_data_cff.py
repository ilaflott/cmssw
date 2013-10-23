
import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.patHeavyIonSequences_cff import *
from CmsHi.JetAnalysis.inclusiveJetAnalyzer_cff import *


akPu4PFmatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akPu4PFJets"),
    matched = cms.InputTag("ak4HiGenJets")
    )

akPu4PFparton = patJetPartonMatch.clone(src = cms.InputTag("akPu4PFJets"),
                                                        matched = cms.InputTag("hiGenParticles")
                                                        )

akPu4PFcorr = patJetCorrFactors.clone(
    useNPV = False,
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),                                                                
    src = cms.InputTag("akPu4PFJets"),
    payload = "AK4PF_hiIterativeTracks"
    )

akPu4PFpatJets = patJets.clone(jetSource = cms.InputTag("akPu4PFJets"),
                                               jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu4PFcorr")),
                                               genJetMatch = cms.InputTag("akPu4PFmatch"),
                                               genPartonMatch = cms.InputTag("akPu4PFparton"),
                                               jetIDMap = cms.InputTag("akPu4PFJetID"),
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

akPu4PFAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akPu4PFpatJets"),
                                                             genjetTag = 'ak4HiGenJets',
                                                             rParam = 0.5,
                                                             matchJets = cms.untracked.bool(True),
                                                             matchTag = 'akPu4PFpatJets',
                                                             pfCandidateLabel = cms.untracked.InputTag('particleFlowTmp'),
                                                             trackTag = cms.InputTag("hiGeneralTracks"),
                                                             fillGenJets = True,
                                                             isMC = True,
                                                             genParticles = cms.untracked.InputTag("hiGenParticles")
                                                             )


akPu4PFJetSequence_mc = cms.Sequence(akPu4PFmatch
                                                  *
                                                  akPu4PFparton
                                                  *
                                                  akPu4PFcorr
                                                  *
                                                  akPu4PFpatJets
                                                  *
                                                  akPu4PFAnalyzer
                                                  )

akPu4PFJetSequence_data = cms.Sequence(akPu4PFcorr
                                                    *
                                                    akPu4PFpatJets
                                                    *
                                                    akPu4PFAnalyzer
                                                    )

akPu4PFJetSequence = cms.Sequence(akPu4PFJetSequence_data)
