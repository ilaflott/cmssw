#ifndef ggHiNtuplizerMuTree_h
#define ggHiNtuplizerMuTree_h

#include "TTree.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

using namespace std;

class ggHiNtuplizerMuTree : public edm::EDAnalyzer {

 public:

   ggHiNtuplizerMuTree(const edm::ParameterSet&);
   virtual ~ggHiNtuplizerMuTree() {};

 private:

   // event routines
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void cleanUp(const edm::Event&, const edm::EventSetup&);

   // fill routines
   void fillGenPileupInfo(const edm::Event&);
   void fillGenParticles (const edm::Event&);
   void fillMuons        (const edm::Event&, const edm::EventSetup&, const reco::Vertex&);   

   // helper routines
   float getGenCalIso(edm::Handle<vector<reco::GenParticle> >&, reco::GenParticleCollection::const_iterator, 
		      float dRMax, bool removeMu, bool removeNu);
   float getGenTrkIso(edm::Handle<vector<reco::GenParticle> >&, reco::GenParticleCollection::const_iterator, 
		      float dRMax);

   // handles to collections of objects
   edm::EDGetTokenT<vector<PileupSummaryInfo> >  genPileupCollection_    ;
   edm::EDGetTokenT<vector<reco::GenParticle> >  genParticlesCollection_ ;
   edm::EDGetTokenT<edm::View<reco::Muon> >      recoMuonsCollection_    ;
   edm::EDGetTokenT<vector<reco::Vertex> >       vtxCollection_          ;
   //edm::EDGetTokenT<edm::View<reco::PFCandidate> >    pfCollection_;
   //edm::EDGetTokenT<reco::BeamSpot>                   beamSpotToken_;

   // switches
   bool doGenParticles_   ;
   bool runOnParticleGun_ ;

   TTree*         tree_ ;

   //// event variables
   UInt_t     run_    ;
   ULong64_t  event_  ;
   UInt_t     lumis_  ;
   Bool_t     isData_ ;

   //// PileupSummaryInfo
   Int_t          nPUInfo_ ;
   vector<int>    nPU_     ;
   vector<int>    puBX_    ;
   vector<float>  puTrue_  ;

   //// reco::GenParticle
   Int_t          nMC_          ;
   vector<int>    mcPID_        ;
   vector<int>    mcStatus_     ;
   vector<float>  mcVtx_x_      ;
   vector<float>  mcVtx_y_      ;
   vector<float>  mcVtx_z_      ;
   vector<float>  mcPt_         ;
   vector<float>  mcEta_        ;
   vector<float>  mcPhi_        ;
   vector<float>  mcE_          ;
   vector<float>  mcEt_         ;
   vector<float>  mcMass_       ;
   vector<int>    mcParentage_  ;
   vector<int>    mcMomPID_     ;
   vector<float>  mcMomPt_      ;
   vector<float>  mcMomEta_     ;
   vector<float>  mcMomPhi_     ;
   vector<float>  mcMomMass_    ;
   vector<int>    mcGMomPID_    ;
   vector<int>    mcIndex_      ;
   vector<float>  mcCalIsoDR03_ ;
   vector<float>  mcCalIsoDR04_ ;
   vector<float>  mcTrkIsoDR03_ ;
   vector<float>  mcTrkIsoDR04_ ;

   //// reco::Muon
   // by order of appearance in fillMuons
   
   //type flags
   vector<int>  muType_  ;
   vector<int>  muIsPF_  ;
   vector<int>  muIsGlb_ ;
   vector<int>  muIsInn_ ;
   vector<int>  muIsSta_ ;
   //vector<int>  muIsCalo_ ;
   //vector<int>  muIsRPC_ ;

   //kinematics and charge
   vector<float>  muPt_     ;
   vector<float>  muEta_    ;
   vector<float>  muPhi_    ;
   vector<int>    muCharge_ ;

   //# muons looped over
   Int_t nMu_ ;

   //selection+quality flags
   vector<int>  muIsSelected_ ;
   vector<int>  muIsLoose_    ;
   vector<int>  muIsMedium_   ;
   vector<int>  muIsTight_    ;
   vector<int>  muIsGood_     ;
   vector<int>  muIsSoft_     ;
   vector<int>  muIsHighPt_   ;

   //track references go here 
   //must figure out valid way to write

   //best track info
   vector<float>  muBestTrkDxy_   ;
   vector<float>  muBestTrkDz_    ;
   vector<float>  muBestTrkPt_    ;
   vector<float>  muBestTrkPtErr_ ;

   //inner track info
   vector<float>  muInnerDxy_    ;
   vector<float>  muInnerDz_     ;
   vector<int>    muTrkLayers_   ;
   vector<int>    muPixelLayers_ ;  
   vector<int>    muPixelHits_   ;
   vector<float>  muValidFrac_   ;
   vector<int>    muTrkQuality_  ;

   //stand alone track info
   //currently none

   //global track info
   vector<float>  muChi2NDF_  ;
   vector<int>    muMuonHits_ ;

   //other muon info
   vector<int>    muStations_     ;
   vector<float>  muSegCompat_    ;
   vector<float>  muTrkKink_      ;
   vector<float>  muChi2LocalPos_ ;

   //isolation info
   vector<float>  muIsoTrk_   ;
   vector<float>  muPFChIso_  ;
   vector<float>  muPFPhoIso_ ;
   vector<float>  muPFNeuIso_ ;
   vector<float>  muPFPUIso_  ;

   //# selected muons
   Int_t nMuSel_ ;

}; // end class declaration

#endif
