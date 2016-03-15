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

   virtual void analyze(const edm::Event&, const edm::EventSetup&);

   void fillMuons        (const edm::Event&, const edm::EventSetup&, math::XYZPoint& pv);
   
   void fillGenParticles (const edm::Event&);
   void fillGenPileupInfo(const edm::Event&);
   
   // Et and pT sums
   float getGenCalIso(edm::Handle<vector<reco::GenParticle> >&, reco::GenParticleCollection::const_iterator, float dRMax, bool removeMu, bool removeNu);
   float getGenTrkIso(edm::Handle<vector<reco::GenParticle> >&, reco::GenParticleCollection::const_iterator, float dRMax);

   // switches
   bool doGenParticles_;
   bool runOnParticleGun_;

   // handles to collections of objects
   edm::EDGetTokenT<vector<PileupSummaryInfo> >       genPileupCollection_;
   edm::EDGetTokenT<vector<reco::GenParticle> >       genParticlesCollection_;
   edm::EDGetTokenT<edm::View<reco::Muon> >           recoMuonsCollection_;
   edm::EDGetTokenT<vector<reco::Vertex> >            vtxCollection_;
   //edm::EDGetTokenT<edm::View<reco::PFCandidate> >    pfCollection_;
   //edm::EDGetTokenT<reco::BeamSpot>                   beamSpotToken_;

   TTree*         tree_;

   // variables associated with tree branches
   UInt_t          run_;
   ULong64_t       event_;
   UInt_t          lumis_;
   Bool_t         isData_;

   // PileupSummaryInfo
   Int_t          nPUInfo_;
   vector<int>    nPU_;
   vector<int>    puBX_;
   vector<float>  puTrue_;

   // reco::GenParticle
   Int_t          nMC_;
   vector<int>    mcPID_;
   vector<int>    mcStatus_;
   vector<float>  mcVtx_x_;
   vector<float>  mcVtx_y_;
   vector<float>  mcVtx_z_;
   vector<float>  mcPt_;
   vector<float>  mcEta_;
   vector<float>  mcPhi_;
   vector<float>  mcE_;
   vector<float>  mcEt_;
   vector<float>  mcMass_;
   vector<int>    mcParentage_;
   vector<int>    mcMomPID_;
   vector<float>  mcMomPt_;
   vector<float>  mcMomEta_;
   vector<float>  mcMomPhi_;
   vector<float>  mcMomMass_;
   vector<int>    mcGMomPID_;
   vector<int>    mcIndex_;
   vector<float>  mcCalIsoDR03_;
   vector<float>  mcCalIsoDR04_;
   vector<float>  mcTrkIsoDR03_;
   vector<float>  mcTrkIsoDR04_;

   // reco::Muon
   Int_t          nMu_;
   vector<float>  muPt_;
   vector<float>  muEta_;
   vector<float>  muPhi_;
   vector<int>    muCharge_;
   vector<int>    muType_;
   vector<int>    muIsGood_;
   vector<float>  muD0_;
   vector<float>  muDz_;
   vector<float>  muChi2NDF_;
   vector<float>  muInnerD0_;
   vector<float>  muInnerDz_;
   vector<int>    muTrkLayers_;
   vector<int>    muPixelLayers_;
   vector<int>    muPixelHits_;
   vector<int>    muMuonHits_;
   vector<int>    muTrkQuality_;
   vector<int>    muStations_;
   vector<float>  muIsoTrk_;
   vector<float>  muPFChIso_;
   vector<float>  muPFPhoIso_;
   vector<float>  muPFNeuIso_;
   vector<float>  muPFPUIso_;

};

#endif
