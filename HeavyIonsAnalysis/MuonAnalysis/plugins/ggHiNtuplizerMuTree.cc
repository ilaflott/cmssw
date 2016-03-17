#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "HeavyIonsAnalysis/MuonAnalysis/interface/GenParticleParentageMuTree.h"

#include "HeavyIonsAnalysis/MuonAnalysis/interface/ggHiNtuplizerMuTree.h"


using namespace std;
//using namespace genpartparentage;

ggHiNtuplizerMuTree::ggHiNtuplizerMuTree(const edm::ParameterSet& ps)//:
{

  // class instance configuration
  doGenParticles_         = ps.getParameter<bool>("doGenParticles");
  runOnParticleGun_       = ps.getParameter<bool>("runOnParticleGun");
  genPileupCollection_    = consumes<vector<PileupSummaryInfo> >   (ps.getParameter<edm::InputTag>("pileupCollection"));
  genParticlesCollection_ = consumes<vector<reco::GenParticle> >   (ps.getParameter<edm::InputTag>("genParticleSrc"));
  recoMuonsCollection_    = consumes<edm::View<reco::Muon> >       (ps.getParameter<edm::InputTag>("recoMuonSrc"));
  vtxCollection_          = consumes<vector<reco::Vertex> >        (ps.getParameter<edm::InputTag>("VtxLabel"));
  //pfCollection_           = consumes<edm::View<reco::PFCandidate> > (ps.getParameter<edm::InputTag>("particleFlowCollection"));
  //beamSpotToken_          = consumes<reco::BeamSpot>(ps.getParameter <edm::InputTag>("beamSpot"));

  // initialize output TTree
  edm::Service<TFileService> fs;

  tree_ = fs->make<TTree>("EventTree", "Event data");

  //event variables
  tree_->Branch("run",    &run_);
  tree_->Branch("event",  &event_);
  tree_->Branch("lumis",  &lumis_);
  tree_->Branch("isData", &isData_);
  
  //mc variables
  if (doGenParticles_) {
    tree_->Branch("nPUInfo",      &nPUInfo_);
    tree_->Branch("nPU",          &nPU_);
    tree_->Branch("puBX",         &puBX_);
    tree_->Branch("puTrue",       &puTrue_);

    tree_->Branch("nMC",          &nMC_);
    tree_->Branch("mcPID",        &mcPID_);
    tree_->Branch("mcStatus",     &mcStatus_);
    tree_->Branch("mcVtx_x",      &mcVtx_x_);
    tree_->Branch("mcVtx_y",      &mcVtx_y_);
    tree_->Branch("mcVtx_z",      &mcVtx_z_);
    tree_->Branch("mcPt",         &mcPt_);
    tree_->Branch("mcEta",        &mcEta_);
    tree_->Branch("mcPhi",        &mcPhi_);
    tree_->Branch("mcE",          &mcE_);
    tree_->Branch("mcEt",         &mcEt_);
    tree_->Branch("mcMass",       &mcMass_);
    tree_->Branch("mcParentage",  &mcParentage_);
    tree_->Branch("mcMomPID",     &mcMomPID_);
    tree_->Branch("mcMomPt",      &mcMomPt_);
    tree_->Branch("mcMomEta",     &mcMomEta_);
    tree_->Branch("mcMomPhi",     &mcMomPhi_);
    tree_->Branch("mcMomMass",    &mcMomMass_);
    tree_->Branch("mcGMomPID",    &mcGMomPID_);
    tree_->Branch("mcIndex",      &mcIndex_);
    tree_->Branch("mcCalIsoDR03", &mcCalIsoDR03_);
    tree_->Branch("mcCalIsoDR04", &mcCalIsoDR04_);
    tree_->Branch("mcTrkIsoDR03", &mcTrkIsoDR03_);
    tree_->Branch("mcTrkIsoDR04", &mcTrkIsoDR04_);
  }
  
  //muon variables
  tree_->Branch("nMu",                   &nMu_);
  tree_->Branch("muPt",                  &muPt_);
  tree_->Branch("muEta",                 &muEta_);
  tree_->Branch("muPhi",                 &muPhi_);
  tree_->Branch("muCharge",              &muCharge_);
  tree_->Branch("muType",                &muType_);
  tree_->Branch("muIsGood",              &muIsGood_);
  tree_->Branch("muD0",                  &muD0_);
  tree_->Branch("muDz",                  &muDz_);
  tree_->Branch("muChi2NDF",             &muChi2NDF_);
  tree_->Branch("muInnerD0",             &muInnerD0_);
  tree_->Branch("muInnerDz",             &muInnerDz_);
  tree_->Branch("muTrkLayers",           &muTrkLayers_);
  tree_->Branch("muPixelLayers",         &muPixelLayers_);
  tree_->Branch("muPixelHits",           &muPixelHits_);
  tree_->Branch("muMuonHits",            &muMuonHits_);
  tree_->Branch("muTrkQuality",          &muTrkQuality_);
  tree_->Branch("muStations",            &muStations_);
  tree_->Branch("muIsoTrk",              &muIsoTrk_);
  tree_->Branch("muPFChIso",             &muPFChIso_);
  tree_->Branch("muPFPhoIso",            &muPFPhoIso_);
  tree_->Branch("muPFNeuIso",            &muPFNeuIso_);
  tree_->Branch("muPFPUIso",             &muPFPUIso_);

} // constructor

// event analyzer

void ggHiNtuplizerMuTree::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  //clean up event stuff
  run_    = e.id().run();
  event_  = e.id().event();
  lumis_  = e.luminosityBlock();
  isData_ = e.isRealData();

  //clean up gen pileup stuff
  nPUInfo_ = 0;
  nPU_                  .clear();
  puBX_                 .clear();
  puTrue_               .clear();
  
  //clean up MC+gen particle stuff
  nMC_ = 0;
  mcPID_                .clear();
  mcStatus_             .clear();
  mcVtx_x_              .clear();
  mcVtx_y_              .clear();
  mcVtx_z_              .clear();
  mcPt_                 .clear();
  mcEta_                .clear();
  mcPhi_                .clear();
  mcE_                  .clear();
  mcEt_                 .clear();
  mcMass_               .clear();
  mcParentage_          .clear();
  mcMomPID_             .clear();
  mcMomPt_              .clear();
  mcMomEta_             .clear();
  mcMomPhi_             .clear();
  mcMomMass_            .clear();
  mcGMomPID_            .clear();
  mcIndex_              .clear();
  mcCalIsoDR03_         .clear();
  mcCalIsoDR04_         .clear();
  mcTrkIsoDR03_         .clear();
  mcTrkIsoDR04_         .clear();

  //clean up muon stuff

  //type flags
  muType_               .clear();
  muIsPF_               .clear();
  muIsGlb_              .clear();
  muIsInn_              .clear();
  muIsSta_              .clear();
  
  //kinematics
  muPt_                 .clear();
  muEta_                .clear();
  muPhi_                .clear();
  muCharge_             .clear();
  
  //muons looped over
  nMu_ = 0;

  muIsSelected_         .clear();

  //mu quality
  muIsLoose_            .clear();
  muIsMedium_           .clear();
  muIsTight_            .clear();
  muIsGood_             .clear();
  //  muIsSoft_               .clear();
  //  muIsHighPt_             .clear();

  //best track info
  muD0_                 .clear();
  muDz_                 .clear();

  //other track info
  muChi2NDF_            .clear();
  muInnerD0_            .clear();
  muInnerDz_            .clear();
  muTrkLayers_          .clear();
  muPixelLayers_        .clear();
  muPixelHits_          .clear();
  muMuonHits_           .clear();
  muTrkQuality_         .clear();
  muStations_           .clear();

  //isolation info
  muIsoTrk_             .clear();
  muPFChIso_            .clear();
  muPFPhoIso_           .clear();
  muPFNeuIso_           .clear();
  muPFPUIso_            .clear();

  //selected muons looped over
  nMuSel_ = 0;
  //done cleaning

  // MC truth
  if (doGenParticles_ && !isData_) {
    fillGenPileupInfo(e);
    fillGenParticles(e);
  }

  //void ggHiNtuplizerMuTree::analyze(const edm::Event& e, const edm::EventSetup& es)

  edm::Handle<vector<reco::Vertex> > vtxHandle;
  e.getByToken(vtxCollection_, vtxHandle);

  reco::Vertex& vtx;
  math::XYZPoint pv(0, 0, 0);  

  // grab best primary vertex coordinates
  for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v)
    if (!v->isFake()) {
      pv.SetXYZ(v->x(), v->y(), v->z());
      vtx = (reco::Vertex&) v;//attempting typecast of reco::vertex iterator, may not work
      break;
    }

  //fillMuons(e, es, pv  ); //pv xyz coordinate version
  fillMuons(e, es, (const) vtx, (const) pv ); //const vtx reference version, may not work
  //does passing the es reference do anything for the routine? consider removing
  //can i typecast the input to the fillmuons routine as a constant? may not work

  tree_->Fill();

} // analyze

// tree fill routines

void ggHiNtuplizerMuTree::fillGenPileupInfo(const edm::Event& e)
{
  // Fills information about pileup from MC truth.

  edm::Handle<vector<PileupSummaryInfo> > genPileupHandle;
  e.getByToken(genPileupCollection_, genPileupHandle);

  for (vector<PileupSummaryInfo>::const_iterator pu = genPileupHandle->begin(); pu != genPileupHandle->end(); ++pu) {
    nPU_   .push_back(pu->getPU_NumInteractions());
    puTrue_.push_back(pu->getTrueNumInteractions());
    puBX_  .push_back(pu->getBunchCrossing());

    nPUInfo_++;
  }

} // fillGenPileupInfo

void ggHiNtuplizerMuTree::fillGenParticles(const edm::Event& e)
{
  // Fills tree branches with generated particle info.

  edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
  e.getByToken(genParticlesCollection_, genParticlesHandle);

  int genIndex = 0;

  // loop over MC particles
  for (vector<reco::GenParticle>::const_iterator p = genParticlesHandle->begin(); p != genParticlesHandle->end(); ++p) {
    genIndex++;

    // skip all primary particles if not particle gun MC
    if (!runOnParticleGun_ && !p->mother()) continue;

    // stable particles with pT > 5 GeV
    bool isStableFast = (p->status() == 1 && p->pt() > 5.0);

    // stable leptons
    bool isStableLepton = (p->status() == 1 && abs(p->pdgId()) >= 11 && abs(p->pdgId()) <= 16);

    // (unstable) Z, W, H, top, bottom
    bool isHeavy = (p->pdgId() == 23 || abs(p->pdgId()) == 24 || p->pdgId() == 25 ||
		    abs(p->pdgId()) == 6 || abs(p->pdgId()) == 5);

    // reduce size of output root file
    if (!isStableFast && !isStableLepton && !isHeavy)
      continue;

    mcPID_   .push_back(p->pdgId());
    mcStatus_.push_back(p->status());
    mcVtx_x_ .push_back(p->vx());
    mcVtx_y_ .push_back(p->vy());
    mcVtx_z_ .push_back(p->vz());
    mcPt_    .push_back(p->pt());
    mcEta_   .push_back(p->eta());
    mcPhi_   .push_back(p->phi());
    mcE_     .push_back(p->energy());
    mcEt_    .push_back(p->et());
    mcMass_  .push_back(p->mass());

    reco::GenParticleRef partRef = reco::GenParticleRef(
      genParticlesHandle, p - genParticlesHandle->begin());
    genpartparentage::GenParticleParentageMuTree particleHistory(partRef);

    mcParentage_.push_back(particleHistory.hasLeptonParent()*16   +
			   particleHistory.hasBosonParent()*8     +
			   particleHistory.hasNonPromptParent()*4 +
			   particleHistory.hasQCDParent()*2       +
			   particleHistory.hasExoticParent());

    int   momPID  = -999;
    float momPt   = -999;
    float momEta  = -999;
    float momPhi  = -999;
    float momMass = -999;
    int   gmomPID = -999;

    if (particleHistory.hasRealParent()) {
      reco::GenParticleRef momRef = particleHistory.parent();

      // mother
      if (momRef.isNonnull() && momRef.isAvailable()) {
	momPID  = momRef->pdgId();
	momPt   = momRef->pt();
	momEta  = momRef->eta();
	momPhi  = momRef->phi();
	momMass = momRef->mass();

	// granny
	genpartparentage::GenParticleParentageMuTree motherParticle(momRef);
	if (motherParticle.hasRealParent()) {
	  reco::GenParticleRef granny = motherParticle.parent();
	  gmomPID = granny->pdgId();
	}
      }
    }

    mcMomPID_ .push_back(momPID);
    mcMomPt_  .push_back(momPt);
    mcMomEta_ .push_back(momEta);
    mcMomPhi_ .push_back(momPhi);
    mcMomMass_.push_back(momMass);
    mcGMomPID_.push_back(gmomPID);

    mcIndex_  .push_back(genIndex - 1);

    mcCalIsoDR03_.push_back(getGenCalIso(genParticlesHandle, p, 0.3, false, false));
    mcCalIsoDR04_.push_back(getGenCalIso(genParticlesHandle, p, 0.4, false, false));
    mcTrkIsoDR03_.push_back(getGenTrkIso(genParticlesHandle, p, 0.3) );
    mcTrkIsoDR04_.push_back(getGenTrkIso(genParticlesHandle, p, 0.4) );

    nMC_++;

  } // gen-level particles loop

}// fillGenParticles

//void ggHiNtuplizerMuTree::fillMuons(const edm::Event& e, const edm::EventSetup& es, 
void ggHiNtuplizerMuTree::fillMuons(const edm::Event& e, const edm::EventSetup& es, const reco::Vertex& vtx, const math::XYZPoint& pv)
{
  // Fills tree branches with reco muons.

  edm::Handle<edm::View<reco::Muon> > recoMuonsHandle;
  e.getByToken(recoMuonsCollection_, recoMuonsHandle);

  for (edm::View<reco::Muon>::const_iterator mu = recoMuonsHandle->begin(); mu != recoMuonsHandle->end(); ++mu) {
    
    // fill type + kinematic variables before seleciton happens for later evaluation
    // type flags
    muType_  .push_back( mu->type()             );// still not sure what this variable means
    muIsPF_  .push_back( (int) mu->isPFMuon()         );
    muIsGlb_ .push_back( (int) mu->isGlobalMuon()     );
    muIsInn_ .push_back( (int) mu->isTrackerMuon()    );
    muIsSta_ .push_back( (int) mu->isStandAloneMuon() );

    //kinematics and charge
    muPt_    .push_back( mu->pt()     );
    muEta_   .push_back( mu->eta()    );
    muPhi_   .push_back( mu->phi()    );
    muCharge_.push_back( mu->charge() );
    
    nMu_++;//iterate total number of muons passed over

    //implement muon selection here
    bool isSelected =  
      ( mu->isPFMuon() || mu->isGlobalMuon() || mu->isTrackerMuon() || mu->isStandAlone() ) // type selection
      && ( mu->pt()>0 ) ; // kinematic selection
                         
    muIsSelected_ .push_back( (int) isSelected );
    
    if( !muIsSelected ) continue; //if it doesn't meet the basic criteria, don't write the more computationally expensive info
    else std::cout << "hello world, I have found a muon" << std::endl;
    
    ////////////////
    // muon flags //
    ////////////////
    
    // main quality selection flags
    muIsLoose_   .push_back( (int) mu->isLooseMuon()  );
    muIsMedium_  .push_back( (int) mu->isMediumMuon() );
    muIsTight_   .push_back( (int) muon::isTightMuon(*mu, vtx)  );//needs a vertex input, see MuonSelectors.h, may not work

    // secondary quality selection flags
    muIsGood_   .push_back( (int) muon::isGoodMuon( *mu, muon::selectionTypeFromString("TMOneStationTight") ) );
    //    muIsSoft_   .push_back( (int) mu->isSoftMuon()   );//needs a vertex input, see MuonSelectors.h
    //    muIsHighPt_ .push_back( (int) mu->isHighPtMuon() );//needs a vertex input, see MuonSelectors.h
      
    /////////////////////////////
    // track and other mu info //
    /////////////////////////////

    // best track impact parameters
    const reco::TrackRef bestTrack = mu->muonBestTrack(); 
      
    muD0_ .push_back( bestTrack->dxy(pv) );
    muDz_ .push_back( bestTrack->dz(pv)  );

    //grab tracks associated with the muon
    const reco::TrackRef innMu = mu->innerTrack();
    const reco::TrackRef staMu = mu->outerTrack();
    const reco::TrackRef glbMu = mu->globalTrack();

    int garbage;// temporary place holder for compilation + cmsRun tests
    
    // if no inner track, fill inner track variables with placeholders
    if ( innMu.isNull() ) 
      {
      muInnerD0_     .push_back(-99);
      muInnerDz_     .push_back(-99);
      muTrkLayers_   .push_back(-99);
      muPixelLayers_ .push_back(-99);
      muPixelHits_   .push_back(-99);
      muTrkQuality_  .push_back(-99);
      } else // innMu specific info
      {
      muInnerD0_     .push_back( innMu->dxy(pv)                                     );
      muInnerDz_     .push_back( innMu->dz(pv)                                      );
      muTrkLayers_   .push_back( innMu->hitPattern().trackerLayersWithMeasurement() );
      muPixelLayers_ .push_back( innMu->hitPattern().pixelLayersWithMeasurement()   );
      muPixelHits_   .push_back( innMu->hitPattern().numberOfValidPixelHits()       );
      muTrkQuality_  .push_back( innMu->quality(reco::TrackBase::highPurity)        );
      }
    
    // if no outer track, fill outer track variables with placeholders
    if ( staMu.isNull() ) 
      {
	garbarge = 0; 
      } else //staMu specific info, if anything
      {
	garbarge = 1; 
      }

    // if no global track, fill global track variables with placeholders
    if ( glbMu.isNull() ) 
      {
      muChi2NDF_  .push_back(-99);
      muMuonHits_ .push_back(-99);
      } else // glbMu specific info
      {
      muChi2NDF_  .push_back( glbMu->normalizedChi2()                     );
      muMuonHits_ .push_back( glbMu->hitPattern().numberOfValidMuonHits() );
      }

    //other info to keep
    muStations_ .push_back(mu->numberOfMatchedStations());
    muIsoTrk_   .push_back(mu->isolationR03().sumPt);
    muPFChIso_  .push_back(mu->pfIsolationR04().sumChargedHadronPt);
    muPFPhoIso_ .push_back(mu->pfIsolationR04().sumPhotonEt);
    muPFNeuIso_ .push_back(mu->pfIsolationR04().sumNeutralHadronEt);
    muPFPUIso_  .push_back(mu->pfIsolationR04().sumPUPt);

    //iterate number of muons that pass selections
    nMuSel_++;
  } // muons loop

} //fillMuons

// end tree fill routines

// helper functions

float ggHiNtuplizerMuTree::getGenCalIso(edm::Handle<vector<reco::GenParticle> > &handle,
				  reco::GenParticleCollection::const_iterator thisPart,
				  float dRMax, bool removeMu, bool removeNu)
{
  // Returns Et sum.

  float etSum = 0;

  for (reco::GenParticleCollection::const_iterator p = handle->begin(); p != handle->end(); ++p) {
    if (p == thisPart) continue;
    if (p->status() != 1) continue;

    // has to come from the same collision
    if (thisPart->collisionId() != p->collisionId())
      continue;

    int pdgCode = abs(p->pdgId());

    // skip muons/neutrinos, if requested
    if (removeMu && pdgCode == 13) continue;
    if (removeNu && (pdgCode == 12 || pdgCode == 14 || pdgCode == 16)) continue;

    // must be within deltaR cone
    float dR = reco::deltaR(thisPart->momentum(), p->momentum());
    if (dR > dRMax) continue;

    etSum += p->et();
  }

  return etSum;
} // getGenCalIso

float ggHiNtuplizerMuTree::getGenTrkIso(edm::Handle<vector<reco::GenParticle> > &handle,
				  reco::GenParticleCollection::const_iterator thisPart, float dRMax)
{
  // Returns pT sum without counting neutral particles.

  float ptSum = 0;

  for (reco::GenParticleCollection::const_iterator p = handle->begin(); p != handle->end(); ++p) {
    if (p == thisPart) continue;
    if (p->status() != 1) continue;
    if (p->charge() == 0) continue;  // do not count neutral particles

    // has to come from the same collision
    if (thisPart->collisionId() != p->collisionId())
      continue;

    // must be within deltaR cone
    float dR = reco::deltaR(thisPart->momentum(), p->momentum());
    if (dR > dRMax) continue;

    ptSum += p->pt();
  }

  return ptSum;
} // getGenTrkIso

// end helper functions

DEFINE_FWK_MODULE(ggHiNtuplizerMuTree);
