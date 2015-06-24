// -*- C++ -*-
//
// Package:    IsoTreeWriterMiniAOD
// Class:      IsoTreeWriterMiniAOD
// 
/**\class IsoTreeWriterMiniAOD IsoTreeWriterMiniAOD.cc NiklasMohr/IsoTreeWriterMiniAOD/src/IsoTreeWriterMiniAOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Niklas Mohr,32 4-C02,+41227676330,
//         Created:  Tue Jan  5 13:23:46 CET 2010
// $Id: IsoTreeWriterMiniAOD.cc,v 1.20 2012/05/31 20:54:58 edelhoff Exp $
//
//

// system include files
#include <memory>
#include "TFile.h"
#include "TTree.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "SuSyAachen/TagAndProbeTreeWriter/interface/LeptonKindFunctor.h"
#include "SuSyAachen/TagAndProbeTreeWriter/interface/IsoTreeEventFootprintExtensionsMiniAOD.h"
#include "SuSyAachen/TagAndProbeTreeWriter/interface/IsoTreeTauExtensionsMiniAOD.h"
#include "SuSyAachen/TagAndProbeTreeWriter/interface/IsolationFunctor.h"
#include "SuSyAachen/TagAndProbeTreeWriter/interface/IsoTreeSecondLeptonExtensionsMiniAOD.h"
#include <SuSyAachen/DiLeptonHistograms/interface/VertexWeightFunctor.h>

#include "DataFormats/Math/interface/deltaR.h"

#include <iostream>

//
// class declaration
//

template< typename T > 
class IsoTreeWriterMiniAOD : public edm::EDAnalyzer {
public:
	explicit IsoTreeWriterMiniAOD(const edm::ParameterSet&);
	~IsoTreeWriterMiniAOD();


private:
	virtual void beginJob() ;
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;

  virtual void fillIso(const edm::Handle< std::vector<T> >&, const edm::Event&);
	virtual void fillExtraVars(const pat::Electron&, const edm::Event& iEvent);
	virtual void fillExtraVars(const pat::Muon&, const edm::Event& iEvent);
	virtual void fillExtraVars(const pat::Tau&, const edm::Event& iEvent);
	//        virtual double calcIso(const T &);
	virtual double calcIso(const pat::Electron &);
	virtual double calcIso(const pat::Muon &);
	virtual double calcIsoMinPt(const pat::Electron &);
        int findNextHardId(const T &lepton, const reco::GenParticleCollection &genParticles);

        float transverseMass(const math::XYZTLorentzVector& p, const math::XYZTLorentzVector& met);

	// ----------member data ---------------------------
	// Switch for debug output
        bool debug_;
	bool mcInfo;
        bool eventFootprintExtensionsActive_;
	bool tauExtensionsActive_;
	bool secondLeptonExtensionsActive_;
        bool genParticleActive_;

	edm::InputTag leptonSrc;
	edm::InputTag jetTag_;
	edm::InputTag metTag_;
	edm::InputTag vertexTag_;
        edm::InputTag genTag_;

	TFile *theFile;

	LeptonKindFunctor fctLeptonKind_;
        IsolationFunctor fctIsolation_;
        VertexWeightFunctor fctVtxWeight_;

	//Trees
	TTree*  treeIso;
	TTree*  treeIsoEvent;
        float pfIso;
        float pfIsoAbs;
        float pfIsoAbsChargedHadrons;
        float pfIsoAbsNeutralHadrons;
        float pfIsoAbsPhotons;
	float mva;
	float iso;
	float isoMinPt;
	float pt;
	float eta;
	float tauDiscr;
	float ht;
	float met;
        float mT;
        float ptJet1;
        float ptJet2; 
	int nLept;
	int nJets;
	int nVertices;
	int leptonKind;
        int hardId;
	
        float weight;

	//Extensions
        IsoTreeEventFootprintExtensionsMiniAOD eventFootprintExtensions_;
	IsoTreeTauExtensionsMiniAOD tauExtensions_;
	IsoTreeSecondLeptonExtensionsMiniAOD secondLeptonExtensions_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
//
// constructors and destructor
//
template< typename T  > 
IsoTreeWriterMiniAOD<T>::IsoTreeWriterMiniAOD(const edm::ParameterSet& iConfig):
  fctIsolation_  (iConfig.getParameter<edm::ParameterSet>("isolationDefinitions")),
  fctVtxWeight_    (iConfig.getParameter<edm::ParameterSet>("vertexWeights") )
{
	//now do what ever initialization is needed
        debug_ = true;
        eventFootprintExtensionsActive_ = true;
	tauExtensionsActive_ = false;
	secondLeptonExtensionsActive_ = false;
	mcInfo = false;
	if( iConfig.existsAs<bool>("useEventFootprintExtensions") )
	  eventFootprintExtensionsActive_ = iConfig.getParameter<bool> ("useEventFootprintExtensions");
	if( iConfig.existsAs<bool>("useTauExtensions") )
	  tauExtensionsActive_ = iConfig.getParameter<bool> ("useTauExtensions");
	if( iConfig.existsAs<bool>("useMcInfo")  && iConfig.getParameter<bool> ("useMcInfo"))
	  mcInfo = true;
	

	//Input collections
	leptonSrc          = iConfig.getParameter<edm::InputTag> ("src");
	jetTag_          = iConfig.getParameter<edm::InputTag> ("jets");
	metTag_          = iConfig.getParameter<edm::InputTag> ("met");
	vertexTag_       = iConfig.getParameter<edm::InputTag> ("vertices");

	// Create the root file
	edm::Service<TFileService> theFile;

	TFileDirectory Tree = theFile->mkdir( "Trees" );
	treeIso = Tree.make<TTree>("Iso","Iso");
	treeIso->Branch("pfIso",&pfIso,"pfIso/F");
	treeIso->Branch("pfIsoAbs",&pfIsoAbs,"pfIsoAbs/F");
	treeIso->Branch("pfIsoAbsChargedHadrons",&pfIsoAbsChargedHadrons,"pfIsoAbsChargedHadrons/F");
	treeIso->Branch("pfIsoAbsNeutralHadrons",&pfIsoAbsNeutralHadrons,"pfIsoAbsNeutralHadrons/F");
	treeIso->Branch("pfIsoAbsPhtotons",&pfIsoAbsPhotons,"pfIsoAbsPhotons/F");
	treeIso->Branch("mva",&mva,"mva/F");
	treeIso->Branch("iso",&iso,"iso/F");
	treeIso->Branch("isoMinPt",&isoMinPt,"isoMinPt/F");
	treeIso->Branch("pt",&pt,"pt/F");
	treeIso->Branch("eta",&eta,"eta/F");
	treeIso->Branch("tauDiscr",&tauDiscr,"tauDiscr/F");
	treeIso->Branch("ht",&ht,"ht/F");
        treeIso->Branch("ptJet1",&ptJet1,"ptJet1/F");
        treeIso->Branch("ptJet2",&ptJet2,"ptJet2/F");
	treeIso->Branch("met",&met,"met/F");
	treeIso->Branch("nLept",&nLept,"nLept/I");
	treeIso->Branch("nVertices",&nVertices,"nVertices/I");
	treeIso->Branch("leptonKind",&leptonKind,"leptonsKind/I");
	treeIso->Branch("mT",&mT,"mT/F");
	
    treeIso->Branch("weight",&weight,"weight/F");

	if(tauExtensionsActive_) tauExtensions_.init(iConfig, *treeIso);
	if(eventFootprintExtensionsActive_) eventFootprintExtensions_.init(iConfig, *treeIso);
	if( iConfig.existsAs<edm::InputTag>("secondLeptonElectronSrc") &&
	    iConfig.existsAs<edm::InputTag>("secondLeptonMuonSrc") &&
	    iConfig.existsAs<edm::InputTag>("secondLeptonTauSrc") &&
	    iConfig.existsAs<double>("secondLeptonMinDeltaR") ){
		secondLeptonExtensionsActive_ = true;
		secondLeptonExtensions_.init(iConfig, *treeIso);
	}
	if( iConfig.existsAs<edm::InputTag>("genSrc")  && mcInfo){
	  genParticleActive_ = true;
	  genTag_ =  iConfig.getParameter<edm::InputTag> ("genSrc");
	  treeIso->Branch("hardId",&hardId,"hardId/I");	  
	}else{
	  genParticleActive_ = false;
	  genTag_ = edm::InputTag("empty");
	}
}


template< typename T  > 
IsoTreeWriterMiniAOD<T>::~IsoTreeWriterMiniAOD()
{
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}

template < typename T >
double IsoTreeWriterMiniAOD<T>::calcIso(const pat::Electron& lepton)
{
	double value = -1.0;
	if (lepton.isEB())
		value = (lepton.dr03HcalTowerSumEt()+lepton.dr03EcalRecHitSumEt()+lepton.dr03TkSumPt())/lepton.pt();
	else
		value = (lepton.dr03HcalTowerSumEt() + std::max(0.0, lepton.dr03EcalRecHitSumEt() - 1.0) + lepton.dr03TkSumPt())/lepton.pt();

	return value;
}

template < typename T >
double IsoTreeWriterMiniAOD<T>::calcIso(const pat::Muon& lepton)
{
  double value = (lepton.isolationR03().hadEt + lepton.isolationR03().emEt + lepton.isolationR03().sumPt) / lepton.pt();
  return value;
}

template < typename T >
double IsoTreeWriterMiniAOD<T>::calcIsoMinPt(const pat::Electron& lepton)
{
	double value = -1.0;
	if (lepton.isEB())
		value = (lepton.dr03HcalTowerSumEt()+lepton.dr03EcalRecHitSumEt()+lepton.dr03TkSumPt())/std::max(lepton.pt(), (double) 20.0);
	else
		value = (lepton.dr03HcalTowerSumEt() + std::max(0.0, lepton.dr03EcalRecHitSumEt() - 1.0) + lepton.dr03TkSumPt())/std::max(lepton.pt(), (double) 20.0);

	return value;
}

template< typename T > 
void IsoTreeWriterMiniAOD<T>::fillExtraVars(const pat::Electron& lepton, const edm::Event& iEvent)
{
	tauDiscr = -1.0;
	mva = -1.;
	iso = calcIso(lepton);
	isoMinPt = calcIsoMinPt(lepton);
        edm::Handle<reco::VertexCollection> vertices;
        iEvent.getByLabel(vertexTag_, vertices);
	if(debug_){
	  std::cout << "****Electron***" 
		    << std::endl
		    << "++++> pt "<< lepton.pt() << " eta "<< lepton.eta()
		    << std::endl << "\t\t"                                
		    << " deltaEtaSuperClusterTrackAtVtx " << lepton.deltaEtaSuperClusterTrackAtVtx()
		    << " deltaPhiSuperClusterTrackAtVtx " << lepton.deltaPhiSuperClusterTrackAtVtx()
		    << std::endl << "\t\t" 
		    << " sigmaIetaIeta " << lepton.full5x5_sigmaIetaIeta() 
		    << " hadronicOverEm "<< lepton.hadronicOverEm()
		    << std::endl << "\t\t"
		    << " abs(1.0/ecalEnergy - eSuperClusterOverP/ecalEnergy) " << fabs(1.0/lepton.ecalEnergy() - lepton.eSuperClusterOverP()/lepton.ecalEnergy())
		    << std::endl << "\t\t"	
		    << "d0 " << lepton.gsfTrack()->dxy(vertices->at(0).position()) << std::endl << "\t\t"
		    << "dZ " << fabs(lepton.gsfTrack()->dz(vertices->at(0).position())) << std::endl << "\t\t"
		    << " passes Conversion " << lepton.passConversionVeto()  << std::endl << "\t\t" 
		    << " missing inner hits " << lepton.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)	 << std::endl << "\t\t"	    
		    << " iso " << fctIsolation_(lepton,"miniIsoEA")/lepton.pt() 
	  	    << std::endl << "\t\t"	
	    //		  << std::endl << "\t\t"
	    //		  << " ch "<<lepton.pfIsolationR03().sumChargedHadronPt <<" neut "<< lepton.pfIsolationR03().sumNeutralHadronEt <<" photo "<<  lepton.pfIsolationR03().sumPhotonEt << " pu " <<lepton.pfIsolationR03().sumPUPt       
		    <<std::endl;  
	}
	

}

template< typename T > 
void IsoTreeWriterMiniAOD<T>::fillExtraVars(const pat::Muon& lepton, const edm::Event& iEvent)
{
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vertexTag_, vertices);
  if(debug_){
    std::cout << "*******Muon***********"
	      << "++++> pt "<< lepton.pt() << " eta "<< lepton.eta() 
	      << std::endl << "\t\t"
	      << " emVetoEt " << lepton.isolationR03().emVetoEt
	      << " hadVetoEt " << lepton.isolationR03().hadVetoEt 
	      << std::endl << "\t\t";


  	     if (lepton.isGlobalMuon()){
	     	 std::cout << "Chi^2 " << lepton.globalTrack()->normalizedChi2() << std::endl << "\t\t"
	      	<< "NumberOfValidMuonHits " << lepton.globalTrack()->hitPattern().numberOfValidMuonHits() << std::endl << "\t\t";
	     }
	     else{std::cout  << "Not Global! " << std::endl << "\t\t"; } 
	    std::cout	
	     << "MatchedStations " << lepton.numberOfMatchedStations() << std::endl << "\t\t";
	     if (lepton.isTrackerMuon()){
	     	std::cout << "ValidPixelHits " << lepton.innerTrack()->hitPattern().numberOfValidPixelHits() << std::endl << "\t\t"
			  << "NLayers " << lepton.track()->hitPattern().trackerLayersWithMeasurement() << std::endl << "\t\t";}
	     else{std::cout << "Not Tracker Muon!" << std::endl << "\t\t";}
	     std::cout << "Is PF " << lepton.isPFMuon() << std::endl << "\t\t"
	      << "d0 " << lepton.dB() << std::endl << "\t\t"
	      << "dZ " << fabs(lepton.muonBestTrack()->dz(vertices->at(0).position())) << std::endl << "\t\t"
		    << " iso " << fctIsolation_(lepton,"miniIsoEA")/lepton.pt() 
	      << std::endl;  
  }
  tauDiscr = -1.0;
  mva = -1.0;
  iso = calcIso(lepton);
  isoMinPt = -1.0;
}

template< typename T > 
void IsoTreeWriterMiniAOD<T>::fillExtraVars(const pat::Tau& lepton, const edm::Event& iEvent)
{
  tauDiscr = -1;//lepton.tauID("byTaNCfrOnePercent");
  mva = -1;//lepton.tauID("byTaNC");
	iso = -1.0;
	isoMinPt = -1.0;
	if(tauExtensionsActive_) tauExtensions_.fill(*treeIso, lepton);
}

template< typename T > 
float IsoTreeWriterMiniAOD<T>::transverseMass(const math::XYZTLorentzVector& p, const math::XYZTLorentzVector& met)
{
  reco::Candidate::LorentzVector otherMet(met.Px(),met.Py(),met.Pz(),met.E());
  reco::Candidate::LorentzVector leptonT(p.Px(),p.Py(),0.,p.E()*sin(p.Theta()));
  reco::Candidate::LorentzVector sumT=leptonT+otherMet;

  return std::sqrt(sumT.M2());
}

//
// member functions
template< typename T > 
void IsoTreeWriterMiniAOD<T>::fillIso(const edm::Handle< std::vector<T> >& leptons, const edm::Event& iEvent)
{
	edm::Handle< std::vector< pat::MET > > mets;
	iEvent.getByLabel(metTag_, mets);

	edm::Handle< std::vector< pat::Jet > > jets;
        iEvent.getByLabel(jetTag_, jets);	
	
	fctIsolation_.init(iEvent);


	if(debug_){
		for(std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end() ; ++it){
		  std::cout << "pt: " << (*it).pt() << " eta: " << (*it).eta() << std::endl;

		  }	
	
	}


	nLept = 0;
	edm::Handle<reco::GenParticleCollection> genParticles;
	if(eventFootprintExtensionsActive_) eventFootprintExtensions_.fill(iEvent);
        if(genParticleActive_) iEvent.getByLabel(genTag_, genParticles);
	for (typename std::vector<T>::const_iterator lep_i = leptons->begin(); lep_i != leptons->end(); ++lep_i){
	        ht = 0.0;
		ptJet1 = 0.0;
		ptJet2 = 0.0;
		for(std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end() ; ++it){
		  ht += (*it).pt();
		  if((*it).pt() > ptJet1 &&  reco::deltaR<const T, const pat::Jet>( *lep_i, *it) > 0.5) {
		      ptJet2 = ptJet1;
		      ptJet1 = (*it).pt();
		  }
		}
		
	        nLept = leptons->size();
		pt = lep_i->pt();
		eta = lep_i->eta();
		mT = transverseMass(lep_i->p4(), mets->front().p4());
		leptonKind = fctLeptonKind_(*lep_i);
		fillExtraVars(*lep_i,iEvent);

		// isolation
		pfIsoAbs = fctIsolation_(*lep_i,"pfIsolation");
		pfIso = fctIsolation_(*lep_i,"pfIsolation") / lep_i->pt();
		pfIsoAbsChargedHadrons = -1.;//fctIsolationChargedHadrons_(*lep_i);
		pfIsoAbsNeutralHadrons = -1.;//fctIsolationNeutralHadrons_(*lep_i);
		pfIsoAbsPhotons = -1.;//fctIsolationPhotons_(*lep_i);

		if(secondLeptonExtensionsActive_) secondLeptonExtensions_.fill<T>(*treeIso, *lep_i, iEvent);
		if(genParticleActive_) hardId = findNextHardId(*lep_i, *genParticles);
		treeIso->Fill();
	}
}


// ------------ method called to for each event  ------------
template< typename T  > 
void IsoTreeWriterMiniAOD<T >::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	// does not work ?
  //	if (iEvent.isRealData())
  //		mcInfo = false;

	//Collection
	edm::Handle< std::vector<T> > leptons;
	iEvent.getByLabel(leptonSrc, leptons);

	edm::Handle< std::vector< pat::Jet > > jets;
	iEvent.getByLabel(jetTag_, jets);

	edm::Handle< std::vector< pat::MET > > mets;
	iEvent.getByLabel(metTag_, mets);

	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByLabel(vertexTag_, vertices);

	met = mets->front().pt();

	// count number of vertices
	nVertices = vertices->size();
    weight = fctVtxWeight_( iEvent );

	//Probes
	//run the TnP
	fillIso(leptons, iEvent);
        std::cout << "end event" << std::endl;

}

template< typename T >
int IsoTreeWriterMiniAOD<T>::findNextHardId(const T &lepton, const reco::GenParticleCollection &genParticles)
{
  int result = 0;
  double  minDeltaR = 1e10;
  for(reco::GenParticleCollection::const_iterator it = genParticles.begin(); it != genParticles.end(); ++it){
    if((*it).status() == 3 && reco::deltaR<const T, const reco::GenParticle>( lepton, *it) < minDeltaR){
      minDeltaR = reco::deltaR<const T, const reco::GenParticle>( lepton, *it);
      result = (*it).pdgId();
    }
  }
  return result;
}


// ------------ method called once each job just before starting event loop  ------------
template< typename T > 
void IsoTreeWriterMiniAOD<T>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template< typename T > 
void IsoTreeWriterMiniAOD<T>::endJob() 
{
}

//define this as a plug-in
typedef IsoTreeWriterMiniAOD< pat::Muon > MuonIsoTreeWriterMiniAOD;
typedef IsoTreeWriterMiniAOD< pat::Electron > ElectronIsoTreeWriterMiniAOD;
typedef IsoTreeWriterMiniAOD< pat::Tau > TauIsoTreeWriterMiniAOD;
DEFINE_FWK_MODULE(MuonIsoTreeWriterMiniAOD);
DEFINE_FWK_MODULE(ElectronIsoTreeWriterMiniAOD);
DEFINE_FWK_MODULE(TauIsoTreeWriterMiniAOD);
