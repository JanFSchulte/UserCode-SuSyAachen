// -*- C++ -*-
//
// Package:    IsoTreeWriter
// Class:      IsoTreeWriter
// 
/**\class IsoTreeWriter IsoTreeWriter.cc NiklasMohr/IsoTreeWriter/src/IsoTreeWriter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Niklas Mohr,32 4-C02,+41227676330,
//         Created:  Tue Jan  5 13:23:46 CET 2010
// $Id: IsoTreeWriter.cc,v 1.7 2011/02/10 15:44:34 edelhoff Exp $
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

#include "SuSyAachen/TagAndProbeTreeWriter/interface/LeptonKindFunctor.h"
#include "SuSyAachen/TagAndProbeTreeWriter/interface/IsoTreeTauExtensions.h"

#include <iostream>

//
// class declaration
//

template< typename T > 
class IsoTreeWriter : public edm::EDAnalyzer {
public:
	explicit IsoTreeWriter(const edm::ParameterSet&);
	~IsoTreeWriter();


private:
	virtual void beginJob() ;
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;

	virtual void fillIso(const edm::Handle< std::vector<T> >&);
	virtual void fillExtraVars(const pat::Electron&);
	virtual void fillExtraVars(const pat::Muon&);
	virtual void fillExtraVars(const pat::Tau&);
	//        virtual double calcIso(const T &);
	virtual double calcIso(const pat::Electron &);
	virtual double calcIsoMinPt(const pat::Electron &);
	virtual double calcPfIso(const T &);

	// ----------member data ---------------------------
	// Switch for debug output
	bool mcInfo;
	bool tauExtensionsActive_;

	edm::InputTag leptonSrc;
	edm::InputTag jetTag_;
	edm::InputTag metTag_;

	TFile *theFile;

	LeptonKindFunctor fctLeptonKind_;

	//Trees
	TTree*  treeIso;
	TTree*  treeIsoEvent;
	float pfIso;
	float mva;
	float iso;
	float isoMinPt;
	float pt;
	float eta;
	float tauDiscr;
	float ht;
	float met;
	int nLept;
	int nJets;
	int leptonKind;

	//Extensions
	IsoTreeTauExtensions tauExtensions_;

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
IsoTreeWriter<T>::IsoTreeWriter(const edm::ParameterSet& iConfig)
{
	//now do what ever initialization is needed
	tauExtensionsActive_ = true;

	//Input collections
	leptonSrc          = iConfig.getParameter<edm::InputTag> ("src");
	jetTag_          = iConfig.getParameter<edm::InputTag> ("jets");
	metTag_          = iConfig.getParameter<edm::InputTag> ("met");

	// Create the root file
	edm::Service<TFileService> theFile;

	TFileDirectory Tree = theFile->mkdir( "Trees" );
	treeIso = Tree.make<TTree>("Iso","Iso");
	treeIso->Branch("pfIso",&pfIso,"pfIso/F");
	treeIso->Branch("mva",&mva,"mva/F");
	treeIso->Branch("iso",&iso,"iso/F");
	treeIso->Branch("isoMinPt",&isoMinPt,"isoMinPt/F");
	treeIso->Branch("pt",&pt,"pt/F");
	treeIso->Branch("eta",&eta,"eta/F");
	treeIso->Branch("tauDiscr",&tauDiscr,"tauDiscr/F");
	treeIso->Branch("ht",&ht,"ht/F");
	treeIso->Branch("met",&met,"met/F");
	treeIso->Branch("nLept",&nLept,"nLept/I");
	treeIso->Branch("leptonKind",&leptonKind,"leptonsKind/I");

	if(tauExtensionsActive_) tauExtensions_.init(iConfig, *treeIso);

}


template< typename T  > 
IsoTreeWriter<T>::~IsoTreeWriter()
{
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}

template < typename T >
double IsoTreeWriter<T>::calcIso(const pat::Electron& lepton)
{
	double value = -1.0;
	if (lepton.eta() <= 1.479)
		value = (lepton.dr03HcalTowerSumEt()+lepton.dr03EcalRecHitSumEt()+lepton.dr03TkSumPt())/lepton.pt();
	else
		value = (lepton.dr03HcalTowerSumEt() + std::max(0.0, lepton.dr03EcalRecHitSumEt() - 1.0) + lepton.dr03TkSumPt())/lepton.pt();

	return value;
}

template < typename T >
double IsoTreeWriter<T>::calcIsoMinPt(const pat::Electron& lepton)
{
	double value = -1.0;
	if (lepton.eta() <= 1.479)
		value = (lepton.dr03HcalTowerSumEt()+lepton.dr03EcalRecHitSumEt()+lepton.dr03TkSumPt())/std::max(lepton.pt(), 20.0);
	else
		value = (lepton.dr03HcalTowerSumEt() + std::max(0.0, lepton.dr03EcalRecHitSumEt() - 1.0) + lepton.dr03TkSumPt())/std::max(lepton.pt(), 20.0);

	return value;
}

template < typename T >
double IsoTreeWriter<T>::calcPfIso(const T& lepton)
{
	double value = (lepton.chargedHadronIso()+lepton.photonIso()+lepton.neutralHadronIso())/lepton.pt();
	return value;
}

template< typename T > 
void IsoTreeWriter<T>::fillExtraVars(const pat::Electron& lepton)
{
	tauDiscr = -1.0;
	mva = lepton.pfCandidateRef()->mva_e_pi();
	iso = calcIso(lepton);
	isoMinPt = calcIsoMinPt(lepton);
}

template< typename T > 
void IsoTreeWriter<T>::fillExtraVars(const pat::Muon& lepton)
{
	tauDiscr = -1.0;
	mva = -1.0;
	iso = -1.0;
	isoMinPt = -1.0;
}

template< typename T > 
void IsoTreeWriter<T>::fillExtraVars(const pat::Tau& lepton)
{
	tauDiscr = lepton.tauID("byTaNCfrOnePercent");
	mva = lepton.tauID("byTaNC");
	iso = -1.0;
	isoMinPt = -1.0;
	if(tauExtensionsActive_) tauExtensions_.fill(*treeIso, lepton);
}


//
// member functions
template< typename T > 
void IsoTreeWriter<T>::fillIso(const edm::Handle< std::vector<T> >& leptons)
{
	nLept = 0;
	for (typename std::vector<T>::const_iterator lep_i = leptons->begin(); lep_i != leptons->end(); ++lep_i){
		nLept = leptons->size();
		pfIso = calcPfIso(*lep_i);
		pt = lep_i->pt();
		eta = lep_i->eta();
		leptonKind = fctLeptonKind_(*lep_i);
		fillExtraVars(*lep_i);
		treeIso->Fill();
	}
}


// ------------ method called to for each event  ------------
template< typename T  > 
void IsoTreeWriter<T >::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	// does not work ?
	if (iEvent.isRealData())
		mcInfo = false;

	//Collection
	edm::Handle< std::vector<T> > leptons;
	iEvent.getByLabel(leptonSrc, leptons);

	edm::Handle< std::vector< pat::Jet > > jets;
	iEvent.getByLabel(jetTag_, jets);

	edm::Handle< std::vector< pat::MET > > mets;
	iEvent.getByLabel(metTag_, mets);

	met = mets->front().pt();

	ht = 0.0;
	for(std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end() ; ++it){
		ht += (*it).pt();
	}

	//Probes
	//run the TnP
	fillIso(leptons);

}


// ------------ method called once each job just before starting event loop  ------------
template< typename T > 
void IsoTreeWriter<T>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template< typename T > 
void IsoTreeWriter<T>::endJob() 
{
}

//define this as a plug-in
typedef IsoTreeWriter< pat::Muon > MuonIsoTreeWriter;
typedef IsoTreeWriter< pat::Electron > ElectronIsoTreeWriter;
typedef IsoTreeWriter< pat::Tau > TauIsoTreeWriter;
DEFINE_FWK_MODULE(MuonIsoTreeWriter);
DEFINE_FWK_MODULE(ElectronIsoTreeWriter);
DEFINE_FWK_MODULE(TauIsoTreeWriter);
