// -*- C++ -*-
//
// Package:    DummyTrees
// Class:      DummyTrees
// 
/**\class DummyTrees DummyTrees.cc SusyAachen/DiLepton/src/DummyTrees.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Niklas Mohr,32 4-C02,+41227676330,
//         Created:  Tue Jan  5 13:23:46 CET 2010
// $Id: DummyTrees.cc,v 1.1 2010/10/11 08:58:50 nmohr Exp $
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
//
// class declaration
//

class DummyTrees : public edm::EDAnalyzer {
   public:
      explicit DummyTrees(const edm::ParameterSet&);
      ~DummyTrees();

    private:
        typedef reco::Candidate::LorentzVector LorentzVector;
        
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        virtual void fillVectors(const edm::Handle< std::vector<pat::Electron> >&,const edm::Handle< std::vector<pat::Muon> >&,const edm::Handle< std::vector<pat::Tau> >&,const edm::Handle< std::vector<pat::Jet> >&,const edm::Handle< std::vector<pat::MET> >&);

        // ----------member data ---------------------------
        // Switch for debug output
        bool mcInfo;
        edm::InputTag elecSrc;
        edm::InputTag muonSrc;
        edm::InputTag tauSrc;
        edm::InputTag jetSrc;
        edm::InputTag metSrc;

        std::string bJetAlgo;
        double cut_bTagDiscriminator; 

        //File
        TFile *theFile;

        //Trees
        TTree*  dummyTree;
        //Objects
        std::vector<LorentzVector> *theObjects;
        std::vector<int> *theIDs;
        int runNr;
        int eventNr;
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
DummyTrees::DummyTrees(const edm::ParameterSet& iConfig)

{
    //now do what ever initialization is needed
    //Debug flag
    mcInfo      = iConfig.getUntrackedParameter<bool> ("mcInfoAvailable",false);
    //Input collections
    elecSrc     = iConfig.getParameter<edm::InputTag> ("electrons");
    muonSrc     = iConfig.getParameter<edm::InputTag> ("muons");
    tauSrc      = iConfig.getParameter<edm::InputTag> ("taus");
    jetSrc      = iConfig.getParameter<edm::InputTag> ("jets");
    metSrc      = iConfig.getParameter<edm::InputTag> ("met");
    bJetAlgo    = iConfig.getUntrackedParameter<std::string> ("bJetAlgo");
    cut_bTagDiscriminator = iConfig.getUntrackedParameter<double> ("bTagDiscriminator");
    
    // Create the root file
    edm::Service<TFileService> theFile;

    TFileDirectory Tree = theFile->mkdir( "Trees" );
    dummyTree = Tree.make<TTree>("Events","Events");

    theObjects = new std::vector<LorentzVector>;
    theIDs = new std::vector<int>;

    dummyTree->Branch("objects","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&theObjects);
    dummyTree->Branch("ids","std::vector<int>",&theIDs);
    dummyTree->Branch("runNr",&runNr,"runNr/I");
    dummyTree->Branch("eventNr",&eventNr,"eventNr/I");

}

DummyTrees::~DummyTrees()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

void DummyTrees::fillVectors(const edm::Handle< std::vector<pat::Electron> >& elecs,const edm::Handle< std::vector<pat::Muon> >& muons,const edm::Handle< std::vector<pat::Tau> >& taus,const edm::Handle< std::vector<pat::Jet> >& jets,const edm::Handle< std::vector<pat::MET> >& mets){
    for (std::vector<pat::Electron>::const_iterator ele_i = elecs->begin(); ele_i != elecs->end(); ++ele_i){
        theObjects->push_back(ele_i->p4());
        theIDs->push_back(11*ele_i->charge());
    }
    for (std::vector<pat::Muon>::const_iterator mu_i = muons->begin(); mu_i != muons->end(); ++mu_i){
        theObjects->push_back(mu_i->p4());
        theIDs->push_back(13*mu_i->charge());
    }
    for (std::vector<pat::Tau>::const_iterator tau_i = taus->begin(); tau_i != taus->end(); ++tau_i){
        theObjects->push_back(tau_i->p4());
        theIDs->push_back(15*tau_i->charge());
    }
    for (std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i){
        theObjects->push_back(jet_i->p4());
        if(jet_i->bDiscriminator(bJetAlgo)>cut_bTagDiscriminator) theIDs->push_back(5);
        else theIDs->push_back(1);
    }
    for (std::vector<pat::MET>::const_iterator met_i = mets->begin(); met_i != mets->end(); ++met_i){
        theObjects->push_back(met_i->p4());
        theIDs->push_back(0);
    }
}

// ------------ method called to for each event  ------------
void DummyTrees::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    if (iEvent.isRealData()) mcInfo = false;
    //Run and event number
    runNr = iEvent.id().run();
    eventNr = iEvent.id().event();

    //Electrons
    edm::Handle< std::vector<pat::Electron> > elecs;
    iEvent.getByLabel(elecSrc, elecs);
  
    //Muons
    edm::Handle< std::vector<pat::Muon> > muons;
    iEvent.getByLabel(muonSrc, muons);
    
    //Taus
    edm::Handle< std::vector<pat::Tau> > taus;
    iEvent.getByLabel(tauSrc, taus);
   
    //Jets
    edm::Handle< std::vector<pat::Jet> > jets;
    iEvent.getByLabel(jetSrc, jets);
    
    //MET
    edm::Handle< std::vector<pat::MET> > mets;
    iEvent.getByLabel(metSrc, mets);

    //Fill vectors
    fillVectors(elecs,muons,taus,jets,mets);
    dummyTree->Fill();

    //Clear contents
    theObjects->clear();
    theIDs->clear();
}


// ------------ method called once each job just before starting event loop  ------------
void DummyTrees::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void DummyTrees::endJob() {
}

DEFINE_FWK_MODULE(DummyTrees);
