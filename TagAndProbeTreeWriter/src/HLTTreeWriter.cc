// -*- C++ -*-
//
// Package:    HLTTreeWriter
// Class:      HLTTreeWriter
// 
/**\class HLTTreeWriter HLTTreeWriter.cc NiklasMohr/HLTTreeWriter/src/HLTTreeWriter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Niklas Mohr,32 4-C02,+41227676330,
//         Created:  Tue Jan  5 13:23:46 CET 2010
// $Id: HLTTreeWriter.cc,v 1.7 2011/02/10 15:44:34 edelhoff Exp $
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

#include <iostream>

//
// class declaration
//

class HLTTreeWriter : public edm::EDAnalyzer {
   public:
      explicit HLTTreeWriter(const edm::ParameterSet&);
      ~HLTTreeWriter();


    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        // ----------member data ---------------------------
        // Switch for debug output
        bool mcInfo;
        bool debug_;
        edm::InputTag offlineJetTag_;
        edm::InputTag offlineMetTag_;
        edm::InputTag onlineJetTag_;
        
        TFile *theFile;

        //Trees
        TTree*  treeHLT;
        
        float onlineHt;
        float onlineMHt;
        float met;
        float ht;
        float mht;
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
HLTTreeWriter::HLTTreeWriter(const edm::ParameterSet& iConfig)
{
	debug_ = false;
    if(debug_) std::cout << "init ";
    //now do what ever initialization is needed
    //Input collections
    offlineJetTag_          = iConfig.getParameter<edm::InputTag> ("jets");
    offlineMetTag_          = iConfig.getParameter<edm::InputTag> ("met");
	onlineJetTag_          = iConfig.getParameter<edm::InputTag> ("onlineJets");

    
    // Create the root file
    edm::Service<TFileService> theFile;

    TFileDirectory Tree = theFile->mkdir( "Trees" );
    treeHLT = Tree.make<TTree>("HLT","HLT"); 
    treeHLT->Branch("onlineHt",&onlineHt,"onlineHt/F");
    treeHLT->Branch("onlineMHt",&onlineMHt,"onlineMHt/F");	
    treeHLT->Branch("ht",&ht,"ht/F");
    treeHLT->Branch("met",&met,"met/F");
    treeHLT->Branch("mht",&met,"mht/F");
    if(debug_) std::cout << "Done!"<< std::endl;
}


HLTTreeWriter::~HLTTreeWriter()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


// ------------ method called to for each event  ------------
void HLTTreeWriter::analyze(const edm::Event& iEvent, const edm::EventSetup&)
{
    if(debug_) std::cout << "analyze ";
    // does not work ?
    if (iEvent.isRealData())
      mcInfo = false;
    if(debug_) std::cout << "handles";
    //Collection
    edm::Handle< std::vector< pat::Jet > > jets;
    iEvent.getByLabel(offlineJetTag_, jets);

    edm::Handle< std::vector< pat::MET > > mets;
    iEvent.getByLabel(offlineMetTag_, mets);
    
    edm::Handle< std::vector< pat::Jet > > onlineJets;
    iEvent.getByLabel(onlineJetTag_, onlineJets);
    if(debug_) std::cout << ". met";
    met = mets->front().pt();
    
    if(debug_) std::cout << ". offline";
    ht = 0.0;
    reco::Particle::LorentzVector vOfflineHt(0.,0.,0.,0.);
    for(std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end() ; ++it){
         ht += (*it).pt();
         mht += (*it).pt();
         vOfflineHt += (*it).p4();
	}	
    mht = vOfflineHt.pt();
    if(debug_) std::cout << ". online";
    onlineHt = 0.0;
    reco::Particle::LorentzVector vOnlineHt(0.,0.,0.,0.);
    for(std::vector<pat::Jet>::const_iterator it = onlineJets->begin(); it != onlineJets->end() ; ++it){
         onlineHt += (*it).pt();
	     vOnlineHt += (*it).p4();
    }
    onlineMHt = vOnlineHt.pt();
    if(debug_) std::cout << ". filling";
    treeHLT->Fill();
    if(debug_) std::cout << ". Done!"<< std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void HLTTreeWriter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void HLTTreeWriter::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTTreeWriter);