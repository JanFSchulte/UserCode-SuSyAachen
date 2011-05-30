// -*- C++ -*-
//
// Package:    HadronicTree
// Class:      HadronicTree
// 
/**\class HadronicTree HadronicTree.cc NiklasMohr/HadronicTree/src/HadronicTree.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Niklas Mohr,32 4-C02,+41227676330,
//         Created:  Tue Jan  5 13:23:46 CET 2010
// $Id: HadronicTree.cc,v 1.3 2011/05/30 17:58:22 nmohr Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <map>
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
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/JetReco/interface/GenJet.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>

#include <SuSyAachen/DiLeptonHistograms/interface/VertexWeightFunctor.h>

#include <iostream>

//
// class declaration
//

class HadronicTree : public edm::EDAnalyzer {
   public:
      explicit HadronicTree(const edm::ParameterSet&);
      ~HadronicTree();


    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        // ----------member data ---------------------------
        // Switch for debug output
        bool mcInfo;
        edm::InputTag genTag_;
        edm::InputTag genJetTag_;
        edm::InputTag jetTag_;
        edm::InputTag metTag_;
        edm::InputTag pfCandTag_;
        edm::InputTag vxTag_;
        std::vector< std::string > triggers_;
        
        TFile *theFile;

        //Trees
        TTree*  tree;
        float sumPt;
        float ht;
        float met;
        float metPhi;
        float leadJetPt;
        float leadJetPhi;
        float genHt;
        float genMet;
        int nJets;
        int genNJets;
        int nVertices;
        int nPrimaryInts;
        
        float weight;
        VertexWeightFunctor fctVtxWeight_;

        std::map< std::string, int*>  hlTrigDecision_;
        std::map< std::string, float*>  hlTrigPrescale_;
        std::map< std::string, float*> l1TrigPrescale_;
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
HadronicTree::HadronicTree(const edm::ParameterSet& iConfig):
fctVtxWeight_    (iConfig.getParameter<edm::ParameterSet>("vertexWeights") )
{
    //now do what ever initialization is needed
    //Input collections
    genTag_        = iConfig.getParameter<edm::InputTag> ("genSrc");
    genJetTag_     = iConfig.getParameter<edm::InputTag> ("genJets");
    jetTag_          = iConfig.getParameter<edm::InputTag> ("jets");
    metTag_          = iConfig.getParameter<edm::InputTag> ("met");
    pfCandTag_          = iConfig.getParameter<edm::InputTag> ("pfCandidates");
    vxTag_          = iConfig.getParameter<edm::InputTag> ("vertices");
    triggers_          = iConfig.getParameter< std::vector<std::string> > ("triggers");
    
    // Create the root file
    edm::Service<TFileService> theFile;

    TFileDirectory Tree = theFile->mkdir( "Trees" );
    tree = Tree.make<TTree>("MC","MC"); 
    tree->Branch("sumPt",&sumPt,"sumPt/F");
    tree->Branch("ht",&ht,"ht/F");
    tree->Branch("nJets",&nJets,"nJets/I");
    tree->Branch("met",&met,"met/F");
    tree->Branch("metPhi",&metPhi,"metPhi/F");
    tree->Branch("leadJetPt",&leadJetPt,"leadJetPt/F");
    tree->Branch("leadJetPhi",&leadJetPhi,"leadJetPhi/F");
    tree->Branch("genHt",&genHt,"genHt/F");
    tree->Branch("genMet",&genMet,"genMet/F");
    tree->Branch("genNJets",&genNJets,"genNJets/I");
    tree->Branch("nVertices",&nVertices,"nVertices/I");
    tree->Branch("nPrimaryInts",&nPrimaryInts,"nPrimaryInts/I");
    tree->Branch("weight",&weight,"weight/F");
    for ( std::vector<std::string>::iterator trig_i = triggers_.begin(); trig_i != triggers_.end(); ++trig_i ) {
        std::string trigPath = *trig_i;
        hlTrigDecision_[trigPath] = new int;
        l1TrigPrescale_[trigPath] = new float;
        hlTrigPrescale_[trigPath] = new float;
        tree->Branch(("hl_"+trigPath).c_str(),hlTrigDecision_[trigPath],("hl_"+trigPath+"/I").c_str());
        tree->Branch(("l1Prescale_"+trigPath).c_str(),l1TrigPrescale_[trigPath],("l1Prescale_"+trigPath+"/F").c_str());
        tree->Branch(("hlPrescale_"+trigPath).c_str(),hlTrigPrescale_[trigPath],("hlPrescale_"+trigPath+"/F").c_str());
    }
}


HadronicTree::~HadronicTree()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions


// ------------ method called to for each event  ------------
void HadronicTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // does not work ?
    if( iEvent.isRealData() ) mcInfo = false;
    else mcInfo = true;

    //Collection
    edm::Handle< reco::VertexCollection > vertices;
    iEvent.getByLabel(vxTag_, vertices);
    
    edm::Handle< std::vector< reco::PFCandidate > > pfs;
    iEvent.getByLabel(pfCandTag_, pfs);

    edm::Handle< std::vector< pat::Jet > > jets;
    iEvent.getByLabel(jetTag_, jets);

    edm::Handle< std::vector< pat::MET > > mets;
    iEvent.getByLabel(metTag_, mets);
    
    met = mets->front().pt();
    metPhi = mets->front().phi();

    sumPt = 0.0;
    ht = 0.0;
    genHt = 0.;
    genMet = 0.;
    genNJets = 0;
    leadJetPt = 0;
    leadJetPhi = 0;
    nJets = jets->size();
    nVertices = vertices->size();
    nPrimaryInts = -1;
    weight = fctVtxWeight_( iEvent );
    if (nJets > 0) {
        leadJetPt = jets->front().pt();
        leadJetPhi = jets->front().phi();
    }
    for(std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end() ; ++it){
         ht += (*it).pt();
    }
    
    for(std::vector<reco::PFCandidate>::const_iterator it = pfs->begin(); it != pfs->end() ; ++it){
         sumPt += (*it).pt();
    }
   
    // TODO Fill dummy variables with real prescales 
    for ( std::vector<std::string>::iterator trig_i = triggers_.begin(); trig_i != triggers_.end(); ++trig_i ) {
        std::string trigPath = *trig_i;
        *(hlTrigDecision_[trigPath]) = 1;
        *(l1TrigPrescale_[trigPath]) = 5.;
        *(hlTrigPrescale_[trigPath]) = 6.;
    }
    /*unsigned int HLTConfigProvider::prescaleValue(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& trigger
    */
    
    if (!mcInfo) {
        tree->Fill();
        return;
    }
    
    edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

    for(std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

        int BX = PVI->getBunchCrossing();
        if(BX == 0) { 
            nPrimaryInts = PVI->getPU_NumInteractions();
            continue;
        }

    }
    
    edm::Handle< std::vector< reco::GenJet > > genJets;
    iEvent.getByLabel(genJetTag_, genJets);

    edm::Handle< std::vector< reco::GenParticle > > genParticles;
    iEvent.getByLabel(genTag_, genParticles);
    
    genNJets = genJets->size();
    for(std::vector<reco::GenJet>::const_iterator it = genJets->begin(); it != genJets->end() ; ++it){
         genHt += (*it).pt();
    }
    
    int pid = 0;
    double metx = 0.;
    double mety = 0.;
    for (std::vector<reco::GenParticle>::const_iterator p_i = genParticles->begin(); p_i != genParticles->end(); ++p_i){
    pid = abs(p_i->pdgId());
    if (p_i->status()==1){
 	    if ( pid == 12 || pid == 13 || pid == 14 || pid == 16 || 
		 pid == 1000022 || pid == 2000012 || pid == 2000014 ||
		 pid == 2000016 || pid == 1000039 || pid == 5000039 ||
		 pid == 4000012 || pid == 9900012 || pid == 9900014 ||
		 pid == 9900016 || pid == 39 ){
	        metx += p_i->px(); //TODO define met using et
		    mety += p_i->py();
	    }
    }
    }
    genMet = sqrt(metx*metx+mety*mety);

    tree->Fill();

 
}


// ------------ method called once each job just before starting event loop  ------------
void HadronicTree::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void HadronicTree::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(HadronicTree);
