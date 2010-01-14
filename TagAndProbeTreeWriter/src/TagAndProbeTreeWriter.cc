// -*- C++ -*-
//
// Package:    TagAndProbeTreeWriter
// Class:      TagAndProbeTreeWriter
// 
/**\class TagAndProbeTreeWriter TagAndProbeTreeWriter.cc NiklasMohr/TagAndProbeTreeWriter/src/TagAndProbeTreeWriter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Niklas Mohr,32 4-C02,+41227676330,
//         Created:  Tue Jan  5 13:23:46 CET 2010
// $Id: TagAndProbeTreeWriter.cc,v 1.1 2010/01/05 15:27:05 nmohr Exp $
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
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
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
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//
// class declaration
//

template< typename T, typename P > 
class TagAndProbeTreeWriter : public edm::EDAnalyzer {
   public:
      explicit TagAndProbeTreeWriter(const edm::ParameterSet&);
      ~TagAndProbeTreeWriter();


    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        virtual void TnP(const edm::Handle< std::vector<T> >&,const edm::Handle< P >&);
        virtual void mcAnalysis(const edm::Handle< std::vector<T> >&,const edm::Handle< std::vector<reco::GenParticle> >&);

        // ----------member data ---------------------------
        // Switch for debug output
        bool mcInfo;
        edm::InputTag mcSrc;
        edm::InputTag tagSrc;
        edm::InputTag probeSrc;
        edm::InputTag jetSrc;

        double cut_Dr;
        double cut_lowInvM;
        double cut_highInvM;
        
        TFile *theFile;

        //Trees
        TTree*  treeTnP;
        TTree*  treeGen;
        TTree*  treeMatch;
        float ptGen;
        float ptMatch;
        float etaGen;
        float etaMatch;
        int pdgIdLepton;
        int motherPdgIdLepton;
        float invM;
        float ptProbe;
        float etaProbe;
        int nMatchProbe;
        int nJets;
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
template< typename T, typename P > 
TagAndProbeTreeWriter<T,P>::TagAndProbeTreeWriter(const edm::ParameterSet& iConfig)

{
    //now do what ever initialization is needed
    //Debug flag
    mcInfo          = iConfig.getUntrackedParameter<bool>   ("mcInfoAvailable",false);
    //Input collections
    mcSrc           = iConfig.getParameter<edm::InputTag> ("mcSource");
    tagSrc          = iConfig.getParameter<edm::InputTag> ("tagSource");
    probeSrc        = iConfig.getParameter<edm::InputTag> ("probeSource");
    jetSrc          = iConfig.getParameter<edm::InputTag> ("jetSource");
    
    cut_Dr          = iConfig.getUntrackedParameter<double> ("cut_TnPDr",0.1);
    cut_lowInvM     = iConfig.getUntrackedParameter<double> ("cut_TnPlowInvM",40.);
    cut_highInvM    = iConfig.getUntrackedParameter<double> ("cut_TnPhighInvM",120.);
    
    // Create the root file
    edm::Service<TFileService> theFile;

    TFileDirectory Tree = theFile->mkdir( "Trees" );
    treeTnP = Tree.make<TTree>("TnP","TnP"); 
    treeTnP->Branch("inv",&invM,"invM/F");
    treeTnP->Branch("pt",&ptProbe,"ptProbe/F");
    treeTnP->Branch("eta",&etaProbe,"etaProbe/F");
    treeTnP->Branch("nMatch",&nMatchProbe,"nMatchProbe/I");
    treeTnP->Branch("nJets",&nJets,"nJets/I");

    if (mcInfo){ 
        treeGen = Tree.make<TTree>("Gen tree", "Gen tree"); 
        treeGen->Branch("pt",&ptGen,"ptGen/F");
        treeGen->Branch("eta",&etaGen,"etaGen/F");
        treeGen->Branch("pdgId",&pdgIdLepton,"pdgIdLepton/I");
        treeGen->Branch("motherPdgId",&motherPdgIdLepton,"motherPdgIdLepton/I");
        treeGen->Branch("nJets",&nJets,"nJets/I");
        treeMatch = Tree.make<TTree>("Match tree", "Match tree tree"); 
        treeMatch->Branch("pt",&ptMatch,"ptMatch/F");
        treeMatch->Branch("eta",&etaMatch,"etaMatch/F");
        treeMatch->Branch("pdgId",&pdgIdLepton,"pdgIdLepton/I");
        treeMatch->Branch("nJets",&nJets,"nJets/I");
    }

}


template< typename T, typename P > 
TagAndProbeTreeWriter<T,P>::~TagAndProbeTreeWriter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions

template< typename T, typename P > 
void TagAndProbeTreeWriter<T,P>::TnP(const edm::Handle< std::vector<T> >& tags, const edm::Handle< P >& probes){
    for (typename std::vector<T>::const_iterator tag_i = tags->begin(); tag_i != tags->end(); ++tag_i){
        for (typename P::const_iterator pb_j = probes->begin(); pb_j != probes->end(); ++pb_j){
            nMatchProbe = 0;
            invM = 0.;
            ptProbe = pb_j->pt();
            etaProbe = pb_j->eta();
            double deltaRTnP = 9999999.;
            double chargesign = 0.;
            for (typename std::vector<T>::const_iterator tag_j = tags->begin(); tag_j != tags->end(); ++tag_j){
                deltaRTnP = reco::deltaR(tag_j->eta(),tag_j->phi(),pb_j->eta(),pb_j->phi());
                chargesign = tag_j->charge()*pb_j->charge(); 
                if (deltaRTnP < cut_Dr && chargesign>0){
                    ++nMatchProbe;
                }
                reco::Particle::LorentzVector pb = reco::Particle::LorentzVector(pb_j->px(),pb_j->py(),pb_j->pz(),pb_j->p());
                invM = (tag_i->p4()+pb).M();
                if (invM > cut_lowInvM && invM < cut_highInvM){treeTnP->Fill();}
            }
        }
    }
}

template< typename T, typename P > 
void TagAndProbeTreeWriter<T,P>::mcAnalysis(const edm::Handle< std::vector<T> >& tags, const edm::Handle< std::vector<reco::GenParticle> >& genParticles){
    for (std::vector<reco::GenParticle>::const_iterator p_i = genParticles->begin(); p_i != genParticles->end(); ++p_i){
        ptGen = p_i->pt();
        etaGen = p_i->eta();
        pdgIdLepton = p_i->pdgId();
        if (p_i->mother()){
            motherPdgIdLepton = p_i->mother()->pdgId();
        }
        else motherPdgIdLepton = 0;
        treeGen->Fill();
 	}
    for (typename std::vector<T>::const_iterator tag_i = tags->begin(); tag_i != tags->end(); ++tag_i){
        if(tag_i->genLepton()){
            pdgIdLepton = tag_i->genLepton()->pdgId();
            ptMatch = tag_i->genLepton()->pt();
            etaMatch = tag_i->genLepton()->eta();
            treeMatch->Fill();
        }
    }
}

// ------------ method called to for each event  ------------
template< typename T, typename P > 
void TagAndProbeTreeWriter<T,P>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    //Tags
    edm::Handle< std::vector<T> > tags;
    iEvent.getByLabel(tagSrc, tags);
  
    //Probes
    edm::Handle< P > probes;
    iEvent.getByLabel(probeSrc, probes);
   
    //Jets
    edm::Handle< std::vector<pat::Jet> > jets;
    iEvent.getByLabel(jetSrc, jets);

    //count the number of jets
    nJets = 0;
    for (std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i){       
        ++nJets;
    }

    //run the TnP
    TnP(tags,probes);

    if (mcInfo){
        //MC gen Particle
        edm::Handle< std::vector<reco::GenParticle> > genParticles;
        iEvent.getByLabel(mcSrc, genParticles);
        mcAnalysis(tags,genParticles);
    }
}


// ------------ method called once each job just before starting event loop  ------------
template< typename T, typename P > 
void TagAndProbeTreeWriter<T,P>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template< typename T, typename P > 
void TagAndProbeTreeWriter<T,P>::endJob() {
}

//define this as a plug-in
typedef TagAndProbeTreeWriter< pat::Muon, reco::TrackCollection > MuonTnPTreeWriter;
typedef TagAndProbeTreeWriter< pat::Electron, reco::CandidateCollection > ElectronTnPTreeWriter;
DEFINE_FWK_MODULE(MuonTnPTreeWriter);
DEFINE_FWK_MODULE(ElectronTnPTreeWriter);
