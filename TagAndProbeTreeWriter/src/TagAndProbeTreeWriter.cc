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
// $Id: TagAndProbeTreeWriter.cc,v 1.8 2010/08/09 16:02:42 nmohr Exp $
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

        virtual void TnP(const edm::Handle< std::vector<T> >&,const edm::Handle< P >&,const edm::Handle< std::vector<T> >&);
        virtual void mcAnalysis(const edm::Handle< std::vector<T> >&,const edm::Handle< std::vector<reco::GenParticle> >&);

        // ----------member data ---------------------------
        // Switch for debug output
        bool mcInfo;
        edm::InputTag mcSrc;
        edm::InputTag tagSrc;
        edm::InputTag probeSrc;
        edm::InputTag passProbeSrc;
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
        int chargeTagProbe;
        float pfIso;
        float mva;
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
    passProbeSrc    = iConfig.getParameter<edm::InputTag> ("passProbeSource");
    jetSrc          = iConfig.getParameter<edm::InputTag> ("jetSource");
    
    cut_Dr          = iConfig.getUntrackedParameter<double> ("cut_TnPDr",0.2);
    cut_lowInvM     = iConfig.getUntrackedParameter<double> ("cut_TnPlowInvM",40.);
    cut_highInvM    = iConfig.getUntrackedParameter<double> ("cut_TnPhighInvM",120.);
    
    // Create the root file
    edm::Service<TFileService> theFile;

    TFileDirectory Tree = theFile->mkdir( "Trees" );
    treeTnP = Tree.make<TTree>("TnP","TnP"); 
    treeTnP->Branch("inv",&invM,"invM/F");
    treeTnP->Branch("pt",&ptProbe,"ptProbe/F");
    treeTnP->Branch("eta",&etaProbe,"etaProbe/F");
    treeTnP->Branch("chargeTP",&chargeTagProbe,"chargeTagProbe/I");
    treeTnP->Branch("nMatch",&nMatchProbe,"nMatchProbe/I");
    treeTnP->Branch("nJets",&nJets,"nJets/I");
    treeTnP->Branch("pfIso",&pfIso,"pfIso/F");
    treeTnP->Branch("mva",&mva,"mva/F");

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
template< class PB > 
void fillExtraVars(const PB *pb_j, float *iso, float *mva){
    *iso = -1.;
    *mva = -1.;
}
template< > 
void fillExtraVars(const pat::Muon *lepton, float *iso, float *mva){
    *iso = (lepton->chargedHadronIso()+lepton->photonIso()+lepton->neutralHadronIso())/lepton->pt();
    *mva = -1.;
}
template< > 
void fillExtraVars(const pat::Electron *lepton, float *iso, float *mva){
    *iso = (lepton->chargedHadronIso()+lepton->photonIso()+lepton->neutralHadronIso())/lepton->pt();
    *mva = lepton->mva();
}

template< typename T, typename P > 
void TagAndProbeTreeWriter<T,P>::TnP(const edm::Handle< std::vector<T> >& tags, const edm::Handle< P >& probes, const edm::Handle< std::vector<T> >& pass_probes){
    for (typename std::vector<T>::const_iterator tag_i = tags->begin(); tag_i != tags->end(); ++tag_i){
        for (typename P::const_iterator pb_j = probes->begin(); pb_j != probes->end(); ++pb_j){
            nMatchProbe = 0;
            invM = 0.;
            ptProbe = pb_j->pt();
            etaProbe = pb_j->eta();
            chargeTagProbe = tag_i->charge()*pb_j->charge(); 
            double deltaRTnP = 9999999.;
            for (typename std::vector<T>::const_iterator tag_j = pass_probes->begin(); tag_j != pass_probes->end(); ++tag_j){
                deltaRTnP = reco::deltaR(tag_j->eta(),tag_j->phi(),pb_j->eta(),pb_j->phi());
                if (deltaRTnP < cut_Dr){
                    ++nMatchProbe;
                }
            }
            reco::Particle::LorentzVector pb = reco::Particle::LorentzVector(pb_j->px(),pb_j->py(),pb_j->pz(),pb_j->p());
            fillExtraVars(&(*pb_j),&pfIso,&mva);

            invM = (tag_i->p4()+pb).M();
            if (invM > cut_lowInvM && invM < cut_highInvM){treeTnP->Fill();}
        }
    }
}

template< typename T, typename P > 
void TagAndProbeTreeWriter<T,P>::mcAnalysis(const edm::Handle< std::vector<T> >& pass_probes, const edm::Handle< std::vector<reco::GenParticle> >& genParticles){
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
    for (typename std::vector<T>::const_iterator tag_i = pass_probes->begin(); tag_i != pass_probes->end(); ++tag_i){
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
    if (iEvent.isRealData()) mcInfo = false;

    //Tags
    edm::Handle< std::vector<T> > tags;
    iEvent.getByLabel(tagSrc, tags);
  
    //Probes
    edm::Handle< P > probes;
    iEvent.getByLabel(probeSrc, probes);
    
    //Passing Probes
    edm::Handle< std::vector< T > > pass_probes;
    iEvent.getByLabel(passProbeSrc, pass_probes);
   
    //Jets
    edm::Handle< std::vector<pat::Jet> > jets;
    iEvent.getByLabel(jetSrc, jets);

    //count the number of jets
    nJets = 0;
    for (std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i){       
        ++nJets;
    }

    //run the TnP
    TnP(tags,probes,pass_probes);

    if (mcInfo){
        //MC gen Particle
        edm::Handle< std::vector<reco::GenParticle> > genParticles;
        iEvent.getByLabel(mcSrc, genParticles);
        mcAnalysis(pass_probes,genParticles);
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
typedef TagAndProbeTreeWriter< pat::Muon, pat::MuonCollection > MuonIsoTnPTreeWriter;
typedef TagAndProbeTreeWriter< pat::Electron, pat::ElectronCollection > ElectronIsoTnPTreeWriter;
DEFINE_FWK_MODULE(MuonTnPTreeWriter);
DEFINE_FWK_MODULE(ElectronTnPTreeWriter);
DEFINE_FWK_MODULE(MuonIsoTnPTreeWriter);
DEFINE_FWK_MODULE(ElectronIsoTnPTreeWriter);
