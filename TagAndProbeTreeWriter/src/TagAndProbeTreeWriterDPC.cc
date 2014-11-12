// -*- C++ -*-
//
// Package:    TagAndProbeTreeWriter
// Class:      TagAndProbeTreeWriterDPC
// 
// Description: TagAndProbeTreeWriter for Dilepton Pair Candidates
//
// Original Author:  sprenger,


// system include files
#include <memory>
#include <math.h>
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
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SuSyAachen/TagAndProbeTreeWriter/interface/DecayLengthFunctor.h"
#include "SuSyAachen/TagAndProbeTreeWriter/interface/IsolationFunctor.h"

//
// class declaration
//

template< typename T, typename P > 
class TagAndProbeTreeWriterDPC : public edm::EDAnalyzer {
public:
  explicit TagAndProbeTreeWriterDPC(const edm::ParameterSet&);
  ~TagAndProbeTreeWriterDPC();


private:
  typedef reco::CompositeCandidate cand;
  typedef edm::View<cand> collection;

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void TnP(const collection& pairs,const edm::Handle< std::vector<T> >&,
		   const edm::Handle< std::vector<pat::Jet> >&, const edm::Handle<reco::VertexCollection>&,
		   const edm::EventSetup& iSetup);
  virtual void mcAnalysis(const edm::Handle< std::vector<T> >&,const edm::Handle< std::vector<reco::GenParticle> >&);
  virtual void fillExtraVars(const reco::Candidate &probe);
  virtual void fillExtraVars(const reco::Track &probe);
  virtual void fillExtraVars(const pat::Electron &probe);
  virtual void fillExtraVars(const pat::Muon &probe);

  virtual void fillDecayLength(const reco::Candidate& tag, const reco::Candidate& probe,
			       const reco::Vertex& primaryVtx, const edm::EventSetup& iSetup);
  virtual void fillDecayLength(const reco::Candidate& tag, const reco::Track& probe,
			       const reco::Vertex& primaryVtx, const edm::EventSetup& iSetup);
  virtual void fillDecayLength(const pat::Electron& tag, const pat::Electron& probe,
			       const reco::Vertex& primaryVtx, const edm::EventSetup& iSetup);

  // ----------member data ---------------------------
  bool mcInfo;
  edm::InputTag mcSrc;
  //edm::InputTag tagSrc;
  //edm::InputTag probeSrc;
  edm::InputTag tnpPairsSrc;
  edm::InputTag passProbeSrc;
  edm::InputTag jetSrc;
  edm::InputTag vertexSrc;

  double cut_Dr;
  double cut_lowInvM;
  double cut_highInvM;
        
  TFile *theFile;

  DecayLengthFunctor fctDecayLength_;
  IsolationFunctor fctIsolation_;

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
  float ptTP;
  float etaProbe;
  float dPt;
  float dRJet;
  int chargeTagProbe;
  float pfIso;
  float pfIsoAbs;
  float pfIsoAbsChargedHadrons;
  float pfIsoAbsNeutralHadrons;
  float pfIsoAbsPhotons;
  float eOverP;
  float fBrem;
  float logSigIetaIeta;
  float deltaEtaIn;
  float mva;
  int nMatchProbe;
  int nJets;
  int nVertices;
  int nLostHits;
  int chargeMethodsProbe;
  int chargeDeviatingProbe;
  float decayLength;
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
TagAndProbeTreeWriterDPC<T,P>::TagAndProbeTreeWriterDPC(const edm::ParameterSet& iConfig):
  fctIsolation_  (iConfig.getParameter<edm::ParameterSet>("isolationDefinitions"))
{
    //now do what ever initialization is needed
    //Debug flag
    mcInfo          = iConfig.getUntrackedParameter<bool>   ("mcInfoAvailable",false);
    //Input collections
    mcSrc           = iConfig.getParameter<edm::InputTag> ("mcSource");
    //tagSrc          = iConfig.getParameter<edm::InputTag> ("tagSource");
    //probeSrc        = iConfig.getParameter<edm::InputTag> ("probeSource");
    tnpPairsSrc     = iConfig.getParameter<edm::InputTag> ("tnpPairsSource");
    passProbeSrc    = iConfig.getParameter<edm::InputTag> ("passProbeSource");
    jetSrc          = iConfig.getParameter<edm::InputTag> ("jetSource");
    vertexSrc       = iConfig.getParameter<edm::InputTag> ("vertexSource");
    
    cut_Dr          = iConfig.getUntrackedParameter<double> ("cut_TnPDr",0.2);
    cut_lowInvM     = iConfig.getUntrackedParameter<double> ("cut_TnPlowInvM",40.);
    cut_highInvM    = iConfig.getUntrackedParameter<double> ("cut_TnPhighInvM",120.);
    
    // Create the root file
    edm::Service<TFileService> theFile;

    TFileDirectory Tree = theFile->mkdir( "Trees" );
    treeTnP = Tree.make<TTree>("TnP","TnP"); 
    treeTnP->Branch("inv",&invM,"invM/F");
    treeTnP->Branch("pt",&ptProbe,"ptProbe/F");
    treeTnP->Branch("ptTP",&ptTP,"ptTP/F");
    treeTnP->Branch("eta",&etaProbe,"etaProbe/F");
    treeTnP->Branch("dPt",&dPt,"dPt/F");
    treeTnP->Branch("dRJet",&dRJet,"dRJet/F");
    treeTnP->Branch("chargeTP",&chargeTagProbe,"chargeTagProbe/I");
    treeTnP->Branch("nMatch",&nMatchProbe,"nMatchProbe/I");
    treeTnP->Branch("nJets",&nJets,"nJets/I");
    treeTnP->Branch("nVertices",&nVertices,"nVertices/I");
    treeTnP->Branch("nLostHits",&nLostHits,"nLostHits/I");
    treeTnP->Branch("pfIso",&pfIso,"pfIso/F");
    treeTnP->Branch("pfIsoAbs",&pfIsoAbs,"pfIsoAbs/F");
    treeTnP->Branch("pfIsoAbsChargedHadrons",&pfIsoAbsChargedHadrons,"pfIsoAbsChargedHadrons/F");
    treeTnP->Branch("pfIsoAbsNeutralHadrons",&pfIsoAbsNeutralHadrons,"pfIsoAbsNeutralHadrons/F");
    treeTnP->Branch("pfIsoAbsPhtotons",&pfIsoAbsPhotons,"pfIsoAbsPhotons/F");
    treeTnP->Branch("eOverP",&eOverP,"eOverP/F");
    treeTnP->Branch("fBrem",&fBrem,"fBrem/F");
    treeTnP->Branch("logSigIetaIeta",&logSigIetaIeta,"logSigIetaIeta/F");
    treeTnP->Branch("mva",&mva,"mva/F");
    treeTnP->Branch("chargeMethods",&chargeMethodsProbe,"chargeMethods/I");
    treeTnP->Branch("chargeDeviatingMethod",&chargeDeviatingProbe,"chargeDeviatingMethod/I");
    treeTnP->Branch("decayLength",&decayLength,"decayLength/F");

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
TagAndProbeTreeWriterDPC<T,P>::~TagAndProbeTreeWriterDPC()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
template< typename T, typename P > 
void TagAndProbeTreeWriterDPC<T,P>::fillExtraVars(const reco::Track& probe){
  mva = -1.;
  eOverP = -1;
  fBrem = -1;
  deltaEtaIn =-1;
  logSigIetaIeta = -1;
  nLostHits = -1;

  // charge methods
  chargeMethodsProbe = -1;
  chargeDeviatingProbe = -1;
}

template< typename T, typename P > 
void TagAndProbeTreeWriterDPC<T,P>::fillExtraVars(const reco::Candidate& probe){
  mva = -1.;
  eOverP = -1;
  fBrem = -1;
  deltaEtaIn =-1;
  logSigIetaIeta = -1;
  nLostHits = -1;

  // charge methods
  chargeMethodsProbe = -1;
  chargeDeviatingProbe = -1;
}

template< typename T, typename P > 
void TagAndProbeTreeWriterDPC<T,P>::fillExtraVars(const pat::Muon& probe){
  mva = -1.;
  eOverP = -1;
  fBrem = -1;
  deltaEtaIn =-1;
  logSigIetaIeta = -1;
  nLostHits = -1;

  // charge methods
  chargeMethodsProbe = -1;
  chargeDeviatingProbe = -1;
}

template< typename T, typename P > 
void TagAndProbeTreeWriterDPC<T,P>::fillExtraVars(const pat::Electron& probe){
  mva = -1.;
  eOverP = probe.eSuperClusterOverP();
  float pin  = probe.trackMomentumAtVtx().R();
  float pout = probe.trackMomentumOut().R();
  fBrem = (pin-pout)/pin;
  logSigIetaIeta = log(probe.sigmaIetaIeta());
  deltaEtaIn = probe.deltaEtaSuperClusterTrackAtVtx();
  nLostHits = probe.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
 
  // charge methods
  chargeMethodsProbe = 2;
  if (probe.isGsfCtfScPixChargeConsistent()){
    chargeMethodsProbe += 1;
    chargeDeviatingProbe = 0;
  }else{
    if (probe.isGsfScPixChargeConsistent())
      chargeDeviatingProbe = 1; //CTF
    else if (probe.isGsfCtfChargeConsistent())
      chargeDeviatingProbe = 2; //SC
    else
      chargeDeviatingProbe = 3; //GSF
  }
}

template< typename T, typename P > 
void TagAndProbeTreeWriterDPC<T,P>::fillDecayLength(const reco::Candidate& tag, const reco::Candidate& probe,
						    const reco::Vertex& primaryVtx, const edm::EventSetup& iSetup)
{
  decayLength = -10.0;
}

template< typename T, typename P > 
void TagAndProbeTreeWriterDPC<T,P>::fillDecayLength(const reco::Candidate& tag, const reco::Track& probe,
						    const reco::Vertex& primaryVtx, const edm::EventSetup& iSetup)
{
  decayLength = -10.0;
}

template< typename T, typename P > 
void TagAndProbeTreeWriterDPC<T,P>::fillDecayLength(const pat::Electron& tag, const pat::Electron& probe,
						    const reco::Vertex& primaryVtx, const edm::EventSetup& iSetup)
{
  decayLength = fctDecayLength_(tag, probe, primaryVtx, iSetup);
}

template< typename T, typename P > 
void TagAndProbeTreeWriterDPC<T,P>::TnP(const collection& pairs, const edm::Handle< std::vector<T> >& pass_probes, const edm::Handle< std::vector<pat::Jet> >& jets, const edm::Handle<reco::VertexCollection>& vertices, const edm::EventSetup& iSetup){
  for(collection::const_iterator it_pair = pairs.begin(); it_pair != pairs.end(); ++it_pair){
    edm::Ref<std::vector <T> > tag = (*it_pair).daughter("tag")->masterRef<edm::Ref< std::vector<T> > >();
    edm::Ref<std::vector <P> > probe = (*it_pair).daughter("probe")->masterRef<edm::Ref< std::vector<P> > >();

    nMatchProbe = 0;
    invM = 0.;
    ptProbe = probe->pt();
    etaProbe = probe->eta();
    chargeTagProbe = tag->charge()*probe->charge();
    dRJet = 9999999.;
    ptTP = (tag->p4() + probe->p4()).Pt();

    // isolation
    pfIsoAbs = fctIsolation_(*probe);
    pfIso = fctIsolation_(*probe) / probe->pt();
    pfIsoAbsChargedHadrons = -1.;//fctIsolationChargedHadrons_(*probe);
    pfIsoAbsNeutralHadrons = -1.;//fctIsolationNeutralHadrons_(*probe);
    pfIsoAbsPhotons = -1.;//fctIsolationPhotons_(*probe);

    // jet information
    for (typename std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i){
      dRJet = std::min(dRJet,(float)reco::deltaR(jet_i->eta(),jet_i->phi(),probe->eta(),probe->phi()));
    }

    // pass probe matches
    double deltaRTnP = 9999999.;
    dPt = 9999999.;

    for (typename std::vector<T>::const_iterator tag_j = pass_probes->begin(); tag_j != pass_probes->end(); ++tag_j){
      deltaRTnP = reco::deltaR(tag_j->eta(),tag_j->phi(),probe->eta(),probe->phi());
      if (deltaRTnP < cut_Dr){
	if (std::abs(probe->pt() - tag_j->pt()) < std::abs(dPt))
	  dPt = probe->pt()-tag_j->pt();
	++nMatchProbe;
      }
    }

    //reco::Particle::LorentzVector pb = reco::Particle::LorentzVector(probe->px(),probe->py(),probe->pz(),probe->p());
    fillExtraVars(*probe);
    fillDecayLength(*tag, *probe, *(vertices->begin()), iSetup);
    //invM = (tag->p4()+pb).M();
    invM = (tag->p4() + probe->p4()).M();
    if (invM > cut_lowInvM && invM < cut_highInvM){
      //std::cout << "DL: " << decayLength << std::endl;
      treeTnP->Fill();
    }
  }
}

template< typename T, typename P > 
void TagAndProbeTreeWriterDPC<T,P>::mcAnalysis(const edm::Handle< std::vector<T> >& pass_probes, const edm::Handle< std::vector<reco::GenParticle> >& genParticles){
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
void TagAndProbeTreeWriterDPC<T,P>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    if (iEvent.isRealData()) mcInfo = false;

    // TnP Pairs
    edm::Handle< collection > tnpPairs;
    iEvent.getByLabel(tnpPairsSrc, tnpPairs);

    //Passing Probes
    edm::Handle< std::vector< T > > pass_probes;
    iEvent.getByLabel(passProbeSrc, pass_probes);
   
    //Jets
    edm::Handle< std::vector<pat::Jet> > jets;
    iEvent.getByLabel(jetSrc, jets);
    
    //Vertices
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByLabel(vertexSrc, vertices);

    //count the number of jets
    nJets = jets->size();    
    nVertices = vertices->size();

    fctIsolation_.init(iEvent);

    //run the TnP
    TnP(*tnpPairs, pass_probes, jets, vertices, iSetup);
   
    if (mcInfo){
        //MC gen Particle
        edm::Handle< std::vector<reco::GenParticle> > genParticles;
        iEvent.getByLabel(mcSrc, genParticles);
        mcAnalysis(pass_probes,genParticles);
    }
}


// ------------ method called once each job just before starting event loop  ------------
template< typename T, typename P > 
void TagAndProbeTreeWriterDPC<T,P>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template< typename T, typename P > 
void TagAndProbeTreeWriterDPC<T,P>::endJob() {
}

//define this as a plug-in
//typedef TagAndProbeTreeWriterDPC< pat::Muon, reco::Track > MuonTnPTreeWriterDPC;
typedef TagAndProbeTreeWriterDPC< pat::Electron, pat::Electron > ElectronIsoTnPTreeWriterDPC;
typedef TagAndProbeTreeWriterDPC< pat::Electron, reco::Candidate > ElectronTnPTreeWriterDPC;
//typedef TagAndProbeTreeWriterDPC< pat::Electron, reco::Track > ElectronTrackTnPTreeWriterDPC;
typedef TagAndProbeTreeWriterDPC< pat::Muon, pat::Muon > MuonIsoTnPTreeWriterDPC;
//DEFINE_FWK_MODULE(MuonTnPTreeWriterDPC);
DEFINE_FWK_MODULE(ElectronTnPTreeWriterDPC);
//DEFINE_FWK_MODULE(ElectronTrackTnPTreeWriterDPC);
DEFINE_FWK_MODULE(MuonIsoTnPTreeWriterDPC);
DEFINE_FWK_MODULE(ElectronIsoTnPTreeWriterDPC);
