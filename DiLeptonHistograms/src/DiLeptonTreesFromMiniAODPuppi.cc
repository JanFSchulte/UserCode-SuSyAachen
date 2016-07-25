// -*- C++ -*-
//
// Package:    Histograms
// Class:      DiLeptonTreesFromMiniAODPuppi
// 
/**\class DiLeptonTreesFromMiniAODPuppi DiLeptonTreesFromMiniAODPuppi.cc brot/DiLeptonTreesFromMiniAODPuppi/src/DiLeptonTreesFromMiniAODPuppi.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  matthias edelhoff
//         Created:  Tue Oct 27 13:50:40 CET 2009
// $Id: DiLeptonTreesFromMiniAODPuppi.cc,v 1.31 2012/09/17 17:38:58 sprenger Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/Lepton.h>

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/MET.h>

#include <DataFormats/Provenance/interface/EventID.h>

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

// For JES
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "RecoJets/JetProducers/interface/PileupJetIdAlgo.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"


#include <SuSyAachen/DiLeptonHistograms/interface/MT2Functor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/WeightFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/PdgIdFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/VertexWeightFunctor.h>
#include <SuSyAachen/TagAndProbeTreeWriter/interface/IsolationFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/TriggerMatchFunctorMiniAOD.h>


//ROOT
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace std;

//
// class decleration
//

class DiLeptonTreesFromMiniAODPuppi : public edm::EDAnalyzer {
public:
  explicit DiLeptonTreesFromMiniAODPuppi(const edm::ParameterSet&);
  ~DiLeptonTreesFromMiniAODPuppi();

private:
  //  typedef reco::Candidate candidate;
  typedef pat::Lepton<reco::Candidate> candidate;
  typedef edm::View<candidate> collection;

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  void initFloatBranch( const std::string &name);
  void initIntBranch( const std::string &name);
  void initTLorentzVectorBranch( const std::string &name);
  template <class aT, class bT> void makeCombinations( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT >&b, const std::vector<pat::PackedCandidate>&pfCands,const std::vector<pat::Electron>&looseElectrons,const std::vector<pat::Muon>&looseMuons, const std::vector<pat::Jet>&jets, const std::vector<pat::Jet>&bJets35, const edm::Event &ev, const pat::MET &patMet, const pat::MET &patMetPuppi, const TLorentzVector &MHT, const TLorentzVector &MHT20, const TLorentzVector &MHTPuppi, const TLorentzVector &MHT20Puppi, const edm::Handle<reco::VertexCollection> &vertices, const float &rho, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties);
  template <class aT> void makeCombinations( const std::string &treeName, const std::vector<aT> &a , const std::vector<pat::PackedCandidate>&pfCands,const std::vector<pat::Electron>&looseElectrons,const std::vector<pat::Muon>&looseMuons,const std::vector<pat::Jet>&jets,const std::vector<pat::Jet>&bJets35,   const edm::Event &ev, const pat::MET &patMet, const pat::MET &patMetPuppi, const TLorentzVector &MHT, const TLorentzVector &MHT20,const TLorentzVector &MHTPuppi, const TLorentzVector &MHT20Puppi, const edm::Handle<reco::VertexCollection> &vertices, const float &rho, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties);
  template<class aT, class bT> void fillTree( const std::string &treeName, const aT &a, const bT &b, const std::vector<pat::PackedCandidate>&pfCands,const std::vector<pat::Electron>&looseElectrons,const std::vector<pat::Muon>&looseMuons,const std::vector<pat::Jet>&jets,const std::vector<pat::Jet>&bJets35,  const pat::MET &patMet, const pat::MET &patMetPuppi, const TLorentzVector &MHT, const TLorentzVector &MHT20, const TLorentzVector &MHTPuppi, const TLorentzVector &MHT20Puppi, const edm::Handle<reco::VertexCollection> &vertices, const float &rho);
  //~ std::pair<double, double> calcPZeta(const TLorentzVector& p1,const TLorentzVector& p2, const TLorentzVector& met);
  void fillPdfUncert(const edm::Handle< std::vector<double> >& weightHandle, const std::string& pdfIdentifier, const std::string& treeName);

  const TLorentzVector getMomentum(const  pat::Electron &e);
  const TLorentzVector getMomentum(const  pat::Muon &mu);
  void fillLeptonIDs(const std::string &treeName, const  pat::Electron &ele1, const  pat::Electron &ele2, const edm::Handle<reco::VertexCollection> &vertices);
  void fillLeptonIDs(const std::string &treeName, const  pat::Muon &mu1, const  pat::Muon &mu2, const edm::Handle<reco::VertexCollection> &vertices);
  void fillLeptonIDs(const std::string &treeName, const  pat::Electron &ele1, const  pat::Muon &mu2, const edm::Handle<reco::VertexCollection> &vertices);
  float topPtWeightBen(double topPt);
  float topPtWeightTOP(double topPt);
  float getIso(const  pat::Electron &e, const std::string &method);
  float getIso(const  pat::Muon &mu, const std::string &method);
  float getDeltaB(const  pat::Electron &e);
  float getDeltaB(const  pat::Muon &mu);
  float transverseMass(const TLorentzVector& p, const TLorentzVector& met);
  float getAEffEle(double eta);
  float getAEffMu(double eta);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  std::string convertInputTag(const edm::InputTag tag);
  

  edm::EDGetTokenT< std::vector< pat::Electron > > 			electronToken_;
  edm::EDGetTokenT< std::vector< pat::Electron > > 			looseElectronToken_;
  edm::EDGetTokenT< std::vector< pat::Muon > > 				muonToken_;
  edm::EDGetTokenT< std::vector< pat::Muon > > 				looseMuonToken_;
  edm::EDGetTokenT< std::vector< pat::Jet > >				jetToken_;
  edm::EDGetTokenT< std::vector< reco::GenJet > >			genJetToken_;
  edm::EDGetTokenT< std::vector< pat::Jet > >				jetPuppiToken_;
  edm::EDGetTokenT< std::vector< pat::Jet > >		 		bJetToken_;
  edm::EDGetTokenT< std::vector< pat::Jet > >		 		bJet35Token_;
  edm::EDGetTokenT< std::vector< pat::Jet > >				bJetPuppiToken_;
  edm::EDGetTokenT< std::vector< pat::MET > > 				metToken_;
  edm::EDGetTokenT< std::vector< pat::MET > >				metPuppiToken_;
  edm::EDGetTokenT<reco::VertexCollection> 					vertexToken_;
  edm::EDGetTokenT< std::vector<pat::PackedCandidate>  > 	pfCandToken_;
  edm::EDGetTokenT< std::vector< reco::GenParticle > >	 	genParticleToken_;
  edm::EDGetTokenT<GenEventInfoProduct>	 					genEventInfoToken_;
  edm::EDGetTokenT<LHEEventProduct>	 						LHEEventToken_;
  edm::EDGetTokenT<double>	 								rhoToken_;
  

  std::map<double, double> electronCorrections_;
  //data
  std::map<std::string, TTree*> trees_;  
  std::map<std::string, std::map< std::string, float*> > floatBranches_; 
  std::map<std::string, std::map< std::string, unsigned int*> > intBranches_; 
  std::map<std::string, std::map< std::string, bool*> > boolBranches_; 
  std::map<std::string, std::map< std::string, TLorentzVector*> > tLorentzVectorBranches_;

  edm::Handle< std::vector< pat::Jet > > jets;

  VertexWeightFunctor fctVtxWeight_;
  VertexWeightFunctor fctVtxWeightUp_;
  VertexWeightFunctor fctVtxWeightDown_;
  IsolationFunctor fctIsolation_;
  TriggerMatchFunctorMiniAOD fctTrigger_;  
  PdgIdFunctor getPdgId_;
  MT2Functor fctMT2_;
 
  bool debug;
  bool writeID_;
  bool writeTrigger_;
  bool metUncert_;  
  bool triggerMatches_;
  std::vector<std::string> triggerNames_;
  bool newLumiBlock_;
  std::map<std::string, Bool_t > triggerDecision_;
  std::map<std::string, int > triggerIndex_;
};

// constructors and destructor
DiLeptonTreesFromMiniAODPuppi::DiLeptonTreesFromMiniAODPuppi(const edm::ParameterSet& iConfig):
  electronToken_		(consumes< std::vector< pat::Electron > > 		(iConfig.getParameter<edm::InputTag>("electrons"))),
  looseElectronToken_	(consumes< std::vector< pat::Electron > > 		(iConfig.getParameter<edm::InputTag>("looseElectrons"))),
  muonToken_			(consumes< std::vector< pat::Muon > >			(iConfig.getParameter<edm::InputTag>("muons"))),
  looseMuonToken_		(consumes< std::vector< pat::Muon > >			(iConfig.getParameter<edm::InputTag>("looseMuons"))),
  jetToken_				(consumes< std::vector< pat::Jet > >			(iConfig.getParameter<edm::InputTag>("jets"))),
  genJetToken_			(consumes< std::vector< reco::GenJet  > >		(iConfig.getParameter<edm::InputTag>("genJets"))),
  jetPuppiToken_		(consumes< std::vector< pat::Jet > >			(iConfig.getParameter<edm::InputTag>("jetsPuppi"))),
  bJetToken_			(consumes< std::vector< pat::Jet > >			(iConfig.getParameter<edm::InputTag>("bJets"))),
  bJet35Token_			(consumes< std::vector< pat::Jet > >			(iConfig.getParameter<edm::InputTag>("bJets35"))),
  bJetPuppiToken_		(consumes< std::vector< pat::Jet > >			(iConfig.getParameter<edm::InputTag>("bJetsPuppi"))),
  metToken_				(consumes< std::vector< pat::MET > >			(iConfig.getParameter<edm::InputTag>("met"))),
  metPuppiToken_		(consumes< std::vector< pat::MET > >			(iConfig.getParameter<edm::InputTag>("metPuppi"))),
  vertexToken_			(consumes<reco::VertexCollection>				(iConfig.getParameter<edm::InputTag>("vertices"))),
  pfCandToken_			(consumes< std::vector<pat::PackedCandidate>  >	(iConfig.getParameter<edm::InputTag>("pfCands"))),
  genParticleToken_		(consumes< std::vector< reco::GenParticle > >	(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genEventInfoToken_	(consumes<GenEventInfoProduct>					(iConfig.getParameter<edm::InputTag>("pdfInfo"))),
  LHEEventToken_		(consumes<LHEEventProduct>						(iConfig.getParameter<edm::InputTag>("LHEInfo"))),
  rhoToken_				(consumes<double>								(iConfig.getParameter<edm::InputTag>("rho"))),
  
  fctVtxWeight_   	(iConfig.getParameter<edm::ParameterSet>("vertexWeights") ,consumesCollector()),
  fctVtxWeightUp_   (iConfig.getParameter<edm::ParameterSet>("vertexWeightsUp") ,consumesCollector()),
  fctVtxWeightDown_ (iConfig.getParameter<edm::ParameterSet>("vertexWeightsDown"),consumesCollector() ),
  fctIsolation_  	(iConfig.getParameter<edm::ParameterSet>("isolationDefinitions"),consumesCollector()),
  fctTrigger_  		(iConfig.getParameter<edm::ParameterSet>("triggerDefinitions"),consumesCollector()),  
  getPdgId_			(iConfig.getParameter< edm::ParameterSet>("pdgIdDefinition"),consumesCollector() ),
  triggerNames_		(iConfig.getParameter< std::vector <std::string> >("triggerNames")),
  newLumiBlock_(true)
{
  debug = false; 
  //~ writeID_ = iConfig.existsAs<edm::InputTag>("baseTrees");
  writeID_ = iConfig.getUntrackedParameter<bool>("writeID");
  writeTrigger_ = iConfig.getUntrackedParameter<bool>("writeTrigger");
  triggerMatches_ = iConfig.existsAs<edm::InputTag>("triggerSummaryTag");  
  metUncert_ = iConfig.existsAs<edm::InputTag>("doMETUncert");  
  
  consumes<std::vector< PileupSummaryInfo > >(edm::InputTag("slimmedAddPileupInfo"));
  
  consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
  consumes<edm::TriggerResults>(edm::InputTag("TriggerResults",""));
  //~ consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"));
  //~ consumes<std::vector<pat::TriggerObjectStandAlone>>(edm::InputTag("selectedPatTrigger"));
	  
  if(iConfig.existsAs<edm::VParameterSet>("electronCorrections")){
    edm::VParameterSet bins = iConfig.getParameter<edm::VParameterSet>("electronCorrections");
    for(edm::VParameterSet::const_iterator it = bins.begin(); it != bins.end(); ++it){
      float absEta = (*it).getParameter<double>("absEta");
      electronCorrections_[absEta] = (*it).getParameter<double>("correction");
    }
  }

  // init trees
  edm::Service<TFileService> file;
  trees_["EE"] = file->make<TTree>("EEDileptonTree", "EE DileponTree");
  trees_["EMu"] = file->make<TTree>("EMuDileptonTree", "EMu DileponTree");
  trees_["MuMu"] = file->make<TTree>("MuMuDileptonTree", "MuMu DileponTree");

  initFloatBranch( "genWeight" );  
  initFloatBranch( "genWeightAbsValue" );    
  initFloatBranch( "weight" );
  initFloatBranch( "weightUp" );
  initFloatBranch( "weightDown" );
  initFloatBranch( "genPtTop1" );
  initFloatBranch( "genPtTop2" );
  initFloatBranch( "chargeProduct" );
  initFloatBranch( "charge1" );
  initFloatBranch( "charge2" );
  initTLorentzVectorBranch( "p4" );
  initTLorentzVectorBranch( "p4Gen" );
  initTLorentzVectorBranch( "lepton1" );
  initTLorentzVectorBranch( "lepton2" );
  initTLorentzVectorBranch( "genLepton1" );
  initTLorentzVectorBranch( "genLepton2" );
  initTLorentzVectorBranch( "jet1" );
  initTLorentzVectorBranch( "jet2" );
  initTLorentzVectorBranch( "jet3" );
  initTLorentzVectorBranch( "jet4" );
  initTLorentzVectorBranch( "jet1Puppi" );
  initTLorentzVectorBranch( "jet2Puppi" );
  initTLorentzVectorBranch( "jet3Puppi" );
  initTLorentzVectorBranch( "jet4Puppi" );
  initTLorentzVectorBranch( "genJet1" );
  initTLorentzVectorBranch( "genJet2" );
  initTLorentzVectorBranch( "genJet3" );
  initTLorentzVectorBranch( "genJet4" );
  initTLorentzVectorBranch( "bJet1" );
  initTLorentzVectorBranch( "bJet2" );
  initTLorentzVectorBranch( "bJet1Puppi" );
  initTLorentzVectorBranch( "bJet2Puppi" );
  initTLorentzVectorBranch( "vMet" );  
  initTLorentzVectorBranch( "vMetPuppi" );  
  initTLorentzVectorBranch( "vMetUncorrected" ); 
  initTLorentzVectorBranch( "vGenMet" ); 
  initTLorentzVectorBranch( "vMHT" );   
  initTLorentzVectorBranch( "vMHT20" );
  initTLorentzVectorBranch( "vMHTPuppi" );   
  initTLorentzVectorBranch( "vMHT20Puppi" );
  initFloatBranch( "rho" );
  initFloatBranch( "pt1" );
  initFloatBranch( "pt2" );
  initFloatBranch( "pt3" );
  initFloatBranch( "HoverE3" );
  initFloatBranch( "deltaEtaSuperCluster3" );
  initFloatBranch( "deltaPhiSuperCluster3" );
  initFloatBranch( "eInvMinusPInv3" );
  initFloatBranch( "sigmaIEtaIEta3" );
  initFloatBranch( "eta1" );
  initFloatBranch( "eta2" );
  initFloatBranch( "miniIsoEffArea1" );
  initFloatBranch( "miniIsoEffArea2" );
  initFloatBranch( "miniIsoPuppi1" );
  initFloatBranch( "miniIsoPuppi2" );
    //~ 
  //~ initFloatBranch( "chargedIso1");
  //~ initFloatBranch( "neutralIso1");
  //~ initFloatBranch( "photonIso1");
  //~ initFloatBranch( "puIso1");
  //~ initFloatBranch( "effectiveArea1");
  //~ 
  //~ initFloatBranch( "chargedIso2");
  //~ initFloatBranch( "neutralIso2");
  //~ initFloatBranch( "photonIso2");
  //~ initFloatBranch( "puIso2");
  //~ initFloatBranch( "effectiveArea2");
  
  initFloatBranch( "mt1" );
  initFloatBranch( "mt2" );
  //~ initFloatBranch( "fakeWeight1" );
  //~ initFloatBranch( "fakeWeight2" );
  initFloatBranch( "deltaPhi" );
  initFloatBranch( "deltaR" );
  initFloatBranch( "angle3D" );
  initFloatBranch( "min_mlb" );
  initFloatBranch( "max_mlb" );
  initFloatBranch( "sumMlb" );
  initFloatBranch( "NLL" );
  initFloatBranch( "MT2" );
  initFloatBranch( "jzbResponse" );
  initFloatBranch( "jzbResponseUncorr" );  
  initFloatBranch( "jzbResponsePuppi" );  
  initFloatBranch( "jzb" );
  initFloatBranch( "jzbUncorr" );  
  initFloatBranch( "jzbPuppi" );  
  initFloatBranch( "ht" );
  initFloatBranch( "htPuppi" );
  initFloatBranch( "htJESUp" );
  initFloatBranch( "htJESDown" );
  initFloatBranch( "ht40" );  
  initFloatBranch( "leptHT" );  
  initFloatBranch( "genLeptHT" );  
  initFloatBranch( "leptHTPuppi" );  
  initFloatBranch( "ET" );  
  initFloatBranch( "genET" );  
  initFloatBranch( "ETPuppi" );
  initFloatBranch( "packedET" );
  initFloatBranch( "genParticleET" );
  initFloatBranch( "mht" );
  initFloatBranch( "mhtPuppi" );
  initFloatBranch( "met" );
  initFloatBranch( "metPuppi" );
  initFloatBranch( "genMet" );
  initFloatBranch( "uncorrectedMet" ); 
  initFloatBranch( "metJESUp" );
  initFloatBranch( "metJESDown" );
  initFloatBranch( "metHadRecoil" );
  initFloatBranch( "metHadRecoil10" );
  initFloatBranch( "metHadRecoil20" );
  initFloatBranch( "metHadRecoil30" );
  initFloatBranch( "metHadRecoilPuppi" );
  initFloatBranch( "metHadRecoil10Puppi" );
  initFloatBranch( "metHadRecoil20Puppi" );
  initFloatBranch( "metHadRecoil30Puppi" );
  initFloatBranch( "genMetHadRecoil" );
  initFloatBranch( "genMetHadRecoil10" );
  initFloatBranch( "genMetHadRecoil20" );
  initFloatBranch( "genMetHadRecoil30" );
  //~ initFloatBranch( "pZeta" );
  //~ initFloatBranch( "pZetaVis" );
  initIntBranch( "nJets" );
  initIntBranch( "nJetsPuppi" );
  initIntBranch( "nGenJets" );
  initIntBranch( "nMatchedJets" );
  initIntBranch( "nMatchedJetsPuppi" );
  //~ initIntBranch( "nJets30" );  
  //initIntBranch( "nJetsNoPULoose" );
  //initIntBranch( "nJetsNoPUMedium" );
  //initIntBranch( "nJetsNoPUTight" ); 
  initIntBranch( "nJets10" );
  initIntBranch( "nJets20" );
  initIntBranch( "nJets30" );
  initIntBranch( "nJetsPuppi10" );
  initIntBranch( "nJetsPuppi20" );
  initIntBranch( "nJetsPuppi30" );
  initIntBranch( "nGenJets10" );     
  initIntBranch( "nGenJets20" );     
  initIntBranch( "nGenJets30" );     
  initIntBranch( "nJetsOld" );
  initIntBranch( "nBJets" );
  initIntBranch( "nBJets35" );
  initIntBranch( "nBJetsPuppi" );
  initIntBranch( "nShiftedJetsJESUp" );
  initIntBranch( "nShiftedJetsJESDown" );
  initIntBranch( "nVertices" );
  initIntBranch( "nGenVertices" );
  initIntBranch( "nLightLeptons" );
  initIntBranch( "nLooseLeptons" );
  initFloatBranch( "jet1pt" );
  initFloatBranch( "jet2pt" );
  initFloatBranch( "jet3pt" );
  initFloatBranch( "jet4pt" );
  initFloatBranch( "jet1ptPuppi" );
  initFloatBranch( "jet2ptPuppi" );
  initFloatBranch( "jet3ptPuppi" );
  initFloatBranch( "jet4ptPuppi" );
  initFloatBranch( "genJet1pt" );
  initFloatBranch( "genJet2pt" );
  initFloatBranch( "genJet3pt" );
  initFloatBranch( "genJet4pt" );
  initFloatBranch( "bjet1pt" );
  initFloatBranch( "bjet2pt" );
  initFloatBranch( "bjet1ptPuppi" );
  initFloatBranch( "bjet2ptPuppi" );
  initFloatBranch( "genHT" );
  initFloatBranch( "genParticleHT" );
  //~ initFloatBranch( "sqrts" ); 
  initIntBranch( "runNr" );
  initIntBranch( "lumiSec" );
  initIntBranch( "eventNr" );
  initIntBranch( "pdgId1" );
  initIntBranch( "pdgId2" );
  initIntBranch( "matched" );
  initIntBranch( "motherPdgId1" );
  initIntBranch( "motherPdgId2" );
  initIntBranch( "grandMotherPdgId1" );
  initIntBranch( "grandMotherPdgId2" );  
  initIntBranch( "isPrompt1" );
  initIntBranch( "isPrompt2" );
  initIntBranch( "isFromTau1" );
  initIntBranch( "isFromTau2" );
  initIntBranch( "isPromptHardProcess1" );
  initIntBranch( "isPromptHardProcess2" );
  initIntBranch( "isFromTauHardProcess1" );
  initIntBranch( "isFromTauHardProcess2" );
  
  if (writeTrigger_){
	  for (const auto& n : triggerNames_){
		  triggerIndex_[n] = -10; //not set and not found
		  triggerDecision_[n] = false;
		  //~ initIntBranch(n.c_str())
		  //~ eventTree_->Branch( n.c_str(), &triggerDecision_[n], (n+"/O").c_str() );
		  
		  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();it != trees_.end(); ++it){
			  if(debug) std::cout << (*it).first <<" - "<< n.c_str() << std::endl;
			  boolBranches_[(*it).first][n.c_str()] = new bool;
			  (*it).second->Branch(n.c_str(), &triggerDecision_[n] ,(n+"/O").c_str() );
  }
	  }
  }
      
      
  if (writeID_){
 
 
	  initFloatBranch( "deltaEtaSuperClusterTrackAtVtx1");
	  initFloatBranch( "deltaPhiSuperClusterTrackAtVtx1");
	  initFloatBranch( "sigmaIetaIeta1");
	  initFloatBranch( "hadronicOverEm1");
	  initFloatBranch( "eOverP1");
	  initFloatBranch( "missingHits1");
	  initFloatBranch( "d01");
	  initFloatBranch( "dZ1");
	  initFloatBranch( "globalMuon1");
	  initFloatBranch( "trackerMuon1");
	  initFloatBranch( "pfMuon1");
	  initFloatBranch( "trackChi21");
	  initFloatBranch( "numberOfValidMuonHits1");
	  initFloatBranch( "numberOfMatchedStations1");
	  initFloatBranch( "numberOfValidPixelHits1");
	  initFloatBranch( "trackerLayersWithMeasurement1");
	  initFloatBranch( "chargedIso1");
	  initFloatBranch( "neutralIso1");
	  initFloatBranch( "photonIso1");
	  initFloatBranch( "puIso1");
	  initFloatBranch( "effectiveArea1");
	  initFloatBranch( "passConversion1");

	  initFloatBranch( "deltaEtaSuperClusterTrackAtVtx2");
	  initFloatBranch( "deltaPhiSuperClusterTrackAtVtx2");
	  initFloatBranch( "sigmaIetaIeta2");
	  initFloatBranch( "hadronicOverEm2");
	  initFloatBranch( "eOverP2");
	  initFloatBranch( "missingHits2");
	  initFloatBranch( "d02");
	  initFloatBranch( "dZ2");
	  initFloatBranch( "globalMuon2");
	  initFloatBranch( "trackerMuon2");
	  initFloatBranch( "pfMuon2");
	  initFloatBranch( "trackChi22");
	  initFloatBranch( "numberOfValidMuonHits2");
	  initFloatBranch( "numberOfMatchedStations2");
	  initFloatBranch( "numberOfValidPixelHits2");
	  initFloatBranch( "trackerLayersWithMeasurement2");
	  initFloatBranch( "chargedIso2");
	  initFloatBranch( "neutralIso2");
	  initFloatBranch( "photonIso2");
	  initFloatBranch( "puIso2");
	  initFloatBranch( "effectiveArea2");
	  initFloatBranch( "passConversion2"); 
  
  }
  
  if (triggerMatches_){
  		
	  initIntBranch( "matchesSingleElectron1" );
	  initIntBranch( "matchesSingleElectron2" );
	  initIntBranch( "matchesSingleMuon1" );
	  initIntBranch( "matchesSingleMuon2" );
	  initIntBranch( "matchesDoubleElectronTrailing1" );
	  initIntBranch( "matchesDoubleElectronTrailing2" );
	  initIntBranch( "matchesDoubleElectronTrailingNonIso1" );
	  initIntBranch( "matchesDoubleElectronTrailingNonIso2" );	  
	  initIntBranch( "matchesDoubleElectronLeading1" );
	  initIntBranch( "matchesDoubleElectronLeading2" );	 
	  initIntBranch( "matchesDoubleElectronLeadingNonIso1" );
	  initIntBranch( "matchesDoubleElectronLeadingNonIso2" );		  
	  initIntBranch( "matchesDoubleMuonLeading1" );
	  initIntBranch( "matchesDoubleMuonLeading2" );
	  initIntBranch( "matchesDoubleMuonLeadingNonIso1" );
	  initIntBranch( "matchesDoubleMuonLeadingNonIso2" );	  
	  initIntBranch( "matchesDoubleMuonLeadingBoth1" );
	  initIntBranch( "matchesDoubleMuonLeadingBoth2" );	  
	  initIntBranch( "matchesDoubleMuonLeadingTk1" );
	  initIntBranch( "matchesDoubleMuonLeadingTk2" );
	  initIntBranch( "matchesEMuLeading1" );
	  initIntBranch( "matchesEMuLeading2" );
	  initIntBranch( "matchesMuELeading1" );
	  initIntBranch( "matchesMuELeading2" );
	  initIntBranch( "matchesDoubleMuonTrailing1" );
	  initIntBranch( "matchesDoubleMuonTrailing2" );
	  initIntBranch( "matchesDoubleMuonTrailingNonIso1" );
	  initIntBranch( "matchesDoubleMuonTrailingNonIso2" );	  
	  initIntBranch( "matchesDoubleMuonTrailingBoth1" );
	  initIntBranch( "matchesDoubleMuonTrailingBoth2" );	  
	  initIntBranch( "matchesDoubleMuonTrailingTk1" );
	  initIntBranch( "matchesDoubleMuonTrailingTk2" );
	  initIntBranch( "matchesEMuTrailing1" );
	  initIntBranch( "matchesEMuTrailing2" );
	  initIntBranch( "matchesMuETrailing1" );
	  initIntBranch( "matchesMuETrailing2" );  
	  initIntBranch( "matchesMuEGMuonNonIso1" );
	  initIntBranch( "matchesMuEGMuonNonIso2" );
	  initIntBranch( "matchesMuEGElectronNonIso1" );
	  initIntBranch( "matchesMuEGElectronNonIso2" );  
  
  }
  
  
  if (metUncert_){
  		
	  initFloatBranch( "metJetEnUp");
	  initFloatBranch( "metJetEnDown");
	  initFloatBranch( "metJetResUp");
	  initFloatBranch( "metJetResDown");
	  initFloatBranch( "metMuonEnUp");
	  initFloatBranch( "metMuonEnDown");
	  initFloatBranch( "metElectronEnUp");
	  initFloatBranch( "metElectronEnDown");
	  initFloatBranch( "metTauEnUp");
	  initFloatBranch( "metTauEnDown");
	  initFloatBranch( "metUnclusteredEnUp");
	  initFloatBranch( "metUnclusteredEnDown");	  	  
  
  
  }  
  
  //~ for ( std::vector<edm::ParameterSet>::iterator susyVar_i = susyVars_.begin(); susyVar_i != susyVars_.end(); ++susyVar_i ) {
    //~ edm::InputTag var = susyVar_i->getParameter<edm::InputTag>( "var" );
    //~ std::string type = susyVar_i->getParameter<std::string>( "type" );
    //~ if(debug) std::cout << var << " of type " << type << std::endl;
    //~ if (type=="int") initIntBranch( convertInputTag(var) );
    //~ else if (type=="float") initFloatBranch( convertInputTag(var) );
    //~ else throw cms::Exception("Unrecognized type") << 
      //~ "Unknown type " << type << " for variable" << var << " found\n";
  //~ }
}

void 
DiLeptonTreesFromMiniAODPuppi::initTLorentzVectorBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    tLorentzVectorBranches_[(*it).first][name] = new TLorentzVector;
    (*it).second->Branch(name.c_str(), "TLorentzVector" ,&tLorentzVectorBranches_[(*it).first][name]);
  }
}

void 
DiLeptonTreesFromMiniAODPuppi::initFloatBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    floatBranches_[(*it).first][name] = new float;
    (*it).second->Branch(name.c_str(), floatBranches_[(*it).first][name], (name+"/F").c_str());
  }
}

void 
DiLeptonTreesFromMiniAODPuppi::initIntBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    intBranches_[(*it).first][name] = new unsigned int;
    (*it).second->Branch(name.c_str(), intBranches_[(*it).first][name], (name+"/I").c_str());
  }
}

DiLeptonTreesFromMiniAODPuppi::~DiLeptonTreesFromMiniAODPuppi()
{ 
  for( std::map<std::string, std::map< std::string, float*> >::const_iterator it = floatBranches_.begin();
       it != floatBranches_.end(); ++it){
    for( std::map< std::string, float*>::const_iterator it2 = (*it).second.begin();
	 it2 != (*it).second.end(); ++it2){
      if(debug)std::cout << "deleting: " << (*it).first << " - "<< (*it2).first << std::endl;
      delete (*it2).second;
    }
  }
  for( std::map<std::string, std::map< std::string, unsigned int*> >::const_iterator it = intBranches_.begin();
       it != intBranches_.end(); ++it){
    for( std::map< std::string, unsigned int*>::const_iterator it2 = (*it).second.begin();
	 it2 != (*it).second.end(); ++it2){
      if(debug) std::cout << "deleting: " << (*it).first << " - "<< (*it2).first << std::endl;
      delete (*it2).second;
    }
  }


}


// member functions
// ------------ method called to for each event  ------------
void
DiLeptonTreesFromMiniAODPuppi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< std::vector< pat::Electron > > electrons;
  iEvent.getByToken(electronToken_, electrons);
  
  edm::Handle< std::vector< pat::Electron > > looseElectrons;
  iEvent.getByToken(looseElectronToken_, looseElectrons);

  edm::Handle< std::vector< pat::Muon > > muons;
  iEvent.getByToken(muonToken_, muons);
  
  edm::Handle< std::vector< pat::Muon > > looseMuons;
  iEvent.getByToken(looseMuonToken_, looseMuons);


  edm::Handle< std::vector<pat::PackedCandidate>  > pfCands;
  //~ iEvent.getByLabel(pfCandTag_, pfCands); 
  iEvent.getByToken(pfCandToken_, pfCands); 

  
  //~ edm::Handle< std::vector< pat::Jet > > jets;
  iEvent.getByToken(jetToken_, jets);

  edm::Handle< std::vector< reco::GenJet  > > genJets;
  iEvent.getByToken(genJetToken_, genJets);

  edm::Handle< std::vector< pat::Jet > > jetsPuppi;
  iEvent.getByToken(jetPuppiToken_, jetsPuppi);

  edm::Handle< std::vector< pat::Jet > > bJets;
  iEvent.getByToken(bJetToken_, bJets);
  
  edm::Handle< std::vector< pat::Jet > > bJets35;
  iEvent.getByToken(bJet35Token_, bJets35);
  
  edm::Handle< std::vector< pat::Jet > > bJetsPuppi;
  iEvent.getByToken(bJetPuppiToken_, bJetsPuppi);
  
  edm::Handle< std::vector< pat::MET > > mets;
  iEvent.getByToken(metToken_, mets);
  
  edm::Handle< std::vector< pat::MET > > metsPuppi;
  iEvent.getByToken(metPuppiToken_, metsPuppi);


  edm::Handle< std::vector< reco::GenParticle > > genParticles;
  iEvent.getByToken(genParticleToken_, genParticles);
  
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::InputTag triggerTag("TriggerResults","","HLT");
  iEvent.getByLabel(triggerTag, triggerBits);
  
  // for each lumiBlock, re-read the trigger indices (rather changes for new run)
   if( triggerIndex_.size() && newLumiBlock_ ) {
      newLumiBlock_=false;
      // set all trigger indeces to -1 as "not available"-flag
      for( auto& it : triggerIndex_ )
        it.second = -1;

      // store the indices of the trigger names that we really find
      const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
      for( unsigned i=0; i<triggerNames.size(); i++ ) {
         for( auto& it : triggerIndex_ ) {
            if( triggerNames.triggerName(i).find( it.first ) == 0 ) {
               it.second = i;
            }
         }
      } // end trigger names
   } // found indices
   
   // set trigger decision
   for( auto& it : triggerIndex_ ) {
      if( it.second != -1 ) {
         triggerDecision_[it.first] = triggerBits->accept( it.second );
      }
   }
   
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);

  getPdgId_.loadGenParticles(iEvent);
  fctIsolation_.init(iEvent);
  fctTrigger_.loadTrigger(iEvent);
  std::map<std::string, int> intEventProperties;
  std::map<std::string, float> floatEventProperties;
  std::map<std::string, TLorentzVector> tLorentzVectorEventProperties;


  edm::Handle<GenEventInfoProduct> genInfoProduct;
  iEvent.getByToken(genEventInfoToken_, genInfoProduct);	
  if (genInfoProduct.isValid()){
  	floatEventProperties["genWeightAbsValue"] = (*genInfoProduct).weight();
   	if ((*genInfoProduct).weight() < 0.0){
   	
   		floatEventProperties["genWeight"] = -1;
   	}
   	else{
   		floatEventProperties["genWeight"] = 1;   	
   	}
  }
  else{
 
 	floatEventProperties["genWeight"] = 1; 
 	floatEventProperties["genWeightAbsValue"] = 1;
  
  }	  
  
	// stolen from https://github.com/Aachen-3A/PxlSkimmer/blob/master/Skimming/src/PxlSkimmer_miniAOD.cc#L590
    edm::Handle<LHEEventProduct> lheInfoHandle;
    iEvent.getByToken(LHEEventToken_ , lheInfoHandle);

    if (lheInfoHandle.isValid()) {
        lhef::HEPEUP lheParticleInfo = lheInfoHandle->hepeup();
        // get the five vector
        // (Px, Py, Pz, E and M in GeV)
        std::vector<lhef::HEPEUP::FiveVector> allParticles = lheParticleInfo.PUP;
        std::vector<int> statusCodes = lheParticleInfo.ISTUP;

        double ht = 0;
        double et = 0;
        for (unsigned int i = 0; i < statusCodes.size(); i++) {
            if (statusCodes[i] == 1) {
				if (abs(lheParticleInfo.IDUP[i]) != 12 && abs(lheParticleInfo.IDUP[i]) != 14 && abs(lheParticleInfo.IDUP[i]) != 16) {
					et += sqrt(pow(allParticles[i][0], 2) + pow(allParticles[i][1], 2));
				}
                if (abs(lheParticleInfo.IDUP[i]) < 11 || abs(lheParticleInfo.IDUP[i]) > 16) {
                    ht += sqrt(pow(allParticles[i][0], 2) + pow(allParticles[i][1], 2));
                }
            }
        }
        floatEventProperties["genParticleET"] = et;
        floatEventProperties["genParticleHT"] = ht;
    }
	
	else {
		floatEventProperties["genParticleET"] = -999.;
		floatEventProperties["genParticleHT"] = -999.;
	}
	
  double packedET = 0;
  for (std::vector<pat::PackedCandidate>::const_iterator it = pfCands->begin(); it != pfCands->end(); it++) {
	  packedET += (*it).pt();
  }
  floatEventProperties["packedET"] = packedET;
	  

      	

	
  intEventProperties["nVertices"] = vertices->size();
  
  int nGenVertices = -1;
 
  
  if (genParticles.isValid()){	  
	  edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
	  iEvent.getByLabel(edm::InputTag("slimmedAddPileupInfo"), PupInfo);
	
	  std::vector<PileupSummaryInfo>::const_iterator PVI;

  
	  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
		
		int BX = PVI->getBunchCrossing();
		
		if(BX == 0) { 
		  nGenVertices = PVI->getTrueNumInteractions();
		  continue;
		}
	  }
  }
  
  intEventProperties["nGenVertices"] = nGenVertices;
    

  edm::Handle<double> rho;
  iEvent.getByToken(rhoToken_,rho);
  floatEventProperties["rho"] = (float)(*rho);
  const float Rho = *rho;


  intEventProperties["nBJets"] = bJets->size();
  intEventProperties["nBJets35"] = bJets35->size();
  intEventProperties["nBJetsPuppi"] = bJetsPuppi->size();
  intEventProperties["nLightLeptons"] = electrons->size() + muons->size();
  intEventProperties["nLooseLeptons"] = looseElectrons->size() + looseMuons->size();
  intEventProperties["runNr"] = iEvent.id().run();
  intEventProperties["lumiSec"] = iEvent.id().luminosityBlock();
  intEventProperties["eventNr"] = iEvent.id().event();
	  



  pat::MET met = mets->front();
  TLorentzVector metVector(mets->front().px(), mets->front().py(), mets->front().pz(), mets->front().energy());
  TLorentzVector uncorrectedMetVector;
  uncorrectedMetVector.SetPtEtaPhiE(mets->front().uncorPt(), 0,	mets->front().uncorPhi(), mets->front().uncorPt());
 
  
  floatEventProperties["met"] = metVector.Pt();
  tLorentzVectorEventProperties["vMet"] = metVector; 
  
  tLorentzVectorEventProperties["vMetUncorrected"] = uncorrectedMetVector;
  floatEventProperties["uncorrectedMet"] = uncorrectedMetVector.Pt();
  
  
  pat::MET metPuppi = metsPuppi->front();
  TLorentzVector metVectorPuppi(metsPuppi->front().px(), metsPuppi->front().py(), metsPuppi->front().pz(), metsPuppi->front().energy());
  
  floatEventProperties["metPuppi"] = metVectorPuppi.Pt();
  tLorentzVectorEventProperties["vMetPuppi"] = metVectorPuppi; 
  

  pat::METCollection const& metsForUncert = *mets;	
  	
  if (metUncert_){

 	floatEventProperties["met"] =  metsForUncert[0].pt();
 	floatEventProperties["metJetEnUp"] = metsForUncert[0].shiftedPt(pat::MET::JetEnUp); 
 	floatEventProperties["metJetEnDown"] = metsForUncert[0].shiftedPt(pat::MET::JetEnDown); 
 	floatEventProperties["metJetResUp"] = metsForUncert[0].shiftedPt(pat::MET::JetResUp);
 	floatEventProperties["metJetResDown"] = metsForUncert[0].shiftedPt(pat::MET::JetResDown);
 	floatEventProperties["metMuonEnUp"] = metsForUncert[0].shiftedPt(pat::MET::MuonEnUp);
 	floatEventProperties["metMuonEnDown"] = metsForUncert[0].shiftedPt(pat::MET::MuonEnDown);
 	floatEventProperties["metElectronEnUp"] = metsForUncert[0].shiftedPt(pat::MET::ElectronEnUp);
 	floatEventProperties["metElectronEnDown"] = metsForUncert[0].shiftedPt(pat::MET::ElectronEnDown);
 	floatEventProperties["metTauEnUp"] = metsForUncert[0].shiftedPt(pat::MET::TauEnUp);
 	floatEventProperties["metTauEnDown"] = metsForUncert[0].shiftedPt(pat::MET::TauEnDown);
 	floatEventProperties["metUnclusteredEnUp"] = metsForUncert[0].shiftedPt(pat::MET::UnclusteredEnUp);
 	floatEventProperties["metUnclusteredEnDown"] = metsForUncert[0].shiftedPt(pat::MET::UnclusteredEnDown);  
  
  }
  
  TLorentzVector genMetVector(0.,0.,0.,0.);
  TLorentzVector vGenParticle(0.,0.,0.,0.);
  
  floatEventProperties["genPtTop1"] = -1;
  floatEventProperties["genPtTop2"] = -1;
  if (genParticles.isValid()){
	  	
	genMetVector.SetPxPyPzE(mets->front().genMET()->px(),mets->front().genMET()->py(),mets->front().genMET()->pz(),mets->front().genMET()->energy());
	
	for (std::vector<reco::GenParticle>::const_iterator itGenParticle = genParticles->begin(); itGenParticle != genParticles->end(); itGenParticle++) {

		if (abs((*itGenParticle).pdgId())== 6){

			if ((*itGenParticle).pdgId()== 6){
				floatEventProperties["genPtTop1"] = (*itGenParticle).pt();
			}
			else if ((*itGenParticle).pdgId()== -6){
				floatEventProperties["genPtTop2"] = (*itGenParticle).pt();
			}

		}


	}

  }
  
  floatEventProperties["genMet"] = genMetVector.Pt();
  tLorentzVectorEventProperties["vGenMet"] = genMetVector;  
  
  TLorentzVector MHT;
  TLorentzVector MHT10;
  TLorentzVector MHT20;  
  TLorentzVector MHT30;  
  TLorentzVector MHTPuppi;
  TLorentzVector MHT10Puppi;  
  TLorentzVector MHT20Puppi;  
  TLorentzVector MHT30Puppi;  
  TLorentzVector genMHT;  
  TLorentzVector genMHT10;  
  TLorentzVector genMHT20;  
  TLorentzVector genMHT30;  
  TLorentzVector tempMHT;
  TLorentzVector tempMHT10;
  TLorentzVector tempMHT20;
  TLorentzVector tempMHT30;
  TLorentzVector tempGenMHT;
  
  floatEventProperties["ht"] = 0.0;
  floatEventProperties["ht40"] = 0.0;
  floatEventProperties["leptHT"] = 0.0;  
  floatEventProperties["genLeptHT"] = 0.0;  
  floatEventProperties["leptHTPuppi"] = 0.0;  
  floatEventProperties["htPuppi"] = 0.0;
  floatEventProperties["htJESUp"] = 0.0;
  floatEventProperties["htJESDown"] = 0.0;
  floatEventProperties["mht"] = 0.0;
  floatEventProperties["mhtPuppi"] = 0.0;

  edm::FileInPath fip("SuSyAachen/DiLeptonHistograms/data/GR_P_V40_AN1::All_Uncertainty_AK5PF.txt");
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(fip.fullPath());

  //edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  //iSetup.get<JetCorrectionsRecord>().get("AK4PF",JetCorParColl); 
  //JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  //JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);


  TLorentzVector changeUp(0,0,0,0);
  TLorentzVector changeDown(0,0,0,0);
  TLorentzVector tempChangeUp;
  TLorentzVector tempChangeDown;
  // loop over jets
  std::vector<pat::Jet> * shiftedJetsJESUp = new std::vector<pat::Jet>(); 
  std::vector<pat::Jet> * shiftedJetsJESDown = new std::vector<pat::Jet>(); 



  for (std::vector<pat::Jet>::const_iterator itJet = jets->begin(); itJet != jets->end(); itJet++) {
	
	

    pat::Jet ajetUp( *itJet );
    pat::Jet ajetDown( *itJet );
    jecUnc->setJetEta(itJet->eta());
    jecUnc->setJetPt(itJet->pt()); // IMPORTANT: the uncertainty is a function of the CORRECTED pt
    double correction = jecUnc->getUncertainty(true);
    //double correction = itJet->relCorrUncert(dir_); 
    
    //use pat::Jet::relCorrUncert
    if (correction==0.0){
	correction = 0.045; 
    }
    //    std::cout << "jet pt = "<<itJet->p4()<<", uncertainty = "<<correction<<std::endl;
    ajetUp.setP4( itJet->p4() * (1.0 + correction) );
    ajetDown.setP4( itJet->p4() * (1.0 - correction) );

    //keep track of the change, in order to correct JEC-corrected MET, too.
    tempChangeUp.SetPxPyPzE((ajetUp.px()-(*itJet).px()), (ajetUp.py()-(*itJet).py()), (ajetUp.pz()-(*itJet).pz()), (ajetUp.energy()-(*itJet).energy()));	
    tempChangeDown.SetPxPyPzE((ajetDown.px()-itJet->px()), (ajetDown.py()-itJet->py()), (ajetDown.pz()-itJet->pz()), (ajetDown.energy()-itJet->energy()));	
    changeUp +=tempChangeUp;
    changeDown += tempChangeDown;

    shiftedJetsJESUp->push_back(ajetUp);
    shiftedJetsJESDown->push_back(ajetDown);
  }
  int nJets=0;
  int nJets10=0;
  int nJets20=0;
  int nJets30=0;
  int matchedJets=0;
  int nGenJets=0;
  int nGenJets10=0;
  int nGenJets20=0;
  int nGenJets30=0;
  int nJetsPuppi=0;
  int nJetsPuppi10=0;
  int nJetsPuppi20=0;
  int nJetsPuppi30=0;
  int matchedJetsPuppi=0;
  //int nJetsNoPULoose = 0;
  //int nJetsNoPUMedium = 0;
  //int nJetsNoPUTight = 0;
  
  //float puJetID = 0.;	
  
  for(std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end() ; ++it){
	if ((*it).pt() >=35.0 && fabs((*it).eta())<2.4){
		nJets++;
		
		//puJetID =  it->userFloat("pileupJetId:fullIdLoose");
		//std::cout << puJetID << endl;
      	//if( PileupJetIdentifier::passJetId( puJetID, PileupJetIdentifier::kLoose )) {
        //   nJetsNoPULoose++;
      	//}
     	//if( PileupJetIdentifier::passJetId( puJetID, PileupJetIdentifier::kMedium )) {
        //   nJetsNoPUMedium++;
        //}
        //if( PileupJetIdentifier::passJetId( puJetID, PileupJetIdentifier::kTight )) {
        //   nJetsNoPUTight++;		
		//}
		tempMHT.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());	
		MHT = MHT + tempMHT;
		floatEventProperties["ht"] += (*it).pt();
		
		if (genParticles.isValid()){
			
			TLorentzVector jetVector( (*it).px(), (*it).py(), (*it).pz(), (*it).energy() );
			//~ 
			for(std::vector<reco::GenJet >::const_iterator itGenJet = genJets->begin(); itGenJet != genJets->end() ; ++itGenJet){
				if ( (*itGenJet).pileup()) continue;
				TLorentzVector genJetVector( (*itGenJet).px(), (*itGenJet).py(), (*itGenJet).pz(), (*itGenJet).energy() );
				if ( fabs((*it).pt()-(*itGenJet).pt())/(*it).pt() < 0.2 && jetVector.DeltaR(genJetVector) < 0.1) {
					matchedJets++;
				}
			}
			//~ if ( (*it).genJet())	matchedJets++;
		}
	}
	
	if ((*it).pt() >=20.0 && fabs((*it).eta())<5.0){
		tempMHT20.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());	
		MHT20 = MHT20 + tempMHT20;
		nJets20++;
	}	
	if ((*it).pt() >=10.0 ){
		tempMHT10.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());
		MHT10 = MHT10 + tempMHT10;
		nJets10++;
	}
	if ((*it).pt() >=30.0 ){
		tempMHT30.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());
		MHT30 = MHT30 + tempMHT30;
		nJets30++;
	}
	
  }
  intEventProperties["nJets"] = nJets;
  intEventProperties["nJets10"] = nJets10;
  intEventProperties["nJets20"] = nJets20;
  intEventProperties["nJets30"] = nJets30;
  intEventProperties["nMatchedJets"] = matchedJets;
  
  floatEventProperties["metHadRecoil"] = (MHT.Px()*metVector.Px()+MHT.Py()*metVector.Py())/sqrt(MHT.Px()*MHT.Px() + MHT.Py()*MHT.Py());
  floatEventProperties["metHadRecoil10"] = (MHT10.Px()*metVector.Px()+MHT10.Py()*metVector.Py())/sqrt(MHT10.Px()*MHT10.Px() + MHT10.Py()*MHT10.Py());
  floatEventProperties["metHadRecoil20"] = (MHT20.Px()*metVector.Px()+MHT20.Py()*metVector.Py())/sqrt(MHT20.Px()*MHT20.Px() + MHT20.Py()*MHT20.Py());
  floatEventProperties["metHadRecoil30"] = (MHT30.Px()*metVector.Px()+MHT30.Py()*metVector.Py())/sqrt(MHT30.Px()*MHT30.Px() + MHT30.Py()*MHT30.Py());
  
  float genHT = 0.;
  
  if (genParticles.isValid()){
	  for(std::vector<reco::GenJet >::const_iterator it = genJets->begin(); it != genJets->end() ; ++it){
		  
		if ((*it).pt() >=35.0 && fabs((*it).eta())<2.4){
			
			bool leptonClean = true;
			TLorentzVector genJetVector( (*it).px(), (*it).py(), (*it).pz(), (*it).energy() );
			
			for(std::vector< pat::Electron >::const_iterator itEle = electrons->begin(); itEle != electrons->end() ; ++itEle){
				
				if (leptonClean == false) break;
				
				if (getIso((*itEle),"miniIsoEA")<0.1) {
					TLorentzVector eleVector( (*itEle).px(), (*itEle).py(), (*itEle).pz(), (*itEle).energy() );
					if (genJetVector.DeltaR( eleVector ) < 0.4){
						leptonClean = false;
					}
				}
			}
			
			for(std::vector< pat::Muon >::const_iterator itMu = muons->begin(); itMu != muons->end() ; ++itMu){
				
				if (leptonClean == false) break;
				
				if (getIso((*itMu),"miniIsoEA")<0.2) {
					TLorentzVector muVector( (*itMu).px(), (*itMu).py(), (*itMu).pz(), (*itMu).energy() );
					if (genJetVector.DeltaR( muVector ) < 0.4){
						leptonClean = false;
					}
				}
			}
			
			if (leptonClean == true) {
				nGenJets++;
				genHT += (*it).pt();
				tempGenMHT.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());	
				genMHT = genMHT + tempGenMHT;
			}
			
		}	
		
	  }
	  
  }
  intEventProperties["nGenJets"] = nGenJets;
  floatEventProperties["genHT"] = genHT;
  floatEventProperties["genMetHadRecoil"] = (genMHT.Px()*genMetVector.Px()+genMHT.Py()*genMetVector.Py())/sqrt(genMHT.Px()*genMHT.Px() + genMHT.Py()*genMHT.Py());
  
  
  
  for(std::vector<pat::Jet>::const_iterator it = jetsPuppi->begin(); it != jetsPuppi->end() ; ++it){
	if ((*it).pt() >=35.0 && fabs((*it).eta())<2.4){
		nJetsPuppi++;
		
		tempMHT.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());	
		MHTPuppi = MHTPuppi + tempMHT;
		floatEventProperties["htPuppi"] += (*it).pt();
		
		if (genParticles.isValid()){
			
			//~ if ( (*it).genJet())	matchedJets++;
			//~ 
			TLorentzVector jetVector( (*it).px(), (*it).py(), (*it).pz(), (*it).energy() );
						
			for(std::vector<reco::GenJet >::const_iterator itGenJet = genJets->begin(); itGenJet != genJets->end() ; ++itGenJet){
				TLorentzVector genJetVector( (*itGenJet).px(), (*itGenJet).py(), (*itGenJet).pz(), (*itGenJet).energy() );
				if ( fabs((*it).pt()-(*itGenJet).pt())/(*it).pt() < 0.2 && jetVector.DeltaR(genJetVector) < 0.1) {
					matchedJetsPuppi++;
				}
			}
		}
	}
	if ((*it).pt() >=10.0){
		tempMHT10.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());
		MHT10Puppi = MHT10Puppi + tempMHT10;
		nJetsPuppi10++;
	}	
	if ((*it).pt() >=20.0 && fabs((*it).eta())<5.0){
		tempMHT20.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());	
		MHT20Puppi = MHT20Puppi + tempMHT20;
		nJetsPuppi20++;
	}	
	
	if ((*it).pt() >=30.0){
		tempMHT30.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());
		MHT30Puppi = MHT30Puppi + tempMHT30;
		nJetsPuppi30++;
	}
	
  }
  intEventProperties["nJetsPuppi"] = nJetsPuppi;
  intEventProperties["nJetsPuppi10"] = nJetsPuppi10;
  intEventProperties["nJetsPuppi20"] = nJetsPuppi20;
  intEventProperties["nJetsPuppi30"] = nJetsPuppi30;
  intEventProperties["nMatchedJetsPuppi"] = matchedJetsPuppi;
  floatEventProperties["metHadRecoilPuppi"] = (MHTPuppi.Px()*metVectorPuppi.Px()+MHTPuppi.Py()*metVectorPuppi.Py())/sqrt(MHTPuppi.Px()*MHTPuppi.Px() + MHTPuppi.Py()*MHTPuppi.Py());
  floatEventProperties["metHadRecoil10Puppi"] = (MHT10Puppi.Px()*metVectorPuppi.Px()+MHT10Puppi.Py()*metVectorPuppi.Py())/sqrt(MHT10Puppi.Px()*MHT10Puppi.Px() + MHT10Puppi.Py()*MHT10Puppi.Py());
  floatEventProperties["metHadRecoil20Puppi"] = (MHT20Puppi.Px()*metVectorPuppi.Px()+MHT20Puppi.Py()*metVectorPuppi.Py())/sqrt(MHT20Puppi.Px()*MHT20Puppi.Px() + MHT20Puppi.Py()*MHT20Puppi.Py());
  floatEventProperties["metHadRecoil30Puppi"] = (MHT30Puppi.Px()*metVectorPuppi.Px()+MHT30Puppi.Py()*metVectorPuppi.Py())/sqrt(MHT30Puppi.Px()*MHT30Puppi.Px() + MHT30Puppi.Py()*MHT30Puppi.Py());
  
  //~ nJets30 = 0;
  //~ 
  //~ for(std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end() ; ++it){
	//~ if ((*it).pt() >=30.0 && fabs((*it).eta())<2.4){
		//~ nJets30++;
	//~ }
  //~ }
  //~ intEventProperties["nJets30"] = nJets30;  
  //~ 
  //~ //intEventProperties["nJetsNoPULoose"] = nJetsNoPULoose;
  //~ //intEventProperties["nJetsNoPUMedium"] = nJetsNoPUMedium;
  //~ //intEventProperties["nJetsNoPUTight"] = nJetsNoPUTight;    
  int nJetsOld=0;
  for(std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end() ; ++it){
	if ((*it).pt() >=40.0){
		nJetsOld++;
		floatEventProperties["leptHT"] += (*it).pt();	
		floatEventProperties["ht40"] += (*it).pt();						
	}
  }
  intEventProperties["nJetsOld"] = nJetsOld;
  
  for(std::vector<pat::Jet>::const_iterator it = jetsPuppi->begin(); it != jetsPuppi->end() ; ++it){
	if ((*it).pt() >=40.0){
		floatEventProperties["leptHTPuppi"] += (*it).pt();					
	}
  }


  TLorentzVector jet1Vector(0.,0.,0.,0.);
  TLorentzVector jet2Vector(0.,0.,0.,0.);
  TLorentzVector jet3Vector(0.,0.,0.,0.);
  TLorentzVector jet4Vector(0.,0.,0.,0.); 


  if (nJets > 0)
    jet1Vector.SetPxPyPzE(jets->at(0).px(),jets->at(0).py(),jets->at(0).pz(),jets->at(0).energy());
  if (nJets > 1)
    jet2Vector.SetPxPyPzE(jets->at(1).px(),jets->at(1).py(),jets->at(1).pz(),jets->at(1).energy());
  if (nJets > 2)
    jet3Vector.SetPxPyPzE(jets->at(2).px(),jets->at(2).py(),jets->at(2).pz(),jets->at(2).energy());
  if (nJets > 3)
    jet4Vector.SetPxPyPzE(jets->at(3).px(),jets->at(3).py(),jets->at(3).pz(),jets->at(3).energy());
  tLorentzVectorEventProperties["jet1"] = jet1Vector;
  tLorentzVectorEventProperties["jet2"] = jet2Vector;
  tLorentzVectorEventProperties["jet3"] = jet3Vector;
  tLorentzVectorEventProperties["jet4"] = jet4Vector;

  TLorentzVector jet1VectorPuppi(0.,0.,0.,0.);
  TLorentzVector jet2VectorPuppi(0.,0.,0.,0.);
  TLorentzVector jet3VectorPuppi(0.,0.,0.,0.);
  TLorentzVector jet4VectorPuppi(0.,0.,0.,0.); 

  if (nJetsPuppi > 0)
    jet1VectorPuppi.SetPxPyPzE(jetsPuppi->at(0).px(),jetsPuppi->at(0).py(),jetsPuppi->at(0).pz(),jetsPuppi->at(0).energy());
  if (nJetsPuppi > 1)
    jet2VectorPuppi.SetPxPyPzE(jetsPuppi->at(1).px(),jetsPuppi->at(1).py(),jetsPuppi->at(1).pz(),jetsPuppi->at(1).energy());
  if (nJetsPuppi > 2)
    jet3VectorPuppi.SetPxPyPzE(jetsPuppi->at(2).px(),jetsPuppi->at(2).py(),jetsPuppi->at(2).pz(),jetsPuppi->at(2).energy());
  if (nJetsPuppi > 3)
    jet4VectorPuppi.SetPxPyPzE(jetsPuppi->at(3).px(),jetsPuppi->at(3).py(),jetsPuppi->at(3).pz(),jetsPuppi->at(3).energy());
  tLorentzVectorEventProperties["jet1Puppi"] = jet1VectorPuppi;
  tLorentzVectorEventProperties["jet2Puppi"] = jet2VectorPuppi;
  tLorentzVectorEventProperties["jet3Puppi"] = jet3VectorPuppi;
  tLorentzVectorEventProperties["jet4Puppi"] = jet4VectorPuppi;

  TLorentzVector genJet1Vector(0.,0.,0.,0.);
  TLorentzVector genJet2Vector(0.,0.,0.,0.);
  TLorentzVector genJet3Vector(0.,0.,0.,0.);
  TLorentzVector genJet4Vector(0.,0.,0.,0.); 

  if (nGenJets > 0)
    genJet1Vector.SetPxPyPzE(genJets->at(0).px(),genJets->at(0).py(),genJets->at(0).pz(),genJets->at(0).energy());
  if (nGenJets > 1)
    genJet2Vector.SetPxPyPzE(genJets->at(1).px(),genJets->at(1).py(),genJets->at(1).pz(),genJets->at(1).energy());
  if (nGenJets > 2)
    genJet3Vector.SetPxPyPzE(genJets->at(2).px(),genJets->at(2).py(),genJets->at(2).pz(),genJets->at(2).energy());
  if (nGenJets > 3)
    genJet4Vector.SetPxPyPzE(genJets->at(3).px(),genJets->at(3).py(),genJets->at(3).pz(),genJets->at(3).energy());
  tLorentzVectorEventProperties["genJet1"] = genJet1Vector;
  tLorentzVectorEventProperties["genJet2"] = genJet2Vector;
  tLorentzVectorEventProperties["genJet3"] = genJet3Vector;
  tLorentzVectorEventProperties["genJet4"] = genJet4Vector;




  int nJetsJESUp=0;
  for(std::vector<pat::Jet>::const_iterator it = shiftedJetsJESUp->begin(); it != shiftedJetsJESUp->end() ; ++it){
	if ((*it).pt() >=35.0){
        	nJetsJESUp++;
		floatEventProperties["htJESUp"] += (*it).pt();
	}
  }
  intEventProperties["nShiftedJetsJESUp"] = nJetsJESUp;
  int nJetsJESDown=0;
  for(std::vector<pat::Jet>::const_iterator it = shiftedJetsJESDown->begin(); it != shiftedJetsJESDown->end() ; ++it){
	if ((*it).pt() >=35.0){	
        	nJetsJESDown++;
		floatEventProperties["htJESDown"] += (*it).pt();
	}
  }
  intEventProperties["nShiftedJetsJESDown"] = nJetsJESDown;
  TLorentzVector metVectorJESUp = metVector + changeUp;
  TLorentzVector metVectorJESDown = metVector + changeDown;
  floatEventProperties["metJESUp"] = metVectorJESUp.Pt();


  floatEventProperties["metJESDown"] = metVectorJESDown.Pt();

//   floatEventProperties["metJESUp"] = 0.;
// 
// 
//   floatEventProperties["metJESDown"] = 0.;


//   floatEventProperties["mht"] = MHT.Pt();

  // Jet pt
  floatEventProperties["jet1pt"] = -1.0;
  floatEventProperties["jet2pt"] = -1.0;
  floatEventProperties["jet3pt"] = -1.0;
  floatEventProperties["jet4pt"] = -1.0;
  if (jets->size() > 0)
    floatEventProperties["jet1pt"] = jets->at(0).pt();
  if (jets->size() > 1)
    floatEventProperties["jet2pt"] = jets->at(1).pt();
  if (jets->size() > 2)
    floatEventProperties["jet3pt"] = jets->at(2).pt();
  if (jets->size() > 3)
    floatEventProperties["jet4pt"] = jets->at(3).pt();

  floatEventProperties["genJet1pt"] = -1.0;
  floatEventProperties["genJet2pt"] = -1.0;
  floatEventProperties["genJet3pt"] = -1.0;
  floatEventProperties["genJet4pt"] = -1.0;
  
  
  if (genParticles.isValid()){
	  int genJet = 0;
	  for(std::vector<reco::GenJet >::const_iterator it = genJets->begin(); it != genJets->end() ; ++it){
		if ((*it).pt() >=10.0 && fabs((*it).eta())<3.0){	
					
			bool leptonClean = true;
			TLorentzVector genJetVector( (*it).px(), (*it).py(), (*it).pz(), (*it).energy() );
			
			for(std::vector< pat::Electron >::const_iterator itEle = electrons->begin(); itEle != electrons->end() ; ++itEle){
				
				if (leptonClean == false) break;
				
				if (getIso((*itEle),"miniIsoEA")<0.1) {
					TLorentzVector eleVector( (*itEle).px(), (*itEle).py(), (*itEle).pz(), (*itEle).energy() );
					if (genJetVector.DeltaR( eleVector ) < 0.4){
						leptonClean = false;
					}
				}
			}
			
			for(std::vector< pat::Muon >::const_iterator itMu = muons->begin(); itMu != muons->end() ; ++itMu){
				
				if (leptonClean == false) break;
				
				if (getIso((*itMu),"miniIsoEA")<0.2) {
					TLorentzVector muVector( (*itMu).px(), (*itMu).py(), (*itMu).pz(), (*itMu).energy() );
					if (genJetVector.DeltaR( muVector ) < 0.4){
						leptonClean = false;
					}
				}
			}
			
			if (leptonClean == true) {
				genJet++;
				nGenJets10++;
				if (genJet == 1) floatEventProperties["genJet1pt"] =(*it).pt();
				if (genJet == 2) floatEventProperties["genJet2pt"] =(*it).pt();
				if (genJet == 3) floatEventProperties["genJet3pt"] =(*it).pt();
				if (genJet == 4) floatEventProperties["genJet4pt"] =(*it).pt();
				
				if ((*it).pt() >=10.0){
					tempMHT10.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());	
					genMHT10 = genMHT10 + tempMHT10;
					nGenJets10++;
				}	
				if ((*it).pt() >=20.0){
					tempMHT20.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());	
					genMHT20 = genMHT20 + tempMHT20;
					nGenJets20++;
				}	
				if ((*it).pt() >=30.0){
					tempMHT30.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());	
					genMHT30 = genMHT30 + tempMHT30;
					nGenJets30++;
				}
			}	
		}	
		
	  }
	  
  }
  intEventProperties["nGenJets10"] = nGenJets10;
  intEventProperties["nGenJets20"] = nGenJets20;
  intEventProperties["nGenJets30"] = nGenJets30;
  
  floatEventProperties["genMetHadRecoil10"] = (genMHT10.Px()*genMetVector.Px()+genMHT10.Py()*genMetVector.Py())/sqrt(genMHT10.Px()*genMHT10.Px() + genMHT10.Py()*genMHT10.Py());
  floatEventProperties["genMetHadRecoil20"] = (genMHT20.Px()*genMetVector.Px()+genMHT20.Py()*genMetVector.Py())/sqrt(genMHT20.Px()*genMHT20.Px() + genMHT20.Py()*genMHT20.Py());
  floatEventProperties["genMetHadRecoil30"] = (genMHT30.Px()*genMetVector.Px()+genMHT30.Py()*genMetVector.Py())/sqrt(genMHT30.Px()*genMHT30.Px() + genMHT30.Py()*genMHT30.Py());
  
  
  floatEventProperties["jet1ptPuppi"] = -1.0;
  floatEventProperties["jet2ptPuppi"] = -1.0;
  floatEventProperties["jet3ptPuppi"] = -1.0;
  floatEventProperties["jet4ptPuppi"] = -1.0;
  if (jetsPuppi->size() > 0)
    floatEventProperties["jet1ptPuppi"] = jetsPuppi->at(0).pt();
  if (jetsPuppi->size() > 1)
    floatEventProperties["jet2ptPuppi"] = jetsPuppi->at(1).pt();
  if (jetsPuppi->size() > 2)
    floatEventProperties["jet3ptPuppi"] = jetsPuppi->at(2).pt();
  if (jetsPuppi->size() > 3)
    floatEventProperties["jet4ptPuppi"] = jetsPuppi->at(3).pt();

  // bjet pt

  TLorentzVector bJet1Vector(0.,0.,0.,0.);
  TLorentzVector bJet2Vector(0.,0.,0.,0.);



  floatEventProperties["bjet1pt"] = -1.0;
  floatEventProperties["bjet2pt"] = -1.0;
  if (bJets->size() > 0){
    floatEventProperties["bjet1pt"] = bJets->at(0).pt();
    bJet1Vector.SetPxPyPzE(bJets->at(0).px(),bJets->at(0).py(),bJets->at(0).pz(),bJets->at(0).energy());
  }
  if (bJets->size() > 1){
    floatEventProperties["bjet2pt"] = bJets->at(1).pt();
    bJet2Vector.SetPxPyPzE(bJets->at(1).px(),bJets->at(1).py(),bJets->at(1).pz(),bJets->at(1).energy());
  }
 
  tLorentzVectorEventProperties["bJet1"] = bJet1Vector;
  tLorentzVectorEventProperties["bJet2"] = bJet2Vector;
  
  
  TLorentzVector bJet1VectorPuppi(0.,0.,0.,0.);
  TLorentzVector bJet2VectorPuppi(0.,0.,0.,0.);

  floatEventProperties["bjet1ptPuppi"] = -1.0;
  floatEventProperties["bjet2ptPuppi"] = -1.0;
  if (bJetsPuppi->size() > 0){
    floatEventProperties["bjet1ptPuppi"] = bJetsPuppi->at(0).pt();
    bJet1VectorPuppi.SetPxPyPzE(bJetsPuppi->at(0).px(),bJetsPuppi->at(0).py(),bJetsPuppi->at(0).pz(),bJetsPuppi->at(0).energy());
  }
  if (bJetsPuppi->size() > 1){
    floatEventProperties["bjet2ptPuppi"] = bJetsPuppi->at(1).pt();
    bJet2VectorPuppi.SetPxPyPzE(bJetsPuppi->at(1).px(),bJetsPuppi->at(1).py(),bJetsPuppi->at(1).pz(),bJetsPuppi->at(1).energy());
  }
 
  tLorentzVectorEventProperties["bJet1Puppi"] = bJet1VectorPuppi;
  tLorentzVectorEventProperties["bJet2Puppi"] = bJet2VectorPuppi;
  
  floatEventProperties["ET"] = 0.0;  
  floatEventProperties["genET"] = 0.0;   
  floatEventProperties["ETPuppi"] = 0.0;
  
  for(std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end() ; ++it){
	if ((*it).pt() >=10.0 && fabs((*it).eta())<3.0){
		floatEventProperties["ET"] += (*it).pt();
	}	
  }
  
  for(std::vector<pat::Jet>::const_iterator it = jetsPuppi->begin(); it != jetsPuppi->end() ; ++it){
	if ((*it).pt() >=10.0 && fabs((*it).eta())<3.0){
		floatEventProperties["ETPuppi"] += (*it).pt();
	}	
  }
  
  if (genParticles.isValid()){
	  for(std::vector<reco::GenJet >::const_iterator it = genJets->begin(); it != genJets->end() ; ++it){
		  if ((*it).pt() >=10.0 && fabs((*it).eta())<3.0){
			  
			  bool leptonClean = true;
			  TLorentzVector genJetVector( (*it).px(), (*it).py(), (*it).pz(), (*it).energy() );
			  
			  for(std::vector< pat::Electron >::const_iterator itEle = electrons->begin(); itEle != electrons->end() ; ++itEle){
				  if (leptonClean == false) break;
				  
				  if (getIso((*itEle),"miniIsoEA")<0.1) {
					  TLorentzVector eleVector( (*itEle).px(), (*itEle).py(), (*itEle).pz(), (*itEle).energy() );
					  if (genJetVector.DeltaR( eleVector ) < 0.4){
						  leptonClean = false;
					  }
				  }
			  }
			  
			  for(std::vector< pat::Muon >::const_iterator itMu = muons->begin(); itMu != muons->end() ; ++itMu){
				  if (leptonClean == false) break;
				  
				  if (getIso((*itMu),"miniIsoEA")<0.2) {
					  TLorentzVector muVector( (*itMu).px(), (*itMu).py(), (*itMu).pz(), (*itMu).energy() );
					  if (genJetVector.DeltaR( muVector ) < 0.4){
						  leptonClean = false;
					  }
				  }
			  }
			  if (leptonClean == true) {
				  floatEventProperties["genET"] += (*it).pt();
				  if ((*it).pt() >=40.0 && fabs((*it).eta())<3.0){
					  floatEventProperties["genLeptHT"] += (*it).pt();
				  }
			  }
		  }
	  }
  }
  
  for(std::vector< pat::Electron >::const_iterator it = electrons->begin(); it != electrons->end() ; ++it){
	  floatEventProperties["ET"] += (*it).pt();
	  floatEventProperties["ETPuppi"] += (*it).pt();
	  if((*it).genLepton() != NULL) {
		  floatEventProperties["genET"] += (*it).genLepton()->pt();
	  }
  }
  for(std::vector< pat::Muon >::const_iterator it = muons->begin(); it != muons->end() ; ++it){
	  floatEventProperties["ET"] += (*it).pt();
	  floatEventProperties["ETPuppi"] += (*it).pt();
	  if((*it).genLepton() != NULL) {
		  floatEventProperties["genET"] += (*it).genLepton()->pt();
	  }
  }
			
  
  
  floatEventProperties["weight"] = fctVtxWeight_( iEvent );
  floatEventProperties["weightUp"] = fctVtxWeightUp_( iEvent );
  floatEventProperties["weightDown"] = fctVtxWeightDown_( iEvent );

  makeCombinations< pat::Electron >("EE", *electrons,*pfCands, *looseElectrons, *looseMuons, *jets, *bJets35, iEvent, met, metPuppi,MHT,MHT20,MHTPuppi,MHT20Puppi,vertices, Rho, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  makeCombinations< pat::Electron, pat::Muon >("EMu", *electrons, *muons,*pfCands, *looseElectrons, *looseMuons, *jets, *bJets35, iEvent, met, metPuppi,MHT,MHT20,MHTPuppi,MHT20Puppi,vertices, Rho, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  makeCombinations< pat::Muon >("MuMu", *muons,*pfCands, *looseElectrons, *looseMuons, *jets, *bJets35, iEvent, met, metPuppi,MHT,MHT20,MHTPuppi,MHT20Puppi,vertices, Rho, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  
  //  if( nMu != 2) std::cout << "-------! "<<nMu<<std::endl;

  delete jecUnc;
  delete shiftedJetsJESUp;
  delete shiftedJetsJESDown;


}

template <class aT, class bT> void 
DiLeptonTreesFromMiniAODPuppi::makeCombinations ( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT> &b,const std::vector<pat::PackedCandidate>&pfCands,const std::vector<pat::Electron>&looseElectrons,const std::vector<pat::Muon>&looseMuons,const std::vector<pat::Jet>&jets,const std::vector<pat::Jet>&bJets35, const edm::Event &ev, const pat::MET &patMet, const pat::MET &patMetPuppi, const TLorentzVector &MHT, const TLorentzVector &MHT20 , const TLorentzVector &MHTPuppi, const TLorentzVector &MHT20Puppi , const edm::Handle<reco::VertexCollection> &vertices, const float &rho , const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties)
{
  TLorentzVector met(patMet.px(), patMet.py(), patMet.pz(), patMet.energy());
  TLorentzVector uncorrectedMet2;
  uncorrectedMet2.SetPtEtaPhiE(patMet.uncorPt(), 0,\
                              patMet.uncorPhi(), patMet.uncorPt());
  *(tLorentzVectorBranches_[treeName]["vMet"]) = met;
  for(std::map<std::string, int>::const_iterator it = intEventProperties.begin(); it != intEventProperties.end(); ++it){
    assert(intBranches_[treeName].find((*it).first) != intBranches_[treeName].end());
    *(intBranches_[treeName][(*it).first]) = (*it).second;
  }
  for(std::map<std::string, float>::const_iterator it = floatEventProperties.begin(); it != floatEventProperties.end(); ++it){
    assert(floatBranches_[treeName].find((*it).first) != floatBranches_[treeName].end());
    *(floatBranches_[treeName][(*it).first]) = (*it).second;
  }

 for(std::map<std::string, TLorentzVector>::const_iterator it = tLorentzVectorEventProperties.begin(); it != tLorentzVectorEventProperties.end(); ++it){
    assert(tLorentzVectorBranches_[treeName].find((*it).first) != tLorentzVectorBranches_[treeName].end());
    *(tLorentzVectorBranches_[treeName][(*it).first]) = (*it).second;
  }
  
  //~ for ( std::vector<edm::ParameterSet>::iterator susyVar_i = susyVars_.begin(); susyVar_i != susyVars_.end(); ++susyVar_i ) {
        //~ edm::InputTag var = susyVar_i->getParameter<edm::InputTag>( "var" );
        //~ std::string type = susyVar_i->getParameter<std::string>( "type" );
        //~ edm::Handle< double > var_;
        //~ ev.getByToken(var, var_);
        //~ if (type=="float") *(floatBranches_[treeName][convertInputTag(var)]) = float(*var_);
        //~ else if (type=="int") *(intBranches_[treeName][convertInputTag(var)]) = int(*var_);
  //~ }
  //std::cout << std::endl;
  for( typename std::vector<aT>::const_iterator itA = a.begin(); itA != a.end(); ++itA){
    for( typename std::vector<bT>::const_iterator itB = b.begin(); itB != b.end(); ++itB){
//      std::cout << treeName <<": "<< fakeRates_(*itA) << std::endl;
//      weight *= fakeRates_();
      fillTree<aT,bT>( treeName, *itA, *itB,pfCands,looseElectrons,looseMuons,jets,bJets35, patMet, patMetPuppi,MHT,MHT20,MHTPuppi,MHT20Puppi,vertices,rho); 
    }
  }
}

template <class aT> void 
DiLeptonTreesFromMiniAODPuppi::makeCombinations ( const std::string &treeName, const std::vector<aT> &a,const std::vector<pat::PackedCandidate>&pfCands,const std::vector<pat::Electron>&looseElectrons,const std::vector<pat::Muon>&looseMuons,const std::vector<pat::Jet>&jets,const std::vector<pat::Jet>&bJets35,const edm::Event &ev, const pat::MET &patMet , const pat::MET &patMetPuppi , const TLorentzVector &MHT, const TLorentzVector &MHT20, const TLorentzVector &MHTPuppi, const TLorentzVector &MHT20Puppi, const edm::Handle<reco::VertexCollection> &vertices, const float &rho , const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties)
{
  TLorentzVector met(patMet.px(), patMet.py(), patMet.pz(), patMet.energy());
  TLorentzVector uncorrectedMet;
  uncorrectedMet.SetPtEtaPhiE(patMet.uncorPt(), 0,\
                              patMet.uncorPhi(), patMet.uncorPt());
  *(tLorentzVectorBranches_[treeName]["vMet"]) = met;
  for(std::map<std::string, int>::const_iterator it = intEventProperties.begin(); it != intEventProperties.end(); ++it){
    assert(intBranches_[treeName].find((*it).first) != intBranches_[treeName].end());
    *(intBranches_[treeName][(*it).first]) = (*it).second;  
  }
  for(std::map<std::string, float>::const_iterator it = floatEventProperties.begin(); it != floatEventProperties.end(); ++it){
    assert(floatBranches_[treeName].find((*it).first) != floatBranches_[treeName].end());
    *(floatBranches_[treeName][(*it).first]) = (*it).second;
  }

 for(std::map<std::string, TLorentzVector>::const_iterator it = tLorentzVectorEventProperties.begin(); it != tLorentzVectorEventProperties.end(); ++it){
    assert(tLorentzVectorBranches_[treeName].find((*it).first) != tLorentzVectorBranches_[treeName].end());
    *(tLorentzVectorBranches_[treeName][(*it).first]) = (*it).second;
  }

  //~ for ( std::vector<edm::ParameterSet>::iterator susyVar_i = susyVars_.begin(); susyVar_i != susyVars_.end(); ++susyVar_i ) {
        //~ edm::InputTag var = susyVar_i->getParameter<edm::InputTag>( "var" );
        //~ std::string type = susyVar_i->getParameter<std::string>( "type" );
        //~ edm::Handle< double > var_;
        //~ ev.getByLabel(var, var_);
        //~ if (type=="int") *(intBranches_[treeName][convertInputTag(var)]) = int(*var_);
        //~ else if (type=="float") *(floatBranches_[treeName][convertInputTag(var)]) = float(*var_);
  //~ }


  for( typename std::vector<aT>::const_iterator itA = a.begin(); itA != a.end(); ++itA){
    for( typename std::vector<aT>::const_iterator itB = a.begin(); itB != itA; ++itB){
      //std::cout << treeName<<"("<<(*itA).pt()<<", "<<(*itA).eta() <<"): "<< fakeRates_(*itA) << std::endl;
//      weight *= fakeRates_();
      fillTree<aT, aT>( treeName, *itA, *itB,pfCands,looseElectrons,looseMuons,jets,bJets35, patMet,patMetPuppi,MHT,MHT20,MHTPuppi,MHT20Puppi,vertices,rho); 
    }
  }
}

template <class aT, class bT> void 
DiLeptonTreesFromMiniAODPuppi::fillTree( const std::string &treeName, const aT& a, const bT& b,const std::vector<pat::PackedCandidate>&pfCands,const std::vector<pat::Electron>&looseElectrons,const std::vector<pat::Muon>&looseMuons,const std::vector<pat::Jet>&jets,const std::vector<pat::Jet>&bJets35, const pat::MET &patMet, const pat::MET &patMetPuppi, const TLorentzVector &MHT, const TLorentzVector &MHT20, const TLorentzVector &MHTPuppi, const TLorentzVector &MHT20Puppi, const edm::Handle<reco::VertexCollection> &vertices,const float &rho)
{
  if(debug) std::cout << treeName << "- pts:"<< a.pt() << " " << b.pt();
  TLorentzVector aVec = getMomentum(a);//( a.px(), a.py(), a.pz(), a.energy() );
  TLorentzVector bVec = getMomentum(b); //( b.px(), b.py(), b.pz(), b.energy() );
  TLorentzVector met(patMet.px(), patMet.py(), patMet.pz(), patMet.energy());
  TLorentzVector metPuppi(patMetPuppi.px(), patMetPuppi.py(), patMetPuppi.pz(), patMetPuppi.energy());
  TLorentzVector uncorrectedMet; 
  uncorrectedMet.SetPtEtaPhiE(patMet.uncorPt(), 0, \
			      patMet.uncorPhi(), patMet.uncorPt());  

  //  std::cout << "met: "<<met.Et()<< ", unCorr met: "<< uncorrectedMet.Et()
  //<< "=> "<< met.Et()* 1./uncorrectedMet.Et()<< " (xCheck: "<< patMet.corSumEt()*1./patMet.uncorrectedPt(pat::MET::uncorrALL) <<")"<<std::endl;

   *(floatBranches_[treeName]["leptHT"]) += aVec.Pt();	
   *(floatBranches_[treeName]["leptHT"]) += bVec.Pt();	

   *(floatBranches_[treeName]["leptHTPuppi"]) += aVec.Pt();	
   *(floatBranches_[treeName]["leptHTPuppi"]) += bVec.Pt();	

  TLorentzVector comb = aVec+bVec;
  TLorentzVector MHT2 = comb + MHT;
  TLorentzVector MHT2Loose = comb + MHT20;
  TLorentzVector MHT2Puppi = comb + MHTPuppi;
  TLorentzVector MHT2LoosePuppi = comb + MHT20Puppi;
  
  //~ float min_mlb_1 = 1.E6;
  //~ float min_mlb_2 = 1.E6;
  //~ 
  //~ TLorentzVector lb (0.,0.,0.,0.);
  //~ 
  //~ int leptonNumber = 0;
  //~ 
  //~ if (bJets35.size() > 1) {
	  //~ int bJetNumber = 0;
	  //~ int bJet1 = 0;
	  //~ 
	  //~ // Get minimal mlB
	  //~ for(std::vector<pat::Jet>::const_iterator it = bJets35.begin(); it != bJets35.end() ; ++it){
		  //~ bJetNumber += 1;
		  //~ TLorentzVector bJetVector1( (*it).px(), (*it).py(), (*it).pz(), (*it).energy() );
		  //~ TLorentzVector lb1 = aVec + bJetVector1;
		  //~ TLorentzVector lb2 = bVec + bJetVector1;
		  //~ if ( lb1.M() < lb2.M()) {
			  //~ if (lb1.M() < min_mlb_1){
				  //~ min_mlb_1 = lb1.M();
				  //~ bJet1 = bJetNumber;
				  //~ leptonNumber = 1;
			  //~ }
		  //~ }
		  //~ else {
			  //~ if (lb2.M() < min_mlb_1){
				  //~ min_mlb_1 = lb2.M();
				  //~ bJet1 = bJetNumber;
				  //~ leptonNumber = 2;
			  //~ }
		  //~ }
	  //~ }
	  //~ 
	  //~ // Get minimal mlb for second lepton
	  //~ bJetNumber = 0;
	  //~ for(std::vector<pat::Jet>::const_iterator it = bJets35.begin(); it != bJets35.end() ; ++it){
		  //~ bJetNumber += 1;
		  //~ if (bJetNumber == bJet1) continue;
		  //~ TLorentzVector bJetVector2( (*it).px(), (*it).py(), (*it).pz(), (*it).energy() );
		  //~ 
		  //~ if (leptonNumber == 1) {
			  //~ lb = bVec + bJetVector2;
		  //~ }
		  //~ else {
			  //~ lb = aVec + bJetVector2;
		  //~ }
		  //~ if (lb.M() < min_mlb_2){
			  //~ min_mlb_2 = lb.M();
		  //~ }
	  //~ }
  //~ }
  //~ 
  //~ if (bJets35.size() == 1) {
	  //~ 
	  //~ // Get minimal mlB
	  //~ TLorentzVector bJetVector( bJets35.at(0).px(), bJets35.at(0).py(), bJets35.at(0).pz(), bJets35.at(0).energy() );
	  //~ TLorentzVector lb1 = aVec + bJetVector;
	  //~ TLorentzVector lb2 = bVec + bJetVector;
		  //~ if ( lb1.M() < lb2.M()) {
			  //~ min_mlb_1 = lb1.M();
			  //~ leptonNumber = 1;
		  //~ }
		  //~ else {
			  //~ min_mlb_1 = lb2.M();
			  //~ leptonNumber = 2;
		  //~ }
	  //~ 
	  //~ // Get minimal mlj for second lepton 
	  //~ for(std::vector<pat::Jet>::const_iterator it = jets.begin(); it != jets.end() ; ++it){
		  //~ if ((*it).pt() >=35.0 && fabs((*it).eta())<2.4){
			  //~ TLorentzVector jetVector( (*it).px(), (*it).py(), (*it).pz(), (*it).energy() );
			  //~ if (jetVector.DeltaR( bJetVector ) < 0.3) continue;
			  //~ 
			  //~ if (leptonNumber == 1) {
				  //~ lb = bVec + jetVector;
			  //~ }
			  //~ else {
				  //~ lb = aVec + jetVector;
			  //~ }
			  //~ if (lb.M() < min_mlb_2){
				  //~ min_mlb_2 = lb.M();
			  //~ }
		  //~ }
	  //~ }
  //~ }
  //~ 
  //~ if (bJets35.size() == 0) {
	  //~ int jetNumber = 0;
	  //~ int jet1 = 0;
	  //~ 
	  //~ // Loop over jets if no b-jets are present
	  //~ for(std::vector<pat::Jet>::const_iterator it = jets.begin(); it != jets.end() ; ++it){
		  //~ 
		  //~ if ((*it).pt() >=35.0 && fabs((*it).eta())<2.4){
			  //~ jetNumber += 1;
			  //~ TLorentzVector jetVector1( (*it).px(), (*it).py(), (*it).pz(), (*it).energy() );
			  //~ TLorentzVector lb1 = aVec + jetVector1;
			  //~ TLorentzVector lb2 = bVec + jetVector1;
			  //~ if ( lb1.M() < lb2.M()) {
				  //~ if (lb1.M() < min_mlb_1){
					  //~ min_mlb_1 = lb1.M();
					  //~ jet1 = jetNumber;
					  //~ leptonNumber = 1;
				  //~ }
			  //~ }
			  //~ else {
				  //~ if (lb2.M() < min_mlb_1){
					  //~ min_mlb_1 = lb2.M();
					  //~ jet1 = jetNumber;
					  //~ leptonNumber = 2;
				  //~ }
			  //~ }
		  //~ }
	  //~ }
	  //~ 
	  //~ // Get minimal mlj for second lepton
	  //~ jetNumber = 0;
	  //~ for(std::vector<pat::Jet>::const_iterator it = jets.begin(); it != jets.end() ; ++it){
		  //~ if ((*it).pt() >=35.0 && fabs((*it).eta())<2.4){
			  //~ jetNumber += 1;
			  //~ if (jetNumber == jet1) continue;
			  //~ TLorentzVector jetVector2( (*it).px(), (*it).py(), (*it).pz(), (*it).energy() );
			  //~ 
			  //~ if (leptonNumber == 1) {
				  //~ lb = bVec + jetVector2;
			  //~ }
			  //~ else {
				  //~ lb = aVec + jetVector2;
			  //~ }
			  //~ if (lb.M() < min_mlb_2){
				  //~ min_mlb_2 = lb.M();
			  //~ }
		  //~ }
	  //~ }
  //~ }
  //~ if (min_mlb_1 < 1.E6 and min_mlb_2 < 1.E6) *(floatBranches_[treeName]["sumMlb"]) = min_mlb_1 + min_mlb_2; 
  //~ else *(floatBranches_[treeName]["sumMlb"]) = 1.E-6;
  //~ 
  //~ if (min_mlb_1 == 1.E6) min_mlb_1 = 0.0001;
  //~ if (min_mlb_2 == 1.E6) min_mlb_2 = 0.0001;
  //~ 
  //~ *(floatBranches_[treeName]["min_mlb"]) = min_mlb_1; 
  //~ *(floatBranches_[treeName]["max_mlb"]) = min_mlb_2;
  //~ 
  // //~ *(floatBranches_[treeName]["sumMlb"]) = min_mlb_1 + min_mlb_2; 
  //~ 
  float mlb_min = 1.E6;
  float mlb_max = 1.E6;
  
  float temp_mlb;
  
  TLorentzVector jet (0.,0.,0.,0.);
  TLorentzVector jet1 (0.,0.,0.,0.);
  TLorentzVector lepton (0.,0.,0.,0.);
  
  TLorentzVector leptList[2] = {aVec,bVec};
  
  int lmin = -1;
  
  std::vector< pat::Jet > jet1Coll;
  std::vector< pat::Jet > jet2Coll;
  
  if (bJets35.size() >= 1) jet1Coll = bJets35;
  else jet1Coll = jets;
  
  if (bJets35.size() >= 2) jet2Coll = bJets35;
  else jet2Coll = jets;
  
  for (int il=0; il < 2; ++il){
	  for (std::size_t ij=0; ij < jet1Coll.size(); ++ij){
		  jet.SetPxPyPzE(jet1Coll.at(ij).px(),jet1Coll.at(ij).py(),jet1Coll.at(ij).pz(),jet1Coll.at(ij).energy());
		  if (jet.Pt() > 35. && abs(jet.Eta()) < 2.4){
			  lepton.SetPxPyPzE(leptList[il].Px(),leptList[il].Py(),leptList[il].Pz(),leptList[il].Energy());
			  
			  temp_mlb = (jet+lepton).M();
			  if (temp_mlb < mlb_min) {
				  mlb_min = temp_mlb;
				  lmin = il;
				  jet1 = jet;
			  }
		  }
	  }
  }
  for (int il=0; il < 2; ++il){
	  if (il == lmin) continue;
	  for (std::size_t ij=0; ij < jet2Coll.size(); ++ij){
		  jet.SetPxPyPzE(jet2Coll.at(ij).px(),jet2Coll.at(ij).py(),jet2Coll.at(ij).pz(),jet2Coll.at(ij).energy());
		  if (jet.Pt() > 35. && abs(jet.Eta()) < 2.4){
			  if (bJets35.size() == 1 && jet.DeltaR (jet1) < 0.1) continue;
			  if ( (bJets35.size() == 0 || bJets35.size() >= 2) && 	jet == jet1) continue;	
			  lepton.SetPxPyPzE(leptList[il].Px(),leptList[il].Py(),leptList[il].Pz(),leptList[il].Energy());
			  
			  temp_mlb = (jet+lepton).M();
			  if (temp_mlb < mlb_max) {
				  mlb_max = temp_mlb;
			  }
		  }
	  }
  }
	 
  if (mlb_min < 1.E6 and mlb_max < 1.E6) *(floatBranches_[treeName]["sumMlb"]) = mlb_min + mlb_max; 
  else *(floatBranches_[treeName]["sumMlb"]) = 1.E-6;
  
  if (mlb_min == 1.E6) mlb_min = 0.0001;
  if (mlb_max == 1.E6) mlb_max = 0.0001;
  
  *(floatBranches_[treeName]["min_mlb"]) = mlb_min; 
  *(floatBranches_[treeName]["max_mlb"]) = mlb_max; 
   //~ 
  
  //~ fctMT2_.mt2_bisect::mt2 mt2_event;
  
  double pa[3] = {aVec.M(),aVec.Px(),aVec.Py()};
  double pb[3] = {bVec.M(),bVec.Px(),bVec.Py()};
  double pmiss[3] = {0.,met.Px(),met.Py()};
  
  fctMT2_.set_momenta(pa,pb,pmiss);
  
  *(floatBranches_[treeName]["MT2"]) = static_cast<float>(fctMT2_.get_mt2()); 
  
  *(floatBranches_[treeName]["pt3"]) = 0.;
  
  if (*(intBranches_[treeName]["nLooseLeptons"]) > 2) {
	  float leptonPts[*(intBranches_[treeName]["nLooseLeptons"])];
	  float HoverEs[*(intBranches_[treeName]["nLooseLeptons"])];
	  float deltaEtaSuperCluster[*(intBranches_[treeName]["nLooseLeptons"])];
	  float deltaPhiSuperCluster[*(intBranches_[treeName]["nLooseLeptons"])];
	  float eInvMinusPInv[*(intBranches_[treeName]["nLooseLeptons"])];
	  float sigmaIEtaIEta[*(intBranches_[treeName]["nLooseLeptons"])];
	  int nLepton = 0;
	  
	  for(std::vector< pat::Electron >::const_iterator it = looseElectrons.begin(); it != looseElectrons.end() ; ++it){
		  TLorentzVector looseLeptonVector( (*it).px(), (*it).py(), (*it).pz(), (*it).energy() );
		  // Check if the loose lepton is the same as a
		  if ( looseLeptonVector.DeltaR( aVec) > 0.05 &&  looseLeptonVector.DeltaR( bVec) > 0.05) {
			   leptonPts[nLepton] = (*it).pt();
			   if ((*it).ecalEnergy() > 0.) eInvMinusPInv[nLepton]  = 1.0/(*it).ecalEnergy() - (*it).eSuperClusterOverP()/(*it).ecalEnergy();
			   else eInvMinusPInv[nLepton] = -1.;
			   HoverEs[nLepton] =(*it).hadronicOverEm();
			   deltaEtaSuperCluster[nLepton] =(*it).deltaEtaSuperClusterTrackAtVtx();
			   deltaPhiSuperCluster[nLepton] =(*it).deltaPhiSuperClusterTrackAtVtx();
			   sigmaIEtaIEta[nLepton] =(*it).full5x5_sigmaIetaIeta();
		   }
		  else {
			  leptonPts[nLepton] = 0.;
			  eInvMinusPInv[nLepton] = -1.;
			  HoverEs[nLepton] =-1.;
			  deltaEtaSuperCluster[nLepton] =-1.;
			  deltaPhiSuperCluster[nLepton] =-1.;
			  sigmaIEtaIEta[nLepton] =-1.;
		   }
		  nLepton++;
	  }
	  for(std::vector< pat::Muon >::const_iterator it = looseMuons.begin(); it != looseMuons.end() ; ++it){
		  TLorentzVector looseLeptonVector( (*it).px(), (*it).py(), (*it).pz(), (*it).energy() );
		  eInvMinusPInv[nLepton] = -1.;
		  HoverEs[nLepton] =-1.;
		  deltaEtaSuperCluster[nLepton] =-1.;
		  deltaPhiSuperCluster[nLepton] =-1.;
		  sigmaIEtaIEta[nLepton] =-1.;
		  if ( looseLeptonVector.DeltaR( aVec) > 0.05 &&  looseLeptonVector.DeltaR( bVec) > 0.05) leptonPts[nLepton] = (*it).pt();
		  else leptonPts[nLepton] = 0.;
		  nLepton++;
	  }
	  
	  float looseLeptonPt = 0.;
	  int thirdLeptonNumber = -1;
	  for (int i = 0; i < nLepton; ++i){
		  if (leptonPts[i] > looseLeptonPt) thirdLeptonNumber = i;
	  }
		  
	  //~ std::sort(leptonPts, leptonPts + nLepton, std::greater<float>());
	  
	  if (thirdLeptonNumber >= 0) {
		  *(floatBranches_[treeName]["pt3"]) = leptonPts[thirdLeptonNumber];
		  *(floatBranches_[treeName]["HoverE3"]) = HoverEs[thirdLeptonNumber];
		  *(floatBranches_[treeName]["deltaEtaSuperCluster3"]) = deltaEtaSuperCluster[thirdLeptonNumber];
		  *(floatBranches_[treeName]["deltaPhiSuperCluster3"]) = deltaPhiSuperCluster[thirdLeptonNumber];
		  *(floatBranches_[treeName]["eInvMinusPInv3"]) = eInvMinusPInv[thirdLeptonNumber];
		  *(floatBranches_[treeName]["sigmaIEtaIEta3"]) = sigmaIEtaIEta[thirdLeptonNumber];
	  }
	  else {
		  *(floatBranches_[treeName]["pt3"]) = -1.;
		  *(floatBranches_[treeName]["HoverE3"]) = -1.;
		  *(floatBranches_[treeName]["deltaEtaSuperCluster3"]) = -1.;
		  *(floatBranches_[treeName]["deltaPhiSuperCluster3"]) = -1.;
		  *(floatBranches_[treeName]["eInvMinusPInv3"]) = -1.;
		  *(floatBranches_[treeName]["sigmaIEtaIEta3"]) = -1.;
	  }
	 
  }
  
		  	  
  
  //~ std::pair<double, double> pZeta = calcPZeta(aVec, bVec, met);
  *(floatBranches_[treeName]["chargeProduct"]) = a.charge()*b.charge();
  *(tLorentzVectorBranches_[treeName]["p4"]) = comb;
  *(tLorentzVectorBranches_[treeName]["lepton1"]) = aVec;
  *(tLorentzVectorBranches_[treeName]["lepton2"]) = bVec;
  *(tLorentzVectorBranches_[treeName]["vMHT"]) = MHT2;
  *(tLorentzVectorBranches_[treeName]["vMHT20"]) = MHT2Loose;
    *(floatBranches_[treeName]["mht"]) = MHT2.Pt();
  *(tLorentzVectorBranches_[treeName]["vMHTPuppi"]) = MHT2Puppi;
  *(tLorentzVectorBranches_[treeName]["vMHT20Puppi"]) = MHT2LoosePuppi;
    *(floatBranches_[treeName]["mhtPuppi"]) = MHT2Puppi.Pt();
  *(floatBranches_[treeName]["pt1"]) = aVec.Pt();
  *(floatBranches_[treeName]["pt2"]) = bVec.Pt();
  *(floatBranches_[treeName]["charge1"]) = a.charge();
  *(floatBranches_[treeName]["charge2"]) = b.charge();
  *(floatBranches_[treeName]["eta1"]) = aVec.Eta();
  *(floatBranches_[treeName]["eta2"]) = bVec.Eta();
  *(floatBranches_[treeName]["miniIsoEffArea1"]) = getIso(a,"miniIsoEA");
  *(floatBranches_[treeName]["miniIsoEffArea2"]) = getIso(b,"miniIsoEA");
  *(floatBranches_[treeName]["miniIsoPuppi1"]) = getIso(a,"miniIsoPuppi");
  *(floatBranches_[treeName]["miniIsoPuppi2"]) = getIso(b,"miniIsoPuppi");
  *(floatBranches_[treeName]["mt1"]) = transverseMass(aVec, met);
  *(floatBranches_[treeName]["mt2"]) = transverseMass(bVec, met);
  //~ *(floatBranches_[treeName]["fakeWeight1"]) = fakeRates_(a);
  //~ *(floatBranches_[treeName]["fakeWeight2"]) = fakeRates_(b);
  *(floatBranches_[treeName]["deltaPhi"]) = aVec.DeltaPhi( bVec );
  *(floatBranches_[treeName]["deltaR"]) = aVec.DeltaR( bVec );
  *(floatBranches_[treeName]["angle3D"]) = aVec.Angle( bVec.Vect() );
  *(floatBranches_[treeName]["jzb"]) = (met+comb).Pt() - comb.Pt();
  *(floatBranches_[treeName]["jzbPuppi"]) = (metPuppi+comb).Pt() - comb.Pt();
  *(floatBranches_[treeName]["jzbUncorr"]) = (uncorrectedMet+comb).Pt() - comb.Pt();  
  *(floatBranches_[treeName]["jzbResponse"]) = (met+comb).Pt() / comb.Pt();
  *(floatBranches_[treeName]["jzbResponsePuppi"]) = (metPuppi+comb).Pt() / comb.Pt();
  *(floatBranches_[treeName]["jzbResponseUncorr"]) = (uncorrectedMet+comb).Pt() / comb.Pt();    
  
   *(floatBranches_[treeName]["NLL"]) = log(*(floatBranches_[treeName]["sumMlb"])) + log(comb.Pt()) + log(met.Pt()) + log (aVec.DeltaPhi( bVec )) ;
   
  //~ *(floatBranches_[treeName]["pZeta"]) = pZeta.first;
  //~ *(floatBranches_[treeName]["pZetaVis"]) = pZeta.second;

  //~ if (jets->size() >= 2){
    //~ TLorentzVector vJet1 = TLorentzVector(jets->at(0).p4().x(), jets->at(0).p4().y(), jets->at(0).p4().z(), jets->at(0).p4().t());
    //~ TLorentzVector vJet2 = TLorentzVector(jets->at(1).p4().x(), jets->at(1).p4().y(), jets->at(1).p4().z(), jets->at(1).p4().t());
    //~ TLorentzVector sub = aVec + bVec + vJet1 + vJet2;
    //~ *(floatBranches_[treeName]["sqrts"]) = std::sqrt(TMath::Power((std::sqrt(sub.M2() + sub.Perp2()) + met.Et()), 2.0) - (sub.Vect().XYvector() + met.Vect().XYvector())*(sub.Vect().XYvector() + met.Vect().XYvector()));
  //~ } else {
    //~ *(floatBranches_[treeName]["sqrts"]) = -1.0;
  //~ }
  
  fillLeptonIDs(treeName, a , b , vertices );

  //~ if(debug) std::cout << "dB1: "<< *(floatBranches_[treeName]["dB1"]) 
		      //~ << "dB2: "<< *(floatBranches_[treeName]["dB2"])<< std::endl;

  int matched = 0;
  //ETH style genMatching
  *(intBranches_[treeName]["pdgId1"]) = -9999;
  *(intBranches_[treeName]["pdgId2"]) = -9999;  
  *(intBranches_[treeName]["motherPdgId1"]) =-9999;
  *(intBranches_[treeName]["motherPdgId2"]) =-9999;  
  *(intBranches_[treeName]["grandMotherPdgId1"]) = -9999;
  *(intBranches_[treeName]["grandMotherPdgId2"]) = -9999; 
  *(intBranches_[treeName]["isPrompt1"]) = 0;
  *(intBranches_[treeName]["isPrompt2"]) = 0;
  *(intBranches_[treeName]["isFromTau1"]) = 0;
  *(intBranches_[treeName]["isFromTau2"]) = 0;
  *(intBranches_[treeName]["isPromptHardProcess1"]) = 0;
  *(intBranches_[treeName]["isPromptHardProcess2"]) = 0;  
  *(intBranches_[treeName]["isFromTauHardProcess1"]) = 0;
  *(intBranches_[treeName]["isFromTauHardProcess2"]) = 0;          
  std::vector<int> pdgIds1 = getPdgId_.operator()<aT>(a); 
  std::vector<int> pdgIds2 = getPdgId_.operator()<bT>(b); 	
  *(intBranches_[treeName]["pdgId1"]) = pdgIds1[0];
  *(intBranches_[treeName]["pdgId2"]) = pdgIds2[0];  
  *(intBranches_[treeName]["motherPdgId1"]) = pdgIds1[1];
  *(intBranches_[treeName]["motherPdgId2"]) = pdgIds2[1];  
  *(intBranches_[treeName]["grandMotherPdgId1"]) = pdgIds1[2];
  *(intBranches_[treeName]["grandMotherPdgId2"]) = pdgIds2[2]; 
  *(intBranches_[treeName]["isPrompt1"]) = pdgIds1[3];
  *(intBranches_[treeName]["isPrompt2"]) = pdgIds2[3];
  *(intBranches_[treeName]["isFromTau1"]) = pdgIds1[4];
  *(intBranches_[treeName]["isFromTau2"]) = pdgIds2[4];
  *(intBranches_[treeName]["isPromptHardProcess1"]) = pdgIds1[5];
  *(intBranches_[treeName]["isPromptHardProcess2"]) = pdgIds2[5]; 
  *(intBranches_[treeName]["isFromTauHardProcess1"]) = pdgIds1[6];
  *(intBranches_[treeName]["isFromTauHardProcess2"]) = pdgIds2[6];

  TLorentzVector genVec( 0., 0., 0., 0. );
  if(a.genLepton() != NULL){
    matched |= 1;
  }
  if(b.genLepton() != NULL){
    matched |= 2;
  }
  TLorentzVector genLepton1;
  TLorentzVector genLepton2;
  if(a.genLepton() != NULL && b.genLepton() != NULL){
      genLepton1.SetPxPyPzE(a.genLepton()->px(),a.genLepton()->py(),a.genLepton()->pz(),a.genLepton()->energy());
      genLepton2.SetPxPyPzE(b.genLepton()->px(),b.genLepton()->py(),b.genLepton()->pz(),b.genLepton()->energy());

      genVec.SetPxPyPzE(a.genLepton()->px()+b.genLepton()->px(),a.genLepton()->py()+b.genLepton()->py(),a.genLepton()->pz()+b.genLepton()->pz(),a.genLepton()->energy()+b.genLepton()->energy());
  }
  if( matched == 3 && pdgIds1[1] == pdgIds2[1] ) {
    matched |= 4;
  }
  *(floatBranches_[treeName]["genLeptHT"]) += genLepton1.Pt();	
  *(floatBranches_[treeName]["genLeptHT"]) += genLepton2.Pt();	
  
  *(intBranches_[treeName]["matched"]) = matched; 
  *(tLorentzVectorBranches_[treeName]["p4Gen"]) = genVec;
  *(tLorentzVectorBranches_[treeName]["genLepton1"]) = genLepton1;
  *(tLorentzVectorBranches_[treeName]["genLepton2"]) = genLepton2;
  if(debug) std::cout << ", matched = "<<matched<<", motherId = "<<pdgIds1[1];
  if(debug) std::cout<<", M = "<< comb.M() <<", chargeProduct = "<< a.charge()*b.charge() <<std::endl;
  
  
  if (triggerMatches_){

		std::map<string,int> triggerMatches1 = fctTrigger_.operator()<aT>(a);   
 		std::map<string,int> triggerMatches2 = fctTrigger_.operator()<bT>(b);  

		*(intBranches_[treeName]["matchesSingleElectron1"]) = triggerMatches1["matchesSingleElectron"];
		*(intBranches_[treeName]["matchesSingleMuon1"]) = triggerMatches1["matchesSingleMuon"];
		*(intBranches_[treeName]["matchesDoubleElectronLeading1"]) = triggerMatches1["matchesDoubleElectronLeading"];
		*(intBranches_[treeName]["matchesDoubleElectronTrailing1"]) = triggerMatches1["matchesDoubleElectronTrailing"];	
		*(intBranches_[treeName]["matchesDoubleMuonLeading1"]) = triggerMatches1["matchesDoubleMuonLeading"];
		*(intBranches_[treeName]["matchesDoubleMuonTrailing1"]) = triggerMatches1["matchesDoubleMuonTrailing"];
		*(intBranches_[treeName]["matchesDoubleMuonLeadingBoth1"]) = triggerMatches1["matchesDoubleMuonLeadingBoth"];
		*(intBranches_[treeName]["matchesDoubleMuonTrailingBoth1"]) = triggerMatches1["matchesDoubleMuonTrailingBoth"];		
		*(intBranches_[treeName]["matchesDoubleMuonLeadingTk1"]) = triggerMatches1["matchesDoubleMuonLeadingTk"];
		*(intBranches_[treeName]["matchesDoubleMuonTrailingTk1"]) = triggerMatches1["matchesDoubleMuonTrailingTk"];
		*(intBranches_[treeName]["matchesMuELeading1"]) = triggerMatches1["matchesMuELeading"];
		*(intBranches_[treeName]["matchesMuETrailing1"]) = triggerMatches1["matchesMuETrailing"];
		*(intBranches_[treeName]["matchesEMuLeading1"]) = triggerMatches1["matchesMuETrailing"];
		*(intBranches_[treeName]["matchesEMuTrailing1"]) = triggerMatches1["matchesEMuTrailing"];
		
		*(intBranches_[treeName]["matchesDoubleElectronLeadingNonIso1"]) = triggerMatches1["matchesDoubleElectronLeadingNonIso"];
		*(intBranches_[treeName]["matchesDoubleElectronTrailingNonIso1"]) = triggerMatches1["matchesDoubleElectronTrailingNonIso"];			
		*(intBranches_[treeName]["matchesDoubleMuonLeadingNonIso1"]) = triggerMatches1["matchesDoubleMuonLeadingNonIso"];
		*(intBranches_[treeName]["matchesDoubleMuonTrailingNonIso1"]) = triggerMatches1["matchesDoubleMuonTrailingNonIso"];		
		*(intBranches_[treeName]["matchesMuEGElectronNonIso1"]) = triggerMatches1["matchesMuEGElectronNonIso"];
		*(intBranches_[treeName]["matchesMuEGMuonNonIso1"]) = triggerMatches1["matchesMuEGMuonNonIso"];
		
		*(intBranches_[treeName]["matchesSingleElectron2"]) = triggerMatches2["matchesSingleElectron"];
		*(intBranches_[treeName]["matchesSingleMuon2"]) = triggerMatches2["matchesSingleMuon"];
		*(intBranches_[treeName]["matchesDoubleElectronLeading2"]) = triggerMatches2["matchesDoubleElectronLeading"];
		*(intBranches_[treeName]["matchesDoubleElectronTrailing2"]) = triggerMatches2["matchesDoubleElectronTrailing"];		
		*(intBranches_[treeName]["matchesDoubleMuonLeading2"]) = triggerMatches2["matchesDoubleMuonLeading"];
		*(intBranches_[treeName]["matchesDoubleMuonTrailing2"]) = triggerMatches2["matchesDoubleMuonTrailing"];
		*(intBranches_[treeName]["matchesDoubleMuonLeadingBoth2"]) = triggerMatches2["matchesDoubleMuonLeadingBoth"];
		*(intBranches_[treeName]["matchesDoubleMuonTrailingBoth2"]) = triggerMatches2["matchesDoubleMuonTrailingBoth"];		
		*(intBranches_[treeName]["matchesDoubleMuonLeadingTk2"]) = triggerMatches2["matchesDoubleMuonLeadingTk"];
		*(intBranches_[treeName]["matchesDoubleMuonTrailingTk2"]) = triggerMatches2["matchesDoubleMuonTrailingTk"];
		*(intBranches_[treeName]["matchesMuELeading2"]) = triggerMatches2["matchesMuELeading"];
		*(intBranches_[treeName]["matchesMuETrailing2"]) = triggerMatches2["matchesMuETrailing"];
		*(intBranches_[treeName]["matchesEMuLeading2"]) = triggerMatches2["matchesEMuLeading"];
		*(intBranches_[treeName]["matchesEMuTrailing2"]) = triggerMatches2["matchesEMuTrailing"];

 		*(intBranches_[treeName]["matchesDoubleElectronLeadingNonIso2"]) = triggerMatches2["matchesDoubleElectronLeadingNonIso"];
		*(intBranches_[treeName]["matchesDoubleElectronTrailingNonIso2"]) = triggerMatches2["matchesDoubleElectronTrailingNonIso"];			
		*(intBranches_[treeName]["matchesDoubleMuonLeadingNonIso2"]) = triggerMatches2["matchesDoubleMuonLeadingNonIso"];
		*(intBranches_[treeName]["matchesDoubleMuonTrailingNonIso2"]) = triggerMatches2["matchesDoubleMuonTrailingNonIso"];		
		*(intBranches_[treeName]["matchesMuEGElectronNonIso2"]) = triggerMatches2["matchesMuEGElectronNonIso"];
		*(intBranches_[treeName]["matchesMuEGMuonNonIso2"]) = triggerMatches2["matchesMuEGMuonNonIso"]; 
  }

  trees_[treeName]->Fill();
}



//from DQM/Physics/src/EwkTauDQM.cc
//~ std::pair<double, double>
//~ DiLeptonTreesFromMiniAODPuppi::calcPZeta(const TLorentzVector& p1,const TLorentzVector& p2, const TLorentzVector& met)
//~ {
	//~ double cosPhi1 = cos(p1.Phi());
	//~ double sinPhi1 = sin(p1.Phi());
	//~ double cosPhi2 = cos(p2.Phi());
	//~ double sinPhi2 = sin(p2.Phi());
	//~ double zetaX = cosPhi1 + cosPhi2;
	//~ double zetaY = sinPhi1 + sinPhi2;
	//~ double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
	//~ if ( zetaR > 0. ) {
		//~ zetaX /= zetaR;
		//~ zetaY /= zetaR;
	//~ }
//~ 
	//~ double pxVis = p1.Px() + p2.Px();
	//~ double pyVis = p1.Py() + p2.Py();
	//~ double pZetaVis = pxVis*zetaX + pyVis*zetaY;
//~ 
	//~ double px = pxVis + met.Px();
	//~ double py = pyVis + met.Py();
	//~ double pZeta = px*zetaX + py*zetaY;
//~ 
	//~ return std::pair<double, double>(pZeta, pZetaVis);
//~ }

void DiLeptonTreesFromMiniAODPuppi::fillPdfUncert(const edm::Handle< std::vector<double> >& weightHandle, const std::string& pdfIdentifier, const std::string& treeName){
     std::string up = "Up";
     std::string down = "Down";
     bool nnpdfFlag = (pdfIdentifier.substr(0,5)=="NNPDF");
     double centralValue = (*weightHandle)[0];
     *(floatBranches_[treeName][pdfIdentifier]) = float(centralValue);
     if(debug) std::cout << "Cen" << treeName << ": " << centralValue << std::endl;
     unsigned int nmembers = weightHandle->size();
     double wminus = 0.;
     double wplus = 0.;
     unsigned int nplus = 0;
     unsigned int nminus = 0;
     for (unsigned int j=1; j<nmembers; j+=2) {
        float wa = ((*weightHandle)[j]-centralValue)/centralValue;
        float wb = ((*weightHandle)[j+1]-centralValue)/centralValue;
	if (nnpdfFlag) {
	  if (wa>0.) {
	    wplus += wa*wa; 
	    nplus++;
	  } else {
	    wminus += wa*wa;
	    nminus++;
	  }
	  if (wb>0.) {
	    wplus += wb*wb; 
	    nplus++;
	  } else {
	    wminus += wb*wb;
	    nminus++;
	  }
	} else {
	  
	  if (wa>wb) {
            if (wa<0.) wa = 0.;
            if (wb>0.) wb = 0.;
            wplus += wa*wa;
            wminus += wb*wb;
	  } else {
            if (wb<0.) wb = 0.;
            if (wa>0.) wa = 0.;
            wplus += wb*wb;
            wminus += wa*wa;
	  }
	}
    }
    if (wplus>0.) wplus = sqrt(wplus);
    if (wminus>0.) wminus = sqrt(wminus);
    if (nnpdfFlag) {
      if (nplus>0) wplus /= sqrt(nplus);
      if (nminus>0) wminus /= sqrt(nminus);
    }
    *(floatBranches_[treeName][pdfIdentifier+down]) = float(wminus);
    *(floatBranches_[treeName][pdfIdentifier+up]) = float(wplus);
}

const TLorentzVector DiLeptonTreesFromMiniAODPuppi::getMomentum(const  pat::Electron &e)
{
  double corr = 1.;
  double lowEdge = 0.;
  for(std::map<double, double>::iterator it = electronCorrections_.begin(); 
      it != electronCorrections_.end(); ++it){
    if(lowEdge <= fabs(e.superCluster()->eta()) && fabs(e.superCluster()->eta()) < (*it).first ){
      corr = (*it).second;
    }
    lowEdge = (*it).second;
  }

  const TLorentzVector result = TLorentzVector(corr*e.px(), corr*e.py(), corr*e.pz(), corr*e.energy());
  if(debug)std::cout << "correction: "<< corr << ", pt = "<< result.Pt()<<std::endl;
  return result;
}

const TLorentzVector DiLeptonTreesFromMiniAODPuppi::getMomentum(const  pat::Muon &mu)
{
  const TLorentzVector result = TLorentzVector(mu.px(), mu.py(), mu.pz(), mu.energy());
  return result;
}

float DiLeptonTreesFromMiniAODPuppi::getIso(const  pat::Electron &e, const std::string &method)
{
  //  if (e.isEE())
  //  return (e.dr03HcalTowerSumEt() + e.dr03EcalRecHitSumEt() + e.dr03TkSumPt())/e.pt();
  // else
  //  return (e.dr03HcalTowerSumEt() + std::max(0.0, e.dr03EcalRecHitSumEt() - 1.0) + e.dr03TkSumPt())/e.pt();

  //  std::cout<<"electron " << (lepton.trackIso()+lepton.ecalIso()+lepton.hcalIso())/lepton.pt() << std::endl;
  // return (lepton.trackIso()+lepton.ecalIso()+lepton.hcalIso())/lepton.pt() ;
    //  return (e.chargedHadronIso() + e.photonIso() + e.neutralHadronIso()) / e.pt();
  return fctIsolation_(e,method)* 1./e.pt();
}

float DiLeptonTreesFromMiniAODPuppi::getIso(const  pat::Muon &mu, const std::string &method)
{
  //  std::cout<<"muon " << (lepton.trackIso()+lepton.ecalIso()+lepton.hcalIso())/lepton.pt() << std::endl;
  //  return (mu.isolationR03().hadEt + mu.isolationR03().emEt + mu.isolationR03().sumPt) / mu.pt();
  //  return (mu.chargedHadronIso() + mu.photonIso() + mu.neutralHadronIso()) / mu.pt();
  return fctIsolation_(mu,method)* 1./mu.pt();
}


void DiLeptonTreesFromMiniAODPuppi::fillLeptonIDs(const std::string &treeName, const  pat::Electron &ele1, const  pat::Electron &ele2, const edm::Handle<reco::VertexCollection> &vertices)
{

  *(floatBranches_[treeName]["effectiveArea1"]) = getAEffEle(ele1.eta());
  *(floatBranches_[treeName]["chargedIso1"]) = ele1.pfIsolationVariables().sumChargedHadronPt;
  *(floatBranches_[treeName]["neutralIso1"]) = ele1.pfIsolationVariables().sumNeutralHadronEt;
  *(floatBranches_[treeName]["photonIso1"]) = ele1.pfIsolationVariables().sumPhotonEt;
  *(floatBranches_[treeName]["puIso1"]) = ele1.pfIsolationVariables().sumPUPt;

  *(floatBranches_[treeName]["effectiveArea2"]) = getAEffEle(ele2.eta());
  *(floatBranches_[treeName]["chargedIso2"]) = ele2.pfIsolationVariables().sumChargedHadronPt;
  *(floatBranches_[treeName]["neutralIso2"]) = ele2.pfIsolationVariables().sumNeutralHadronEt;
  *(floatBranches_[treeName]["photonIso2"]) = ele2.pfIsolationVariables().sumPhotonEt;
  *(floatBranches_[treeName]["puIso2"]) = ele2.pfIsolationVariables().sumPUPt;
  
  if (writeID_){
  
  *(floatBranches_[treeName]["deltaEtaSuperClusterTrackAtVtx1"]) = ele1.deltaEtaSuperClusterTrackAtVtx();
  *(floatBranches_[treeName]["deltaPhiSuperClusterTrackAtVtx1"]) = ele1.deltaPhiSuperClusterTrackAtVtx();
  *(floatBranches_[treeName]["sigmaIetaIeta1"]) = ele1.sigmaIetaIeta();
  *(floatBranches_[treeName]["hadronicOverEm1"]) = ele1.hadronicOverEm();
  *(floatBranches_[treeName]["eOverP1"]) = abs(1.0/ele1.ecalEnergy() - ele1.eSuperClusterOverP()/ele1.ecalEnergy());
  *(floatBranches_[treeName]["missingHits1"]) = ele1.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
  *(floatBranches_[treeName]["effectiveArea1"]) = getAEffEle(ele1.eta());
  *(floatBranches_[treeName]["passConversion1"]) = ele1.passConversionVeto();
  *(floatBranches_[treeName]["d01"]) = ele1.gsfTrack()->dxy(vertices->at(0).position());
  *(floatBranches_[treeName]["dZ1"]) = fabs(ele1.gsfTrack()->dz(vertices->at(0).position())); 

  *(floatBranches_[treeName]["globalMuon1"]) = -999.;
  *(floatBranches_[treeName]["trackerMuon1"]) = -999.;
  *(floatBranches_[treeName]["pfMuon1"]) = -999.;
  *(floatBranches_[treeName]["trackChi21"]) = -999.;
  *(floatBranches_[treeName]["numberOfValidMuonHits1"]) = -999.;
  *(floatBranches_[treeName]["numberOfMatchedStations1"]) = -999.;
  *(floatBranches_[treeName]["numberOfValidPixelHits1"]) = -999.;
  *(floatBranches_[treeName]["trackerLayersWithMeasurement1"]) = -999.;

  *(floatBranches_[treeName]["deltaEtaSuperClusterTrackAtVtx2"]) = ele2.deltaEtaSuperClusterTrackAtVtx();
  *(floatBranches_[treeName]["deltaPhiSuperClusterTrackAtVtx2"]) = ele2.deltaPhiSuperClusterTrackAtVtx();
  *(floatBranches_[treeName]["sigmaIetaIeta2"]) = ele2.sigmaIetaIeta();
  *(floatBranches_[treeName]["hadronicOverEm2"]) = ele2.hadronicOverEm();
  *(floatBranches_[treeName]["eOverP2"]) = abs(1.0/ele2.ecalEnergy() - ele2.eSuperClusterOverP()/ele2.ecalEnergy());
  *(floatBranches_[treeName]["missingHits2"]) = ele2.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
  *(floatBranches_[treeName]["effectiveArea2"]) = getAEffEle(ele2.eta());
  *(floatBranches_[treeName]["passConversion2"]) = ele2.passConversionVeto();
  *(floatBranches_[treeName]["d02"]) = ele2.gsfTrack()->dxy(vertices->at(0).position());
  *(floatBranches_[treeName]["dZ2"]) = fabs(ele2.gsfTrack()->dz(vertices->at(0).position()));
  
   *(floatBranches_[treeName]["globalMuon2"]) = -999.;
  *(floatBranches_[treeName]["trackerMuon2"]) = -999.;
  *(floatBranches_[treeName]["pfMuon2"]) = -999.;
  *(floatBranches_[treeName]["trackChi22"]) = -999.;
  *(floatBranches_[treeName]["numberOfValidMuonHits2"]) = -999.;
  *(floatBranches_[treeName]["numberOfMatchedStations2"]) = -999.;
  *(floatBranches_[treeName]["numberOfValidPixelHits2"]) = -999.;
  *(floatBranches_[treeName]["trackerLayersWithMeasurement2"]) = -999.; 
  
  }

}


void DiLeptonTreesFromMiniAODPuppi::fillLeptonIDs(const std::string &treeName, const  pat::Electron &ele1, const  pat::Muon &mu2, const edm::Handle<reco::VertexCollection> &vertices)
{
	
  *(floatBranches_[treeName]["effectiveArea1"]) = getAEffEle(ele1.eta());
  *(floatBranches_[treeName]["chargedIso1"]) = ele1.pfIsolationVariables().sumChargedHadronPt;
  *(floatBranches_[treeName]["neutralIso1"]) = ele1.pfIsolationVariables().sumNeutralHadronEt;
  *(floatBranches_[treeName]["photonIso1"]) = ele1.pfIsolationVariables().sumPhotonEt;
  *(floatBranches_[treeName]["puIso1"]) = ele1.pfIsolationVariables().sumPUPt;

  *(floatBranches_[treeName]["effectiveArea2"]) = getAEffMu(mu2.eta());
  *(floatBranches_[treeName]["chargedIso2"]) = mu2.pfIsolationR03().sumChargedHadronPt;
  *(floatBranches_[treeName]["neutralIso2"]) = mu2.pfIsolationR03().sumNeutralHadronEt;
  *(floatBranches_[treeName]["photonIso2"]) = mu2.pfIsolationR03().sumPhotonEt;
  *(floatBranches_[treeName]["puIso2"]) = mu2.pfIsolationR03().sumPUPt;
  
  if (writeID_){
  
  *(floatBranches_[treeName]["deltaEtaSuperClusterTrackAtVtx1"]) = ele1.deltaEtaSuperClusterTrackAtVtx();
  *(floatBranches_[treeName]["deltaPhiSuperClusterTrackAtVtx1"]) = ele1.deltaPhiSuperClusterTrackAtVtx();
  *(floatBranches_[treeName]["sigmaIetaIeta1"]) = ele1.sigmaIetaIeta();
  *(floatBranches_[treeName]["hadronicOverEm1"]) = ele1.hadronicOverEm();
  *(floatBranches_[treeName]["eOverP1"]) = abs(1.0/ele1.ecalEnergy() - ele1.eSuperClusterOverP()/ele1.ecalEnergy());
  *(floatBranches_[treeName]["missingHits1"]) = ele1.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
  *(floatBranches_[treeName]["effectiveArea1"]) = getAEffEle(ele1.eta());
  *(floatBranches_[treeName]["passConversion1"]) = ele1.passConversionVeto();
  *(floatBranches_[treeName]["d01"]) = ele1.gsfTrack()->dxy(vertices->at(0).position());
  *(floatBranches_[treeName]["dZ1"]) = fabs(ele1.gsfTrack()->dz(vertices->at(0).position())); 
 
  *(floatBranches_[treeName]["globalMuon1"]) = -999.;
  *(floatBranches_[treeName]["trackerMuon1"]) = -999.;
  *(floatBranches_[treeName]["pfMuon1"]) = -999.;
  *(floatBranches_[treeName]["trackChi21"]) = -999.;
  *(floatBranches_[treeName]["numberOfValidMuonHits1"]) = -999.;
  *(floatBranches_[treeName]["numberOfMatchedStations1"]) = -999.;
  *(floatBranches_[treeName]["numberOfValidPixelHits1"]) = -999.;
  *(floatBranches_[treeName]["trackerLayersWithMeasurement1"]) = -999.;
  
  *(floatBranches_[treeName]["deltaEtaSuperClusterTrackAtVtx2"]) = -999.;
  *(floatBranches_[treeName]["deltaPhiSuperClusterTrackAtVtx2"]) = -999.;
  *(floatBranches_[treeName]["sigmaIetaIeta2"]) = -999.;
  *(floatBranches_[treeName]["hadronicOverEm2"]) = -999.;
  *(floatBranches_[treeName]["eOverP2"]) = -999.;
  *(floatBranches_[treeName]["missingHits2"]) = -999.;
  *(floatBranches_[treeName]["effectiveArea2"]) = getAEffMu(mu2.eta());
  *(floatBranches_[treeName]["passConversion2"]) = -999.;
  
  *(floatBranches_[treeName]["d02"]) = mu2.dB();
  *(floatBranches_[treeName]["dZ2"]) = abs(mu2.muonBestTrack()->dz(vertices->at(0).position())); 
   if (mu2.isGlobalMuon()){
	  *(floatBranches_[treeName]["trackChi22"]) = mu2.globalTrack()->normalizedChi2();
	  *(floatBranches_[treeName]["numberOfValidMuonHits2"]) = mu2.globalTrack()->hitPattern().numberOfValidMuonHits();
  }
  else{
	 *(floatBranches_[treeName]["trackChi22"]) = -999.;
         *(floatBranches_[treeName]["numberOfValidMuonHits2"]) = -999.;
  }
  *(floatBranches_[treeName]["numberOfMatchedStations2"]) =mu2.numberOfMatchedStations();
  if (mu2.isTrackerMuon()){
  	*(floatBranches_[treeName]["numberOfValidPixelHits2"]) = mu2.innerTrack()->hitPattern().numberOfValidPixelHits();
  	*(floatBranches_[treeName]["trackerLayersWithMeasurement2"]) = mu2.track()->hitPattern().trackerLayersWithMeasurement();
  }
  else{
  	*(floatBranches_[treeName]["numberOfValidPixelHits2"]) = -999.;
  	*(floatBranches_[treeName]["trackerLayersWithMeasurement2"]) = -999.;
  }    
  
  }

}

void DiLeptonTreesFromMiniAODPuppi::fillLeptonIDs(const std::string &treeName, const  pat::Muon &mu1, const  pat::Muon &mu2, const edm::Handle<reco::VertexCollection> &vertices)
{

  *(floatBranches_[treeName]["effectiveArea1"]) = getAEffMu(mu1.eta());
  *(floatBranches_[treeName]["chargedIso1"]) = mu1.pfIsolationR03().sumChargedHadronPt;
  *(floatBranches_[treeName]["neutralIso1"]) = mu1.pfIsolationR03().sumNeutralHadronEt;
  *(floatBranches_[treeName]["photonIso1"]) = mu1.pfIsolationR03().sumPhotonEt;
  *(floatBranches_[treeName]["puIso1"]) = mu1.pfIsolationR03().sumPUPt;

  
  *(floatBranches_[treeName]["effectiveArea2"]) = getAEffMu(mu2.eta());
  *(floatBranches_[treeName]["chargedIso2"]) = mu2.pfIsolationR03().sumChargedHadronPt;
  *(floatBranches_[treeName]["neutralIso2"]) = mu2.pfIsolationR03().sumNeutralHadronEt;
  *(floatBranches_[treeName]["photonIso2"]) = mu2.pfIsolationR03().sumPhotonEt;
  *(floatBranches_[treeName]["puIso2"]) = mu2.pfIsolationR03().sumPUPt;
  
  if (writeID_){
  *(floatBranches_[treeName]["deltaEtaSuperClusterTrackAtVtx1"]) = -999.;
  *(floatBranches_[treeName]["deltaPhiSuperClusterTrackAtVtx1"]) = -999.;
  *(floatBranches_[treeName]["sigmaIetaIeta1"]) = -999.;
  *(floatBranches_[treeName]["hadronicOverEm1"]) = -999.;
  *(floatBranches_[treeName]["eOverP1"]) = -999.;
  *(floatBranches_[treeName]["missingHits1"]) = -999.;
  *(floatBranches_[treeName]["effectiveArea1"]) = getAEffMu(mu1.eta());
  *(floatBranches_[treeName]["passConversion1"]) = -999.;

  *(floatBranches_[treeName]["d01"]) = mu1.dB();
  *(floatBranches_[treeName]["dZ1"]) = abs(mu1.muonBestTrack()->dz(vertices->at(0).position()));


  *(floatBranches_[treeName]["globalMuon1"]) = mu1.isGlobalMuon();
  *(floatBranches_[treeName]["trackerMuon1"]) = mu1.isTrackerMuon();
  *(floatBranches_[treeName]["pfMuon1"]) = mu1.isPFMuon();
  if (mu1.isGlobalMuon()){
	  *(floatBranches_[treeName]["trackChi21"]) = mu1.globalTrack()->normalizedChi2();
	  *(floatBranches_[treeName]["numberOfValidMuonHits1"]) = mu1.globalTrack()->hitPattern().numberOfValidMuonHits();
  }
  else{
	 *(floatBranches_[treeName]["trackChi21"]) = -999.;
         *(floatBranches_[treeName]["numberOfValidMuonHits1"]) = -999.;
  }
  *(floatBranches_[treeName]["numberOfMatchedStations1"]) =mu1.numberOfMatchedStations();
  if (mu1.isTrackerMuon()){
  	*(floatBranches_[treeName]["numberOfValidPixelHits1"]) = mu1.innerTrack()->hitPattern().numberOfValidPixelHits();
  	*(floatBranches_[treeName]["trackerLayersWithMeasurement1"]) = mu1.track()->hitPattern().trackerLayersWithMeasurement();
  }
  else{
  	*(floatBranches_[treeName]["numberOfValidPixelHits1"]) = -999.;
  	*(floatBranches_[treeName]["trackerLayersWithMeasurement1"]) = -999.;
  }

  *(floatBranches_[treeName]["deltaEtaSuperClusterTrackAtVtx2"]) = -999.;
  *(floatBranches_[treeName]["deltaPhiSuperClusterTrackAtVtx2"]) = -999.;
  *(floatBranches_[treeName]["sigmaIetaIeta2"]) = -999.;
  *(floatBranches_[treeName]["hadronicOverEm2"]) = -999.;
  *(floatBranches_[treeName]["eOverP2"]) = -999.;
  *(floatBranches_[treeName]["missingHits2"]) = -999.;
  *(floatBranches_[treeName]["effectiveArea2"]) = getAEffMu(mu2.eta());
  *(floatBranches_[treeName]["passConversion2"]) = -999.;

  *(floatBranches_[treeName]["d02"]) = mu2.dB();
  *(floatBranches_[treeName]["dZ2"]) = abs(mu2.muonBestTrack()->dz(vertices->at(0).position()));



  *(floatBranches_[treeName]["globalMuon2"]) = mu2.isGlobalMuon();
  *(floatBranches_[treeName]["trackerMuon2"]) = mu2.isTrackerMuon();
  *(floatBranches_[treeName]["pfMuon2"]) = mu2.isPFMuon();
  if (mu2.isGlobalMuon()){
	  *(floatBranches_[treeName]["trackChi22"]) = mu2.globalTrack()->normalizedChi2();
	  *(floatBranches_[treeName]["numberOfValidMuonHits2"]) = mu2.globalTrack()->hitPattern().numberOfValidMuonHits();
  }
  else{
	 *(floatBranches_[treeName]["trackChi22"]) = -999.;
         *(floatBranches_[treeName]["numberOfValidMuonHits2"]) = -999.;
  }
  *(floatBranches_[treeName]["numberOfMatchedStations2"]) =mu2.numberOfMatchedStations();
  if (mu2.isTrackerMuon()){
  	*(floatBranches_[treeName]["numberOfValidPixelHits2"]) = mu2.innerTrack()->hitPattern().numberOfValidPixelHits();
  	*(floatBranches_[treeName]["trackerLayersWithMeasurement2"]) = mu2.track()->hitPattern().trackerLayersWithMeasurement();
  }
  else{
  	*(floatBranches_[treeName]["numberOfValidPixelHits2"]) = -999.;
  	*(floatBranches_[treeName]["trackerLayersWithMeasurement2"]) = -999.;
  }
  
  
  }


}


float DiLeptonTreesFromMiniAODPuppi::topPtWeightBen(double topPt){
  if( topPt<0 ) return 1;

  float p0 = 1.18246e+00;
  float p1 = 4.63312e+02;
  float p2 = 2.10061e-06;

  if( topPt>p1 ) topPt = p1;

  float result = p0 + p2 * topPt * ( topPt - 2 * p1 );
  return result;
}

float DiLeptonTreesFromMiniAODPuppi::topPtWeightTOP(double topPt){
  if( topPt<0 ) return 1;

  float p0 = 0.156;
  float p1 = 0.00137;

  float result = exp(p0 - p1 * topPt);
  return result;
}


float DiLeptonTreesFromMiniAODPuppi::getDeltaB(const  pat::Electron &e)
{
  float result = e.dB(pat::Electron::PV3D);
  return result;
}

float DiLeptonTreesFromMiniAODPuppi::getDeltaB(const  pat::Muon &mu)
{
  float result = mu.dB(pat::Muon::PV3D);
  return result;
}

float DiLeptonTreesFromMiniAODPuppi::getAEffEle(double eta)
{
    //from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/src/EGammaCutBasedEleId.cc
    // but gamma + neutral hadrons values...
    double etaAbs = fabs(eta);
    double AEff = 0.1013;
    if (etaAbs > 0.8 && etaAbs <= 1.3) AEff = 0.0988;
    if (etaAbs > 1.3 && etaAbs <= 2.0) AEff = 0.0572;
    if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.0842;
    if (etaAbs > 2.2) AEff = 0.153;
    return AEff;
}

float DiLeptonTreesFromMiniAODPuppi::getAEffMu(double eta)
{
    //from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/src/EGammaCutBasedEleId.cc
    // but gamma + neutral hadrons values...
    double etaAbs = fabs(eta);
    double AEff = 0.0913;
    if (etaAbs > 0.8 && etaAbs <= 1.3) AEff = 0.0765;
    if (etaAbs > 1.3 && etaAbs <= 2.0) AEff = 0.0546;
    if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.0728;
    if (etaAbs > 2.2) AEff = 0.1177;
    return AEff;
}


float DiLeptonTreesFromMiniAODPuppi::transverseMass(const TLorentzVector& p, const TLorentzVector& met)
{
  reco::Candidate::LorentzVector otherMet(met.Px(),met.Py(),met.Pz(),met.E());
  reco::Candidate::LorentzVector leptonT(p.Px(),p.Py(),0.,p.E()*sin(p.Theta()));
  reco::Candidate::LorentzVector sumT=leptonT+otherMet;

  return std::sqrt(sumT.M2());
}

// ------------ method called when starting to processes a luminosity block  ------------

void
DiLeptonTreesFromMiniAODPuppi::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
   newLumiBlock_=true;
}

std::string DiLeptonTreesFromMiniAODPuppi::convertInputTag(const edm::InputTag tag)
{
  std::string result = tag.label();
  if(tag.instance().length() > 0)
    result = tag.instance();
  //  std::cerr << "'"<<tag.label() << "', '"<< tag.instance()<<"' = '"<< result<<"'"<<std::endl;
  return result;
}


// ------------ Method called once each job just before starting event loop  ------------
void 
DiLeptonTreesFromMiniAODPuppi::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiLeptonTreesFromMiniAODPuppi::endJob() 
{
}



//define this as a plug-in
DEFINE_FWK_MODULE(DiLeptonTreesFromMiniAODPuppi);