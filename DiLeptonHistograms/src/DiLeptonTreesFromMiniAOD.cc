// -*- C++ -*-
//
// Package:    Histograms
// Class:      DiLeptonTreesFromMiniAOD
// 
/**\class DiLeptonTreesFromMiniAOD DiLeptonTreesFromMiniAOD.cc brot/DiLeptonTreesFromMiniAOD/src/DiLeptonTreesFromMiniAOD.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  matthias edelhoff
//         Created:  Tue Oct 27 13:50:40 CET 2009
// $Id: DiLeptonTreesFromMiniAOD.cc,v 1.31 2012/09/17 17:38:58 sprenger Exp $
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
//~ #include <SuSyAachen/DiLeptonHistograms/interface/LeptonFullSimScaleFactorMapFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/BTagCalibrationStandalone.h>
#include <SuSyAachen/DiLeptonHistograms/interface/BTagEffMapFunctor.h>


//ROOT
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace std;

//
// class decleration
//

class DiLeptonTreesFromMiniAOD : public edm::EDAnalyzer {
public:
  explicit DiLeptonTreesFromMiniAOD(const edm::ParameterSet&);
  ~DiLeptonTreesFromMiniAOD();

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
  template<class aT, class bT> void fillTree( const std::string &treeName, const aT &a, const bT &b, const std::vector<pat::PackedCandidate>&pfCands,const std::vector<pat::Electron>&looseElectrons,const std::vector<pat::Muon>&looseMuons,const std::vector<pat::Jet>&jets,const std::vector<pat::Jet>&shiftedJetsJESUp,const std::vector<pat::Jet>&shiftedJetsJESDown,const std::vector<pat::Jet>&bJets35,const std::vector<pat::Jet>&shiftedJetsBJESUp,const std::vector<pat::Jet>&shiftedBJetsJESDown,  const pat::MET &patMet, const TLorentzVector &MHT, const edm::Handle<reco::VertexCollection> &vertices, const float &rho, const std::map<std::string, long> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties, const bool &isMC);
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
  float getIso(const  pat::PackedCandidate &track, const std::string &method);
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
  edm::EDGetTokenT< std::vector< pat::Jet > >		 		bJetToken_;
  edm::EDGetTokenT< std::vector< pat::Jet > >		 		bJet35Token_;
  edm::EDGetTokenT< std::vector< pat::MET > > 				metToken_;
  edm::EDGetTokenT<reco::VertexCollection> 					vertexToken_;
  edm::EDGetTokenT< std::vector<pat::PackedCandidate>  > 	pfCandToken_;
  edm::EDGetTokenT< std::vector< reco::GenParticle > >	 	genParticleToken_;
  edm::EDGetTokenT<GenEventInfoProduct>	 					genEventInfoToken_;
  edm::EDGetTokenT<LHEEventProduct>	 						LHEEventToken_;
  edm::EDGetTokenT<double>	 								rhoToken_;
  
  edm::EDGetTokenT<bool>	 								badPFMuonFilterToken_;
  edm::EDGetTokenT<bool>	 								badChargedCandidateFilterToken_;
  edm::EDGetTokenT<bool>	 								cloneGlobalMuonFilterToken_;
  edm::EDGetTokenT<bool>	 								badGlobalMuonFilterToken_;
  
  edm::EDGetTokenT<edm::TriggerResults>	 					metFilterToken_;
  

  std::map<double, double> electronCorrections_;
  //data
  std::map<std::string, TTree*> trees_;  
  std::map<std::string, std::map< std::string, float*> > floatBranches_; 
  std::map<std::string, std::map< std::string, long*> > intBranches_; 
  std::map<std::string, std::map< std::string, bool*> > boolBranches_; 
  std::map<std::string, std::map< std::string, TLorentzVector*> > tLorentzVectorBranches_;

  edm::Handle< std::vector< pat::Jet > > jets;
  
  BTagEffMapFunctor fctBTagEff_;
  
  BTagCalibration fctBTagCalibFullSim_;
  BTagCalibrationReader fctBTagCalibReaderFullSimBJets_;
  BTagCalibrationReader fctBTagCalibReaderFullSimCJets_;
  BTagCalibrationReader fctBTagCalibReaderFullSimLightJets_;
  BTagCalibrationReader fctBTagCalibReaderFullSimBJetsUp_;
  BTagCalibrationReader fctBTagCalibReaderFullSimCJetsUp_;
  BTagCalibrationReader fctBTagCalibReaderFullSimLightJetsUp_;
  BTagCalibrationReader fctBTagCalibReaderFullSimBJetsDown_;
  BTagCalibrationReader fctBTagCalibReaderFullSimCJetsDown_;
  BTagCalibrationReader fctBTagCalibReaderFullSimLightJetsDown_;
  
  //~ LeptonFullSimScaleFactorMapFunctor fctLeptonFullSimScaleFactors_;

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
  bool storeMetFilters_;
  std::vector<std::string> metFilterNames_;
  std::vector<std::string> triggerNames_;
  bool newLumiBlock_;
  std::map<std::string, Bool_t > triggerDecision_;
  std::map<std::string, int > triggerIndex_;
};

// constructors and destructor
DiLeptonTreesFromMiniAOD::DiLeptonTreesFromMiniAOD(const edm::ParameterSet& iConfig):
  electronToken_		(consumes< std::vector< pat::Electron > > 		(iConfig.getParameter<edm::InputTag>("electrons"))),
  looseElectronToken_	(consumes< std::vector< pat::Electron > > 		(iConfig.getParameter<edm::InputTag>("looseElectrons"))),
  muonToken_			(consumes< std::vector< pat::Muon > >			(iConfig.getParameter<edm::InputTag>("muons"))),
  looseMuonToken_		(consumes< std::vector< pat::Muon > >			(iConfig.getParameter<edm::InputTag>("looseMuons"))),
  jetToken_				(consumes< std::vector< pat::Jet > >			(iConfig.getParameter<edm::InputTag>("jets"))),
  genJetToken_			(consumes< std::vector< reco::GenJet  > >		(iConfig.getParameter<edm::InputTag>("genJets"))),
  bJetToken_			(consumes< std::vector< pat::Jet > >			(iConfig.getParameter<edm::InputTag>("bJets"))),
  bJet35Token_			(consumes< std::vector< pat::Jet > >			(iConfig.getParameter<edm::InputTag>("bJets35"))),
  metToken_				(consumes< std::vector< pat::MET > >			(iConfig.getParameter<edm::InputTag>("met"))),
  vertexToken_			(consumes<reco::VertexCollection>				(iConfig.getParameter<edm::InputTag>("vertices"))),
  pfCandToken_			(consumes< std::vector<pat::PackedCandidate>  >	(iConfig.getParameter<edm::InputTag>("pfCands"))),
  genParticleToken_		(consumes< std::vector< reco::GenParticle > >	(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genEventInfoToken_	(consumes<GenEventInfoProduct>					(iConfig.getParameter<edm::InputTag>("pdfInfo"))),
  LHEEventToken_		(consumes<LHEEventProduct>						(iConfig.getParameter<edm::InputTag>("LHEInfo"))),
  rhoToken_				(consumes<double>								(iConfig.getParameter<edm::InputTag>("rho"))),
  
  badPFMuonFilterToken_				(consumes<bool>						(iConfig.getParameter<edm::InputTag>("badPFMuonFilter"))),
  badChargedCandidateFilterToken_	(consumes<bool>						(iConfig.getParameter<edm::InputTag>("badChargedCandidateFilter"))),
  cloneGlobalMuonFilterToken_		(consumes<bool>						(iConfig.getParameter<edm::InputTag>("cloneGlobalMuonFilter"))),
  badGlobalMuonFilterToken_			(consumes<bool>						(iConfig.getParameter<edm::InputTag>("badGlobalMuonFilter"))),
  metFilterToken_					(consumes<edm::TriggerResults>		(edm::InputTag("TriggerResults",""))),
	
  fctBTagEff_    (iConfig.getParameter<edm::ParameterSet>("bTagEfficiencies") ),
  
  fctBTagCalibFullSim_    (iConfig.getParameter<edm::ParameterSet>("BTagCalibration").getParameter<std::string>("CSVFullSimTagger"),iConfig.getParameter<edm::ParameterSet>("BTagCalibration").getParameter<std::string>("CSVFullSimFileName") ),
  fctBTagCalibReaderFullSimBJets_    (&fctBTagCalibFullSim_,BTagEntry::OP_MEDIUM,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_bJets"),"central" ),
  fctBTagCalibReaderFullSimCJets_    (&fctBTagCalibFullSim_,BTagEntry::OP_MEDIUM,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_cJets"),"central" ),
  fctBTagCalibReaderFullSimLightJets_    (&fctBTagCalibFullSim_,BTagEntry::OP_MEDIUM,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_lightJets"),"central" ),
  fctBTagCalibReaderFullSimBJetsUp_    (&fctBTagCalibFullSim_,BTagEntry::OP_MEDIUM,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_bJets"),"up" ),
  fctBTagCalibReaderFullSimCJetsUp_    (&fctBTagCalibFullSim_,BTagEntry::OP_MEDIUM,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_cJets"),"up" ),
  fctBTagCalibReaderFullSimLightJetsUp_    (&fctBTagCalibFullSim_,BTagEntry::OP_MEDIUM,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_lightJets"),"up" ),
  fctBTagCalibReaderFullSimBJetsDown_    (&fctBTagCalibFullSim_,BTagEntry::OP_MEDIUM,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_bJets"),"down" ),
  fctBTagCalibReaderFullSimCJetsDown_    (&fctBTagCalibFullSim_,BTagEntry::OP_MEDIUM,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_cJets"),"down" ),
  fctBTagCalibReaderFullSimLightJetsDown_    (&fctBTagCalibFullSim_,BTagEntry::OP_MEDIUM,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_lightJets"),"down" ),
  
  //~ fctLeptonFullSimScaleFactors_ (iConfig.getParameter<edm::ParameterSet>("LeptonFullSimScaleFactors") ),
  
  fctVtxWeight_   	(iConfig.getParameter<edm::ParameterSet>("vertexWeights") ,consumesCollector()),
  fctVtxWeightUp_   (iConfig.getParameter<edm::ParameterSet>("vertexWeightsUp") ,consumesCollector()),
  fctVtxWeightDown_ (iConfig.getParameter<edm::ParameterSet>("vertexWeightsDown"),consumesCollector() ),
  fctIsolation_  	(iConfig.getParameter<edm::ParameterSet>("isolationDefinitions"),consumesCollector()),
  fctTrigger_  		(iConfig.getParameter<edm::ParameterSet>("triggerDefinitions"),consumesCollector()),  
  getPdgId_			(iConfig.getParameter< edm::ParameterSet>("pdgIdDefinition"),consumesCollector() ),
  metFilterNames_	(iConfig.getUntrackedParameter< std::vector <std::string> >("metFilterNames")),
  triggerNames_		(iConfig.getUntrackedParameter< std::vector <std::string> >("triggerNames")),
  newLumiBlock_(true)
{
  debug = false; 
  //~ writeID_ = iConfig.existsAs<edm::InputTag>("baseTrees");
  writeID_ = iConfig.getUntrackedParameter<bool>("writeID");
  writeTrigger_ = iConfig.getUntrackedParameter<bool>("writeTrigger");
  storeMetFilters_ = iConfig.getUntrackedParameter<bool>("storeMetFilters");
  triggerMatches_ = iConfig.existsAs<edm::InputTag>("triggerSummaryTag");  
  metUncert_ = iConfig.getUntrackedParameter<bool>("doMETUncert");  
  
  consumes<std::vector< PileupSummaryInfo > >(edm::InputTag("slimmedAddPileupInfo"));
  
  consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
  consumes<edm::TriggerResults>(edm::InputTag("TriggerResults",""));
  
  //~ consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"));
  //~ consumes<std::vector<pat::TriggerObjectStandAlone>>(edm::InputTag("selectedPatTrigger"));
	  

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
  initFloatBranch( "bTagWeight" );
  initFloatBranch( "bTagWeightErrHeavy" );
  initFloatBranch( "bTagWeightErrLight" );
  //~ initFloatBranch( "leptonFullSimScaleFactor1" );
  //~ initFloatBranch( "leptonFullSimScaleFactor2" );
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
  initTLorentzVectorBranch( "genJet1" );
  initTLorentzVectorBranch( "genJet2" );
  initTLorentzVectorBranch( "bJet1" );
  initTLorentzVectorBranch( "bJet2" );
  initTLorentzVectorBranch( "vMet" );   
  //~ initTLorentzVectorBranch( "vMetUncorrected" ); 
  initTLorentzVectorBranch( "vGenMet" );  
  initFloatBranch( "rho" );
  initFloatBranch( "pt1" );
  initFloatBranch( "pt2" );
  initFloatBranch( "ptErr1" );
  initFloatBranch( "ptErr2" );
  initFloatBranch( "pt3" );
  initFloatBranch( "ptTrack3" );
  initFloatBranch( "eta1" );
  initFloatBranch( "eta2" );
  initFloatBranch( "miniIsoEffArea1" );
  initFloatBranch( "miniIsoEffArea2" );

  initFloatBranch( "mt1" );
  initFloatBranch( "mt2" );
  //~ initFloatBranch( "fakeWeight1" );
  //~ initFloatBranch( "fakeWeight2" );
  initFloatBranch( "deltaPhi" );
  initFloatBranch( "deltaR" );
  //~ initFloatBranch( "angle3D" );
  initFloatBranch( "min_mlb" );
  initFloatBranch( "max_mlb" );
  initFloatBranch( "sumMlb" );
  initFloatBranch( "sumMlbJESUp" );
  initFloatBranch( "sumMlbJESDown" );
  initFloatBranch( "MT2" );
  initFloatBranch( "deltaPhiJetMet1" );
  initFloatBranch( "deltaPhiJetMet2" );
  initFloatBranch( "jzbResponse" );
  initFloatBranch( "jzbResponseUncorr" ); 
  initFloatBranch( "jzb" );
  initFloatBranch( "jzbUncorr" );  
  initFloatBranch( "ht" );
  initFloatBranch( "htJESUp" );
  initFloatBranch( "htJESDown" );
  //~ initFloatBranch( "ht40" );  
  //~ initFloatBranch( "leptHT" );  
  //~ initFloatBranch( "genLeptHT" );  
  initFloatBranch( "mht" );
  initFloatBranch( "met" );
  initFloatBranch( "caloMet" );
  initFloatBranch( "genMet" );
  initFloatBranch( "uncorrectedMet" ); 
  initFloatBranch( "metJESUp" );
  initFloatBranch( "metJESDown" );
  //~ initFloatBranch( "pZeta" );
  //~ initFloatBranch( "pZetaVis" );
  initIntBranch( "nJets" );
  initIntBranch( "nGenJets" );
  initIntBranch( "nBadMuonJets" );
  //~ initIntBranch( "nJets30" );  
  //initIntBranch( "nJetsNoPULoose" );
  //initIntBranch( "nJetsNoPUMedium" );
  //initIntBranch( "nJetsNoPUTight" );     
  //~ initIntBranch( "nJetsOld" );
  initIntBranch( "nBJets" );
  initIntBranch( "nBJets35" );
  initIntBranch( "nShiftedJetsJESUp" );
  initIntBranch( "nShiftedJetsJESDown" );
  initIntBranch( "nVertices" );
  initIntBranch( "nGenVertices" );
  initIntBranch( "nLightLeptons" );
  initIntBranch( "nLooseLeptons" );
  initIntBranch( "nIsoTracks" );
  initFloatBranch( "jet1pt" );
  initFloatBranch( "jet2pt" );
  initFloatBranch( "jet3pt" );
  initFloatBranch( "jet4pt" );
  initFloatBranch( "genJet1pt" );
  initFloatBranch( "genJet2pt" );
  initFloatBranch( "genJet3pt" );
  initFloatBranch( "genJet4pt" );
  initFloatBranch( "bjet1pt" );
  initFloatBranch( "bjet2pt" );
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
  
  if (storeMetFilters_){
  
	  initIntBranch( "badPFMuonFilter" );
	  initIntBranch( "badChargedCandidateFilter" );
	  initIntBranch( "cloneGlobalMuonFilter" );
	  initIntBranch( "badGlobalMuonFilter" );
	  
	  for (const auto& n : metFilterNames_){
		  initIntBranch( n.c_str() );
	  }
	  
  }
  
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
	  
	  initFloatBranch( "chargedIso1");
	  initFloatBranch( "neutralIso1");
	  initFloatBranch( "photonIso1");
	  initFloatBranch( "puIso1");
	  initFloatBranch( "effectiveArea1");
	  
	  initFloatBranch( "chargedIso2");
	  initFloatBranch( "neutralIso2");
	  initFloatBranch( "photonIso2");
	  initFloatBranch( "puIso2");
	  initFloatBranch( "effectiveArea2");
 
 
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
  
}

void 
DiLeptonTreesFromMiniAOD::initTLorentzVectorBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    tLorentzVectorBranches_[(*it).first][name] = new TLorentzVector;
    (*it).second->Branch(name.c_str(), "TLorentzVector" ,&tLorentzVectorBranches_[(*it).first][name]);
  }
}

void 
DiLeptonTreesFromMiniAOD::initFloatBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    floatBranches_[(*it).first][name] = new float;
    (*it).second->Branch(name.c_str(), floatBranches_[(*it).first][name], (name+"/F").c_str());
  }
}

void 
DiLeptonTreesFromMiniAOD::initIntBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    intBranches_[(*it).first][name] = new long;
    (*it).second->Branch(name.c_str(), intBranches_[(*it).first][name], (name+"/L").c_str());
  }
}

DiLeptonTreesFromMiniAOD::~DiLeptonTreesFromMiniAOD()
{ 
  for( std::map<std::string, std::map< std::string, float*> >::const_iterator it = floatBranches_.begin();
       it != floatBranches_.end(); ++it){
    for( std::map< std::string, float*>::const_iterator it2 = (*it).second.begin();
	 it2 != (*it).second.end(); ++it2){
      if(debug)std::cout << "deleting: " << (*it).first << " - "<< (*it2).first << std::endl;
      delete (*it2).second;
    }
  }
  for( std::map<std::string, std::map< std::string, long*> >::const_iterator it = intBranches_.begin();
       it != intBranches_.end(); ++it){
    for( std::map< std::string, long*>::const_iterator it2 = (*it).second.begin();
	 it2 != (*it).second.end(); ++it2){
      if(debug) std::cout << "deleting: " << (*it).first << " - "<< (*it2).first << std::endl;
      delete (*it2).second;
    }
  }


}


// member functions
// ------------ method called to for each event  ------------
void
DiLeptonTreesFromMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  iEvent.getByToken(pfCandToken_, pfCands); 

  
  iEvent.getByToken(jetToken_, jets);

  edm::Handle< std::vector< reco::GenJet  > > genJets;
  iEvent.getByToken(genJetToken_, genJets);

  edm::Handle< std::vector< pat::Jet > > bJets;
  iEvent.getByToken(bJetToken_, bJets);
  
  edm::Handle< std::vector< pat::Jet > > bJets35;
  iEvent.getByToken(bJet35Token_, bJets35);
  
  edm::Handle< std::vector< pat::MET > > mets;
  iEvent.getByToken(metToken_, mets);


  edm::Handle< std::vector< reco::GenParticle > > genParticles;
  iEvent.getByToken(genParticleToken_, genParticles);
  
  
  bool isMC = genParticles.isValid();
  
  if (writeTrigger_){ 
	  edm::Handle<edm::TriggerResults> triggerBits;
	  edm::InputTag triggerTag("TriggerResults","","HLT");
	  iEvent.getByLabel(triggerTag, triggerBits);
	  
	  // for each lumiBlock, re-read the trigger indices (rather changes for new run)
	   if( triggerIndex_.size() && newLumiBlock_ ) {
	      newLumiBlock_=false;
	      // set all trigger indices to -1 as "not available"-flag
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
  }
   
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);

  getPdgId_.loadGenParticles(iEvent);
  fctIsolation_.init(iEvent);
  fctTrigger_.loadTrigger(iEvent);
  std::map<std::string, long> intEventProperties;
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
        for (unsigned int i = 0; i < statusCodes.size(); i++) {
            if (statusCodes[i] == 1) {
                if (abs(lheParticleInfo.IDUP[i]) < 11 || abs(lheParticleInfo.IDUP[i]) > 16) {
                    ht += sqrt(pow(allParticles[i][0], 2) + pow(allParticles[i][1], 2));
                }
            }
        }
        floatEventProperties["genParticleHT"] = ht;
    }
	
	else {
		floatEventProperties["genParticleHT"] = -999.;
	}      	

	
  intEventProperties["nVertices"] = vertices->size();
  
  int nGenVertices = -1;
 
  
  if (isMC){	  
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
  
  edm::Handle<bool> iFilterBadPFMuon;
  iEvent.getByToken(badPFMuonFilterToken_, iFilterBadPFMuon);
  
  edm::Handle<bool> iFilterBadChargedCandidate;
  iEvent.getByToken(badChargedCandidateFilterToken_, iFilterBadChargedCandidate);
  
  edm::Handle<bool> iFilterCloneGlobalMuon;
  iEvent.getByToken(cloneGlobalMuonFilterToken_, iFilterCloneGlobalMuon);
  
  edm::Handle<bool> iFilterBadGlobalMuon;
  iEvent.getByToken(badGlobalMuonFilterToken_, iFilterBadGlobalMuon);
  
  edm::Handle<edm::TriggerResults> metFilterBits;
  iEvent.getByToken(metFilterToken_, metFilterBits);
  const edm::TriggerNames &allFilterNames = iEvent.triggerNames(*metFilterBits);
  
  
  if (storeMetFilters_){ 
	    
	  if (*iFilterBadPFMuon) intEventProperties["badPFMuonFilter"] = 1;
	  else intEventProperties["badPFMuonFilter"] = 0;  
	  
	  if (*iFilterBadChargedCandidate) intEventProperties["badChargedCandidateFilter"] = 1;
	  else intEventProperties["badChargedCandidateFilter"] = 0;  
	  
	  if (*iFilterCloneGlobalMuon) intEventProperties["cloneGlobalMuonFilter"] = 1;
	  else intEventProperties["cloneGlobalMuonFilter"] = 0;  
	  
	  if (*iFilterBadGlobalMuon) intEventProperties["badGlobalMuonFilter"] = 1;
	  else intEventProperties["badGlobalMuonFilter"] = 0;
	  
	  
	  for (std::string const &filterName : metFilterNames_){
		  const unsigned index = allFilterNames.triggerIndex(filterName);
		  if (index >= allFilterNames.size()) std::cerr << "MET filter " << filterName << "not found!" << std::endl;
		  if (metFilterBits->accept(index)) intEventProperties[filterName] = 1;
		  else intEventProperties[filterName] = 0;

	  }  
  }
 else{ 
	    
	  if (!*iFilterBadPFMuon) return;
	  if (!*iFilterBadChargedCandidate) return;
	  if (!*iFilterCloneGlobalMuon) return;
	  if (!*iFilterBadGlobalMuon) return;	  
	  
	  for (std::string const &filterName : metFilterNames_){
		  const unsigned index = allFilterNames.triggerIndex(filterName);
		  if (index >= allFilterNames.size()) std::cerr << "MET filter " << filterName << "not found!" << std::endl;
		  if (!metFilterBits->accept(index)) return;
	  }  
  }


  intEventProperties["nBJets"] = bJets->size();
  intEventProperties["nBJets35"] = bJets35->size();
  intEventProperties["nLightLeptons"] = electrons->size() + muons->size();
  intEventProperties["nLooseLeptons"] = looseElectrons->size() + looseMuons->size();
  intEventProperties["runNr"] = iEvent.id().run();
  intEventProperties["lumiSec"] = iEvent.id().luminosityBlock();
  intEventProperties["eventNr"] = iEvent.id().event();
	  



  pat::MET met = mets->front();
  TLorentzVector metVector(met.px(), met.py(), met.pz(), met.energy());
  TLorentzVector uncorrectedMetVector;
  uncorrectedMetVector.SetPtEtaPhiE(met.uncorPt(), 0,	met.uncorPhi(), met.uncorPt());
  
  
  floatEventProperties["met"] = metVector.Pt();
  tLorentzVectorEventProperties["vMet"] = metVector; 
  
  //~ tLorentzVectorEventProperties["vMetUncorrected"] = uncorrectedMetVector;
  floatEventProperties["uncorrectedMet"] = uncorrectedMetVector.Pt();
  
  floatEventProperties["caloMet"] = met.caloMETPt();
  
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
  if (isMC){
	  	
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
  TLorentzVector genMHT;  
  TLorentzVector tempMHT;
  
  floatEventProperties["ht"] = 0.0;
  floatEventProperties["htJESUp"] = 0.0;
  floatEventProperties["htJESDown"] = 0.0;
  floatEventProperties["mht"] = 0.0;

  //~ edm::FileInPath fip("SuSyAachen/DiLeptonHistograms/data/GR_P_V40_AN1::All_Uncertainty_AK5PF.txt");
  //~ JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(fip.fullPath());

  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl); 
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);
  //~ JetCorrectionUncertainty jecUnc(JetCorPar);


  TLorentzVector JESChangeUp(0,0,0,0);
  TLorentzVector JESChangeDown(0,0,0,0);
  TLorentzVector tempJESChangeUp;
  TLorentzVector tempJESChangeDown;
  // loop over jets
  std::vector<pat::Jet> * shiftedJetsJESUp = new std::vector<pat::Jet>(); 
  std::vector<pat::Jet> * shiftedJetsJESDown = new std::vector<pat::Jet>(); 
  std::vector<pat::Jet> * shiftedBJetsJESUp = new std::vector<pat::Jet>(); 
  std::vector<pat::Jet> * shiftedBJetsJESDown = new std::vector<pat::Jet>(); 


  // Shift jets by JEC
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
    tempJESChangeUp.SetPxPyPzE((ajetUp.px()-(*itJet).px()), (ajetUp.py()-(*itJet).py()), (ajetUp.pz()-(*itJet).pz()), (ajetUp.energy()-(*itJet).energy()));	
    tempJESChangeDown.SetPxPyPzE((ajetDown.px()-itJet->px()), (ajetDown.py()-itJet->py()), (ajetDown.pz()-itJet->pz()), (ajetDown.energy()-itJet->energy()));	
    JESChangeUp +=tempJESChangeUp;
    JESChangeDown += tempJESChangeDown;

    shiftedJetsJESUp->push_back(ajetUp);
    shiftedJetsJESDown->push_back(ajetDown);
  }

  // Shift b-jets by JEC
  for (std::vector<pat::Jet>::const_iterator itBJet = bJets35->begin(); itBJet != bJets35->end(); itBJet++) {
	
    pat::Jet ajetUp( *itBJet );
    pat::Jet ajetDown( *itBJet );
    jecUnc->setJetEta(itBJet->eta());
    jecUnc->setJetPt(itBJet->pt()); // IMPORTANT: the uncertainty is a function of the CORRECTED pt
    double correction = jecUnc->getUncertainty(true);
    //double correction = itJet->relCorrUncert(dir_); 
    
    //use pat::Jet::relCorrUncert
    if (correction==0.0){
	correction = 0.045; 
    }
    //    std::cout << "jet pt = "<<itJet->p4()<<", uncertainty = "<<correction<<std::endl;
    ajetUp.setP4( itBJet->p4() * (1.0 + correction) );
    ajetDown.setP4( itBJet->p4() * (1.0 - correction) );

    shiftedBJetsJESUp->push_back(ajetUp);
    shiftedBJetsJESDown->push_back(ajetDown);
  }
  
  int nJets=0;
  int nBadMuonJets=0;
  int nGenJets=0;

  
  for(std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end() ; ++it){
	if ((*it).pt() >=35.0 && fabs((*it).eta())<2.4){
		nJets++;
		

		tempMHT.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());	
		MHT = MHT + tempMHT;
		floatEventProperties["ht"] += (*it).pt();
		
		if ((*it).pt() >=200.0 && (*it).muonEnergyFraction() > 0.5 && fabs(tempMHT.DeltaPhi( metVector )) > M_PI - 0.4 ){
			nBadMuonJets++;
		}
		
	}
	
	
  }
  intEventProperties["nJets"] = nJets;
  intEventProperties["nBadMuonJets"] = nBadMuonJets;
    
  float genHT = 0.;
  
  if (isMC){
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
			}
			
		}	
		
	  }
	  
  }
  intEventProperties["nGenJets"] = nGenJets;
  floatEventProperties["genHT"] = genHT;
 


  TLorentzVector jet1Vector(0.,0.,0.,0.);
  TLorentzVector jet2Vector(0.,0.,0.,0.); 
  
  floatEventProperties["deltaPhiJetMet1"] = -9999.;
  floatEventProperties["deltaPhiJetMet2"] = -9999.;


  if (nJets > 0)
    jet1Vector.SetPxPyPzE(jets->at(0).px(),jets->at(0).py(),jets->at(0).pz(),jets->at(0).energy());
    floatEventProperties["deltaPhiJetMet1"] = fabs(jet1Vector.DeltaPhi( metVector ));
  if (nJets > 1)
    jet2Vector.SetPxPyPzE(jets->at(1).px(),jets->at(1).py(),jets->at(1).pz(),jets->at(1).energy());
    floatEventProperties["deltaPhiJetMet2"] = fabs(jet2Vector.DeltaPhi( metVector ));
  
  tLorentzVectorEventProperties["jet1"] = jet1Vector;
  tLorentzVectorEventProperties["jet2"] = jet2Vector;
  

  TLorentzVector genJet1Vector(0.,0.,0.,0.);
  TLorentzVector genJet2Vector(0.,0.,0.,0.);

  if (nGenJets > 0)
    genJet1Vector.SetPxPyPzE(genJets->at(0).px(),genJets->at(0).py(),genJets->at(0).pz(),genJets->at(0).energy());
  if (nGenJets > 1)
    genJet2Vector.SetPxPyPzE(genJets->at(1).px(),genJets->at(1).py(),genJets->at(1).pz(),genJets->at(1).energy());
  tLorentzVectorEventProperties["genJet1"] = genJet1Vector;
  tLorentzVectorEventProperties["genJet2"] = genJet2Vector;


  int nJetsJESUp=0;
  for(std::vector<pat::Jet>::const_iterator it = shiftedJetsJESUp->begin(); it != shiftedJetsJESUp->end() ; ++it){
	if ((*it).pt() >=35.0 && fabs((*it).eta())<2.4){
        	nJetsJESUp++;
		floatEventProperties["htJESUp"] += (*it).pt();
	}
  }
  intEventProperties["nShiftedJetsJESUp"] = nJetsJESUp;
  int nJetsJESDown=0;
  for(std::vector<pat::Jet>::const_iterator it = shiftedJetsJESDown->begin(); it != shiftedJetsJESDown->end() ; ++it){
	if ((*it).pt() >=35.0 && fabs((*it).eta())<2.4){	
        	nJetsJESDown++;
		floatEventProperties["htJESDown"] += (*it).pt();
	}
  }
  intEventProperties["nShiftedJetsJESDown"] = nJetsJESDown;
  TLorentzVector metVectorJESUp = metVector + JESChangeUp;
  TLorentzVector metVectorJESDown = metVector + JESChangeDown;
  floatEventProperties["metJESUp"] = metVectorJESUp.Pt();


  floatEventProperties["metJESDown"] = metVectorJESDown.Pt();



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
    
  
  if (isMC){
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
				if (genJet == 1) floatEventProperties["genJet1pt"] =(*it).pt();
				if (genJet == 2) floatEventProperties["genJet2pt"] =(*it).pt();
				if (genJet == 3) floatEventProperties["genJet3pt"] =(*it).pt();
				if (genJet == 4) floatEventProperties["genJet4pt"] =(*it).pt();
				
			}	
		}	
		
	  }
	  
  }
  
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
  
   floatEventProperties["bTagWeight"] = -9999.;
  floatEventProperties["bTagWeightErrHeavy"] = -9999.;
  floatEventProperties["bTagWeightErrLight"] = -9999.;
  
  // Calculate btagging weights
  if (genParticles.isValid()){
	  
	  float P_MC = 1.;
	  float P_Data = 1.;
	  float err1 = 0.;
	  float err2 = 0.;
	  float err3 = 0.;
	  float err4 = 0.;
	  
	  for(std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end() ; ++it){
		  if ((*it).pt() >=35.0 && fabs((*it).eta())<2.4){
			  
			  int jetFlavor;
			  float eff;
			  float SF, SF_up, SF_down;
			  float SF_err;
			  
			  jetFlavor = (*it).hadronFlavour();
			  if (jetFlavor == 5)
			  {
				  SF = fctBTagCalibReaderFullSimBJets_.eval(BTagEntry::FLAV_B, (*it).eta(), (*it).pt());
				  SF_up = fctBTagCalibReaderFullSimBJetsUp_.eval(BTagEntry::FLAV_B, (*it).eta(), (*it).pt());
				  SF_down = fctBTagCalibReaderFullSimBJetsDown_.eval(BTagEntry::FLAV_B, (*it).eta(), (*it).pt());
				  
				 
			  }
			  else if (jetFlavor == 4)
			  {
				  SF = fctBTagCalibReaderFullSimCJets_.eval(BTagEntry::FLAV_C, (*it).eta(), (*it).pt());
				  SF_up = fctBTagCalibReaderFullSimCJetsUp_.eval(BTagEntry::FLAV_C, (*it).eta(), (*it).pt());
				  SF_down = fctBTagCalibReaderFullSimCJetsDown_.eval(BTagEntry::FLAV_C, (*it).eta(), (*it).pt());

			  }
			  else
			  {
				  SF = fctBTagCalibReaderFullSimLightJets_.eval(BTagEntry::FLAV_UDSG, (*it).eta(), (*it).pt());
				  SF_up = fctBTagCalibReaderFullSimLightJetsUp_.eval(BTagEntry::FLAV_UDSG, (*it).eta(), (*it).pt());
				  SF_down = fctBTagCalibReaderFullSimLightJetsDown_.eval(BTagEntry::FLAV_UDSG, (*it).eta(), (*it).pt());

			  }
			  
			 			  
			  eff = fctBTagEff_(jetFlavor, (*it).pt(), fabs((*it).eta()));
			  
			  SF_err = std::max(fabs(SF_up-SF),fabs(SF_down-SF));
				  
			  // check if jet is btagged
			  bool tagged = false;
			  if (bJets->size() > 0){
				  for (unsigned int i = 0; i < bJets->size(); ++i){
					  if ((*it).pt()==bJets->at(i).pt()){
						  tagged = true;
						  break;
					  }
				  }
			  }
			  
			  if (tagged){
				  P_MC = P_MC * eff;
				  P_Data = P_Data * SF * eff;
				  if (jetFlavor == 5 || jetFlavor == 4) err1 += SF_err/SF;
				  else err3 += SF_err/SF;
			  }
			  else{
				  P_MC = P_MC * (1 - eff);
				  P_Data = P_Data  * (1 - SF * eff);
				  if (jetFlavor == 5 || jetFlavor == 4) err2 += (-eff*SF_err) / (1-eff*SF);
				  else err4 += (-eff*SF_err) / (1-eff*SF);
			  }
		  }
	  }
	  if (P_MC > 0. && P_Data > 0.){
		  floatEventProperties["bTagWeight"] = P_Data/P_MC;
		  floatEventProperties["bTagWeightErrHeavy"] = (err1+err2)*P_Data/P_MC;
		  floatEventProperties["bTagWeightErrLight"] = (err3+err4)*P_Data/P_MC;
	  }
	  else{
		  floatEventProperties["bTagWeight"] = 0.;
		  floatEventProperties["bTagWeightErrHeavy"] = 0.;
		  floatEventProperties["bTagWeightErrLight"] = 0.;
	  }
  }
  
  
  
  	
  floatEventProperties["weight"] = fctVtxWeight_( iEvent );
  floatEventProperties["weightUp"] = fctVtxWeightUp_( iEvent );
  floatEventProperties["weightDown"] = fctVtxWeightDown_( iEvent );
  
  float pt1 = 0.;
  float pt2 = 0.;
  
  std::string leptonFlavor1 = "None";
  std::string leptonFlavor2 = "None";
  int leptonNr1=-1;
  int leptonNr2=-1;
  
  for(std::size_t it = 0; it < (*electrons).size() ; ++it){
	  if ((*electrons).at(it).pt() > pt1){
		  pt1 = (*electrons).at(it).pt();
		  leptonFlavor1 = "Ele";
		  leptonNr1 = it;
	  }
  }
  for(std::size_t it = 0; it < (*muons).size() ; ++it){
	  if ((*muons).at(it).pt() > pt1){
		  pt1 = (*muons).at(it).pt();
		  leptonFlavor1 = "Mu";
		  leptonNr1 = it;
	  }
  }
  for(std::size_t it = 0; it < (*electrons).size() ; ++it){
	  if ((*electrons).at(it).pt() < pt1 && (*electrons).at(it).pt() > pt2){
		  pt2 = (*electrons).at(it).pt();
		  leptonFlavor2 = "Ele";
		  leptonNr2 = it;
	  }
  }
  
  for(std::size_t it = 0; it < (*muons).size() ; ++it){
	  if ((*muons).at(it).pt() < pt1 && (*muons).at(it).pt() > pt2){
		  pt2 = (*muons).at(it).pt();
		  leptonFlavor2 = "Mu";
		  leptonNr2 = it;
	  }
  } 
 
  
  if (leptonFlavor1 == "Ele" && leptonFlavor2 == "Ele") {
	  fillTree<pat::Electron, pat::Electron>( "EE", (*electrons).at(leptonNr1), (*electrons).at(leptonNr2),*pfCands,*looseElectrons,*looseMuons,*jets,*shiftedJetsJESUp,*shiftedJetsJESDown,*bJets35,*shiftedBJetsJESUp,*shiftedBJetsJESDown, met,MHT,vertices,Rho, intEventProperties, floatEventProperties,tLorentzVectorEventProperties,isMC); 
	  
  }
  else if (leptonFlavor1 == "Mu" && leptonFlavor2 == "Mu") {
	  fillTree<pat::Muon, pat::Muon>( "MuMu", (*muons).at(leptonNr1), (*muons).at(leptonNr2),*pfCands,*looseElectrons,*looseMuons,*jets,*shiftedJetsJESUp,*shiftedJetsJESDown,*bJets35,*shiftedBJetsJESUp,*shiftedBJetsJESDown, met,MHT,vertices,Rho, intEventProperties, floatEventProperties,tLorentzVectorEventProperties,isMC); 
  }
  else if (leptonFlavor1 == "Ele" && leptonFlavor2 == "Mu") {
	  fillTree<pat::Electron, pat::Muon>( "EMu", (*electrons).at(leptonNr1), (*muons).at(leptonNr2),*pfCands,*looseElectrons,*looseMuons,*jets,*shiftedJetsJESUp,*shiftedJetsJESDown,*bJets35,*shiftedBJetsJESUp,*shiftedBJetsJESDown, met,MHT,vertices,Rho, intEventProperties, floatEventProperties,tLorentzVectorEventProperties,isMC); 
  }
  // Change ordering for Mu E events, in such a way that the electron is always the first lepton (required by some tools that select the lepton flavor)
  else if (leptonFlavor1 == "Mu" && leptonFlavor2 == "Ele") {
	  fillTree<pat::Electron, pat::Muon>( "EMu", (*electrons).at(leptonNr2), (*muons).at(leptonNr1),*pfCands,*looseElectrons,*looseMuons,*jets,*shiftedJetsJESUp,*shiftedJetsJESDown,*bJets35,*shiftedBJetsJESUp,*shiftedBJetsJESDown, met,MHT,vertices,Rho, intEventProperties, floatEventProperties,tLorentzVectorEventProperties,isMC); 
  }
	  

  delete jecUnc;
  delete shiftedJetsJESUp;
  delete shiftedJetsJESDown;
  delete shiftedBJetsJESUp;
  delete shiftedBJetsJESDown;


}


template <class aT, class bT> void 
DiLeptonTreesFromMiniAOD::fillTree( const std::string &treeName, const aT& a, const bT& b,const std::vector<pat::PackedCandidate>&pfCands,const std::vector<pat::Electron>&looseElectrons,const std::vector<pat::Muon>&looseMuons,const std::vector<pat::Jet>&jets,const std::vector<pat::Jet>&shiftedJetsUp,const std::vector<pat::Jet>&shiftedJetsDown,const std::vector<pat::Jet>&bJets35,const std::vector<pat::Jet>&shiftedBJetsUp,const std::vector<pat::Jet>&shiftedBJetsDown, const pat::MET &patMet,const TLorentzVector &MHT,const edm::Handle<reco::VertexCollection> &vertices,const float &rho, const std::map<std::string, long> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties, const bool &isMC)
{

  for(std::map<std::string, long>::const_iterator it = intEventProperties.begin(); it != intEventProperties.end(); ++it){
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
  
  if(debug) std::cout << treeName << "- pts:"<< a.pt() << " " << b.pt();
  TLorentzVector aVec = getMomentum(a);//( a.px(), a.py(), a.pz(), a.energy() );
  TLorentzVector bVec = getMomentum(b); //( b.px(), b.py(), b.pz(), b.energy() );
  TLorentzVector met(patMet.px(), patMet.py(), patMet.pz(), patMet.energy());
  TLorentzVector uncorrectedMet; 
  uncorrectedMet.SetPtEtaPhiE(patMet.uncorPt(), 0, \
			      patMet.uncorPhi(), patMet.uncorPt());
			      

  //  std::cout << "met: "<<met.Et()<< ", unCorr met: "<< uncorrectedMet.Et()
  //<< "=> "<< met.Et()* 1./uncorrectedMet.Et()<< " (xCheck: "<< patMet.corSumEt()*1./patMet.uncorrectedPt(pat::MET::uncorrALL) <<")"<<std::endl;

   //~ *(floatBranches_[treeName]["leptHT"]) += aVec.Pt();	
   //~ *(floatBranches_[treeName]["leptHT"]) += bVec.Pt();	
   

  TLorentzVector comb = aVec+bVec;
  TLorentzVector MHT2 = comb + MHT;
  
  
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
  
  // Calculate sum Mlb
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
  

  // Calculate sum Mlb with JES shifted up
  mlb_min = 1.E6;
  mlb_max = 1.E6;
  jet.SetPxPyPzE(0.,0.,0.,0.);
  
  if (shiftedBJetsUp.size() >= 1) jet1Coll = shiftedBJetsUp;
  else jet1Coll = shiftedJetsUp;
  
  if (shiftedBJetsUp.size() >= 2) jet2Coll = shiftedBJetsUp;
  else jet2Coll = shiftedJetsUp;
  
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
			  if (shiftedBJetsUp.size() == 1 && jet.DeltaR (jet1) < 0.1) continue;
			  if ( (shiftedBJetsUp.size() == 0 || shiftedBJetsUp.size() >= 2) && 	jet == jet1) continue;	
			  lepton.SetPxPyPzE(leptList[il].Px(),leptList[il].Py(),leptList[il].Pz(),leptList[il].Energy());
			  
			  temp_mlb = (jet+lepton).M();
			  if (temp_mlb < mlb_max) {
				  mlb_max = temp_mlb;
			  }
		  }
	  }
  }
	 
  if (mlb_min < 1.E6 and mlb_max < 1.E6) *(floatBranches_[treeName]["sumMlbJESUp"]) = mlb_min + mlb_max; 
  else *(floatBranches_[treeName]["sumMlbJESUp"]) = 1.E-6;
  

  // Calculate sum Mlb with JES shifted down
  mlb_min = 1.E6;
  mlb_max = 1.E6;
  jet.SetPxPyPzE(0.,0.,0.,0.);
  
  if (shiftedBJetsDown.size() >= 1) jet1Coll = shiftedBJetsDown;
  else jet1Coll = shiftedJetsDown;  

  if (shiftedBJetsDown.size() >= 2) jet2Coll = shiftedBJetsDown;
  else jet2Coll = shiftedJetsDown;
  
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
			  if (shiftedBJetsDown.size() == 1 && jet.DeltaR (jet1) < 0.1) continue;
			  if ( (shiftedBJetsDown.size() == 0 || shiftedBJetsDown.size() >= 2) && 	jet == jet1) continue;	
			  lepton.SetPxPyPzE(leptList[il].Px(),leptList[il].Py(),leptList[il].Pz(),leptList[il].Energy());
			  
			  temp_mlb = (jet+lepton).M();
			  if (temp_mlb < mlb_max) {
				  mlb_max = temp_mlb;
			  }
		  }
	  }
  }
	 
  if (mlb_min < 1.E6 and mlb_max < 1.E6) *(floatBranches_[treeName]["sumMlbJESDown"]) = mlb_min + mlb_max; 
  else *(floatBranches_[treeName]["sumMlbJESDown"]) = 1.E-6;

  
  double pa[3] = {aVec.M(),aVec.Px(),aVec.Py()};
  double pb[3] = {bVec.M(),bVec.Px(),bVec.Py()};
  
  double pmiss[3] = {0.,met.Px(),met.Py()};
  fctMT2_.set_momenta(pa,pb,pmiss);
  *(floatBranches_[treeName]["MT2"]) = static_cast<float>(fctMT2_.get_mt2()); 
  
  *(floatBranches_[treeName]["pt3"]) = 0.;
  
  if (*(intBranches_[treeName]["nLooseLeptons"]) > 2) {
	  float leptonPts[*(intBranches_[treeName]["nLooseLeptons"])];
	  int nLepton = 0;
	  
	  for(std::vector< pat::Electron >::const_iterator it = looseElectrons.begin(); it != looseElectrons.end() ; ++it){
		  TLorentzVector looseLeptonVector( (*it).px(), (*it).py(), (*it).pz(), (*it).energy() );
		  // Check if the loose lepton is the same as a
		  if ( looseLeptonVector.DeltaR( aVec) > 0.05 &&  looseLeptonVector.DeltaR( bVec) > 0.05)leptonPts[nLepton] = (*it).pt();
		  else leptonPts[nLepton] = 0.;
		  nLepton++;
	  }
	  for(std::vector< pat::Muon >::const_iterator it = looseMuons.begin(); it != looseMuons.end() ; ++it){
		  TLorentzVector looseLeptonVector( (*it).px(), (*it).py(), (*it).pz(), (*it).energy() );
		  if ( looseLeptonVector.DeltaR( aVec) > 0.05 &&  looseLeptonVector.DeltaR( bVec) > 0.05) leptonPts[nLepton] = (*it).pt();
		  else leptonPts[nLepton] = 0.;
		  nLepton++;
	  }
		  
	  std::sort(leptonPts, leptonPts + nLepton, std::greater<float>());
	  *(floatBranches_[treeName]["pt3"]) = leptonPts[0];	 
  }
  
  int nIsoTracks = 0;
  double absIso = 0.;
  
  std::vector < float > trackPts;
  
  for(std::vector<pat::PackedCandidate>::const_iterator it = pfCands.begin(); it != pfCands.end() ; ++it){
	
	  if ((*it).charge()==0) continue;
	  
	  if ( (*it).pt() < 5. ) continue;
	  if ( abs((*it).dz()) > 0.1 ) continue;
	  if ( (*it).fromPV() <= 1 ) continue;
	  
	  absIso = getIso(*it,"trackIso");
	  if (absIso >= min(0.2 * (*it).pt(),8.0) ) continue;
	  
	  if ( deltaR(a,(*it)) < 0.01 || deltaR(b,(*it)) < 0.01 ) continue;
	  
	  if ( abs((*it).pdgId()) == 11 || abs((*it).pdgId()) == 13 || ( (*it).pt() > 10. && abs((*it).pdgId()) == 211  && absIso < 0.1 * (*it).pt() ) ){
		  nIsoTracks++;
		  trackPts.push_back((*it).pt());
	  }
	  
  }
  *(intBranches_[treeName]["nIsoTracks"]) = nIsoTracks;
  if (nIsoTracks > 0) *(floatBranches_[treeName]["ptTrack3"]) = *std::max_element(std::begin(trackPts),std::end(trackPts));
  else *(floatBranches_[treeName]["ptTrack3"]) = 0.;
  
		  	  
  
  //~ std::pair<double, double> pZeta = calcPZeta(aVec, bVec, met);
  *(floatBranches_[treeName]["chargeProduct"]) = a.charge()*b.charge();
  *(tLorentzVectorBranches_[treeName]["p4"]) = comb;
  *(tLorentzVectorBranches_[treeName]["lepton1"]) = aVec;
  *(tLorentzVectorBranches_[treeName]["lepton2"]) = bVec;
    *(floatBranches_[treeName]["mht"]) = MHT2.Pt();
  *(floatBranches_[treeName]["pt1"]) = aVec.Pt();
  *(floatBranches_[treeName]["pt2"]) = bVec.Pt();
  *(floatBranches_[treeName]["ptErr1"]) = a.bestTrack()->ptError()/aVec.Pt();
  *(floatBranches_[treeName]["ptErr2"]) = b.bestTrack()->ptError()/bVec.Pt();
  *(floatBranches_[treeName]["charge1"]) = a.charge();
  *(floatBranches_[treeName]["charge2"]) = b.charge();
  *(floatBranches_[treeName]["eta1"]) = aVec.Eta();
  *(floatBranches_[treeName]["eta2"]) = bVec.Eta();
   //~ if (isMC){	
	  //~ *(floatBranches_[treeName]["leptonFullSimScaleFactor1"]) = fctLeptonFullSimScaleFactors_(a, a.pt(), a.eta() );
	  //~ *(floatBranches_[treeName]["leptonFullSimScaleFactor2"]) = fctLeptonFullSimScaleFactors_(b, b.pt(), b.eta() );
  //~ }
   //~ else{	
	  //~ *(floatBranches_[treeName]["leptonFullSimScaleFactor1"]) = 1.;
	  //~ *(floatBranches_[treeName]["leptonFullSimScaleFactor2"]) = 1.;
  //~ }
  *(floatBranches_[treeName]["miniIsoEffArea1"]) = getIso(a,"miniIsoEA");
  *(floatBranches_[treeName]["miniIsoEffArea2"]) = getIso(b,"miniIsoEA");
  *(floatBranches_[treeName]["mt1"]) = transverseMass(aVec, met);
  *(floatBranches_[treeName]["mt2"]) = transverseMass(bVec, met);
  //~ *(floatBranches_[treeName]["fakeWeight1"]) = fakeRates_(a);
  //~ *(floatBranches_[treeName]["fakeWeight2"]) = fakeRates_(b);
  *(floatBranches_[treeName]["deltaPhi"]) = aVec.DeltaPhi( bVec );
  *(floatBranches_[treeName]["deltaR"]) = aVec.DeltaR( bVec );
  //~ *(floatBranches_[treeName]["angle3D"]) = aVec.Angle( bVec.Vect() );
  *(floatBranches_[treeName]["jzb"]) = (met+comb).Pt() - comb.Pt();
  *(floatBranches_[treeName]["jzbUncorr"]) = (uncorrectedMet+comb).Pt() - comb.Pt();  
  *(floatBranches_[treeName]["jzbResponse"]) = (met+comb).Pt() / comb.Pt();
  *(floatBranches_[treeName]["jzbResponseUncorr"]) = (uncorrectedMet+comb).Pt() / comb.Pt();    
 
  if (writeID_){
	  fillLeptonIDs(treeName, a , b , vertices );
  }


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
  //~ *(floatBranches_[treeName]["genLeptHT"]) += genLepton1.Pt();	
  //~ *(floatBranches_[treeName]["genLeptHT"]) += genLepton2.Pt();	
  
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


void DiLeptonTreesFromMiniAOD::fillPdfUncert(const edm::Handle< std::vector<double> >& weightHandle, const std::string& pdfIdentifier, const std::string& treeName){
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

const TLorentzVector DiLeptonTreesFromMiniAOD::getMomentum(const  pat::Electron &e)
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

const TLorentzVector DiLeptonTreesFromMiniAOD::getMomentum(const  pat::Muon &mu)
{
  const TLorentzVector result = TLorentzVector(mu.px(), mu.py(), mu.pz(), mu.energy());
  return result;
}

float DiLeptonTreesFromMiniAOD::getIso(const  pat::Electron &e, const std::string &method)
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

float DiLeptonTreesFromMiniAOD::getIso(const  pat::Muon &mu, const std::string &method)
{
  //  std::cout<<"muon " << (lepton.trackIso()+lepton.ecalIso()+lepton.hcalIso())/lepton.pt() << std::endl;
  //  return (mu.isolationR03().hadEt + mu.isolationR03().emEt + mu.isolationR03().sumPt) / mu.pt();
  //  return (mu.chargedHadronIso() + mu.photonIso() + mu.neutralHadronIso()) / mu.pt();
  return fctIsolation_(mu,method)* 1./mu.pt();
}

float DiLeptonTreesFromMiniAOD::getIso(const  pat::PackedCandidate &track, const std::string &method)
{
  return fctIsolation_(track,method);
}

void DiLeptonTreesFromMiniAOD::fillLeptonIDs(const std::string &treeName, const  pat::Electron &ele1, const  pat::Electron &ele2, const edm::Handle<reco::VertexCollection> &vertices)
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
  
  
  *(floatBranches_[treeName]["deltaEtaSuperClusterTrackAtVtx1"]) = ele1.deltaEtaSuperClusterTrackAtVtx();
  *(floatBranches_[treeName]["deltaPhiSuperClusterTrackAtVtx1"]) = ele1.deltaPhiSuperClusterTrackAtVtx();
  *(floatBranches_[treeName]["sigmaIetaIeta1"]) = ele1.sigmaIetaIeta();
  *(floatBranches_[treeName]["hadronicOverEm1"]) = ele1.hadronicOverEm();
  *(floatBranches_[treeName]["eOverP1"]) = abs(1.0/ele1.ecalEnergy() - ele1.eSuperClusterOverP()/ele1.ecalEnergy());
  *(floatBranches_[treeName]["missingHits1"]) = ele1.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
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


void DiLeptonTreesFromMiniAOD::fillLeptonIDs(const std::string &treeName, const  pat::Electron &ele1, const  pat::Muon &mu2, const edm::Handle<reco::VertexCollection> &vertices)
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

void DiLeptonTreesFromMiniAOD::fillLeptonIDs(const std::string &treeName, const  pat::Muon &mu1, const  pat::Muon &mu2, const edm::Handle<reco::VertexCollection> &vertices)
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


float DiLeptonTreesFromMiniAOD::topPtWeightBen(double topPt){
  if( topPt<0 ) return 1;

  float p0 = 1.18246e+00;
  float p1 = 4.63312e+02;
  float p2 = 2.10061e-06;

  if( topPt>p1 ) topPt = p1;

  float result = p0 + p2 * topPt * ( topPt - 2 * p1 );
  return result;
}

float DiLeptonTreesFromMiniAOD::topPtWeightTOP(double topPt){
  if( topPt<0 ) return 1;

  float p0 = 0.156;
  float p1 = 0.00137;

  float result = exp(p0 - p1 * topPt);
  return result;
}

float DiLeptonTreesFromMiniAOD::getAEffEle(double eta)
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

float DiLeptonTreesFromMiniAOD::getAEffMu(double eta)
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


float DiLeptonTreesFromMiniAOD::transverseMass(const TLorentzVector& p, const TLorentzVector& met)
{
  reco::Candidate::LorentzVector otherMet(met.Px(),met.Py(),met.Pz(),met.E());
  reco::Candidate::LorentzVector leptonT(p.Px(),p.Py(),0.,p.E()*sin(p.Theta()));
  reco::Candidate::LorentzVector sumT=leptonT+otherMet;

  return std::sqrt(sumT.M2());
}

// ------------ method called when starting to processes a luminosity block  ------------

void
DiLeptonTreesFromMiniAOD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
   newLumiBlock_=true;
}

std::string DiLeptonTreesFromMiniAOD::convertInputTag(const edm::InputTag tag)
{
  std::string result = tag.label();
  if(tag.instance().length() > 0)
    result = tag.instance();
  //  std::cerr << "'"<<tag.label() << "', '"<< tag.instance()<<"' = '"<< result<<"'"<<std::endl;
  return result;
}


// ------------ Method called once each job just before starting event loop  ------------
void 
DiLeptonTreesFromMiniAOD::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiLeptonTreesFromMiniAOD::endJob() 
{
}



//define this as a plug-in
DEFINE_FWK_MODULE(DiLeptonTreesFromMiniAOD);
