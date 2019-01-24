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
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/Lepton.h>
#include <DataFormats/PatCandidates/interface/IsolatedTrack.h>

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
#include <SuSyAachen/DiLeptonHistograms/interface/IsolationFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/TriggerMatchFunctorMiniAOD.h>
#include <SuSyAachen/DiLeptonHistograms/interface/LeptonFullSimScaleFactorMapFunctor.h>
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


class DiLeptonTreesFromMiniAOD : public edm::one::EDAnalyzer<edm::one::WatchLuminosityBlocks,edm::one::SharedResources> {
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
  void initLongIntBranch( const std::string &name);
  void initTLorentzVectorBranch( const std::string &name);
  template<class aT, class bT> void fillTree( const std::string &treeName, const aT &a, const bT &b, const std::vector<pat::IsolatedTrack>&isoTracks,const std::vector<pat::Electron>&looseElectrons,const std::vector<pat::Muon>&looseMuons,const std::vector<pat::Jet>&jets,const std::vector<pat::Jet>&shiftedJetsJESUp,const std::vector<pat::Jet>&shiftedJetsJESDown,const std::vector<pat::Jet>&bJets35,const std::vector<pat::Jet>&shiftedJetsBJESUp,const std::vector<pat::Jet>&shiftedBJetsJESDown,  const pat::MET &patMet, const TLorentzVector &MHT, const edm::Handle<reco::VertexCollection> &vertices, const float &rho, const std::map<std::string, int> &intEventProperties, const std::map<std::string, unsigned long> &longIntEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties, const bool &isMC);
  void sumMlb(TLorentzVector &lepton1, TLorentzVector &lepton2, const std::vector<pat::Jet> &jets, const std::vector<pat::Jet> &bjets, float &result_sum_mlb, float &result_mlb_min, float &result_mlb_max);
  const TLorentzVector getMomentum(const  pat::Electron &e);
  const TLorentzVector getMomentum(const  pat::Muon &mu);
  float getIso(const  pat::Electron &e, const std::string &method);
  float getIso(const  pat::Muon &mu, const std::string &method);
  float getIso(const  pat::PackedCandidate &track, const std::string &method);
  float transverseMass(const TLorentzVector& p, const TLorentzVector& met);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;  

  edm::EDGetTokenT< std::vector< pat::Electron > >      electronToken_;
  edm::EDGetTokenT< std::vector< pat::Electron > >      looseElectronToken_;
  edm::EDGetTokenT< std::vector< pat::Muon > >        muonToken_;
  edm::EDGetTokenT< std::vector< pat::Muon > >        looseMuonToken_;
  edm::EDGetTokenT< std::vector< pat::Jet > >       jetToken_;
  edm::EDGetTokenT< std::vector< reco::GenJet > >     genJetToken_;
  edm::EDGetTokenT< std::vector< pat::Jet > >       bJetToken_;
  edm::EDGetTokenT< std::vector< pat::Jet > >       bJet35Token_;
  edm::EDGetTokenT< std::vector< pat::MET > >         metToken_;
  edm::EDGetTokenT< std::vector< pat::MET > >         metRawToken_;
  edm::EDGetTokenT<reco::VertexCollection>          vertexToken_;
  edm::EDGetTokenT< std::vector< pat::IsolatedTrack >  >  isoTrackToken_;
  edm::EDGetTokenT< std::vector< reco::GenParticle > >    genParticleToken_;
  edm::EDGetTokenT<GenEventInfoProduct>           genEventInfoToken_;
  edm::EDGetTokenT<LHEEventProduct>             LHEEventToken_;
  edm::EDGetTokenT<double>                  rhoToken_;
  
  edm::EDGetTokenT< double > prefweight_token;
  edm::EDGetTokenT< double > prefweightup_token;
  edm::EDGetTokenT< double > prefweightdown_token;
  
  edm::EDGetTokenT<edm::TriggerResults>           metFilterToken_;
  

  std::map<double, double> electronCorrections_;
  //data
  std::map<std::string, TTree*> trees_;  
  std::map<std::string, std::map< std::string, float*> > floatBranches_; 
  std::map<std::string, std::map< std::string, int*> > intBranches_; 
  std::map<std::string, std::map< std::string, unsigned long*> > longIntBranches_; 
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
  
  LeptonFullSimScaleFactorMapFunctor fctLeptonFullSimScaleFactors_;

  VertexWeightFunctor fctVtxWeight_;
  VertexWeightFunctor fctVtxWeightUp_;
  VertexWeightFunctor fctVtxWeightDown_;
  IsolationFunctor fctIsolation_; 
  TriggerMatchFunctorMiniAOD fctTrigger_;  
  PdgIdFunctor getPdgId_;
  MT2Functor fctMT2_;
 
  bool debug;
  bool writeTrigger_;
  bool metUncert_;  
  bool storeMetFilters_;
  std::vector<std::string> metFilterNames_;
  std::vector<std::string> eeTriggerNames_;
  std::vector<std::string> emTriggerNames_;
  std::vector<std::string> mmTriggerNames_;
  std::vector<std::string> htTriggerNames_;
  std::vector<std::string> metTriggerNames_;
  bool eeNewLumiBlock_;
  bool emNewLumiBlock_;
  bool mmNewLumiBlock_;
  bool htNewLumiBlock_;
  bool metNewLumiBlock_;
  std::map<std::string, Bool_t > eeTriggerDecision_;
  std::map<std::string, Bool_t > emTriggerDecision_;
  std::map<std::string, Bool_t > mmTriggerDecision_;
  std::map<std::string, Bool_t > htTriggerDecision_;
  std::map<std::string, Bool_t > metTriggerDecision_;
  std::map<std::string, int > eeTriggerIndex_;
  std::map<std::string, int > emTriggerIndex_;
  std::map<std::string, int > mmTriggerIndex_;
  std::map<std::string, int > htTriggerIndex_;
  std::map<std::string, int > metTriggerIndex_;
  
};

// constructors and destructor
DiLeptonTreesFromMiniAOD::DiLeptonTreesFromMiniAOD(const edm::ParameterSet& iConfig):
  electronToken_      (consumes< std::vector< pat::Electron > >     (iConfig.getParameter<edm::InputTag>("electrons"))),
  looseElectronToken_   (consumes< std::vector< pat::Electron > >     (iConfig.getParameter<edm::InputTag>("looseElectrons"))),
  muonToken_        (consumes< std::vector< pat::Muon > >     (iConfig.getParameter<edm::InputTag>("muons"))),
  looseMuonToken_     (consumes< std::vector< pat::Muon > >     (iConfig.getParameter<edm::InputTag>("looseMuons"))),
  jetToken_         (consumes< std::vector< pat::Jet > >      (iConfig.getParameter<edm::InputTag>("jets"))),
  genJetToken_        (consumes< std::vector< reco::GenJet  > >   (iConfig.getParameter<edm::InputTag>("genJets"))),
  bJetToken_        (consumes< std::vector< pat::Jet > >      (iConfig.getParameter<edm::InputTag>("bJets"))),
  bJet35Token_        (consumes< std::vector< pat::Jet > >      (iConfig.getParameter<edm::InputTag>("bJets35"))),
  metToken_         (consumes< std::vector< pat::MET > >      (iConfig.getParameter<edm::InputTag>("met"))),
  metRawToken_         (consumes< std::vector< pat::MET > >      (iConfig.getParameter<edm::InputTag>("met_normal"))),
  vertexToken_        (consumes<reco::VertexCollection>       (iConfig.getParameter<edm::InputTag>("vertices"))),
  isoTrackToken_        (consumes< std::vector<pat::IsolatedTrack>  > (iConfig.getParameter<edm::InputTag>("isoTracks"))),
  genParticleToken_     (consumes< std::vector< reco::GenParticle > > (iConfig.getParameter<edm::InputTag>("genParticles"))),
  genEventInfoToken_    (consumes<GenEventInfoProduct>          (iConfig.getParameter<edm::InputTag>("pdfInfo"))),
  LHEEventToken_      (consumes<LHEEventProduct>            (iConfig.getParameter<edm::InputTag>("LHEInfo"))),
  rhoToken_         (consumes<double>               (iConfig.getParameter<edm::InputTag>("rho"))),
  
  metFilterToken_ (consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","",edm::InputTag::kSkipCurrentProcess))),
  
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
  
  fctLeptonFullSimScaleFactors_ (iConfig.getParameter<edm::ParameterSet>("LeptonFullSimScaleFactors") ),
  
  fctVtxWeight_     (iConfig.getParameter<edm::ParameterSet>("vertexWeights") ,consumesCollector()),
  fctVtxWeightUp_   (iConfig.getParameter<edm::ParameterSet>("vertexWeightsUp") ,consumesCollector()),
  fctVtxWeightDown_ (iConfig.getParameter<edm::ParameterSet>("vertexWeightsDown"),consumesCollector() ),
  fctIsolation_   (iConfig.getParameter<edm::ParameterSet>("isolationDefinitions"),consumesCollector()),
  fctTrigger_     (iConfig.getParameter<edm::ParameterSet>("triggerDefinitions"),consumesCollector()),  
  getPdgId_     (iConfig.getParameter< edm::ParameterSet>("pdgIdDefinition"),consumesCollector() ),
  metFilterNames_ (iConfig.getUntrackedParameter< std::vector <std::string> >("metFilterNames")),
  eeTriggerNames_   (iConfig.getUntrackedParameter< std::vector <std::string> >("eeTriggerNames")),
  emTriggerNames_   (iConfig.getUntrackedParameter< std::vector <std::string> >("emTriggerNames")),
  mmTriggerNames_   (iConfig.getUntrackedParameter< std::vector <std::string> >("mmTriggerNames")),
  htTriggerNames_   (iConfig.getUntrackedParameter< std::vector <std::string> >("htTriggerNames")),
  metTriggerNames_   (iConfig.getUntrackedParameter< std::vector <std::string> >("metTriggerNames")),
  eeNewLumiBlock_(true),
  emNewLumiBlock_(true),
  mmNewLumiBlock_(true),
  htNewLumiBlock_(true),
  metNewLumiBlock_(true)

{
  usesResource("TFileService");
  debug = false; 
  writeTrigger_ = iConfig.getUntrackedParameter<bool>("writeTrigger");
  storeMetFilters_ = iConfig.getUntrackedParameter<bool>("storeMetFilters");
  metUncert_ = iConfig.getUntrackedParameter<bool>("doMETUncert");  
  
  consumes<std::vector< PileupSummaryInfo > >(edm::InputTag("slimmedAddPileupInfo"));
  
  prefweight_token = consumes< double >(edm::InputTag("prefiringweight:NonPrefiringProb"));
  prefweightup_token = consumes< double >(edm::InputTag("prefiringweight:NonPrefiringProbUp"));
  prefweightdown_token = consumes< double >(edm::InputTag("prefiringweight:NonPrefiringProbDown"));
  
  consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
  consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","",edm::InputTag::kSkipCurrentProcess));
  
  //~ consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"));
  //~ consumes<std::vector<pat::TriggerObjectStandAlone>>(edm::InputTag("selectedPatTrigger"));

  // init trees
  edm::Service<TFileService> file;
  trees_["EE"] = file->make<TTree>("EEDileptonTree", "EE DileponTree");
  trees_["EMu"] = file->make<TTree>("EMuDileptonTree", "EMu DileponTree");
  trees_["MuMu"] = file->make<TTree>("MuMuDileptonTree", "MuMu DileponTree");
  
  initFloatBranch( "mll" );  
  initFloatBranch( "genWeight" );  
  initFloatBranch( "prefireWeight" );   
  initFloatBranch( "genWeightAbsValue" );    
  initFloatBranch( "weight" );
  initFloatBranch( "weightUp" );
  initFloatBranch( "weightDown" );
  initFloatBranch( "bTagWeight" );
  initFloatBranch( "bTagWeightErrHeavy" );
  initFloatBranch( "bTagWeightErrLight" );
  initFloatBranch( "leptonFullSimScaleFactor1" );
  initFloatBranch( "leptonFullSimScaleFactor2" );
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
  initFloatBranch( "deltaPhi" );
  initFloatBranch( "deltaR" );
  initFloatBranch( "min_mlb" );
  initFloatBranch( "max_mlb" );
  initFloatBranch( "sumMlb" );
  initFloatBranch( "sumMlbJESUp" );
  initFloatBranch( "sumMlbJESDown" );
  initFloatBranch( "MT2" );
  initFloatBranch( "deltaPhiJetMet1" );
  initFloatBranch( "deltaPhiJetMet2" );
  initFloatBranch( "ht" );
  initFloatBranch( "htJESUp" );
  initFloatBranch( "htJESDown" ); 
  initFloatBranch( "mht" );
  initFloatBranch( "met" );
  initFloatBranch( "met_raw" );
  initFloatBranch( "caloMet" );
  initFloatBranch( "genMet" );
  initFloatBranch( "uncorrectedMet" ); 
  initFloatBranch( "metJESUp" );
  initFloatBranch( "metJESDown" );
  initIntBranch( "nJets" );
  initIntBranch( "nGenJets" );
  initIntBranch( "nBadMuonJets" );
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
  initIntBranch( "runNr" );
  initIntBranch( "lumiSec" );
  initLongIntBranch( "eventNr" );
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
  
  initIntBranch( "triggerSummary" );
  initIntBranch( "triggerSummaryHT" );
  initIntBranch( "triggerSummaryMET" );
  initIntBranch( "metFilterSummary" );
  
  if (storeMetFilters_){
  
    for (const auto& n : metFilterNames_){
      initIntBranch( n.c_str() );
    }
    
  }
  
  if (writeTrigger_){
    for (const auto& n : eeTriggerNames_){
      eeTriggerIndex_[n] = -10; //not set and not found
      eeTriggerDecision_[n] = false;

      for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();it != trees_.end(); ++it){
        if(debug) std::cout << (*it).first <<" - "<< n.c_str() << std::endl;
        boolBranches_[(*it).first][n.c_str()] = new bool;
        (*it).second->Branch(n.c_str(), &eeTriggerDecision_[n] ,(n+"/O").c_str() );
      }
    }
    for (const auto& n : emTriggerNames_){
      emTriggerIndex_[n] = -10; //not set and not found
      emTriggerDecision_[n] = false;

      for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();it != trees_.end(); ++it){
        if(debug) std::cout << (*it).first <<" - "<< n.c_str() << std::endl;
        boolBranches_[(*it).first][n.c_str()] = new bool;
        (*it).second->Branch(n.c_str(), &emTriggerDecision_[n] ,(n+"/O").c_str() );
      }
    }
    for (const auto& n : mmTriggerNames_){
      mmTriggerIndex_[n] = -10; //not set and not found
      mmTriggerDecision_[n] = false;

      for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();it != trees_.end(); ++it){
        if(debug) std::cout << (*it).first <<" - "<< n.c_str() << std::endl;
        boolBranches_[(*it).first][n.c_str()] = new bool;
        (*it).second->Branch(n.c_str(), &mmTriggerDecision_[n] ,(n+"/O").c_str() );
      }
    }
    for (const auto& n : htTriggerNames_){
      htTriggerIndex_[n] = -10; //not set and not found
      htTriggerDecision_[n] = false;

      for( const auto& it : trees_){
        if(debug) std::cout << it.first <<" - "<< n.c_str() << std::endl;
        boolBranches_[it.first][n.c_str()] = new bool;
        it.second->Branch(n.c_str(), &htTriggerDecision_[n] ,(n+"/O").c_str() );
      }
    }
    for (const auto& n : metTriggerNames_){
      metTriggerIndex_[n] = -10; //not set and not found
      metTriggerDecision_[n] = false;

      for( const auto& it : trees_){
        if(debug) std::cout << it.first <<" - "<< n.c_str() << std::endl;
        boolBranches_[it.first][n.c_str()] = new bool;
        it.second->Branch(n.c_str(), &metTriggerDecision_[n] ,(n+"/O").c_str() );
      }
    }
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
  for( const auto& it : trees_){
    if(debug) std::cout << it.first <<" - "<< name << std::endl;
    tLorentzVectorBranches_[it.first][name] = new TLorentzVector;
    it.second->Branch(name.c_str(), "TLorentzVector" ,&tLorentzVectorBranches_[it.first][name]);
  }
}

void 
DiLeptonTreesFromMiniAOD::initFloatBranch(const std::string &name)
{
  for( const auto& it : trees_){
    if(debug) std::cout << it.first <<" - "<< name << std::endl;
    floatBranches_[it.first][name] = new float;
    it.second->Branch(name.c_str(), floatBranches_[it.first][name], (name+"/F").c_str());
  }
}

void 
DiLeptonTreesFromMiniAOD::initIntBranch(const std::string &name)
{
  for( const auto& it : trees_){
    if(debug) std::cout << it.first <<" - "<< name << std::endl;
    intBranches_[it.first][name] = new int;
    it.second->Branch(name.c_str(), intBranches_[it.first][name], (name+"/I").c_str());
  }
}

void 
DiLeptonTreesFromMiniAOD::initLongIntBranch(const std::string &name)
{
  for( const auto& it : trees_){
    if(debug) std::cout << it.first <<" - "<< name << std::endl;
    longIntBranches_[it.first][name] = new unsigned long;
    it.second->Branch(name.c_str(), longIntBranches_[it.first][name], (name+"/l").c_str());
  }
}

DiLeptonTreesFromMiniAOD::~DiLeptonTreesFromMiniAOD()
{ 
  for( const auto& it: tLorentzVectorBranches_){
    for( const auto& it2 : it.second){
      if(debug)std::cout << "deleting: " << it.first << " - "<< it2.first << std::endl;
      delete it2.second;
    }
  }
  for( const auto& it: floatBranches_){
    for( const auto& it2 : it.second){
      if(debug)std::cout << "deleting: " << it.first << " - "<< it2.first << std::endl;
      delete it2.second;
    }
  }
  
  for( const auto& it: intBranches_){
    for( const auto& it2 : it.second){
      if(debug)std::cout << "deleting: " << it.first << " - "<< it2.first << std::endl;
      delete it2.second;
    }
  }
  
  for( const auto& it: longIntBranches_){
    for( const auto& it2 : it.second){
      if(debug)std::cout << "deleting: " << it.first << " - "<< it2.first << std::endl;
      delete it2.second;
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


  edm::Handle< std::vector<pat::IsolatedTrack>  > isoTracks;
  iEvent.getByToken(isoTrackToken_, isoTracks); 

  
  iEvent.getByToken(jetToken_, jets);

  edm::Handle< std::vector< reco::GenJet  > > genJets;
  iEvent.getByToken(genJetToken_, genJets);

  edm::Handle< std::vector< pat::Jet > > bJets;
  iEvent.getByToken(bJetToken_, bJets);
  
  edm::Handle< std::vector< pat::Jet > > bJets35;
  iEvent.getByToken(bJet35Token_, bJets35);
  
  edm::Handle< std::vector< pat::MET > > mets;
  iEvent.getByToken(metToken_, mets);
  
  edm::Handle< std::vector< pat::MET > > mets_raw;
  iEvent.getByToken(metRawToken_, mets_raw);

  edm::Handle< std::vector< reco::GenParticle > > genParticles;
  iEvent.getByToken(genParticleToken_, genParticles);
  
  
  bool isMC = genParticles.isValid();
  
   
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);

  getPdgId_.loadGenParticles(iEvent);
  fctIsolation_.init(iEvent);
  fctTrigger_.loadTrigger(iEvent);
  std::map<std::string, int> intEventProperties;
  std::map<std::string, unsigned long> longIntEventProperties;
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
  
  edm::Handle< double > theprefweight;
  iEvent.getByToken(prefweight_token, theprefweight ) ;
  double _prefiringweight =(*theprefweight);
  if (isMC){
    floatEventProperties["prefireWeight"] = _prefiringweight;
  }else{
    floatEventProperties["prefireWeight"] = 1.0;
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
  
  edm::Handle<edm::TriggerResults> metFilterBits;
  iEvent.getByToken(metFilterToken_, metFilterBits);
  const edm::TriggerNames &allFilterNames = iEvent.triggerNames(*metFilterBits);
  
  
  int metFilterSummary = 1;
  if (storeMetFilters_){ 
    
    for (std::string const &filterName : metFilterNames_){
      const unsigned index = allFilterNames.triggerIndex(filterName);
      if (index >= allFilterNames.size()) std::cerr << "MET filter " << filterName << "not found!" << std::endl;
      if (metFilterBits->accept(index)){
        intEventProperties[filterName] = 1;
      }else{
        intEventProperties[filterName] = 0;
        metFilterSummary = 0;
      }
    }  
  }else{ 
    
    for (std::string const &filterName : metFilterNames_){
      const unsigned index = allFilterNames.triggerIndex(filterName);
      if (index >= allFilterNames.size()) std::cerr << "MET filter " << filterName << "not found!" << std::endl;
      if (!metFilterBits->accept(index)) return;
    }  
  }
  intEventProperties["metFilterSummary"] = metFilterSummary; 


  intEventProperties["nBJets"] = bJets->size();
  intEventProperties["nBJets35"] = bJets35->size();
  intEventProperties["nLightLeptons"] = electrons->size() + muons->size();
  intEventProperties["nLooseLeptons"] = looseElectrons->size() + looseMuons->size();
  intEventProperties["runNr"] = iEvent.id().run();
  intEventProperties["lumiSec"] = iEvent.id().luminosityBlock();
  longIntEventProperties["eventNr"] = iEvent.id().event();
    

  
  pat::MET met = mets->front();
  TLorentzVector metVector(met.px(), met.py(), met.pz(), met.energy());
  
  pat::MET met_raw = mets_raw->front();
  TLorentzVector metRawVector(met_raw.px(), met_raw.py(), met_raw.pz(), met_raw.energy());
  
  TLorentzVector uncorrectedMetVector;
  uncorrectedMetVector.SetPtEtaPhiE(met.uncorPt(), 0, met.uncorPhi(), met.uncorPt());
  
  
  floatEventProperties["met"] = metVector.Pt();
  floatEventProperties["met_raw"] = metRawVector.Pt();
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
  std::vector<reco::GenJet> genJetsCleaned;
  
  floatEventProperties["genJet1pt"] = -1.0;
  floatEventProperties["genJet2pt"] = -1.0;
  floatEventProperties["genJet3pt"] = -1.0;
  floatEventProperties["genJet4pt"] = -1.0;
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
          genJetsCleaned.push_back(*it);
          if (nGenJets == 1) floatEventProperties["genJet1pt"] =(*it).pt();
          if (nGenJets == 2) floatEventProperties["genJet2pt"] =(*it).pt();
          if (nGenJets == 3) floatEventProperties["genJet3pt"] =(*it).pt();
          if (nGenJets == 4) floatEventProperties["genJet4pt"] =(*it).pt();
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
    genJet1Vector.SetPxPyPzE(genJetsCleaned.at(0).px(),genJetsCleaned.at(0).py(),genJetsCleaned.at(0).pz(),genJetsCleaned.at(0).energy());
  if (nGenJets > 1)
    genJet2Vector.SetPxPyPzE(genJetsCleaned.at(1).px(),genJetsCleaned.at(1).py(),genJetsCleaned.at(1).pz(),genJetsCleaned.at(1).energy());
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
  
  floatEventProperties["bTagWeight"] = 1.;
  floatEventProperties["bTagWeightErrHeavy"] = 0.;
  floatEventProperties["bTagWeightErrLight"] = 0.;
  
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
        float temp_pt = std::min(599.,(*it).pt());
        
        jetFlavor = (*it).hadronFlavour();
        if (jetFlavor == 5)
        {
          SF = fctBTagCalibReaderFullSimBJets_.eval(BTagEntry::FLAV_B, (*it).eta(), temp_pt);
          SF_up = fctBTagCalibReaderFullSimBJetsUp_.eval(BTagEntry::FLAV_B, (*it).eta(), temp_pt);
          SF_down = fctBTagCalibReaderFullSimBJetsDown_.eval(BTagEntry::FLAV_B, (*it).eta(), temp_pt);
          
         
        }
        else if (jetFlavor == 4)
        {
          SF = fctBTagCalibReaderFullSimCJets_.eval(BTagEntry::FLAV_C, (*it).eta(), temp_pt);
          SF_up = fctBTagCalibReaderFullSimCJetsUp_.eval(BTagEntry::FLAV_C, (*it).eta(), temp_pt);
          SF_down = fctBTagCalibReaderFullSimCJetsDown_.eval(BTagEntry::FLAV_C, (*it).eta(), temp_pt);

        }
        else
        {
          SF = fctBTagCalibReaderFullSimLightJets_.eval(BTagEntry::FLAV_UDSG, (*it).eta(), temp_pt);
          SF_up = fctBTagCalibReaderFullSimLightJetsUp_.eval(BTagEntry::FLAV_UDSG, (*it).eta(), temp_pt);
          SF_down = fctBTagCalibReaderFullSimLightJetsDown_.eval(BTagEntry::FLAV_UDSG, (*it).eta(), temp_pt);

        }
        
              
        eff = fctBTagEff_(jetFlavor, temp_pt, fabs((*it).eta()));
        
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
 
  int eeTriggerSummary = 0;
  int emTriggerSummary = 0;
  int mmTriggerSummary = 0;
  int htTriggerSummary = 0;
  int metTriggerSummary = 0;
  
  
  if (writeTrigger_){ 
    edm::Handle<edm::TriggerResults> triggerBits;
    edm::InputTag triggerTag("TriggerResults","","HLT");
    iEvent.getByLabel(triggerTag, triggerBits);
    
    // for each lumiBlock, re-read the trigger indices (rather changes for new run)
    if( eeTriggerIndex_.size() && eeNewLumiBlock_ ) {
      eeNewLumiBlock_=false;
      // set all trigger indices to -1 as "not available"-flag
      for( auto& it : eeTriggerIndex_ )
        it.second = -1;
  
      // store the indices of the trigger names that we really find
      const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
      for( unsigned i=0; i<triggerNames.size(); i++ ) {
          for( auto& it : eeTriggerIndex_ ) {
            if( triggerNames.triggerName(i).find( it.first ) == 0 ) {
                it.second = i;
            }
          }
      } // end trigger names
    } // found indices
     
    // for each lumiBlock, re-read the trigger indices (rather changes for new run)
    if( emTriggerIndex_.size() && emNewLumiBlock_ ) {
      emNewLumiBlock_=false;
      // set all trigger indices to -1 as "not available"-flag
      for( auto& it : emTriggerIndex_ )
        it.second = -1;
  
      // store the indices of the trigger names that we really find
      const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
      for( unsigned i=0; i<triggerNames.size(); i++ ) {
          for( auto& it : emTriggerIndex_ ) {
            if( triggerNames.triggerName(i).find( it.first ) == 0 ) {
                it.second = i;
            }
          }
      } // end trigger names
    } // found indices
     
    // for each lumiBlock, re-read the trigger indices (rather changes for new run)
    if( mmTriggerIndex_.size() && mmNewLumiBlock_ ) {
      mmNewLumiBlock_=false;
      // set all trigger indices to -1 as "not available"-flag
      for( auto& it : mmTriggerIndex_ )
        it.second = -1;
  
      // store the indices of the trigger names that we really find
      const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
      for( unsigned i=0; i<triggerNames.size(); i++ ) {
          for( auto& it : mmTriggerIndex_ ) {
            if( triggerNames.triggerName(i).find( it.first ) == 0 ) {
                it.second = i;
            }
          }
      } // end trigger names
    } // found indices
    
    // for each lumiBlock, re-read the trigger indices (rather changes for new run)
    if( htTriggerIndex_.size() && htNewLumiBlock_ ) {
      htNewLumiBlock_=false;
      // set all trigger indices to -1 as "not available"-flag
      for( auto& it : htTriggerIndex_ )
        it.second = -1;
  
      // store the indices of the trigger names that we really find
      const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
      for( unsigned i=0; i<triggerNames.size(); i++ ) {
          for( auto& it : htTriggerIndex_ ) {
            if( triggerNames.triggerName(i).find( it.first ) == 0 ) {
                it.second = i;
            }
          }
      } // end trigger names
    } // found indices
    
    if( metTriggerIndex_.size() && metNewLumiBlock_ ) {
      metNewLumiBlock_=false;
      // set all trigger indices to -1 as "not available"-flag
      for( auto& it : metTriggerIndex_ )
        it.second = -1;
  
      // store the indices of the trigger names that we really find
      const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
      for( unsigned i=0; i<triggerNames.size(); i++ ) {
          for( auto& it : metTriggerIndex_ ) {
            if( triggerNames.triggerName(i).find( it.first ) == 0 ) {
                it.second = i;
            }
          }
      } // end trigger names
    } // found indices
    
    // set trigger decision
    for( auto& it : eeTriggerIndex_ ) {
      if( it.second != -1 ) {
          eeTriggerDecision_[it.first] = triggerBits->accept( it.second );
          if (eeTriggerDecision_[it.first]){
            eeTriggerSummary = 1;
          }
      }
    }
    for( auto& it : emTriggerIndex_ ) {
      if( it.second != -1 ) {
          emTriggerDecision_[it.first] = triggerBits->accept( it.second );
          if (emTriggerDecision_[it.first]){
            emTriggerSummary = 1;
          }
      }
    }
    for( auto& it : mmTriggerIndex_ ) {
      if( it.second != -1 ) {
          mmTriggerDecision_[it.first] = triggerBits->accept( it.second );
          if (mmTriggerDecision_[it.first]){
            if (it.first.compare("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") == 0){
              if (iEvent.id().run() <= 299329){
                mmTriggerSummary = 1;
              }
            }else{
              mmTriggerSummary = 1;
            }
          }
      }
    }
    for( auto& it : htTriggerIndex_ ) {
      if( it.second != -1 ) {
          htTriggerDecision_[it.first] = triggerBits->accept( it.second );
          if (htTriggerDecision_[it.first]){
            htTriggerSummary = 1;
          }
      }
    }
    for( auto& it : metTriggerIndex_ ) {
      if( it.second != -1 ) {
          metTriggerDecision_[it.first] = triggerBits->accept( it.second );
          if (metTriggerDecision_[it.first]){
            metTriggerSummary = 1;
          }
      }
    }
  }
 
 
  
 
  intEventProperties["triggerSummaryHT"] = htTriggerSummary;
  intEventProperties["triggerSummaryMET"] = metTriggerSummary;
  
  
  if (leptonFlavor1 == "Ele" && leptonFlavor2 == "Ele") {
    intEventProperties["triggerSummary"] = eeTriggerSummary;
    fillTree<pat::Electron, pat::Electron>( "EE", (*electrons).at(leptonNr1), (*electrons).at(leptonNr2),*isoTracks,*looseElectrons,*looseMuons,*jets,*shiftedJetsJESUp,*shiftedJetsJESDown,*bJets35,*shiftedBJetsJESUp,*shiftedBJetsJESDown, met,MHT,vertices,Rho, intEventProperties, longIntEventProperties, floatEventProperties,tLorentzVectorEventProperties,isMC); 
  }
  else if (leptonFlavor1 == "Mu" && leptonFlavor2 == "Mu") {
    intEventProperties["triggerSummary"] = mmTriggerSummary;
    fillTree<pat::Muon, pat::Muon>( "MuMu", (*muons).at(leptonNr1), (*muons).at(leptonNr2),*isoTracks,*looseElectrons,*looseMuons,*jets,*shiftedJetsJESUp,*shiftedJetsJESDown,*bJets35,*shiftedBJetsJESUp,*shiftedBJetsJESDown, met,MHT,vertices,Rho, intEventProperties, longIntEventProperties, floatEventProperties,tLorentzVectorEventProperties,isMC); 
  }
  else if (leptonFlavor1 == "Ele" && leptonFlavor2 == "Mu") {
    intEventProperties["triggerSummary"] = emTriggerSummary;
    fillTree<pat::Electron, pat::Muon>( "EMu", (*electrons).at(leptonNr1), (*muons).at(leptonNr2),*isoTracks,*looseElectrons,*looseMuons,*jets,*shiftedJetsJESUp,*shiftedJetsJESDown,*bJets35,*shiftedBJetsJESUp,*shiftedBJetsJESDown, met,MHT,vertices,Rho, intEventProperties, longIntEventProperties, floatEventProperties,tLorentzVectorEventProperties,isMC); 
  }
  // Change ordering for Mu E events, in such a way that the electron is always the first lepton (required by some tools that select the lepton flavor)
  else if (leptonFlavor1 == "Mu" && leptonFlavor2 == "Ele") {
    intEventProperties["triggerSummary"] = emTriggerSummary;
    fillTree<pat::Electron, pat::Muon>( "EMu", (*electrons).at(leptonNr2), (*muons).at(leptonNr1),*isoTracks,*looseElectrons,*looseMuons,*jets,*shiftedJetsJESUp,*shiftedJetsJESDown,*bJets35,*shiftedBJetsJESUp,*shiftedBJetsJESDown, met,MHT,vertices,Rho, intEventProperties, longIntEventProperties, floatEventProperties,tLorentzVectorEventProperties,isMC); 
  }
    

  delete jecUnc;
  delete shiftedJetsJESUp;
  delete shiftedJetsJESDown;
  delete shiftedBJetsJESUp;
  delete shiftedBJetsJESDown;


}


template <class aT, class bT> void 
DiLeptonTreesFromMiniAOD::fillTree( const std::string &treeName, const aT& a, const bT& b,const std::vector<pat::IsolatedTrack>&isoTracks,const std::vector<pat::Electron>&looseElectrons,const std::vector<pat::Muon>&looseMuons,const std::vector<pat::Jet>&jets,const std::vector<pat::Jet>&shiftedJetsUp,const std::vector<pat::Jet>&shiftedJetsDown,const std::vector<pat::Jet>&bJets35,const std::vector<pat::Jet>&shiftedBJetsUp,const std::vector<pat::Jet>&shiftedBJetsDown, const pat::MET &patMet,const TLorentzVector &MHT,const edm::Handle<reco::VertexCollection> &vertices,const float &rho, const std::map<std::string, int> &intEventProperties, const std::map<std::string, unsigned long> &longIntEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties, const bool &isMC)
{

  for(const auto& it : intEventProperties){
    assert(intBranches_[treeName].find(it.first) != intBranches_[treeName].end());
    *(intBranches_[treeName][it.first]) = it.second;
  }
  for(const auto& it : longIntEventProperties){
    assert(longIntBranches_[treeName].find(it.first) != longIntBranches_[treeName].end());
    *(longIntBranches_[treeName][it.first]) = it.second;
  }
  for(const auto& it : floatEventProperties){
    assert(floatBranches_[treeName].find(it.first) != floatBranches_[treeName].end());
    *(floatBranches_[treeName][it.first]) = it.second;
  }
 for(const auto& it : tLorentzVectorEventProperties){
    assert(tLorentzVectorBranches_[treeName].find(it.first) != tLorentzVectorBranches_[treeName].end());
    *(tLorentzVectorBranches_[treeName][it.first]) = it.second;
  }   
  
  if(debug) std::cout << treeName << "- pts:"<< a.pt() << " " << b.pt();
  TLorentzVector aVec = getMomentum(a);//( a.px(), a.py(), a.pz(), a.energy() );
  TLorentzVector bVec = getMomentum(b); //( b.px(), b.py(), b.pz(), b.energy() );
  TLorentzVector met(patMet.px(), patMet.py(), patMet.pz(), patMet.energy());
  TLorentzVector uncorrectedMet; 
  uncorrectedMet.SetPtEtaPhiE(patMet.uncorPt(), 0, \
            patMet.uncorPhi(), patMet.uncorPt());
    
   

  TLorentzVector comb = aVec+bVec;
  TLorentzVector MHT2 = comb + MHT;
  
  float sum_mlb;
  float mlb_min;
  float mlb_max;
  
  sumMlb(aVec, bVec, jets, bJets35, sum_mlb, mlb_min, mlb_max);
  *(floatBranches_[treeName]["sumMlb"]) = sum_mlb; 
  *(floatBranches_[treeName]["min_mlb"]) = mlb_min; 
  *(floatBranches_[treeName]["max_mlb"]) = mlb_max; 
  
  sum_mlb = 0;
  mlb_min = 0;
  mlb_max = 0;
  sumMlb(aVec, bVec, shiftedJetsUp, shiftedBJetsUp, sum_mlb, mlb_min, mlb_max);
  *(floatBranches_[treeName]["sumMlbJESUp"]) = sum_mlb;
  
  sum_mlb = 0;
  mlb_min = 0;
  mlb_max = 0;
  sumMlb(aVec, bVec, shiftedJetsDown, shiftedBJetsDown, sum_mlb, mlb_min, mlb_max);
  *(floatBranches_[treeName]["sumMlbJESDown"]) = sum_mlb;

  
  double pa[3] = {aVec.M(),aVec.Px(),aVec.Py()};
  double pb[3] = {bVec.M(),bVec.Px(),bVec.Py()};
  
  double pmiss[3] = {0.,met.Px(),met.Py()};
  fctMT2_.set_mn(0.);
  fctMT2_.set_momenta(pa,pb,pmiss,*(longIntBranches_[treeName]["eventNr"]));
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
  //double absIso = 0.;
  
  std::vector < float > trackPts;
  
  for(auto track : isoTracks){
    if (not track.fromPV()) continue;
    if (not track.packedCandRef().isNonnull()) continue;
    if (abs(track.pdgId()) == 11 || abs(track.pdgId()) == 13){
      if (track.pt() < 5 or abs(track.eta()) > 2.4) continue;
      if (abs(track.dz())  > 0.1) continue;
      if (abs(track.dxy()) > 0.2) continue;
      if (track.pfIsolationDR03().chargedHadronIso() > 5) continue;
      if (track.pfIsolationDR03().chargedHadronIso()/track.pt() < 0.2) continue;
    }else{
      if (track.pt() < 10 or abs(track.eta()) > 2.4) continue;
      if (abs(track.dz())  > 0.1) continue;
      if (abs(track.dxy()) > 0.2) continue;
      if (track.pfIsolationDR03().chargedHadronIso() > 5) continue;
      if (track.pfIsolationDR03().chargedHadronIso()/track.pt() < 0.2) continue;
    }
    nIsoTracks++;
    trackPts.push_back(track.pt());
    
  }
  *(intBranches_[treeName]["nIsoTracks"]) = nIsoTracks;
  if (nIsoTracks > 0) *(floatBranches_[treeName]["ptTrack3"]) = *std::max_element(std::begin(trackPts),std::end(trackPts));
  else *(floatBranches_[treeName]["ptTrack3"]) = 0.;
  
          
  
  //~ std::pair<double, double> pZeta = calcPZeta(aVec, bVec, met);
  *(floatBranches_[treeName]["chargeProduct"]) = a.charge()*b.charge();
  *(floatBranches_[treeName]["mll"]) = comb.M();
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
  if (isMC){ 
    *(floatBranches_[treeName]["leptonFullSimScaleFactor1"]) = fctLeptonFullSimScaleFactors_(a, a.pt(), a.eta() );
    *(floatBranches_[treeName]["leptonFullSimScaleFactor2"]) = fctLeptonFullSimScaleFactors_(b, b.pt(), b.eta() );
  }
  else{  
    *(floatBranches_[treeName]["leptonFullSimScaleFactor1"]) = 1.;
    *(floatBranches_[treeName]["leptonFullSimScaleFactor2"]) = 1.;
  }
  *(floatBranches_[treeName]["miniIsoEffArea1"]) = getIso(a,"miniIsoEA");
  *(floatBranches_[treeName]["miniIsoEffArea2"]) = getIso(b,"miniIsoEA");
  *(floatBranches_[treeName]["mt1"]) = transverseMass(aVec, met);
  *(floatBranches_[treeName]["mt2"]) = transverseMass(bVec, met);
  //~ *(floatBranches_[treeName]["fakeWeight1"]) = fakeRates_(a);
  //~ *(floatBranches_[treeName]["fakeWeight2"]) = fakeRates_(b);
  *(floatBranches_[treeName]["deltaPhi"]) = aVec.DeltaPhi( bVec );
  *(floatBranches_[treeName]["deltaR"]) = aVec.DeltaR( bVec );
  //~ *(floatBranches_[treeName]["angle3D"]) = aVec.Angle( bVec.Vect() );   

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
  


  trees_[treeName]->Fill();
}

void 
DiLeptonTreesFromMiniAOD::sumMlb(TLorentzVector &lepton1, TLorentzVector &lepton2, const std::vector<pat::Jet> &jets, const std::vector<pat::Jet> &bjets, float &result_sum_mlb, float &result_mlb_min, float &result_mlb_max){
  float mlb_min = 1.E6;
  float mlb_max = 1.E6;
  
  float temp_mlb;
  
  TLorentzVector jet (0.,0.,0.,0.);
  TLorentzVector jet1 (0.,0.,0.,0.);
  TLorentzVector lepton (0.,0.,0.,0.);
  
  TLorentzVector leptList[2] = {lepton1,lepton2};
  
  int lmin = -1;
  
  std::vector< pat::Jet > jet1Coll;
  std::vector< pat::Jet > jet2Coll;
  
  // Calculate sum Mlb
  if (bjets.size() >= 1) jet1Coll = bjets;
  else jet1Coll = jets;
  
  if (bjets.size() >= 2) jet2Coll = bjets;
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
        if (bjets.size() == 1 && jet.DeltaR (jet1) < 0.1) continue;
        if ( (bjets.size() == 0 || bjets.size() >= 2) &&  jet == jet1) continue;  
        lepton.SetPxPyPzE(leptList[il].Px(),leptList[il].Py(),leptList[il].Pz(),leptList[il].Energy());
        
        temp_mlb = (jet+lepton).M();
        if (temp_mlb < mlb_max) {
          mlb_max = temp_mlb;
        }
      }
    }
  }
   
  if (mlb_min < 1.E6 and mlb_max < 1.E6) result_sum_mlb = mlb_min + mlb_max; 
  else result_sum_mlb = 1.E-6;
  
  if (mlb_min == 1.E6) mlb_min = 0.0001;
  if (mlb_max == 1.E6) mlb_max = 0.0001;
  
  result_mlb_min = mlb_min; 
  result_mlb_max = mlb_max; 
  
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
  // needs to be produced first before using it
  //corr = e.userFloat("ecalTrkEnergyPostCorr") / e.energy();
  
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
  return fctIsolation_(e,method)* 1./getMomentum(e).Pt();
}

float DiLeptonTreesFromMiniAOD::getIso(const  pat::Muon &mu, const std::string &method)
{
  return fctIsolation_(mu,method)* 1./getMomentum(mu).Pt();
}

float DiLeptonTreesFromMiniAOD::getIso(const  pat::PackedCandidate &track, const std::string &method)
{
  return fctIsolation_(track,method);
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
  eeNewLumiBlock_=true;
  emNewLumiBlock_=true;
  mmNewLumiBlock_=true;
  htNewLumiBlock_=true;
}

void 
DiLeptonTreesFromMiniAOD::endLuminosityBlock(edm::LuminosityBlock const& iEvent, edm::EventSetup const&)
{
  
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
