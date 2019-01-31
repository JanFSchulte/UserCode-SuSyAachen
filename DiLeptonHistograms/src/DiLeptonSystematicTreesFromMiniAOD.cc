// -*- C++ -*-
//
// Package:    Histograms
// Class:      DiLeptonSystematicTreesFromMiniAOD
// 
/**\class DiLeptonSystematicTreesFromMiniAOD DiLeptonSystematicTreesFromMiniAOD.cc brot/DiLeptonSystematicTreesFromMiniAOD/src/DiLeptonSystematicTreesFromMiniAOD.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  matthias edelhoff
//         Created:  Tue Oct 27 13:50:40 CET 2009
// $Id: DiLeptonSystematicTreesFromMiniAOD.cc,v 1.31 2012/09/17 17:38:58 sprenger Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"

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
#include <SuSyAachen/DiLeptonHistograms/interface/PdgIdFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/VertexWeightFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/IsolationFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/BTagCalibrationStandalone.h>
#include <SuSyAachen/DiLeptonHistograms/interface/BTagEffMapFunctor.h>
//~ #include <SuSyAachen/DiLeptonHistograms/interface/LeptonFullSimScaleFactorMapFunctor.h>
//~ #include <SuSyAachen/DiLeptonHistograms/interface/LeptonFastSimScaleFactorMapFunctor.h>


//ROOT
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace std;

//
// class decleration
//

class DiLeptonSystematicTreesFromMiniAOD : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit DiLeptonSystematicTreesFromMiniAOD(const edm::ParameterSet&);
  ~DiLeptonSystematicTreesFromMiniAOD();

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
  template<class aT, class bT> void fillTree( const std::string &treeName, const aT &a, const bT &b, const std::vector<pat::PackedCandidate>&pfCands,const std::vector<pat::Electron>&looseElectrons,const std::vector<pat::Muon>&looseMuons,const std::vector<pat::Jet>&jets,const std::vector<pat::Jet>&shiftedJetsJESUp,const std::vector<pat::Jet>&shiftedJetsJESDown,const std::vector<pat::Jet>&bJets35,const std::vector<pat::Jet>&shiftedJetsBJESUp,const std::vector<pat::Jet>&shiftedBJetsJESDown,  const pat::MET &patMet, const TLorentzVector &MHT, const edm::Handle<reco::VertexCollection> &vertices, const float &rho, const std::map<std::string, long> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties, const TLorentzVector &JESChangeUp, const TLorentzVector &JESChangeDown);
  void sumMlb(TLorentzVector &lepton1, TLorentzVector &lepton2, const std::vector<pat::Jet> &jets, const std::vector<pat::Jet> &bjets, float &result_sum_mlb, float &result_mlb_min, float &result_mlb_max);
  const TLorentzVector getMomentum(const  pat::Electron &e);
  const TLorentzVector getMomentum(const  pat::Muon &mu);
  float getIso(const  pat::Electron &e, const std::string &method);
  float getIso(const  pat::Muon &mu, const std::string &method);
  float getIso(const  pat::PackedCandidate &track, const std::string &method);
  float transverseMass(const TLorentzVector& p, const TLorentzVector& met);
  edm::EDGetTokenT< std::vector< pat::Electron > >      electronToken_;
  edm::EDGetTokenT< std::vector< pat::Electron > >      looseElectronToken_;
  edm::EDGetTokenT< std::vector< pat::Muon > >        muonToken_;
  edm::EDGetTokenT< std::vector< pat::Muon > >        looseMuonToken_;
  edm::EDGetTokenT< std::vector< pat::Jet > >       jetToken_;
  edm::EDGetTokenT< std::vector< reco::GenJet > >     genJetToken_;
  edm::EDGetTokenT< std::vector< pat::Jet > >       bJetToken_;
  edm::EDGetTokenT< std::vector< pat::Jet > >       bJet35Token_;
  edm::EDGetTokenT< std::vector< pat::MET > >         metToken_;
  edm::EDGetTokenT<reco::VertexCollection>          vertexToken_;
  edm::EDGetTokenT< std::vector<pat::PackedCandidate>  >  pfCandToken_;
  edm::EDGetTokenT< std::vector< reco::GenParticle > >    genParticleToken_;
  edm::EDGetTokenT<GenEventInfoProduct>           genEventInfoToken_;
  edm::EDGetTokenT<LHEEventProduct>             LHEEventToken_;
  edm::EDGetTokenT<double>                  rhoToken_;
  
  edm::EDGetTokenT<edm::TriggerResults>           metFilterToken_;


  std::map<double, double> electronCorrections_;
  //data
  std::map<std::string, TTree*> trees_;  
  std::map<std::string, std::map< std::string, float*> > floatBranches_; 
  std::map<std::string, std::map< std::string, int*> > intBranches_; 
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
  
  BTagCalibration fctBTagCalibFastSim_;
  BTagCalibrationReader fctBTagCalibReaderFastSim_;
  BTagCalibrationReader fctBTagCalibReaderFastSimUp_;
  BTagCalibrationReader fctBTagCalibReaderFastSimDown_;
  
  //~ LeptonFullSimScaleFactorMapFunctor fctLeptonFullSimScaleFactors_;
  //~ LeptonFastSimScaleFactorMapFunctor fctLeptonFastSimScaleFactors_;

  VertexWeightFunctor fctVtxWeight_;
  VertexWeightFunctor fctVtxWeightUp_;
  VertexWeightFunctor fctVtxWeightDown_;
  IsolationFunctor fctIsolation_;
  PdgIdFunctor getPdgId_;
  MT2Functor fctMT2_;
 

  bool debug;
  bool metUncert_; 
  bool storeMetFilters_;
  std::vector<std::string> metFilterNames_;
  
};

// constructors and destructor
DiLeptonSystematicTreesFromMiniAOD::DiLeptonSystematicTreesFromMiniAOD(const edm::ParameterSet& iConfig):
  electronToken_    (consumes< std::vector< pat::Electron > >     (iConfig.getParameter<edm::InputTag>("electrons"))),
  looseElectronToken_ (consumes< std::vector< pat::Electron > >     (iConfig.getParameter<edm::InputTag>("looseElectrons"))),
  muonToken_      (consumes< std::vector< pat::Muon > >     (iConfig.getParameter<edm::InputTag>("muons"))),
  looseMuonToken_   (consumes< std::vector< pat::Muon > >     (iConfig.getParameter<edm::InputTag>("looseMuons"))),
  jetToken_       (consumes< std::vector< pat::Jet > >      (iConfig.getParameter<edm::InputTag>("jets"))),
  genJetToken_      (consumes< std::vector< reco::GenJet  > >   (iConfig.getParameter<edm::InputTag>("genJets"))),
  bJetToken_      (consumes< std::vector< pat::Jet > >      (iConfig.getParameter<edm::InputTag>("bJets"))),
  bJet35Token_      (consumes< std::vector< pat::Jet > >      (iConfig.getParameter<edm::InputTag>("bJets35"))),
  metToken_       (consumes< std::vector< pat::MET > >      (iConfig.getParameter<edm::InputTag>("met"))),
  vertexToken_      (consumes<reco::VertexCollection>       (iConfig.getParameter<edm::InputTag>("vertices"))),
  pfCandToken_      (consumes< std::vector<pat::PackedCandidate>  > (iConfig.getParameter<edm::InputTag>("pfCands"))),
  genParticleToken_   (consumes< std::vector< reco::GenParticle > > (iConfig.getParameter<edm::InputTag>("genParticles"))),
  genEventInfoToken_  (consumes<GenEventInfoProduct>          (iConfig.getParameter<edm::InputTag>("pdfInfo"))),
  LHEEventToken_    (consumes<LHEEventProduct>            (iConfig.getParameter<edm::InputTag>("LHEInfo"))),
  rhoToken_       (consumes<double>               (iConfig.getParameter<edm::InputTag>("rho"))),
  
  metFilterToken_         (consumes<edm::TriggerResults>    (edm::InputTag("TriggerResults",""))),

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
  
  fctBTagCalibFastSim_    (iConfig.getParameter<edm::ParameterSet>("BTagCalibration").getParameter<std::string>("CSVFastSimTagger"),iConfig.getParameter<edm::ParameterSet>("BTagCalibration").getParameter<std::string>("CSVFastSimFileName") ),
  fctBTagCalibReaderFastSim_    (&fctBTagCalibFastSim_,BTagEntry::OP_MEDIUM,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_fastSim"),"central" ),
  fctBTagCalibReaderFastSimUp_    (&fctBTagCalibFastSim_,BTagEntry::OP_MEDIUM,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_fastSim"),"up" ),
  fctBTagCalibReaderFastSimDown_    (&fctBTagCalibFastSim_,BTagEntry::OP_MEDIUM,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_fastSim"),"down" ),
  
  //~ fctLeptonFullSimScaleFactors_ (iConfig.getParameter<edm::ParameterSet>("LeptonFullSimScaleFactors") ),
  //~ fctLeptonFastSimScaleFactors_ (iConfig.getParameter<edm::ParameterSet>("LeptonFastSimScaleFactors") ),
  
  fctVtxWeight_    (iConfig.getParameter<edm::ParameterSet>("vertexWeights") ,consumesCollector()),
  fctVtxWeightUp_    (iConfig.getParameter<edm::ParameterSet>("vertexWeightsUp") ,consumesCollector()),
  fctVtxWeightDown_    (iConfig.getParameter<edm::ParameterSet>("vertexWeightsDown") ,consumesCollector()),
  fctIsolation_  (iConfig.getParameter<edm::ParameterSet>("isolationDefinitions"),consumesCollector()), 
  getPdgId_( iConfig.getParameter< edm::ParameterSet>("pdgIdDefinition"),consumesCollector() ),
  metFilterNames_ (iConfig.getUntrackedParameter< std::vector <std::string> >("metFilterNames"))
 
  
{
  usesResource("TFileService");
  debug = false;
  metUncert_ = iConfig.getUntrackedParameter<bool>("doMETUncert");
  storeMetFilters_ = iConfig.getUntrackedParameter<bool>("storeMetFilters");  
  // read config
  
  consumes<std::vector< PileupSummaryInfo > >(edm::InputTag("slimmedAddPileupInfo"));
  
  consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
  consumes<edm::TriggerResults>(edm::InputTag("TriggerResults",""));
  
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
  initFloatBranch( "weight" );
  initFloatBranch( "weightUp" );
  initFloatBranch( "weightDown" );
  initFloatBranch( "bTagWeight" );
  initFloatBranch( "bTagWeightErrHeavy" );
  initFloatBranch( "bTagWeightErrLight" );
  //~ initFloatBranch( "leptonFullSimScaleFactor1" );
  //~ initFloatBranch( "leptonFullSimScaleFactor2" );
  //~ initFloatBranch( "leptonFastSimScaleFactor1" );
  //~ initFloatBranch( "leptonFastSimScaleFactor2" );
  initFloatBranch( "genPtDiSbottom" );
  initFloatBranch( "chargeProduct" );
  initFloatBranch( "mll" );
  initFloatBranch( "charge1" );
  initFloatBranch( "charge2" );
  
  initFloatBranch( "mGluino" );
  initFloatBranch( "mLightSquarks" );
  initFloatBranch( "mSbottom" );
  initFloatBranch( "mStop" );
  initFloatBranch( "mSlepton" );
  initFloatBranch( "mSneutrino" );
  initFloatBranch( "mNeutralino1" );
  initFloatBranch( "mNeutralino2" );
  initFloatBranch( "mNeutralino3" );
  initFloatBranch( "mNeutralino4" );
  initFloatBranch( "mChargino1" );
  initFloatBranch( "mChargino2" );
  initFloatBranch( "mGravitino" );
  
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
  initFloatBranch( "deltaPhi" );
  initFloatBranch( "deltaR" );
  initFloatBranch( "min_mlb" );
  initFloatBranch( "max_mlb" );
  initFloatBranch( "sumMlb" );
  initFloatBranch( "sumMlbJESUp" );
  initFloatBranch( "sumMlbJESDown" );
  initFloatBranch( "MT2" );
  initFloatBranch( "MT2JESUp" );
  initFloatBranch( "MT2JESDown" );
  initFloatBranch( "deltaPhiJetMet1" );
  initFloatBranch( "deltaPhiJetMet2" );
  initFloatBranch( "deltaPhiJetMet1JESUp" );
  initFloatBranch( "deltaPhiJetMet2JESUp" );
  initFloatBranch( "deltaPhiJetMet1JESDown" );
  initFloatBranch( "deltaPhiJetMet2JESDown" ); 
  initFloatBranch( "ht" );
  initFloatBranch( "htJESUp" );
  initFloatBranch( "htJESDown" );
  initFloatBranch( "mht" );
  initFloatBranch( "met" );
  initFloatBranch( "uncorrectedMet" ); 
  initFloatBranch( "genMet" );
  initFloatBranch( "genMetNeutrinos" );
  initFloatBranch( "genMetPromptNeutrinos" );
  initFloatBranch( "caloMet" );
  initFloatBranch( "metJESUp" );
  initFloatBranch( "metJESDown" );
  initIntBranch( "nJets" );
  initIntBranch( "nGenJets" );
  initIntBranch( "nBadMuonJets" );
  initIntBranch( "nUnmatchedJets" );
  initIntBranch( "nISRJets" );
  initFloatBranch( "ISRCorrection" );
  initFloatBranch( "ISRUncertainty" );
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
  
  initFloatBranch( "scaleWeight1" );
  initFloatBranch( "scaleWeight2" );
  initFloatBranch( "scaleWeight3" );
  initFloatBranch( "scaleWeight4" );
  initFloatBranch( "scaleWeight5" );
  initFloatBranch( "scaleWeight6" );
  initFloatBranch( "scaleWeight7" );
  initFloatBranch( "scaleWeight8" );
  
  if (storeMetFilters_){
    
    for (const auto& n : metFilterNames_){
      initIntBranch( n.c_str() );
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
DiLeptonSystematicTreesFromMiniAOD::initFloatBranch(const std::string &name)
{
  for( const auto& it : trees_){
    if(debug) std::cout << it.first <<" - "<< name << std::endl;
    floatBranches_[it.first][name] = new float;
    it.second->Branch(name.c_str(), floatBranches_[it.first][name], (name+"/F").c_str());
  }
}

void 
DiLeptonSystematicTreesFromMiniAOD::initIntBranch(const std::string &name)
{
  for( const auto& it : trees_){
    if(debug) std::cout << it.first <<" - "<< name << std::endl;
    intBranches_[it.first][name] = new int;
    it.second->Branch(name.c_str(), intBranches_[it.first][name], (name+"/I").c_str());
  }
}

void 
DiLeptonSystematicTreesFromMiniAOD::initTLorentzVectorBranch(const std::string &name)
{
  for( const auto& it : trees_){
    if(debug) std::cout << it.first <<" - "<< name << std::endl;
    tLorentzVectorBranches_[it.first][name] = new TLorentzVector;
    it.second->Branch(name.c_str(), "TLorentzVector" ,&tLorentzVectorBranches_[it.first][name]);
  }
}

DiLeptonSystematicTreesFromMiniAOD::~DiLeptonSystematicTreesFromMiniAOD()
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


}


// member functions
// ------------ method called to for each event  ------------
void
DiLeptonSystematicTreesFromMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  
  //~ edm::Handle< std::vector< pat::Jet > > jets;
  iEvent.getByToken(jetToken_, jets);
  
  edm::Handle< std::vector< reco::GenJet > > genJets;
  iEvent.getByToken(genJetToken_, genJets);

  edm::Handle< std::vector< pat::Jet > > bJets;
  iEvent.getByToken(bJetToken_, bJets);
  
  edm::Handle< std::vector< pat::Jet > > bJets35;
  iEvent.getByToken(bJet35Token_, bJets35);
  
  edm::Handle< std::vector< pat::MET > > mets;
  iEvent.getByToken(metToken_, mets);


  edm::Handle< std::vector< reco::GenParticle > > genParticles;
  iEvent.getByToken(genParticleToken_, genParticles);


  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);

  getPdgId_.loadGenParticles(iEvent);
  fctIsolation_.init(iEvent);
  std::map<std::string, long> intEventProperties;
  std::map<std::string, float> floatEventProperties;
  std::map<std::string, TLorentzVector> tLorentzVectorEventProperties;


  edm::Handle<GenEventInfoProduct> genInfoProduct;
  iEvent.getByToken(genEventInfoToken_, genInfoProduct);  
  if (genInfoProduct.isValid()){
    if ((*genInfoProduct).weight() < 0.0){
    
      floatEventProperties["genWeight"] = -1;
    }
    else{
      floatEventProperties["genWeight"] = 1;    
    }
  }
  else{
 
  floatEventProperties["genWeight"] = 1; 
  
  } 
  
  if (genParticles.isValid()) {
      edm::Handle<GenEventInfoProduct> GenEventInfoHandle;
      iEvent.getByLabel("generator", GenEventInfoHandle);
      auto weightsize = GenEventInfoHandle->weights().size();
      if (weightsize < 2) {   // for old scans
         edm::Handle<LHEEventProduct> LHEEventProductHandle;
         iEvent.getByToken(LHEEventToken_, LHEEventProductHandle);
         if (LHEEventProductHandle.isValid()) {               
             
             floatEventProperties["scaleWeight1"] = LHEEventProductHandle->weights()[1].wgt/LHEEventProductHandle->originalXWGTUP();
             floatEventProperties["scaleWeight2"] = LHEEventProductHandle->weights()[2].wgt/LHEEventProductHandle->originalXWGTUP();
             floatEventProperties["scaleWeight3"] = LHEEventProductHandle->weights()[3].wgt/LHEEventProductHandle->originalXWGTUP();
             floatEventProperties["scaleWeight4"] = LHEEventProductHandle->weights()[4].wgt/LHEEventProductHandle->originalXWGTUP();
             floatEventProperties["scaleWeight5"] = LHEEventProductHandle->weights()[5].wgt/LHEEventProductHandle->originalXWGTUP();
             floatEventProperties["scaleWeight6"] = LHEEventProductHandle->weights()[6].wgt/LHEEventProductHandle->originalXWGTUP();
             floatEventProperties["scaleWeight7"] = LHEEventProductHandle->weights()[7].wgt/LHEEventProductHandle->originalXWGTUP();
             floatEventProperties["scaleWeight8"] = LHEEventProductHandle->weights()[8].wgt/LHEEventProductHandle->originalXWGTUP();
             
         }
      } else { // for SMS scans
          floatEventProperties["scaleWeight1"] = GenEventInfoHandle->weights()[2]/GenEventInfoHandle->weights()[1];
          floatEventProperties["scaleWeight2"] = GenEventInfoHandle->weights()[3]/GenEventInfoHandle->weights()[1];
          floatEventProperties["scaleWeight3"] = GenEventInfoHandle->weights()[4]/GenEventInfoHandle->weights()[1];
          floatEventProperties["scaleWeight4"] = GenEventInfoHandle->weights()[5]/GenEventInfoHandle->weights()[1];
          floatEventProperties["scaleWeight5"] = GenEventInfoHandle->weights()[6]/GenEventInfoHandle->weights()[1];
          floatEventProperties["scaleWeight6"] = GenEventInfoHandle->weights()[7]/GenEventInfoHandle->weights()[1];
          floatEventProperties["scaleWeight7"] = GenEventInfoHandle->weights()[8]/GenEventInfoHandle->weights()[1];
          floatEventProperties["scaleWeight8"] = GenEventInfoHandle->weights()[9]/GenEventInfoHandle->weights()[1];

      }
   }
   else {
     floatEventProperties["scaleWeight1"] = 1.;
     floatEventProperties["scaleWeight2"] = 1.;
     floatEventProperties["scaleWeight3"] = 1.;
     floatEventProperties["scaleWeight4"] = 1.;
     floatEventProperties["scaleWeight5"] = 1.;
     floatEventProperties["scaleWeight6"] = 1.;
     floatEventProperties["scaleWeight7"] = 1.;
     floatEventProperties["scaleWeight8"] = 1.; 
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
        floatEventProperties["genParticleHT"]  = ht;
    }
  
  else floatEventProperties["genParticleHT"]  = -999.;
        

  
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
  
  
  edm::Handle<edm::TriggerResults> metFilterBits;
  iEvent.getByToken(metFilterToken_, metFilterBits);
  const edm::TriggerNames &allFilterNames = iEvent.triggerNames(*metFilterBits);
  
  
  if (storeMetFilters_){ 
          
    for (std::string const &filterName : metFilterNames_){
      const unsigned index = allFilterNames.triggerIndex(filterName);
      if (index >= allFilterNames.size()) std::cerr << "MET filter " << filterName << "not found!" << std::endl;
      if (metFilterBits->accept(index)) intEventProperties[filterName] = 1;
      else intEventProperties[filterName] = 0;

    }  
  }
 else{ 
    
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
  uncorrectedMetVector.SetPtEtaPhiE(met.uncorPt(), 0, met.uncorPhi(), met.uncorPt());
  
  floatEventProperties["met"] = metVector.Pt();
  tLorentzVectorEventProperties["vMet"] = metVector; 
  
  //~ tLorentzVectorEventProperties["vMetUncorrected"] = uncorrectedMetVector;
  floatEventProperties["uncorrectedMet"] = uncorrectedMetVector.Pt();
  
  floatEventProperties["caloMet"] = met.caloMETPt();
  

  pat::METCollection const& metsForUncert = *mets;  
  floatEventProperties["met"] =  metsForUncert[0].pt(); 
  if (metUncert_){

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
  TLorentzVector genMetNeutrinosVector(0.,0.,0.,0.);
  TLorentzVector genMetPromptNeutrinosVector(0.,0.,0.,0.);
  TLorentzVector vGenParticle(0.,0.,0.,0.);
  TLorentzVector tempParticle(0.,0.,0.,0.);
  TLorentzVector genNu1Vector(0.,0.,0.,0.);
  TLorentzVector genNu2Vector(0.,0.,0.,0.);
  int nuNum = 1;
  
  TLorentzVector genSbottom1(0.,0.,0.,0.);
  TLorentzVector genSbottom2(0.,0.,0.,0.);
  TLorentzVector genDiSbottom(0.,0.,0.,0.);
  
  floatEventProperties["mGluino"] = -1.;
  floatEventProperties["mLightSquarks"] = -1.;
  floatEventProperties["mStop"] = -1.;
  floatEventProperties["mSbottom"] = -1.;
  floatEventProperties["mSlepton"] = -1.;
  floatEventProperties["mSneutrino"] = -1.;
  floatEventProperties["mNeutralino1"] = -1.;
  floatEventProperties["mNeutralino2"] = -1.;
  floatEventProperties["mNeutralino3"] = -1.;
  floatEventProperties["mNeutralino4"] = -1.;
  floatEventProperties["mChargino1"] = -1.;
  floatEventProperties["mChargino2"] = -1.;
  floatEventProperties["mGravitino"] = -1.;
  //~ floatEventProperties["genPtDiSbottom"] = -1.;
  
  if (genParticles.isValid()){
  
  genMetVector.SetPxPyPzE(mets->front().genMET()->px(),mets->front().genMET()->py(),mets->front().genMET()->pz(),mets->front().genMET()->energy());
    
  for (std::vector<reco::GenParticle>::const_iterator itGenParticle = genParticles->begin(); itGenParticle != genParticles->end(); itGenParticle++) { 
    if (abs((*itGenParticle).pdgId())== 12 || abs((*itGenParticle).pdgId())== 14 || abs((*itGenParticle).pdgId())== 16){
        tempParticle.SetPxPyPzE((*itGenParticle).px(),(*itGenParticle).py(),0, 0);
        genMetNeutrinosVector = genMetNeutrinosVector + tempParticle;

        //if (abs((*itGenParticle).mother()->pdgId()) == 24){
        if ((*itGenParticle).isPromptFinalState()){
          genMetPromptNeutrinosVector = genMetPromptNeutrinosVector + tempParticle;
          if (nuNum == 1){
            genNu1Vector.SetPxPyPzE((*itGenParticle).px(),(*itGenParticle).py(),(*itGenParticle).pz(), (*itGenParticle).energy());
            nuNum += 1;
          }else if (nuNum == 2){
            genNu2Vector.SetPxPyPzE((*itGenParticle).px(),(*itGenParticle).py(),(*itGenParticle).pz(), (*itGenParticle).energy());
            nuNum += 1;
          }
        }
      }
    
    
    if (abs((*itGenParticle).pdgId())== 1000001 || abs((*itGenParticle).pdgId())== 1000002 || abs((*itGenParticle).pdgId())== 1000003 || abs((*itGenParticle).pdgId())== 1000004 || abs((*itGenParticle).pdgId())== 2000001 || abs((*itGenParticle).pdgId())== 2000002 || abs((*itGenParticle).pdgId())== 2000003 || abs((*itGenParticle).pdgId())== 2000004){
      floatEventProperties["mLightSquarks"] = (*itGenParticle).mass();
    }
    if (abs((*itGenParticle).pdgId())== 1000005 || abs((*itGenParticle).pdgId())== 2000005){
      floatEventProperties["mSbottom"] = (*itGenParticle).mass();
      //~ if ((*itGenParticle).pdgId()== 1000005){
        //~ genSbottom1.SetPxPyPzE((*itGenParticle).px(), (*itGenParticle).py(), (*itGenParticle).pz(), (*itGenParticle).energy());       
      //~ }
      //~ else if ((*itGenParticle).pdgId()== -1000005){
        //~ genSbottom2.SetPxPyPzE((*itGenParticle).px(), (*itGenParticle).py(), (*itGenParticle).pz(), (*itGenParticle).energy());
      //~ }
    }
    if (abs((*itGenParticle).pdgId())== 1000006 || abs((*itGenParticle).pdgId())== 2000006){
      floatEventProperties["mStop"] = (*itGenParticle).mass();
    }
    if (abs((*itGenParticle).pdgId())== 1000011 || abs((*itGenParticle).pdgId())== 1000013 || abs((*itGenParticle).pdgId())== 1000015 || abs((*itGenParticle).pdgId())== 2000011 || abs((*itGenParticle).pdgId())== 2000013 || abs((*itGenParticle).pdgId())== 2000015){
      floatEventProperties["mSlepton"] = (*itGenParticle).mass();
    }
    if (abs((*itGenParticle).pdgId())== 1000012 || abs((*itGenParticle).pdgId())== 1000014 || abs((*itGenParticle).pdgId())== 1000016){
      floatEventProperties["mSneutrino"] = (*itGenParticle).mass();
    }
    if (abs((*itGenParticle).pdgId())== 1000021 ){
      floatEventProperties["mGluino"] = (*itGenParticle).mass();
    }
    if (abs((*itGenParticle).pdgId())== 1000022 ){
      floatEventProperties["mNeutralino1"] = (*itGenParticle).mass();
    }
    if (abs((*itGenParticle).pdgId())== 1000023 ){
      floatEventProperties["mNeutralino2"] = (*itGenParticle).mass();
    }
    if (abs((*itGenParticle).pdgId())== 1000024 ){
      floatEventProperties["mChargino1"] = (*itGenParticle).mass();
    }
    if (abs((*itGenParticle).pdgId())== 1000025 ){
      floatEventProperties["mNeutralino3"] = (*itGenParticle).mass();
    }
    if (abs((*itGenParticle).pdgId())== 1000035 ){
      floatEventProperties["mNeutralino4"] = (*itGenParticle).mass();
    }
    if (abs((*itGenParticle).pdgId())== 1000037 ){
      floatEventProperties["mChargino2"] = (*itGenParticle).mass();
    }
    if (abs((*itGenParticle).pdgId())== 1000039 ){
      floatEventProperties["mGravitino"] = (*itGenParticle).mass();
    }


  }

  }
  
  floatEventProperties["genMet"] = genMetVector.Pt();
  floatEventProperties["genMetNeutrinos"] = genMetNeutrinosVector.Pt();
  floatEventProperties["genMetPromptNeutrinos"] = genMetPromptNeutrinosVector.Pt();
  tLorentzVectorEventProperties["vGenMet"] = genMetVector;  
  tLorentzVectorEventProperties["vGenMetNeutrinos"] = genMetNeutrinosVector;  
  tLorentzVectorEventProperties["vGenMetPromptNeutrinos"] = genMetPromptNeutrinosVector;    
  
  tLorentzVectorEventProperties["genNu1"] = genNu1Vector;    
  tLorentzVectorEventProperties["genNu2"] = genNu2Vector;   
  
  //~ genDiSbottom = genSbottom1 + genSbottom2;
  //~ if (genDiSbottom.Pt() > 1){
    //~ floatEventProperties["genPtDiSbottom"] = genDiSbottom.Pt();
  //~ }
  
  //~ if (genDiSbottom.Pt() < 400) {
  //~ ISRUncertainty = 0.;
  //~ }
  //~ else {
  //~ if (genDiSbottom.Pt() < 600){
    //~ ISRUncertainty = 0.15;
  //~ }
  //~ else {
    //~ ISRUncertainty = 0.3;
  //~ }
  //~ }
  //~ floatEventProperties["ISRUncertainty"] = ISRUncertainty;


  TLorentzVector MHT;
  TLorentzVector tempMHT;
  
  floatEventProperties["ht"] = 0.0;
  floatEventProperties["htJESUp"] = 0.0;
  floatEventProperties["htJESDown"] = 0.0;
  floatEventProperties["mht"] = 0.0;

  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl); 
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);


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
  int nUnmatchedJets=0;
  int nISRJets=0;
  int nGenJets=0;  
  int nBadMuonJets=0;
  
  for(std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end() ; ++it){
  if ((*it).pt() >=35.0 && fabs((*it).eta())<2.4){
    nJets++;
    
    bool matchedJet = false;
    for (std::vector<reco::GenParticle>::const_iterator itGenParticle = genParticles->begin(); itGenParticle != genParticles->end(); itGenParticle++) { 
      if (matchedJet) break;
      if ( abs((*itGenParticle).pdgId()) > 5 || (*itGenParticle).status() != 23) continue;
      int momid =  abs((*itGenParticle).mother()->pdgId());
      if(!(momid == 6 || momid == 23 || momid == 24 || momid == 25 || momid > 1e6)) continue;
      //check against daughter in case of hard initial splitting
      for (size_t idau(0); idau < (*itGenParticle).numberOfDaughters(); idau++) {
        float dR = deltaR((*it), (*itGenParticle).daughter(idau)->p4() );
        if (dR < 0.3) {
          matchedJet = true;
          break;
        }
      }
    }
    if (!matchedJet) nISRJets++;
    
    tempMHT.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy()); 
    MHT = MHT + tempMHT;
    floatEventProperties["ht"] += (*it).pt();
    
    if ((*it).pt() >=200.0 && (*it).muonEnergyFraction() > 0.5 && fabs(tempMHT.DeltaPhi( metVector )) > M_PI - 0.4 ){
      nBadMuonJets++;
    }
  }
  
  if ((*it).pt() >=20.0 && fabs((*it).eta())<2.5){
    
    if (genParticles.isValid()){
      
      bool matched = false;
      TLorentzVector jetVector( (*it).px(), (*it).py(), (*it).pz(), (*it).energy() );
      
      for(std::vector<reco::GenJet >::const_iterator itGenJet = genJets->begin(); itGenJet != genJets->end() ; ++itGenJet){
        TLorentzVector genJetVector( (*itGenJet).px(), (*itGenJet).py(), (*itGenJet).pz(), (*itGenJet).energy() );
        
        if ( jetVector.DeltaR(genJetVector) < 0.3) {
          matched = true;
          break;
        }
      }
      if (matched == false && (*it).chargedHadronEnergy()/(*it).energy() < 0.1) {
        nUnmatchedJets+=1;
      }
    }
  }
  
  }
  intEventProperties["nJets"] = nJets;
  intEventProperties["nUnmatchedJets"] = nUnmatchedJets;
  intEventProperties["nISRJets"] = nISRJets;
  intEventProperties["nBadMuonJets"] = nBadMuonJets;
  
  floatEventProperties["genJet1pt"] = -1.0;
  floatEventProperties["genJet2pt"] = -1.0;
  floatEventProperties["genJet3pt"] = -1.0;
  floatEventProperties["genJet4pt"] = -1.0;
  
  floatEventProperties["ISRCorrection"] = 1.;
  floatEventProperties["ISRUncertainty"] = 0.;
  
  if (nISRJets==1) {
    floatEventProperties["ISRCorrection"] = 0.920;
    floatEventProperties["ISRUncertainty"] = 0.040;   
  }
  if (nISRJets==2) {
    floatEventProperties["ISRCorrection"] = 0.821;
    floatEventProperties["ISRUncertainty"] = 0.090;   
  }
  if (nISRJets==3) {
    floatEventProperties["ISRCorrection"] = 0.715;
    floatEventProperties["ISRUncertainty"] = 0.143;   
  }
  if (nISRJets==4) {
    floatEventProperties["ISRCorrection"] = 0.662;
    floatEventProperties["ISRUncertainty"] = 0.169;   
  }
  if (nISRJets==5) {
    floatEventProperties["ISRCorrection"] = 0.561;
    floatEventProperties["ISRUncertainty"] = 0.219;   
  }
  if (nISRJets>=6) {
    floatEventProperties["ISRCorrection"] = 0.511;
    floatEventProperties["ISRUncertainty"] = 0.244;   
  }

  float genHT = 0.;
  std::vector<reco::GenJet> genJetsCleaned;
  
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
  
  TLorentzVector metVectorJESUp = metVector + JESChangeUp;
  TLorentzVector metVectorJESDown = metVector + JESChangeDown;
  floatEventProperties["metJESUp"] = metVectorJESUp.Pt();
  floatEventProperties["metJESDown"] = metVectorJESDown.Pt();

  int nJetsJESUp=0;
  TLorentzVector jet1VectorJESUp(0.,0.,0.,0.);
  TLorentzVector jet2VectorJESUp(0.,0.,0.,0.);
  floatEventProperties["deltaPhiJetMet1JESUp"] = -9999.;
  floatEventProperties["deltaPhiJetMet2JESUp"] = -9999.;
  
  for(std::vector<pat::Jet>::const_iterator it = shiftedJetsJESUp->begin(); it != shiftedJetsJESUp->end() ; ++it){
  if ((*it).pt() >=35.0 && fabs((*it).eta())<2.4){
          nJetsJESUp++;
    floatEventProperties["htJESUp"] += (*it).pt();
  }
  }
  intEventProperties["nShiftedJetsJESUp"] = nJetsJESUp;
  
  if (nJetsJESUp > 0)
    jet1VectorJESUp.SetPxPyPzE(shiftedJetsJESUp->at(0).px(),shiftedJetsJESUp->at(0).py(),shiftedJetsJESUp->at(0).pz(),shiftedJetsJESUp->at(0).energy()); 
  floatEventProperties["deltaPhiJetMet1JESUp"] = fabs(jet1VectorJESUp.DeltaPhi( metVectorJESUp ));
  if (nJetsJESUp > 1)
    jet2VectorJESUp.SetPxPyPzE(shiftedJetsJESUp->at(1).px(),shiftedJetsJESUp->at(1).py(),shiftedJetsJESUp->at(1).pz(),shiftedJetsJESUp->at(1).energy());
  floatEventProperties["deltaPhiJetMet2JESUp"] = fabs(jet2VectorJESUp.DeltaPhi( metVectorJESUp ));
    
  
  int nJetsJESDown=0;
  TLorentzVector jet1VectorJESDown(0.,0.,0.,0.);
  TLorentzVector jet2VectorJESDown(0.,0.,0.,0.);
  floatEventProperties["deltaPhiJetMet1JESDown"] = -9999.;
  floatEventProperties["deltaPhiJetMet2JESDown"] = -9999.;
  for(std::vector<pat::Jet>::const_iterator it = shiftedJetsJESDown->begin(); it != shiftedJetsJESDown->end() ; ++it){
  if ((*it).pt() >=35.0 && fabs((*it).eta())<2.4){  
          nJetsJESDown++;
    floatEventProperties["htJESDown"] += (*it).pt();
  }
  }
  intEventProperties["nShiftedJetsJESDown"] = nJetsJESDown;
  
  if (nJetsJESDown > 0)
    jet1VectorJESDown.SetPxPyPzE(shiftedJetsJESDown->at(0).px(),shiftedJetsJESDown->at(0).py(),shiftedJetsJESDown->at(0).pz(),shiftedJetsJESDown->at(0).energy()); 
  floatEventProperties["deltaPhiJetMet1JESDown"] = fabs(jet1VectorJESDown.DeltaPhi( metVectorJESDown ));
  if (nJetsJESDown > 1)
    jet2VectorJESDown.SetPxPyPzE(shiftedJetsJESDown->at(1).px(),shiftedJetsJESDown->at(1).py(),shiftedJetsJESDown->at(1).pz(),shiftedJetsJESDown->at(1).energy());
  floatEventProperties["deltaPhiJetMet2JESDown"] = fabs(jet2VectorJESDown.DeltaPhi( metVectorJESDown ));
  

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
    float temp_pt;
    
    for(std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end() ; ++it){
      if ((*it).pt() >=35.0 && fabs((*it).eta())<2.4){
        
        int jetFlavor;
        float eff;
        float SF, SF_up, SF_down;
        float SF_FullSim, SF_FullSim_up, SF_FullSim_down;
        float SF_FastSim, SF_FastSim_up, SF_FastSim_down;
        float SF_err;
        temp_pt = std::min(599.,(*it).pt());
        
        jetFlavor = (*it).hadronFlavour();
        if (jetFlavor == 5)
        {
          SF_FullSim = fctBTagCalibReaderFullSimBJets_.eval(BTagEntry::FLAV_B, (*it).eta(), temp_pt);
          SF_FullSim_up = fctBTagCalibReaderFullSimBJetsUp_.eval(BTagEntry::FLAV_B, (*it).eta(), temp_pt);
          SF_FullSim_down = fctBTagCalibReaderFullSimBJetsDown_.eval(BTagEntry::FLAV_B, (*it).eta(), temp_pt);
          
          SF_FastSim = fctBTagCalibReaderFastSim_.eval(BTagEntry::FLAV_B, (*it).eta(), temp_pt);
          SF_FastSim_up = fctBTagCalibReaderFastSimUp_.eval(BTagEntry::FLAV_B, (*it).eta(), temp_pt);
          SF_FastSim_down = fctBTagCalibReaderFastSimDown_.eval(BTagEntry::FLAV_B, (*it).eta(), temp_pt);
          
         
        }
        else if (jetFlavor == 4)
        {
          SF_FullSim = fctBTagCalibReaderFullSimCJets_.eval(BTagEntry::FLAV_C, (*it).eta(), temp_pt);
          SF_FullSim_up = fctBTagCalibReaderFullSimCJetsUp_.eval(BTagEntry::FLAV_C, (*it).eta(), temp_pt);
          SF_FullSim_down = fctBTagCalibReaderFullSimCJetsDown_.eval(BTagEntry::FLAV_C, (*it).eta(), temp_pt);
          
          SF_FastSim = fctBTagCalibReaderFastSim_.eval(BTagEntry::FLAV_C, (*it).eta(), temp_pt);
          SF_FastSim_up = fctBTagCalibReaderFastSimUp_.eval(BTagEntry::FLAV_C, (*it).eta(), temp_pt);
          SF_FastSim_down = fctBTagCalibReaderFastSimDown_.eval(BTagEntry::FLAV_C, (*it).eta(), temp_pt);

        }
        else
        {
          SF_FullSim = fctBTagCalibReaderFullSimLightJets_.eval(BTagEntry::FLAV_UDSG, (*it).eta(), temp_pt);
          SF_FullSim_up = fctBTagCalibReaderFullSimLightJetsUp_.eval(BTagEntry::FLAV_UDSG, (*it).eta(), temp_pt);
          SF_FullSim_down = fctBTagCalibReaderFullSimLightJetsDown_.eval(BTagEntry::FLAV_UDSG, (*it).eta(), temp_pt);
          
          SF_FastSim = fctBTagCalibReaderFastSim_.eval(BTagEntry::FLAV_UDSG, (*it).eta(), temp_pt);
          SF_FastSim_up = fctBTagCalibReaderFastSimUp_.eval(BTagEntry::FLAV_UDSG, (*it).eta(), temp_pt);
          SF_FastSim_down = fctBTagCalibReaderFastSimDown_.eval(BTagEntry::FLAV_UDSG, (*it).eta(), temp_pt);

        }
        
        SF = SF_FullSim * SF_FastSim;
        SF_up = SF_FullSim_up * SF_FastSim_up;
        SF_down = SF_FullSim_down * SF_FastSim_down;
        
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
  
 
  
  if (leptonFlavor1 == "Ele" && leptonFlavor2 == "Ele") {
    fillTree<pat::Electron, pat::Electron>( "EE", (*electrons).at(leptonNr1), (*electrons).at(leptonNr2),*pfCands,*looseElectrons,*looseMuons,*jets,*shiftedJetsJESUp,*shiftedJetsJESDown,*bJets35,*shiftedBJetsJESUp,*shiftedBJetsJESDown, met,MHT,vertices,Rho, intEventProperties, floatEventProperties,tLorentzVectorEventProperties,JESChangeUp,JESChangeDown); 
  }
  else if (leptonFlavor1 == "Mu" && leptonFlavor2 == "Mu") {
    fillTree<pat::Muon, pat::Muon>( "MuMu", (*muons).at(leptonNr1), (*muons).at(leptonNr2),*pfCands,*looseElectrons,*looseMuons,*jets,*shiftedJetsJESUp,*shiftedJetsJESDown,*bJets35,*shiftedBJetsJESUp,*shiftedBJetsJESDown, met,MHT,vertices,Rho, intEventProperties, floatEventProperties,tLorentzVectorEventProperties,JESChangeUp,JESChangeDown); 
  }
  else if (leptonFlavor1 == "Ele" && leptonFlavor2 == "Mu") {
    fillTree<pat::Electron, pat::Muon>( "EMu", (*electrons).at(leptonNr1), (*muons).at(leptonNr2),*pfCands,*looseElectrons,*looseMuons,*jets,*shiftedJetsJESUp,*shiftedJetsJESDown,*bJets35,*shiftedBJetsJESUp,*shiftedBJetsJESDown, met,MHT,vertices,Rho, intEventProperties, floatEventProperties,tLorentzVectorEventProperties,JESChangeUp,JESChangeDown); 
  }
  // Change ordering for Mu E events, in such a way that the electron is always the first lepton (required by some tools that select the lepton flavor)
  else if (leptonFlavor1 == "Mu" && leptonFlavor2 == "Ele") {
    fillTree<pat::Electron, pat::Muon>( "EMu", (*electrons).at(leptonNr2), (*muons).at(leptonNr1),*pfCands,*looseElectrons,*looseMuons,*jets,*shiftedJetsJESUp,*shiftedJetsJESDown,*bJets35,*shiftedBJetsJESUp,*shiftedBJetsJESDown, met,MHT,vertices,Rho, intEventProperties, floatEventProperties,tLorentzVectorEventProperties,JESChangeUp,JESChangeDown); 
  }
  
  delete jecUnc;
  delete shiftedJetsJESUp;
  delete shiftedJetsJESDown;
  delete shiftedBJetsJESUp;
  delete shiftedBJetsJESDown;

}


template <class aT, class bT> void 
DiLeptonSystematicTreesFromMiniAOD::fillTree( const std::string &treeName, const aT& a, const bT& b,const std::vector<pat::PackedCandidate>&pfCands,const std::vector<pat::Electron>&looseElectrons,const std::vector<pat::Muon>&looseMuons,const std::vector<pat::Jet>&jets,const std::vector<pat::Jet>&shiftedJetsUp,const std::vector<pat::Jet>&shiftedJetsDown,const std::vector<pat::Jet>&bJets35,const std::vector<pat::Jet>&shiftedBJetsUp,const std::vector<pat::Jet>&shiftedBJetsDown, const pat::MET &patMet,const TLorentzVector &MHT,const edm::Handle<reco::VertexCollection> &vertices,const float &rho, const std::map<std::string, long> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties, const TLorentzVector &JESChangeUp, const TLorentzVector &JESChangeDown)
{

  for(const auto& it : intEventProperties){
    assert(intBranches_[treeName].find(it.first) != intBranches_[treeName].end());
    *(intBranches_[treeName][it.first]) = it.second;
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
            
  TLorentzVector metJESUp = met+JESChangeUp;
  TLorentzVector metJESDown = met+JESChangeDown;  

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
  fctMT2_.set_momenta(pa,pb,pmiss);
  *(floatBranches_[treeName]["MT2"]) = static_cast<float>(fctMT2_.get_mt2()); 
  
  double pmissJESUp[3] = {0.,metJESUp.Px(),metJESUp.Py()};
  fctMT2_.set_momenta(pa,pb,pmissJESUp);
  *(floatBranches_[treeName]["MT2JESUp"]) = static_cast<float>(fctMT2_.get_mt2()); 
  
  double pmissJESDown[3] = {0.,metJESDown.Px(),metJESDown.Py()};
  fctMT2_.set_momenta(pa,pb,pmissJESDown);
  *(floatBranches_[treeName]["MT2JESDown"]) = static_cast<float>(fctMT2_.get_mt2()); 
  
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
    
    if ( abs((*it).pdgId()) == 11 || abs((*it).pdgId()) == 13 || ( abs((*it).pdgId()) == 211  && absIso < 0.1 * (*it).pt() ) ){
      nIsoTracks++;
      trackPts.push_back((*it).pt());
    }
    
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
  //~ *(floatBranches_[treeName]["leptonFullSimScaleFactor1"]) = fctLeptonFullSimScaleFactors_(a, a.pt(), a.eta() );
  //~ *(floatBranches_[treeName]["leptonFullSimScaleFactor2"]) = fctLeptonFullSimScaleFactors_(b, b.pt(), b.eta() );
  //~ *(floatBranches_[treeName]["leptonFastSimScaleFactor1"]) = fctLeptonFastSimScaleFactors_(a, a.pt(), fabs(a.eta()));
  //~ *(floatBranches_[treeName]["leptonFastSimScaleFactor2"]) = fctLeptonFastSimScaleFactors_(b, b.pt(), fabs(b.eta()) );
  *(floatBranches_[treeName]["miniIsoEffArea1"]) = getIso(a,"miniIsoEA");
  *(floatBranches_[treeName]["miniIsoEffArea2"]) = getIso(b,"miniIsoEA");
  *(floatBranches_[treeName]["mt1"]) = transverseMass(aVec, met);
  *(floatBranches_[treeName]["mt2"]) = transverseMass(bVec, met);
  *(floatBranches_[treeName]["deltaPhi"]) = aVec.DeltaPhi( bVec );
  *(floatBranches_[treeName]["deltaR"]) = aVec.DeltaR( bVec );

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
DiLeptonSystematicTreesFromMiniAOD::sumMlb(TLorentzVector &lepton1, TLorentzVector &lepton2, const std::vector<pat::Jet> &jets, const std::vector<pat::Jet> &bjets, float &result_sum_mlb, float &result_mlb_min, float &result_mlb_max){
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


const TLorentzVector DiLeptonSystematicTreesFromMiniAOD::getMomentum(const  pat::Electron &e)
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

const TLorentzVector DiLeptonSystematicTreesFromMiniAOD::getMomentum(const  pat::Muon &mu)
{
  const TLorentzVector result = TLorentzVector(mu.px(), mu.py(), mu.pz(), mu.energy());
  return result;
}

float DiLeptonSystematicTreesFromMiniAOD::getIso(const  pat::Electron &e, const std::string &method)
{
  return fctIsolation_(e,method)* 1./e.pt();
}

float DiLeptonSystematicTreesFromMiniAOD::getIso(const  pat::Muon &mu, const std::string &method)
{
  return fctIsolation_(mu,method)* 1./mu.pt();
}

float DiLeptonSystematicTreesFromMiniAOD::getIso(const  pat::PackedCandidate &track, const std::string &method)
{
  return fctIsolation_(track,method);
}


float DiLeptonSystematicTreesFromMiniAOD::transverseMass(const TLorentzVector& p, const TLorentzVector& met)
{
  reco::Candidate::LorentzVector otherMet(met.Px(),met.Py(),met.Pz(),met.E());
  reco::Candidate::LorentzVector leptonT(p.Px(),p.Py(),0.,p.E()*sin(p.Theta()));
  reco::Candidate::LorentzVector sumT=leptonT+otherMet;

  return std::sqrt(sumT.M2());
}

// ------------ Method called once each job just before starting event loop  ------------
void 
DiLeptonSystematicTreesFromMiniAOD::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiLeptonSystematicTreesFromMiniAOD::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiLeptonSystematicTreesFromMiniAOD);
