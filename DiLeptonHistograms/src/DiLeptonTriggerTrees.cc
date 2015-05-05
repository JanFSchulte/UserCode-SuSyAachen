// -*- C++ -*-
//
// Package:    Histograms
// Class:      DiLeptonTriggerTrees
// 
/**\class DiLeptonTriggerTrees DiLeptonTriggerTrees.cc brot/DiLeptonTriggerTrees/src/DiLeptonTriggerTrees.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  matthias edelhoff
//         Created:  Tue Oct 27 13:50:40 CET 2009
// $Id: DiLeptonTriggerTrees.cc,v 1.31 2012/09/17 17:38:58 sprenger Exp $
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
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/Lepton.h>

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/MET.h>

#include <DataFormats/Provenance/interface/EventID.h>

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include <SuSyAachen/DiLeptonHistograms/interface/WeightFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/PdgIdFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/VertexWeightFunctor.h>
#include <SuSyAachen/TagAndProbeTreeWriter/interface/IsolationFunctor.h>

//ROOT
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace std;

//
// class decleration
//

class DiLeptonTriggerTrees : public edm::EDAnalyzer {
public:
  explicit DiLeptonTriggerTrees(const edm::ParameterSet&);
  ~DiLeptonTriggerTrees();

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
  template <class aT, class bT> void makeCombinations( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT >&b, const std::vector<reco::PFCandidate>&pfCands, const edm::Event &ev, const pat::MET &patMet, const TLorentzVector &MHT, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties,const edm::Handle<trigger::TriggerEvent> triggerEvent);
  template <class aT> void makeCombinations( const std::string &treeName, const std::vector<aT> &a , const std::vector<reco::PFCandidate>&pfCands,  const edm::Event &ev, const pat::MET &patMet, const TLorentzVector &MHT, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties,const edm::Handle<trigger::TriggerEvent> triggerEvent);
  template<class aT, class bT> void fillTree( const std::string &treeName, const aT &a, const bT &b, const std::vector<reco::PFCandidate>&pfCands, const pat::MET &patMet, const TLorentzVector &MHT,const edm::Handle<trigger::TriggerEvent> triggerEvent);
  int getLeptonPdgId( const reco::GenParticle &p);
  int getMotherPdgId( const reco::GenParticle &p);
  std::pair<double, double> calcPZeta(const TLorentzVector& p1,const TLorentzVector& p2, const TLorentzVector& met);
  void fillPdfUncert(const edm::Handle< std::vector<double> >& weightHandle, const std::string& pdfIdentifier, const std::string& treeName);

  const TLorentzVector getMomentum(const  pat::Electron &e);
  const TLorentzVector getMomentum(const  pat::Muon &mu);
  const TLorentzVector getMomentum(const  pat::Tau &tau);
  float getPhiAtECAL(const  pat::Electron &e, const std::vector<reco::PFCandidate>&pfCands);
  float getPhiAtECAL(const  pat::Muon &mu, const std::vector<reco::PFCandidate>&pfCands);
  float getPhiAtECAL(const  pat::Tau &tau, const std::vector<reco::PFCandidate>&pfCands);
  float getIso(const  pat::Electron &e, const std::string &method);
  float getIso(const  pat::Muon &mu, const std::string &method);
  float getIso(const  pat::Tau &tau, const std::string &method);
  float getDeltaB(const  pat::Electron &e);
  float getDeltaB(const  pat::Muon &mu);
  float getDeltaB(const  pat::Tau &tau);
  float transverseMass(const TLorentzVector& p, const TLorentzVector& met);
  std::string convertInputTag(const edm::InputTag tag);

  edm::InputTag eTag_;
  edm::InputTag muTag_;
  edm::InputTag tauTag_;
  edm::InputTag jetTag_;
  edm::InputTag jet2Tag_;
  edm::InputTag bJetTag_;
  edm::InputTag metTag_;
  edm::InputTag type1metTag_;
  edm::InputTag tcmetTag_;
  edm::InputTag caloMetTag_;
  edm::InputTag genMetTrueTag_;
  edm::InputTag vertexTag_;
  edm::InputTag pfCandTag_;
  edm::InputTag genParticleTag_;
  edm::InputTag triggerSummary_;
  std::vector<edm::ParameterSet> susyVars_;
  std::vector<edm::InputTag> pdfs_;
  std::string tauId_;


  std::map<double, double> electronCorrections_;
  //data
  std::map<std::string, TTree*> trees_;  
  std::map<std::string, std::map< std::string, float*> > floatBranches_; 
  std::map<std::string, std::map< std::string, int*> > intBranches_; 
  std::map<std::string, std::map< std::string, TLorentzVector*> > tLorentzVectorBranches_;

  edm::Handle< std::vector< pat::Jet > > jets_;

  WeightFunctor fakeRates_;
  WeightFunctor efficiencies_;
  VertexWeightFunctor fctVtxWeight_;
  IsolationFunctor fctIsolation_;
  PdgIdFunctor getPdgId_;

  bool debug;
  bool useJets2_;
  bool useTaus_;
};

// constructors and destructor
DiLeptonTriggerTrees::DiLeptonTriggerTrees(const edm::ParameterSet& iConfig):
  fctVtxWeight_    (iConfig.getParameter<edm::ParameterSet>("vertexWeights") ),
  fctIsolation_  (iConfig.getParameter<edm::ParameterSet>("isolationDefinitions")),
  getPdgId_( iConfig.getParameter< edm::ParameterSet>("pdgIdDefinition") )
{
  debug = false;
  useTaus_ = iConfig.existsAs<edm::InputTag>("taus");
  useJets2_ = iConfig.existsAs<edm::InputTag>("jets2");
  
  // read config
  eTag_ = iConfig.getParameter<edm::InputTag>("electrons");
  muTag_ = iConfig.getParameter<edm::InputTag>("muons");
  if(useTaus_)tauTag_ = iConfig.getParameter<edm::InputTag>("taus");
  jetTag_ = iConfig.getParameter<edm::InputTag>("jets");
  if(useJets2_) jet2Tag_ = iConfig.getParameter<edm::InputTag>("jets2");
  bJetTag_ = iConfig.getParameter<edm::InputTag>("bJets");
  metTag_ = iConfig.getParameter<edm::InputTag>("met");
  type1metTag_ = iConfig.getParameter<edm::InputTag>("type1met");
  tcmetTag_ = iConfig.getParameter<edm::InputTag>("tcmet");
  caloMetTag_ = iConfig.getParameter<edm::InputTag>("caloMet");
  genMetTrueTag_ = iConfig.getParameter<edm::InputTag>("genMetTrue");
  vertexTag_ = iConfig.getParameter<edm::InputTag>("vertices");
  pfCandTag_ = iConfig.getParameter<edm::InputTag>("pfCands");
  susyVars_ = iConfig.getParameter< std::vector<edm::ParameterSet> >("susyVars");
  genParticleTag_ = iConfig.getParameter<edm::InputTag>("genParticles");
  pdfs_ = iConfig.getParameter<std::vector<edm::InputTag> > ("pdfWeightTags");
  triggerSummary_ = iConfig.getParameter<edm::InputTag> ("triggerSummaryTag");

  tauId_ = iConfig.getParameter<std::string >("tauId");
  fakeRates_.SetSource(iConfig,"fakeRates");// TODO use these and add mcInfo flag to choose right rates...
  efficiencies_.SetSource(iConfig,"efficiencies");// TODO use these and add mcInfo flag to choose right rates...
  
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
  if(useTaus_){
    trees_["ETau"] = file->make<TTree>("ETauDileptonTree", "ETau DileponTree");
    trees_["MuTau"] = file->make<TTree>("MuTauDileptonTree", "MuTau DileponTree");
    trees_["TauTau"] = file->make<TTree>("TauTauDileptonTree", "TauTau DileponTree");
  }
  initFloatBranch( "weight" );
  initFloatBranch( "chargeProduct" );
  initFloatBranch( "charge1" );
  initFloatBranch( "charge2" );
  initTLorentzVectorBranch( "p4" );
  initTLorentzVectorBranch( "p4Gen" );
  initTLorentzVectorBranch( "vMet" );
  initIntBranch( "matchesSingleElectron1" );
  initIntBranch( "matchesSingleElectron2" );
  initIntBranch( "matchesSingleMuon1" );
  initIntBranch( "matchesSingleMuon2" );
  initIntBranch( "matchesDoubleElectron1" );
  initIntBranch( "matchesDoubleElectron2" );
  initIntBranch( "matchesDoubleMuonLeading1" );
  initIntBranch( "matchesDoubleMuonLeading2" );
  initIntBranch( "matchesDoubleMuonLeadingTk1" );
  initIntBranch( "matchesDoubleMuonLeadingTk2" );
  initIntBranch( "matchesEMuLeading1" );
  initIntBranch( "matchesEMuLeading2" );
  initIntBranch( "matchesMuELeading1" );
  initIntBranch( "matchesMuELeading2" );
  initIntBranch( "matchesDoubleMuonTrailing1" );
  initIntBranch( "matchesDoubleMuonTrailing2" );
  initIntBranch( "matchesDoubleMuonTrailingTk1" );
  initIntBranch( "matchesDoubleMuonTrailingTk2" );
  initIntBranch( "matchesEMuTrailing1" );
  initIntBranch( "matchesEMuTrailing2" );
  initIntBranch( "matchesMuETrailing1" );
  initIntBranch( "matchesMuETrailing2" );
  initFloatBranch( "pt1" );
  initFloatBranch( "pt2" );
  initFloatBranch( "eta1" );
  initFloatBranch( "eta2" );
  initFloatBranch( "isoEffArea1" );
  initFloatBranch( "isoEffArea2" );
  initFloatBranch( "isoPF1" );
  initFloatBranch( "isoPF2" );
  initFloatBranch( "dB1" );
  initFloatBranch( "dB2" );
  initFloatBranch( "mt1" );
  initFloatBranch( "mt2" );
  initFloatBranch( "eff1" );
  initFloatBranch( "eff2" );
  initFloatBranch( "fakeWeight1" );
  initFloatBranch( "fakeWeight2" );
  initFloatBranch( "deltaPhi" );
  initFloatBranch( "deltaR" );
  initFloatBranch( "jzb" );
  initFloatBranch( "ht" );
  initFloatBranch( "htJESUp" );
  initFloatBranch( "htJESDown" );
  initFloatBranch( "mht" );
  initFloatBranch( "met" );
  initFloatBranch( "metJESUp" );
  initFloatBranch( "metJESDown" );
  initFloatBranch( "type1Met" );
  initFloatBranch( "genMetTrue" );
  initIntBranch( "nJets" );
  initIntBranch( "nBJets" );
  initIntBranch( "nShiftedJetsJESUp" );
  initIntBranch( "nShiftedJetsJESDown" );
  initIntBranch( "nVertices" );
  initIntBranch( "nLightLeptons" );
  initFloatBranch( "jet1pt" );
  initFloatBranch( "jet2pt" );
  initFloatBranch( "jet3pt" );
  initFloatBranch( "jet4pt" );
  initIntBranch( "runNr" );
  initIntBranch( "lumiSec" );
  initIntBranch( "eventNr" );
  initIntBranch( "pdgId1" );
  initIntBranch( "pdgId2" );
  initIntBranch( "matched" );
  initIntBranch( "motherPdgId1" );
  initIntBranch( "motherPdgId2" );
  if(useJets2_) {
    initFloatBranch( "ht2" );
    initIntBranch( "nJets2" );    
  }
  for ( std::vector<edm::ParameterSet>::iterator susyVar_i = susyVars_.begin(); susyVar_i != susyVars_.end(); ++susyVar_i ) {
    edm::InputTag var = susyVar_i->getParameter<edm::InputTag>( "var" );
    std::string type = susyVar_i->getParameter<std::string>( "type" );
    if(debug) std::cout << var << " of type " << type << std::endl;
    if (type=="int") initIntBranch( convertInputTag(var) );
    else if (type=="float") initFloatBranch( convertInputTag(var) );
    else throw cms::Exception("Unrecognized type") << 
      "Unknown type " << type << " for variable" << var << " found\n";
  }
  for ( std::vector<edm::InputTag>::iterator pdf_i = pdfs_.begin(); pdf_i != pdfs_.end(); ++pdf_i ) {
     std::string pdfIdentifier = (*pdf_i).instance();
     std::string up = "Up";
     std::string down = "Down";
     initFloatBranch( pdfIdentifier );
     initFloatBranch( pdfIdentifier+up );
     initFloatBranch( pdfIdentifier+down );
  }
}

void 
DiLeptonTriggerTrees::initTLorentzVectorBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    tLorentzVectorBranches_[(*it).first][name] = new TLorentzVector;
    (*it).second->Branch(name.c_str(), "TLorentzVector" ,&tLorentzVectorBranches_[(*it).first][name]);
  }
}

void 
DiLeptonTriggerTrees::initFloatBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    floatBranches_[(*it).first][name] = new float;
    (*it).second->Branch(name.c_str(), floatBranches_[(*it).first][name], (name+"/F").c_str());
  }
}

void 
DiLeptonTriggerTrees::initIntBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    intBranches_[(*it).first][name] = new int;
    (*it).second->Branch(name.c_str(), intBranches_[(*it).first][name], (name+"/I").c_str());
  }
}

DiLeptonTriggerTrees::~DiLeptonTriggerTrees()
{ 
  for( std::map<std::string, std::map< std::string, float*> >::const_iterator it = floatBranches_.begin();
       it != floatBranches_.end(); ++it){
    for( std::map< std::string, float*>::const_iterator it2 = (*it).second.begin();
	 it2 != (*it).second.end(); ++it2){
      if(debug)std::cout << "deleting: " << (*it).first << " - "<< (*it2).first << std::endl;
      delete (*it2).second;
    }
  }
  for( std::map<std::string, std::map< std::string, int*> >::const_iterator it = intBranches_.begin();
       it != intBranches_.end(); ++it){
    for( std::map< std::string, int*>::const_iterator it2 = (*it).second.begin();
	 it2 != (*it).second.end(); ++it2){
      if(debug) std::cout << "deleting: " << (*it).first << " - "<< (*it2).first << std::endl;
      delete (*it2).second;
    }
  }


}

// member functions
// ------------ method called to for each event  ------------
void
DiLeptonTriggerTrees::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< std::vector< pat::Electron > > electrons;
  iEvent.getByLabel(eTag_, electrons);

  edm::Handle< std::vector< pat::Muon > > muons;
  iEvent.getByLabel(muTag_, muons);

  edm::Handle< std::vector<reco::PFCandidate>  > pfCands;
  iEvent.getByLabel(pfCandTag_, pfCands); 


  edm::Handle< std::vector< pat::Tau > > taus;
  if(useTaus_)
    iEvent.getByLabel(tauTag_, taus);
  
  //edm::Handle< std::vector< pat::Jet > > jets;
  iEvent.getByLabel(jetTag_, jets_);




  edm::Handle< std::vector< pat::Jet > > bJets;
  iEvent.getByLabel(bJetTag_, bJets);
  
  edm::Handle< std::vector< pat::MET > > mets;
  iEvent.getByLabel(metTag_, mets);


  edm::Handle< std::vector< pat::MET > > tcmets;
  iEvent.getByLabel(tcmetTag_, tcmets);


  edm::Handle< std::vector< pat::MET > > calomets;
  iEvent.getByLabel(caloMetTag_, calomets);


  edm::Handle< std::vector< reco::GenMET > > genMetTrues;
  iEvent.getByLabel(genMetTrueTag_, genMetTrues);

  edm::Handle< std::vector< reco::GenParticle > > genParticles;
  iEvent.getByLabel(genParticleTag_, genParticles);

  edm::Handle< std::vector< reco::PFMET > > type1mets;
  iEvent.getByLabel(type1metTag_, type1mets);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vertexTag_, vertices);

  getPdgId_.loadGenParticles(iEvent);
  fctIsolation_.init(iEvent);

  std::map<std::string, int> intEventProperties;
  std::map<std::string, float> floatEventProperties;
  std::map<std::string, TLorentzVector> tLorentzVectorEventProperties;

	
  intEventProperties["nVertices"] = vertices->size();



  intEventProperties["nBJets"] = bJets->size();
  intEventProperties["nLightLeptons"] = electrons->size() + muons->size();
  intEventProperties["runNr"] = iEvent.id().run();
  intEventProperties["lumiSec"] = iEvent.id().luminosityBlock();
  intEventProperties["eventNr"] = iEvent.id().event();



  pat::MET met = mets->front();
  TLorentzVector metVector(mets->front().px(), mets->front().py(), mets->front().pz(), mets->front().energy());
  //TLorentzVector uncorrectedMet;
  //  uncorrectedMet.SetPtEtaPhiE(mets->front().uncorrectedPt(), 0,	
  //  mets->front().uncorrectedPhi(), mets->front().uncorrectedPt());
  floatEventProperties["met"] = metVector.Pt();



  reco::MET type1met = type1mets->front();
  TLorentzVector type1metVector(type1mets->front().px(), type1mets->front().py(), type1mets->front().pz(), type1mets->front().energy());

  floatEventProperties["type1Met"] = type1metVector.Pt();

  pat::MET tcmet = tcmets->front();
  TLorentzVector tcmetVector(tcmets->front().px(), tcmets->front().py(), tcmets->front().pz(), tcmets->front().energy());



  if (!genMetTrues.isValid()){
	floatEventProperties["genMetTrue"] = -1.;

  }
  else{
	reco::GenMET genMetTrue = genMetTrues->front();
  	TLorentzVector genMetTrueVector(genMetTrues->front().px(), genMetTrues->front().py(), genMetTrues->front().pz(), genMetTrues->front().energy());
	floatEventProperties["genMetTrue"] = genMetTrueVector.Pt();
// 


  }





  TLorentzVector MHT;
  TLorentzVector tempMHT;
  TLorentzVector leadingJetMomentum;
  TLorentzVector subLeadingJetMomentum;
  floatEventProperties["ht"] = 0.0;
  floatEventProperties["htJESUp"] = 0.0;
  floatEventProperties["htJESDown"] = 0.0;
  floatEventProperties["mht"] = 0.0;

  edm::FileInPath fip("SuSyAachen/DiLeptonHistograms/data/GR_P_V40_AN1::All_Uncertainty_AK5PF.txt");
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(fip.fullPath());



  TLorentzVector changeUp(0,0,0,0);
  TLorentzVector changeDown(0,0,0,0);
  TLorentzVector tempChangeUp;
  TLorentzVector tempChangeDown;
  // loop over jets
  std::vector<pat::Jet> * shiftedJetsJESUp = new std::vector<pat::Jet>(); 
  std::vector<pat::Jet> * shiftedJetsJESDown = new std::vector<pat::Jet>(); 

  for (std::vector<pat::Jet>::const_iterator itJet = jets_->begin(); itJet != jets_->end(); itJet++) {

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
  for(std::vector<pat::Jet>::const_iterator it = jets_->begin(); it != jets_->end() ; ++it){
	if ((*it).pt() >=40.0){
		nJets++;
		tempMHT.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());	
		MHT = MHT + tempMHT;
		floatEventProperties["ht"] += (*it).pt();
		if((*it).pt() > leadingJetMomentum.Pt()){
		leadingJetMomentum.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());
		}
		if((*it).pt() < leadingJetMomentum.Pt() && ((*it).pt() > subLeadingJetMomentum.Pt())){
		subLeadingJetMomentum.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());
		}
	}
  }
  intEventProperties["nJets"] = nJets;











  int nJetsJESUp=0;
  for(std::vector<pat::Jet>::const_iterator it = shiftedJetsJESUp->begin(); it != shiftedJetsJESUp->end() ; ++it){
	if ((*it).pt() >=40.0){
        	nJetsJESUp++;
		floatEventProperties["htJESUp"] += (*it).pt();
	}
  }
  intEventProperties["nShiftedJetsJESUp"] = nJetsJESUp;
  int nJetsJESDown=0;
  for(std::vector<pat::Jet>::const_iterator it = shiftedJetsJESDown->begin(); it != shiftedJetsJESDown->end() ; ++it){
	if ((*it).pt() >=40.0){	
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
  if (jets_->size() > 0)
    floatEventProperties["jet1pt"] = jets_->at(0).pt();
  if (jets_->size() > 1)
    floatEventProperties["jet2pt"] = jets_->at(1).pt();
  if (jets_->size() > 2)
    floatEventProperties["jet3pt"] = jets_->at(2).pt();
  if (jets_->size() > 3)
    floatEventProperties["jet4pt"] = jets_->at(3).pt();

  // bjet pt

 



 
    	

  floatEventProperties["weight"] = fctVtxWeight_( iEvent );
  if(useJets2_) {
    edm::Handle< std::vector< pat::Jet > > jets2;
    iEvent.getByLabel(jet2Tag_, jets2);
    intEventProperties["nJets2"] = jets2->size();
    floatEventProperties["ht2"] = 0.0;
    for(std::vector<pat::Jet>::const_iterator it = jets2->begin(); it != jets2->end() ; ++it){
      floatEventProperties["ht2"] += (*it).pt();
    }
  }






  edm::Handle<trigger::TriggerEvent> triggerEvent;
  iEvent.getByLabel( triggerSummary_ , triggerEvent);







  makeCombinations< pat::Electron >("EE", *electrons,*pfCands, iEvent, met,MHT, intEventProperties, floatEventProperties,tLorentzVectorEventProperties,triggerEvent);
  makeCombinations< pat::Electron, pat::Muon >("EMu", *electrons, *muons,*pfCands, iEvent, met,MHT, intEventProperties, floatEventProperties,tLorentzVectorEventProperties,triggerEvent);
  makeCombinations< pat::Muon >("MuMu", *muons,*pfCands, iEvent, met,MHT, intEventProperties, floatEventProperties,tLorentzVectorEventProperties,triggerEvent);
  if(useTaus_){
    makeCombinations< pat::Electron, pat::Tau >("ETau", *electrons, *taus,*pfCands, iEvent, met,MHT, intEventProperties, floatEventProperties,tLorentzVectorEventProperties,triggerEvent);
    makeCombinations< pat::Muon, pat::Tau>("MuTau", *muons, *taus,*pfCands, iEvent, met,MHT, intEventProperties, floatEventProperties,tLorentzVectorEventProperties,triggerEvent);
    makeCombinations< pat::Tau >("TauTau", *taus,*pfCands, iEvent, met,MHT, intEventProperties, floatEventProperties,tLorentzVectorEventProperties,triggerEvent);
  }
  //  if( nMu != 2) std::cout << "-------! "<<nMu<<std::endl;

  delete jecUnc;
  delete shiftedJetsJESUp;
  delete shiftedJetsJESDown;


}

template <class aT, class bT> void 
DiLeptonTriggerTrees::makeCombinations ( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT> &b,const std::vector<reco::PFCandidate>&pfCands, const edm::Event &ev, const pat::MET &patMet, const TLorentzVector &MHT, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties,const edm::Handle<trigger::TriggerEvent> triggerEvent)
{
  TLorentzVector met(patMet.px(), patMet.py(), patMet.pz(), patMet.energy());
  TLorentzVector uncorrectedMet;
  uncorrectedMet.SetPtEtaPhiE(patMet.uncorrectedPt(), 0,\
                              patMet.uncorrectedPhi(), patMet.uncorrectedPt());
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


  for ( std::vector<edm::ParameterSet>::iterator susyVar_i = susyVars_.begin(); susyVar_i != susyVars_.end(); ++susyVar_i ) {
        edm::InputTag var = susyVar_i->getParameter<edm::InputTag>( "var" );
        std::string type = susyVar_i->getParameter<std::string>( "type" );
        edm::Handle< double > var_;
        ev.getByLabel(var, var_);
        if (type=="float") *(floatBranches_[treeName][convertInputTag(var)]) = float(*var_);
        else if (type=="int") *(intBranches_[treeName][convertInputTag(var)]) = int(*var_);
  }
  //std::cout << std::endl;
  for ( std::vector<edm::InputTag>::iterator pdf_i = pdfs_.begin(); pdf_i != pdfs_.end(); ++pdf_i ) {
     const std::string pdfIdentifier = (*pdf_i).instance();
     edm::Handle<std::vector<double> > weightHandle;
     ev.getByLabel((*pdf_i), weightHandle);
     fillPdfUncert(weightHandle,pdfIdentifier,treeName);
  }
  for( typename std::vector<aT>::const_iterator itA = a.begin(); itA != a.end(); ++itA){
    for( typename std::vector<bT>::const_iterator itB = b.begin(); itB != b.end(); ++itB){
//      std::cout << treeName <<": "<< fakeRates_(*itA) << std::endl;
//      weight *= fakeRates_();
      fillTree<aT,bT>( treeName, *itA, *itB,pfCands, patMet,MHT,triggerEvent); 
    }
  }
}

template <class aT> void 
DiLeptonTriggerTrees::makeCombinations ( const std::string &treeName, const std::vector<aT> &a,const std::vector<reco::PFCandidate>&pfCands, const edm::Event &ev, const pat::MET &patMet , const TLorentzVector &MHT, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties,const edm::Handle<trigger::TriggerEvent> triggerEvent)
{
  TLorentzVector met(patMet.px(), patMet.py(), patMet.pz(), patMet.energy());
  TLorentzVector uncorrectedMet;
  uncorrectedMet.SetPtEtaPhiE(patMet.uncorrectedPt(), 0,\
                              patMet.uncorrectedPhi(), patMet.uncorrectedPt());
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

  for ( std::vector<edm::ParameterSet>::iterator susyVar_i = susyVars_.begin(); susyVar_i != susyVars_.end(); ++susyVar_i ) {
        edm::InputTag var = susyVar_i->getParameter<edm::InputTag>( "var" );
        std::string type = susyVar_i->getParameter<std::string>( "type" );
        edm::Handle< double > var_;
        ev.getByLabel(var, var_);
        if (type=="int") *(intBranches_[treeName][convertInputTag(var)]) = int(*var_);
        else if (type=="float") *(floatBranches_[treeName][convertInputTag(var)]) = float(*var_);
  }
  for ( std::vector<edm::InputTag>::iterator pdf_i = pdfs_.begin(); pdf_i != pdfs_.end(); ++pdf_i ) {
     std::string pdfIdentifier = (*pdf_i).instance();
     edm::Handle<std::vector<double> > weightHandle;
     ev.getByLabel((*pdf_i), weightHandle);
     fillPdfUncert(weightHandle,pdfIdentifier,treeName);
  }

  for( typename std::vector<aT>::const_iterator itA = a.begin(); itA != a.end(); ++itA){
    for( typename std::vector<aT>::const_iterator itB = a.begin(); itB != itA; ++itB){
      //std::cout << treeName<<"("<<(*itA).pt()<<", "<<(*itA).eta() <<"): "<< fakeRates_(*itA) << std::endl;
//      weight *= fakeRates_();
      fillTree<aT, aT>( treeName, *itA, *itB,pfCands, patMet,MHT,triggerEvent); 
    }
  }
}

template <class aT, class bT> void 
DiLeptonTriggerTrees::fillTree( const std::string &treeName, const aT& a, const bT& b,const std::vector<reco::PFCandidate>&pfCands, const pat::MET &patMet, const TLorentzVector &MHT,edm::Handle<trigger::TriggerEvent> triggerEvent)
{
  if(debug) std::cout << treeName << "- pts:"<< a.pt() << " " << b.pt();
  TLorentzVector aVec = getMomentum(a);//( a.px(), a.py(), a.pz(), a.energy() );
  TLorentzVector bVec = getMomentum(b); //( b.px(), b.py(), b.pz(), b.energy() );
  TLorentzVector met(patMet.px(), patMet.py(), patMet.pz(), patMet.energy());
  TLorentzVector uncorrectedMet; 
  uncorrectedMet.SetPtEtaPhiE(patMet.uncorrectedPt(pat::MET::uncorrALL), 0, \
			      patMet.uncorrectedPhi(pat::MET::uncorrALL), patMet.uncorrectedPt(pat::MET::uncorrALL));  

  //  std::cout << "met: "<<met.Et()<< ", unCorr met: "<< uncorrectedMet.Et()
  //<< "=> "<< met.Et()* 1./uncorrectedMet.Et()<< " (xCheck: "<< patMet.corSumEt()*1./patMet.uncorrectedPt(pat::MET::uncorrALL) <<")"<<std::endl;


  TLorentzVector comb = aVec+bVec;
  TLorentzVector MHT2 = comb + MHT;
  *(floatBranches_[treeName]["chargeProduct"]) = a.charge()*b.charge();
  *(tLorentzVectorBranches_[treeName]["p4"]) = comb;
  *(floatBranches_[treeName]["mht"]) = MHT2.Pt();
  *(floatBranches_[treeName]["pt1"]) = aVec.Pt();
  *(floatBranches_[treeName]["pt2"]) = bVec.Pt();
  *(floatBranches_[treeName]["charge1"]) = a.charge();
  *(floatBranches_[treeName]["charge2"]) = b.charge();
  *(floatBranches_[treeName]["eta1"]) = aVec.Eta();
  *(floatBranches_[treeName]["eta2"]) = bVec.Eta();
  *(floatBranches_[treeName]["isoEffArea1"]) = getIso(a,"effectiveArea");
  *(floatBranches_[treeName]["isoEffArea2"]) = getIso(b,"effectiveArea");
  *(floatBranches_[treeName]["isoPF1"]) = getIso(a,"deltaBeta");
  *(floatBranches_[treeName]["isoPF2"]) = getIso(b,"deltaBeta");
  *(floatBranches_[treeName]["dB1"]) = getDeltaB(a);
  *(floatBranches_[treeName]["dB2"]) = getDeltaB(b);
  *(floatBranches_[treeName]["mt1"]) = transverseMass(aVec, met);
  *(floatBranches_[treeName]["mt2"]) = transverseMass(bVec, met);
  *(floatBranches_[treeName]["eff1"]) = efficiencies_(a);
  *(floatBranches_[treeName]["eff2"]) = efficiencies_(b);
  *(floatBranches_[treeName]["fakeWeight1"]) = fakeRates_(a);
  *(floatBranches_[treeName]["fakeWeight2"]) = fakeRates_(b);
  *(floatBranches_[treeName]["deltaPhi"]) = aVec.DeltaPhi( bVec );
  *(floatBranches_[treeName]["deltaR"]) = aVec.DeltaR( bVec );
  *(floatBranches_[treeName]["jzb"]) = (uncorrectedMet+comb).Pt() - comb.Pt();




  if(debug) std::cout << "dB1: "<< *(floatBranches_[treeName]["dB1"]) 
		      << "dB2: "<< *(floatBranches_[treeName]["dB2"])<< std::endl;


  *(intBranches_[treeName]["matchesSingleElectron1"]) = 0;
  *(intBranches_[treeName]["matchesSingleElectron2"]) = 0;
  *(intBranches_[treeName]["matchesSingleMuon1"]) = 0;
  *(intBranches_[treeName]["matchesSingleMuon2"]) = 0;
  *(intBranches_[treeName]["matchesDoubleElectron1"]) = 0;
  *(intBranches_[treeName]["matchesDoubleElectron2"]) = 0;

  *(intBranches_[treeName]["matchesDoubleMuonLeading1"]) = 0;
  *(intBranches_[treeName]["matchesDoubleMuonLeading2"]) = 0;
  *(intBranches_[treeName]["matchesDoubleMuonTrailing1"]) = 0;
  *(intBranches_[treeName]["matchesDoubleMuonTrailing2"]) = 0;
  *(intBranches_[treeName]["matchesDoubleMuonLeadingTk1"]) = 0;
  *(intBranches_[treeName]["matchesDoubleMuonLeadingTk2"]) = 0;
  *(intBranches_[treeName]["matchesDoubleMuonTrailingTk1"]) = 0;
  *(intBranches_[treeName]["matchesDoubleMuonTrailingTk2"]) = 0;
  *(intBranches_[treeName]["matchesEMuLeading1"]) = 0;
  *(intBranches_[treeName]["matchesEMuLeading2"]) = 0;
  *(intBranches_[treeName]["matchesEMuTrailing1"]) = 0;
  *(intBranches_[treeName]["matchesEMuTrailing2"]) = 0;
  *(intBranches_[treeName]["matchesMuELeading1"]) = 0;
  *(intBranches_[treeName]["matchesMuELeading2"]) = 0;
  *(intBranches_[treeName]["matchesMuETrailing1"]) = 0;
  *(intBranches_[treeName]["matchesMuETrailing2"]) = 0;

 const trigger::TriggerObjectCollection & triggerObjects = triggerEvent -> getObjects();

 size_t nFilters       = triggerEvent -> sizeFilters();
 for (size_t iFilter = 0; iFilter < nFilters; ++iFilter) {
      TString          name = triggerEvent -> filterTag ( iFilter ).label();

//       if (name.Contains("hltEle27WP80TrackIsoFilter")){
	const trigger::Keys& keys = triggerEvent -> filterKeys( iFilter );
	const trigger::Vids& vids = triggerEvent -> filterIds ( iFilter );
	int nKeys = (int) keys.size();
	int nVids = (int) vids.size();
	assert(nKeys == nVids);
	vector<TLorentzVector> triggerObjectP4s;
	vector<int>            triggerObjectIds;

	for (int iTriggerObject = 0; iTriggerObject < nKeys; ++iTriggerObject ) {
	
			// Get the object ID and key
			int                id  = vids[iTriggerObject];
			trigger::size_type key = keys[iTriggerObject];
	
			// Get the trigger object from the key
			const trigger::TriggerObject & triggerObject = triggerObjects[key];
	
			// Store the trigger object as a TLorentzVector (borrowed from S. Harper)
			TLorentzVector p4;
			p4.SetPxPyPzE(triggerObject.px  (),
				triggerObject.py (),
				triggerObject.pz (),
				triggerObject.energy() );
	
			triggerObjectP4s.push_back ( p4 ) ;
			triggerObjectIds.push_back ( id ) ;
	
	} // end loop over keys/trigger objects passing filters
		
	if ( nKeys > 0 ) {
			for (int iFilterObject = 0; iFilterObject < nKeys; ++iFilterObject) {
				if (!triggerObjectP4s[iFilterObject].Pt() == 0) {
					if (sqrt(pow(aVec.Eta()-triggerObjectP4s[iFilterObject].Eta(),2)+pow(aVec.Phi()-triggerObjectP4s[iFilterObject].Phi(),2)) < 0.2){
						if (name.Contains("hltEle27WP80TrackIsoFilter")){
							*(intBranches_[treeName]["matchesSingleElectron1"]) = 1;
						}
						if (name.Contains("hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15")){
							*(intBranches_[treeName]["matchesSingleMuon1"]) = 1;
						}
						if (name.Contains("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ")){
							*(intBranches_[treeName]["matchesDoubleElectron1"]) = 1;
						}
						if  (name.Contains("hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17") || name.Contains("hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17")){
							*(intBranches_[treeName]["matchesDoubleMuonLeading1"]) = 1;
						}
						if  (name.Contains("hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8") || name.Contains("hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8")){
							*(intBranches_[treeName]["matchesDoubleMuonTrailing1"]) = 1;
						}
						if          (name.Contains("hltL3fL1sMu10MuOpenL1f0L2f10L3Filtered17")|| name.Contains("hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17")){
							*(intBranches_[treeName]["matchesDoubleMuonLeadingTk1"]) = 1;
						}
						if (name.Contains("hltDiMuonGlbFiltered17TrkFiltered8")){
							*(intBranches_[treeName]["matchesDoubleMuonTrailingTk1"]) = 1;
						}
						if (name.Contains("hltL1Mu12EG7L3MuFiltered17")){
							*(intBranches_[treeName]["matchesMuELeading1"]) = 1;
						}
						if (name.Contains("hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter")){
							*(intBranches_[treeName]["matchesMuETrailing1"]) = 1;
						}
						if (name.Contains("hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter")){
							*(intBranches_[treeName]["matchesEMuLeading1"]) = 1;
						}
						if (name.Contains("hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8")){
							*(intBranches_[treeName]["matchesEMuTrailing1"]) = 1;
						}
					}
					if (sqrt(pow(bVec.Eta()-triggerObjectP4s[iFilterObject].Eta(),2)+pow(bVec.Phi()-triggerObjectP4s[iFilterObject].Phi(),2)) < 0.2){
						if (name.Contains("hltEle27WP80TrackIsoFilter")){
							*(intBranches_[treeName]["matchesSingleElectron2"]) = 1;
						}
						if (name.Contains("hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15")){
							*(intBranches_[treeName]["matchesSingleMuon2"]) = 1;
						}
						if (name.Contains("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ")){
							*(intBranches_[treeName]["matchesDoubleElectron2"]) = 1;
						}
						if (name.Contains("hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17") || name.Contains("hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17")){
							*(intBranches_[treeName]["matchesDoubleMuonLeading2"]) = 1;
						}
						if (name.Contains("hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8")|| name.Contains("hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8")){
							*(intBranches_[treeName]["matchesDoubleMuonTrailing2"]) = 1;
						}
						if (name.Contains("hltL3fL1sMu10MuOpenL1f0L2f10L3Filtered17")|| name.Contains("hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17")){
							*(intBranches_[treeName]["matchesDoubleMuonLeadingTk2"]) = 1;
						}
						if (name.Contains("hltDiMuonGlbFiltered17TrkFiltered8")){
							*(intBranches_[treeName]["matchesDoubleMuonTrailingTk2"]) = 1;
						}
						if (name.Contains("hltL1Mu12EG7L3MuFiltered17")){
							*(intBranches_[treeName]["matchesMuELeading2"]) = 1;
						}
						if (name.Contains("hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter")){
							*(intBranches_[treeName]["matchesMuETrailing2"]) = 1;
						}
						if (name.Contains("hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter")){
							*(intBranches_[treeName]["matchesEMuLeading2"]) = 1;
						}
						if (name.Contains("hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8")){
							*(intBranches_[treeName]["matchesEMuTrailing2"]) = 1;
						}
					}
				}
			}
	}


//       }
  







    }




  int matched = 0;
  int pdgId1 = 0;
  int pdgId2 = 0;
  int aMother = -99999;
  int bMother = -99999;

  pdgId1 = getPdgId_.operator()<aT>(a);
  pdgId2 = getPdgId_.operator()<bT>(b);

  TLorentzVector genVec( 0., 0., 0., 0. );
  if(a.genLepton() != NULL){
    matched |= 1;
 // getLeptonPdgId(*(a.genLepton()));
    aMother = getMotherPdgId(*(a.genLepton()));
  }
  if(b.genLepton() != NULL){
    matched |= 2;
//getLeptonPdgId(*(b.genLepton()));
    bMother = getMotherPdgId(*(b.genLepton()));
  }
  TLorentzVector genLepton1;
  TLorentzVector genLepton2;
  if(a.genLepton() != NULL && b.genLepton() != NULL){
      genLepton1.SetPxPyPzE(a.genLepton()->px(),a.genLepton()->py(),a.genLepton()->pz(),a.genLepton()->energy());
      genLepton2.SetPxPyPzE(b.genLepton()->px(),b.genLepton()->py(),b.genLepton()->pz(),b.genLepton()->energy());

      genVec.SetPxPyPzE(a.genLepton()->px()+b.genLepton()->px(),a.genLepton()->py()+b.genLepton()->py(),a.genLepton()->pz()+b.genLepton()->pz(),a.genLepton()->energy()+b.genLepton()->energy());
  }
  if( matched == 3 && aMother == bMother ) {
    matched |= 4;
  }
  *(intBranches_[treeName]["pdgId1"]) = pdgId1;
  *(intBranches_[treeName]["pdgId2"]) = pdgId2;
  *(intBranches_[treeName]["matched"]) = matched;
  *(intBranches_[treeName]["motherPdgId1"]) = aMother;
  *(intBranches_[treeName]["motherPdgId2"]) = bMother;
  *(tLorentzVectorBranches_[treeName]["p4Gen"]) = genVec;
  if(debug) std::cout << ", matched = "<<matched<<", motherId = "<<aMother;
  if(debug) std::cout<<", M = "<< comb.M() <<", chargeProduct = "<< a.charge()*b.charge() <<std::endl;
  
  trees_[treeName]->Fill();
}

int 
DiLeptonTriggerTrees::getMotherPdgId( const reco::GenParticle &p)
{
  int result = -9999;
  if(p.mother() != NULL){
    if(p.status() == 3)
      result = p.mother()->pdgId();
    else if(p.mother()->mother() != NULL)
      result = p.mother()->mother()->pdgId();
  }
  return result;
}

int 
DiLeptonTriggerTrees::getLeptonPdgId( const reco::GenParticle &p)
{
  int result = -9999;
  if(p.status() == 3)
     result = p.pdgId();
   else if(p.mother() != NULL){
      if(abs(p.mother()->pdgId()) == 11 || abs(p.mother()->pdgId()) == 13 || abs(p.mother()->pdgId()) == 15)
        result = p.mother()->pdgId();
      else
        result = p.pdgId();
   }
  return result;
}

//from DQM/Physics/src/EwkTauDQM.cc
std::pair<double, double>
DiLeptonTriggerTrees::calcPZeta(const TLorentzVector& p1,const TLorentzVector& p2, const TLorentzVector& met)
{
	double cosPhi1 = cos(p1.Phi());
	double sinPhi1 = sin(p1.Phi());
	double cosPhi2 = cos(p2.Phi());
	double sinPhi2 = sin(p2.Phi());
	double zetaX = cosPhi1 + cosPhi2;
	double zetaY = sinPhi1 + sinPhi2;
	double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
	if ( zetaR > 0. ) {
		zetaX /= zetaR;
		zetaY /= zetaR;
	}

	double pxVis = p1.Px() + p2.Px();
	double pyVis = p1.Py() + p2.Py();
	double pZetaVis = pxVis*zetaX + pyVis*zetaY;

	double px = pxVis + met.Px();
	double py = pyVis + met.Py();
	double pZeta = px*zetaX + py*zetaY;

	return std::pair<double, double>(pZeta, pZetaVis);
}

void DiLeptonTriggerTrees::fillPdfUncert(const edm::Handle< std::vector<double> >& weightHandle, const std::string& pdfIdentifier, const std::string& treeName){
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

const TLorentzVector DiLeptonTriggerTrees::getMomentum(const  pat::Electron &e)
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

const TLorentzVector DiLeptonTriggerTrees::getMomentum(const  pat::Muon &mu)
{
  const TLorentzVector result = TLorentzVector(mu.px(), mu.py(), mu.pz(), mu.energy());
  return result;
}

const TLorentzVector DiLeptonTriggerTrees::getMomentum(const  pat::Tau &tau)
{
  const TLorentzVector result = TLorentzVector(tau.px(), tau.py(), tau.pz(), tau.energy());
  return result;
}

float DiLeptonTriggerTrees::getIso(const  pat::Electron &e, const std::string &method)
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

float DiLeptonTriggerTrees::getIso(const  pat::Muon &mu, const std::string &method)
{
  //  std::cout<<"muon " << (lepton.trackIso()+lepton.ecalIso()+lepton.hcalIso())/lepton.pt() << std::endl;
  //  return (mu.isolationR03().hadEt + mu.isolationR03().emEt + mu.isolationR03().sumPt) / mu.pt();
  //  return (mu.chargedHadronIso() + mu.photonIso() + mu.neutralHadronIso()) / mu.pt();
  return fctIsolation_(mu,method)* 1./mu.pt();
}

float DiLeptonTriggerTrees::getIso(const  pat::Tau &tau, const std::string &method)
{
  float result = fctIsolation_(tau,method);
  if(tau.tauID(tauId_) < 0.5)
    result *= -1.0;
  return result;
}


float DiLeptonTriggerTrees::getPhiAtECAL(const  pat::Tau &tau, const std::vector<reco::PFCandidate>&pfCands)
{

  float result = -99.0;
  return result;
}


float DiLeptonTriggerTrees::getPhiAtECAL(const  pat::Muon &mu, const std::vector<reco::PFCandidate>&pfCands)
{
  float minDist = 0.5; 
  float deltaR = 0.;
  float result = -99.0;	
	for(std::vector<reco::PFCandidate>::const_iterator it = pfCands.begin(); it != pfCands.end() ; ++it){
		
		deltaR = sqrt(pow((it->eta()-mu.eta()),2)+pow((it->phi()-mu.phi()),2));
		if (deltaR < minDist && abs(it->pdgId()) == 13 && abs(mu.pt()-it->pt())/mu.pt() < 0.1){
			minDist = deltaR;
			result = it->positionAtECALEntrance().Phi();
			
		}
	}

  return result;
}

float DiLeptonTriggerTrees::getPhiAtECAL(const  pat::Electron &e, const std::vector<reco::PFCandidate>&pfCands)
{
  float minDist = 0.5; 
  float deltaR = 0.;
  float result = -99.0;	
	for(std::vector<reco::PFCandidate>::const_iterator it = pfCands.begin(); it != pfCands.end() ; ++it){
		
		deltaR = sqrt(pow((it->eta()-e.eta()),2)+pow((it->phi()-e.phi()),2));
		if (deltaR < minDist && abs(it->pdgId()) == 11 && abs(e.pt()-it->pt())/e.pt() < 0.1){
			minDist = deltaR;
			result = it->positionAtECALEntrance().Phi();
			
		}
	}

  return result;

}

float DiLeptonTriggerTrees::getDeltaB(const  pat::Electron &e)
{
  float result = e.dB(pat::Electron::PV3D);
  return result;
}

float DiLeptonTriggerTrees::getDeltaB(const  pat::Muon &mu)
{
  float result = mu.dB(pat::Muon::PV3D);
  return result;
}

float DiLeptonTriggerTrees::getDeltaB(const  pat::Tau &tau)
{
  float result = -1; // not available for pat::Tau could use ip of leading ch. Hadr if needed.
  return result;
}


float DiLeptonTriggerTrees::transverseMass(const TLorentzVector& p, const TLorentzVector& met)
{
  reco::Candidate::LorentzVector otherMet(met.Px(),met.Py(),met.Pz(),met.E());
  reco::Candidate::LorentzVector leptonT(p.Px(),p.Py(),0.,p.E()*sin(p.Theta()));
  reco::Candidate::LorentzVector sumT=leptonT+otherMet;

  return std::sqrt(sumT.M2());
}

std::string DiLeptonTriggerTrees::convertInputTag(const edm::InputTag tag)
{
  std::string result = tag.label();
  if(tag.instance().length() > 0)
    result = tag.instance();
  //  std::cerr << "'"<<tag.label() << "', '"<< tag.instance()<<"' = '"<< result<<"'"<<std::endl;
  return result;
}

// ------------ Method called once each job just before starting event loop  ------------
void 
DiLeptonTriggerTrees::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiLeptonTriggerTrees::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiLeptonTriggerTrees);
