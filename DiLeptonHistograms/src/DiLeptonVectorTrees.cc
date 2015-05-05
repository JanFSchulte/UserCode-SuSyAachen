// -*- C++ -*-
//
// Package:    Histograms
// Class:      DiLeptonVectorTrees
// 
/**\class DiLeptonVectorTrees DiLeptonVectorTrees.cc brot/DiLeptonVectorTrees/src/DiLeptonVectorTrees.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  matthias edelhoff
//         Created:  Tue Oct 27 13:50:40 CET 2009
// $Id: DiLeptonVectorTrees.cc,v 1.31 2012/09/17 17:38:58 sprenger Exp $
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

#include <SuSyAachen/DiLeptonHistograms/interface/WeightFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/PdgIdFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/VertexWeightFunctor.h>
#include <SuSyAachen/TagAndProbeTreeWriter/interface/IsolationFunctor.h>

//ROOT
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
//
// class decleration
//

class DiLeptonVectorTrees : public edm::EDAnalyzer {
public:
  explicit DiLeptonVectorTrees(const edm::ParameterSet&);
  ~DiLeptonVectorTrees();

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
  template <class aT, class bT> void makeCombinations( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT >&b, const edm::Event &ev, const pat::MET &patMet, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties);
  template <class aT> void makeCombinations( const std::string &treeName, const std::vector<aT> &a, const edm::Event &ev, const pat::MET &patMet, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties,  const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties);
  template<class aT, class bT> void fillTree( const std::string &treeName, const aT &a, const bT &b, const pat::MET &patMet);
  int getLeptonPdgId( const reco::GenParticle &p);
  int getMotherPdgId( const reco::GenParticle &p);
  std::pair<double, double> calcPZeta(const TLorentzVector& p1,const TLorentzVector& p2, const TLorentzVector& met);
  void fillPdfUncert(const edm::Handle< std::vector<double> >& weightHandle, const std::string& pdfIdentifier, const std::string& treeName);

  const TLorentzVector getMomentum(const  pat::Electron &e);
  const TLorentzVector getMomentum(const  pat::Muon &mu);
  const TLorentzVector getMomentum(const  pat::Tau &tau);
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
  edm::InputTag vertexTag_;
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
DiLeptonVectorTrees::DiLeptonVectorTrees(const edm::ParameterSet& iConfig):
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
  vertexTag_ = iConfig.getParameter<edm::InputTag>("vertices");
  susyVars_ = iConfig.getParameter< std::vector<edm::ParameterSet> >("susyVars");
  pdfs_ = iConfig.getParameter<std::vector<edm::InputTag> > ("pdfWeightTags");

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
  trees_["EE"] = file->make<TTree>("EEDileptonVectorTree", "EE DileptonVectorTree");
  trees_["EMu"] = file->make<TTree>("EMuDileptonVectorTree", "EMu DileptonVectorTree");
  trees_["MuMu"] = file->make<TTree>("MuMuDileptonVectorTree", "MuMu DileptonVectorTree");
  if(useTaus_){
    trees_["ETau"] = file->make<TTree>("ETauDileptonVectorTree", "ETau DileptonTree");
    trees_["MuTau"] = file->make<TTree>("MuTauDileptonVectorTree", "MuTau DileptonVectorTree");
    trees_["TauTau"] = file->make<TTree>("TauTauDileptonVectorTree", "TauTau DileptonVectorTree");
  }
  initFloatBranch( "weight" );
  initFloatBranch( "chargeProduct" );
  initTLorentzVectorBranch( "p4" );
  initTLorentzVectorBranch( "vMet" );
  initTLorentzVectorBranch( "p4Gen" );
  initTLorentzVectorBranch( "lepton1" );
  initTLorentzVectorBranch( "lepton2" );
  initTLorentzVectorBranch( "jet1" );
  initTLorentzVectorBranch( "jet2" );
  initTLorentzVectorBranch( "jet3" );
  initTLorentzVectorBranch( "jet4" );
  initTLorentzVectorBranch( "bJet1" );
  initTLorentzVectorBranch( "bJet2" );
  initTLorentzVectorBranch( "bJet3" );
  initTLorentzVectorBranch( "bJet4" );
  initTLorentzVectorBranch( "vMetType1" );
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
  initFloatBranch( "jzb" );
  initFloatBranch( "ht" );
  initFloatBranch( "met" );
  initFloatBranch( "type1met" );
  initFloatBranch( "tcmet" );
  initFloatBranch( "pZeta" );
  initFloatBranch( "pZetaVis" );
  initIntBranch( "nJets" );
  initIntBranch( "nBJets" );
  initIntBranch( "nVertices" );
  initIntBranch( "nLightLeptons" );
  initFloatBranch( "sqrts" );
  initIntBranch( "runNr" );
  initIntBranch( "lumiSec" );
  initIntBranch( "eventNr" );
  initIntBranch( "pdgId1" );
  initIntBranch( "pdgId2" );
  initIntBranch( "matched" );
  initIntBranch( "motherPdgId" );
  initIntBranch( "IntConvCand" );
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
DiLeptonVectorTrees::initTLorentzVectorBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    tLorentzVectorBranches_[(*it).first][name] = new TLorentzVector;
    (*it).second->Branch(name.c_str(), "TLorentzVector" ,&tLorentzVectorBranches_[(*it).first][name]);
  }
}

void 
DiLeptonVectorTrees::initFloatBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    floatBranches_[(*it).first][name] = new float;
    (*it).second->Branch(name.c_str(), floatBranches_[(*it).first][name], (name+"/F").c_str());
  }
}

void 
DiLeptonVectorTrees::initIntBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    intBranches_[(*it).first][name] = new int;
    (*it).second->Branch(name.c_str(), intBranches_[(*it).first][name], (name+"/I").c_str());
  }
}

DiLeptonVectorTrees::~DiLeptonVectorTrees()
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
DiLeptonVectorTrees::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< std::vector< pat::Electron > > electrons;
  iEvent.getByLabel(eTag_, electrons);

  edm::Handle< std::vector< pat::Muon > > muons;
  iEvent.getByLabel(muTag_, muons);

  
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
  intEventProperties["nJets"] = jets_->size();
  intEventProperties["nBJets"] = bJets->size();
  intEventProperties["nLightLeptons"] = electrons->size() + muons->size();
  intEventProperties["runNr"] = iEvent.id().run();
  intEventProperties["lumiSec"] = iEvent.id().luminosityBlock();
  intEventProperties["eventNr"] = iEvent.id().event();

  intEventProperties["IntConvCand"] = 0;
  if (electrons->size()>2 || muons->size()>2 || (muons->size()+electrons->size()) >2) intEventProperties["IntConvCand"] = 1;

  pat::MET met = mets->front();
  TLorentzVector metVector(mets->front().px(), mets->front().py(), mets->front().pz(), mets->front().energy());
  //TLorentzVector uncorrectedMet;
  //  uncorrectedMet.SetPtEtaPhiE(mets->front().uncorrectedPt(), 0,	
  //  mets->front().uncorrectedPhi(), mets->front().uncorrectedPt());
  floatEventProperties["met"] = metVector.Pt();


  reco::MET type1met = type1mets->front();
  TLorentzVector type1metVector(type1mets->front().px(), type1mets->front().py(), type1mets->front().pz(), type1mets->front().energy());
  tLorentzVectorEventProperties["vMetType1"] = type1metVector;
  floatEventProperties["type1met"] = type1metVector.Pt();

  pat::MET tcmet = tcmets->front();
  TLorentzVector tcmetVector(tcmets->front().px(), tcmets->front().py(), tcmets->front().pz(), tcmets->front().energy());

  floatEventProperties["tcmet"] = tcmetVector.Pt();


  floatEventProperties["ht"] = 0.0;
  for(std::vector<pat::Jet>::const_iterator it = jets_->begin(); it != jets_->end() ; ++it){
        floatEventProperties["ht"] += (*it).pt();
	
  }

  // Jets
  TLorentzVector jet1Vector(0.,0.,0.,0.);
  TLorentzVector jet2Vector(0.,0.,0.,0.);
  TLorentzVector jet3Vector(0.,0.,0.,0.);
  TLorentzVector jet4Vector(0.,0.,0.,0.); 


  if (jets_->size() > 0)
    jet1Vector.SetPxPyPzE(jets_->at(0).px(),jets_->at(0).py(),jets_->at(0).pz(),jets_->at(0).energy());
  if (jets_->size() > 1)
    jet2Vector.SetPxPyPzE(jets_->at(1).px(),jets_->at(1).py(),jets_->at(1).pz(),jets_->at(1).energy());
  if (jets_->size() > 2)
    jet3Vector.SetPxPyPzE(jets_->at(2).px(),jets_->at(2).py(),jets_->at(2).pz(),jets_->at(2).energy());
  if (jets_->size() > 3)
    jet4Vector.SetPxPyPzE(jets_->at(3).px(),jets_->at(3).py(),jets_->at(3).pz(),jets_->at(3).energy());
  tLorentzVectorEventProperties["jet1"] = jet1Vector;
  tLorentzVectorEventProperties["jet2"] = jet2Vector;
  tLorentzVectorEventProperties["jet3"] = jet3Vector;
  tLorentzVectorEventProperties["jet4"] = jet4Vector;

  // bjets

  TLorentzVector bJet1Vector(0.,0.,0.,0.);
  TLorentzVector bJet2Vector(0.,0.,0.,0.);
  TLorentzVector bJet3Vector(0.,0.,0.,0.);
  TLorentzVector bJet4Vector(0.,0.,0.,0.); 


  if (bJets->size() > 0)
    bJet1Vector.SetPxPyPzE(bJets->at(0).px(),bJets->at(0).py(),bJets->at(0).pz(),bJets->at(0).energy());
  if (bJets->size() > 1)
    bJet2Vector.SetPxPyPzE(bJets->at(1).px(),bJets->at(1).py(),bJets->at(1).pz(),bJets->at(1).energy());
  if (bJets->size() > 2)
    bJet3Vector.SetPxPyPzE(bJets->at(2).px(),bJets->at(2).py(),bJets->at(2).pz(),bJets->at(2).energy());
  if (bJets->size() > 3)
    bJet4Vector.SetPxPyPzE(bJets->at(3).px(),bJets->at(3).py(),bJets->at(3).pz(),bJets->at(3).energy());
  tLorentzVectorEventProperties["bJet1"] = bJet1Vector;
  tLorentzVectorEventProperties["bJet2"] = bJet2Vector;
  tLorentzVectorEventProperties["bJet3"] = bJet3Vector;
  tLorentzVectorEventProperties["bJet4"] = bJet4Vector;


    	

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

  makeCombinations< pat::Electron >("EE", *electrons, iEvent, met, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  makeCombinations< pat::Electron, pat::Muon >("EMu", *electrons, *muons, iEvent, met, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  makeCombinations< pat::Muon >("MuMu", *muons, iEvent, met, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  if(useTaus_){
    makeCombinations< pat::Electron, pat::Tau >("ETau", *electrons, *taus, iEvent, met, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
    makeCombinations< pat::Muon, pat::Tau>("MuTau", *muons, *taus, iEvent, met, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
    makeCombinations< pat::Tau >("TauTau", *taus, iEvent, met, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  }
  //  if( nMu != 2) std::cout << "-------! "<<nMu<<std::endl;
}

template <class aT, class bT> void 
DiLeptonVectorTrees::makeCombinations ( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT> &b, const edm::Event &ev, const pat::MET &patMet, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties)
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
      fillTree<aT,bT>( treeName, *itA, *itB, patMet); 
    }
  }
}

template <class aT> void 
DiLeptonVectorTrees::makeCombinations ( const std::string &treeName, const std::vector<aT> &a, const edm::Event &ev, const pat::MET &patMet, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties)
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
      fillTree<aT, aT>( treeName, *itA, *itB, patMet); 
    }
  }
}

template <class aT, class bT> void 
DiLeptonVectorTrees::fillTree( const std::string &treeName, const aT& a, const bT& b, const pat::MET &patMet)
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
  std::pair<double, double> pZeta = calcPZeta( aVec, bVec, met);
  *(floatBranches_[treeName]["chargeProduct"]) = a.charge()*b.charge();
  *(tLorentzVectorBranches_[treeName]["p4"]) = comb;
  *(tLorentzVectorBranches_[treeName]["lepton1"]) = aVec;
  *(tLorentzVectorBranches_[treeName]["lepton2"]) = bVec;
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
  *(floatBranches_[treeName]["jzb"]) = (uncorrectedMet+comb).Pt() - comb.Pt();
  *(floatBranches_[treeName]["pZeta"]) = pZeta.first;
  *(floatBranches_[treeName]["pZetaVis"]) = pZeta.second;

  if (jets_->size() >= 2){
    TLorentzVector vJet1 = TLorentzVector(jets_->at(0).p4().x(), jets_->at(0).p4().y(), jets_->at(0).p4().z(), jets_->at(0).p4().t());
    TLorentzVector vJet2 = TLorentzVector(jets_->at(0).p4().x(), jets_->at(0).p4().y(), jets_->at(0).p4().z(), jets_->at(0).p4().t());
    TLorentzVector sub = aVec + bVec + vJet1 + vJet2;
    *(floatBranches_[treeName]["sqrts"]) = std::sqrt(TMath::Power((std::sqrt(sub.M2() + sub.Perp2()) + met.Et()), 2.0) - (sub.Vect() + met.Vect()).Mag2());
  } else {
    *(floatBranches_[treeName]["sqrts"]) = -1.0;
  }

  if(debug) std::cout << "dB1: "<< *(floatBranches_[treeName]["dB1"]) 
		      << "dB2: "<< *(floatBranches_[treeName]["dB2"])<< std::endl;

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
  if(a.genLepton() != NULL && b.genLepton() != NULL){
      genVec.SetPxPyPzE(a.genLepton()->px()+b.genLepton()->px(),a.genLepton()->py()+b.genLepton()->py(),a.genLepton()->pz()+b.genLepton()->pz(),a.genLepton()->energy()+b.genLepton()->energy());
  }
  if( matched == 3 && aMother == bMother ) {
    matched |= 4;
  }
  *(intBranches_[treeName]["pdgId1"]) = pdgId1;
  *(intBranches_[treeName]["pdgId2"]) = pdgId2;
  *(intBranches_[treeName]["matched"]) = matched;
  *(intBranches_[treeName]["motherPdgId"]) = aMother;
  *(tLorentzVectorBranches_[treeName]["p4Gen"]) = genVec;
  if(debug) std::cout << ", matched = "<<matched<<", motherId = "<<aMother;
  if(debug) std::cout<<", M = "<< comb.M() <<", chargeProduct = "<< a.charge()*b.charge() <<std::endl;
  
  trees_[treeName]->Fill();
}

int 
DiLeptonVectorTrees::getMotherPdgId( const reco::GenParticle &p)
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
DiLeptonVectorTrees::getLeptonPdgId( const reco::GenParticle &p)
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
DiLeptonVectorTrees::calcPZeta(const TLorentzVector& p1,const TLorentzVector& p2, const TLorentzVector& met)
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

void DiLeptonVectorTrees::fillPdfUncert(const edm::Handle< std::vector<double> >& weightHandle, const std::string& pdfIdentifier, const std::string& treeName){
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

const TLorentzVector DiLeptonVectorTrees::getMomentum(const  pat::Electron &e)
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

const TLorentzVector DiLeptonVectorTrees::getMomentum(const  pat::Muon &mu)
{
  const TLorentzVector result = TLorentzVector(mu.px(), mu.py(), mu.pz(), mu.energy());
  return result;
}

const TLorentzVector DiLeptonVectorTrees::getMomentum(const  pat::Tau &tau)
{
  const TLorentzVector result = TLorentzVector(tau.px(), tau.py(), tau.pz(), tau.energy());
  return result;
}

float DiLeptonVectorTrees::getIso(const  pat::Electron &e, const std::string &method)
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

float DiLeptonVectorTrees::getIso(const  pat::Muon &mu, const std::string &method)
{
  //  std::cout<<"muon " << (lepton.trackIso()+lepton.ecalIso()+lepton.hcalIso())/lepton.pt() << std::endl;
  //  return (mu.isolationR03().hadEt + mu.isolationR03().emEt + mu.isolationR03().sumPt) / mu.pt();
  //  return (mu.chargedHadronIso() + mu.photonIso() + mu.neutralHadronIso()) / mu.pt();
  return fctIsolation_(mu,method)* 1./mu.pt();
}

float DiLeptonVectorTrees::getIso(const  pat::Tau &tau, const std::string &method)
{
  float result = fctIsolation_(tau,method);
  if(tau.tauID(tauId_) < 0.5)
    result *= -1.0;
  return result;
}

float DiLeptonVectorTrees::getDeltaB(const  pat::Electron &e)
{
  float result = e.dB(pat::Electron::PV3D);
  return result;
}

float DiLeptonVectorTrees::getDeltaB(const  pat::Muon &mu)
{
  float result = mu.dB(pat::Muon::PV3D);
  return result;
}

float DiLeptonVectorTrees::getDeltaB(const  pat::Tau &tau)
{
  float result = -1; // not available for pat::Tau could use ip of leading ch. Hadr if needed.
  return result;
}


float DiLeptonVectorTrees::transverseMass(const TLorentzVector& p, const TLorentzVector& met)
{
  reco::Candidate::LorentzVector otherMet(met.Px(),met.Py(),met.Pz(),met.E());
  reco::Candidate::LorentzVector leptonT(p.Px(),p.Py(),0.,p.E()*sin(p.Theta()));
  reco::Candidate::LorentzVector sumT=leptonT+otherMet;

  return std::sqrt(sumT.M2());
}

std::string DiLeptonVectorTrees::convertInputTag(const edm::InputTag tag)
{
  std::string result = tag.label();
  if(tag.instance().length() > 0)
    result = tag.instance();
  //  std::cerr << "'"<<tag.label() << "', '"<< tag.instance()<<"' = '"<< result<<"'"<<std::endl;
  return result;
}

// ------------ Method called once each job just before starting event loop  ------------
void 
DiLeptonVectorTrees::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiLeptonVectorTrees::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiLeptonVectorTrees);
