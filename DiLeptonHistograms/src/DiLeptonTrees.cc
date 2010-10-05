// -*- C++ -*-
//
// Package:    Histograms
// Class:      DiLeptonTrees
// 
/**\class DiLeptonTrees DiLeptonTrees.cc brot/DiLeptonTrees/src/DiLeptonTrees.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  matthias edelhoff
//         Created:  Tue Oct 27 13:50:40 CET 2009
// $Id: DiLeptonTrees.cc,v 1.3 2010/09/30 18:58:36 edelhoff Exp $
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

#include <DataFormats/Provenance/interface/EventID.h>

//ROOT
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

//
// class decleration
//

class DiLeptonTrees : public edm::EDAnalyzer {
public:
  explicit DiLeptonTrees(const edm::ParameterSet&);
  ~DiLeptonTrees();

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
  template <class aT, class bT> void makeCombinations( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT >&b, const edm::EventID &id, double &weight);
  template <class aT> void makeCombinations( const std::string &treeName, const std::vector<aT> &a, const edm::EventID &id, double &weight);
  template<class aT, class bT> void fillTree( const std::string &treeName, const aT &a, const bT &b, double &weight);
  int getMotherPdgId( const reco::GenParticle &p);

  edm::InputTag eTag_;
  edm::InputTag muTag_;
  edm::InputTag tauTag_;

  //data
  std::map<std::string, TTree*> trees_;  
  std::map<std::string, std::map< std::string, float*> > floatBranches_; 
  std::map<std::string, std::map< std::string, int*> > intBranches_; 
  std::map<std::string, std::map< std::string, TLorentzVector*> > tLorentzVectorBranches_; 

  bool debug;
};

// constructors and destructor
DiLeptonTrees::DiLeptonTrees(const edm::ParameterSet& iConfig)
{
  debug = false;

  // read config
  eTag_ = iConfig.getParameter<edm::InputTag>("electrons");
  muTag_ = iConfig.getParameter<edm::InputTag>("muons");
  tauTag_ = iConfig.getParameter<edm::InputTag>("taus");

  // init trees
  edm::Service<TFileService> file;
  trees_["EE"] = file->make<TTree>("EEDileptonTree", "EE DileponTree");
  trees_["EMu"] = file->make<TTree>("EMuDileptonTree", "EMu DileponTree");
  trees_["MuMu"] = file->make<TTree>("MuMuDileptonTree", "MuMu DileponTree");
  trees_["ETau"] = file->make<TTree>("ETauDileptonTree", "ETau DileponTree");
  trees_["MuTau"] = file->make<TTree>("MuTauDileptonTree", "MuTau DileponTree");
  trees_["TauTau"] = file->make<TTree>("TauTauDileptonTree", "TauTau DileponTree");
  initFloatBranch( "weight" );
  initFloatBranch( "chargeProduct" );
  initTLorentzVectorBranch( "p4" );
  initFloatBranch( "deltaPhi" );
  initFloatBranch( "deltaR" );
  initIntBranch( "runNr" );
  initIntBranch( "eventNr" );
  initIntBranch( "matched" );
  initIntBranch( "motherPdgId" );
}

void 
DiLeptonTrees::initTLorentzVectorBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    tLorentzVectorBranches_[(*it).first][name] = new TLorentzVector;
    (*it).second->Branch(name.c_str(), "TLorentzVector" ,&tLorentzVectorBranches_[(*it).first][name]);
  }
}

void 
DiLeptonTrees::initFloatBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    floatBranches_[(*it).first][name] = new float;
    (*it).second->Branch(name.c_str(), floatBranches_[(*it).first][name], (name+"/F").c_str());
  }
}

void 
DiLeptonTrees::initIntBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    intBranches_[(*it).first][name] = new int;
    (*it).second->Branch(name.c_str(), intBranches_[(*it).first][name], (name+"/I").c_str());
  }
}

DiLeptonTrees::~DiLeptonTrees()
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
DiLeptonTrees::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< std::vector< pat::Electron > > electrons;
  iEvent.getByLabel(eTag_, electrons);

  edm::Handle< std::vector< pat::Muon > > muons;
  iEvent.getByLabel(muTag_, muons);
  
  edm::Handle< std::vector< pat::Tau > > taus;
  iEvent.getByLabel(tauTag_, taus);

  double weight = 1.;// TODO get weight...

  makeCombinations< pat::Electron >("EE", *electrons, iEvent.id(), weight);
  makeCombinations< pat::Electron, pat::Muon >("EMu", *electrons, *muons, iEvent.id(), weight);
  makeCombinations< pat::Muon >("MuMu", *muons, iEvent.id(), weight);
  makeCombinations< pat::Electron, pat::Tau >("ETau", *electrons, *taus, iEvent.id(), weight);
  makeCombinations< pat::Muon, pat::Tau>("MuTau", *muons, *taus, iEvent.id(), weight);
  makeCombinations< pat::Tau >("TauTau", *taus, iEvent.id(), weight);

  //  if( nMu != 2) std::cout << "-------! "<<nMu<<std::endl;
}

template <class aT, class bT> void 
DiLeptonTrees::makeCombinations ( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT> &b, const edm::EventID &id, double &weight)
{
  *(intBranches_[treeName]["runNr"]) = id.run();
  *(intBranches_[treeName]["eventNr"]) = id.event();
  for( typename std::vector<aT>::const_iterator itA = a.begin(); itA != a.end(); ++itA){
    for( typename std::vector<bT>::const_iterator itB = b.begin(); itB != b.end(); ++itB){
      fillTree<aT,bT>( treeName, *itA, *itB, weight); 
    }
  }
}

template <class aT> void 
DiLeptonTrees::makeCombinations ( const std::string &treeName, const std::vector<aT> &a, const edm::EventID &id, double &weight)
{
  *(intBranches_[treeName]["runNr"]) = id.run();
  *(intBranches_[treeName]["eventNr"]) = id.event();
  for( typename std::vector<aT>::const_iterator itA = a.begin(); itA != a.end(); ++itA){
    for( typename std::vector<aT>::const_iterator itB = a.begin(); itB != itA; ++itB){
      fillTree<aT, aT>( treeName, *itA, *itB, weight); 
    }
  }
}

template <class aT, class bT> void 
DiLeptonTrees::fillTree( const std::string &treeName, const aT& a, const bT& b, double &weight)
{
  if(debug) std::cout << treeName << "- pts:"<< a.pt() << " " << b.pt();
  TLorentzVector aVec( a.px(), a.py(), a.pz(), a.energy() );
  TLorentzVector bVec( b.px(), b.py(), b.pz(), b.energy() );
  TLorentzVector comb = aVec+bVec;
  *(floatBranches_[treeName]["weight"]) = weight;
  *(floatBranches_[treeName]["chargeProduct"]) = a.charge()*b.charge();
  *(tLorentzVectorBranches_[treeName]["p4"]) = comb;
  *(floatBranches_[treeName]["deltaPhi"]) = aVec.DeltaPhi( bVec );
  *(floatBranches_[treeName]["deltaR"]) = aVec.DeltaR( bVec );
  int matched = 0;
  int aMother = -99999;
  int bMother = -99999;

  if(a.genLepton() != NULL){
    matched |= 1;
    aMother = getMotherPdgId(*(a.genLepton()));
  }
  if(b.genLepton() != NULL){
    matched |= 2;
    bMother = getMotherPdgId(*(b.genLepton()));
  }
  if( matched == 3 && aMother == bMother ) {
    matched |= 4;
  }
  *(intBranches_[treeName]["matched"]) = matched;
  *(intBranches_[treeName]["motherPdgId"]) = aMother;
  if(debug) std::cout << ", matched = "<<matched<<", motherId = "<<aMother;
  if(debug) std::cout<<", M = "<< comb.M() <<", chargeProduct = "<< a.charge()*b.charge() <<std::endl;
  
  trees_[treeName]->Fill();
}

int 
DiLeptonTrees::getMotherPdgId( const reco::GenParticle &p)
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

// ------------ method called once each job just before starting event loop  ------------
void 
DiLeptonTrees::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiLeptonTrees::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiLeptonTrees);
