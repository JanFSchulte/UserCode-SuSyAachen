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
// $Id: DiLeptonTrees.cc,v 1.10 2011/01/17 10:57:12 nmohr Exp $
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
  template <class aT, class bT> void makeCombinations( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT >&b, const edm::Event &ev, double &ht, const TLorentzVector &met, double &weight);
  template <class aT> void makeCombinations( const std::string &treeName, const std::vector<aT> &a, const edm::Event &ev, double &ht, const TLorentzVector &met, double &weight);
  template<class aT, class bT> void fillTree( const std::string &treeName, const aT &a, const bT &b, double &ht, const TLorentzVector &met, double &weight);
  int getMotherPdgId( const reco::GenParticle &p);
  std::pair<double, double> calcPZeta(const TLorentzVector& p1,const TLorentzVector& p2, const TLorentzVector& met);
  void fillPdfUncert(const edm::Handle< std::vector<double> >& weightHandle, const std::string& pdfIdentifier, const std::string& treeName);

  edm::InputTag eTag_;
  edm::InputTag muTag_;
  edm::InputTag tauTag_;
  edm::InputTag jetTag_;
  edm::InputTag metTag_;
  std::vector<edm::ParameterSet> susyVars_;
  std::vector<edm::InputTag> pdfs_;

  //data
  std::map<std::string, TTree*> trees_;  
  std::map<std::string, std::map< std::string, float*> > floatBranches_; 
  std::map<std::string, std::map< std::string, int*> > intBranches_; 
  std::map<std::string, std::map< std::string, TLorentzVector*> > tLorentzVectorBranches_;

  WeightFunctor fakeRates_;

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
  jetTag_ = iConfig.getParameter<edm::InputTag>("jets");
  metTag_ = iConfig.getParameter<edm::InputTag>("met");
  susyVars_ = iConfig.getParameter< std::vector<edm::ParameterSet> >("susyVars");
  pdfs_ = iConfig.getParameter<std::vector<edm::InputTag> > ("pdfWeightTags");

  fakeRates_.SetSource(iConfig,"fakeRates");// TODO use these and add mcInfo flag to choose right rates...

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
  initFloatBranch( "pt1" );
  initFloatBranch( "pt2" );
  initFloatBranch( "eta1" );
  initFloatBranch( "eta2" );
  initFloatBranch( "deltaPhi" );
  initFloatBranch( "deltaR" );
  initFloatBranch( "jzb" );
  initFloatBranch( "ht" );
  initFloatBranch( "met" );
  initFloatBranch( "pZeta" );
  initFloatBranch( "pZetaVis" );
  initIntBranch( "runNr" );
  initIntBranch( "lumiSec" );
  initIntBranch( "eventNr" );
  initIntBranch( "matched" );
  initIntBranch( "motherPdgId" );
  for ( std::vector<edm::ParameterSet>::iterator susyVar_i = susyVars_.begin(); susyVar_i != susyVars_.end(); ++susyVar_i ) {
        std::string var = susyVar_i->getParameter<std::string>( "var" );
        std::string type = susyVar_i->getParameter<std::string>( "type" );
        if(debug) std::cout << var << " of type " << type << std::endl;
        if (type=="int") initIntBranch( var );
        else if (type=="float") initFloatBranch( var );
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
  
  edm::Handle< std::vector< pat::Jet > > jets;
  iEvent.getByLabel(jetTag_, jets);
  
  edm::Handle< std::vector< pat::MET > > mets;
  iEvent.getByLabel(metTag_, mets);

  TLorentzVector met(mets->front().px(), mets->front().py(), mets->front().pz(), mets->front().energy());

  double ht = 0.0;
  for(std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end() ; ++it){
        ht += (*it).pt();
  }
  double weight = 1.;// TODO get weight...
  
  makeCombinations< pat::Electron >("EE", *electrons, iEvent, ht, met, weight);
  makeCombinations< pat::Electron, pat::Muon >("EMu", *electrons, *muons, iEvent, ht,met, weight);
  makeCombinations< pat::Muon >("MuMu", *muons, iEvent, ht, met , weight);
  makeCombinations< pat::Electron, pat::Tau >("ETau", *electrons, *taus, iEvent, ht, met , weight);
  makeCombinations< pat::Muon, pat::Tau>("MuTau", *muons, *taus, iEvent, ht, met, weight);
  makeCombinations< pat::Tau >("TauTau", *taus, iEvent, ht, met, weight);

  //  if( nMu != 2) std::cout << "-------! "<<nMu<<std::endl;
}

template <class aT, class bT> void 
DiLeptonTrees::makeCombinations ( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT> &b, const edm::Event &ev, double &ht, const TLorentzVector &met, double &weight)
{
  *(intBranches_[treeName]["runNr"]) = ev.id().run();
  *(intBranches_[treeName]["lumiSec"]) = ev.id().luminosityBlock();
  *(intBranches_[treeName]["eventNr"]) = ev.id().event();
  for ( std::vector<edm::ParameterSet>::iterator susyVar_i = susyVars_.begin(); susyVar_i != susyVars_.end(); ++susyVar_i ) {
        std::string var = susyVar_i->getParameter<std::string>( "var" );
        std::string type = susyVar_i->getParameter<std::string>( "type" );
        edm::Handle< double > var_;
        ev.getByLabel(var, var_);
        if (type=="float") *(floatBranches_[treeName][var]) = float(*var_);
        else if (type=="int") *(intBranches_[treeName][var]) = int(*var_);
  }
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
      fillTree<aT,bT>( treeName, *itA, *itB, ht, met, weight); 
    }
  }
}

template <class aT> void 
DiLeptonTrees::makeCombinations ( const std::string &treeName, const std::vector<aT> &a, const edm::Event &ev, double &ht, const TLorentzVector &met, double &weight)
{
  *(intBranches_[treeName]["runNr"]) = ev.id().run();
  *(intBranches_[treeName]["lumiSec"]) = ev.id().luminosityBlock();
  *(intBranches_[treeName]["eventNr"]) = ev.id().event();
  for ( std::vector<edm::ParameterSet>::iterator susyVar_i = susyVars_.begin(); susyVar_i != susyVars_.end(); ++susyVar_i ) {
        std::string var = susyVar_i->getParameter<std::string>( "var" );
        std::string type = susyVar_i->getParameter<std::string>( "type" );
        edm::Handle< double > var_;
        ev.getByLabel(var, var_);
        if (type=="int") *(intBranches_[treeName][var]) = int(*var_);
        else if (type=="float") *(floatBranches_[treeName][var]) = float(*var_);
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
      fillTree<aT, aT>( treeName, *itA, *itB, ht, met, weight); 
    }
  }
}

template <class aT, class bT> void 
DiLeptonTrees::fillTree( const std::string &treeName, const aT& a, const bT& b, double &ht, const TLorentzVector &met, double &weight)
{
  if(debug) std::cout << treeName << "- pts:"<< a.pt() << " " << b.pt();
  TLorentzVector aVec( a.px(), a.py(), a.pz(), a.energy() );
  TLorentzVector bVec( b.px(), b.py(), b.pz(), b.energy() );
  TLorentzVector comb = aVec+bVec;
  std::pair<double, double> pZeta = calcPZeta(a.p(), b.p(), met);
  *(floatBranches_[treeName]["weight"]) = weight;
  *(floatBranches_[treeName]["chargeProduct"]) = a.charge()*b.charge();
  *(tLorentzVectorBranches_[treeName]["p4"]) = comb;
  *(floatBranches_[treeName]["pt1"]) = aVec.Pt();
  *(floatBranches_[treeName]["pt2"]) = bVec.Pt();
  *(floatBranches_[treeName]["eta1"]) = aVec.Eta();
  *(floatBranches_[treeName]["eta2"]) = bVec.Eta();
  *(floatBranches_[treeName]["deltaPhi"]) = aVec.DeltaPhi( bVec );
  *(floatBranches_[treeName]["deltaR"]) = aVec.DeltaR( bVec );
  *(floatBranches_[treeName]["jzb"]) = (met+comb).Pt() - comb.Pt();
  *(floatBranches_[treeName]["ht"]) = ht;
  *(floatBranches_[treeName]["met"]) = (met).Pt();
  *(floatBranches_[treeName]["pZeta"]) = pZeta.first;
  *(floatBranches_[treeName]["pZetaVis"]) = pZeta.second;

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

//from DQM/Physics/src/EwkTauDQM.cc
std::pair<double, double>
DiLeptonTrees::calcPZeta(const TLorentzVector& p1,const TLorentzVector& p2, const TLorentzVector& met)
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

void DiLeptonTrees::fillPdfUncert(const edm::Handle< std::vector<double> >& weightHandle, const std::string& pdfIdentifier, const std::string& treeName){
     std::string up = "Up";
     std::string down = "Down";
     double centralValue = (*weightHandle)[0];
     *(floatBranches_[treeName][pdfIdentifier]) = float(centralValue);
     if(debug) std::cout << "Cen" << treeName << ": " << centralValue << std::endl;
     unsigned int nmembers = weightHandle->size();
     double wminus = 0.;
     double wplus = 0.;
     for (unsigned int j=1; j<nmembers; j+=2) {
        float wa = ((*weightHandle)[j]-centralValue)/centralValue;
        float wb = ((*weightHandle)[j+1]-centralValue)/centralValue;
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
    if (wplus>0.) wplus = sqrt(wplus);
    if (wminus>0.) wminus = sqrt(wminus);
    *(floatBranches_[treeName][pdfIdentifier+down]) = float(wminus);
    *(floatBranches_[treeName][pdfIdentifier+up]) = float(wplus);
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
