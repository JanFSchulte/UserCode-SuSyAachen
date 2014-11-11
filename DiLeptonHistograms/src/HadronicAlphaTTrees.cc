// -*- C++ -*-
//
// Package:    Histograms
// Class:      HadronicAlphaTTrees
// 
/**\class HadronicAlphaTTrees HadronicAlphaTTrees.cc brot/HadronicAlphaTTrees/src/HadronicAlphaTTrees.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  matthias edelhoff
//         Created:  Tue Oct 27 13:50:40 CET 2009
// $Id: HadronicAlphaTTrees.cc,v 1.31 2012/09/17 17:38:58 sprenger Exp $
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

class HadronicAlphaTTrees : public edm::EDAnalyzer {
public:
  explicit HadronicAlphaTTrees(const edm::ParameterSet&);
  ~HadronicAlphaTTrees();

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
  template <class aT, class bT> void makeCombinations( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT >&b, const std::vector<reco::PFCandidate>&pfCands, const edm::Event &ev, const pat::MET &patMet, const TLorentzVector &MHT, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties);
  template <class aT> void makeCombinations( const std::string &treeName, const std::vector<aT> &a , const std::vector<reco::PFCandidate>&pfCands,  const edm::Event &ev, const pat::MET &patMet, const TLorentzVector &MHT, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties);
  template<class aT, class bT> void fillTree( const std::string &treeName, const aT &a, const bT &b, const std::vector<reco::PFCandidate>&pfCands, const pat::MET &patMet, const TLorentzVector &MHT);
  //~ template<class aT> void fillLepton( const std::string &treeName, const aT &a, const std::string &leptonFlavorNumber);
  //~ template<class aT> void fillLepton( const std::string &treeName, const aT &a, const std::string &leptonFlavorNumber, const  std::map<std::string, float> &floatEventProperties);
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
  float topPtWeightBen(double topPt);
  float topPtWeightTOP(double topPt);
  float getId(const  pat::Electron &e);
  float getId(const  pat::Muon &mu);
  float getId(const  pat::Tau &tau);
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
  edm::InputTag vertexTagUp_;
  edm::InputTag vertexTagDown_;
  edm::InputTag vertexTagBlockA_;
  edm::InputTag vertexTagBlockAUp_;
  edm::InputTag vertexTagBlockADown_;
  edm::InputTag vertexTagBlockB_;
  edm::InputTag vertexTagBlockBUp_;
  edm::InputTag vertexTagBlockBDown_;
  edm::InputTag pfCandTag_;
  edm::InputTag genParticleTag_;
  edm::InputTag rhoTag_;
  std::vector<edm::ParameterSet> susyVars_;
  std::vector<edm::InputTag> pdfs_;
  std::string tauId_;


  std::map<double, double> electronCorrections_;
  //data
  std::map<std::string, TTree*> trees_;  
  std::map<std::string, std::map< std::string, float*> > floatBranches_; 
  std::map<std::string, std::map< std::string, unsigned int*> > intBranches_; 
  std::map<std::string, std::map< std::string, TLorentzVector*> > tLorentzVectorBranches_;

  edm::Handle< std::vector< pat::Jet > > jets_;

  WeightFunctor fakeRates_;
  WeightFunctor efficiencies_;
  VertexWeightFunctor fctVtxWeight_;
  VertexWeightFunctor fctVtxWeightBlockA_;
  VertexWeightFunctor fctVtxWeightBlockB_;
  VertexWeightFunctor fctVtxWeightUp_;
  VertexWeightFunctor fctVtxWeightBlockAUp_;
  VertexWeightFunctor fctVtxWeightBlockBUp_;
  VertexWeightFunctor fctVtxWeightDown_;
  VertexWeightFunctor fctVtxWeightBlockADown_;
  VertexWeightFunctor fctVtxWeightBlockBDown_;
  IsolationFunctor fctIsolation_;
  PdgIdFunctor getPdgId_;

  bool debug;
  bool useJets2_;
  bool useTaus_;
};

// constructors and destructor
HadronicAlphaTTrees::HadronicAlphaTTrees(const edm::ParameterSet& iConfig):
  fctVtxWeight_    (iConfig.getParameter<edm::ParameterSet>("vertexWeights") ),
  fctVtxWeightBlockA_    (iConfig.getParameter<edm::ParameterSet>("vertexWeightsBlockA") ),
  fctVtxWeightBlockB_    (iConfig.getParameter<edm::ParameterSet>("vertexWeightsBlockB") ),
  fctVtxWeightUp_    (iConfig.getParameter<edm::ParameterSet>("vertexWeightsUp") ),
  fctVtxWeightBlockAUp_    (iConfig.getParameter<edm::ParameterSet>("vertexWeightsBlockAUp") ),
  fctVtxWeightBlockBUp_    (iConfig.getParameter<edm::ParameterSet>("vertexWeightsBlockBUp") ),
  fctVtxWeightDown_    (iConfig.getParameter<edm::ParameterSet>("vertexWeightsDown") ),
  fctVtxWeightBlockADown_    (iConfig.getParameter<edm::ParameterSet>("vertexWeightsBlockADown") ),
  fctVtxWeightBlockBDown_    (iConfig.getParameter<edm::ParameterSet>("vertexWeightsBlockBDown") ),
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
  rhoTag_ = iConfig.getParameter<edm::InputTag>("rho");
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
  trees_["HadronicAlphaT"] = file->make<TTree>("HadronicAlphaTTree", "HadronicAlphaTTree");
  //~ if(useTaus_){
    //~ trees_["ETau"] = file->make<TTree>("ETauDileptonTree", "ETau DileponTree");
    //~ trees_["MuTau"] = file->make<TTree>("MuTauDileptonTree", "MuTau DileponTree");
    //~ trees_["TauTau"] = file->make<TTree>("TauTauDileptonTree", "TauTau DileponTree");
  //~ }
  initFloatBranch( "weight" );
  initFloatBranch( "weightBlockA" );
  initFloatBranch( "weightBlockB" );
  initFloatBranch( "weightUp" );
  initFloatBranch( "weightBlockAUp" );
  initFloatBranch( "weightBlockBUp" );
  initFloatBranch( "weightDown" );
  initFloatBranch( "weightBlockADown" );
  initFloatBranch( "weightBlockBDown" );
  //~ initFloatBranch( "topWeightBen" );
  //~ initFloatBranch( "topWeightTOP" );
  //~ initFloatBranch( "genPtTop1" );
  //~ initFloatBranch( "genPtTop2" );
  //~ initFloatBranch( "chargeProduct" );
  initTLorentzVectorBranch( "vMet" );
  initTLorentzVectorBranch( "vGenMet" );
  initTLorentzVectorBranch( "vGenMetNeutrinos" );
  initTLorentzVectorBranch( "vHt" );
  initTLorentzVectorBranch( "vHt37" );
  initTLorentzVectorBranch( "vHt43" );
  initTLorentzVectorBranch( "jet1" );
  initTLorentzVectorBranch( "jet2" );
  initTLorentzVectorBranch( "jet3" );
  initTLorentzVectorBranch( "jet4" );
  initTLorentzVectorBranch( "jet5" );
  initTLorentzVectorBranch( "jet6" );
  initTLorentzVectorBranch( "jet7" );
  initTLorentzVectorBranch( "jet8" );
  initTLorentzVectorBranch( "jet9" );
  initTLorentzVectorBranch( "jet10" );
  initTLorentzVectorBranch( "jet11" );
  initTLorentzVectorBranch( "jet1_37" );
  initTLorentzVectorBranch( "jet2_37" );
  initTLorentzVectorBranch( "jet3_37" );
  initTLorentzVectorBranch( "jet4_37" );
  initTLorentzVectorBranch( "jet5_37" );
  initTLorentzVectorBranch( "jet6_37" );
  initTLorentzVectorBranch( "jet1_43" );
  initTLorentzVectorBranch( "jet2_43" );
  initTLorentzVectorBranch( "jet3_43" );
  initTLorentzVectorBranch( "jet4_43" );
  initTLorentzVectorBranch( "jet5_43" );
  initTLorentzVectorBranch( "jet6_43" );
  initTLorentzVectorBranch( "bJet1" );
  initTLorentzVectorBranch( "bJet2" );
  initTLorentzVectorBranch( "bJet3" );
  initTLorentzVectorBranch( "bJet4" );
  initTLorentzVectorBranch( "bJet5" );
  initTLorentzVectorBranch( "bJet6" );
  initTLorentzVectorBranch( "bJet1_37" );
  initTLorentzVectorBranch( "bJet2_37" );
  initTLorentzVectorBranch( "bJet3_37" );
  initTLorentzVectorBranch( "bJet4_37" );
  initTLorentzVectorBranch( "bJet5_37" );
  initTLorentzVectorBranch( "bJet6_37" );
  initTLorentzVectorBranch( "bJet1_43" );
  initTLorentzVectorBranch( "bJet2_43" );
  initTLorentzVectorBranch( "bJet3_43" );
  initTLorentzVectorBranch( "bJet4_43" );
  initTLorentzVectorBranch( "bJet5_43" );
  initTLorentzVectorBranch( "bJet6_43" );
  initTLorentzVectorBranch( "vMetType1" );
  initFloatBranch( "rho" );
  //~ initFloatBranch( "pt1" );
  //~ initFloatBranch( "pt2" );
  //~ initFloatBranch( "phiAtECALEntrance1");
  //~ initFloatBranch( "phiAtECALEntrance2");
  //~ initFloatBranch( "eta1" );
  //~ initFloatBranch( "eta2" );
  //~ initFloatBranch( "mt1" );
  //~ initFloatBranch( "mt2" );
  initFloatBranch( "mbb" );
  //~ initFloatBranch( "eff1" );
  //~ initFloatBranch( "eff2" );
  //~ initFloatBranch( "fakeWeight1" );
  //~ initFloatBranch( "fakeWeight2" );
  //~ initFloatBranch( "deltaPhiJetMET" );
  //~ initFloatBranch( "deltaPhiSecondJetMET" );
  //~ initFloatBranch( "deltaPhi" );
  //~ initFloatBranch( "deltaPhiLeptonMET1");
  //~ initFloatBranch( "deltaPhiLeptonMET2");
  //~ initFloatBranch( "deltaR" );
  //~ initFloatBranch( "jzb" );
  initFloatBranch( "alphaT" );
  initFloatBranch( "alphaT37" );
  initFloatBranch( "alphaT43" );
  initFloatBranch( "absvHt" );
  initFloatBranch( "absvHt37" );
  initFloatBranch( "absvHt43" );
  initFloatBranch( "deltaHt" );
  initFloatBranch( "deltaHt37" );
  initFloatBranch( "deltaHt43" );
  initFloatBranch( "ht" );
  initFloatBranch( "ht37" );
  initFloatBranch( "ht43" );
  initFloatBranch( "htJESUp" );
  initFloatBranch( "htJESDown" );
  //~ initFloatBranch( "mht" );
  initFloatBranch( "met" );
  initFloatBranch( "metJESUp" );
  initFloatBranch( "metJESDown" );
  initFloatBranch( "metChargedEMEtFraction" );
  initFloatBranch( "metChargedHadEtFraction" );
  initFloatBranch( "metMuonEtFraction" );
  initFloatBranch( "metNeutralEMFraction" );
  initFloatBranch( "metNeutralHadEtFraction" );
  initFloatBranch( "metType6EtFraction" );
  initFloatBranch( "metType7EtFraction" );
  initFloatBranch( "type1Met" );
  initFloatBranch( "tcMet" );
  initFloatBranch( "caloMet" );
  initFloatBranch( "genMetTrue" );
  //~ initFloatBranch( "pZeta" );
  //~ initFloatBranch( "pZetaVis" );
  initIntBranch( "nJets" );
  initIntBranch( "nJets37" );
  initIntBranch( "nJets43" );
  initIntBranch( "nBJets" );
  initIntBranch( "nBJets37" );
  initIntBranch( "nBJets43" );
  initIntBranch( "nShiftedJetsJESUp" );
  initIntBranch( "nShiftedJetsJESDown" );
  initIntBranch( "nVertices" );
  initIntBranch( "nLightLeptons" );
  initIntBranch( "nGenSUSYLeptons" );
  initIntBranch( "nGenSUSYLightLeptons" );
  initIntBranch( "nGenSUSYTaus" );
  initIntBranch( "nGenSUSYNeutrinos" );
  initIntBranch( "nElectrons" );
  initIntBranch( "nMuons" );
  //~ initFloatBranch( "jet1pt" );
  //~ initFloatBranch( "jet2pt" );
  //~ initFloatBranch( "jet3pt" );
  //~ initFloatBranch( "jet4pt" );
  //~ initFloatBranch( "bjet1pt" );
  //~ initFloatBranch( "bjet2pt" );
  //~ initFloatBranch( "bjet3pt" );
  //~ initFloatBranch( "bjet4pt" );
  //~ initFloatBranch( "sqrts" );
  initIntBranch( "runNr" );
  initIntBranch( "lumiSec" );
  initIntBranch( "eventNr" );
  //~ initIntBranch( "pdgId1" );
  //~ initIntBranch( "pdgId2" );
  //~ initIntBranch( "matched" );
  //~ initIntBranch( "motherPdgId1" );
  //~ initIntBranch( "motherPdgId2" );
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
HadronicAlphaTTrees::initTLorentzVectorBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    tLorentzVectorBranches_[(*it).first][name] = new TLorentzVector;
    (*it).second->Branch(name.c_str(), "TLorentzVector" ,&tLorentzVectorBranches_[(*it).first][name]);
  }
}

void 
HadronicAlphaTTrees::initFloatBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    floatBranches_[(*it).first][name] = new float;
    (*it).second->Branch(name.c_str(), floatBranches_[(*it).first][name], (name+"/F").c_str());
  }
}

void 
HadronicAlphaTTrees::initIntBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    intBranches_[(*it).first][name] = new unsigned int;
    (*it).second->Branch(name.c_str(), intBranches_[(*it).first][name], (name+"/I").c_str());
  }
}

HadronicAlphaTTrees::~HadronicAlphaTTrees()
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
HadronicAlphaTTrees::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  edm::Handle<double> rho;
  iEvent.getByLabel(rhoTag_,rho);
  floatEventProperties["rho"] = (float)(*rho);


  //~ intEventProperties["nBJets"] = bJets->size();
  intEventProperties["nLightLeptons"] = electrons->size() + muons->size();
  intEventProperties["nElectrons"] = electrons->size();
  intEventProperties["nMuons"] = muons->size();
  intEventProperties["runNr"] = iEvent.id().run();
  intEventProperties["lumiSec"] = iEvent.id().luminosityBlock();
  intEventProperties["eventNr"] = iEvent.id().event();



  pat::MET met = mets->front();
  TLorentzVector metVector(mets->front().px(), mets->front().py(), mets->front().pz(), mets->front().energy());
  //TLorentzVector uncorrectedMet;
  //  uncorrectedMet.SetPtEtaPhiE(mets->front().uncorrectedPt(), 0,	
  //  mets->front().uncorrectedPhi(), mets->front().uncorrectedPt());
  floatEventProperties["met"] = metVector.Pt();
  floatEventProperties["metChargedEMEtFraction"]= met.ChargedEMEtFraction();
  floatEventProperties["metChargedHadEtFraction"]= met.ChargedHadEtFraction();
  floatEventProperties["metMuonEtFraction"]= met.MuonEtFraction();
  floatEventProperties["metNeutralEMFraction"]= met.NeutralEMFraction();
  floatEventProperties["metNeutralHadEtFraction"]= met.NeutralHadEtFraction();
  floatEventProperties["metType6EtFraction"]= met.Type6EtFraction();
  floatEventProperties["metType7EtFraction"]= met.Type7EtFraction();
  
  tLorentzVectorEventProperties["vMet"] = metVector;



  reco::MET type1met = type1mets->front();
  TLorentzVector type1metVector(type1mets->front().px(), type1mets->front().py(), type1mets->front().pz(), type1mets->front().energy());

  floatEventProperties["type1Met"] = type1metVector.Pt();
  tLorentzVectorEventProperties["vMetType1"] = type1metVector;
  pat::MET tcmet = tcmets->front();
  TLorentzVector tcmetVector(tcmets->front().px(), tcmets->front().py(), tcmets->front().pz(), tcmets->front().energy());

  floatEventProperties["tcMet"] = tcmetVector.Pt();

  pat::MET calomet = calomets->front();
  TLorentzVector calometVector(calomets->front().px(), calomets->front().py(), calomets->front().pz(), calomets->front().energy());

  floatEventProperties["caloMet"] = calometVector.Pt();


  if (!genMetTrues.isValid()){
	floatEventProperties["genMetTrue"] = -1.;
	tLorentzVectorEventProperties["vGenMet"] = TLorentzVector(0.,0.,0.,0.);

  }
  else{
	reco::GenMET genMetTrue = genMetTrues->front();
  	TLorentzVector genMetTrueVector(genMetTrues->front().px(), genMetTrues->front().py(), genMetTrues->front().pz(), genMetTrues->front().energy());
	floatEventProperties["genMetTrue"] = genMetTrueVector.Pt();
	tLorentzVectorEventProperties["vGenMet"] = genMetTrueVector;
// 


  }

  TLorentzVector genNeutrino1(0.,0.,0.,0.);
  TLorentzVector genNeutrino2(0.,0.,0.,0.);
  bool foundOne = false;
  int nGenLeptons = 0;
  int nGenLightLeptons = 0;
  int nGenTaus = 0;
  int nGenNeutrinos = 0;
  if (genParticles.isValid()){
	
	for (std::vector<reco::GenParticle>::const_iterator itGenParticle = genParticles->begin(); itGenParticle != genParticles->end(); itGenParticle++) {
		if ((abs((*itGenParticle).pdgId())==12 || abs((*itGenParticle).pdgId())==14) && abs(getMotherPdgId(*itGenParticle))==24 ){
			
			if (!foundOne){
				genNeutrino1.SetPxPyPzE((*itGenParticle).px(), (*itGenParticle).py(), (*itGenParticle).pz(), (*itGenParticle).energy());
				foundOne=true;
			}
			else{
				genNeutrino2.SetPxPyPzE((*itGenParticle).px(), (*itGenParticle).py(), (*itGenParticle).pz(), (*itGenParticle).energy());
			}				

		}
		if ((abs((*itGenParticle).pdgId())==12 || abs((*itGenParticle).pdgId())==14 || abs((*itGenParticle).pdgId())==16) && ((*itGenParticle).status() == 1 || (*itGenParticle).status() == 2) && abs(getMotherPdgId(*itGenParticle))==1000023 ){
			nGenNeutrinos++;


		}
		if ((abs((*itGenParticle).pdgId())==11 || abs((*itGenParticle).pdgId())==13 || abs((*itGenParticle).pdgId())==15) && ((*itGenParticle).status() == 1 || (*itGenParticle).status() == 2) && abs(getMotherPdgId(*itGenParticle))==1000023 ){
			nGenLeptons++;


		}
		if ((abs((*itGenParticle).pdgId())==11 || abs((*itGenParticle).pdgId())==13) && ((*itGenParticle).status() == 1 || (*itGenParticle).status() == 2) && abs(getMotherPdgId(*itGenParticle))==1000023 ){
			nGenLightLeptons++;

		}
		if ( abs((*itGenParticle).pdgId())==15 && ((*itGenParticle).status() == 1 || (*itGenParticle).status() == 2) && abs(getMotherPdgId(*itGenParticle))==1000023 ){
			nGenTaus++;
		}
	}
  }

  intEventProperties["nGenSUSYLeptons"] = nGenLeptons;
  intEventProperties["nGenSUSYLightLeptons"] = nGenLightLeptons;
  intEventProperties["nGenSUSYTaus"] = nGenTaus;
  intEventProperties["nGenSUSYNeutrinos"] = nGenNeutrinos;
  tLorentzVectorEventProperties["vGenMetNeutrinos"] = genNeutrino1+genNeutrino2;

  //~ TLorentzVector MHT;
  //~ TLorentzVector tempMHT;
  //~ TLorentzVector leadingJetMomentum;
  //~ TLorentzVector subLeadingJetMomentum;
  floatEventProperties["alphaT"] = 0.0;
  floatEventProperties["absvHt"] = 0.0;
  floatEventProperties["deltaHt"] = 0.0;
  floatEventProperties["ht"] = 0.0;
  floatEventProperties["htJESUp"] = 0.0;
  floatEventProperties["htJESDown"] = 0.0;
  //~ floatEventProperties["mht"] = 0.0;

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
  int nJets37=0;  
  int nJets43=0;  
  TLorentzVector htVector(0.,0.,0.,0.);
  TLorentzVector htVector37(0.,0.,0.,0.);
  TLorentzVector htVector43(0.,0.,0.,0.);
  
  float px = 0.;
  float py = 0.;
  float pz = 0.;
  float energy = 0.;
  
  for(std::vector<pat::Jet>::const_iterator it = jets_->begin(); it != jets_->end() ; ++it){
	if ((*it).et() >=50.0){
		nJets++;
		floatEventProperties["ht"] += (*it).et();
		px += (*it).px();
		py += (*it).py();
		pz += (*it).pz();
		energy += (*it).energy();
	}
  }
  intEventProperties["nJets"] = nJets;
  htVector.SetPxPyPzE(px,py,pz,energy);
  tLorentzVectorEventProperties["vHt"] = htVector;
  floatEventProperties["absvHt"] = htVector.Pt();


  TLorentzVector jet1Vector(0.,0.,0.,0.);
  TLorentzVector jet2Vector(0.,0.,0.,0.);
  TLorentzVector jet3Vector(0.,0.,0.,0.);
  TLorentzVector jet4Vector(0.,0.,0.,0.); 
  TLorentzVector jet5Vector(0.,0.,0.,0.); 
  TLorentzVector jet6Vector(0.,0.,0.,0.); 
  TLorentzVector jet7Vector(0.,0.,0.,0.); 
  TLorentzVector jet8Vector(0.,0.,0.,0.); 
  TLorentzVector jet9Vector(0.,0.,0.,0.); 
  TLorentzVector jet10Vector(0.,0.,0.,0.); 
  TLorentzVector jet11Vector(0.,0.,0.,0.); 


  if (nJets > 0)
    jet1Vector.SetPxPyPzE(jets_->at(0).px(),jets_->at(0).py(),jets_->at(0).pz(),jets_->at(0).energy());
  if (nJets > 1)
    jet2Vector.SetPxPyPzE(jets_->at(1).px(),jets_->at(1).py(),jets_->at(1).pz(),jets_->at(1).energy());
  if (nJets > 2)
    jet3Vector.SetPxPyPzE(jets_->at(2).px(),jets_->at(2).py(),jets_->at(2).pz(),jets_->at(2).energy());
  if (nJets > 3)
    jet4Vector.SetPxPyPzE(jets_->at(3).px(),jets_->at(3).py(),jets_->at(3).pz(),jets_->at(3).energy());
  if (nJets > 4)
    jet5Vector.SetPxPyPzE(jets_->at(4).px(),jets_->at(4).py(),jets_->at(4).pz(),jets_->at(4).energy());
  if (nJets > 5)
    jet6Vector.SetPxPyPzE(jets_->at(5).px(),jets_->at(5).py(),jets_->at(5).pz(),jets_->at(5).energy());
  if (nJets > 6)
    jet7Vector.SetPxPyPzE(jets_->at(6).px(),jets_->at(6).py(),jets_->at(6).pz(),jets_->at(6).energy());
  if (nJets > 7)
    jet8Vector.SetPxPyPzE(jets_->at(7).px(),jets_->at(7).py(),jets_->at(7).pz(),jets_->at(7).energy());
  if (nJets > 8)
    jet9Vector.SetPxPyPzE(jets_->at(8).px(),jets_->at(8).py(),jets_->at(8).pz(),jets_->at(8).energy());
  if (nJets > 9)
    jet10Vector.SetPxPyPzE(jets_->at(9).px(),jets_->at(9).py(),jets_->at(9).pz(),jets_->at(9).energy());
  if (nJets > 10)
    jet11Vector.SetPxPyPzE(jets_->at(10).px(),jets_->at(10).py(),jets_->at(10).pz(),jets_->at(10).energy());
    
  tLorentzVectorEventProperties["jet1"] = jet1Vector;
  tLorentzVectorEventProperties["jet2"] = jet2Vector;
  tLorentzVectorEventProperties["jet3"] = jet3Vector;
  tLorentzVectorEventProperties["jet4"] = jet4Vector;
  tLorentzVectorEventProperties["jet5"] = jet5Vector;
  tLorentzVectorEventProperties["jet6"] = jet6Vector;
  tLorentzVectorEventProperties["jet7"] = jet7Vector;
  tLorentzVectorEventProperties["jet8"] = jet8Vector;
  tLorentzVectorEventProperties["jet9"] = jet9Vector;
  tLorentzVectorEventProperties["jet10"] = jet10Vector;
  tLorentzVectorEventProperties["jet11"] = jet11Vector;
  
  
  ////~ first bin 275 < HT < 325 GeV, jet ET > 37 cut
  
  px = 0.;
  py = 0.;
  pz = 0.;
  energy = 0.;
  
  for(std::vector<pat::Jet>::const_iterator it = jets_->begin(); it != jets_->end() ; ++it){
	if ((*it).et() >=37.0){
		nJets37++;
		floatEventProperties["ht37"] += (*it).et();
		px += (*it).px();
		py += (*it).py();
		pz += (*it).pz();
		energy += (*it).energy();
	}
  }
  intEventProperties["nJets37"] = nJets37;
  htVector37.SetPxPyPzE(px,py,pz,energy);
  tLorentzVectorEventProperties["vHt37"] = htVector37;
  floatEventProperties["absvHt37"] = htVector37.Pt();


  TLorentzVector jet1Vector37(0.,0.,0.,0.);
  TLorentzVector jet2Vector37(0.,0.,0.,0.);
  TLorentzVector jet3Vector37(0.,0.,0.,0.);
  TLorentzVector jet4Vector37(0.,0.,0.,0.); 
  TLorentzVector jet5Vector37(0.,0.,0.,0.); 
  TLorentzVector jet6Vector37(0.,0.,0.,0.); 


  if (nJets37 > 0)
    jet1Vector37.SetPxPyPzE(jets_->at(0).px(),jets_->at(0).py(),jets_->at(0).pz(),jets_->at(0).energy());
  if (nJets37 > 1)
    jet2Vector37.SetPxPyPzE(jets_->at(1).px(),jets_->at(1).py(),jets_->at(1).pz(),jets_->at(1).energy());
  if (nJets37 > 2)
    jet3Vector37.SetPxPyPzE(jets_->at(2).px(),jets_->at(2).py(),jets_->at(2).pz(),jets_->at(2).energy());
  if (nJets37 > 3)
    jet4Vector37.SetPxPyPzE(jets_->at(3).px(),jets_->at(3).py(),jets_->at(3).pz(),jets_->at(3).energy());
  if (nJets37 > 4)
    jet5Vector37.SetPxPyPzE(jets_->at(4).px(),jets_->at(4).py(),jets_->at(4).pz(),jets_->at(4).energy());
  if (nJets37 > 5)
    jet6Vector37.SetPxPyPzE(jets_->at(5).px(),jets_->at(5).py(),jets_->at(5).pz(),jets_->at(5).energy());
    
  tLorentzVectorEventProperties["jet1_37"] = jet1Vector37;
  tLorentzVectorEventProperties["jet2_37"] = jet2Vector37;
  tLorentzVectorEventProperties["jet3_37"] = jet3Vector37;
  tLorentzVectorEventProperties["jet4_37"] = jet4Vector37;
  tLorentzVectorEventProperties["jet5_37"] = jet5Vector37;
  tLorentzVectorEventProperties["jet6_37"] = jet6Vector37;
  
  ////~ second bin 325 < HT < 375 GeV, jet ET > 43 cut
  
  px = 0.;
  py = 0.;
  pz = 0.;
  energy = 0.;
  
  for(std::vector<pat::Jet>::const_iterator it = jets_->begin(); it != jets_->end() ; ++it){
	if ((*it).et() >=43.0){
		nJets43++;
		floatEventProperties["ht43"] += (*it).et();
		px += (*it).px();
		py += (*it).py();
		pz += (*it).pz();
		energy += (*it).energy();
	}
  }
  intEventProperties["nJets43"] = nJets43;
  htVector43.SetPxPyPzE(px,py,pz,energy);
  tLorentzVectorEventProperties["vHt43"] = htVector43;
  floatEventProperties["absvHt43"] = htVector43.Pt();


  TLorentzVector jet1Vector43(0.,0.,0.,0.);
  TLorentzVector jet2Vector43(0.,0.,0.,0.);
  TLorentzVector jet3Vector43(0.,0.,0.,0.);
  TLorentzVector jet4Vector43(0.,0.,0.,0.); 
  TLorentzVector jet5Vector43(0.,0.,0.,0.); 
  TLorentzVector jet6Vector43(0.,0.,0.,0.); 


  if (nJets43 > 0)
    jet1Vector43.SetPxPyPzE(jets_->at(0).px(),jets_->at(0).py(),jets_->at(0).pz(),jets_->at(0).energy());
  if (nJets43 > 1)
    jet2Vector43.SetPxPyPzE(jets_->at(1).px(),jets_->at(1).py(),jets_->at(1).pz(),jets_->at(1).energy());
  if (nJets43 > 2)
    jet3Vector43.SetPxPyPzE(jets_->at(2).px(),jets_->at(2).py(),jets_->at(2).pz(),jets_->at(2).energy());
  if (nJets43 > 3)
    jet4Vector43.SetPxPyPzE(jets_->at(3).px(),jets_->at(3).py(),jets_->at(3).pz(),jets_->at(3).energy());
  if (nJets43 > 4)
    jet5Vector43.SetPxPyPzE(jets_->at(4).px(),jets_->at(4).py(),jets_->at(4).pz(),jets_->at(4).energy());
  if (nJets43 > 5)
    jet6Vector43.SetPxPyPzE(jets_->at(5).px(),jets_->at(5).py(),jets_->at(5).pz(),jets_->at(5).energy());
    
  tLorentzVectorEventProperties["jet1_43"] = jet1Vector43;
  tLorentzVectorEventProperties["jet2_43"] = jet2Vector43;
  tLorentzVectorEventProperties["jet3_43"] = jet3Vector43;
  tLorentzVectorEventProperties["jet4_43"] = jet4Vector43;
  tLorentzVectorEventProperties["jet5_43"] = jet5Vector43;
  tLorentzVectorEventProperties["jet6_43"] = jet6Vector43;
  
  
  float deltaHT = 0.;
  if (nJets > 1){
	  deltaHT = 1000.;
  }
  if (nJets == 2){
	  deltaHT = abs(jet1Vector.Et() - jet2Vector.Et());
  }
  if (nJets == 3){
	  for (int n1 = 0; n1 < nJets; n1++){
		  for (int n2 = 0; n2 < nJets; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						  if (deltaHT > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et())){
							  deltaHT = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et());
						  }
					  }
				  }
			  }
		  }
	  }
  }
  if (nJets == 4){
	  for (int n1 = 0; n1 < nJets; n1++){
		  for (int n2 = 0; n2 < nJets; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								if (deltaHT > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et())){
									deltaHT = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et());
								}
								if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et())){
									deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et());
								}
								
							}
						}
					}
				}
			}
		}
	 }
  }
  if (nJets == 5){
	  for (int n1 = 0; n1 < nJets; n1++){
		  for (int n2 = 0; n2 < nJets; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										if (deltaHT > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et())){
											deltaHT = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et());
										}
										if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et())){
											deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et());
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

  if (nJets == 6){
	  for (int n1 = 0; n1 < nJets; n1++){
		  for (int n2 = 0; n2 < nJets; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												if (deltaHT > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et())){
													deltaHT = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et());
												}
												if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et())){
													deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et());
												}
												if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et())){
													deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et());
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
  if (nJets == 7){
	  for (int n1 = 0; n1 < nJets; n1++){
		  for (int n2 = 0; n2 < nJets; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												for (int n7 = 0; n7 < nJets; n7++){
													if ((n7 != n1) && (n7 != n2) && (n7 != n3) && (n7 != n4) && (n7 != n5) && (n7 != n6)){
														if (deltaHT > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et())){
															deltaHT = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et());
														}
														if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et())){
															deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et());
														}
														if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et())){
															deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et());
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
  if (nJets == 8){
	  for (int n1 = 0; n1 < nJets; n1++){
		  for (int n2 = 0; n2 < nJets; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												for (int n7 = 0; n7 < nJets; n7++){
													if ((n7 != n1) && (n7 != n2) && (n7 != n3) && (n7 != n4) && (n7 != n5) && (n7 != n6)){
														for (int n8 = 0; n8 < nJets; n8++){
															if ((n8 != n1) && (n8 != n2) && (n8 != n3) && (n8 != n4) && (n8 != n5) && (n8 != n6) && (n8 != n7)){
																if (deltaHT > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et())){
																	deltaHT = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et());
																}
																if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et())){
																	deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et());
																}
																if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et())){
																	deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et());
																}
																if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et())){
																	deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et());
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
  if (nJets == 9){
	  for (int n1 = 0; n1 < nJets; n1++){
		  for (int n2 = 0; n2 < nJets; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												for (int n7 = 0; n7 < nJets; n7++){
													if ((n7 != n1) && (n7 != n2) && (n7 != n3) && (n7 != n4) && (n7 != n5) && (n7 != n6)){
														for (int n8 = 0; n8 < nJets; n8++){
															if ((n8 != n1) && (n8 != n2) && (n8 != n3) && (n8 != n4) && (n8 != n5) && (n8 != n6) && (n8 != n7)){
																for (int n9 = 0; n9 < nJets; n9++){
																	if ((n9 != n1) && (n9 != n2) && (n9 != n3) && (n9 != n4) && (n9 != n5) && (n9 != n6) && (n9 != n7) && (n9 != n8)){
																		if (deltaHT > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et())){
																			deltaHT = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et());
																		}
																		if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et())){
																			deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et());
																		}
																		if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et())){
																			deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et());
																		}
																		if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et())){
																			deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et());
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
  if (nJets == 10){
	  for (int n1 = 0; n1 < nJets; n1++){
		  for (int n2 = 0; n2 < nJets; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												for (int n7 = 0; n7 < nJets; n7++){
													if ((n7 != n1) && (n7 != n2) && (n7 != n3) && (n7 != n4) && (n7 != n5) && (n7 != n6)){
														for (int n8 = 0; n8 < nJets; n8++){
															if ((n8 != n1) && (n8 != n2) && (n8 != n3) && (n8 != n4) && (n8 != n5) && (n8 != n6) && (n8 != n7)){
																for (int n9 = 0; n9 < nJets; n9++){
																	if ((n9 != n1) && (n9 != n2) && (n9 != n3) && (n9 != n4) && (n9 != n5) && (n9 != n6) && (n9 != n7) && (n9 != n8)){
																		for (int n10 = 0; n10 < nJets; n10++){
																			if ((n10 != n1) && (n10 != n2) && (n10 != n3) && (n10 != n4) && (n10 != n5) && (n10 != n6) && (n10 != n7) && (n10 != n8) && (n10 != n9)){
																				if (deltaHT > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et())){
																					deltaHT = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et());
																				}
																				if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et())){
																					deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et());
																				}
																				if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et())){
																					deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et());
																				}
																				if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et())){
																					deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et());
																				}
																				if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() + jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et())){
																					deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() + jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et());
																				}
																			}
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
  if (nJets == 11){
	  for (int n1 = 0; n1 < nJets; n1++){
		  for (int n2 = 0; n2 < nJets; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												for (int n7 = 0; n7 < nJets; n7++){
													if ((n7 != n1) && (n7 != n2) && (n7 != n3) && (n7 != n4) && (n7 != n5) && (n7 != n6)){
														for (int n8 = 0; n8 < nJets; n8++){
															if ((n8 != n1) && (n8 != n2) && (n8 != n3) && (n8 != n4) && (n8 != n5) && (n8 != n6) && (n8 != n7)){
																for (int n9 = 0; n9 < nJets; n9++){
																	if ((n9 != n1) && (n9 != n2) && (n9 != n3) && (n9 != n4) && (n9 != n5) && (n9 != n6) && (n9 != n7) && (n9 != n8)){
																		for (int n10 = 0; n10 < nJets; n10++){
																			if ((n10 != n1) && (n10 != n2) && (n10 != n3) && (n10 != n4) && (n10 != n5) && (n10 != n6) && (n10 != n7) && (n10 != n8) && (n10 != n9)){
																				for (int n11 = 0; n11 < nJets; n11++){
																					if ((n11 != n1) && (n11 != n2) && (n11 != n3) && (n11 != n4) && (n11 != n5) && (n11 != n6) && (n11 != n7) && (n11 != n8) && (n11 != n9) && (n11 != n10)){
																						if (deltaHT > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et())){
																							deltaHT = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et());
																						}
																						if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et())){
																							deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et());
																						}
																						if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et())){
																							deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et());
																						}
																						if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et())){
																							deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et());
																						}
																						if (deltaHT > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() + jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et())){
																							deltaHT = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() + jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et());
																						}
																					}
																				}
																			}
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
										
  if (nJets > 1){
	  floatEventProperties["deltaHt"] = deltaHT;
	  floatEventProperties["alphaT"] = 1./2. * (1 - floatEventProperties["deltaHt"] / floatEventProperties["ht"]) / sqrt( 1. - (floatEventProperties["absvHt"] / floatEventProperties["ht"]) * (floatEventProperties["absvHt"] / floatEventProperties["ht"]));
  }
  else{
	  floatEventProperties["deltaHt"] = 0.;
	  floatEventProperties["alphaT"] = 0.;
  }
  
  float deltaHT37 = 0.;
  if (nJets37 > 1){
	  deltaHT37 = 1000.;
  }
  if (nJets37 == 2){
	  deltaHT37 = abs(jet1Vector37.Et() - jet2Vector37.Et());
  }
  if (nJets37 == 3){
	  for (int n1 = 0; n1 < nJets37; n1++){
		  for (int n2 = 0; n2 < nJets37; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets37; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						  if (deltaHT37 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et())){
							  deltaHT37 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et());
						  }
					  }
				  }
			  }
		  }
	  }
  }
  if (nJets37 == 4){
	  for (int n1 = 0; n1 < nJets37; n1++){
		  for (int n2 = 0; n2 < nJets37; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets37; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets37; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								if (deltaHT37 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et())){
									deltaHT37 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et());
								}
								if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et())){
									deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et());
								}
								
							}
						}
					}
				}
			}
		}
	 }
  }
  if (nJets37 == 5){
	  for (int n1 = 0; n1 < nJets37; n1++){
		  for (int n2 = 0; n2 < nJets37; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets37; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets37; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets37; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										if (deltaHT37 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et())){
											deltaHT37 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et());
										}
										if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et())){
											deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et());
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

  if (nJets37 == 6){
	  for (int n1 = 0; n1 < nJets37; n1++){
		  for (int n2 = 0; n2 < nJets37; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets37; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets37; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets37; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets37; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												if (deltaHT37 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et())){
													deltaHT37 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et());
												}
												if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et())){
													deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et());
												}
												if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et())){
													deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et());
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

 if (nJets37 == 7){
	  for (int n1 = 0; n1 < nJets37; n1++){
		  for (int n2 = 0; n2 < nJets37; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets37; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets37; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets37; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets37; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												for (int n7 = 0; n7 < nJets37; n7++){
													if ((n7 != n1) && (n7 != n2) && (n7 != n3) && (n7 != n4) && (n7 != n5) && (n7 != n6)){
														if (deltaHT37 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et())){
															deltaHT37 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et());
														}
														if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et())){
															deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et());
														}
														if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et())){
															deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et());
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
  if (nJets37 == 8){
	  for (int n1 = 0; n1 < nJets37; n1++){
		  for (int n2 = 0; n2 < nJets37; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets37; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets37; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets37; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets37; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												for (int n7 = 0; n7 < nJets37; n7++){
													if ((n7 != n1) && (n7 != n2) && (n7 != n3) && (n7 != n4) && (n7 != n5) && (n7 != n6)){
														for (int n8 = 0; n8 < nJets37; n8++){
															if ((n8 != n1) && (n8 != n2) && (n8 != n3) && (n8 != n4) && (n8 != n5) && (n8 != n6) && (n8 != n7)){
																if (deltaHT37 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et())){
																	deltaHT37 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et());
																}
																if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et())){
																	deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et());
																}
																if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et())){
																	deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et());
																}
																if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et())){
																	deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et());
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
  if (nJets37 == 9){
	  for (int n1 = 0; n1 < nJets37; n1++){
		  for (int n2 = 0; n2 < nJets37; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets37; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets37; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets37; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets37; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												for (int n7 = 0; n7 < nJets37; n7++){
													if ((n7 != n1) && (n7 != n2) && (n7 != n3) && (n7 != n4) && (n7 != n5) && (n7 != n6)){
														for (int n8 = 0; n8 < nJets37; n8++){
															if ((n8 != n1) && (n8 != n2) && (n8 != n3) && (n8 != n4) && (n8 != n5) && (n8 != n6) && (n8 != n7)){
																for (int n9 = 0; n9 < nJets37; n9++){
																	if ((n9 != n1) && (n9 != n2) && (n9 != n3) && (n9 != n4) && (n9 != n5) && (n9 != n6) && (n9 != n7) && (n9 != n8)){
																		if (deltaHT37 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et())){
																			deltaHT37 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et());
																		}
																		if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et())){
																			deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et());
																		}
																		if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et())){
																			deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et());
																		}
																		if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et())){
																			deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et());
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
  if (nJets37 == 10){
	  for (int n1 = 0; n1 < nJets37; n1++){
		  for (int n2 = 0; n2 < nJets37; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets37; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets37; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets37; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets37; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												for (int n7 = 0; n7 < nJets37; n7++){
													if ((n7 != n1) && (n7 != n2) && (n7 != n3) && (n7 != n4) && (n7 != n5) && (n7 != n6)){
														for (int n8 = 0; n8 < nJets37; n8++){
															if ((n8 != n1) && (n8 != n2) && (n8 != n3) && (n8 != n4) && (n8 != n5) && (n8 != n6) && (n8 != n7)){
																for (int n9 = 0; n9 < nJets37; n9++){
																	if ((n9 != n1) && (n9 != n2) && (n9 != n3) && (n9 != n4) && (n9 != n5) && (n9 != n6) && (n9 != n7) && (n9 != n8)){
																		for (int n10 = 0; n10 < nJets37; n10++){
																			if ((n10 != n1) && (n10 != n2) && (n10 != n3) && (n10 != n4) && (n10 != n5) && (n10 != n6) && (n10 != n7) && (n10 != n8) && (n10 != n9)){
																				if (deltaHT37 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et())){
																					deltaHT37 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et());
																				}
																				if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et())){
																					deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et());
																				}
																				if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et())){
																					deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et());
																				}
																				if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et())){
																					deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et());
																				}
																				if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() + jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et())){
																					deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() + jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et());
																				}
																			}
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
  if (nJets37 == 11){
	  for (int n1 = 0; n1 < nJets37; n1++){
		  for (int n2 = 0; n2 < nJets37; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets37; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets37; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets37; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets37; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												for (int n7 = 0; n7 < nJets37; n7++){
													if ((n7 != n1) && (n7 != n2) && (n7 != n3) && (n7 != n4) && (n7 != n5) && (n7 != n6)){
														for (int n8 = 0; n8 < nJets37; n8++){
															if ((n8 != n1) && (n8 != n2) && (n8 != n3) && (n8 != n4) && (n8 != n5) && (n8 != n6) && (n8 != n7)){
																for (int n9 = 0; n9 < nJets37; n9++){
																	if ((n9 != n1) && (n9 != n2) && (n9 != n3) && (n9 != n4) && (n9 != n5) && (n9 != n6) && (n9 != n7) && (n9 != n8)){
																		for (int n10 = 0; n10 < nJets37; n10++){
																			if ((n10 != n1) && (n10 != n2) && (n10 != n3) && (n10 != n4) && (n10 != n5) && (n10 != n6) && (n10 != n7) && (n10 != n8) && (n10 != n9)){
																				for (int n11 = 0; n11 < nJets37; n11++){
																					if ((n11 != n1) && (n11 != n2) && (n11 != n3) && (n11 != n4) && (n11 != n5) && (n11 != n6) && (n11 != n7) && (n11 != n8) && (n11 != n9) && (n11 != n10)){
																						if (deltaHT37 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et())){
																							deltaHT37 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et());
																						}
																						if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et())){
																							deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et());
																						}
																						if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et())){
																							deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et());
																						}
																						if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et())){
																							deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et());
																						}
																						if (deltaHT37 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() + jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et())){
																							deltaHT37 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() + jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et());
																						}
																					}
																				}
																			}
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

  if (nJets37 > 1){
	  floatEventProperties["deltaHt37"] = deltaHT37;
	  floatEventProperties["alphaT37"] = 1./2. * (1 - floatEventProperties["deltaHt37"] / floatEventProperties["ht37"]) / sqrt( 1. - (floatEventProperties["absvHt37"] / floatEventProperties["ht37"]) * (floatEventProperties["absvHt37"] / floatEventProperties["ht37"]));
  }
  else{
	  floatEventProperties["deltaHt37"] = 0.;
	  floatEventProperties["alphaT37"] = 0.;
  }

  float deltaHT43 = 0.;
  if (nJets43 > 1){
	  deltaHT43 = 1000.;
  }
  if (nJets43 == 2){
	  deltaHT43 = abs(jet1Vector43.Et() - jet2Vector43.Et());
  }
  if (nJets43 == 3){
	  for (int n1 = 0; n1 < nJets43; n1++){
		  for (int n2 = 0; n2 < nJets43; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets43; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						  if (deltaHT43 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et())){
							  deltaHT43 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et());
						  }
					  }
				  }
			  }
		  }
	  }
  }
  if (nJets43 == 4){
	  for (int n1 = 0; n1 < nJets43; n1++){
		  for (int n2 = 0; n2 < nJets43; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets43; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets43; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								if (deltaHT43 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et())){
									deltaHT43 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et());
								}
								if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et())){
									deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et());
								}
								
							}
						}
					}
				}
			}
		}
	 }
  }
  if (nJets43 == 5){
	  for (int n1 = 0; n1 < nJets43; n1++){
		  for (int n2 = 0; n2 < nJets43; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets43; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets43; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets43; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										if (deltaHT43 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et())){
											deltaHT43 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et());
										}
										if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et())){
											deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et());
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

  if (nJets43 == 6){
	  for (int n1 = 0; n1 < nJets43; n1++){
		  for (int n2 = 0; n2 < nJets43; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets43; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets43; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets43; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets43; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												if (deltaHT43 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et())){
													deltaHT43 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et());
												}
												if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et())){
													deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et());
												}
												if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et())){
													deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et());
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

 if (nJets43 == 7){
	  for (int n1 = 0; n1 < nJets43; n1++){
		  for (int n2 = 0; n2 < nJets43; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets43; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets43; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets43; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets43; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												for (int n7 = 0; n7 < nJets43; n7++){
													if ((n7 != n1) && (n7 != n2) && (n7 != n3) && (n7 != n4) && (n7 != n5) && (n7 != n6)){
														if (deltaHT43 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et())){
															deltaHT43 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et());
														}
														if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et())){
															deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et());
														}
														if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et())){
															deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et());
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
  if (nJets43 == 8){
	  for (int n1 = 0; n1 < nJets43; n1++){
		  for (int n2 = 0; n2 < nJets43; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets43; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets43; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets43; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets43; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												for (int n7 = 0; n7 < nJets43; n7++){
													if ((n7 != n1) && (n7 != n2) && (n7 != n3) && (n7 != n4) && (n7 != n5) && (n7 != n6)){
														for (int n8 = 0; n8 < nJets43; n8++){
															if ((n8 != n1) && (n8 != n2) && (n8 != n3) && (n8 != n4) && (n8 != n5) && (n8 != n6) && (n8 != n7)){
																if (deltaHT43 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et())){
																	deltaHT43 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et());
																}
																if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et())){
																	deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et());
																}
																if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et())){
																	deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et());
																}
																if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et())){
																	deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et());
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
  if (nJets43 == 9){
	  for (int n1 = 0; n1 < nJets43; n1++){
		  for (int n2 = 0; n2 < nJets43; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets43; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets43; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets43; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets43; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												for (int n7 = 0; n7 < nJets43; n7++){
													if ((n7 != n1) && (n7 != n2) && (n7 != n3) && (n7 != n4) && (n7 != n5) && (n7 != n6)){
														for (int n8 = 0; n8 < nJets43; n8++){
															if ((n8 != n1) && (n8 != n2) && (n8 != n3) && (n8 != n4) && (n8 != n5) && (n8 != n6) && (n8 != n7)){
																for (int n9 = 0; n9 < nJets43; n9++){
																	if ((n9 != n1) && (n9 != n2) && (n9 != n3) && (n9 != n4) && (n9 != n5) && (n9 != n6) && (n9 != n7) && (n9 != n8)){
																		if (deltaHT43 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et())){
																			deltaHT43 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et());
																		}
																		if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et())){
																			deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et());
																		}
																		if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et())){
																			deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et());
																		}
																		if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et())){
																			deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et());
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
  if (nJets43 == 10){
	  for (int n1 = 0; n1 < nJets43; n1++){
		  for (int n2 = 0; n2 < nJets43; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets43; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets43; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets43; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets43; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												for (int n7 = 0; n7 < nJets43; n7++){
													if ((n7 != n1) && (n7 != n2) && (n7 != n3) && (n7 != n4) && (n7 != n5) && (n7 != n6)){
														for (int n8 = 0; n8 < nJets43; n8++){
															if ((n8 != n1) && (n8 != n2) && (n8 != n3) && (n8 != n4) && (n8 != n5) && (n8 != n6) && (n8 != n7)){
																for (int n9 = 0; n9 < nJets43; n9++){
																	if ((n9 != n1) && (n9 != n2) && (n9 != n3) && (n9 != n4) && (n9 != n5) && (n9 != n6) && (n9 != n7) && (n9 != n8)){
																		for (int n10 = 0; n10 < nJets43; n10++){
																			if ((n10 != n1) && (n10 != n2) && (n10 != n3) && (n10 != n4) && (n10 != n5) && (n10 != n6) && (n10 != n7) && (n10 != n8) && (n10 != n9)){
																				if (deltaHT43 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et())){
																					deltaHT43 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et());
																				}
																				if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et())){
																					deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et());
																				}
																				if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et())){
																					deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et());
																				}
																				if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et())){
																					deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et());
																				}
																				if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() + jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et())){
																					deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() + jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et());
																				}
																			}
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
  if (nJets43 == 11){
	  for (int n1 = 0; n1 < nJets43; n1++){
		  for (int n2 = 0; n2 < nJets43; n2++){
			  if (n2 != n1){
				  for (int n3 = 0; n3 < nJets43; n3++){
					  if ((n3 != n1) && (n3 != n2)){
						for (int n4 = 0; n4 < nJets43; n4++){
							if ((n4 != n1) && (n4 != n2) && (n4 != n3)){
								for (int n5 = 0; n5 < nJets43; n5++){
									if ((n5 != n1) && (n5 != n2) && (n5 != n3) && (n5 != n4)){
										for (int n6 = 0; n6 < nJets43; n6++){
											if ((n6 != n1) && (n6 != n2) && (n6 != n3) && (n6 != n4) && (n6 != n5)){
												for (int n7 = 0; n7 < nJets43; n7++){
													if ((n7 != n1) && (n7 != n2) && (n7 != n3) && (n7 != n4) && (n7 != n5) && (n7 != n6)){
														for (int n8 = 0; n8 < nJets43; n8++){
															if ((n8 != n1) && (n8 != n2) && (n8 != n3) && (n8 != n4) && (n8 != n5) && (n8 != n6) && (n8 != n7)){
																for (int n9 = 0; n9 < nJets43; n9++){
																	if ((n9 != n1) && (n9 != n2) && (n9 != n3) && (n9 != n4) && (n9 != n5) && (n9 != n6) && (n9 != n7) && (n9 != n8)){
																		for (int n10 = 0; n10 < nJets43; n10++){
																			if ((n10 != n1) && (n10 != n2) && (n10 != n3) && (n10 != n4) && (n10 != n5) && (n10 != n6) && (n10 != n7) && (n10 != n8) && (n10 != n9)){
																				for (int n11 = 0; n11 < nJets43; n11++){
																					if ((n11 != n1) && (n11 != n2) && (n11 != n3) && (n11 != n4) && (n11 != n5) && (n11 != n6) && (n11 != n7) && (n11 != n8) && (n11 != n9) && (n11 != n10)){
																						if (deltaHT43 > abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et())){
																							deltaHT43 = abs(jets_->at(n1).et() - jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et());
																						}
																						if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et())){
																							deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() - jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et());
																						}
																						if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et())){
																							deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() - jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et());
																						}
																						if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et())){
																							deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() - jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et());
																						}
																						if (deltaHT43 > abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() + jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et())){
																							deltaHT43 = abs(jets_->at(n1).et() + jets_->at(n2).et() + jets_->at(n3).et() + jets_->at(n4).et() + jets_->at(n5).et() - jets_->at(n6).et() - jets_->at(n7).et() - jets_->at(n8).et() - jets_->at(n9).et() - jets_->at(n10).et() - jets_->at(n11).et());
																						}
																					}
																				}
																			}
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

  if (nJets43 > 1){
	  floatEventProperties["deltaHt43"] = deltaHT43;
	  floatEventProperties["alphaT43"] = 1./2. * (1 - floatEventProperties["deltaHt43"] / floatEventProperties["ht43"]) / sqrt( 1. - (floatEventProperties["absvHt43"] / floatEventProperties["ht43"]) * (floatEventProperties["absvHt43"] / floatEventProperties["ht43"]));
  }
  else{
	  floatEventProperties["deltaHt43"] = 0.;
	  floatEventProperties["alphaT43"] = 0.;
  }

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
  //~ floatEventProperties["deltaPhiJetMET"] = leadingJetMomentum.DeltaPhi(metVector);
  //~ floatEventProperties["deltaPhiSecondJetMET"] = 999.;
  //~ floatEventProperties["deltaPhiSecondJetMET"] = subLeadingJetMomentum.DeltaPhi(metVector);

  // Jet pt
  // bjet pt
  
  int nBJets = 0;

  TLorentzVector bJet1Vector(0.,0.,0.,0.);
  TLorentzVector bJet2Vector(0.,0.,0.,0.);
  TLorentzVector bJet3Vector(0.,0.,0.,0.);
  TLorentzVector bJet4Vector(0.,0.,0.,0.); 
  TLorentzVector bJet5Vector(0.,0.,0.,0.); 
  TLorentzVector bJet6Vector(0.,0.,0.,0.); 


  if ((bJets->size() > 0) && (bJets->at(0).et() > 50)){
	nBJets++;
    bJet1Vector.SetPxPyPzE(bJets->at(0).px(),bJets->at(0).py(),bJets->at(0).pz(),bJets->at(0).energy());
  }
  if ((bJets->size() > 1) && (bJets->at(1).et() > 50)){
	nBJets++;
    bJet2Vector.SetPxPyPzE(bJets->at(1).px(),bJets->at(1).py(),bJets->at(1).pz(),bJets->at(1).energy());
  }
  if ((bJets->size() > 2) && (bJets->at(2).et() > 50)){
	nBJets++;
    bJet3Vector.SetPxPyPzE(bJets->at(2).px(),bJets->at(2).py(),bJets->at(2).pz(),bJets->at(2).energy());
  }
  if ((bJets->size() > 3) && (bJets->at(3).et() > 50)){
	nBJets++;
    bJet4Vector.SetPxPyPzE(bJets->at(3).px(),bJets->at(3).py(),bJets->at(3).pz(),bJets->at(3).energy());
  }
  if ((bJets->size() > 4) && (bJets->at(4).et() > 50)){
	nBJets++;
    bJet5Vector.SetPxPyPzE(bJets->at(4).px(),bJets->at(4).py(),bJets->at(4).pz(),bJets->at(4).energy());
  }
  if ((bJets->size() > 5) && (bJets->at(5).et() > 50)){
	nBJets++;
    bJet6Vector.SetPxPyPzE(bJets->at(5).px(),bJets->at(5).py(),bJets->at(5).pz(),bJets->at(5).energy());
  }
  
  intEventProperties["nBJets"] = nBJets;
  tLorentzVectorEventProperties["bJet1"] = bJet1Vector;
  tLorentzVectorEventProperties["bJet2"] = bJet2Vector;
  tLorentzVectorEventProperties["bJet3"] = bJet3Vector;
  tLorentzVectorEventProperties["bJet4"] = bJet4Vector;
  tLorentzVectorEventProperties["bJet5"] = bJet5Vector;
  tLorentzVectorEventProperties["bJet6"] = bJet6Vector;

  TLorentzVector leadingbjet;
  TLorentzVector subleadingbjet;
  TLorentzVector bjetsum;
  floatEventProperties["mbb"] = -1.0;
  if (bJets->size() > 1) {

	for(std::vector<pat::Jet>::const_iterator it = bJets->begin(); it != bJets->end() ; ++it){

		if((*it).pt() > leadingbjet.Pt()){
		leadingbjet.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());
		}
		else if ((*it).pt() < leadingbjet.Pt() && (*it).pt() > subleadingbjet.Pt() ){
		subleadingbjet.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());
		}

	}

	bjetsum = leadingbjet+subleadingbjet;
	floatEventProperties["mbb"] = bjetsum.M();		
	
	
	

   }

  int nBJets37 = 0;   
  TLorentzVector bJet1Vector37(0.,0.,0.,0.);
  TLorentzVector bJet2Vector37(0.,0.,0.,0.);
  TLorentzVector bJet3Vector37(0.,0.,0.,0.);
  TLorentzVector bJet4Vector37(0.,0.,0.,0.); 
  TLorentzVector bJet5Vector37(0.,0.,0.,0.); 
  TLorentzVector bJet6Vector37(0.,0.,0.,0.); 


  if ((bJets->size() > 0) && (bJets->at(0).et() > 37)){
	nBJets37++;
    bJet1Vector37.SetPxPyPzE(bJets->at(0).px(),bJets->at(0).py(),bJets->at(0).pz(),bJets->at(0).energy());
  }
  if ((bJets->size() > 1) && (bJets->at(1).et() > 37)){
	nBJets37++;
    bJet2Vector37.SetPxPyPzE(bJets->at(1).px(),bJets->at(1).py(),bJets->at(1).pz(),bJets->at(1).energy());
  }
  if ((bJets->size() > 2) && (bJets->at(2).et() > 37)){
	nBJets37++;
    bJet3Vector37.SetPxPyPzE(bJets->at(2).px(),bJets->at(2).py(),bJets->at(2).pz(),bJets->at(2).energy());
  }
  if ((bJets->size() > 3) && (bJets->at(3).et() > 37)){
	nBJets37++;
    bJet4Vector37.SetPxPyPzE(bJets->at(3).px(),bJets->at(3).py(),bJets->at(3).pz(),bJets->at(3).energy());
  }
  if ((bJets->size() > 4) && (bJets->at(4).et() > 37)){
	nBJets37++;
    bJet4Vector37.SetPxPyPzE(bJets->at(4).px(),bJets->at(4).py(),bJets->at(4).pz(),bJets->at(4).energy());
  }
  if ((bJets->size() > 5) && (bJets->at(5).et() > 37)){
	nBJets37++;
    bJet5Vector37.SetPxPyPzE(bJets->at(5).px(),bJets->at(5).py(),bJets->at(5).pz(),bJets->at(5).energy());
  }

  intEventProperties["nBJets37"] = nBJets37;
  tLorentzVectorEventProperties["bJet1_37"] = bJet1Vector37;
  tLorentzVectorEventProperties["bJet2_37"] = bJet2Vector37;
  tLorentzVectorEventProperties["bJet3_37"] = bJet3Vector37;
  tLorentzVectorEventProperties["bJet4_37"] = bJet4Vector37;
  tLorentzVectorEventProperties["bJet5_37"] = bJet5Vector37;
  tLorentzVectorEventProperties["bJet6_37"] = bJet6Vector37;
  
  int nBJets43 = 0;
  TLorentzVector bJet1Vector43(0.,0.,0.,0.);
  TLorentzVector bJet2Vector43(0.,0.,0.,0.);
  TLorentzVector bJet3Vector43(0.,0.,0.,0.);
  TLorentzVector bJet4Vector43(0.,0.,0.,0.); 
  TLorentzVector bJet5Vector43(0.,0.,0.,0.); 
  TLorentzVector bJet6Vector43(0.,0.,0.,0.); 


  if ((bJets->size() > 0) && (bJets->at(0).et() > 43)){  
	nBJets43++;
    bJet1Vector43.SetPxPyPzE(bJets->at(0).px(),bJets->at(0).py(),bJets->at(0).pz(),bJets->at(0).energy());
  }
  if ((bJets->size() > 1) && (bJets->at(1).et() > 43)){  
	nBJets43++;
    bJet2Vector43.SetPxPyPzE(bJets->at(1).px(),bJets->at(1).py(),bJets->at(1).pz(),bJets->at(1).energy());
  }
  if ((bJets->size() > 2) && (bJets->at(2).et() > 43)){  
	nBJets43++;
    bJet3Vector43.SetPxPyPzE(bJets->at(2).px(),bJets->at(2).py(),bJets->at(2).pz(),bJets->at(2).energy());
  }
  if ((bJets->size() > 3) && (bJets->at(3).et() > 43)){  
	nBJets43++,
    bJet4Vector43.SetPxPyPzE(bJets->at(3).px(),bJets->at(3).py(),bJets->at(3).pz(),bJets->at(3).energy());
  }
  if ((bJets->size() > 4) && (bJets->at(4).et() > 43)){  
	nBJets43++;
    bJet4Vector43.SetPxPyPzE(bJets->at(4).px(),bJets->at(4).py(),bJets->at(4).pz(),bJets->at(4).energy());
  }
  if ((bJets->size() > 5) && (bJets->at(5).et() > 43)){  
	nBJets43++;
    bJet5Vector43.SetPxPyPzE(bJets->at(5).px(),bJets->at(5).py(),bJets->at(5).pz(),bJets->at(5).energy());
  }

  intEventProperties["nBJets43"] = nBJets43;
  tLorentzVectorEventProperties["bJet1_43"] = bJet1Vector43;
  tLorentzVectorEventProperties["bJet2_43"] = bJet2Vector43;
  tLorentzVectorEventProperties["bJet3_43"] = bJet3Vector43;
  tLorentzVectorEventProperties["bJet4_43"] = bJet4Vector43;
  tLorentzVectorEventProperties["bJet5_43"] = bJet5Vector43;
  tLorentzVectorEventProperties["bJet6_43"] = bJet6Vector43;
   

  floatEventProperties["weight"] = fctVtxWeight_( iEvent );
  floatEventProperties["weightBlockA"] = fctVtxWeightBlockA_( iEvent );
  floatEventProperties["weightBlockB"] = fctVtxWeightBlockB_( iEvent );
  floatEventProperties["weightUp"] = fctVtxWeightUp_( iEvent );
  floatEventProperties["weightBlockAUp"] = fctVtxWeightBlockAUp_( iEvent );
  floatEventProperties["weightBlockBUp"] = fctVtxWeightBlockBUp_( iEvent );
  floatEventProperties["weightDown"] = fctVtxWeightDown_( iEvent );
  floatEventProperties["weightBlockADown"] = fctVtxWeightBlockADown_( iEvent );
  floatEventProperties["weightBlockBDown"] = fctVtxWeightBlockBDown_( iEvent );
  
  if(useJets2_) {
    edm::Handle< std::vector< pat::Jet > > jets2;
    iEvent.getByLabel(jet2Tag_, jets2);
    intEventProperties["nJets2"] = jets2->size();
    floatEventProperties["ht2"] = 0.0;
    for(std::vector<pat::Jet>::const_iterator it = jets2->begin(); it != jets2->end() ; ++it){
      floatEventProperties["ht2"] += (*it).pt();
    }
  }

  //~ makeCombinations< pat::Electron >("EE", *electrons,*pfCands, iEvent, met,MHT, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  //~ makeCombinations< pat::Electron, pat::Muon >("EMu", *electrons, *muons,*pfCands, iEvent, met,MHT, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  //~ makeCombinations< pat::Muon >("MuMu", *muons,*pfCands, iEvent, met,MHT, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  //~ if(useTaus_){
    //~ makeCombinations< pat::Electron, pat::Tau >("ETau", *electrons, *taus,*pfCands, iEvent, met,MHT, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
    //~ makeCombinations< pat::Muon, pat::Tau>("MuTau", *muons, *taus,*pfCands, iEvent, met,MHT, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
    //~ makeCombinations< pat::Tau >("TauTau", *taus,*pfCands, iEvent, met,MHT, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  //~ }
  //  if( nMu != 2) std::cout << "-------! "<<nMu<<std::endl;
  
  const std::string treeName = "HadronicAlphaT";
  
  for(std::map<std::string, int>::const_iterator it = intEventProperties.begin(); it != intEventProperties.end(); ++it){
	//~ std::cout << (*it).first<< endl;
    assert(intBranches_[treeName].find((*it).first) != intBranches_[treeName].end());
    *(intBranches_[treeName][(*it).first]) = (*it).second;
  }
  for(std::map<std::string, float>::const_iterator it = floatEventProperties.begin(); it != floatEventProperties.end(); ++it){
	//~ std::cout << (*it).first<< endl;
    assert(floatBranches_[treeName].find((*it).first) != floatBranches_[treeName].end());
    *(floatBranches_[treeName][(*it).first]) = (*it).second;
  }

  for(std::map<std::string, TLorentzVector>::const_iterator it = tLorentzVectorEventProperties.begin(); it != tLorentzVectorEventProperties.end(); ++it){
	//~ std::cout << (*it).first<< endl;
    assert(tLorentzVectorBranches_[treeName].find((*it).first) != tLorentzVectorBranches_[treeName].end());
    *(tLorentzVectorBranches_[treeName][(*it).first]) = (*it).second;
  }
  
  trees_[treeName]->Fill();
	  

  delete jecUnc;
  delete shiftedJetsJESUp;
  delete shiftedJetsJESDown;


}

template <class aT, class bT> void 
HadronicAlphaTTrees::makeCombinations ( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT> &b,const std::vector<reco::PFCandidate>&pfCands, const edm::Event &ev, const pat::MET &patMet, const TLorentzVector &MHT, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties)
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
      fillTree<aT,bT>( treeName, *itA, *itB,pfCands, patMet,MHT); 
    }
  }
}

template <class aT> void 
HadronicAlphaTTrees::makeCombinations ( const std::string &treeName, const std::vector<aT> &a,const std::vector<reco::PFCandidate>&pfCands, const edm::Event &ev, const pat::MET &patMet , const TLorentzVector &MHT, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties)
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
      fillTree<aT, aT>( treeName, *itA, *itB,pfCands, patMet,MHT); 
    }
  }
}
  
//~ template <class aT> void 
//~ HadronicAlphaTTrees::fillLepton( const std::string &treeName, const aT& a, const std::string &leptonFlavorNumber,const  std::map<std::string, float> &floatEventProperties)
//~ {
  //~ if(debug) std::cout << treeName << "- pts:"<< a.pt();

  //~ std::cout << treeName << " charge"+leptonFlavorNumber+": " << a.charge() << endl;
  //~ std::cout << treeName << " id"+leptonFlavorNumber+": " << getId(a) << endl;
  //~ std::cout << treeName << " dB"+leptonFlavorNumber+": " << getDeltaB(a) << endl;

  //~ *(floatBranches_[treeName]["charge"+leptonFlavorNumber]) = a.charge();
  //~ *(floatBranches_[treeName]["id"+leptonFlavorNumber]) = getId(a);
  //~ *(floatBranches_[treeName]["dB"+leptonFlavorNumber]) = getDeltaB(a);
//~ 
//~ 
  //~ TLorentzVector genLepton;
  //~ if(a.genLepton() != NULL){
      //~ genLepton.SetPxPyPzE(a.genLepton()->px(),a.genLepton()->py(),a.genLepton()->pz(),a.genLepton()->energy());
  //~ }
//~ 
  //~ *(tLorentzVectorBranches_[treeName]["gen"+leptonFlavorNumber]) = genLepton;
//~ }

template <class aT, class bT> void 
HadronicAlphaTTrees::fillTree( const std::string &treeName, const aT& a, const bT& b,const std::vector<reco::PFCandidate>&pfCands, const pat::MET &patMet, const TLorentzVector &MHT)
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

  *(floatBranches_[treeName]["phiAtECALEntrance1"]) = getPhiAtECAL(a,pfCands);
  *(floatBranches_[treeName]["phiAtECALEntrance2"]) = getPhiAtECAL(b,pfCands); 

  TLorentzVector comb = aVec+bVec;
  TLorentzVector MHT2 = comb + MHT;
  std::pair<double, double> pZeta = calcPZeta(a.p(), b.p(), met);
  *(floatBranches_[treeName]["chargeProduct"]) = a.charge()*b.charge();
  *(tLorentzVectorBranches_[treeName]["p4"]) = comb;
  *(tLorentzVectorBranches_[treeName]["lepton1"]) = aVec;
  *(tLorentzVectorBranches_[treeName]["lepton2"]) = bVec;
  *(floatBranches_[treeName]["mht"]) = MHT2.Pt();
  *(floatBranches_[treeName]["pt1"]) = aVec.Pt();
  *(floatBranches_[treeName]["pt2"]) = bVec.Pt();
  *(floatBranches_[treeName]["charge1"]) = a.charge();
  *(floatBranches_[treeName]["charge2"]) = b.charge();
  *(floatBranches_[treeName]["eta1"]) = aVec.Eta();
  *(floatBranches_[treeName]["eta2"]) = bVec.Eta();
  *(floatBranches_[treeName]["id1"]) = getId(a);
  *(floatBranches_[treeName]["id2"]) = getId(b);
  *(floatBranches_[treeName]["dB1"]) = getDeltaB(a);
  *(floatBranches_[treeName]["dB2"]) = getDeltaB(b);
  *(floatBranches_[treeName]["mt1"]) = transverseMass(aVec, met);
  *(floatBranches_[treeName]["mt2"]) = transverseMass(bVec, met);
  *(floatBranches_[treeName]["eff1"]) = efficiencies_(a);
  *(floatBranches_[treeName]["eff2"]) = efficiencies_(b);
  *(floatBranches_[treeName]["fakeWeight1"]) = fakeRates_(a);
  *(floatBranches_[treeName]["fakeWeight2"]) = fakeRates_(b);
  *(floatBranches_[treeName]["deltaPhi"]) = aVec.DeltaPhi( bVec );
  *(floatBranches_[treeName]["deltaPhiLeptonMET1"]) = aVec.DeltaPhi( met ); 
  *(floatBranches_[treeName]["deltaPhiLeptonMET2"]) = bVec.DeltaPhi( met );
  *(floatBranches_[treeName]["deltaR"]) = aVec.DeltaR( bVec );
  *(floatBranches_[treeName]["jzb"]) = (uncorrectedMet+comb).Pt() - comb.Pt();
  *(floatBranches_[treeName]["pZeta"]) = pZeta.first;
  *(floatBranches_[treeName]["pZetaVis"]) = pZeta.second;

  if (jets_->size() >= 2){
    TLorentzVector vJet1 = TLorentzVector(jets_->at(0).p4().x(), jets_->at(0).p4().y(), jets_->at(0).p4().z(), jets_->at(0).p4().t());
    TLorentzVector vJet2 = TLorentzVector(jets_->at(1).p4().x(), jets_->at(1).p4().y(), jets_->at(1).p4().z(), jets_->at(1).p4().t());
    TLorentzVector sub = aVec + bVec + vJet1 + vJet2;
    *(floatBranches_[treeName]["sqrts"]) = std::sqrt(TMath::Power((std::sqrt(sub.M2() + sub.Perp2()) + met.Et()), 2.0) - (sub.Vect().XYvector() + met.Vect().XYvector())*(sub.Vect().XYvector() + met.Vect().XYvector()));
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
  *(tLorentzVectorBranches_[treeName]["genLepton1"]) = genLepton1;
  *(tLorentzVectorBranches_[treeName]["genLepton2"]) = genLepton2;
  if(debug) std::cout << ", matched = "<<matched<<", motherId = "<<aMother;
  if(debug) std::cout<<", M = "<< comb.M() <<", chargeProduct = "<< a.charge()*b.charge() <<std::endl;
  
  trees_[treeName]->Fill();
}

int 
HadronicAlphaTTrees::getMotherPdgId( const reco::GenParticle &p)
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
HadronicAlphaTTrees::getLeptonPdgId( const reco::GenParticle &p)
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
HadronicAlphaTTrees::calcPZeta(const TLorentzVector& p1,const TLorentzVector& p2, const TLorentzVector& met)
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

void HadronicAlphaTTrees::fillPdfUncert(const edm::Handle< std::vector<double> >& weightHandle, const std::string& pdfIdentifier, const std::string& treeName){
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

const TLorentzVector HadronicAlphaTTrees::getMomentum(const  pat::Electron &e)
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

const TLorentzVector HadronicAlphaTTrees::getMomentum(const  pat::Muon &mu)
{
  const TLorentzVector result = TLorentzVector(mu.px(), mu.py(), mu.pz(), mu.energy());
  return result;
}

const TLorentzVector HadronicAlphaTTrees::getMomentum(const  pat::Tau &tau)
{
  const TLorentzVector result = TLorentzVector(tau.px(), tau.py(), tau.pz(), tau.energy());
  return result;
}

float HadronicAlphaTTrees::getId(const  pat::Electron &e)
{
  //  if (e.isEE())
  //  return (e.dr03HcalTowerSumEt() + e.dr03EcalRecHitSumEt() + e.dr03TkSumPt())/e.pt();
  // else
  //  return (e.dr03HcalTowerSumEt() + std::max(0.0, e.dr03EcalRecHitSumEt() - 1.0) + e.dr03TkSumPt())/e.pt();

  //  std::cout<<"electron " << (lepton.trackIso()+lepton.ecalIso()+lepton.hcalIso())/lepton.pt() << std::endl;
  // return (lepton.trackIso()+lepton.ecalIso()+lepton.hcalIso())/lepton.pt() ;
    //  return (e.chargedHadronIso() + e.photonIso() + e.neutralHadronIso()) / e.pt();
  return fctIsolation_(e)* 1./e.pt();
}

float HadronicAlphaTTrees::getId(const  pat::Muon &mu)
{
  //  std::cout<<"muon " << (lepton.trackIso()+lepton.ecalIso()+lepton.hcalIso())/lepton.pt() << std::endl;
  //  return (mu.isolationR03().hadEt + mu.isolationR03().emEt + mu.isolationR03().sumPt) / mu.pt();
  //  return (mu.chargedHadronIso() + mu.photonIso() + mu.neutralHadronIso()) / mu.pt();
  return fctIsolation_(mu)* 1./mu.pt();
}

float HadronicAlphaTTrees::getId(const  pat::Tau &tau)
{
  float result = fctIsolation_(tau);
  if(tau.tauID(tauId_) < 0.5)
    result *= -1.0;
  return result;
}


float HadronicAlphaTTrees::getPhiAtECAL(const  pat::Tau &tau, const std::vector<reco::PFCandidate>&pfCands)
{

  float result = -99.0;
  return result;
}


float HadronicAlphaTTrees::getPhiAtECAL(const  pat::Muon &mu, const std::vector<reco::PFCandidate>&pfCands)
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

float HadronicAlphaTTrees::getPhiAtECAL(const  pat::Electron &e, const std::vector<reco::PFCandidate>&pfCands)
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


float HadronicAlphaTTrees::topPtWeightBen(double topPt){
  if( topPt<0 ) return 1;

  float p0 = 1.18246e+00;
  float p1 = 4.63312e+02;
  float p2 = 2.10061e-06;

  if( topPt>p1 ) topPt = p1;

  float result = p0 + p2 * topPt * ( topPt - 2 * p1 );
  return result;
}

float HadronicAlphaTTrees::topPtWeightTOP(double topPt){
  if( topPt<0 ) return 1;

  float p0 = 0.156;
  float p1 = 0.00137;

  float result = exp(p0 - p1 * topPt);
  return result;
}


float HadronicAlphaTTrees::getDeltaB(const  pat::Electron &e)
{
  float result = e.dB(pat::Electron::PV3D);
  return result;
}

float HadronicAlphaTTrees::getDeltaB(const  pat::Muon &mu)
{
  float result = mu.dB(pat::Muon::PV3D);
  return result;
}

float HadronicAlphaTTrees::getDeltaB(const  pat::Tau &tau)
{
  float result = -1; // not available for pat::Tau could use ip of leading ch. Hadr if needed.
  return result;
}


float HadronicAlphaTTrees::transverseMass(const TLorentzVector& p, const TLorentzVector& met)
{
  reco::Candidate::LorentzVector otherMet(met.Px(),met.Py(),met.Pz(),met.E());
  reco::Candidate::LorentzVector leptonT(p.Px(),p.Py(),0.,p.E()*sin(p.Theta()));
  reco::Candidate::LorentzVector sumT=leptonT+otherMet;

  return std::sqrt(sumT.M2());
}

std::string HadronicAlphaTTrees::convertInputTag(const edm::InputTag tag)
{
  std::string result = tag.label();
  if(tag.instance().length() > 0)
    result = tag.instance();
  //  std::cerr << "'"<<tag.label() << "', '"<< tag.instance()<<"' = '"<< result<<"'"<<std::endl;
  return result;
}

// ------------ Method called once each job just before starting event loop  ------------
void 
HadronicAlphaTTrees::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HadronicAlphaTTrees::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(HadronicAlphaTTrees);
