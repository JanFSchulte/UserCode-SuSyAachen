// -*- C++ -*-
//
// Package:    Histograms
// Class:      MultiLeptonTrees
// 
/**\class MultiLeptonTrees MultiLeptonTrees.cc brot/MultiLeptonTrees/src/MultiLeptonTrees.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  matthias edelhoff
//         Created:  Tue Oct 27 13:50:40 CET 2009
// $Id: MultiLeptonTrees.cc,v 1.31 2012/09/17 17:38:58 sprenger Exp $
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

class MultiLeptonTrees : public edm::EDAnalyzer {
public:
  explicit MultiLeptonTrees(const edm::ParameterSet&);
  ~MultiLeptonTrees();

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
MultiLeptonTrees::MultiLeptonTrees(const edm::ParameterSet& iConfig):
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
  trees_["MultiLepton"] = file->make<TTree>("MultiLeptonTree", "MultiLeptonTree");
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
  initFloatBranch( "chargeElectron1" );
  initFloatBranch( "chargeElectron2" );
  initFloatBranch( "chargeElectron3" );
  initFloatBranch( "chargeElectron4" );
  initFloatBranch( "chargeMuon1" );
  initFloatBranch( "chargeMuon2" );
  initFloatBranch( "chargeMuon3" );
  initFloatBranch( "chargeMuon4" );
  //~ initTLorentzVectorBranch( "p4" );
  initTLorentzVectorBranch( "vMet" );
  initTLorentzVectorBranch( "vGenMet" );
  //~ initTLorentzVectorBranch( "vGenMetNeutrinos" );
  //~ initTLorentzVectorBranch( "p4Gen" );
  //~ initTLorentzVectorBranch( "lepton1" );
  //~ initTLorentzVectorBranch( "lepton2" );
  //~ initTLorentzVectorBranch( "lepton3" );
  //~ initTLorentzVectorBranch( "lepton4" );
  initTLorentzVectorBranch( "electron1" );
  initTLorentzVectorBranch( "electron2" );
  initTLorentzVectorBranch( "electron3" );
  initTLorentzVectorBranch( "electron4" );
  initTLorentzVectorBranch( "muon1" );
  initTLorentzVectorBranch( "muon2" );
  initTLorentzVectorBranch( "muon3" );
  initTLorentzVectorBranch( "muon4" );
  //~ initTLorentzVectorBranch( "genLepton1" );
  //~ initTLorentzVectorBranch( "genLepton2" );
  //~ initTLorentzVectorBranch( "genLepton3" );
  //~ initTLorentzVectorBranch( "genLepton4" );
  //~ initTLorentzVectorBranch( "genElectron1" );
  //~ initTLorentzVectorBranch( "genElectron2" );
  //~ initTLorentzVectorBranch( "genElectron3" );
  //~ initTLorentzVectorBranch( "genElectron4" );
  //~ initTLorentzVectorBranch( "genMuon1" );
  //~ initTLorentzVectorBranch( "genMuon2" );
  //~ initTLorentzVectorBranch( "genMuon3" );
  //~ initTLorentzVectorBranch( "genMuon4" );
  //~ initTLorentzVectorBranch( "p4Electrons12" );
  //~ initTLorentzVectorBranch( "p4Electrons13" );
  //~ initTLorentzVectorBranch( "p4Electrons14" );
  //~ initTLorentzVectorBranch( "p4Electrons23" );
  //~ initTLorentzVectorBranch( "p4Electrons24" );
  //~ initTLorentzVectorBranch( "p4Electrons34" );
  //~ initTLorentzVectorBranch( "p4Muons12" );
  //~ initTLorentzVectorBranch( "p4Muons13" );
  //~ initTLorentzVectorBranch( "p4Muons14" );
  //~ initTLorentzVectorBranch( "p4Muons23" );
  //~ initTLorentzVectorBranch( "p4Muons24" );
  //~ initTLorentzVectorBranch( "p4Muons34" );
  //~ initTLorentzVectorBranch( "jet1" );
  //~ initTLorentzVectorBranch( "jet2" );
  //~ initTLorentzVectorBranch( "jet3" );
  //~ initTLorentzVectorBranch( "jet4" );
  //~ initTLorentzVectorBranch( "bJet1" );
  //~ initTLorentzVectorBranch( "bJet2" );
  //~ initTLorentzVectorBranch( "bJet3" );
  //~ initTLorentzVectorBranch( "bJet4" );
  initTLorentzVectorBranch( "vMetType1" );
  initFloatBranch( "Electrons12Mass" );
  initFloatBranch( "Electrons13Mass" );
  initFloatBranch( "Electrons14Mass" );
  initFloatBranch( "Electrons23Mass" );
  initFloatBranch( "Electrons24Mass" );
  initFloatBranch( "Electrons34Mass" );
  initFloatBranch( "Muons12Mass" );
  initFloatBranch( "Muons13Mass" );
  initFloatBranch( "Muons14Mass" );
  initFloatBranch( "Muons23Mass" );
  initFloatBranch( "Muons24Mass" );
  initFloatBranch( "Muons34Mass" );
  initFloatBranch( "rho" );
  //~ initFloatBranch( "pt1" );
  //~ initFloatBranch( "pt2" );
  //~ initFloatBranch( "phiAtECALEntrance1");
  //~ initFloatBranch( "phiAtECALEntrance2");
  //~ initFloatBranch( "eta1" );
  //~ initFloatBranch( "eta2" );
  initFloatBranch( "idElectron1" );
  initFloatBranch( "idElectron2" );
  initFloatBranch( "idElectron3" );
  initFloatBranch( "idElectron4" );
  initFloatBranch( "idMuon1" );
  initFloatBranch( "idMuon2" );
  initFloatBranch( "idMuon3" );
  initFloatBranch( "idMuon4" );
  initFloatBranch( "dBElectron1" );
  initFloatBranch( "dBElectron2" );
  initFloatBranch( "dBElectron3" );
  initFloatBranch( "dBElectron4" );
  initFloatBranch( "dBMuon1" );
  initFloatBranch( "dBMuon2" );
  initFloatBranch( "dBMuon3" );
  initFloatBranch( "dBMuon4" );
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
  initFloatBranch( "ht" );
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
  initIntBranch( "nBJets" );
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
  initIntBranch( "nOSSF" );
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
MultiLeptonTrees::initTLorentzVectorBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    tLorentzVectorBranches_[(*it).first][name] = new TLorentzVector;
    (*it).second->Branch(name.c_str(), "TLorentzVector" ,&tLorentzVectorBranches_[(*it).first][name]);
  }
}

void 
MultiLeptonTrees::initFloatBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    floatBranches_[(*it).first][name] = new float;
    (*it).second->Branch(name.c_str(), floatBranches_[(*it).first][name], (name+"/F").c_str());
  }
}

void 
MultiLeptonTrees::initIntBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    intBranches_[(*it).first][name] = new unsigned int;
    (*it).second->Branch(name.c_str(), intBranches_[(*it).first][name], (name+"/I").c_str());
  }
}

MultiLeptonTrees::~MultiLeptonTrees()
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
MultiLeptonTrees::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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


  intEventProperties["nBJets"] = bJets->size();
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


  int nGenLeptons = 0;
  int nGenLightLeptons = 0;
  int nGenTaus = 0;
  int nGenNeutrinos = 0;
  if (genParticles.isValid()){
	
	for (std::vector<reco::GenParticle>::const_iterator itGenParticle = genParticles->begin(); itGenParticle != genParticles->end(); itGenParticle++) {
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

  //~ TLorentzVector MHT;
  //~ TLorentzVector tempMHT;
  //~ TLorentzVector leadingJetMomentum;
  //~ TLorentzVector subLeadingJetMomentum;
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
  for(std::vector<pat::Jet>::const_iterator it = jets_->begin(); it != jets_->end() ; ++it){
	if ((*it).pt() >=30.0){
		nJets++;
		//~ tempMHT.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());	
		//~ MHT = MHT + tempMHT;
		floatEventProperties["ht"] += (*it).pt();
		//~ if((*it).pt() > leadingJetMomentum.Pt()){
		//~ leadingJetMomentum.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());
		//~ }
		//~ if((*it).pt() < leadingJetMomentum.Pt() && ((*it).pt() > subLeadingJetMomentum.Pt())){
		//~ subLeadingJetMomentum.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());
		//~ }
	}
  }
  intEventProperties["nJets"] = nJets;


  //~ TLorentzVector jet1Vector(0.,0.,0.,0.);
  //~ TLorentzVector jet2Vector(0.,0.,0.,0.);
  //~ TLorentzVector jet3Vector(0.,0.,0.,0.);
  //~ TLorentzVector jet4Vector(0.,0.,0.,0.); 


  //~ if (nJets > 0)
    //~ jet1Vector.SetPxPyPzE(jets_->at(0).px(),jets_->at(0).py(),jets_->at(0).pz(),jets_->at(0).energy());
  //~ if (nJets > 1)
    //~ jet2Vector.SetPxPyPzE(jets_->at(1).px(),jets_->at(1).py(),jets_->at(1).pz(),jets_->at(1).energy());
  //~ if (nJets > 2)
    //~ jet3Vector.SetPxPyPzE(jets_->at(2).px(),jets_->at(2).py(),jets_->at(2).pz(),jets_->at(2).energy());
  //~ if (nJets > 3)
    //~ jet4Vector.SetPxPyPzE(jets_->at(3).px(),jets_->at(3).py(),jets_->at(3).pz(),jets_->at(3).energy());
  //~ tLorentzVectorEventProperties["jet1"] = jet1Vector;
  //~ tLorentzVectorEventProperties["jet2"] = jet2Vector;
  //~ tLorentzVectorEventProperties["jet3"] = jet3Vector;
  //~ tLorentzVectorEventProperties["jet4"] = jet4Vector;








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

  //~ TLorentzVector bJet1Vector(0.,0.,0.,0.);
  //~ TLorentzVector bJet2Vector(0.,0.,0.,0.);
  //~ TLorentzVector bJet3Vector(0.,0.,0.,0.);
  //~ TLorentzVector bJet4Vector(0.,0.,0.,0.); 


  //~ if (bJets->size() > 0){
    //~ bJet1Vector.SetPxPyPzE(bJets->at(0).px(),bJets->at(0).py(),bJets->at(0).pz(),bJets->at(0).energy());
  //~ }
  //~ if (bJets->size() > 1){
    //~ bJet2Vector.SetPxPyPzE(bJets->at(1).px(),bJets->at(1).py(),bJets->at(1).pz(),bJets->at(1).energy());
  //~ }
  //~ if (bJets->size() > 2){
    //~ bJet3Vector.SetPxPyPzE(bJets->at(2).px(),bJets->at(2).py(),bJets->at(2).pz(),bJets->at(2).energy());
  //~ }
  //~ if (bJets->size() > 3){
    //~ bJet4Vector.SetPxPyPzE(bJets->at(3).px(),bJets->at(3).py(),bJets->at(3).pz(),bJets->at(3).energy());
  //~ }
//~ 
  //~ tLorentzVectorEventProperties["bJet1"] = bJet1Vector;
  //~ tLorentzVectorEventProperties["bJet2"] = bJet2Vector;
  //~ tLorentzVectorEventProperties["bJet3"] = bJet3Vector;
  //~ tLorentzVectorEventProperties["bJet4"] = bJet4Vector;

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
   
  for (int i = 0; i < 4; i++){
	  std::string number = static_cast<ostringstream*>( &(ostringstream() << i+1) )->str();
	  floatEventProperties["chargeElectron"+number] = 0.;
	  floatEventProperties["idElectron"+number] = 0.;
	  floatEventProperties["dBElectron"+number] = 0.;
	  floatEventProperties["chargeMuon"+number] = 0.;
	  floatEventProperties["idMuon"+number] = 0.;
	  floatEventProperties["dBMuon"+number] = 0.;
  }
  
  int n_OSSF_pairs = 0;
  
  TLorentzVector electron1Vector(0.,0.,0.,0.);
  TLorentzVector electron2Vector(0.,0.,0.,0.);
  TLorentzVector electron3Vector(0.,0.,0.,0.);
  TLorentzVector electron4Vector(0.,0.,0.,0.);
  
  //~ TLorentzVector genElectron1Vector(0.,0.,0.,0.);
  //~ TLorentzVector genElectron2Vector(0.,0.,0.,0.);
  //~ TLorentzVector genElectron3Vector(0.,0.,0.,0.);
  //~ TLorentzVector genElectron4Vector(0.,0.,0.,0.);
  
  TLorentzVector electronComb12Vector(0.,0.,0.,0.);  
  TLorentzVector electronComb13Vector(0.,0.,0.,0.);  
  TLorentzVector electronComb14Vector(0.,0.,0.,0.);  
  TLorentzVector electronComb23Vector(0.,0.,0.,0.);  
  TLorentzVector electronComb24Vector(0.,0.,0.,0.);  
  TLorentzVector electronComb34Vector(0.,0.,0.,0.); 

  if (electrons->size() > 0){
    electron1Vector.SetPxPyPzE(electrons->at(0).px(),electrons->at(0).py(),electrons->at(0).pz(),electrons->at(0).energy());
    //~ fillLepton< pat::Electron>( "MultiLepton", electrons->at(0), "Electron1", floatEventProperties);
    
    floatEventProperties["chargeElectron1"] = electrons->at(0).charge();
    floatEventProperties["idElectron1"] = getId(electrons->at(0));
    floatEventProperties["dBElectron1"] = getDeltaB(electrons->at(0));
    
    TLorentzVector genLepton(0.,0.,0.,0.);
	//~ if(electrons->at(0).genLepton() != NULL){
      //~ genElectron1Vector.SetPxPyPzE(electrons->at(0).genLepton()->px(),electrons->at(0).genLepton()->py(),electrons->at(0).genLepton()->pz(),electrons->at(0).genLepton()->energy());
    //~ }
  }
  if (electrons->size() > 1){
    electron2Vector.SetPxPyPzE(electrons->at(1).px(),electrons->at(1).py(),electrons->at(1).pz(),electrons->at(1).energy());
    //~ fillLepton< pat::Electron>( "MultiLepton", electrons->at(1), "Electron2", floatEventProperties);
    
    floatEventProperties["chargeElectron2"] = electrons->at(1).charge();
    floatEventProperties["idElectron2"] = getId(electrons->at(1));
    floatEventProperties["dBElectron2"] = getDeltaB(electrons->at(1));
    
    if (electrons->at(0).charge() * electrons->at(1).charge() == -1.){
		n_OSSF_pairs++;
		electronComb12Vector.SetPxPyPzE(electrons->at(0).px()+electrons->at(1).px(),electrons->at(0).py()+electrons->at(1).py(),electrons->at(0).pz()+electrons->at(1).pz(),electrons->at(0).energy()+electrons->at(1).energy());
	}
    
    TLorentzVector genLepton(0.,0.,0.,0.);
	//~ if(electrons->at(1).genLepton() != NULL){
      //~ genElectron2Vector.SetPxPyPzE(electrons->at(1).genLepton()->px(),electrons->at(1).genLepton()->py(),electrons->at(1).genLepton()->pz(),electrons->at(1).genLepton()->energy());
    //~ }
  }
  if (electrons->size() > 2){
    electron3Vector.SetPxPyPzE(electrons->at(2).px(),electrons->at(2).py(),electrons->at(2).pz(),electrons->at(2).energy());
    //~ fillLepton< pat::Electron>( "MultiLepton", electrons->at(2), "Electron3", floatEventProperties);
    
    floatEventProperties["chargeElectron3"] = electrons->at(2).charge();
    floatEventProperties["idElectron3"] = getId(electrons->at(2));
    floatEventProperties["dBElectron3"] = getDeltaB(electrons->at(2));
    
    if (electrons->at(0).charge()* electrons->at(2).charge()== -1.){
		n_OSSF_pairs++; 
		electronComb13Vector.SetPxPyPzE(electrons->at(0).px()+electrons->at(2).px(),electrons->at(0).py()+electrons->at(2).py(),electrons->at(0).pz()+electrons->at(2).pz(),electrons->at(0).energy()+electrons->at(2).energy());
	}
    
    if (electrons->at(1).charge()* electrons->at(2).charge()== -1.){
		n_OSSF_pairs++;
		electronComb23Vector.SetPxPyPzE(electrons->at(1).px()+electrons->at(2).px(),electrons->at(1).py()+electrons->at(2).py(),electrons->at(1).pz()+electrons->at(2).pz(),electrons->at(1).energy()+electrons->at(2).energy());
	}
    
    TLorentzVector genLepton(0.,0.,0.,0.);
	//~ if(electrons->at(2).genLepton() != NULL){
      //~ genElectron3Vector.SetPxPyPzE(electrons->at(2).genLepton()->px(),electrons->at(2).genLepton()->py(),electrons->at(2).genLepton()->pz(),electrons->at(2).genLepton()->energy());
    //~ }
  }
  if (electrons->size() > 3){
    electron4Vector.SetPxPyPzE(electrons->at(3).px(),electrons->at(3).py(),electrons->at(3).pz(),electrons->at(3).energy());
    //~ fillLepton< pat::Electron>( "MultiLepton", electrons->at(3), "Electron4", floatEventProperties);
    
    floatEventProperties["chargeElectron4"] = electrons->at(3).charge();
    floatEventProperties["idElectron4"] = getId(electrons->at(3));
    floatEventProperties["dBElectron4"] = getDeltaB(electrons->at(3));
    
     if (electrons->at(0).charge()* electrons->at(3).charge()== -1.){
		 n_OSSF_pairs++;
		 electronComb14Vector.SetPxPyPzE(electrons->at(0).px()+electrons->at(3).px(),electrons->at(0).py()+electrons->at(3).py(),electrons->at(0).pz()+electrons->at(3).pz(),electrons->at(0).energy()+electrons->at(3).energy());
	}
    
    if (electrons->at(1).charge()* electrons->at(3).charge()== -1.){
		n_OSSF_pairs++;
		electronComb24Vector.SetPxPyPzE(electrons->at(1).px()+electrons->at(3).px(),electrons->at(1).py()+electrons->at(3).py(),electrons->at(1).pz()+electrons->at(3).pz(),electrons->at(1).energy()+electrons->at(3).energy());
	}
    
    if (electrons->at(2).charge()* electrons->at(3).charge()== -1.){
		n_OSSF_pairs++;
		electronComb34Vector.SetPxPyPzE(electrons->at(2).px()+electrons->at(3).px(),electrons->at(2).py()+electrons->at(3).py(),electrons->at(2).pz()+electrons->at(3).pz(),electrons->at(2).energy()+electrons->at(3).energy());
	}
    
    TLorentzVector genLepton(0.,0.,0.,0.);
	//~ if(electrons->at(3).genLepton() != NULL){
      //~ genElectron4Vector.SetPxPyPzE(electrons->at(3).genLepton()->px(),electrons->at(3).genLepton()->py(),electrons->at(3).genLepton()->pz(),electrons->at(3).genLepton()->energy());
    //~ }
  }


  tLorentzVectorEventProperties["electron1"] = electron1Vector;
  tLorentzVectorEventProperties["electron2"] = electron2Vector;
  tLorentzVectorEventProperties["electron3"] = electron3Vector;
  tLorentzVectorEventProperties["electron4"] = electron4Vector;
  
  //~ tLorentzVectorEventProperties["genElectron1"] = genElectron1Vector;
  //~ tLorentzVectorEventProperties["genElectron2"] = genElectron2Vector;
  //~ tLorentzVectorEventProperties["genElectron3"] = genElectron3Vector;
  //~ tLorentzVectorEventProperties["genElectron4"] = genElectron4Vector;
  
  //~ tLorentzVectorEventProperties["p4Electrons12"] = electronComb12Vector;
  //~ tLorentzVectorEventProperties["p4Electrons13"] = electronComb13Vector;
  //~ tLorentzVectorEventProperties["p4Electrons14"] = electronComb14Vector;
  //~ tLorentzVectorEventProperties["p4Electrons23"] = electronComb23Vector;
  //~ tLorentzVectorEventProperties["p4Electrons24"] = electronComb24Vector;
  //~ tLorentzVectorEventProperties["p4Electrons34"] = electronComb34Vector;

  floatEventProperties["Electrons12Mass"] = electronComb12Vector.M();
  floatEventProperties["Electrons13Mass"] = electronComb13Vector.M();
  floatEventProperties["Electrons14Mass"] = electronComb14Vector.M();
  floatEventProperties["Electrons23Mass"] = electronComb23Vector.M();
  floatEventProperties["Electrons24Mass"] = electronComb24Vector.M();
  floatEventProperties["Electrons34Mass"] = electronComb34Vector.M();
    
  TLorentzVector muon1Vector(0.,0.,0.,0.);
  TLorentzVector muon2Vector(0.,0.,0.,0.);
  TLorentzVector muon3Vector(0.,0.,0.,0.);
  TLorentzVector muon4Vector(0.,0.,0.,0.); 
    
  //~ TLorentzVector genMuon1Vector(0.,0.,0.,0.);
  //~ TLorentzVector genMuon2Vector(0.,0.,0.,0.);
  //~ TLorentzVector genMuon3Vector(0.,0.,0.,0.);
  //~ TLorentzVector genMuon4Vector(0.,0.,0.,0.);
  
  TLorentzVector muonComb12Vector(0.,0.,0.,0.);  
  TLorentzVector muonComb13Vector(0.,0.,0.,0.);  
  TLorentzVector muonComb14Vector(0.,0.,0.,0.);  
  TLorentzVector muonComb23Vector(0.,0.,0.,0.);  
  TLorentzVector muonComb24Vector(0.,0.,0.,0.);  
  TLorentzVector muonComb34Vector(0.,0.,0.,0.);  
  

  if (muons->size() > 0){
    muon1Vector.SetPxPyPzE(muons->at(0).px(),muons->at(0).py(),muons->at(0).pz(),muons->at(0).energy());
    //~ fillLepton< pat::Muon>( "MultiLepton", muons->at(0), "Muon1", floatEventProperties);
    
    floatEventProperties["chargeMuon1"] = muons->at(0).charge();
    floatEventProperties["idMuon1"] = getId(muons->at(0));
    floatEventProperties["dBMuon1"] = getDeltaB(muons->at(0));
    
	//~ if(muons->at(0).genLepton() != NULL){
      //~ genMuon1Vector.SetPxPyPzE(muons->at(0).genLepton()->px(),muons->at(0).genLepton()->py(),muons->at(0).genLepton()->pz(),muons->at(0).genLepton()->energy());
    //~ }
  }
  if (muons->size() > 1){
    muon2Vector.SetPxPyPzE(muons->at(1).px(),muons->at(1).py(),muons->at(1).pz(),muons->at(1).energy());
    //~ fillLepton< pat::Muon>( "MultiLepton", muons->at(1), "Muon2", floatEventProperties);
    
    floatEventProperties["chargeMuon2"] = muons->at(1).charge();
    floatEventProperties["idMuon2"] = getId(muons->at(1));
    floatEventProperties["dBMuon2"] = getDeltaB(muons->at(1));
    
	if (muons->at(0).charge()* muons->at(1).charge()== -1.){
		n_OSSF_pairs++;
		muonComb12Vector.SetPxPyPzE(muons->at(0).px()+muons->at(1).px(),muons->at(0).py()+muons->at(1).py(),muons->at(0).pz()+muons->at(1).pz(),muons->at(0).energy()+muons->at(1).energy());
	}
    
	//~ if(muons->at(1).genLepton() != NULL){
      //~ genMuon2Vector.SetPxPyPzE(muons->at(1).genLepton()->px(),muons->at(1).genLepton()->py(),muons->at(1).genLepton()->pz(),muons->at(1).genLepton()->energy());
    //~ }
  }
  if (muons->size() > 2){
    muon3Vector.SetPxPyPzE(muons->at(2).px(),muons->at(2).py(),muons->at(2).pz(),muons->at(2).energy());
    //~ fillLepton< pat::Muon>( "MultiLepton", muons->at(2), "Muon3", floatEventProperties);
    
    floatEventProperties["chargeMuon3"] = muons->at(2).charge();
    floatEventProperties["idMuon3"] = getId(muons->at(2));
    floatEventProperties["dBMuon3"] = getDeltaB(muons->at(2));
       
	if (muons->at(0).charge()* muons->at(2).charge()== -1.){
		n_OSSF_pairs++;
		muonComb13Vector.SetPxPyPzE(muons->at(0).px()+muons->at(2).px(),muons->at(0).py()+muons->at(2).py(),muons->at(0).pz()+muons->at(2).pz(),muons->at(0).energy()+muons->at(2).energy());
	}
       
	if (muons->at(1).charge()* muons->at(2).charge()== -1.){
		n_OSSF_pairs++;
		muonComb23Vector.SetPxPyPzE(muons->at(1).px()+muons->at(2).px(),muons->at(1).py()+muons->at(2).py(),muons->at(1).pz()+muons->at(2).pz(),muons->at(1).energy()+muons->at(2).energy());
	}
    
	//~ if(muons->at(2).genLepton() != NULL){
      //~ genMuon3Vector.SetPxPyPzE(muons->at(2).genLepton()->px(),muons->at(2).genLepton()->py(),muons->at(2).genLepton()->pz(),muons->at(2).genLepton()->energy());
    //~ }
  }
  if (muons->size() > 3){
    muon4Vector.SetPxPyPzE(muons->at(3).px(),muons->at(3).py(),muons->at(3).pz(),muons->at(3).energy());
    //~ fillLepton< pat::Muon>( "MultiLepton", muons->at(3), "Muon4", floatEventProperties);
    
    floatEventProperties["chargeMuon4"] = muons->at(3).charge();
    floatEventProperties["idMuon4"] = getId(muons->at(3));
    floatEventProperties["dBMuon4"] = getDeltaB(muons->at(3));
           
	if (muons->at(0).charge()* muons->at(3).charge()== -1.){
		n_OSSF_pairs++;
		muonComb14Vector.SetPxPyPzE(muons->at(0).px()+muons->at(3).px(),muons->at(0).py()+muons->at(3).py(),muons->at(0).pz()+muons->at(3).pz(),muons->at(0).energy()+muons->at(3).energy());
	}
       
	if (muons->at(1).charge()* muons->at(3).charge()== -1.){
		n_OSSF_pairs++;
		muonComb24Vector.SetPxPyPzE(muons->at(1).px()+muons->at(3).px(),muons->at(1).py()+muons->at(3).py(),muons->at(1).pz()+muons->at(3).pz(),muons->at(1).energy()+muons->at(3).energy());
	}
       
	if (muons->at(2).charge()* muons->at(3).charge()== -1.){
		n_OSSF_pairs++;
		muonComb34Vector.SetPxPyPzE(muons->at(2).px()+muons->at(3).px(),muons->at(2).py()+muons->at(3).py(),muons->at(2).pz()+muons->at(3).pz(),muons->at(2).energy()+muons->at(3).energy());
	}
    
	//~ if(muons->at(3).genLepton() != NULL){
      //~ genMuon4Vector.SetPxPyPzE(muons->at(3).genLepton()->px(),muons->at(3).genLepton()->py(),muons->at(3).genLepton()->pz(),muons->at(3).genLepton()->energy());
    //~ }
  }

  tLorentzVectorEventProperties["muon1"] = muon1Vector;
  tLorentzVectorEventProperties["muon2"] = muon2Vector;
  tLorentzVectorEventProperties["muon3"] = muon3Vector;
  tLorentzVectorEventProperties["muon4"] = muon4Vector;
  
  //~ tLorentzVectorEventProperties["genMuon1"] = genMuon1Vector;
  //~ tLorentzVectorEventProperties["genMuon2"] = genMuon2Vector;
  //~ tLorentzVectorEventProperties["genMuon3"] = genMuon3Vector;
  //~ tLorentzVectorEventProperties["genMuon4"] = genMuon4Vector;
   
  //~ tLorentzVectorEventProperties["p4Muons12"] = muonComb12Vector;
  //~ tLorentzVectorEventProperties["p4Muons13"] = muonComb13Vector;
  //~ tLorentzVectorEventProperties["p4Muons14"] = muonComb14Vector;
  //~ tLorentzVectorEventProperties["p4Muons23"] = muonComb23Vector;
  //~ tLorentzVectorEventProperties["p4Muons24"] = muonComb24Vector;
  //~ tLorentzVectorEventProperties["p4Muons34"] = muonComb34Vector;

  floatEventProperties["Muons12Mass"] = muonComb12Vector.M();
  floatEventProperties["Muons13Mass"] = muonComb13Vector.M();
  floatEventProperties["Muons14Mass"] = muonComb14Vector.M();
  floatEventProperties["Muons23Mass"] = muonComb23Vector.M();
  floatEventProperties["Muons24Mass"] = muonComb24Vector.M();
  floatEventProperties["Muons34Mass"] = muonComb34Vector.M();

  
  if (n_OSSF_pairs == 0){
	  intEventProperties["nOSSF"] = 0;
  }
  if (n_OSSF_pairs == 1){
	  intEventProperties["nOSSF"] = 1;
  }
  if (n_OSSF_pairs > 1){
	  intEventProperties["nOSSF"] = 2;
  }

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
  
  const std::string treeName = "MultiLepton";
  
  for(std::map<std::string, int>::const_iterator it = intEventProperties.begin(); it != intEventProperties.end(); ++it){
    assert(intBranches_[treeName].find((*it).first) != intBranches_[treeName].end());
    *(intBranches_[treeName][(*it).first]) = (*it).second;
  }
  for(std::map<std::string, float>::const_iterator it = floatEventProperties.begin(); it != floatEventProperties.end(); ++it){
	//~ std::cout << (*it).first<< endl;
    assert(floatBranches_[treeName].find((*it).first) != floatBranches_[treeName].end());
    *(floatBranches_[treeName][(*it).first]) = (*it).second;
  }

  for(std::map<std::string, TLorentzVector>::const_iterator it = tLorentzVectorEventProperties.begin(); it != tLorentzVectorEventProperties.end(); ++it){
    assert(tLorentzVectorBranches_[treeName].find((*it).first) != tLorentzVectorBranches_[treeName].end());
    *(tLorentzVectorBranches_[treeName][(*it).first]) = (*it).second;
  }
  
  trees_[treeName]->Fill();
	  

  delete jecUnc;
  delete shiftedJetsJESUp;
  delete shiftedJetsJESDown;


}

template <class aT, class bT> void 
MultiLeptonTrees::makeCombinations ( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT> &b,const std::vector<reco::PFCandidate>&pfCands, const edm::Event &ev, const pat::MET &patMet, const TLorentzVector &MHT, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties)
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
MultiLeptonTrees::makeCombinations ( const std::string &treeName, const std::vector<aT> &a,const std::vector<reco::PFCandidate>&pfCands, const edm::Event &ev, const pat::MET &patMet , const TLorentzVector &MHT, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties)
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
//~ MultiLeptonTrees::fillLepton( const std::string &treeName, const aT& a, const std::string &leptonFlavorNumber,const  std::map<std::string, float> &floatEventProperties)
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
MultiLeptonTrees::fillTree( const std::string &treeName, const aT& a, const bT& b,const std::vector<reco::PFCandidate>&pfCands, const pat::MET &patMet, const TLorentzVector &MHT)
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
MultiLeptonTrees::getMotherPdgId( const reco::GenParticle &p)
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
MultiLeptonTrees::getLeptonPdgId( const reco::GenParticle &p)
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
MultiLeptonTrees::calcPZeta(const TLorentzVector& p1,const TLorentzVector& p2, const TLorentzVector& met)
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

void MultiLeptonTrees::fillPdfUncert(const edm::Handle< std::vector<double> >& weightHandle, const std::string& pdfIdentifier, const std::string& treeName){
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

const TLorentzVector MultiLeptonTrees::getMomentum(const  pat::Electron &e)
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

const TLorentzVector MultiLeptonTrees::getMomentum(const  pat::Muon &mu)
{
  const TLorentzVector result = TLorentzVector(mu.px(), mu.py(), mu.pz(), mu.energy());
  return result;
}

const TLorentzVector MultiLeptonTrees::getMomentum(const  pat::Tau &tau)
{
  const TLorentzVector result = TLorentzVector(tau.px(), tau.py(), tau.pz(), tau.energy());
  return result;
}

float MultiLeptonTrees::getId(const  pat::Electron &e)
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

float MultiLeptonTrees::getId(const  pat::Muon &mu)
{
  //  std::cout<<"muon " << (lepton.trackIso()+lepton.ecalIso()+lepton.hcalIso())/lepton.pt() << std::endl;
  //  return (mu.isolationR03().hadEt + mu.isolationR03().emEt + mu.isolationR03().sumPt) / mu.pt();
  //  return (mu.chargedHadronIso() + mu.photonIso() + mu.neutralHadronIso()) / mu.pt();
  return fctIsolation_(mu)* 1./mu.pt();
}

float MultiLeptonTrees::getId(const  pat::Tau &tau)
{
  float result = fctIsolation_(tau);
  if(tau.tauID(tauId_) < 0.5)
    result *= -1.0;
  return result;
}


float MultiLeptonTrees::getPhiAtECAL(const  pat::Tau &tau, const std::vector<reco::PFCandidate>&pfCands)
{

  float result = -99.0;
  return result;
}


float MultiLeptonTrees::getPhiAtECAL(const  pat::Muon &mu, const std::vector<reco::PFCandidate>&pfCands)
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

float MultiLeptonTrees::getPhiAtECAL(const  pat::Electron &e, const std::vector<reco::PFCandidate>&pfCands)
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


float MultiLeptonTrees::topPtWeightBen(double topPt){
  if( topPt<0 ) return 1;

  float p0 = 1.18246e+00;
  float p1 = 4.63312e+02;
  float p2 = 2.10061e-06;

  if( topPt>p1 ) topPt = p1;

  float result = p0 + p2 * topPt * ( topPt - 2 * p1 );
  return result;
}

float MultiLeptonTrees::topPtWeightTOP(double topPt){
  if( topPt<0 ) return 1;

  float p0 = 0.156;
  float p1 = 0.00137;

  float result = exp(p0 - p1 * topPt);
  return result;
}


float MultiLeptonTrees::getDeltaB(const  pat::Electron &e)
{
  float result = e.dB(pat::Electron::PV3D);
  return result;
}

float MultiLeptonTrees::getDeltaB(const  pat::Muon &mu)
{
  float result = mu.dB(pat::Muon::PV3D);
  return result;
}

float MultiLeptonTrees::getDeltaB(const  pat::Tau &tau)
{
  float result = -1; // not available for pat::Tau could use ip of leading ch. Hadr if needed.
  return result;
}


float MultiLeptonTrees::transverseMass(const TLorentzVector& p, const TLorentzVector& met)
{
  reco::Candidate::LorentzVector otherMet(met.Px(),met.Py(),met.Pz(),met.E());
  reco::Candidate::LorentzVector leptonT(p.Px(),p.Py(),0.,p.E()*sin(p.Theta()));
  reco::Candidate::LorentzVector sumT=leptonT+otherMet;

  return std::sqrt(sumT.M2());
}

std::string MultiLeptonTrees::convertInputTag(const edm::InputTag tag)
{
  std::string result = tag.label();
  if(tag.instance().length() > 0)
    result = tag.instance();
  //  std::cerr << "'"<<tag.label() << "', '"<< tag.instance()<<"' = '"<< result<<"'"<<std::endl;
  return result;
}

// ------------ Method called once each job just before starting event loop  ------------
void 
MultiLeptonTrees::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MultiLeptonTrees::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(MultiLeptonTrees);
