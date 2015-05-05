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
// $Id: DiLeptonTrees.cc,v 1.31 2012/09/17 17:38:58 sprenger Exp $
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
  template <class aT, class bT> void makeCombinations( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT >&b, const std::vector<reco::PFCandidate>&pfCands, const std::vector<reco::GenParticle>&genParticles, const edm::Event &ev, const pat::MET &patMet, const TLorentzVector &MHT, const edm::Handle<reco::VertexCollection> &vertices, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties);
  template <class aT> void makeCombinations( const std::string &treeName, const std::vector<aT> &a , const std::vector<reco::PFCandidate>&pfCands, const std::vector<reco::GenParticle>&genParticles,  const edm::Event &ev, const pat::MET &patMet, const TLorentzVector &MHT, const edm::Handle<reco::VertexCollection> &vertices, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties);
  template<class aT, class bT> void fillTree( const std::string &treeName, const aT &a, const bT &b, const std::vector<reco::PFCandidate>&pfCands, const std::vector<reco::GenParticle>&genParticles, const pat::MET &patMet, const TLorentzVector &MHT, const edm::Handle<reco::VertexCollection> &vertices);
  int getLeptonPdgId( const reco::GenParticle &p);
  int getMotherPdgId( const reco::GenParticle &p);
  //~ std::pair<double, double> calcPZeta(const TLorentzVector& p1,const TLorentzVector& p2, const TLorentzVector& met);
  void fillPdfUncert(const edm::Handle< std::vector<double> >& weightHandle, const std::string& pdfIdentifier, const std::string& treeName);

  const TLorentzVector getMomentum(const  pat::Electron &e);
  const TLorentzVector getMomentum(const  pat::Muon &mu);
  const TLorentzVector getMomentum(const  pat::Tau &tau);
  void fillLeptonIDs(const std::string &treeName, const  pat::Electron &ele1, const  pat::Electron &ele2, const edm::Handle<reco::VertexCollection> &vertices);
  void fillLeptonIDs(const std::string &treeName, const  pat::Muon &mu1, const  pat::Muon &mu2, const edm::Handle<reco::VertexCollection> &vertices);
  void fillLeptonIDs(const std::string &treeName, const  pat::Electron &ele1, const  pat::Muon &mu2, const edm::Handle<reco::VertexCollection> &vertices);
  void fillLeptonIDs(const std::string &treeName, const  pat::Electron &ele1, const  pat::Tau &tau2, const edm::Handle<reco::VertexCollection> &vertices);
  void fillLeptonIDs(const std::string &treeName, const  pat::Muon &mu1, const  pat::Tau &tau2, const edm::Handle<reco::VertexCollection> &vertices);
  void fillLeptonIDs(const std::string &treeName, const  pat::Tau &tau1, const  pat::Tau &tau2, const edm::Handle<reco::VertexCollection> &vertices);
  template <class aT> std::vector<int>  getGenMatches(const aT &l, const std::vector<reco::GenParticle>&genParticles);
  float topPtWeightBen(double topPt);
  float topPtWeightTOP(double topPt);
  float getIso(const  pat::Electron &e, const std::string &method);
  float getIso(const  pat::Muon &mu, const std::string &method);
  float getIso(const  pat::Tau &tau, const std::string &method);
  float getDeltaB(const  pat::Electron &e);
  float getDeltaB(const  pat::Muon &mu);
  float getDeltaB(const  pat::Tau &tau);
  float transverseMass(const TLorentzVector& p, const TLorentzVector& met);
  float getAEffEle(double eta);
  float getAEffMu(double eta);
  std::string convertInputTag(const edm::InputTag tag);

  edm::InputTag eTag_;
  edm::InputTag muTag_;
  edm::InputTag tauTag_;
  edm::InputTag jetTag_;
  edm::InputTag jet2Tag_;
  edm::InputTag bJetTag_;
  edm::InputTag metTag_;
  //~ edm::InputTag type1metTag_;
  edm::InputTag vertexTag_;
  edm::InputTag vertexTagUp_;
  edm::InputTag vertexTagDown_;
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
  VertexWeightFunctor fctVtxWeightUp_;
  VertexWeightFunctor fctVtxWeightDown_;
  IsolationFunctor fctIsolation_;
  PdgIdFunctor getPdgId_;

  bool debug;
  bool useJets2_;
  bool useTaus_;
};

// constructors and destructor
DiLeptonTrees::DiLeptonTrees(const edm::ParameterSet& iConfig):
  fctVtxWeight_    (iConfig.getParameter<edm::ParameterSet>("vertexWeights") ),
  fctVtxWeightUp_    (iConfig.getParameter<edm::ParameterSet>("vertexWeightsUp") ),
  fctVtxWeightDown_    (iConfig.getParameter<edm::ParameterSet>("vertexWeightsDown") ),
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
  //~ type1metTag_ = iConfig.getParameter<edm::InputTag>("type1met");
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
  trees_["EE"] = file->make<TTree>("EEDileptonTree", "EE DileponTree");
  trees_["EMu"] = file->make<TTree>("EMuDileptonTree", "EMu DileponTree");
  trees_["MuMu"] = file->make<TTree>("MuMuDileptonTree", "MuMu DileponTree");
  if(useTaus_){
    trees_["ETau"] = file->make<TTree>("ETauDileptonTree", "ETau DileponTree");
    trees_["MuTau"] = file->make<TTree>("MuTauDileptonTree", "MuTau DileponTree");
    trees_["TauTau"] = file->make<TTree>("TauTauDileptonTree", "TauTau DileponTree");
  }
  initFloatBranch( "weight" );
  initFloatBranch( "weightUp" );
  initFloatBranch( "weightDown" );
  initFloatBranch( "genPtTop1" );
  initFloatBranch( "genPtTop2" );
  initFloatBranch( "chargeProduct" );
  initFloatBranch( "charge1" );
  initFloatBranch( "charge2" );
  initTLorentzVectorBranch( "p4" );
  initTLorentzVectorBranch( "vMet" );
  initTLorentzVectorBranch( "p4Gen" );
  initTLorentzVectorBranch( "lepton1" );
  initTLorentzVectorBranch( "lepton2" );
  initTLorentzVectorBranch( "genLepton1" );
  initTLorentzVectorBranch( "genLepton2" );
  initTLorentzVectorBranch( "jet1" );
  initTLorentzVectorBranch( "jet2" );
  initTLorentzVectorBranch( "jet3" );
  initTLorentzVectorBranch( "jet4" );
  initTLorentzVectorBranch( "bJet1" );
  initTLorentzVectorBranch( "bJet2" );
  initTLorentzVectorBranch( "vMetUncorrected" );
  initFloatBranch( "rho" );
  initFloatBranch( "pt1" );
  initFloatBranch( "pt2" );
  initFloatBranch( "eta1" );
  initFloatBranch( "eta2" );
  initFloatBranch( "isoEffArea1" );
  initFloatBranch( "isoEffArea2" );
  initFloatBranch( "isoPF1" );
  initFloatBranch( "isoPF2" );
    
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
  
  initFloatBranch( "dB1" );
  initFloatBranch( "dB2" );
  initFloatBranch( "mt1" );
  initFloatBranch( "mt2" );
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
  initFloatBranch( "uncorrectedMet" );
  initFloatBranch( "metJESUp" );
  initFloatBranch( "metJESDown" );
  //~ initFloatBranch( "pZeta" );
  //~ initFloatBranch( "pZetaVis" );
  initIntBranch( "nJets" );
  initIntBranch( "nJetsOld" );
  initIntBranch( "nBJets" );
  initIntBranch( "nShiftedJetsJESUp" );
  initIntBranch( "nShiftedJetsJESDown" );
  initIntBranch( "nVertices" );
  initIntBranch( "nLightLeptons" );
  initFloatBranch( "jet1pt" );
  initFloatBranch( "jet2pt" );
  initFloatBranch( "jet3pt" );
  initFloatBranch( "jet4pt" );
  initFloatBranch( "bjet1pt" );
  initFloatBranch( "bjet2pt" );
  initFloatBranch( "sqrts" );
  initIntBranch( "runNr" );
  initIntBranch( "lumiSec" );
  initIntBranch( "eventNr" );
  initIntBranch( "pdgId1" );
  initIntBranch( "pdgId2" );
  initIntBranch( "matched" );
  initIntBranch( "motherPdgId1" );
  initIntBranch( "motherPdgId2" );
  initIntBranch( "pdgIdETH1" );
  initIntBranch( "pdgIdETH2" );
  initIntBranch( "motherPdgIdETH1" );
  initIntBranch( "motherPdgIdETH2" ); 
  initIntBranch( "grandMotherPdgIdETH1" );
  initIntBranch( "grandMotherPdgIdETH2" );
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
    intBranches_[(*it).first][name] = new unsigned int;
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
DiLeptonTrees::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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


  edm::Handle< std::vector< reco::GenParticle > > genParticles;
  iEvent.getByLabel(genParticleTag_, genParticles);

  //~ edm::Handle< std::vector< reco::PFMET > > type1mets;
  //~ iEvent.getByLabel(type1metTag_, type1mets);

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
  intEventProperties["runNr"] = iEvent.id().run();
  intEventProperties["lumiSec"] = iEvent.id().luminosityBlock();
  intEventProperties["eventNr"] = iEvent.id().event();

  pat::MET met = mets->front();
  TLorentzVector metVector(mets->front().px(), mets->front().py(), mets->front().pz(), mets->front().energy());
  TLorentzVector uncorrectedMetVector;
  uncorrectedMetVector.SetPtEtaPhiE(mets->front().uncorrectedPt(), 0,	mets->front().uncorrectedPhi(), mets->front().uncorrectedPt());
  
  floatEventProperties["uncorrectedMet"] =uncorrectedMetVector.Pt();
  tLorentzVectorEventProperties["vMetUncorrected"] = uncorrectedMetVector;

  //~ reco::MET type1met = type1mets->front();
  //~ TLorentzVector type1metVector(type1mets->front().px(), type1mets->front().py(), type1mets->front().pz(), type1mets->front().energy());

  floatEventProperties["met"] = metVector.Pt();
  tLorentzVectorEventProperties["vMet"] = metVector;

  floatEventProperties["genPtTop1"] = -1;
  floatEventProperties["genPtTop2"] = -1;
  if (genParticles.isValid()){
	
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
	if ((*it).pt() >=35.0 and abs((*it).pt())<2.4){
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

  int nJetsOld=0;
  for(std::vector<pat::Jet>::const_iterator it = jets_->begin(); it != jets_->end() ; ++it){
	if ((*it).pt() >=40.0){
		nJetsOld++;
	}
  }
  intEventProperties["nJetsOld"] = nJetsOld;


  TLorentzVector jet1Vector(0.,0.,0.,0.);
  TLorentzVector jet2Vector(0.,0.,0.,0.);
  TLorentzVector jet3Vector(0.,0.,0.,0.);
  TLorentzVector jet4Vector(0.,0.,0.,0.); 


  if (nJets > 0)
    jet1Vector.SetPxPyPzE(jets_->at(0).px(),jets_->at(0).py(),jets_->at(0).pz(),jets_->at(0).energy());
  if (nJets > 1)
    jet2Vector.SetPxPyPzE(jets_->at(1).px(),jets_->at(1).py(),jets_->at(1).pz(),jets_->at(1).energy());
  if (nJets > 2)
    jet3Vector.SetPxPyPzE(jets_->at(2).px(),jets_->at(2).py(),jets_->at(2).pz(),jets_->at(2).energy());
  if (nJets > 3)
    jet4Vector.SetPxPyPzE(jets_->at(3).px(),jets_->at(3).py(),jets_->at(3).pz(),jets_->at(3).energy());
  tLorentzVectorEventProperties["jet1"] = jet1Vector;
  tLorentzVectorEventProperties["jet2"] = jet2Vector;
  tLorentzVectorEventProperties["jet3"] = jet3Vector;
  tLorentzVectorEventProperties["jet4"] = jet4Vector;




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

    	
  floatEventProperties["weight"] = fctVtxWeight_( iEvent );
  floatEventProperties["weightUp"] = fctVtxWeightUp_( iEvent );
  floatEventProperties["weightDown"] = fctVtxWeightDown_( iEvent );
  
  if(useJets2_) {
    edm::Handle< std::vector< pat::Jet > > jets2;
    iEvent.getByLabel(jet2Tag_, jets2);
    intEventProperties["nJets2"] = jets2->size();
    floatEventProperties["ht2"] = 0.0;
    for(std::vector<pat::Jet>::const_iterator it = jets2->begin(); it != jets2->end() ; ++it){
      floatEventProperties["ht2"] += (*it).pt();
    }
  }

  makeCombinations< pat::Electron >("EE", *electrons,*pfCands,*genParticles, iEvent, met,MHT,vertices, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  makeCombinations< pat::Electron, pat::Muon >("EMu", *electrons, *muons,*pfCands,*genParticles, iEvent, met,MHT,vertices, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  makeCombinations< pat::Muon >("MuMu", *muons,*pfCands,*genParticles, iEvent, met,MHT,vertices, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  if(useTaus_){
    makeCombinations< pat::Electron, pat::Tau >("ETau", *electrons, *taus,*pfCands,*genParticles, iEvent, met,MHT,vertices, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
    makeCombinations< pat::Muon, pat::Tau>("MuTau", *muons, *taus,*pfCands,*genParticles, iEvent, met,MHT,vertices, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
    makeCombinations< pat::Tau >("TauTau", *taus,*pfCands,*genParticles, iEvent, met,MHT,vertices, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  }
  //  if( nMu != 2) std::cout << "-------! "<<nMu<<std::endl;

  delete jecUnc;
  delete shiftedJetsJESUp;
  delete shiftedJetsJESDown;


}

template <class aT, class bT> void 
DiLeptonTrees::makeCombinations ( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT> &b,const std::vector<reco::PFCandidate>&pfCands, const std::vector<reco::GenParticle>&genParticles, const edm::Event &ev, const pat::MET &patMet, const TLorentzVector &MHT, const edm::Handle<reco::VertexCollection> &vertices, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties)
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
      fillTree<aT,bT>( treeName, *itA, *itB,pfCands,genParticles, patMet,MHT,vertices); 
    }
  }
}

template <class aT> void 
DiLeptonTrees::makeCombinations ( const std::string &treeName, const std::vector<aT> &a,const std::vector<reco::PFCandidate>&pfCands, const std::vector<reco::GenParticle>&genParticles, const edm::Event &ev, const pat::MET &patMet , const TLorentzVector &MHT, const edm::Handle<reco::VertexCollection> &vertices, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties)
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
      fillTree<aT, aT>( treeName, *itA, *itB,pfCands,genParticles, patMet,MHT,vertices); 
    }
  }
}

template <class aT, class bT> void 
DiLeptonTrees::fillTree( const std::string &treeName, const aT& a, const bT& b,const std::vector<reco::PFCandidate>&pfCands, const std::vector<reco::GenParticle>&genParticles, const pat::MET &patMet, const TLorentzVector &MHT, const edm::Handle<reco::VertexCollection> &vertices)
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
  //~ std::pair<double, double> pZeta = calcPZeta(aVec, bVec, met);
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
  *(floatBranches_[treeName]["isoEffArea1"]) = getIso(a,"effectiveArea");
  *(floatBranches_[treeName]["isoEffArea2"]) = getIso(b,"effectiveArea");
  *(floatBranches_[treeName]["isoPF1"]) = getIso(a,"deltaBeta");
  *(floatBranches_[treeName]["isoPF2"]) = getIso(b,"deltaBeta");
  *(floatBranches_[treeName]["dB1"]) = getDeltaB(a);
  *(floatBranches_[treeName]["dB2"]) = getDeltaB(b);
  *(floatBranches_[treeName]["mt1"]) = transverseMass(aVec, met);
  *(floatBranches_[treeName]["mt2"]) = transverseMass(bVec, met);
  *(floatBranches_[treeName]["fakeWeight1"]) = fakeRates_(a);
  *(floatBranches_[treeName]["fakeWeight2"]) = fakeRates_(b);
  *(floatBranches_[treeName]["deltaPhi"]) = aVec.DeltaPhi( bVec );
  *(floatBranches_[treeName]["deltaR"]) = aVec.DeltaR( bVec );
  *(floatBranches_[treeName]["jzb"]) = (uncorrectedMet+comb).Pt() - comb.Pt();
  //~ *(floatBranches_[treeName]["pZeta"]) = pZeta.first;
  //~ *(floatBranches_[treeName]["pZetaVis"]) = pZeta.second;

  if (jets_->size() >= 2){
    TLorentzVector vJet1 = TLorentzVector(jets_->at(0).p4().x(), jets_->at(0).p4().y(), jets_->at(0).p4().z(), jets_->at(0).p4().t());
    TLorentzVector vJet2 = TLorentzVector(jets_->at(1).p4().x(), jets_->at(1).p4().y(), jets_->at(1).p4().z(), jets_->at(1).p4().t());
    TLorentzVector sub = aVec + bVec + vJet1 + vJet2;
    *(floatBranches_[treeName]["sqrts"]) = std::sqrt(TMath::Power((std::sqrt(sub.M2() + sub.Perp2()) + met.Et()), 2.0) - (sub.Vect().XYvector() + met.Vect().XYvector())*(sub.Vect().XYvector() + met.Vect().XYvector()));
  } else {
    *(floatBranches_[treeName]["sqrts"]) = -1.0;
  }
  
  fillLeptonIDs(treeName, a , b , vertices );

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
  
  //ETH style genMatching
  *(intBranches_[treeName]["pdgIdETH1"]) = -9999;
  *(intBranches_[treeName]["pdgIdETH2"]) = -9999;  
  *(intBranches_[treeName]["motherPdgIdETH1"]) =-9999;
  *(intBranches_[treeName]["motherPdgIdETH2"]) =-9999;  
  *(intBranches_[treeName]["grandMotherPdgIdETH1"]) = -9999;
  *(intBranches_[treeName]["grandMotherPdgIdETH2"]) = -9999;   
  if (genParticles.size() > 0){

	std::vector<int> pdgIds1 = getGenMatches(a,genParticles); 
	std::vector<int> pdgIds2 = getGenMatches(b,genParticles); 	
	*(intBranches_[treeName]["pdgIdETH1"]) = pdgIds1[0];
	*(intBranches_[treeName]["pdgIdETH2"]) = pdgIds2[0];  
	*(intBranches_[treeName]["motherPdgIdETH1"]) = pdgIds1[1];
	*(intBranches_[treeName]["motherPdgIdETH2"]) = pdgIds2[1];  
	*(intBranches_[treeName]["grandMotherPdgIdETH1"]) = pdgIds1[2];
	*(intBranches_[treeName]["grandMotherPdgIdETH2"]) = pdgIds2[2]; 
  }
  
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

int 
DiLeptonTrees::getLeptonPdgId( const reco::GenParticle &p)
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

template <class aT> std::vector<int>  
DiLeptonTrees::getGenMatches(const aT  &l, const std::vector<reco::GenParticle>&genParticles)
{
	double minDR = 999.;
	double deltaR = 0.;
	const reco::GenParticle * matchedGenPart = new reco::GenParticle();
	const reco::GenParticle * matchedGenPartMother = new reco::GenParticle();
	const reco::GenParticle * matchedGenPartGrandMother = new reco::GenParticle();	
	bool matched = false;

	int pdgId = -9999;
	int motherPdgId = -9999;
	int grandMotherPdgId = -9999;
	
	 std::vector<int> res;
	
	for (std::vector<reco::GenParticle>::const_iterator itGenParticle = genParticles.begin(); itGenParticle != genParticles.end(); itGenParticle++) {
	    
	    if (itGenParticle->status() != 1) continue;
	  
	    deltaR = reco::deltaR(l.eta(),l.phi(), itGenParticle->eta(), itGenParticle->phi());
	    if (deltaR > 0.1) continue;
	    
	    double ndpt = fabs(l.pt() - itGenParticle->pt())/itGenParticle->pt();
	    if(ndpt > 2.) continue;
	    
	    if(deltaR > minDR) continue;
	    
	    minDR = deltaR;
	    
	    matched = true;
	    matchedGenPart = &(*itGenParticle);
	    

	}
	if (matched == true){
	    
	  pdgId = matchedGenPart->pdgId();
	  
	  if (!matchedGenPart->mother()){
	      motherPdgId = -9999;
	      grandMotherPdgId = -9999;
	  }
	  else{
	      matchedGenPartMother = static_cast<const reco::GenParticle*> (matchedGenPart->mother());
	      
	      if (matchedGenPartMother->pdgId() != pdgId){
		motherPdgId = matchedGenPartMother->pdgId();
	      }
	      else{
		

		int tempid = matchedGenPartMother->pdgId();
		int loop_counter = 0;
		while(tempid == pdgId){
		  loop_counter++;
		  if(loop_counter>=10){
		    matchedGenPartMother = static_cast<const reco::GenParticle*>(matchedGenPart->mother());
		    break;
		  }
		  matchedGenPart = static_cast<const reco::GenParticle*>(matchedGenPartMother->mother());
		  matchedGenPartMother = static_cast<const reco::GenParticle*>(matchedGenPartMother->mother());		  
		  tempid = matchedGenPartMother->pdgId();
		}
		
	      }
	      motherPdgId = matchedGenPartMother->pdgId();
		
	  }
	  if (!matchedGenPartMother->mother()){
	      grandMotherPdgId = -9999;
	  }	  
	  else{
	      matchedGenPartGrandMother = static_cast<const reco::GenParticle*> (matchedGenPartMother->mother());
	      
	      if (matchedGenPartGrandMother->pdgId() != motherPdgId){
		grandMotherPdgId = matchedGenPartGrandMother->pdgId();
	      }
	      else{
		
		
		int tempid = matchedGenPartGrandMother->pdgId();
		int loop_counter = 0;
		while(tempid == motherPdgId){
		  std::cout << tempid << std::cout << endl;
		  loop_counter++;
		  if(loop_counter>=10){
		    matchedGenPartGrandMother = static_cast<const reco::GenParticle*>(matchedGenPartMother->mother());
		    break;
		  }
		  matchedGenPartMother = static_cast<const reco::GenParticle*>(matchedGenPartMother->mother());
		  matchedGenPartGrandMother = static_cast<const reco::GenParticle*>(matchedGenPartMother->mother());		  
		  tempid = matchedGenPartGrandMother->pdgId();
		}
		
	      }
	      grandMotherPdgId = matchedGenPartGrandMother->pdgId();
		
	  }	      
	}
	res.push_back(pdgId);
	res.push_back(motherPdgId);
	res.push_back(grandMotherPdgId);

	return res;

}

//from DQM/Physics/src/EwkTauDQM.cc
//~ std::pair<double, double>
//~ DiLeptonTrees::calcPZeta(const TLorentzVector& p1,const TLorentzVector& p2, const TLorentzVector& met)
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

void DiLeptonTrees::fillPdfUncert(const edm::Handle< std::vector<double> >& weightHandle, const std::string& pdfIdentifier, const std::string& treeName){
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

const TLorentzVector DiLeptonTrees::getMomentum(const  pat::Electron &e)
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

const TLorentzVector DiLeptonTrees::getMomentum(const  pat::Muon &mu)
{
  const TLorentzVector result = TLorentzVector(mu.px(), mu.py(), mu.pz(), mu.energy());
  return result;
}

const TLorentzVector DiLeptonTrees::getMomentum(const  pat::Tau &tau)
{
  const TLorentzVector result = TLorentzVector(tau.px(), tau.py(), tau.pz(), tau.energy());
  return result;
}

float DiLeptonTrees::getIso(const  pat::Electron &e, const std::string &method)
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

float DiLeptonTrees::getIso(const  pat::Muon &mu, const std::string &method)
{
  //  std::cout<<"muon " << (lepton.trackIso()+lepton.ecalIso()+lepton.hcalIso())/lepton.pt() << std::endl;
  //  return (mu.isolationR03().hadEt + mu.isolationR03().emEt + mu.isolationR03().sumPt) / mu.pt();
  //  return (mu.chargedHadronIso() + mu.photonIso() + mu.neutralHadronIso()) / mu.pt();
  return fctIsolation_(mu,method)* 1./mu.pt();
}

float DiLeptonTrees::getIso(const  pat::Tau &tau, const std::string &method)
{
  float result = fctIsolation_(tau,method);
  if(tau.tauID(tauId_) < 0.5)
    result *= -1.0;
  return result;
}

void DiLeptonTrees::fillLeptonIDs(const std::string &treeName, const  pat::Electron &ele1, const  pat::Electron &ele2, const edm::Handle<reco::VertexCollection> &vertices)
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

}


void DiLeptonTrees::fillLeptonIDs(const std::string &treeName, const  pat::Electron &ele1, const  pat::Tau &tau2, const edm::Handle<reco::VertexCollection> &vertices)
{

  *(floatBranches_[treeName]["effectiveArea1"]) = getAEffEle(ele1.eta());
  *(floatBranches_[treeName]["chargedIso1"]) = ele1.pfIsolationVariables().sumChargedHadronPt;
  *(floatBranches_[treeName]["neutralIso1"]) = ele1.pfIsolationVariables().sumNeutralHadronEt;
  *(floatBranches_[treeName]["photonIso1"]) = ele1.pfIsolationVariables().sumPhotonEt;
  *(floatBranches_[treeName]["puIso1"]) = ele1.pfIsolationVariables().sumPUPt;

  *(floatBranches_[treeName]["effectiveArea2"]) = -999.;
  *(floatBranches_[treeName]["chargedIso2"]) = -999.;
  *(floatBranches_[treeName]["neutralIso2"]) = -999.;
  *(floatBranches_[treeName]["photonIso2"]) = -999.;
  *(floatBranches_[treeName]["puIso2"]) = -999.;

}


void DiLeptonTrees::fillLeptonIDs(const std::string &treeName, const  pat::Tau &tau1, const  pat::Tau &tau2, const edm::Handle<reco::VertexCollection> &vertices)
{

  *(floatBranches_[treeName]["effectiveArea1"]) = -999.;
  *(floatBranches_[treeName]["chargedIso1"]) = -999.;
  *(floatBranches_[treeName]["neutralIso1"]) = -999.;
  *(floatBranches_[treeName]["photonIso1"]) = -999.;
  *(floatBranches_[treeName]["puIso1"]) = -999.;

  *(floatBranches_[treeName]["effectiveArea2"]) = -999.;
  *(floatBranches_[treeName]["chargedIso2"]) = -999.;
  *(floatBranches_[treeName]["neutralIso2"]) = -999.;
  *(floatBranches_[treeName]["photonIso2"]) = -999.;
  *(floatBranches_[treeName]["puIso2"]) = -999.;

}



void DiLeptonTrees::fillLeptonIDs(const std::string &treeName, const  pat::Electron &ele1, const  pat::Muon &mu2, const edm::Handle<reco::VertexCollection> &vertices)
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

}

void DiLeptonTrees::fillLeptonIDs(const std::string &treeName, const  pat::Muon &mu1, const  pat::Muon &mu2, const edm::Handle<reco::VertexCollection> &vertices)
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


}


void DiLeptonTrees::fillLeptonIDs(const std::string &treeName, const  pat::Muon &mu1, const  pat::Tau &tau2, const edm::Handle<reco::VertexCollection> &vertices)
{

  *(floatBranches_[treeName]["effectiveArea1"]) = getAEffMu(mu1.eta());
  *(floatBranches_[treeName]["chargedIso1"]) = mu1.pfIsolationR03().sumChargedHadronPt;
  *(floatBranches_[treeName]["neutralIso1"]) = mu1.pfIsolationR03().sumNeutralHadronEt;
  *(floatBranches_[treeName]["photonIso1"]) = mu1.pfIsolationR03().sumPhotonEt;
  *(floatBranches_[treeName]["puIso1"]) = mu1.pfIsolationR03().sumPUPt;

  *(floatBranches_[treeName]["effectiveArea2"]) = -999.;
  *(floatBranches_[treeName]["chargedIso2"]) = -999.;
  *(floatBranches_[treeName]["neutralIso2"]) = -999.;
  *(floatBranches_[treeName]["photonIso2"]) = -999.;
  *(floatBranches_[treeName]["puIso2"]) = -999.;

}

float DiLeptonTrees::topPtWeightBen(double topPt){
  if( topPt<0 ) return 1;

  float p0 = 1.18246e+00;
  float p1 = 4.63312e+02;
  float p2 = 2.10061e-06;

  if( topPt>p1 ) topPt = p1;

  float result = p0 + p2 * topPt * ( topPt - 2 * p1 );
  return result;
}

float DiLeptonTrees::topPtWeightTOP(double topPt){
  if( topPt<0 ) return 1;

  float p0 = 0.156;
  float p1 = 0.00137;

  float result = exp(p0 - p1 * topPt);
  return result;
}


float DiLeptonTrees::getDeltaB(const  pat::Electron &e)
{
  float result = e.dB(pat::Electron::PV3D);
  return result;
}

float DiLeptonTrees::getDeltaB(const  pat::Muon &mu)
{
  float result = mu.dB(pat::Muon::PV3D);
  return result;
}

float DiLeptonTrees::getDeltaB(const  pat::Tau &tau)
{
  float result = -1; // not available for pat::Tau could use ip of leading ch. Hadr if needed.
  return result;
}

float DiLeptonTrees::getAEffEle(double eta)
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

float DiLeptonTrees::getAEffMu(double eta)
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


float DiLeptonTrees::transverseMass(const TLorentzVector& p, const TLorentzVector& met)
{
  reco::Candidate::LorentzVector otherMet(met.Px(),met.Py(),met.Pz(),met.E());
  reco::Candidate::LorentzVector leptonT(p.Px(),p.Py(),0.,p.E()*sin(p.Theta()));
  reco::Candidate::LorentzVector sumT=leptonT+otherMet;

  return std::sqrt(sumT.M2());
}

std::string DiLeptonTrees::convertInputTag(const edm::InputTag tag)
{
  std::string result = tag.label();
  if(tag.instance().length() > 0)
    result = tag.instance();
  //  std::cerr << "'"<<tag.label() << "', '"<< tag.instance()<<"' = '"<< result<<"'"<<std::endl;
  return result;
}

// ------------ Method called once each job just before starting event loop  ------------
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
