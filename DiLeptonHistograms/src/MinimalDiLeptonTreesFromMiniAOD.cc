// -*- C++ -*-
//
// Package:    Histograms
// Class:      MinimalDiLeptonTreesFromMiniAOD
// 
/**\class MinimalDiLeptonTreesFromMiniAOD MinimalDiLeptonTreesFromMiniAOD.cc brot/MinimalDiLeptonTreesFromMiniAOD/src/MinimalDiLeptonTreesFromMiniAOD.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  matthias edelhoff
//         Created:  Tue Oct 27 13:50:40 CET 2009
// $Id: MinimalDiLeptonTreesFromMiniAOD.cc,v 1.31 2012/09/17 17:38:58 sprenger Exp $
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


#include "RecoJets/JetProducers/interface/PileupJetIdAlgo.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"


#include <SuSyAachen/DiLeptonHistograms/interface/WeightFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/PdgIdFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/VertexWeightFunctor.h>
#include <SuSyAachen/TagAndProbeTreeWriter/interface/IsolationFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/TriggerMatchFunctorMiniAOD.h>


//ROOT
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace std;

//
// class decleration
//

class MinimalDiLeptonTreesFromMiniAOD : public edm::EDAnalyzer {
public:
  explicit MinimalDiLeptonTreesFromMiniAOD(const edm::ParameterSet&);
  ~MinimalDiLeptonTreesFromMiniAOD();

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
  template <class aT, class bT> void makeCombinations( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT >&b, const std::vector<pat::PackedCandidate>&pfCands, const edm::Event &ev, const pat::MET &patMet, const TLorentzVector &MHT, const edm::Handle<reco::VertexCollection> &vertices, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties);
  template <class aT> void makeCombinations( const std::string &treeName, const std::vector<aT> &a , const std::vector<pat::PackedCandidate>&pfCands,   const edm::Event &ev, const pat::MET &patMet, const TLorentzVector &MHT, const edm::Handle<reco::VertexCollection> &vertices, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties);
  template<class aT, class bT> void fillTree( const std::string &treeName, const aT &a, const bT &b, const std::vector<pat::PackedCandidate>&pfCands,  const pat::MET &patMet, const TLorentzVector &MHT, const edm::Handle<reco::VertexCollection> &vertices);
  //~ std::pair<double, double> calcPZeta(const TLorentzVector& p1,const TLorentzVector& p2, const TLorentzVector& met);
 
  const TLorentzVector getMomentum(const  pat::Electron &e);
  const TLorentzVector getMomentum(const  pat::Muon &mu);

  float transverseMass(const TLorentzVector& p, const TLorentzVector& met);
  

  edm::EDGetTokenT< std::vector< pat::Electron > > 			electronToken_;
  edm::EDGetTokenT< std::vector< pat::Muon > > 				muonToken_;
  edm::EDGetTokenT< std::vector< pat::Jet > >				jetToken_;
  edm::EDGetTokenT< std::vector< pat::Jet > >		 		bJetToken_;
  edm::EDGetTokenT< std::vector< pat::MET > > 				metToken_;
  edm::EDGetTokenT<reco::VertexCollection> 					vertexToken_;
  edm::EDGetTokenT< std::vector<pat::PackedCandidate>  > 	pfCandToken_;
  edm::EDGetTokenT< std::vector< reco::GenParticle > >	 	genParticleToken_;
  edm::EDGetTokenT<GenEventInfoProduct>	 					genEventInfoToken_;
  edm::EDGetTokenT<LHEEventProduct>	 						LHEEventToken_;
  
  //data
  std::map<std::string, TTree*> trees_;  
  std::map<std::string, std::map< std::string, float*> > floatBranches_; 
  std::map<std::string, std::map< std::string, unsigned int*> > intBranches_; 
  std::map<std::string, std::map< std::string, TLorentzVector*> > tLorentzVectorBranches_;

  edm::Handle< std::vector< pat::Jet > > jets;

  VertexWeightFunctor fctVtxWeight_;
  PdgIdFunctor getPdgId_;
 
  bool debug;
  bool writeID_;
  bool metUncert_;  
  bool triggerMatches_;
};

// constructors and destructor
MinimalDiLeptonTreesFromMiniAOD::MinimalDiLeptonTreesFromMiniAOD(const edm::ParameterSet& iConfig):
  electronToken_	(consumes< std::vector< pat::Electron > > 		(iConfig.getParameter<edm::InputTag>("electrons"))),
  muonToken_		(consumes< std::vector< pat::Muon > >			(iConfig.getParameter<edm::InputTag>("muons"))),
  jetToken_			(consumes< std::vector< pat::Jet > >			(iConfig.getParameter<edm::InputTag>("jets"))),
  bJetToken_		(consumes< std::vector< pat::Jet > >			(iConfig.getParameter<edm::InputTag>("bJets"))),
  metToken_			(consumes< std::vector< pat::MET > >			(iConfig.getParameter<edm::InputTag>("met"))),
  vertexToken_		(consumes<reco::VertexCollection>				(iConfig.getParameter<edm::InputTag>("vertices"))),
  pfCandToken_		(consumes< std::vector<pat::PackedCandidate>  >	(iConfig.getParameter<edm::InputTag>("pfCands"))),
  genParticleToken_	(consumes< std::vector< reco::GenParticle > >	(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genEventInfoToken_(consumes<GenEventInfoProduct>					(iConfig.getParameter<edm::InputTag>("pdfInfo"))),
  LHEEventToken_	(consumes<LHEEventProduct>						(iConfig.getParameter<edm::InputTag>("LHEInfo"))),
  
  //~ susyVars_			(consumes< std:vector<edm::InputTag> >			(iConfig.getParameter<std::vector<edm::InputTag> > ("pdfWeightTags"))),
  
  fctVtxWeight_    (iConfig.getParameter<edm::ParameterSet>("vertexWeights") ,consumesCollector()),
  getPdgId_( iConfig.getParameter< edm::ParameterSet>("pdgIdDefinition"),consumesCollector() )
{
  consumes<std::vector< PileupSummaryInfo > >(edm::InputTag("slimmedAddPileupInfo"));


  // init trees
  edm::Service<TFileService> file;
  trees_["EE"] = file->make<TTree>("EEDileptonTree", "EE DileptonTree");
  trees_["EMu"] = file->make<TTree>("EMuDileptonTree", "EMu DileptonTree");
  trees_["MuMu"] = file->make<TTree>("MuMuDileptonTree", "MuMu DileptonTree");

  initFloatBranch( "genWeight" );     
  initFloatBranch( "weight" );
  initFloatBranch( "chargeProduct" );
  initFloatBranch( "charge1" );
  initFloatBranch( "charge2" );
  initTLorentzVectorBranch( "leptonPair" );
  initTLorentzVectorBranch( "genLeptonPair" );
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
  initTLorentzVectorBranch( "vMet" );  
  initTLorentzVectorBranch( "vMetUncorrected" ); 
  initTLorentzVectorBranch( "vGenMet" ); 
  initFloatBranch( "pt1" );
  initFloatBranch( "pt2" );
  initFloatBranch( "eta1" );
  initFloatBranch( "eta2" );
  initFloatBranch( "mt1" );
  initFloatBranch( "mt2" );
  initFloatBranch( "deltaPhi" );
  initFloatBranch( "deltaR" );
  initFloatBranch( "ht" );
  initFloatBranch( "genHT" );
  initFloatBranch( "mht" );
  initFloatBranch( "met" );
  initFloatBranch( "uncorrectedMet" ); 
  initFloatBranch( "genMet" );
  initIntBranch( "nJets" );
  initIntBranch( "nBJets" );
  initIntBranch( "nVertices" );
  initIntBranch( "nGenVertices" );
  initIntBranch( "nLightLeptons" );
  initIntBranch( "runNr" );
  initIntBranch( "lumiSec" );
  initIntBranch( "eventNr" );
  initIntBranch( "pdgId1" );
  initIntBranch( "pdgId2" );
  initIntBranch( "motherPdgId1" );
  initIntBranch( "motherPdgId2" ); 
      
}

void 
MinimalDiLeptonTreesFromMiniAOD::initTLorentzVectorBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    tLorentzVectorBranches_[(*it).first][name] = new TLorentzVector;
    (*it).second->Branch(name.c_str(), "TLorentzVector" ,&tLorentzVectorBranches_[(*it).first][name]);
  }
}

void 
MinimalDiLeptonTreesFromMiniAOD::initFloatBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    floatBranches_[(*it).first][name] = new float;
    (*it).second->Branch(name.c_str(), floatBranches_[(*it).first][name], (name+"/F").c_str());
  }
}

void 
MinimalDiLeptonTreesFromMiniAOD::initIntBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    intBranches_[(*it).first][name] = new unsigned int;
    (*it).second->Branch(name.c_str(), intBranches_[(*it).first][name], (name+"/I").c_str());
  }
}

MinimalDiLeptonTreesFromMiniAOD::~MinimalDiLeptonTreesFromMiniAOD()
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
MinimalDiLeptonTreesFromMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< std::vector< pat::Electron > > electrons;
  iEvent.getByToken(electronToken_, electrons);

  edm::Handle< std::vector< pat::Muon > > muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle< std::vector<pat::PackedCandidate>  > pfCands;
  iEvent.getByToken(pfCandToken_, pfCands); 

  iEvent.getByToken(jetToken_, jets);

  edm::Handle< std::vector< pat::Jet > > bJets;
  iEvent.getByToken(bJetToken_, bJets);
  
  edm::Handle< std::vector< pat::MET > > mets;
  iEvent.getByToken(metToken_, mets);


  edm::Handle< std::vector< reco::GenParticle > > genParticles;
  iEvent.getByToken(genParticleToken_, genParticles);


  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);

  getPdgId_.loadGenParticles(iEvent);
  std::map<std::string, int> intEventProperties;
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
        floatEventProperties["genHT"] = ht;
    }
	
	else floatEventProperties["genHT"] = -999.;
      	

	
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
    

  intEventProperties["nBJets"] = bJets->size();
  intEventProperties["nLightLeptons"] = electrons->size() + muons->size();
  intEventProperties["runNr"] = iEvent.id().run();
  intEventProperties["lumiSec"] = iEvent.id().luminosityBlock();
  intEventProperties["eventNr"] = iEvent.id().event();



  pat::MET met = mets->front();
  TLorentzVector metVector(mets->front().px(), mets->front().py(), mets->front().pz(), mets->front().energy());
  TLorentzVector uncorrectedMetVector;
  uncorrectedMetVector.SetPtEtaPhiE(mets->front().uncorPt(), 0,	mets->front().uncorPhi(), mets->front().uncorPt());
  
  floatEventProperties["met"] = metVector.Pt();
  tLorentzVectorEventProperties["vMet"] = metVector; 
  
  tLorentzVectorEventProperties["vMetUncorrected"] = uncorrectedMetVector;
  floatEventProperties["uncorrectedMet"] = uncorrectedMetVector.Pt();
  	
  
  TLorentzVector genMetVector(0.,0.,0.,0.);
  
  if (genParticles.isValid()){
	genMetVector.SetPxPyPzE(mets->front().genMET()->px(),mets->front().genMET()->py(),mets->front().genMET()->pz(),mets->front().genMET()->energy());
  }
  
  floatEventProperties["genMet"] = genMetVector.Pt();
  tLorentzVectorEventProperties["vGenMet"] = genMetVector;  
  
  TLorentzVector MHT;  
  TLorentzVector tempMHT;
  
  floatEventProperties["ht"] = 0.0;
  floatEventProperties["mht"] = 0.0;

	
  int nJets=0;

  
  for(std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end() ; ++it){
	if ((*it).pt() >=35.0 && fabs((*it).eta())<2.4){
		nJets++;

		tempMHT.SetPxPyPzE((*it).px(), (*it).py(), (*it).pz(), (*it).energy());	
		MHT = MHT + tempMHT;
		floatEventProperties["ht"] += (*it).pt();
	}		
	
  }
  intEventProperties["nJets"] = nJets;
  

  TLorentzVector jet1Vector(0.,0.,0.,0.);
  TLorentzVector jet2Vector(0.,0.,0.,0.);
  TLorentzVector jet3Vector(0.,0.,0.,0.);
  TLorentzVector jet4Vector(0.,0.,0.,0.); 


  if (nJets > 0)
    jet1Vector.SetPxPyPzE(jets->at(0).px(),jets->at(0).py(),jets->at(0).pz(),jets->at(0).energy());
  if (nJets > 1)
    jet2Vector.SetPxPyPzE(jets->at(1).px(),jets->at(1).py(),jets->at(1).pz(),jets->at(1).energy());
  if (nJets > 2)
    jet3Vector.SetPxPyPzE(jets->at(2).px(),jets->at(2).py(),jets->at(2).pz(),jets->at(2).energy());
  if (nJets > 3)
    jet4Vector.SetPxPyPzE(jets->at(3).px(),jets->at(3).py(),jets->at(3).pz(),jets->at(3).energy());
  tLorentzVectorEventProperties["jet1"] = jet1Vector;
  tLorentzVectorEventProperties["jet2"] = jet2Vector;
  tLorentzVectorEventProperties["jet3"] = jet1Vector;
  tLorentzVectorEventProperties["jet4"] = jet2Vector;


  // bjets

  TLorentzVector bJet1Vector(0.,0.,0.,0.);
  TLorentzVector bJet2Vector(0.,0.,0.,0.);

  if (bJets->size() > 0){
    bJet1Vector.SetPxPyPzE(bJets->at(0).px(),bJets->at(0).py(),bJets->at(0).pz(),bJets->at(0).energy());
  }
  if (bJets->size() > 1){
    bJet2Vector.SetPxPyPzE(bJets->at(1).px(),bJets->at(1).py(),bJets->at(1).pz(),bJets->at(1).energy());
  }
 
  tLorentzVectorEventProperties["bJet1"] = bJet1Vector;
  tLorentzVectorEventProperties["bJet2"] = bJet2Vector;
  
 
  floatEventProperties["weight"] = fctVtxWeight_( iEvent );

  makeCombinations< pat::Electron >("EE", *electrons,*pfCands, iEvent, met, MHT,vertices, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  makeCombinations< pat::Electron, pat::Muon >("EMu", *electrons, *muons,*pfCands, iEvent, met, MHT,vertices, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);
  makeCombinations< pat::Muon >("MuMu", *muons,*pfCands, iEvent, met,MHT,vertices, intEventProperties, floatEventProperties,tLorentzVectorEventProperties);


}

template <class aT, class bT> void 
MinimalDiLeptonTreesFromMiniAOD::makeCombinations ( const std::string &treeName, const std::vector<aT> &a, const std::vector<bT> &b,const std::vector<pat::PackedCandidate>&pfCands, const edm::Event &ev, const pat::MET &patMet, const TLorentzVector &MHT, const edm::Handle<reco::VertexCollection> &vertices, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties)
{
  TLorentzVector met(patMet.px(), patMet.py(), patMet.pz(), patMet.energy());
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
  
  for( typename std::vector<aT>::const_iterator itA = a.begin(); itA != a.end(); ++itA){
    for( typename std::vector<bT>::const_iterator itB = b.begin(); itB != b.end(); ++itB){
      fillTree<aT,bT>( treeName, *itA, *itB,pfCands, patMet,MHT,vertices); 
    }
  }
}

template <class aT> void 
MinimalDiLeptonTreesFromMiniAOD::makeCombinations ( const std::string &treeName, const std::vector<aT> &a,const std::vector<pat::PackedCandidate>&pfCands,const edm::Event &ev, const pat::MET &patMet , const TLorentzVector &MHT, const edm::Handle<reco::VertexCollection> &vertices, const std::map<std::string, int> &intEventProperties, const  std::map<std::string, float> &floatEventProperties, const  std::map<std::string, TLorentzVector> &tLorentzVectorEventProperties)
{
  TLorentzVector met(patMet.px(), patMet.py(), patMet.pz(), patMet.energy());
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

  for( typename std::vector<aT>::const_iterator itA = a.begin(); itA != a.end(); ++itA){
    for( typename std::vector<aT>::const_iterator itB = a.begin(); itB != itA; ++itB){
      fillTree<aT, aT>( treeName, *itA, *itB,pfCands, patMet,MHT,vertices); 
    }
  }
}

template <class aT, class bT> void 
MinimalDiLeptonTreesFromMiniAOD::fillTree( const std::string &treeName, const aT& a, const bT& b,const std::vector<pat::PackedCandidate>&pfCands, const pat::MET &patMet, const TLorentzVector &MHT, const edm::Handle<reco::VertexCollection> &vertices)
{
  if(debug) std::cout << treeName << "- pts:"<< a.pt() << " " << b.pt();
  TLorentzVector aVec( a.px(), a.py(), a.pz(), a.energy() );
  TLorentzVector bVec( b.px(), b.py(), b.pz(), b.energy() );
  TLorentzVector met(patMet.px(), patMet.py(), patMet.pz(), patMet.energy());
  TLorentzVector uncorrectedMet; 
  uncorrectedMet.SetPtEtaPhiE(patMet.uncorPt(), 0, \
			      patMet.uncorPhi(), patMet.uncorPt());  

  TLorentzVector comb = aVec+bVec;
  TLorentzVector MHT2 = comb + MHT;
  
  *(floatBranches_[treeName]["chargeProduct"]) = a.charge()*b.charge();
  *(tLorentzVectorBranches_[treeName]["leptonPair"]) = comb;
  *(tLorentzVectorBranches_[treeName]["lepton1"]) = aVec;
  *(tLorentzVectorBranches_[treeName]["lepton2"]) = bVec;
  *(floatBranches_[treeName]["mht"]) = MHT2.Pt();
  *(floatBranches_[treeName]["pt1"]) = aVec.Pt();
  *(floatBranches_[treeName]["pt2"]) = bVec.Pt();
  *(floatBranches_[treeName]["charge1"]) = a.charge();
  *(floatBranches_[treeName]["charge2"]) = b.charge();
  *(floatBranches_[treeName]["eta1"]) = aVec.Eta();
  *(floatBranches_[treeName]["eta2"]) = bVec.Eta();
  *(floatBranches_[treeName]["mt1"]) = transverseMass(aVec, met);
  *(floatBranches_[treeName]["mt2"]) = transverseMass(bVec, met);
  *(floatBranches_[treeName]["deltaPhi"]) = aVec.DeltaPhi( bVec );
  *(floatBranches_[treeName]["deltaR"]) = aVec.DeltaR( bVec );

  int matched = 0;
  *(intBranches_[treeName]["pdgId1"]) = -9999;
  *(intBranches_[treeName]["pdgId2"]) = -9999;  
  *(intBranches_[treeName]["motherPdgId1"]) =-9999;
  *(intBranches_[treeName]["motherPdgId2"]) =-9999;         
  std::vector<int> pdgIds1 = getPdgId_.operator()<aT>(a); 
  std::vector<int> pdgIds2 = getPdgId_.operator()<bT>(b); 	
  *(intBranches_[treeName]["pdgId1"]) = pdgIds1[0];
  *(intBranches_[treeName]["pdgId2"]) = pdgIds2[0];  
  *(intBranches_[treeName]["motherPdgId1"]) = pdgIds1[1];
  *(intBranches_[treeName]["motherPdgId2"]) = pdgIds2[1];  

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
  *(tLorentzVectorBranches_[treeName]["genLeptonPair"]) = genVec;
  *(tLorentzVectorBranches_[treeName]["genLepton1"]) = genLepton1;
  *(tLorentzVectorBranches_[treeName]["genLepton2"]) = genLepton2;
  if(debug) std::cout << ", matched = "<<matched<<", motherId = "<<pdgIds1[1];
  if(debug) std::cout<<", M = "<< comb.M() <<", chargeProduct = "<< a.charge()*b.charge() <<std::endl;  
  
  trees_[treeName]->Fill();
}


float MinimalDiLeptonTreesFromMiniAOD::transverseMass(const TLorentzVector& p, const TLorentzVector& met)
{
  reco::Candidate::LorentzVector otherMet(met.Px(),met.Py(),met.Pz(),met.E());
  reco::Candidate::LorentzVector leptonT(p.Px(),p.Py(),0.,p.E()*sin(p.Theta()));
  reco::Candidate::LorentzVector sumT=leptonT+otherMet;

  return std::sqrt(sumT.M2());
}


// ------------ Method called once each job just before starting event loop  ------------
void 
MinimalDiLeptonTreesFromMiniAOD::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MinimalDiLeptonTreesFromMiniAOD::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(MinimalDiLeptonTreesFromMiniAOD);
