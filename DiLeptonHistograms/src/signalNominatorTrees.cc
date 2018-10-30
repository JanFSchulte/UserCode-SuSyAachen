// -*- C++ -*-
//
// Package:    Histograms
// Class:      signalNominatorTrees
// 
/**\class signalNominatorTrees signalNominatorTrees.cc brot/signalNominatorTrees/src/signalNominatorTrees.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  matthias edelhoff
//         Created:  Tue Oct 27 13:50:40 CET 2009
// $Id: nominatorTrees.cc,v 1.31 2012/09/17 17:38:58 sprenger Exp $
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

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"


#include <SuSyAachen/DiLeptonHistograms/interface/WeightFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/PdgIdFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/VertexWeightFunctor.h>

//ROOT
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace std;

//
// class decleration
//

class signalNominatorTrees : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit signalNominatorTrees(const edm::ParameterSet&);
  ~signalNominatorTrees();

private:
  typedef pat::Lepton<reco::Candidate> candidate;
  typedef edm::View<candidate> collection;

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  void initFloatBranch( const std::string &name);
  void initIntBranch( const std::string &name);
  

  edm::EDGetTokenT< std::vector< pat::Electron > >      electronToken_;
  edm::EDGetTokenT< std::vector< pat::Muon > >        muonToken_;
  edm::EDGetTokenT< std::vector< pat::Jet > >       jetToken_;  
  edm::EDGetTokenT< std::vector< reco::GenParticle > >    genParticleToken_;
  edm::EDGetTokenT<reco::VertexCollection>          vertexToken_;
  edm::EDGetTokenT<LHEEventProduct>             LHEEventToken_;
  
  //~ 
  edm::Handle< std::vector< pat::Jet > > jets;
  
  PdgIdFunctor getPdgId_;


  //data
  std::map<std::string, TTree*> trees_;  
  std::map<std::string, std::map< std::string, float*> > floatBranches_; 
  std::map<std::string, std::map< std::string, unsigned int*> > intBranches_;
  
  

  bool debug;
};

// constructors and destructor
signalNominatorTrees::signalNominatorTrees(const edm::ParameterSet& iConfig):
  electronToken_    (consumes< std::vector< pat::Electron > >     (iConfig.getParameter<edm::InputTag>("electrons"))),
  muonToken_      (consumes< std::vector< pat::Muon > >     (iConfig.getParameter<edm::InputTag>("muons"))),
  jetToken_       (consumes< std::vector< pat::Jet > >      (iConfig.getParameter<edm::InputTag>("jets"))),
  genParticleToken_   (consumes< std::vector< reco::GenParticle > > (iConfig.getParameter<edm::InputTag>("genParticles"))),
  vertexToken_      (consumes<reco::VertexCollection>       (iConfig.getParameter<edm::InputTag>("vertices"))), 
  LHEEventToken_    (consumes<LHEEventProduct>            (iConfig.getParameter<edm::InputTag>("LHEInfo"))),
  
  getPdgId_( iConfig.getParameter< edm::ParameterSet>("pdgIdDefinition") , consumesCollector() )
{
  usesResource("TFileService");
  debug = false;
  consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  consumes<std::vector< PileupSummaryInfo > >(edm::InputTag("slimmedAddPileupInfo"));

  // init trees
  edm::Service<TFileService> file;
  trees_["Tree"] = file->make<TTree>("Tree", "Tree");
 
  initFloatBranch( "mSbottom" );
  initFloatBranch( "mNeutralino2" );
  initFloatBranch( "ISRCorrection" );
  initFloatBranch( "ISRUncertainty" );
  initIntBranch( "nISRJets" );
  initIntBranch( "nVertices" );
  initIntBranch( "nGenVertices" );
  
  initFloatBranch( "scaleWeight1" );
  initFloatBranch( "scaleWeight2" );
  initFloatBranch( "scaleWeight3" );
  initFloatBranch( "scaleWeight4" );
  initFloatBranch( "scaleWeight5" );
  initFloatBranch( "scaleWeight6" );
  initFloatBranch( "scaleWeight7" );
  initFloatBranch( "scaleWeight8" );
  

}

void
signalNominatorTrees::initFloatBranch(const std::string &name)
{
  for( const auto& it : trees_){
    if(debug) std::cout << it.first <<" - "<< name << std::endl;
    floatBranches_[it.first][name] = new float;
    it.second->Branch(name.c_str(), floatBranches_[it.first][name], (name+"/F").c_str());
  }
}

void 
signalNominatorTrees::initIntBranch(const std::string &name)
{
  for( const auto& it : trees_){
    if(debug) std::cout << it.first <<" - "<< name << std::endl;
    intBranches_[it.first][name] = new unsigned int;
    it.second->Branch(name.c_str(), intBranches_[it.first][name], (name+"/I").c_str());
  }
}

signalNominatorTrees::~signalNominatorTrees()
{ 
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
signalNominatorTrees::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 
  std::map<std::string, int> intEventProperties;
  std::map<std::string, float> floatEventProperties;
  std::map<std::string, TLorentzVector> tLorentzVectorEventProperties;
  
  edm::Handle< std::vector< reco::GenParticle > > genParticles;
  iEvent.getByToken(genParticleToken_, genParticles);
  
  edm::Handle< std::vector< pat::Electron > > electrons;
  iEvent.getByToken(electronToken_, electrons);

  edm::Handle< std::vector< pat::Muon > > muons;
  iEvent.getByToken(muonToken_, muons);
  
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);
  
  iEvent.getByToken(jetToken_, jets);
  
  getPdgId_.loadGenParticles(iEvent);

 
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

  
  floatEventProperties["mSbottom"] = -1;
  floatEventProperties["mNeutralino2"] = -1;
  
  if (genParticles.isValid()){
  
  for (std::vector<reco::GenParticle>::const_iterator itGenParticle = genParticles->begin(); itGenParticle != genParticles->end(); itGenParticle++) {

    
    if ((*itGenParticle).pdgId()== 1000005){
      floatEventProperties["mSbottom"] = (*itGenParticle).mass();
    }
    if ((*itGenParticle).pdgId()== 1000023 ){
      floatEventProperties["mNeutralino2"] = (*itGenParticle).mass();
    }


  }

  }
  
  int nISRJets = 0;
  
  for(std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end() ; ++it){
  if ((*it).pt() >=35.0 && fabs((*it).eta())<2.4){    
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
  }
  
  }
  intEventProperties["nISRJets"] = nISRJets;
  
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

  for(const auto& it : intEventProperties){
    assert(intBranches_["Tree"].find(it.first) != intBranches_["Tree"].end());
    *(intBranches_["Tree"][it.first]) = it.second;
  }
  for(const auto& it : floatEventProperties){
    assert(floatBranches_["Tree"].find(it.first) != floatBranches_["Tree"].end());
    *(floatBranches_["Tree"][it.first]) = it.second;
  }
  
  trees_["Tree"]->Fill();
}



// ------------ Method called once each job just before starting event loop  ------------
void 
signalNominatorTrees::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
signalNominatorTrees::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(signalNominatorTrees);
