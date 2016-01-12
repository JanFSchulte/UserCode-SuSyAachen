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

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"


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

class signalNominatorTrees : public edm::EDAnalyzer {
public:
  explicit signalNominatorTrees(const edm::ParameterSet&);
  ~signalNominatorTrees();

private:
  //  typedef reco::Candidate candidate;
  typedef pat::Lepton<reco::Candidate> candidate;
  typedef edm::View<candidate> collection;

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  void initFloatBranch( const std::string &name);
  void initIntBranch( const std::string &name);

  
  edm::InputTag genParticleTag_;
  edm::InputTag LHEEventTag_;  
  
  PdgIdFunctor getPdgId_;


  //data
  std::map<std::string, TTree*> trees_;  
  std::map<std::string, std::map< std::string, float*> > floatBranches_; 
  std::map<std::string, std::map< std::string, unsigned int*> > intBranches_; 




  bool debug;
};

// constructors and destructor
signalNominatorTrees::signalNominatorTrees(const edm::ParameterSet& iConfig): 
  getPdgId_( iConfig.getParameter< edm::ParameterSet>("pdgIdDefinition") )
{
  debug = false;
  
  
  // read config
  genParticleTag_ = iConfig.getParameter<edm::InputTag>("genParticles");
  //~ LHEEventTag_ = iConfig.getParameter<edm::InputTag>("LHEInfo");



  // init trees
  edm::Service<TFileService> file;
  trees_["Tree"] = file->make<TTree>("Tree", "Tree");
 
  initFloatBranch( "mSbottom" );
  initFloatBranch( "mNeutralino2" );

}



void 
signalNominatorTrees::initFloatBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    floatBranches_[(*it).first][name] = new float;
    (*it).second->Branch(name.c_str(), floatBranches_[(*it).first][name], (name+"/F").c_str());
  }
}

void 
signalNominatorTrees::initIntBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    intBranches_[(*it).first][name] = new unsigned int;
    (*it).second->Branch(name.c_str(), intBranches_[(*it).first][name], (name+"/I").c_str());
  }
}

signalNominatorTrees::~signalNominatorTrees()
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
signalNominatorTrees::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 
  std::map<std::string, int> intEventProperties;
  std::map<std::string, float> floatEventProperties;
  std::map<std::string, TLorentzVector> tLorentzVectorEventProperties;
  
  edm::Handle< std::vector< reco::GenParticle > > genParticles;
  iEvent.getByLabel(genParticleTag_, genParticles);
  
  getPdgId_.loadGenParticles(iEvent);

  //~ floatEventProperties["mSbottomLHE"] = -1.;
  //~ floatEventProperties["mNeutralino2LHE"] = -1.;
  //~ 
  //~ edm::Handle<LHEEventProduct> lheInfoHandle;
  //~ iEvent.getByLabel(LHEEventTag_ , lheInfoHandle);
//~ 
  //~ if (lheInfoHandle.isValid()) {
        //~ lhef::HEPEUP lheParticleInfo = lheInfoHandle->hepeup();
        //~ // get the five vector
        //~ // (Px, Py, Pz, E and M in GeV)
        //~ std::vector<lhef::HEPEUP::FiveVector> allParticles = lheParticleInfo.PUP;
        //~ std::vector<int> statusCodes = lheParticleInfo.ISTUP;
//~ 
        //~ for (unsigned int i = 0; i < statusCodes.size(); i++) {
            //~ if (statusCodes[i] == 1) {
                //~ if (abs(lheParticleInfo.IDUP[i]) == 1000005) {
                    //~ floatEventProperties["mSbottomLHE"] = allParticles[i][4];
                //~ }
                //~ if (abs(lheParticleInfo.IDUP[i]) == 1000023) {
                    //~ floatEventProperties["mNeutralino2LHE"] = allParticles[i][4];
                //~ }
            //~ }
        //~ }
  //~ }
  //~ 
  
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

  
  for(std::map<std::string, int>::const_iterator it = intEventProperties.begin(); it != intEventProperties.end(); ++it){
    assert(intBranches_["Tree"].find((*it).first) != intBranches_["Tree"].end());
    *(intBranches_["Tree"][(*it).first]) = (*it).second;
  }
  for(std::map<std::string, float>::const_iterator it = floatEventProperties.begin(); it != floatEventProperties.end(); ++it){
    assert(floatBranches_["Tree"].find((*it).first) != floatBranches_["Tree"].end());
    *(floatBranches_["Tree"][(*it).first]) = (*it).second;
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
