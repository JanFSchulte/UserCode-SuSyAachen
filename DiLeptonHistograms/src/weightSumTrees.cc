// -*- C++ -*-
//
// Package:    Histograms
// Class:      weightSumTrees
// 
/**\class weightSumTrees weightSumTrees.cc brot/weightSumTrees/src/weightSumTrees.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  matthias edelhoff
//         Created:  Tue Oct 27 13:50:40 CET 2009
// $Id: weightSumTrees.cc,v 1.4 2018/10/17 11:56:00 teroerde Exp $
//
//


// system include files
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <DataFormats/Provenance/interface/EventID.h>

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//ROOT
#include "TTree.h"
#include "TFile.h"

using namespace std;

//
// class decleration
//

class weightSumTrees : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit weightSumTrees(const edm::ParameterSet&);
  ~weightSumTrees();

private:  
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  void initFloatBranch( const std::string &name);
 
  edm::EDGetTokenT<GenEventInfoProduct>           genEventInfoToken_;

  //data
  std::map<std::string, TTree*> trees_;  
  std::map<std::string, std::map< std::string, float*> > floatBranches_; 




  bool debug;
};

// constructors and destructor
weightSumTrees::weightSumTrees(const edm::ParameterSet& iConfig):
  genEventInfoToken_    (consumes<GenEventInfoProduct>          (iConfig.getParameter<edm::InputTag>("genInfo")))
{
  usesResource("TFileService");
  debug = false;

  // init trees
  edm::Service<TFileService> file;
  trees_["Tree"] = file->make<TTree>("Tree", "Tree");
 
  initFloatBranch( "genWeight" );  
  initFloatBranch( "genWeightAbsValue" );  

}



void 
weightSumTrees::initFloatBranch(const std::string &name)
{
  for( const auto& it : trees_){
    if(debug) std::cout << it.first <<" - "<< name << std::endl;
    floatBranches_[it.first][name] = new float;
    it.second->Branch(name.c_str(), floatBranches_[it.first][name], (name+"/F").c_str());
  }
}
 
weightSumTrees::~weightSumTrees()
{ 
  for( const auto& it: floatBranches_){
    for( const auto& it2 : it.second){
      if(debug)std::cout << "deleting: " << it.first << " - "<< it2.first << std::endl;
      delete it2.second;
    }
  }


}

// member functions
// ------------ method called to for each event  ------------
void
weightSumTrees::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 
  std::map<std::string, int> intEventProperties;
  std::map<std::string, float> floatEventProperties;

  edm::Handle<GenEventInfoProduct> genInfoProduct;
  iEvent.getByToken(genEventInfoToken_, genInfoProduct);  
  if (genInfoProduct.isValid()){
    floatEventProperties["genWeightAbsValue"] = (*genInfoProduct).weight();
    if ((*genInfoProduct).weight() < 0.0){
    
      floatEventProperties["genWeight"] = -1;
    }
    else{
      floatEventProperties["genWeight"] = 1;    
    }
  }
  else{
 
  floatEventProperties["genWeight"] = 1; 
  floatEventProperties["genWeightAbsValue"] = 1;
  
  }  
  
  for(const auto& it : floatEventProperties){
    assert(floatBranches_["Tree"].find(it.first) != floatBranches_["Tree"].end());
    *(floatBranches_["Tree"][it.first]) = it.second;
  }
  trees_["Tree"]->Fill();
}





// ------------ Method called once each job just before starting event loop  ------------
void 
weightSumTrees::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
weightSumTrees::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(weightSumTrees);
