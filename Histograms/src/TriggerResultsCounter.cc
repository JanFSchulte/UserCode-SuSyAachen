// -*- C++ -*-
//
// Package:    Histograms
// Class:      TriggerResultsCounter
// 
/**\class TriggerResultsCounter TriggerResultsCounter.cc brot/TriggerResultsCounter/src/TriggerResultsCounter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  matthias edelhoff
//         Created:  Tue Oct 27 13:50:40 CET 2009
// $Id: TriggerResultsCounter.cc,v 1.3 2010/04/01 09:36:57 sprenger Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
//#include <FWCore/Framework/interface/TriggerReport.h>



//ROOT
#include "TH1.h"
//
// class decleration
//

class TriggerResultsCounter : public edm::EDAnalyzer {
public:
  explicit TriggerResultsCounter(const edm::ParameterSet&);
  ~TriggerResultsCounter();

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  void addTriggerSet( const edm::ParameterSet& tSet);
  edm::InputTag triggerTag_;
  std::string prefix_;

  //histos
  std::map<std::string, TH1F* > count_; 
  std::map<std::string, std::vector<std::string> > triggerNames_; 
};

// constructors and destructor
TriggerResultsCounter::TriggerResultsCounter(const edm::ParameterSet& iConfig)
{
  // read config
  //  subDir_ = iConfig.getParameter<std::string> ("subDir");
  triggerTag_ = iConfig.getParameter<edm::InputTag>("triggerTag");
  prefix_ = iConfig.getParameter<std::string> ("prefix");
  // init histos
  //  fs->mkdir(subDir_);
  std::cout << triggerTag_<<" "<<prefix_<<std::endl;
  std::vector<edm::ParameterSet> triggerSets = iConfig.getParameter<std::vector<edm::ParameterSet> >("count");
  for( std::vector<edm::ParameterSet>::const_iterator iSet = triggerSets.begin();
       iSet != triggerSets.end(); ++iSet){
    std::cout << iSet->dump() <<std::endl;
    addTriggerSet(*iSet);
  }
  //  setBinLabels(count_);
}

void TriggerResultsCounter::addTriggerSet( const edm::ParameterSet& tSet)
{

  edm::Service<TFileService> fs;
  std::string name = tSet.getParameter<std::string>("name");
  triggerNames_[name] = tSet.getParameter<std::vector<std::string> >("triggerNames");
  count_[name] = fs->make<TH1F>(name.c_str() , name.c_str() , triggerNames_[name].size() , 0 ,  triggerNames_[name].size());
  for(unsigned int i = 0; i< triggerNames_[name].size(); ++i){
    count_[name]->GetXaxis()->SetBinLabel(i+1, triggerNames_[name].at(i).c_str());
  }
}


TriggerResultsCounter::~TriggerResultsCounter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// member functions
// ------------ method called to for each event  ------------
void
TriggerResultsCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //  edm::Handle<edm::TriggerReport> triggerReport;
  //iEvent.getByLabel(triggerTag_, triggerReport);

  //for(std::vector<edm::WorkerSummary>::const_iterator iSummary = triggerReport->workerSummaries.begin();
  //    iSummary != triggerReport->workerSummaries.end(); ++iSummary){
  //  std::cout << (*iSummary).moduleLabel <<" "<<(*iSummary).timesPassed <<", " ;
  //}
  //std::cout << std::endl;

  // get hold of TriggerResults Object
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByLabel(triggerTag_, triggerResults);

  if (triggerResults.failedToGet()) {
    std::cout<<" could not get: "<<triggerTag_<<std::endl;
    return;
  }

  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
  //  triggerNames.triggerIndex(HLTPathsByName_[trig]);
  for(unsigned int i = 0; i < triggerNames.size(); ++i){
    for(std::map<std::string, std::vector<std::string> >::iterator iNameVec = triggerNames_.begin();
	iNameVec != triggerNames_.end(); ++iNameVec){
      for(std::vector<std::string>::iterator iName = iNameVec->second.begin();
	  iName != iNameVec->second.end(); ++iName){
	if( triggerNames.triggerName(i) == prefix_+(*iName) && triggerResults->accept(i)){
	  count_[iNameVec->first]->Fill((*iName).c_str(), 1.0);
	}
      }
    }
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
TriggerResultsCounter::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerResultsCounter::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerResultsCounter);
