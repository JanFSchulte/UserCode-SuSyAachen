// -*- C++ -*-
//
// Package:    DiLeptonFilter
// Class:      DiLeptonFilter
// 
/**\class DiLeptonFilter DiLeptonFilter.cc SuSyAachen/Skimming/src/DiLeptonFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Matthias Edelhoff
//         Created:  Mon Nov 16 11:26:19 CET 2009
// $Id: DiLeptonFilter.cc,v 1.3 2009/11/18 10:58:25 edelhoff Exp $
//
//


// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"


#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"

//include "CommonTools/UtilAlgos/interface/DeltaR.h"


//
// class declaration
//

class DiLeptonFilter : public edm::EDFilter {
public:
  explicit DiLeptonFilter(const edm::ParameterSet&);
  ~DiLeptonFilter();
  
private:
  typedef reco::Candidate cand;
  typedef edm::View<cand> collection;
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  bool inclusiveFilter( const collection &primary, const collection &secondary);
  bool different(const cand &a, const cand &b);

  // ----------member data ---------------------------
  edm::InputTag primaryTag_;
  edm::InputTag secondaryTag_;

  bool sameSign_;
  bool matching_;
  std::string method_;

  double minDR_;
  double minDpt_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DiLeptonFilter::DiLeptonFilter(const edm::ParameterSet& iConfig)
{
  primaryTag_ = iConfig.getParameter<edm::InputTag> ("primarySrc");
  secondaryTag_ = iConfig.getParameter<edm::InputTag> ("secondarySrc");
  sameSign_ = iConfig.getParameter<bool> ("sameSign");
  matching_ = iConfig.getParameter<bool> ("matching");
  method_ = iConfig.getParameter<std::string> ("method");
  minDR_ = iConfig.getParameter<double> ("minDR");
  minDpt_ = iConfig.getParameter<double> ("minDpt");

}


DiLeptonFilter::~DiLeptonFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
DiLeptonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool result = false;
  
  edm::Handle< collection > primary;
  iEvent.getByLabel(primaryTag_, primary);

  edm::Handle< collection > secondary;
  iEvent.getByLabel(secondaryTag_, secondary);

  if(method_ == "inclusive")
    result = inclusiveFilter( *primary, *secondary);
  else if (method_ == "exclusive"){
    int overlap = 0;
    for(collection::const_iterator iPrime = primary->begin(); iPrime != primary->end() ; ++iPrime)
      for(collection::const_iterator iSec = secondary->begin(); iSec != secondary->end() ; ++iSec)
	if(!different( (*iPrime), (*iSec)))
	  overlap++;
    result = inclusiveFilter( *primary, *secondary);
    result = result && (( primary->size() + secondary->size() - overlap) == 2);
  }else
    throw new cms::Exception("unknown method: "+method_);
  
  return result;
}

bool DiLeptonFilter::inclusiveFilter( const collection &primary, const collection &secondary)
{
  bool result = false;
  int countPrime = 0;
  for(collection::const_iterator iPrime = primary.begin(); iPrime != primary.end() ; ++iPrime, ++countPrime){
    int countSec = 0;
    for(collection::const_iterator iSec = secondary.begin(); iSec != secondary.end() ; ++iSec, ++countSec){
      if( ( !(primaryTag_ == secondaryTag_) || countPrime != countSec)
	  &&( different( (*iPrime), (*iSec) ) ) 
	  &&( fabs((*iPrime).charge()) == 1 && fabs((*iSec).charge()) == 1 ) ){
	if( sameSign_ ){
	  result = result || ( (*iPrime).charge() == (*iSec).charge() );
	}
	else
	  result = result || ( (*iPrime).charge() != (*iSec).charge() );
	//	std::cout << (*iPrime).charge() << (*iSec).charge() <<"| ";
      }else{
	//	std::cout << "O| ";
      }
    }
    //    std::cout << "   ";
  }
  //  std::cout << primary.size()<<" result: "<< result;
  return result;
}

bool DiLeptonFilter::different(const cand &a, const cand &b)
{
  if(matching_)
    return &a != &b;
  else
    return reco::deltaR<const cand, const cand>( a, b) > minDR_ 
      && fabs( a.pt() - b.pt() ) > minDpt_ ;
}

// ------------ method called once each job just before starting event loop  ------------
void 
DiLeptonFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiLeptonFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiLeptonFilter);
