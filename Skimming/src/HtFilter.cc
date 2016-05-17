// -*- C++ -*-
//
// Package:    HtFilter
// Class:      HtFilter
// 
/**\class HtFilter HtFilter.cc SuSyAachen/Skimming/src/HtFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Matthias Edelhoff
//         Created:  Mon Nov 16 11:26:19 CET 2009
// $Id: HtFilter.cc,v 1.2 2010/06/28 13:33:29 nmohr Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"


#include "DataFormats/Candidate/interface/Candidate.h"

//
// class declaration
//

class HtFilter : public edm::EDFilter {
public:
  explicit HtFilter(const edm::ParameterSet&);
  ~HtFilter();
  
private:
  typedef reco::Candidate cand;
  typedef edm::View<cand> collection;
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  
  //~ edm::InputTag inputTag_;
  edm::EDGetTokenT< collection > inputToken_;

  double minHT_;
  double maxHT_;

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
HtFilter::HtFilter(const edm::ParameterSet& iConfig):
  inputToken_(consumes< collection >(iConfig.getParameter<edm::InputTag>("src")))
{
  //~ inputTag_ = iConfig.getParameter<edm::InputTag> ("src");
  minHT_ = iConfig.getParameter<double> ("minHT");
  maxHT_ = iConfig.getParameter<double> ("maxHT");
}


HtFilter::~HtFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HtFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle< collection > candidates;
  //~ iEvent.getByLabel(inputTag_, candidates);
  iEvent.getByToken(inputToken_, candidates);

  bool result = false;
  double ht = 0.0;
  for(collection::const_iterator it = candidates->begin(); it != candidates->end() ; ++it){
    ht += (*it).pt();
  }
  if (maxHT_ > 0. ) result = (ht >=minHT_ && ht <= maxHT_);
  else result = ht >=minHT_;

  return result;
}

// ------------ method called once each job just before starting event loop  ------------
void 
HtFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HtFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HtFilter);
