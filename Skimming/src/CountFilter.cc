// -*- C++ -*-
//
// Package:    CountFilter
// Class:      CountFilter
// 
/**\class CountFilter CountFilter.cc SuSyAachen/Skimming/src/CountFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Niklas Mohr
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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

//
// class declaration
//

template< typename T >
class CountFilter : public edm::EDFilter {
public:
  explicit CountFilter(const edm::ParameterSet&);
  ~CountFilter();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  
  edm::EDGetTokenT< T > inputToken_;

  unsigned int minN_;
  unsigned int maxN_;

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
template< typename T >
CountFilter<T>::CountFilter(const edm::ParameterSet& iConfig):
  inputToken_(consumes<T>(iConfig.getParameter<edm::InputTag>("src")))
{
  minN_ = iConfig.getParameter<unsigned int> ("minNumber");
  maxN_ = iConfig.getParameter<unsigned int> ("maxNumber");
}


template< typename T >
CountFilter<T>::~CountFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
template<typename T >
bool
CountFilter<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle< T > candidates;
  iEvent.getByToken(inputToken_, candidates);

  unsigned int number = candidates->size();
  bool result = (number >= minN_ && number <= maxN_);

  return result;
}

// ------------ method called once each job just before starting event loop  ------------
template<typename T >
void 
CountFilter<T>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template<typename T >
void 
CountFilter<T>::endJob() {
}

typedef CountFilter< reco::CandidateCollection > CandCountFilter;
typedef CountFilter< reco::VertexCollection > nVertexCountFilter;


//define this as a plug-in
DEFINE_FWK_MODULE(CandCountFilter);
DEFINE_FWK_MODULE(nVertexCountFilter);
