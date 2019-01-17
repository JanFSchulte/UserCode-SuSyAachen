#ifndef eleLooseMVAIDSelector_h
#define eleLooseMVAIDSelector_h

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/ValueMap.h"


//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"

//STL
#include <vector>

// OBSOLETE, currently only eleMVAIDSelector is used, even for loose electrons!

template<typename T, typename collectionType, typename containerType>
struct eleLooseMVAIDSelector {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  eleLooseMVAIDSelector ( const edm::ParameterSet & cfg, edm::ConsumesCollector iC ):
    idMapToken_(iC.consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>( "idMapSource" )))  { }
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {
    
    edm::Handle<edm::ValueMap<float> > id_values;
    ev.getByToken(idMapToken_, id_values);


    selected_.clear();
    for(typename collection::const_iterator it = col.product()->begin(); it != col.product()->end(); ++it ){
    const edm::Ptr<pat::Electron> elPtr(col, it - col->begin() );
    float mvaValue  = (*id_values)[ elPtr ];

    float workingPoint = -9999.;
    float eta = fabs((*it).superCluster()->eta();
    if (eta < 0.8) workingPoint =  -0.70;
    else if (eta < 1.479) workingPoint = -0.83;
    else workingPoint = -0.92;
    if (mvaValue > workingPoint){
      selected_.push_back( & (*it) );
      }
    }
  }

    

  size_t size() const { return selected_.size(); }
private:
  container selected_;
  edm::EDGetTokenT<edm::ValueMap<float> > idMapToken_;
};

#endif
