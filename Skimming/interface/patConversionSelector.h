#ifndef patConversionSelector_h
#define patConversionSelector_h

//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

//STL
#include <vector>

template<typename collectionType, typename containerType>
struct patConversionSelector {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  patConversionSelector ( const edm::ParameterSet & cfg ):
   cutLow_( cfg.getParameter<int>( "minLostHits") ), 
   cutHigh_( cfg.getParameter<int>( "maxLostHits" ) ) { }
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {
    
    selected_.clear();
    for(typename collection::const_iterator it = col.product()->begin(); it != col.product()->end(); ++it ){
        if(it->gsfTrack()->trackerExpectedHitsInner().numberOfHits() >= cutLow_ && it->gsfTrack()->trackerExpectedHitsInner().numberOfHits() <= cutHigh_) selected_.push_back( & (*it) );
    }
  }
  
  size_t size() const { return selected_.size(); }
private:
  container selected_;
  int cutLow_;
  int cutHigh_;
};

#endif
