#ifndef d0Selector_h
#define d0Selector_h

//DataFormats
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

//STL
#include <vector>

template<typename collectionType, typename containerType>
struct d0Selector {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  d0Selector ( const edm::ParameterSet & cfg ):
    d0Min_( cfg.getParameter<double>( "d0Min" ) ),
    beamSpotSrc_( cfg.getParameter<edm::InputTag>( "beamSpotSource" ) )  { }
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {
    
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    ev.getByLabel(beamSpotSrc_, beamSpotHandle);
    
    selected_.clear();
    double d0 = 0.0;
    for(typename collection::const_iterator it = col.product()->begin(); 
	 it != col.product()->end(); ++it ){
      
      d0 = calcD0( *it, beamSpotHandle );
      if ( d0 < d0Min_ )
	selected_.push_back( & (*it) );
    }
  }
  // fast hack: this should be specialized
  double calcD0( pat::Electron p, edm::Handle<reco::BeamSpot> beamSpotHandle)
  {    
    return fabs( p.gsfTrack()->dxy( beamSpotHandle->position() ));  
  }
  // fast hack: this should be the normal one
  double calcD0( pat::Muon p, edm::Handle<reco::BeamSpot> beamSpotHandle )
  {
    return fabs( p.track()->dxy( beamSpotHandle->position() ));  
  }

  
  size_t size() const { return selected_.size(); }
private:
  container selected_;
  double d0Min_;
  edm::InputTag beamSpotSrc_;
};

#endif
