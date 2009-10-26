// this is broken and not usefull since electrons and muons do need different treatment anyway :/


#ifndef d0Selector_h
#define d0Selector_h

//DataFormats
#include "DataFormats/PatCandidates/interface/Muon.h"
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
  //  typedef pat::MuonCollection collection;
  //  typedef std::vector<const pat::Muon *> container;
  typedef typename collection::const_iterator const_iterator;
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
    for( const_iterator it = col.product()->begin(); 
	 it != col.product()->end(); ++it ){
      d0 = fabs( (*it).innerTrack()->dxy( beamSpotHandle->position() ));
      if ( d0 < d0Min_ )
	selected_.push_back( & (*it) );
    }
  }
  
  size_t size() const { return selected_.size(); }
private:
  container selected_;
  double d0Min_;
  edm::InputTag beamSpotSrc_;
};

#endif
