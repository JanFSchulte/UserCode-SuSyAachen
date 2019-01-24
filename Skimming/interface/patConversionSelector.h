#ifndef patConversionSelector_h
#define patConversionSelector_h

//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"

//STL
#include <vector>

template<typename collectionType, typename containerType>
struct patConversionSelector {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  patConversionSelector ( const edm::ParameterSet & cfg, edm::ConsumesCollector ):
   cutLow_( cfg.getParameter<int>( "minLostHits") ), 
   cutHigh_( cfg.getParameter<int>( "maxLostHits" ) ),
   conv_( cfg.getParameter<bool>( "convInfo") ), 
   cutDist_( cfg.getParameter<double>( "maxDistance") ), 
   cutCot_( cfg.getParameter<double>( "maxCotangentTheta" ) ) { }
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {

    edm::Handle<reco::TrackCollection> tracks_;
    ConversionFinder cf;
    if(conv_){
        ev.getByLabel("generalTracks", tracks_);
    }
    
    selected_.clear();
    //double dcot = 0;
    //double dist = 0;
    for(typename collection::const_iterator it = col.product()->begin(); it != col.product()->end(); ++it ){
        if(it->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) >= cutLow_ && it->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= cutHigh_  ) selected_.push_back( & (*it) );
    }
  }
  
  size_t size() const { return selected_.size(); }
private:
  container selected_;
  int cutLow_;
  int cutHigh_;
  bool conv_;
  double cutDist_;
  double cutCot_;
};

#endif
