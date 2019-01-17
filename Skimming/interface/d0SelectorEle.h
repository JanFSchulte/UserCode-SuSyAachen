#ifndef d0SelectorEle_h
#define d0SelectorEle_h

//DataFormats
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"

//STL
#include <vector>

template<typename T, typename collectionType, typename containerType>
struct d0SelectorEle {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  d0SelectorEle ( const edm::ParameterSet & cfg, edm::ConsumesCollector iC ):
    d0Min_( cfg.getParameter<double>( "d0Min") ),
    d0MaxEB_( cfg.getParameter<double>( "d0MaxEB" ) ),
    d0MaxEE_( cfg.getParameter<double>( "d0MaxEE" ) ),
    dZMin_( cfg.getParameter<double>( "dZMin") ),
    dZMaxEB_( cfg.getParameter<double>( "dZMaxEB" ) ),
    dZMaxEE_( cfg.getParameter<double>( "dZMaxEE" ) ),
    SIP3DMin_( cfg.getParameter<double>( "SIP3DMin" ) ),
    SIP3DMax_( cfg.getParameter<double>( "SIP3DMax" ) ),
    beamSpotToken_(iC.consumes<T>(cfg.getParameter<edm::InputTag>( "beamSpotSource" )))  { }
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {
    
    edm::Handle<T> beamSpotHandle;
    ev.getByToken(beamSpotToken_, beamSpotHandle);
    
    point_ = getPoint(beamSpotHandle);

    selected_.clear();
    std::pair<double, double> dS (0.,0.);
    double SIP3D = 0.;
    for(typename collection::const_iterator it = col.product()->begin(); 
   it != col.product()->end(); ++it ){
      
      dS = calcDs( *it, point_);
      SIP3D = calcSIP3D( *it);
      if ( dS.first >= d0Min_ && dS.first < d0MaxEB_ && dS.second >= dZMin_ && dS.second < dZMaxEB_ && SIP3D < SIP3DMax_ && SIP3D >= SIP3DMin_ && (*it).isEB() )
  selected_.push_back( & (*it) );
      if ( dS.first >= d0Min_ && dS.first < d0MaxEE_ && dS.second >= dZMin_ && dS.second < dZMaxEE_  && SIP3D < SIP3DMax_ && SIP3D >= SIP3DMin_&& (*it).isEE() )
  selected_.push_back( & (*it) );
    }
  }
  
  // fast hack: this should be specialized
  std::pair<double, double> calcDs( pat::Electron p, math::XYZPoint vx)
  {    
    return std::make_pair(std::abs( p.gsfTrack()->dxy( vx )), std::abs( p.gsfTrack()->dz( vx )));  
    //return std::make_pair(std::abs( p.gsfTrack()->dxy( vx )), std::abs( p.dB(pat::Electron::PVDZ)));  
  }
  
  double calcSIP3D( pat::Electron p)
  {
    return std::abs(p.dB(pat::Electron::PV3D)/p.edB(pat::Electron::PV3D));  
  }
  
  // fast hack: this should be specialized
  math::XYZPoint getPoint( edm::Handle<reco::BeamSpot> beamSpotHandle)
  { 
      return beamSpotHandle->position();  
  }
  // fast hack: this should be the normal one
  math::XYZPoint getPoint( edm::Handle<reco::VertexCollection> vertexHandle )
  {
      math::XYZPoint pv;
      for (reco::VertexCollection::const_iterator it = vertexHandle->begin(); it != vertexHandle->end(); ++it) {
            pv = it->position();
            break;
      }
      return pv;
  }
    

  size_t size() const { return selected_.size(); }
private:
  container selected_;
  double d0Min_;
  double d0MaxEB_;
  double d0MaxEE_;
  double dZMin_;
  double dZMaxEB_;
  double dZMaxEE_;
  double SIP3DMin_;
  double SIP3DMax_;
  math::XYZPoint point_;
  edm::EDGetTokenT<T> beamSpotToken_;
};

#endif
