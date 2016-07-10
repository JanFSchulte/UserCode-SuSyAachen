#ifndef eleTriggerSelector_h
#define eleTriggerSelector_h

#include "DataFormats/PatCandidates/interface/Electron.h"

//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"

//STL
#include <vector>

template<typename T, typename collectionType, typename containerType>
struct eleTriggerSelector {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  eleTriggerSelector ( const edm::ParameterSet & cfg, edm::ConsumesCollector iC ):
    HoverEEB_( cfg.getParameter<double>( "HoverEEB") ),
    HoverEEE_( cfg.getParameter<double>( "HoverEEE") ),
    deltaEtaEB_( cfg.getParameter<double>( "deltaEtaEB" ) ),
    deltaEtaEE_( cfg.getParameter<double>( "deltaEtaEE" ) ),
    deltaPhiEB_( cfg.getParameter<double>( "deltaPhiEB" ) ),
    deltaPhiEE_( cfg.getParameter<double>( "deltaPhiEE" ) ),
    eInvMinusPInvEB_( cfg.getParameter<double>( "eInvMinusPInvEB" ) ),
    eInvMinusPInvEE_( cfg.getParameter<double>( "eInvMinusPInvEE" ) ),
    sigmaIEtaIEtaEB_( cfg.getParameter<double>( "sigmaIEtaIEtaEB" ) ),
    sigmaIEtaIEtaEE_( cfg.getParameter<double>( "sigmaIEtaIEtaEE" ) )   
    
  { }
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {
	  
	  selected_.clear();
	  
	  for(typename collection::const_iterator it = col.product()->begin(); it != col.product()->end(); ++it ){
		  
		  etasc = (*it).superCluster()->eta();
		  if ( (*it).ecalEnergy() > 0.) eInvMinusPInv = 1.0/(*it).ecalEnergy() - (*it).eSuperClusterOverP()/(*it).ecalEnergy();
		  else eInvMinusPInv = 99999.;
		  
		  if ( eInvMinusPInv <= -0.05) continue; 
		  if (abs(etasc) > 1.479) {
			  if ( (*it).hadronicOverEm() >= HoverEEE_) continue;
			  if ( abs((*it).deltaEtaSuperClusterTrackAtVtx()) >= deltaEtaEE_) continue;
			  if ( abs((*it).deltaPhiSuperClusterTrackAtVtx()) >= deltaPhiEE_) continue;
			  if ( eInvMinusPInv >= eInvMinusPInvEE_) continue;
			  if ( (*it).full5x5_sigmaIetaIeta() >= sigmaIEtaIEtaEE_) continue;
		  }
		  else {
			  if ( (*it).hadronicOverEm() >= HoverEEB_) continue;
			  if ( abs((*it).deltaEtaSuperClusterTrackAtVtx()) >= deltaEtaEB_) continue;
			  if ( abs((*it).deltaPhiSuperClusterTrackAtVtx()) >= deltaPhiEB_) continue;
			  if ( eInvMinusPInv >= eInvMinusPInvEB_) continue;
			  if ( (*it).full5x5_sigmaIetaIeta() >= sigmaIEtaIEtaEB_) continue;
		  }
		  //~ if (abs(etasc) > 1.479) {
			  //~ if ( (*it).hadronicOverEm() >= 0.07) continue;
			  //~ if ( abs((*it).deltaEtaSuperClusterTrackAtVtx()) >= 0.008) continue;
			  //~ if ( abs((*it).deltaPhiSuperClusterTrackAtVtx()) >= 0.07) continue;
			  //~ if ( eInvMinusPInv >= 0.005) continue;
			  //~ if ( (*it).full5x5_sigmaIetaIeta() >= 0.03) continue;
		  //~ }
		  //~ else {
			  //~ if ( (*it).hadronicOverEm() >= 0.1) continue;
			  //~ if ( abs((*it).deltaEtaSuperClusterTrackAtVtx()) >= 0.01) continue;
			  //~ if ( abs((*it).deltaPhiSuperClusterTrackAtVtx()) >= 0.04) continue;
			  //~ if ( eInvMinusPInv >= 0.01) continue;
			  //~ if ( (*it).full5x5_sigmaIetaIeta() >= 0.011) continue;
		  //~ }
		  
		  selected_.push_back( & (*it) );
	  }
  }
  
  size_t size() const { return selected_.size(); }
private:
  container selected_;
  double etasc;
  double eInvMinusPInv;
  double HoverEEB_;
  double HoverEEE_;
  double deltaEtaEB_;
  double deltaEtaEE_;
  double deltaPhiEB_;
  double deltaPhiEE_;
  double eInvMinusPInvEB_;
  double eInvMinusPInvEE_;
  double sigmaIEtaIEtaEB_;
  double sigmaIEtaIEtaEE_;
};

#endif
