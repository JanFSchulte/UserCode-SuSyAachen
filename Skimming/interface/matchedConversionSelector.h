#ifndef matchedConversionSelector_h
#define matchedConversionSelector_h

//DataFormats
//include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"

//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//STL
#include <vector>

template<typename convType, typename bsType, typename collectionType, typename containerType>
struct matchedConversionSelector {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  matchedConversionSelector ( const edm::ParameterSet & cfg, edm::ConsumesCollector ):
    conversionsSrc_( cfg.getParameter<edm::InputTag>( "conversionsSource" ) ), 
    beamspotSrc_( cfg.getParameter<edm::InputTag>( "beamspotSource" ) )  { }
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {

    edm::Handle<convType> conversions_h;
    ev.getByLabel(conversionsSrc_, conversions_h);
    // stupod way egamma uses cmssw
    //    convType conversions = *(conversions_h.product());
    
    edm::Handle<bsType> beamspot_h;
    ev.getByLabel(beamspotSrc_, beamspot_h);
    bsType beamspot = *(beamspot_h.product());

    selected_.clear();
    for(typename collection::const_iterator it = col.product()->begin(); 
	it != col.product()->end(); ++it ){
      const pat::Electron *e = &(*(it));

      if(! ConversionTools::hasMatchedConversion(*(dynamic_cast<const reco::GsfElectron*>(e)), conversions_h, beamspot.position()))
	selected_.push_back( & (*it) );

    }
  }

  size_t size() const { return selected_.size(); }
  
  private:
  container selected_;
  edm::InputTag conversionsSrc_;
  edm::InputTag beamspotSrc_;
};

#endif
