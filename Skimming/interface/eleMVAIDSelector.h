#ifndef eleMVAIDSelector_h
#define eleMVAIDSelector_h

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/ValueMap.h"


//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"

//STL
#include <vector>

template<typename T, typename collectionType, typename containerType>
struct eleMVAIDSelector {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  eleMVAIDSelector ( const edm::ParameterSet & cfg, edm::ConsumesCollector iC):
    workingPointCentralBarrel_( cfg.getParameter<double>( "workingPointCentralBarrel") ),
    workingPointOuterBarrel_( cfg.getParameter<double>( "workingPointOuterBarrel") ),
    workingPointEndcap_( cfg.getParameter<double>( "workingPointEndcap") ),
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

		//std::cout << "pt: " << (*it).pt() << " eta: " << (*it).eta() << " MVA value: " << mvaValue << endl;
		float workingPoint = -9999.;
		if (fabs((*it).eta()) < 0.8) workingPoint =  workingPointCentralBarrel_;
		else if (fabs((*it).eta()) > 0.8 && fabs((*it).eta()) < 1.479) workingPoint = workingPointOuterBarrel_;
		else workingPoint = workingPointEndcap_;
		if (mvaValue > workingPoint){
			selected_.push_back( & (*it) );
    	}
    }
  }

    

  size_t size() const { return selected_.size(); }
private:
  container selected_;
  float workingPointCentralBarrel_;
  float workingPointOuterBarrel_;
  float workingPointEndcap_;
  edm::EDGetTokenT<edm::ValueMap<float> > idMapToken_;
};

#endif
