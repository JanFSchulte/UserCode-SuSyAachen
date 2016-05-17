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

template<typename T, typename collectionType, typename containerType>
struct eleLooseMVAIDSelector {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  eleLooseMVAIDSelector ( const edm::ParameterSet & cfg, edm::ConsumesCollector iC ):
    //~ idMapSrc_( cfg.getParameter<edm::InputTag>( "idMapSource" ) )  { }
    idMapToken_(iC.consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>( "idMapSource" )))  { }
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {
    
    edm::Handle<edm::ValueMap<float> > id_values;
    //~ ev.getByLabel(idMapSrc_, id_values);
    ev.getByToken(idMapToken_, id_values);


    selected_.clear();
	//std::cout << "--------------------------" << std::endl;
	//std::cout << ev.id().event() << endl;
    //~ double etaValue;
    for(typename collection::const_iterator it = col.product()->begin(); it != col.product()->end(); ++it ){
		const edm::Ptr<pat::Electron> elPtr(col, it - col->begin() );
		float mvaValue  = (*id_values)[ elPtr ];

		//std::cout << "pt: " << (*it).pt() << " eta: " << (*it).eta() << " MVA value: " << mvaValue << endl;
		float workingPoint = -9999.;
		if (fabs((*it).eta()) < 0.8) workingPoint =  0.35;
		else if (fabs((*it).eta()) > 0.8 && fabs((*it).eta()) < 1.479) workingPoint = 0.20;
		else workingPoint = -0.52;
		if (mvaValue > workingPoint){
			//std::cout << "pushing back" << std::endl;
			selected_.push_back( & (*it) );
    	}
    }
  }

    

  size_t size() const { return selected_.size(); }
private:
  container selected_;
  //~ edm::InputTag idMapSrc_;
  edm::EDGetTokenT<edm::ValueMap<float> > idMapToken_;
};

#endif
