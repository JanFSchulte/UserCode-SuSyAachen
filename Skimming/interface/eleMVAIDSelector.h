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
    workingPointCentralBarrelHighPt_( cfg.getParameter<double>( "workingPointCentralBarrelHighPt") ),
    workingPointCentralBarrelLowPt_( cfg.getParameter<double>( "workingPointCentralBarrelLowPt") ),
    workingPointOuterBarrelHighPt_( cfg.getParameter<double>( "workingPointOuterBarrelHighPt") ),
    workingPointOuterBarrelLowPt_( cfg.getParameter<double>( "workingPointOuterBarrelLowPt") ),
    workingPointEndcapHighPt_( cfg.getParameter<double>( "workingPointEndcapHighPt") ),
    workingPointEndcapLowPt_( cfg.getParameter<double>( "workingPointEndcapLowPt") ),
    idMapToken_(iC.consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>( "idMapSource" )))  { }
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {
    
    edm::Handle<edm::ValueMap<float> > id_values;
    ev.getByToken(idMapToken_, id_values);
    
    float slope;
    //~ int eventNr;


    selected_.clear();
    for(typename collection::const_iterator it = col.product()->begin(); it != col.product()->end(); ++it ){
		const edm::Ptr<pat::Electron> elPtr(col, it - col->begin() );
		float mvaValue  = (*id_values)[ elPtr ];

		//std::cout << "pt: " << (*it).pt() << " eta: " << (*it).eta() << " MVA value: " << mvaValue << endl;
		float workingPoint = -9999.;
		if (fabs((*it).eta()) < 0.8) {
			slope = (workingPointCentralBarrelHighPt_ - workingPointCentralBarrelLowPt_) / 10.;
			workingPoint =  std::min(workingPointCentralBarrelLowPt_, std::max(workingPointCentralBarrelHighPt_, workingPointCentralBarrelLowPt_ + slope * ( static_cast<float>((*it).pt() - 15.))));
		}
		else if (fabs((*it).eta()) > 0.8 && fabs((*it).eta()) < 1.479) {
			slope = (workingPointOuterBarrelHighPt_ - workingPointOuterBarrelLowPt_) / 10.;
			workingPoint =  std::min(workingPointOuterBarrelLowPt_, std::max(workingPointOuterBarrelHighPt_, workingPointOuterBarrelLowPt_ + slope * ( static_cast<float>((*it).pt() - 15.))));
		}
		else {
			slope = (workingPointEndcapHighPt_ - workingPointEndcapLowPt_) / 10.;
			workingPoint =  std::min(workingPointEndcapLowPt_, std::max(workingPointEndcapHighPt_, workingPointEndcapLowPt_ + slope * ( static_cast<float>((*it).pt() - 15.))));
		}
		if (mvaValue > workingPoint){
			selected_.push_back( & (*it) );
			//~ eventNr = ev.id().event();
	    	//~ eventNr = (eventNr>0?eventNr:eventNr+4294967296);
	    	//~ if (eventNr == 11933833 || eventNr == 15994148 || eventNr == 16656578 || eventNr == 15817262 || eventNr == 13134895 || eventNr == 12128012 ||  eventNr == 13210877 ||  eventNr == 11933854){
				//~ std::cout << "eventNr: " << eventNr << endl;
				//~ std::cout << "pt: " << (*it).pt() << " eta: " << (*it).eta() << " MVA value: " << mvaValue << " working Point: " << workingPoint << endl;
			//~ }
    	}
    	
    }
  }

    

  size_t size() const { return selected_.size(); }
private:
  container selected_;
  float workingPointCentralBarrelHighPt_;
  float workingPointCentralBarrelLowPt_;
  float workingPointOuterBarrelHighPt_;
  float workingPointOuterBarrelLowPt_;
  float workingPointEndcapHighPt_;
  float workingPointEndcapLowPt_;
  edm::EDGetTokenT<edm::ValueMap<float> > idMapToken_;
};

#endif
