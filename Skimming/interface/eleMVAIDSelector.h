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
    workingPointCentralBarrelLowPtLinear_( cfg.getParameter<double>( "workingPointCentralBarrelLowPtLinear") ),
    workingPointOuterBarrelHighPt_( cfg.getParameter<double>( "workingPointOuterBarrelHighPt") ),
    workingPointOuterBarrelLowPt_( cfg.getParameter<double>( "workingPointOuterBarrelLowPt") ),
    workingPointOuterBarrelLowPtLinear_( cfg.getParameter<double>( "workingPointOuterBarrelLowPtLinear") ),
    workingPointEndcapHighPt_( cfg.getParameter<double>( "workingPointEndcapHighPt") ),
    workingPointEndcapLowPt_( cfg.getParameter<double>( "workingPointEndcapLowPt") ),
    workingPointEndcapLowPtLinear_( cfg.getParameter<double>( "workingPointEndcapLowPtLinear") ),
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
        float pt = (*it).pt();
        float eta = fabs((*it).eta());
        if (eta < 0.8) {
            if (pt < 25){
                workingPoint = workingPointCentralBarrelLowPt_ + workingPointCentralBarrelLowPtLinear_ * (pt-10);
            }else{
                workingPoint = workingPointCentralBarrelHighPt_;
            }
        }
        else if (eta > 0.8 && eta < 1.479) {
            if (pt < 25){
                workingPoint = workingPointOuterBarrelLowPt_ + workingPointOuterBarrelLowPtLinear_ * (pt-10);
            }else{
                workingPoint = workingPointOuterBarrelHighPt_;
            }
        }
        else {
            if (pt < 25){
                workingPoint = workingPointEndcapLowPt_ + workingPointEndcapLowPtLinear_ * (pt-10);
            }else{
                workingPoint = workingPointEndcapHighPt_;
            }
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
  float workingPointCentralBarrelLowPtLinear_;
  float workingPointOuterBarrelHighPt_;
  float workingPointOuterBarrelLowPt_;
  float workingPointOuterBarrelLowPtLinear_;
  float workingPointEndcapHighPt_;
  float workingPointEndcapLowPt_;
  float workingPointEndcapLowPtLinear_;
  edm::EDGetTokenT<edm::ValueMap<float> > idMapToken_;
};

#endif
