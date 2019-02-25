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
    lowPtHighPtCutOff_( cfg.getParameter<double>( "lowPtHighPtCutOff") ), 
    lowPtLinearSubtraction_( cfg.getParameter<double>( "lowPtLinearSubtraction") ), 
    idMapName_(cfg.getParameter<std::string>( "idMapSource" ) )  { }
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {
    selected_.clear();
    for(typename collection::const_iterator it = col.product()->begin(); it != col.product()->end(); ++it ){
        float mvaValue = (*it).userFloat(idMapName_);
        //for (const auto &userF : it->userFloatNames()){
            //std::cout << userF << std::endl;
        //}
        
        //std::cout << "pt: " << (*it).pt() << " eta: " << (*it).eta() << " MVA value: " << mvaValue << endl;
        float workingPoint = -9999.;
        float pt = (*it).pt();
        float eta = fabs((*it).superCluster()->eta());
        
        if (eta < 0.8) {
            if (pt < lowPtHighPtCutOff_){
                workingPoint = workingPointCentralBarrelLowPt_ + workingPointCentralBarrelLowPtLinear_ * (pt-lowPtLinearSubtraction_);
            }else{
                workingPoint = workingPointCentralBarrelHighPt_;
            }
        }
        else if (eta > 0.8 && eta < 1.479) {
            if (pt < lowPtHighPtCutOff_){
                workingPoint = workingPointOuterBarrelLowPt_ + workingPointOuterBarrelLowPtLinear_ * (pt-lowPtLinearSubtraction_);
            }else{
                workingPoint = workingPointOuterBarrelHighPt_;
            }
        }
        else {
            if (pt < lowPtHighPtCutOff_){
                workingPoint = workingPointEndcapLowPt_ + workingPointEndcapLowPtLinear_ * (pt-lowPtLinearSubtraction_);
            }else{
                workingPoint = workingPointEndcapHighPt_;
            }
        }
        if (mvaValue > workingPoint){
            selected_.push_back( & (*it) );
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
  float lowPtHighPtCutOff_;
  float lowPtLinearSubtraction_;
  std::string idMapName_;
};

#endif
