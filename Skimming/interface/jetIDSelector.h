#ifndef jetIDSelector_h
#define jetIDSelector_h

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/ValueMap.h"


//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"

//STL
#include <vector>

template<typename T, typename collectionType, typename containerType>
struct jetIDSelector {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  jetIDSelector ( const edm::ParameterSet & cfg, edm::ConsumesCollector ):
  year_( cfg.getParameter<int>( "year") ) {}
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {
    


    selected_.clear();

    for(typename collection::const_iterator it = col.product()->begin(); it != col.product()->end(); ++it ){


    bool isPass = true;
    auto NHF  = it->neutralHadronEnergyFraction();
    auto NEMF = it->neutralEmEnergyFraction();
    auto CHF  = it->chargedHadronEnergyFraction();
    //auto MUF  = it->muonEnergyFraction();
    auto CEMF = it->chargedEmEnergyFraction();
    auto NumConst = it->chargedMultiplicity()+it->neutralMultiplicity();
    //auto NumNeutralParticles =it->neutralMultiplicity();
    auto CHM      = it->chargedMultiplicity();
    
    if (year_ == 2016){
      // Using loose Jet ID https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
      isPass =  (NHF < 0.99) && 
                (NEMF < 0.99) && 
                (NumConst > 1) && 
                (CHF > 0) && 
                (CHM > 0) && 
                (CEMF < 0.99);
    }
    
    if (year_ == 2017){
      // Tight Jet ID now recommended https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
      isPass =  (NHF < 0.90) && 
                (NEMF < 0.90) && 
                (NumConst > 1) && 
                (CHF > 0) && 
                (CHM > 0);
    }
    
    if (year_ == 2018){
      // Only 1 Jet ID available https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
      isPass =  (NHF < 0.90) && 
                (NEMF < 0.90) && 
                (NumConst > 1) && 
                (CHF > 0) && 
                (CHM > 0);
    }
    if (isPass){
      
      selected_.push_back( & (*it) );
      }
    }
  }

    

  size_t size() const { return selected_.size(); }
private:
  container selected_;
  int year_;
};

#endif
