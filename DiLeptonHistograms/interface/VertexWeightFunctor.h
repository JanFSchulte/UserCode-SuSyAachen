/*
 * VertexWeightFunctor.h
 *
 *  Created on: 29.05.2011
 *      Author: nmohr
 */

#ifndef VERTEXWEIGHTFUNCTOR_H_
#define VERTEXWEIGHTFUNCTOR_H_

#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <PhysicsTools/Utilities/interface/LumiReWeighting.h>
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>

#include "FWCore/Framework/interface/ConsumesCollector.h"

#include <iostream>

class VertexWeightFunctor
{
public:
  VertexWeightFunctor( ){}
  VertexWeightFunctor( edm::ParameterSet const & params , edm::ConsumesCollector iC ){
            std::string mcFile_ = params.getParameter<std::string>("mcFile");
            std::string mcName_ = params.getParameter<std::string>("mcName");
            std::string dataFile_ = params.getParameter<std::string>("dataFile");
            std::string dataName_ = params.getParameter<std::string>("dataName");
            std::string outputName_ = "WeightOutputs";
            int verbosity_  = params.getParameter<int>("verbosity");
            doWeight_ = params.getParameter<bool>("doWeight");
            iC.consumes<std::vector< PileupSummaryInfo > >(edm::InputTag("slimmedAddPileupInfo"));
            /*std::cout << mcFile_ << std::endl;
            std::cout << mcName_ << std::endl;
            std::cout << dataFile_ << std::endl;
            std::cout << dataName_ << std::endl;*/
            if (!verbosity_){
              std::cout.setstate(std::ios_base::failbit);
            }
            LumiWeights_ = edm::LumiReWeighting(mcFile_, dataFile_, mcName_, dataName_);
            if (!verbosity_){
              std::cout.clear();
            }
  }
  
  const double operator()(const edm::Event &iEvent){
    //std::cout << LumiWeights_.weight( iEvent ) << std::endl
    doWeight_ = !(iEvent.isRealData());
    int nTruePV = getTrueNPV(iEvent);
    return operator()(nTruePV);
  }
  const double operator()(const int nVertices){
    //std::cout << LumiWeights_.weight( nVertices ) << std::endl;
    if (!(doWeight_)) return 1.;
    else{  
      return LumiWeights_.weight( nVertices );
    }
  }
  
  int getTrueNPV(const edm::Event &iEvent){
    // from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupMCReweightingUtilities
    edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel(edm::InputTag("slimmedAddPileupInfo"), PupInfo);

    std::vector<PileupSummaryInfo>::const_iterator PVI;

    float Tnpv = -1;
    if(doWeight_){
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
  
  int BX = PVI->getBunchCrossing();
  
  if(BX == 0) { 
    Tnpv = PVI->getTrueNumInteractions();
    continue;
  }
      }
    }
    return Tnpv;
  }

private:

  edm::LumiReWeighting LumiWeights_;
  bool doWeight_;
};

#endif /* VERTEXWEIGHTFUNCTOR_H_ */
