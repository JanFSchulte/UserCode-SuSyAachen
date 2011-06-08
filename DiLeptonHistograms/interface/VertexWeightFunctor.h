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

#include <iostream>

class VertexWeightFunctor
{
public:
  VertexWeightFunctor( ){}
  VertexWeightFunctor( edm::ParameterSet const & params ){
            std::string mcFile_ = params.getParameter<std::string>("mcFile");
            std::string mcName_ = params.getParameter<std::string>("mcName");
            std::string dataFile_ = params.getParameter<std::string>("dataFile");
            std::string dataName_ = params.getParameter<std::string>("dataName");
            doWeight_ = params.getParameter<bool>("doWeight");
            /*std::cout << mcFile_ << std::endl;
            std::cout << mcName_ << std::endl;
            std::cout << dataFile_ << std::endl;
            std::cout << dataName_ << std::endl;*/
            LumiWeights_ = edm::LumiReWeighting(mcFile_, dataFile_, mcName_, dataName_);
  }
  const double operator()(const edm::Event &iEvent){
    //std::cout << LumiWeights_.weight( iEvent ) << std::endl;
    if (iEvent.isRealData() || !(doWeight_)) return 1.;
    else return LumiWeights_.weight( iEvent );
  }
  const double operator()(const int nVertices){
    //std::cout << LumiWeights_.weight( nVertices ) << std::endl;
    if (!doWeight_) return 1.;
    else return LumiWeights_.weight( nVertices );
  }

private:
  edm::LumiReWeighting LumiWeights_;
  bool doWeight_;
};

#endif /* VERTEXWEIGHTFUNCTOR_H_ */
