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
#include <PhysicsTools/Utilities/interface/Lumi3DReWeighting.h>

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
            std::string mc3DFile_ = params.getParameter<std::string>("mc3DFile");
            std::string mc3DName_ = params.getParameter<std::string>("mc3DName");
            std::string data3DFile_ = params.getParameter<std::string>("data3DFile");
            std::string data3DName_ = params.getParameter<std::string>("data3DName");
            doWeight_ = params.getParameter<bool>("doWeight");
            doWeight3D_ = params.getParameter<bool>("doWeight3D");
            fractionRunA_ = params.getParameter<double>("fractionRunA");
            fractionRunB_ = params.getParameter<double>("fractionRunB");
            /*std::cout << mcFile_ << std::endl;
            std::cout << mcName_ << std::endl;
            std::cout << dataFile_ << std::endl;
            std::cout << dataName_ << std::endl;*/
            LumiWeights_ = edm::LumiReWeighting(mcFile_, dataFile_, mcName_, dataName_);
            if (doWeight3D_) {
                LumiWeights3D_ = edm::Lumi3DReWeighting(mc3DFile_, data3DFile_, mc3DName_, data3DName_);
                LumiWeights3D_.weight3D_init(1.);
            }
    
  }
  const double operator()(const edm::Event &iEvent){
    //std::cout << LumiWeights_.weight( iEvent ) << std::endl;
    if (iEvent.isRealData() || !(doWeight_)) return 1.;
    else {
        if (doWeight3D_) return (fractionRunA_*LumiWeights_.weight( iEvent )+fractionRunB_*LumiWeights3D_.weight3D( iEvent ));
        else return LumiWeights_.weight( iEvent ); 
    }
  }
  const double operator()(const int nVertices){
    //std::cout << LumiWeights_.weight( nVertices ) << std::endl;
    if (!doWeight_) return 1.;
    else return LumiWeights_.weight( nVertices );
  }

private:
  edm::LumiReWeighting LumiWeights_;
  edm::Lumi3DReWeighting LumiWeights3D_;
  bool doWeight_;
  bool doWeight3D_;
  double fractionRunA_;
  double fractionRunB_;
};

#endif /* VERTEXWEIGHTFUNCTOR_H_ */
