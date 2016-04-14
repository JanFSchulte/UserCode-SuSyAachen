/*
 * LeptonScaleFactorMapFunctor.h
 *
 *  Created on: 04.12.2015
 *      Author: cschomak
 */

#ifndef LEPTONFASTSIMSCALEFACTORMAPFUNCTOR_H_
#define LEPTONFASTSIMSCALEFACTORMAPFUNCTOR_H_

#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <DataFormats/PatCandidates/interface/Lepton.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Tau.h>

#include <iostream>

#include "TFile.h"
#include "TH2.h"
#include "TH3.h"

class LeptonFastSimScaleFactorMapFunctor
{
public:
  LeptonFastSimScaleFactorMapFunctor( ){}
  LeptonFastSimScaleFactorMapFunctor( edm::ParameterSet const & params ){
            

            
            std::string FastSimScaleFactorFile_mu_ = params.getParameter<std::string>("FastSimScaleFactorFile_mu");
            std::string FastSimScaleFactorFile_ele_ = params.getParameter<std::string>("FastSimScaleFactorFile_ele");
            
            std::string FastSimScaleFactorHisto_mu_ = params.getParameter<std::string>("FastSimScaleFactorHisto_mu");
            std::string FastSimScaleFactorHisto_ele_ = params.getParameter<std::string>("FastSimScaleFactorHisto_ele");
            
            
			
			getFastSimScaleFactorHistos(
								FastSimScaleFactorFile_mu_, 
								FastSimScaleFactorFile_ele_, 
								FastSimScaleFactorHisto_mu_, 
								FastSimScaleFactorHisto_ele_
											);
    
  }
  
  
  
  const double operator()(const  pat::Electron &ele, double pt, double absEta, double nVertices){

	  float result = 1.0;
	  
	  result *= FastSimScaleFactorHisto_ele_->GetBinContent(FastSimScaleFactorHisto_ele_->GetXaxis()->FindBin(pt),FastSimScaleFactorHisto_ele_->GetYaxis()->FindBin(absEta),FastSimScaleFactorHisto_ele_->GetZaxis()->FindBin(nVertices));
	  
	  return result;
  }
  
  const double operator()(const  pat::Muon &mu, double pt, double absEta, double nVertices){

	  float result = 1.0;
	  result *= FastSimScaleFactorHisto_mu_->GetBinContent(FastSimScaleFactorHisto_mu_->GetXaxis()->FindBin(pt),FastSimScaleFactorHisto_mu_->GetYaxis()->FindBin(absEta),FastSimScaleFactorHisto_mu_->GetZaxis()->FindBin(nVertices));
	  
	  return result;
  }
  
  const double operator()(const  pat::Tau &tau, double pt, double absEta, double nVertices){

	  float result = 1.0;
	  return result;
  }
  
  
  
  

private:

  bool useFastSim_;
  
  void getFastSimScaleFactorHistos(
							std::string FastSimScaleFactorFile_mu, 
							std::string FastSimScaleFactorFile_ele, 
							std::string FastSimScaleFactorHisto_mu, 
							std::string FastSimScaleFactorHisto_ele
							);

  
  TH3F * FastSimScaleFactorHisto_mu_;
  TH3F * FastSimScaleFactorHisto_ele_;
};

#endif /* LEPTONFASTSIMSCALEFACTORMAPFUNCTOR_H_ */
