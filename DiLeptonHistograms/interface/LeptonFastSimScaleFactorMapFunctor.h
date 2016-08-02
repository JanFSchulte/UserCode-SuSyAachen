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
            

            
            std::string FastSimScaleFactorFile_mu_Id_ = params.getParameter<std::string>("FastSimScaleFactorFile_mu_Id");
            std::string FastSimScaleFactorFile_mu_Iso_ = params.getParameter<std::string>("FastSimScaleFactorFile_mu_Iso");
            std::string FastSimScaleFactorFile_mu_Impact_ = params.getParameter<std::string>("FastSimScaleFactorFile_mu_Impact");
            std::string FastSimScaleFactorFile_ele_ = params.getParameter<std::string>("FastSimScaleFactorFile_ele");
            
            std::string FastSimScaleFactorHisto_mu_Id_ = params.getParameter<std::string>("FastSimScaleFactorHisto_mu_Id");
            std::string FastSimScaleFactorHisto_mu_Iso_ = params.getParameter<std::string>("FastSimScaleFactorHisto_mu_Iso");
            std::string FastSimScaleFactorHisto_mu_Impact_ = params.getParameter<std::string>("FastSimScaleFactorHisto_mu_Impact");
            std::string FastSimScaleFactorHisto_ele_ = params.getParameter<std::string>("FastSimScaleFactorHisto_ele");
            
            
			
			getFastSimScaleFactorHistos(
								FastSimScaleFactorFile_mu_Id_, 
								FastSimScaleFactorFile_mu_Iso_, 
								FastSimScaleFactorFile_mu_Impact_, 
								FastSimScaleFactorFile_ele_, 
								
								FastSimScaleFactorHisto_mu_Id_, 
								FastSimScaleFactorHisto_mu_Iso_, 
								FastSimScaleFactorHisto_mu_Impact_, 
								FastSimScaleFactorHisto_ele_
											);
    
  }
  
  
  
  //~ const double operator()(const  pat::Electron &ele, double pt, double absEta, double nVertices){
  const double operator()(const  pat::Electron &ele, double pt, double absEta){

	  float result = 1.0;
	  
	  if(pt > 200.) pt= 199.;
	  
	  //~ result *= FastSimScaleFactorHisto_ele_->GetBinContent(FastSimScaleFactorHisto_ele_->GetXaxis()->FindBin(pt),FastSimScaleFactorHisto_ele_->GetYaxis()->FindBin(absEta),FastSimScaleFactorHisto_ele_->GetZaxis()->FindBin(nVertices));
	  result *= FastSimScaleFactorHisto_ele_->GetBinContent(FastSimScaleFactorHisto_ele_->GetXaxis()->FindBin(pt),FastSimScaleFactorHisto_ele_->GetYaxis()->FindBin(absEta));
	  
	  return result;
  }
  
  //~ const double operator()(const  pat::Muon &mu, double pt, double absEta, double nVertices){
  const double operator()(const  pat::Muon &mu, double pt, double absEta){

	  float result = 1.0;
	  
	  if(pt > 200.) pt= 199.;
	  
	  //~ result *= FastSimScaleFactorHisto_mu_Id_->GetBinContent(FastSimScaleFactorHisto_mu_Id_->GetXaxis()->FindBin(pt),FastSimScaleFactorHisto_mu_Id_->GetYaxis()->FindBin(absEta),FastSimScaleFactorHisto_mu_->GetZaxis()->FindBin(nVertices));
	  result *= FastSimScaleFactorHisto_mu_Id_->GetBinContent(FastSimScaleFactorHisto_mu_Id_->GetXaxis()->FindBin(pt),FastSimScaleFactorHisto_mu_Id_->GetYaxis()->FindBin(absEta));
	  result *= FastSimScaleFactorHisto_mu_Iso_->GetBinContent(FastSimScaleFactorHisto_mu_Iso_->GetXaxis()->FindBin(pt),FastSimScaleFactorHisto_mu_Iso_->GetYaxis()->FindBin(absEta));
	  result *= FastSimScaleFactorHisto_mu_Impact_->GetBinContent(FastSimScaleFactorHisto_mu_Impact_->GetXaxis()->FindBin(pt),FastSimScaleFactorHisto_mu_Impact_->GetYaxis()->FindBin(absEta));
	  
	  return result;
  }
  
  //~ const double operator()(const  pat::Tau &tau, double pt, double absEta, double nVertices){
  const double operator()(const  pat::Tau &tau, double pt, double absEta){

	  float result = 1.0;
	  return result;
  }
  
  
  
  

private:

  bool useFastSim_;
  
  void getFastSimScaleFactorHistos(
							std::string FastSimScaleFactorFile_mu_Id, 
							std::string FastSimScaleFactorFile_mu_Iso, 
							std::string FastSimScaleFactorFile_mu_Impact, 
							std::string FastSimScaleFactorFile_ele,
							 
							std::string FastSimScaleFactorHisto_mu_Id, 
							std::string FastSimScaleFactorHisto_mu_Iso, 
							std::string FastSimScaleFactorHisto_mu_Impact, 
							std::string FastSimScaleFactorHisto_ele
							);

  
  //~ TH3F * FastSimScaleFactorHisto_mu_;
  //~ TH3F * FastSimScaleFactorHisto_ele_;
//~ };
  TH2F * FastSimScaleFactorHisto_mu_Id_;
  TH2F * FastSimScaleFactorHisto_mu_Iso_;
  TH2F * FastSimScaleFactorHisto_mu_Impact_;
  TH2F * FastSimScaleFactorHisto_ele_;
};

#endif /* LEPTONFASTSIMSCALEFACTORMAPFUNCTOR_H_ */
