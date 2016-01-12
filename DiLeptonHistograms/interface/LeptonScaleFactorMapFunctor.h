/*
 * LeptonScaleFactorMapFunctor.h
 *
 *  Created on: 04.12.2015
 *      Author: cschomak
 */

#ifndef LEPTONSCALEFACTORMAPFUNCTOR_H_
#define LEPTONSCALEFACTORMAPFUNCTOR_H_

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

class LeptonScaleFactorMapFunctor
{
public:
  LeptonScaleFactorMapFunctor( ){}
  LeptonScaleFactorMapFunctor( edm::ParameterSet const & params ){
            
            std::string dataMCScaleFactorFile_mu_ID_ = params.getParameter<std::string>("dataMCScaleFactorFile_mu_ID");
            std::string dataMCScaleFactorFile_mu_Iso_ = params.getParameter<std::string>("dataMCScaleFactorFile_mu_Iso");
            std::string dataMCScaleFactorFile_ele_ = params.getParameter<std::string>("dataMCScaleFactorFile_ele");
            
            std::string dataMCScaleFactorHisto_mu_ID_ = params.getParameter<std::string>("dataMCScaleFactorHisto_mu_ID");
            std::string dataMCScaleFactorHisto_mu_Iso_ = params.getParameter<std::string>("dataMCScaleFactorHisto_mu_Iso");
            std::string dataMCScaleFactorHisto_ele_ID_ = params.getParameter<std::string>("dataMCScaleFactorHisto_ele_ID");
            std::string dataMCScaleFactorHisto_ele_Iso_ = params.getParameter<std::string>("dataMCScaleFactorHisto_ele_Iso");
            
            std::string FastSimScaleFactorFile_mu_ = params.getParameter<std::string>("FastSimScaleFactorFile_mu");
            std::string FastSimScaleFactorFile_ele_ = params.getParameter<std::string>("FastSimScaleFactorFile_ele");
            
            std::string FastSimScaleFactorHisto_mu_ = params.getParameter<std::string>("FastSimScaleFactorHisto_mu");
            std::string FastSimScaleFactorHisto_ele_ = params.getParameter<std::string>("FastSimScaleFactorHisto_ele");
            
            
			
			getScaleFactorHistos(dataMCScaleFactorFile_mu_ID_, 
								dataMCScaleFactorFile_mu_Iso_, 
								dataMCScaleFactorFile_ele_, 
								dataMCScaleFactorHisto_mu_ID_, 
								dataMCScaleFactorHisto_mu_Iso_, 
								dataMCScaleFactorHisto_ele_ID_, 
								dataMCScaleFactorHisto_ele_Iso_, 
								FastSimScaleFactorFile_mu_, 
								FastSimScaleFactorFile_ele_, 
								FastSimScaleFactorHisto_mu_, 
								FastSimScaleFactorHisto_ele_
											);
    
  }
  
  
  
  const double operator()(const  pat::Electron &ele, double pt, double absEta, double nVertices){

	  float result = 1.0;
	  
	  result *= dataMCScaleFactorHisto_ele_ID_->GetBinContent(dataMCScaleFactorHisto_ele_ID_->GetXaxis()->FindBin(pt),dataMCScaleFactorHisto_ele_ID_->GetYaxis()->FindBin(absEta));
	  result *= dataMCScaleFactorHisto_ele_Iso_->GetBinContent(dataMCScaleFactorHisto_ele_Iso_->GetXaxis()->FindBin(pt),dataMCScaleFactorHisto_ele_Iso_->GetYaxis()->FindBin(absEta));
	  
	  result *= FastSimScaleFactorHisto_ele_->GetBinContent(FastSimScaleFactorHisto_ele_->GetXaxis()->FindBin(pt),FastSimScaleFactorHisto_ele_->GetYaxis()->FindBin(absEta),FastSimScaleFactorHisto_ele_->GetZaxis()->FindBin(nVertices));
	  
	  return result;
  }
  
  const double operator()(const  pat::Muon &mu, double pt, double absEta, double nVertices){

	  float result = 1.0;
	  
	  result *= dataMCScaleFactorHisto_mu_ID_->GetBinContent(dataMCScaleFactorHisto_mu_ID_->GetXaxis()->FindBin(pt),dataMCScaleFactorHisto_mu_ID_->GetYaxis()->FindBin(absEta));
	  result *= dataMCScaleFactorHisto_mu_Iso_->GetBinContent(dataMCScaleFactorHisto_mu_Iso_->GetXaxis()->FindBin(pt),dataMCScaleFactorHisto_mu_Iso_->GetYaxis()->FindBin(absEta));
	  
	  result *= FastSimScaleFactorHisto_mu_->GetBinContent(FastSimScaleFactorHisto_mu_->GetXaxis()->FindBin(pt),FastSimScaleFactorHisto_mu_->GetYaxis()->FindBin(absEta),FastSimScaleFactorHisto_mu_->GetZaxis()->FindBin(nVertices));
	  
	  return result;
  }
  
  const double operator()(const  pat::Tau &tau, double pt, double absEta, double nVertices){

	  float result = 1.0;
	  return result;
  }
  
  
  
  

private:

  bool useFastSim_;
  
  void getScaleFactorHistos(std::string dataMCScaleFactorFile_mu_ID, 
							std::string dataMCScaleFactorFile_mu_Iso, 
							std::string dataMCScaleFactorFile_ele, 
							std::string dataMCScaleFactorHisto_mu_ID, 
							std::string dataMCScaleFactorHisto_mu_Iso, 
							std::string dataMCScaleFactorHisto_ele_ID, 
							std::string dataMCScaleFactorHisto_ele_Iso, 
							std::string FastSimScaleFactorFile_mu, 
							std::string FastSimScaleFactorFile_ele, 
							std::string FastSimScaleFactorHisto_mu, 
							std::string FastSimScaleFactorHisto_ele
							);

  TH2F * dataMCScaleFactorHisto_mu_ID_;
  TH2F * dataMCScaleFactorHisto_mu_Iso_;
  TH2F * dataMCScaleFactorHisto_ele_ID_;
  TH2F * dataMCScaleFactorHisto_ele_Iso_;
  
  TH3F * FastSimScaleFactorHisto_mu_;
  TH3F * FastSimScaleFactorHisto_ele_;
};

#endif /* LEPTONSCALEFACTORMAPFUNCTOR_H_ */
