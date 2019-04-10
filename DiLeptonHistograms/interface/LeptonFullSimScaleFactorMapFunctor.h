/*
 * LeptonFullSimScaleFactorMapFunctor.h
 *
 *  Created on: 04.12.2015
 *      Author: cschomak
 */

#ifndef LEPTONFULLSIMSCALEFACTORMAPFUNCTOR_H_
#define LEPTONFULLSIMSCALEFACTORMAPFUNCTOR_H_

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

class LeptonFullSimScaleFactorMapFunctor
{
public:
  LeptonFullSimScaleFactorMapFunctor( ){}
  LeptonFullSimScaleFactorMapFunctor( edm::ParameterSet const & params ){
            
            std::string dataMCScaleFactorFile_mu_ID_ = params.getParameter<std::string>("dataMCScaleFactorFile_mu_ID");
            std::string dataMCScaleFactorFile_mu_Iso_ = params.getParameter<std::string>("dataMCScaleFactorFile_mu_Iso");
            //std::string dataMCScaleFactorFile_mu_Impact_ = params.getParameter<std::string>("dataMCScaleFactorFile_mu_Impact");
            //std::string dataMCScaleFactorFile_mu_Track_ = params.getParameter<std::string>("dataMCScaleFactorFile_mu_Track");
            //std::string dataMCScaleFactorFile_mu_SIP3D_ = params.getParameter<std::string>("dataMCScaleFactorFile_mu_SIP3D");
            
            std::string dataMCScaleFactorFile_ele_ = params.getParameter<std::string>("dataMCScaleFactorFile_ele");
            std::string dataMCScaleFactorFile_ele_Reco_ = params.getParameter<std::string>("dataMCScaleFactorFile_ele_Reco");
            //std::string dataMCScaleFactorFile_ele_Track_ = params.getParameter<std::string>("dataMCScaleFactorFile_ele_Track");
            
            std::string dataMCScaleFactorHisto_mu_ID_ = params.getParameter<std::string>("dataMCScaleFactorHisto_mu_ID");
            std::string dataMCScaleFactorHisto_mu_Iso_ = params.getParameter<std::string>("dataMCScaleFactorHisto_mu_Iso");
            //std::string dataMCScaleFactorHisto_mu_Impact_ = params.getParameter<std::string>("dataMCScaleFactorHisto_mu_Impact");
            //std::string dataMCScaleFactorHisto_mu_Track_ = params.getParameter<std::string>("dataMCScaleFactorHisto_mu_Track");
            //std::string dataMCScaleFactorHisto_mu_SIP3D_ = params.getParameter<std::string>("dataMCScaleFactorHisto_mu_SIP3D");
            
            std::string dataMCScaleFactorHisto_ele_Reco_ = params.getParameter<std::string>("dataMCScaleFactorHisto_ele_Reco");
            std::string dataMCScaleFactorHisto_ele_ID_ = params.getParameter<std::string>("dataMCScaleFactorHisto_ele_ID");
            std::string dataMCScaleFactorHisto_ele_Iso_ = params.getParameter<std::string>("dataMCScaleFactorHisto_ele_Iso");
            std::string dataMCScaleFactorHisto_ele_ConvHit_ = params.getParameter<std::string>("dataMCScaleFactorHisto_ele_ConvHit");
            //std::string dataMCScaleFactorHisto_ele_Track_ = params.getParameter<std::string>("dataMCScaleFactorHisto_ele_Track");
            
            
      
      getFullSimScaleFactorHistos(dataMCScaleFactorFile_mu_ID_, 
                dataMCScaleFactorFile_mu_Iso_, 
                //dataMCScaleFactorFile_mu_Impact_, 
                //~ dataMCScaleFactorFile_mu_Track_, 
                //dataMCScaleFactorFile_mu_SIP3D_, 
                
                dataMCScaleFactorFile_ele_, 
                dataMCScaleFactorFile_ele_Reco_, 
                //dataMCScaleFactorFile_ele_Track_, 
                
                dataMCScaleFactorHisto_mu_ID_, 
                dataMCScaleFactorHisto_mu_Iso_, 
                //dataMCScaleFactorHisto_mu_Impact_, 
                //~ dataMCScaleFactorHisto_mu_Track_, 
                //dataMCScaleFactorHisto_mu_SIP3D_, 
                
                dataMCScaleFactorHisto_ele_Reco_, 
                dataMCScaleFactorHisto_ele_ID_, 
                dataMCScaleFactorHisto_ele_Iso_,
                dataMCScaleFactorHisto_ele_ConvHit_
                //dataMCScaleFactorHisto_ele_Track_
                      );
    
  }
  
  
  
  //~ const double operator()(const  pat::Electron &ele, double pt, double eta, double nVertices){
  const double operator()(const  pat::Electron &ele, double pt, double eta){

    float result = 1.0;
    //double tempPtTrack;
        
    //if (pt > 500.)tempPtTrack = 499.;
    //else if (pt < 25.)tempPtTrack = 26.;
    //else tempPtTrack = pt;
    
    eta = ele.superCluster()->eta();
    if(pt > 200.) pt= 199.;
    
    result *= dataMCScaleFactorHisto_ele_Reco_->GetBinContent(dataMCScaleFactorHisto_ele_Reco_->GetXaxis()->FindBin(eta), dataMCScaleFactorHisto_ele_Reco_->GetYaxis()->FindBin(pt));
    result *= dataMCScaleFactorHisto_ele_ID_->GetBinContent(dataMCScaleFactorHisto_ele_ID_->GetXaxis()->FindBin(eta), dataMCScaleFactorHisto_ele_ID_->GetYaxis()->FindBin(pt));
    result *= dataMCScaleFactorHisto_ele_Iso_->GetBinContent(dataMCScaleFactorHisto_ele_Iso_->GetXaxis()->FindBin(eta), dataMCScaleFactorHisto_ele_Iso_->GetYaxis()->FindBin(pt));
    result *= dataMCScaleFactorHisto_ele_ConvHit_->GetBinContent(dataMCScaleFactorHisto_ele_ConvHit_->GetXaxis()->FindBin(eta),dataMCScaleFactorHisto_ele_ConvHit_->GetYaxis()->FindBin(pt));
    //result *=  dataMCScaleFactorHisto_ele_Track_->GetBinContent(dataMCScaleFactorHisto_ele_Track_->GetXaxis()->FindBin(eta),dataMCScaleFactorHisto_ele_Track_->GetYaxis()->FindBin(tempPtTrack)); 
    
    return result;
  }
  
  //~ const double operator()(const  pat::Muon &mu, double pt, double eta, double nVertices){
  const double operator()(const  pat::Muon &mu, double pt, double eta){

    float result = 1.0;
    
    if(pt > 120.) pt= 119.;
    
    result *= dataMCScaleFactorHisto_mu_ID_->GetBinContent(dataMCScaleFactorHisto_mu_ID_->GetXaxis()->FindBin(pt),dataMCScaleFactorHisto_mu_ID_->GetYaxis()->FindBin(fabs(eta)));
    result *= dataMCScaleFactorHisto_mu_Iso_->GetBinContent(dataMCScaleFactorHisto_mu_Iso_->GetXaxis()->FindBin(pt),dataMCScaleFactorHisto_mu_Iso_->GetYaxis()->FindBin(fabs(eta)));
    //result *= dataMCScaleFactorHisto_mu_Impact_->GetBinContent(dataMCScaleFactorHisto_mu_Impact_->GetXaxis()->FindBin(pt),dataMCScaleFactorHisto_mu_Impact_->GetYaxis()->FindBin(fabs(eta)));
    //result *= dataMCScaleFactorHisto_mu_SIP3D_->GetBinContent(dataMCScaleFactorHisto_mu_SIP3D_->GetXaxis()->FindBin(pt),dataMCScaleFactorHisto_mu_SIP3D_->GetYaxis()->FindBin(fabs(eta)));
    //result *= dataMCScaleFactorHisto_mu_Track_->GetBinContent(dataMCScaleFactorHisto_mu_Track_->GetXaxis()->FindBin(eta));
    
    return result;
  }
 
  

private:

  bool useFastSim_;
  
  void getFullSimScaleFactorHistos(std::string dataMCScaleFactorFile_mu_ID, 
              std::string dataMCScaleFactorFile_mu_Iso, 
              //std::string dataMCScaleFactorFile_mu_Impact, 
              //std::string dataMCScaleFactorFile_mu_Track, 
              //std::string dataMCScaleFactorFile_mu_SIP3D, 
              
              std::string dataMCScaleFactorFile_ele, 
              std::string dataMCScaleFactorFile_ele_Reco, 
              //std::string dataMCScaleFactorFile_ele_Track, 
              
              std::string dataMCScaleFactorHisto_mu_ID, 
              std::string dataMCScaleFactorHisto_mu_Iso, 
              //std::string dataMCScaleFactorHisto_mu_Impact, 
              //std::string dataMCScaleFactorHisto_mu_Track, 
              //std::string dataMCScaleFactorHisto_mu_SIP3D, 
              
              std::string dataMCScaleFactorHisto_ele_Reco, 
              std::string dataMCScaleFactorHisto_ele_ID, 
              std::string dataMCScaleFactorHisto_ele_Iso,
              std::string dataMCScaleFactorHisto_ele_ConvHit
              //std::string dataMCScaleFactorHisto_ele_Track
              );

  TH2F * dataMCScaleFactorHisto_mu_ID_;
  TH2F * dataMCScaleFactorHisto_mu_Iso_;
  //TH2F * dataMCScaleFactorHisto_mu_Impact_;
  //TH2F * dataMCScaleFactorHisto_mu_SIP3D_;
  //TH1F * dataMCScaleFactorHisto_mu_Track_;
  
  TH2F * dataMCScaleFactorHisto_ele_Reco_;
  TH2F * dataMCScaleFactorHisto_ele_ID_;
  TH2F * dataMCScaleFactorHisto_ele_Iso_;
  TH2F * dataMCScaleFactorHisto_ele_ConvHit_;
  //TH2F * dataMCScaleFactorHisto_ele_Track_;
};

#endif /* LEPTONFULLSIMSCALEFACTORMAPFUNCTOR_H_ */
