/*
 * LeptonScaleFactorMapFunctor.h
 *
 *  Created on: 18.10.2019
 *      Author: teroerde
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
            

            
    std::vector<edm::ParameterSet> electronSFNames = params.getParameter< std::vector<edm::ParameterSet> >("electronScaleFactors");
    std::vector<edm::ParameterSet> muonSFNames = params.getParameter< std::vector<edm::ParameterSet> >("muonScaleFactors");
    
    electronPtThreshold = params.getParameter< double >("electronPtThreshold");
    muonPtThreshold = params.getParameter< double >("muonPtThreshold");
    
    TH2F* hist;
    for (auto sf : electronSFNames){
      std::string fileName = sf.getParameter< std::string >("fileName");
      TFile file(fileName.c_str());
      std::string histName = sf.getParameter< std::string >("histName");
      int order = sf.getParameter< int >("order");
      hist = (TH2F*)file.Get(histName.c_str());
      electronScaleFactorHists.push_back(hist);
      electronOrder.push_back(order);
      
    }
    for (auto sf : muonSFNames){
      std::string fileName = sf.getParameter< std::string >("fileName");
      TFile file(fileName.c_str());
      std::string histName = sf.getParameter< std::string >("histName");
      int order = sf.getParameter< int >("order");
      hist = (TH2F*)file.Get(histName.c_str());
      muonScaleFactorHists.push_back(hist);
      muonOrder.push_back(order);
      
    }
    hist = nullptr;
            
    
  }
  
  
  
  //~ const double operator()(const  pat::Electron &ele, double pt, double absEta, double nVertices){
  const std::pair<double, double> operator()(const  pat::Electron &ele, double pt, double eta){

    float result = 1.0;
    float error = 0.0;
    eta = ele.superCluster()->eta();
    if(pt > electronPtThreshold ) pt= electronPtThreshold-1.0;
  
    float firstVar = 0;
    float secondVar = 0;
  
    for(unsigned int i=0; i < electronScaleFactorHists.size(); i++){
      if (electronOrder[i] == 1){
        firstVar = eta;
        secondVar = pt;
      }else if (electronOrder[i] == 11){
        firstVar = fabs(eta);
        secondVar = pt;
      }else if (electronOrder[i] == 0){
        firstVar = pt;
        secondVar = eta;
      }else if (electronOrder[i] == 10){
        firstVar = pt;
        secondVar = fabs(eta);
      }
      result *= electronScaleFactorHists[i]->GetBinContent(electronScaleFactorHists[i]->GetXaxis()->FindBin(firstVar), electronScaleFactorHists[i]->GetYaxis()->FindBin(secondVar) );
      error += std::pow(electronScaleFactorHists[i]->GetBinError(electronScaleFactorHists[i]->GetXaxis()->FindBin(firstVar), electronScaleFactorHists[i]->GetYaxis()->FindBin(secondVar) ), 2);
      
    }
    //error += std::pow(0.01, 2);
    error = std::sqrt(error);
    

    return std::pair<double, double>(result, error);
  }
  
  //~ const double operator()(const  pat::Muon &mu, double pt, double absEta, double nVertices){
  const std::pair<double, double> operator()(const  pat::Muon &mu, double pt, double eta){

    float result = 1.0;
    float error = 0.0;
    if(pt > muonPtThreshold) {
      pt= muonPtThreshold-1.0;
    }

    float firstVar  = 0;
    float secondVar = 0;
    
    for(unsigned int i=0; i < muonScaleFactorHists.size(); i++){
      if (muonOrder[i] == 1){
        firstVar = eta;
        secondVar = pt;
      }else if (muonOrder[i] == 11){
        firstVar = fabs(eta);
        secondVar = pt;
      }else if (muonOrder[i] == 0){
        firstVar = pt;
        secondVar = eta;
      }else if (muonOrder[i] == 10){
        firstVar = pt;
        secondVar = fabs(eta);
      }
      result *= muonScaleFactorHists[i]->GetBinContent(muonScaleFactorHists[i]->GetXaxis()->FindBin(firstVar), muonScaleFactorHists[i]->GetYaxis()->FindBin(fabs(secondVar)) );
    }
    error = 0.03;
    
    
    
    return std::pair<double, double>(result, error);
  }

  

private:

  bool useFastSim_;

  double electronPtThreshold;
  double muonPtThreshold;
  
  std::vector<TH2F*> electronScaleFactorHists;
  std::vector<int> electronOrder;
  std::vector<TH2F*> muonScaleFactorHists;
  std::vector<int> muonOrder;
};

#endif /* LEPTONFASTSIMSCALEFACTORMAPFUNCTOR_H_ */
