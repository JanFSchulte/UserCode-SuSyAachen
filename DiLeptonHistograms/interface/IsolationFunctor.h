/*
 * IsolationFunctor.h
 *
 *  Created on: 20.04.2011
 *      Author: sprenger
 */

#ifndef ISOLATIONFUNCTOR_H_
#define ISOLATIONFUNCTOR_H_

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include "DataFormats/Candidate/interface/Candidate.h"


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"

#include <iostream>

using namespace std;

// base functor class
//====================
class IsolationFunctor
{
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT< std::vector<pat::PackedCandidate>  > candidateToken_;
  
  double rhoIso_;
  edm::Handle< std::vector<pat::PackedCandidate>  > pfCands;  
  std::vector<double> electronEta, muonEta;
  std::vector<double> electronValue, muonValue;
  
public:
  IsolationFunctor(const edm::ParameterSet & cfg, edm::ConsumesCollector && iC  ):
    rhoToken_(iC.consumes<double>(cfg.getParameter<edm::InputTag>("rhoSource"))), 
    candidateToken_(iC.consumes< std::vector<pat::PackedCandidate>  >(cfg.getParameter<edm::InputTag>("candSource"))), 
    rhoIso_(0.) 
    {
    electronEta = cfg.getParameter< std::vector<double>  >("effAreaElectronEta");
    electronValue = cfg.getParameter< std::vector<double>  >("effAreaElectronValue");
    muonEta = cfg.getParameter< std::vector<double>  >("effAreaMuonEta");
    muonValue = cfg.getParameter< std::vector<double>  >("effAreaMuonValue");
  
    }
  const void init(const edm::Event &ev)
  {
    edm::Handle<double> rhoIso_h;     
    ev.getByToken(rhoToken_, rhoIso_h);
    ev.getByToken(candidateToken_, pfCands);        
    rhoIso_ = *(rhoIso_h.product());
    
  }

  template<class T> const double operator()(const T& lepton, const std::string& method)
    { return GetIsolation(lepton,method);}

private:
  const double GetIsolation(const pat::Electron& lepton, const std::string& method)
  {
    double caloIso = 0.;
    caloIso +=  lepton.pfIsolationVariables().sumNeutralHadronEt;
    caloIso += lepton.pfIsolationVariables().sumPhotonEt;
    if (method == "effectiveArea"){
    caloIso -= GetAEff(lepton) * rhoIso_;
  }
  else if (method == "deltaBeta"){
    caloIso -= 0.5* lepton.pfIsolationVariables().sumPUPt;
  }
  


    double iso = 0.;
    iso += lepton.pfIsolationVariables().sumChargedHadronPt;

    if (caloIso > 0)
      iso += caloIso;

  if (method == "miniIsoEA"){
    
    iso = GetMiniIsolation(lepton, *pfCands, "effectiveArea", rhoIso_);
  
  }

  else if (method == "miniIsoDB"){
    
    iso = GetMiniIsolation(lepton, *pfCands, "deltaBeta", rhoIso_);
  
  }
  
  else if (method == "miniIsoPFWeight"){
    
    iso = GetMiniIsolation(lepton, *pfCands, "pfWeight", rhoIso_);
  
  } 
  else if (method == "miniIsoPuppi"){
    
    iso = GetMiniIsolation(lepton, *pfCands, "Puppi", rhoIso_);
  
  } 


    return iso;
  }

  const double GetIsolation(const pat::Muon& lepton, const std::string& method)
  {
    double caloIso = 0.;
    caloIso += lepton.pfIsolationR03().sumNeutralHadronEt;
    caloIso += lepton.pfIsolationR03().sumPhotonEt;
    if (method == "effectiveArea"){
    caloIso -= GetAEff(lepton) * rhoIso_;
  }
  else if (method == "deltaBeta"){
    caloIso -= 0.5* lepton.pfIsolationR03().sumPUPt;
  }

    double iso = 0.;
    iso += lepton.pfIsolationR03().sumChargedHadronPt;
    if(caloIso > 0)
      iso+=caloIso;

  
  if (method == "miniIsoEA"){
    
    iso = GetMiniIsolation(lepton, *pfCands, "effectiveArea", rhoIso_);
  
  }

  else if (method == "miniIsoDB"){
    
    iso = GetMiniIsolation(lepton, *pfCands, "deltaBeta", rhoIso_);
  
  }
  
  else if (method == "miniIsoPFWeight"){
    
    iso = GetMiniIsolation(lepton, *pfCands, "pfWeight", rhoIso_);
  
  }
  else if (method == "miniIsoPuppi"){
    
    iso = GetMiniIsolation(lepton, *pfCands, "Puppi", rhoIso_);
  
  }



    return iso;
  }

  const double GetIsolation(const pat::PackedCandidate& track, const std::string& method)
  {
  double iso = 0.;
  if (method == "trackIso"){
    
    iso = GetTrackIsolation(track, *pfCands);
  
  }


    return iso;
  }


template<class T>
double GetMiniIsolation(const T& lepton, const std::vector<pat::PackedCandidate> &pfCands, const std::string &method, const double &rho)
{
  if (lepton.pt()<5.) return 99999.;

  double r_iso_min = 0.05;
  double r_iso_max = 0.2;
  double kt_scale = 10;
  double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/lepton.pt()));
  
  float iso = 0.;

  
  iso = lepton.miniPFIsolation().neutralHadronIso() + lepton.miniPFIsolation().photonIso();
  
  if (method == "effectiveArea"){
    iso -= GetAEff(lepton) * rho * pow(r_iso/0.3,2);
  }
  if (method == "deltaBeta"){
    iso -= 0.5*lepton.miniPFIsolation().puChargedHadronIso();
  }

  if (iso>0) iso += lepton.miniPFIsolation().chargedHadronIso();
  else iso = lepton.miniPFIsolation().chargedHadronIso();

  return iso;
 } 



const double GetTrackIsolation(const pat::PackedCandidate &track, const std::vector<pat::PackedCandidate> &pfCands)
{
  double r_cone = 0.3;
  
  float iso = 0.;
  
  for (std::vector<pat::PackedCandidate>::const_iterator itPFC = pfCands.begin(); itPFC != pfCands.end(); itPFC++) {
    if ((*itPFC).charge()==0) continue; // skip neutrals  
    
    double dr = deltaR((*itPFC), track);
    
    if (dr > r_cone) continue;
    if (dr == 0. && track.pt()==(*itPFC).pt()) continue;
    
    // Only consider charged pions from PV  
    if ((*itPFC).pt() > 0. && (*itPFC).fromPV()>1){
      if (abs((*itPFC).pdgId())==211) iso += (*itPFC).pt();
    }
  }
    

  return iso;
 } 




  const double GetIsolation(const reco::Candidate& lepton, const std::string& method){return -1.;}


const double GetAEff(const pat::Muon& lepton){

  double etaAbs = fabs(lepton.eta());
  double effArea = 0;
  for (unsigned int i = 0; i < muonEta.size(); i++){
    if (etaAbs > muonEta[i]){
      effArea = muonValue[i];
    }else{
      break;
    }
  }
  return effArea;
}


const double GetAEff(const pat::Electron& lepton){

  double etaAbs = fabs(lepton.eta());
  double effArea = 0;
  for (unsigned int i = 0; i < electronEta.size(); i++){
    if (etaAbs > electronEta[i]){
      effArea = electronValue[i];
    }else{
      break;
    }
  }
  return effArea;



}

bool isNH( long pdgid ){

  const long id = abs( pdgid );

  //     pdgId = cms.vint32(111,130,310,2112),
  if( id == 111 ) return true ; 
  if( id == 130 ) return true ; 
  if( id == 310 ) return true ; 
  if( id == 2112 ) return true ; 

  return false;
}

bool isCH( long pdgid ){

  const long id = abs( pdgid );

  //  pdgId = cms.vint32(211,-211,321,-321,999211,2212,-2212),
  if( id == 211    ) return true ; 
  if( id == 321    ) return true ; 
  if( id == 999211 ) return true ; 
  if( id == 2212   ) return true ; 

  return false;

}

bool isPH( long pdgid ){

  const long id = abs( pdgid );
  if( id == 22 ) return true ; 
  return false;
}

  
};


#endif /* ISOLATIONFUNCTOR_H_ */
