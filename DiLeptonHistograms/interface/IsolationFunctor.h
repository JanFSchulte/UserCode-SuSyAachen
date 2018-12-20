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
public:
  IsolationFunctor(const edm::ParameterSet & cfg, edm::ConsumesCollector && iC  ):
    rhoToken_(iC.consumes<double>(cfg.getParameter<edm::InputTag>("rhoSource"))), 
    candidateToken_(iC.consumes< std::vector<pat::PackedCandidate>  >(cfg.getParameter<edm::InputTag>("candSource"))), 
    rhoIso_(0.) {}

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
  const virtual double GetIsolation(const pat::Electron& lepton, const std::string& method)
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

  const virtual double GetIsolation(const pat::Muon& lepton, const std::string& method)
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

  const virtual double GetIsolation(const pat::PackedCandidate& track, const std::string& method)
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
  
  //double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
  //if (abs(lepton.pdgId()) == 13) {
      //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;  
  //} 
  //else{
    //if (fabs(lepton.eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
  //}
  
  //double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
  //if(lepton.isElectron()) {
  //  if (fabs(lepton.eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
  //} else if(lepton.isMuon()) {
  //  deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;  
  //}
  //
  //
  //double iso_nh(0.); double iso_ch(0.); 
  //double iso_ph(0.); double iso_pu(0.);
  //double ptThresh = 0.5;
  //if (lepton.isElectron()){
  //  ptThresh = 0.0;
  //}
  
  //
  //for (std::vector<pat::PackedCandidate>::const_iterator itPFC = pfCands.begin(); itPFC != pfCands.end(); itPFC++) {
  //  if (abs((*itPFC).pdgId())<7) continue;
  //  double dr = deltaR((*itPFC), lepton);
  //  if (dr > r_iso) continue;
  //    
  //    //////////////////  NEUTRALS  /////////////////////////
  //  if ((*itPFC).charge()==0){
  //    if ((*itPFC).pt()>ptThresh) {
  //      double wpf(1.);
  //        /////////// PHOTONS ////////////
  //      if (abs((*itPFC).pdgId())==22) {
  //        if(dr < deadcone_ph) continue;
  //        iso_ph += wpf*(*itPFC).pt();
  //    /////////// NEUTRAL HADRONS ////////////
  //      } else if (abs((*itPFC).pdgId())==130) {
  //        if(dr < deadcone_nh) continue;
  //        iso_nh += wpf*(*itPFC).pt();
  //      }
  //    }
  //      //////////////////  CHARGED from PV  /////////////////////////
  //  } else if ((*itPFC).fromPV()>1){
  //    if (abs((*itPFC).pdgId())==211) {
  //      if(dr < deadcone_ch) continue;
  //      iso_ch += (*itPFC).pt();
  //    }
  //      //////////////////  CHARGED from PU  /////////////////////////
  //  } else {
  //    if ((*itPFC).pt()>ptThresh){
  //      if(dr < deadcone_pu) continue;
  //      iso_pu += (*itPFC).pt();
  //    }
  //  }
  //}
  //
  //iso = iso_ph + iso_nh;
  //
  ////std::cout << iso << " " << rho << " " << GetAEff(lepton) << " " << r_iso << " " << pow(r_iso/0.3,2) << std::endl;
  //if (method == "deltaBeta"){
  //  iso -= 0.5*iso_pu;
  //  }
  //if (method == "effectiveArea"){
  //iso -= GetAEff(lepton) * rho * pow(r_iso/0.3,2);
  //}
  //if (iso>0) iso += iso_ch;
  //else iso = iso_ch;
  
  iso = lepton.miniPFIsolation().neutralHadronIso() + lepton.miniPFIsolation().photonIso();
  iso -= GetAEff(lepton) * rho * pow(r_iso/0.3,2);

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




  const virtual double GetIsolation(const reco::Candidate& lepton, const std::string& method){return -1.;}


// AEffs updated to recommendations on https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF on 12/10/2015


//# Fall17 version of muon effective areas for dR=0.3,
//# neutral+photon PF isolation component, for usage with mini-isolation
//# |eta| min   |eta| max   effective area
//0.0000         0.8000        0.0566
//0.8000         1.3000        0.0562
//1.3000         2.0000        0.0363
//2.0000         2.2000        0.0119
//2.2000 2.4000 0.0064

const double GetAEff(const pat::Muon& lepton){

  double etaAbs = fabs(lepton.eta());
  
    double AEff = 0.0566;
    if (etaAbs > 0.8 && etaAbs <= 1.3) AEff =  0.0562;
    if (etaAbs > 1.3 && etaAbs <= 2.0) AEff =  0.0363;
    if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.0119;
    if (etaAbs > 2.2) AEff = 0.0064;
    return AEff;
}

//#  The effective areas are based on 90% efficient contours
//#
//# |eta| min   |eta| max   effective area
//0.000        1.000       0.1440
//1.000        1.479       0.1562
//1.479        2.000       0.1032
//2.000        2.200       0.0859
//2.200        2.300       0.1116
//2.300        2.400       0.1321
//2.400 2.500 0.1654

const double GetAEff(const pat::Electron& lepton){

  double etaAbs = fabs(lepton.eta());

    double AEff = 0.1440;
    if (etaAbs > 1.0 && etaAbs <= 1.479) AEff = 0.1562;
    if (etaAbs > 1.479 && etaAbs <= 2.0) AEff = 0.1032;
    if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.0859;
    if (etaAbs > 2.2 && etaAbs <= 2.3) AEff = 0.1116;
    if (etaAbs > 2.3 && etaAbs <= 2.4) AEff = 0.1321;        
    if (etaAbs > 2.4) AEff = 0.1654;
    return AEff;



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

  


  //~ edm::InputTag rhoSrc_;
  //~ edm::InputTag candidateSrc_;  
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT< std::vector<pat::PackedCandidate>  > candidateToken_;
  double rhoIso_;
  edm::Handle< std::vector<pat::PackedCandidate>  > pfCands;  
};


#endif /* ISOLATIONFUNCTOR_H_ */
