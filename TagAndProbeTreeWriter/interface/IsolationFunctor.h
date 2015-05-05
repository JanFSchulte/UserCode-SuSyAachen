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

#include <iostream>


// base functor class
//====================
class IsolationFunctor
{
public:
  IsolationFunctor(const edm::ParameterSet & cfg):
    rhoSrc_( cfg.getParameter<edm::InputTag>( "rhoSource" ) ),
    rhoIso_(0.) {}

  const void init(const edm::Event &ev )
  {
    edm::Handle<double> rhoIso_h;
    ev.getByLabel(rhoSrc_, rhoIso_h);
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
		caloIso -= GetAEffEle(lepton.eta()) * rhoIso_;
	}
	else if (method == "deltaBeta"){
		caloIso -= 0.5* lepton.pfIsolationVariables().sumPUPt;
	}

    double iso = 0.;
    iso += lepton.pfIsolationVariables().sumChargedHadronPt;

    if (caloIso > 0)
      iso += caloIso;

    return iso;
  }

  const virtual double GetIsolation(const pat::Muon& lepton, const std::string& method)
  {
    double caloIso = 0.;
    caloIso += lepton.pfIsolationR03().sumNeutralHadronEt;
    caloIso += lepton.pfIsolationR03().sumPhotonEt;
    if (method == "effectiveArea"){
		caloIso -= GetAEffMu(lepton.eta()) * rhoIso_;
	}
	else if (method == "deltaBeta"){
		caloIso -= 0.5* lepton.pfIsolationR03().sumPUPt;
	}

    double iso = 0.;
    iso += lepton.pfIsolationR03().sumChargedHadronPt;
    if(caloIso > 0)
      iso+=caloIso;

    return iso;
  }

  const virtual double GetIsolation(const pat::Tau& lepton, const std::string& method)
  {
    double value = lepton.chargedHadronIso()+lepton.photonIso()+lepton.neutralHadronIso();
    return value;
  }

template<class T>
int GetMiniIsolation(const T& lepton, const std::string& method)
{
  
  if (lepton.pt()<5.) return 99999.;

  double r_iso_min = 0.05;
  double r_iso_max = 0.2;
  double kt_scale = 10;

  double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
  if (fabs(e.eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
  else {
      //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
  }

  double iso_nh(0.); double iso_ch(0.); 
  double iso_ph(0.); double iso_pu(0.);
  double ptThresh = 0;
  double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/lepton.pt()));
  for (std::vector<pat::PackedCandidate>::const_iterator itPFC = pfCands.begin(); itPFC != pfCands.end(); itPFC++) {
    if (abs((*itPFC).pdgId())<7) continue;

    double dr = deltaR((*itPFC), lepton);
    if (dr > r_iso) continue;
      
      //////////////////  NEUTRALS  /////////////////////////
    if ((*itPFC).charge()==0){
      if ((*itPFC).pt()>ptThresh) {
        double wpf(1.);
        if (method == "pfWeight"){
          double wpv(1.), wpu(1.);
          for (std::vector<pat::PackedCandidate>::const_iterator itJPFC = pfCands.begin(); itJPFC != pfCands.end(); itJPFC++) {
            double jdr = deltaR2((*itPFC), (*itJPFC));
            if ((*itPFC).charge()!=0 || jdr<0.00001) continue;
            double jpt = (*itJPFC).pt();
            double weight = jpt*jpt/jdr;
            if ((*itJPFC).fromPV()>1 and abs((*itJPFC).charge())>0 and weight>1) wpv *= weight;
            else if (weight>1) wpu *= weight;
          }
          wpv = 0.5*log(wpv);
          wpu = 0.5*log(wpu);
          wpf = wpv/(wpv+wpu);
        }
          /////////// PHOTONS ////////////
        if (abs((*itPFC).pdgId())==22) {
          if(dr < deadcone_ph) continue;
          iso_ph += wpf*(*itPFC).pt();
	    /////////// NEUTRAL HADRONS ////////////
        } else if (abs((*itPFC).pdgId())==130) {
          if(dr < deadcone_nh) continue;
          iso_nh += wpf*(*itPFC).pt();
        }
      }
        //////////////////  CHARGED from PV  /////////////////////////
    } else if ((*itPFC).fromPV()>1){
      if (abs((*itPFC).pdgId())==211) {
        if(dr < deadcone_ch) continue;
        iso_ch += (*itPFC).pt();
      }
        //////////////////  CHARGED from PU  /////////////////////////
    } else {
      if ((*itPFC).pt()>ptThresh){
        if(dr < deadcone_pu) continue;
        iso_pu += (*itPFC).pt();
      }
    }
  }
    
  float iso = 0.;
  iso = iso_ph + iso_nh;
  if (method == "deltaBeta"){
	  iso -= 0.5*iso_pu;
	  }
  if (method == "effectiveArea"){
	iso -= getAEffEle(e.eta()) * rho;
	}
  if (iso>0) iso += iso_ch;
  else iso = iso_ch;
  iso = iso*1./lepton.pt();
  return iso;
 } 



  const virtual double GetIsolation(const reco::Candidate& lepton, const std::string& method){return -1.;}
  
  const double GetAEffEle(double eta)
  {
    //from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/src/EGammaCutBasedEleId.c 
    // but gamma + neutral hadrons values...
    double etaAbs = fabs(eta);
    //~ double AEff = 0.1;
    //~ if (etaAbs > 1.0 && etaAbs <= 1.479) AEff = 0.12;
    //~ if (etaAbs > 1.479 && etaAbs <= 2.0) AEff = 0.085;
    //~ if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.11;
    //~ if (etaAbs > 2.2 && etaAbs <= 2.3) AEff = 0.12;
    //~ if (etaAbs > 2.3 && etaAbs <= 2.4) AEff = 0.12;
    //~ if (etaAbs > 2.4) AEff = 0.13;
    double AEff = 0.1013;
    if (etaAbs > 0.8 && etaAbs <= 1.3) AEff = 0.0988;
    if (etaAbs > 1.3 && etaAbs <= 2.0) AEff = 0.0572;
    if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.0842;
    if (etaAbs > 2.2) AEff = 0.153;
    return AEff;
  }

  const double GetAEffMu(double eta)
  {
    //from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/src/EGammaCutBasedEleId.c 
    // but gamma + neutral hadrons values...
    double etaAbs = fabs(eta);
    double AEff = 0.0913;
    if (etaAbs > 0.8 && etaAbs <= 1.3) AEff = 0.0765;
    if (etaAbs > 1.3 && etaAbs <= 2.0) AEff = 0.0546;
    if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.0728;
    if (etaAbs > 2.2) AEff = 0.1177;
    return AEff;
  }

  edm::InputTag rhoSrc_;
  double rhoIso_;
};


#endif /* ISOLATIONFUNCTOR_H_ */
