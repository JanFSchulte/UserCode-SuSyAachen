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
