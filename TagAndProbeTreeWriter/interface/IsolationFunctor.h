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
    //~ rhoSrc_( cfg.getParameter<edm::InputTag>( "rhoSource" ) ),
    //~ candidateSrc_( cfg.getParameter<edm::InputTag>( "candSource" ) ),    
    rhoIso_(0.) {}

  const void init(const edm::Event &ev)
  {
    edm::Handle<double> rhoIso_h;
    //~ ev.getByLabel(rhoSrc_, rhoIso_h);
    //~ ev.getByLabel(candidateSrc_, pfCands);        
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


template<class T>
double GetMiniIsolation(const T& lepton, const std::vector<pat::PackedCandidate> &pfCands, const std::string &method, const double &rho)
{
  if (lepton.pt()<5.) return 99999.;

  double r_iso_min = 0.05;
  double r_iso_max = 0.2;
  double kt_scale = 10;
  
  float iso = 0.;

  double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
  

  
  if (abs(lepton.pdgId()) == 13) {
      deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;  
  } 
  else{
  	if (fabs(lepton.eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
  }

  double iso_nh(0.); double iso_ch(0.); 
  double iso_ph(0.); double iso_pu(0.);
  double ptThresh = 0;
  if (abs(lepton.pdgId()) == 13) {
		ptThresh = 0.5;
  }
  double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/lepton.pt()));
  
  if (method == "Puppi"){
	  
	  double val_PuppiWithLep    [3]= {0,0,0} ;
	  double val_PuppiWithoutLep [3]= {0,0,0} ;
	
	  enum particleType{
	      CH = 0, 
	      NH = 1,
	      PH = 2,
	      OTHER = 100000
	  };
	  
	  for (std::vector<pat::PackedCandidate>::const_iterator cand = pfCands.begin(); cand != pfCands.end(); cand++) {
		  if (abs((*cand).pdgId())<7) continue;
		  double dr = deltaR((*cand), lepton);
		  if (dr > r_iso) continue;
		
		  // check particleTyple (CH/NH/PH or other). remove 'other'.
		  const particleType pType =
		  isCH( abs((*cand).pdgId()) ) ? CH :
		  isNH( abs((*cand).pdgId()) ) ? NH :
		  isPH( abs((*cand).pdgId()) ) ? PH : OTHER ;
		  if( pType == OTHER ){
			if( abs((*cand).pdgId()) != 1 && abs((*cand).pdgId()) != 2 
			&& abs( (*cand).pdgId() ) != 11
			&& abs( (*cand).pdgId() ) != 13){
				std::cout <<"candidate with PDGID = " << abs( (*cand).pdgId() ) << " is not CH/NH/PH/e/mu or 1/2 (and this is removed from isolation calculation)"  << std::endl ; 
			}
		  continue ;
		  }

		  // Check particleType dependent DR cut (remove overlapped candiadte)
		  // The threshold values were taken from 'MuonPFIsolationSequence_cff.py'.
		  if( pType == CH && dr < deadcone_ch ) continue ;
		  if( pType == NH && dr < deadcone_nh ) continue ;
		  if( pType == PH && dr < deadcone_ph ) continue ;

		  // The candidate passed all the selection.
		  // Now, add its PT to the variable with weight.
		
		  val_PuppiWithLep   [ pType ] += (*cand).pt() * (*cand).puppiWeight() ;
		  val_PuppiWithoutLep[ pType ] += (*cand).pt() * (*cand).puppiWeightNoLep();
	  }
	  
	  const double reliso_Puppi_withLep    = ( val_PuppiWithLep   [CH] + val_PuppiWithLep   [NH] + val_PuppiWithLep   [PH] ) ;
	  const double reliso_Puppi_withoutlep = ( val_PuppiWithoutLep[CH] + val_PuppiWithoutLep[NH] + val_PuppiWithoutLep[PH] ) ;
	  
	  iso = 0.5 * ( reliso_Puppi_withLep + reliso_Puppi_withoutlep);
	  
  }  
	  
  else {
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
	  
	  iso = iso_ph + iso_nh;
	  
	  //std::cout << iso << " " << rho << " " << GetAEff(lepton) << " " << r_iso << " " << pow(r_iso/0.3,2) << std::endl;
	  if (method == "deltaBeta"){
		  iso -= 0.5*iso_pu;
		  }
	  if (method == "effectiveArea"){
		iso -= GetAEff(lepton) * rho * pow(r_iso/0.3,2);
		}
	  if (iso>0) iso += iso_ch;
	  else iso = iso_ch;
  }

  return iso;
 } 




  const virtual double GetIsolation(const reco::Candidate& lepton, const std::string& method){return -1.;}


// AEffs updated to recommendations on https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF on 12/10/2015

const double GetAEff(const pat::Muon& lepton){

	double etaAbs = fabs(lepton.eta());
	
    double AEff = 0.0735;
    if (etaAbs > 0.8 && etaAbs <= 1.3) AEff =  0.0619;
    if (etaAbs > 1.3 && etaAbs <= 2.0) AEff =  0.0465;
    if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.0433;
    if (etaAbs > 2.2) AEff = 0.0577;
    return AEff;
}


const double GetAEff(const pat::Electron& lepton){

	double etaAbs = fabs(lepton.eta());

    double AEff = 0.1752;
    if (etaAbs > 1.0 && etaAbs <= 1.479) AEff = 0.1862;
    if (etaAbs > 1.479 && etaAbs <= 2.0) AEff = 0.1411;
    if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.1534;
    if (etaAbs > 2.2 && etaAbs <= 2.3) AEff = 0.1903;
    if (etaAbs > 2.3 && etaAbs <= 2.4) AEff = 0.2243;        
    if (etaAbs > 2.4) AEff = 0.2687;
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
