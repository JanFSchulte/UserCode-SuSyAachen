/*
 * TriggerMatchFunctorMiniAOD.h
 *
 *  Created on: 06.05.2015
 *      Author: jschulte
 */

#ifndef TriggerMatchFunctorMiniAOD_H_
#define TriggerMatchFunctorMiniAOD_H_

#include <vector>
#include <iostream>

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
//#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
//#include <DataFormats/PatCandidates/interface/Electron.h>
//#include <DataFormats/PatCandidates/interface/Muon.h>
//#include <DataFormats/PatCandidates/interface/Tau.h>

#include "TLorentzVector.h"

using namespace std;

class TriggerMatchFunctorMiniAOD
{
public:
 TriggerMatchFunctorMiniAOD(edm::ParameterSet const & params):triggerObjects_(){
    hasTrigger_ = false;
    bitsTag_ = params.getParameter<edm::InputTag> ("bits");
    //prescaleTag_ = params.getParameter<edm::InputTag> ("prescales");  
    objectTag_ = params.getParameter<edm::InputTag> ("objects");       
  }
  
  void loadTrigger( const edm::Event& ev){
    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    //edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    hasTrigger_ = ev.getByLabel(bitsTag_, triggerBits);
    //hasTrigger_ = ev.getByLabel(prescaleTag_, triggerPrescales);
    hasTrigger_ = ev.getByLabel(objectTag_, triggerObjects);        
    if(hasTrigger_) triggerObjects_ = triggerObjects;
    if(hasTrigger_) triggerBits_ = triggerBits;  
    if(hasTrigger_) names_ = ev.triggerNames(*triggerBits);
      
  }
  
  bool hasTrigger(){return hasTrigger_;}
  template<class T> std::map<string,int> operator()(const T& lepton);
  
 private:
  bool hasTrigger_;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::Handle<edm::TriggerResults> triggerBits_;
  edm::InputTag bitsTag_; 
  //edm::InputTag prescaleTag_;
  edm::InputTag objectTag_;
  edm::TriggerNames names_;     
  
};

template<class T>
std::map<string,int>
TriggerMatchFunctorMiniAOD::operator()(const T& lepton)
{

	 TLorentzVector aVec = TLorentzVector(lepton.px(), lepton.py(), lepton.pz(), lepton.energy());;	
	 std::map<string,int> res;
	 
	 res["matchesSingleElectron"] = 0;
	 res["matchesSingleMuon"] = 0;
	 res["matchesDoubleElectronLeading"] = 0;
	 res["matchesDoubleElectronTrailing"] = 0;	 
	 res["matchesDoubleMuonLeading"] = 0;
	 res["matchesDoubleMuonTrailing"] = 0;
	 res["matchesDoubleMuonLeadingTk"] = 0;
	 res["matchesDoubleMuonTrailingTk"] = 0;
	 res["matchesDoubleMuonLeadingBoth"] = 0;
	 res["matchesDoubleMuonTrailingBoth"] = 0;	 
	 res["matchesMuELeading"] = 0;
	 res["matchesMuETrailing"] = 0;
	 res["matchesEMuLeading"] = 0;
	 res["matchesEMuTrailing"] = 0;
	 
	 
	 bool passesTrailingMuonPt = false;
	 bool passesMuonIso = false;
	 bool passesLeadingMuonPt = false;

	 bool passesTrailingTkMuonPt = false;
	 bool passesMuonTkIso = false;
	 bool passesLeadingTkMuonPt = false;
	 	 
    for (pat::TriggerObjectStandAlone obj : *triggerObjects_) { // note: not "const &" since we want to call unpackPathNames

        obj.unpackPathNames(names_);
	 	size_t nFilters       = obj.filterLabels().size();
	 	for (size_t iFilter = 0; iFilter < nFilters; ++iFilter) {
		  TString          name = obj.filterLabels()[iFilter];

			if (!obj.pt() == 0) {
				if (sqrt(pow(aVec.Eta()-obj.eta(),2)+pow(aVec.Phi()-obj.phi(),2)) < 0.2){

					// corresponds to menu for Phys14 samples. Update when appropriate		
					if (name.Contains("hltEle27WP85GsfTrackIsoFilter")){
						res["matchesSingleElectron"] = 1;
					}
					if (name.Contains("hltL3crIsoL1sMu20Eta2p1L1f0L2f20QL3f24QL3crIsoRhoFiltered0p15IterTrk02")){
						res["matchesSingleMuon"] = 1;
					}
					if (name.Contains("hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg1Filter")){
						res["matchesDoubleElectronLeading"] = 1;
					}
					if (name.Contains("hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg2Filter")){
						res["matchesDoubleElectronTrailing"] = 1;						
					}
					
					if (name.Contains("hltDiMuonGlb17Trk8DzFiltered0p2RelTrkIsoFiltered0p4")){
						passesMuonTkIso = true;
					}
					if (name.Contains("hltDiMuonGlb17Glb8DzFiltered0p2RelTrkIsoFiltered0p4")){
						passesMuonIso = true;
					}
					if (name.Contains("hltDiMuonGlbFiltered17TrkFiltered8")){
						passesTrailingTkMuonPt = true;
					}
					if (name.Contains("hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8")){
						passesTrailingMuonPt = true;
					}
					if (name.Contains("hltL3fL1sDoubleMu103p5L1f0L2f10L3Filtered17")){
						passesLeadingTkMuonPt = true;
					}
					if (name.Contains("hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17")){
						passesLeadingMuonPt = true;
					}															

					
										
					if (name.Contains("hltL1Mu12EG7L3IsoMuFiltered23")){
						res["matchesMuELeading"] = 1;
					}
					if (name.Contains("hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter")){
						res["matchesMuETrailing"] = 1;
					}
					if (name.Contains("hltL1sL1Mu5EG20ORL1Mu5IsoEG18L3IsoFiltered8")){
						res["matchesEMuLeading"] = 1;
					}
					if (name.Contains("hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter")){
						res["matchesEMuTrailing"] = 1;
					}

				}
			}


		}

		if  ((passesLeadingTkMuonPt || passesLeadingMuonPt) && (passesMuonIso || passesMuonTkIso)){
			res["matchesDoubleMuonLeadingBoth"] = 1;
		}
		if  ((passesTrailingTkMuonPt || passesTrailingMuonPt) && (passesMuonIso || passesMuonTkIso)){
			res["matchesDoubleMuonTrailingBoth"] = 1;
		}
		if  (passesLeadingTkMuonPt && passesMuonTkIso){
			res["matchesDoubleMuonLeadingTk"] = 1;
		}
		if  (passesTrailingTkMuonPt && passesMuonTkIso){
			res["matchesDoubleMuonTrailingTk"] = 1;
		}
		if  (passesLeadingMuonPt && passesMuonIso){
			res["matchesDoubleMuonLeading"] = 1;
		}
		if  (passesTrailingMuonPt && passesMuonIso){
			res["matchesDoubleMuonTrailing"] = 1;
		}		

    }
    return res;
}


#endif /* TriggerMatchFunctorMiniAOD_H_ */
