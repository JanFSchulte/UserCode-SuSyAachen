/*
 * TriggerMatchFunctor.h
 *
 *  Created on: 06.05.2015
 *      Author: jschulte
 */

#ifndef TriggerMatchFunctor_H_
#define TriggerMatchFunctor_H_

#include <vector>
#include <iostream>

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"

//#include <DataFormats/PatCandidates/interface/Electron.h>
//#include <DataFormats/PatCandidates/interface/Muon.h>
//#include <DataFormats/PatCandidates/interface/Tau.h>

#include "TLorentzVector.h"

using namespace std;

class TriggerMatchFunctor
{
public:
 TriggerMatchFunctor(edm::ParameterSet const & params):triggerEvent_(){
    hasTrigger_ = false;
    triggerTag_ = params.getParameter<edm::InputTag> ("triggerSrc");
  }
  
  void loadTrigger( const edm::Event& ev){
    edm::Handle<trigger::TriggerEvent> triggerEvent;
    hasTrigger_ = ev.getByLabel(triggerTag_, triggerEvent);
    if(hasTrigger_) triggerEvent_ = *triggerEvent;
  }
  
  bool hasTrigger(){return hasTrigger_;}
  template<class T> std::map<string,int> operator()(const T& lepton);
  
 private:
  bool hasTrigger_;
  trigger::TriggerEvent triggerEvent_; 
  edm::InputTag triggerTag_;  
  
};

template<class T>
std::map<string,int>
TriggerMatchFunctor::operator()(const T& lepton)
{

	 TLorentzVector aVec = TLorentzVector(lepton.px(), lepton.py(), lepton.pz(), lepton.energy());;	
	 std::map<string,int> res;
	 
	 res["matchesSingleElectron"] = 0;
	 res["matchesSingleMuon"] = 0;
	 res["matchesDoubleElectron"] = 0;
	 res["matchesDoubleMuonLeading"] = 0;
	 res["matchesDoubleMuonTrailing"] = 0;
	 res["matchesDoubleMuonLeadingTk"] = 0;
	 res["matchesDoubleMuonTrailingTk"] = 0;
	 res["matchesMuELeading"] = 0;
	 res["matchesMuETrailing"] = 0;
	 res["matchesEMuLeading"] = 0;
	 res["matchesEMuTrailing"] = 0;
	 
	 	 
	 const trigger::TriggerObjectCollection & triggerObjects = triggerEvent_.getObjects();

	 size_t nFilters       = triggerEvent_.sizeFilters();
	 for (size_t iFilter = 0; iFilter < nFilters; ++iFilter) {
		  TString          name = triggerEvent_.filterTag ( iFilter ).label();

	//       if (name.Contains("hltEle27WP80TrackIsoFilter")){
		const trigger::Keys& keys = triggerEvent_.filterKeys( iFilter );
		const trigger::Vids& vids = triggerEvent_.filterIds ( iFilter );
		int nKeys = (int) keys.size();
		int nVids = (int) vids.size();
		assert(nKeys == nVids);
		vector<TLorentzVector> triggerObjectP4s;
		vector<int>            triggerObjectIds;

		for (int iTriggerObject = 0; iTriggerObject < nKeys; ++iTriggerObject ) {
	
				// Get the object ID and key
				int                id  = vids[iTriggerObject];
				trigger::size_type key = keys[iTriggerObject];
	
				// Get the trigger object from the key
				const trigger::TriggerObject & triggerObject = triggerObjects[key];
	
				// Store the trigger object as a TLorentzVector (borrowed from S. Harper)
				TLorentzVector p4;
				p4.SetPxPyPzE(triggerObject.px  (),
					triggerObject.py (),
					triggerObject.pz (),
					triggerObject.energy() );
	
				triggerObjectP4s.push_back ( p4 ) ;
				triggerObjectIds.push_back ( id ) ;
	
		} // end loop over keys/trigger objects passing filters
		
		if ( nKeys > 0 ) {
				for (int iFilterObject = 0; iFilterObject < nKeys; ++iFilterObject) {
					if (!triggerObjectP4s[iFilterObject].Pt() == 0) {
						if (sqrt(pow(aVec.Eta()-triggerObjectP4s[iFilterObject].Eta(),2)+pow(aVec.Phi()-triggerObjectP4s[iFilterObject].Phi(),2)) < 0.2){
							if (name.Contains("hltEle27WP80TrackIsoFilter")){
								res["matchesSingleElectron"] = 1;
							}
							if (name.Contains("hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15")){
								res["matchesSingleMuon"] = 1;
							}
							if (name.Contains("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ")){
								res["matchesDoubleElectron"] = 1;
							}
							if  (name.Contains("hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17") || name.Contains("hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17")){
								res["matchesDoubleMuonLeading"] = 1;
							}
							if  (name.Contains("hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8") || name.Contains("hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8")){
								res["matchesDoubleMuonTrailing"] = 1;
							}
							if          (name.Contains("hltL3fL1sMu10MuOpenL1f0L2f10L3Filtered17")|| name.Contains("hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17")){
								res["matchesDoubleMuonLeadingTk"] = 1;
							}
							if (name.Contains("hltDiMuonGlbFiltered17TrkFiltered8")){
								res["matchesDoubleMuonTrailingTk"] = 1;
							}
							if (name.Contains("hltL1Mu12EG7L3MuFiltered17")){
								res["matchesMuELeading"] = 1;
							}
							if (name.Contains("hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter")){
								res["matchesMuETrailing"] = 1;
							}
							if (name.Contains("hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter")){
								res["matchesEMuLeading"] = 1;
							}
							if (name.Contains("hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8")){
								res["matchesEMuTrailing"] = 1;
							}
						}

					}
				}
			}


		}


	 return res;
}


#endif /* TriggerMatchFunctor_H_ */
