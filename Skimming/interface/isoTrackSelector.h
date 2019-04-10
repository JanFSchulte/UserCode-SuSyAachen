#ifndef isoTrackSelector_h
#define isoTrackSelector_h

// Data formats

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include <DataFormats/PatCandidates/interface/IsolatedTrack.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/RefToBase.h"

//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"

//STL
#include <vector>

bool closeToLepton(const std::vector<const reco::LeafCandidate*>& leptons, const reco::LeafCandidate& track){
  for (const auto& lep: leptons){
    double dPhi = (*lep).phi() - track.phi();
    double dEta = (*lep).eta() - track.eta();
    double dR2 = dPhi*dPhi + dEta*dEta;
    if (dR2 < 0.2*0.2){
      return true;
    }
  }
  return false;
}


template<typename T, typename collectionType, typename containerType>
struct isoTrackSelector {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  isoTrackSelector ( const edm::ParameterSet & cfg, edm::ConsumesCollector iC):
    allElectronToken_(iC.consumes<std::vector< pat::Electron >>(cfg.getParameter<edm::InputTag>("electronSource"))),
    allMuonToken_   (iC.consumes<std::vector< pat::Muon >>(cfg.getParameter<edm::InputTag>("muonSource"))),
    pfCandToken_(iC.consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("pfCandSource"))){ }
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {
    
    edm::Handle< std::vector< pat::Electron > > allElectrons;
    ev.getByToken(allElectronToken_, allElectrons);
    
    edm::Handle< std::vector< pat::Muon > > allMuons;
    ev.getByToken(allMuonToken_, allMuons);
    
    edm::Handle<pat::PackedCandidateCollection> pfCands;
    ev.getByToken(pfCandToken_, pfCands); 
    
    // remove tracks from last event
    selected_.clear();
    
    std::vector<const reco::LeafCandidate*> selectedLeps;
    
    
    // add muons to isoTrack collection
    for (const auto & lep : *allMuons) {
      if (!lep.passed(reco::Muon::CutBasedIdLoose)) continue; // nanoAOD muons have to pass loose cuts
      if (lep.pt() > 3 && lep.track().isNonnull()){
        selectedLeps.push_back( &(lep) );
      }
      if (lep.pt() < 5 or abs(lep.eta()) > 2.4) continue;
      if (abs(lep.dB(pat::Muon::PVDZ))  > 0.1) continue;
      if (abs(lep.dB(pat::Muon::PV2D)) > 0.2) continue;
      if (lep.pfIsolationR03().sumChargedHadronPt > 5.0) continue;
      if (lep.pfIsolationR03().sumChargedHadronPt/lep.pt() > 0.2) continue;
      
      
      selected_.push_back( &(lep) );
    }
    
    // add electrons to isoTrack collection
    for (const auto & lep : *allElectrons) {
      if (lep.pt() > 5){
        selectedLeps.push_back( &(lep) );
      }
      if (not lep.pfCandidateRef().isNonnull()) continue;
      if (lep.pt() < 5 or abs(lep.eta()) > 2.4) continue;
      if (abs(lep.dB(pat::Electron::PVDZ))  > 0.1) continue;
      if (abs(lep.dB(pat::Electron::PV2D)) > 0.2) continue;
      if (lep.pfIsolationVariables().sumChargedHadronPt > 5.0) continue;
      if (lep.pfIsolationVariables().sumChargedHadronPt/lep.pt() > 0.2) continue;

      
      selected_.push_back( &(lep) );
    }

    // add others to isoTrack collection
    for(typename collection::const_iterator it = col.product()->begin(); it != col.product()->end(); ++it ){
      if (not (*it).fromPV()) continue;
      if (not (*it).packedCandRef().isNonnull()) continue;
      if (not ((*it).packedCandRef().id() == pfCands.id())) continue;
      if (closeToLepton(selectedLeps, (*it))) continue; // remove tracks matched to closeby leptons to emulate nanoAOD cleaning
      if (abs((*it).pdgId()) == 11 || abs((*it).pdgId()) == 13) continue; // take leptons from the lepton collections, ignore in this collection
      if ((*it).pt() < 10 or abs((*it).eta()) > 2.4) continue;
      if (abs((*it).dz())  > 0.1) continue;
      if (abs((*it).dxy()) > 0.2) continue;
      if ((*it).pfIsolationDR03().chargedHadronIso() > 5.0) continue;
      if ((*it).pfIsolationDR03().chargedHadronIso()/(*it).pt() > 0.2) continue;
      selected_.push_back( &(*it) );
    }
  }

  size_t size() const { return selected_.size(); }
private:
  container selected_;
  edm::EDGetTokenT< std::vector< pat::Electron > > allElectronToken_;
  edm::EDGetTokenT< std::vector< pat::Muon > >  allMuonToken_;
  edm::EDGetTokenT< pat::PackedCandidateCollection > pfCandToken_;
};




#endif
