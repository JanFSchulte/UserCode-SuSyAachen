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
    
    selected_.clear();
    std::vector<reco::CandidatePtr> leptonPfCands;
    
    // fill pointers for pfCand matching
    for (const auto & lep : *allMuons) {
      for (unsigned int i = 0, n = lep.numberOfSourceCandidatePtrs(); i < n; ++i) { // 0, 1, or 2
          auto ptr = lep.sourceCandidatePtr(i);
          if (ptr.isNonnull()) leptonPfCands.push_back(ptr);
      }
    }
    
    for (const auto & lep : *allElectrons) {
      for (unsigned int i = 0, n = lep.numberOfSourceCandidatePtrs(); i < n; ++i) { // 0, 1, or 2
          auto ptr = lep.sourceCandidatePtr(i);
          if (ptr.isNonnull()) leptonPfCands.push_back(ptr);
      }
    }
    
    // has to be sorted for binary search later
    std::sort(leptonPfCands.begin(), leptonPfCands.end());
    
    
    // add muons to isoTrack collection
    for (const auto & lep : *allMuons) {
      reco::CandidatePtr pfCand(edm::refToPtr(lep.pfCandidateRef()));
      float absIso,absEta,pt,dxy,dz;
      absIso = lep.pfIsolationR03().sumChargedHadronPt;
      absEta = abs(lep.eta());
      pt = lep.pt();
      dxy = lep.dB(pat::Muon::PV2D);
      dz = lep.dB(pat::Muon::PVDZ);
      if (!lep.passed(reco::Muon::CutBasedIdLoose)) continue; // nanoAOD muons have to pass loose cuts
      if (pt < 5 or absEta > 2.4) continue;
      if (abs(dz)  > 0.1) continue;
      if (abs(dxy) > 0.2) continue;
      if (absIso > 5.0) continue;
      if (absIso/pt > 0.2) continue;
      
      int* thePdgId = new int(lep.pdgId()); // currently, only pdgId is stored. have not checked if this memory leaks
      selected_.push_back(thePdgId);
    }
    
    // add electrons to isoTrack collection
    for (const auto & lep : *allElectrons) {
      reco::CandidatePtr pfCand(edm::refToPtr(lep.pfCandidateRef())); 
      if (!std::binary_search(leptonPfCands.begin(), leptonPfCands.end(), pfCand)) continue; // need to match electron to a pfCandidate to match nanoAOD
      float absIso,absEta,pt,dxy,dz;
      absIso = lep.pfIsolationVariables().sumChargedHadronPt;
      absEta = abs(lep.eta());
      pt = lep.pt();
      dxy = lep.dB(pat::Electron::PV2D);
      dz = lep.dB(pat::Electron::PVDZ);
      if (pt < 5 or absEta > 2.4) continue;
      if (abs(dz)  > 0.1) continue;
      if (abs(dxy) > 0.2) continue;
      if (absIso > 5.0) continue;
      if (absIso/pt > 0.2) continue;
      
      int* thePdgId = new int(lep.pdgId()); // currently, only pdgId is stored. have not checked if this memory leaks
      selected_.push_back(thePdgId);
      
    }

    // add others to isoTrack collection
    for(typename collection::const_iterator it = col.product()->begin(); it != col.product()->end(); ++it ){
      if (not (*it).fromPV()) continue;
      if (not (*it).packedCandRef().isNonnull()) continue;
      if (not ((*it).packedCandRef().id() == pfCands.id())) continue;
      reco::CandidatePtr pfCand(edm::refToPtr((*it).packedCandRef()));
      if (std::binary_search(leptonPfCands.begin(), leptonPfCands.end(), pfCand)) continue; // pfCandidates are cleaned from nanoAOD, so they have to be removed
      
      if (abs((*it).pdgId()) == 11 || abs((*it).pdgId()) == 13) continue; // take leptons from the lepton collections

      if ((*it).pt() < 10 or abs((*it).eta()) > 2.4) continue;
      if (abs((*it).dz())  > 0.1) continue;
      if (abs((*it).dxy()) > 0.2) continue;
      if ((*it).pfIsolationDR03().chargedHadronIso() > 5.0) continue;
      if ((*it).pfIsolationDR03().chargedHadronIso()/(*it).pt() > 0.2) continue;
        
      int* thePdgId = new int((*it).pdgId()); // currently, only pdgId is stored. have not checked if this memory leaks
      selected_.push_back(thePdgId);
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
