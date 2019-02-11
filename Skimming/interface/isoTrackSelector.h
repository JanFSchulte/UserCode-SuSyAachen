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

//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"

//STL
#include <vector>

bool matchTrackToLepton(const std::vector<pat::Electron> &lepColl,const pat::IsolatedTrack &track, pat::Electron &matched){
  float closestDR=10;
  bool isMatched = false;
  for (auto const &lep: lepColl){
    float dEta = lep.p4().eta()-track.eta();
    float dPhi = lep.p4().phi()-track.phi();
    float dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
    float ptAbsDiff = abs(track.pt()-lep.pt())/track.pt();
    if (dR < 0.1 && dR < closestDR && ptAbsDiff < 0.1){
      matched = lep;
      closestDR = dR;
      isMatched = true;
    }
  }
  return isMatched;
}

bool matchTrackToLepton(const std::vector<pat::Muon> &lepColl, const pat::IsolatedTrack &track, pat::Muon &matched){
  float closestDR=10;
  bool isMatched = false;
  for (auto const &lep: lepColl){
    float dEta = lep.p4().eta()-track.eta();
    float dPhi = lep.p4().phi()-track.phi();
    float dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
    float ptAbsDiff = abs(track.pt()-lep.pt())/track.pt();
    if (dR < 0.1 && dR < closestDR && ptAbsDiff < 0.1){
      matched = lep;
      closestDR = dR;
      isMatched = true;
    }
  }
  return isMatched;
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
    
    std::vector<reco::CandidatePtr> leptonPfCands;
    
    for (const auto & lep : *allMuons) {
      for (unsigned int i = 0, n = lep.numberOfSourceCandidatePtrs(); i < n; ++i) {
          auto ptr = lep.sourceCandidatePtr(i);
          if (ptr.isNonnull()) leptonPfCands.push_back(ptr);
      }
    }
    
    for (const auto & lep : *allElectrons) {
      for (unsigned int i = 0, n = lep.numberOfSourceCandidatePtrs(); i < n; ++i) {
          auto ptr = lep.sourceCandidatePtr(i);
          if (ptr.isNonnull()) leptonPfCands.push_back(ptr);
      }
    }
    std::sort(leptonPfCands.begin(), leptonPfCands.end());

    selected_.clear();
    for(typename collection::const_iterator it = col.product()->begin(); it != col.product()->end(); ++it ){
        if (not (*it).fromPV()) continue;
        if (not (*it).packedCandRef().isNonnull()) continue;
        if (not ((*it).packedCandRef().id() == pfCands.id())) continue;
        reco::CandidatePtr pfCand(edm::refToPtr((*it).packedCandRef()));

        
        if (abs((*it).pdgId()) == 11 || abs((*it).pdgId()) == 13){
          float absIso,absEta,pt,dxy,dz;
          if (abs((*it).pdgId()) == 11){
            pat::Electron lep;
            matchTrackToLepton(*allElectrons, (*it), lep);
            if (lep.pt() > 1 && std::binary_search(leptonPfCands.begin(), leptonPfCands.end(), pfCand)){
              absIso = lep.pfIsolationVariables().sumChargedHadronPt;
              absEta = abs(lep.eta());
              pt = lep.pt();
              dxy = lep.dB(pat::Electron::PV2D);
              dz = lep.dB(pat::Electron::PVDZ);
            }else{
              absIso = (*it).pfIsolationDR03().chargedHadronIso();
              absEta = abs((*it).eta());
              pt = (*it).pt();
              dxy = (*it).dxy();
              dz = (*it).dz();
            }
          }else{
            pat::Muon lep;
            matchTrackToLepton(*allMuons, (*it), lep);
            if (lep.pt() > 1 && std::binary_search(leptonPfCands.begin(), leptonPfCands.end(), pfCand)){
              absIso = lep.pfIsolationR03().sumChargedHadronPt;
              absEta = abs(lep.eta());
              pt = lep.pt();
              dxy = lep.dB(pat::Muon::PV2D);
              dz = lep.dB(pat::Muon::PVDZ);
            }else{
              absIso = (*it).pfIsolationDR03().chargedHadronIso();
              absEta = abs((*it).eta());
              pt = (*it).pt();
              dxy = (*it).dxy();
              dz = (*it).dz();
            }
          }
          if (pt < 5 or absEta > 2.4) continue;
          if (abs(dz)  > 0.1) continue;
          if (abs(dxy) > 0.2) continue;
          if (absIso > 5.0) continue;
          if (absIso/pt > 0.2) continue;
          //std::cout << ev.id().event() << " pdgId " << (*it).pdgId() << " pt " << pt << " eta " << absEta << " ch_iso " << absIso << " relIso " << absIso/pt << std::endl;
        }else{
          if ((*it).pt() < 10 or abs((*it).eta()) > 2.4) continue;
          if (abs((*it).dz())  > 0.1) continue;
          if (abs((*it).dxy()) > 0.2) continue;
          if ((*it).pfIsolationDR03().chargedHadronIso() > 5.0) continue;
          if ((*it).pfIsolationDR03().chargedHadronIso()/(*it).pt() > 0.2) continue;
          // to reproduce cleaning in nanoAOD
          if (std::binary_search(leptonPfCands.begin(), leptonPfCands.end(), pfCand)) continue;
  
          //std::cout << ev.id().event() << " pdgId " << (*it).pdgId() << " pt " << (*it).pt() << " eta " << (*it).eta() << " phi " << (*it).phi() << std::endl;
          //std::cout << "absIso " << (*it).pfIsolationDR03().chargedHadronIso() << " relIso " << (*it).pfIsolationDR03().chargedHadronIso()/(*it).pt() << std::endl;
        }
        
            selected_.push_back( & (*it) );

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
