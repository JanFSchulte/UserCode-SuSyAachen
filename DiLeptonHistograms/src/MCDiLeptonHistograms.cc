// -*- C++ -*-
//
// Package:    Histograms
// Class:      MCDiLeptonHistograms
// 
/**\class MCDiLeptonHistograms MCDiLeptonHistograms.cc brot/MCDiLeptonHistograms/src/MCDiLeptonHistograms.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  matthias edelhoff
//         Created:  Tue Oct 27 13:50:40 CET 2009
// $Id: MCDiLeptonHistograms.cc,v 1.2 2010/04/01 09:36:57 sprenger Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/JetReco/interface/GenJet.h>

//ROOT
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

//
// class decleration
//

class MCDiLeptonHistograms : public edm::EDAnalyzer {
public:
  explicit MCDiLeptonHistograms(const edm::ParameterSet&);
  ~MCDiLeptonHistograms();

private:
  typedef edm::View<reco::Candidate> collection;

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  void initInvMass(TFileDirectory &file);
  void fillInvMass( const reco::Candidate &mother, const collection &electrons,\
		    const collection &muons, const collection &taus, double &weight);

  void initCommon(const std::string &name, TFileDirectory &file);
  template < class colT> void fillCommon(const std::string &name,  const colT &particles, double &weight);

  bool matches( const reco::Candidate& aCand, const reco::Candidate& bCand);
  bool contains ( collection col, const reco::Candidate& a);

  edm::InputTag genParticleTag_;

  edm::InputTag electronTag_;
  edm::InputTag muonTag_;
  edm::InputTag tauTag_;

  edm::InputTag tauJetTag_;

  //histos
  std::map<std::string, TH1F* > invMassHistos_; 
  std::map<std::string, TH1F* > commonHistos_; 
  std::map<std::string, TH2F* > commonHistos2D_; 


};

// constructors and destructor
MCDiLeptonHistograms::MCDiLeptonHistograms(const edm::ParameterSet& iConfig)
{
  // read config
  genParticleTag_ = iConfig.getParameter<edm::InputTag>("genParticles");

  electronTag_ = iConfig.getParameter<edm::InputTag>("electrons");
  muonTag_ = iConfig.getParameter<edm::InputTag>("muons");
  tauTag_ = iConfig.getParameter<edm::InputTag>("taus");

  tauJetTag_ = iConfig.getParameter<edm::InputTag>("tauJets");

  // init histos
  edm::Service<TFileDirectory> file;
  initInvMass( *file);
  initCommon( "Electron", *file);
  initCommon( "Muon", *file);
  initCommon( "Tau", *file);
  initCommon( "TauJets", *file);
}

MCDiLeptonHistograms::~MCDiLeptonHistograms()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// member functions
// ------------ method called to for each event  ------------
void
MCDiLeptonHistograms::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< collection > genParticles;
  iEvent.getByLabel(genParticleTag_, genParticles);

  edm::Handle<collection> electrons;
  iEvent.getByLabel(electronTag_, electrons);

  edm::Handle<collection> muons;
  iEvent.getByLabel(muonTag_, muons);

  edm::Handle<collection> taus;
  iEvent.getByLabel(tauTag_, taus);

  edm::Handle<std::vector<reco::GenJet> > tauJets;
  iEvent.getByLabel(tauJetTag_, tauJets);


  double weight = 1.;// TODO get weight...
  //int nMu = 0;
  //    std::cout << "++++++"<< std::endl;
  for( collection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it){
    //if(contains(*muons, *it))  nMu++; 

    fillInvMass(*it, *electrons, *muons, *taus, weight); 
  }
  //  if( nMu != 2) std::cout << "-------! "<<nMu<<std::endl;

  fillCommon<collection>("Electron", *electrons, weight);
  fillCommon<collection>("Muon", *muons, weight);
  fillCommon<collection>("Tau", *taus, weight);
  fillCommon< std::vector<reco::GenJet> >("TauJets", *tauJets, weight);
  
  //for(std::vector<edm::WorkerSummary>::const_iterator iSummary = triggerReport->workerSummaries.begin();
  //    iSummary != triggerReport->workerSummaries.end(); ++iSummary){
  //  std::cout << (*iSummary).moduleLabel <<" "<<(*iSummary).timesPassed <<", " ;
  //}
  
}

void 
MCDiLeptonHistograms::initCommon(const std::string &name, TFileDirectory &file)
{
  TFileDirectory dir = file.mkdir(("Common "+name).c_str());
  commonHistos_[name+"Count"] = dir.make<TH1F>( (name +" count").c_str(), (name + " count").c_str(), 100, 0.0, 100.0);
  commonHistos_[name+"Pt"] = dir.make<TH1F>( (name +" pt").c_str(), (name + " pt").c_str(), 1000, 0.0, 1000.0);
  commonHistos_[name+"Eta"] = dir.make<TH1F>( (name + " eta").c_str(), (name +" eta").c_str(), 250, -2.5, 2.5);
  commonHistos_[name+"Phi"] = dir.make<TH1F>( (name + " phi").c_str(), (name +" phi").c_str(), 250, -3.2, 3.2);
  commonHistos_[name+"Charge"] = dir.make<TH1F>( (name +" charge").c_str(), (name + " charge").c_str(), 16, -7.5, 8.5);

  commonHistos2D_[name+"EtaVsPt"] = dir.make<TH2F>( (name + "Eta-Pt").c_str(), (name + "Eta - Pt").c_str(), 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
}

template < class colT> void 
MCDiLeptonHistograms::fillCommon(const std::string &name,  const colT &particles, double &weight)
{
  commonHistos_[name+"Count"]->Fill(particles.size());
  for (typename colT::const_iterator it = particles.begin(); it != particles.end(); ++it){

    commonHistos_[name+"Pt"]->Fill( (*it).pt() );
    commonHistos_[name+"Eta"]->Fill( (*it).eta() );
    commonHistos_[name+"Phi"]->Fill( (*it).phi() );
    commonHistos_[name+"Charge"]->Fill( (*it).charge() );
    
    commonHistos2D_[name+"EtaVsPt"]->Fill( (*it).pt(), (*it).eta() );
  }
  
}


void 
MCDiLeptonHistograms::initInvMass(TFileDirectory &file)
{
  TFileDirectory InvMass = file.mkdir("Invariant Mass");
  invMassHistos_["FromNeutralinoInvMass"] = InvMass.make<TH1F>( "Invariant mass of OS signal decays", "Invariant mass of signal decays", 300, 0, 300);
  invMassHistos_["FromNeutralinoInvMassE"] = InvMass.make<TH1F>( "Invariant mass of OS signal decays into Electrons", "Invariant mass of signal decays into Electrons", 300, 0, 300);
  invMassHistos_["FromNeutralinoInvMassMu"] = InvMass.make<TH1F>( "Invariant mass of OS signal decays into Muons", "Invariant mass of signal decays into Muons", 300, 0, 300);
  invMassHistos_["FromNeutralinoInvMassTau"] = InvMass.make<TH1F>( "Invariant mass of OS signal decays into Taus", "Invariant mass of signal decays Taus", 300, 0, 300);
  
  invMassHistos_["FromZInvMass"] = InvMass.make<TH1F>( "Invariant mass of Z decays", "Invariant mass of Z decays", 300, 0, 300);
  invMassHistos_["FromZInvMassE"] = InvMass.make<TH1F>( "Invariant mass of Z decays into Electrons", "Invariant mass of Z decays into Electrons", 300, 0, 300);
  invMassHistos_["FromZInvMassMu"] = InvMass.make<TH1F>( "Invariant mass of Z decays into Muons", "Invariant mass of Z decays into Muons", 300, 0, 300);  
  invMassHistos_["FromZInvMassTau"] = InvMass.make<TH1F>( "Invariant mass of Z decays into Taus", "Invariant mass of Z decays into Taus", 300, 0, 300);
}

void 
MCDiLeptonHistograms::fillInvMass( const reco::Candidate &mother, const  collection &electrons, \
				   const collection &muons, const collection &taus, double &weight)
{
  std::vector<const reco::Candidate *> chi_1;
  std::vector<const reco::Candidate *> selectedElectrons;
  std::vector<const reco::Candidate *> selectedMuons;
  std::vector<const reco::Candidate *> selectedTaus;
  for(size_t j = 0; j < mother.numberOfDaughters(); ++j){
    const reco::Candidate *d = mother.daughter(j);
    //std::cout << d->pdgId() << std::endl;
    //if( abs(d->pdgId()) == 1000022){chi_1.push_back(d);}
    if( abs(d->pdgId()) == 11 && contains( electrons, *d)){selectedElectrons.push_back(d); }
    if( abs(d->pdgId()) == 13 && contains( muons, *d)){ selectedMuons.push_back(d); }
    if( abs(d->pdgId()) == 15 && contains( taus, *d)){selectedTaus.push_back(d);}
  }
    
  // OS signal decay
  if(abs(mother.pdgId()) == 1000023){
    if (chi_1.size()==1 && selectedElectrons.size()==2){
      double m = (selectedElectrons[0]->p4()+selectedElectrons[1]->p4()).M();
      invMassHistos_["FromNeutralinoInvMass"]->Fill(m,weight);
      invMassHistos_["FromNeutralinoInvMassE"]->Fill(m,weight);
     }
    if (chi_1.size()==1 && selectedMuons.size()==2){
      double m = (selectedMuons[0]->p4()+selectedMuons[1]->p4()).M();
      invMassHistos_["FromNeutralinoInvMass"]->Fill(m,weight);
      invMassHistos_["FromNeutralinoInvMassMu"]->Fill(m,weight);
    }
    if (chi_1.size()==1 && selectedTaus.size()==2){
      double m = (selectedTaus[0]->p4()+selectedTaus[1]->p4()).M();
      invMassHistos_["FromNeutralinoInvMass"]->Fill(m,weight);
      invMassHistos_["FromNeutralinoInvMassTau"]->Fill(m,weight);
    }
  }
  if ( abs(mother.pdgId())== 22 || abs(mother.pdgId()) == 23){
    if (selectedElectrons.size()==2 && selectedElectrons[0]->charge()*selectedElectrons[1]->charge()<0){
      double m = (selectedElectrons[0]->p4()+selectedElectrons[1]->p4()).M();
      invMassHistos_["FromZInvMass"]->Fill(m,weight);
      invMassHistos_["FromZInvMassE"]->Fill(m,weight);
      //std::cout << selectedElectrons[0].pdgId() << std::endl;
      //std::cout << selectedElectrons[1].pdgId() << std::endl;
    }
    if (selectedMuons.size()==2 && selectedMuons[0]->charge()*selectedMuons[1]->charge()<0){
      double m = (selectedMuons[0]->p4()+selectedMuons[1]->p4()).M();
      invMassHistos_["FromZInvMass"]->Fill(m,weight);
      invMassHistos_["FromZInvMassMu"]->Fill(m,weight);
      //std::cout << selectedMuons[0].pdgId() << std::endl;
      //std::cout << selectedMuons[1].pdgId() << std::endl;
    }
    if (selectedTaus.size()==2 && selectedTaus[0]->charge()*selectedTaus[1]->charge()<0){
      double m = (selectedTaus[0]->p4()+selectedTaus[1]->p4()).M();
      invMassHistos_["FromZInvMass"]->Fill(m,weight);
      invMassHistos_["FromZInvMassTau"]->Fill(m,weight);
      //std::cout << selectedTaus[0].pdgId() << std::endl;
      //std::cout << selectedTaus[1].pdgId() << std::endl;
    }
  }
}

bool 
MCDiLeptonHistograms::matches( const reco::Candidate& aCand, const reco::Candidate& bCand)
{
  double eps = 10e-10;
  bool result = fabs( aCand.pt() - bCand.pt()) < eps	
    && fabs( aCand.eta() - bCand.eta()) < eps 
    && fabs( aCand.phi() - bCand.phi()) < eps			
    && fabs( aCand.status() - bCand.status()) < eps 
    && fabs(aCand.pdgId() - bCand.pdgId()) < eps;
//  std::cout << aCand.pt() - bCand.pt()
//	    <<" "<< aCand.eta() - bCand.eta()
//	    <<" "<< aCand.phi() - bCand.phi()
//	    <<" "<< aCand.status() - bCand.status()
//	    <<" "<< aCand.pdgId() - bCand.pdgId()<< std::endl;
  if( !result ){
    for(size_t j = 0; j < aCand.numberOfDaughters(); ++j){
      result |= matches( *(aCand.daughter(j)), bCand );      
    }
  }
  return result;
}

bool 
MCDiLeptonHistograms::contains ( collection col, const reco::Candidate& a)
{
  bool result = false;
  for (collection::const_iterator it = col.begin(); it != col.end(); ++it){
    
    result |= matches(a, *it);
  }
  return result;
}

// ------------ method called once each job just before starting event loop  ------------
void 
MCDiLeptonHistograms::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MCDiLeptonHistograms::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(MCDiLeptonHistograms);
