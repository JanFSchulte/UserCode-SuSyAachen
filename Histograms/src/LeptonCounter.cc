// -*- C++ -*-
//
// Package:    Histograms
// Class:      LeptonCounter
// 
/**\class LeptonCounter LeptonCounter.cc brot/LeptonCounter/src/LeptonCounter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  matthias edelhoff
//         Created:  Tue Oct 27 13:50:40 CET 2009
// $Id: LeptonCounter.cc,v 1.4 2010/04/01 09:36:57 sprenger Exp $
//
//


// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <DataFormats/Candidate/interface/Candidate.h>
//#include <DataFormats/PatCandidates/interface/Electron.h>

//ROOT
#include "TH1.h"
//
// class decleration
//

class LeptonCounter : public edm::EDAnalyzer {
public:
  explicit LeptonCounter(const edm::ParameterSet&);
  ~LeptonCounter();

  
private:
  typedef edm::View<reco::Candidate> collection;
  enum leptonCombi { ee, mumu, tautau, emu, etau, mutau };
  enum lepton {e, mu, tau};
  
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void fillInclusive( const collection &electrons, const collection &muons, const collection &taus );
  
  void countByCharge( const collection &input, int &nPlus, int &nMinus);
  void fillCombination(int combination, const collection &input);
  void fillCombination(int combination, const collection &input1xs, const collection &input2);
  void setBinLabels(TH1F* histo);
  // ----------member data ---------------------------
  //names
  std::string subDir_;
  std::string method_;

  edm::InputTag electronTag_;
  edm::InputTag muonTag_;
  edm::InputTag tauTag_;

  //histos
  TH1F* sameSignLeptons_; 
  TH1F* oppositeSignLeptons_; 

  TH1F* leptonCount_[3];  
};

// constructors and destructor
LeptonCounter::LeptonCounter(const edm::ParameterSet& iConfig)
{
  // read config
  subDir_ = iConfig.getParameter<std::string> ("subDir");
  method_ = iConfig.getParameter<std::string> ("method");

  electronTag_ = iConfig.getParameter<edm::InputTag> ("electronSource");
  muonTag_ = iConfig.getParameter<edm::InputTag> ("muonSource");
  tauTag_ = iConfig.getParameter<edm::InputTag> ("tauSource");

  // init histos
  edm::Service<TFileService> fs;
  //  fs->mkdir(subDir_);
  sameSignLeptons_ = fs->make<TH1F>("SSFlavor" , "Same Sign Flavor Combination" , 6 , -0.5 , 5.5 );
  setBinLabels(sameSignLeptons_);
  oppositeSignLeptons_ = fs->make<TH1F>("OSFlavor" , "Oposite Sign Flavor Combination" , 6 , -0.5 , 5.5 );
  setBinLabels(oppositeSignLeptons_);

  leptonCount_[e] = fs->make<TH1F>("NumElectrons" , "Number of Electrons" , 21 , -0.5 , 20.5 );
  leptonCount_[mu] = fs->make<TH1F>("NumMuons" , "Number of Muons" , 21 , -0.5 , 20.5 );
  leptonCount_[tau] = fs->make<TH1F>("NumTaus" , "Number of Taus" , 21 , -0.5 , 20.5 );
}


LeptonCounter::~LeptonCounter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// member functions
// ------------ method called to for each event  ------------
void
LeptonCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //Electrons
  edm::Handle< collection > electrons;
  iEvent.getByLabel(electronTag_, electrons);
 
  //Muons
  edm::Handle< collection > muons;
  iEvent.getByLabel(muonTag_, muons);
  
  //Taus
  edm::Handle< collection > taus;
  iEvent.getByLabel(tauTag_, taus);

  if(method_ == "inclusive") fillInclusive( *electrons, *muons, *taus);
  if(method_ == "exclusive") {
    int sumLeptons = electrons->size() + muons->size() + taus->size(); 
    if( sumLeptons == 2){
      fillInclusive( *electrons, *muons, *taus);
    }
  }

  leptonCount_[e]->Fill(electrons->size());
  leptonCount_[mu]->Fill(muons->size());
  leptonCount_[tau]->Fill(taus->size());

}

void LeptonCounter::fillInclusive( const collection &electrons, const collection &muons, const collection &taus )
{
  fillCombination(ee, electrons);
  fillCombination(mumu, muons);
  fillCombination(tautau, taus);

  fillCombination(emu, electrons, muons);
  fillCombination(etau, electrons, taus);
  fillCombination(mutau, muons, taus);
}

void LeptonCounter::fillCombination(int combination, const collection &input)
{
  int nPlus = 0;
  int nMinus = 0;
  countByCharge(input, nPlus, nMinus);
  if( nPlus >= 1 && nMinus >= 1) {
    oppositeSignLeptons_->Fill( combination );
  }
  if( nPlus >= 2 || nMinus >= 2) {
    sameSignLeptons_->Fill( combination );
    // if (combination == 2)std::cout << "count"<<std::endl;
  }
  //  if (combination == 2) std::cout << nPlus << ":"<<nMinus<< std::endl;
}

void LeptonCounter::fillCombination(int combination, const collection &input1, const collection &input2)
{
  int n1Plus = 0;
  int n1Minus = 0;
  int n2Plus = 0;
  int n2Minus = 0;
  countByCharge(input1, n1Plus, n1Minus);
  countByCharge(input2, n2Plus, n2Minus);
  if( (n1Plus >= 1 && n2Minus >= 1) || (n2Plus >= 1 && n1Minus >= 1) ) {
    oppositeSignLeptons_->Fill( combination );
  }
  if( (n1Plus >= 1 && n2Plus >=1) || (n1Minus >= 1 && n2Minus >=1)) {
    sameSignLeptons_->Fill( combination );
  }
}

void LeptonCounter::countByCharge( const collection &input, int &nPlus, int &nMinus)
{
  nPlus = nMinus = 0;
  for(collection::const_iterator it = input.begin(); it != input.end() ; ++it) {
    if((*it).charge() == 1.) nPlus++;
    if((*it).charge() == -1.) nMinus++;
  }
}

void LeptonCounter::setBinLabels(TH1F* histo)
{
  // bin counting in root starts at 1 *sigh*
  histo->GetXaxis()->SetBinLabel(ee+1,"ee");
  histo->GetXaxis()->SetBinLabel(mumu+1,"#mu#mu");
  histo->GetXaxis()->SetBinLabel(tautau+1,"#tau#tau");
  histo->GetXaxis()->SetBinLabel(emu+1,"e#mu");
  histo->GetXaxis()->SetBinLabel(etau+1,"e#tau");
  histo->GetXaxis()->SetBinLabel(mutau+1,"#mu#tau");
}

// ------------ method called once each job just before starting event loop  ------------
void 
LeptonCounter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LeptonCounter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(LeptonCounter);
