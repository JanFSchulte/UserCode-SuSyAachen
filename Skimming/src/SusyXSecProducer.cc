// -*- C++ -*-
//
// Package:    SusyXSecProducer
// Class:      SusyXSecProducer
// 
/**\class SusyXSecProducer SusyXSecProducer.cc SuSyAachen/Skimming/src/SusyXSecProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Niklas Mohr
//         Created:  Wed Aug 18 15:37:34 CEST 2010
// $Id: SusyXSecProducer.cc,v 1.1 2011/01/07 15:47:27 nmohr Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


//
// class declaration
//

class CMSSMStruct{
    /* 
#                   0 ng     neutralino/chargino + gluino                     #
#                   1 ns     neutralino/chargino + squark                     #
#                   2 nn     neutralino/chargino pair combinations            #
#                   3 ll     slepton pair combinations                        #
#                   4 sb     squark-antisquark                                #
#                   5 ss     squark-squark                                    #
#                   6 tb     stop-antistop                                    #
#                   7 bb     sbottom-antisbottom                              #
#                   8 gg     gluino pair                                      #
#                   9 sg     squark + gluino                                  #
*/

    std::map<std::string,double> map_; 
    std::vector<double> val_;
    public:
        void set_values(std::vector<std::string>,edm::ParameterSet);
        void set_values(std::map<std::string,double>);
        double getXSec(unsigned int index);
        bool operator==(CMSSMStruct);
};

void CMSSMStruct::set_values(std::vector<std::string> vars, edm::ParameterSet values) {
    for( std::vector<std::string>::iterator var_i = vars.begin(); var_i != vars.end(); ++var_i ) {
        double val = values.getParameter<double>( *var_i );
        map_[*var_i] = val;
    } 
    val_.push_back(values.getParameter<double>("ng"));
    val_.push_back(values.getParameter<double>("ns"));
    val_.push_back(values.getParameter<double>("nn"));
    val_.push_back(values.getParameter<double>("ll"));
    val_.push_back(values.getParameter<double>("sb"));
    val_.push_back(values.getParameter<double>("ss"));
    val_.push_back(values.getParameter<double>("tb"));
    val_.push_back(values.getParameter<double>("bb"));
    val_.push_back(values.getParameter<double>("gg"));
    val_.push_back(values.getParameter<double>("sg"));

}

bool CMSSMStruct::operator==(CMSSMStruct a)
{
  return (*this).map_ == a.map_;
}


void CMSSMStruct::set_values(std::map<std::string,double> map) {
    map_ = map;
    for (int i = 0; i < 10; ++i) val_.push_back(0.);
}


double CMSSMStruct::getXSec(unsigned int index) {
    if (index < val_.size() ) return val_[index];
    else return 0.;
}


class SusyXSecProducer : public edm::EDProducer {
   public:
      explicit SusyXSecProducer(const edm::ParameterSet&);
      ~SusyXSecProducer();

   private:
      edm::InputTag src;
      std::vector<edm::ParameterSet> susyVars_;
      std::vector<edm::ParameterSet> scan_;
      bool debug;
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual int getSUSYSubprocess(const edm::Handle< std::vector<reco::GenParticle> >& );
      virtual int sFinalState(int, int);
      virtual double susyNLOXSec(int, CMSSMStruct);
      
      // ----------member data ---------------------------
      std::vector<CMSSMStruct> theNLO_;
      std::vector<std::string> theVars_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
SusyXSecProducer::SusyXSecProducer(const edm::ParameterSet& iConfig)
{ 
   produces< int > ( ); //"susyScanProductionProcess" );
   produces< double > ( );//"susyScanNLOCrossSection" );

   debug = false;
   
   src = iConfig.getParameter<edm::InputTag> ("src");
   susyVars_ = iConfig.getParameter< std::vector<edm::ParameterSet> >("susyVars");
   scan_ = iConfig.getParameter< std::vector<edm::ParameterSet> >("crossSections");
   for ( std::vector<edm::ParameterSet>::iterator susyVar_i = susyVars_.begin(); susyVar_i != susyVars_.end(); ++susyVar_i ) {
        std::string var = susyVar_i->getParameter<std::string>( "var" );
        theVars_.push_back(var);
   }
   for ( std::vector<edm::ParameterSet>::iterator scan_i = scan_.begin(); scan_i != scan_.end(); ++scan_i ) {
        CMSSMStruct point;
        point.set_values(theVars_,*scan_i);
        theNLO_.push_back(point);
   }
}


SusyXSecProducer::~SusyXSecProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

int SusyXSecProducer::getSUSYSubprocess(const edm::Handle< std::vector<reco::GenParticle> >& genParticles)
{
  std::vector<int> interactions;

  for (std::vector<reco::GenParticle>::const_iterator p_i = genParticles->begin(); p_i != genParticles->end(); ++p_i){
    if (p_i->status() != 3) continue;
    int ID = abs(p_i->pdgId());
    int mID = 999999999;
    if(p_i->mother()) mID = abs(p_i->mother()->pdgId());
    if(debug) std::cout << "Particle with status 3: "<< ID << " has mother:" << mID << std::endl;
    if (ID > 1000000 && ID < 2000016 && (mID < 7 || mID ==21 || mID == 22 || mID == 23 || mID == 24)) { // kept the bosons in case of screw ups

     // Check the mother of the SM is a proton
       bool isProton = false;       
       for (std::vector<reco::GenParticle>::const_iterator p_k = genParticles->begin(); p_k < p_i; ++p_k){
            int pID = abs(p_k->pdgId());
            int mpID = 999999999;
            if(p_k->mother()) mpID = abs(p_k->mother()->pdgId());
         if ((pID < 7 || pID ==21 || pID == 22 || pID == 23 || pID == 24)&&(mpID == 2212)) isProton = true; 
       }
      if (!isProton) continue;
      if(debug) std::cout << "Accepted susy particles for final state: "<< ID << " has mother:" << mID << std::endl;
    

      interactions.push_back(p_i->pdgId());
      if (interactions.size() > 2){ 
         std::cout << std::setw(4) << std::left << " WARNING mcSUSYkfactor: Something is wrong with "
         << std::setw(10) << std::left << ID << " "
         << std::setw(4) << std::right << p_i->status() << " "
         << std::setw(10) << std::left << mID
         << " using k=1 " << std::endl;
      }
    }
  }

  if (interactions.size() == 2)  return sFinalState(interactions[0],interactions[1]);
  else return -1;
}


int SusyXSecProducer::sFinalState(int ipart1, int ipart2) {
/* 
#                   0 ng     neutralino/chargino + gluino                     #
#                   1 ns     neutralino/chargino + squark                     #
#                   2 nn     neutralino/chargino pair combinations            #
#                   3 ll     slepton pair combinations                        #
#                   4 sb     squark-antisquark                                #
#                   5 ss     squark-squark                                    #
#                   6 tb     stop-antistop                                    #
#                   7 bb     sbottom-antisbottom                              #
#                   8 gg     gluino pair                                      #
#                   9 sg     squark + gluino                                  #
*/

 int index = -1;

 // ng case
 if ((abs(ipart1) > 1000021 && abs(ipart1) < 1000039 && abs(ipart2) == 1000021) ||
      (abs(ipart2) > 1000021 && abs(ipart2) < 1000039 && abs(ipart1) == 1000021)) index = 0;
// ns case 
 else if ((abs(ipart1) > 1000021 && abs(ipart1) < 1000039 && (abs(ipart2) < 1000007 || (abs(ipart2) > 2000000 && abs(ipart2) < 2000007))) ||
    (abs(ipart2) > 1000021 && abs(ipart2) < 1000039 && (abs(ipart1) < 1000007 || (abs(ipart1) > 2000000 && abs(ipart1) < 2000007)))) index = 1;
// nn case 
 else if (abs(ipart1) > 1000021 && abs(ipart1) < 1000039 &&  abs(ipart2) > 1000021 && abs(ipart2) < 1000039) index = 2;
// ll case
 else if (((abs(ipart1) > 1000010 && abs(ipart1) < 1000017) || (abs(ipart1) > 2000010 && abs(ipart1) < 2000016)) &&
         ((abs(ipart2) > 1000010 && abs(ipart2) < 1000017) || (abs(ipart2) > 2000010 && abs(ipart2) < 2000016))) index = 3;
// tb case
 else if ((abs(ipart1) == 1000006 || abs(ipart1) == 2000006) && (abs(ipart2) == 1000006 || abs(ipart2) == 2000006)) index = 6;
// bb case
 else if ((abs(ipart1) == 1000005 || abs(ipart1) == 2000005) && (abs(ipart2) == 1000005 || abs(ipart2) == 2000005)) index = 7;
// sb case
 else if ((ipart1 * ipart2 < 0) &&
         (abs(ipart1) < 1000007 || (abs(ipart1) > 2000000 && abs(ipart1) < 2000007)) &&
         (abs(ipart2) < 1000007 || (abs(ipart2) > 2000000 && abs(ipart2) < 2000007))) index = 4;
// ss case
 else if ((ipart1 * ipart2 > 0) &&
         (abs(ipart1) < 1000007 || (abs(ipart1) > 2000000 && abs(ipart1) < 2000007)) &&
         (abs(ipart2) < 1000007 || (abs(ipart2) > 2000000 && abs(ipart2) < 2000007))) index = 5;
// gg case 
  else if (ipart1 == 1000021 && ipart2 == 1000021) index = 8;
// sg case 
 else if ((abs(ipart1) == 1000021 && (abs(ipart2) < 1000007 || (abs(ipart2) > 2000000 && abs(ipart2) < 2000007))) ||
    (abs(ipart2) == 1000021 && (abs(ipart1) < 1000007 || (abs(ipart1) > 2000000 && abs(ipart1) < 2000007)))) index = 9;
 else
   std::cout << "Warning mcSUSYkfactor: index not defined " << ipart1 << "  " << ipart2 << std::endl;

  return index;
}

double SusyXSecProducer::susyNLOXSec(int index, CMSSMStruct point) {

  if (index < 0) return 0;

  std::vector<CMSSMStruct>::iterator it;
  it = find (theNLO_.begin(), theNLO_.end(), point);

  if (it == theNLO_.end()) return 0.;
  else return (*it).getXSec(index);
    
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
SusyXSecProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::map<std::string,double> map_;
    for( std::vector<edm::ParameterSet>::iterator susyVar_i = susyVars_.begin(); susyVar_i != susyVars_.end(); ++susyVar_i ) {
        std::string var = susyVar_i->getParameter<std::string>( "var" );
        edm::Handle< double > var_;
        iEvent.getByLabel("seqSUSYPARS",var, var_);
        map_[var] = (*var_);
    } 
    CMSSMStruct point;
    point.set_values(map_);

    edm::Handle< std::vector<reco::GenParticle> > genParticles;
    iEvent.getByLabel(src, genParticles);
    
    int index = getSUSYSubprocess(genParticles);
    double xSec = susyNLOXSec(index,point);
    if(debug){
        std::cout << "I am alive" << std::endl;
        std::cout << index << std::endl;
        std::cout << xSec << std::endl;
    }
   
    std::auto_ptr< int > theIndex ( new int(index) );
    std::auto_ptr< double > theXSec ( new double(xSec) );
   
    iEvent.put( theIndex ); //, "susyScanProductionProcess" );
    iEvent.put( theXSec );   //, "susyScanNLOCrossSection" );
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
SusyXSecProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SusyXSecProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(SusyXSecProducer);
