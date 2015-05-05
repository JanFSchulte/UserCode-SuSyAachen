#ifndef effectiveAreaIsolationSelector_h
#define effectiveAreaIsolationSelector_h

//DataFormats
//include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"

//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "SuSyAachen/TagAndProbeTreeWriter/interface/IsolationFunctor.h"

//STL
#include <vector>

template<typename rhoType, typename collectionType, typename containerType>
struct effectiveAreaIsolationSelector {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  effectiveAreaIsolationSelector ( const edm::ParameterSet & cfg, edm::ConsumesCollector ):
    fctIsolation_ (cfg.getParameter<edm::ParameterSet>("isolationDefinitions")),
    isoMin_( cfg.getParameter<double>( "isoMin") ),
    isoMax_( cfg.getParameter<double>( "isoMax" ) ),
    isoMaxEE_(0.1),
    rhoSrc_( cfg.getParameter<edm::InputTag>( "rhoSource" ) )  { }
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {
    
    fctIsolation_.init(ev);

    edm::Handle<rhoType> rhoIso_h;
    ev.getByLabel(rhoSrc_, rhoIso_h);
    double rhoIso = *(rhoIso_h.product());

    selected_.clear();
    int i = 0;
    for(typename collection::const_iterator it = col.product()->begin(); 
	it != col.product()->end(); ++it ){
      
      double caloIso = 0.;
      //      caloIso += (*it).pfIsolationVariables().neutralHadronIso;
      //      caloIso += (*it).pfIsolationVariables().photonIso;
      caloIso +=  (*it).neutralHadronIso();
      caloIso += (*it).photonIso();

      caloIso -= getAEff((*it).eta()) * rhoIso;
  
      double iso = 0.;
      //      iso += (*it).pfIsolationVariables().chargedHadronIso;
      iso += (*it).chargedHadronIso();

      if (caloIso > 0)
	iso += caloIso;
      iso /= (*it).pt();
      if (isoMin_ < 0 &&  isoMax_ < 0 && col.product()->size() == 2){
	/*	std::cerr << i 
		  << "\t eventNr: "<<ev.id().event()
                  << "pt "<<(*it).pt() << std::endl
		  << "neut "<< (*it).pfIsolationVariables().neutralHadronIso 
                  << ", phot "<<(*it).pfIsolationVariables().photonIso 
                  << ", char "<<(*it).pfIsolationVariables().chargedHadronIso
		  << std::endl
		  << "neut "<< (*it).neutralHadronIso() 
                  << ", phot "<<(*it).photonIso()
                  << ", char "<<(*it).chargedHadronIso() 
                  << std::endl;*/
	/*
	std::cerr //<<<  col.product()->size() << std::endl
		  << i<< ";"
		  << ev.id().luminosityBlock() << ";"<< ev.id().run() << ";"<< ev.id().event()<< ";"
		  << (*it).pt() << ";" 
		  << (*it).pfIsolationVariables().neutralHadronIso << ";"
		  << (*it).pfIsolationVariables().photonIso << ";"
		  << (*it).pfIsolationVariables().chargedHadronIso << ";"
		  << getAEff((*it).eta()) * rhoIso 
		  << std::endl;
	*/
	i++;
      }

      bool passesIsolation=false;
      if(iso < isoMaxEE_ && (*it).pt() < 20. && (*it).isEE()) passesIsolation=true;
      if(iso < isoMax_   && (*it).pt() > 20. && (*it).isEE()) passesIsolation=true;
      if(iso < isoMax_   &&                    !(*it).isEE()) passesIsolation=true;
      if( isoMax_ < 0)  passesIsolation=true;
      
      if( (iso > isoMin_ || isoMin_ < 0) && passesIsolation){
	
	/*std::cout << "++picked+++> pt "<< (*it).pt()<<", iso " << iso <<" isEB "<< (*it).isEB()<<std::endl;
	std::cout << "   neutr "<< (*it).pfIsolationVariables().neutralHadronIso 
		  << " gamma "<< (*it).pfIsolationVariables().photonIso 
		  << " cha "<< (*it).pfIsolationVariables().chargedHadronIso 
		  << " pu "<<getAEff((*it).eta()) * rhoIso << std::endl;*/
	  selected_.push_back( & (*it) );
      }
    }
  }

  size_t size() const { return selected_.size(); }
  
  private:
  double getAEff(double eta)
  {
    //from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/src/EGammaCutBasedEleId.cc
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
    
    container selected_;
    IsolationFunctor fctIsolation_;
    double isoMin_;
    double isoMax_;
    double isoMaxEE_;
    edm::InputTag rhoSrc_;
};

#endif
