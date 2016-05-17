// -*- C++ -*-
//
// Package:    SignatureFilter
// Class:      SignatureFilter
//
/**\class SignatureFilter SignatureFilter.cc SuSyAachen/Skimming/src/SignatureFilter.cc

 Description: Filters for a configurable signature given the genParticle collection

 Implementation:
     vectors for IDs and status are given. any combination of both has to exist.
     if the daughter VPSet exists and is filled the AND f all declared particles as to exist.
     Again the definition can incorporate a arbitrary combination of IDs and status

     Example usage:
		IDs = cms.vint32( 23 ),
		status = cms.vint32(3),
		daughters = cms.VPSet(
			cms.PSet(
				IDs = cms.vint32( 15 ),
				status = cms.vint32(3)
			),
			cms.PSet(
				IDs = cms.vint32( -15 ),
				status = cms.vint32(3)
			),
		)

 */
//
// Original Author:  Matthais Edelhoff
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//
// class declaration
//

class SignatureFilter : public edm::EDFilter {
public:
	explicit SignatureFilter(const edm::ParameterSet&);
	~SignatureFilter();

private:
	virtual void beginJob() ;
	virtual bool filter(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;


	class ParticleTemplate{
	public:
		ParticleTemplate( const edm::ParameterSet& pSet );
		ParticleTemplate( const std::vector< int> &ids, const std::vector<int> &stauts):ids_(ids), status_(stauts){};

		bool match(const reco::GenParticle &particle);
		void addDaughter(const edm::ParameterSet &pSet){ daughters_.push_back(ParticleTemplate(pSet) ); };
	private:
		std::vector<int> ids_;
		std::vector<int> status_;
		std::vector<ParticleTemplate> daughters_;

	};
	// ----------member data ---------------------------

	edm::EDGetTokenT< reco::GenParticleCollection > inputToken_;
	//instead of a single mother we could go for a vector of ParticleTemplates and read a VPSet then the AND of mother would prevail
	ParticleTemplate mother_;

	bool debug_;

};

// constructors and destructor
SignatureFilter::SignatureFilter(const edm::ParameterSet& iConfig):
		inputToken_(consumes< reco::GenParticleCollection >(iConfig.getParameter<edm::InputTag>("src"))),
		mother_(iConfig)
{
	debug_ = false;

}

SignatureFilter::~SignatureFilter()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SignatureFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	bool result = false;
	edm::Handle< reco::GenParticleCollection > genParticles;
	iEvent.getByToken(inputToken_, genParticles);

	for(reco::GenParticleCollection::const_iterator it = genParticles->begin();
			it != genParticles->end(); ++it){
		result |= mother_.match(*it);
	}
	return result;
}

// ------------ method called once each job just before starting event loop  ------------
void
SignatureFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
SignatureFilter::endJob() {
}
//#######################  SignatureFilter::ParticleTemplate  #############
SignatureFilter::ParticleTemplate::ParticleTemplate(const edm::ParameterSet& pSet):ids_(pSet.getParameter< std::vector<int> > ("IDs")),
		  status_(pSet.getParameter< std::vector<int> > ("status"))
{
	if(pSet.exists("daughters")){
		std::vector<edm::ParameterSet> daughterSets = pSet.getParameter<std::vector<edm::ParameterSet> >("daughters");
		for( std::vector<edm::ParameterSet>::const_iterator iSet = daughterSets.begin();
				iSet != daughterSets.end(); ++iSet){
			//if(debug_) std::cout << iSet->dump() <<std::endl;
			this->addDaughter(*iSet);
		}
	}
}

bool
SignatureFilter::ParticleTemplate::match(const reco::GenParticle &particle)
{
	bool id = false;
	bool status = false;
	bool daughtersFound = true;

	for(std::vector<int>::iterator it = ids_.begin();
			it != ids_.end(); ++it)
		id |= ((*it) == particle.pdgId());

	for(std::vector<int>::iterator it = status_.begin();
				it != status_.end(); ++it)
		status |= ((*it) == particle.status());

	for (std::vector<ParticleTemplate>::iterator it = daughters_.begin();
			it != daughters_.end(); ++it){
		bool match = false;
		reco::GenParticleRefVector daughters = particle.daughterRefVector();
		for(reco::GenParticleRefVector::const_iterator itDaughter = daughters.begin();
				itDaughter != daughters.end(); ++itDaughter)
			match |= (*it).match(*(*itDaughter) );
		daughtersFound &= match;
	}
	return id && status && daughtersFound;
}


//define this as a plug-in
DEFINE_FWK_MODULE(SignatureFilter);
