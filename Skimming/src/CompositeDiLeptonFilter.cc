// -*- C++ -*-
//
// Package:    CompositeDiLeptonFilter
// Class:      CompositeDiLeptonFilter
// 
/**\class CompositeDiLeptonFilter CompositeDiLeptonFilter.cc SuSyAachen/Skimming/src/CompositeDiLeptonFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
 */
//
// Original Author:  Matthias Edelhoff
//         Created:  Mon Nov 16 11:26:19 CET 2009
// $Id: CompositeDiLeptonFilter.cc,v 1.6 2010/05/14 11:10:08 edelhoff Exp $
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

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

//#include "DataFormats/Math/interface/deltaR.h"

//include "CommonTools/UtilAlgos/interface/DeltaR.h"


//
// class declaration
//

class CompositeDiLeptonFilter : public edm::EDFilter {
public:
	explicit CompositeDiLeptonFilter(const edm::ParameterSet&);
	~CompositeDiLeptonFilter();

private:
	typedef reco::CompositeCandidate cand;
	typedef edm::View<cand> collection;
	virtual void beginJob() ;
	virtual bool filter(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;

	bool filterChannel(const collection &col, const std::string mode);

	// ----------member data ---------------------------
	std::vector<std::string> channels_;
	std::map<std::string, std::string> modes_;
	std::map<std::string, edm::InputTag> inputTags_;

	bool sameSign_;
	bool debug_;

	class ChargeCombination
	{
	public:
		ChargeCombination(const collection& col, bool debug=false);
		bool isOppositeSign(){return absChargeSum_ == 0;};
		bool operator<(const ChargeCombination& cc){return metric_ < cc.metric();}
		double metric() const {return metric_;}
	private:
		double calcMetric(const cand& theCand) const {return fabs(theCand.daughter(0)->pt()) + fabs(theCand.daughter(1)->pt());}
		double metric_;
		double absChargeSum_;
		bool debug_;
	};

};

CompositeDiLeptonFilter::CompositeDiLeptonFilter(const edm::ParameterSet& iConfig):
		channels_( iConfig.getParameter< std::vector<std::string> > ("channels") ),
		sameSign_(iConfig.getParameter<bool> ("sameSign")),
		debug_(true)
{
	for(std::vector<std::string>::const_iterator it = channels_.begin(); it != channels_.end(); ++it){
		modes_[*it] = iConfig.getParameter<std::string> (*it);
		inputTags_[*it] = iConfig.getParameter<edm::InputTag> (*it+"Src");
	}
}


CompositeDiLeptonFilter::~CompositeDiLeptonFilter()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
CompositeDiLeptonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	bool result = true;
	bool useInclusive = false;
	bool inclusiveResult = false;
	std::list<ChargeCombination> leadingChargeCombinations;

	for(std::vector<std::string>::const_iterator it = channels_.begin(); it != channels_.end(); ++it){
		if(debug_) std::cout << *it <<": ";
		edm::Handle< collection > channelHandle;
		iEvent.getByLabel(inputTags_[*it], channelHandle);
		if(modes_[*it] == "inclusive"){
			useInclusive = true;
			inclusiveResult |= filterChannel(*channelHandle, modes_[*it]);
		}
		else
			result &= filterChannel(*channelHandle, modes_[*it]);
		if(modes_[*it] == "once" || modes_[*it] == "require" || modes_[*it] == "inclusive" ){
			if(debug_) std::cout << "(";
			leadingChargeCombinations.push_back( ChargeCombination(*channelHandle, debug_) );
			if(debug_) std::cout << ") ";
		}
	}
	leadingChargeCombinations.sort();
	if(useInclusive)
		result &= inclusiveResult;
	result &= ( sameSign_ && !leadingChargeCombinations.front().isOppositeSign() ) ||
			  (!sameSign_ &&  leadingChargeCombinations.front().isOppositeSign() );
	if(debug_) std::cout << "-> "<< result << std::endl;
	return result;
}

bool
CompositeDiLeptonFilter::filterChannel(const collection &col, const std::string mode)
{
	bool result = false;
	if(debug_) std::cout << mode;
	if(mode == "reject") result = col.size() == 0;
	else if(mode == "require" || mode == "inclusive") result = col.size() > 0;
	else if(mode == "allow") result = true;
	else if(mode == "once") result = col.size() == 1;
	if(debug_) std::cout <<" " << result <<" ["<<col.size()<<"], ";
	return result;
}

// ------------ method called once each job just before starting event loop  ------------
void 
CompositeDiLeptonFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CompositeDiLeptonFilter::endJob() {
}


CompositeDiLeptonFilter::ChargeCombination::ChargeCombination(const collection& col, bool debug):
		debug_(debug)
{
	metric_ = -1;
	absChargeSum_ = -1;
	for(collection::const_iterator it = col.begin(); it != col.end(); ++it){
		if( metric_ < calcMetric(*it))
			metric_ = calcMetric(*it);
			absChargeSum_ = fabs((*it).daughter(0)->charge() + (*it).daughter(1)->charge());
		if(debug_) std::cout << metric_<<", "<<absChargeSum_;
	}
}
//define this as a plug-in
DEFINE_FWK_MODULE(CompositeDiLeptonFilter);
