/*
 * WeightFunctor.h
 *
 *  Created on: 21.12.2010
 *      Author: heron
 */

#ifndef WEIGHTFUNCTOR_H_
#define WEIGHTFUNCTOR_H_

#include <vector>
#include <iostream>

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Tau.h>


class WeightFunctor
{
public:
	WeightFunctor();
	void SetSource( const edm::ParameterSet& p, const std::string name);
	template<class T> double operator()(const T& lepton, const std::string type);
	double operator()(const pat::Electron& lepton){return this->operator()(lepton, "electrons");};
	double operator()(const pat::Muon& lepton){return this->operator()(lepton, "muons");};
	double operator()(const pat::Tau& lepton){return this->operator()(lepton, "taus");};
	bool isUseable(){return initialized;};


private:
	bool isInBin(const edm::ParameterSet& p, const std::string name, double value);

	edm::ParameterSet	src_;
	bool initialized;
};

template<class T>
double
WeightFunctor::operator()(const T& lepton, const std::string type)
{
	double result = 1.0;
	std::map<std::string, double> binVariables;
	binVariables["pt"] = lepton.pt();
	binVariables["eta"] = lepton.eta();
	if(src_.existsAs<edm::VParameterSet>(type)){
		edm::VParameterSet bins = src_.getParameter<edm::VParameterSet>(type);
		for(edm::VParameterSet::const_iterator it = bins.begin(); it != bins.end(); ++it){
			bool inBin = true;
			for(std::map<std::string, double>::const_iterator binName = binVariables.begin(); binName != binVariables.end(); ++binName)
				inBin &= isInBin(*it, (*binName).first, (*binName).second);

			if(inBin)
				result = (*it).getParameter<double>("weight");
		}
	}

	return result;
}


#endif /* WEIGHTFUNCTOR_H_ */
