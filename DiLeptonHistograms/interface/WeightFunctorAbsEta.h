/*
 * WeightFunctorAbsEta.h
 *
 *  Created on: 21.12.2010
 *      Author: heron
 */

#ifndef WEIGHTFUNCTORABSETA_H_
#define WEIGHTFUNCTORABSETA_H_

#include <vector>
#include <iostream>

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Tau.h>


class WeightFunctorAbsEta
{
public:
	WeightFunctorAbsEta();
	void SetSource( const edm::ParameterSet& p, const std::string name, const bool mc=false);
	template<class T> double operator()(const T& lepton, const std::string type);
	double operator()(const pat::Electron& lepton){return this->operator()(lepton, "electrons");};
	double operator()(const pat::Muon& lepton){return this->operator()(lepton, "muons");};
	double operator()(const pat::Tau& lepton){return this->operator()(lepton, "taus");};
	void looseNotTight(){looseNotTight_ = true;}
	void plain(){looseNotTight_ = false;}
	bool isUseable(){return initialized_;};


private:
	bool isInBin(const edm::ParameterSet& p, const std::string name, double value);

	edm::ParameterSet	src_;
	bool initialized_;
	bool looseNotTight_;
};

template<class T>
double
WeightFunctorAbsEta::operator()(const T& lepton, const std::string type)
{
	double result = 1.0;
	std::map<std::string, double> binVariables;
	binVariables["pt"] = lepton.pt();
	binVariables["eta"] = fabs(lepton.eta());
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
	if(looseNotTight_) result = result/(1. - result);
	return result;
}


#endif /* WEIGHTFUNCTORABSETA_H_ */
