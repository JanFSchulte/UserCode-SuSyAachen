/*
 * WeightFunctorFastSim.h
 *
 *  Created on: 21.12.2010
 *      Author: heron
 */

#ifndef WEIGHTFUNCTORFASTSIM_H_
#define WEIGHTFUNCTORFASTSIM_H_

#include <vector>
#include <iostream>

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Tau.h>


class WeightFunctorFastSim
{
public:
	WeightFunctorFastSim();
	void SetSource( const edm::ParameterSet& p, const std::string name, const bool mc=false);


// 	double operator()(const pat::Muon& lepton,const pat::Muon& lepton2,const double nJets){return this->operator()(lepton,lepton2,nJets, "dimuons");};
// 	double operator()(const pat::Tau& lepton){return this->operator()(lepton, "taus");};

	template<class T> double operator()(const T& lepton, const T& lepton2, unsigned int& nJets,const std::string mode, const std::string type);
	double operator()(const pat::Electron& lepton,const pat::Electron& lepton2,unsigned int& nJets,const std::string mode){return this->operator()(lepton,lepton2,nJets,mode, "dielectrons");};
// 	double operator()(const pat::Electron& lepton,const pat::Muon& lepton2,const double nJets){return this->operator()(lepton,lepton2,nJets, "electronmuons");};
	double operator()(const pat::Muon& lepton,const pat::Muon& lepton2,unsigned int& nJets,const std::string mode){return this->operator()(lepton,lepton2,nJets,mode, "dimuons");};
	double operator()(const pat::Tau& lepton,const pat::Tau& lepton2,unsigned int& nJets,const std::string mode){return this->operator()(lepton,lepton2,nJets,mode, "ditaus");};

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
WeightFunctorFastSim::operator()(const T& lepton, const T& lepton2,unsigned int& nJets,const std::string mode, const std::string type)
{
	double result = 1.0;

	std::map<std::string, double> binVariables;
	binVariables["pt1"] = lepton.pt();
	binVariables["eta1"] = fabs(lepton.eta());
	binVariables["pt2"] = lepton2.pt();
	binVariables["eta2"] = fabs(lepton2.eta());
	binVariables["nJets"] = nJets;
	if(src_.existsAs<edm::VParameterSet>(type)){
		edm::VParameterSet bins = src_.getParameter<edm::VParameterSet>(type);
		for(edm::VParameterSet::const_iterator it = bins.begin(); it != bins.end(); ++it){
			bool inBin = true;
			for(std::map<std::string, double>::const_iterator binName = binVariables.begin(); binName != binVariables.end(); ++binName)
				inBin &= isInBin(*it, (*binName).first, (*binName).second);

			if(inBin){
				if (mode.compare("Up") == 0){
					result = (*it).getParameter<double>("weight")+(*it).getParameter<double>("weight")*(*it).getParameter<double>("uncert");

				}
				else if (mode.compare("Down") == 0){
					result = (*it).getParameter<double>("weight")-(*it).getParameter<double>("weight")*(*it).getParameter<double>("uncert");
				}
				else{
					result = (*it).getParameter<double>("weight");
				}
			}
		}
	}
	if(looseNotTight_) result = result/(1. - result);
	return result;
}



#endif /* WEIGHTFUNCTORFASTSIM_H_ */
