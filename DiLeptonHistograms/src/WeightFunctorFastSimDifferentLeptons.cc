#include "SuSyAachen/DiLeptonHistograms/interface/WeightFunctorFastSimDifferentLeptons.h"

WeightFunctorFastSimDifferentLeptons::WeightFunctorFastSimDifferentLeptons()
{
	initialized_ = false;
	bool looseNotTight_ = false;
}

void
WeightFunctorFastSimDifferentLeptons::SetSource( const edm::ParameterSet& p, const std::string rawName, const bool mc )
{
	std::string name = rawName;
	if(p.existsAs<edm::ParameterSet>(name)){
		initialized_ = true;
		src_ = edm::ParameterSet(p.getParameter< edm::ParameterSet >(name));
	}
}

bool
WeightFunctorFastSimDifferentLeptons::isInBin(const edm::ParameterSet& p, const std::string name, double value)
{
	bool result = false;
	if(p.existsAs<double>(name+"Min") && p.existsAs<double>(name+"Max")){
		double min = p.getParameter<double>(name+"Min");
		double max = p.getParameter<double>(name+"Max");
		result = value >= min && value < max;
	}
	return result;
}

