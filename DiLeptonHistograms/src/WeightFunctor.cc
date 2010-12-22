#include "SuSyAachen/DiLeptonHistograms/interface/WeightFunctor.h"

WeightFunctor::WeightFunctor()
{
	initialized = false;
}

void
WeightFunctor::SetSource( const edm::ParameterSet& p, const std::string name )
{
	if(p.existsAs<edm::ParameterSet>(name)){
		initialized = true;
		src_ = edm::ParameterSet(p.getParameter< edm::ParameterSet >(name));
	}
}

bool
WeightFunctor::isInBin(const edm::ParameterSet& p, const std::string name, double value)
{
	bool result = false;
	if(p.existsAs<double>(name+"Min") && p.existsAs<double>(name+"Max")){
		double min = p.getParameter<double>(name+"Min");
		double max = p.getParameter<double>(name+"Max");
		result = value >= min && value < max;
	}
	return result;
}

