/*
 * LeptonKindFunctor.h
 *
 *  Created on: 03.02.2011
 *      Author: sprenger
 */

#ifndef LEPTONKINDFUNCTOR_H_
#define LEPTONKINDFUNCTOR_H_

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Tau.h>

#include <iostream>


class LeptonKindFunctor
{
public:
	LeptonKindFunctor();
	template<class T> const int operator()(const T& lepton);
	
private:
	template<class T> const int GetLeptKind(const T * lepton);
	const int PromptCategory(const reco::Candidate * genParticle);
};


template<class T>
const int LeptonKindFunctor::operator()(const T& lepton)
{
  return GetLeptKind(&lepton);
}


template<class T> 
const int LeptonKindFunctor::GetLeptKind(const T * lepton)
{
    int value = 0;
    if(lepton->hasUserInt("classByHitsGlb")){
        if(lepton->userInt("classByHitsGlb") == 4 || (lepton->userInt("classByHitsGlb") == 3 && lepton->userInt("classByHitsGlb:flav") == 15) ) value = 3; //Prompt
        if(lepton->userInt("classByHitsGlb") == 3 && lepton->userInt("classByHitsGlb:flav") != 15) value = 4; //Heavy
        if(lepton->userInt("classByHitsGlb") == 2) value = 5; //Light
        if(lepton->userInt("classByHitsGlb") <= 1 || lepton->userInt("classByHitsGlb") > 4) value = 2; //Fake
    }
    else if(lepton->genLepton()){
        const reco::Candidate * genLept = lepton->genLepton();
        value = PromptCategory(genLept);
    }
    else value=2;
    //std::cout << " leptKindResult: "<< value <<std::endl;
    return value;
}


#endif /* LEPTONKINDFUNCTOR_H_ */
