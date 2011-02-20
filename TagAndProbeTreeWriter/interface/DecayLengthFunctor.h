/*
 * DecayLenghtFunctor.h
 *
 *  Created on: 17.02.2011
 *      Author: sprenger
 */

#ifndef DECAYLENGTHFUNCTOR_H_
#define DECAYLENGTHFUNCTOR_H_

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <iostream>


class DecayLengthFunctor
{
public:
	DecayLengthFunctor();
	//template<class T> const int operator()(const T& lepton);
	const double operator()(const pat::Electron& e1, const pat::Electron& e2,
				const reco::Vertex& primaryVtx, const edm::EventSetup& iSetup);
	
private:
	const double GetDecayLength(const pat::Electron& e1, const pat::Electron& e2,
				    const reco::Vertex& primaryVtx, const edm::EventSetup& iSetup);
};

#endif /* DECAYLENGTHFUNCTOR_H_ */
