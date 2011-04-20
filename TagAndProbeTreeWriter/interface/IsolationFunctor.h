/*
 * IsolationFunctor.h
 *
 *  Created on: 20.04.2011
 *      Author: sprenger
 */

#ifndef ISOLATIONFUNCTOR_H_
#define ISOLATIONFUNCTOR_H_

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include "DataFormats/Candidate/interface/Candidate.h"

#include <iostream>


// base functor class
//====================
class IsolationFunctor
{
public:
  IsolationFunctor();
  template<class T> const double operator()(const T& lepton);
	
private:
  const virtual double GetIsolation(const pat::Electron& lepton);
  const virtual double GetIsolation(const pat::Muon& lepton);
  const virtual double GetIsolation(const pat::Tau& lepton);
  const virtual double GetIsolation(const reco::Candidate& lepton);
};


template<class T>
const double IsolationFunctor::operator()(const T& lepton)
{
  return GetIsolation(lepton);
}


// subclasses
//============
// charged hadron isolation
class ChargedHadronIsolationFunctor : public IsolationFunctor
{
public:
	ChargedHadronIsolationFunctor();

private:
	const double GetIsolation(const pat::Electron& lepton);
	const double GetIsolation(const pat::Muon& lepton);
	const double GetIsolation(const pat::Tau& lepton);
};

// neutral hadron isolation
class NeutralHadronIsolationFunctor : public IsolationFunctor
{
public:
	NeutralHadronIsolationFunctor();
	
private:
	const double GetIsolation(const pat::Electron& lepton);
	const double GetIsolation(const pat::Muon& lepton);
	const double GetIsolation(const pat::Tau& lepton);
};

// photon isolation
class PhotonIsolationFunctor : public IsolationFunctor
{
public:
	PhotonIsolationFunctor();
	
private:
	const double GetIsolation(const pat::Electron& lepton);
	const double GetIsolation(const pat::Muon& lepton);
	const double GetIsolation(const pat::Tau& lepton);
};

#endif /* ISOLATIONFUNCTOR_H_ */
