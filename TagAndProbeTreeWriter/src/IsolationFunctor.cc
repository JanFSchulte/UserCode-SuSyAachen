#include "SuSyAachen/TagAndProbeTreeWriter/interface/IsolationFunctor.h"

// #include "DataFormats/Math/interface/LorentzVector.h"
// #include "TVector3.h"

#include <iostream>


// base functor class
//====================
IsolationFunctor::IsolationFunctor()
{
}

const double IsolationFunctor::GetIsolation(const pat::Electron& lepton)
{
  double value = lepton.chargedHadronIso()+lepton.photonIso()+lepton.neutralHadronIso();
  return value;
}

const double IsolationFunctor::GetIsolation(const pat::Muon& lepton)
{
  double value = lepton.chargedHadronIso()+lepton.photonIso()+lepton.neutralHadronIso();
  return value;
}

const double IsolationFunctor::GetIsolation(const pat::Tau& lepton)
{
  double value = lepton.chargedHadronIso()+lepton.photonIso()+lepton.neutralHadronIso();
  return value;
}

const double IsolationFunctor::GetIsolation(const reco::Candidate& lepton)
{
  double value = -1.0;
  return value;
}


// subclasses
//============
// charged hadron isolation
ChargedHadronIsolationFunctor::ChargedHadronIsolationFunctor()
{
}

const double ChargedHadronIsolationFunctor::GetIsolation(const pat::Electron& lepton)
{
  double value = lepton.chargedHadronIso();
  return value;
}

const double ChargedHadronIsolationFunctor::GetIsolation(const pat::Muon& lepton)
{
  double value = lepton.chargedHadronIso();
  return value;
}

const double ChargedHadronIsolationFunctor::GetIsolation(const pat::Tau& lepton)
{
  double value = lepton.chargedHadronIso();
  return value;
}

// neutral hadron isolation
NeutralHadronIsolationFunctor::NeutralHadronIsolationFunctor()
{
}

const double NeutralHadronIsolationFunctor::GetIsolation(const pat::Electron& lepton)
{
  double value = lepton.neutralHadronIso();
  return value;
}

const double NeutralHadronIsolationFunctor::GetIsolation(const pat::Muon& lepton)
{
  double value = lepton.neutralHadronIso();
  return value;
}

const double NeutralHadronIsolationFunctor::GetIsolation(const pat::Tau& lepton)
{
  double value = lepton.neutralHadronIso();
  return value;
}

// photon isolation
PhotonIsolationFunctor::PhotonIsolationFunctor()
{
}

const double PhotonIsolationFunctor::GetIsolation(const pat::Electron& lepton)
{
  double value = lepton.photonIso();
  return value;
}

const double PhotonIsolationFunctor::GetIsolation(const pat::Muon& lepton)
{
  double value = lepton.photonIso();
  return value;
}

const double PhotonIsolationFunctor::GetIsolation(const pat::Tau& lepton)
{
  double value = lepton.photonIso();
  return value;
}
