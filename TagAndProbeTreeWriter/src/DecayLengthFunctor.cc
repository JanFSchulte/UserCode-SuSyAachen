#include "SuSyAachen/TagAndProbeTreeWriter/interface/DecayLengthFunctor.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TVector3.h"

#include <iostream>

DecayLengthFunctor::DecayLengthFunctor()
{
}

const double DecayLengthFunctor::operator()(const pat::Electron& e1, const pat::Electron& e2,
					    const reco::Vertex& primaryVtx, const edm::EventSetup& iSetup)
{
  return GetDecayLength(e1, e2, primaryVtx, iSetup);
}

const double DecayLengthFunctor::GetDecayLength(const pat::Electron& e1, const pat::Electron& e2,
						const reco::Vertex& primaryVtx, const edm::EventSetup& iSetup)
{
  double value = -10.0;
  std::vector<TransientVertex> pvs;
  math::XYZTLorentzVectorD jpsi = e1.p4() + e2.p4();
  // to be consistent with TnP Analyzer
  // math::XYZTLorentzVectorD e2P4 = math::XYZTLorentzVectorD(e2.px(), e2.py(), e2.pz(), e2.p());
  // math::XYZTLorentzVectorD jpsi = e1.p4() + e2P4;

  // for debugging reasons slightly different return values if quality cuts fail
  // --> removing sanity checks, can do all those on the final tree
  /*
  if (jpsi.M() > 3.4)
    return -10.1;
  if (jpsi.M() < 2.8)
    return -10.2;
  if (e1.charge() + e2.charge()!=0)
    return -10.3;
  */
    
  // what is this for?
  /*
  if (  it->trackRef().id() == pvtracks.id() && it2->trackRef().id() == pvtracks.id()) { 
    TrackCollection elecLess;
    elecLess.reserve(pvtracks->size());
    for (size_t i = 0, n = pvtracks->size(); i < n; ++i) {
      if (i == it->trackRef().key()) continue;
      if (i == it2->trackRef().key()) continue;
      elecLess.push_back((*pvtracks)[i]);
    }
    pvs = revertex.makeVertices(elecLess, *pvbeamspot, iSetup) ;
    if (!pvs.empty()) {
      Vertex elecLessPV = Vertex(pvs.front());
      thePrimaryV = elecLessPV;
    }
  }
  */

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter;

  std::vector<reco::TransientTrack> t_tks;
  t_tks.push_back(theTTBuilder->build(*e1.gsfTrack()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
  t_tks.push_back(theTTBuilder->build(*e2.gsfTrack())); // otherwise the vertex will have transient refs inside.
  TransientVertex myVertex = vtxFitter.vertex(t_tks);
  if (myVertex.isValid()) {
    TVector3 vtx;
    TVector3 pvtx;
    VertexDistanceXY vdistXY;
	  
    vtx.SetXYZ(myVertex.position().x(),myVertex.position().y(),0);
    TVector3 pperp(jpsi.px(), jpsi.py(), 0);
    AlgebraicVector vpperp(3);
    vpperp[0] = pperp.x();
    vpperp[1] = pperp.y();
    vpperp[2] = 0.;
	  
    // lifetime using PV
    pvtx.SetXYZ(primaryVtx.position().x(),primaryVtx.position().y(),0);
    TVector3 vdiff = vtx - pvtx;
    double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
    Measurement1D distXY = vdistXY.distance(reco::Vertex(myVertex), primaryVtx);
    double ctauPV = distXY.value()*cosAlpha*3.09688/pperp.Perp();
    // GlobalError v1e = (reco::Vertex(myVertex)).error();
    // GlobalError v2e = primaryVtx.error();
    // AlgebraicSymMatrix vXYe = v1e.matrix()+ v2e.matrix();
    // double ctauErrPV = sqrt(vXYe.similarity(vpperp))*3.09688/(pperp.Perp2());

    //std::cout << "DecayLength: " << ctauPV << std::endl;
    value = ctauPV;
  }

  return value;

  //GENERATOR LEVEL
  /*
  Handle<GenParticleCollection> GEN;
  iEvent.getByLabel("genParticles",GEN);
  math::XYZPoint momVx(0,0,0);
  math::XYZPoint dauVx(0,0,0);
  bool isMom=false;
  bool isDau=false;
  for(GenParticleCollection::const_iterator it = GEN->begin(); it != GEN->end(); ++it) {
    for (uint uu=0;uu<it->daughterRefVector().size();uu++){
      if ((*it->daughterRefVector()[uu]).pdgId()==443 && (it->pdgId()%10000)!=443) {
	isMom=true;
	momVx=it->vertex();
      }
    }
    
    if (it->pdgId()==443 && it->status()==2) { 
      dauVx=it->vertex();
      isDau=true;
    }
  }
  if (isMom && isDau){
    math::XYZPoint DeltaVx(momVx.x()-dauVx.x(), momVx.y()-dauVx.y(),momVx.z()-dauVx.z());
    DeltaZ->Fill(DeltaVx.z());
    DeltaRho->Fill(DeltaVx.Rho());
    DeltaR->Fill(DeltaVx.R());
  }
  */
  //END GENERATOR

  /*
  //RECONSTRUCTION
  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter;

  //VERTICES
  double value = 0.0;
  return value;
  Handle<VertexCollection> priVtxs;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS", priVtxs);

  //BEAM SPOT
  Handle<BeamSpot> theBeamSpot;
  iEvent.getByLabel("offlineBeamSpot",theBeamSpot);

  BeamSpot bs = *theBeamSpot;
  Vertex theBeamSpotV(bs.position(), bs.covariance3D());
  Vertex thePrimaryV=( priVtxs->begin() != priVtxs->end() )? 
    Vertex(*(priVtxs->begin())): Vertex(bs.position(), bs.covariance3D());

  VertexReProducer revertex(priVtxs, iEvent);
  Handle<TrackCollection> pvtracks;   
  iEvent.getByLabel(revertex.inputTracks(),   pvtracks);

  Handle<BeamSpot>        pvbeamspot; 
  iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);

  if (pvbeamspot.id() != theBeamSpot.id()) 
    cout << "The BeamSpot used for PV reco is not the same used in this analyzer."<<endl;



  if (conf_.getParameter<bool>("isMuon")){
    
    Handle< View<reco::Muon> > muons;
    iEvent.getByLabel("muons",muons);
    
    
    // JPsi candidates only from muons
    for(View<reco::Muon>::const_iterator it = muons->begin(), itend = muons->end(); it != itend; ++it){
      // both must pass low quality
      if (it->track().isNull()) continue;
      for(View<reco::Muon>::const_iterator it2 = it+1; it2 != itend;++it2){
	if (it2->track().isNull()) continue;
	
	pat::CompositeCandidate myCand;
	vector<TransientVertex> pvs;
	
	
	// ---- define and set candidate's 4momentum  ----  
	LorentzVector jpsi = it->p4() + it2->p4();
	if (jpsi.M()>3.4) continue;
	if (jpsi.M()<2.8) continue;
	if (it->charge()+it2->charge()!=0)continue;

	if (it->track().isNonnull() && it2->track().isNonnull()) {
	  
	  
	  if (  it->track().id() == pvtracks.id() && it2->track().id() == pvtracks.id()) { 
	    
	    TrackCollection muonLess;
	    muonLess.reserve(pvtracks->size());
	    for (size_t i = 0, n = pvtracks->size(); i < n; ++i) {
	      if (i == it->track().key()) continue;
	      if (i == it2->track().key()) continue;
	      muonLess.push_back((*pvtracks)[i]);
	    }
	    pvs = revertex.makeVertices(muonLess, *pvbeamspot, iSetup) ;
	    if (!pvs.empty()) {
	      Vertex muonLessPV = Vertex(pvs.front());
	      thePrimaryV = muonLessPV;
	    }
	  }
	}
	vector<TransientTrack> t_tks;
	t_tks.push_back(theTTBuilder->build(*it->track()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
	t_tks.push_back(theTTBuilder->build(*it2->track())); // otherwise the vertex will have transient refs inside.
	TransientVertex myVertex = vtxFitter.vertex(t_tks);
	if (myVertex.isValid()) {
	  
	  
	  TVector3 vtx;
	  TVector3 pvtx;
	  VertexDistanceXY vdistXY;
	
	  vtx.SetXYZ(myVertex.position().x(),myVertex.position().y(),0);
	  TVector3 pperp(jpsi.px(), jpsi.py(), 0);
	  AlgebraicVector vpperp(3);
	  vpperp[0] = pperp.x();
	  vpperp[1] = pperp.y();
	  vpperp[2] = 0.;
	  
	  
	  
        


	  // lifetime using PV
	  pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
	  TVector3 vdiff = vtx - pvtx;
	  double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
	  Measurement1D distXY = vdistXY.distance(Vertex(myVertex), thePrimaryV);
	  double ctauPV = distXY.value()*cosAlpha*3.09688/pperp.Perp();
	  GlobalError v1e = (Vertex(myVertex)).error();
	  GlobalError v2e = thePrimaryV.error();
	  AlgebraicSymMatrix vXYe = v1e.matrix()+ v2e.matrix();
	  double ctauErrPV = sqrt(vXYe.similarity(vpperp))*3.09688/(pperp.Perp2());
	  PseudoDL->Fill(ctauPV);
	  PseudoDLSign->Fill(ctauPV/ctauErrPV);
	  //	myCand.addUserFloat("ppdlPV",ctauPV);
	  //	  cout<<"PSEUDO MUON "<< ctauPV <<"+/- "<<ctauErrPV<<endl;
	  
	}
      }
    }
  }
  //END MUONS
  //START ELECTRONS
  if (!(conf_.getParameter<bool>("isMuon"))){

    Handle< View<reco::PFCandidate> > electrons;
    iEvent.getByLabel("particleFlow","electrons",electrons);
    
    // JPsi candidates only from muons
    for(View<reco::PFCandidate>::const_iterator it = electrons->begin(), itend = electrons->end(); it != itend; ++it){
      // both must pass low quality
      if (it->trackRef().isNull()) continue;
      for(View<reco::PFCandidate>::const_iterator it2 = it+1; it2 != itend;++it2){
	if (it2->trackRef().isNull()) continue;
	
	pat::CompositeCandidate myCand;
	vector<TransientVertex> pvs;
	

	LorentzVector jpsi = it->p4() + it2->p4();

	if (jpsi.M()>3.4) continue;
	if (jpsi.M()<2.8) continue;
	if (it->charge()+it2->charge()!=0)continue;

	if (it->trackRef().isNonnull() && it2->trackRef().isNonnull()) {
	  //	  cout<<"PL "<<it->trackRef()->pt()<<" "<<it->gsfTrackRef()->pt()<<endl;
	  if (  it->trackRef().id() == pvtracks.id() && it2->trackRef().id() == pvtracks.id()) { 
	    
	    TrackCollection elecLess;
	    elecLess.reserve(pvtracks->size());
	    for (size_t i = 0, n = pvtracks->size(); i < n; ++i) {
	      if (i == it->trackRef().key()) continue;
	      if (i == it2->trackRef().key()) continue;
	      elecLess.push_back((*pvtracks)[i]);
	    }
	    pvs = revertex.makeVertices(elecLess, *pvbeamspot, iSetup) ;
	    if (!pvs.empty()) {
	      Vertex elecLessPV = Vertex(pvs.front());
	      thePrimaryV = elecLessPV;
	    }
	  }
	}
	vector<TransientTrack> t_tks;
	t_tks.push_back(theTTBuilder->build(*it->trackRef()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
	t_tks.push_back(theTTBuilder->build(*it2->trackRef())); // otherwise the vertex will have transient refs inside.
	TransientVertex myVertex = vtxFitter.vertex(t_tks);
	if (myVertex.isValid()) {
	  
	  
	  TVector3 vtx;
	  TVector3 pvtx;
	  VertexDistanceXY vdistXY;
	  
	  vtx.SetXYZ(myVertex.position().x(),myVertex.position().y(),0);
	  TVector3 pperp(jpsi.px(), jpsi.py(), 0);
	  AlgebraicVector vpperp(3);
	  vpperp[0] = pperp.x();
	  vpperp[1] = pperp.y();
	  vpperp[2] = 0.;
	  
	  // lifetime using PV
	  pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
	  TVector3 vdiff = vtx - pvtx;
	  double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
	  Measurement1D distXY = vdistXY.distance(Vertex(myVertex), thePrimaryV);
	  double ctauPV = distXY.value()*cosAlpha*3.09688/pperp.Perp();
	  GlobalError v1e = (Vertex(myVertex)).error();
	  GlobalError v2e = thePrimaryV.error();
	  AlgebraicSymMatrix vXYe = v1e.matrix()+ v2e.matrix();
	  double ctauErrPV = sqrt(vXYe.similarity(vpperp))*3.09688/(pperp.Perp2());
	  PseudoDL->Fill(ctauPV);
	  PseudoDLSign->Fill(ctauPV/ctauErrPV);
	  //	  cout<<"PSEUDO ELEC "<< ctauPV <<"+/- "<<ctauErrPV<<endl;
	  
	}
	//	}
      }
    }
  }
  */
}



