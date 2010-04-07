#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/UtilAlgos/interface/ObjectCountFilter.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
//#include "PhysicsTools/UtilAlgos/interface/MinNumberSelector.h"
//#include "PhysicsTools/UtilAlgos/interface/MaxNumberSelector.h"
//#include "PhysicsTools/UtilAlgos/interface/AndSelector.h"
#include "CommonTools/UtilAlgos/interface/AnySelector.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

// define your producer name
typedef ObjectCountFilter< pat::JetCollection, 
			   StringCutObjectSelector<pat::Jet> 
			   > PATJetCountFilter;

typedef ObjectCountFilter< edm::View<reco::Candidate>, 
			   StringCutObjectSelector<reco::Candidate> 
			   > CandViewJetCountFilter;


// declare the module as plugin
DEFINE_FWK_MODULE( PATJetCountFilter );
DEFINE_FWK_MODULE( CandViewJetCountFilter );
