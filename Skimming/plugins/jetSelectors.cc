#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/UtilAlgos/interface/ObjectCountFilter.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/AnySelector.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "SuSyAachen/Skimming/interface/jetIDSelector.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"

// define your producer name
typedef ObjectCountFilter< pat::JetCollection, 
			   StringCutObjectSelector<pat::Jet> >::type PATJetCountFilter;

typedef ObjectCountFilter< edm::View<reco::Candidate>, 
			   StringCutObjectSelector<reco::Candidate> >::type CandViewJetCountFilter;

typedef ObjectSelector< jetIDSelector <double, pat::JetCollection, std::vector<const pat::Jet *> > > PATPFJetIDSelector;
typedef SingleObjectSelector< pat::JetCollection, JetIDSelectionFunctor, pat::JetCollection > PATJetIDSelector;

// declare the module as plugin
DEFINE_FWK_MODULE( PATJetCountFilter );
DEFINE_FWK_MODULE( CandViewJetCountFilter );
DEFINE_FWK_MODULE( PATPFJetIDSelector );
DEFINE_FWK_MODULE( PATJetIDSelector );

