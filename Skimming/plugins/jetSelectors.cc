#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/UtilAlgos/interface/ObjectCountFilter.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/AnySelector.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "SuSyAachen/Skimming/interface/jetIDSelector.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"

typedef ObjectSelector< jetIDSelector <double, pat::JetCollection, std::vector<const pat::Jet *> > > PATPFJetIDSelector;
// declare the module as plugin
DEFINE_FWK_MODULE( PATPFJetIDSelector );


