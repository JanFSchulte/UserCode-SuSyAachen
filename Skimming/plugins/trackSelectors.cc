#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SingleElementCollectionSelector.h"

#include "DataFormats/TrackReco/interface/Track.h"

// define your producer name
typedef SingleObjectSelector< std::vector<reco::Track>,
                   StringCutObjectSelector<reco::Track>
                   > GenericTrackSelector;

// declare the module as plugin
DEFINE_FWK_MODULE( GenericTrackSelector );
