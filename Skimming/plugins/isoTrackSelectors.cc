#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/UtilAlgos/interface/AnySelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "SuSyAachen/Skimming/interface/isoTrackSelector.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// define your producer name

typedef ObjectSelector< isoTrackSelector<double, pat::IsolatedTrackCollection, std::vector<const int*> >, std::vector<int> > PATIsoTrackSelector;


// declare the module as plugin
DEFINE_FWK_MODULE( PATIsoTrackSelector );


