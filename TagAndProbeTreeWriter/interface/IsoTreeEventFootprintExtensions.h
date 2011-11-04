#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TTree.h"

class IsoTreeEventFootprintExtensions
{
 public:
  IsoTreeEventFootprintExtensions(){debug_=false;};
  void init(const edm::ParameterSet& iConfig, TTree& tree);
  void fill(const edm::Event& iEvent);
 private:
    bool debug_;

    int runNr_;
    int lumiSec_;
    int eventNr_;
};

void IsoTreeEventFootprintExtensions::init(const edm::ParameterSet& iConfig, TTree& tree)
{
  if(debug_) std::cout << "init EventFootprint";
  runNr_ = 0;
  lumiSec_ = 0;
  eventNr_ = 0;

  tree.Branch("runNr",&runNr_,"runNr/I");
  tree.Branch("lumiSec",&lumiSec_,"lumiSec/I");
  tree.Branch("eventNr",&eventNr_,"eventNr/I");

  if(debug_) std::cout << "Done!"<< std::endl;
}

void IsoTreeEventFootprintExtensions::fill(const edm::Event& iEvent)
{
  if(debug_) std::cout << "analyze EventFootprint";
  runNr_ = iEvent.id().run();
  lumiSec_ = iEvent.id().luminosityBlock();
  eventNr_ = iEvent.id().event();

  if(debug_) std::cout << "Done!"<< std::endl;
}

