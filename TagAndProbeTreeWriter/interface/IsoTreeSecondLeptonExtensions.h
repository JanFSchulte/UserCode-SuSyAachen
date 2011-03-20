#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "TTree.h"
#include <DataFormats/Math/interface/deltaR.h>

#include "SuSyAachen/TagAndProbeTreeWriter/interface/LeptonKindFunctor.h"

class IsoTreeSecondLeptonExtensions
{
public:
	IsoTreeSecondLeptonExtensions(){debug_=false;};
	void init(const edm::ParameterSet& iConfig, TTree& tree);
	template< typename T > void fill(TTree& tree, const T& lepton, const edm::Event& iEvent);
private:
	int trueDecayMode( const reco::Candidate& p);
	template< typename T, typename U > int getDileptonCount(const T& lepton, const std::vector<U>& secondLeptons, const bool oppositeSign, const int leptonKind);

    bool debug_;
    double minDeltaR_;
    edm::InputTag electronSrc_;
    edm::InputTag muonSrc_;
    edm::InputTag tauSrc_;

    LeptonKindFunctor fctLeptonKind_;

    int nElectronOS;
    int nElectronSS;
    int nMuonOS;
    int nMuonSS;
    int nTauOS;
    int nTauSS;
    int nPromptElectronOS;
    int nPromptElectronSS;
    int nPromptMuonOS;
    int nPromptMuonSS;
    int nPromptTauOS;
    int nPromptTauSS;
};

void IsoTreeSecondLeptonExtensions::init(const edm::ParameterSet& iConfig, TTree& tree)
{
    if(debug_) std::cout << "init ";
    nElectronOS = 0;
    nElectronSS = 0;
    nMuonOS = 0;
    nMuonSS = 0;
    nTauOS = 0;
    nTauSS = 0;
    tree.Branch("nElectronOS",&nElectronOS,"nElectronOS/I");
    tree.Branch("nElectronSS",&nElectronSS,"nElectronSS/I");
    tree.Branch("nMuonOS",&nMuonOS,"nMuonOS/I");
    tree.Branch("nMuonSS",&nMuonSS,"nMuonSS/I");
    tree.Branch("nTauOS",&nTauOS,"nTauOS/I");
    tree.Branch("nTauSS",&nTauSS,"nTauSS/I");
    tree.Branch("nPromptElectronOS",&nPromptElectronOS,"nPromptElectronOS/I");
    tree.Branch("nPromptElectronSS",&nPromptElectronSS,"nPromptElectronSS/I");
    tree.Branch("nPromptMuonOS",&nPromptMuonOS,"nPromptMuonOS/I");
    tree.Branch("nPromptMuonSS",&nPromptMuonSS,"nPromptMuonSS/I");
    tree.Branch("nPromptTauOS",&nPromptTauOS,"nPromptTauOS/I");
    tree.Branch("nPromptTauSS",&nPromptTauSS,"nPromptTauSS/I");

    if(debug_) std::cout << " InputTags ";
    electronSrc_ = iConfig.getParameter<edm::InputTag> ("secondLeptonElectronSrc");
    muonSrc_ = iConfig.getParameter<edm::InputTag> ("secondLeptonMuonSrc");
    tauSrc_ = iConfig.getParameter<edm::InputTag> ("secondLeptonTauSrc");
	minDeltaR_ = iConfig.getParameter<double> ("secondLeptonMinDeltaR");
    if(debug_) std::cout << "Done!"<< std::endl;
}

template< typename T >
void IsoTreeSecondLeptonExtensions::fill(TTree& tree, const T& lepton, const edm::Event& iEvent)
{
	if(debug_) std::cout << "analyze ";
	edm::Handle< std::vector<pat::Electron> > electrons;
	iEvent.getByLabel(electronSrc_, electrons);
	edm::Handle< std::vector<pat::Muon> > muons;
	iEvent.getByLabel(muonSrc_, muons);
	edm::Handle< std::vector<pat::Tau> > taus;
	iEvent.getByLabel(tauSrc_, taus);
	nElectronOS = getDileptonCount<T, pat::Electron>(lepton, *electrons, true, -1);
	nElectronSS = getDileptonCount<T, pat::Electron>(lepton, *electrons, false, -1);
	nMuonOS = getDileptonCount<T, pat::Muon>(lepton, *muons, true, -1);
	nMuonSS = getDileptonCount<T, pat::Muon>(lepton, *muons, false, -1);
	nTauOS = getDileptonCount<T, pat::Tau>(lepton, *taus, true, -1);
	nTauSS = getDileptonCount<T, pat::Tau>(lepton, *taus, false, -1);
	nPromptElectronOS = getDileptonCount<T, pat::Electron>(lepton, *electrons, true, 3);
	nPromptElectronSS = getDileptonCount<T, pat::Electron>(lepton, *electrons, false, 3);
	nPromptMuonOS = getDileptonCount<T, pat::Muon>(lepton, *muons, true, 3);
	nPromptMuonSS = getDileptonCount<T, pat::Muon>(lepton, *muons, false, 3);
	nPromptTauOS = getDileptonCount<T, pat::Tau>(lepton, *taus, true, 3);
	nPromptTauSS = getDileptonCount<T, pat::Tau>(lepton, *taus, false, 3);

	if(debug_) std::cout << "e:"<<nElectronOS<<","<<nElectronSS
			<< "mu:"<<nMuonOS<<","<<nMuonSS
			<< "tau:"<<nTauOS<<","<<nTauSS;

//    tree.Fill();
    if(debug_) std::cout << "Done!"<< std::endl;
}

template< typename T, typename U >
int IsoTreeSecondLeptonExtensions::getDileptonCount(const T& lepton, const std::vector<U>& secondLeptons, const bool oppositeSign, const int leptonKind)
{
	int result = 0;
	for(typename std::vector<U>::const_iterator it = secondLeptons.begin(); it != secondLeptons.end(); ++it ){
		if(deltaR<T,U>(lepton, *it) > minDeltaR_){
			if( ( ( oppositeSign && lepton.charge()* (*it).charge() < 0) ||
			      (!oppositeSign && lepton.charge()* (*it).charge() > 0) )
			    && (leptonKind == fctLeptonKind_(*it) || leptonKind == -1) )
				result++;
		}else{
			if(debug_) std::cout << "same("<<deltaR<T,U>(lepton, *it)<<":"<<lepton.p4()<<"=="<<(*it).p4()<<") ";
		}
	}
	return result;
}
