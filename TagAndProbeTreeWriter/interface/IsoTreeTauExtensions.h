#include "DataFormats/PatCandidates/interface/Tau.h"
#include "TTree.h"

class IsoTreeTauExtensions
{
public:
	IsoTreeTauExtensions(){debug_=false;};
	void init(const edm::ParameterSet& iConfig, TTree& tree);
	void fill(TTree& tree, const pat::Tau& tau);
private:
	int trueDecayMode( const reco::Candidate& p);

    bool debug_;
    int decayMode;
    float invMass;
    float mvaElectron;
    float electronPreIDOutput;

    float isolationPFGammaCandsEtSum;
    float isolationPFChargedHadrCandsPtSum;
    float neutralHadronIso;
    float photonIso;
    float chargedHadronIso;
    float hcalIso;
    float ecalIso;
    float trackIso;

    int decayModeTruth;
    float ptTruth;
    float etaTruth;
    float phiTruth;

    std::map<std::string, float> discriminators_;
};

void IsoTreeTauExtensions::init(const edm::ParameterSet& iConfig, TTree& tree)
{
    if(debug_) std::cout << "init ";
    tree.Branch("decayMode",&decayMode,"decayMode/I");
    tree.Branch("invMass",&invMass,"invMass/F");
    tree.Branch("mvaElectron",&mvaElectron,"mvaElectron/F");
    tree.Branch("electronPreIDOutput",&electronPreIDOutput,"electronPreIDOutput/F");
    tree.Branch("isolationPFGammaCandsEtSum",&isolationPFGammaCandsEtSum,"isolationPFGammaCandsEtSum/F");
    tree.Branch("isolationPFChargedHadrCandsPtSum",&isolationPFChargedHadrCandsPtSum,"isolationPFChargedHadrCandsPtSum/F");
    tree.Branch("neutralHadronIso",&neutralHadronIso,"neutralHadronIso/F");
    tree.Branch("photonIso",&photonIso,"photonIso/F");
    tree.Branch("chargedHadronIso",&chargedHadronIso,"chargedHadronIso/F");
    tree.Branch("hcalIso",&chargedHadronIso,"chargedHadronIso/F");
    tree.Branch("ecalIso",&ecalIso,"ecalIso/F");
    tree.Branch("trackIso",&trackIso,"trackIso/F");
    tree.Branch("electronPreIDOutput",&electronPreIDOutput,"electronPreIDOutput/F");

    tree.Branch("decayModeTruth",&decayModeTruth,"decayModeTruth/I");
    tree.Branch("ptTruth",&ptTruth,"ptTruth/F");
    tree.Branch("etaTruth",&etaTruth,"etaTruth/F");
    tree.Branch("phiTruth",&phiTruth,"phiTruth/F");

    if(debug_) std::cout << "Done!"<< std::endl;
}

void IsoTreeTauExtensions::fill(TTree& tree, const pat::Tau& tau)
{
	if(debug_) std::cout << "analyze ";
    decayMode = tau.decayMode();
    reco::Particle::LorentzVector pSignal(0.,0.,0.,0.);
    reco::PFCandidateRefVector signalPFCands = tau.signalPFCands();
    for(reco::PFCandidateRefVector::const_iterator it = signalPFCands.begin();
    		it != signalPFCands.end(); ++it){
    	pSignal += (*it)->p4();
    }
    invMass = pSignal.M();
    mvaElectron = -2.;
    if(tau.leadPFChargedHadrCand().isNonnull()) mvaElectron = tau.leadPFChargedHadrCand()->mva_e_pi();
    electronPreIDOutput = tau.electronPreIDOutput();
    isolationPFGammaCandsEtSum=tau.isolationPFGammaCandsEtSum();
    isolationPFChargedHadrCandsPtSum=tau.isolationPFChargedHadrCandsPtSum();
    neutralHadronIso=tau.neutralHadronIso();
    photonIso=tau.photonIso();
    chargedHadronIso=tau.chargedHadronIso();
    hcalIso=tau.hcalIso();
    ecalIso=tau.ecalIso();
    trackIso=tau.trackIso();

    decayModeTruth = -20;
    ptTruth = -1.;
    etaTruth = -10.;
    phiTruth = -1.;
    if(debug_) std::cout << "dm("<<decayMode;
    if(tau.genJet() != NULL){
    	assert(tau.genLepton()); //broken genLepton found while looking for true decay mode
    	decayModeTruth = trueDecayMode(*dynamic_cast<const reco::Candidate*>(tau.genLepton()));
    	if(debug_) std::cout <<","<<decayModeTruth;
    	ptTruth = tau.genJet()->pt();
    	etaTruth = tau.genJet()->eta();
    	phiTruth = tau.genJet()->phi();
    }
    if(debug_) std::cout <<")";

	std::vector< pat::Tau::IdPair> tauIds = tau.tauIDs();
	for(std::vector< pat::Tau::IdPair>::iterator it = tauIds.begin(); it != tauIds.end(); ++it){
		bool init = discriminators_.find((*it).first) == discriminators_.end();
		discriminators_[(*it).first] = (*it).second;
		if(init){
			tree.Branch((*it).first.c_str(), &discriminators_[(*it).first], ((*it).first+std::string("/F")).c_str());
			if(debug_) std::cout <<"adding: " <<(*it).first<< " = "<<(*it).second<<", ";
		}
	}
	tree.Fill();
    if(debug_) std::cout << "Done!"<< std::endl;
}

int IsoTreeTauExtensions::trueDecayMode( const reco::Candidate& p)
{
	int result = 0;
	if(p.pdgId() == 111)
		result = 1;
	else if(abs(p.pdgId()) == 211)
		result = 5;
	else if(abs(p.pdgId()) == 12 || abs(p.pdgId()) == 14 || abs(p.pdgId()) == 16)
		result = 0;
	else if(p.status() == 1)
		result = -100;
	else{
		if(abs(p.pdgId()) == 15)
			result = -5;
		if(debug_) std::cout <<"(";
		for(reco::Candidate::const_iterator it = p.begin();
				it != p.end(); ++it){
			//assert(it);//broken candidate found while looking for true decay mode!
			int dm = trueDecayMode(*it);
			result += dm;
			if(debug_) std::cout << dm <<"+";
		}
		if(debug_) std::cout <<")";
	}
	return result;
}

