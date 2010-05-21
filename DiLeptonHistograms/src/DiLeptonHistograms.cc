/** \class DiLeptonHistograms
 *
 *  
 *  This class is an EDAnalyzer for PAT 
 *  Layer 0 and Layer 1 output
 *
 *  $Date: 2010/05/21 15:26:53 $
 *  $Revision: 1.13 $ for CMSSW 3_6_X
 *
 *  \author: Niklas Mohr -- niklas.mohr@cern.ch
 *  
 */

#include "SuSyAachen/DiLeptonHistograms/interface/DiLeptonHistograms.h"


//Constructor
DiLeptonHistograms::DiLeptonHistograms(const edm::ParameterSet &iConfig)
{
    //now do what ever initialization is needed
    //Debug flag
    debug             = iConfig.getUntrackedParameter<bool>   ("debug");

    //Monte carlo information
    mcInfo            = iConfig.getUntrackedParameter<bool>   ("mcInfo",false);

    //tree information for unbinned fit
    treeInfo          = iConfig.getUntrackedParameter<bool>   ("treeInfo",false);
    effInfo           = iConfig.getUntrackedParameter<bool>   ("effInfo",false);

    // how many tauDiscriminators to take
    maxTauDiscriminators_ = 20; //TODO read from config
    //Input collections
    mcSrc             = iConfig.getParameter<edm::InputTag> ("mcSource");
    beamSpotSrc        = iConfig.getParameter<edm::InputTag> ("beamSpotSource");
    primaryVertexSrc  = iConfig.getParameter<edm::InputTag> ("primaryVertexSource");
    muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
    electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
    tauSrc            = iConfig.getParameter<edm::InputTag> ("tauSource");
    metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
    jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
    trackSrc          = iConfig.getParameter<edm::InputTag> ("trackSource");
 
    externalWeight   = iConfig.getUntrackedParameter<double>   ("external_Weight",1.);
    Signal_Analysis  = iConfig.getUntrackedParameter<bool>   ("Signal_Analysis",false);

    //Cuts
    cut_MuonPt     = iConfig.getUntrackedParameter<double> ("acc_MuonPt");
    cut_MuonEta    = iConfig.getUntrackedParameter<double> ("acc_MuonEta");
   
    cut_ElectronPt  = iConfig.getUntrackedParameter<double> ("acc_ElectronPt");
    cut_ElectronEta  = iConfig.getUntrackedParameter<double> ("acc_ElectronEta");

    bJetAlgo  = iConfig.getUntrackedParameter<std::string> ("user_bJetAlgo");
    muon_fname = iConfig.getUntrackedParameter<std::string>("muon_fname","SuSyAachen/DiLeptonHistograms/data/efficiency_Muon.data");
    electron_fname = iConfig.getUntrackedParameter<std::string>("electron_fname","SuSyAachen/DiLeptonHistograms/data/efficiency_Electron.data");
    cut_bTagDiscriminator = iConfig.getUntrackedParameter<double> ("user_bTagDiscriminator");
    
    electron_Scale  = iConfig.getUntrackedParameter<double>   ("electron_Scale",0.);
    eff_Electron_Scale  = iConfig.getUntrackedParameter<double>   ("eff_Electron_Scale",0.);
    eff_Muon_Scale  = iConfig.getUntrackedParameter<double>   ("eff_Muon_Scale",0.);
    jet_Scale  = iConfig.getUntrackedParameter<double>   ("jet_Scale",0.);

    //initialize global counters
    numTotEvents = 0;
    numTotElectrons = 0;
    numTotMuons = 0;
    numTotTaus = 0;
    numTotJets = 0;

    const int nHistos=5;

    // Create the root file
    edm::Service<TFileService> theFile;

    // book histograms for multiplicities of leptons and jets
    hLeptonMult = new TH1F * [nHistos];
    hLightLeptonMult = new TH1F * [nHistos];
    hElectronMult = new TH1F * [nHistos];
    hMuonMult = new TH1F * [nHistos];
    hTauMult = new TH1F * [nHistos];
    hJetMult = new TH1F * [nHistos];
    hbJetMult = new TH1F * [nHistos];

    //histograms for lepton isolation cuts
    hElectronIso = new TH1F * [nHistos];
    hElectronTrackIso = new TH1F * [nHistos];
    hElectronEcalIso = new TH1F * [nHistos];
    hElectronHcalIso = new TH1F * [nHistos];
    hElectronPfIso = new TH1F * [nHistos];
    hElectronChargedHadronIso = new TH1F * [nHistos];
    hElectronPhotonIso = new TH1F * [nHistos];
    hElectronNeutralHadronIso = new TH1F * [nHistos];
    hMuonIso = new TH1F * [nHistos];
    hMuonTrackIso = new TH1F * [nHistos];
    hMuonEcalIso = new TH1F * [nHistos];
    hMuonHcalIso = new TH1F * [nHistos];
    hMuonPfIso = new TH1F * [nHistos];
    hMuonChargedHadronIso = new TH1F * [nHistos];
    hMuonPhotonIso = new TH1F * [nHistos];
    hMuonNeutralHadronIso = new TH1F * [nHistos];
    hTauIso = new TH1F * [nHistos];
    hTauTrackIso = new TH1F * [nHistos];
    hTauEcalIso = new TH1F * [nHistos];
    hTauHcalIso = new TH1F * [nHistos];
    hTauPfIso = new TH1F * [nHistos];
    hTauChargedHadronIso = new TH1F * [nHistos];
    hTauPhotonIso = new TH1F * [nHistos];
    hTauNeutralHadronIso = new TH1F * [nHistos];
    
    //histograms for the invariant mass of the leptons
    hInvMSFOS = new TH1F * [nHistos];
    hInvMOFOS = new TH1F * [nHistos];
    hInvMass = new TH1F * [nHistos];
    hDileptonPt = new TH1F * [nHistos];
    hJZB = new TH1F * [nHistos];
    hInvMElectron = new TH1F * [nHistos];
    hInvMElectronSS = new TH1F * [nHistos];
    hInvMMuon = new TH1F * [nHistos];
    hInvMMuonSS = new TH1F * [nHistos];
    hInvMassMC = new TH1F * [nHistos];
    hInvMassZMC = new TH1F * [nHistos];
    hInvMTau = new TH1F * [nHistos];
    hInvMTauSS = new TH1F * [nHistos];
    
    hMuonResolution = new TH1F * [nHistos];
    hElectronResolution = new TH1F * [nHistos];

    //Lepton corr plots
    hDeltaR = new TH1F * [nHistos];
    hDeltaPhi = new TH1F * [nHistos];
    hDeltaPhiMET = new TH1F * [nHistos];
  
    //bbll histos
    hInvMbbllSFOS = new TH1F * [nHistos];
    hInvMbbllOFOS = new TH1F * [nHistos];

    //2D histos
    h2dMETInvMassSFOS = new TH2F * [nHistos];
    h2dHTInvMassSFOS = new TH2F * [nHistos];
    h2dEtJetsInvMassSFOS = new TH2F * [nHistos];
    h2dEtJet4InvMassSFOS = new TH2F * [nHistos];
    h2dIsoInvMassSFOS = new TH2F * [nHistos];

    h2dMETInvMassOFOS = new TH2F * [nHistos];
    h2dHTInvMassOFOS = new TH2F * [nHistos];
    h2dEtJetsInvMassOFOS = new TH2F * [nHistos];
    h2dEtJet4InvMassOFOS = new TH2F * [nHistos];
    h2dIsoInvMassOFOS = new TH2F * [nHistos];

    //muon histograms
    hMuonPt = new TH1F * [nHistos];
    hMuonCharge = new TH1F * [nHistos];
    hMuonSumPt = new TH1F * [nHistos];
    hMuon1Pt = new TH1F * [nHistos];
    hMuon2Pt = new TH1F * [nHistos];
    hMuonEta = new TH1F * [nHistos];
    hMuon1Eta = new TH1F * [nHistos];
    hMuon2Eta = new TH1F * [nHistos];
    hMuonPhi = new TH1F * [nHistos];
    hGenMuonPt = new TH1F * [nHistos];
    hGenMuonEta = new TH1F * [nHistos];
    hMuonChi2 = new TH1F * [nHistos];
    hMuond0Pv = new TH1F * [nHistos];
    hMuond0Bs = new TH1F * [nHistos];
    hMuond0SigPv = new TH1F * [nHistos];
    hMuond0SigBs = new TH1F * [nHistos];
    hMuonnHits = new TH1F * [nHistos];
    hMuonEtaPhi = new TH2F * [nHistos];
    hMuonPfIsod0Pv = new TH2F * [nHistos];
    hMuonPfIsod0Bs = new TH2F * [nHistos];

    //electron histograms
    hElectronPt = new TH1F * [nHistos];
    hElectronCharge = new TH1F * [nHistos];
    hElectronSumPt = new TH1F * [nHistos];
    hElectron1Pt = new TH1F * [nHistos];
    hElectron2Pt = new TH1F * [nHistos];
    hElectronEta = new TH1F * [nHistos];
    hElectron1Eta = new TH1F * [nHistos];
    hElectron2Eta = new TH1F * [nHistos];
    hElectronPhi = new TH1F * [nHistos];
    hElectronLostHits = new TH1F * [nHistos];
    hElectrond0Pv = new TH1F * [nHistos];
    hElectrond0Bs = new TH1F * [nHistos];
    hElectrond0SigPv = new TH1F * [nHistos];
    hElectrond0SigBs = new TH1F * [nHistos];
    hGenElectronPt = new TH1F * [nHistos];
    hGenElectronEta = new TH1F * [nHistos];
    hElectronMva = new TH1F * [nHistos];
    hElectronEoverP = new TH1F * [nHistos];
    hElectronfBrem = new TH1F * [nHistos];
    hElectronHoverE = new TH1F * [nHistos];
    hElectrondeltaPhiIn = new TH1F * [nHistos];
    hElectrondeltaEtaIn = new TH1F * [nHistos];
    hElectrone25Maxoe55 = new TH1F * [nHistos];
    hElectrone15oe55 = new TH1F * [nHistos];
    hElectronsigmaEtaEta = new TH1F * [nHistos];
    hElectronsigmaIetaIeta = new TH1F * [nHistos];
    hElectronsigmaIetaIeta = new TH1F * [nHistos];
    hElectronEtaPhi = new TH2F * [nHistos];
    hElectronPfIsod0Pv = new TH2F * [nHistos];
    hElectronPfIsod0Bs = new TH2F * [nHistos];
    
    //Tau histograms
    hTauPt = new TH1F * [nHistos];
    hTauCharge = new TH1F * [nHistos];
    hTauSumPt = new TH1F * [nHistos];
    hTau1Pt = new TH1F * [nHistos];
    hTau2Pt = new TH1F * [nHistos];
    hTauEta = new TH1F * [nHistos];
    hTau1Eta = new TH1F * [nHistos];
    hTau2Eta = new TH1F * [nHistos];
    hTauPhi = new TH1F * [nHistos];
    hTauEtaPhi = new TH2F * [nHistos];

    hTauDiscriminators = new TH1F * [nHistos];
    hGenTauPt = new TH1F * [nHistos];
    hGenTauVisPt = new TH1F * [nHistos];

    h2dMuonEtaPt = new TH2F * [nHistos];
    h2dMatchedMuonEtaPt = new TH2F * [nHistos];
    h2dGenMuonEtaPt = new TH2F * [nHistos];
    h2dElectronEtaPt = new TH2F * [nHistos];
    h2dMatchedElectronEtaPt = new TH2F * [nHistos];
    h2dGenElectronEtaPt = new TH2F * [nHistos];
    h2dTauEtaPt = new TH2F * [nHistos];
    h2dMatchedTauEtaPt = new TH2F * [nHistos];
  
    //histograms for Missing ET
    hMissingET = new TH1F * [nHistos];
    hMissingETmc =  new TH1F * [nHistos];
    hEtSum = new TH1F * [nHistos];
    halphaT = new TH1F * [nHistos];
    hHT = new TH1F * [nHistos];
    hMHT = new TH1F * [nHistos];
    h2dMETEtSumJets = new TH2F * [nHistos];
    h2dMETHT = new TH2F * [nHistos];
    h2dMETMHT = new TH2F * [nHistos];
    
    hJet1METPhi = new TH1F * [nHistos];
    hJet1EtaMET = new TH2F * [nHistos];

    //histograms for jets
    hJetEt = new TH1F * [nHistos];
    hJetEta = new TH1F * [nHistos];
    hJetPhi = new TH1F * [nHistos];
    hEtJet1 = new TH1F * [nHistos];
    hEtJet2 = new TH1F * [nHistos];
    hEtJet3 = new TH1F * [nHistos];
    hEtJet4 = new TH1F * [nHistos];
    hJet1Eta = new TH1F * [nHistos];
    hJet2Eta = new TH1F * [nHistos];
    hJet3Eta = new TH1F * [nHistos];
    hJet4Eta = new TH1F * [nHistos];
    hJetCHF = new TH1F * [nHistos];
    hJetNHF = new TH1F * [nHistos];
    hJetCEF = new TH1F * [nHistos];
    hJetNEF = new TH1F * [nHistos];
    hJetCM = new TH1F * [nHistos];
    hJetNConst = new TH1F * [nHistos];
    
    h2dJetEtaPhi = new TH2F * [nHistos];
    h2dJetEtaPt = new TH2F * [nHistos];
    h2dMatchedJetEtaPt = new TH2F * [nHistos];
  
    hTrigger = new TH1F * [nHistos];
    hWeight = new TH1F * [nHistos];

    TFileDirectory General = theFile->mkdir( "General" );
    TFileDirectory Effcor = theFile->mkdir( "Efficiency corrected" );
    TFileDirectory Unmatched = theFile->mkdir( "Unmatched" );
    TFileDirectory Promt = theFile->mkdir( "Promt" );
    TFileDirectory Decay = theFile->mkdir( "Decay" );

    //Trees for unbinned maximum likelihood fit
    if (treeInfo){
        TFileDirectory Tree = theFile->mkdir( "Trees" );
        treeOFOS = Tree.make<TTree>("OFOS tree", "OFOS tree"); 
        treeOFOS->Branch("inv",&invMOFOS,"invMOFOS/F");
        treeOFOS->Branch("weight",&invweight,"invweight/F");
        treeSFOS = Tree.make<TTree>("SFOS tree", "SFOS tree"); 
        treeSFOS->Branch("inv",&invMSFOS,"invMSFOS/F");
        treeSFOS->Branch("weight",&invweight,"invweight/F");
        treeElec = Tree.make<TTree>("Electron tree", "Electron tree"); 
        treeElec->Branch("inv",&invMElec,"invMElec/F");
        treeElec->Branch("weight",&invweight,"invweight/F");
        treeMuon = Tree.make<TTree>("Muon tree", "Muon tree"); 
        treeMuon->Branch("inv",&invMMuon,"invMMuon/F");
        treeMuon->Branch("weight",&invweight,"invweight/F");
    }
   
    general = 0;
    effcor = 1;
    unmatched = 2;
    promt = 3;
    decay = 4;
    tauInitialized_ =  new bool[nHistos];
    for(int i =0; i< nHistos; ++i){
      tauInitialized_[i] = false;
    }

    InitHisto(&General,general);
    if (effInfo) InitHisto(&Effcor,effcor);
    if (mcInfo){
        InitHisto(&Unmatched,unmatched);
        InitHisto(&Promt,promt);
        InitHisto(&Decay,decay);
    }
   
    //Read the efficiencies from the files 
    ReadEfficiency();
}

//Initialize all histos including their boundaries
void inline DiLeptonHistograms::InitHisto(TFileDirectory *theFile, const int process)
{
    //Multiplicity plots
    TFileDirectory Multiplicity = theFile->mkdir("Multiplicity"); 
    hLeptonMult[process] = Multiplicity.make<TH1F>( "LeptonMultiplicity", "Multiplicity of electrons + muons + taus", 15, -0.5, 14.5);
    hLightLeptonMult[process] = Multiplicity.make<TH1F>( "LightLeptonMultiplicity", "Multiplicity of electrons + muons", 15, -0.5, 14.5);
    hElectronMult[process] = Multiplicity.make<TH1F>( "ElectronMultiplicity", "Multiplicity of electrons", 10, -0.5, 9.5);
    hMuonMult[process] = Multiplicity.make<TH1F>( "MuonMultiplicity", "Multiplicity of muons", 10, -0.5, 9.5);
    hTauMult[process] = Multiplicity.make<TH1F>( "TauMultiplicity", "Multiplicity of taus", 10, -0.5, 9.5);
    hJetMult[process] = Multiplicity.make<TH1F>( "JetMultiplicity", "Multiplicity of jets", 30, -0.5, 29.5);
    hbJetMult[process] = Multiplicity.make<TH1F>( "bJetMultiplicity", "Multiplicity of b jets", 15, -0.5, 14.5);

        
    TFileDirectory InvMass = theFile->mkdir("Invariant Mass"); 
    //histograms for the invariant mass of the leptons
    hInvMSFOS[process] = InvMass.make<TH1F>( "Invariant mass of SFOS lepton pairs", "Invariant mass of SFOS lepton pairs", 300, 0, 300);
    hInvMOFOS[process] = InvMass.make<TH1F>( "Invariant mass of OFOS lepton pairs", "Invariant mass of OFOS lepton pairs", 300, 0, 300);
    hInvMass[process] = InvMass.make<TH1F>( "Invariant mass of Fsubtracted OS lepton pairs", "Invariant mass of flavor subtracted opposite sign lepton pairs", 300, 0, 300);
    hDileptonPt[process] = InvMass.make<TH1F>( "pt of lepton pairs", "pt of lepton pairs", 500, 0, 500);
    hJZB[process] = InvMass.make<TH1F>( "JZB", "JZB of lepton pairs", 1000, -500, 500);
    hInvMElectron[process] = InvMass.make<TH1F>( "Invariant mass of OS electron pairs", "Invariant mass of opposite sign electron pairs", 300, 0, 300);
    hInvMElectronSS[process] = InvMass.make<TH1F>( "Invariant mass of SS electron pairs", "Invariant mass of same sign electron pairs", 300, 0, 300);
    hInvMMuon[process] = InvMass.make<TH1F>( "Invariant mass of OS muon pairs", "Invariant mass of opposite sign muon pairs", 300, 0, 300);
    hInvMMuonSS[process] = InvMass.make<TH1F>( "Invariant mass of SS muon pairs", "Invariant mass of same sign muon pairs", 300, 0, 300);
    hInvMTau[process] = InvMass.make<TH1F>( "Invariant mass of OS tau pairs", "Invariant mass of opposite sign tau pairs", 300, 0, 300);
    hInvMTauSS[process] = InvMass.make<TH1F>( "Invariant mass of SS tau pairs", "Invariant mass of same sign tau pairs", 300, 0, 300);
    hInvMassMC[process] = InvMass.make<TH1F>( "Invariant mass of signal decays", "Invariant mass of signal decays", 300, 0, 300);
    hInvMassZMC[process] = InvMass.make<TH1F>( "Invariant mass of Z decays", "Invariant mass of Z decays", 300, 0, 300);
 
    hDeltaR[process] = InvMass.make<TH1F>( "delta R", "delta R", 350, 0., 7.);
    hDeltaPhi[process] = InvMass.make<TH1F>( "delta Phi", "delta phi", 350, -3.5, 3.5);
    hDeltaPhiMET[process] = InvMass.make<TH1F>( "delta Phi MET", "delta phi MET", 350, -3.5, 3.5);
    
    //histograms for bbll inv mass
    hInvMbbllSFOS[process] = InvMass.make<TH1F>( "Invariant mass of SFOS bbll", "Invariant mass of SFOS bbll", 1500, 0, 1500);
    hInvMbbllOFOS[process] = InvMass.make<TH1F>( "Invariant mass of OFOS bbll", "Invariant mass of OFOS bbll", 1500, 0, 1500);

    h2dMETInvMassSFOS[process] = InvMass.make<TH2F>( "MET - Invariant mass SFOS", "MET - Invariant mass SFOS", 300, 0., 300., 1000, 0., 1000.0);
    h2dHTInvMassSFOS[process] = InvMass.make<TH2F>( "HT - Invariant mass SFOS", "HT - Invariant mass SFOS", 300, 0., 300., 2000, 0., 2000.0);
    h2dEtJetsInvMassSFOS[process] = InvMass.make<TH2F>( "SumJets - Invariant mass SFOS", "SumJets - Invariant mass SFOS", 300, 0., 300., 2000, 0., 2000.0);
    h2dEtJet4InvMassSFOS[process] = InvMass.make<TH2F>( "PtJet4 - Invariant mass SFOS", "PtJet4 - Invariant mass SFOS", 300, 0., 300., 2000, 0., 1000.0);
    h2dIsoInvMassSFOS[process] = InvMass.make<TH2F>( "Iso - Invariant mass SFOS", "Iso - Invariant mass SFOS", 300, 0., 300., 600, 0., 6.0);

    h2dMETInvMassOFOS[process] = InvMass.make<TH2F>( "MET - Invariant mass OFOS", "MET - Invariant mass OFOS", 300, 0., 300., 1000, 0., 1000.0);
    h2dHTInvMassOFOS[process] = InvMass.make<TH2F>( "HT - Invariant mass OFOS", "HT - Invariant mass OFOS", 300, 0., 300., 2000, 0., 2000.0);
    h2dEtJetsInvMassOFOS[process] = InvMass.make<TH2F>( "SumJets - Invariant mass OFOS", "SumJets - Invariant mass OFOS", 300, 0., 300., 2000, 0., 2000.0);
    h2dEtJet4InvMassOFOS[process] = InvMass.make<TH2F>( "PtJet4 - Invariant mass OFOS", "PtJet4 - Invariant mass OFOS", 300, 0., 300., 2000, 0., 1000.0);
    h2dIsoInvMassOFOS[process] = InvMass.make<TH2F>( "Iso - Invariant mass OFOS", "Iso - Invariant mass OFOS", 300, 0., 300., 600, 0., 6.0);

    TFileDirectory Muons = theFile->mkdir("Muons"); 
    //muon histograms
    hMuonPt[process] = Muons.make<TH1F>( "muon pt", "muon pt", 1000, 0.0, 1000.0);
    hMuonCharge[process] = Muons.make<TH1F>( "muon charge", "muon charge", 8, -2.0, 2.0);
    hMuonSumPt[process] = Muons.make<TH1F>( "sum muon pt", "sum muon pt", 1000, 0.0, 1000.0);
    hMuon1Pt[process] = Muons.make<TH1F>( "muon 1 pt", "pt of first muon", 1000, 0.0, 1000.0);
    hMuon2Pt[process] = Muons.make<TH1F>( "muon 2 pt", "pt of second muon", 1000, 0.0, 1000.0);
    hMuonEta[process] = Muons.make<TH1F>( "muon eta", "muon eta", 250, -2.5, 2.5);
    hMuon1Eta[process] = Muons.make<TH1F>( "muon 1 eta", "muon 1 eta", 250, -2.5, 2.5);
    hMuon2Eta[process] = Muons.make<TH1F>( "muon 2 eta", "muon 2 eta", 250, -2.5, 2.5);
    hMuonPhi[process] = Muons.make<TH1F>( "muon phi", "muon phi", 350, -3.5, 3.5);
    //Muon quality variables
    hMuonChi2[process] = Muons.make<TH1F>( "muon track chi2 per dof", "muon track chi2 / dof", 200, 0.0, 20.0);
    hMuond0Pv[process] = Muons.make<TH1F>( "muon track d0 pv", "muon track d0 pv", 400, -0.2, 0.2);
    hMuond0Bs[process] = Muons.make<TH1F>( "muon track d0 bs", "muon track d0 bs", 400, -0.2, 0.2);
    hMuond0SigPv[process] = Muons.make<TH1F>( "muon track d0 significance pv", "muon track d0 significance pv", 100, -50., 50.0);
    hMuond0SigBs[process] = Muons.make<TH1F>( "muon track d0 significance bs", "muon track d0 significance bs", 100, -50., 50.0);
    hMuonnHits[process] = Muons.make<TH1F>( "muon track hits", "muon track number of valid hits", 50, 0.0, 50.0);
    hMuonPfIsod0Pv[process] = Muons.make<TH2F>( "muon pf iso d0 pv", "muon pf iso d0 pv", 300, 0.0 , 3.0, 200, -0.2, 0.2);
    hMuonPfIsod0Bs[process] = Muons.make<TH2F>( "muon pf iso d0 bs", "muon pf iso d0 bs", 300, 0.0 , 3.0, 200, -0.2, 0.2);
    hMuonEtaPhi[process] = Muons.make<TH2F>( "muon eta phi", "muon eta phi", 250, -2.5 , 2.5, 350, -3.5, 3.5);
    //histograms for lepton isolation cuts
    hMuonIso[process] = Muons.make<TH1F>( "muon iso", "Isolation of muons", 300, 0.0, 3.0);
    hMuonTrackIso[process] = Muons.make<TH1F>( "muon track iso", "Isolation of muons in tracker", 1000, 0.0, 10.0);
    hMuonEcalIso[process] = Muons.make<TH1F>( "muon ecal iso", "Isolation of muons in electromagnetic calorimeter", 1000, 0.0, 10.0);
    hMuonHcalIso[process] = Muons.make<TH1F>( "muon hcal iso", "Isolation of muons in hadronic calorimeter", 1000, 0.0, 10.0);
    hMuonPfIso[process] = Muons.make<TH1F>( "muon pf iso", "Isolation of pf muons", 300, 0.0, 3.0);
    hMuonChargedHadronIso[process] = Muons.make<TH1F>( "muon charged hadron iso", "Isolation of muons from charged hadrons", 1000, 0.0, 10.0);
    hMuonPhotonIso[process] = Muons.make<TH1F>( "muon photon iso", "Isolation of muons from photons", 1000, 0.0, 10.0);
    hMuonNeutralHadronIso[process] = Muons.make<TH1F>( "muon neutral hadron iso", "Isolation of muons from neutral hadrons", 1000, 0.0, 10.0);
    hGenMuonPt[process] = Muons.make<TH1F>( "generator muon pt", "Generator muon pt", 1000, 0.0, 1000.0);
    hGenMuonEta[process] = Muons.make<TH1F>( "generator muon eta", "Generator muon eta", 250, -2.5, 2.5);

    hMuonResolution[process] = Muons.make<TH1F>( "muon resolution", "Resolution of muons", 1000, -5.0, 5.0);
    
    TFileDirectory Electrons = theFile->mkdir("Electrons"); 
    //electron histograms
    hElectronPt[process] = Electrons.make<TH1F>( "electron pt", "electron pt", 1000, 0.0, 1000.0);
    hElectronCharge[process] = Electrons.make<TH1F>( "electron charge", "electron charge", 8, -2.0, 2.0);
    hElectronSumPt[process] = Electrons.make<TH1F>( "sum electron pt", "sum electron pt", 1000, 0.0, 1000.0);
    hElectron1Pt[process] = Electrons.make<TH1F>( "electron 1 pt", "pt of first electron", 1000, 0.0, 1000.0);
    hElectron2Pt[process] = Electrons.make<TH1F>( "electron 2 pt", "pt of second electron", 1000, 0.0, 1000.0);
    hElectronEta[process] = Electrons.make<TH1F>( "electron eta", "electron eta", 250, -2.5, 2.5);
    hElectron1Eta[process] = Electrons.make<TH1F>( "electron 1 eta", "electron 1 eta", 250, -2.5, 2.5);
    hElectron2Eta[process] = Electrons.make<TH1F>( "electron 2 eta", "electron 2 eta", 250, -2.5, 2.5);
    hElectronPhi[process] = Electrons.make<TH1F>( "electron phi", "electron phi", 350, -3.5, 3.5);
    hElectronLostHits[process] = Electrons.make<TH1F>( "electron lost hits", "electron lost hits", 15, -0.5, 14.5);
    hElectrond0Pv[process] = Electrons.make<TH1F>( "electron track d0 pv", "electron track d0 pv", 400, -0.2, 0.2);
    hElectrond0Bs[process] = Electrons.make<TH1F>( "electron track d0 bs", "electron track d0 bs", 400, -0.2, 0.2);
    hElectrond0SigPv[process] = Electrons.make<TH1F>( "electron track d0 significance pv", "electron track d0 significance pv", 100, -50., 50.0);
    hElectrond0SigBs[process] = Electrons.make<TH1F>( "electron track d0 significance bs", "electron track d0 significance bs", 100, -50., 50.0);
    //histograms for lepton isolation cuts
    hElectronIso[process] = Electrons.make<TH1F>( "electron iso", "Isolation of electrons", 300, 0.0, 3.0);
    hElectronTrackIso[process] = Electrons.make<TH1F>( "electron track iso", "Isolation of electrons in tracker", 1000, 0.0, 10.0); 
    hElectronEcalIso[process] = Electrons.make<TH1F>( "electron ecal iso", "Isolation of electrons in electromagnetic calorimeter", 1000, 0.0, 10.0); 
    hElectronHcalIso[process] = Electrons.make<TH1F>( "electron hcal iso", "Isolation of electrons in hadronic calorimeter", 1000, 0.0, 10.0); 
    hElectronPfIso[process] = Electrons.make<TH1F>( "electron pf iso", "Isolation of pf electrons", 300, 0.0, 3.0);
    hElectronChargedHadronIso[process] = Electrons.make<TH1F>( "electron charged hadron iso", "Isolation of electrons from charged hadrons", 1000, 0.0, 10.0);
    hElectronPhotonIso[process] = Electrons.make<TH1F>( "electron photon iso", "Isolation of electron from photons", 1000, 0.0, 10.0);
    hElectronNeutralHadronIso[process] = Electrons.make<TH1F>( "electron neutral hadron iso", "Isolation of electrons from neutral hadrons", 1000, 0.0, 10.0);
    hGenElectronPt[process] = Electrons.make<TH1F>( "generator electron pt", "Generator electron pt", 1000, 0.0, 1000.0);
    hGenElectronEta[process] = Electrons.make<TH1F>( "generator electron eta", "Generator electron eta", 250, -2.5, 2.5);
    //electron variables
    hElectronMva[process] = Electrons.make<TH1F>( "electron mva", "electron mva", 200, -1., 1.);
    hElectronEoverP[process] = Electrons.make<TH1F>( "electron E over P", "electron E over P", 250, 0.0, 2.5);
    hElectronfBrem[process] = Electrons.make<TH1F>( "electron fBrem", "electron fBrem", 110, 0.0, 1.1);
    hElectronHoverE[process] = Electrons.make<TH1F>( "electron H over E", "electron H over E", 200, 0.0, 0.2);
    hElectrondeltaPhiIn[process] = Electrons.make<TH1F>( "electron deltaPhiIn", "electron deltaPhiIn", 200, -0.1, 0.1);
    hElectrondeltaEtaIn[process] = Electrons.make<TH1F>( "electron deltaEtaIn", "electron deltaEtaIn", 400, -0.2, 0.2);
    hElectrone25Maxoe55[process] = Electrons.make<TH1F>( "electron e25Maxoe55", "electron e25Maxoe55", 100, -1, 2);
    hElectrone15oe55[process] = Electrons.make<TH1F>( "electron e15oe55", "electron e15oe55", 100, -1, 2);
    hElectronsigmaEtaEta[process] = Electrons.make<TH1F>( "electron sigma eta eta", "electron sigma eta eta", 100, 0, 0.04);
    hElectronsigmaIetaIeta[process] = Electrons.make<TH1F>( "electron sigma Ieta Ieta", "electron sigma Ieta Ieta", 100, 0, 0.04);
    hElectronPfIsod0Pv[process] = Electrons.make<TH2F>( "electron pf iso d0 pv", "electron pf iso d0 pv", 300, 0.0 , 3.0, 200, -0.2, 0.2);
    hElectronPfIsod0Bs[process] = Electrons.make<TH2F>( "electron pf iso d0 bs", "electron pf iso d0 bs", 300, 0.0 , 3.0, 200, -0.2, 0.2);
    hElectronEtaPhi[process] = Electrons.make<TH2F>( "electron eta phi", "electron eta phi", 250, -2.5 , 2.5, 350, -3.5, 3.5);

    hElectronResolution[process] = Electrons.make<TH1F>( "electron resolution", "Resolution of electrons", 1000, -5.0, 5.0);
    
    TFileDirectory Taus = theFile->mkdir("Taus"); 
    //tau histograms
    hTauPt[process] = Taus.make<TH1F>( "tau pt", "tau pt", 1000, 0.0, 1000.0);
    hTauCharge[process] = Taus.make<TH1F>( "tau charge", "tau charge", 8, -2.0, 2.0);
    hTauSumPt[process] = Taus.make<TH1F>( "sum tau pt", "sum tau pt", 1000, 0.0, 1000.0);
    hTau1Pt[process] = Taus.make<TH1F>( "tau 1 pt", "pt of first tau", 1000, 0.0, 1000.0);
    hTau2Pt[process] = Taus.make<TH1F>( "tau 2 pt", "pt of second tau", 1000, 0.0, 1000.0);
    hTauEta[process] = Taus.make<TH1F>( "tau eta", "tau eta", 250, -2.5, 2.5);
    hTau1Eta[process] = Taus.make<TH1F>( "tau 1 eta", "tau 1 eta", 250, -2.5, 2.5);
    hTau2Eta[process] = Taus.make<TH1F>( "tau 2 eta", "tau 2 eta", 250, -2.5, 2.5);
    hTauPhi[process] = Taus.make<TH1F>( "tau phi", "tau phi", 350, -3.5, 3.5);
    hTauEtaPhi[process] = Taus.make<TH2F>( "tau eta phi", "tau eta phi", 250, -2.5 , 2.5, 350, -3.5, 3.5);
    hTauIso[process] = Taus.make<TH1F>( "tau iso", "Isolation of taus", 300, 0.0, 3.0);
    hTauTrackIso[process] = Taus.make<TH1F>( "tau track iso", "Isolation of taus in tracker", 1000, 0.0, 10.0);
    hTauEcalIso[process] = Taus.make<TH1F>( "tau ecal iso", "Isolation of taus in electromagnetc calorimeter", 1000, 0.0, 10.0);
    hTauHcalIso[process] = Taus.make<TH1F>( "tau hcal iso", "Isolation of taus in hadronic calorimeter", 1000, 0.0, 10.0);
    hTauPfIso[process] = Taus.make<TH1F>( "tau pf iso", "Isolation of pf taus", 300, 0.0, 3.0);
    hTauChargedHadronIso[process] = Taus.make<TH1F>( "tau charged hadron iso", "Isolation of taus from charged hadrons", 1000, 0.0, 10.0);
    hTauPhotonIso[process] = Taus.make<TH1F>( "tau photon iso", "Isolation of taus from photons", 1000, 0.0, 10.0);
    hTauNeutralHadronIso[process] = Taus.make<TH1F>( "tau neutral hadron iso", "Isolation of taus from neutral hadrons", 1000, 0.0, 10.0);
    hGenTauPt[process] = Taus.make<TH1F>( "generator tau pt", "matched gen tau pt", 1000, 0.0, 1000.0);
    hGenTauVisPt[process] = Taus.make<TH1F>( "generator tau vis pt", "matched gen tau visible pt", 1000, 0.0, 1000.0);
    hTauDiscriminators[process] = Taus.make<TH1F>( "tau discriminators", "Tau ID discriminators", maxTauDiscriminators_+1, 0.0, maxTauDiscriminators_+1);
    
    h2dMuonEtaPt[process] = Muons.make<TH2F>( "muon eta pt", "muon eta pt", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dMatchedMuonEtaPt[process] = Muons.make<TH2F>( "muon matched eta pt", "muon matched eta pt", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dGenMuonEtaPt[process] = Muons.make<TH2F>( "muon generator eta pt", "muon generator eta pt", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dElectronEtaPt[process] = Electrons.make<TH2F>( "electron eta pt", "electron eta pt", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dMatchedElectronEtaPt[process] = Electrons.make<TH2F>( "electron matched eta pt", "electron matched eta pt", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dGenElectronEtaPt[process] = Electrons.make<TH2F>( "electron generator eta pt", "electron generator eta pt", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dTauEtaPt[process] = Taus.make<TH2F>( "tau eta pt", "tau eta pt", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dMatchedTauEtaPt[process] = Taus.make<TH2F>( "tau matched eta pt", "tau matched eta pt", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
  
    //histograms for Missing ET
    TFileDirectory MET = theFile->mkdir("MET"); 
    hMissingET[process] = MET.make<TH1F>( "MET", "Missing transverse energy", 1000, 0.0, 1000.0);
    hMissingETmc[process] =  MET.make<TH1F>( "MET MC", "Missing transverse energy MC", 1000, 0.0, 1000.0);
    halphaT[process] = MET.make<TH1F>( "alphaT", "alphaT", 100, 0.0, 1.0);
    hEtSum[process] = MET.make<TH1F>( "ETsum", "Transverse energy sum ET", 2000, 0.0, 2000.0);
    hHT[process] = MET.make<TH1F>( "HT", "Transverse energy sum HT", 2000, 0.0, 2000.0);
    hMHT[process] = MET.make<TH1F>( "MHT", "Transverse energy sum HT + MET", 2000, 0.0, 2000.0);
    h2dMETEtSumJets[process] = MET.make<TH2F>( "MET - sum4Jets", "MET - sum4Jets", 2000, 0.0, 2000.0, 1000, 0.0, 1000.0);
    h2dMETHT[process] = MET.make<TH2F>( "MET - HT", "MET - HT", 2000, 0.0, 2000.0, 1000, 0.0, 1000.0);
    h2dMETMHT[process] = MET.make<TH2F>( "MET - MHT", "MET - MHT", 2000, 0.0, 2000.0, 1000, 0.0, 1000.0);
    hJet1METPhi[process] = MET.make<TH1F>( "jet 1 MET delta phi", "jet 1 delta phi", 300, -3., 3.);
    hJet1EtaMET[process] = MET.make<TH2F>( "jet 1 eta MET", "jet 1 eta MET", 1000, 0.0, 1000.0, 250, -2.5, 2.5);

    //histograms for jets
    TFileDirectory Jets = theFile->mkdir("Jets"); 
    hJetEt[process] = Jets.make<TH1F>( "jet pt", "jet pt", 1500, 0.0, 1500.0);
    hJetEta[process] = Jets.make<TH1F>( "jet eta", "jet eta", 300, -3., 3.);
    hJetPhi[process] = Jets.make<TH1F>( "jet phi", "jet phi", 350, -3.5, 3.5);
    hEtJet1[process] = Jets.make<TH1F>( "pt first jet", "Pt spectrum of the 1st jet", 1500, 0.0, 1500.0);
    hEtJet2[process] = Jets.make<TH1F>( "pt second jet", "Pt spectrum of the 2nd jet", 1500, 0.0, 1500.0);
    hEtJet3[process] = Jets.make<TH1F>( "pt third jet", "Pt spectrum of the 3rd jet", 1500, 0.0, 1500.0);
    hEtJet4[process] = Jets.make<TH1F>( "pt fourth jet", "Pt spectrum of the 4th jet", 1500, 0.0, 1500.0);
    hJet1Eta[process] = Jets.make<TH1F>( "jet 1 eta", "jet 1 eta", 300, -3., 3.);
    hJet2Eta[process] = Jets.make<TH1F>( "jet 2 eta", "jet 2 eta", 300, -3., 3.);
    hJet3Eta[process] = Jets.make<TH1F>( "jet 3 eta", "jet 3 eta", 300, -3., 3.);
    hJet4Eta[process] = Jets.make<TH1F>( "jet 4 eta", "jet 4 eta", 300, -3., 3.);
    hJetCHF[process] = Jets.make<TH1F>( "jet charged hadron fraction", "jet charged hadron fraction", 100, 0., 1.);
    hJetNHF[process] = Jets.make<TH1F>( "jet neutral hadron fraction", "jet neutral hadron fraction", 100, 0., 1.);
    hJetCEF[process] = Jets.make<TH1F>( "jet charged electromagnetic fraction", "jet charged electromagnetic fraction", 100, 0., 1.);
    hJetNEF[process] = Jets.make<TH1F>( "jet neutral electromagnetic fraction", "jet neutral electromagnetic fraction", 100, 0., 1.);
    hJetCM[process] = Jets.make<TH1F>( "jet charged multiplity", "jet charged multiplicity", 40, -0.5, 39.5);
    hJetNConst[process] = Jets.make<TH1F>( "jet number of constituents", "jet number of constituents", 40, -0.5, 39.5);
    h2dJetEtaPhi[process] = Jets.make<TH2F>( "jet eta phi", "jet eta phi", 350, -3.5 , 3.5, 350, -3.5, 3.5);
    h2dJetEtaPt[process] = Jets.make<TH2F>( "jet eta pt", "jet eta pt", 1000, 0.0, 1000.0, 300, -3.0, 3.0); 
    h2dMatchedJetEtaPt[process] = Jets.make<TH2F>( "jet matched eta pt", "jet matched eta pt", 1000, 0.0, 1000.0, 300, -3.0, 3.0); 
 
    hTrigger[process] = theFile->make<TH1F>( "Trigger paths", "Trigger paths", 160, 0, 160);
    hWeight[process] = theFile->make<TH1F>( "Weights", "Weights", 1000, 0, 1000);
}

//Destructor
DiLeptonHistograms::~DiLeptonHistograms()
{
    PrintStatistics();
    if (debug) std::cout << "************* Finished analysis" << std::endl;
} 

void DiLeptonHistograms::ReadEfficiency(){
    edm::LogPrint("Summary")  << "Reading muon efficiencies from file " << muon_fname << "\n"; 
    std::ifstream muon_listfile(edm::FileInPath(muon_fname).fullPath().c_str());
    edm::LogPrint("Summary")  << "Reading electron efficiencies from file " << electron_fname << "\n";
    std::ifstream electron_listfile(edm::FileInPath(electron_fname).fullPath().c_str());
   
    float muon_eff=0.;
    float electron_eff=0.;
    int muon_nEta;
    int electron_nEta;
    int muon_nPt;
    int electron_nPt;
    muon_listfile >> muon_nEta;
    muon_listfile >> muon_nPt;
    electron_listfile >> electron_nEta;
    electron_listfile >> electron_nPt;
    if (muon_nEta != nEtaBins && electron_nEta != nEtaBins && muon_nPt != nPtBins && electron_nPt != nPtBins) { 
        std::cout << " *** ERROR -> Check data files : nEta bins " 
        << muon_nEta << " instead of " << nEtaBins << std::endl;
    }
    else{
        for (int i=0; i<nEtaBins; i++) {
            for (int j=0; j<nPtBins; j++) {
                muon_listfile >> muon_eff;
                electron_listfile >> electron_eff;
                Muon_Eff[i][j]=muon_eff;
                Electron_Eff[i][j]=electron_eff;
                if (debug) std::cout << "Cat: " << i << "," << j << " muon: " << Muon_Eff[i][j] << " electron: " << Electron_Eff[i][j] << std::endl; 
            }
        }
    }
}    

//calculate the isolation of a pat lepton relative sum in cone: (tracker+ecal+hcal)/lept.pt
template < class T > 
double CalcIso(const T & lepton)
{
    //double cut_ConeSize = 0.3;
    //reco::isodeposit::Direction candDir(lepton.eta(), lepton.phi());
    double value = (lepton.trackIso()+lepton.ecalIso()+lepton.hcalIso())/lepton.pt();
    /*double value =  (lepton.trackerIsoDeposit()->depositWithin(cut_ConeSize)+
                    lepton.ecalIsoDeposit()->depositWithin(cut_ConeSize)+
                    lepton.hcalIsoDeposit()->depositWithin(cut_ConeSize))/lepton.pt();*/
    return value;
}

template < class T > 
double CalcPfIso(const T & lepton)
{
    double value = (lepton.chargedHadronIso()+lepton.photonIso()+lepton.neutralHadronIso()/3.)/lepton.pt();
    return value;
}

template < class T > 
int GetLeptKind(const T * lepton)
{
    int value = 0;
    if(lepton->genLepton()){
        if(lepton->genLepton()->status()==1 && lepton->genLepton()->numberOfMothers()==1){
            const reco::Candidate * mom = lepton->genLepton()->mother();
            //Check if lepton is promt (itself,tau,Z,W,SUSY)
            if(mom->pdgId()==lepton->genLepton()->pdgId()||abs(mom->pdgId())==15||abs(mom->pdgId())==23||abs(mom->pdgId())==24||abs(mom->pdgId())>1000000){
                value=3;
            }
            else {
                //LogPrint("Lepton") << "Non promt: " << mom->pdgId();
                value=4;
            }
        }
        else value = 4;
    }
    else value=2;
    return value;
}


//Find the combination of objects with the closest distance in pt
bool DiLeptonHistograms::FindMinCombo(std::vector<const reco::Candidate*> objects,
				std::vector<unsigned int> & lista,
				std::vector<unsigned int> & listb) {
 
    unsigned int n = objects.size();
    //if (n>Combinations::nmax) {LogPrint("alhpaCalc") << "FindMinDptCombo: " << n << " too big"; return false;}

    lista.clear(); listb.clear(); // clears the lists, just to be sure

    double mindiff = 1000000000., diff = 0.;
    // Strip indices from the relevant combination set
    std::vector<unsigned int> la; //!< Temporary list a for calculating best combo
    std::vector<unsigned int> lb; //!< Temporary list b for calculating best combo

    for (unsigned int r=0; r<Combinations::rmax[n]; r++) {
        for (unsigned int j=0; j<Combinations::jmax[n]; j++) {
        //std::cout << r << " " << j << std::endl;
        // populate list a
        for (unsigned int ia=0; ia<(n-1); ia++) {
	    if ((n==2) && (Combinations::la2[r][j][ia]!=-1) ) la.push_back( Combinations::la2[r][j][ia] );	
	    if ((n==3) && (Combinations::la3[r][j][ia]!=-1) ) la.push_back( Combinations::la3[r][j][ia] );
	    if ((n==4) && (Combinations::la4[r][j][ia]!=-1) ) la.push_back( Combinations::la4[r][j][ia] );
	    if ((n==5) && (Combinations::la5[r][j][ia]!=-1) ) la.push_back( Combinations::la5[r][j][ia] );
	    if ((n==6) && (Combinations::la6[r][j][ia]!=-1) ) la.push_back( Combinations::la6[r][j][ia] );
	    if ((n==7) && (Combinations::la7[r][j][ia]!=-1) ) la.push_back( Combinations::la7[r][j][ia] );
	//if (lista[ia] < 0) { lista.pop_back(); break; }
        }
        // populate list b
        for (unsigned int ib=0; ib<(n-1); ib++) {
	    if ((n==2) && (Combinations::lb2[r][j][ib]!=-1) ) lb.push_back( Combinations::lb2[r][j][ib] );
	    if ((n==3) && (Combinations::lb3[r][j][ib]!=-1) ) lb.push_back( Combinations::lb3[r][j][ib] );
	    if ((n==4) && (Combinations::lb4[r][j][ib]!=-1) ) lb.push_back( Combinations::lb4[r][j][ib] );
	    if ((n==5) && (Combinations::lb5[r][j][ib]!=-1) ) lb.push_back( Combinations::lb5[r][j][ib] );
	    if ((n==6) && (Combinations::lb6[r][j][ib]!=-1) ) lb.push_back( Combinations::lb6[r][j][ib] );
	    if ((n==7) && (Combinations::lb7[r][j][ib]!=-1) ) lb.push_back( Combinations::lb7[r][j][ib] );
	    //if (listb[ib] < 0) { listb.pop_back(); break; }
        }


        if ( la.size()==0 || lb.size()==0 ) break;

        double JetaEt = 0., JetbEt = 0.;
        for (std::vector<unsigned int>::iterator ia=la.begin();ia!=la.end();++ia) {
	//std::cout << (*ia) << " ";
	JetaEt += objects[ (*ia) ]->pt();
        }
        //std::cout << ", ";
        for (std::vector<unsigned int>::iterator ib=lb.begin();ib!=lb.end();++ib) {
	    JetbEt += objects[ (*ib) ]->pt();
	    //std::cout << (*ib) << " ";
        }
        //std::cout << std::endl;
        diff = fabs(JetaEt - JetbEt);
        //std::cout << "Difference in Et is " << diff << std::endl;
        if (diff < mindiff) { mindiff = diff; lista = la; listb = lb; }
        la.clear(); lb.clear();
        } // end of permutation list loop
    } // end of combination list loop
  
    return true;
}


double DiLeptonHistograms::CalcalphaT(std::vector<const reco::Candidate*> objects,
			  std::vector<unsigned int> la,
			  std::vector<unsigned int> lb) {

    // Check we have enough jets for the lists supplied.
    if ((la[la.size()-1]>objects.size())||(lb[lb.size()-1]>objects.size())) {return 0.;}
    
    double jetaEt = 0., jetbEt = 0.; // sums for Et
    // Loop over jet list a to get pseudo-jet a
    for ( std::vector<unsigned int>::iterator ia = la.begin(); ia != la.end(); ++ia ) {
      jetaEt += objects[(*ia)]->pt();
    }
      
    // Loop over jet list b to get pseudo-jet b
    for ( std::vector<unsigned int>::iterator ib = lb.begin(); ib != lb.end(); ++ib ) {
      jetbEt += objects[(*ib)]->pt();
    }
    double ptSum = 0., pxSum = 0., pySum = 0.;
    for (unsigned int i=0; i < objects.size(); ++i){
        ptSum += objects[i]->pt();
        pxSum += objects[i]->px();
        pySum += objects[i]->py();
    }
    double MT = sqrt(ptSum*ptSum-pxSum*pxSum-pySum*pySum);
    double alphaT = std::min(jetaEt,jetbEt)/MT;
    return alphaT; 
}


//Calculate the muon effificiency
double DiLeptonHistograms::getMuonWeight(const pat::Muon* muon)
{
    int catEta = 0;
    int catPt = 0;
    double lowPt = 2;
    for(int i=0; i<nEtaBins; ++i){
        if(boundEta[i]<muon->eta() && muon->eta()<boundEta[i+1]){catEta=i;}
    }
    for(int j=0; j<nPtBins; ++j){
        if(boundPt[j]<muon->pt() && muon->pt()<boundPt[j+1]){catPt=j;}
    }
    float eff = Muon_Eff[catEta][catPt];
    if (debug) std::cout << "Muon eta: " << muon->eta() << " pt: " << muon->pt() << " category: " << catEta << "," << catPt << " Eff: " << eff << std::endl;
    if(eff!=0.) return 1/(eff+eff_Muon_Scale*eff+lowPt/muon->pt()*eff);
    else return 0;
}

//Calculate the electron effificiency
double DiLeptonHistograms::getElectronWeight(const pat::Electron* electron)
{
    int catEta = 0;
    int catPt = 0;
    double lowPt = 2;
    for(int i=0; i<nEtaBins; ++i){
        if(boundEta[i]<electron->eta() && electron->eta()<boundEta[i+1]){catEta=i;}
    }
    for(int j=0; j<nPtBins; ++j){
        if(boundPt[j]<electron->pt() && electron->pt()<boundPt[j+1]){catPt=j;}
    }
    float eff = Electron_Eff[catEta][catPt];
    if (debug) std::cout << "Electron eta: " << electron->eta() << " pt: " << electron->pt() << " category: " << catEta << "," << catPt << " Eff: " << eff << std::endl;
    if(eff!=0.) return 1/(eff+eff_Electron_Scale*eff+lowPt/electron->pt()*eff);
    else return 0;
}


//Filling of all histograms and calculation of kinematics
//main part
void DiLeptonHistograms::Analysis(const edm::Handle< std::vector<pat::Muon> >& muons, const edm::Handle< std::vector<pat::Electron> >& electrons, const edm::Handle< std::vector<pat::Tau> >& taus, const edm::Handle< std::vector<pat::Jet> >& jets, const edm::Handle< std::vector<pat::MET> >& met, double weight, const int process){
    hWeight[process]->Fill(weight,weight);
 
    double MET=0;
    pat::MET meti;
    for (std::vector<pat::MET>::const_iterator met_i = met->begin(); met_i != met->end(); ++met_i){
        if (debug) std::cout <<"MET et = "<< met_i->et() << std::endl;
        MET = met_i->et();
        meti = *met_i;
    }

    hMissingET[process]->Fill(MET,weight);
  
    int n_Jet=0;
    int n_bJet=0;
    std::vector< pat::Jet > bJets;
    float et4Jets=0;
    float etFourthJet=0;
    float HT=0;

    std::vector<const reco::Candidate*> objects;
    math::PtEtaPhiMLorentzVector JPt(0., 0., 0., 0.);

    for (std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i){       
        if (debug) std::cout <<"Jet et = "<< jet_i->pt() << std::endl;
        ++numTotJets;
	    ++n_Jet;
        JPt += jet_i->p4();
        objects.push_back(static_cast<const reco::Candidate*>( &(*jet_i) ));
	    //Plots of leading Jets
   	    JetMonitor(&(*jet_i),n_Jet,weight,general); 
	    if(n_Jet==1){
            hJet1METPhi[process]->Fill(reco::deltaPhi(jet_i->phi(),meti.phi()),weight);
            hJet1EtaMET[process]->Fill(MET,jet_i->eta(),weight);
	    }
	    if(n_Jet<5){
		et4Jets+=jet_i->pt();
		etFourthJet=jet_i->pt();}
	    //Jet base plots
	    if(jet_i->bDiscriminator(bJetAlgo)>cut_bTagDiscriminator){
	        ++n_bJet;
		bJets.push_back( *jet_i );
	    }
	    HT += jet_i->pt();
    }

    hJetMult[process]->Fill(n_Jet,weight);
    hbJetMult[process]->Fill(n_bJet,weight);
    if (debug) std::cout <<" Number of Jets = " << n_Jet << std::endl;

    //Sum of four leading jets
    hEtSum[process]->Fill(et4Jets,weight);
    h2dMETEtSumJets[process]->Fill(et4Jets,MET,weight);
    //Sum of jet et
    hHT[process]->Fill(HT,weight);
    h2dMETHT[process]->Fill(HT,MET,weight);
    //Sum of jet et + MET
    hMHT[process]->Fill(HT+MET,weight);
    h2dMETMHT[process]->Fill(HT+MET,MET,weight);
    
    //Inv mass variables
    invMOFOS = 0; 
    invMSFOS = 0; 
    invMMuon = 0; 
    invMElec = 0; 
    invweight = 1.; 
    double inv = 0;
    double res = 0;
    double weightcorr = 0;
    double iso = 0;
    double dileptonPt = 0;
    float deltaPhiMET = 0;
    //Muon histograms
    int n_Muons = 0;
    float muonPt = 0.;
    for (std::vector<pat::Muon>::const_iterator mu_i = muons->begin(); mu_i != muons->end(); ++mu_i){
        if (debug) std::cout <<"mu eta = "<< mu_i->eta() << std::endl;
        ++numTotMuons;
        ++n_Muons;
   	    MuonMonitor(&(*mu_i),n_Muons,weight,general); 
        if(mcInfo){MuonMonitor(&(*mu_i),n_Muons,weight,GetLeptKind(&(*mu_i)));}   
	    //Clean and isolated muons
        objects.push_back(static_cast<const reco::Candidate*>( &(*mu_i) ));
        muonPt += mu_i->pt();
   	    if (effInfo) MuonMonitor(&(*mu_i),n_Muons,weight,effcor); 
	
        //Invariant mass plots 
	    //Muon pairs
   	    for (std::vector<pat::Muon>::const_iterator mu_j = muons->begin(); mu_j != muons->end(); ++mu_j){
       	    if( mu_i->charge()==+1 && mu_j->charge()==-1){
	        if (debug) std::cout <<"Invariant Mass mu+ mu-: " << (mu_i->p4()+mu_j->p4()).M() <<std::endl;
                inv = (mu_i->p4()+mu_j->p4()).M();
                iso = CalcIso(*mu_i)+CalcIso(*mu_j);
                deltaPhiMET = reco::deltaPhi((mu_i->p4()+mu_j->p4()).phi(),meti.phi());
                weightcorr = getMuonWeight(&(*mu_i))*getMuonWeight(&(*mu_j))*weight;
                hDeltaPhiMET[general]->Fill(deltaPhiMET,weight);
                if (mcInfo){
                    if (mu_i->genLepton()&&mu_j->genLepton()){
                        res = (mu_i->genLepton()->p4()+mu_j->genLepton()->p4()).M()-inv;
                    }
                }
                hMuonResolution[general]->Fill(res,weight);
                MuonInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weight,general);
                if (effInfo) MuonInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weightcorr,effcor);
                hDeltaPhiMET[general]->Fill(deltaPhiMET,weight);
                if (effInfo) hDeltaPhiMET[effcor]->Fill(deltaPhiMET,weightcorr);
                dileptonPt = (mu_i->p4()+mu_j->p4()).pt();
                hDileptonPt[general]->Fill(dileptonPt,weight);
                hJZB[general]->Fill(JPt.pt()-dileptonPt,weight);
                //LeptCorrMonitor(*mu_i,*mu_j,meti, weightcorr,effcor);
                invMSFOS = (mu_i->p4()+mu_j->p4()).M();
                invMMuon = (mu_i->p4()+mu_j->p4()).M();
                invweight = weightcorr;
                if (treeInfo) treeSFOS->Fill();
                if (treeInfo) treeMuon->Fill();
	            if(n_bJet>=2){
		            for(std::vector<pat::Jet>::const_iterator bjet_i = bJets.begin(); bjet_i!=bJets.end(); ++bjet_i){
                        for(std::vector<pat::Jet>::const_iterator bjet_j = bJets.begin(); bjet_j!=bJets.end();++bjet_j){
		                    if(bjet_i<bjet_j){hInvMbbllSFOS[process]->Fill((bjet_i->p4()+bjet_j->p4()+mu_i->p4()+mu_j->p4()).M(),weight);}
		                }
		            }
	            }
            }
       	    if( mu_i->charge()==mu_j->charge()&&mu_i<mu_j){hInvMMuonSS[process]->Fill((mu_i->p4()+mu_j->p4()).M(),weight);}
        }
            //Wrong pairings for different flavour subtraction
            for (std::vector<pat::Electron>::const_iterator ele_j = electrons->begin(); ele_j != electrons->end(); ++ele_j){
       	        if(( mu_i->charge()==+1 && ele_j->charge()==-1)|(mu_i->charge()==-1 && ele_j->charge()==+1)){
	                if (debug) std::cout <<"Invariant Mass mu+ e-, mu- e+: " << (mu_i->p4()+ele_j->p4()).M() <<std::endl;
                    inv = (mu_i->p4()+ele_j->p4()).M();
                    weightcorr = getMuonWeight(&(*mu_i))*getElectronWeight(&(*ele_j))*weight;
                    iso = CalcIso(*mu_i)+CalcIso(*ele_j);
                    OFOSInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weight,general);
                    if (effInfo) OFOSInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weightcorr,effcor);
                    dileptonPt = (mu_i->p4()+ele_j->p4()).pt();
                    hDileptonPt[general]->Fill(dileptonPt,weight);
                    hJZB[general]->Fill(JPt.pt()-dileptonPt,weight);
                    invMOFOS = (mu_i->p4()+ele_j->p4()).M();
                    invweight = weightcorr;
                    if (treeInfo) treeOFOS->Fill();
                }
	            if(n_bJet>=2){
	                for(std::vector<pat::Jet>::const_iterator bjet_i = bJets.begin(); bjet_i!=bJets.end(); ++bjet_i){
                        for(std::vector<pat::Jet>::const_iterator bjet_j = bJets.begin(); bjet_j!=bJets.end();++bjet_j){
		                    if(bjet_i<bjet_j){hInvMbbllOFOS[process]->Fill((bjet_i->p4()+bjet_j->p4()+mu_i->p4()+ele_j->p4()).M(),weight);}
		                }
	                }
	            }
            }
    }            
    
    //Loop over electrons 
    int n_Electrons = 0;
    float elePt = 0.;
    for (std::vector<pat::Electron>::const_iterator ele_i = electrons->begin(); ele_i != electrons->end(); ++ele_i){
        if (debug) std::cout <<"ele eta = "<< ele_i->eta() << std::endl;
        ++numTotElectrons;
	    ++n_Electrons;
   	    ElectronMonitor(&(*ele_i),n_Electrons,weight,general); 
        if(mcInfo){ElectronMonitor(&(*ele_i),n_Electrons,weight,GetLeptKind(&(*ele_i)));}  
        elePt += ele_i->pt();
        objects.push_back(static_cast<const reco::Candidate*>( &(*ele_i) ));
   	    if (effInfo) ElectronMonitor(&(*ele_i),n_Electrons,weight,effcor);
        
        //Invariant mass plots
	    //Electron pairs
   	    for (std::vector<pat::Electron>::const_iterator ele_j = electrons->begin(); ele_j != electrons->end(); ++ele_j){
       	    if( ele_i->charge()==-1 && ele_j->charge()==+1){
	        if (debug) std::cout <<"Invariant Mass e+ e-: " << (ele_i->p4()+ele_j->p4()).M() <<std::endl;
                inv = (ele_i->p4()+ele_j->p4()).M();
                if(electron_Scale<-0.00001||electron_Scale>0.00001){
                    reco::Particle::LorentzVector pb_i = reco::Particle::LorentzVector((ele_i->px()+ele_i->px()*electron_Scale),(ele_i->py()+ele_i->py()*electron_Scale),(ele_i->pz()+ele_i->pz()*electron_Scale),(ele_i->energy()+ele_i->energy()*electron_Scale));
                    reco::Particle::LorentzVector pb_j = reco::Particle::LorentzVector((ele_j->px()+ele_j->px()*electron_Scale),(ele_j->py()+ele_j->py()*electron_Scale),(ele_j->pz()+ele_j->pz()*electron_Scale),(ele_j->energy()+ele_j->energy()*electron_Scale));
                    inv = (pb_i+pb_j).M();
                }
                weightcorr = getElectronWeight(&(*ele_i))*getElectronWeight(&(*ele_j))*weight;
                iso = CalcIso(*ele_i)+CalcIso(*ele_j);
                ElectronInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weight,general);
                if (mcInfo){
                    if (ele_i->genLepton()&&ele_j->genLepton()){
                            res = (ele_i->genLepton()->p4()+ele_j->genLepton()->p4()).M()-inv;
                    }
                }
                hElectronResolution[general]->Fill(res,weight);
                if (effInfo) ElectronInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weightcorr,effcor);
                dileptonPt = (ele_i->p4()+ele_j->p4()).pt();
                hDileptonPt[general]->Fill(dileptonPt,weight);
                hJZB[general]->Fill(JPt.pt()-dileptonPt,weight);
                invMSFOS = inv;
                invMElec = inv;
                invweight = weightcorr;
                if (treeInfo) treeSFOS->Fill();
                if (treeInfo) treeElec->Fill();
                if(n_bJet>=2){
		            for(std::vector<pat::Jet>::const_iterator bjet_i = bJets.begin(); bjet_i!=bJets.end(); ++bjet_i){
                	    for(std::vector<pat::Jet>::const_iterator bjet_j = bJets.begin(); bjet_j!=bJets.end();++bjet_j){
			                if(bjet_i<bjet_j){hInvMbbllSFOS[process]->Fill((bjet_i->p4()+bjet_j->p4()+ele_i->p4()+ele_j->p4()).M(),weight);}
			            }
                    }
		        }
            }        
       	    if( ele_i->charge()==ele_j->charge()&&ele_i<ele_j){hInvMElectronSS[process]->Fill((ele_i->p4()+ele_j->p4()).M(),weight);}   
        }
    }
    
    //Loop over taus 
    int n_Taus = 0;
    float tauPt = 0.;
    for (std::vector<pat::Tau>::const_iterator tau_i = taus->begin(); tau_i != taus->end(); ++tau_i){
        if (debug) std::cout <<"tau eta = "<< tau_i->eta() << std::endl;
        tauPt += tau_i->pt();
        ++numTotTaus;
	    ++n_Taus;
   	    TauMonitor(&(*tau_i),n_Taus,weight,general); 
        if(mcInfo){TauMonitor(&(*tau_i),n_Taus,weight,GetLeptKind(&(*tau_i)));}  
        
        //Invariant mass plots
	    //Tau pairs
   	    for (std::vector<pat::Tau>::const_iterator tau_j = taus->begin(); tau_j != taus->end(); ++tau_j){
       	    if( tau_i->charge()*tau_j->charge() < 0&&tau_i<tau_j){
	        if (debug) std::cout <<"Invariant Mass tau+ tau-: " << (tau_i->p4()+tau_j->p4()).M() <<std::endl;
                inv = (tau_i->p4()+tau_j->p4()).M();
                iso = CalcIso(*tau_i)+CalcIso(*tau_j);
                TauInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weight,general);
            }
       	    if( tau_i->charge()*tau_j->charge()>0&&tau_i<tau_j){hInvMTauSS[process]->Fill((tau_i->p4()+tau_j->p4()).M(),weight);}   
        }
    }
 
    //Global alphaT from all objects
    //double alphaT = 0; 
    //if(FindMinCombo(objects,mlaMinDpt,mlbMinDpt)){alphaT =  CalcalphaT(objects,mlaMinDpt,mlbMinDpt);}
    //halphaT[process]->Fill(alphaT,weight);
    //Muon multiplicity 
    hMuonMult[general]->Fill(n_Muons,weight);
    hMuonSumPt[general]->Fill(muonPt,weight);
    //Electron multiplicity 
    hElectronMult[general]->Fill(n_Electrons,weight);
    hElectronSumPt[general]->Fill(elePt,weight);
    //Tau multipilcity
    hTauMult[general]->Fill(n_Taus,weight);
    hTauSumPt[general]->Fill(tauPt,weight);
    //Lepton multiplicity
    hLightLeptonMult[general]->Fill(n_Electrons+n_Muons, weight);
    hLeptonMult[general]->Fill(n_Electrons+n_Muons+n_Taus, weight);
}

//Fill all muon inv mass related quantities
void DiLeptonHistograms::MuonInvMonitor(const double inv,const double MET,const double HT,
        const double et4Jets, const double etFourthJet, const double iso, double weight, const int process){
    hInvMass[process]->Fill(inv,weight);
    hInvMMuon[process]->Fill(inv,weight);
    hInvMSFOS[process]->Fill(inv,weight);
    h2dMETInvMassSFOS[process]->Fill(inv,MET,weight);
    h2dEtJetsInvMassSFOS[process]->Fill(inv,et4Jets,weight);
    h2dHTInvMassSFOS[process]->Fill(inv,HT,weight);
    h2dEtJet4InvMassSFOS[process]->Fill(inv,etFourthJet,weight);
    h2dIsoInvMassSFOS[process]->Fill(inv,iso,weight);
}

//Fill all electron inv mass related quantities
void DiLeptonHistograms::ElectronInvMonitor(const double inv,const double MET,const double HT,
        const double et4Jets, const double etFourthJet, const double iso, double weight, const int process){
    hInvMass[process]->Fill(inv,weight);
    hInvMElectron[process]->Fill(inv,weight);
    hInvMSFOS[process]->Fill(inv,weight);
    h2dMETInvMassSFOS[process]->Fill(inv,MET,weight);
    h2dEtJetsInvMassSFOS[process]->Fill(inv,et4Jets,weight);
    h2dHTInvMassSFOS[process]->Fill(inv,HT,weight);
    h2dEtJet4InvMassSFOS[process]->Fill(inv,etFourthJet,weight);
    h2dIsoInvMassSFOS[process]->Fill(inv,iso,weight);
}

void DiLeptonHistograms::TauInvMonitor(const double inv,const double MET,const double HT,
        const double et4Jets, const double etFourthJet, const double iso, double weight, const int process){
    hInvMTau[process]->Fill(inv,weight);
}

//Fill all OFOS inv mass related quantities
void DiLeptonHistograms::OFOSInvMonitor(const double inv,const double MET,const double HT,
        const double et4Jets, const double etFourthJet, const double iso, double weight, const int process){
    hInvMass[process]->Fill(inv,-weight);
    hInvMOFOS[process]->Fill(inv,weight);
    h2dMETInvMassOFOS[process]->Fill(inv,MET,weight);
    h2dEtJetsInvMassOFOS[process]->Fill(inv,et4Jets,weight);
    h2dHTInvMassOFOS[process]->Fill(inv,HT,weight);
    h2dEtJet4InvMassOFOS[process]->Fill(inv,etFourthJet,weight);
    h2dIsoInvMassOFOS[process]->Fill(inv,iso,weight);
}

//Fill all the trigger related quantities
void DiLeptonHistograms::TriggerMonitor(const edm::Handle< edm::TriggerResults>& trigger, const edm::TriggerNames & names, double weight, const int process){
    //const int nFilters = trigger->size();
    //int nAccept = 0;

    std::vector< std::string > hlNames;
    for (edm::TriggerNames::Strings::const_iterator j = names.triggerNames().begin(); j !=names.triggerNames().end(); ++j ) { 
        hlNames.push_back(*j);
    }
    hTrigger[process]->Fill(hlNames[0].c_str(),weight);
    //LogPrint("Trigger") << "Event" << "\n";
    for (unsigned int i=0; i<trigger->size()-1; ++i){
        hTrigger[process]->Fill(hlNames[i].c_str(),0);
   	if((*trigger)[i].accept()){
	    //++nAccept;
   	    hTrigger[process]->Fill(hlNames[i].c_str(),weight);
	    //LogPrint("Event")  << i << " : " << hlNames[i] ;
	}
    }
    //LogPrint("Event") << "Number of triggers: " << nFilters << "\n"
    //		     << "Number of accepted triggers: " << nAccept;
}

//Fill all electron related quantities
void DiLeptonHistograms::ElectronMonitor(const pat::Electron* electron,const int n_Electron, double weight, const int process){
    if(process==effcor){weight=getElectronWeight(electron);}
    //Electron base plot
    hElectronPt[process]->Fill(electron->pt(),weight);
    hElectronCharge[process]->Fill(electron->charge(),weight);
    if(n_Electron == 1){
        hElectron1Pt[process]->Fill(electron->pt(),weight);
        hElectron1Eta[process]->Fill(electron->eta(),weight);
    }
    if(n_Electron == 2){
        hElectron2Pt[process]->Fill(electron->pt(),weight);
        hElectron2Eta[process]->Fill(electron->eta(),weight);
    }
    hElectronEta[process]->Fill(electron->eta(),weight);
    hElectronPhi[process]->Fill(electron->phi(),weight);
    hElectronEtaPhi[process]->Fill(electron->eta(),electron->phi(),weight);
    //Electron isolation
    double IsoValue = CalcIso(*electron);
    double pfIsoValue = CalcPfIso(*electron);
    hElectronIso[process]->Fill(IsoValue,weight);
    hElectronTrackIso[process]->Fill(electron->trackIso(),weight);
    hElectronEcalIso[process]->Fill(electron->ecalIso(),weight);
    hElectronHcalIso[process]->Fill(electron->hcalIso(),weight);
    hElectronPfIso[process]->Fill(pfIsoValue,weight);
    hElectronChargedHadronIso[process]->Fill(electron->chargedHadronIso(),weight);
    hElectronPhotonIso[process]->Fill(electron->photonIso(),weight);
    hElectronNeutralHadronIso[process]->Fill(electron->neutralHadronIso(),weight);
    h2dElectronEtaPt[process]->Fill(electron->pt(),electron->eta(),weight);
  	
    double eOverP = electron->eSuperClusterOverP();
    //double eSeed = electron->superCluster()->seed()->energy();
    double pin  = electron->trackMomentumAtVtx().R();   
    //double eSeedOverPin = eSeed/pin; 
    double pout = electron->trackMomentumOut().R(); 
    double fBrem = (pin-pout)/pin;
    double e25Maxoe55 = electron->e2x5Max()/electron->e5x5();  
    double e15oe55 = electron->e1x5()/electron->e5x5();  
    double sigmaEtaEta = electron->sigmaEtaEta();  
    double sigmaIetaIeta = electron->sigmaIetaIeta();  
      
    double hOverE = electron->hadronicOverEm();
    //double sigmaee = sqrt((shapeRef)->covEtaEta());
    double deltaPhiIn = electron->deltaPhiSuperClusterTrackAtVtx();
    double deltaEtaIn = electron->deltaEtaSuperClusterTrackAtVtx();
  
    hElectronMva[process]->Fill(electron->mva(),weight); 
    hElectronEoverP[process]->Fill(eOverP,weight); 
    hElectronfBrem[process]->Fill(fBrem,weight); 
    hElectronHoverE[process]->Fill(hOverE,weight); 
    hElectrondeltaPhiIn[process]->Fill(deltaPhiIn,weight); 
    hElectrondeltaEtaIn[process]->Fill(deltaEtaIn,weight);
    hElectrone25Maxoe55[process]->Fill(e25Maxoe55,weight);
    hElectrone15oe55[process]->Fill(e15oe55,weight);
    hElectronsigmaEtaEta[process]->Fill(sigmaEtaEta,weight);
    hElectronsigmaIetaIeta[process]->Fill(sigmaIetaIeta,weight);
    if(!electron->gsfTrack().isNull()){
        //impact parameter
        int lostHits = electron->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
        double d0tobs = electron->gsfTrack()->dxy(bs);
        double d0topv = electron->gsfTrack()->dxy(pv);
        double d0Err = electron->gsfTrack()->dxyError();
        hElectronLostHits[process]->Fill(lostHits,weight);
        hElectrond0Pv[process]->Fill(d0topv,weight);
        hElectrond0SigPv[process]->Fill(d0topv/d0Err,weight);
        hElectronPfIsod0Pv[process]->Fill(pfIsoValue,d0topv,weight);
        hElectrond0Bs[process]->Fill(d0tobs,weight);
        hElectrond0SigBs[process]->Fill(d0tobs/d0Err,weight);
        hElectronPfIsod0Bs[process]->Fill(pfIsoValue,d0tobs,weight);
    }
    if (mcInfo){
        if(electron->genLepton()){
            h2dMatchedElectronEtaPt[process]->Fill(electron->pt(),electron->eta(),weight);
        }
    } 
}


//Fill all muon related quantities
void DiLeptonHistograms::MuonMonitor(const pat::Muon* muon,const int n_Muon, double weight, const int process){
    if(process==effcor){weight=getMuonWeight(muon);}
    h2dMuonEtaPt[process]->Fill(muon->pt(),muon->eta(),weight);
    //Muon base plots
    hMuonPt[process]->Fill(muon->pt(),weight);
    hMuonCharge[process]->Fill(muon->charge(),weight);
	
    if(n_Muon == 1){
        hMuon1Pt[process]->Fill(muon->pt(),weight);
        hMuon1Eta[process]->Fill(muon->eta(),weight);
    }
    if(n_Muon == 2){
        hMuon2Pt[process]->Fill(muon->pt(),weight);
        hMuon2Eta[process]->Fill(muon->eta(),weight);
    }
    hMuonEta[process]->Fill(muon->eta(),weight);
    hMuonPhi[process]->Fill(muon->phi(),weight);
    hMuonEtaPhi[process]->Fill(muon->eta(),muon->phi(),weight);
    
    //Muon isolation
    double IsoValue = CalcIso(*muon);
    double pfIsoValue = CalcPfIso(*muon);
    hMuonIso[process]->Fill(IsoValue,weight);
    hMuonTrackIso[process]->Fill(muon->trackIso(),weight);
    hMuonEcalIso[process]->Fill(muon->ecalIso(),weight);
    hMuonHcalIso[process]->Fill(muon->hcalIso(),weight);
    hMuonPfIso[process]->Fill(pfIsoValue,weight);
    hMuonChargedHadronIso[process]->Fill(muon->chargedHadronIso(),weight);
    hMuonPhotonIso[process]->Fill(muon->photonIso(),weight);
    hMuonNeutralHadronIso[process]->Fill(muon->neutralHadronIso(),weight);

    if(!muon->innerTrack().isNull()){
        double d0tobs = muon->innerTrack()->dxy(bs);
        double d0topv = muon->innerTrack()->dxy(pv);
        double d0Err = muon->innerTrack()->dxyError();
    	hMuonChi2[process]->Fill(muon->innerTrack()->normalizedChi2(),weight);
	    hMuonnHits[process]->Fill(muon->innerTrack()->numberOfValidHits(),weight);
        hMuond0Pv[process]->Fill(d0topv,weight);
	    hMuond0SigPv[process]->Fill(d0topv/d0Err,weight);
        hMuonPfIsod0Pv[process]->Fill(pfIsoValue,d0topv,weight);
        hMuond0Bs[process]->Fill(d0tobs,weight);
	    hMuond0SigBs[process]->Fill(d0tobs/d0Err,weight);
        hMuonPfIsod0Bs[process]->Fill(pfIsoValue,d0tobs,weight);
    }
    if (mcInfo){
        if(muon->genLepton()){
            h2dMatchedMuonEtaPt[process]->Fill(muon->pt(),muon->eta(),weight);
        }
    }
}

void DiLeptonHistograms::JetMonitor(const pat::Jet* jet,const int n_Jet, double weight, const int process){
	if(n_Jet==1){
        hEtJet1[process]->Fill(jet->pt(),weight);
        hJet1Eta[process]->Fill(jet->eta(),weight);
    }
	if(n_Jet==2){
        hEtJet2[process]->Fill(jet->pt(),weight);
        hJet2Eta[process]->Fill(jet->eta(),weight);
	}
	if(n_Jet==3){
        hEtJet3[process]->Fill(jet->pt(),weight);
        hJet3Eta[process]->Fill(jet->eta(),weight);
	}
	if(n_Jet==4){
        hEtJet4[process]->Fill(jet->pt(),weight);
        hJet4Eta[process]->Fill(jet->eta(),weight);
	}
	hJetEt[process]->Fill(jet->pt(),weight);
	hJetEta[process]->Fill(jet->eta(),weight);
	hJetPhi[process]->Fill(jet->phi(),weight);
	hJetCHF[process]->Fill(jet->chargedHadronEnergyFraction(),weight);
	hJetNHF[process]->Fill(jet->neutralHadronEnergyFraction(),weight);
	hJetCEF[process]->Fill(jet->chargedEmEnergyFraction(),weight);
	hJetNEF[process]->Fill(jet->neutralEmEnergyFraction(),weight);
	hJetCM[process]->Fill(jet->chargedMultiplicity(),weight);
	hJetNConst[process]->Fill(jet->numberOfDaughters(),weight);
    h2dJetEtaPt[process]->Fill(jet->pt(),jet->eta(),weight);
    h2dJetEtaPhi[process]->Fill(jet->eta(),jet->phi(),weight);
    if (mcInfo){
        if(jet->genJet()){
            h2dMatchedJetEtaPt[process]->Fill(jet->pt(),jet->eta(),weight);
        }
    }
}

void DiLeptonHistograms::InitTauHistos( const pat::Tau& tau, const int process)
{
  std::vector< pat::Tau::IdPair  > tauIds = tau.tauIDs();
  unsigned int binNr = 1;
  hTauDiscriminators[process]->GetXaxis()->SetBinLabel(binNr,"None");
  assert( maxTauDiscriminators_ > tauIds.size());// std::cerr << "maxTauDiscriminators too small: "<< tauIds.size() << std::endl;
  for(std::vector< pat::Tau::IdPair  >::iterator it = tauIds.begin(); it != tauIds.end() && binNr <= maxTauDiscriminators_; ++it){
    binNr++;
    hTauDiscriminators[process]->GetXaxis()->SetBinLabel(binNr,(*it).first.c_str());
  }
  tauInitialized_[process] = true;
}

//Fill all tau related quantities
void DiLeptonHistograms::TauMonitor(const pat::Tau* tau,const int n_Tau, double weight, const int process){
    if(!tauInitialized_[process]) InitTauHistos(*tau, process);
    std::vector< pat::Tau::IdPair  > tauIds = tau->tauIDs();
    hTauDiscriminators[process]->Fill("None", weight);
    for(std::vector< pat::Tau::IdPair  >::iterator it = tauIds.begin(); it != tauIds.end(); ++it){
      if((*it).second > 0.5 ) // TODO make this configurable
      	hTauDiscriminators[process]->Fill( (*it).first.c_str(), weight);
    }
    if(mcInfo && tau->genJet() != 0 && tau->genParticle() != 0){
      hGenTauPt[process]->Fill(tau->genParticle()->pt(), weight);
      hGenTauVisPt[process]->Fill(tau->genJet()->pt(), weight);
    }
	
    h2dTauEtaPt[process]->Fill(tau->pt(),tau->eta(),weight);
    //Tau base plots
    hTauPt[process]->Fill(tau->pt(),weight);
    hTauCharge[process]->Fill(tau->charge(),weight);
	
    if(n_Tau == 1){
        hTau1Pt[process]->Fill(tau->pt(),weight);
        hTau1Eta[process]->Fill(tau->eta(),weight);
    }
    if(n_Tau == 2){
        hTau2Pt[process]->Fill(tau->pt(),weight);
        hTau2Eta[process]->Fill(tau->eta(),weight);
    }
    hTauEta[process]->Fill(tau->eta(),weight);
    hTauPhi[process]->Fill(tau->phi(),weight);
    hTauEtaPhi[process]->Fill(tau->eta(),tau->phi(),weight);
    
    //Tau isolation
    double IsoValue = CalcIso(*tau);
    double pfIsoValue = CalcPfIso(*tau);
    hTauIso[process]->Fill(IsoValue,weight);
    hTauTrackIso[process]->Fill(tau->trackIso(),weight);
    hTauEcalIso[process]->Fill(tau->ecalIso(),weight);
    hTauHcalIso[process]->Fill(tau->hcalIso(),weight);
    hTauPfIso[process]->Fill(pfIsoValue,weight);
    hTauChargedHadronIso[process]->Fill(tau->chargedHadronIso(),weight);
    hTauPhotonIso[process]->Fill(tau->photonIso(),weight);
    hTauNeutralHadronIso[process]->Fill(tau->neutralHadronIso(),weight);

    if (mcInfo){
        if(tau->genLepton()){
            h2dMatchedTauEtaPt[process]->Fill(tau->pt(),tau->eta(),weight);
        }
    }
}
 
//Event loop
//Gets all collections and calls Analysis
void DiLeptonHistograms::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
    double weight = externalWeight;

    //edm::LogPrint("Evt")  << "Run = " << iEvent.id().run() << ", Event = " << iEvent.id().event();
    //std::cout << "Run = " << iEvent.id().run() << ", Event = " << iEvent.id().event() << std::endl;

    //Beam spot
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByLabel(beamSpotSrc, beamSpotHandle);
    bs = beamSpotHandle->position();

    edm::Handle<reco::VertexCollection> Vertices;
    iEvent.getByLabel(primaryVertexSrc, Vertices);
    for (reco::VertexCollection::const_iterator it = Vertices->begin(); it != Vertices->end(); ++it) {
        pv = it->position();
        break;
    }
  
    // retrieve the PAT-objects
    //Muons
    edm::Handle< std::vector<pat::Muon> > muons;
    iEvent.getByLabel(muonSrc, muons);
  
    //Electrons
    edm::Handle< std::vector<pat::Electron> > electrons;
    iEvent.getByLabel(electronSrc, electrons);
    
    //Taus
    edm::Handle< std::vector<pat::Tau> > taus;
    iEvent.getByLabel(tauSrc, taus);
   
    //Jets
    edm::Handle< std::vector<pat::Jet> > jets;
    iEvent.getByLabel(jetSrc, jets);
   
    //MET
    edm::Handle< std::vector<pat::MET> > met;
    iEvent.getByLabel(metSrc, met);
  
    //Trigger 
    edm::Handle< edm::TriggerResults > trigger;
    iEvent.getByLabel("TriggerResults", trigger);
    const edm::TriggerNames & triggerNames = iEvent.triggerNames(*trigger);

    bool signal = false;
    ++numTotEvents;
    if (mcInfo){
        //MC gen Particle
        edm::Handle< std::vector<reco::GenParticle> > genParticles;
        iEvent.getByLabel(mcSrc, genParticles);
        signal = MCAnalysis(muons,electrons,genParticles,weight,general);
    }
    if (signal&&Signal_Analysis){   
        TriggerMonitor(trigger,triggerNames,weight,general);
        Analysis(muons, electrons, taus, jets, met, weight, general);
    }
    if(!Signal_Analysis){
        TriggerMonitor(trigger,triggerNames,weight,general);
        Analysis(muons, electrons, taus, jets, met, weight, general);
    }
     
}


//MC analysis of leptons
bool DiLeptonHistograms::MCAnalysis(const edm::Handle< std::vector<pat::Muon> >& muons, const edm::Handle< std::vector<pat::Electron> >& electrons, const edm::Handle< std::vector<reco::GenParticle> >& genParticles, double weight, const int process){
    bool signal = false;
    int pid = 0;
    float metx = 0;
    float mety = 0;
    for (std::vector<reco::GenParticle>::const_iterator p_i = genParticles->begin(); p_i != genParticles->end(); ++p_i){
        pid = abs(p_i->pdgId());
        if (pid==1000023){
            //std::cout << "Found neutralino with 3 daughters" << std::endl;
            std::vector<const reco::Candidate *> mc_electrons;
            std::vector<const reco::Candidate *> mc_muons;
            std::vector<const reco::Candidate *> mc_chi1;
            for(size_t j = 0; j < p_i->numberOfDaughters(); ++j){
                const reco::Candidate * d = p_i->daughter(j);
                //std::cout << d->pdgId() << std::endl;
 	        if (d->pt()>cut_MuonPt&&fabs(d->eta())<cut_MuonEta&& abs(d->pdgId())==11){mc_electrons.push_back(d);}     
 	        if (d->pt()>cut_MuonPt&&fabs(d->eta())<cut_MuonEta&& abs(d->pdgId())==13){mc_muons.push_back(d);}     
                if( abs(d->pdgId()) == 1000022){mc_chi1.push_back(d);}
            }
            if (mc_chi1.size()==1 && mc_electrons.size()==2){
		hInvMassMC[process]->Fill((mc_electrons[0]->p4()+mc_electrons[1]->p4()).M(),weight);
                signal=true;
            }
            if (mc_chi1.size()==1 && mc_muons.size()==2){
		hInvMassMC[process]->Fill((mc_muons[0]->p4()+mc_muons[1]->p4()).M(),weight);
                signal=true;
            }
        }
        if (pid== 22 || pid == 23){
            std::vector<const reco::Candidate *> mc_electrons;
            std::vector<const reco::Candidate *> mc_muons;
            for(size_t j = 0; j < p_i->numberOfDaughters(); ++j){
                //if (j == 0){std::cout << "Found Z with daughters" << std::endl;}
                const reco::Candidate * d = p_i->daughter(j);
                //std::cout << d->pdgId() << std::endl;
                //std::cout << d->status() << std::endl;
 	        if (d->pt()>cut_MuonPt&&fabs(d->eta())<cut_MuonEta&& abs(d->pdgId())==11){mc_electrons.push_back(d);}     
 	        if (d->pt()>cut_MuonPt&&fabs(d->eta())<cut_MuonEta&& abs(d->pdgId())==13){mc_muons.push_back(d);}     
            }
            if (mc_electrons.size()==2&&mc_electrons[0]->pdgId()*mc_electrons[1]->pdgId()<0){
		hInvMassZMC[process]->Fill((mc_electrons[0]->p4()+mc_electrons[1]->p4()).M(),weight);
                //std::cout << mc_electrons[0]->pdgId() << std::endl;
                //std::cout << mc_electrons[1]->pdgId() << std::endl;
            }
            if (mc_muons.size()==2&&mc_muons[0]->pdgId()*mc_muons[1]->pdgId()<0){
		hInvMassZMC[process]->Fill((mc_muons[0]->p4()+mc_muons[1]->p4()).M(),weight);
                //std::cout << mc_muons[0]->pdgId() << std::endl;
                //std::cout << mc_muons[1]->pdgId() << std::endl;
            }
        }
        if (p_i->status()==1){
	    //Muons (13) with status 1
 	    if (p_i->pt()>cut_MuonPt&&fabs(p_i->eta())<cut_MuonEta&&pid==13){
		    hGenMuonPt[process]->Fill(p_i->pt(),weight);
		    hGenMuonEta[process]->Fill(p_i->eta(),weight);
            h2dGenMuonEtaPt[process]->Fill(p_i->pt(),p_i->eta(),weight);
 	    }
	    //Electron (11) with status 1
 	    if (p_i->pt()>cut_ElectronPt&&fabs(p_i->eta())<cut_ElectronEta&&abs(p_i->pdgId())==11&&p_i->status()==1){
	        hGenElectronPt[process]->Fill(p_i->pt(),weight);
		    hGenElectronEta[process]->Fill(p_i->eta(),weight);
    		h2dGenElectronEtaPt[process]->Fill(p_i->pt(),p_i->eta(),weight);
            }
 	    if ( pid == 12 || pid == 13 || pid == 14 || pid == 16 || 
		 pid == 1000022 || pid == 2000012 || pid == 2000014 ||
		 pid == 2000016 || pid == 1000039 || pid == 5000039 ||
		 pid == 4000012 || pid == 9900012 || pid == 9900014 ||
		 pid == 9900016 || pid == 39 ){
	        metx += p_i->px(); //TODO define met using et
		    mety += p_i->py();
	    }
        }
    }
    hMissingETmc[process]->Fill(sqrt(metx*metx+mety*mety),weight);

    return signal;
}

//PrintStatics
void DiLeptonHistograms::PrintStatistics(void)
{ 
    edm::LogPrint("Summary")  << "Total number of events processed = " << numTotEvents << "\n"
			  << "Accepted good objects: " << "\n"
 			  << "Total number of electrons = " << numTotElectrons 
				 << " per event = " << (float)numTotElectrons / (float)numTotEvents << "\n"
			  << "Total number of muons     = " << numTotMuons
				<< " per event = " << (float)numTotMuons / (float)numTotEvents << "\n"
			  << "Total number of taus     = " << numTotTaus
				<< " per event = " << (float)numTotTaus / (float)numTotEvents << "\n"
			  << "Total number of jets      = " << numTotJets
			    << " per event = " << (float)numTotJets / (float)numTotEvents << "\n";
}

DEFINE_FWK_MODULE(DiLeptonHistograms);
