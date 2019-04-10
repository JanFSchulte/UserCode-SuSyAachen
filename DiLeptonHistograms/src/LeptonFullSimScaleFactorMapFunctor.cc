#include "SuSyAachen/DiLeptonHistograms/interface/LeptonFullSimScaleFactorMapFunctor.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH3.h"
#include "TH3F.h"
#include "TFile.h"
#include <string>

#include <DataFormats/PatCandidates/interface/Lepton.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>


void
LeptonFullSimScaleFactorMapFunctor::getFullSimScaleFactorHistos(std::string dataMCScaleFactorFile_mu_ID, 
                                                                                        std::string dataMCScaleFactorFile_mu_Iso, 
                                                                                        //std::string dataMCScaleFactorFile_mu_Impact, 
                                                                                        //std::string dataMCScaleFactorFile_mu_Track, 
                                                                                        //std::string dataMCScaleFactorFile_mu_SIP3D, 
                                                                                        
                                                                                        std::string dataMCScaleFactorFile_ele, 
                                                                                        std::string dataMCScaleFactorFile_ele_Reco, 
                                                                                        //std::string dataMCScaleFactorFile_ele_Track, 
                                                                                        
                                                                                        std::string dataMCScaleFactorHisto_mu_ID, 
                                                                                        std::string dataMCScaleFactorHisto_mu_Iso, 
                                                                                        //std::string dataMCScaleFactorHisto_mu_Impact, 
                                                                                        //std::string dataMCScaleFactorHisto_mu_Track, 
                                                                                        //std::string dataMCScaleFactorHisto_mu_SIP3D, 
                                                                                        
                                                                                        std::string dataMCScaleFactorHisto_ele_Reco, 
                                                                                        std::string dataMCScaleFactorHisto_ele_ID, 
                                                                                        std::string dataMCScaleFactorHisto_ele_Iso,
                                                                                        std::string dataMCScaleFactorHisto_ele_ConvHit
                                                                                        //std::string dataMCScaleFactorHisto_ele_Track
                                                                                        )
{
        TFile dataMCScaleFactorFile_mu_ID_(dataMCScaleFactorFile_mu_ID.c_str());
        TFile dataMCScaleFactorFile_mu_Iso_(dataMCScaleFactorFile_mu_Iso.c_str());
        //TFile dataMCScaleFactorFile_mu_Impact_(dataMCScaleFactorFile_mu_Impact.c_str());
        //TFile dataMCScaleFactorFile_mu_Track_(dataMCScaleFactorFile_mu_Track.c_str());
        //TFile dataMCScaleFactorFile_mu_SIP3D_(dataMCScaleFactorFile_mu_SIP3D.c_str());
        
        TFile dataMCScaleFactorFile_ele_(dataMCScaleFactorFile_ele.c_str());
        TFile dataMCScaleFactorFile_ele_Reco_(dataMCScaleFactorFile_ele_Reco.c_str());
        //TFile dataMCScaleFactorFile_ele_Track_(dataMCScaleFactorFile_ele_Track.c_str());
        
        dataMCScaleFactorHisto_mu_ID_ = (TH2F*)dataMCScaleFactorFile_mu_ID_.Get( dataMCScaleFactorHisto_mu_ID.c_str() );
        dataMCScaleFactorHisto_mu_Iso_ = (TH2F*)dataMCScaleFactorFile_mu_Iso_.Get( dataMCScaleFactorHisto_mu_Iso.c_str() );
        //dataMCScaleFactorHisto_mu_Impact_ = (TH2F*)dataMCScaleFactorFile_mu_Impact_.Get( dataMCScaleFactorHisto_mu_Impact.c_str() );
        //dataMCScaleFactorHisto_mu_SIP3D_ = (TH2F*)dataMCScaleFactorFile_mu_SIP3D_.Get( dataMCScaleFactorHisto_mu_SIP3D.c_str() );
        //dataMCScaleFactorHisto_mu_Track_ = (TH1F*)dataMCScaleFactorFile_mu_Track_.Get( dataMCScaleFactorHisto_mu_Track.c_str() );
        
        dataMCScaleFactorHisto_ele_Reco_ = (TH2F*)dataMCScaleFactorFile_ele_Reco_.Get( dataMCScaleFactorHisto_ele_Reco.c_str() );
        dataMCScaleFactorHisto_ele_ID_ = (TH2F*)dataMCScaleFactorFile_ele_.Get( dataMCScaleFactorHisto_ele_ID.c_str() );
        dataMCScaleFactorHisto_ele_Iso_ = (TH2F*)dataMCScaleFactorFile_ele_.Get( dataMCScaleFactorHisto_ele_Iso.c_str() );
        dataMCScaleFactorHisto_ele_ConvHit_ = (TH2F*)dataMCScaleFactorFile_ele_.Get( dataMCScaleFactorHisto_ele_ConvHit.c_str() );
        //dataMCScaleFactorHisto_ele_Track_ = (TH2F*)dataMCScaleFactorFile_ele_Track_.Get( dataMCScaleFactorHisto_ele_Track.c_str() );
}

