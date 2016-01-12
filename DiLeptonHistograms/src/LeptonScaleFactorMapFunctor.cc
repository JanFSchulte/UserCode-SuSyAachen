#include "SuSyAachen/DiLeptonHistograms/interface/LeptonScaleFactorMapFunctor.h"
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
LeptonScaleFactorMapFunctor::getScaleFactorHistos(std::string dataMCScaleFactorFile_mu_ID, 
											std::string dataMCScaleFactorFile_mu_Iso, 
											std::string dataMCScaleFactorFile_ele, 
											std::string dataMCScaleFactorHisto_mu_ID, 
											std::string dataMCScaleFactorHisto_mu_Iso, 
											std::string dataMCScaleFactorHisto_ele_ID, 
											std::string dataMCScaleFactorHisto_ele_Iso, 
											std::string FastSimScaleFactorFile_mu, 
											std::string FastSimScaleFactorFile_ele, 
											std::string FastSimScaleFactorHisto_mu, 
											std::string FastSimScaleFactorHisto_ele
											)
{
	TFile dataMCScaleFactorFile_mu_ID_(dataMCScaleFactorFile_mu_ID.c_str());
	TFile dataMCScaleFactorFile_mu_Iso_(dataMCScaleFactorFile_mu_Iso.c_str());
	TFile dataMCScaleFactorFile_ele_(dataMCScaleFactorFile_ele.c_str());
	
	TFile FastSimScaleFactorFile_mu_(FastSimScaleFactorFile_mu.c_str());
	TFile FastSimScaleFactorFile_ele_(FastSimScaleFactorFile_ele.c_str());
	
	dataMCScaleFactorHisto_mu_ID_ = (TH2F*)dataMCScaleFactorFile_mu_ID_.Get( dataMCScaleFactorHisto_mu_ID.c_str() );
	dataMCScaleFactorHisto_mu_Iso_ = (TH2F*)dataMCScaleFactorFile_mu_Iso_.Get( dataMCScaleFactorHisto_mu_Iso.c_str() );
	dataMCScaleFactorHisto_ele_ID_ = (TH2F*)dataMCScaleFactorFile_ele_.Get( dataMCScaleFactorHisto_ele_ID.c_str() );
	dataMCScaleFactorHisto_ele_Iso_ = (TH2F*)dataMCScaleFactorFile_ele_.Get( dataMCScaleFactorHisto_ele_Iso.c_str() );
	
	FastSimScaleFactorHisto_mu_ = (TH3F*)FastSimScaleFactorFile_mu_.Get( FastSimScaleFactorHisto_mu.c_str() );
	FastSimScaleFactorHisto_ele_ = (TH3F*)FastSimScaleFactorFile_ele_.Get( FastSimScaleFactorHisto_ele.c_str() );
}

