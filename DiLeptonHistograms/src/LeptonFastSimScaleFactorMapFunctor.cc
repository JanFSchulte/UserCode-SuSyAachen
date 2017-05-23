#include "SuSyAachen/DiLeptonHistograms/interface/LeptonFastSimScaleFactorMapFunctor.h"
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
LeptonFastSimScaleFactorMapFunctor::getFastSimScaleFactorHistos(
											std::string FastSimScaleFactorFile_mu_Id, 
											std::string FastSimScaleFactorFile_mu_Iso, 
											std::string FastSimScaleFactorFile_mu_Impact, 
											std::string FastSimScaleFactorFile_mu_SIP3D,
											 
											std::string FastSimScaleFactorFile_ele_Id, 
											std::string FastSimScaleFactorFile_ele_Iso, 
											std::string FastSimScaleFactorFile_ele_ConvHit, 
											
											std::string FastSimScaleFactorHisto_mu_Id, 
											std::string FastSimScaleFactorHisto_mu_Iso, 
											std::string FastSimScaleFactorHisto_mu_Impact, 
											std::string FastSimScaleFactorHisto_mu_SIP3D, 
											
											std::string FastSimScaleFactorHisto_ele_Id,
											std::string FastSimScaleFactorHisto_ele_Iso,
											std::string FastSimScaleFactorHisto_ele_ConvHit
											)
{
	
	TFile FastSimScaleFactorFile_mu_Id_(FastSimScaleFactorFile_mu_Id.c_str());
	TFile FastSimScaleFactorFile_mu_Iso_(FastSimScaleFactorFile_mu_Iso.c_str());
	TFile FastSimScaleFactorFile_mu_Impact_(FastSimScaleFactorFile_mu_Impact.c_str());
	TFile FastSimScaleFactorFile_mu_SIP3D_(FastSimScaleFactorFile_mu_SIP3D.c_str());
	
	TFile FastSimScaleFactorFile_ele_Id_(FastSimScaleFactorFile_ele_Id.c_str());
	TFile FastSimScaleFactorFile_ele_Iso_(FastSimScaleFactorFile_ele_Iso.c_str());
	TFile FastSimScaleFactorFile_ele_ConvHit_(FastSimScaleFactorFile_ele_ConvHit.c_str());
	
	FastSimScaleFactorHisto_mu_Id_ = (TH2D*)FastSimScaleFactorFile_mu_Id_.Get( FastSimScaleFactorHisto_mu_Id.c_str() );
	FastSimScaleFactorHisto_mu_Iso_ = (TH2D*)FastSimScaleFactorFile_mu_Iso_.Get( FastSimScaleFactorHisto_mu_Iso.c_str() );
	FastSimScaleFactorHisto_mu_Impact_ = (TH2D*)FastSimScaleFactorFile_mu_Impact_.Get( FastSimScaleFactorHisto_mu_Impact.c_str() );
	FastSimScaleFactorHisto_mu_SIP3D_ = (TH2D*)FastSimScaleFactorFile_mu_SIP3D_.Get( FastSimScaleFactorHisto_mu_SIP3D.c_str() );
	
	FastSimScaleFactorHisto_ele_Id_ = (TH2D*)FastSimScaleFactorFile_ele_Id_.Get( FastSimScaleFactorHisto_ele_Id.c_str() );
	FastSimScaleFactorHisto_ele_Iso_ = (TH2D*)FastSimScaleFactorFile_ele_Iso_.Get( FastSimScaleFactorHisto_ele_Iso.c_str() );
	FastSimScaleFactorHisto_ele_ConvHit_ = (TH2D*)FastSimScaleFactorFile_ele_ConvHit_.Get( FastSimScaleFactorHisto_ele_ConvHit.c_str() );
}

