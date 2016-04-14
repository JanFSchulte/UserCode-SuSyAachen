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
											std::string FastSimScaleFactorFile_mu, 
											std::string FastSimScaleFactorFile_ele, 
											std::string FastSimScaleFactorHisto_mu, 
											std::string FastSimScaleFactorHisto_ele
											)
{
	
	TFile FastSimScaleFactorFile_mu_(FastSimScaleFactorFile_mu.c_str());
	TFile FastSimScaleFactorFile_ele_(FastSimScaleFactorFile_ele.c_str());
	
	FastSimScaleFactorHisto_mu_ = (TH3F*)FastSimScaleFactorFile_mu_.Get( FastSimScaleFactorHisto_mu.c_str() );
	FastSimScaleFactorHisto_ele_ = (TH3F*)FastSimScaleFactorFile_ele_.Get( FastSimScaleFactorHisto_ele.c_str() );
}

