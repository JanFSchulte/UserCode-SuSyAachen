import FWCore.ParameterSet.Config as cms
# central place to edit used JEC so it is always consistent
# don't forget to also add the DB files to additional import files for crab


# whether to load from GT or from DB file
dataUseDB = {}
dataUseDB["2016"] = True
dataUseDB["2017"] = True
dataUseDB["2018"] = False

# db name, can be empty string if not used
dataDB = {}
dataDB["2016"] = "Summer16_07Aug2017All_V11_DATA"
dataDB["2017"] = "Fall17_17Nov2017_V32_94X_DATA"
dataDB["2018"] = "Autumn18_RunABCD_V8_DATA"

# whether to load from GT or from DB file
mcUseDB = {}
mcUseDB["2016"] = True
mcUseDB["2017"] = True
mcUseDB["2018"] = False

# db name, can be empty string if not used
mcDB = {}
mcDB["2016"] = "Summer16_07Aug2017_V11_MC"
mcDB["2017"] = "Fall17_17Nov2017_V32_94X_MC"
mcDB["2018"] = "Autumn18_V8_MC"
