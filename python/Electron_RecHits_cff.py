import FWCore.ParameterSet.Config as cms

def setupElectronNtuplizer(process, dataFormat, isMC):
    if process=="AOD":
        if isMC: process.load("EM_GNN_ID.EM_Skimmer.Electron_RecHits_AODSIM_cfi")
        else: process.load("EM_GNN_ID.EM_Skimmer.Electron_RecHits_AOD_cfi")
    elif process=="MINIAOD":
        if isMC: process.load("EM_GNN_ID.EM_Skimmer.Electron_RecHits_MINIAOD_cfi")
        else: process.load("EM_GNN_ID.EM_Skimmer.Electron_RecHits_MINIAODSIM_cfi")
