from CRABClient.UserUtilities import config
import FWCore.ParameterSet.Config as cms

config = config()

config.General.requestName = ''
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'

config.JobType.outputFiles= [
    'nTupleMC.root'
]

config.JobType.psetName = 'Photon_RecHit_AOD_cfg.py'

config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 2750
#config.JobType.sendPythonFolder = True

#config.Data.inputDataset  = '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-RECOSIMstep_94X_mc2017_realistic_v10-v1/MINIAODSIM'
config.Data.inputDataset  = '/DoublePhoton_FlatPt-5To300/RunIIWinter19PFCalibDR-2018ConditionsFlatPU0to70ECALGT_105X_upgrade2018_realistic_IdealEcalIC_v4-v1/AODSIM'
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'FileBased' #'LumiBased'
config.Data.unitsPerJob   = 10 #30000
config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples/rchatter/2018_UL_Photon_MC_v3'
config.Data.publication   = False

# GRID
config.Site.storageSite   =  'T2_CH_CERN'
#config.Site.whitelist     = ['T3_CH_CERN_CAF']
#config.Site.blacklist     = ['T1_US_FNAL','T2_UA_KIPT','T2_UK_London_Brunel','T2_CH_CSCS','T2_US_*']
#config.Site.ignoreGlobalBlacklist = True
