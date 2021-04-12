import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2500)

#options = VarParsing.VarParsing('standard')
options = VarParsing.VarParsing('analysis')

#options.register('inputFile',
#        "~/",
#        VarParsing.VarParsing.multiplicity.singleton,
#        VarParsing.VarParsing.varType.string,
#        "File containing a list of the EXACT location of the output file  (default = ~/)"
#        )


options.parseArguments()
#options.inputFile = 'root://eoscms//' + options.inputFile
#print(options.inputFile)
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
     #       'root://cms-xrd-global.cern.ch///store/mc/RunIISummer19UL18RECO/DoublePhoton_FlatPt-5To300/AODSIM/FlatPU0to70RAW_106X_upgrade2018_realistic_v11_L1v1_ext1-v2/00000/6F7BDF5D-2A6A-7342-993F-E8DBABF920C8.root'
   #         'root://cms-xrd-global.cern.ch///store/mc/RunIISummer19UL18RECO/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/250000/35BAE18A-E8EA-254E-B693-9212D3725212.root'
     #       'root://cms-xrd-global.cern.ch///store/mc/RunIISummer19UL18RECO/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/250000/26B3A9F9-290D-6346-8862-2BE71645C57D.root'
                options.inputFiles
                )
                            )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.AOD
switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff','RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1p1_cff','RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff']

for idmod in my_id_modules:
        setupAllVIDIdsInModule(process, idmod, setupVIDPhotonSelection)

process.nTuplelize = cms.EDAnalyzer('Photon_RefinedRecHit_NTuplizer',
        rhoFastJet = cms.InputTag("fixedGridRhoFastjetAll"),
        photons = cms.InputTag("gedPhotons"),
        genParticles = cms.InputTag("genParticles"),
        #MVA Based Id
        eleMediumIdMap = cms.InputTag(""),
        eleTightIdMap = cms.InputTag("")
	)


process.TFileService = cms.Service("TFileService",
     fileName = cms.string("nTupleMC.root"),
#     fileName = cms.string("Tree_Gamma_ABCD.root"),
      closeFileFast = cms.untracked.bool(True)
  )


process.p = cms.Path(process.egmPhotonIDSequence*process.nTuplelize)

