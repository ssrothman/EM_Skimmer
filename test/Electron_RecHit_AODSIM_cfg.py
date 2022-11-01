import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2500)

options = VarParsing.VarParsing('standard')

options.register('inputFile',
        "~/",
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.string,
        "File containing a list of the EXACT location of the output file  (default = ~/)"
        )


options.parseArguments()
options.inputFile = 'root://eoscms//' + options.inputFile
print(options.inputFile)
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            #'/store/mc/RunIISummer20UL18RECO/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/260005/C9B7018E-9032-E84C-9865-A7141C895350.root'
            #'root://cms-xrd-global.cern.ch:1094//store/mc/RunIISummer20UL18RECO/GluGluHToGG_M-125_TuneCP5_13TeV-powheg-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/40000/418EDE70-2DF3-714A-8436-71F580A9ED86.root'
            'root://cms-xrd-global.cern.ch:1094//store/mc/RunIISummer20UL18RECO/DoubleElectron_Pt-1To300-gun/AODSIM/FlatPU0to70EdalIdealGT_EdalIdealGT_106X_upgrade2018_realistic_v11_L1v1_EcalIdealIC-v2/270000/10D911D6-79F1-0642-8B23-0A3FC578DB82.root'
#                options.inputFile
                ),
                inputCommands=cms.untracked.vstring(
                'drop recoTrackExtrasedmAssociation_muonReducedTrackExtras_*_*'
                )
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.AOD
switchOnVIDElectronIdProducer(process, dataFormat)


# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff']


for idmod in my_id_modules:
        setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)


process.nTuplelize = cms.EDAnalyzer('Electron_RefinedRecHit_NTuplizer',
        rhoFastJet = cms.InputTag("fixedGridRhoAll"),
        electrons = cms.InputTag("gedGsfElectrons","","RECO"),
        photons = cms.InputTag("photons"),
        genParticles = cms.InputTag("genParticles"),
        refinedCluster = cms.bool(False),
        miniAODRun = cms.bool(False),
        isMC = cms.bool(True), 
        #MVA Based Id
	eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose"),
        eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium"),
        eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight")
        )

process.TFileService = cms.Service("TFileService",
     fileName = cms.string("ElectronRecHits_ntuple.root"),
     closeFileFast = cms.untracked.bool(True)
  )

process.p = cms.Path(process.egmGsfElectronIDSequence*process.nTuplelize)
