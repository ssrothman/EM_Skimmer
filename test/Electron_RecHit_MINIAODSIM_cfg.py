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
            '/store/mc/RunIISummer20UL18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/260000/8D4B008D-BE14-DD47-8949-C4640525DA1F.root'
#                options.inputFile
                ),
                inputCommands=cms.untracked.vstring(
                    'keep *',
                    'drop recoTrackExtrasedmAssociation_muonReducedTrackExtras_*_*'
                )
                            )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)


# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff']


for idmod in my_id_modules:
        setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)

process.nTuplelize = cms.EDAnalyzer('Electron_RefinedRecHit_MiniAOD_NTuplizer',
        ebRecHits = cms.InputTag("reducedEgamma","reducedEBRecHits","PAT"),
        eeRecHits = cms.InputTag("reducedEgamma","reducedEERecHits","PAT"),
        esRecHits = cms.InputTag("reducedEgamma","reducedESRecHits","PAT"),
        rhoFastJet = cms.InputTag("fixedGridRhoAll"),
        electrons = cms.InputTag("slimmedElectrons"),
        photons = cms.InputTag("slimmedPhotons"),
        genParticles = cms.InputTag("packedGenParticles"),
        refinedCluster = cms.bool(True),
        isMC = cms.bool(True), 
        miniAODRun = cms.bool(True),
        #MVA Based Id
	eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose"),
        eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium"),
        eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight"),
        egammaDeltaRMatch = cms.double(0.03),
        )


process.TFileService = cms.Service("TFileService",
     fileName = cms.string("ElectronRecHits_ntuple.root"),
     closeFileFast = cms.untracked.bool(True)
  )


process.p = cms.Path(process.egmGsfElectronIDSequence*process.nTuplelize)

