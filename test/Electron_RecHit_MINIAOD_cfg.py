import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
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
            '/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/270002/0C306E94-C57B-4347-A6A0-58C4340C6E2F.root'
#                options.inputFile
                ),
                inputCommands=cms.untracked.vstring(
                    'keep *',
                    'drop recoTrackExtrasedmAssociation_muonReducedTrackExtras_*_*'
                )
                            )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v27', '')

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)


# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff']


for idmod in my_id_modules:
        setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)

process.nTuplelize = cms.EDAnalyzer('Electron_RefinedRecHit_MiniAOD_NTuplizer',
        ebRecHits = cms.InputTag("reducedEgamma","reducedEBRecHits","RECO"),
        eeRecHits = cms.InputTag("reducedEgamma","reducedEERecHits","RECO"),
        esRecHits = cms.InputTag("reducedEgamma","reducedESRecHits","RECO"),
        rhoFastJet = cms.InputTag("fixedGridRhoAll"),
        electrons = cms.InputTag("slimmedElectrons"),
        photons = cms.InputTag("slimmedPhotons"),
        genParticles = cms.InputTag("packedGenParticles"),
        refinedCluster = cms.bool(False),
        isMC = cms.bool(False), 
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

