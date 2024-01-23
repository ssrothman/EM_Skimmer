import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

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
                                #'/store/mc/RunIISummer20UL18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/260000/8D4B008D-BE14-DD47-8949-C4640525DA1F.root'
                                #'file:EEDF16FC-5BD5-2D49-9DF8-DA84826DF752.root'
                                '/store/mc/RunIISummer20UL18MiniAODv2/QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2560000/04CAD391-974F-104E-97C5-8FCB48DD8874.root'
                                )
                            )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoEgamma/PhotonIdentification/photonIDValueMapProducer_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v27', '')

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = [
#'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff',
#'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1p1_cff',
'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff',
'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff']

for idmod in my_id_modules:
        setupAllVIDIdsInModule(process, idmod, setupVIDPhotonSelection)

process.nTuplelize = cms.EDAnalyzer('Photon_RefinedRecHit_MiniAOD_NTuplizer',
        ebRecHits = cms.InputTag("reducedEgamma","reducedEBRecHits","PAT"),
        eeRecHits = cms.InputTag("reducedEgamma","reducedEERecHits","PAT"),
        esRecHits = cms.InputTag("reducedEgamma","reducedESRecHits","PAT"),
        debug = cms.bool(False),
        isMC = cms.bool(True),
        useOuterHits = cms.bool(False),
        rhoFastJet = cms.InputTag("fixedGridRhoFastjetAll"),
        photons = cms.InputTag("slimmedPhotons"),
        genParticles = cms.InputTag("packedGenParticles"),
        #MVA Based Id
        eleMediumIdMap = cms.InputTag(""),
        eleTightIdMap = cms.InputTag(""),
        photonLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Fall17-94X-V2-loose"),
        photonMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Fall17-94X-V2-medium"),
        photonTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Fall17-94X-V2-tight"),
        photonMVAwp80IdMap = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIFall17-v2-wp80"),
        photonMVAwp90IdMap = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIFall17-v2-wp90"),
        #Calo clusters
        ebNeighbourXtalMap = cms.FileInPath("EM_ID_GNN/EM_Skimmer/data/EB_xtal_dR0p3_map.root"),
        eeNeighbourXtalMap = cms.FileInPath("EM_ID_GNN/EM_Skimmer/data/EE_xtal_dR0p3_map.root")
	)


process.TFileService = cms.Service("TFileService",
     fileName = cms.string("GammaRecHits_ntuple.root"),
#     fileName = cms.string("Tree_Gamma_ABCD.root"),
      closeFileFast = cms.untracked.bool(True)
  )


process.p = cms.Path(process.egmPhotonIDSequence*process.nTuplelize)

