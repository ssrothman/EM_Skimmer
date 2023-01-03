import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2500)

options = VarParsing.VarParsing('analysis')
options.parseArguments()
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            #'root://cms-xrd-global.cern.ch//store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/70000/3DB4D471-74C3-5E46-87FA-B6867DBC7AED.root',
            'file:EEDF16FC-5BD5-2D49-9DF8-DA84826DF752.root'
            )
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v27', '')

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.AOD
switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff','RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1p1_cff','RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff']

for idmod in my_id_modules:
        setupAllVIDIdsInModule(process, idmod, setupVIDPhotonSelection)

process.nTuplelize = cms.EDAnalyzer('Photon_RefinedRecHit_MiniAOD_NTuplizer',
        ebRecHits = cms.InputTag("reducedEgamma","reducedEBRecHits","RECO"),
        eeRecHits = cms.InputTag("reducedEgamma","reducedEERecHits","RECO"),
        esRecHits = cms.InputTag("reducedEgamma","reducedESRecHits","RECO"),
        debug = cms.bool(False),
        isMC = cms.bool(False),
        useOuterHits = cms.bool(False),
        rhoFastJet = cms.InputTag("fixedGridRhoFastjetAll"),
        photons = cms.InputTag("slimmedPhotons"),
        genParticles = cms.InputTag("packedGenParticles"),
        #MVA Based Id
        eleMediumIdMap = cms.InputTag(""),
        eleTightIdMap = cms.InputTag(""),
        #Calo clusters
        ebNeighbourXtalMap = cms.FileInPath("EM_ID_GNN/EM_Skimmer/data/EB_xtal_dR0p3_map.root"),
        eeNeighbourXtalMap = cms.FileInPath("EM_ID_GNN/EM_Skimmer/data/EE_xtal_dR0p3_map.root")
	)


process.TFileService = cms.Service("TFileService",
     fileName = cms.string("GammaRecHits_ntuple.root"),
     closeFileFast = cms.untracked.bool(True)
  )

process.p = cms.Path(process.egmPhotonIDSequence*process.nTuplelize)
