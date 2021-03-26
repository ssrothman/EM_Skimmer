import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

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
#            'root://cms-xrd-global.cern.ch///store/mc/RunIIWinter19PFCalibDR/DoublePhoton_FlatPt-5To300/AODSIM/2018ConditionsFlatPU0to70ECALGT_105X_upgrade2018_realistic_IdealEcalIC_v4-v1/60000/AD91349A-854D-DD47-B82A-B77A907D1A9C.root'
	    'root://cms-xrd-global.cern.ch///store/mc/RunIISummer19UL18RECO/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/00000/CE71D346-5BDC-9348-8887-872D3C69264D.root'
#                options.inputFile
                )
                            )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.AOD
switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff']

for idmod in my_id_modules:
        setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)

process.nTuplelize = cms.EDAnalyzer('Photon_RefinedRecHit_NTuplizer',
        rhoFastJet = cms.InputTag("fixedGridRhoFastjetAll"),
        photons = cms.InputTag("gedPhotons"),
        genParticles = cms.InputTag("genParticles"),
        #Cut Based Id
        eleMediumIdMap = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIFall17-v1-wp90"),
        eleTightIdMap = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIFall17-v1-wp90")
	)


process.TFileService = cms.Service("TFileService",
     fileName = cms.string("nTupleMC.root"),
#     fileName = cms.string("Tree_Gamma_ABCD.root"),
      closeFileFast = cms.untracked.bool(True)
  )


process.p = cms.Path(process.egmGsfElectronIDSequence * process.nTuplelize)

