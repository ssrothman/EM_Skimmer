import FWCore.ParameterSet.Config as cms

nTuplize = cms.EDAnalyzer('Electron_RefinedRecHit_NTuplizer',
        isMC = cms.bool(False),
        miniAODRun = cms.bool(False),
        refinedCluster = cms.bool(True), # use refined clusters by default 
        rhoFastJet = cms.InputTag("fixedGridRhoAll"),
        electrons = cms.InputTag("gedGsfElectrons"),
        genParticles = cms.InputTag(""),
        #MVA Based Id
        eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID_Fall17_iso_V2_wp90"),
        eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID_Fall17_iso_V2_wp80")
	)

#AOD
nTuplelize = cms.EDAnalyzer('Electron_RefinedRecHit_NTuplizer',
        rhoFastJet = cms.InputTag("fixedGridRhoAll"),
        photons = cms.InputTag("photons"),
        electrons = cms.InputTag("gedGsfElectrons"),
        genParticles = cms.InputTag("genParticles"),
        refinedCluster = cms.bool(False),
        isMC = cms.bool(False), 
        #MVA Based Id
	eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose"),
        eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium"),
        eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight")
        )

#AODSIM
nTuplelize = cms.EDAnalyzer('Electron_RefinedRecHit_NTuplizer',
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

#MINIAOD
nTuplelize = cms.EDAnalyzer('Electron_RefinedRecHit_MiniAOD_NTuplizer',
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

#MINIAODSIM
nTuplelize = cms.EDAnalyzer('Electron_RefinedRecHit_MiniAOD_NTuplizer',
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
