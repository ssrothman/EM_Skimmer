// -*- C++ -*-
//
// Package:    Electron_GNN_Regression/Photon_RefinedRecHit_NTuplizer
// Class:      Photon_RefinedRecHit_NTuplizer
//
/**\class Photon_RefinedRecHit_NTuplizer Photon_RefinedRecHit_NTuplizer.cc Electron_GNN_Regression/Photon_RefinedRecHit_NTuplizer/plugins/Photon_RefinedRecHit_NTuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rajdeep Mohan Chatterjee
//         Created:  Fri, 21 Feb 2020 11:38:58 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include "TTree.h"
#include "Math/VectorUtil.h"
#include "TFile.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
//
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class Photon_RefinedRecHit_NTuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Photon_RefinedRecHit_NTuplizer(const edm::ParameterSet&);
      ~Photon_RefinedRecHit_NTuplizer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

     std::vector<float> Hit_ES_Eta[2];
     std::vector<float> Hit_ES_Phi[2];
     std::vector<float> Hit_ES_X[2];
     std::vector<float> Hit_ES_Y[2];
     std::vector<float> Hit_ES_Z[2];
     std::vector<float> ES_RecHitEn[2];


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


//   cluster tools
      EcalClusterLazyTools *clustertools;
      noZS::EcalClusterLazyTools *clustertools_NoZS;

//   Identify if the SC lies in EB OR EE based on its seed
     bool isEB = 0;
     bool isEE = 0; // !isEB not sufficient since later will try to include the preshower as well


//     bool GetGenMatchType(const reco::Photon& Photon, const reco::GenParticle& GenColl, int pdgId, double dRThresh);
// Get the hits from the ES
//     std::vector<GlobalPoint> GetESPlaneRecHits(const reco::SuperCluster& sc, unsigned int planeIndex) const;
     void GetESPlaneRecHits(const reco::SuperCluster& sc, const CaloGeometry* &geo, unsigned int phonum, unsigned int planeIndex);

//   clear the vectors 
     void ClearTreeVectors();
// ----------member data ---------------------------
     TTree* T;
     // Variables for Run info.
     int run;
     int event;
     int lumi;
     // Electron variables
     int nPhotons_;
     Float_t rho;
     std::vector<float> iEta[2];
     std::vector<float> iPhi[2];
     std::vector<float> Hit_Eta[2];
     std::vector<float> Hit_Phi[2];
     std::vector<float> Hit_X[2];
     std::vector<float> Hit_Y[2];
     std::vector<float> Hit_Z[2];


     std::vector<float> RecHitFrac[2];
     std::vector<float> RecHitEn[2];

     std::vector<float> Pho_pt_;
     std::vector<float> Pho_eta_;
     std::vector<float> Pho_phi_;
     std::vector<float> Pho_energy_;
     std::vector<float> Pho_ecal_mustache_energy_;

     std::vector<float> Pho_R9;
     std::vector<float> Pho_S4;
     std::vector<float> Pho_SigIEIE;
     std::vector<float> Pho_SigIPhiIPhi;
     std::vector<float> Pho_SCEtaW;
     std::vector<float> Pho_SCPhiW;
     std::vector<float> Pho_CovIEtaIEta;
     std::vector<float> Pho_CovIEtaIPhi;
     std::vector<float> Pho_ESSigRR;
     std::vector<float> Pho_SCRawE;
     std::vector<float> Pho_SC_ESEnByRawE;
     std::vector<float> Pho_HadOverEm;

     std::vector<float> Pho_Gen_Pt;
     std::vector<float> Pho_Gen_Eta;
     std::vector<float> Pho_Gen_Phi;
     std::vector<float> Pho_Gen_E;


     std::vector<int> passLooseId_;
     std::vector<int> passMediumId_;
     std::vector<int> passTightId_;
     std::vector<int> passMVAMediumId_;

     std::vector<int> isTrue_;

      // -----------------Handles--------------------------
      edm::Handle<double> rhoHandle;
      edm::Handle<EcalRecHitCollection> EBRechitsHandle;
      edm::Handle<EcalRecHitCollection> EERechitsHandle;
      edm::Handle<EcalRecHitCollection> ESRechitsHandle;
      edm::Handle<edm::View<reco::Photon> > photons;
      edm::Handle<edm::View<reco::GenParticle> > genParticles;
      edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
      edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
      //---------------- Input Tags-----------------------
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<EcalRecHitCollection> recHitCollectionEBToken_;
      edm::EDGetTokenT<EcalRecHitCollection> recHitCollectionEEToken_;
      edm::EDGetTokenT<EcalRecHitCollection> recHitCollectionESToken_;
      edm::EDGetToken photonsToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Photon_RefinedRecHit_NTuplizer::Photon_RefinedRecHit_NTuplizer(const edm::ParameterSet& iConfig):
   rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastJet"))),
   recHitCollectionEBToken_(consumes<EcalRecHitCollection>(edm::InputTag("reducedEcalRecHitsEB"))),
   recHitCollectionEEToken_(consumes<EcalRecHitCollection>(edm::InputTag("reducedEcalRecHitsEE"))),
   recHitCollectionESToken_(consumes<EcalRecHitCollection>(edm::InputTag("reducedEcalRecHitsES"))),
   eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
   eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap")))
{
   //now do what ever initialization is needed
   photonsToken_ = mayConsume<edm::View<reco::Photon> >(iConfig.getParameter<edm::InputTag>("photons"));
   genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"));
   usesResource("TFileService");
}


Photon_RefinedRecHit_NTuplizer::~Photon_RefinedRecHit_NTuplizer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)




}


//
// member functions
//

// ------------ method called for each event  ------------
void
Photon_RefinedRecHit_NTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;

   iEvent.getByToken(rhoToken_, rhoHandle);
   iEvent.getByToken(recHitCollectionEBToken_, EBRechitsHandle);
   iEvent.getByToken(recHitCollectionEEToken_, EERechitsHandle);
   iEvent.getByToken(recHitCollectionESToken_, ESRechitsHandle);
   iEvent.getByToken(photonsToken_, photons);
   iEvent.getByToken(genParticlesToken_, genParticles);
   iEvent.getByToken(eleMediumIdMapToken_, medium_id_decisions);
   iEvent.getByToken(eleTightIdMapToken_ , tight_id_decisions);


   ESHandle<CaloGeometry> pG;
   iSetup.get<CaloGeometryRecord>().get(pG);
   const CaloGeometry* geo = pG.product();
   const CaloSubdetectorGeometry* ecalEBGeom = static_cast<const CaloSubdetectorGeometry*>(geo->getSubdetectorGeometry(DetId::Ecal, EcalBarrel));
   const CaloSubdetectorGeometry* ecalEEGeom = static_cast<const CaloSubdetectorGeometry*>(geo->getSubdetectorGeometry(DetId::Ecal, EcalEndcap));

/////////////////////////////PU related variables////////////////////////////////////////////////////
   rho = *rhoHandle;
///////////////////////////////////////////////////////////////////////////////////////////////////////


   clustertools_NoZS = new noZS::EcalClusterLazyTools(iEvent, iSetup, recHitCollectionEBToken_, recHitCollectionEEToken_);

//Clear all vectors to be written to the tree
   ClearTreeVectors();
   run=0;  event=0;  lumi=0;

///////////////////////////Fill Electron/Photon related stuff/////////////////////////////////////////////////////
   nPhotons_ = 0;
   const EcalRecHitCollection *recHitsEB = clustertools_NoZS->getEcalEBRecHitCollection();
   const EcalRecHitCollection *recHitsEE = clustertools_NoZS->getEcalEERecHitCollection();
   for (size_t i = 0; i < photons->size(); ++i){
	if (nPhotons_ == 2) break;
	const auto pho = photons->ptrAt(i);
        if( pho->pt() < 10 ) continue;
	const SuperClusterRef& sc = pho->superCluster();
	std::vector< std::pair<DetId, float> > hitsAndFractions = sc->hitsAndFractions();
	isEB = ((*sc->seed()).hitsAndFractions().at(0).first.subdetId() == EcalBarrel);
	isEE = ((*sc->seed()).hitsAndFractions().at(0).first.subdetId() == EcalEndcap);
	EBDetId* DidEB;
	EEDetId* DidEE;
	EcalRecHitCollection::const_iterator oneHit;
        for (const auto&  detitr : hitsAndFractions) {
		if(isEB){
			DidEB = new EBDetId(detitr.first.rawId());
			DetId Did   = detitr.first.rawId();
                        shared_ptr<const CaloCellGeometry> geom = ecalEBGeom->getGeometry(Did);
			oneHit = recHitsEB->find( (detitr.first) ) ;
			iEta[nPhotons_].push_back(DidEB->ieta());
			iPhi[nPhotons_].push_back(DidEB->iphi());
			Hit_Eta[nPhotons_].push_back(geom->etaPos());
			Hit_Phi[nPhotons_].push_back(geom->phiPos());
			Hit_X[nPhotons_].push_back(geom->getPosition().x());
			Hit_Y[nPhotons_].push_back(geom->getPosition().y());
			Hit_Z[nPhotons_].push_back(geom->getPosition().z());
		}
		else if(isEE){
			DidEE = new EEDetId(detitr.first.rawId());
			DetId Did   = detitr.first.rawId();
			shared_ptr<const CaloCellGeometry> geom = ecalEEGeom->getGeometry(Did);
                        oneHit = recHitsEE->find( (detitr.first) ) ;
                        iEta[nPhotons_].push_back(DidEE->ix());
                        iPhi[nPhotons_].push_back(DidEE->iy());
			Hit_Eta[nPhotons_].push_back(geom->etaPos());
                        Hit_Phi[nPhotons_].push_back(geom->phiPos());
                        Hit_X[nPhotons_].push_back(geom->getPosition().x());
                        Hit_Y[nPhotons_].push_back(geom->getPosition().y());
                        Hit_Z[nPhotons_].push_back(geom->getPosition().z());
		}
	
		RecHitEn[nPhotons_].push_back(oneHit->energy());
		RecHitFrac[nPhotons_].push_back(detitr.second);
	}  

	if(isEE){
		GetESPlaneRecHits(*sc, geo, nPhotons_, 1);     
 		GetESPlaneRecHits(*sc, geo, nPhotons_, 2);
	}

        nPhotons_++;
        Pho_pt_.push_back( pho->pt() );
        Pho_eta_.push_back( pho->superCluster()->eta() );
        Pho_phi_.push_back( pho->superCluster()->phi() );
        Pho_energy_.push_back( pho->energy() );
	Pho_ecal_mustache_energy_.push_back( sc->energy() );
        Pho_R9.push_back(pho->full5x5_r9());
        Pho_SigIEIE.push_back(pho->full5x5_sigmaIetaIeta());
        Pho_SigIPhiIPhi.push_back(pho->full5x5_showerShapeVariables().sigmaIphiIphi);
        Pho_SCEtaW.push_back(pho->superCluster()->etaWidth());
        Pho_SCPhiW.push_back(pho->superCluster()->phiWidth());
	Pho_HadOverEm.push_back(pho->hadronicOverEm());
        const CaloClusterPtr seed_clu = pho->superCluster()->seed();
//        if (!seed_clu) continue;
//        Pho_CovIEtaIEta.push_back(clustertools_NoZS->localCovariances(*seed_clu)[0]);
//        Pho_CovIEtaIPhi.push_back(clustertools_NoZS->localCovariances(*seed_clu)[1]);
//	Pho_ESSigRR.push_back(clustertools->eseffsirir( *(pho->superCluster()) ) );
	Pho_SCRawE.push_back(pho->superCluster()->rawEnergy());
        Pho_SC_ESEnByRawE.push_back( (pho->superCluster()->preshowerEnergy())/(pho->superCluster()->rawEnergy()) );
//        Pho_S4.push_back(clustertools_NoZS->e2x2( *seed_clu ) / clustertools_NoZS->e5x5( *seed_clu ) );

////// Look up and save the ID decisions
//        bool isPassMedium = (*medium_id_decisions)[pho];
//        bool isPassTight  = (*tight_id_decisions)[pho];
//        passMediumId_.push_back( (int) isPassMedium);
//        passTightId_.push_back ( (int) isPassTight );

   }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////// Gen Stuff hardcaded for status 1 photons for now /////////////////////////////////////

  for(edm::View<GenParticle>::const_iterator part = genParticles->begin(); part != genParticles->end(); ++part){
        if( part->status()==1  && abs(part->pdgId())==22 ){
                Pho_Gen_Pt.push_back(part->pt());
                Pho_Gen_Eta.push_back(part->eta());
                Pho_Gen_Phi.push_back(part->phi());
                Pho_Gen_E.push_back(part->energy());
        }
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////Run, event, lumi//////////////////////////////////
   run=iEvent.id().run();
   event=iEvent.id().event();
   lumi=iEvent.luminosityBlock();
///////////////////////////////////////////////////////////////////////////



   T->Fill(); // Write out the events
   delete clustertools_NoZS;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
Photon_RefinedRecHit_NTuplizer::beginJob()
{
        edm::Service<TFileService> fs;
        T=fs->make<TTree>("T","MyTuple");
	T->Branch("iEtaPho1"  ,  &(iEta[0]));
        T->Branch("iPhiPho1"  ,  &(iPhi[0]));
	T->Branch("Hit_ES_Eta_Pho1"  ,  &(Hit_ES_Eta[0]));
        T->Branch("Hit_ES_Phi_Pho1"  ,  &(Hit_ES_Phi[0]));
        T->Branch("Hit_ES_X_Pho1"  ,  &(Hit_ES_X[0]));
        T->Branch("Hit_ES_Y_Pho1"  ,  &(Hit_ES_Y[0]));
        T->Branch("Hit_ES_Z_Pho1"  ,  &(Hit_ES_Z[0]));
        T->Branch("ES_RecHitEnPho1"  ,  &(ES_RecHitEn[0]));

	T->Branch("Hit_Eta_Pho1"  ,  &(Hit_Eta[0]));
        T->Branch("Hit_Phi_Pho1"  ,  &(Hit_Phi[0]));
	T->Branch("Hit_X_Pho1"  ,  &(Hit_X[0]));
        T->Branch("Hit_Y_Pho1"  ,  &(Hit_Y[0]));
	T->Branch("Hit_Z_Pho1"  ,  &(Hit_Z[0]));
        T->Branch("RecHitEnPho1"  ,  &(RecHitEn[0]));
	T->Branch("RecHitFracPho1"  ,  &(RecHitFrac[0]));
	T->Branch("iEtaPho2"  ,  &(iEta[1]));
        T->Branch("iPhiPho2"  ,  &(iPhi[1]));
	T->Branch("Hit_ES_Eta_Pho2"  ,  &(Hit_ES_Eta[1]));
        T->Branch("Hit_ES_Phi_Pho2"  ,  &(Hit_ES_Phi[1]));
        T->Branch("Hit_ES_X_Pho2"  ,  &(Hit_ES_X[1]));
        T->Branch("Hit_ES_Y_Pho2"  ,  &(Hit_ES_Y[1]));
        T->Branch("Hit_ES_Z_Pho2"  ,  &(Hit_ES_Z[1]));
	T->Branch("ES_RecHitEnPho2"  ,  &(ES_RecHitEn[1]));

	T->Branch("Hit_Eta_Pho2"  ,  &(Hit_Eta[1]));
        T->Branch("Hit_Phi_Pho2"  ,  &(Hit_Phi[1]));
	T->Branch("Hit_X_Pho2"  ,  &(Hit_X[1]));
        T->Branch("Hit_Y_Pho2"  ,  &(Hit_Y[1]));
        T->Branch("Hit_Z_Pho2"  ,  &(Hit_Z[1]));
        T->Branch("RecHitEnPho2"  ,  &(RecHitEn[1]));
        T->Branch("RecHitFracPho2"  ,  &(RecHitFrac[1]));

        T->Branch("nPhotons",  &nPhotons_ , "nPho/I");
        T->Branch("pt"  ,  &Pho_pt_);
        T->Branch("eta" ,  &Pho_eta_ );
        T->Branch("phi" ,  &Pho_phi_ );
        T->Branch("energy", &Pho_energy_);
	T->Branch("energy_ecal_mustache", &Pho_ecal_mustache_energy_);

        T->Branch("passMediumId" ,  &passMediumId_ );
        T->Branch("passTightId"  ,  &passTightId_ );
        T->Branch("passMVAMediumId", &passMVAMediumId_);

        T->Branch("Pho_R9"  ,  &Pho_R9);
        T->Branch("Pho_S4"  ,  &Pho_S4);
        T->Branch("Pho_SigIEIE"  ,  &Pho_SigIEIE);
        T->Branch("Pho_SigIPhiIPhi" , &Pho_SigIPhiIPhi);
        T->Branch("Pho_SCEtaW"  ,  &Pho_SCEtaW);
        T->Branch("Pho_SCPhiW"  ,  &Pho_SCPhiW);
        T->Branch("Pho_CovIEtaIEta"  ,  &Pho_CovIEtaIEta);
        T->Branch("Pho_CovIEtaIPhi"  ,  &Pho_CovIEtaIPhi);
        T->Branch("Pho_ESSigRR"  ,  &Pho_ESSigRR);
        T->Branch("Pho_SCRawE"  ,  &Pho_SCRawE);
        T->Branch("Pho_SC_ESEnByRawE"  ,  &Pho_SC_ESEnByRawE);
	T->Branch("Pho_HadOverEm"  ,  &Pho_HadOverEm);

        T->Branch("Pho_Gen_Pt" , &Pho_Gen_Pt);
        T->Branch("Pho_Gen_Eta" , &Pho_Gen_Eta);
        T->Branch("Pho_Gen_Phi" , &Pho_Gen_Phi);
        T->Branch("Pho_Gen_E" , &Pho_Gen_E);

        T->Branch("rho", &rho, "rho/F");

        T->Branch("run",&run,"run/I");
        T->Branch("event",&event,"event/I");
        T->Branch("lumi",&lumi,"lumi/I");

}

// ------------ method called once each job just after ending the event loop  ------------
void
Photon_RefinedRecHit_NTuplizer::endJob()
{
}

/*
//Evaluate if the gen particle dR matched to a reco photon is also a photon
bool Photon_RefinedRecHit_NTuplizer::GetGenMatchType(const reco::Photon& Photon, const reco::GenParticle& GenColl, int pdgId, double dRThresh){
	

}
*/


// Extract the rechits of the SC from the ES layers
void Photon_RefinedRecHit_NTuplizer::GetESPlaneRecHits(const reco::SuperCluster& sc, const CaloGeometry* &geo, unsigned int phonum, unsigned int planeIndex) {
        double RawenergyPlane = 0.;
        double pfRawenergyPlane = 0.;
//      if(!_ESRechitsHandle.isValid())
//              return RawenergyPlane;
              
//        if (!sc.preshowerClusters().isAvailable()) //protection for miniAOD
//                break;

	int NumHits = 0;

	const CaloSubdetectorGeometry* ecalESGeom = static_cast<const CaloSubdetectorGeometry*>(geo->getSubdetectorGeometry(DetId::Ecal, EcalPreshower));


        for(auto iES = sc.preshowerClustersBegin(); iES != sc.preshowerClustersEnd(); ++iES) {//loop over preshower clusters
                const std::vector< std::pair<DetId, float> > hits = (*iES)->hitsAndFractions();
                for(std::vector<std::pair<DetId, float> >::const_iterator rh = hits.begin(); rh != hits.end(); ++rh) { // loop over recHits of the cluster
                        //      std::cout << "print = " << (*iES)->printHitAndFraction(iCount);
                        //      ++iCount;
                        for(ESRecHitCollection::const_iterator esItr = ESRechitsHandle->begin(); esItr != ESRechitsHandle->end(); ++esItr) {//loop over ES rechits to find the one in the cluster
                                ESDetId rhid = ESDetId(esItr->id());
                                if(rhid == (*rh).first) { // found ESrechit
//                                        std::cout << " ES energy = " << esItr->energy() << " pf energy = " << (*rh).second << std::endl;
                                        if((int) rhid.plane() == (int) planeIndex) {
						std::shared_ptr<const CaloCellGeometry> geom = ecalESGeom->getGeometry(rhid);
						Hit_ES_Eta[phonum].push_back( geom->etaPos() );
                                                Hit_ES_Phi[phonum].push_back( geom->phiPos() );
						Hit_ES_X[phonum].push_back( geom->getPosition().x() );
						Hit_ES_Y[phonum].push_back( geom->getPosition().y() );
						Hit_ES_Z[phonum].push_back( geom->getPosition().z() ) ;
						ES_RecHitEn[phonum].push_back(esItr->energy());
//						std::cout << "Preshower" <<std::setprecision(4) << " Eta = " <<geom->etaPos() << " : " <<" Phi = "<< geom->phiPos() << " 3D position" << geom->getPosition().z() << std::endl;
                                                RawenergyPlane += esItr->energy();
                                                pfRawenergyPlane += rh->second;
						NumHits++;
                                        }
                                        break;
                                }
                        }
                }

//		std::cout<<std::endl<<" Number of ES hits in plane 1 = "<<NumHits<<std::endl;
        }

 //       return RawenergyPlane;
}


//Clear tree vectors each time analyze method is called
void Photon_RefinedRecHit_NTuplizer::ClearTreeVectors()
{
	nPhotons_ = 0;
	iEta[0].clear();
	iPhi[0].clear();

	
	Hit_ES_Eta[0].clear();
        Hit_ES_Phi[0].clear();
        Hit_ES_X[0].clear();
        Hit_ES_Y[0].clear();
        Hit_ES_Z[0].clear();
	ES_RecHitEn[0].clear();


        Hit_Eta[0].clear();
	Hit_Phi[0].clear();
	Hit_X[0].clear();
	Hit_Y[0].clear();
	Hit_Z[0].clear();
	RecHitEn[0].clear();
	RecHitFrac[0].clear();
	iEta[1].clear();
        iPhi[1].clear();

	Hit_ES_Eta[1].clear();
        Hit_ES_Phi[1].clear();
        Hit_ES_X[1].clear();
        Hit_ES_Y[1].clear();
        Hit_ES_Z[1].clear();
	ES_RecHitEn[1].clear();

	Hit_Eta[1].clear();
        Hit_Phi[1].clear();
        Hit_X[1].clear();
        Hit_Y[1].clear();
        Hit_Z[1].clear();
        RecHitEn[1].clear();
        RecHitFrac[1].clear();
	Pho_pt_.clear();
	Pho_eta_.clear();
	Pho_phi_.clear();
	Pho_energy_.clear();
	Pho_ecal_mustache_energy_.clear();
	Pho_R9.clear();
	Pho_S4.clear();
	Pho_SigIEIE.clear();
	Pho_SigIPhiIPhi.clear();
	Pho_SCEtaW.clear();
	Pho_SCPhiW.clear();
	Pho_CovIEtaIEta.clear();
	Pho_CovIEtaIPhi.clear();
	Pho_ESSigRR.clear();
	Pho_SCRawE.clear();
	Pho_SC_ESEnByRawE.clear();
	Pho_HadOverEm.clear();

	Pho_Gen_Pt.clear();
	Pho_Gen_Eta.clear();
	Pho_Gen_Phi.clear();
	Pho_Gen_E.clear();

	passMediumId_.clear();
	passTightId_ .clear();
	passMVAMediumId_.clear();

	isTrue_.clear();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Photon_RefinedRecHit_NTuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Photon_RefinedRecHit_NTuplizer);
