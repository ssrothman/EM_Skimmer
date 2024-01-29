#include "../interface/Photon_RefinedRecHit_MiniAOD_NTuplizer.h"
#include "DataFormats/Math/interface/deltaR.h"

//using reco::TrackCollection;

using namespace std;
using namespace edm;
using namespace reco;
using namespace pat;

Photon_RefinedRecHit_MiniAOD_NTuplizer::Photon_RefinedRecHit_MiniAOD_NTuplizer(const edm::ParameterSet& iConfig):
   recHitCollectionEBToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("ebRecHits"))),
   recHitCollectionEEToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("eeRecHits"))),
   recHitCollectionESToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("esRecHits"))),
   rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastJet"))),
   eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
   eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
   photonLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonLooseIdMap"))),
   photonMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonMediumIdMap"))),
   photonTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonTightIdMap"))),
   photonMVAwp80IdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonMVAwp80IdMap"))),
   photonMVAwp90IdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonMVAwp90IdMap")))
{
   //now do what ever initialization is needed
   DEBUG = iConfig.getParameter<bool>("debug");
   isMC_ = iConfig.getParameter<bool>("isMC");
   useOuterHits_ = iConfig.getParameter<bool>("useOuterHits");
   ebNeighbourXtalMap_ = iConfig.getParameter<edm::FileInPath>("ebNeighbourXtalMap");
   eeNeighbourXtalMap_ = iConfig.getParameter<edm::FileInPath>("eeNeighbourXtalMap");
   photonsToken_ = mayConsume<std::vector<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photons"));
   genParticlesToken_ = mayConsume<std::vector<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"));
   usesResource("TFileService");
   //if (useOuterHits_) readLUTables(); // read look up tables (dR = 0.3 for now)
   if (useOuterHits_) readROOTLUTables(); // read look up tables (dR = 0.3 for now)

}


Photon_RefinedRecHit_MiniAOD_NTuplizer::~Photon_RefinedRecHit_MiniAOD_NTuplizer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// Method to load look up tables
void Photon_RefinedRecHit_MiniAOD_NTuplizer::readROOTLUTables(){
   int id;
   vector<int>* nids;

   nids = 0;

   TFile* file_eb = new TFile(ebNeighbourXtalMap_.fullPath().c_str(), "read");
   TFile* file_ee = new TFile(eeNeighbourXtalMap_.fullPath().c_str(), "read");

   TTree* tree_eb = (TTree*)file_eb->Get("xtal_map");
   TTree* tree_ee = (TTree*)file_ee->Get("xtal_map");

   tree_eb->SetBranchAddress("id", &id);
   tree_eb->SetBranchAddress("nids", &nids);

   for (int ientry=0; ientry<tree_eb->GetEntriesFast(); ientry++){
      tree_eb->GetEntry(ientry);
      ebnxtals.insert({id, nids[0]});
   }
   
   tree_ee->SetBranchAddress("id", &id);
   tree_ee->SetBranchAddress("nids", &nids);

   for (int ientry=0; ientry<tree_ee->GetEntriesFast(); ientry++){
      tree_ee->GetEntry(ientry);
      eenxtals.insert({id, nids[0]});
   }

   file_eb->Close();
   file_ee->Close();

}

// ------------ method called for each event  ------------
   void
Photon_RefinedRecHit_MiniAOD_NTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   if (DEBUG) cout << "getting handles"<< endl;
   iEvent.getByToken(rhoToken_, rhoHandle);
   iEvent.getByToken(recHitCollectionEBToken_, EBRechitsHandle);
   iEvent.getByToken(recHitCollectionEEToken_, EERechitsHandle);
   iEvent.getByToken(recHitCollectionESToken_, ESRechitsHandle);
   iEvent.getByToken(photonsToken_, photons);
   iEvent.getByToken(genParticlesToken_, genParticles);
   iEvent.getByToken(eleMediumIdMapToken_, medium_id_decisions);
   iEvent.getByToken(eleTightIdMapToken_ , tight_id_decisions);
   iEvent.getByToken(photonLooseIdMapToken_, photon_loose_id_decisions);
   iEvent.getByToken(photonMediumIdMapToken_, photon_medium_id_decisions);
   iEvent.getByToken(photonTightIdMapToken_, photon_tight_id_decisions);
   iEvent.getByToken(photonMVAwp80IdMapToken_, photon_mvawp80_id_decisions);
   iEvent.getByToken(photonMVAwp90IdMapToken_, photon_mvawp90_id_decisions);

   if (DEBUG) cout << "getting geometry"<< endl;
   ESHandle<CaloGeometry> pG;
   iSetup.get<CaloGeometryRecord>().get(pG);
   iSetup.get<EcalPedestalsRcd>().get(_ped);
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

   if (DEBUG) cout << "getting cluters"<< endl;
   ///////////////////////////Fill Electron/Photon related stuff/////////////////////////////////////////////////////
   nPhotons_ = 0;
   const EcalRecHitCollection *recHitsEB = clustertools_NoZS->getEcalEBRecHitCollection();
   const EcalRecHitCollection *recHitsEE = clustertools_NoZS->getEcalEERecHitCollection();

   //cout<<"================ Event ==================="<<endl;

   if (DEBUG) cout << "getting photons"<< endl;
   for (size_t i = 0; i < photons->size(); ++i){
      //cout<<"----------- Photon ------------"<<endl;
      if (nPhotons_ == 2) break;
      
      const auto pho = photons->at(i);
      const edm::Ptr<pat::Photon> phoPtr(photons, i);

      if( pho.pt() < 10 ) continue;
   
      if (DEBUG) cout << "getting photon superclusters"<< endl;
      const SuperClusterRef& sc = pho.superCluster(); 
      //const SuperClusterRef& sc = pho->parentSuperCluster(); // mustache cluster
      std::vector< std::pair<DetId, float> > hitsAndFractions = sc->hitsAndFractions();
      isEB = ((*sc->seed()).hitsAndFractions().at(0).first.subdetId() == EcalBarrel);
      isEE = ((*sc->seed()).hitsAndFractions().at(0).first.subdetId() == EcalEndcap);

      EBDetId* DidEB;
      EEDetId* DidEE;

      EcalRecHitCollection::const_iterator oneHit;

      std::vector<int> EBHitsAndFractions;
      std::vector<int> EEHitsAndFractions;

      //cout<<"*** supercluster hits ***"<<endl;
      
      if (DEBUG) cout<<"Hits and Fractions"<<endl;

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

            // Catalog the hits EB
            EBHitsAndFractions.push_back(DidEB->numberByEtaPhi());
            if (DEBUG) cout<<*DidEB<<endl;

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

            // Catalog the hits in EE
            EEHitsAndFractions.push_back(DidEE->hashedIndex());
            if (DEBUG) cout<<*DidEE<<endl;
         }

         RecHitEn[nPhotons_].push_back(oneHit->energy());
         RecHitFrac[nPhotons_].push_back(detitr.second);
         if(oneHit->checkFlag(EcalRecHit::kGood))	RecHitQuality[nPhotons_].push_back(1);
         else RecHitQuality[nPhotons_].push_back(0);

         for (int iflag=0; iflag<EcalRecHit::kHasSwitchToGain1+1; iflag++){
            RecHitFlag_container[iflag][nPhotons_].push_back(oneHit->checkFlag(iflag));
         }

         if (DEBUG) cout<<endl<<" Reco Flags = "<<oneHit->recoFlag()<<endl;

         if(oneHit->checkFlag(EcalRecHit::kHasSwitchToGain6)) 		RecHitGain[nPhotons_].push_back(6);
         else if(oneHit->checkFlag(EcalRecHit::kHasSwitchToGain1)) RecHitGain[nPhotons_].push_back(1);
         else RecHitGain[nPhotons_].push_back(12);
         HitNoise[nPhotons_].push_back(_ped->find(detitr.first)->rms(1));
      }
      
      CaloClusterPtr seed_clu;

      // Get the cluster seed
      if (DEBUG) cout<<"cluster seeds"<<endl;
      try {
         CaloClusterPtr seed_clu = sc->seed();
      }
      catch (...){
         if (DEBUG) cout<<"Invalid cluster seed pointer!"<<endl;
      }

      
      math::XYZPoint seed_pos;
      DetId seedRawId;

      if (!seed_clu.isNull()){
         seed_pos = seed_clu->position();
         seedRawId = (seed_clu->seed()).rawId();
      }
      else {
         seed_pos = math::XYZPoint(-999.,-999.,-999.);
         seedRawId = DetId();
      }

      if (DEBUG) cout<<seed_pos<<endl;

      float seed_eta=-99.;
      float seed_phi=-99.;

      // Add the rechits within (dR==0.3 for now)
      // EB region outer hits
      if (useOuterHits_ && isEB){

         DidEB = new EBDetId(seedRawId);
         shared_ptr<const CaloCellGeometry> geom = ecalEBGeom->getGeometry(seedRawId);
         seed_eta = geom->etaPos();
         seed_phi = geom->phiPos();

         // Get list of all neighboutring xtals to the seed
         auto uncheckedXtals = ebnxtals[DidEB->numberByEtaPhi()];

         for (unsigned int ixtal=0; ixtal<uncheckedXtals.size(); ixtal++){

            if (std::find(EBHitsAndFractions.begin(),
                     EBHitsAndFractions.end(),
                     uncheckedXtals[ixtal])!=EBHitsAndFractions.end()) continue;

            EBDetId detitr = EBDetId::unhashIndex(uncheckedXtals[ixtal]); // covert from hash to EBDetId
            DetId Did   = detitr.rawId();
            DidEB = new EBDetId(Did);

            shared_ptr<const CaloCellGeometry> geom = ecalEBGeom->getGeometry(Did);
            oneHit = recHitsEB->find(detitr) ;

            if (oneHit==recHitsEB->end()) continue;

            //cout<<"detitr = "<<detitr<<endl;
            //float tmpdr = reco::deltaR(seed_ieta, seed_iphi, geom->etaPos(), geom->phiPos());
            //if (tmpdr>0.3) cout<<"--------------!!!-------- dR = "<<tmpdr<<endl;

            iEta[nPhotons_].push_back(DidEB->ieta());
            iPhi[nPhotons_].push_back(DidEB->iphi());
            Hit_Eta[nPhotons_].push_back(geom->etaPos());
            Hit_Phi[nPhotons_].push_back(geom->phiPos());
            Hit_X[nPhotons_].push_back(geom->getPosition().x());
            Hit_Y[nPhotons_].push_back(geom->getPosition().y());
            Hit_Z[nPhotons_].push_back(geom->getPosition().z());

            RecHitEn[nPhotons_].push_back(oneHit->energy());
            RecHitFrac[nPhotons_].push_back(-1); // fraction does not apply to hits outside the cluster

            if(oneHit->checkFlag(EcalRecHit::kGood))	RecHitQuality[nPhotons_].push_back(1);
            else RecHitQuality[nPhotons_].push_back(0);

            for (int iflag=0; iflag<EcalRecHit::kHasSwitchToGain1+1; iflag++){
               RecHitFlag_container[iflag][nPhotons_].push_back(oneHit->checkFlag(iflag));
            }

            if (DEBUG) cout<<endl<<" (Outer Hits EB) Reco Flags = "<<oneHit->recoFlag()<<endl;

            if(oneHit->checkFlag(EcalRecHit::kHasSwitchToGain6)) 		RecHitGain[nPhotons_].push_back(6);
            else if(oneHit->checkFlag(EcalRecHit::kHasSwitchToGain1)) RecHitGain[nPhotons_].push_back(1);
            else RecHitGain[nPhotons_].push_back(12);
            HitNoise[nPhotons_].push_back(_ped->find(detitr)->rms(1));
         }
      }    
      else {
         if (DEBUG) cout << "didn't use eb external hits" <<endl;
      }

      // EE region outer hits
      if (useOuterHits_ && isEE){

         DidEE = new EEDetId(seedRawId);
         int eextalid = DidEE->hashedIndex(); 
         shared_ptr<const CaloCellGeometry> geom = ecalEEGeom->getGeometry(seedRawId);
         seed_eta = geom->etaPos();
         seed_phi = geom->phiPos();

         // Get list of all neighbouring xtals to the seed
         auto uncheckedXtals = eenxtals[eextalid];

         for (unsigned int ixtal=0; ixtal<uncheckedXtals.size(); ixtal++){

            if (std::find(EEHitsAndFractions.begin(),
                     EEHitsAndFractions.end(),
                     uncheckedXtals[ixtal])!=EEHitsAndFractions.end()) continue;

            EEDetId detitr = EEDetId::unhashIndex(uncheckedXtals[ixtal]); // covert from hash to EEDetId
            DetId Did   = detitr.rawId();
            DidEE = new EEDetId(Did);

            shared_ptr<const CaloCellGeometry> geom = ecalEEGeom->getGeometry(Did);
            oneHit = recHitsEE->find(detitr) ;

            if (oneHit==recHitsEE->end()) continue;
            if (DEBUG) cout<<"detitr = "<<detitr<<", xtal = "<<uncheckedXtals[ixtal]<<endl;
            
            //float tmpdr = reco::deltaR(seed_eta, seed_phi, geom->etaPos(), geom->phiPos());
            //if (tmpdr>0.3) cout<<"--------------!!!-------- dR = "<<tmpdr<<endl;

            iEta[nPhotons_].push_back(DidEE->ix());
            iPhi[nPhotons_].push_back(DidEE->iy());
            Hit_Eta[nPhotons_].push_back(geom->etaPos());
            Hit_Phi[nPhotons_].push_back(geom->phiPos());
            Hit_X[nPhotons_].push_back(geom->getPosition().x());
            Hit_Y[nPhotons_].push_back(geom->getPosition().y());
            Hit_Z[nPhotons_].push_back(geom->getPosition().z());

            RecHitEn[nPhotons_].push_back(oneHit->energy());
            RecHitFrac[nPhotons_].push_back(-1);

            if(oneHit->checkFlag(EcalRecHit::kGood))	RecHitQuality[nPhotons_].push_back(1);
            else RecHitQuality[nPhotons_].push_back(0);

            for (int iflag=0; iflag<EcalRecHit::kHasSwitchToGain1+1; iflag++){
               RecHitFlag_container[iflag][nPhotons_].push_back(oneHit->checkFlag(iflag));
            }

            if (DEBUG) cout<<endl<<" (Outer Hits EE) Reco Flags = "<<oneHit->recoFlag()<<endl;

            if(oneHit->checkFlag(EcalRecHit::kHasSwitchToGain6)) 		RecHitGain[nPhotons_].push_back(6);
            else if(oneHit->checkFlag(EcalRecHit::kHasSwitchToGain1)) RecHitGain[nPhotons_].push_back(1);
            else RecHitGain[nPhotons_].push_back(12);
            HitNoise[nPhotons_].push_back(_ped->find(detitr)->rms(1));
         }
      }
      
      else {
         if (DEBUG) cout << "didn't use ee external hits" <<endl;
      }

      // Add preshower hits
      if (DEBUG) cout<<"Getting ES plane hits"<<endl;
      if(isEE){
         GetESPlaneRecHits(*sc, geo, nPhotons_, 1);     
         GetESPlaneRecHits(*sc, geo, nPhotons_, 2);
      }

      if (DEBUG) cout << "Filling Photon Variables" << endl;
      if (DEBUG) cout<< sc->eta()<<endl;

      nPhotons_++;
      Pho_pt_.push_back( pho.pt() );
      Pho_eta_.push_back( sc->eta() );
      Pho_phi_.push_back( sc->phi() );
      Pho_energy_.push_back( pho.energy() );
      Pho_ecal_mustache_energy_.push_back( sc->energy() );
      Pho_R9.push_back(pho.full5x5_r9());
      Pho_SigIEIE.push_back(pho.full5x5_sigmaIetaIeta());
      Pho_SigIPhiIPhi.push_back(pho.full5x5_showerShapeVariables().sigmaIphiIphi);
      Pho_SCEtaW.push_back(sc->etaWidth());
      Pho_SCPhiW.push_back(sc->phiWidth());
      Pho_HadOverEm.push_back(pho.hadronicOverEm());

      // if (!seed_clu) continue;
      // Pho_CovIEtaIEta.push_back(clustertools_NoZS->localCovariances(*seed_clu)[0]);
      // Pho_CovIEtaIPhi.push_back(clustertools_NoZS->localCovariances(*seed_clu)[1]);
      //	Pho_ESSigRR.push_back(clustertools->eseffsirir( *(sc) ) );
      Pho_SCRawE.push_back(sc->rawEnergy());
      Pho_SC_ESEnByRawE.push_back( (sc->preshowerEnergy())/(sc->rawEnergy()) );
      // Pho_S4.push_back(clustertools_NoZS->e2x2( *seed_clu ) / clustertools_NoZS->e5x5( *seed_clu ) );

      // Fill Isolation variables
      Pho_PFChIso.push_back(pho.photonIso());
      Pho_PFChPVIso.push_back(pho.chargedHadronPFPVIso());
      Pho_PFPhoIso.push_back(pho.chargedHadronIso());
      Pho_PFNeuIso.push_back(pho.neutralHadronIso());
      Pho_PFChWorstVetoIso.push_back(pho.chargedHadronWorstVtxGeomVetoIso());
      Pho_PFChWorstIso.push_back(pho.chargedHadronWorstVtxIso());
      Pho_EcalPFClusterIso.push_back(pho.ecalPFClusterIso());
      Pho_HcalPFClusterIso.push_back(pho.hcalPFClusterIso());

      Pho_cluster_seed_x.push_back(seed_pos.x());
      Pho_cluster_seed_y.push_back(seed_pos.y());
      Pho_cluster_seed_z.push_back(seed_pos.z());

      Pho_cluster_seed_eta.push_back(seed_eta);
      Pho_cluster_seed_phi.push_back(seed_phi);

      Pho_CorrectedEnergy.push_back(pho.getCorrectedEnergy(pho.getCandidateP4type()));
      Pho_CorrectedEnergyError.push_back(pho.getCorrectedEnergyError(pho.getCandidateP4type())); // Error in corrected energy

      bool gammaLoosePass  = (*photon_loose_id_decisions)[phoPtr];
      bool gammaMediumPass = (*photon_medium_id_decisions)[phoPtr];
      bool gammaTightPass = (*photon_tight_id_decisions)[phoPtr];
      bool gammaMVAwp80Pass = (*photon_mvawp80_id_decisions)[phoPtr];
      bool gammaMVAwp90Pass = (*photon_mvawp90_id_decisions)[phoPtr];

      Pho_loose_id.push_back(gammaLoosePass);
      Pho_medium_id.push_back(gammaMediumPass);
      Pho_tight_id.push_back(gammaTightPass);
      Pho_mvawp80_id.push_back(gammaMVAwp80Pass);
      Pho_mvawp90_id.push_back(gammaMVAwp90Pass);

   }
   
   if (DEBUG) cout << "Done Filling Photon Variables" << endl;

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////////////// Gen Stuff hardcaded for status 1 photons for now /////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
   if (isMC_){
      for(std::vector<pat::PackedGenParticle>::const_iterator part = genParticles->begin(); part != genParticles->end(); ++part){
         if( part->status()==1  && abs(part->pdgId())==22 ){
            Pho_Gen_Pt.push_back(part->pt());
            Pho_Gen_Eta.push_back(part->eta());
            Pho_Gen_Phi.push_back(part->phi());
            Pho_Gen_E.push_back(part->energy());
         }
      }

      //Do gen-matching 
      Pho_GenIdx.resize(Pho_pt_.size(), -1);
      Pho_DR.resize(Pho_pt_.size(), -1);
      
      std::vector<bool> taken(Pho_Gen_Pt.size(), false);
      for(unsigned iReco=0; iReco<Pho_pt_.size(); ++iReco){
          std::vector<double> distances(Pho_Gen_Pt.size(), 99999);
          for(unsigned iGen=0; iGen<Pho_Gen_Pt.size(); ++iGen){
            if(taken[iGen]){
              continue;
            }
            double dist = reco::deltaR(Pho_Gen_Eta[iGen], Pho_Gen_Phi[iGen], Pho_eta_[iReco], Pho_phi_[iReco]);
            distances[iGen] = dist;
          }
          int minIndex = std::min_element(distances.begin(), distances.end()) - distances.begin();
          printf("The minIndex is %d\n", minIndex);
          Pho_GenIdx[iReco] = minIndex;
          Pho_DR[iReco] = distances[minIndex];
          taken[minIndex] = true;
      }
   }

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////// Run, event, lumi //////////////////////////////////

   run=iEvent.id().run();
   event=iEvent.id().event();
   lumi=iEvent.luminosityBlock();

   T->Fill(); // Write out the events
   delete clustertools_NoZS;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
   void
Photon_RefinedRecHit_MiniAOD_NTuplizer::beginJob()
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
   T->Branch("RecHitGain1"  ,  &(RecHitGain[0]));
   T->Branch("RecHitQuality1", &(RecHitQuality[0]));
   T->Branch("HitNoisePho1", &(HitNoise[0]));

   T->Branch("RecHitFlag_kGood_pho1", &(RecHitFlag_kGood[0]));
   T->Branch("RecHitFlag_kPoorReco_pho1", &(RecHitFlag_kPoorReco[0]));
   T->Branch("RecHitFlag_kOutOfTime_pho1", &(RecHitFlag_kOutOfTime[0]));
   T->Branch("RecHitFlag_kFaultyHardware_pho1", &(RecHitFlag_kFaultyHardware[0]));
   T->Branch("RecHitFlag_kNoisy_pho1", &(RecHitFlag_kNoisy[0]));
   T->Branch("RecHitFlag_kPoorCalib_pho1", &(RecHitFlag_kPoorCalib[0]));
   T->Branch("RecHitFlag_kSaturated_pho1", &(RecHitFlag_kSaturated[0]));
   T->Branch("RecHitFlag_kLeadingEdgeRecovered_pho1", &(RecHitFlag_kLeadingEdgeRecovered[0]));
   T->Branch("RecHitFlag_kNeighboursRecovered_pho1", &(RecHitFlag_kNeighboursRecovered[0]));
   T->Branch("RecHitFlag_kTowerRecovered_pho1", &(RecHitFlag_kTowerRecovered[0]));
   T->Branch("RecHitFlag_kDead_pho1", &(RecHitFlag_kDead[0]));
   T->Branch("RecHitFlag_kKilled_pho1", &(RecHitFlag_kKilled[0]));
   T->Branch("RecHitFlag_kTPSaturated_pho1", &(RecHitFlag_kTPSaturated[0]));
   T->Branch("RecHitFlag_kL1SpikeFlag_pho1", &(RecHitFlag_kL1SpikeFlag[0]));
   T->Branch("RecHitFlag_kWeird_pho1", &(RecHitFlag_kWeird[0]));
   T->Branch("RecHitFlag_kDiWeird_pho1", &(RecHitFlag_kDiWeird[0]));
   T->Branch("RecHitFlag_kHasSwitchToGain6_pho1", &(RecHitFlag_kHasSwitchToGain6[0]));
   T->Branch("RecHitFlag_kHasSwitchToGain1_pho1", &(RecHitFlag_kHasSwitchToGain1[0]));

   T->Branch("RecHitFlag_kESGood_pho1", &(RecHitFlag_kESGood[0]));
   T->Branch("RecHitFlag_kESDead_pho1", &(RecHitFlag_kESDead[0]));
   T->Branch("RecHitFlag_kESHot_pho1", &(RecHitFlag_kESHot[0]));
   T->Branch("RecHitFlag_kESPassBX_pho1", &(RecHitFlag_kESPassBX[0]));
   T->Branch("RecHitFlag_kESTwoGoodRatios_pho1", &(RecHitFlag_kESTwoGoodRatios[0]));
   T->Branch("RecHitFlag_kESBadRatioFor12_pho1", &(RecHitFlag_kESBadRatioFor12[0]));
   T->Branch("RecHitFlag_kESBadRatioFor23Upper_pho1", &(RecHitFlag_kESBadRatioFor23Upper[0]));
   T->Branch("RecHitFlag_kESBadRatioFor23Lower_pho1", &(RecHitFlag_kESBadRatioFor23Lower[0]));
   T->Branch("RecHitFlag_kESTS1Largest_pho1", &(RecHitFlag_kESTS1Largest[0]));
   T->Branch("RecHitFlag_kESTS3Largest_pho1", &(RecHitFlag_kESTS3Largest[0]));
   T->Branch("RecHitFlag_kESTS3Negative_pho1", &(RecHitFlag_kESTS3Negative[0]));
   T->Branch("RecHitFlag_kESSaturated_pho1", &(RecHitFlag_kESSaturated[0]));
   T->Branch("RecHitFlag_kESTS2Saturated_pho1", &(RecHitFlag_kESTS2Saturated[0]));
   T->Branch("RecHitFlag_kESTS3Saturated_pho1", &(RecHitFlag_kESTS3Saturated[0]));
   T->Branch("RecHitFlag_kESTS13Sigmas_pho1", &(RecHitFlag_kESTS13Sigmas[0]));
   T->Branch("RecHitFlag_kESTS15Sigmas_pho1", &(RecHitFlag_kESTS15Sigmas[0]));

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
   T->Branch("RecHitGain2"  ,  &(RecHitGain[1]));
   T->Branch("RecHitQuality2", &(RecHitQuality[1]));
   T->Branch("HitNoisePho2", &(HitNoise[1]));

   T->Branch("RecHitFlag_kGood_pho2", &(RecHitFlag_kGood[1]));
   T->Branch("RecHitFlag_kPoorReco_pho2", &(RecHitFlag_kPoorReco[1]));
   T->Branch("RecHitFlag_kOutOfTime_pho2", &(RecHitFlag_kOutOfTime[1]));
   T->Branch("RecHitFlag_kFaultyHardware_pho2", &(RecHitFlag_kFaultyHardware[1]));
   T->Branch("RecHitFlag_kNoisy_pho2", &(RecHitFlag_kNoisy[1]));
   T->Branch("RecHitFlag_kPoorCalib_pho2", &(RecHitFlag_kPoorCalib[1]));
   T->Branch("RecHitFlag_kSaturated_pho2", &(RecHitFlag_kSaturated[1]));
   T->Branch("RecHitFlag_kLeadingEdgeRecovered_pho2", &(RecHitFlag_kLeadingEdgeRecovered[1]));
   T->Branch("RecHitFlag_kNeighboursRecovered_pho2", &(RecHitFlag_kNeighboursRecovered[1]));
   T->Branch("RecHitFlag_kTowerRecovered_pho2", &(RecHitFlag_kTowerRecovered[1]));
   T->Branch("RecHitFlag_kDead_pho2", &(RecHitFlag_kDead[1]));
   T->Branch("RecHitFlag_kKilled_pho2", &(RecHitFlag_kKilled[1]));
   T->Branch("RecHitFlag_kTPSaturated_pho2", &(RecHitFlag_kTPSaturated[1]));
   T->Branch("RecHitFlag_kL1SpikeFlag_pho2", &(RecHitFlag_kL1SpikeFlag[1]));
   T->Branch("RecHitFlag_kWeird_pho2", &(RecHitFlag_kWeird[1]));
   T->Branch("RecHitFlag_kDiWeird_pho2", &(RecHitFlag_kDiWeird[1]));
   T->Branch("RecHitFlag_kHasSwitchToGain6_pho2", &(RecHitFlag_kHasSwitchToGain6[1]));
   T->Branch("RecHitFlag_kHasSwitchToGain1_pho2", &(RecHitFlag_kHasSwitchToGain1[1]));

   T->Branch("RecHitFlag_kESGood_pho2", &(RecHitFlag_kESGood[1]));
   T->Branch("RecHitFlag_kESDead_pho2", &(RecHitFlag_kESDead[1]));
   T->Branch("RecHitFlag_kESHot_pho2", &(RecHitFlag_kESHot[1]));
   T->Branch("RecHitFlag_kESPassBX_pho2", &(RecHitFlag_kESPassBX[1]));
   T->Branch("RecHitFlag_kESTwoGoodRatios_pho2", &(RecHitFlag_kESTwoGoodRatios[1]));
   T->Branch("RecHitFlag_kESBadRatioFor12_pho2", &(RecHitFlag_kESBadRatioFor12[1]));
   T->Branch("RecHitFlag_kESBadRatioFor23Upper_pho2", &(RecHitFlag_kESBadRatioFor23Upper[1]));
   T->Branch("RecHitFlag_kESBadRatioFor23Lower_pho2", &(RecHitFlag_kESBadRatioFor23Lower[1]));
   T->Branch("RecHitFlag_kESTS1Largest_pho2", &(RecHitFlag_kESTS1Largest[1]));
   T->Branch("RecHitFlag_kESTS3Largest_pho2", &(RecHitFlag_kESTS3Largest[1]));
   T->Branch("RecHitFlag_kESTS3Negative_pho2", &(RecHitFlag_kESTS3Negative[1]));
   T->Branch("RecHitFlag_kESSaturated_pho2", &(RecHitFlag_kESSaturated[1]));
   T->Branch("RecHitFlag_kESTS2Saturated_pho2", &(RecHitFlag_kESTS2Saturated[1]));
   T->Branch("RecHitFlag_kESTS3Saturated_pho2", &(RecHitFlag_kESTS3Saturated[1]));
   T->Branch("RecHitFlag_kESTS13Sigmas_pho2", &(RecHitFlag_kESTS13Sigmas[1]));
   T->Branch("RecHitFlag_kESTS15Sigmas_pho2", &(RecHitFlag_kESTS15Sigmas[1]));

   T->Branch("nPhotons",  &nPhotons_ , "nPho/I");
   T->Branch("pt"  ,  &Pho_pt_);
   T->Branch("eta" ,  &Pho_eta_ );
   T->Branch("phi" ,  &Pho_phi_ );

   T->Branch("Pho_cluster_seed_x", &Pho_cluster_seed_x);
   T->Branch("Pho_cluster_seed_y", &Pho_cluster_seed_y);
   T->Branch("Pho_cluster_seed_z", &Pho_cluster_seed_z);
   T->Branch("Pho_cluster_seed_eta", &Pho_cluster_seed_eta);
   T->Branch("Pho_cluster_seed_phi", &Pho_cluster_seed_phi);

   T->Branch("Pho_loose_id", &Pho_loose_id);
   T->Branch("Pho_medium_id", &Pho_medium_id);
   T->Branch("Pho_tight_id", &Pho_tight_id);
   T->Branch("Pho_mvawp80_id", &Pho_mvawp80_id);
   T->Branch("Pho_mvawp90_id", &Pho_mvawp90_id);

   T->Branch("energy", &Pho_energy_);
   T->Branch("energy_ecal_mustache", &Pho_ecal_mustache_energy_);

   T->Branch("Pho_R9"  ,  &Pho_R9);
   //T->Branch("Pho_S4"  ,  &Pho_S4);
   T->Branch("Pho_SigIEIE"  ,  &Pho_SigIEIE);
   T->Branch("Pho_SigIPhiIPhi" , &Pho_SigIPhiIPhi);
   T->Branch("Pho_SCEtaW"  ,  &Pho_SCEtaW);
   T->Branch("Pho_SCPhiW"  ,  &Pho_SCPhiW);
   //T->Branch("Pho_CovIEtaIEta"  ,  &Pho_CovIEtaIEta);
   //T->Branch("Pho_CovIEtaIPhi"  ,  &Pho_CovIEtaIPhi);
   //T->Branch("Pho_ESSigRR"  ,  &Pho_ESSigRR);
   T->Branch("Pho_SCRawE"  ,  &Pho_SCRawE);
   T->Branch("Pho_SC_ESEnByRawE"  ,  &Pho_SC_ESEnByRawE);
   T->Branch("Pho_HadOverEm"  ,  &Pho_HadOverEm);

   T->Branch("Pho_PFChIso"	,	&Pho_PFChIso);
   T->Branch("Pho_PFChPVIso"	,	&Pho_PFChPVIso);
   T->Branch("Pho_PFPhoIso"	,	&Pho_PFPhoIso);
   T->Branch("Pho_PFNeuIso"	,	&Pho_PFNeuIso);
   T->Branch("Pho_PFChWorstVetoIso"	,	&Pho_PFChWorstVetoIso);
   T->Branch("Pho_PFChWorstIso"	,	&Pho_PFChWorstIso);
   T->Branch("Pho_EcalPFClusterIso"	,	&Pho_EcalPFClusterIso);
   T->Branch("Pho_HcalPFClusterIso"	,	 &Pho_HcalPFClusterIso);

   T->Branch("Pho_CorrectedEnergy", &Pho_CorrectedEnergy); // Add corrected energy for photon 
   T->Branch("Pho_CorrectedEnergyError", &Pho_CorrectedEnergyError); //Error in energy correction

   if (isMC_){
      T->Branch("Pho_Gen_Pt" , &Pho_Gen_Pt);
      T->Branch("Pho_Gen_Eta" , &Pho_Gen_Eta);
      T->Branch("Pho_Gen_Phi" , &Pho_Gen_Phi);
      T->Branch("Pho_Gen_E" , &Pho_Gen_E);

      T->Branch("Pho_GenIdx", &Pho_GenIdx);
      T->Branch("Pho_DR", &Pho_DR);
   }

   T->Branch("rho", &rho, "rho/F");

   T->Branch("run",&run,"run/I");
   T->Branch("event",&event,"event/I");
   T->Branch("lumi",&lumi,"lumi/I");

}

// ------------ method called once each job just after ending the event loop  ------------
   void
Photon_RefinedRecHit_MiniAOD_NTuplizer::endJob()
{
}

/*
//Evaluate if the gen particle dR matched to a reco photon is also a photon
bool Photon_RefinedRecHit_MiniAOD_NTuplizer::GetGenMatchType(const reco::Photon& Photon, const reco::GenParticle& GenColl, int pdgId, double dRThresh){


}
*/


// Extract the rechits of the SC from the ES layers
void Photon_RefinedRecHit_MiniAOD_NTuplizer::GetESPlaneRecHits(const reco::SuperCluster& sc, const CaloGeometry* &geo, unsigned int phonum, unsigned int planeIndex) {
   double RawenergyPlane = 0.;
   double pfRawenergyPlane = 0.;
   // if(!_ESRechitsHandle.isValid())
   //    return RawenergyPlane;

   if (!sc.preshowerClusters().isAvailable()) //protection for miniAOD
       return;

   int NumHits = 0;

   const CaloSubdetectorGeometry* ecalESGeom = static_cast<const CaloSubdetectorGeometry*>(geo->getSubdetectorGeometry(DetId::Ecal, EcalPreshower));


   for(auto iES = sc.preshowerClustersBegin(); iES != sc.preshowerClustersEnd(); ++iES) {//loop over preshower clusters
      const std::vector< std::pair<DetId, float> > hits = (*iES)->hitsAndFractions();
      for(std::vector<std::pair<DetId, float> >::const_iterator rh = hits.begin(); rh != hits.end(); ++rh) { // loop over recHits of the cluster
         //std::cout << "print = " << (*iES)->printHitAndFraction(iCount);
         //++iCount;
         for(ESRecHitCollection::const_iterator esItr = ESRechitsHandle->begin(); esItr != ESRechitsHandle->end(); ++esItr) {
            //loop over ES rechits to find the one in the cluster
            ESDetId rhid = ESDetId(esItr->id());
            if(rhid == (*rh).first) { // found ESrechit
               //std::cout << " ES energy = " << esItr->energy() << " pf energy = " << (*rh).second << std::endl;
               if((int) rhid.plane() == (int) planeIndex) {
                  std::shared_ptr<const CaloCellGeometry> geom = ecalESGeom->getGeometry(rhid);
                  Hit_ES_Eta[phonum].push_back( geom->etaPos() );
                  Hit_ES_Phi[phonum].push_back( geom->phiPos() );
                  Hit_ES_X[phonum].push_back( geom->getPosition().x() );
                  Hit_ES_Y[phonum].push_back( geom->getPosition().y() );
                  Hit_ES_Z[phonum].push_back( geom->getPosition().z() ) ;
                  ES_RecHitEn[phonum].push_back(esItr->energy());

                  for (int iflag=EcalRecHit::kESGood; iflag<EcalRecHit::kESTS15Sigmas+1; iflag++){
                     bool check_bit = esItr->checkFlag(iflag);
                     RecHitESFlag_container[iflag][phonum].push_back(check_bit);

                     if (DEBUG) cout<< "ES Flag: "<<iflag<<endl;
                  }
                  //std::cout << "Preshower" <<std::setprecision(4) << " Eta = " <<geom->etaPos() << " : " <<" Phi = "<< geom->phiPos() << " 3D position" << geom->getPosition().z() << std::endl;
                  RawenergyPlane += esItr->energy();
                  pfRawenergyPlane += rh->second;
                  NumHits++;
               }
               break;
            }
         }
      }

      //std::cout<<std::endl<<" Number of ES hits in plane 1 = "<<NumHits<<std::endl;
   }

   //return RawenergyPlane;
}


//Clear tree vectors each time analyze method is called
void Photon_RefinedRecHit_MiniAOD_NTuplizer::ClearTreeVectors()
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
   RecHitGain[0].clear();
   RecHitQuality[0].clear();
   HitNoise[0].clear();
   iEta[1].clear();
   iPhi[1].clear();

   RecHitFlag_kGood[0].clear();
   RecHitFlag_kPoorReco[0].clear();
   RecHitFlag_kOutOfTime[0].clear();
   RecHitFlag_kFaultyHardware[0].clear();
   RecHitFlag_kNoisy[0].clear();
   RecHitFlag_kPoorCalib[0].clear();
   RecHitFlag_kSaturated[0].clear();
   RecHitFlag_kLeadingEdgeRecovered[0].clear();
   RecHitFlag_kNeighboursRecovered[0].clear();
   RecHitFlag_kTowerRecovered[0].clear();
   RecHitFlag_kDead[0].clear();
   RecHitFlag_kKilled[0].clear();
   RecHitFlag_kTPSaturated[0].clear();
   RecHitFlag_kL1SpikeFlag[0].clear();
   RecHitFlag_kWeird[0].clear();
   RecHitFlag_kDiWeird[0].clear();
   RecHitFlag_kHasSwitchToGain6[0].clear();
   RecHitFlag_kHasSwitchToGain1[0].clear();

   RecHitFlag_kESGood[0].clear();
   RecHitFlag_kESDead[0].clear();
   RecHitFlag_kESHot[0].clear();
   RecHitFlag_kESPassBX[0].clear();
   RecHitFlag_kESTwoGoodRatios[0].clear();
   RecHitFlag_kESBadRatioFor12[0].clear();
   RecHitFlag_kESBadRatioFor23Upper[0].clear();
   RecHitFlag_kESBadRatioFor23Lower[0].clear();
   RecHitFlag_kESTS1Largest[0].clear();
   RecHitFlag_kESTS3Largest[0].clear();
   RecHitFlag_kESTS3Negative[0].clear();
   RecHitFlag_kESSaturated[0].clear();
   RecHitFlag_kESTS2Saturated[0].clear();
   RecHitFlag_kESTS3Saturated[0].clear();
   RecHitFlag_kESTS13Sigmas[0].clear();
   RecHitFlag_kESTS15Sigmas[0].clear();

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
   RecHitGain[1].clear();
   RecHitQuality[1].clear();
   HitNoise[1].clear();

   RecHitFlag_kGood[1].clear();
   RecHitFlag_kPoorReco[1].clear();
   RecHitFlag_kOutOfTime[1].clear();
   RecHitFlag_kFaultyHardware[1].clear();
   RecHitFlag_kNoisy[1].clear();
   RecHitFlag_kPoorCalib[1].clear();
   RecHitFlag_kSaturated[1].clear();
   RecHitFlag_kLeadingEdgeRecovered[1].clear();
   RecHitFlag_kNeighboursRecovered[1].clear();
   RecHitFlag_kTowerRecovered[1].clear();
   RecHitFlag_kDead[1].clear();
   RecHitFlag_kKilled[1].clear();
   RecHitFlag_kTPSaturated[1].clear();
   RecHitFlag_kL1SpikeFlag[1].clear();
   RecHitFlag_kWeird[1].clear();
   RecHitFlag_kDiWeird[1].clear();
   RecHitFlag_kHasSwitchToGain6[1].clear();
   RecHitFlag_kHasSwitchToGain1[1].clear();

   RecHitFlag_kESGood[1].clear();
   RecHitFlag_kESDead[1].clear();
   RecHitFlag_kESHot[1].clear();
   RecHitFlag_kESPassBX[1].clear();
   RecHitFlag_kESTwoGoodRatios[1].clear();
   RecHitFlag_kESBadRatioFor12[1].clear();
   RecHitFlag_kESBadRatioFor23Upper[1].clear();
   RecHitFlag_kESBadRatioFor23Lower[1].clear();
   RecHitFlag_kESTS1Largest[1].clear();
   RecHitFlag_kESTS3Largest[1].clear();
   RecHitFlag_kESTS3Negative[1].clear();
   RecHitFlag_kESSaturated[1].clear();
   RecHitFlag_kESTS2Saturated[1].clear();
   RecHitFlag_kESTS3Saturated[1].clear();
   RecHitFlag_kESTS13Sigmas[1].clear();
   RecHitFlag_kESTS15Sigmas[1].clear();

   Pho_pt_.clear();
   Pho_eta_.clear();
   Pho_phi_.clear();

   Pho_cluster_seed_x.clear();
   Pho_cluster_seed_y.clear();
   Pho_cluster_seed_z.clear();
   Pho_cluster_seed_eta.clear();
   Pho_cluster_seed_phi.clear();

   Pho_loose_id.clear();
   Pho_medium_id.clear();
   Pho_tight_id.clear();
   Pho_mvawp80_id.clear();
   Pho_mvawp90_id.clear();

   Pho_energy_.clear();
   Pho_ecal_mustache_energy_.clear();
   Pho_R9.clear();
   //Pho_S4.clear();
   Pho_SigIEIE.clear();
   Pho_SigIPhiIPhi.clear();
   Pho_SCEtaW.clear();
   Pho_SCPhiW.clear();
   //Pho_CovIEtaIEta.clear();
   //Pho_CovIEtaIPhi.clear();
   //Pho_ESSigRR.clear();
   Pho_SCRawE.clear();
   Pho_SC_ESEnByRawE.clear();
   Pho_HadOverEm.clear();

   Pho_PFChIso.clear();
   Pho_PFChPVIso.clear();
   Pho_PFPhoIso.clear();
   Pho_PFNeuIso.clear();
   Pho_PFChWorstVetoIso.clear();
   Pho_PFChWorstIso.clear();
   Pho_EcalPFClusterIso.clear();
   Pho_HcalPFClusterIso.clear();

   Pho_CorrectedEnergy.clear();
   Pho_CorrectedEnergyError.clear();

   if (isMC_){
      Pho_Gen_Pt.clear();
      Pho_Gen_Eta.clear();
      Pho_Gen_Phi.clear();
      Pho_Gen_E.clear();

      Pho_GenIdx.clear();
      Pho_DR.clear();
   }

   isTrue_.clear();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Photon_RefinedRecHit_MiniAOD_NTuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(Photon_RefinedRecHit_MiniAOD_NTuplizer);
