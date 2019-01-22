#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CustoTnP/Analyzer/src/PATUtilities.h"
#include "CustoTnP/Analyzer/src/TriggerUtilities.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "TLorentzVector.h"

//~
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTSuperLayer.h"
#include "DataFormats/DTRecHit/interface/DTSLRecSegment2D.h"
#include "RecoLocalMuon/DTSegment/src/DTSegmentUpdator.h"
#include "RecoLocalMuon/DTSegment/src/DTSegmentCleaner.h"
#include "RecoLocalMuon/DTSegment/src/DTHitPairForFit.h"

#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/CSCRecHit/interface/CSCRangeMapAccessor.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/MuonReco/interface/MuonShower.h"


class CustoTnPLeptonProducer : public edm::EDProducer {
public:
  explicit CustoTnPLeptonProducer(const edm::ParameterSet&);

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  pat::Muon*     cloneAndSwitchMuonTrack     (const pat::Muon&, const edm::Event& event)     const;

  void embedTriggerMatch(pat::Muon*, std::string, const pat::TriggerObjectStandAloneCollection&, std::vector<int>&);
  // void embedTriggerMatch_or(pat::Muon*, const std::string&, const pat::TriggerObjectStandAloneCollection&, const pat::TriggerObjectStandAloneCollection&, std::vector<int>&, std::vector<int>&);

  //~
  void embedExpectedMatchedStations(pat::Muon*, float);

  std::pair<pat::Muon*,     int> doLepton(const edm::Event&, const pat::Muon&,     const reco::CandidateBaseRef&);

  template <typename T> edm::OrphanHandle<std::vector<T> > doLeptons(edm::Event&, const edm::InputTag&, const edm::InputTag&, const std::string&);

  edm::InputTag muon_src;
  edm::InputTag muon_view_src;
  StringCutObjectSelector<pat::Muon> muon_selector;
  std::string muon_track_for_momentum;
  std::string muon_track_for_momentum_primary;
  std::vector<std::string> muon_tracks_for_momentum;
  edm::InputTag muon_photon_match_src;
  edm::Handle<reco::CandViewMatchMap> muon_photon_match_map;

  double trigger_match_max_dR;
  std::vector<std::string> trigger_filters;
  std::vector<std::string> trigger_path_names;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigger_summary_src_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

  std::vector<pat::TriggerObjectStandAloneCollection> vec_L3_muons;
  std::vector<std::vector<int>> vec_L3_muons_matched;

  edm::EDGetTokenT<l1t::MuonBxCollection> l1_src_;

  //@
  bool isAOD;
  edm::InputTag reco_muon_src;
  edm::InputTag muonshower_src;
  edm::InputTag dtseg_src;
  edm::InputTag cscseg_src;

  std::pair<bool, reco::MuonRef> getMuonRef(const edm::Event&, pat::Muon*);
  std::vector<int> countDTsegs(const edm::Event&, reco::MuonRef);
  std::vector<int> countCSCsegs(const edm::Event&, reco::MuonRef);
  void embedShowerInfo(const edm::Event&, pat::Muon*, reco::MuonRef);
};

CustoTnPLeptonProducer::CustoTnPLeptonProducer(const edm::ParameterSet& cfg)
  : muon_src(cfg.getParameter<edm::InputTag>("muon_src")),
    muon_view_src(cfg.getParameter<edm::InputTag>("muon_src")),
    muon_selector(cfg.getParameter<std::string>("muon_cuts")),
    muon_track_for_momentum(cfg.getParameter<std::string>("muon_track_for_momentum")),
    muon_track_for_momentum_primary(muon_track_for_momentum),
    muon_photon_match_src(cfg.getParameter<edm::InputTag>("muon_photon_match_src")),

    trigger_match_max_dR(cfg.getParameter<double>("trigger_match_max_dR")),
    trigger_filters(cfg.getParameter<std::vector<std::string>>("trigger_filters")),
    trigger_path_names(cfg.getParameter<std::vector<std::string>>("trigger_path_names")),

    triggerBits_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("bits"))),
    trigger_summary_src_(consumes<pat::TriggerObjectStandAloneCollection>(cfg.getParameter<edm::InputTag>("trigger_summary"))),
    triggerPrescales_(consumes<pat::PackedTriggerPrescales>(cfg.getParameter<edm::InputTag>("prescales"))),
    l1_src_(consumes<l1t::MuonBxCollection>(cfg.getParameter<edm::InputTag>("l1"))),

    isAOD(cfg.getParameter<bool>("isAOD"))
{
  consumes<edm::View<reco::Candidate>>(muon_view_src);
  consumes<pat::MuonCollection>(muon_src);
  consumes<reco::CandViewMatchMap >(muon_photon_match_src);

  if (cfg.existsAs<std::vector<std::string> >("muon_tracks_for_momentum"))
    muon_tracks_for_momentum = cfg.getParameter<std::vector<std::string> >("muon_tracks_for_momentum");

  if (muon_tracks_for_momentum.size())
    for (size_t i = 0; i < muon_tracks_for_momentum.size(); ++i)
      produces<pat::MuonCollection>(muon_tracks_for_momentum[i]);

  if(isAOD) {
    reco_muon_src  = (cfg.getParameter<edm::InputTag>("reco_muon_src"));
    muonshower_src = (cfg.getParameter<edm::InputTag>("muonshower_src"));
    dtseg_src      = (cfg.getParameter<edm::InputTag>("dtseg_src"));
    cscseg_src     = (cfg.getParameter<edm::InputTag>("cscseg_src"));

    consumes<std::vector< reco::Muon >>(reco_muon_src);
    consumes<edm::ValueMap<reco::MuonShower>>(muonshower_src);
    consumes<DTRecSegment4DCollection>(dtseg_src);
    consumes<CSCSegmentCollection>(cscseg_src);
  }

  produces<pat::MuonCollection>("muons");
}


pat::Muon* CustoTnPLeptonProducer::cloneAndSwitchMuonTrack(const pat::Muon& muon, const edm::Event& event) const {

  // Muon mass to make a four-vector out of the new track.

  pat::Muon* mu = muon.clone();

  // Start with null track/invalid type before we find the right one.
  reco::TrackRef newTrack;
  newTrack = muon.tunePMuonBestTrack();
  patmuon::TrackType type = patmuon::nTrackTypes;

  // If the muon has the track embedded using the UserData mechanism,
  // take it from there first. Otherwise, try to get the track the
  // standard way.

  if (muon.pt() > 100.) {
	mu->addUserInt("hasTeVMuons", 1);
  }
  else{

	 mu->addUserInt("hasTeVMuons", 0); 
  }
  if (muon.hasUserData(muon_track_for_momentum))
    newTrack = patmuon::userDataTrack(muon, muon_track_for_momentum);
  else {
 
	type = patmuon::trackNameToType(muon_track_for_momentum);
    	newTrack = patmuon::trackByType(*mu, type);
  }
  
  // If we didn't find the appropriate track, indicate failure by a
  // null pointer.
  
  if (newTrack.isNull()){
    std::cout << "No TuneP" << std::endl;
    return 0;
    
  }
  
  
  
  static const double mass = 0.10566;
  
  reco::Particle::Point vtx(newTrack->vx(), newTrack->vy(), newTrack->vz());
  reco::Particle::LorentzVector p4;

  //////////   Comment following lines to apply pt bias correction /////
    const double p = newTrack->p();  
    p4.SetXYZT(newTrack->px(), newTrack->py(), newTrack->pz(), sqrt(p*p + mass*mass));  
  //////////   Comment previous lines to apply pt bias correction ----->  Uncomment following lines /////


  ///////// uncomment following lines to apply pt bias correction -----> comment previous lines /////////
    //  float phi = newTrack->phi()*TMath::RadToDeg();

    //  float mupt = GeneralizedEndpoint().GeneralizedEndpointPt(newTrack->pt(),newTrack->charge(),newTrack->eta(),phi,-1,1); //for DATA
    //  float mupt = GeneralizedEndpoint().GeneralizedEndpointPt(newTrack->pt(),newTrack->charge(),newTrack->eta(),phi,0,1);  // for MC


    //	float px = mupt*TMath::Cos(newTrack->phi());
    //	float py = mupt*TMath::Sin(newTrack->phi());
    //	float pz = mupt*TMath::SinH(newTrack->eta());
    //	float p = mupt*TMath::CosH(newTrack->eta());
    //	p4.SetXYZT(px, py, pz, sqrt(p*p + mass*mass));

    // 	std::cout<<"my definition = "<<mupt<<std::endl;
  /////// uncomment previous lines to apply pt bias correction /////////


  mu->setP4(p4);  

  mu->setCharge(newTrack->charge());

  mu->setVertex(vtx);

  mu->addUserInt("trackUsedForMomentum", type);
  
  return mu;
}


void CustoTnPLeptonProducer::embedTriggerMatch(pat::Muon* new_mu, std::string ex, const pat::TriggerObjectStandAloneCollection& L3, std::vector<int>& L3_matched) {
  
  int best = -1;
  float defaultpTvalue = -1.;
  float best_dR = trigger_match_max_dR;
  //std::cout << "size of L3 collection: " << L3.size() << std::endl;
  for (size_t i = 0; i < L3.size(); ++i) {
    // Skip those already used.
    if (L3_matched[i])
      continue;

    const float dR = reco::deltaR(L3[i], *new_mu);
    if (dR < best_dR) {
      best = int(i);
      best_dR = dR;
    }
  }

  // if (best < 0)
  //  return;
  if (ex.length()>0) ex += "_";
  if (best >= 0) {
    const pat::TriggerObjectStandAlone& L3_mu = L3[best];
    L3_matched[best] = 1;
    
    int id = L3_mu.pdgId();
    new_mu->addUserFloat(ex + "TriggerMatchCharge", -id/abs(id));
    new_mu->addUserFloat(ex + "TriggerMatchPt",     L3_mu.pt());
    new_mu->addUserFloat(ex + "TriggerMatchEta",    L3_mu.eta());
    new_mu->addUserFloat(ex + "TriggerMatchPhi",    L3_mu.phi());
  }
  else {
    new_mu->addUserFloat(ex + "TriggerMatchPt",    defaultpTvalue);
  }

}

//~
void CustoTnPLeptonProducer::embedExpectedMatchedStations(pat::Muon* new_mu, float minDistanceFromEdge = 10.) {
  unsigned int stationMask = 0;
  for( auto& chamberMatch : new_mu->matches() )
  {
    if (chamberMatch.detector()!=MuonSubdetId::DT && chamberMatch.detector()!=MuonSubdetId::CSC) continue;
    float edgeX = chamberMatch.edgeX;
    float edgeY = chamberMatch.edgeY;
    // check we if the trajectory is well within the acceptance
    /*
    if(chamberMatch.detector()==MuonSubdetId::DT) {
      DTChamberId ID( chamberMatch.id.rawId() );
      std::cout << "\tChamberId: " << ID << std::endl;
    }
    if(chamberMatch.detector()==MuonSubdetId::CSC) {
      CSCDetId    ID( chamberMatch.id.rawId() );
      std::cout << "\tChamberId: " << ID << std::endl;
    }
    std::cout << "\t\tedgeX: " << edgeX << "\tedgeY: " << edgeY << "\t" << (edgeX<0 && fabs(edgeX)>fabs(minDistanceFromEdge) && edgeY<0 && fabs(edgeY)>fabs(minDistanceFromEdge)) <<std::endl;
    */
    if(edgeX<0 && fabs(edgeX)>fabs(minDistanceFromEdge) && edgeY<0 && fabs(edgeY)>fabs(minDistanceFromEdge))
      stationMask |= 1<<( (chamberMatch.station()-1)+4*(chamberMatch.detector()-1) );
  }
  unsigned int n = 0;
  for(unsigned int i=0; i<8; ++i)
    if (stationMask&(1<<i)) n++;

  std::string var_temp = "expectedNnumberOfMatchedStations"+std::to_string( int(minDistanceFromEdge) );
  new_mu->addUserInt(var_temp, n);
}

//@
std::pair<bool, reco::MuonRef> CustoTnPLeptonProducer::getMuonRef(const edm::Event& event, pat::Muon* new_mu) {

  edm::Handle< std::vector< reco::Muon > > recoMuons;
  event.getByLabel(reco_muon_src, recoMuons);

  bool isMatched = false;
  reco::MuonRef matchedMuRef = reco::MuonRef(recoMuons, 0);

  int imucount = 0;
  for (std::vector<reco::Muon>::const_iterator imu = recoMuons->begin(); imu != recoMuons->end(); imu++) {
    if(imu->globalTrack()->pt() == new_mu->globalTrack()->pt() &&
       imu->globalTrack()->eta() == new_mu->globalTrack()->eta() &&
       imu->globalTrack()->phi() == new_mu->globalTrack()->phi() ) {
      isMatched = true;
      matchedMuRef = reco::MuonRef(recoMuons, imucount);
      break;
    }
    imucount++;
  }

  return make_pair(isMatched, matchedMuRef);
}

std::vector<int> CustoTnPLeptonProducer::countDTsegs(const edm::Event& event, reco::MuonRef muon) {
  double DTCut = 30.;

  std::vector<int> stations={0,0,0,0};
  std::vector<int> removed={0,0,0,0};

  edm::Handle<DTRecSegment4DCollection> dtSegments;
  event.getByLabel(dtseg_src, dtSegments);

  if (!ShutUp) std::cout << std::endl << " *** DT Segment search" << std::endl;
  for (const auto &ch : muon->matches()) {
    if( ch.detector() != MuonSubdetId::DT )  continue;
    DTChamberId DTid( ch.id.rawId() );
    if (!ShutUp) std::cout << "   DT chamber in station " << ch.station() << "  DTChamberId:" << DTid << "  local position: (" << ch.x << ", " << ch.y << ", 0)" << std::endl;

    int nsegs_temp = 0;
    std::vector<float> nsegs_x_temp, nsegs_y_temp;

    for (auto seg = dtSegments->begin(); seg!=dtSegments->end(); ++seg) {
      DTChamberId myChamber((*seg).geographicalId().rawId());
      if (!(DTid==myChamber))  continue;
      LocalPoint posLocalSeg = seg->localPosition();
      if (!ShutUp) std::cout << "     Found segment in station " << ch.station() << "  DTChamberId:" << myChamber << "  local position: " << posLocalSeg << std::endl;

      /*if( ( (posLocalSeg.x()==0 && posLocalSeg.y()!=0) && (fabs(posLocalSeg.y()-ch.y)<DTCut) ) ||
          ( (posLocalSeg.x()!=0 && posLocalSeg.y()==0) && (fabs(posLocalSeg.x()-ch.x)<DTCut) ) ||
          ( (posLocalSeg.x()!=0 && posLocalSeg.y()!=0) && (sqrt( (posLocalSeg.x()-ch.x)*(posLocalSeg.x()-ch.x) + (posLocalSeg.y()-ch.y)*(posLocalSeg.y()-ch.y) )<DTCut) )
        ) {
        nsegs_temp++;
      }*/

      if( ( posLocalSeg.x()!=0 && ch.x!=0 ) && (fabs(posLocalSeg.x()-ch.x)<DTCut) ) {
        bool found = false;
        for( auto prev_x : nsegs_x_temp) {
          if( fabs(prev_x-posLocalSeg.x()) < 0.1 ) {
            found = true;
            break;
          }
        }
        if( !found )  nsegs_x_temp.push_back(posLocalSeg.x());
      }

      if( ( posLocalSeg.y()!=0 && ch.y!=0 ) && (fabs(posLocalSeg.y()-ch.y)<DTCut) ) {
        bool found = false;
        for( auto prev_y : nsegs_y_temp) {
          if( fabs(prev_y-posLocalSeg.y()) < 0.1 ) {
            found = true;
            break;
          }
        }
        if( !found )  nsegs_y_temp.push_back(posLocalSeg.y());
      }

    }

    nsegs_temp = (int)std::max(nsegs_x_temp.size(), nsegs_y_temp.size());

    //--- subtract best matched segment from given muon
    bool isBestMatched = false;
    for(std::vector<reco::MuonSegmentMatch>::const_iterator matseg = ch.segmentMatches.begin(); matseg != ch.segmentMatches.end(); matseg++) {
      if( matseg->isMask(reco::MuonSegmentMatch::BestInChamberByDR) ) {
        if (!ShutUp) std::cout << "     Found BestInChamberByDR in station " << ch.station() << "  DTChamberId:" << DTid << "  local position: (" << matseg->x << ", " << matseg->y << ", 0)" << std::endl;
        removed[ch.station()-1]++;
        isBestMatched = true;
        break;
      }
    }
    if (!ShutUp) std::cout << std::endl;

    if(isBestMatched) nsegs_temp = nsegs_temp-1;
    if(nsegs_temp>0)  stations[ch.station()-1] += nsegs_temp;

  }

  if (!ShutUp) {
    std::cout << " DT Shower pattern: ";
    int DTSum = 0;
    for (int i=0;i<4;i++) {
      std::cout << stations[i] << " ";
      DTSum += stations[i];
    }
    std::cout << std::endl;
    std::cout << " DTSum = " << DTSum << std::endl;

    std::cout << " DT removed pattern: ";
    for (int i=0;i<4;i++) {
      std::cout << removed[i] << " ";
    }
    std::cout << std::endl;
  }

  return stations;
}

std::vector<int> CustoTnPLeptonProducer::countCSCsegs(const edm::Event& event, reco::MuonRef muon) {
  double CSCCut = 30.;

  std::vector<int> stations={0,0,0,0};
  std::vector<int> removed={0,0,0,0};

  edm::Handle<CSCSegmentCollection> cscSegments;
  event.getByLabel(cscseg_src, cscSegments);

  if (!ShutUp) std::cout << std::endl << " *** CSC Segment search" << std::endl;
  for (const auto &ch : muon->matches()) {
    if( ch.detector() != MuonSubdetId::CSC )  continue;
    CSCDetId CSCid( ch.id.rawId() );
    if (!ShutUp) std::cout << "   CSC chamber in station " << ch.station() << "  CSCDetId:" << CSCid << "  local position: (" << ch.x << ", " << ch.y << ", 0)" << std::endl;

    int nsegs_temp = 0;

    for (auto seg = cscSegments->begin(); seg!=cscSegments->end(); ++seg) {
      CSCDetId myChamber((*seg).geographicalId().rawId());
      if (!(CSCid==myChamber))  continue;
      LocalPoint posLocalSeg = seg->localPosition();
      if (!ShutUp) std::cout << "     Found segment in station " << ch.station() << "  CSCDetId:" << myChamber << "  local position: " << posLocalSeg << std::endl;

      if( (posLocalSeg.x()!=0 && posLocalSeg.y()!=0) && (sqrt( (posLocalSeg.x()-ch.x)*(posLocalSeg.x()-ch.x) + (posLocalSeg.y()-ch.y)*(posLocalSeg.y()-ch.y) )<CSCCut) )  { 
        nsegs_temp++;
      }
    }

    //--- subtract best matched segment from given muon
    bool isBestMatched = false;
    for(std::vector<reco::MuonSegmentMatch>::const_iterator matseg = ch.segmentMatches.begin(); matseg != ch.segmentMatches.end(); matseg++) {
      if( matseg->isMask(reco::MuonSegmentMatch::BestInChamberByDR) ) {
        if (!ShutUp) std::cout << "     Found BestInChamberByDR in station " << ch.station() << "  CSCDetId:" << CSCid << "  local position: (" << matseg->x << ", " << matseg->y << ", 0)" << std::endl;
        removed[ch.station()-1]++;
        isBestMatched = true;
        break;
      }
    }
    if (!ShutUp) std::cout << std::endl;

    if(isBestMatched) nsegs_temp = nsegs_temp-1;
    if(nsegs_temp>0)  stations[ch.station()-1] += nsegs_temp;

    /*for(std::vector<reco::MuonSegmentMatch>::const_iterator seg = ch.segmentMatches.begin(); seg != ch.segmentMatches.end(); seg++) {
      if( seg->isMask(reco::MuonSegmentMatch::BestInChamberByDR) ) {
        if (!ShutUp) std::cout << "     Found BestInChamberByDR in station " << ch.station() << "  local position: (" << seg->x << ", " << seg->y << ", 0)" << std::endl;
        removed[ch.station()-1]++;
      }
      else {
        if (!ShutUp) std::cout << "     Found segment in station " << ch.station() << "  local position: (" << seg->x << ", " << seg->y << ", 0)" << std::endl;
        if( (seg->x!=0 && seg->y!=0) && (sqrt( (seg->x-ch.x)*(seg->x-ch.x) + (seg->y-ch.y)*(seg->y-ch.y) )<CSCCut) )  { stations[ch.station()-1]++; }
      }
    }*/
  }

  if (!ShutUp) {
    std::cout << " CSC Shower pattern: ";
    int CSCSum = 0;
    for (int i=0;i<4;i++) {
      std::cout << stations[i] << " ";
      CSCSum += stations[i];
    }
    std::cout << std::endl;
    std::cout << " CSCSum = " << CSCSum << std::endl;

    std::cout << " CSC removed pattern: ";
    for (int i=0;i<4;i++) {
      std::cout << removed[i] << " ";
    }
    std::cout << std::endl;
  }

  return stations;
}

void CustoTnPLeptonProducer::embedShowerInfo(const edm::Event& event, pat::Muon* new_mu, reco::MuonRef MuRef) {
  edm::Handle<edm::ValueMap<reco::MuonShower> > muonShowerInformationValueMap;
  event.getByLabel(muonshower_src, muonShowerInformationValueMap);

  reco::MuonShower muonShowerInformation = (*muonShowerInformationValueMap)[MuRef];
  std::vector<int> vec_DTsegs = countDTsegs(event, MuRef);
  std::vector<int> vec_CSCsegs = countDTsegs(event, MuRef);

  for(int i=0; i<4; ++i) {
    int nsegs_DT_temp  = vec_DTsegs[i];
    int nsegs_CSC_temp = vec_CSCsegs[i];
    int nhits_temp = (muonShowerInformation.nStationHits).at(i);

    std::string var_segs_DT_temp  = "nSegsDT"+std::to_string( int(i+1) );
    std::string var_segs_CSC_temp = "nSegsCSC"+std::to_string( int(i+1) );
    std::string var_hits_temp     = "nHits"+std::to_string( int(i+1) );

    new_mu->addUserInt(var_segs_DT_temp,  nsegs_DT_temp);
    new_mu->addUserInt(var_segs_CSC_temp, nsegs_CSC_temp);
    new_mu->addUserInt(var_hits_temp,     nhits_temp);
  }

}

std::pair<pat::Muon*,int> CustoTnPLeptonProducer::doLepton(const edm::Event& event, const pat::Muon& mu, const reco::CandidateBaseRef& cand) {
  // Failure is indicated by a null pointer as the first member of the
  // pair.

  // To use one of the refit tracks, we have to have a global muon.
      
  if (!mu.isGlobalMuon())
    return std::make_pair((pat::Muon*)(0), -1);

  // Copy the input muon, and switch its p4/charge/vtx out for that of
  // the selected refit track.
  
  pat::Muon* new_mu = cloneAndSwitchMuonTrack(mu, event);

  if (new_mu == 0){
    return std::make_pair(new_mu, -1);
    //std::cout << "Warning" << std::endl;
  }  

  // Simply store the photon four-vector for now in the muon as a
  // userData.
  if (muon_photon_match_map.isValid()) {
    const reco::CandViewMatchMap& mm = *muon_photon_match_map;
    if (mm.find(cand) != mm.end()) {
      new_mu->addUserData<reco::Particle::LorentzVector>("photon_p4", mm[cand]->p4());
      new_mu->addUserInt("photon_index", mm[cand].key());
    }    
  }

  //--- Trig match
  for(unsigned i_f=0; i_f<trigger_filters.size(); ++i_f)
    embedTriggerMatch(new_mu, trigger_path_names[i_f], vec_L3_muons[i_f], vec_L3_muons_matched[i_f]);

  //--- L1 match
  edm::Handle<l1t::MuonBxCollection> l1_src;
  event.getByToken(l1_src_, l1_src);

  float the_dr  = 999.;
  float the_pt  = -1.;
  float the_eta = -999.;
  float the_phi = -999.;
  int the_q     = -1;
  for(int ibx = l1_src->getFirstBX(); ibx<=l1_src->getLastBX(); ++ibx) {
    if(ibx != 0) continue;
    for(auto it=l1_src->begin(ibx); it!=l1_src->end(ibx); it++) {
      l1t::MuonRef ref_L1Mu(l1_src, distance(l1_src->begin(l1_src->getFirstBX()), it) );
      if( ref_L1Mu->pt() < 22. ) // L1 pT>22
        continue;
      float dr_temp = reco::deltaR(*ref_L1Mu, *new_mu);
      if( dr_temp > 0.3 )
        continue;
      if( dr_temp < the_dr ) {
        the_dr = dr_temp;
        the_pt = ref_L1Mu->pt();
        the_eta = ref_L1Mu->eta();
        the_phi = ref_L1Mu->phi();
        the_q = ref_L1Mu->hwQual();
      }
    }
  }
  if(the_pt > 0) {
    new_mu->addUserFloat("L122MatchPt",     the_pt);
    new_mu->addUserFloat("L122MatchEta",    the_eta);
    new_mu->addUserFloat("L122MatchPhi",    the_phi);
    new_mu->addUserInt("L122MatchQ",        the_q);
  }
  else
    new_mu->addUserFloat("L122MatchPt",     the_pt);

  //~
  // std::cout << "Muon pT= " << new_mu->pt() << "\teta= " << new_mu->eta() << "\tphi= " << new_mu->phi() << std::endl;
  // embedExpectedMatchedStations(new_mu, 5.);
  embedExpectedMatchedStations(new_mu, 10.);
  // embedExpectedMatchedStations(new_mu, 15.);
  // embedExpectedMatchedStations(new_mu, 20.);
  // std::cout << "\texpectedNnumberOfMatchedStations: " << new_mu->userInt("expectedNnumberOfMatchedStations") << std::endl;

  //@
  if(isAOD) {
    std::pair<bool, reco::MuonRef> pair_muRef = getMuonRef(event, new_mu);
    if( pair_muRef.first ) {
      embedShowerInfo(event, new_mu, pair_muRef.second);
      // std::cout << "\tembedShowerInfo: nHits1=" << new_mu->userInt("nHits1") << std::endl;
      // std::cout << "\tembedShowerInfo: nHits2=" << new_mu->userInt("nHits2") << std::endl;
      // std::cout << "\tembedShowerInfo: nHits3=" << new_mu->userInt("nHits3") << std::endl;
      // std::cout << "\tembedShowerInfo: nHits4=" << new_mu->userInt("nHits4") << std::endl;
    }
    else {
      std::cout << "Reco muon ref not found..." << "  Event ID = " << event.id() << endl;
    }
  }

  // Evaluate cuts here with string object selector, and any code that
  // cannot be done in the string object selector (none so far).
  int cutFor = muon_selector(*new_mu) ? 0 : 1;

  return std::make_pair(new_mu, cutFor);
}

template <typename T>
edm::OrphanHandle<std::vector<T> > CustoTnPLeptonProducer::doLeptons(edm::Event& event, const edm::InputTag& src, const edm::InputTag& view_src, const std::string& instance_label) {
  typedef std::vector<T> TCollection;
  edm::Handle<TCollection> leptons; 
  event.getByLabel(src, leptons); 

  static std::map<std::string, bool> warned;
  if (!leptons.isValid()) {
    if (!warned[instance_label]) {
      edm::LogWarning("LeptonsNotFound") << src << " for " << instance_label << " not found, not producing anything -- not warning any more either.";
      warned[instance_label] = true;
    }
    return edm::OrphanHandle<std::vector<T> >();
  }

  edm::Handle<reco::CandidateView> lepton_view;
  event.getByLabel(view_src, lepton_view);

  std::unique_ptr<TCollection> new_leptons(new TCollection);

  for (size_t i = 0; i < leptons->size(); ++i) {
    std::pair<T*,int> res = doLepton(event, leptons->at(i), lepton_view->refAt(i));
    if (res.first == 0)
      continue;
    res.first->addUserInt("cutFor", res.second);
    new_leptons->push_back(*res.first);
    delete res.first;
  }

  return event.put(std::move(new_leptons), instance_label);
}

void CustoTnPLeptonProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  // Grab the match map between PAT photons and PAT muons so we can
  // embed the photon candidates later.
  //std::cout << event.id() << std::endl;
  event.getByLabel(muon_photon_match_src, muon_photon_match_map);
  static bool warned = false;
  if (!warned && !muon_photon_match_map.isValid()) {
    edm::LogWarning("PhotonsNotFound") << muon_photon_match_src << " for photons not found, not trying to embed their matches to muons -- not warning any more either.";
    warned = true;
  }

  // Store the L3 muons for trigger match embedding, and clear the
  // vector of matched flags that indicate whether the particular L3
  // muon has been used in a match already. This means matching
  // ambiguities are resolved by original sort order of the
  // candidates; no attempt is done to find a global best
  // matching. (This is how it was done in our configuration of the
  // PATTrigger matcher previously, so why not.) We do this for both
  // the main path and the pres caled path.
  
  // CustoTnPTriggerPathsAndFilters pandf(event);
  // if (!pandf.valid)
  //   throw cms::Exception("CustoTnPLeptonProducer") << "could not determine the HLT path and filter names for this event\n";

  //-- for HLT match
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> trigger_summary_src;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

  event.getByToken(triggerBits_, triggerBits);
  event.getByToken(trigger_summary_src_, trigger_summary_src);
  event.getByToken(triggerPrescales_, triggerPrescales);

  const edm::TriggerNames &names = event.triggerNames(*triggerBits);

  vec_L3_muons.clear();
  for(unsigned i_f=0; i_f<trigger_filters.size(); ++i_f) {
    vec_L3_muons.push_back({});
    vec_L3_muons.back().clear();
  }

  for(pat::TriggerObjectStandAlone obj : *trigger_summary_src) {
    obj.unpackPathNames(names);
    obj.unpackFilterLabels(event, *triggerBits); // for 2017~
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h) {

      for(unsigned i_f=0; i_f<trigger_filters.size(); ++i_f) {
        if (obj.filterLabels()[h] == trigger_filters[i_f]){ 
          vec_L3_muons[i_f].push_back(obj);
        }
      }

    }
  }

  vec_L3_muons_matched.clear();
  for(unsigned i_f=0; i_f<trigger_filters.size(); ++i_f) {
    vec_L3_muons_matched.push_back({});
    vec_L3_muons_matched.back().clear();
    vec_L3_muons_matched.back().resize(vec_L3_muons[i_f].size(), 0);
  }

  // Using the main choice for momentum assignment, make the primary
  // collection of muons, which will have branch name
  // e.g. leptons:muons.
  muon_track_for_momentum = muon_track_for_momentum_primary;
  edm::OrphanHandle<pat::MuonCollection> muons = doLeptons<pat::Muon>(event, muon_src, muon_view_src, "muons");

  // Now make secondary collections of muons using the momentum
  // assignments specified. They will come out as e.g. leptons:tpfms,
  // leptons:picky, ...
  for (size_t i = 0; i < muon_tracks_for_momentum.size(); ++i) {

    vec_L3_muons_matched.clear();
    for(unsigned i_f=0; i_f<trigger_filters.size(); ++i_f) {
      vec_L3_muons_matched.push_back({});
      vec_L3_muons_matched.back().clear();
      vec_L3_muons_matched.back().resize(vec_L3_muons[i_f].size(), 0);
    }

    muon_track_for_momentum = muon_tracks_for_momentum[i];
    doLeptons<pat::Muon>(event, muon_src, muon_view_src, muon_track_for_momentum);
  }

}

DEFINE_FWK_MODULE(CustoTnPLeptonProducer);
