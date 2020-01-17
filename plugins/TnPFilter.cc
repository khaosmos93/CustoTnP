// #include "TH1F.h"
// #include "TH2F.h"
// #include "TProfile.h"
// #include "TTree.h"
// #include "TMath.h"
// #include "TRandom3.h"
// #include "TString.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CustoTnP/Analyzer/src/DileptonUtilities.h"
#include "CustoTnP/Analyzer/src/GeneralUtilities.h"
#include "CustoTnP/Analyzer/src/PATUtilities.h"
#include "CustoTnP/Analyzer/src/ToConcrete.h"
#include "CustoTnP/Analyzer/src/TrackUtilities.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//using namespace std;
void CopyVectorToArray( std::vector<double>& vec, double*& arr ) {
  int nBin = (int)vec.size()-1; // -- # bins = # bin edges-1 -- //
  arr = new Double_t[nBin+1]; // -- dynamic allocation -- //
  for(int i=0; i<nBin+1; i++)
    arr[i] = vec[i];
}

class TnPFilter : public edm::EDFilter {
 public:
  explicit TnPFilter(const edm::ParameterSet&);

 private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  void getBSandPV(edm::Event&);

  std::vector<int> calcNShowers(const reco::CandidateBaseRef&, bool);

  edm::InputTag dilepton_src;
  edm::InputTag beamspot_src;
  edm::InputTag vertex_src;
  const bool use_bs_and_pv;
  const reco::BeamSpot* beamspot;
  const reco::Vertex*   vertex;
  int                   nVtx;

  int shower_tag;    // 1: Digis, 2: Segments
  std::vector<int>  threshold_b;
  std::vector<int>  threshold_e;

  StringCutObjectSelector<pat::Muon> passing_probe_selector;

  double minMass;
  double maxMass;
  double min_nVtx;
  double max_nVtx;
  double probe_pt_min;
  double passing_probe_dpt_over_pt_max;
  double passing_probe_dz_max;

  bool isAOD;
  const bool ShutUp;
};

TnPFilter::TnPFilter(const edm::ParameterSet& cfg)
  : dilepton_src(cfg.getParameter<edm::InputTag>("dilepton_src")),
    beamspot_src(cfg.getParameter<edm::InputTag>("beamspot_src")),
    vertex_src(cfg.getParameter<edm::InputTag>("vertex_src")),
    use_bs_and_pv(cfg.getParameter<bool>("use_bs_and_pv")),
    beamspot(0),
    vertex(0),
    nVtx(0),

    shower_tag(cfg.getParameter<int>("shower_tag")),
    threshold_b(cfg.getParameter<std::vector<int>>("threshold_b")),
    threshold_e(cfg.getParameter<std::vector<int>>("threshold_e")),

    passing_probe_selector(cfg.getParameter<std::string>("passing_probe_cut")),

    minMass(cfg.getParameter<double>("minMass")),
    maxMass(cfg.getParameter<double>("maxMass")),

    min_nVtx(cfg.getParameter<double>("min_nVtx")),
    max_nVtx(cfg.getParameter<double>("max_nVtx")),

    probe_pt_min(cfg.getParameter<double>("probe_pt_min")),
    passing_probe_dpt_over_pt_max(cfg.getParameter<double>("passing_probe_dpt_over_pt_max")),
    passing_probe_dz_max(cfg.getParameter<double>("passing_probe_dz_max")),

    isAOD(cfg.getParameter<bool>("isAOD")),

    ShutUp(cfg.getParameter<bool>("ShutUp"))
{
  consumes<pat::CompositeCandidateCollection>(dilepton_src);
  consumes<reco::BeamSpot>(beamspot_src);
  consumes<reco::VertexCollection>(vertex_src);
  mayConsume<GenEventInfoProduct>(edm::InputTag("generator"));
}

void TnPFilter::getBSandPV(edm::Event& event) {
  // We store these as bare pointers. Should find better way, but
  // don't want to pass them around everywhere...
  edm::Handle<reco::BeamSpot> hbs;
  event.getByLabel(beamspot_src, hbs);
  beamspot = hbs.isValid() ? &*hbs : 0; // nice and fragile

  edm::Handle<reco::VertexCollection> vertices;
  event.getByLabel(vertex_src, vertices);
  vertex = 0;
  int vertex_count = 0;
  for (reco::VertexCollection::const_iterator it = vertices->begin(), ite = vertices->end(); it != ite; ++it) {
    if (it->ndof() > 4 && fabs(it->z()) <= 24 && fabs(it->position().rho()) <= 2) {
      if (vertex == 0)
        vertex = &*it;
      ++vertex_count;
    }
  }
  nVtx = vertex_count;
}

std::vector<int> TnPFilter::calcNShowers(
                                      const reco::CandidateBaseRef& mu,
                                      bool verbos = false )
{

  int etaCat = -1;
  if( fabs(mu->eta()) < 0.9 )
    etaCat = 1;
  else if( fabs(mu->eta()) >= 1.2 )
    etaCat = 2;

  if( etaCat<0 )
    return {-999};

  const pat::Muon* muPat = toConcretePtr<pat::Muon>(mu);
  int st1, st2, st3, st4;
  if(shower_tag == 1) {  // Digis
    st1 = (etaCat==1) ? muPat->userInt("nDigisDT1") : muPat->userInt("nDigisCSC1");
    st2 = (etaCat==1) ? muPat->userInt("nDigisDT2") : muPat->userInt("nDigisCSC2");
    st3 = (etaCat==1) ? muPat->userInt("nDigisDT3") : muPat->userInt("nDigisCSC3");
    st4 = (etaCat==1) ? muPat->userInt("nDigisDT4") : muPat->userInt("nDigisCSC4");
  }
  else if( shower_tag == 2) {  // Segments
    st1 = (etaCat==1) ? muPat->userInt("nSegsDT1") : muPat->userInt("nSegsCSC1");
    st2 = (etaCat==1) ? muPat->userInt("nSegsDT2") : muPat->userInt("nSegsCSC2");
    st3 = (etaCat==1) ? muPat->userInt("nSegsDT3") : muPat->userInt("nSegsCSC3");
    st4 = (etaCat==1) ? muPat->userInt("nSegsDT4") : muPat->userInt("nSegsCSC4");
  }
  else {
    std::cout << "WARNING: TnPFilter::calcNShowers: wrong shower_tag" << std::endl;
    return {-999};
  }

  int threshold_st1 = -999;
  int threshold_st2 = -999;
  int threshold_st3 = -999;
  int threshold_st4 = -999;

  if( etaCat==1 ) {
    threshold_st1 = threshold_b[0];
    threshold_st2 = threshold_b[1];
    threshold_st3 = threshold_b[2];
    threshold_st4 = threshold_b[3];
  }
  else if( etaCat==2 ) {
    threshold_st1 = threshold_e[0];
    threshold_st2 = threshold_e[1];
    threshold_st3 = threshold_e[2];
    threshold_st4 = threshold_e[3];
  }

  bool is_st1 = (st1 >= threshold_st1);
  bool is_st2 = (st2 >= threshold_st2);
  bool is_st3 = (st3 >= threshold_st3);
  bool is_st4 = (st4 >= threshold_st4);

  std::vector<int> out_vec = {};

  out_vec.push_back( (int)( !is_st1 && !is_st2 && !is_st3 && !is_st4 ) );

  out_vec.push_back( (int)(  is_st1 && !is_st2 && !is_st3 && !is_st4 ) );
  out_vec.push_back( (int)( !is_st1 &&  is_st2 && !is_st3 && !is_st4 ) );
  out_vec.push_back( (int)( !is_st1 && !is_st2 &&  is_st3 && !is_st4 ) );
  out_vec.push_back( (int)( !is_st1 && !is_st2 && !is_st3 &&  is_st4 ) );

  out_vec.push_back( (int)(  is_st1 &&  is_st2 && !is_st3 && !is_st4 ) );
  out_vec.push_back( (int)(  is_st1 && !is_st2 &&  is_st3 && !is_st4 ) );
  out_vec.push_back( (int)(  is_st1 && !is_st2 && !is_st3 &&  is_st4 ) );
  out_vec.push_back( (int)( !is_st1 &&  is_st2 &&  is_st3 && !is_st4 ) );
  out_vec.push_back( (int)( !is_st1 &&  is_st2 && !is_st3 &&  is_st4 ) );
  out_vec.push_back( (int)( !is_st1 && !is_st2 &&  is_st3 &&  is_st4 ) );

  out_vec.push_back( (int)(  is_st1 &&  is_st2 &&  is_st3 && !is_st4 ) );
  out_vec.push_back( (int)(  is_st1 &&  is_st2 && !is_st3 &&  is_st4 ) );
  out_vec.push_back( (int)(  is_st1 && !is_st2 &&  is_st3 &&  is_st4 ) );
  out_vec.push_back( (int)( !is_st1 &&  is_st2 &&  is_st3 &&  is_st4 ) );

  out_vec.push_back( (int)(  is_st1 &&  is_st2 &&  is_st3 &&  is_st4 ) );

  if(verbos) {
    std::cout << "\neta = " << mu->eta() << "  => etaCat: " << etaCat << std::endl;
    std::cout << "\tst1= " << st1 << std::endl;
    std::cout << "\tst2= " << st2 << std::endl;
    std::cout << "\tst3= " << st3 << std::endl;
    std::cout << "\tst4= " << st4 << std::endl;
    std::cout << "\t--> nShowers= ";
    for(int is=0; is<16; ++is) {
      std::cout << out_vec[is] << ", ";
    }
    std::cout << std::endl;
  }

  return out_vec;
}

bool TnPFilter::filter(edm::Event& event, const edm::EventSetup& setup) {

  bool pass_tnp = false;

  if (use_bs_and_pv)
    getBSandPV(event);

  if( !(min_nVtx <= nVtx && max_nVtx >= nVtx ) )
    return false;

  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  event.getByLabel(dilepton_src, dileptons);

  if( !dileptons.isValid() ) {
    std::cout << "TnPFilter::filter : !dileptons.isValid() ---> return" << std::endl;
    return false;
  }

  pat::CompositeCandidateCollection::const_iterator dil = dileptons->begin(), dile = dileptons->end();
  for ( ; dil != dile; ++dil) {

    if( !( dil->hasUserInt("isTag0Probe1") && dil->hasUserInt("isTag1Probe0") ) ) {
      std::cout << "TnPFilter::filter : no isTag0Probe1 or isTag1Probe0 in dil" << std::endl;
      continue;
    }

    // double dil_mass = dil->mass();
    // double dil_vertex_mass = dil->userFloat("vertexM");
    float lep0_dpt_over_pt = dil->userFloat("lep0_dpt_over_pt");
    float lep1_dpt_over_pt = dil->userFloat("lep1_dpt_over_pt");

    const reco::CandidateBaseRef& lep0 = dileptonDaughter(*dil, 0);
    const reco::CandidateBaseRef& lep1 = dileptonDaughter(*dil, 1);

    if(lep0.isNonnull() && lep1.isNonnull()) {
      const pat::Muon* mu0 = toConcretePtr<pat::Muon>(lep0);
      const pat::Muon* mu1 = toConcretePtr<pat::Muon>(lep1);
      if(mu0 && mu1) {

        //---- Tag0 and Probe1
        if( dil->userInt("isTag0Probe1") ) {
          // const reco::CandidateBaseRef& TagMu = lep0;
          const reco::CandidateBaseRef& ProbeMu = lep1;

          // int probe_nshowers = -999;
          // if(isAOD) {
          //   std::vector<int> vec_showers = calcNShowers(ProbeMu);

          //   if(vec_showers.size()==16) {
          //     for(int is=0; is<(int)vec_showers.size(); ++is) {
          //       if(vec_showers[is]) {
          //         probe_nshowers = is;
          //         break;
          //       }
          //     }
          //   }
          // }

          if(!ShutUp)  std::cout << "TnPFilter::filter : Tag0 and Probe1" << std::endl;
          if(!ShutUp)  std::cout << "                                              pT=" << ProbeMu->pt() << std::endl;
          if(!ShutUp)  std::cout << "                                             eta=" << ProbeMu->eta() << std::endl;
          if(!ShutUp)  std::cout << "                                             phi=" << ProbeMu->phi() << std::endl;
          if(!ShutUp)  std::cout << "                                            nVtx=" << nVtx << std::endl;

          bool isPassingProbe = ( passing_probe_selector(*mu1) && (ProbeMu->pt() > probe_pt_min) && (lep1_dpt_over_pt < passing_probe_dpt_over_pt_max) && (fabs(mu1->innerTrack()->dz( vertex->position() )) < passing_probe_dz_max) );
          if(!ShutUp)  std::cout << "                                  isPassingProbe = " << isPassingProbe << std::endl;

          pass_tnp = (pass_tnp || isPassingProbe);
        } // Tag0 and Probe1

        //---- Tag1 and Probe0
        if( dil->userInt("isTag1Probe0") ) {
          // const reco::CandidateBaseRef& TagMu = lep1;
          const reco::CandidateBaseRef& ProbeMu = lep0;

          // int probe_nshowers = -999;
          // if(isAOD) {
          //   std::vector<int> vec_showers = calcNShowers(ProbeMu);

          //   if(vec_showers.size()==16) {
          //     for(int is=0; is<(int)vec_showers.size(); ++is) {
          //       if(vec_showers[is]) {
          //         probe_nshowers = is;
          //         break;
          //       }
          //     }
          //   }
          // }

          if(!ShutUp)  std::cout << "TnPFilter::filter : Tag1 and Probe0" << std::endl;
          if(!ShutUp)  std::cout << "                                              pT=" << ProbeMu->pt() << std::endl;
          if(!ShutUp)  std::cout << "                                             eta=" << ProbeMu->eta() << std::endl;
          if(!ShutUp)  std::cout << "                                             phi=" << ProbeMu->phi() << std::endl;
          if(!ShutUp)  std::cout << "                                            nVtx=" << nVtx << std::endl;

          bool isPassingProbe = ( passing_probe_selector(*mu0) && (ProbeMu->pt() > probe_pt_min) && (lep0_dpt_over_pt < passing_probe_dpt_over_pt_max) && (fabs(mu0->innerTrack()->dz( vertex->position() )) < passing_probe_dz_max) );
          if(!ShutUp)  std::cout << "                                  isPassingProbe = " << isPassingProbe << std::endl;

          pass_tnp = (pass_tnp || isPassingProbe);
        } // Tag1 and Probe0

        if(!ShutUp)  std::cout << std::endl;
      } // if(mu0 && mu1)
    } // if(lep0.isNonnull() && lep1.isNonnull())

  } // for ( ; dil != dile; ++dil)

  return pass_tnp;
}

DEFINE_FWK_MODULE(TnPFilter);
