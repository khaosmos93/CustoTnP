#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"

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
#include "FWCore/Framework/interface/EDAnalyzer.h"
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

class CustoTnPHistosForTnP : public edm::EDAnalyzer {
 public:
  explicit CustoTnPHistosForTnP(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  void getBSandPV(const edm::Event&);
  void fillTnPControlHistos(const pat::CompositeCandidate&, const reco::CandidateBaseRef&, const reco::CandidateBaseRef&, float, int, bool );

  typedef std::pair< std::vector<TH1F*>, std::vector<TH1F*> > BinHistos;
  void fillTnPBinHistos( double, double, bool, double, std::vector<double>& , BinHistos&, bool );

  typedef std::vector< TH2F* > BinHistos2D;
  void fillTnPBinHistos2D( double, double, bool, double, BinHistos2D&, bool );

  std::vector<int> calcNShowers(const reco::CandidateBaseRef&, int, int, int, int, int, int, int, int, bool);

  edm::InputTag dilepton_src;
  edm::InputTag beamspot_src;
  edm::InputTag vertex_src;
  const bool use_bs_and_pv;
  const reco::BeamSpot* beamspot;
  const reco::Vertex*   vertex;
  int                   nVtx;

  StringCutObjectSelector<pat::Muon> passing_probe_selector;
  StringCutObjectSelector<pat::Muon> comparison_probe_selector;

  double minMass;
  double maxMass;
  double probe_pt_min;
  double passing_probe_dpt_over_pt_max;
  double passing_probe_dz_max;
  double comparison_probe_dpt_over_pt_max;
  double comparison_probe_dz_max;

  bool   _usePrescaleWeight;
  int    _prescaleWeight;
  double _eventWeight;
  bool   _useMadgraphWeight;
  double _madgraphWeight;

  edm::InputTag pileup_src;
  double _pileupWeight;
  std::vector<double> vec_PileUpWeight;

  double _totalWeight;

  bool isAOD;
  bool useBinHistos2D;
  const bool ShutUp;


  // -- Histograms -- //
  TH1F* NBeamSpot;
  TH1F* NVertices;
  TH1F* NDileptons;
  TH1F* WeightMadGraph;

  TH1F *TagPt;
  TH1F *TagEta;
  TH1F *TagPhi;

  typedef std::vector< TH1F* > ProbeHistos;
  typedef std::vector< TH2F* > ProbeHistos2D;


  std::vector< TH1F* > make_probe_histos(TString name, std::vector<double>& vec_x_bins) {
    edm::Service<TFileService> fs;

    std::vector< TH1F* > vec_h = {};

    int n_x_bins = ( (int)vec_x_bins.size() ) - 1;
    double *arr_x_bins;
    CopyVectorToArray(vec_x_bins, arr_x_bins);

    TH1F* Probe = fs->make<TH1F>("Probe"+name,     "", n_x_bins, arr_x_bins );  Probe->Sumw2();
    TH1F* Pass  = fs->make<TH1F>("ProbePass"+name, "", n_x_bins, arr_x_bins );  Pass->Sumw2();
    TH1F* Fail  = fs->make<TH1F>("ProbeFail"+name, "", n_x_bins, arr_x_bins );  Fail->Sumw2();

    vec_h.push_back( Probe );
    vec_h.push_back( Pass );
    vec_h.push_back( Fail );

    return vec_h;
  }

  std::vector< TH1F* > make_probe_histos(TString name, int n_x_bins, double x_min, double x_max) {
    edm::Service<TFileService> fs;

    std::vector< TH1F* > vec_h = {};

    TH1F* Probe = fs->make<TH1F>("Probe"+name,     "", n_x_bins, x_min, x_max );  Probe->Sumw2();
    TH1F* Pass  = fs->make<TH1F>("ProbePass"+name, "", n_x_bins, x_min, x_max );  Pass->Sumw2();
    TH1F* Fail  = fs->make<TH1F>("ProbeFail"+name, "", n_x_bins, x_min, x_max );  Fail->Sumw2();

    vec_h.push_back( Probe );
    vec_h.push_back( Pass );
    vec_h.push_back( Fail );

    return vec_h;
  }

  std::vector< TH2F* > make_probe_histos_2D(TString name, std::vector<double>& vec_x_bins, int n_y_bins, double y_min, double y_max) {
    edm::Service<TFileService> fs;

    std::vector< TH2F* > vec_h = {};

    int n_x_bins = ( (int)vec_x_bins.size() ) - 1;
    double *arr_x_bins;
    CopyVectorToArray(vec_x_bins, arr_x_bins);

    TH2F* Probe = fs->make<TH2F>("Probe"+name,     "", n_x_bins, arr_x_bins, n_y_bins, y_min, y_max );  Probe->Sumw2();
    TH2F* Pass  = fs->make<TH2F>("ProbePass"+name, "", n_x_bins, arr_x_bins, n_y_bins, y_min, y_max );  Pass->Sumw2();
    TH2F* Fail  = fs->make<TH2F>("ProbeFail"+name, "", n_x_bins, arr_x_bins, n_y_bins, y_min, y_max );  Fail->Sumw2();

    vec_h.push_back( Probe );
    vec_h.push_back( Pass );
    vec_h.push_back( Fail );

    return vec_h;
  }

  std::vector< TH2F* > make_probe_histos_2D(TString name, int n_x_bins, double x_min, double x_max, int n_y_bins, double y_min, double y_max) {
    edm::Service<TFileService> fs;

    std::vector< TH2F* > vec_h = {};

    TH2F* Probe = fs->make<TH2F>("Probe"+name,     "", n_x_bins, x_min, x_max, n_y_bins, y_min, y_max );  Probe->Sumw2();
    TH2F* Pass  = fs->make<TH2F>("ProbePass"+name, "", n_x_bins, x_min, x_max, n_y_bins, y_min, y_max );  Pass->Sumw2();
    TH2F* Fail  = fs->make<TH2F>("ProbeFail"+name, "", n_x_bins, x_min, x_max, n_y_bins, y_min, y_max );  Fail->Sumw2();

    vec_h.push_back( Probe );
    vec_h.push_back( Pass );
    vec_h.push_back( Fail );

    return vec_h;
  }

  ProbeHistos   ProbePt;
  ProbeHistos   ProbeEta;
  ProbeHistos   ProbePhi;
  ProbeHistos   ProbeNVertices;

  ProbeHistos2D ProbeEtaPhi;
  ProbeHistos2D ProbeEtaPt;
  ProbeHistos2D ProbeEtaDptPt;
  ProbeHistos2D ProbeEtaShower;
  ProbeHistos2D ProbePtShowerB;
  ProbeHistos2D ProbePtShowerE;
  ProbeHistos2D ProbePShowerB;
  ProbeHistos2D ProbePShowerE;

  ProbeHistos2D ProbeEtaHitsSt1;
  ProbeHistos2D ProbePtHitsSt1B;
  ProbeHistos2D ProbePtHitsSt1E;
  ProbeHistos2D ProbePHitsSt1B;
  ProbeHistos2D ProbePHitsSt1E;

  ProbeHistos2D ProbeEtaHitsSt2;
  ProbeHistos2D ProbePtHitsSt2B;
  ProbeHistos2D ProbePtHitsSt2E;
  ProbeHistos2D ProbePHitsSt2B;
  ProbeHistos2D ProbePHitsSt2E;

  ProbeHistos2D ProbeEtaHitsSt3;
  ProbeHistos2D ProbePtHitsSt3B;
  ProbeHistos2D ProbePtHitsSt3E;
  ProbeHistos2D ProbePHitsSt3B;
  ProbeHistos2D ProbePHitsSt3E;

  ProbeHistos2D ProbeEtaHitsSt4;
  ProbeHistos2D ProbePtHitsSt4B;
  ProbeHistos2D ProbePtHitsSt4E;
  ProbeHistos2D ProbePHitsSt4B;
  ProbeHistos2D ProbePHitsSt4E;

  ProbeHistos   PairNoPtMass;
  ProbeHistos   PairMass;
  ProbeHistos   PairPt;
  ProbeHistos   PairEta;
  ProbeHistos   PairRap;

  // -- bin distributions for passing and failing probes -- //
  std::vector<double> vec_PtBins;
  std::vector<double> vec_AbsPBins;
  std::vector<double> vec_EtaBins;
  std::vector<double> vec_PhiBins;
  std::vector<double> vec_VtxBins;
  std::vector<double> vec_ShowerBins;

  TH1F *hEffTemplatePt;
  TH1F *hEffTemplateAbsP;
  TH1F *hEffTemplateEta;
  TH1F *hEffTemplatePhi;
  TH1F *hEffTemplateVtx;
  TH1F *hEffTemplateShower;

  std::pair< std::vector<TH1F*>, std::vector<TH1F*> > make_bin_histos(TString name, std::vector<double>& vec_bins) {
    edm::Service<TFileService> fs;

    int nbins = ( (int)vec_bins.size() ) - 1;

    unsigned nMassBin = (unsigned)(maxMass - minMass);

    std::vector<TH1F*> vec_pass;
    std::vector<TH1F*> vec_fail;

    TString histNameBase = TString::Format("h%s", name.Data());

    for(int i=0; i<nbins; ++i) {
      TString binInfo = TString::Format("%02dbin", i);

      TString histNamePass = histNameBase + "Pass_" + binInfo;
      TH1F* hTempPass = fs->make<TH1F>(histNamePass, "", nMassBin, minMass, maxMass);  hTempPass->Sumw2();
      vec_pass.push_back( hTempPass );

      TString histNameFail = histNameBase + "Fail_" + binInfo;
      TH1F* hTempFail = fs->make<TH1F>(histNameFail, "", nMassBin, minMass, maxMass);  hTempFail->Sumw2();
      vec_fail.push_back( hTempFail );
    }

    return std::make_pair(vec_pass, vec_fail);
  }

  std::vector< TH2F* > make_bin_histos_2D(TString name, std::vector<double>& vec_bins) {
    edm::Service<TFileService> fs;

    std::vector< TH2F* > vec_h2d = {};

    Double_t mass_bins_for_2D[13] = {0, 60, 120, 200, 300, 400, 800, 1400, 2300, 3500, 4500, 6000, 1e6};
    unsigned nMassBin = (unsigned)(12);

    int nbins = ( (int)vec_bins.size() ) - 1;
    double *arr_bins;
    CopyVectorToArray(vec_bins, arr_bins);

    TString histNameBase = TString::Format("h2d%s", name.Data());

    TH2F* h2d_pass = fs->make<TH2F>(histNameBase + "Pass", "", nMassBin, mass_bins_for_2D, nbins, arr_bins );  h2d_pass->Sumw2();
    TH2F* h2d_fail = fs->make<TH2F>(histNameBase + "Fail", "", nMassBin, mass_bins_for_2D, nbins, arr_bins );  h2d_fail->Sumw2();
    TH2F* h2d_sum  = fs->make<TH2F>(histNameBase + "Sum",  "", nMassBin, mass_bins_for_2D, nbins, arr_bins );  h2d_sum->Sumw2();
    TH2F* h2d_sqs  = fs->make<TH2F>(histNameBase + "Sqs",  "", nMassBin, mass_bins_for_2D, nbins, arr_bins );  h2d_sqs->Sumw2();

    vec_h2d.push_back( h2d_pass );
    vec_h2d.push_back( h2d_fail );
    vec_h2d.push_back( h2d_sum );
    vec_h2d.push_back( h2d_sqs );

    return vec_h2d;
  }

  // -- Bin histograms -- //
  BinHistos Pt;
  BinHistos PtB;  // 0.0 < |eta| < 0.9
  BinHistos PtO;  // 0.9 < |eta| < 1.2
  BinHistos PtE;  // 1.2 < |eta| < 2.1
  BinHistos PtF;  // 2.1 < |eta| < 2.4

  BinHistos AbsP;
  BinHistos AbsPB;  // 0.0 < |eta| < 0.9
  BinHistos AbsPO;  // 0.9 < |eta| < 1.2
  BinHistos AbsPE;  // 1.2 < |eta| < 2.1
  BinHistos AbsPF;  // 2.1 < |eta| < 2.4

  BinHistos Eta;
  BinHistos Phi;
  BinHistos Vtx;

  BinHistos ShowerB;
  BinHistos ShowerBPt200;
  BinHistos ShowerBPt400;
  BinHistos ShowerBP200;
  BinHistos ShowerBP400;
  BinHistos ShowerE;
  BinHistos ShowerEPt200;
  BinHistos ShowerEPt400;
  BinHistos ShowerEP200;
  BinHistos ShowerEP400;

  BinHistos2D Pt2D;
  BinHistos2D PtB2D;  // 0.0 < |eta| < 0.9
  BinHistos2D PtO2D;  // 0.9 < |eta| < 1.2
  BinHistos2D PtE2D;  // 1.2 < |eta| < 2.1
  BinHistos2D PtF2D;  // 2.1 < |eta| < 2.4

  BinHistos2D AbsP2D;
  BinHistos2D AbsPB2D;  // 0.0 < |eta| < 0.9
  BinHistos2D AbsPO2D;  // 0.9 < |eta| < 1.2
  BinHistos2D AbsPE2D;  // 1.2 < |eta| < 2.1
  BinHistos2D AbsPF2D;  // 2.1 < |eta| < 2.4

  BinHistos2D Eta2D;
  BinHistos2D Phi2D;
  BinHistos2D Vtx2D;

  BinHistos2D ShowerB2D;
  BinHistos2D ShowerBPt2002D;
  BinHistos2D ShowerBPt4002D;
  BinHistos2D ShowerBP2002D;
  BinHistos2D ShowerBP4002D;
  BinHistos2D ShowerE2D;
  BinHistos2D ShowerEPt2002D;
  BinHistos2D ShowerEPt4002D;
  BinHistos2D ShowerEP2002D;
  BinHistos2D ShowerEP4002D;

  bool IsRealData;
  int RunNum;
  int LumiBlockNum;
  unsigned long long EventNum;
  double Mass;
  double VertexMass;
  double Probe_Pt;
  double Probe_Eta;
  TTree* comparison_tree;
};

CustoTnPHistosForTnP::CustoTnPHistosForTnP(const edm::ParameterSet& cfg)
  : dilepton_src(cfg.getParameter<edm::InputTag>("dilepton_src")),
    beamspot_src(cfg.getParameter<edm::InputTag>("beamspot_src")),
    vertex_src(cfg.getParameter<edm::InputTag>("vertex_src")),
    use_bs_and_pv(cfg.getParameter<bool>("use_bs_and_pv")),
    beamspot(0),
    vertex(0),
    nVtx(0),

    passing_probe_selector(cfg.getParameter<std::string>("passing_probe_cut")),
    comparison_probe_selector(cfg.getParameter<std::string>("comparison_probe_cut")),

    minMass(cfg.getParameter<double>("minMass")),
    maxMass(cfg.getParameter<double>("maxMass")),

    probe_pt_min(cfg.getParameter<double>("probe_pt_min")),
    passing_probe_dpt_over_pt_max(cfg.getParameter<double>("passing_probe_dpt_over_pt_max")),
    passing_probe_dz_max(cfg.getParameter<double>("passing_probe_dz_max")),
    comparison_probe_dpt_over_pt_max(cfg.getParameter<double>("comparison_probe_dpt_over_pt_max")),
    comparison_probe_dz_max(cfg.getParameter<double>("comparison_probe_dz_max")),

    _usePrescaleWeight(cfg.getUntrackedParameter<bool>("usePrescaleWeight",false)),
    _prescaleWeight(1),
    _useMadgraphWeight(cfg.getParameter<bool>("useMadgraphWeight")),
    _madgraphWeight(1.),

    pileup_src(cfg.getParameter<edm::InputTag>("pileup_src")),
    _pileupWeight(1.),
    vec_PileUpWeight(cfg.getParameter<std::vector<double>>("vec_PileUpWeight")),

    _totalWeight(1.),

    isAOD(cfg.getParameter<bool>("isAOD")),
    useBinHistos2D(cfg.getParameter<bool>("useBinHistos2D")),

    ShutUp(cfg.getParameter<bool>("ShutUp")),

    vec_PtBins(cfg.getParameter<std::vector<double>>("vec_PtBins")),
    vec_AbsPBins(cfg.getParameter<std::vector<double>>("vec_AbsPBins")),
    vec_EtaBins(cfg.getParameter<std::vector<double>>("vec_EtaBins")),
    vec_PhiBins(cfg.getParameter<std::vector<double>>("vec_PhiBins")),
    vec_VtxBins(cfg.getParameter<std::vector<double>>("vec_VtxBins")),
    vec_ShowerBins(cfg.getParameter<std::vector<double>>("vec_ShowerBins"))
{
  consumes<pat::CompositeCandidateCollection>(dilepton_src);
  consumes<reco::BeamSpot>(beamspot_src);
  consumes<reco::VertexCollection>(vertex_src);
  consumes<std::vector<PileupSummaryInfo> >(pileup_src);
  mayConsume<GenEventInfoProduct>(edm::InputTag("generator"));

  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  NBeamSpot = fs->make<TH1F>("NBeamSpot", "# beamspots/event",  2, 0,  2);
  NVertices = fs->make<TH1F>("NVertices", "# vertices/event",  200, 0, 200);

  // Dilepton multiplicity.
  NDileptons = fs->make<TH1F>("NDileptons", "# dileptons/event", 10, 0, 10);

  // Generator weights
  WeightMadGraph = fs->make<TH1F>("WeightMadGraph", "weight per event", 4, -2,2);

  // Tag
  TagPt = fs->make<TH1F>("TagPt", "Tag pT", 10000, 0, 10000);
  TagEta = fs->make<TH1F>("TagEta", "Tag #eta",    96, -4.8, 4.8);
  TagPhi = fs->make<TH1F>("TagPhi", "Tag #phi", 41, -TMath::Pi(), TMath::Pi());

  std::vector<double> eta_bins_for_2D = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.0, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4};

  // Probe
  ProbePt        = make_probe_histos("Pt",              10000, 0, 10000);
  ProbeEta       = make_probe_histos("Eta",             96, -4.8, 4.8);
  ProbePhi       = make_probe_histos("Phi",             41, -TMath::Pi(), TMath::Pi());
  ProbeNVertices = make_probe_histos("NVertices",       200, 0, 200);
  ProbeEtaPhi    = make_probe_histos_2D("EtaPhi",       eta_bins_for_2D, 41, -TMath::Pi(), TMath::Pi());
  ProbeEtaPt     = make_probe_histos_2D("EtaPt",        eta_bins_for_2D, 10000, 0, 10000);
  ProbeEtaDptPt  = make_probe_histos_2D("EtaDptPt",     eta_bins_for_2D, 1000, 0, 2);
  if(isAOD) {
  ProbeEtaShower = make_probe_histos_2D("EtaShower",    eta_bins_for_2D, 18, -1.5, 16.5);
  ProbePtShowerB = make_probe_histos_2D("PtShowerB",    100, 0, 10000, 18, -1.5, 16.5);
  ProbePtShowerE = make_probe_histos_2D("PtShowerE",    100, 0, 10000, 18, -1.5, 16.5);
  ProbePShowerB  = make_probe_histos_2D("PShowerB",     100, 0, 10000, 18, -1.5, 16.5);
  ProbePShowerE  = make_probe_histos_2D("PShowerE",     100, 0, 10000, 18, -1.5, 16.5);

  ProbeEtaHitsSt1 = make_probe_histos_2D("EtaHitsSt1",    eta_bins_for_2D, 100, 0, 200);
  ProbePtHitsSt1B = make_probe_histos_2D("PtHitsSt1B",    100, 0, 10000, 100, 0, 200);
  ProbePtHitsSt1E = make_probe_histos_2D("PtHitsSt1E",    100, 0, 10000, 100, 0, 200);
  ProbePHitsSt1B  = make_probe_histos_2D("PHitsSt1B",     100, 0, 10000, 100, 0, 200);
  ProbePHitsSt1E  = make_probe_histos_2D("PHitsSt1E",     100, 0, 10000, 100, 0, 200);

  ProbeEtaHitsSt2 = make_probe_histos_2D("EtaHitsSt2",    eta_bins_for_2D, 100, 0, 200);
  ProbePtHitsSt2B = make_probe_histos_2D("PtHitsSt2B",    100, 0, 10000, 100, 0, 200);
  ProbePtHitsSt2E = make_probe_histos_2D("PtHitsSt2E",    100, 0, 10000, 100, 0, 200);
  ProbePHitsSt2B  = make_probe_histos_2D("PHitsSt2B",     100, 0, 10000, 100, 0, 200);
  ProbePHitsSt2E  = make_probe_histos_2D("PHitsSt2E",     100, 0, 10000, 100, 0, 200);

  ProbeEtaHitsSt3 = make_probe_histos_2D("EtaHitsSt3",    eta_bins_for_2D, 100, 0, 200);
  ProbePtHitsSt3B = make_probe_histos_2D("PtHitsSt3B",    100, 0, 10000, 100, 0, 200);
  ProbePtHitsSt3E = make_probe_histos_2D("PtHitsSt3E",    100, 0, 10000, 100, 0, 200);
  ProbePHitsSt3B  = make_probe_histos_2D("PHitsSt3B",     100, 0, 10000, 100, 0, 200);
  ProbePHitsSt3E  = make_probe_histos_2D("PHitsSt3E",     100, 0, 10000, 100, 0, 200);

  ProbeEtaHitsSt4 = make_probe_histos_2D("EtaHitsSt4",    eta_bins_for_2D, 100, 0, 200);
  ProbePtHitsSt4B = make_probe_histos_2D("PtHitsSt4B",    100, 0, 10000, 100, 0, 200);
  ProbePtHitsSt4E = make_probe_histos_2D("PtHitsSt4E",    100, 0, 10000, 100, 0, 200);
  ProbePHitsSt4B  = make_probe_histos_2D("PHitsSt4B",     100, 0, 10000, 100, 0, 200);
  ProbePHitsSt4E  = make_probe_histos_2D("PHitsSt4E",     100, 0, 10000, 100, 0, 200);
  }

  // TnP pair
  PairNoPtMass   = make_probe_histos("PairNoPtMass",    10000, 0, 10000);
  PairMass       = make_probe_histos("PairMass",        10000, 0, 10000);
  PairPt         = make_probe_histos("PairPt",          10000, 0, 10000);
  PairEta        = make_probe_histos("PairEta",         96, -4.8, 4.8);
  PairRap        = make_probe_histos("PairRap",         96, -4.8, 4.8);

  // Bin histograms
  if(useBinHistos2D) {
    Pt2D      = make_bin_histos_2D("Pt",  vec_PtBins);
    PtB2D     = make_bin_histos_2D("PtB", vec_PtBins);
    PtO2D     = make_bin_histos_2D("PtO", vec_PtBins);
    PtE2D     = make_bin_histos_2D("PtE", vec_PtBins);
    PtF2D     = make_bin_histos_2D("PtF", vec_PtBins);

    AbsP2D    = make_bin_histos_2D("AbsP",  vec_AbsPBins);
    AbsPB2D   = make_bin_histos_2D("AbsPB", vec_AbsPBins);
    AbsPO2D   = make_bin_histos_2D("AbsPO", vec_AbsPBins);
    AbsPE2D   = make_bin_histos_2D("AbsPE", vec_AbsPBins);
    AbsPF2D   = make_bin_histos_2D("AbsPF", vec_AbsPBins);

    Eta2D     = make_bin_histos_2D("Eta", vec_EtaBins);
    Phi2D     = make_bin_histos_2D("Phi", vec_PhiBins);
    Vtx2D     = make_bin_histos_2D("Vtx", vec_VtxBins);

    if(isAOD) {
      ShowerB2D      = make_bin_histos_2D("ShowerB", vec_ShowerBins);
      ShowerBPt2002D = make_bin_histos_2D("ShowerBPt200", vec_ShowerBins);
      ShowerBPt4002D = make_bin_histos_2D("ShowerBPt400", vec_ShowerBins);
      ShowerBP2002D  = make_bin_histos_2D("ShowerBP200", vec_ShowerBins);
      ShowerBP4002D  = make_bin_histos_2D("ShowerBP400", vec_ShowerBins);
      ShowerE2D      = make_bin_histos_2D("ShowerE", vec_ShowerBins);
      ShowerEPt2002D = make_bin_histos_2D("ShowerEPt200", vec_ShowerBins);
      ShowerEPt4002D = make_bin_histos_2D("ShowerEPt400", vec_ShowerBins);
      ShowerEP2002D  = make_bin_histos_2D("ShowerEP200", vec_ShowerBins);
      ShowerEP4002D  = make_bin_histos_2D("ShowerEP400", vec_ShowerBins);
    }
  }

  else {
    Pt      = make_bin_histos("Pt",  vec_PtBins);
    PtB     = make_bin_histos("PtB", vec_PtBins);
    PtO     = make_bin_histos("PtO", vec_PtBins);
    PtE     = make_bin_histos("PtE", vec_PtBins);
    PtF     = make_bin_histos("PtF", vec_PtBins);

    AbsP    = make_bin_histos("AbsP",  vec_AbsPBins);
    AbsPB   = make_bin_histos("AbsPB", vec_AbsPBins);
    AbsPO   = make_bin_histos("AbsPO", vec_AbsPBins);
    AbsPE   = make_bin_histos("AbsPE", vec_AbsPBins);
    AbsPF   = make_bin_histos("AbsPF", vec_AbsPBins);

    Eta     = make_bin_histos("Eta", vec_EtaBins);
    Phi     = make_bin_histos("Phi", vec_PhiBins);
    Vtx     = make_bin_histos("Vtx", vec_VtxBins);
    if(isAOD) {
      ShowerB      = make_bin_histos("ShowerB", vec_ShowerBins);
      ShowerBPt200 = make_bin_histos("ShowerBPt200", vec_ShowerBins);
      ShowerBPt400 = make_bin_histos("ShowerBPt400", vec_ShowerBins);
      ShowerBP200  = make_bin_histos("ShowerBP200", vec_ShowerBins);
      ShowerBP400  = make_bin_histos("ShowerBP400", vec_ShowerBins);
      ShowerE      = make_bin_histos("ShowerE", vec_ShowerBins);
      ShowerEPt200 = make_bin_histos("ShowerEPt200", vec_ShowerBins);
      ShowerEPt400 = make_bin_histos("ShowerEPt400", vec_ShowerBins);
      ShowerEP200  = make_bin_histos("ShowerEP200", vec_ShowerBins);
      ShowerEP400  = make_bin_histos("ShowerEP400", vec_ShowerBins);
    }
  }

  // Templet histograms
  double *arr_PtBins;
  double *arr_AbsPBins;
  double *arr_EtaBins;
  double *arr_PhiBins;
  double *arr_VtxBins;
  double *arr_ShowerBins;
  CopyVectorToArray(vec_PtBins, arr_PtBins);
  CopyVectorToArray(vec_AbsPBins, arr_AbsPBins);
  CopyVectorToArray(vec_EtaBins, arr_EtaBins);
  CopyVectorToArray(vec_PhiBins, arr_PhiBins);
  CopyVectorToArray(vec_VtxBins, arr_VtxBins);
  CopyVectorToArray(vec_ShowerBins, arr_ShowerBins);
  hEffTemplatePt     = fs->make<TH1F>("hEffTemplatePt", "", (int)vec_PtBins.size()-1, arr_PtBins);
  hEffTemplateAbsP   = fs->make<TH1F>("hEffTemplateAbsP", "", (int)vec_AbsPBins.size()-1, arr_AbsPBins);
  hEffTemplateEta    = fs->make<TH1F>("hEffTemplateEta", "", (int)vec_EtaBins.size()-1, arr_EtaBins);
  hEffTemplatePhi    = fs->make<TH1F>("hEffTemplatePhi", "", (int)vec_PhiBins.size()-1, arr_PhiBins);
  hEffTemplateVtx    = fs->make<TH1F>("hEffTemplateVtx", "", (int)vec_VtxBins.size()-1, arr_VtxBins);
  hEffTemplateShower = fs->make<TH1F>("hEffTemplateShower", "", (int)vec_ShowerBins.size()-1, arr_ShowerBins);

  comparison_tree = fs->make<TTree>("t", "");
  comparison_tree->Branch("IsRealData", &IsRealData, "IsRealData/O");
  comparison_tree->Branch("RunNum",&RunNum,"RunNum/I");
  comparison_tree->Branch("LumiBlockNum",&LumiBlockNum,"LumiBlockNum/I");
  comparison_tree->Branch("EventNum",&EventNum,"EventNum/l");
  comparison_tree->Branch("Mass",&Mass,"Mass/D");
  comparison_tree->Branch("VertexMass",&VertexMass,"VertexMass/D");
  comparison_tree->Branch("Probe_Pt",&Probe_Pt,"Probe_Pt/D");
  comparison_tree->Branch("Probe_Eta",&Probe_Eta,"Probe_Eta/D");
}

void CustoTnPHistosForTnP::getBSandPV(const edm::Event& event) {
  // We store these as bare pointers. Should find better way, but
  // don't want to pass them around everywhere...
  edm::Handle<reco::BeamSpot> hbs;
  event.getByLabel(beamspot_src, hbs);
  beamspot = hbs.isValid() ? &*hbs : 0; // nice and fragile
  NBeamSpot->Fill(beamspot != 0);

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
  NVertices->Fill(vertex_count, _totalWeight );
}

std::vector<int> CustoTnPHistosForTnP::calcNShowers(
                                      const reco::CandidateBaseRef& mu,
                                      int min_barrel_st1 = 26,
                                      int min_barrel_st2 = 26,
                                      int min_barrel_st3 = 26,
                                      int min_barrel_st4 = 18,
                                      int min_endcap_st1 = 30,
                                      int min_endcap_st2 = 18,
                                      int min_endcap_st3 = 18,
                                      int min_endcap_st4 = 18,
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
  int st1 = (etaCat==1) ? muPat->userInt("nHits1")/2 : muPat->userInt("nHits1");
  int st2 = (etaCat==1) ? muPat->userInt("nHits2")/2 : muPat->userInt("nHits2");
  int st3 = (etaCat==1) ? muPat->userInt("nHits3")/2 : muPat->userInt("nHits3");
  int st4 = (etaCat==1) ? muPat->userInt("nHits4")/2 : muPat->userInt("nHits4");

  int threshold_st1 = -999;
  int threshold_st2 = -999;
  int threshold_st3 = -999;
  int threshold_st4 = -999;

  if( etaCat==1 ) {
    threshold_st1 = min_barrel_st1;
    threshold_st2 = min_barrel_st2;
    threshold_st3 = min_barrel_st3;
    threshold_st4 = min_barrel_st4;
  }
  else if( etaCat==2 ) {
    threshold_st1 = min_endcap_st1;
    threshold_st2 = min_endcap_st2;
    threshold_st3 = min_endcap_st3;
    threshold_st4 = min_endcap_st4;
  }

  bool is_st1 = (st1 > threshold_st1);
  bool is_st2 = (st2 > threshold_st2);
  bool is_st3 = (st3 > threshold_st3);
  bool is_st4 = (st4 > threshold_st4);

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

void CustoTnPHistosForTnP::fillTnPControlHistos(const pat::CompositeCandidate& dil,
                                                 const reco::CandidateBaseRef& TagMu,
                                                 const reco::CandidateBaseRef& ProbeMu,
                                                 float probe_dpt_over_pt,
                                                 int   probe_nshowers,
                                                 bool  isPassNoPt ) {

  TagPt->Fill( TagMu->pt(), _totalWeight );
  TagEta->Fill( TagMu->eta(), _totalWeight );
  TagPhi->Fill( TagMu->phi(), _totalWeight );

  for(int i=0; i<3; ++i) {
    bool doFill = false;
    if( i==0 )
      doFill = true;
    else if( i==1 && isPassNoPt )
      doFill = true;
    else if( i==2 && !isPassNoPt )
      doFill = true;

    if(doFill) {
      ProbePt[i]->Fill( ProbeMu->pt(), _totalWeight );
      PairNoPtMass[i]->Fill( dil.mass(), _totalWeight );

      if(ProbeMu->pt() > probe_pt_min) {

        ProbeEta[i]->Fill( ProbeMu->eta(), _totalWeight );
        ProbePhi[i]->Fill( ProbeMu->phi(), _totalWeight );
        ProbeNVertices[i]->Fill( nVtx, _totalWeight );
        ProbeEtaPhi[i]->Fill( ProbeMu->eta(), ProbeMu->phi(), _totalWeight );
        ProbeEtaPt[i]->Fill( ProbeMu->eta(), ProbeMu->pt(), _totalWeight );
        ProbeEtaDptPt[i]->Fill( ProbeMu->eta(), probe_dpt_over_pt, _totalWeight );
        if(isAOD && probe_nshowers>-1) {

          const pat::Muon* muPat = toConcretePtr<pat::Muon>(ProbeMu);

          int nSt1 = (fabs(ProbeMu->eta())<0.9) ? muPat->userInt("nHits1")/2 : muPat->userInt("nHits1");
          int nSt2 = (fabs(ProbeMu->eta())<0.9) ? muPat->userInt("nHits2")/2 : muPat->userInt("nHits2");
          int nSt3 = (fabs(ProbeMu->eta())<0.9) ? muPat->userInt("nHits3")/2 : muPat->userInt("nHits3");
          int nSt4 = (fabs(ProbeMu->eta())<0.9) ? muPat->userInt("nHits4")/2 : muPat->userInt("nHits4");

          ProbeEtaShower[i]->Fill( ProbeMu->eta(), probe_nshowers, _totalWeight );
          ProbeEtaHitsSt1[i]->Fill( ProbeMu->eta(), nSt1, _totalWeight );
          ProbeEtaHitsSt2[i]->Fill( ProbeMu->eta(), nSt2, _totalWeight );
          ProbeEtaHitsSt3[i]->Fill( ProbeMu->eta(), nSt3, _totalWeight );
          ProbeEtaHitsSt4[i]->Fill( ProbeMu->eta(), nSt4, _totalWeight );

          if(fabs(ProbeMu->eta())<0.9) {
            ProbePtShowerB[i]->Fill( ProbeMu->pt(), probe_nshowers, _totalWeight );
            ProbePtHitsSt1B[i]->Fill( ProbeMu->pt(), nSt1, _totalWeight );
            ProbePtHitsSt2B[i]->Fill( ProbeMu->pt(), nSt2, _totalWeight );
            ProbePtHitsSt3B[i]->Fill( ProbeMu->pt(), nSt3, _totalWeight );
            ProbePtHitsSt4B[i]->Fill( ProbeMu->pt(), nSt4, _totalWeight );

            ProbePShowerB[i]->Fill( ProbeMu->p(), probe_nshowers, _totalWeight );
            ProbePHitsSt1B[i]->Fill( ProbeMu->p(), nSt1, _totalWeight );
            ProbePHitsSt2B[i]->Fill( ProbeMu->p(), nSt2, _totalWeight );
            ProbePHitsSt3B[i]->Fill( ProbeMu->p(), nSt3, _totalWeight );
            ProbePHitsSt4B[i]->Fill( ProbeMu->p(), nSt4, _totalWeight );
          }
          else if(fabs(ProbeMu->eta())>=1.2) {
            ProbePtShowerE[i]->Fill( ProbeMu->pt(), probe_nshowers, _totalWeight );
            ProbePtHitsSt1E[i]->Fill( ProbeMu->pt(), nSt1, _totalWeight );
            ProbePtHitsSt2E[i]->Fill( ProbeMu->pt(), nSt2, _totalWeight );
            ProbePtHitsSt3E[i]->Fill( ProbeMu->pt(), nSt3, _totalWeight );
            ProbePtHitsSt4E[i]->Fill( ProbeMu->pt(), nSt4, _totalWeight );

            ProbePShowerE[i]->Fill( ProbeMu->p(), probe_nshowers, _totalWeight );
            ProbePHitsSt1E[i]->Fill( ProbeMu->p(), nSt1, _totalWeight );
            ProbePHitsSt2E[i]->Fill( ProbeMu->p(), nSt2, _totalWeight );
            ProbePHitsSt3E[i]->Fill( ProbeMu->p(), nSt3, _totalWeight );
            ProbePHitsSt4E[i]->Fill( ProbeMu->p(), nSt4, _totalWeight );
          }
        }

        PairMass[i]->Fill( dil.mass(), _totalWeight );
        PairPt[i]->Fill( dil.pt(), _totalWeight );
        PairEta[i]->Fill( dil.eta(), _totalWeight );
        PairRap[i]->Fill( dil.rapidity(), _totalWeight );

      }
    }
  }

}

void CustoTnPHistosForTnP::fillTnPBinHistos(  double               dil_mass,
                                              double               probe_pt,
                                              bool                 isPassNoPt,
                                              double               binValue,
                                              std::vector<double>& vec_bins,
                                              BinHistos&           Histos,
                                              bool hasPtCut = true ) {

  if( hasPtCut && (probe_pt < probe_pt_min) )
    return;

  bool isPass = false;

  if( hasPtCut )
    isPass = ( isPassNoPt && (probe_pt > probe_pt_min) );
  else
    isPass = isPassNoPt;

  int nbins = ( (int)vec_bins.size() ) - 1;

  for(int i=0; i<nbins; ++i) {
    if( binValue >= vec_bins[i] && binValue < vec_bins[i+1] ) {

      if(isPass)
        (Histos.first)[i]->Fill( dil_mass, _totalWeight );
      else
        (Histos.second)[i]->Fill( dil_mass, _totalWeight );

      break;
    }
  }

}

void CustoTnPHistosForTnP::fillTnPBinHistos2D(  double               dil_mass,
                                                double               probe_pt,
                                                bool                 isPassNoPt,
                                                double               binValue,
                                                BinHistos2D&         Histos,
                                                bool hasPtCut = true ) {

  if( hasPtCut && (probe_pt < probe_pt_min) )
    return;

  bool isPass = false;

  if( hasPtCut )
    isPass = ( isPassNoPt && (probe_pt > probe_pt_min) );
  else
    isPass = isPassNoPt;

  if(isPass)
    Histos[0]->Fill( dil_mass, binValue, _totalWeight );
  else
    Histos[1]->Fill( dil_mass, binValue, _totalWeight );

  Histos[2]->Fill( dil_mass, binValue, binValue*_totalWeight );
  Histos[3]->Fill( dil_mass, binValue, binValue*binValue*_totalWeight );

}


void CustoTnPHistosForTnP::analyze(const edm::Event& event, const edm::EventSetup& setup) {

  IsRealData = false;
  RunNum = -999;
  LumiBlockNum = -999;
  EventNum = 0;
  Mass = -999.;
  VertexMass = -999.;
  Probe_Pt = -999.;
  Probe_Eta = -999.;

  //---- Prescales : Not using for now...
    //  edm::Handle<int> hltPrescale;
    //  edm::Handle<int> l1Prescale;
    //  event.getByLabel(edm::InputTag("getPrescales","HLTPrescale","CustoTnPAnalysis"), hltPrescale);
    //  event.getByLabel(edm::InputTag("getPrescales","L1Prescale","CustoTnPAnalysis"), l1Prescale);
    if (_usePrescaleWeight) {
      edm::Handle<int> totalPrescale;
      event.getByLabel(edm::InputTag("getPrescales","TotalPrescale","CustoTnPAnalysis"), totalPrescale);
      _prescaleWeight = *totalPrescale;
    }
    //  std::cout<<*hltPrescale<<std::endl;
    //  std::cout<<l1Prescale<<std::endl;
    //  std::cout<<totalPrescale<<std::endl;

  //---- Generator weights
  if (_useMadgraphWeight) {
    _eventWeight = 1.;
    _madgraphWeight = 1.;

    edm::Handle<GenEventInfoProduct> gen_ev_info;
    event.getByLabel(edm::InputTag("generator"), gen_ev_info);
    if (gen_ev_info.isValid()){
      _eventWeight = gen_ev_info->weight();
      _madgraphWeight = ( _eventWeight > 0 ) ? 1.0 : -1.0;
    }
    WeightMadGraph->Fill( _madgraphWeight );
  }

  //---- Get PileUp Weights
  int thePU = -1;
  _pileupWeight = 1.;
  if( !event.isRealData() ){

    edm::Handle<std::vector< PileupSummaryInfo > > pileup;
    if ( event.getByLabel(pileup_src, pileup)){
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      for(PVI = pileup->begin(); PVI != pileup->end(); ++PVI)
      {
        if(PVI->getBunchCrossing()==0){
          thePU = PVI->getTrueNumInteractions();
          continue;
        }
      }

      int nPileUpBin = (int)vec_PileUpWeight.size();
      for(int iPU=0; iPU<nPileUpBin; ++iPU) {
        if(thePU == iPU) {
          _pileupWeight = vec_PileUpWeight[iPU];
        }
      }
    }
    else
      edm::LogError("") << "PU collection not found !!!";

    if(!ShutUp)  std::cout << "CustoTnPHistosForTnP::analyze : PU = " << thePU << "  PU weight = " << _pileupWeight << std::endl;
  }

  _totalWeight = 1.;
  if( _madgraphWeight != 1. || _pileupWeight != 1. )
    _totalWeight = _madgraphWeight * _pileupWeight;  // prescale weight is not considered
  if(!ShutUp)  std::cout << "CustoTnPHistosForTnP::analyze : Total weight = " << _totalWeight << std::endl;

  if (use_bs_and_pv)
    getBSandPV(event);

  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  event.getByLabel(dilepton_src, dileptons);

  if( !dileptons.isValid() ) {
    std::cout << "CustoTnPHistosForTnP::analyze : !dileptons.isValid() ---> return" << std::endl;
    return;
  }

  NDileptons->Fill(dileptons->size(), _totalWeight );

  pat::CompositeCandidateCollection::const_iterator dil = dileptons->begin(), dile = dileptons->end();
  for ( ; dil != dile; ++dil) {

    if( !( dil->hasUserInt("isTag0Probe1") && dil->hasUserInt("isTag1Probe0") ) ) {
      std::cout << "CustoTnPHistosForTnP::analyze : no isTag0Probe1 or isTag1Probe0 in dil" << std::endl;
      continue;
    }

    double dil_mass = dil->mass();
    double dil_vertex_mass = dil->userFloat("vertexM");
    float lep0_dpt_over_pt = dil->userFloat("lep0_dpt_over_pt");
    float lep1_dpt_over_pt = dil->userFloat("lep1_dpt_over_pt");

    const reco::CandidateBaseRef& lep0 = dileptonDaughter(*dil, 0);
    const reco::CandidateBaseRef& lep1 = dileptonDaughter(*dil, 1);

    if(lep0.isNonnull() && lep1.isNonnull()) {
      const pat::Muon* mu0 = toConcretePtr<pat::Muon>(lep0);
      const pat::Muon* mu1 = toConcretePtr<pat::Muon>(lep1);
      if(mu0 && mu1) {

        bool isAleadyFilled = false;

        //---- Tag0 and Probe1
        if( dil->userInt("isTag0Probe1") ) {
          const reco::CandidateBaseRef& TagMu = lep0;
          const reco::CandidateBaseRef& ProbeMu = lep1;

          int probe_nshowers = -999;
          if(isAOD) {
            std::vector<int> vec_showers = calcNShowers(ProbeMu);

            if(vec_showers.size()==16) {
              for(int is=0; is<(int)vec_showers.size(); ++is) {
                if(vec_showers[is]) {
                  probe_nshowers = is;
                  break;
                }
              }
            }

          }

          if(!ShutUp)  std::cout << "CustoTnPHistosForTnP::analyze : Tag0 and Probe1" << std::endl;
          if(!ShutUp)  std::cout << "                                              pT=" << ProbeMu->pt() << std::endl;
          if(!ShutUp)  std::cout << "                                             eta=" << ProbeMu->eta() << std::endl;
          if(!ShutUp)  std::cout << "                                             phi=" << ProbeMu->phi() << std::endl;
          if(!ShutUp)  std::cout << "                                            nVtx=" << nVtx << std::endl;

          bool isPassingProbe = ( passing_probe_selector(*mu1) && (lep1_dpt_over_pt < passing_probe_dpt_over_pt_max) && (fabs(mu1->innerTrack()->dz( vertex->position() )) < passing_probe_dz_max) );
          if(!ShutUp)  std::cout << "                                  isPassingProbe = " << isPassingProbe << std::endl;

          fillTnPControlHistos(*dil, TagMu, ProbeMu, lep1_dpt_over_pt, probe_nshowers, isPassingProbe);

          if(useBinHistos2D) {
            fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), Pt2D, false);
            if( fabs(ProbeMu->eta())<0.9 )
              fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), PtB2D, false);
            else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
              fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), PtO2D, false);
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.1 )
              fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), PtE2D, false);
            else if( fabs(ProbeMu->eta())>=2.1 && fabs(ProbeMu->eta())<2.4 )
              fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), PtF2D, false);

            fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsP2D, true);
            if( fabs(ProbeMu->eta())<0.9 )
              fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsPB2D, true);
            else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
              fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsPO2D, true);
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.1 )
              fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsPE2D, true);
            else if( fabs(ProbeMu->eta())>=2.1 && fabs(ProbeMu->eta())<2.4 )
              fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsPF2D, true);

            fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->eta(), Eta2D, true);
            fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->phi(), Phi2D, true);
            fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, nVtx, Vtx2D, true);

            if(isAOD && probe_nshowers>-1) {
              if( fabs(ProbeMu->eta())<0.9 ) {
                fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerB2D, true);
                if(ProbeMu->pt()>=200 && ProbeMu->pt()<400)
                  fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerBPt2002D, true);
                if(ProbeMu->pt()>=400)
                  fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerBPt4002D, true);
                if(ProbeMu->p()>=200 && ProbeMu->p()<400)
                  fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerBP2002D, true);
                if(ProbeMu->p()>=400)
                  fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerBP4002D, true);
              }
              else if( fabs(ProbeMu->eta())>=1.2 ) {
                fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerE2D, true);
                if(ProbeMu->pt()>=200 && ProbeMu->pt()<400)
                  fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerEPt2002D, true);
                if(ProbeMu->pt()>=400)
                  fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerEPt4002D, true);
                if(ProbeMu->p()>=200 && ProbeMu->p()<400)
                  fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerEP2002D, true);
                if(ProbeMu->p()>=400)
                  fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerEP4002D, true);
              }
            }
          }
          else {
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, Pt, false);
            if( fabs(ProbeMu->eta())<0.9 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, PtB, false);
            else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, PtO, false);
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.1 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, PtE, false);
            else if( fabs(ProbeMu->eta())>=2.1 && fabs(ProbeMu->eta())<2.4 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, PtF, false);

            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsP, true);
            if( fabs(ProbeMu->eta())<0.9 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsPB, true);
            else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsPO, true);
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.1 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsPE, true);
            else if( fabs(ProbeMu->eta())>=2.1 && fabs(ProbeMu->eta())<2.4 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsPF, true);

            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->eta(), vec_EtaBins, Eta, true);
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->phi(), vec_PhiBins, Phi, true);
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, nVtx, vec_VtxBins, Vtx, true);

            if(isAOD && probe_nshowers>-1) {
              if( fabs(ProbeMu->eta())<0.9 ) {
                fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerB, true);
                if(ProbeMu->pt()>=200 && ProbeMu->pt()<400)
                  fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerBPt200, true);
                if(ProbeMu->pt()>=400)
                  fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerBPt400, true);
                if(ProbeMu->p()>=200 && ProbeMu->p()<400)
                  fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerBP200, true);
                if(ProbeMu->p()>=400)
                  fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerBP400, true);
              }
              else if( fabs(ProbeMu->eta())>=1.2 ) {
                fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerE, true);
                if(ProbeMu->pt()>=200 && ProbeMu->pt()<400)
                  fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerEPt200, true);
                if(ProbeMu->pt()>=400)
                  fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerEPt400, true);
                if(ProbeMu->p()>=200 && ProbeMu->p()<400)
                  fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerEP200, true);
                if(ProbeMu->p()>=400)
                  fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerEP400, true);
              }
            }
          }

          if( ProbeMu->pt() > probe_pt_min ) {
            bool isComparisonProbe = ( comparison_probe_selector(*mu1) && (lep1_dpt_over_pt < comparison_probe_dpt_over_pt_max) && (fabs(mu1->innerTrack()->dz( vertex->position() )) < comparison_probe_dz_max) );
            if(!ShutUp)  std::cout << "                               isComparisonProbe = " << isComparisonProbe << std::endl;
            if( isPassingProbe && !isComparisonProbe ) {
              IsRealData = event.isRealData();
              RunNum = event.id().run();
              LumiBlockNum = event.id().luminosityBlock();
              EventNum = event.id().event();
              Mass = dil_mass;
              VertexMass = dil_vertex_mass;
              Probe_Pt = ProbeMu->pt();
              Probe_Eta = ProbeMu->eta();
              comparison_tree->Fill();
              isAleadyFilled = true;
            }
          }
        } // Tag0 and Probe1

        //---- Tag1 and Probe0
        if( dil->userInt("isTag1Probe0") ) {
          const reco::CandidateBaseRef& TagMu = lep1;
          const reco::CandidateBaseRef& ProbeMu = lep0;

          int probe_nshowers = -999;
          if(isAOD) {
            std::vector<int> vec_showers = calcNShowers(ProbeMu);

            if(vec_showers.size()==16) {
              for(int is=0; is<(int)vec_showers.size(); ++is) {
                if(vec_showers[is]) {
                  probe_nshowers = is;
                  break;
                }
              }
            }

          }

          if(!ShutUp)  std::cout << "CustoTnPHistosForTnP::analyze : Tag1 and Probe0" << std::endl;
          if(!ShutUp)  std::cout << "                                              pT=" << ProbeMu->pt() << std::endl;
          if(!ShutUp)  std::cout << "                                             eta=" << ProbeMu->eta() << std::endl;
          if(!ShutUp)  std::cout << "                                             phi=" << ProbeMu->phi() << std::endl;
          if(!ShutUp)  std::cout << "                                            nVtx=" << nVtx << std::endl;

          bool isPassingProbe = ( passing_probe_selector(*mu0) && (lep0_dpt_over_pt < passing_probe_dpt_over_pt_max) && (fabs(mu0->innerTrack()->dz( vertex->position() )) < passing_probe_dz_max) );
          if(!ShutUp)  std::cout << "                                  isPassingProbe = " << isPassingProbe << std::endl;

          fillTnPControlHistos(*dil, TagMu, ProbeMu, lep0_dpt_over_pt, probe_nshowers, isPassingProbe);

          if(useBinHistos2D) {
            fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), Pt2D, false);
            if( fabs(ProbeMu->eta())<0.9 )
              fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), PtB2D, false);
            else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
              fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), PtO2D, false);
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.1 )
              fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), PtE2D, false);
            else if( fabs(ProbeMu->eta())>=2.1 && fabs(ProbeMu->eta())<2.4 )
              fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), PtF2D, false);

            fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsP2D, true);
            if( fabs(ProbeMu->eta())<0.9 )
              fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsPB2D, true);
            else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
              fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsPO2D, true);
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.1 )
              fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsPE2D, true);
            else if( fabs(ProbeMu->eta())>=2.1 && fabs(ProbeMu->eta())<2.4 )
              fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsPF2D, true);

            fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->eta(), Eta2D, true);
            fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->phi(), Phi2D, true);
            fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, nVtx, Vtx2D, true);

            if(isAOD && probe_nshowers>-1) {
              if( fabs(ProbeMu->eta())<0.9 ) {
                fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerB2D, true);
                if(ProbeMu->pt()>=200 && ProbeMu->pt()<400)
                  fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerBPt2002D, true);
                if(ProbeMu->pt()>=400)
                  fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerBPt4002D, true);
                if(ProbeMu->p()>=200 && ProbeMu->p()<400)
                  fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerBP2002D, true);
                if(ProbeMu->p()>=400)
                  fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerBP4002D, true);
              }
              else if( fabs(ProbeMu->eta())>=1.2 ) {
                fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerE2D, true);
                if(ProbeMu->pt()>=200 && ProbeMu->pt()<400)
                  fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerEPt2002D, true);
                if(ProbeMu->pt()>=400)
                  fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerEPt4002D, true);
                if(ProbeMu->p()>=200 && ProbeMu->p()<400)
                  fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerEP2002D, true);
                if(ProbeMu->p()>=400)
                  fillTnPBinHistos2D(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerEP4002D, true);
              }
            }
          }
          else {
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, Pt, false);
            if( fabs(ProbeMu->eta())<0.9 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, PtB, false);
            else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, PtO, false);
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.1 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, PtE, false);
            else if( fabs(ProbeMu->eta())>=2.1 && fabs(ProbeMu->eta())<2.4 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, PtF, false);

            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsP, true);
            if( fabs(ProbeMu->eta())<0.9 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsPB, true);
            else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsPO, true);
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.1 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsPE, true);
            else if( fabs(ProbeMu->eta())>=2.1 && fabs(ProbeMu->eta())<2.4 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsPF, true);

            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->eta(), vec_EtaBins, Eta, true);
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->phi(), vec_PhiBins, Phi, true);
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, nVtx, vec_VtxBins, Vtx, true);

            if(isAOD && probe_nshowers>-1) {
              if( fabs(ProbeMu->eta())<0.9 ) {
                fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerB, true);
                if(ProbeMu->pt()>=200 && ProbeMu->pt()<400)
                  fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerBPt200, true);
                if(ProbeMu->pt()>=400)
                  fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerBPt400, true);
                if(ProbeMu->p()>=200 && ProbeMu->p()<400)
                  fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerBP200, true);
                if(ProbeMu->p()>=400)
                  fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerBP400, true);
              }
              else if( fabs(ProbeMu->eta())>=1.2 ) {
                fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerE, true);
                if(ProbeMu->pt()>=200 && ProbeMu->pt()<400)
                  fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerEPt200, true);
                if(ProbeMu->pt()>=400)
                  fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerEPt400, true);
                if(ProbeMu->p()>=200 && ProbeMu->p()<400)
                  fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerEP200, true);
                if(ProbeMu->p()>=400)
                  fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerEP400, true);
              }
            }
          }

          if( ProbeMu->pt() > probe_pt_min ) {
            bool isComparisonProbe = ( comparison_probe_selector(*mu0) && (lep0_dpt_over_pt < comparison_probe_dpt_over_pt_max) && (fabs(mu0->innerTrack()->dz( vertex->position() )) < comparison_probe_dz_max) );
            if(!ShutUp)  std::cout << "                               isComparisonProbe = " << isComparisonProbe << std::endl;
            if( !isAleadyFilled && isPassingProbe && !isComparisonProbe ) {
              IsRealData = event.isRealData();
              RunNum = event.id().run();
              LumiBlockNum = event.id().luminosityBlock();
              EventNum = event.id().event();
              Mass = dil_mass;
              VertexMass = dil_vertex_mass;
              Probe_Pt = ProbeMu->pt();
              Probe_Eta = ProbeMu->eta();
              comparison_tree->Fill();
            }
          }
        } // Tag1 and Probe0

        if(!ShutUp)  std::cout << std::endl;
      } // if(mu0 && mu1)
    } // if(lep0.isNonnull() && lep1.isNonnull())


  } // for ( ; dil != dile; ++dil)

}

DEFINE_FWK_MODULE(CustoTnPHistosForTnP);
