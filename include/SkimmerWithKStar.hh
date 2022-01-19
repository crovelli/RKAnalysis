#ifndef SkimmerWithKStar_h
#define SkimmerWithKStar_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <vector>
#include <string>
#include <TLorentzVector.h>

#include "./BParkBase.h"

using namespace std;

class SkimmerWithKStar : public BParkBase{
public:

  //! constructor
  SkimmerWithKStar(TChain *tree=0, int isMC=0 );
  //! destructor
  virtual ~SkimmerWithKStar();
  //! loop over events
  void Loop();
  void PrepareOutputs(std::string filename);             
  float DeltaR(float x, float y,float z,float w);
  float DeltaPhi(float x, float y);
  float _ptB;
  float _svProb;
  float _cos2D;
  void SetCuts(float ptB, float svProb, float cos2D){
    _ptB=ptB; _svProb=svProb; _cos2D=cos2D; 
  }

private:
  
  // Analysis methods
  bool isMcB( int myB );
  bool isMcBOtto( int myB );
  bool isMcEleFromJPsi (int myEle );
  void bookOutputTree();
  void bookOutputHistos();
  void SetNvtxWeights(std::string nvtxWeightFile);
  float GetNvtxWeight(float nvtx);
  void SetAnBdtWeights();
  float GetAnBdtWeight(float bdt, int pfpf);

  // to compute weights for pileup
  std::vector<Double_t> nvtxweights_;
  std::vector<Double_t> nvtxlowedge_;

  // to compute systematics
  std::vector<Double_t> anweight_pflpt_common_;
  std::vector<Double_t> anweight_pfpf_common_;
  std::vector<Double_t> lowedge_;

  // settings
  bool dosyst_;
  bool donvtxreweight_;
  int sampleID;
  string nvtxWFileName_;
  float lumiWeight_;

  // ---- outputs
  TFile* outFile_;
  TTree* outTree_;
  TH1F* h_entries;
  TH1F* h_selection;

  // dataset name
  std::string _datasetName;      

  //---output tree branches variables
  int   theRun;
  int   theEvent;
  int   nvtx;
  int   theSampleID;
  float rho;
  float pu_weight;
  int   hlt7;
  int   hlt9;
  float trg_muon_pt;
  // 
  int selectedBSize;
  //
  vector <float> p4Trk={};   
  vector <float> p4TrkKStar={};   
  //
  vector <float> bXYOrder={};   
  vector <float> bPtOrder={};   
  vector <float> b2DOrder={};   
  //
  vector <float> tag_pt={};   
  vector <float> tag_eta={};
  vector <float> tag_phi={};
  vector <bool> tag_isPF={};
  vector <bool> tag_isPFOverlap={};
  vector <bool> tag_isLowPt={};
  vector <float> tag_id={};
  vector <float> tag_drm={}; 
  //
  vector <bool>  tag_matchMcFromJPsi={};
  vector <bool>  tag_matchMc={};
  vector <float> tag_ptMc={};
  //
  vector <float> probe_pt={}; 
  vector <float> probe_eta ={};
  vector <float> probe_phi={};
  vector <bool> probe_isPF ={};
  vector <bool> probe_isPFOverlap ={};
  vector <bool> probe_isLowPt ={};
  vector <float> probe_id={};        
  vector <float> probe_drm={}; 
  //
  vector <bool>  probe_matchMcFromJPsi={};
  vector <bool>  probe_matchMc={};
  vector <float> probe_ptMc={};
  //
  vector <float> K_pt={}; 
  vector <float> K_eta={}; 
  vector <float> K_phi={}; 
  //
  vector <float> mll_fullfit={};
  vector <float> mll_raw={};
  //
  vector <float> fit_Bmass={};   
  vector <float> fit_Bpt={};   
  vector <float> fit_Bcos2D={};   
  vector <float> fit_Bsvprob={};
  vector <float> fit_Bxysig={};
  vector <float> fit_Bi3dsig={};
  vector <bool> bmatchMC={};   
  //
  vector <float> LKdz_vec={};
  vector <float> L1L2dr_vec={};
  vector <float> L1L2dr_raw_vec={};
  vector <float> b_iso04_rel_vec={};
  vector <float> eleDR_vec={};
  vector <float> k_DCASig_vec={};
  vector <float> k_dzTrg_vec={};
  vector <float> k_svip2d_vec={};
  vector <float> k_svip3d_vec={};
  vector <float> LKdr_vec={};
  vector <float> LKdr_raw_vec={};
  vector <float> BTrkdxy2_vec={};
  vector <float> BBDPhi_vec={};
  vector <float> BBDPhi_raw_vec={};
  vector <float> ptAsym_vec={};
  vector <float> ptAsym_raw_vec={};
  vector <float> k_iso04_rel_vec={};
  vector <float> tag_iso04_rel_vec={};
  vector <float> probe_iso04_rel_vec={};
  vector <float> k_dxy_sig_vec={}; 
  vector <float> tag_dxy_sig_vec={};
  vector <float> tag_dzTrg_vec={};
  vector <float> probe_dxy_sig_vec={};
  vector <float> probe_dzTrg_vec={};
  vector <float> llkDR_vec={};

  //vector <float> analysisBdtC={};
  //vector <float> analysisBdtCWithSyst={};

  vector <float> Dmass_vec={};
  vector <float> Dmass_flip_vec={};
  vector <float> Dmass_ll_flip_vec={};
  vector <float> Dmass_ll_vec={};
  vector <float> Dmass_mukhad_vec={};
  vector <float> Dmass_mukmu_vec={};
  vector <float> KLmassD0_vec={};
  //
  int selectedPairsSize;
};

#endif
