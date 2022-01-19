#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TClass.h>

#include <cmath>

#include "../FastForest/include/fastforest.h"
#include "../include/SkimmerWithKStar.hh"    

using namespace std;

float SkimmerWithKStar::DeltaR(float eta1,float  phi1, float eta2, float phi2){
  float PI=3.1415972;
  float deta=eta1-eta2;
  float dphi=phi1-phi2;
  if(dphi>PI){
    dphi-=2.0*PI;
  } else if(dphi<=-PI){
    dphi+=2.0*PI;
  }
  return TMath::Sqrt(deta*deta+dphi*dphi);
}

float SkimmerWithKStar::DeltaPhi(float  phi1, float phi2){
  float PI=3.1415972;
  float dphi=phi1-phi2;
  if(dphi>PI){
    dphi-=2.0*PI;
  } else if(dphi<=-PI){
    dphi+=2.0*PI;
  }
  return dphi;
}

SkimmerWithKStar::SkimmerWithKStar(TChain *tree, int isMC )     
  : BParkBase((TTree*) tree, isMC ) {        

  sampleID = isMC;         

  // To be set by hand - chiara
  donvtxreweight_ = 0;
  nvtxWFileName_ = "/afs/cern.ch/user/c/crovelli/public/bphys/march/nvtxWeights2018ALL.root";

  // To be set by hand - chiara
  dosyst_ = 0;
}

SkimmerWithKStar::~SkimmerWithKStar() {

  // output
  outFile_->cd();
  h_entries -> Write();
  h_selection -> Write();
  outTree_->Write();
  outFile_->Close();
}     

void SkimmerWithKStar::Loop() {

  if (fChain == 0) return;

  /*
  // Load Analysis MVA weights for common BDT
  std::string bdtfileCommonPFPF = "models/xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.2_0.txt";
  std::vector<std::string> featCommonPFPF = {"f0","f1","f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13","f14","f15","f16"};
  const auto fastForestCommonPFPF = fastforest::load_txt(bdtfileCommonPFPF.c_str(), featCommonPFPF);
  // Load Analysis MVA weights PFLP
  std::string bdtfileCommonPFLP = "models/xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.2_0.txt";
  std::vector<std::string> featCommonPFLP = {"f0","f1","f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13","f14","f15","f16"};
  const auto fastForestCommonPFLP = fastforest::load_txt(bdtfileCommonPFLP.c_str(), featCommonPFLP);
  */

  // Loop over events
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  cout << "entries : " <<  nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%500==0) cout << jentry << endl;

    // To keep track of the total number of events 
    h_entries->Fill(5);
    
    // Event info     
    theRun   = run;
    theEvent = event;
    theSampleID = sampleID;

    // # Vertices
    nvtx = PV_npvs;
    
    // Energy density
    rho = fixedGridRhoFastjetAll;
    
    // PU weight (for MC only and if requested)
    pu_weight = 1.;
    if (sampleID>0 && donvtxreweight_==1) {    
      pu_weight = GetNvtxWeight(nvtx);         
    }

    // other weights for the dataset
    float perEveW = 1.;
    if (sampleID>0) perEveW = Generator_weight;

    // Events breakdown  
    h_selection->Fill(0.,perEveW);
    

    // ----------------------------------------------------

    // PV must be there
    bool goodPV = false;
    if (PV_x>-999) goodPV = true;
    if (!goodPV) continue;
    h_selection->Fill(1.,perEveW);

    // Trigger 
    int iHLT_Mu7_IP4 = (int)HLT_Mu7_IP4;
    int iHLT_Mu9_IP6 = (int)HLT_Mu9_IP6;
    hlt7 = iHLT_Mu7_IP4;
    hlt9 = iHLT_Mu9_IP6;
    
    // Triggering muons
    bool nTriggerMuon=0;
    int idTrgMu=-1;
    float ptTrgMu=0.;
    trg_muon_pt=0.;
    if(nMuon>0){
      for (u_int iM=0; iM<nMuon; iM++) {
	if(Muon_isTriggering[iM]) {
	  nTriggerMuon=nTriggerMuon+1;
	  if(Muon_pt[iM]>ptTrgMu) {
	    ptTrgMu=Muon_pt[iM];
	    idTrgMu=iM;
	  }
	}
      }
    }
    trg_muon_pt=ptTrgMu;


    // Needed for calculation of BBDPhi
    TLorentzVector trgmuon_vec(0,0,0,0);
    TLorentzVector sum_track_vec(0,0,0,0);
    float sum_track=0;

    if(nTriggerMuon>0){
      trgmuon_vec.SetPtEtaPhiM(Muon_pt[idTrgMu], Muon_eta[idTrgMu],Muon_phi[idTrgMu],0.105);
      sum_track_vec.SetPtEtaPhiM(Muon_pt[idTrgMu], Muon_eta[idTrgMu],Muon_phi[idTrgMu],0.139);
      sum_track=Muon_pt[idTrgMu];
      
      for (int itrk=0; itrk<nProbeTracks; itrk++){
	TLorentzVector track_vec(0,0,0,0);
	track_vec.SetPtEtaPhiM(ProbeTracks_pt[itrk],ProbeTracks_eta[itrk],ProbeTracks_phi[itrk],0.139);
	if (trgmuon_vec.DeltaR(track_vec)<0.4){ 
	  sum_track_vec=sum_track_vec+track_vec;
	  sum_track=sum_track+ProbeTracks_pt[itrk];
	}
      }
    }


    // B candidates
    if (nBToKEE<=0) continue;
    h_selection->Fill(2.,perEveW);


    // order by pt
    u_int iBGd[1000]={0};
    float ptBG[1000]={0.};
    int n=nBToKEE;
    if(n>1000) n=1000; // max 1000 triplets
    for (u_int iB=0; iB<n; iB++) {
      int kmin=0;
      if(iB==0) {
	ptBG[0]=BToKEE_fit_pt[iB];
	iBGd[0]=iB;
      } else {
	if(BToKEE_fit_pt[iB]>ptBG[n-1]){ // it is higher than the last 
	  if(BToKEE_fit_pt[iB]>ptBG[0]) { 
	    kmin=0; // it is higher than any registered 
	  } else {
	    for(u_int k=0; k<n-1; k++){
	      if(BToKEE_fit_pt[iB]>ptBG[k+1] && BToKEE_fit_pt[iB] <= ptBG[k] ){
		kmin=k+1; // it is in the middle 
	      }
	    }
	  }
	  int nmax=n-2;
	  if(n-2 <=0 ) nmax=0;
	  for(int k=nmax; k>=kmin; k--){
	    ptBG[k+1]=ptBG[k];
	    iBGd[k+1]=iBGd[k];
	  }
	  ptBG[kmin]=BToKEE_fit_pt[iB];
	  iBGd[kmin]=iB;
	}
      }
    }
    
    // order by the cos2D    
    u_int iBGd2D[1000]={0};
    float c2BG[1000]={0.};    
    for (u_int iB=0; iB<n; iB++) {
      int kmin=0;
      if(iB==0) {
	c2BG[0]=BToKEE_fit_cos2D[iB];
	iBGd2D[0]=iB;
      } else {
	if(BToKEE_fit_cos2D[iB]>c2BG[n-1]){ // it is higher than the last 
	  if(BToKEE_fit_cos2D[iB]>c2BG[0]) { 
	    kmin=0; // it is higher than any registered 
	  } else {
	    for(u_int k=0; k<n-1; k++){
	      if(BToKEE_fit_cos2D[iB]>c2BG[k+1] && BToKEE_fit_cos2D[iB] <= c2BG[k] ){
		kmin=k+1; // it is in the middle 
	      }
	    }
	  }
	  int nmax=n-2;
	  if(n-2 <=0 ) nmax=0;
	  for(int k=nmax; k>=kmin; k--){
	    c2BG[k+1]=c2BG[k];
	    iBGd2D[k+1]=iBGd2D[k];
	  }
	  c2BG[kmin]=BToKEE_fit_cos2D[iB];
	  iBGd2D[kmin]=iB;
	}
      }
    }
    
    // order by the Lxy/sigma_Lxy
    u_int iBGdxy[1000]={0};
    float xyBG[1000]={0.};
    for (u_int iB=0; iB<n; iB++) {
      float b_xySig = BToKEE_l_xy[iB]/BToKEE_l_xy_unc[iB];        
      int kmin=0;
      if(iB==0) {
	xyBG[0]=b_xySig;
	iBGdxy[0]=iB;
      } else {
	if(b_xySig>xyBG[n-1]){ // it is higher than the last 
	  if(b_xySig>xyBG[0]) { 
	    kmin=0; // it is higher than any registered 
	  } else {
	    for(u_int k=0; k<n-1; k++){
	      if(b_xySig>xyBG[k+1] && b_xySig <= xyBG[k] ){
		kmin=k+1; // it is in the middle 
	      }
	    }
	  }
	  int nmax=n-2;
	  if(n-2 <=0 ) nmax=0;
	  for(int k=nmax; k>=kmin; k--){
	    xyBG[k+1]=xyBG[k];
	    iBGdxy[k+1]=iBGdxy[k];
	  }
	  xyBG[kmin]=b_xySig;
	  iBGdxy[kmin]=iB;
	}
      }
    }


    // Minimal Bcandidate requirements
    vector<int> goodBs;
    vector<int> goodBorder;
    vector<int> goodBorder2D;
    vector<int> goodBorderXY;

    for (u_int iB=0; iB<n; iB++) {

      // preparing variables for selection
      int ele1_idx = BToKEE_l1Idx[iB];
      int ele2_idx = BToKEE_l2Idx[iB];
      int k_idx    = BToKEE_kIdx[iB];
      
      float ele1_pt = BToKEE_fit_l1_pt[iB];
      float ele2_pt = BToKEE_fit_l2_pt[iB];
      float k_pt    = BToKEE_fit_k_pt[iB];
      
      float ele1_eta = BToKEE_fit_l1_eta[iB];     
      float ele2_eta = BToKEE_fit_l2_eta[iB];    
      float k_eta    = BToKEE_fit_k_eta[iB];    

      float ele1_phi = BToKEE_fit_l1_phi[iB];     
      float ele2_phi = BToKEE_fit_l2_phi[iB];    
      float k_phi    = BToKEE_fit_k_phi[iB];    

      // Electron ID, with swap if needed
      float L1id=-1000.;
      float L2id=-1000.;
      if( Electron_isPF[ele1_idx]==1 && Electron_isPF[ele2_idx]==1 ) { // PFPF: L1=leading, L2=subleading
	L1id=Electron_pfmvaId[ele1_idx];
	L2id=Electron_pfmvaId[ele2_idx];
      }
      if( Electron_isPF[ele1_idx]==1 && Electron_isPF[ele2_idx]==0 ) { // PFLP: L1=PF, L2=LP
	L1id=Electron_pfmvaId[ele1_idx];
	L2id=Electron_mvaId[ele2_idx];
      }
      if( Electron_isPF[ele1_idx]==0 && Electron_isPF[ele2_idx]==1 ) { // PFLP: L1=PF, L2=LP
	L1id=Electron_pfmvaId[ele2_idx];
	L2id=Electron_mvaId[ele1_idx];
      }

      // swap ele1/ele2 if needed: 
      // PFPF: L1=leading, L2=subleading  => no swap needed
      // PFLP, L1=PF, L2=LP               => no swap needed  
      // PFLP, L1=LP, L2=PF               => !! swap needed !!  
      if( !Electron_isPF[ele1_idx] && Electron_isPF[ele2_idx]) { 
	ele1_idx = BToKEE_l2Idx[iB];
	ele2_idx = BToKEE_l1Idx[iB];
	ele1_pt  = BToKEE_fit_l2_pt[iB];
	ele2_pt  = BToKEE_fit_l1_pt[iB];
	ele1_eta = BToKEE_fit_l2_eta[iB];
	ele2_eta = BToKEE_fit_l1_eta[iB];
	ele1_phi = BToKEE_fit_l2_phi[iB];
	ele2_phi = BToKEE_fit_l1_phi[iB];
      }

      float theXySig = BToKEE_l_xy[iB]/BToKEE_l_xy_unc[iB];        

      // BTrkdxy2
      float vx_b=BToKEE_vtx_x[iB];
      float vy_b=BToKEE_vtx_y[iB];
      float vz_b=BToKEE_vtx_z[iB];
      float distance_b_trk1=10;         
      float distance_b_trk2=10;         
      float mean_b_trk=0;      
      int imean_b_trk=0 ;
      float BTrkdxy2=0;
      for (int itrk=0; itrk<nProbeTracks; itrk++){
	if(itrk==k_idx) continue;
	float vx_trk=ProbeTracks_vx[itrk];
	float vy_trk=ProbeTracks_vy[itrk];
	float vz_trk=ProbeTracks_vz[itrk];
	float trk_eta=ProbeTracks_eta[itrk]; 
	float trk_phi=ProbeTracks_phi[itrk];
	
	if( (DeltaR(ele1_eta,ele1_phi,trk_eta,trk_phi)<0.03) || (DeltaR(ele2_eta,ele2_phi,trk_eta,trk_phi)<0.3) || ( DeltaR(k_eta,k_phi,trk_eta,trk_phi)<0.03)) continue;
	
	float distance_b_trk=TMath::Sqrt( pow((vx_b-vx_trk),2) + pow((vy_b-vy_trk),2) ) ;
	if (DeltaR(BToKEE_fit_eta[iB],BToKEE_fit_phi[iB],trk_eta,trk_phi)<0.4){
	  mean_b_trk=mean_b_trk+TMath::Sqrt( pow((vx_b-vx_trk),2) + pow((vy_b-vy_trk),2) );
	  imean_b_trk=imean_b_trk+1;
	}

	if(distance_b_trk<distance_b_trk2 && distance_b_trk<distance_b_trk1){
	  distance_b_trk2=distance_b_trk1;
	  distance_b_trk1=distance_b_trk;
	} else if(distance_b_trk<distance_b_trk2) {
	  distance_b_trk2=distance_b_trk;
	} 
      }      
      if(imean_b_trk>0) {
	mean_b_trk=mean_b_trk/imean_b_trk;
      } else {
	mean_b_trk=10.;
      }
      BTrkdxy2=distance_b_trk2;


      // Preselection: for PF-LPT remove B candidates which are also PF-PF
      if (Electron_isLowPt[ele1_idx] && Electron_isPFoverlap[ele1_idx]) continue;
      if (Electron_isLowPt[ele2_idx] && Electron_isPFoverlap[ele2_idx]) continue;

      // Preselection proposed by George & Otto on top of nanoAOD presel
      bool commonPresel = false;
      bool etaPresel = fabs(ele1_eta)<2.4 && fabs(ele2_eta)<2.4 && fabs(k_eta)<2.4;
      if (Electron_isPF[ele1_idx]==1 && Electron_isPF[ele2_idx]==1) { // PFPF
	commonPresel = (BToKEE_svprob[iB]>0.000001) && (BToKEE_fit_cos2D[iB]>0.85) && (theXySig>3.) && (k_pt>0.5) && (BTrkdxy2>0.000001) && (ele1_pt)>2.0 && (ele2_pt)>2.0 && (L1id)>-1.5 && (L2id)>-3;
      } else {
	commonPresel = (BToKEE_svprob[iB]>0.000001) && (BToKEE_fit_cos2D[iB]>0.95) && (fabs(BToKEE_k_svip3d[iB])<0.06) && (k_pt>0.5) && (BTrkdxy2>0.000001) && (ele1_pt)>2.0 && (ele2_pt>1.0) && (L1id>-2.0) && (L2id>0);
      }
      if (!etaPresel)    continue;
      if (!commonPresel) continue;


      // Further B selection (de-activated for now)
      /*
	bool vtxFitSel = BToKEE_fit_pt[iB]>_ptB && BToKEE_svprob[iB]>_svProb && BToKEE_fit_cos2D[iB]>_cos2D;
	bool ele1Sel = fabs(ele1_eta)<2.4 && ele1_pt>0.7;
	bool ele2Sel = fabs(ele2_eta)<2.4 && ele2_pt>0.7;
	bool kSel = k_pt>0.7 && fabs(k_eta)<2.4; 
	bool additionalSel = BToKEE_fit_mass[iB]>4.5 && BToKEE_fit_mass[iB]<6.0;
	bool isBsel = vtxFitSel && ele1Sel && ele2Sel && kSel && additionalSel;
      */


      // -------------------------------------------------------------------------
      // Minimal selection to reduce the ntuples size - chiara
      /*
      bool extraEle1Sel = fabs(ele1_eta)<2.4 && ele1_pt>1.0;   // questo non era applicato x produzione ntuple
      bool extraEle2Sel = fabs(ele2_eta)<2.4 && ele2_pt>1.0;   // questo non era applicato x produzione ntuple
      bool extraKSel = k_pt>0.7 && fabs(k_eta)<2.4;
      bool extraPFLPT = true;
      if (Electron_isLowPt[ele1_idx]==1 && Electron_isLowPt[ele2_idx]==1) extraPFLPT = false;     
      bool extraOverlap = true;
      if (Electron_isPFoverlap[ele1_idx]==1 || Electron_isPFoverlap[ele2_idx]==1) extraOverlap = false; 
      bool extraConv = true;
      if (Electron_convVeto[ele1_idx]==0 || Electron_convVeto[ele2_idx]==0) extraConv = false;    
      bool extraMll = true;
      if (BToKEE_mll_raw[iB]>3.4 || BToKEE_mll_raw[iB]<2.6) extraMll = false;
      bool extraIdLPt = true;
      if (Electron_mvaId[ele1_idx]<-3 || Electron_mvaId[ele2_idx]<-3) extraIdLPt = false; // questo non era applicato x produzione ntuple,
                                                                                          // ma c'e' un taglio a -2.5 che credo venga dai nani
      bool extraIdPF = true;
      if (Electron_pfmvaId[ele1_idx]<-3 || Electron_pfmvaId[ele2_idx]<-3) extraIdPF = false;  // applicato
      bool isExtraSel = extraEle1Sel && extraEle2Sel && extraKSel && extraPFLPT && extraOverlap && extraConv && extraMll && extraIdLPt && extraIdPF;
      if (!isExtraSel) continue;
      */


      // -------------------------------------------------------------------------
      // Order Bs with different algos
      int i2D=n;
      for (u_int m=0; m<n; m++){
	if(iBGd2D[m]==iB) {
	  i2D=m; 
	  break;
	}
      }
      int ik=n;
      for (u_int m=0; m<n; m++){
	if(iBGd[m]==iB) {
	  ik=m; 
	  break;
	}
      }
      int ixy=n;
      for (u_int m=0; m<n; m++){
	if(iBGdxy[m]==iB) {
	    ixy=m; 
	    break;
	}
      }      
      goodBs.push_back(iB);
      goodBorder.push_back(ik);
      goodBorder2D.push_back(i2D);
      goodBorderXY.push_back(ixy);
    }

    if (goodBs.size()>0) h_selection->Fill(3.,perEveW);
    
    selectedBSize = goodBs.size();
    
    // At least one good B candidate
    if (goodBs.size()<=0) continue;
    h_selection->Fill(4.,perEveW);

    
    // ------------------------------------------------------------------
    for (u_int iB=0; iB<goodBs.size(); iB++) {

      // B-infos for further cuts offline
      int thisB = goodBs[iB];
      float thisBmass   = BToKEE_fit_mass[thisB];
      float thisBpt     = BToKEE_fit_pt[thisB];
      float thisBcos    = BToKEE_fit_cos2D[thisB];
      float thisBsvprob = BToKEE_svprob[thisB];
      float thisBxysig  = BToKEE_l_xy[thisB]/BToKEE_l_xy_unc[thisB];
      float thisBi3dsig = BToKEE_k_svip3d[thisB]/BToKEE_k_svip3d_err[thisB]; 

      // MC truth
      bool isThisAMcB = -1;
      if (sampleID>0) isThisAMcB = isMcBOtto(thisB);            

      // Index
      int k_idx    = BToKEE_kIdx[thisB];  
      int ele1_idx = BToKEE_l1Idx[thisB];
      int ele2_idx = BToKEE_l2Idx[thisB];

      // Kine
      float k_pt     = BToKEE_fit_k_pt[thisB];   
      float ele1_pt  = BToKEE_fit_l1_pt[thisB];
      float ele2_pt  = BToKEE_fit_l2_pt[thisB];
      float k_eta    = BToKEE_fit_k_eta[thisB];   
      float ele1_eta = BToKEE_fit_l1_eta[thisB];
      float ele2_eta = BToKEE_fit_l2_eta[thisB];
      float k_phi    = BToKEE_fit_k_phi[thisB];   
      float ele1_phi = BToKEE_fit_l1_phi[thisB];
      float ele2_phi = BToKEE_fit_l2_phi[thisB];

      // Kine raw
      float k_pt_raw     = ProbeTracks_pt[k_idx];   
      float ele1_pt_raw  = Electron_pt[ele1_idx];
      float ele2_pt_raw  = Electron_pt[ele2_idx];  
      float k_eta_raw    = ProbeTracks_eta[k_idx];   
      float ele1_eta_raw = Electron_eta[ele1_idx];
      float ele2_eta_raw = Electron_eta[ele2_idx];  
      float k_phi_raw    = ProbeTracks_phi[k_idx];   
      float ele1_phi_raw = Electron_phi[ele1_idx];
      float ele2_phi_raw = Electron_phi[ele2_idx];  


      // Ele-ID and isolation, with swap if needed
      float L1id=-1000.;
      float L2id=-1000.;
      float L1iso04Rel=-1000;
      float L2iso04Rel=-1000;
      if( Electron_isPF[ele1_idx]==1 && Electron_isPF[ele2_idx]==1 ) { 
	L1id=Electron_pfmvaId[ele1_idx];
	L2id=Electron_pfmvaId[ele2_idx];
	L1iso04Rel=BToKEE_l1_iso04[thisB]/BToKEE_fit_l1_pt[thisB];
	L2iso04Rel=BToKEE_l2_iso04[thisB]/BToKEE_fit_l2_pt[thisB];
      }
      if( Electron_isPF[ele1_idx]==1 && Electron_isPF[ele2_idx]==0 ) { 
	L1id=Electron_pfmvaId[ele1_idx];
	L2id=Electron_mvaId[ele2_idx];
	L1iso04Rel=BToKEE_l1_iso04[thisB]/BToKEE_fit_l1_pt[thisB];
	L2iso04Rel=BToKEE_l2_iso04[thisB]/BToKEE_fit_l2_pt[thisB];
      }
      if( Electron_isPF[ele1_idx]==0 && Electron_isPF[ele2_idx]==1 ) { 
	L1id=Electron_pfmvaId[ele2_idx];
	L2id=Electron_mvaId[ele1_idx];
	L1iso04Rel=BToKEE_l2_iso04[thisB]/BToKEE_fit_l2_pt[thisB];
	L2iso04Rel=BToKEE_l1_iso04[thisB]/BToKEE_fit_l1_pt[thisB];
      }
      
      // swap ele1/ele2 if needed: 
      // PFPF: L1=leading, L2=subleading  => no swap needed
      // PFLP, L1=PF, L2=LP               => no swap needed  
      // PFLP, L1=LP, L2=PF               => !! swap needed !!  
      // This swap will be automatically propagated to all variables computed in the following
      if( Electron_isPF[ele1_idx]==0 && Electron_isPF[ele2_idx]==1 ) { 
	ele1_idx = BToKEE_l2Idx[thisB];
	ele2_idx = BToKEE_l1Idx[thisB];
	ele1_pt  = BToKEE_fit_l2_pt[thisB];
	ele2_pt  = BToKEE_fit_l1_pt[thisB];
	ele1_eta = BToKEE_fit_l2_eta[thisB];
	ele2_eta = BToKEE_fit_l1_eta[thisB];
	ele1_phi = BToKEE_fit_l2_phi[thisB];
	ele2_phi = BToKEE_fit_l1_phi[thisB];
	// 
	ele1_pt_raw  = Electron_pt[ele1_idx];   // swap done at idx level
	ele2_pt_raw  = Electron_pt[ele2_idx];  
	ele1_eta_raw = Electron_eta[ele1_idx];
	ele2_eta_raw = Electron_eta[ele2_idx];  
	ele1_phi_raw = Electron_phi[ele1_idx];
	ele2_phi_raw = Electron_phi[ele2_idx];  
      }

      // Vectors
      TLorentzVector kTLV(0,0,0,0);
      kTLV.SetPtEtaPhiM(k_pt,k_eta,k_phi,0.493);
      TLorentzVector ele1TLV(0,0,0,0);
      ele1TLV.SetPtEtaPhiM(ele1_pt,ele1_eta,ele1_phi,0.000511);
      TLorentzVector ele2TLV(0,0,0,0);
      ele2TLV.SetPtEtaPhiM(ele2_pt,ele2_eta,ele2_phi,0.000511);

      // Vertices
      float kvx = ProbeTracks_vx[k_idx];
      float kvy = ProbeTracks_vy[k_idx];
      float kvz = ProbeTracks_vz[k_idx];
      
      float l1vx = Electron_vx[ele1_idx];
      float l1vy = Electron_vy[ele1_idx];
      float l1vz = Electron_vz[ele1_idx];
      
      float l2vx = Electron_vx[ele2_idx];
      float l2vy = Electron_vy[ele2_idx];
      float l2vz = Electron_vz[ele2_idx];

      // check if it exists another triplet with the same electrons and check the 
      // highest 2ndary vertex probability for those
      int iPbest=0;
      float prob4trk=0; 
      for (int iP=0; iP<nBToKEE; iP++){
	if(iP!=thisB){
	  if( (BToKEE_l1Idx[iP]==ele1_idx && BToKEE_l2Idx[iP]==ele2_idx) || (BToKEE_l1Idx[iP]==ele2_idx && BToKEE_l2Idx[iP]==ele1_idx) || 
	      (ele1_pt==Electron_pt[BToKEE_l1Idx[iP]] && Electron_isPF[ele1_idx]==Electron_isPF[BToKEE_l1Idx[iP]] 
	       && ele2_pt== Electron_pt[BToKEE_l2Idx[iP]] && Electron_isPF[ele2_idx]==Electron_isPF[BToKEE_l2Idx[iP]]) || 
	      (ele2_pt==Electron_pt[BToKEE_l1Idx[iP]] && Electron_isPF[ele2_idx]==Electron_isPF[BToKEE_l1Idx[iP]] &&  
	       ele1_pt==Electron_pt[BToKEE_l2Idx[iP]] && Electron_isPF[ele1_idx]==Electron_isPF[BToKEE_l2Idx[iP]] ) ){
	    if(BToKEE_svprob[iP]>prob4trk){
	      iPbest=iP;
	      prob4trk=BToKEE_svprob[iP];
	    }
	  }
	}
      }
      
      int iPbestKStar=0;
      float prob4trkKStar=0;
      if(sampleID>1 ){  // there is the KStar part in the nanoAOD
	for (int iP=0; iP<nBToKsEE; iP++){
	  if( 
	     (( BToKsEE_l1_idx[iP]==ele1_idx &&BToKsEE_l2_idx[iP]==ele2_idx)||
	      (BToKsEE_l1_idx[iP]==ele2_idx &&BToKsEE_l2_idx[iP]==ele1_idx)|| 
	      (ele1_pt== Electron_pt[BToKsEE_l1_idx[iP]] && 
	       Electron_isPF[ele1_idx]==Electron_isPF[BToKsEE_l1_idx[iP]] 
	       && ele2_pt== Electron_pt[BToKsEE_l2_idx[iP]] && 
	       Electron_isPF[ele2_idx]==Electron_isPF[BToKsEE_l2_idx[iP]])|| 
	      (ele2_pt== Electron_pt[BToKsEE_l1_idx[iP]] && 
	       Electron_isPF[ele2_idx]==Electron_isPF[BToKsEE_l1_idx[iP]] &&  
	       ele1_pt== Electron_pt[BToKsEE_l2_idx[iP]] && 
	       Electron_isPF[ele1_idx]==Electron_isPF[BToKsEE_l2_idx[iP]] )) && 
	     (((BToKEE_kIdx[thisB]==Kstar_trk1_idx[BToKsEE_kstar_idx[iP]] || 
		BToKEE_kIdx[thisB]==Kstar_trk2_idx[BToKsEE_kstar_idx[iP]])) || 
	      (ProbeTracks_pt[BToKEE_kIdx[thisB]]==ProbeTracks_pt[Kstar_trk1_idx[BToKsEE_kstar_idx[iP]]] 
	       || ProbeTracks_pt[BToKEE_kIdx[thisB]]==ProbeTracks_pt[Kstar_trk2_idx[BToKsEE_kstar_idx[iP]]]))){
	    if(BToKsEE_svprob[iP]>prob4trkKStar){
	      iPbestKStar=iP;
	      prob4trkKStar=BToKsEE_svprob[iP];
	    }
	  }
	}
      }

      // check if not close to a muon 
      float dr1m=100.;
      float dr2m=100.;
      if(nMuon>0){
	for (u_int iM=0; iM<nMuon; iM++) {
	  float dr1=DeltaR(ele1_eta, ele1_phi, Muon_eta[iM], Muon_phi[iM]);
	  if(dr1<dr1m) dr1m=dr1; 	  
	  float dr2=DeltaR(ele2_eta, ele2_phi, Muon_eta[iM], Muon_phi[iM]);
	  if(dr2<dr2m) dr2m=dr2; 	  
	}
      }
    
      float x=goodBorder[iB];
      bPtOrder.push_back(x);
      float y=goodBorder2D[iB];
      b2DOrder.push_back(y);
      float z=goodBorderXY[iB];
      bXYOrder.push_back(z);

      p4Trk.push_back(prob4trk);
      p4TrkKStar.push_back(prob4trkKStar);

      fit_Bmass.push_back(thisBmass);
      fit_Bpt.push_back(thisBpt);
      fit_Bcos2D.push_back(thisBcos);
      fit_Bsvprob.push_back(thisBsvprob);
      fit_Bxysig.push_back(thisBxysig);
      fit_Bi3dsig.push_back(thisBi3dsig);
      if(sampleID>0) bmatchMC.push_back(isThisAMcB);
      // 
      tag_pt.push_back(ele1_pt);
      tag_eta.push_back(ele1_eta);
      tag_phi.push_back(ele1_phi);
      tag_isPF.push_back(Electron_isPF[ele1_idx]);           
      tag_isPFOverlap.push_back(Electron_isPFoverlap[ele1_idx]); 
      tag_isLowPt.push_back(Electron_isLowPt[ele1_idx]);
      tag_id.push_back(L1id);
      tag_drm.push_back(dr1m);
      //
      probe_pt.push_back(ele2_pt);
      probe_eta.push_back(ele2_eta);
      probe_phi.push_back(ele2_phi);
      probe_isPF.push_back(Electron_isPF[ele2_idx]);
      probe_isPFOverlap.push_back(Electron_isPFoverlap[ele2_idx]); 
      probe_isLowPt.push_back(Electron_isLowPt[ele2_idx]); 
      probe_id.push_back(L2id);
      probe_drm.push_back(dr2m);
      // 
      K_pt.push_back(k_pt);
      K_eta.push_back(k_eta);
      K_phi.push_back(k_phi);
      // 
      mll_fullfit.push_back(BToKEE_mll_fullfit[thisB]);
      mll_raw.push_back(BToKEE_mll_raw[thisB]);
      
      // Vertex distances
      float l1kDz = fabs(kvz-l1vz);
      float l2kDz = fabs(kvz-l2vz);
      float LKdz  = l1kDz;
      if (l2kDz<LKdz) LKdz = l2kDz;
	
      // DeltaR
      float l1kDr = DeltaR(ele1_eta, ele1_phi, k_eta, k_phi);
      float l2kDr = DeltaR(ele2_eta, ele2_phi, k_eta, k_phi);
      float LKdr  = l1kDr;
      if (l2kDr<LKdr) LKdr = l2kDr;
      float L1L2dr = DeltaR(ele1_eta, ele1_phi, ele2_eta, ele2_phi);

      float l1kDr_raw = DeltaR(ele1_eta_raw, ele1_phi_raw, k_eta_raw, k_phi_raw);
      float l2kDr_raw = DeltaR(ele2_eta_raw, ele2_phi_raw, k_eta_raw, k_phi_raw);
      float LKdr_raw  = l1kDr_raw;
      if (l2kDr_raw<LKdr_raw) LKdr_raw = l2kDr_raw;
      float L1L2dr_raw = DeltaR(ele1_eta_raw, ele1_phi_raw, ele2_eta_raw, ele2_phi_raw);
	
      // Isolation
      float Kiso = BToKEE_k_iso04[thisB];

      // Vtx
      float vx_b=BToKEE_vtx_x[thisB];
      float vy_b=BToKEE_vtx_y[thisB];
      float vz_b=BToKEE_vtx_z[thisB];
      TVector3 b_vtx;
      b_vtx.SetXYZ(vx_b,vy_b,vz_b);
  
      // BTrkdxy2
      float distance_b_trk1=10;         
      float distance_b_trk2=10;         
      float mean_b_trk=0;      
      int imean_b_trk=0 ;
      float BTrkdxy2=0;
      for (int itrk=0; itrk<nProbeTracks; itrk++){
	if(itrk==k_idx) continue;
	float vx_trk=ProbeTracks_vx[itrk];
	float vy_trk=ProbeTracks_vy[itrk];
	float vz_trk=ProbeTracks_vz[itrk];
	float trk_eta=ProbeTracks_eta[itrk]; 
	float trk_phi=ProbeTracks_phi[itrk];
	
	if( (DeltaR(ele1_eta,ele1_phi,trk_eta,trk_phi)<0.03) || (DeltaR(ele2_eta,ele2_phi,trk_eta,trk_phi)<0.3) || ( DeltaR(k_eta,k_phi,trk_eta,trk_phi)<0.03)) continue;
	
	float distance_b_trk=TMath::Sqrt( pow((vx_b-vx_trk),2) + pow((vy_b-vy_trk),2) ) ;
	if (DeltaR(BToKEE_fit_eta[thisB],BToKEE_fit_phi[thisB],trk_eta,trk_phi)<0.4){
	  mean_b_trk=mean_b_trk+TMath::Sqrt( pow((vx_b-vx_trk),2) + pow((vy_b-vy_trk),2) );
	  imean_b_trk=imean_b_trk+1;
	}

	if(distance_b_trk<distance_b_trk2 && distance_b_trk<distance_b_trk1){
	  distance_b_trk2=distance_b_trk1;
	  distance_b_trk1=distance_b_trk;
	} else if(distance_b_trk<distance_b_trk2) {
	  distance_b_trk2=distance_b_trk;
	} 
      }      
      if(imean_b_trk>0) {
	mean_b_trk=mean_b_trk/imean_b_trk;
      } else {
	mean_b_trk=10.;
      }
      BTrkdxy2=distance_b_trk2;
            
      float BBDPhi=DeltaPhi(BToKEE_fit_phi[thisB],sum_track_vec.Phi());

      // p-asym
      TLorentzVector dll=ele1TLV+ele2TLV;
      TVector3 diele_p3 = dll.Vect();
      TVector3 k_p3 = kTLV.Vect();
      TVector3 pv2sv_p3(PV_x-BToKEE_vtx_x[thisB], PV_y-BToKEE_vtx_y[thisB], PV_z-BToKEE_vtx_z[thisB]);
      float BToKEE_ptAsym = ( (diele_p3.Cross(pv2sv_p3)).Mag() - (k_p3.Cross(pv2sv_p3)).Mag() ) 
			  / ( (diele_p3.Cross(pv2sv_p3)).Mag() + (k_p3.Cross(pv2sv_p3)).Mag() ); 

      float K_iso04_rel = BToKEE_k_iso04[thisB] /BToKEE_fit_k_pt[thisB];

      // Tree : save variables for old George's BDT and Common BDTs
      LKdz_vec.push_back(LKdz);
      L1L2dr_vec.push_back(L1L2dr);
      L1L2dr_raw_vec.push_back(L1L2dr_raw);
      LKdr_vec.push_back(LKdr);
      LKdr_raw_vec.push_back(LKdr_raw);
      BBDPhi_vec.push_back(BBDPhi);
      BTrkdxy2_vec.push_back(BTrkdxy2);
      ptAsym_vec.push_back(BToKEE_ptAsym);     
      k_iso04_rel_vec.push_back(K_iso04_rel);
      tag_iso04_rel_vec.push_back(L1iso04Rel);
      probe_iso04_rel_vec.push_back(L2iso04Rel);


      // -----------------------------------------------------------
      // Common analysis BDT
      float L1_normpt = ele1_pt/BToKEE_fit_mass[thisB];
      float L2_normpt = ele2_pt/BToKEE_fit_mass[thisB];
      float K_normpt  = k_pt/BToKEE_fit_mass[thisB];

      /*
      float scoreBdtCommon=0.;
      if(Electron_isPF[ele1_idx] && Electron_isPF[ele2_idx]){
	std::vector<float> vecBdtCommon = { thisBsvprob, thisBxysig, thisBcos, L1_normpt, L2_normpt, K_normpt, LKdz, L1L2dr_raw, LKdr_raw, L1id, L2id, L2iso04Rel, K_iso04_rel, BBDPhi, BTrkdxy2, BToKEE_ptAsym, thisBi3dsig };
	scoreBdtCommon = fastForestCommonPFPF(vecBdtCommon.data());
	analysisBdtC.push_back(scoreBdtCommon);

	float scoreBdtCommonWithSyst = scoreBdtCommon;
	if (dosyst_) scoreBdtCommonWithSyst = GetAnBdtWeight(scoreBdtCommon, 1);   
	analysisBdtCWithSyst.push_back(scoreBdtCommonWithSyst);  

      } else {
	std::vector<float> vecBdtCommon = { thisBsvprob, thisBxysig, thisBcos, L1_normpt, L2_normpt, K_normpt, LKdz, L1L2dr_raw, LKdr_raw, L1id, L2id, L2iso04Rel, K_iso04_rel, BBDPhi, BTrkdxy2, BToKEE_ptAsym, thisBi3dsig };
	scoreBdtCommon = fastForestCommonPFLP(vecBdtCommon.data());
	analysisBdtC.push_back(scoreBdtCommon);

	float scoreBdtCommonWithSyst = scoreBdtCommon;
	if (dosyst_) scoreBdtCommonWithSyst = GetAnBdtWeight(scoreBdtCommon, 0);
	analysisBdtCWithSyst.push_back(scoreBdtCommonWithSyst);  
      }
      */
      
      // ----------------------------------------------------------------
      // Variables for old analysis BDT by Otto, PFPF
      float BToKEE_l1_dxy_sig=(Electron_dxy[ele1_idx]) /Electron_dxyErr[ele1_idx];
      float BToKEE_l2_dxy_sig=(Electron_dxy[ele2_idx]) /Electron_dxyErr[ele2_idx];
      float BToKEE_k_dxy_sig=ProbeTracks_dxyS[k_idx];
      float BToKEE_l1_dzTrg=Electron_dzTrg[ele1_idx];
      float BToKEE_l2_dzTrg=Electron_dzTrg[ele2_idx];
      float BToKEE_k_dzTrg=ProbeTracks_dzTrg[k_idx];
      // 
      float BToKEE_eleDR= DeltaR(ele1_eta,ele1_phi,ele2_eta,ele2_phi);
      float BToKEE_llkDR=dll.DeltaR(kTLV);
      float BToKEE_b_iso04_rel =BToKEE_b_iso04[thisB]/BToKEE_fit_pt[thisB];
      //
      // D-related variables
      float BToKEE_k_DCASig = ProbeTracks_DCASig[k_idx];
      // 
      TLorentzVector l1_pihypoTLV(0,0,0,0);
      l1_pihypoTLV.SetPtEtaPhiM(ele1_pt,ele1_eta,ele1_phi,0.139);
      TLorentzVector l2_pihypoTLV(0,0,0,0);
      l2_pihypoTLV.SetPtEtaPhiM(ele2_pt,ele2_eta,ele2_phi,0.139);
      TLorentzVector l1_khypoTLV(0,0,0,0);
      l1_khypoTLV.SetPtEtaPhiM(ele1_pt,ele1_eta,ele1_phi,0.493);
      TLorentzVector l2_khypoTLV(0,0,0,0);
      l2_khypoTLV.SetPtEtaPhiM(ele2_pt,ele2_eta,ele2_phi,0.493);
      TLorentzVector k_pihypoTLV(0,0,0,0);
      k_pihypoTLV.SetPtEtaPhiM(k_pt,k_eta,k_phi,0.139);
      float BToKEE_Dmass_l1=(l1_pihypoTLV + kTLV).M();
      float BToKEE_Dmass_l2=(l2_pihypoTLV + kTLV).M();
      float BToKEE_Dmass_flip_l1=(l1_khypoTLV + k_pihypoTLV).M();
      float BToKEE_Dmass_flip_l2=(l2_khypoTLV + k_pihypoTLV).M();
      float BToKEE_Dmass=0.;              
      if(ProbeTracks_charge[k_idx]*Electron_charge[ele1_idx]<0) BToKEE_Dmass=BToKEE_Dmass_l1;
      if(ProbeTracks_charge[k_idx]*Electron_charge[ele2_idx]<0) BToKEE_Dmass=BToKEE_Dmass_l2;
      float BToKEE_Dmass_flip=0;          
      if(ProbeTracks_charge[k_idx]*Electron_charge[ele1_idx]<0) BToKEE_Dmass_flip=BToKEE_Dmass_flip_l1;       
      if(ProbeTracks_charge[k_idx]*Electron_charge[ele2_idx]<0) BToKEE_Dmass_flip=BToKEE_Dmass_flip_l2;
      float BToKEE_Dmass_ll=(l1_khypoTLV + l2_pihypoTLV).M();
      float BToKEE_Dmass_ll_flip=(l1_pihypoTLV + l2_khypoTLV).M();      

      // This is the variable to be used for anti-D0 cut
      float k_opp_l_mass_n1 = 0.;
      if (ProbeTracks_charge[k_idx]*Electron_charge[ele1_idx]<0) {
	k_opp_l_mass_n1 = BToKEE_Dmass_flip_l1;
      } else {
	k_opp_l_mass_n1 = BToKEE_Dmass_flip_l2;
      }
      float k_opp_l_mass_n2 = 0.;
      if (ProbeTracks_charge[k_idx]*Electron_charge[ele1_idx]<0) {
	k_opp_l_mass_n2 = BToKEE_Dmass_l1;
      } else {
	k_opp_l_mass_n2 = BToKEE_Dmass_l2;
      }
      float k_opp_l_mass = k_opp_l_mass_n1;
      if (k_opp_l_mass_n2<k_opp_l_mass_n1) k_opp_l_mass = k_opp_l_mass_n2;


      float m_muk_hadron_mass=-99.;  
      float m_muk_muon_mass=-99.;           // this is George/Otto's KMumassJpsi
      TLorentzVector trgmuon_khypoTLV(0,0,0,0);
      if(nTriggerMuon>0){
	trgmuon_khypoTLV.SetPtEtaPhiM(Muon_pt[idTrgMu], Muon_eta[idTrgMu],Muon_phi[idTrgMu],0.493);
	float kcharge  = ProbeTracks_charge[k_idx];  
	float mucharge = Muon_charge[idTrgMu];
	m_muk_hadron_mass=0.;
	m_muk_muon_mass=0.;
	if (kcharge!=mucharge) {
	  m_muk_hadron_mass=(trgmuon_khypoTLV+k_pihypoTLV).M();

	  TLorentzVector k_muhypoTLV(0,0,0,0);
	  k_muhypoTLV.SetPtEtaPhiM(k_pihypoTLV.Pt(),k_pihypoTLV.Eta(),k_pihypoTLV.Phi(),0.105);
	  TLorentzVector trgmuon_muhypoTLV(0,0,0,0);
	  trgmuon_muhypoTLV.SetPtEtaPhiM(trgmuon_khypoTLV.Pt(),trgmuon_khypoTLV.Eta(),trgmuon_khypoTLV.Phi(),0.105);
	  m_muk_muon_mass=(k_muhypoTLV+trgmuon_muhypoTLV).M();
	}
      }


      // Tree: save further variables for old Otto's BDT
      eleDR_vec.push_back(BToKEE_eleDR);
      k_DCASig_vec.push_back(BToKEE_k_DCASig);
      k_svip3d_vec.push_back(BToKEE_k_svip3d[thisB]);
      // Not used anymore - start -------------
      b_iso04_rel_vec.push_back(BToKEE_b_iso04_rel);
      k_dzTrg_vec.push_back(BToKEE_k_dzTrg);
      k_svip2d_vec.push_back(BToKEE_k_svip2d[thisB]);
      k_dxy_sig_vec.push_back(BToKEE_k_dxy_sig);
      tag_dxy_sig_vec.push_back(BToKEE_l1_dxy_sig);
      tag_dzTrg_vec.push_back(BToKEE_l1_dzTrg);
      probe_dxy_sig_vec.push_back(BToKEE_l2_dxy_sig);
      probe_dzTrg_vec.push_back(BToKEE_l2_dzTrg);
      llkDR_vec.push_back(BToKEE_llkDR);
      // Not used anymore - end -------------
      Dmass_vec.push_back(BToKEE_Dmass);
      Dmass_flip_vec.push_back(BToKEE_Dmass_flip);
      Dmass_ll_vec.push_back(BToKEE_Dmass_ll);
      Dmass_ll_flip_vec.push_back(BToKEE_Dmass_ll_flip);
      Dmass_mukhad_vec.push_back(m_muk_hadron_mass);  
      Dmass_mukmu_vec.push_back(m_muk_muon_mass);          
      KLmassD0_vec.push_back(k_opp_l_mass);
      
      if (sampleID>0) {     // MC      
	bool isTagAMcEleFromJPsi   = isMcEleFromJPsi(ele1_idx);
	bool isProbeAMcEleFromJPsi = isMcEleFromJPsi(ele2_idx);
	bool isTagAMcEle   = (Electron_genPartIdx[ele1_idx]>-0.5);
	bool isProbeAMcEle = (Electron_genPartIdx[ele2_idx]>-0.5);

	if(isTagAMcEle){
	  tag_ptMc.push_back(GenPart_pt[Electron_genPartIdx[ele1_idx]]);
	}else{
	  tag_ptMc.push_back(0.);
	}
	if(isProbeAMcEle){
	  probe_ptMc.push_back(GenPart_pt[Electron_genPartIdx[ele2_idx]]);
	}else{
	  probe_ptMc.push_back(0.);
	}
	tag_matchMc.push_back(isTagAMcEle);
	probe_matchMc.push_back(isProbeAMcEle);
	tag_matchMcFromJPsi.push_back(isTagAMcEleFromJPsi);
	probe_matchMcFromJPsi.push_back(isProbeAMcEleFromJPsi);
	
      } else {
	tag_matchMcFromJPsi.push_back(0);  
	probe_matchMcFromJPsi.push_back(0);  
	tag_matchMc.push_back(0);  
	probe_matchMc.push_back(0);  
	tag_ptMc.push_back(0.);
	probe_ptMc.push_back(0.);
      }
      
    } // Loop over good Bs
      

    // At least one tag and one probe
    selectedPairsSize = tag_pt.size();
    if (selectedPairsSize<=0) continue;
    h_selection->Fill(5.,perEveW);

    
    // ------------------------------------------------------
    // Minimal selection to reduce the ntuples size - chiara
    // if (hlt9==0) continue;
    // ------------------------------------------------------

    // Filling the output tree
    outTree_->Fill();
    
    // Cleaning all vectors used for the selection
    goodBs.clear();
    
    // Cleaning all vectors used for the output tree, ready for a new entry
    bPtOrder.clear();
    b2DOrder.clear();
    bXYOrder.clear();
    p4Trk.clear();
    p4TrkKStar.clear();
    tag_pt.clear();  
    tag_eta.clear();  
    tag_phi.clear();  
    tag_isPF.clear();  
    tag_isPFOverlap.clear();
    tag_isLowPt.clear();  
    tag_id.clear();
    if(sampleID>0){
      tag_ptMc.clear();
      probe_ptMc.clear();
      tag_matchMcFromJPsi.clear();  
      tag_matchMc.clear();
      bmatchMC.clear();
      probe_matchMcFromJPsi.clear();  
      probe_matchMc.clear();  
    }  
    mll_fullfit.clear();
    mll_raw.clear();

    fit_Bmass.clear();
    fit_Bpt.clear();
    fit_Bcos2D.clear();
    fit_Bsvprob.clear();    
    fit_Bxysig.clear();    
    fit_Bi3dsig.clear();

    K_pt.clear();
    K_eta.clear();
    K_phi.clear();
    probe_pt.clear();  
    probe_drm.clear();  
    tag_drm.clear();  
    probe_eta.clear();  
    probe_phi.clear();  
    probe_isPF.clear();  
    probe_isPFOverlap.clear();
    probe_isLowPt.clear();  
    probe_id.clear();  

    //analysisBdtC.clear();
    //analysisBdtCWithSyst.clear(); 

    LKdz_vec.clear();
    L1L2dr_vec.clear();
    L1L2dr_raw_vec.clear();
    LKdr_vec.clear();
    LKdr_raw_vec.clear();
    BBDPhi_vec.clear();
    BTrkdxy2_vec.clear(); 

    Dmass_vec.clear();
    Dmass_flip_vec.clear();
    Dmass_ll_vec.clear();
    Dmass_ll_flip_vec.clear();
    Dmass_mukhad_vec.clear();
    Dmass_mukmu_vec.clear();
    KLmassD0_vec.clear();

    b_iso04_rel_vec.clear();
    eleDR_vec.clear();
    k_DCASig_vec.clear();
    k_dzTrg_vec.clear();
    k_dxy_sig_vec.clear();
    k_iso04_rel_vec.clear();
    k_svip2d_vec.clear();
    k_svip3d_vec.clear();
    tag_dxy_sig_vec.clear();
    tag_dzTrg_vec.clear();
    tag_iso04_rel_vec.clear();
    probe_dxy_sig_vec.clear();
    probe_dzTrg_vec.clear();
    probe_iso04_rel_vec.clear();
    llkDR_vec.clear();
    ptAsym_vec.clear();
  }
}

void SkimmerWithKStar::SetNvtxWeights(std::string nvtxWeightFile) {

  if (nvtxWeightFile == "") {
    std::cout << "you need a weights file to use this function" << std::endl;
    return;
  }
  std::cout << "PU REWEIGHTING Based on #vertices:: Using file " << nvtxWeightFile << std::endl;
  
  TFile *f_nvtx = new TFile(nvtxWeightFile.c_str(),"READ");
  f_nvtx->cd();
  
  TH1F *nvtxweights = 0;
  TH1F *mc_nvtx = 0;
  mc_nvtx     = (TH1F*) f_nvtx->Get("mcNvtx");
  nvtxweights = (TH1F*) f_nvtx->Get("weights");
  
  if (!nvtxweights || !mc_nvtx) {
    std::cout << "weights histograms not found in file " << nvtxWeightFile << std::endl;
    return;
  }
  TH1F* weightedNvtx= (TH1F*)mc_nvtx->Clone("weightedNvtx");
  weightedNvtx->Multiply(nvtxweights);
  
  // Rescaling weights in order to preserve same integral of events     
  TH1F* weights = (TH1F*)nvtxweights->Clone("rescaledWeights");
  weights->Scale( mc_nvtx->Integral() / weightedNvtx->Integral() );
  
  float sumNvtxweights=0.;
  for (int i = 0; i<nvtxweights->GetNbinsX(); i++) {
    float weight=1.;
    weight=weights->GetBinContent(i+1);
    sumNvtxweights+=weight;
    nvtxweights_.push_back(weight);
    float lowedge=weights->GetBinLowEdge(i+1);
    nvtxlowedge_.push_back(lowedge);
  }
}

float SkimmerWithKStar::GetNvtxWeight(float nvtx) {

  int thesize   = nvtxlowedge_.size();
  int thesizem1 = nvtxlowedge_.size()-1;
  float weight=1;

  if (sampleID!=0 && thesize>0 && donvtxreweight_) {
    for (int i = 0; i<thesizem1; i++) {   
      if (nvtxlowedge_[i]<=nvtx && nvtxlowedge_[i+1]>nvtx) { 
	weight = nvtxweights_[i];
      }
    }
    if (nvtxlowedge_[thesizem1]<=nvtx) { 
      weight = nvtxweights_[thesizem1];
    }
  }

  return weight;
}

bool SkimmerWithKStar::isMcB( int theB ) {
  
  // taking index
  int ele1_idx = BToKEE_l1Idx[theB];
  int ele2_idx = BToKEE_l2Idx[theB];
  int k_idx    = BToKEE_kIdx[theB];

  // Gen tree
  int k_genPartIdx      = ProbeTracks_genPartIdx[k_idx];  
  int k_genMotherIdx    = GenPart_genPartIdxMother[k_genPartIdx];
  int k_genGMotherIdx   = GenPart_genPartIdxMother[k_genMotherIdx];
  int k_genPdgId        = GenPart_pdgId[k_genPartIdx];
  int k_genMotherPdgId  = GenPart_pdgId[k_genMotherIdx];
  int k_genGMotherPdgId = GenPart_pdgId[k_genGMotherIdx];

  int ele1_genPartIdx      = Electron_genPartIdx[ele1_idx];  
  int ele1_genMotherIdx    = GenPart_genPartIdxMother[ele1_genPartIdx];
  int ele1_genGMotherIdx   = GenPart_genPartIdxMother[ele1_genMotherIdx];
  int ele1_genPdgId        = GenPart_pdgId[ele1_genPartIdx];
  int ele1_genMotherPdgId  = GenPart_pdgId[ele1_genMotherIdx];
  int ele1_genGMotherPdgId = GenPart_pdgId[ele1_genGMotherIdx];

  int ele2_genPartIdx      = Electron_genPartIdx[ele2_idx];  
  int ele2_genMotherIdx    = GenPart_genPartIdxMother[ele2_genPartIdx];
  int ele2_genGMotherIdx   = GenPart_genPartIdxMother[ele2_genMotherIdx];
  int ele2_genPdgId        = GenPart_pdgId[ele2_genPartIdx];
  int ele2_genMotherPdgId  = GenPart_pdgId[ele2_genMotherIdx];
  int ele2_genGMotherPdgId = GenPart_pdgId[ele2_genGMotherIdx];

  // B -> K J/psi(ll) at gen level
  bool okMatch = (ele1_genPartIdx>-0.5 && ele2_genPartIdx>-0.5 && k_genPartIdx>-0.5);
  bool RK_res1 = abs(ele1_genMotherPdgId)==443 && abs(k_genMotherPdgId)==521;
  bool RK_res2 = (ele1_genMotherPdgId==ele2_genMotherPdgId) && (k_genMotherPdgId==ele1_genGMotherPdgId) && (k_genMotherPdgId==ele2_genGMotherPdgId);

  bool RK_res3 = abs(ele1_genMotherPdgId)==521 && abs(k_genMotherPdgId)==521;
  bool RK_res4 = (ele1_genMotherPdgId==ele2_genMotherPdgId) && (k_genMotherPdgId==ele1_genMotherPdgId) && (k_genMotherPdgId==ele2_genMotherPdgId);

  bool RK_res = okMatch  && (( RK_res1 && RK_res2) || (RK_res3 &&RK_res4)) ;

  return RK_res;
}

bool SkimmerWithKStar::isMcBOtto( int theB ) {
  
  // taking index
  int ele1_idx = BToKEE_l1Idx[theB];
  int ele2_idx = BToKEE_l2Idx[theB];
  int k_idx    = BToKEE_kIdx[theB];

  // Gen tree
  int k_genPartIdx      = ProbeTracks_genPartIdx[k_idx];  
  int k_genMotherIdx    = GenPart_genPartIdxMother[k_genPartIdx];
  int k_genGMotherIdx   = GenPart_genPartIdxMother[k_genMotherIdx];
  int k_genPdgId        = GenPart_pdgId[k_genPartIdx];
  int k_genMotherPdgId  = GenPart_pdgId[k_genMotherIdx];
  int k_genGMotherPdgId = GenPart_pdgId[k_genGMotherIdx];

  int ele1_genPartIdx      = Electron_genPartIdx[ele1_idx];  
  int ele1_genMotherIdx    = GenPart_genPartIdxMother[ele1_genPartIdx];
  int ele1_genGMotherIdx   = GenPart_genPartIdxMother[ele1_genMotherIdx];
  int ele1_genPdgId        = GenPart_pdgId[ele1_genPartIdx];
  int ele1_genMotherPdgId  = GenPart_pdgId[ele1_genMotherIdx];
  int ele1_genGMotherPdgId = GenPart_pdgId[ele1_genGMotherIdx];

  int ele2_genPartIdx      = Electron_genPartIdx[ele2_idx];  
  int ele2_genMotherIdx    = GenPart_genPartIdxMother[ele2_genPartIdx];
  int ele2_genGMotherIdx   = GenPart_genPartIdxMother[ele2_genMotherIdx];
  int ele2_genPdgId        = GenPart_pdgId[ele2_genPartIdx];
  int ele2_genMotherPdgId  = GenPart_pdgId[ele2_genMotherIdx];
  int ele2_genGMotherPdgId = GenPart_pdgId[ele2_genGMotherIdx];

  // Reco TVector3
  float reco_ele1_pt  = BToKEE_fit_l1_pt[theB];
  float reco_ele2_pt  = BToKEE_fit_l2_pt[theB];
  float reco_k_pt     = BToKEE_fit_k_pt[theB];
  float reco_ele1_eta = BToKEE_fit_l1_eta[theB];
  float reco_ele2_eta = BToKEE_fit_l2_eta[theB];
  float reco_k_eta    = BToKEE_fit_k_eta[theB];
  float reco_ele1_phi = BToKEE_fit_l1_phi[theB];
  float reco_ele2_phi = BToKEE_fit_l2_phi[theB];
  float reco_k_phi    = BToKEE_fit_k_phi[theB];
  // 
  TVector3 reco_ele1(0,0,0);
  TVector3 reco_ele2(0,0,0);
  TVector3 reco_k(0,0,0);
  reco_ele1.SetPtEtaPhi(reco_ele1_pt, reco_ele1_eta, reco_ele1_phi);
  reco_ele2.SetPtEtaPhi(reco_ele2_pt, reco_ele2_eta, reco_ele2_phi);
  reco_k.SetPtEtaPhi(reco_k_pt,reco_k_eta,reco_k_phi);

  // Gen TVector3  
  float gen_ele1_pt  = GenPart_pt[ele1_genPartIdx];
  float gen_ele2_pt  = GenPart_pt[ele2_genPartIdx];      
  float gen_k_pt     = GenPart_pt[k_genPartIdx];      
  float gen_ele1_eta = GenPart_eta[ele1_genPartIdx];
  float gen_ele2_eta = GenPart_eta[ele2_genPartIdx];      
  float gen_k_eta    = GenPart_eta[k_genPartIdx];      
  float gen_ele1_phi = GenPart_phi[ele1_genPartIdx];
  float gen_ele2_phi = GenPart_phi[ele2_genPartIdx];      
  float gen_k_phi    = GenPart_phi[k_genPartIdx];      
  // 
  TVector3 gen_ele1(0,0,0);
  TVector3 gen_ele2(0,0,0);
  TVector3 gen_k(0,0,0);
  gen_ele1.SetPtEtaPhi(gen_ele1_pt, gen_ele1_eta, gen_ele1_phi);
  gen_ele2.SetPtEtaPhi(gen_ele2_pt, gen_ele2_eta, gen_ele2_phi);
  gen_k.SetPtEtaPhi(gen_k_pt,gen_k_eta,gen_k_phi);

  // B -> K J/psi(ll) at gen level
  bool okMatch = (ele1_genPartIdx>-0.5 && ele2_genPartIdx>-0.5 && k_genPartIdx>-0.5);
  if (!okMatch) return okMatch;

  float dR1 = gen_ele1.DeltaR(reco_ele1);
  float dR2 = gen_ele2.DeltaR(reco_ele2);
  float dRK = gen_k.DeltaR(reco_k);
  bool RK_dR = (dR1<0.3) && (dR2<0.3) && (dRK<0.3);

  bool RK_res1 = abs(ele1_genMotherPdgId)==443 && abs(k_genMotherPdgId)==521;
  bool RK_res2 = (ele1_genMotherPdgId==ele2_genMotherPdgId) && (k_genMotherPdgId==ele1_genGMotherPdgId) && (k_genMotherPdgId==ele2_genGMotherPdgId);

  bool RK_res3 = abs(ele1_genMotherPdgId)==521 && abs(k_genMotherPdgId)==521;
  bool RK_res4 = (ele1_genMotherPdgId==ele2_genMotherPdgId) && (k_genMotherPdgId==ele1_genMotherPdgId) && (k_genMotherPdgId==ele2_genMotherPdgId);

  bool RK_res = okMatch && RK_dR && (( RK_res1 && RK_res2) || (RK_res3 &&RK_res4)) ;

  return RK_res;
}

bool SkimmerWithKStar::isMcEleFromJPsi( int ele_idx ) {

  // Gen tree
  int ele_genPartIdx      = Electron_genPartIdx[ele_idx];  
  int ele_genMotherIdx    = GenPart_genPartIdxMother[ele_genPartIdx];
  int ele_genGMotherIdx   = GenPart_genPartIdxMother[ele_genMotherIdx];
  int ele_genPdgId        = GenPart_pdgId[ele_genPartIdx];
  int ele_genMotherPdgId  = GenPart_pdgId[ele_genMotherIdx];
  int ele_genGMotherPdgId = GenPart_pdgId[ele_genGMotherIdx];

  // B -> K J/psi(ll) at gen level
  bool okMatch = (ele_genPartIdx>-0.5) && (abs(ele_genMotherPdgId)==443);

  return okMatch;
}

void SkimmerWithKStar::PrepareOutputs(std::string filename) 
{
  _datasetName=filename;

  std::string outname = _datasetName+".root";    
  cout << "output: " << outname << endl;
  outFile_ = new TFile(outname.c_str(),"RECREATE");

  bookOutputTree();
  bookOutputHistos();

  // loading weights for pileup if needed
  if (donvtxreweight_) SetNvtxWeights(nvtxWFileName_);

  // loading weights for systematics if needed
  if (dosyst_) SetAnBdtWeights();  
};


void SkimmerWithKStar::bookOutputTree() 
{
  outTree_ = new TTree("TaPtree", "TaPtree");
  
  cout << "Booking tree" << endl;
  
  outTree_->Branch("theRun", &theRun, "theRun/I");    
  outTree_->Branch("theEvent", &theEvent, "theEvent/I");    
  outTree_->Branch("nvtx", &nvtx, "nvtx/I");    
  outTree_->Branch("sampleID", &sampleID, "sampleID/I");    
  outTree_->Branch("rho", &rho, "rho/F");    

  outTree_->Branch("hlt7", &hlt7, "hlt7/I");
  outTree_->Branch("hlt9", &hlt9, "hlt9/I");
  outTree_->Branch("trg_muon_pt", &trg_muon_pt, "trg_muon_pt/F");

  outTree_->Branch("selectedBSize",  &selectedBSize,  "selectedBSize/I");   
  
  outTree_->Branch("tag_pt", "std::vector<float>", &tag_pt);  
  outTree_->Branch("tag_eta", "std::vector<float>", &tag_eta);  
  outTree_->Branch("tag_phi", "std::vector<float>", &tag_phi);  
  outTree_->Branch("tag_isPF", "std::vector<bool>", &tag_isPF);  
  outTree_->Branch("tag_isPFOverlap", "std::vector<bool>", &tag_isPFOverlap);  
  outTree_->Branch("tag_isLowPt", "std::vector<bool>", &tag_isLowPt);  
  outTree_->Branch("tag_id", "std::vector<float>", &tag_id);  
  outTree_->Branch("tag_drm", "std::vector<float>", &tag_drm);  

  if(sampleID>0){
    outTree_->Branch("tag_ptMc", "std::vector<float>", &tag_ptMc);
    outTree_->Branch("tag_matchMcFromJPsi", "std::vector<bool>", &tag_matchMcFromJPsi);  
    outTree_->Branch("tag_matchMc", "std::vector<bool>", &tag_matchMc);  
  }

  outTree_->Branch("mll_fullfit", "std::vector<float>", &mll_fullfit);  
  outTree_->Branch("mll_raw", "std::vector<float>", &mll_raw);  

  outTree_->Branch("fit_Bmass",   "std::vector<float>", &fit_Bmass);  
  outTree_->Branch("fit_Bpt",     "std::vector<float>", &fit_Bpt);  
  outTree_->Branch("fit_Bcos2D",  "std::vector<float>", &fit_Bcos2D);  
  outTree_->Branch("fit_Bsvprob", "std::vector<float>", &fit_Bsvprob);  
  outTree_->Branch("fit_Bxysig",  "std::vector<float>", &fit_Bxysig);  
  outTree_->Branch("fit_Bi3dsig", "std::vector<float>", &fit_Bi3dsig);

  outTree_->Branch("bmatchMC", "std::vector<bool>", &bmatchMC);  
  outTree_->Branch("K_pt", "std::vector<float>", &K_pt);  
  outTree_->Branch("K_eta", "std::vector<float>", &K_eta);  
  outTree_->Branch("K_phi", "std::vector<float>", &K_phi);  

  outTree_->Branch("probe_pt", "std::vector<float>", &probe_pt);  
  outTree_->Branch("probe_eta", "std::vector<float>", &probe_eta);  
  outTree_->Branch("probe_phi", "std::vector<float>", &probe_phi);  
  outTree_->Branch("probe_isPF", "std::vector<bool>", &probe_isPF);  
  outTree_->Branch("probe_isPFOverlap", "std::vector<bool>", &probe_isPFOverlap);  
  outTree_->Branch("probe_isLowPt", "std::vector<bool>", &probe_isLowPt);  
  outTree_->Branch("probe_id", "std::vector<float>", &probe_id);  

  if(sampleID>0){
    outTree_->Branch("probe_ptMc", "std::vector<float>", &probe_ptMc);
    outTree_->Branch("probe_matchMcFromJPsi", "std::vector<bool>", &probe_matchMcFromJPsi);  
    outTree_->Branch("probe_matchMc", "std::vector<bool>", &probe_matchMc);  
  }
  outTree_->Branch("probe_drm", "std::vector<float>", &probe_drm);  

  outTree_->Branch("bPtOrder", "std::vector<float>", &bPtOrder);  
  outTree_->Branch("b2DOrder", "std::vector<float>", &b2DOrder);  
  outTree_->Branch("bXYOrder", "std::vector<float>", &bXYOrder);  
  outTree_->Branch("p4Trk", "std::vector<float>", &p4Trk);  
  if(sampleID>1){
    outTree_->Branch("p4TrkKStar", "std::vector<float>", &p4TrkKStar);  
  }
  // outTree_->Branch("analysisBdtC", "std::vector<float>", &analysisBdtC);
  // outTree_->Branch("analysisBdtCWithSyst", "std::vector<float>", &analysisBdtCWithSyst);
  
  outTree_->Branch("LKdz", "std::vector<float>", &LKdz_vec);
  outTree_->Branch("L1L2dr", "std::vector<float>", &L1L2dr_vec);
  outTree_->Branch("L1L2dr_raw", "std::vector<float>", &L1L2dr_raw_vec);
  outTree_->Branch("LKdr", "std::vector<float>", &LKdr_vec);
  outTree_->Branch("LKdr_raw", "std::vector<float>", &LKdr_raw_vec);
  outTree_->Branch("BBDPhi", "std::vector<float>", &BBDPhi_vec);
  outTree_->Branch("BTrkdxy2", "std::vector<float>", &BTrkdxy2_vec);

  outTree_->Branch("Dmass",      "std::vector<float>", &Dmass_vec);
  outTree_->Branch("Dmass_flip", "std::vector<float>", &Dmass_flip_vec);
  outTree_->Branch("Dmass_ll",   "std::vector<float>", &Dmass_ll_vec);
  outTree_->Branch("Dmass_ll_flip", "std::vector<float>", &Dmass_ll_flip_vec);
  outTree_->Branch("Dmass_mukhad",  "std::vector<float>", &Dmass_mukhad_vec);   
  outTree_->Branch("Dmass_mukmu",   "std::vector<float>", &Dmass_mukmu_vec);
  outTree_->Branch("KLmassD0", "std::vector<float>", &KLmassD0_vec);

  outTree_->Branch("B_iso04_rel", "std::vector<float>", &b_iso04_rel_vec);
  outTree_->Branch("eleDR", "std::vector<float>", &eleDR_vec);
  outTree_->Branch("k_DCASig", "std::vector<float>", &k_DCASig_vec);
  outTree_->Branch("k_dzTrg", "std::vector<float>", &k_dzTrg_vec);
  outTree_->Branch("k_dxy_sig", "std::vector<float>", &k_dxy_sig_vec);
  outTree_->Branch("k_iso04_rel", "std::vector<float>", &k_iso04_rel_vec);
  outTree_->Branch("k_svip2d", "std::vector<float>", &k_svip2d_vec);
  outTree_->Branch("k_svip3d", "std::vector<float>", &k_svip3d_vec);
  outTree_->Branch("tag_dxy_sig", "std::vector<float>", &tag_dxy_sig_vec);
  outTree_->Branch("tag_dzTrg", "std::vector<float>", &tag_dzTrg_vec);
  outTree_->Branch("tag_iso04_rel", "std::vector<float>", &tag_iso04_rel_vec);
  outTree_->Branch("probe_dxy_sig", "std::vector<float>", &probe_dxy_sig_vec);
  outTree_->Branch("probe_dzTrg", "std::vector<float>", &probe_dzTrg_vec);
  outTree_->Branch("probe_iso04_rel", "std::vector<float>", &probe_iso04_rel_vec);

  outTree_->Branch("llkDR", "std::vector<float>", &llkDR_vec);
  outTree_->Branch("ptAsym", "std::vector<float>", &ptAsym_vec);

  outTree_->Branch("selectedPairsSize",  &selectedPairsSize,  "selectedPairsSize/I");   
}

void SkimmerWithKStar::bookOutputHistos() 
{
  cout << "Booking histos" << endl;
  //
  h_entries   = new TH1F("h_entries",  "Number of entries",   3,  3.5, 6.5);
  h_selection = new TH1F("h_selection","Selection breakdown", 6, -0.5, 5.5);
}

void SkimmerWithKStar::SetAnBdtWeights() {

  TFile *f_pflpt_commonRatio = new TFile("models/anBdtSystematics_PFLPT_common.root","READ");
  TFile *f_pfpf_commonRatio  = new TFile("models/anBdtSystematics_PFPF_common.root","READ"); 

  TH1F *anweight_pflpt_common = (TH1F*)f_pflpt_commonRatio->Get("ratio_commonBdt");
  TH1F *anweight_pfpf_common  = (TH1F*)f_pfpf_commonRatio->Get("ratio_commonBdt");
  
  for (int i = 0; i<anweight_pflpt_common->GetNbinsX(); i++) {

    float weight1=anweight_pflpt_common->GetBinContent(i+1);
    anweight_pflpt_common_.push_back(weight1);             

    float weight2=anweight_pfpf_common->GetBinContent(i+1);
    anweight_pfpf_common_.push_back(weight2);              

    float lowedge=anweight_pflpt_common->GetBinLowEdge(i+1); 
    lowedge_.push_back(lowedge);
  }
}

float SkimmerWithKStar::GetAnBdtWeight(float bdt, int pfpf ) {

  int thesizem1 = lowedge_.size()-1;

  float weight=1;
  for (int i = 0; i<thesizem1; i++) {   
    if (lowedge_[i]<=bdt && lowedge_[i+1]>bdt) { 
      if (pfpf==1) weight = anweight_pfpf_common_[i];  
      if (pfpf==0) weight = anweight_pflpt_common_[i];  
    }
  }
  
  return weight;
}

int main(int argc, char **argv){
  if (argc <2 ) {
    std::cerr << "Usage: " << argv[0] << " [file_list.txt] [outputFile] [0=data, 1=MC noK*, 2=MC with K*] [nmin] [nmax] [ptB] [svprob] [co2D]" << std::endl;
    return 1;
  }
  char inputFileName[500];
  char outputFileName[500];

  strcpy(inputFileName,argv[1]);
  if (argc <= 2 ) strcpy(outputFileName,argv[1]);
  else strcpy(outputFileName,argv[2]);
  int isMC=0;
  isMC=atoi(argv[3]);
  int nmin=0;
  int nmax=10000;

  if(argc>4){
    nmin=atoi(argv[4]);
    nmax=atoi(argv[5]);
  }
  float ptB=0.;
  float svprob =0.;
  float cos2D =0.;

  if(argc>6){
    ptB=atof(argv[6]);
    svprob=atof(argv[7]);
    cos2D=atof(argv[8]);
  }


  std::cout<<" your request: isMC="<<isMC<<" nmin="<<nmin<<" nmax="<<nmax<<" cuts ptb="<<ptB<<" svprob="<<svprob<<" cos2D="<<cos2D<<std::endl;
  
  // -------------------------
  // Loading the file
  TChain *theChain = new TChain("Events");
  char Buffer[5000];
  char MyRootFile[10000];
  std::cout << "input: " << inputFileName << std::endl;
  ifstream *inputFile = new ifstream(inputFileName);
  char tmpFileName[256];
  vector<string> filesToRemove;
  int nfiles=0;
  int i=0;
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,500);
    if(i>=nmin&& i<nmax){
      nfiles=nfiles+1;
      if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
	{
	  sscanf(Buffer,"%s",MyRootFile);
	  theChain->Add(TString(MyRootFile));
	  std::cout << "chaining " << MyRootFile << std::endl;
	}
    }
    i=i+1;
  }
  inputFile->close();
  delete inputFile;
  std::cout<<"we will process "<<nfiles<<" files"<< std::endl;
  
  SkimmerWithKStar tnp(theChain, isMC);
  std::cout<<"we are here 1"<< std::endl;

  tnp.PrepareOutputs(outputFileName);   
  std::cout<<"we are here 2 "<< std::endl;

  tnp.SetCuts(ptB,svprob,cos2D);
  std::cout<<"we are here 3 "<< std::endl;

  tnp.Loop();
  std::cout<<"we are here 4"<< std::endl;

  return 0;

}

