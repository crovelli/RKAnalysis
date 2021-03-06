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

  // To be set by hand   
  donvtxreweight_ = 0;
  nvtxWFileName_ = "/afs/cern.ch/user/c/crovelli/public/bphys/march/nvtxWeights_run2018D.root"; 
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

  // Load Analysis MVA weights for G.
  std::string bdtfile = "models/model_PFPF_george.txt";
  std::vector<std::string> feat = {"f0","f1", "f2","f3","f4","f5","f6","f7","f8","f9","f10","f11"};
  const auto fastForest = fastforest::load_txt(bdtfile.c_str(), feat);
  // Load Analysis MVA weights PFLP
  std::string bdtfilePFLP = "models/model_PFLP_george.txt";
  std::vector<std::string> featPFLP = {"f0","f1", "f2","f3","f4","f5","f6","f7","f8","f9","f10","f11"};
  const auto fastForestPFLP = fastforest::load_txt(bdtfilePFLP.c_str(), featPFLP);

  // Load Analysis MVA weights for Otto
  std::string bdtfileOttoPFPF = "models/otto_model.txt";
  std::vector<std::string> featOttoPFPF = {"f0","f1", "f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13","f14","f15","f16","f17","f18","f19","f20","f21", "f22","f23","f24","f25"};
  const auto fastForestOttoPFPF = fastforest::load_txt(bdtfileOttoPFPF.c_str(), featOttoPFPF);
  // PF+LP
  std::string bdtfileOttoPFLP = "models/otto_model_pflp.txt";
  std::vector<std::string> featOttoPFLP = {"f0","f1", "f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13","f14","f15","f16","f17","f18","f19","f20","f21", "f22","f23","f24","f25","f26","f27"};
  const auto fastForestOttoPFLP = fastforest::load_txt(bdtfileOttoPFLP.c_str(), featOttoPFLP);




  // Loop over events
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  cout << "entries : " <<  nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%10000==0) cout << jentry << endl;

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
    perEveW = 1.;
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
    int iHLT_Mu9_IP6 = (int)HLT_Mu9_IP6;
    hlt9 = iHLT_Mu9_IP6;
    h_selection->Fill(2.,perEveW);
    
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

    // calculation of BBDPhi
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
    h_selection->Fill(3.,perEveW);


    u_int iBGd[1000]={0};
    float ptBG[1000]={0.};
    int n=nBToKEE;
    if(n>1000) n=1000; // max 1000 triplets

    // order by pt
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

      // preparing variables
      int ele1_idx = BToKEE_l1Idx[iB];
      int ele2_idx = BToKEE_l2Idx[iB];
      
      float ele1_pt = BToKEE_fit_l1_pt[iB];
      float ele2_pt = BToKEE_fit_l2_pt[iB];
      float k_pt    = BToKEE_fit_k_pt[iB];
      
      float ele1_eta = BToKEE_fit_l1_eta[iB];     
      float ele2_eta = BToKEE_fit_l2_eta[iB];    
      float k_eta    = BToKEE_fit_k_eta[iB];    
      
      // B selection (standard cut)
      bool vtxFitSel = BToKEE_fit_pt[iB]>_ptB && BToKEE_svprob[iB]>_svProb && BToKEE_fit_cos2D[iB]>_cos2D;
      bool ele1Sel = fabs(ele1_eta)<2.4 && ele1_pt>0.7;
      bool ele2Sel = fabs(ele2_eta)<2.4 && ele2_pt>0.7;
      bool kSel = k_pt>0.7 && fabs(k_eta)<2.4; 
      bool additionalSel = BToKEE_fit_mass[iB]>4.5 && BToKEE_fit_mass[iB]<6.0;
      bool isBsel = vtxFitSel && ele1Sel && ele2Sel && kSel && additionalSel;
	
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

    if (goodBs.size()>0) h_selection->Fill(4.,perEveW);
    
    selectedBSize = goodBs.size();
    
    // At least one good B candidate
    if (goodBs.size()<=0) continue;
    h_selection->Fill(5.,perEveW);

    
    // ------------------------------------------------------------------
    for (u_int iB=0; iB<goodBs.size(); iB++) {
      
      // B-infos for further cuts offline
      int thisB = goodBs[iB];
      float thisBmass   = BToKEE_fit_mass[thisB];
      float thisBpt     = BToKEE_fit_pt[thisB];
      float thisBcos    = BToKEE_fit_cos2D[thisB];
      float thisBsvprob = BToKEE_svprob[thisB];
      float thisBxysig  = BToKEE_l_xy[thisB]/BToKEE_l_xy_unc[thisB];
      bool isThisAMcB = -1;
      if (sampleID>0) isThisAMcB = isMcB(thisB);      
      
      // K-infos for further cuts offline
      float thisKpt = BToKEE_fit_k_pt[thisB];   
      
      // Electrons
      int ele1_idx = BToKEE_l1Idx[thisB];
      int ele2_idx = BToKEE_l2Idx[thisB];
      int k_idx    = BToKEE_kIdx[thisB];  
      
      float k_pt     = BToKEE_fit_k_pt[thisB];   
      float ele1_pt  = BToKEE_fit_l1_pt[thisB];
      float ele2_pt  = BToKEE_fit_l2_pt[thisB];
      float k_eta    = BToKEE_fit_k_eta[thisB];   
      float ele1_eta = BToKEE_fit_l1_eta[thisB];
      float ele2_eta = BToKEE_fit_l2_eta[thisB];
      float k_phi    = BToKEE_fit_k_phi[thisB];   
      float ele1_phi = BToKEE_fit_l1_phi[thisB];
      float ele2_phi = BToKEE_fit_l2_phi[thisB];
      
      TLorentzVector kTLV(0,0,0,0);
      kTLV.SetPtEtaPhiM(k_pt,k_eta,k_phi,0.493677);

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

      // BToKEE_svprob
      // check if it exists another triplet with the same electrons and check the 
      // highest 2ndary vertex probability for those

      int iPbest=0;
      float prob4trk=0; 
      for (int iP=0; iP<nBToKEE; iP++){
	if(iP!=thisB){
	  if(( BToKEE_l1Idx[iP]==ele1_idx &&BToKEE_l2Idx[iP]==ele2_idx)||(BToKEE_l1Idx[iP]==ele2_idx &&BToKEE_l2Idx[iP]==ele1_idx)|| 
	     (ele1_pt== Electron_pt[BToKEE_l1Idx[iP]] && Electron_isPF[ele1_idx]==Electron_isPF[BToKEE_l1Idx[iP]] 
	      && ele2_pt== Electron_pt[BToKEE_l2Idx[iP]] && Electron_isPF[ele2_idx]==Electron_isPF[BToKEE_l2Idx[iP]])|| 
	     (ele2_pt== Electron_pt[BToKEE_l1Idx[iP]] && Electron_isPF[ele2_idx]==Electron_isPF[BToKEE_l1Idx[iP]] &&  
	      ele1_pt== Electron_pt[BToKEE_l2Idx[iP]] && Electron_isPF[ele1_idx]==Electron_isPF[BToKEE_l2Idx[iP]] ) ){
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
	  if(dr1<dr1m) {
	    dr1m=dr1; 	  
	  }
	  float dr2=DeltaR(ele2_eta, ele2_phi, Muon_eta[iM], Muon_phi[iM]);
	  if(dr2<dr2m) {
	    dr2m=dr2; 	  
	  }
	}
      }

      
      TLorentzVector ele1TLV(0,0,0,0);
      ele1TLV.SetPtEtaPhiM(ele1_pt,ele1_eta,ele1_phi,0.000511);
      TLorentzVector ele2TLV(0,0,0,0);
      ele2TLV.SetPtEtaPhiM(ele2_pt,ele2_eta,ele2_phi,0.000511);



      float x=goodBorder[iB];
      bPtOrder.push_back(x);
      float y=goodBorder2D[iB];
      b2DOrder.push_back(y);
      float z=goodBorderXY[iB];
      bXYOrder.push_back(z);
      p4Trk.push_back(prob4trk);
      p4TrkKStar.push_back(prob4trkKStar);
      
      tag_pt.push_back(Electron_pt[ele1_idx]);
      tag_eta.push_back(Electron_eta[ele1_idx]);
      tag_phi.push_back(Electron_phi[ele1_idx]);
      tag_isPF.push_back(Electron_isPF[ele1_idx]);           
      tag_isPFOverlap.push_back(Electron_isPFoverlap[ele1_idx]); 
      tag_isLowPt.push_back(Electron_isLowPt[ele1_idx]);
      tag_mvaId.push_back(Electron_mvaId[ele1_idx]);
      tag_pfmvaId.push_back(Electron_pfmvaId[ele1_idx]);
      tag_convveto.push_back(Electron_convVeto[ele1_idx]);
      //
      tag_drm.push_back(dr1m);
      probe_drm.push_back(dr2m);
      //
      fit_Bmass.push_back(thisBmass);
      fit_Bpt.push_back(thisBpt);
      fit_Bcos2D.push_back(thisBcos);
      fit_Bsvprob.push_back(thisBsvprob);
      fit_Bxysig.push_back(thisBxysig);
      if(sampleID>0){
	bmatchMC.push_back(isThisAMcB);
      }
      Kpt.push_back(thisKpt);
      K_eta.push_back(k_eta);
      K_phi.push_back(k_phi);
      probe_pt.push_back(Electron_pt[ele2_idx]);
      probe_eta.push_back(Electron_eta[ele2_idx]);
      probe_phi.push_back(Electron_phi[ele2_idx]);
      probe_isPF.push_back(Electron_isPF[ele2_idx]);
      probe_isPFOverlap.push_back(Electron_isPFoverlap[ele2_idx]); 
      probe_isLowPt.push_back(Electron_isLowPt[ele2_idx]); 
      probe_mvaId.push_back(Electron_mvaId[ele2_idx]);
      probe_pfmvaId.push_back(Electron_pfmvaId[ele2_idx]);
      probe_dxySig.push_back(Electron_dxy[ele2_idx]/Electron_dxyErr[ele2_idx]);  
      probe_dzSig.push_back(Electron_dz[ele2_idx]/Electron_dzErr[ele2_idx]);  
      probe_fBrem.push_back(Electron_fBrem[ele2_idx]);  
      probe_convveto.push_back(Electron_convVeto[ele2_idx]);
      mll_fullfit.push_back(BToKEE_mll_fullfit[thisB]);
      mll_raw.push_back(BToKEE_mll_raw[thisB]);
      fit_mass.push_back(BToKEE_fit_mass[thisB]); 
      
      // Vertex distances
      float l1kDz = fabs(kvz-l1vz);
      float l2kDz = fabs(kvz-l2vz);
      float LKdz = l1kDz;
      if (l2kDz<LKdz) LKdz = l2kDz;
      
      // Vectors
      TVector3 l1v3;
      TVector3 l2v3;
      TVector3 kv3;
      l1v3.SetPtEtaPhi(ele1_pt,ele1_eta,ele1_phi);
      l2v3.SetPtEtaPhi(ele2_pt,ele2_eta,ele2_phi);
      kv3.SetPtEtaPhi(k_pt,k_eta,k_phi);
	
      // DeltaR
      float l1kDr = l1v3.DeltaR(kv3);
      float l2kDr = l2v3.DeltaR(kv3);
      float LKdr = l1kDr;
      if (l2kDr<LKdr) LKdr = l2kDr;
      float L1L2dr = l1v3.DeltaR(l2v3);
      
      // Isolation
      float L1iso = BToKEE_l1_iso04[thisB];
      float L2iso = BToKEE_l2_iso04[thisB];
      float Kiso = BToKEE_k_iso04[thisB];
      
      
      float distance_b_trk1=10;         
      float distance_b_trk2=10;         
      float mean_b_trk=0;      
      int imean_b_trk=0 ;
      float vx_b=BToKEE_vtx_x[thisB];
      float vy_b=BToKEE_vtx_y[thisB];

      float BTrkdxy2=0;
      
      for (int itrk=0; itrk<nProbeTracks; itrk++){
	if(itrk==k_idx) continue;
	float vx_trk=ProbeTracks_vx[itrk];
	float vy_trk=ProbeTracks_vy[itrk];
	float vz_trk=ProbeTracks_vz[itrk];
	float trk_eta=ProbeTracks_eta[itrk]; 
	float trk_phi=ProbeTracks_phi[itrk];
	
	if(( DeltaR(ele1_eta,ele1_phi,trk_eta,trk_phi)<0.03)|| 
	   (DeltaR(ele2_eta,ele2_phi,trk_eta,trk_phi)<0.3)|| ( DeltaR(k_eta,k_phi,trk_eta,trk_phi)<0.03)) continue;
	
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
      
      // Analysis BDT George
      
      float L2id=Electron_mvaId[ele2_idx];
      if( Electron_isPF[ele2_idx]) L2id=Electron_pfmvaId[ele2_idx];

      LKdz_vec.push_back(LKdz);
      L1L2dr_vec.push_back(L1L2dr);
      LKdr_vec.push_back(LKdr);
      Kiso_vec.push_back(Kiso);
      BBDPhi_vec.push_back(BBDPhi);
      BTrkdxy2_vec.push_back(BTrkdxy2);


      float scoreBdt=0.;
      if(Electron_isPF[ele1_idx]&&Electron_isPF[ele2_idx]){
	std::vector<float> vecBdt = {thisBsvprob, thisBxysig, ele2_pt, thisKpt, thisBcos, LKdz,  L1L2dr, LKdr, L2id, Kiso, BBDPhi, BTrkdxy2 };
	scoreBdt = fastForest(vecBdt.data());
	analysisBdtG.push_back(scoreBdt);
      } else {
	// change here with the other BDT by George 
	std::vector<float> vecBdt = {thisBsvprob, thisBxysig, ele2_pt, thisKpt, thisBcos, LKdz,  L1L2dr, LKdr, L2id, Kiso, BBDPhi, BTrkdxy2 };
	scoreBdt = fastForestPFLP(vecBdt.data());
	analysisBdtG.push_back(scoreBdt);
      }


      // Analysis BDT Otto PFPF

      float BToKEE_fit_l1_normpt=BToKEE_fit_l1_pt[thisB]/BToKEE_fit_mass[thisB];
      float BToKEE_fit_l2_normpt=BToKEE_fit_l2_pt[thisB]/BToKEE_fit_mass[thisB];
      float BToKEE_l1_dxy_sig=(Electron_dxy[ele1_idx]) /Electron_dxyErr[ele1_idx];
      float BToKEE_l2_dxy_sig=(Electron_dxy[ele2_idx]) /Electron_dxyErr[ele2_idx];
      float BToKEE_fit_k_normpt=BToKEE_fit_k_pt[thisB] /BToKEE_fit_mass[thisB];
      float BToKEE_k_DCASig=ProbeTracks_DCASig[k_idx];
      float BToKEE_k_dxy_sig=ProbeTracks_dxyS[k_idx];
      float BToKEE_fit_normpt=BToKEE_fit_pt[thisB] /BToKEE_fit_mass[thisB];
      float BToKEE_l_xy_sig = (BToKEE_l_xy[thisB]) /BToKEE_l_xy_unc[thisB];
      float BToKEE_eleDR= DeltaR(ele1_eta,ele1_phi,ele2_eta,ele2_phi);
      TLorentzVector dll=ele1TLV+ele2TLV;
      
      
      float BToKEE_llkDR=dll.DeltaR(kTLV);
      float BToKEE_l1_iso04_rel=BToKEE_l1_iso04[thisB]/BToKEE_fit_l1_pt[thisB];
      float BToKEE_l2_iso04_rel=BToKEE_l2_iso04[thisB]/BToKEE_fit_l2_pt[thisB];
      float BToKEE_k_iso04_rel = BToKEE_k_iso04[thisB] / BToKEE_fit_k_pt[thisB];
      float BToKEE_b_iso04_rel =BToKEE_b_iso04[thisB]/BToKEE_fit_pt[thisB];
      TVector3 diele_p3 = dll.Vect();
      TVector3 k_p3 = kTLV.Vect();
      TVector3 pv2sv_p3(PV_x-BToKEE_vtx_x[thisB], PV_y-BToKEE_vtx_y[thisB], PV_z-BToKEE_vtx_z[thisB]);
      
      float BToKEE_ptAsym = ( (diele_p3.Cross(pv2sv_p3)).Mag() - (k_p3.Cross(pv2sv_p3)).Mag() ) 
			  / ( (diele_p3.Cross(pv2sv_p3)).Mag() + (k_p3.Cross(pv2sv_p3)).Mag() );
      
      float BToKEE_l1_dzTrg=Electron_dzTrg[ele1_idx];
      float BToKEE_l2_dzTrg=Electron_dzTrg[ele2_idx];
      float BToKEE_k_dzTrg=ProbeTracks_dzTrg[k_idx];
      
      float BToKEE_l1_pfmvaId_lowPt=Electron_pfmvaId[ele1_idx];
      if(ele1_pt>5) BToKEE_l1_pfmvaId_lowPt=20;

      float BToKEE_l2_pfmvaId_lowPt=Electron_pfmvaId[ele2_idx];
      if(ele2_pt>5) BToKEE_l2_pfmvaId_lowPt=20;
      
      float BToKEE_l1_pfmvaId_highPt=Electron_pfmvaId[ele1_idx];
      if(ele1_pt<=5) BToKEE_l1_pfmvaId_highPt=20;
      
      float BToKEE_l2_pfmvaId_highPt=Electron_pfmvaId[ele2_idx];
      if(ele2_pt<=5) BToKEE_l2_pfmvaId_highPt=20;

      float scoreBdtOtto=0.;
      if(Electron_isPF[ele1_idx]&&Electron_isPF[ele2_idx]){
	std::vector<float> vecBdtOtto = {BToKEE_b_iso04_rel, BToKEE_eleDR, thisBcos, BToKEE_fit_k_normpt,
					 BToKEE_fit_l1_normpt,BToKEE_fit_l2_normpt, BToKEE_fit_normpt,
					 BToKEE_k_DCASig, BToKEE_k_dzTrg, BToKEE_k_iso04_rel, BToKEE_k_svip2d[thisB], BToKEE_k_svip3d[thisB],
					  BToKEE_l1_dxy_sig, BToKEE_l1_dzTrg, BToKEE_l1_iso04_rel, 
					 //l1_mvaId,
					 BToKEE_l1_pfmvaId_highPt, BToKEE_l1_pfmvaId_lowPt,
					 BToKEE_l2_dxy_sig,BToKEE_l2_dzTrg,BToKEE_l2_iso04_rel,
					 // l2_mvaId,
					 BToKEE_l2_pfmvaId_highPt, BToKEE_l2_pfmvaId_lowPt, 
					 thisBxysig,BToKEE_llkDR, BToKEE_ptAsym, thisBsvprob};
	
	scoreBdtOtto = fastForestOttoPFPF(vecBdtOtto.data());
	analysisBdtO.push_back(scoreBdtOtto);
	/*	if(jentry<10) std::cout<<"1 "<<"0:"<< BToKEE_b_iso04_rel<<" 1:"<< BToKEE_eleDR<<" 2:"<< thisBcos<<" 3:"<< BToKEE_fit_k_normpt<<" 4:"<<
	BToKEE_fit_l1_normpt<<" 5:"<<BToKEE_fit_l2_normpt<<" 6:"<< BToKEE_fit_normpt<<" 7:"<<
	BToKEE_k_DCASig<<" 8:"<< BToKEE_k_dzTrg<<" 9:"<< BToKEE_k_iso04_rel<<" 10:"<< BToKEE_k_svip2d[thisB]<<" 11:"<< BToKEE_k_svip3d[thisB]<<" 12:"<<
	thisBxysig<<" 13:"<< BToKEE_l1_dxy_sig<<" 14:"<< BToKEE_l1_dzTrg<<" 15:"<< BToKEE_l1_iso04_rel<<" 16:"<<BToKEE_l1_pfmvaId_highPt<<" 17:"<< BToKEE_l1_pfmvaId_lowPt<<" 18:"<<
	BToKEE_l2_dxy_sig<<" 19:"<<BToKEE_l2_dzTrg<<" 20:"<<BToKEE_l2_iso04_rel<<" 21:"<<BToKEE_l2_pfmvaId_highPt<<" 22:"<< BToKEE_l2_pfmvaId_lowPt<<" 23:"<<
	BToKEE_llkDR<<" 24:"<< BToKEE_ptAsym<<" 25:"<< thisBsvprob<<std::endl;
	*/ 

      } else {
	std::vector<float> vecBdtOtto = {BToKEE_b_iso04_rel, BToKEE_eleDR, thisBcos, BToKEE_fit_k_normpt,
					 BToKEE_fit_l1_normpt,BToKEE_fit_l2_normpt, BToKEE_fit_normpt,
					 BToKEE_k_DCASig, BToKEE_k_dzTrg, BToKEE_k_iso04_rel, 
					 BToKEE_k_svip2d[thisB], BToKEE_k_svip3d[thisB],
					 BToKEE_l1_dxy_sig, BToKEE_l1_dzTrg, BToKEE_l1_iso04_rel, 
					 Electron_mvaId[ele1_idx],
					 BToKEE_l1_pfmvaId_highPt, BToKEE_l1_pfmvaId_lowPt,
					 BToKEE_l2_dxy_sig,BToKEE_l2_dzTrg,BToKEE_l2_iso04_rel,
					 Electron_mvaId[ele2_idx],
					 BToKEE_l2_pfmvaId_highPt, BToKEE_l2_pfmvaId_lowPt, 
					 thisBxysig,BToKEE_llkDR, BToKEE_ptAsym, thisBsvprob};
	
	scoreBdtOtto = fastForestOttoPFLP(vecBdtOtto.data());
	analysisBdtO.push_back(scoreBdtOtto);

	/* if(jentry<10)  std::cout<<"1 "<<"0:"<< BToKEE_b_iso04_rel<<" 1:"<< BToKEE_eleDR<<" 2:"<< thisBcos<<" 3:"<< BToKEE_fit_k_normpt<<" 4:"<<
	BToKEE_fit_l1_normpt<<" 5:"<<BToKEE_fit_l2_normpt<<" 6:"<< BToKEE_fit_normpt<<" 7:"<<
	BToKEE_k_DCASig<<" 8:"<< BToKEE_k_dzTrg<<" 9:"<< BToKEE_k_iso04_rel<<" 10:"<< BToKEE_k_svip2d[thisB]<<" 11:"<< BToKEE_k_svip3d[thisB]<<" 12:"<<
	thisBxysig<<" 13:"<< BToKEE_l1_dxy_sig<<" 14:"<< BToKEE_l1_dzTrg<<" 15:"<< BToKEE_l1_iso04_rel<<" 16:"<<
	Electron_mvaId[ele1_idx]<<" 17:"<<BToKEE_l1_pfmvaId_highPt<<" 18:"<< BToKEE_l1_pfmvaId_lowPt<<" 19:"<<
	BToKEE_l2_dxy_sig<<" 20:"<<BToKEE_l2_dzTrg<<" 21:"<<BToKEE_l2_iso04_rel<<" 22:"<<
	Electron_mvaId[ele2_idx]<<" 23:"<<BToKEE_l2_pfmvaId_highPt<<" 24:"<< BToKEE_l2_pfmvaId_lowPt<<" 25:"<<
	BToKEE_llkDR<<" 26:"<< BToKEE_ptAsym<<" 27:"<< thisBsvprob<<std::endl;
	*/
      }

      



      b_iso04_rel_vec.push_back(BToKEE_b_iso04_rel);
      eleDR_vec.push_back(BToKEE_eleDR);
      k_DCASig_vec.push_back(BToKEE_k_DCASig);
      k_dzTrg_vec.push_back(BToKEE_k_dzTrg);
      k_dxy_sig_vec.push_back(BToKEE_k_dxy_sig);
      k_iso04_rel_vec.push_back(BToKEE_k_iso04_rel);
      k_svip2d_vec.push_back(BToKEE_k_svip2d[thisB]);
      k_svip3d_vec.push_back(BToKEE_k_svip3d[thisB]);
      tag_dxy_sig_vec.push_back(BToKEE_l1_dxy_sig);
      tag_dzTrg_vec.push_back(BToKEE_l1_dzTrg);
      tag_iso04_rel_vec.push_back(BToKEE_l1_iso04_rel);
      probe_dxy_sig_vec.push_back(BToKEE_l2_dxy_sig);
      probe_dzTrg_vec.push_back(BToKEE_l2_dzTrg);
      probe_iso04_rel_vec.push_back(BToKEE_l2_iso04_rel);
      llkDR_vec.push_back(BToKEE_llkDR);
      ptAsym_vec.push_back(BToKEE_ptAsym);

      TLorentzVector l1_pihypoTLV(0,0,0,0);
      l1_pihypoTLV.SetPtEtaPhiM(ele1_pt,ele1_eta,ele1_phi,0.135);
      TLorentzVector l2_pihypoTLV(0,0,0,0);
      l2_pihypoTLV.SetPtEtaPhiM(ele2_pt,ele2_eta,ele2_phi,0.135);
      TLorentzVector l1_khypoTLV(0,0,0,0);
      l1_khypoTLV.SetPtEtaPhiM(ele1_pt,ele1_eta,ele1_phi,0.493677);
      TLorentzVector l2_khypoTLV(0,0,0,0);
      l2_khypoTLV.SetPtEtaPhiM(ele2_pt,ele2_eta,ele2_phi,0.493677);
      TLorentzVector k_pihypoTLV(0,0,0,0);
      k_pihypoTLV.SetPtEtaPhiM(k_pt,k_eta,k_phi,0.135);
      // kTLV exists already
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

      
      Dmass_vec.push_back(BToKEE_Dmass);
      Dmass_flip_vec.push_back(BToKEE_Dmass_flip);
      Dmass_ll_vec.push_back(BToKEE_Dmass_ll);
      Dmass_ll_flip_vec.push_back(BToKEE_Dmass_ll_flip);

      
      
      if (sampleID>0) {     // MC      
	bool isTagAMcEleFromJPsi   = isMcEleFromJPsi(ele1_idx);
	bool isProbeAMcEleFromJPsi = isMcEleFromJPsi(ele2_idx);
	//
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
    h_selection->Fill(6.,perEveW);

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
    tag_mvaId.clear();  
    tag_pfmvaId.clear();  
    tag_convveto.clear();
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
    fit_mass.clear(); 

    fit_Bpt.clear();
    fit_Bcos2D.clear();
    fit_Bsvprob.clear();    
    fit_Bxysig.clear();    

    Kpt.clear();
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
    probe_mvaId.clear();  
    probe_pfmvaId.clear();  
    probe_dxySig.clear();  
    probe_dzSig.clear();  
    probe_convveto.clear();

    analysisBdtG.clear();
    analysisBdtO.clear();
    LKdz_vec.clear();
    L1L2dr_vec.clear();
    LKdr_vec.clear();
    Kiso_vec.clear();
    BBDPhi_vec.clear();
    BTrkdxy2_vec.clear(); 
    Dmass_vec.clear();
    Dmass_flip_vec.clear();
    Dmass_ll_vec.clear();
    Dmass_ll_flip_vec.clear();

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

  outTree_->Branch("hlt9", &hlt9, "hlt9/I");
  outTree_->Branch("trg_muon_pt", &trg_muon_pt, "trg_muon_pt/F");

  outTree_->Branch("selectedBSize",  &selectedBSize,  "selectedBSize/I");   
  
  outTree_->Branch("tag_pt", "std::vector<float>", &tag_pt);  
  outTree_->Branch("tag_eta", "std::vector<float>", &tag_eta);  
  outTree_->Branch("tag_phi", "std::vector<float>", &tag_phi);  
  outTree_->Branch("tag_isPF", "std::vector<bool>", &tag_isPF);  
  outTree_->Branch("tag_isPFOverlap", "std::vector<bool>", &tag_isPFOverlap);  
  outTree_->Branch("tag_isLowPt", "std::vector<bool>", &tag_isLowPt);  
  outTree_->Branch("tag_mvaId", "std::vector<float>", &tag_mvaId);  
  outTree_->Branch("tag_pfmvaId", "std::vector<float>", &tag_pfmvaId);  
  outTree_->Branch("tag_convveto", "std::vector<bool>", &tag_convveto);  

  if(sampleID>0){
    outTree_->Branch("tag_ptMc", "std::vector<float>", &tag_ptMc);
    outTree_->Branch("probe_ptMc", "std::vector<float>", &probe_ptMc);
    outTree_->Branch("tag_matchMcFromJPsi", "std::vector<bool>", &tag_matchMcFromJPsi);  
    outTree_->Branch("tag_matchMc", "std::vector<bool>", &tag_matchMc);  
  }
  outTree_->Branch("mll_fullfit", "std::vector<float>", &mll_fullfit);  
  outTree_->Branch("mll_raw", "std::vector<float>", &mll_raw);  
  outTree_->Branch("fit_mass", "std::vector<float>", &fit_mass);  

  outTree_->Branch("fit_Bpt",     "std::vector<float>", &fit_Bpt);  
  outTree_->Branch("fit_Bcos2D",  "std::vector<float>", &fit_Bcos2D);  
  outTree_->Branch("fit_Bsvprob", "std::vector<float>", &fit_Bsvprob);  
  outTree_->Branch("fit_Bxysig",  "std::vector<float>", &fit_Bxysig);  
  outTree_->Branch("bmatchMC", "std::vector<bool>", &bmatchMC);  
  outTree_->Branch("K_pt", "std::vector<float>", &Kpt);  
  outTree_->Branch("K_eta", "std::vector<float>", &K_eta);  
  outTree_->Branch("K_phi", "std::vector<float>", &K_phi);  

  outTree_->Branch("probe_pt", "std::vector<float>", &probe_pt);  
  outTree_->Branch("probe_eta", "std::vector<float>", &probe_eta);  
  outTree_->Branch("probe_phi", "std::vector<float>", &probe_phi);  
  outTree_->Branch("probe_isPF", "std::vector<bool>", &probe_isPF);  
  outTree_->Branch("probe_isPFOverlap", "std::vector<bool>", &probe_isPFOverlap);  
  outTree_->Branch("probe_isLowPt", "std::vector<bool>", &probe_isLowPt);  
  outTree_->Branch("probe_mvaId", "std::vector<float>", &probe_mvaId);  
  outTree_->Branch("probe_pfmvaId", "std::vector<float>", &probe_pfmvaId);  
  outTree_->Branch("probe_dxySig", "std::vector<float>", &probe_dxySig);  
  outTree_->Branch("probe_dzSig", "std::vector<float>", &probe_dzSig);  
  outTree_->Branch("probe_convveto", "std::vector<bool>", &probe_convveto);  

  if(sampleID>0){
    outTree_->Branch("probe_matchMcFromJPsi", "std::vector<bool>", &probe_matchMcFromJPsi);  
    outTree_->Branch("probe_matchMc", "std::vector<bool>", &probe_matchMc);  
  }
  outTree_->Branch("probe_drm", "std::vector<float>", &probe_drm);  
  outTree_->Branch("tag_drm", "std::vector<float>", &tag_drm);  
  outTree_->Branch("bPtOrder", "std::vector<float>", &bPtOrder);  
  outTree_->Branch("b2DOrder", "std::vector<float>", &b2DOrder);  
  outTree_->Branch("bXYOrder", "std::vector<float>", &bXYOrder);  
  outTree_->Branch("p4Trk", "std::vector<float>", &p4Trk);  
  if(sampleID>1){
    outTree_->Branch("p4TrkKStar", "std::vector<float>", &p4TrkKStar);  
  }
  outTree_->Branch("analysisBdtG", "std::vector<float>", &analysisBdtG);
  outTree_->Branch("analysisBdtO", "std::vector<float>", &analysisBdtO);
  
  outTree_->Branch("LKdz", "std::vector<float>", &LKdz_vec);
  outTree_->Branch("L1L2dr", "std::vector<float>", &L1L2dr_vec);
  outTree_->Branch("LKdr", "std::vector<float>", &LKdr_vec);
  outTree_->Branch("Kiso", "std::vector<float>", &Kiso_vec);
  outTree_->Branch("BBDPhi", "std::vector<float>", &BBDPhi_vec);
  outTree_->Branch("BTrkdxy2", "std::vector<float>", &BTrkdxy2_vec);

  outTree_->Branch("Dmass", "std::vector<float>", &Dmass_vec);
  outTree_->Branch("Dmass_flip", "std::vector<float>", &Dmass_flip_vec);
  outTree_->Branch("Dmass_ll", "std::vector<float>", &Dmass_ll_vec);
  outTree_->Branch("Dmass_ll_flip", "std::vector<float>", &Dmass_ll_flip_vec);

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
  h_selection = new TH1F("h_selection","Selection breakdown", 8, -0.5, 7.5);
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

  tnp.Loop();
  std::cout<<"we are here 3"<< std::endl;

  return 0;

}

