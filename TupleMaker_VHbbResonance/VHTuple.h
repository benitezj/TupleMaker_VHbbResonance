#ifndef VHTUPLE
#define VHTUPLE

#include "RecoTuple.h"

class VHTuple : public RecoTuple {
  
public:
  VHTuple();

  ///Define here the block to hold the VH resonance quantities
  const static int MAXVH=200;
  int nvh=0;
  int vh_charge[MAXVH];
  float vh_m[MAXVH];//final mass variable
  float vh_m1[MAXVH];
  float vh_m2[MAXVH];
  float vh_m3[MAXVH];
  float vh_m4[MAXVH];
  float vh_E[MAXVH];
  float vh_p[MAXVH];
  float vh_pt[MAXVH];
  float vh_phi[MAXVH];
  float vh_eta[MAXVH];
  float vh_dR[MAXVH];
  int   vh_v[MAXVH];//index of V
  int   vh_h[MAXVH];//index of H
  float vh_dPhiHMET[MAXVH];
  float vh_h_truthdR[MAXVH];//deltaR w.r.t generated Higgs 

  //Lepton selection flags
  int eve_passVHMultiLepVeto=0;

  ////These are the selected candidates
  int vh=0;//index of VH resonance candidate

  ///Jets with b-tag requirements
  int vh_nbjet=0;//number of b-tagged jets
  int vh_bjet1=0;//leading b-tagged jet
  int vh_bjet2=0;//subleading b-tagged jet

  ///Jets with standard jet selections
  int vh_njet=0;//signal candidate jets
  int vh_jet1=0;//leading signal jet
  int vh_jet2=0;
  int vh_jet3=0;

  int vh_njetveto=0;//for central jet veto

  int vh_njetfw=0;//for forward jet veto
  int vh_jetfw=0;//leading forward jet


  ///This is needed for Fat Jet analysis
  int vh_nfatjet=0;
  int vh_fatjet1=0;

  //additional bjets not overlapping with the selected leptons and Higgs
  //Note that this is specific to the V-H combination, should be made part of the vh block or the vh block should be removed
  int vh_naddtrkbjet=0;


  virtual void DefineBranches(TTree * tr); 

private:
  
};
#endif 
