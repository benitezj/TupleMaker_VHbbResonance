#ifndef BASETUPLE
#define BASETUPLE

using namespace std;

#include "TTree.h"

class BaseTuple {
 public:
  BaseTuple();
  virtual ~BaseTuple();


  // Event level properties
  ULong64_t eve_num=0;
  int eve_run=0;
  int eve_lb=0;
  int eve_isMC=0;
  float eve_mu=0;
  float eve_mc_w=0;
  int eve_mc_num=0;
  int eve_mc_chan=0;

  //selection flags
  int eve_isCleanEvent=0;
  int eve_passTrig=0;
  int eve_passTrigMatch=0;
  
  //trigger weight
  float eve_trig_w=1.;
  float eve_trig_w_EL_EFF_Trigger_TotalCorrUncertainty__1up=1.;
  float eve_trig_w_EL_EFF_Trigger_TotalCorrUncertainty__1down=1.;
  float eve_trig_w_MUON_EFF_TrigSystUncertainty__1up=1.;
  float eve_trig_w_MUON_EFF_TrigSystUncertainty__1down=1.;
  float eve_trig_w_MUON_EFF_TrigStatUncertainty__1up=1.;
  float eve_trig_w_MUON_EFF_TrigStatUncertainty__1down=1.;

  ////pile-up weights
  const static int MAXMCPUWEIGHTS=50;
  int nmcpuw=0;
  float eve_mc_puw[MAXMCPUWEIGHTS];//one for each data/mc period

  const static int MAXTRUTH=300;
  int ntruth=0;
  float truth_m[MAXTRUTH];
  float truth_p[MAXTRUTH];
  float truth_pt[MAXTRUTH];
  float truth_eta[MAXTRUTH];
  float truth_phi[MAXTRUTH];
  //int truth_charge[MAXTRUTH];
  int truth_status[MAXTRUTH];
  int truth_pdg[MAXTRUTH];
  
  //Vector boson
  int truth_V=0;//Z or W needed for combining pT slices

  //Specific channels
  int truth_Zmm=0;
  int truth_Zmm_m1=0;
  int truth_Zmm_m2=0;

  int truth_Zee=0;
  int truth_Zee_e1=0;
  int truth_Zee_e2=0;

  int truth_Ztt=0;
  int truth_Ztt_t1=0;
  int truth_Ztt_t2=0;

  int truth_Wmv=0;
  int truth_Wmv_m=0;
  int truth_Wmv_v=0;

  int truth_Wev=0;
  int truth_Wev_e=0;
  int truth_Wev_v=0;

  int truth_Wtv=0;
  int truth_Wtv_t=0;
  int truth_Wtv_v=0;
  
  int truth_top=0;
  int truth_top_mm_m1=0;
  int truth_top_mm_m2=0;
  int truth_top_ee_e1=0;
  int truth_top_ee_e2=0;
  int truth_top_em_e=0;
  int truth_top_em_m=0;

  int truth_ttbar=0;  
  int truth_ttbar_t1=0;
  int truth_ttbar_t2=0;


  int truth_H0bb=0;
  int truth_H0bb_b1=0;
  int truth_H0bb_b2=0;

  //these composites are defined in a special way:
  //A new truth particle will be created from the V and H
  //and will be added at the end of the truth list
  int truth_ZHmmbb=0;
  int truth_ZHeebb=0;
  int truth_WHmvbb=0;
  int truth_WHevbb=0;


  ///VH resonance analysis
  int truth_RToVH=0;
  int truth_RToVH_V=0;
  int truth_RToVH_H=0;
  int truth_RToVH_H_b1=0;
  int truth_RToVH_H_b2=0;
  float truth_RToVH_H_dR=0;
  int truth_RToVH_H_decay=0;

  virtual void DefineBranches(TTree * tr); 
  
private:
  
};
#endif // _HBBTUPLENEW_HPP
