#ifndef RECOTUPLE
#define RECOTUPLE
#include "xAODJet/Jet.h"
#include "BaseTuple.h"

class RecoTuple : public BaseTuple {

public:
  RecoTuple();
  ~RecoTuple();


  //electron SF's
  float ele_reco_w=1.0;//effSFReco  * effSFReco
  float ele_id_w=1.0;// effSFmediumLH * effSFlooseLH *
  float ele_iso_w=1.0;//effSFIsoLooseTrackOnlyMediumLH * effSFIsoLooseTrackOnlyLooseLH

  float ele_reco_w_EL_EFF_Reco_TotalCorrUncertainty__1up=1.0;
  float ele_reco_w_EL_EFF_Reco_TotalCorrUncertainty__1down=1.0;
  float ele_id_w_EL_EFF_ID_TotalCorrUncertainty__1up=1.0;
  float ele_id_w_EL_EFF_ID_TotalCorrUncertainty__1down=1.0;
  float ele_iso_w_EL_EFF_Iso_TotalCorrUncertainty__1up=1.0;
  float ele_iso_w_EL_EFF_Iso_TotalCorrUncertainty__1down=1.0;

  
  //muon SF's
  float muo_ttv_w=1.0;//vertex association weight
  float muo_id_w=1.0;//looseEffSF * mediumEffSF
  float muo_iso_w=1.0;//looseTrackOnlyIsoSF * looseTrackOnlyIsoSF

  float muo_ttv_w_MUON_TTVA_STAT__1up=1.0;
  float muo_ttv_w_MUON_TTVA_STAT__1down=1.0;
  float muo_ttv_w_MUON_TTVA_SYS__1up=1.0;
  float muo_ttv_w_MUON_TTVA_SYS__1down=1.0;
  float muo_id_w_MUON_EFF_STAT__1up=1.0;
  float muo_id_w_MUON_EFF_STAT__1down=1.0;
  float muo_id_w_MUON_EFF_SYS__1up=1.0;
  float muo_id_w_MUON_EFF_SYS__1down=1.0;
  float muo_iso_w_MUON_ISO_STAT__1up=1.0;
  float muo_iso_w_MUON_ISO_STAT__1down=1.0;
  float muo_iso_w_MUON_ISO_SYS__1up=1.0;
  float muo_iso_w_MUON_ISO_SYS__1down=1.0;


  //b-tag weight
  float eve_btag_w = 1.0;
  float eve_btag_w_FT_EFF_Eigen_B_0__1down = 1.0;
  float eve_btag_w_FT_EFF_Eigen_B_0__1up = 1.0;
  float eve_btag_w_FT_EFF_Eigen_B_1__1down = 1.0;
  float eve_btag_w_FT_EFF_Eigen_B_1__1up = 1.0;
  float eve_btag_w_FT_EFF_Eigen_B_2__1down = 1.0;
  float eve_btag_w_FT_EFF_Eigen_B_2__1up = 1.0;
  float eve_btag_w_FT_EFF_Eigen_C_0__1down = 1.0;
  float eve_btag_w_FT_EFF_Eigen_C_0__1up = 1.0;
  float eve_btag_w_FT_EFF_Eigen_C_1__1down = 1.0;
  float eve_btag_w_FT_EFF_Eigen_C_1__1up = 1.0;
  float eve_btag_w_FT_EFF_Eigen_C_2__1down = 1.0;
  float eve_btag_w_FT_EFF_Eigen_C_2__1up = 1.0;
  float eve_btag_w_FT_EFF_Eigen_C_3__1down = 1.0;
  float eve_btag_w_FT_EFF_Eigen_C_3__1up = 1.0;
  float eve_btag_w_FT_EFF_Eigen_Light_0__1down = 1.0;
  float eve_btag_w_FT_EFF_Eigen_Light_0__1up = 1.0;
  float eve_btag_w_FT_EFF_Eigen_Light_1__1down = 1.0;
  float eve_btag_w_FT_EFF_Eigen_Light_1__1up = 1.0;
  float eve_btag_w_FT_EFF_Eigen_Light_2__1down = 1.0;
  float eve_btag_w_FT_EFF_Eigen_Light_2__1up = 1.0;
  float eve_btag_w_FT_EFF_Eigen_Light_3__1down = 1.0;
  float eve_btag_w_FT_EFF_Eigen_Light_3__1up = 1.0;
  float eve_btag_w_FT_EFF_Eigen_Light_4__1down = 1.0;
  float eve_btag_w_FT_EFF_Eigen_Light_4__1up = 1.0;
  float eve_btag_w_FT_EFF_extrapolation_from_charm__1down = 1.0;
  float eve_btag_w_FT_EFF_extrapolation_from_charm__1up = 1.0;
  float eve_btag_w_FT_EFF_extrapolation__1down = 1.0;
  float eve_btag_w_FT_EFF_extrapolation__1up = 1.0;

  //vertices
  const static int MAXVERTICES=100;
  int nvtx=0;
  float vtx_x[MAXVERTICES];
  float vtx_y[MAXVERTICES];
  float vtx_z[MAXVERTICES];

  //MET
  float met_sumet=0;
  float met_pt=0;
  float met_eta=0;
  float met_phi=0;
  float met_px=0;
  float met_py=0;

  //Electrons
  const static int MAXELECTRONS=30;
  int nele=0;
  float ele_E[MAXELECTRONS];
  float ele_m[MAXELECTRONS];
  float ele_pt[MAXELECTRONS];
  float ele_phi[MAXELECTRONS];
  float ele_eta[MAXELECTRONS];
  float ele_cluseta[MAXELECTRONS];
  float ele_z0[MAXELECTRONS];
  float ele_z0sintheta[MAXELECTRONS];
  float ele_d0[MAXELECTRONS];
  float ele_d0sig[MAXELECTRONS];
  int ele_sharedtrk[MAXELECTRONS];
  int ele_isFromPV[MAXELECTRONS];
  int ele_charge[MAXELECTRONS];
  float ele_ptiso[MAXELECTRONS];
  float ele_etiso[MAXELECTRONS];
  int ele_idVHl[MAXELECTRONS];
  int ele_idZHl[MAXELECTRONS];
  int ele_idWHl[MAXELECTRONS];
  int ele_passor[MAXELECTRONS];
  float ele_trkpt[MAXELECTRONS];
  int ele_isoWP1[MAXELECTRONS];
  int ele_isoWP2[MAXELECTRONS];
  int ele_isoWP3[MAXELECTRONS];
  int ele_isoWP4[MAXELECTRONS];
  int ele_isoWP5[MAXELECTRONS];
  int ele_isoWP6[MAXELECTRONS];

  int ele_isLooseIDElectron[MAXELECTRONS];
  int ele_isMediumIDElectron[MAXELECTRONS];
  int ele_isTightIDElectron[MAXELECTRONS];

  int ele_matchHLT_e24_lhmedium_L1EM18VH[MAXELECTRONS];
  int ele_matchHLT_e24_lhmedium_L1EM20VH[MAXELECTRONS];
  int ele_matchHLT_e26_lhtight_iloose[MAXELECTRONS];
  int ele_matchHLT_e60_lhmedium[MAXELECTRONS];
  int ele_matchHLT_e120_lhloose[MAXELECTRONS];


  //Muons
  const static int MAXMUONS=20;
  int nmuo=0;
  float muo_E[MAXMUONS];
  float muo_m[MAXMUONS];
  float muo_pt[MAXMUONS];
  float muo_phi[MAXMUONS];
  float muo_eta[MAXMUONS];
  float muo_z0[MAXMUONS];
  float muo_z0sintheta[MAXMUONS];
  float muo_d0[MAXMUONS];
  float muo_d0sig[MAXMUONS];
  int muo_isFromPV[MAXMUONS];
  int muo_charge[MAXMUONS];
  float muo_ptiso[MAXMUONS];
  float muo_etiso[MAXMUONS];
  int muo_quality[MAXMUONS];
  int muo_idVHl[MAXMUONS];
  int muo_idZH[MAXMUONS];
  int muo_idWH[MAXMUONS];
  int muo_passor[MAXMUONS];
  float muo_trkpt[MAXMUONS];
  int muo_isoWP1[MAXMUONS];
  int muo_isoWP2[MAXMUONS];
  int muo_isoWP3[MAXMUONS];
  int muo_isoWP4[MAXMUONS];
  int muo_isoWP5[MAXMUONS];
  int muo_isoWP6[MAXMUONS];

  int muo_isLooseIDMuon[MAXMUONS];
  int muo_isMediumIDMuon[MAXMUONS];
  int muo_isTightIDMuon[MAXMUONS];

  int muo_matchHLT_mu20_iloose_L1MU15[MAXMUONS];
  int muo_matchHLT_mu26_imedium[MAXMUONS];
  int muo_matchHLT_mu50[MAXMUONS];


  //Taus
  const static int MAXTAUS=30;
  int ntau=0;
  float tau_E[MAXTAUS];
  float tau_m[MAXTAUS];
  float tau_pt[MAXTAUS];
  float tau_phi[MAXTAUS];
  float tau_eta[MAXTAUS];
  float tau_z0[MAXTAUS];
  float tau_d0[MAXTAUS];
  int tau_charge[MAXTAUS];
  int tau_passor[MAXTAUS];

  //Jets
  const static int MAXJETS=50;
  int njet=0;
  int jet_charge[MAXJETS];
  float jet_E[MAXJETS];
  float jet_m[MAXJETS];
  float jet_pt[MAXJETS];
  float jet_phi[MAXJETS];
  float jet_eta[MAXJETS];
  int jet_ntrk[MAXJETS];
  int jet_NumTrkPt500PV[MAXJETS];
  float jet_SumPtTrkPt500PV[MAXJETS];
  float jet_massvx[MAXJETS];
  float jet_normdist[MAXJETS];
  float jet_jvf[MAXJETS];
  float jet_jvt[MAXJETS];

  int jet_good[MAXJETS];
  int jet_passor[MAXJETS];

  int jet_btag[MAXJETS];
  float jet_mv1[MAXJETS];
  float jet_sv1ip3d[MAXJETS];
  float jet_mv2c00[MAXJETS];
  float jet_mv2c10[MAXJETS];
  float jet_mv2c20[MAXJETS];
  float jet_mvb[MAXJETS];


  //Jets
  int nfatjet=0;
  float fatjet_m[MAXJETS];
  float fatjet_pt[MAXJETS];
  float fatjet_phi[MAXJETS];
  float fatjet_eta[MAXJETS];
  int   fatjet_ltrkjet[MAXJETS];
  int   fatjet_sltrkjet[MAXJETS];
  int   fatjet_ntrkjet[MAXJETS];
  int   fatjet_ntrkjetb[MAXJETS];
  int   fatjet_XbbL2b[MAXJETS];
  int   fatjet_XbbL2bmH[MAXJETS];
  int   fatjet_XbbnTag[MAXJETS];
  float fatjet_m_cor[MAXJETS];
  float fatjet_pt_cor[MAXJETS];
  float fatjet_phi_cor[MAXJETS];
  float fatjet_eta_cor[MAXJETS];
  const xAOD::Jet* fatjet_address[MAXJETS];//needed at run time to check track jet links


  //track Jets
  int ntrkjet=0;
  float trkjet_m[MAXJETS];
  float trkjet_pt[MAXJETS];
  float trkjet_phi[MAXJETS];
  float trkjet_eta[MAXJETS];
  int   trkjet_numConstituents[MAXJETS];
  int   trkjet_btag[MAXJETS];
  float trkjet_sv1ip3d[MAXJETS];
  float trkjet_mv2c00[MAXJETS];
  float trkjet_mv2c10[MAXJETS];
  float trkjet_mv2c20[MAXJETS];
  const xAOD::Jet* trkjet_address[MAXJETS];//needed at run time to check link to fat jet

  //electron pairs
  const static int MAXEE=100;
  int nee=0;
  int ee_charge[MAXEE];
  float ee_m[MAXEE];
  float ee_E[MAXEE];
  float ee_p[MAXEE];
  float ee_pt[MAXEE];
  float ee_phi[MAXEE];
  float ee_eta[MAXEE];
  float ee_dR[MAXEE];
  int   ee_leg1[MAXEE];
  int   ee_leg2[MAXEE];

  //muon pairs
  const static int MAXMM=100;
  int nmm=0;
  int mm_charge[MAXMM];
  float mm_m[MAXMM];
  float mm_E[MAXMM];
  float mm_p[MAXMM];
  float mm_pt[MAXMM];
  float mm_phi[MAXMM];
  float mm_eta[MAXMM];
  float mm_dR[MAXMM];
  int   mm_leg1[MAXMM];
  int   mm_leg2[MAXMM];

  //e-mu pairs
  const static int MAXEM=100;
  int nem=0;
  int em_charge[MAXEM];
  float em_m[MAXEM];
  float em_E[MAXEM];
  float em_p[MAXEM];
  float em_pt[MAXEM];
  float em_phi[MAXEM];
  float em_eta[MAXEM];
  float em_dR[MAXEM];
  int   em_leg1[MAXEM];
  int   em_leg2[MAXEM];

  //tau pairs
  const static int MAXTT=100;
  int ntt=0;
  int tt_charge[MAXTT];
  float tt_m[MAXTT];
  float tt_E[MAXTT];
  float tt_p[MAXTT];
  float tt_pt[MAXTT];
  float tt_phi[MAXTT];
  float tt_eta[MAXTT];
  float tt_dR[MAXTT];
  int   tt_leg1[MAXTT];
  int   tt_leg2[MAXTT];

  //jet pairs
  const static int MAXJJ=300;
  int njj=0;
  int jj_charge[MAXJJ];
  float jj_m[MAXJJ];
  float jj_E[MAXJJ];
  float jj_p[MAXJJ];
  float jj_pt[MAXJJ];
  float jj_phi[MAXJJ];
  float jj_eta[MAXJJ];
  float jj_dR[MAXJJ];
  int   jj_leg1[MAXJJ];
  int   jj_leg2[MAXJJ];


  //e-MET pairs
  const static int MAXEMET=MAXELECTRONS;
  int neMet=0;
  int eMet_charge[MAXEMET];
  float eMet_m[MAXEMET];
  float eMet_E[MAXEMET];
  float eMet_p[MAXEMET];
  float eMet_pt[MAXEMET];
  float eMet_phi[MAXEMET];
  float eMet_eta[MAXEMET];
  float eMet_mT[MAXEMET];
  int eMet_ele[MAXEMET];//link to the electron

  //mu-MET pairs
  const static int MAXMUMET=MAXMUONS;
  int nmuMet=0;
  int muMet_charge[MAXMUMET];
  float muMet_m[MAXMUMET];
  float muMet_E[MAXMUMET];
  float muMet_p[MAXMUMET];
  float muMet_pt[MAXMUMET];
  float muMet_phi[MAXMUMET];
  float muMet_eta[MAXMUMET];
  float muMet_mT[MAXMUMET];
  int muMet_mu[MAXMUMET];//link to the muon
  

  virtual void DefineBranches(TTree * tr); 
  
private:

  
};
#endif // _HBBTUPLENEW_HPP
