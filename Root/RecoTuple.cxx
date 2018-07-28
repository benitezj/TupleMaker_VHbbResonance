#include "TupleMaker_VHbbResonance/RecoTuple.h"

RecoTuple::RecoTuple():
  BaseTuple()
{
}

RecoTuple::~RecoTuple() {;}

void RecoTuple::DefineBranches(TTree * tr) {
  BaseTuple::DefineBranches(tr);

  //scale factors
  tr->Branch("ele_reco_w",&ele_reco_w,"ele_reco_w/F");
  tr->Branch("ele_id_w",&ele_id_w,"ele_id_w/F");
  tr->Branch("ele_iso_w",&ele_iso_w,"ele_iso_w/F");

  tr->Branch("ele_reco_w_EL_EFF_Reco_TotalCorrUncertainty__1up",&ele_reco_w_EL_EFF_Reco_TotalCorrUncertainty__1up,"ele_reco_w_EL_EFF_Reco_TotalCorrUncertainty__1up/F");
  tr->Branch("ele_reco_w_EL_EFF_Reco_TotalCorrUncertainty__1down",&ele_reco_w_EL_EFF_Reco_TotalCorrUncertainty__1down,"ele_reco_w_EL_EFF_Reco_TotalCorrUncertainty__1down/F");
  tr->Branch("ele_id_w_EL_EFF_ID_TotalCorrUncertainty__1up",&ele_id_w_EL_EFF_ID_TotalCorrUncertainty__1up,"ele_id_w_EL_EFF_ID_TotalCorrUncertainty__1up/F");
  tr->Branch("ele_id_w_EL_EFF_ID_TotalCorrUncertainty__1down",&ele_id_w_EL_EFF_ID_TotalCorrUncertainty__1down,"ele_id_w_EL_EFF_ID_TotalCorrUncertainty__1down/F");
  tr->Branch("ele_iso_w_EL_EFF_Iso_TotalCorrUncertainty__1up",&ele_iso_w_EL_EFF_Iso_TotalCorrUncertainty__1up,"ele_iso_w_EL_EFF_Iso_TotalCorrUncertainty__1up/F");
  tr->Branch("ele_iso_w_EL_EFF_Iso_TotalCorrUncertainty__1down",&ele_iso_w_EL_EFF_Iso_TotalCorrUncertainty__1down,"ele_iso_w_EL_EFF_Iso_TotalCorrUncertainty__1down/F");


  tr->Branch("muo_ttv_w",&muo_ttv_w,"muo_ttv_w/F");
  tr->Branch("muo_id_w",&muo_id_w,"muo_id_w/F");
  tr->Branch("muo_iso_w",&muo_iso_w,"muo_iso_w/F");

  tr->Branch("muo_ttv_w_MUON_TTVA_STAT__1up",&muo_ttv_w_MUON_TTVA_STAT__1up,"muo_ttv_w_MUON_TTVA_STAT__1up/F");
  tr->Branch("muo_ttv_w_MUON_TTVA_STAT__1down",&muo_ttv_w_MUON_TTVA_STAT__1down,"muo_ttv_w_MUON_TTVA_STAT__1down/F");
  tr->Branch("muo_ttv_w_MUON_TTVA_SYS__1up",&muo_ttv_w_MUON_TTVA_SYS__1up,"muo_ttv_w_MUON_TTVA_SYS__1up/F");
  tr->Branch("muo_ttv_w_MUON_TTVA_SYS__1down",&muo_ttv_w_MUON_TTVA_SYS__1down,"muo_ttv_w_MUON_TTVA_SYS__1down/F");
  tr->Branch("muo_id_w_MUON_EFF_STAT__1up",&muo_id_w_MUON_EFF_STAT__1up,"muo_id_w_MUON_EFF_STAT__1up/F");
  tr->Branch("muo_id_w_MUON_EFF_STAT__1down",&muo_id_w_MUON_EFF_STAT__1down,"muo_id_w_MUON_EFF_STAT__1down/F");
  tr->Branch("muo_id_w_MUON_EFF_SYS__1up",&muo_id_w_MUON_EFF_SYS__1up,"muo_id_w_MUON_EFF_SYS__1up/F");
  tr->Branch("muo_id_w_MUON_EFF_SYS__1down",&muo_id_w_MUON_EFF_SYS__1down,"muo_id_w_MUON_EFF_SYS__1down/F");
  tr->Branch("muo_iso_w_MUON_ISO_STAT__1up",&muo_iso_w_MUON_ISO_STAT__1up,"muo_iso_w_MUON_ISO_STAT__1up/F");
  tr->Branch("muo_iso_w_MUON_ISO_STAT__1down",&muo_iso_w_MUON_ISO_STAT__1down,"muo_iso_w_MUON_ISO_STAT__1down/F");
  tr->Branch("muo_iso_w_MUON_ISO_SYS__1up",&muo_iso_w_MUON_ISO_SYS__1up,"muo_iso_w_MUON_ISO_SYS__1up/F");
  tr->Branch("muo_iso_w_MUON_ISO_SYS__1down",&muo_iso_w_MUON_ISO_SYS__1down,"muo_iso_w_MUON_ISO_SYS__1down/F");
  

  tr->Branch("eve_btag_w",&eve_btag_w,"eve_btag_w/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_B_0__1down",&eve_btag_w_FT_EFF_Eigen_B_0__1down,"eve_btag_w_FT_EFF_Eigen_B_0__1down/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_B_0__1up",&eve_btag_w_FT_EFF_Eigen_B_0__1up,"eve_btag_w_FT_EFF_Eigen_B_0__1up/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_B_1__1down",&eve_btag_w_FT_EFF_Eigen_B_1__1down,"eve_btag_w_FT_EFF_Eigen_B_1__1down/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_B_1__1up",&eve_btag_w_FT_EFF_Eigen_B_1__1up,"eve_btag_w_FT_EFF_Eigen_B_1__1up/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_B_2__1down",&eve_btag_w_FT_EFF_Eigen_B_2__1down,"eve_btag_w_FT_EFF_Eigen_B_2__1down/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_B_2__1up",&eve_btag_w_FT_EFF_Eigen_B_2__1up,"eve_btag_w_FT_EFF_Eigen_B_2__1up/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_C_0__1down",&eve_btag_w_FT_EFF_Eigen_C_0__1down,"eve_btag_w_FT_EFF_Eigen_C_0__1down/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_C_0__1up",&eve_btag_w_FT_EFF_Eigen_C_0__1up,"eve_btag_w_FT_EFF_Eigen_C_0__1up/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_C_1__1down",&eve_btag_w_FT_EFF_Eigen_C_1__1down,"eve_btag_w_FT_EFF_Eigen_C_1__1down/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_C_1__1up",&eve_btag_w_FT_EFF_Eigen_C_1__1up,"eve_btag_w_FT_EFF_Eigen_C_1__1up/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_C_2__1down",&eve_btag_w_FT_EFF_Eigen_C_2__1down,"eve_btag_w_FT_EFF_Eigen_C_2__1down/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_C_2__1up",&eve_btag_w_FT_EFF_Eigen_C_2__1up,"eve_btag_w_FT_EFF_Eigen_C_2__1up/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_C_3__1down",&eve_btag_w_FT_EFF_Eigen_C_3__1down,"eve_btag_w_FT_EFF_Eigen_C_3__1down/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_C_3__1up",&eve_btag_w_FT_EFF_Eigen_C_3__1up,"eve_btag_w_FT_EFF_Eigen_C_3__1up/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_Light_0__1down",&eve_btag_w_FT_EFF_Eigen_Light_0__1down,"eve_btag_w_FT_EFF_Eigen_Light_0__1down/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_Light_0__1up",&eve_btag_w_FT_EFF_Eigen_Light_0__1up,"eve_btag_w_FT_EFF_Eigen_Light_0__1up/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_Light_1__1down",&eve_btag_w_FT_EFF_Eigen_Light_1__1down,"eve_btag_w_FT_EFF_Eigen_Light_1__1down/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_Light_1__1up",&eve_btag_w_FT_EFF_Eigen_Light_1__1up,"eve_btag_w_FT_EFF_Eigen_Light_1__1up/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_Light_2__1down",&eve_btag_w_FT_EFF_Eigen_Light_2__1down,"eve_btag_w_FT_EFF_Eigen_Light_2__1down/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_Light_2__1up",&eve_btag_w_FT_EFF_Eigen_Light_2__1up,"eve_btag_w_FT_EFF_Eigen_Light_2__1up/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_Light_3__1down",&eve_btag_w_FT_EFF_Eigen_Light_3__1down,"eve_btag_w_FT_EFF_Eigen_Light_3__1down/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_Light_3__1up",&eve_btag_w_FT_EFF_Eigen_Light_3__1up,"eve_btag_w_FT_EFF_Eigen_Light_3__1up/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_Light_4__1down",&eve_btag_w_FT_EFF_Eigen_Light_4__1down,"eve_btag_w_FT_EFF_Eigen_Light_4__1down/F");
  tr->Branch("eve_btag_w_FT_EFF_Eigen_Light_4__1up",&eve_btag_w_FT_EFF_Eigen_Light_4__1up,"eve_btag_w_FT_EFF_Eigen_Light_4__1up/F");       
  tr->Branch("eve_btag_w_FT_EFF_extrapolation_from_charm__1down",&eve_btag_w_FT_EFF_extrapolation_from_charm__1down,"eve_btag_w_FT_EFF_extrapolation_from_charm__1down/F");
  tr->Branch("eve_btag_w_FT_EFF_extrapolation_from_charm__1up",&eve_btag_w_FT_EFF_extrapolation_from_charm__1up,"eve_btag_w_FT_EFF_extrapolation_from_charm__1up/F");
  tr->Branch("eve_btag_w_FT_EFF_extrapolation__1down",&eve_btag_w_FT_EFF_extrapolation__1down,"eve_btag_w_FT_EFF_extrapolation__1down/F");
  tr->Branch("eve_btag_w_FT_EFF_extrapolation__1up",&eve_btag_w_FT_EFF_extrapolation__1up,"eve_btag_w_FT_EFF_extrapolation__1up/F");



  // vertices
  tr->Branch("nvtx",&nvtx,"nvtx/I");
  tr->Branch("vtx_x",vtx_x,"vtx_x[nvtx]/F");
  tr->Branch("vtx_y",vtx_y,"vtx_y[nvtx]/F");
  tr->Branch("vtx_z",vtx_z,"vtx_z[nvtx]/F");

  //met
  tr->Branch("met_sumet",&met_sumet,"met_sumet/F");
  tr->Branch("met_pt",&met_pt,"met_pt/F");
  tr->Branch("met_eta",&met_eta,"met_eta/F");
  tr->Branch("met_phi",&met_phi,"met_phi/F");
  tr->Branch("met_px",&met_px,"met_px/F");
  tr->Branch("met_py",&met_py,"met_py/F");

  // electrons
  tr->Branch("nele",&nele,"nele/I");
  tr->Branch("ele_charge",ele_charge,"ele_charge[nele]/I");
  tr->Branch("ele_E",ele_E,"ele_E[nele]/F");
  tr->Branch("ele_m", ele_m, "ele_m[nele]/F");
  tr->Branch("ele_pt", ele_pt, "ele_pt[nele]/F");
  tr->Branch("ele_phi",ele_phi,"ele_phi[nele]/F");
  tr->Branch("ele_eta",ele_eta,"ele_eta[nele]/F");
  tr->Branch("ele_cluseta",ele_cluseta,"ele_cluseta[nele]/F");
  tr->Branch("ele_z0",ele_z0,"ele_z0[nele]/F");
  tr->Branch("ele_z0sintheta",ele_z0sintheta,"ele_z0sintheta[nele]/F");
  tr->Branch("ele_d0",ele_d0,"ele_d0[nele]/F");
  tr->Branch("ele_d0sig",ele_d0sig,"ele_d0sig[nele]/F");
  tr->Branch("ele_sharedtrk",ele_sharedtrk,"ele_sharedtrk[nele]/I");
  tr->Branch("ele_isFromPV",ele_isFromPV,"ele_isFromPV[nele]/I");
  tr->Branch("ele_ptiso",ele_ptiso,"ele_ptiso[nele]/F");
  tr->Branch("ele_etiso",ele_etiso,"ele_etiso[nele]/F");
  tr->Branch("ele_idVHl",ele_idVHl,"ele_idVHl[nele]/I");
  tr->Branch("ele_idZHl",ele_idZHl,"ele_idZHl[nele]/I");
  tr->Branch("ele_idWHl",ele_idWHl,"ele_idWHl[nele]/I");
  tr->Branch("ele_passor",ele_passor,"ele_passor[nele]/I");
  tr->Branch("ele_trkpt", ele_trkpt, "ele_trkpt[nele]/F");
  
  tr->Branch("ele_isoWP1",ele_isoWP1,"ele_isoWP1[nele]/I");
  tr->Branch("ele_isoWP2",ele_isoWP2,"ele_isoWP2[nele]/I");
  tr->Branch("ele_isoWP3",ele_isoWP3,"ele_isoWP3[nele]/I");
  tr->Branch("ele_isoWP4",ele_isoWP4,"ele_isoWP4[nele]/I");
  tr->Branch("ele_isoWP5",ele_isoWP5,"ele_isoWP5[nele]/I");
  tr->Branch("ele_isoWP6",ele_isoWP6,"ele_isoWP6[nele]/I");

  tr->Branch("ele_isLooseIDElectron",ele_isLooseIDElectron,"ele_isLooseIDElectron[nele]/I");
  tr->Branch("ele_isMediumIDElectron",ele_isMediumIDElectron,"ele_isMediumIDElectron[nele]/I");
  tr->Branch("ele_isTightIDElectron",ele_isTightIDElectron,"ele_isTightIDElectron[nele]/I");

  tr->Branch("ele_matchHLT_e24_lhmedium_L1EM18VH",ele_matchHLT_e24_lhmedium_L1EM18VH,"ele_matchHLT_e24_lhmedium_L1EM18VH[nele]/I");
  tr->Branch("ele_matchHLT_e24_lhmedium_L1EM20VH",ele_matchHLT_e24_lhmedium_L1EM20VH,"ele_matchHLT_e24_lhmedium_L1EM20VH[nele]/I");
  tr->Branch("ele_matchHLT_e26_lhtight_iloose",ele_matchHLT_e26_lhtight_iloose,"ele_matchHLT_e26_lhtight_iloose[nele]/I");
  tr->Branch("ele_matchHLT_e60_lhmedium",ele_matchHLT_e60_lhmedium,"ele_matchHLT_e60_lhmedium[nele]/I");
  tr->Branch("ele_matchHLT_e120_lhloose",ele_matchHLT_e120_lhloose,"ele_matchHLT_e120_lhloose[nele]/I");

  // muons
  tr->Branch("nmuo",&nmuo,"nmuo/I");
  tr->Branch("muo_charge",muo_charge,"muo_charge[nmuo]/I");
  tr->Branch("muo_E",muo_E,"muo_E[nmuo]/F");
  tr->Branch("muo_m",muo_m,"muo_m[nmuo]/F");
  tr->Branch("muo_pt" ,muo_pt, "muo_pt[nmuo]/F");
  tr->Branch("muo_phi",muo_phi,"muo_phi[nmuo]/F");
  tr->Branch("muo_eta",muo_eta,"muo_eta[nmuo]/F");
  tr->Branch("muo_z0",muo_z0,"muo_z0[nmuo]/F");
  tr->Branch("muo_z0sintheta",muo_z0sintheta,"muo_z0sintheta[nmuo]/F");
  tr->Branch("muo_d0",muo_d0,"muo_d0[nmuo]/F");
  tr->Branch("muo_d0sig",muo_d0sig,"muo_d0sig[nmuo]/F");
  tr->Branch("muo_isFromPV",muo_isFromPV,"muo_isFromPV[nmuo]/I");
  tr->Branch("muo_ptiso",muo_ptiso,"muo_ptiso[nmuo]/F");
  tr->Branch("muo_etiso",muo_etiso,"muo_etiso[nmuo]/F");
  tr->Branch("muo_quality",muo_quality,"muo_quality[nmuo]/I");
  tr->Branch("muo_idVHl",muo_idVHl,"muo_idVHl[nmuo]/I");
  tr->Branch("muo_idZH",muo_idZH,"muo_idZH[nmuo]/I");
  tr->Branch("muo_idWH",muo_idWH,"muo_idWH[nmuo]/I");
  tr->Branch("muo_passor",muo_passor,"muo_passor[nmuo]/I");
  tr->Branch("muo_trkpt" ,muo_trkpt, "muo_trkpt[nmuo]/F");

  tr->Branch("muo_isoWP1",muo_isoWP1,"muo_isoWP1[nmuo]/I");
  tr->Branch("muo_isoWP2",muo_isoWP2,"muo_isoWP2[nmuo]/I");
  tr->Branch("muo_isoWP3",muo_isoWP3,"muo_isoWP3[nmuo]/I");
  tr->Branch("muo_isoWP4",muo_isoWP4,"muo_isoWP4[nmuo]/I");
  tr->Branch("muo_isoWP5",muo_isoWP5,"muo_isoWP5[nmuo]/I");
  tr->Branch("muo_isoWP6",muo_isoWP6,"muo_isoWP6[nmuo]/I");

  tr->Branch("muo_isLooseIDMuon",muo_isLooseIDMuon,"muo_isLooseIDMuon[nmuo]/I");
  tr->Branch("muo_isMediumIDMuon",muo_isMediumIDMuon,"muo_isMediumIDMuon[nmuo]/I");
  tr->Branch("muo_isTightIDMuon",muo_isTightIDMuon,"muo_isTightIDMuon[nmuo]/I");

  tr->Branch("muo_matchHLT_mu20_iloose_L1MU15",muo_matchHLT_mu20_iloose_L1MU15,"muo_matchHLT_mu20_iloose_L1MU15[nmuo]/I");
  tr->Branch("muo_matchHLT_mu26_imedium",muo_matchHLT_mu26_imedium,"muo_matchHLT_mu26_imedium[nmuo]/I");
  tr->Branch("muo_matchHLT_mu50",muo_matchHLT_mu50,"muo_matchHLT_mu50[nmuo]/I");


  // taus
  tr->Branch("ntau",&ntau,"ntau/I");
  tr->Branch("tau_charge",tau_charge,"tau_charge[ntau]/I");
  tr->Branch("tau_E",  tau_E,  "tau_E[ntau]/F");
  tr->Branch("tau_m", tau_m, "tau_m[ntau]/F");
  tr->Branch("tau_pt", tau_pt, "tau_pt[ntau]/F");
  tr->Branch("tau_phi",tau_phi,"tau_phi[ntau]/F");
  tr->Branch("tau_eta",tau_eta,"tau_eta[ntau]/F");
  tr->Branch("tau_passor",tau_passor,"tau_passor[ntau]/I");

  // jets
  tr->Branch("njet",&njet,"njet/I");
  tr->Branch("jet_charge",jet_charge,"jet_charge[njet]/I");
  tr->Branch("jet_E",jet_E,"jet_E[njet]/F");
  tr->Branch("jet_m",jet_m,"jet_m[njet]/F");
  tr->Branch("jet_pt", jet_pt,"jet_pt[njet]/F");
  tr->Branch("jet_phi",jet_phi,"jet_phi[njet]/F");
  tr->Branch("jet_eta",jet_eta,"jet_eta[njet]/F");

  tr->Branch("jet_ntrk",jet_ntrk,"jet_ntrk[njet]/I");
  tr->Branch("jet_NumTrkPt500PV",jet_NumTrkPt500PV,"jet_NumTrkPt500PV[njet]/I");
  tr->Branch("jet_SumPtTrkPt500PV",jet_SumPtTrkPt500PV,"jet_SumPtTrkPt500PV[njet]/F");
  tr->Branch("jet_massvx",jet_massvx,"jet_massvx[njet]/F");
  tr->Branch("jet_jvf",jet_jvf,"jet_jvf[njet]/F");
  tr->Branch("jet_jvt",jet_jvt,"jet_jvt[njet]/F");
  tr->Branch("jet_normdist",jet_normdist,"jet_normdist[njet]/F");

  tr->Branch("jet_good",jet_good,"jet_good[njet]/I");
  tr->Branch("jet_passor",jet_passor,"jet_passor[njet]/I");

  tr->Branch("jet_btag",jet_btag,"jet_btag[njet]/I");
  tr->Branch("jet_mv1",jet_mv1,"jet_mv1[njet]/F");
  tr->Branch("jet_sv1ip3d",jet_sv1ip3d,"jet_sv1ip3d[njet]/F");
  tr->Branch("jet_mv2c00",jet_mv2c00,"jet_mv2c00[njet]/F");
  tr->Branch("jet_mv2c10",jet_mv2c10,"jet_mv2c10[njet]/F");
  tr->Branch("jet_mv2c20",jet_mv2c20,"jet_mv2c20[njet]/F");
  tr->Branch("jet_mvb",jet_mvb,"jet_mvb[njet]/F");

  // fat jets
  tr->Branch("nfatjet",&nfatjet,"nfatjet/I");
  tr->Branch("fatjet_m",fatjet_m,"fatjet_m[nfatjet]/F");
  tr->Branch("fatjet_pt", fatjet_pt,"fatjet_pt[nfatjet]/F");
  tr->Branch("fatjet_phi",fatjet_phi,"fatjet_phi[nfatjet]/F");
  tr->Branch("fatjet_eta",fatjet_eta,"fatjet_eta[nfatjet]/F");
  tr->Branch("fatjet_ltrkjet",fatjet_ltrkjet,"fatjet_ltrkjet[nfatjet]/I");
  tr->Branch("fatjet_sltrkjet",fatjet_sltrkjet,"fatjet_sltrkjet[nfatjet]/I");
  tr->Branch("fatjet_ntrkjet",fatjet_ntrkjet,"fatjet_ntrkjet[nfatjet]/I");
  tr->Branch("fatjet_ntrkjetb",fatjet_ntrkjetb,"fatjet_ntrkjetb[nfatjet]/I");
  tr->Branch("fatjet_XbbL2b",fatjet_XbbL2b,"fatjet_XbbL2b[nfatjet]/I");
  tr->Branch("fatjet_XbbL2bmH",fatjet_XbbL2bmH,"fatjet_XbbL2bmH[nfatjet]/I");
  tr->Branch("fatjet_XbbnTag",fatjet_XbbnTag,"fatjet_XbbnTag[nfatjet]/I");
  tr->Branch("fatjet_m_cor",fatjet_m_cor,"fatjet_m_cor[nfatjet]/F");
  tr->Branch("fatjet_pt_cor", fatjet_pt_cor,"fatjet_pt_cor[nfatjet]/F");
  tr->Branch("fatjet_phi_cor",fatjet_phi_cor,"fatjet_phi_cor[nfatjet]/F");
  tr->Branch("fatjet_eta_cor",fatjet_eta_cor,"fatjet_eta_cor[nfatjet]/F");

  // fat jets
  tr->Branch("ntrkjet",&ntrkjet,"ntrkjet/I");
  tr->Branch("trkjet_m",trkjet_m,"trkjet_m[ntrkjet]/F");
  tr->Branch("trkjet_pt", trkjet_pt,"trkjet_pt[ntrkjet]/F");
  tr->Branch("trkjet_phi",trkjet_phi,"trkjet_phi[ntrkjet]/F");
  tr->Branch("trkjet_eta",trkjet_eta,"trkjet_eta[ntrkjet]/F");
  tr->Branch("trkjet_numConstituents",trkjet_numConstituents,"trkjet_numConstituents[ntrkjet]/I");
  tr->Branch("trkjet_btag",trkjet_btag,"trkjet_btag[ntrkjet]/I");
  tr->Branch("trkjet_sv1ip3d",trkjet_sv1ip3d,"trkjet_sv1ip3d[ntrkjet]/F");
  tr->Branch("trkjet_mv2c00",trkjet_mv2c00,"trkjet_mv2c00[ntrkjet]/F");
  tr->Branch("trkjet_mv2c10",trkjet_mv2c10,"trkjet_mv2c10[ntrkjet]/F");
  tr->Branch("trkjet_mv2c20",trkjet_mv2c20,"trkjet_mv2c20[ntrkjet]/F");

  //electron pairs
  tr->Branch("nee",&nee,"nee/I");
  tr->Branch("ee_charge",ee_charge,"ee_charge[nee]/I");
  tr->Branch("ee_m",ee_m,"ee_m[nee]/F");
  tr->Branch("ee_E",ee_E,"ee_E[nee]/F");
  tr->Branch("ee_p" ,ee_p, "ee_p[nee]/F");
  tr->Branch("ee_pt" ,ee_pt, "ee_pt[nee]/F");
  tr->Branch("ee_phi",ee_phi,"ee_phi[nee]/F");
  tr->Branch("ee_eta",ee_eta,"ee_eta[nee]/F");
  tr->Branch("ee_dR",ee_dR,"ee_dR[nee]/F");
  tr->Branch("ee_leg1",ee_leg1,"ee_leg1[nee]/I");
  tr->Branch("ee_leg2",ee_leg2,"ee_leg2[nee]/I");

  //muon pairs
  tr->Branch("nmm",&nmm,"nmm/I");
  tr->Branch("mm_charge",mm_charge,"mm_charge[nmm]/I");
  tr->Branch("mm_m",mm_m,"mm_m[nmm]/F");
  tr->Branch("mm_E",mm_E,"mm_E[nmm]/F");
  tr->Branch("mm_p" ,mm_p, "mm_p[nmm]/F");
  tr->Branch("mm_pt" ,mm_pt, "mm_pt[nmm]/F");
  tr->Branch("mm_phi",mm_phi,"mm_phi[nmm]/F");
  tr->Branch("mm_eta",mm_eta,"mm_eta[nmm]/F");
  tr->Branch("mm_dR",mm_dR,"mm_dR[nmm]/F");
  tr->Branch("mm_leg1",mm_leg1,"mm_leg1[nmm]/I");
  tr->Branch("mm_leg2",mm_leg2,"mm_leg2[nmm]/I");

  //e-mu pairs
  tr->Branch("nem",&nem,"nem/I");
  tr->Branch("em_charge",em_charge,"em_charge[nem]/I");
  tr->Branch("em_m",em_m,"em_m[nem]/F");
  tr->Branch("em_E",em_E,"em_E[nem]/F");
  tr->Branch("em_p" ,em_p, "em_p[nem]/F");
  tr->Branch("em_pt" ,em_pt, "em_pt[nem]/F");
  tr->Branch("em_phi",em_phi,"em_phi[nem]/F");
  tr->Branch("em_eta",em_eta,"em_eta[nem]/F");
  tr->Branch("em_dR",em_dR,"em_dR[nem]/F");
  tr->Branch("em_leg1",em_leg1,"em_leg1[nem]/I");
  tr->Branch("em_leg2",em_leg2,"em_leg2[nem]/I");


  //tau pairs
  tr->Branch("ntt",&ntt,"ntt/I");
  tr->Branch("tt_charge",tt_charge,"tt_charge[ntt]/I");
  tr->Branch("tt_m",tt_m,"tt_m[ntt]/F");
  tr->Branch("tt_E",tt_E,"tt_E[ntt]/F");
  tr->Branch("tt_p" ,tt_p, "tt_p[ntt]/F");
  tr->Branch("tt_pt" ,tt_pt, "tt_pt[ntt]/F");
  tr->Branch("tt_phi",tt_phi,"tt_phi[ntt]/F");
  tr->Branch("tt_eta",tt_eta,"tt_eta[ntt]/F");
  tr->Branch("tt_dR",tt_dR,"tt_dR[ntt]/F");
  tr->Branch("tt_leg1",tt_leg1,"tt_leg1[ntt]/I");
  tr->Branch("tt_leg2",tt_leg2,"tt_leg2[ntt]/I");


  //jet pairs
  tr->Branch("njj",&njj,"njj/I");
  tr->Branch("jj_charge",jj_charge,"jj_charge[njj]/I");
  tr->Branch("jj_m",jj_m,"jj_m[njj]/F");
  tr->Branch("jj_E",jj_E,"jj_E[njj]/F");
  tr->Branch("jj_p" ,jj_pt, "jj_p[njj]/F");
  tr->Branch("jj_pt" ,jj_pt, "jj_pt[njj]/F");
  tr->Branch("jj_phi",jj_phi,"jj_phi[njj]/F");
  tr->Branch("jj_eta",jj_eta,"jj_eta[njj]/F");
  tr->Branch("jj_dR",jj_dR,"jj_dR[njj]/F");
  tr->Branch("jj_leg1",jj_leg1,"jj_leg1[njj]/I");
  tr->Branch("jj_leg2",jj_leg2,"jj_leg2[njj]/I");


  //e-Met
  tr->Branch("neMet",&neMet,"neMet/I");
  tr->Branch("eMet_charge",eMet_charge,"eMet_charge[neMet]/I");
  tr->Branch("eMet_m",eMet_m,"eMet_m[neMet]/F");
  tr->Branch("eMet_E",eMet_E,"eMet_E[neMet]/F");
  tr->Branch("eMet_p" ,eMet_p, "eMet_p[neMet]/F");
  tr->Branch("eMet_pt" ,eMet_pt, "eMet_pt[neMet]/F");
  tr->Branch("eMet_phi",eMet_phi,"eMet_phi[neMet]/F");
  tr->Branch("eMet_eta",eMet_eta,"eMet_eta[neMet]/F");
  tr->Branch("eMet_mT",eMet_mT,"eMet_mT[neMet]/F");
  tr->Branch("eMet_ele",eMet_ele,"eMet_ele[neMet]/I");

  //mu-Met
  tr->Branch("nmuMet",&nmuMet,"nmuMet/I");
  tr->Branch("muMet_charge",muMet_charge,"muMet_charge[nmuMet]/I");
  tr->Branch("muMet_m",muMet_m,"muMet_m[nmuMet]/F");
  tr->Branch("muMet_E",muMet_E,"muMet_E[nmuMet]/F");
  tr->Branch("muMet_p" ,muMet_p, "muMet_p[nmuMet]/F");
  tr->Branch("muMet_pt" ,muMet_pt, "muMet_pt[nmuMet]/F");
  tr->Branch("muMet_phi",muMet_phi,"muMet_phi[nmuMet]/F");
  tr->Branch("muMet_eta",muMet_eta,"muMet_eta[nmuMet]/F");
  tr->Branch("muMet_mT",muMet_mT,"muMet_mT[nmuMet]/F");
  tr->Branch("muMet_mu",muMet_mu,"muMet_mu[nmuMet]/I");


}



