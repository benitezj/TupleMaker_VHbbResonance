#include "TupleMaker_VHbbResonance/ZHemJTupleMaker.h"

// this is needed to distribute the algorithm to the workers
ClassImp(ZHemJTupleMaker)


ZHemJTupleMaker :: ZHemJTupleMaker ()
{
  ///Do not use this constructor
}

ZHemJTupleMaker :: ZHemJTupleMaker (std::string configPath) :
  ZHllJTupleMaker(configPath),
  m_tuple(0)
{
}

EL::StatusCode ZHemJTupleMaker :: initialize ()
{
  //this is the top level, ntuple must be defined here
  m_tuple = new ZHemJTuple();
  ZHllJTupleMaker::m_tuple = m_tuple;
    
  return ZHllJTupleMaker::initialize(); 
}



EL::StatusCode ZHemJTupleMaker :: processEMJ(){
  //Note: Muon array must be filled prior
  m_tuple->nvh=0;

  int ivh=0;
  for(int z = 0; z < m_tuple->nem ; z++){
    for(int h = 0; h < m_tuple->nfatjet ; h++){
      if( ivh < VHTuple::MAXVH) {

	TLorentzVector zp4;
	zp4.SetPtEtaPhiM(m_tuple->em_pt[z],m_tuple->em_eta[z],m_tuple->em_phi[z],m_tuple->em_m[z]);
	TLorentzVector hp4;
	hp4.SetPtEtaPhiM(m_tuple->fatjet_pt_cor[h],m_tuple->fatjet_eta_cor[h],m_tuple->fatjet_phi_cor[h],m_tuple->fatjet_m_cor[h]);
	TLorentzVector vhp4 = zp4 + hp4;

	m_tuple->vh_charge[ivh] = 0;//m_tuple->em_charge[z] + m_tuple->fatjet_charge[h];
	m_tuple->vh_E[ivh] = vhp4.E();
	m_tuple->vh_p[ivh] = vhp4.P();
	m_tuple->vh_pt[ivh] = vhp4.Pt(); 
	m_tuple->vh_phi[ivh] = vhp4.Phi();
	m_tuple->vh_eta[ivh] = vhp4.Eta();
	m_tuple->vh_dR[ivh] = zp4.DeltaR(hp4);
          
	m_tuple->vh_v[ivh] = z;
	m_tuple->vh_h[ivh] = h;

	m_tuple->vh_dPhiHMET[ivh]= m_tuple->fatjet_phi[h] - m_tuple->met_phi;
	
	///Apply Higgs and Z mass constraint
	TLorentzVector hp4Corr;
        hp4Corr.SetPtEtaPhiM(m_tuple->fatjet_pt_cor[h] * m_HMass / m_tuple->fatjet_m_cor[h], m_tuple->fatjet_eta_cor[h],m_tuple->fatjet_phi_cor[h], m_HMass);
	TLorentzVector Zp4Corr;
	Zp4Corr.SetPtEtaPhiM(m_tuple->em_pt[z] * m_ZMass / m_tuple->em_m[z] ,m_tuple->em_eta[z],m_tuple->em_phi[z], m_ZMass);
 
	m_tuple->vh_m[ivh] = vhp4.M();
	m_tuple->vh_m1[ivh] = vhp4.M();//uncorrected          
	m_tuple->vh_m2[ivh] = (zp4 + hp4Corr).M();//only Higgs corrected
	m_tuple->vh_m3[ivh] = (Zp4Corr + hp4).M();//only Z corrected
	m_tuple->vh_m4[ivh] = (Zp4Corr + hp4Corr).M();//both corrected
	

	++ivh;
      }else{
	Info( "processEMJ()" , " Trying to build too many em-fatjet pairs, max is %i" , VHTuple::MAXVH);
      }
    }
  }
  m_tuple->nvh=ivh;
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ZHemJTupleMaker :: processEvent(){
  
  if( ZHllJTupleMaker :: processEvent()  == EL::StatusCode::FAILURE) {
    //this should never happen, there are no selection applied in base classes
    cout<<"Failed ZHllJTupleMaker :: processEvent() !"<<endl;
    return EL::StatusCode::FAILURE;
  }
  
  //remove events without a di-lepton
  if(m_tuple->nem==0) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_NoDiLepton");

  //Select the e-mu pair with the highest sum pT
  int pair=-1;
  for (int z = 0; z < m_tuple->nem ; z++){
    float sumpt= m_tuple->ele_pt[m_tuple->em_leg1[z]] + m_tuple->muo_pt[m_tuple->em_leg2[z]];
    if(pair==-1 ||
       sumpt > m_tuple->ele_pt[m_tuple->em_leg1[pair]] + m_tuple->muo_pt[m_tuple->em_leg2[pair]])
      pair=z;
  }
  if(pair==-1){//should never happen
    cout<<"Failed finding leading e-mu pair !"<<endl;
    return EL::StatusCode::FAILURE;
  }

  const int l1=m_tuple->em_leg1[pair];
  const int l2=m_tuple->em_leg2[pair];
  float l1pt=m_tuple->ele_pt[m_tuple->em_leg1[pair]];
  float l2pt=m_tuple->muo_pt[m_tuple->em_leg2[pair]];

  l1p4.SetPtEtaPhiM(m_tuple->ele_pt[l1],
  		    m_tuple->ele_eta[l1],
  		    m_tuple->ele_phi[l1],
  		    m_tuple->ele_m[l1]);
 
  l2p4.SetPtEtaPhiM(m_tuple->muo_pt[l2],
  		    m_tuple->muo_eta[l2],
  		    m_tuple->muo_phi[l2],
  		    m_tuple->muo_m[l2]);


  ///Preselection
  if( m_tuple->ele_idVHl[l1] ==0 && m_tuple->muo_idVHl[l2] ==0 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonPresel");

  //if electron is the leading lepton
  if(l1pt>l2pt){
    
    ///Selections on the electron
    if(fabs(m_tuple->ele_cluseta[l1]) > 2.47 )  return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_leadleptonEta");
    
    if( m_tuple->ele_pt[l1] < 25000 )  return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_leadleptonPt");
    
    if( m_tuple->ele_isLooseIDElectron[l1] == 0 )  return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_leadleptonID");
    
    if( m_tuple->ele_isFromPV[l1] != 1 ) return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_leadleptonPV");


    //Selections on the Muon
    if( fabs(m_tuple->muo_eta[l2] ) > 2.5 )  return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_subleadleptonEta");

    if( m_tuple->muo_pt[l2] < 25000)  return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_subleadleptonPt");
  
    if( m_tuple->muo_isLooseIDMuon[l2] == 0 )  return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_subleadleptonID");

    if( m_tuple->muo_isFromPV[l2] != 1 ) return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_subleadleptonPV");
    

  }

  //if muon is the leading lepton
  if(l1pt<l2pt){
    
    //selections on the  Muon
    if(fabs(m_tuple->muo_eta[l2] ) > 2.5)  return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_leadleptonEta");
    
    if( m_tuple->muo_pt[l2] < 25000 )  return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_leadleptonPt");

    if( m_tuple->muo_isLooseIDMuon[l2] == 0 )  return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_leadleptonID");

    if( m_tuple->muo_isFromPV[l2] != 1 ) return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_leadleptonPV");


    //Selections on the electron
    if( fabs(m_tuple->ele_cluseta[l1] ) > 2.47 )  return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_subleadleptonEta");

    if( m_tuple->ele_pt[l1] < 25000 )  return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_subleadleptonPt");

    if( m_tuple->ele_isLooseIDElectron[l1] == 0)  return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_subleadleptonID");

    if( m_tuple->ele_isFromPV[l1] != 1 ) return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_subleadleptonPV");

  }


  ///////Isolation
  if( leptonIsoOption==1 && ( m_tuple->ele_isoWP1[l1] != 1 || m_tuple->muo_isoWP1[l2] != 1 ) ) return EL::StatusCode::FAILURE;
  if( leptonIsoOption==2 && ( m_tuple->ele_isoWP2[l1] != 1 || m_tuple->muo_isoWP2[l2] != 1 ) ) return EL::StatusCode::FAILURE;
  if( leptonIsoOption==3 && ( m_tuple->ele_isoWP3[l1] != 1 || m_tuple->muo_isoWP3[l2] != 1 ) ) return EL::StatusCode::FAILURE;
  if( leptonIsoOption==4 && ( m_tuple->ele_isoWP4[l1] != 1 || m_tuple->muo_isoWP4[l2] != 1 ) ) return EL::StatusCode::FAILURE;
  if( leptonIsoOption==5 && ( m_tuple->ele_isoWP5[l1] != 1 || m_tuple->muo_isoWP5[l2] != 1 ) ) return EL::StatusCode::FAILURE;
  if( leptonIsoOption==6 && ( m_tuple->ele_isoWP6[l1] != 1 || m_tuple->muo_isoWP6[l2] != 1 ) ) return EL::StatusCode::FAILURE;

  if(leptonIsoOption == 100){
    //check electron
    if( (m_tuple->ele_ptiso[l1]-m_tuple->muo_trkpt[l2]*(m_tuple->em_dR[pair]<0.2))/m_tuple->ele_pt[l1] > 0.1 )
      return EL::StatusCode::FAILURE;

    //check muon
    if( (m_tuple->muo_ptiso[l2]-m_tuple->ele_trkpt[l1]*(m_tuple->em_dR[pair]<0.2))/m_tuple->muo_pt[l2] > 0.1 ) 
      return EL::StatusCode::FAILURE;
  }
  incrementCounter("eventCounter_LeptonIsolation");

  //Find the b-jets and fill extra jets
  if( selectHiggs() == EL::StatusCode::FAILURE ) 
    return EL::StatusCode::FAILURE;


  /////////////////////////////////////////////////
  //////No more selections below here, just flags
  /////////////////////////////////////////////////

  ///Create the llJ candidates
  if( processEMJ()  == EL::StatusCode::FAILURE) 
    return EL::StatusCode::FAILURE;


  //Find the selected VH : Note leptons and Jets must be pT ordered
  m_tuple->vh=-1;
  for(int i=0;i<m_tuple->nvh;i++){
    if(m_tuple->em_leg1[m_tuple->vh_v[i]] == l1 &&
       m_tuple->em_leg2[m_tuple->vh_v[i]] == l2 &&
       m_tuple->vh_h[i] == higgsJet
       ) 
      m_tuple->vh=i;
  }  
  if( m_tuple->vh==-1) {
    cout<<"processEvent: something is wrong with the VH logic"<<endl; 
    return EL::StatusCode::FAILURE;
  }


  
  ////////////////////////////////////
  ///check Trigger
  ////////////////////////////////////
  m_tuple->eve_passTrig=0;
  m_tuple->eve_passTrigMatch=0;
  
  //no trigger requirement
  if(m_triggerPaths.size()==0){
    m_tuple->eve_passTrig=1;
    m_tuple->eve_passTrigMatch=1;
  }
  
  if(m_triggerPaths.size() > 0 ){
    if(Props::passHLT_e24_lhmedium_L1EM18VH.get(m_eventInfo)  
       || Props::passHLT_e60_lhmedium.get(m_eventInfo) 
       || Props::passHLT_e120_lhloose.get(m_eventInfo)
       || Props::passHLT_mu20_iloose_L1MU15.get(m_eventInfo) 
       || Props::passHLT_mu50.get(m_eventInfo) 
       )  m_tuple->eve_passTrig=1;
    
    if( m_tuple->ele_matchHLT_e24_lhmedium_L1EM18VH[l1] == 1 
	|| m_tuple->ele_matchHLT_e60_lhmedium[l1] == 1  
	|| m_tuple->ele_matchHLT_e120_lhloose[l1] == 1  
	|| m_tuple->muo_matchHLT_mu20_iloose_L1MU15[l2] == 1 
	|| m_tuple->muo_matchHLT_mu50[l2] == 1
	)  m_tuple->eve_passTrigMatch=1;
     
  }



  if( m_tuple->eve_isMC ){
    
    /////////////////////////////////////
    ///Trigger Weight, to use tool need to find the xAOD objects
    /////////////////////////////////////
    const xAOD::Muon * xAODMuon = NULL;
    const xAOD::Electron * xAODElectron = NULL;
    for( xAOD::ElectronContainer::const_iterator iter = m_Electrons->begin();iter!= m_Electrons->end(); ++iter) {
      if((*iter)->pt() == m_tuple->ele_pt[l1])
	xAODElectron = *iter ;
    }
    for( xAOD::MuonContainer::const_iterator iter = m_Muons->begin();iter!= m_Muons->end(); ++iter) {
      if((*iter)->pt() == m_tuple->muo_pt[l2])
	xAODMuon = *iter ;
    }
    if( xAODElectron == NULL || xAODMuon == NULL){
      cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
      return EL::StatusCode::FAILURE;
    } 
  
    if(l1pt<l2pt){
      //Muon is higher pT
 
      double trig_w=1.;
      m_triggerTool->getTriggerDecision(m_eventInfo,trig_w,0,0,xAODMuon,0,0,0,m_randRun,"Nominal");
      m_tuple->eve_trig_w = trig_w;

      //ttv weight (Muons only)
      m_tuple->muo_ttv_w = Props::TTVAEffSF.get(xAODMuon);

      //reco weight (Electrons only)
      m_tuple->ele_reco_w = Props::effSFReco.get(xAODElectron);

      ///Id Weight
      m_tuple->ele_id_w = Props::effSFlooseLH.get(xAODElectron);
      m_tuple->muo_id_w = Props::looseEffSF.get(xAODMuon);
    
      ///Iso Weight
      m_tuple->ele_iso_w = Props::effSFIsoLooseTrackOnlyLooseLH.get(xAODElectron);
      m_tuple->muo_iso_w = Props::looseTrackOnlyIsoSF.get(xAODMuon);



      if(  m_currentSyst.compare("Nominal") ==0 ) {
	//systematics not implemented

	m_tuple->eve_trig_w_EL_EFF_Trigger_TotalCorrUncertainty__1up = trig_w;
	m_tuple->eve_trig_w_EL_EFF_Trigger_TotalCorrUncertainty__1down = trig_w;

      	m_tuple->eve_trig_w_MUON_EFF_TrigSystUncertainty__1up = trig_w;
	m_tuple->eve_trig_w_MUON_EFF_TrigSystUncertainty__1down = trig_w;
	m_tuple->eve_trig_w_MUON_EFF_TrigStatUncertainty__1up = trig_w;
	m_tuple->eve_trig_w_MUON_EFF_TrigStatUncertainty__1down = trig_w;

	m_tuple->ele_reco_w_EL_EFF_Reco_TotalCorrUncertainty__1up= m_tuple->ele_reco_w ;
	m_tuple->ele_reco_w_EL_EFF_Reco_TotalCorrUncertainty__1down= m_tuple->ele_reco_w ;

	m_tuple->ele_id_w_EL_EFF_ID_TotalCorrUncertainty__1up= m_tuple->ele_id_w ;
	m_tuple->ele_id_w_EL_EFF_ID_TotalCorrUncertainty__1down= m_tuple->ele_id_w ;

	m_tuple->ele_iso_w_EL_EFF_Iso_TotalCorrUncertainty__1up= m_tuple->ele_iso_w ;
	m_tuple->ele_iso_w_EL_EFF_Iso_TotalCorrUncertainty__1down= m_tuple->ele_iso_w ;

	m_tuple->muo_ttv_w_MUON_TTVA_STAT__1up=m_tuple->muo_ttv_w ;
	m_tuple->muo_ttv_w_MUON_TTVA_STAT__1down=m_tuple->muo_ttv_w ;
	m_tuple->muo_ttv_w_MUON_TTVA_SYS__1up=m_tuple->muo_ttv_w ;
	m_tuple->muo_ttv_w_MUON_TTVA_SYS__1down=m_tuple->muo_ttv_w ;

	m_tuple->muo_id_w_MUON_EFF_STAT__1up=m_tuple->muo_id_w ;
	m_tuple->muo_id_w_MUON_EFF_STAT__1down=m_tuple->muo_id_w ;
	m_tuple->muo_id_w_MUON_EFF_SYS__1up=m_tuple->muo_id_w ;
	m_tuple->muo_id_w_MUON_EFF_SYS__1down=m_tuple->muo_id_w ;

	m_tuple->muo_iso_w_MUON_ISO_STAT__1up= m_tuple->muo_iso_w ;
	m_tuple->muo_iso_w_MUON_ISO_STAT__1down= m_tuple->muo_iso_w ;
	m_tuple->muo_iso_w_MUON_ISO_SYS__1up= m_tuple->muo_iso_w ;
	m_tuple->muo_iso_w_MUON_ISO_SYS__1down= m_tuple->muo_iso_w ;

      

      }
    }else {
      //Electron is higher pT 
      double trig_w=1.;
      trig_w =  Props::trigSFlooseLHIsoLooseTrackOnly.get(xAODElectron);
      m_tuple->eve_trig_w = trig_w;

      //ttv weight (Muons only)
      m_tuple->muo_ttv_w = Props::TTVAEffSF.get(xAODMuon)  ;

      //reco weight (Electrons only)
      m_tuple->ele_reco_w = Props::effSFReco.get(xAODElectron);

      ///Id Weight
      m_tuple->ele_id_w = Props::effSFlooseLH.get(xAODElectron);
      m_tuple->muo_id_w = Props::looseEffSF.get(xAODMuon);
    
      ///Iso Weight
      m_tuple->ele_iso_w = Props::effSFIsoLooseTrackOnlyLooseLH.get(xAODElectron);
      m_tuple->muo_iso_w = Props::looseTrackOnlyIsoSF.get(xAODMuon);


      if( m_currentSyst.compare("Nominal") ==0 ) {
	//systematics not implemented
	m_tuple->eve_trig_w_MUON_EFF_TrigSystUncertainty__1up = trig_w;
	m_tuple->eve_trig_w_MUON_EFF_TrigSystUncertainty__1down = trig_w;
	m_tuple->eve_trig_w_MUON_EFF_TrigStatUncertainty__1up = trig_w;
	m_tuple->eve_trig_w_MUON_EFF_TrigStatUncertainty__1down = trig_w;

	m_tuple->eve_trig_w_EL_EFF_Trigger_TotalCorrUncertainty__1up = trig_w;
	m_tuple->eve_trig_w_EL_EFF_Trigger_TotalCorrUncertainty__1down = trig_w;

	m_tuple->ele_reco_w_EL_EFF_Reco_TotalCorrUncertainty__1up= m_tuple->ele_reco_w ;
	m_tuple->ele_reco_w_EL_EFF_Reco_TotalCorrUncertainty__1down= m_tuple->ele_reco_w ;

	m_tuple->ele_id_w_EL_EFF_ID_TotalCorrUncertainty__1up= m_tuple->ele_id_w ;
	m_tuple->ele_id_w_EL_EFF_ID_TotalCorrUncertainty__1down= m_tuple->ele_id_w ;

	m_tuple->ele_iso_w_EL_EFF_Iso_TotalCorrUncertainty__1up= m_tuple->ele_iso_w ;
	m_tuple->ele_iso_w_EL_EFF_Iso_TotalCorrUncertainty__1down= m_tuple->ele_iso_w ;

	m_tuple->muo_ttv_w_MUON_TTVA_STAT__1up=m_tuple->muo_ttv_w ;
	m_tuple->muo_ttv_w_MUON_TTVA_STAT__1down=m_tuple->muo_ttv_w ;
	m_tuple->muo_ttv_w_MUON_TTVA_SYS__1up=m_tuple->muo_ttv_w ;
	m_tuple->muo_ttv_w_MUON_TTVA_SYS__1down=m_tuple->muo_ttv_w ;

	m_tuple->muo_id_w_MUON_EFF_STAT__1up=m_tuple->muo_id_w ;
	m_tuple->muo_id_w_MUON_EFF_STAT__1down=m_tuple->muo_id_w ;
	m_tuple->muo_id_w_MUON_EFF_SYS__1up=m_tuple->muo_id_w ;
	m_tuple->muo_id_w_MUON_EFF_SYS__1down=m_tuple->muo_id_w ;

	m_tuple->muo_iso_w_MUON_ISO_STAT__1up= m_tuple->muo_iso_w ;
	m_tuple->muo_iso_w_MUON_ISO_STAT__1down= m_tuple->muo_iso_w ;
	m_tuple->muo_iso_w_MUON_ISO_SYS__1up= m_tuple->muo_iso_w ;
	m_tuple->muo_iso_w_MUON_ISO_SYS__1down= m_tuple->muo_iso_w ;


      }

    }
  }

  return EL::StatusCode::SUCCESS;
}

