#include "TupleMaker_VHbbResonance/ZHeeJTupleMaker.h"

// this is needed to distribute the algorithm to the workers
ClassImp(ZHeeJTupleMaker)


ZHeeJTupleMaker :: ZHeeJTupleMaker ()
{
  ///Do not use this constructor
}

ZHeeJTupleMaker :: ZHeeJTupleMaker (std::string configPath) :
  ZHllJTupleMaker(configPath),
  m_tuple(0)
{
}

EL::StatusCode ZHeeJTupleMaker :: initialize ()
{
  if( m_debug) std::cout<<" ZHeeJTupleMaker :: initialize "<<std::endl;

  //this is the top level, ntuple must be defined here
  m_tuple = new ZHeeJTuple();
  ZHllJTupleMaker::m_tuple = m_tuple;
    
  return ZHllJTupleMaker::initialize(); 
}



EL::StatusCode ZHeeJTupleMaker :: processEEJ(){
  //Note: Muon array must be filled prior
  m_tuple->nvh=0;

  int ivh=0;
  for(int z = 0; z < m_tuple->nee ; z++){
    for(int h = 0; h < m_tuple->nfatjet ; h++){
      if( ivh < VHTuple::MAXVH) {

	TLorentzVector zp4;
	zp4.SetPtEtaPhiM(m_tuple->ee_pt[z],m_tuple->ee_eta[z],m_tuple->ee_phi[z],m_tuple->ee_m[z]);
	TLorentzVector hp4;
	hp4.SetPtEtaPhiM(m_tuple->fatjet_pt_cor[h],m_tuple->fatjet_eta_cor[h],m_tuple->fatjet_phi_cor[h],m_tuple->fatjet_m_cor[h]);
	TLorentzVector vhp4 = zp4 + hp4;

	m_tuple->vh_charge[ivh] = 0;//m_tuple->ee_charge[z] + m_tuple->fatjet_charge[h];
	m_tuple->vh_E[ivh] = vhp4.E();
	m_tuple->vh_p[ivh] = vhp4.P();
	m_tuple->vh_pt[ivh] = vhp4.Pt(); 
	m_tuple->vh_phi[ivh] = vhp4.Phi();
	m_tuple->vh_eta[ivh] = vhp4.Eta();
	m_tuple->vh_dR[ivh] = zp4.DeltaR(hp4);
          
	m_tuple->vh_v[ivh] = z;
	m_tuple->vh_h[ivh] = h;

	m_tuple->vh_dPhiHMET[ivh]= m_tuple->fatjet_phi[h] - m_tuple->met_phi;
	
	TLorentzVector hp4Corr;
        hp4Corr.SetPtEtaPhiM(m_tuple->fatjet_pt_cor[h] * m_HMass / m_tuple->fatjet_m_cor[h], m_tuple->fatjet_eta_cor[h],m_tuple->fatjet_phi_cor[h], m_HMass);
	TLorentzVector Zp4Corr;
	Zp4Corr.SetPtEtaPhiM(m_tuple->ee_pt[z] * m_ZMass / m_tuple->ee_m[z] ,m_tuple->ee_eta[z],m_tuple->ee_phi[z], m_ZMass);


	m_tuple->vh_m[ivh] = vhp4.M();
	m_tuple->vh_m1[ivh] = vhp4.M();//uncorrected          
	m_tuple->vh_m2[ivh] = (zp4 + hp4Corr).M();//only Higgs corrected
	m_tuple->vh_m3[ivh] = (Zp4Corr + hp4).M();//only Z corrected
	m_tuple->vh_m4[ivh] = (Zp4Corr + hp4Corr).M();//both corrected

	++ivh;
      }else{
	Info( "processEEJ()" , " Trying to build too many ee-fatjet pairs, max is %i" , VHTuple::MAXVH);
      }
    }
  }
  m_tuple->nvh=ivh;
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ZHeeJTupleMaker :: processEvent(){
  
  if( ZHllJTupleMaker :: processEvent()  == EL::StatusCode::FAILURE) {
    //this should never happen, there are no selection applied in base classes
    cout<<"Failed ZHllJTupleMaker :: processEvent() !"<<endl;
    return EL::StatusCode::FAILURE;
  }


  if(m_debug && m_currentSyst.compare("Nominal") ==0 ){
    PrintEventParticles();
  }

  
  //remove events without a di-lepton
  if(m_tuple->nee==0) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_NoDiLepton");


  //Select the pair with the highest sum pT
  int pair=-1;
  for (int z = 0; z < m_tuple->nee ; z++){
    if( pair == -1 || 
	(m_tuple->ele_pt[m_tuple->ee_leg1[z]] + m_tuple->ele_pt[m_tuple->ee_leg2[z]]) 
	> (m_tuple->ele_pt[m_tuple->ee_leg1[pair]] + m_tuple->ele_pt[m_tuple->ee_leg2[pair]]) )
      pair = z;
  }
  if(pair==-1){//should never happen
    cout<<"Failed finding leading pair !"<<endl;
    return EL::StatusCode::FAILURE;
  }

  const int l1=m_tuple->ee_leg1[pair];//this is already the leading
  const int l2=m_tuple->ee_leg2[pair];

  l1p4.SetPtEtaPhiM(m_tuple->ele_pt[l1],
  		    m_tuple->ele_eta[l1],
  		    m_tuple->ele_phi[l1],
  		    m_tuple->ele_m[l1]);
  l2p4.SetPtEtaPhiM(m_tuple->ele_pt[l2],
  		    m_tuple->ele_eta[l2],
  		    m_tuple->ele_phi[l2],
  		    m_tuple->ele_m[l2]);


  ///Preselection
  if( m_tuple->ele_idVHl[l1] ==0 && m_tuple->ele_idVHl[l2] ==0 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonPresel");

  ///Leading Lepton
  if( fabs(m_tuple->ele_cluseta[l1] ) > 2.47 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonEta");

  //if( 1.37 < fabs(m_tuple->ele_cluseta[l1])  && fabs(m_tuple->ele_cluseta[l1]) < 1.52 ) return EL::StatusCode::FAILURE;
  //incrementCounter("eventCounter_leadleptonCrackVeto");
  
  if( m_tuple->ele_pt[l1] < 25000 )return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonPt");
  
  if( m_tuple->ele_isLooseIDElectron[l1] == 0 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonID");

  if( m_tuple->ele_isFromPV[l1] != 1 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonPV");
  

  ///////////Subleading lepton
  if( fabs(m_tuple->ele_cluseta[l2] ) > 2.47 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_subleadleptonEta");
 
  //if( 1.37 < fabs(m_tuple->ele_cluseta[l2])  && fabs(m_tuple->ele_cluseta[l2]) < 1.52 ) return EL::StatusCode::FAILURE;
  //incrementCounter("eventCounter_subleadleptonCrackVeto");
 
  if( m_tuple->ele_pt[l2] < 25000 )return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_subleadleptonPt");

  if( m_tuple->ele_isLooseIDElectron[l2] == 0 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_subleadleptonID");

  if( m_tuple->ele_isFromPV[l2] != 1 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_subleadleptonPV");


  ///////Isolation
  if( leptonIsoOption==1 && (m_tuple->ele_isoWP1[l1] != 1 || m_tuple->ele_isoWP1[l2] != 1  ))  return EL::StatusCode::FAILURE;
  if( leptonIsoOption==2 && (m_tuple->ele_isoWP2[l1] != 1 || m_tuple->ele_isoWP2[l2] != 1  ))  return EL::StatusCode::FAILURE;
  if( leptonIsoOption==3 && (m_tuple->ele_isoWP3[l1] != 1 || m_tuple->ele_isoWP3[l2] != 1  ))  return EL::StatusCode::FAILURE;
  if( leptonIsoOption==4 && (m_tuple->ele_isoWP4[l1] != 1 || m_tuple->ele_isoWP4[l2] != 1  ))  return EL::StatusCode::FAILURE;
  if( leptonIsoOption==5 && (m_tuple->ele_isoWP5[l1] != 1 || m_tuple->ele_isoWP5[l2] != 1  ))  return EL::StatusCode::FAILURE;
  if( leptonIsoOption==6 && (m_tuple->ele_isoWP6[l1] != 1 || m_tuple->ele_isoWP6[l2] != 1  ))  return EL::StatusCode::FAILURE;
  if(leptonIsoOption == 100){
    if( (m_tuple->ele_ptiso[l1]-m_tuple->ele_trkpt[l2]*(m_tuple->ee_dR[pair]<0.2))/m_tuple->ele_pt[l1] > 0.1 ) 
      return EL::StatusCode::FAILURE;
    if( (m_tuple->ele_ptiso[l2]-m_tuple->ele_trkpt[l1]*(m_tuple->ee_dR[pair]<0.2))/m_tuple->ele_pt[l2] > 0.1 )
      return EL::StatusCode::FAILURE;
  }
  incrementCounter("eventCounter_leptonIso");



  //Find the b-jets
  if( selectHiggs() == EL::StatusCode::FAILURE ) 
    return EL::StatusCode::FAILURE;


  /////////////////////////////////////////////////
  //////No more selections below here, just flags
  /////////////////////////////////////////////////
  ///Create the llJ candidates and save the indexes 

  if( processEEJ()  == EL::StatusCode::FAILURE) 
    return EL::StatusCode::FAILURE;


  //Find the selected VH : Note leptons and Jets must be pT ordered
  m_tuple->vh=-1;
  for(int i=0;i<m_tuple->nvh;i++){
    if(m_tuple->ee_leg1[m_tuple->vh_v[i]] == l1 &&
       m_tuple->ee_leg2[m_tuple->vh_v[i]] == l2 &&
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
  
  //Data: e24_lhmedium_L1EM20VH || e60_lhmedium || e120_lhloose
  //MC:   e24_lhmedium_L1EM18VH || e60_lhmedium || e120_lhloose
  ///For now use MC trigger on data too
  if(m_triggerPaths.size() > 0 ){    
    if( Props::passHLT_e24_lhmedium_L1EM18VH.get(m_eventInfo)
	|| Props::passHLT_e60_lhmedium.get(m_eventInfo)
	|| Props::passHLT_e120_lhloose.get(m_eventInfo)
	) m_tuple->eve_passTrig=1;
    
    if(  m_tuple->ele_matchHLT_e24_lhmedium_L1EM18VH[l1] == 1 
	|| m_tuple->ele_matchHLT_e24_lhmedium_L1EM18VH[l2] == 1 
	|| m_tuple->ele_matchHLT_e60_lhmedium[l1] == 1  
	|| m_tuple->ele_matchHLT_e60_lhmedium[l2] == 1  
	|| m_tuple->ele_matchHLT_e120_lhloose[l1] == 1  
	|| m_tuple->ele_matchHLT_e120_lhloose[l2] == 1  
	) m_tuple->eve_passTrigMatch=1;
  }

  

  ////Print Event Info
  if(m_debug && m_currentSyst.compare("Nominal") ==0 ){
    PrintEventFlags();//clean event, trigger, trigger match
    PrintVHFlags();
    PrintVH(m_tuple->vh);
  }



  if( m_tuple->eve_isMC ){

    /////////////////////////////////////
    ///For weights find the xAOD objects
    /////////////////////////////////////
    const xAOD::Electron * xAODElectron1 = NULL;
    const xAOD::Electron * xAODElectron2 = NULL;
    for( xAOD::ElectronContainer::const_iterator iter = m_Electrons->begin();iter!= m_Electrons->end(); ++iter) {
      if((*iter)->pt() == m_tuple->ele_pt[l1])
	xAODElectron1 = *iter ;
      if((*iter)->pt() == m_tuple->ele_pt[l2])
	xAODElectron2 = *iter ;
    }
    if(xAODElectron1 == NULL || xAODElectron2 == NULL){
      cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
      return EL::StatusCode::FAILURE;
    } 


    /////////////////////////////////////
    ///Trigger Weight, to use tool need to find the xAOD objects
    /////////////////////////////////////
    ///Trigger tools is crashing
    //double trig_w=1.;
    //m_triggerTool->getTriggerDecision(m_eventInfo,trig_w,xAODElectron1,0,0,0,0,0,276073,"Nominal");
    //m_tuple->eve_trig_w = trig_w;
    m_tuple->eve_trig_w =  Props::trigSFlooseLHIsoLooseTrackOnly.get(xAODElectron1);

    /////////////////////////////////////
    ///reco Weight
    /////////////////////////////////////
    m_tuple->ele_reco_w = Props::effSFReco.get(xAODElectron1)  * Props::effSFReco.get(xAODElectron2);

    /////////////////////////////////////
    ///Id Weight
    /////////////////////////////////////
    m_tuple->ele_id_w = Props::effSFlooseLH.get(xAODElectron1)  * Props::effSFlooseLH.get(xAODElectron2);

    /////////////////////////////////////
    ///Iso Weight, 
    /////////////////////////////////////
    m_tuple->ele_iso_w = Props::effSFIsoLooseTrackOnlyLooseLH.get(xAODElectron1)  * Props::effSFIsoLooseTrackOnlyLooseLH.get(xAODElectron2);
  

    if( m_currentSyst.compare("Nominal") ==0 ){
      
      ///Muon trigger systematics set to nominal for e+e- channel
      m_tuple->eve_trig_w_MUON_EFF_TrigSystUncertainty__1up = m_tuple->eve_trig_w ;
      m_tuple->eve_trig_w_MUON_EFF_TrigSystUncertainty__1down = m_tuple->eve_trig_w ;
      m_tuple->eve_trig_w_MUON_EFF_TrigStatUncertainty__1up = m_tuple->eve_trig_w ;
      m_tuple->eve_trig_w_MUON_EFF_TrigStatUncertainty__1down = m_tuple->eve_trig_w ;

      //m_triggerTool->getTriggerDecision(m_eventInfo,trig_w,xAODElectron1,0,0,0,0,0,276073,"EL_EFF_Trigger_TotalCorrUncertainty__1up");
      //m_tuple->eve_trig_w_EL_EFF_Trigger_TotalCorrUncertainty__1up = trig_w;
      //m_triggerTool->getTriggerDecision(m_eventInfo,trig_w,xAODElectron1,0,0,0,0,0,276073,"EL_EFF_Trigger_TotalCorrUncertainty__1down");
      //m_tuple->eve_trig_w_EL_EFF_Trigger_TotalCorrUncertainty__1down = trig_w;

    
      const xAOD::ElectronContainer * ElectronsVar = 0;
      const xAOD::Electron * xAODElectron1Var = NULL;
      const xAOD::Electron * xAODElectron2Var = NULL;

      //EL_EFF_Trigger_TotalCorrUncertainty__1up
      if( !m_event->retrieve( ElectronsVar, "Electrons___EL_EFF_Trigger_TotalCorrUncertainty__1up" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Electrons___EL_EFF_Trigger_TotalCorrUncertainty__1up");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::ElectronContainer::const_iterator iter = ElectronsVar->begin();iter!= ElectronsVar->end(); ++iter) {
	if((*iter)->pt() == xAODElectron1->pt())
	  xAODElectron1Var = *iter ;
	if((*iter)->pt() == xAODElectron2->pt())
	  xAODElectron2Var = *iter ;
      }
      if(xAODElectron1Var == NULL || xAODElectron2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->eve_trig_w_EL_EFF_Trigger_TotalCorrUncertainty__1up = Props::trigSFlooseLHIsoLooseTrackOnly.get(xAODElectron1Var);
      xAODElectron1Var = NULL; xAODElectron2Var = NULL;
      

      //EL_EFF_Trigger_TotalCorrUncertainty__1down
      if( !m_event->retrieve( ElectronsVar, "Electrons___EL_EFF_Trigger_TotalCorrUncertainty__1down" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Electrons___EL_EFF_Trigger_TotalCorrUncertainty__1down");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::ElectronContainer::const_iterator iter = ElectronsVar->begin();iter!= ElectronsVar->end(); ++iter) {
	if((*iter)->pt() == xAODElectron1->pt())
	  xAODElectron1Var = *iter ;
	if((*iter)->pt() == xAODElectron2->pt())
	  xAODElectron2Var = *iter ;
      }
      if(xAODElectron1Var == NULL || xAODElectron2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->eve_trig_w_EL_EFF_Trigger_TotalCorrUncertainty__1down = Props::trigSFlooseLHIsoLooseTrackOnly.get(xAODElectron1Var);
      xAODElectron1Var = NULL; xAODElectron2Var = NULL;


      ///reco weight EL_EFF_Reco_TotalCorrUncertainty__1up
      if( !m_event->retrieve( ElectronsVar, "Electrons___EL_EFF_Reco_TotalCorrUncertainty__1up" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Electrons___EL_EFF_Reco_TotalCorrUncertainty__1up");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::ElectronContainer::const_iterator iter = ElectronsVar->begin();iter!= ElectronsVar->end(); ++iter) {
	if((*iter)->pt() == xAODElectron1->pt())
	  xAODElectron1Var = *iter ;
	if((*iter)->pt() == xAODElectron2->pt())
	  xAODElectron2Var = *iter ;
      }
      if(xAODElectron1Var == NULL || xAODElectron2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->ele_reco_w_EL_EFF_Reco_TotalCorrUncertainty__1up = Props::effSFReco.get(xAODElectron1Var)  * Props::effSFReco.get(xAODElectron2Var);
      xAODElectron1Var = NULL; xAODElectron2Var = NULL;

      ///reco weight EL_EFF_Reco_TotalCorrUncertainty__1down
      if( !m_event->retrieve( ElectronsVar, "Electrons___EL_EFF_Reco_TotalCorrUncertainty__1down" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Electrons___EL_EFF_Reco_TotalCorrUncertainty__1down");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::ElectronContainer::const_iterator iter = ElectronsVar->begin();iter!= ElectronsVar->end(); ++iter) {
	if((*iter)->pt() == xAODElectron1->pt())
	  xAODElectron1Var = *iter ;
	if((*iter)->pt() == xAODElectron2->pt())
	  xAODElectron2Var = *iter ;
      }
      if(xAODElectron1Var == NULL || xAODElectron2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->ele_reco_w_EL_EFF_Reco_TotalCorrUncertainty__1down = Props::effSFReco.get(xAODElectron1Var)  * Props::effSFReco.get(xAODElectron2Var);
      xAODElectron1Var = NULL; xAODElectron2Var = NULL;


      ///id weight EL_EFF_ID_TotalCorrUncertainty__1up
      if( !m_event->retrieve( ElectronsVar, "Electrons___EL_EFF_ID_TotalCorrUncertainty__1up" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Electrons___EL_EFF_ID_TotalCorrUncertainty__1up");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::ElectronContainer::const_iterator iter = ElectronsVar->begin();iter!= ElectronsVar->end(); ++iter) {
	if((*iter)->pt() == xAODElectron1->pt())
	  xAODElectron1Var = *iter ;
	if((*iter)->pt() == xAODElectron2->pt())
	  xAODElectron2Var = *iter ;
      }
      if(xAODElectron1Var == NULL || xAODElectron2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->ele_id_w_EL_EFF_ID_TotalCorrUncertainty__1up = Props::effSFlooseLH.get(xAODElectron1Var)  * Props::effSFlooseLH.get(xAODElectron2Var);
      xAODElectron1Var = NULL; xAODElectron2Var = NULL;

      ///id weight EL_EFF_ID_TotalCorrUncertainty__1down
      if( !m_event->retrieve( ElectronsVar, "Electrons___EL_EFF_ID_TotalCorrUncertainty__1down" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Electrons___EL_EFF_ID_TotalCorrUncertainty__1down");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::ElectronContainer::const_iterator iter = ElectronsVar->begin();iter!= ElectronsVar->end(); ++iter) {
	if((*iter)->pt() == xAODElectron1->pt())
	  xAODElectron1Var = *iter ;
	if((*iter)->pt() == xAODElectron2->pt())
	  xAODElectron2Var = *iter ;
      }
      if(xAODElectron1Var == NULL || xAODElectron2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->ele_id_w_EL_EFF_ID_TotalCorrUncertainty__1down = Props::effSFlooseLH.get(xAODElectron1Var)  * Props::effSFlooseLH.get(xAODElectron2Var);
      xAODElectron1Var = NULL; xAODElectron2Var = NULL;


      ///iso weight EL_EFF_Iso_TotalCorrUncertainty__1up
      if( !m_event->retrieve( ElectronsVar, "Electrons___EL_EFF_Iso_TotalCorrUncertainty__1up" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Electrons___EL_EFF_Iso_TotalCorrUncertainty__1up");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::ElectronContainer::const_iterator iter = ElectronsVar->begin();iter!= ElectronsVar->end(); ++iter) {
	if((*iter)->pt() == xAODElectron1->pt())
	  xAODElectron1Var = *iter ;
	if((*iter)->pt() == xAODElectron2->pt())
	  xAODElectron2Var = *iter ;
      }
      if(xAODElectron1Var == NULL || xAODElectron2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->ele_iso_w_EL_EFF_Iso_TotalCorrUncertainty__1up = Props::effSFIsoLooseTrackOnlyLooseLH.get(xAODElectron1Var)  * Props::effSFIsoLooseTrackOnlyLooseLH.get(xAODElectron2Var);
      xAODElectron1Var = NULL; xAODElectron2Var = NULL;

      ///iso weight EL_EFF_Iso_TotalCorrUncertainty__1down
      if( !m_event->retrieve( ElectronsVar, "Electrons___EL_EFF_Iso_TotalCorrUncertainty__1down" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Electrons___EL_EFF_Iso_TotalCorrUncertainty__1down");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::ElectronContainer::const_iterator iter = ElectronsVar->begin();iter!= ElectronsVar->end(); ++iter) {
	if((*iter)->pt() == xAODElectron1->pt())
	  xAODElectron1Var = *iter ;
	if((*iter)->pt() == xAODElectron2->pt())
	  xAODElectron2Var = *iter ;
      }
      if(xAODElectron1Var == NULL || xAODElectron2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->ele_iso_w_EL_EFF_Iso_TotalCorrUncertainty__1down = Props::effSFIsoLooseTrackOnlyLooseLH.get(xAODElectron1Var)  * Props::effSFIsoLooseTrackOnlyLooseLH.get(xAODElectron2Var);
      xAODElectron1Var = NULL; xAODElectron2Var = NULL;

    }
  }


  return EL::StatusCode::SUCCESS;
}

