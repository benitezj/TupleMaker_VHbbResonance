#include "TupleMaker_VHbbResonance/ZHmmJTupleMaker.h"

// this is needed to distribute the algorithm to the workers
ClassImp(ZHmmJTupleMaker)

ZHmmJTupleMaker :: ZHmmJTupleMaker ()
{
  ///Do not use this constructor
}

ZHmmJTupleMaker :: ZHmmJTupleMaker (std::string configPath) :
  ZHllJTupleMaker(configPath),
  m_tuple(0)
{
}

EL::StatusCode ZHmmJTupleMaker :: initialize ()
{
  //this is the top level, ntuple must be defined here
  m_tuple = new ZHmmJTuple();
  ZHllJTupleMaker::m_tuple = m_tuple;
    
  return ZHllJTupleMaker::initialize(); 
}



EL::StatusCode ZHmmJTupleMaker :: processMMJ(){
  //Note: Muon array must be filled prior
  m_tuple->nvh=0;

  int ivh=0;
  for(int z = 0; z < m_tuple->nmm ; z++){
    for(int h = 0; h < m_tuple->nfatjet ; h++){
      if( ivh < VHTuple::MAXVH) {

	TLorentzVector zp4;
	zp4.SetPtEtaPhiM(m_tuple->mm_pt[z],m_tuple->mm_eta[z],m_tuple->mm_phi[z],m_tuple->mm_m[z]);

	TLorentzVector hp4;
	//hp4.SetPtEtaPhiM(m_tuple->fatjet_pt[h],m_tuple->fatjet_eta[h],m_tuple->fatjet_phi[h],m_tuple->fatjet_m[h]);
	hp4.SetPtEtaPhiM(m_tuple->fatjet_pt_cor[h],m_tuple->fatjet_eta_cor[h],m_tuple->fatjet_phi_cor[h],m_tuple->fatjet_m_cor[h]);

	TLorentzVector vhp4 = zp4 + hp4;

	m_tuple->vh_charge[ivh] = 0;//m_tuple->mm_charge[z] + m_tuple->fatjet_charge[h];
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
	Zp4Corr.SetPtEtaPhiM(m_tuple->mm_pt[z] * m_ZMass / m_tuple->mm_m[z] ,m_tuple->mm_eta[z],m_tuple->mm_phi[z], m_ZMass);

	m_tuple->vh_m[ivh] = (Zp4Corr + hp4).M();//only Z corrected
	m_tuple->vh_m1[ivh] = vhp4.M();//uncorrected          
	m_tuple->vh_m2[ivh] = (zp4 + hp4Corr).M();//only Higgs corrected
	m_tuple->vh_m3[ivh] = (Zp4Corr + hp4).M();//only Z corrected
	m_tuple->vh_m4[ivh] = (Zp4Corr + hp4Corr).M();//both corrected


	++ivh;
      }else{
	Info( "processMMJ()" , " Trying to build too many mm-fatjet pairs, max is %i" , VHTuple::MAXVH);
      }
    }
  }
  m_tuple->nvh=ivh;
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ZHmmJTupleMaker :: processEvent(){

  if( ZHllJTupleMaker :: processEvent()  == EL::StatusCode::FAILURE) {
    //this should never happen, there are no selection applied in base classes
    cout<<"Failed ZHllJTupleMaker :: processEvent() !"<<endl;
    return EL::StatusCode::FAILURE;
  }

  if(m_debug && m_currentSyst.compare("Nominal") ==0 ){
    PrintEventParticles();
  }
  
  //remove events without a di-lepton
  if(m_tuple->nmm==0) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_NoDiLepton");


  //Select the pair with the highest sum pT
  int pair=-1;
  for (int z = 0; z < m_tuple->nmm ; z++){
    float sumpt= m_tuple->muo_pt[m_tuple->mm_leg1[z]] + m_tuple->muo_pt[m_tuple->mm_leg2[z]] ;
    if( pair == -1 || 
	sumpt > m_tuple->muo_pt[m_tuple->mm_leg1[pair]] + m_tuple->muo_pt[m_tuple->mm_leg2[pair]] )
      pair = z;
  }
  if(pair==-1){//should never happen
    cout<<"Failed finding leading pair !"<<endl;
    return EL::StatusCode::FAILURE;
  }

  const int l1=m_tuple->mm_leg1[pair];//this is already the leading
  const int l2=m_tuple->mm_leg2[pair];

  l1p4.SetPtEtaPhiM(m_tuple->muo_pt[l1],
  		    m_tuple->muo_eta[l1],
  		    m_tuple->muo_phi[l1],
  		    m_tuple->muo_m[l1]);
  l2p4.SetPtEtaPhiM(m_tuple->muo_pt[l2],
  		    m_tuple->muo_eta[l2],
  		    m_tuple->muo_phi[l2],
  		    m_tuple->muo_m[l2]);


  //Re-Apply pre-selection:
  if( m_tuple->muo_idVHl[l1] ==0 && m_tuple->muo_idVHl[l2] ==0 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonPresel");

  ///Leading Lepton
  if( fabs(m_tuple->muo_eta[l1] ) > 2.5 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonEta");
  
  if( m_tuple->muo_pt[l1] < 25000 )return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonPt");
  
  if( m_tuple->muo_isLooseIDMuon[l1] == 0 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonID");

  if( m_tuple->muo_isFromPV[l1] != 1 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonPV");


  ///////////Subleading lepton
  if( fabs(m_tuple->muo_eta[l2] ) > 2.5 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_subleadleptonEta");

  if( m_tuple->muo_pt[l2] < 25000 )return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_subleadleptonPt");

  if( m_tuple->muo_isLooseIDMuon[l2] == 0 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_subleadleptonID");

  if( m_tuple->muo_isFromPV[l2] != 1 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_subleadleptonPV");

  ///////Isolation
  if( leptonIsoOption==1 && (m_tuple->muo_isoWP1[l1] != 1 || m_tuple->muo_isoWP1[l2] != 1  ))  return EL::StatusCode::FAILURE;
  if( leptonIsoOption==2 && (m_tuple->muo_isoWP2[l1] != 1 || m_tuple->muo_isoWP2[l2] != 1  ))  return EL::StatusCode::FAILURE;
  if( leptonIsoOption==3 && (m_tuple->muo_isoWP3[l1] != 1 || m_tuple->muo_isoWP3[l2] != 1  ))  return EL::StatusCode::FAILURE;
  if( leptonIsoOption==4 && (m_tuple->muo_isoWP4[l1] != 1 || m_tuple->muo_isoWP4[l2] != 1  ))  return EL::StatusCode::FAILURE;
  if( leptonIsoOption==5 && (m_tuple->muo_isoWP5[l1] != 1 || m_tuple->muo_isoWP5[l2] != 1  ))  return EL::StatusCode::FAILURE;
  if( leptonIsoOption==6 && (m_tuple->muo_isoWP6[l1] != 1 || m_tuple->muo_isoWP6[l2] != 1  ))  return EL::StatusCode::FAILURE;
  if(leptonIsoOption == 100){
    if( (m_tuple->muo_ptiso[l1]-m_tuple->muo_trkpt[l2]*(m_tuple->mm_dR[pair]<0.2))/m_tuple->muo_pt[l1] > 0.1 ) 
      return EL::StatusCode::FAILURE;
    if( (m_tuple->muo_ptiso[l2]-m_tuple->muo_trkpt[l1]*(m_tuple->mm_dR[pair]<0.2))/m_tuple->muo_pt[l2] > 0.1 )
      return EL::StatusCode::FAILURE;
  }
  incrementCounter("eventCounter_leptonIso");
  
 
  //select the fat jet
  if( selectHiggs() == EL::StatusCode::FAILURE ) 
    return EL::StatusCode::FAILURE;


  /////////////////////////////////////////////////
  //////No more selections below here, just flags
  /////////////////////////////////////////////////
  ///Create the llJ candidates and save the indexes 
  if( processMMJ()  == EL::StatusCode::FAILURE) 
    return EL::StatusCode::FAILURE;

  m_tuple->vh=-1;
  for(int i=0;i<m_tuple->nvh;i++){
    if(m_tuple->mm_leg1[m_tuple->vh_v[i]] == l1 &&
       m_tuple->mm_leg2[m_tuple->vh_v[i]] == l2 &&
       m_tuple->vh_h[i] == higgsJet
       ) 
      m_tuple->vh=i;
  }
  if( m_tuple->vh==-1){ 
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
  

  //EOYE: HLT_mu20_iloose_L1MU15  HLT_mu50
  if(m_triggerPaths.size() > 0 ){
    if( Props::passHLT_mu20_iloose_L1MU15.get(m_eventInfo) 
	|| Props::passHLT_mu50.get(m_eventInfo) 
	) m_tuple->eve_passTrig=1;

    if( m_tuple->muo_matchHLT_mu20_iloose_L1MU15[l1] == 1
	|| m_tuple->muo_matchHLT_mu20_iloose_L1MU15[l2] == 1 
	|| m_tuple->muo_matchHLT_mu50[l1] == 1 
	|| m_tuple->muo_matchHLT_mu50[l2] == 1
	) m_tuple->eve_passTrigMatch=1;
  }


  ////Print Event Info
  if(m_debug && m_currentSyst.compare("Nominal") ==0 ){
    PrintEventFlags();
    PrintVHFlags();
    PrintVH(m_tuple->vh);
  }
  


  ////////////////BELOW JUST FILL WEIGHTS
  if( m_tuple->eve_isMC ){

    ///For Weights need the xAOD objects
    const xAOD::Muon * xAODMuon1 = NULL;
    const xAOD::Muon * xAODMuon2 = NULL;
    for( xAOD::MuonContainer::const_iterator iter = m_Muons->begin();iter!= m_Muons->end(); ++iter) {
      if((*iter)->pt() == m_tuple->muo_pt[l1])
	xAODMuon1 = *iter ;
      if((*iter)->pt() == m_tuple->muo_pt[l2])
	xAODMuon2 = *iter ;
    }
    if(xAODMuon1 == NULL || xAODMuon2 == NULL){
      cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
      return EL::StatusCode::FAILURE;
    } 



    /////////////////////////////////////
    ///Trigger Weight
    /////////////////////////////////////
    ///Note: Use scale factor for Period D (run # 276073 - 276954)
    // with proper ilumicale use   //
    double trig_w=1.;//need a double for the tool
    m_triggerTool->getTriggerDecision(m_eventInfo,trig_w,0,0,xAODMuon1,0,0,0,m_randRun,"Nominal");
    m_tuple->eve_trig_w = trig_w;

    /////////////////////////////////////
    ////vertex weight
    //////////////////////////////////
    m_tuple->muo_ttv_w = Props::TTVAEffSF.get(xAODMuon1)  * Props::TTVAEffSF.get(xAODMuon2);

    /////////////////////////////////////
    ///Id Weight
    /////////////////////////////////////
    m_tuple->muo_id_w = Props::looseEffSF.get(xAODMuon1)  * Props::looseEffSF.get(xAODMuon2);

    /////////////////////////////////////
    ///Iso Weight
    /////////////////////////////////////
    m_tuple->muo_iso_w = Props::looseTrackOnlyIsoSF.get(xAODMuon1)  * Props::looseTrackOnlyIsoSF.get(xAODMuon2);


    ///Fill weight variations only for Nominal tuple
    if( m_currentSyst.compare("Nominal") ==0 ) {


      ////////////////////////////////////////////////////////////
      ////////////Trigger Weight Variations///////////////////////
      ////////////////////////////////////////////////////////////

      ///set Electron variations to nominal
      m_tuple->eve_trig_w_EL_EFF_Trigger_TotalCorrUncertainty__1up =  m_tuple->eve_trig_w;
      m_tuple->eve_trig_w_EL_EFF_Trigger_TotalCorrUncertainty__1down =  m_tuple->eve_trig_w;

      ///muon variations retrive from corresponding branch
      m_triggerTool->getTriggerDecision(m_eventInfo,trig_w,0,0,xAODMuon1,0,0,0,m_randRun,"MUON_EFF_TrigSystUncertainty__1up");
      m_tuple->eve_trig_w_MUON_EFF_TrigSystUncertainty__1up = trig_w;
      m_triggerTool->getTriggerDecision(m_eventInfo,trig_w,0,0,xAODMuon1,0,0,0,m_randRun,"MUON_EFF_TrigSystUncertainty__1down");
      m_tuple->eve_trig_w_MUON_EFF_TrigSystUncertainty__1down = trig_w;
      m_triggerTool->getTriggerDecision(m_eventInfo,trig_w,0,0,xAODMuon1,0,0,0,m_randRun,"MUON_EFF_TrigStatUncertainty__1up");
      m_tuple->eve_trig_w_MUON_EFF_TrigStatUncertainty__1up = trig_w;
      m_triggerTool->getTriggerDecision(m_eventInfo,trig_w,0,0,xAODMuon1,0,0,0,m_randRun,"MUON_EFF_TrigStatUncertainty__1down");
      m_tuple->eve_trig_w_MUON_EFF_TrigStatUncertainty__1down = trig_w;

      /////this is how to use the CP tool directly
      // if(m_pileupreweighting->getRandomRunNumber(*m_eventInfo)!=0){    
      //   m_trig_sfmuon->setRunNumber(m_pileupreweighting->getRandomRunNumber(*m_eventInfo));
      //   CP::SystematicVariation testSys("Nominal");
      //   CP::SystematicSet shiftSet(testSys.name());
      //   if(m_trig_sfmuon->applySystematicVariation(shiftSet) != CP::SystematicCode::Ok){
      //     cout<<"m_trig_sfmuon->applySystematicVariation failed."<<endl;
      //   }else{
      //     ConstDataVector<xAOD::MuonContainer> selectedMuons(SG::VIEW_ELEMENTS);
      //     selectedMuons.push_back(xAODMuon1);
      //     m_trig_sfmuon->getTriggerScaleFactor(*selectedMuons.asDataVector(),trig_w,"HLT_mu26_imedium_OR_HLT_mu50");
      //   }
      // }



      ///////////////////////////////////////////////////////////////////////
      //for vertex, id , iso variations need to get from a corresponding particle collection
      /////////////////////////////////////////////
      //1. get the muon collection    
      //2. match the muon
      //3. get the SF
      const xAOD::MuonContainer * MuonsVar = 0;
      const xAOD::Muon * xAODMuon1Var = NULL;
      const xAOD::Muon * xAODMuon2Var = NULL;

      /////////////////////////////////////
      //// track-to-vertex scale factor
      /////////////////////////////////////
      //MUON_TTVA_STAT__1up
      if( !m_event->retrieve( MuonsVar, "Muons___MUON_TTVA_STAT__1up" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Muons___MUON_TTVA_STAT__1up");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::MuonContainer::const_iterator iter = MuonsVar->begin();iter!= MuonsVar->end(); ++iter) {
	if((*iter)->pt() == xAODMuon1->pt())
	  xAODMuon1Var = *iter ;
	if((*iter)->pt() == xAODMuon2->pt())
	  xAODMuon2Var = *iter ;
      }
      if(xAODMuon1Var == NULL || xAODMuon2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->muo_ttv_w_MUON_TTVA_STAT__1up = Props::TTVAEffSF.get(xAODMuon1Var)  * Props::TTVAEffSF.get(xAODMuon2Var);
      xAODMuon1Var = NULL; xAODMuon2Var = NULL;


      //MUON_TTVA_STAT__1down
      if( !m_event->retrieve( MuonsVar, "Muons___MUON_TTVA_STAT__1down" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Muons___MUON_TTVA_STAT__1down");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::MuonContainer::const_iterator iter = MuonsVar->begin();iter!= MuonsVar->end(); ++iter) {
	if((*iter)->pt() == xAODMuon1->pt())
	  xAODMuon1Var = *iter ;
	if((*iter)->pt() == xAODMuon2->pt())
	  xAODMuon2Var = *iter ;
      }
      if(xAODMuon1Var == NULL || xAODMuon2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->muo_ttv_w_MUON_TTVA_STAT__1down = Props::TTVAEffSF.get(xAODMuon1Var)  * Props::TTVAEffSF.get(xAODMuon2Var);
      xAODMuon1Var = NULL; xAODMuon2Var = NULL;

      //MUON_TTVA_SYS__1up
      if( !m_event->retrieve( MuonsVar, "Muons___MUON_TTVA_SYS__1up" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Muons___MUON_TTVA_SYS__1up");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::MuonContainer::const_iterator iter = MuonsVar->begin();iter!= MuonsVar->end(); ++iter) {
	if((*iter)->pt() == xAODMuon1->pt())
	  xAODMuon1Var = *iter ;
	if((*iter)->pt() == xAODMuon2->pt())
	  xAODMuon2Var = *iter ;
      }
      if(xAODMuon1Var == NULL || xAODMuon2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->muo_ttv_w_MUON_TTVA_SYS__1up = Props::TTVAEffSF.get(xAODMuon1Var)  * Props::TTVAEffSF.get(xAODMuon2Var);
      xAODMuon1Var = NULL; xAODMuon2Var = NULL;

      //MUON_TTVA_SYS__1down 
      if( !m_event->retrieve( MuonsVar, "Muons___MUON_TTVA_SYS__1down" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Muons___MUON_TTVA_SYS__1down");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::MuonContainer::const_iterator iter = MuonsVar->begin();iter!= MuonsVar->end(); ++iter) {
	if((*iter)->pt() == xAODMuon1->pt())
	  xAODMuon1Var = *iter ;
	if((*iter)->pt() == xAODMuon2->pt())
	  xAODMuon2Var = *iter ;
      }
      if(xAODMuon1Var == NULL || xAODMuon2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->muo_ttv_w_MUON_TTVA_SYS__1down = Props::TTVAEffSF.get(xAODMuon1Var)  * Props::TTVAEffSF.get(xAODMuon2Var);
      xAODMuon1Var = NULL; xAODMuon2Var = NULL;



      /////////////////////////////////////
      //ID
      /////////////////////////////////////
      ///id weight MUON_EFF_STAT__1up
      if( !m_event->retrieve( MuonsVar, "Muons___MUON_EFF_STAT__1up" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Muons___MUON_EFF_STAT__1up");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::MuonContainer::const_iterator iter = MuonsVar->begin();iter!= MuonsVar->end(); ++iter) {
	if((*iter)->pt() == xAODMuon1->pt())
	  xAODMuon1Var = *iter ;
	if((*iter)->pt() == xAODMuon2->pt())
	  xAODMuon2Var = *iter ;
      }
      if(xAODMuon1Var == NULL || xAODMuon2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->muo_id_w_MUON_EFF_STAT__1up = Props::looseEffSF.get(xAODMuon1Var)  * Props::looseEffSF.get(xAODMuon2Var);
      xAODMuon1Var = NULL; xAODMuon2Var = NULL;

      ///id weight MUON_EFF_STAT__1down
      if( !m_event->retrieve( MuonsVar, "Muons___MUON_EFF_STAT__1down" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Muons___MUON_EFF_STAT__1down");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::MuonContainer::const_iterator iter = MuonsVar->begin();iter!= MuonsVar->end(); ++iter) {
	if((*iter)->pt() == xAODMuon1->pt())
	  xAODMuon1Var = *iter ;
	if((*iter)->pt() == xAODMuon2->pt())
	  xAODMuon2Var = *iter ;
      }
      if(xAODMuon1Var == NULL || xAODMuon2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->muo_id_w_MUON_EFF_STAT__1down = Props::looseEffSF.get(xAODMuon1Var)  * Props::looseEffSF.get(xAODMuon2Var);
      xAODMuon1Var = NULL; xAODMuon2Var = NULL;

      ///id weight MUON_EFF_SYS__1up
      if( !m_event->retrieve( MuonsVar, "Muons___MUON_EFF_SYS__1up" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Muons___MUON_EFF_SYS__1up");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::MuonContainer::const_iterator iter = MuonsVar->begin();iter!= MuonsVar->end(); ++iter) {
	if((*iter)->pt() == xAODMuon1->pt())
	  xAODMuon1Var = *iter ;
	if((*iter)->pt() == xAODMuon2->pt())
	  xAODMuon2Var = *iter ;
      }
      if(xAODMuon1Var == NULL || xAODMuon2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->muo_id_w_MUON_EFF_SYS__1up = Props::looseEffSF.get(xAODMuon1Var)  * Props::looseEffSF.get(xAODMuon2Var);
      xAODMuon1Var = NULL; xAODMuon2Var = NULL;

      ///id weight MUON_EFF_SYS__1down
      if( !m_event->retrieve( MuonsVar, "Muons___MUON_EFF_SYS__1down" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Muons___MUON_EFF_SYS__1down");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::MuonContainer::const_iterator iter = MuonsVar->begin();iter!= MuonsVar->end(); ++iter) {
	if((*iter)->pt() == xAODMuon1->pt())
	  xAODMuon1Var = *iter ;
	if((*iter)->pt() == xAODMuon2->pt())
	  xAODMuon2Var = *iter ;
      }
      if(xAODMuon1Var == NULL || xAODMuon2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->muo_id_w_MUON_EFF_SYS__1down = Props::looseEffSF.get(xAODMuon1Var)  * Props::looseEffSF.get(xAODMuon2Var);
      xAODMuon1Var = NULL; xAODMuon2Var = NULL;

      /////////////////////////////////////
      ///Iso 
      /////////////////////////////////////
      ///iso weight MUON_ISO_STAT__1up
      if( !m_event->retrieve( MuonsVar, "Muons___MUON_ISO_STAT__1up" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Muons___MUON_ISO_STAT__1up");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::MuonContainer::const_iterator iter = MuonsVar->begin();iter!= MuonsVar->end(); ++iter) {
	if((*iter)->pt() == xAODMuon1->pt())
	  xAODMuon1Var = *iter ;
	if((*iter)->pt() == xAODMuon2->pt())
	  xAODMuon2Var = *iter ;
      }
      if(xAODMuon1Var == NULL || xAODMuon2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->muo_iso_w_MUON_ISO_STAT__1up = Props::looseTrackOnlyIsoSF.get(xAODMuon1Var)  * Props::looseTrackOnlyIsoSF.get(xAODMuon2Var);
      xAODMuon1Var = NULL; xAODMuon2Var = NULL;

      ///iso weight MUON_ISO_STAT__1down
      if( !m_event->retrieve( MuonsVar, "Muons___MUON_ISO_STAT__1down" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Muons___MUON_ISO_STAT__1down");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::MuonContainer::const_iterator iter = MuonsVar->begin();iter!= MuonsVar->end(); ++iter) {
	if((*iter)->pt() == xAODMuon1->pt())
	  xAODMuon1Var = *iter ;
	if((*iter)->pt() == xAODMuon2->pt())
	  xAODMuon2Var = *iter ;
      }
      if(xAODMuon1Var == NULL || xAODMuon2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->muo_iso_w_MUON_ISO_STAT__1down = Props::looseTrackOnlyIsoSF.get(xAODMuon1Var)  * Props::looseTrackOnlyIsoSF.get(xAODMuon2Var);
      xAODMuon1Var = NULL; xAODMuon2Var = NULL;


      ///iso weight MUON_ISO_SYS__1up
      if( !m_event->retrieve( MuonsVar, "Muons___MUON_ISO_SYS__1up" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Muons___MUON_ISO_SYS__1up");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::MuonContainer::const_iterator iter = MuonsVar->begin();iter!= MuonsVar->end(); ++iter) {
	if((*iter)->pt() == xAODMuon1->pt())
	  xAODMuon1Var = *iter ;
	if((*iter)->pt() == xAODMuon2->pt())
	  xAODMuon2Var = *iter ;
      }
      if(xAODMuon1Var == NULL || xAODMuon2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->muo_iso_w_MUON_ISO_SYS__1up = Props::looseTrackOnlyIsoSF.get(xAODMuon1Var)  * Props::looseTrackOnlyIsoSF.get(xAODMuon2Var);
      xAODMuon1Var = NULL; xAODMuon2Var = NULL;

      ///iso weight MUON_ISO_SYS__1down
      if( !m_event->retrieve( MuonsVar, "Muons___MUON_ISO_SYS__1down" ).isSuccess() ) {
	Error("execute", "Failed to retrieve %s contained. Exiting.", "Muons___MUON_ISO_SYS__1down");
	return EL::StatusCode::FAILURE;
      }
      for( xAOD::MuonContainer::const_iterator iter = MuonsVar->begin();iter!= MuonsVar->end(); ++iter) {
	if((*iter)->pt() == xAODMuon1->pt())
	  xAODMuon1Var = *iter ;
	if((*iter)->pt() == xAODMuon2->pt())
	  xAODMuon2Var = *iter ;
      }
      if(xAODMuon1Var == NULL || xAODMuon2Var == NULL){
	cout<<"processEvent: something is wrong xAOD muon not found for trigger weight"<<endl;
	return EL::StatusCode::FAILURE;
      } 
      m_tuple->muo_iso_w_MUON_ISO_SYS__1down = Props::looseTrackOnlyIsoSF.get(xAODMuon1Var)  * Props::looseTrackOnlyIsoSF.get(xAODMuon2Var);
      xAODMuon1Var = NULL; xAODMuon2Var = NULL;


    }
  }


  return EL::StatusCode::SUCCESS;
}



