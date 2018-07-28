#include "TupleMaker_VHbbResonance/ZeeTupleMaker.h"

// this is needed to distribute the algorithm to the workers
ClassImp(ZeeTupleMaker)


ZeeTupleMaker :: ZeeTupleMaker ()
{
  ///Do not use this constructor
}

ZeeTupleMaker :: ZeeTupleMaker (std::string configPath) :
  ZllTupleMaker(configPath),
  m_tuple(0)
{
  
}

EL::StatusCode ZeeTupleMaker :: initialize ()
{
  
  m_tuple = new ZeeTuple();
  ZllTupleMaker::m_tuple = m_tuple;

  return ZllTupleMaker::initialize(); 
}


EL::StatusCode ZeeTupleMaker :: processEvent ()
{
  m_tuple->zll_ll=-1;

  if( ZllTupleMaker :: processEvent()  == EL::StatusCode::FAILURE) {
    //this should never happen, there are no selection applied in base classes
    cout<<"Failed ZllTupleMaker :: processEvent() !"<<endl;
    return EL::StatusCode::FAILURE;
  }

  
  //remove events without a di-lepton
  if(m_tuple->nee==0) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_NoDiLepton");


  //Select the pair with the highest sum pT
  float maxsumpt=0.;
  for (int z = 0; z < m_tuple->nee ; z++){
    float sumpt= m_tuple->ele_pt[m_tuple->ee_leg1[z]] + m_tuple->ele_pt[m_tuple->ee_leg2[z]] ;
    if( maxsumpt == 0. ||  sumpt > maxsumpt ){
      m_tuple->zll_ll = z;
      maxsumpt = sumpt;
    }
  }
  if(m_tuple->zll_ll ==-1 ){//should never happen
    cout<<"Failed finding leading pair !"<<endl;
    return EL::StatusCode::FAILURE;
  }

  m_tuple->zll_l1 = m_tuple->ee_leg1[m_tuple->zll_ll];//this is already the leading
  m_tuple->zll_l2 = m_tuple->ee_leg2[m_tuple->zll_ll];


  l1p4.SetPtEtaPhiM(m_tuple->ele_pt[ m_tuple->zll_l1],
                    m_tuple->ele_eta[ m_tuple->zll_l1],
                    m_tuple->ele_phi[ m_tuple->zll_l1],
                    m_tuple->ele_m[ m_tuple->zll_l1]);

  l2p4.SetPtEtaPhiM(m_tuple->ele_pt[ m_tuple->zll_l2],
                    m_tuple->ele_eta[ m_tuple->zll_l2],
                    m_tuple->ele_phi[ m_tuple->zll_l2],
                    m_tuple->ele_m[ m_tuple->zll_l2]);


  ///Leading Lepton
  if( m_tuple->ele_idVHl[m_tuple->zll_l1] == 0 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonPresel");

  if( fabs(m_tuple->ele_eta[m_tuple->zll_l1] ) > 2.47 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonEta");

  if( fabs(m_tuple->ele_eta[m_tuple->zll_l1] ) > 1.37 && fabs(m_tuple->ele_eta[m_tuple->zll_l1] ) < 1.52) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonCrackVeto");
  
  if( m_tuple->ele_pt[m_tuple->zll_l1] < 25000 )return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonPt");
  
  if( m_tuple->ele_isMediumIDElectron[m_tuple->zll_l1] == 0 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonID");


  ///////////Subleading lepton
  //if( m_tuple->ele_idVHl[m_tuple->zll_l2] == 0 ) return EL::StatusCode::FAILURE;
  //incrementCounter("eventCounter_subleadleptonPresel");
  //this cut is not needed because only leading lepton is required at preselection

  if( fabs(m_tuple->ele_eta[m_tuple->zll_l2] ) > 2.5 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_subleadleptonEta");

  if( fabs(m_tuple->ele_eta[m_tuple->zll_l2] ) > 1.37 && fabs(m_tuple->ele_eta[m_tuple->zll_l2] ) < 1.52) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonCrackVeto");

  if( m_tuple->ele_pt[m_tuple->zll_l2] < 10000 )return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_subleadleptonPt");

  if( m_tuple->ele_isLooseIDElectron[m_tuple->zll_l2] == 0 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_subleadleptonID");



  ///////Isolation
  if( leptonIsoOption==1 && (m_tuple->ele_isoWP1[m_tuple->zll_l1] != 1 || m_tuple->ele_isoWP1[m_tuple->zll_l2] != 1  ))
    return EL::StatusCode::FAILURE;
  if( leptonIsoOption==2 && (m_tuple->ele_isoWP2[m_tuple->zll_l1] != 1 || m_tuple->ele_isoWP2[m_tuple->zll_l2] != 1  )) 
    return EL::StatusCode::FAILURE;
  if( leptonIsoOption==3 && (m_tuple->ele_isoWP3[m_tuple->zll_l1] != 1 || m_tuple->ele_isoWP3[m_tuple->zll_l2] != 1  )) 
    return EL::StatusCode::FAILURE;
  if( leptonIsoOption==4 && (m_tuple->ele_isoWP4[m_tuple->zll_l1] != 1 || m_tuple->ele_isoWP4[m_tuple->zll_l2] != 1  )) 
    return EL::StatusCode::FAILURE;
  if( leptonIsoOption==5 && (m_tuple->ele_isoWP5[m_tuple->zll_l1] != 1 || m_tuple->ele_isoWP5[m_tuple->zll_l2] != 1  )) 
    return EL::StatusCode::FAILURE;
  if( leptonIsoOption==6 && (m_tuple->ele_isoWP6[m_tuple->zll_l1] != 1 || m_tuple->ele_isoWP6[m_tuple->zll_l2] != 1  )) 
    return EL::StatusCode::FAILURE;
  if(leptonIsoOption == 100){
    if( (m_tuple->ele_ptiso[m_tuple->zll_l1]-m_tuple->ele_trkpt[m_tuple->zll_l2]*(m_tuple->mm_dR[m_tuple->zll_ll]<0.2))/m_tuple->ele_pt[m_tuple->zll_l1] > 0.1 ) 
      return EL::StatusCode::FAILURE;
    if( (m_tuple->ele_ptiso[m_tuple->zll_l2]-m_tuple->ele_trkpt[m_tuple->zll_l1]*(m_tuple->mm_dR[m_tuple->zll_ll]<0.2))/m_tuple->ele_pt[m_tuple->zll_l2] > 0.1 )
      return EL::StatusCode::FAILURE;
  }
  incrementCounter("eventCounter_leptonIso");


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
    if( Props::passHLT_e60_lhmedium.get(m_eventInfo) )
      m_tuple->eve_passTrig=1;
    
    if( m_tuple->ele_matchHLT_e60_lhmedium[m_tuple->zll_l1] == 1  
        ||  m_tuple->ele_matchHLT_e60_lhmedium[m_tuple->zll_l2] == 1  )
      m_tuple->eve_passTrigMatch=1;
  }
  /////////////////////////////////////



  ///
  fillJets();
  
  return EL::StatusCode::SUCCESS;
}

