#include "TupleMaker_VHbbResonance/ZmmTupleMaker.h"


// this is needed to distribute the algorithm to the workers
ClassImp(ZmmTupleMaker)


ZmmTupleMaker :: ZmmTupleMaker ()
{
  ///Do not use this constructor
}

ZmmTupleMaker :: ZmmTupleMaker (std::string configPath) :
  ZllTupleMaker(configPath),
  m_tuple(0)
{
  

}

EL::StatusCode ZmmTupleMaker :: initialize ()
{

  //this is the top level, ntuple must be defined here
  m_tuple = new ZmmTuple();
  ZllTupleMaker::m_tuple = m_tuple;
    
  return ZllTupleMaker::initialize(); 
}

EL::StatusCode ZmmTupleMaker :: processEvent(){

  m_tuple->zll_ll=-1;

  if( ZllTupleMaker :: processEvent()  == EL::StatusCode::FAILURE) {
    //this should never happen, there are no selection applied in base classes
    cout<<"Failed ZllTupleMaker :: processEvent() !"<<endl;
    return EL::StatusCode::FAILURE;
  }

  
  //remove events without a di-lepton
  if(m_tuple->nmm==0) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_NoDiLepton");


  //Select the pair with the highest sum pT
  float maxsumpt=0.;
  for (int z = 0; z < m_tuple->nmm ; z++){
    float sumpt= m_tuple->muo_pt[m_tuple->mm_leg1[z]] + m_tuple->muo_pt[m_tuple->mm_leg2[z]] ;
    if( maxsumpt == 0. ||  sumpt > maxsumpt ){
      m_tuple->zll_ll = z;
      maxsumpt = sumpt;
    }
  }
  if(m_tuple->zll_ll ==-1 ){//should never happen
    cout<<"Failed finding leading pair !"<<endl;
    return EL::StatusCode::FAILURE;
  }

  m_tuple->zll_l1 = m_tuple->mm_leg1[m_tuple->zll_ll];//this is already the leading
  m_tuple->zll_l2 = m_tuple->mm_leg2[m_tuple->zll_ll];

  l1p4.SetPtEtaPhiM(m_tuple->muo_pt[ m_tuple->zll_l1],
                    m_tuple->muo_eta[ m_tuple->zll_l1],
                    m_tuple->muo_phi[ m_tuple->zll_l1],
                    m_tuple->muo_m[ m_tuple->zll_l1]);

  l2p4.SetPtEtaPhiM(m_tuple->muo_pt[ m_tuple->zll_l2],
                    m_tuple->muo_eta[ m_tuple->zll_l2],
                    m_tuple->muo_phi[ m_tuple->zll_l2],
                    m_tuple->muo_m[ m_tuple->zll_l2]);


  ///Leading Lepton
  if( m_tuple->muo_idVHl[m_tuple->zll_l1] == 0 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonPresel");

  if( fabs(m_tuple->muo_eta[m_tuple->zll_l1] ) > 2.5 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonEta");
  
  if( m_tuple->muo_pt[m_tuple->zll_l1] < 25000 )return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonPt");
  
  if( m_tuple->muo_isMediumIDMuon[m_tuple->zll_l1] == 0 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_leadleptonID");


  ///////////Subleading lepton
  //if( m_tuple->muo_idVHl[m_tuple->zll_l2] == 0 ) return EL::StatusCode::FAILURE;
  //incrementCounter("eventCounter_subleadleptonPresel");
  //this cut is not needed because only leading lepton is required at preselection

  if( fabs(m_tuple->muo_eta[m_tuple->zll_l2] ) > 2.5 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_subleadleptonEta");

  if( m_tuple->muo_pt[m_tuple->zll_l2] < 10000 )return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_subleadleptonPt");

  if( m_tuple->muo_isLooseIDMuon[m_tuple->zll_l2] == 0 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_subleadleptonID");



  ///////Isolation
  if( leptonIsoOption==1 && (m_tuple->muo_isoWP1[m_tuple->zll_l1] != 1 || m_tuple->muo_isoWP1[m_tuple->zll_l2] != 1  ))
    return EL::StatusCode::FAILURE;
  if( leptonIsoOption==2 && (m_tuple->muo_isoWP2[m_tuple->zll_l1] != 1 || m_tuple->muo_isoWP2[m_tuple->zll_l2] != 1  )) 
    return EL::StatusCode::FAILURE;
  if( leptonIsoOption==3 && (m_tuple->muo_isoWP3[m_tuple->zll_l1] != 1 || m_tuple->muo_isoWP3[m_tuple->zll_l2] != 1  )) 
    return EL::StatusCode::FAILURE;
  if( leptonIsoOption==4 && (m_tuple->muo_isoWP4[m_tuple->zll_l1] != 1 || m_tuple->muo_isoWP4[m_tuple->zll_l2] != 1  )) 
    return EL::StatusCode::FAILURE;
  if( leptonIsoOption==5 && (m_tuple->muo_isoWP5[m_tuple->zll_l1] != 1 || m_tuple->muo_isoWP5[m_tuple->zll_l2] != 1  )) 
    return EL::StatusCode::FAILURE;
  if( leptonIsoOption==6 && (m_tuple->muo_isoWP6[m_tuple->zll_l1] != 1 || m_tuple->muo_isoWP6[m_tuple->zll_l2] != 1  )) 
    return EL::StatusCode::FAILURE;
  if(leptonIsoOption == 100){
    if( (m_tuple->muo_ptiso[m_tuple->zll_l1]-m_tuple->muo_trkpt[m_tuple->zll_l2]*(m_tuple->mm_dR[m_tuple->zll_ll]<0.2))/m_tuple->muo_pt[m_tuple->zll_l1] > 0.1 ) 
      return EL::StatusCode::FAILURE;
    if( (m_tuple->muo_ptiso[m_tuple->zll_l2]-m_tuple->muo_trkpt[m_tuple->zll_l1]*(m_tuple->mm_dR[m_tuple->zll_ll]<0.2))/m_tuple->muo_pt[m_tuple->zll_l2] > 0.1 )
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
    if( Props::passHLT_mu50.get(m_eventInfo) )
      m_tuple->eve_passTrig=1;
    
    if( m_tuple->muo_matchHLT_mu50[m_tuple->zll_l1] == 1  
	||  m_tuple->muo_matchHLT_mu50[m_tuple->zll_l2] == 1  )
      m_tuple->eve_passTrigMatch=1;
  }
  /////////////////////////////////////


  
  //Fill the jets
  fillJets();

  
  return EL::StatusCode::SUCCESS;
}
