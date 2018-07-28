#include "TupleMaker_VHbbResonance/WenuTupleMaker.h"

// this is needed to distribute the algorithm to the workers
ClassImp(WenuTupleMaker)


WenuTupleMaker :: WenuTupleMaker ()
{
  ///Do not use this constructor
}

WenuTupleMaker :: WenuTupleMaker (std::string configPath) :
  WlnuTupleMaker(configPath),
  m_tuple(0)
{
}

EL::StatusCode WenuTupleMaker :: initialize ()
{
  m_tuple = new WenuTuple();
  WlnuTupleMaker::m_tuple = m_tuple ;

  return WlnuTupleMaker::initialize(); 
}

EL::StatusCode WenuTupleMaker :: processEvent ()
{

  if( WlnuTupleMaker :: processEvent() ==  EL::StatusCode::FAILURE) 
    return EL::StatusCode::FAILURE;

  m_tuple->wlnu_lMet=-1;

  std::vector<int> list;
  std::vector<int> listsel;
  
  //init loop
  for(int i=0;i<m_tuple->neMet;i++){
    list.push_back(i);
  }
  if(list.size()==0) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_neMet");

  //leg1 id cut
  for (std::vector<int>::iterator it = list.begin() ; it != list.end(); ++it){
    if( m_tuple->ele_idVHl[*it] == false ) continue;
    listsel.push_back(*it);
  }
  if(listsel.size()==0) return EL::StatusCode::FAILURE;
  list=listsel;
  listsel.clear();
  incrementCounter("eventCounter_leg1_id");
 
  ///leg1 pT cut
  for(std::vector<int>::iterator it = list.begin() ; it != list.end(); ++it){
    if( m_tuple->ele_pt[*it] < 25000 ) continue;
    listsel.push_back(*it);
  }
  if(listsel.size()==0) return EL::StatusCode::FAILURE;
  list=listsel;
  listsel.clear();
  incrementCounter("eventCounter_leg1_pT");

  //leg1 eta cut
  for (std::vector<int>::iterator it = list.begin() ; it != list.end(); ++it){
    if( fabs(m_tuple->ele_eta[*it]) > 2.5 ) continue;
    listsel.push_back(*it);
  }
  if(listsel.size()==0) return EL::StatusCode::FAILURE;
  list=listsel;
  listsel.clear();
  incrementCounter("eventCounter_leg1_eta");

  //leg1 pt iso cut
  for (std::vector<int>::iterator it = list.begin() ; it != list.end(); ++it){
    if(leptonIsoOption == 1 && m_tuple->ele_ptiso[*it]/m_tuple->ele_pt[*it] > 0.1 ) continue;
    listsel.push_back(*it);
  }
  if(listsel.size()==0) return EL::StatusCode::FAILURE;
  list=listsel;
  listsel.clear();
  incrementCounter("eventCounter_leg1_ptiso");

  //leg1 et iso cut
  for (std::vector<int>::iterator it = list.begin() ; it != list.end(); ++it){
    if(leptonIsoOption == 1 && m_tuple->ele_etiso[*it]/m_tuple->ele_pt[*it] > 0.07 ) continue;
    listsel.push_back(*it);
  }
  if(listsel.size()==0) return EL::StatusCode::FAILURE;
  list=listsel;
  listsel.clear();
  incrementCounter("eventCounter_leg1_etiso");


  //remove events with 0 or >1 candidate
  if(list.size() !=1 )  return EL::StatusCode::FAILURE; 
  m_tuple->wlnu_lMet = list[0];
  incrementCounter("eventCounter_singleLepton");


  //fill 4-vectors for later
  l1p4.SetPtEtaPhiM(m_tuple->ele_pt[m_tuple->wlnu_lMet],
		    m_tuple->ele_eta[m_tuple->wlnu_lMet],
		    m_tuple->ele_phi[m_tuple->wlnu_lMet],
		    m_tuple->ele_m[m_tuple->wlnu_lMet]);
  
  
  //Fill the jets
  fillJets();

  return EL::StatusCode::SUCCESS;
}
