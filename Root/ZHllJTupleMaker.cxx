#include "TupleMaker_VHbbResonance/ZHllJTupleMaker.h"

// this is needed to distribute the algorithm to the workers
ClassImp(ZHllJTupleMaker)


ZHllJTupleMaker :: ZHllJTupleMaker ()
{
  ///Do not use this constructor
}

ZHllJTupleMaker :: ZHllJTupleMaker (std::string configPath) :
  VHTupleMaker(configPath),
  m_tuple(0)
{

}

EL::StatusCode ZHllJTupleMaker :: initialize ()
{
  if( m_debug) std::cout<<" ZHllJTupleMaker :: initialize "<<std::endl;

  if(!m_tuple){
    //tuple must be created in a fully defined TupleMaker
    cout<<"ZHllJTupleMaker :: initialize() : Error, m_tuple should have been defined at this point"<<endl; 
    return EL::StatusCode::FAILURE;
  }
  VHTupleMaker::m_tuple = m_tuple;
  
  return VHTupleMaker::initialize();
}


EL::StatusCode ZHllJTupleMaker :: processEvent(){
  if( VHTupleMaker :: processEvent() == EL::StatusCode::FAILURE)
    return EL::StatusCode::FAILURE;
  
  //set the multi-lepton flag
  m_tuple->eve_passVHMultiLepVeto = 1;
  if( countLooseLeptons() > 2 )   m_tuple->eve_passVHMultiLepVeto = 0;
  
  return EL::StatusCode::SUCCESS;
}


bool ZHllJTupleMaker :: jetOverlap (int jetindex){
  if(deltaR(m_tuple->jet_eta[jetindex],m_tuple->jet_phi[jetindex],l1p4.Eta(),l1p4.Phi()) < 0.4) return 1;
  if(deltaR(m_tuple->jet_eta[jetindex],m_tuple->jet_phi[jetindex],l2p4.Eta(),l2p4.Phi()) < 0.4) return 1;
  return 0; 
}

bool ZHllJTupleMaker :: fatJetOverlap (int jetindex){
  if(deltaR(m_tuple->fatjet_eta[jetindex],m_tuple->fatjet_phi[jetindex],l1p4.Eta(),l1p4.Phi()) < 1.2) return 1;
  if(deltaR(m_tuple->fatjet_eta[jetindex],m_tuple->fatjet_phi[jetindex],l2p4.Eta(),l2p4.Phi()) < 1.2) return 1;
  return 0; 
}


bool ZHllJTupleMaker :: trkJetOverlap (int jetindex){
  if(deltaR(m_tuple->trkjet_eta[jetindex],m_tuple->trkjet_phi[jetindex],l1p4.Eta(),l1p4.Phi()) < 0.3) return 1;
  if(deltaR(m_tuple->trkjet_eta[jetindex],m_tuple->trkjet_phi[jetindex],l2p4.Eta(),l2p4.Phi()) < 0.3) return 1;
  return 0; 
}
