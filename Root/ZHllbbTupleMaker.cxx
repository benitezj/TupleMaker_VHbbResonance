#include "TupleMaker_VHbbResonance/ZHllbbTupleMaker.h"

// this is needed to distribute the algorithm to the workers
ClassImp(ZHllbbTupleMaker)


ZHllbbTupleMaker :: ZHllbbTupleMaker ()
{
  ///Do not use this constructor
}

ZHllbbTupleMaker :: ZHllbbTupleMaker (std::string configPath) :
  VHTupleMaker(configPath),
  m_tuple(0)
{

}

EL::StatusCode ZHllbbTupleMaker :: initialize ()
{
  if(!m_tuple){
    //tuple must be created in a fully defined TupleMaker
    cout<<"ZHllbbTupleMaker :: initialize() : Error, m_tuple should have been defined at this point"<<endl; 
    return EL::StatusCode::FAILURE;
  }
  VHTupleMaker::m_tuple = m_tuple;
  
  return VHTupleMaker::initialize();
}


EL::StatusCode ZHllbbTupleMaker :: processEvent(){
  return VHTupleMaker :: processEvent();
}


bool ZHllbbTupleMaker :: jetOverlap (int jetindex){
  if(deltaR(m_tuple->jet_eta[jetindex],m_tuple->jet_phi[jetindex],l1p4.Eta(),l1p4.Phi()) < 0.4) return 1;
  if(deltaR(m_tuple->jet_eta[jetindex],m_tuple->jet_phi[jetindex],l2p4.Eta(),l2p4.Phi()) < 0.4) return 1;
  return 0; 
}

bool ZHllbbTupleMaker :: fatJetOverlap (int jetindex){
  if(deltaR(m_tuple->fatjet_eta[jetindex],m_tuple->fatjet_phi[jetindex],l1p4.Eta(),l1p4.Phi()) < 1.0) return 1;
  if(deltaR(m_tuple->fatjet_eta[jetindex],m_tuple->fatjet_phi[jetindex],l2p4.Eta(),l2p4.Phi()) < 1.0) return 1;
  return 0; 
}

bool ZHllbbTupleMaker :: trkJetOverlap (int jetindex){
  if(deltaR(m_tuple->trkjet_eta[jetindex],m_tuple->trkjet_phi[jetindex],l1p4.Eta(),l1p4.Phi()) < 0.3) return 1;
  if(deltaR(m_tuple->trkjet_eta[jetindex],m_tuple->trkjet_phi[jetindex],l2p4.Eta(),l2p4.Phi()) < 0.3) return 1;
  return 0; 
}

