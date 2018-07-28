#ifndef ZHMMJTUPLEMAKER_H
#define ZHMMJTUPLEMAKER_H

#include "TupleMaker_VHbbResonance/ZHllJTupleMaker.h"
#include "TupleMaker_VHbbResonance/ZHmmJTuple.h"

class ZHmmJTupleMaker : public ZHllJTupleMaker
{
public:
  ZHmmJTupleMaker ();//Default Do not use.
  ZHmmJTupleMaker (std::string configPath);

  virtual EL::StatusCode initialize ();

  EL::StatusCode processMMJ();  


  void PrintVH(int i){
    std::cout<<"PrintVH: pT"<<m_tuple->vh_pt[i]<<" eta="<<m_tuple->vh_eta[i]<<" phi="<<m_tuple->vh_phi[i]<<" M="<<m_tuple->vh_m[i]<<" dR="<<m_tuple->vh_dR[i]<<std::endl;
    PrintMM(m_tuple->vh_v[i]);
    PrintFatJet(m_tuple->vh_h[i]);
  }

  ClassDef(ZHmmJTupleMaker, 1);

private:
  
  ZHmmJTuple * m_tuple =0 ; //!

  virtual EL::StatusCode processEvent();
 
};

#endif
