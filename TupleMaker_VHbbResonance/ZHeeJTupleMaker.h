#ifndef ZHEEJTUPLEMAKER_H
#define ZHEEJTUPLEMAKER_H

#include "TupleMaker_VHbbResonance/ZHllJTupleMaker.h"
#include "TupleMaker_VHbbResonance/ZHeeJTuple.h"

class ZHeeJTupleMaker : public ZHllJTupleMaker
{
public:
  ZHeeJTupleMaker ();//Default Do not use.
  ZHeeJTupleMaker (std::string configPath);

  virtual EL::StatusCode initialize ();

  EL::StatusCode processEEJ();  
  

  void PrintVH(int i){
    std::cout<<"PrintVH: pT"<<m_tuple->vh_pt[i]<<" eta="<<m_tuple->vh_eta[i]<<" phi="<<m_tuple->vh_phi[i]<<" M="<<m_tuple->vh_m[i]<<" dR="<<m_tuple->vh_dR[i]<<std::endl;
    PrintEE(m_tuple->vh_v[i]);
    PrintFatJet(m_tuple->vh_h[i]);
  }

  ClassDef(ZHeeJTupleMaker, 1);

private:
  
  ZHeeJTuple * m_tuple =0 ; //!

  virtual EL::StatusCode processEvent();
 
};

#endif
