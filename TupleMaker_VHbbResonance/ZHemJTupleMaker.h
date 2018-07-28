#ifndef ZHEMJTUPLEMAKER_H
#define ZHEMJTUPLEMAKER_H

#include "TupleMaker_VHbbResonance/ZHllJTupleMaker.h"
#include "TupleMaker_VHbbResonance/ZHemJTuple.h"

class ZHemJTupleMaker : public ZHllJTupleMaker
{
public:
  ZHemJTupleMaker ();//Default Do not use.
  ZHemJTupleMaker (std::string configPath);

  virtual EL::StatusCode initialize ();

  EL::StatusCode processEMJ();  
  

  ClassDef(ZHemJTupleMaker, 1);

private:
  
  ZHemJTuple * m_tuple =0 ; //!

  virtual EL::StatusCode processEvent();
 
};

#endif
