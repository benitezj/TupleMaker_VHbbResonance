#ifndef ZHEEBBTUPLEMAKER_H
#define ZHEEBBTUPLEMAKER_H

#include "TupleMaker_VHbbResonance/ZHllbbTupleMaker.h"
#include "TupleMaker_VHbbResonance/ZHeebbTuple.h"

class ZHeebbTupleMaker : public ZHllbbTupleMaker
{
public:
  ZHeebbTupleMaker ();//Default Do not use.
  ZHeebbTupleMaker (std::string configPath);

  virtual EL::StatusCode initialize ();

  EL::StatusCode processEEBB();  
  

  ClassDef(ZHeebbTupleMaker, 1);

private:
  
  ZHeebbTuple * m_tuple =0 ; //!

  virtual EL::StatusCode processEvent();
 
};

#endif
