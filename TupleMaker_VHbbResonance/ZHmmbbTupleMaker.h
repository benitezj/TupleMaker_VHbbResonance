#ifndef ZHMMBBTUPLEMAKER_H
#define ZHMMBBTUPLEMAKER_H

#include "TupleMaker_VHbbResonance/ZHllbbTupleMaker.h"
#include "TupleMaker_VHbbResonance/ZHmmbbTuple.h"

class ZHmmbbTupleMaker : public ZHllbbTupleMaker
{
public:
  ZHmmbbTupleMaker ();//Default Do not use.
  ZHmmbbTupleMaker (std::string configPath);

  virtual EL::StatusCode initialize ();

  EL::StatusCode processMMBB();  
  

  ClassDef(ZHmmbbTupleMaker, 1);

private:
  
  ZHmmbbTuple * m_tuple =0 ; //!

  virtual EL::StatusCode processEvent();
 
};

#endif
