#ifndef ZEETUPLEMAKER_H
#define ZEETUPLEMAKER_H

#include "TupleMaker_VHbbResonance/ZllTupleMaker.h"
#include "TupleMaker_VHbbResonance/ZeeTuple.h"

class ZeeTupleMaker : public ZllTupleMaker
{
public:
  ZeeTupleMaker ();//Default Do not use.
  ZeeTupleMaker (std::string configPath);

  virtual EL::StatusCode initialize ();

  ClassDef(ZeeTupleMaker, 1);

private:
  ZeeTuple * m_tuple = 0; //!

  virtual EL::StatusCode processEvent();

};

#endif
