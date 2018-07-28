#ifndef WLNUTUPLEMAKER_H
#define WLNUTUPLEMAKER_H

#include "TupleMaker_VHbbResonance/RecoTupleMaker.h"
#include "TupleMaker_VHbbResonance/WlnuTuple.h"

class WlnuTupleMaker : public RecoTupleMaker
{
public:
  WlnuTupleMaker ();//Default Do not use.
  WlnuTupleMaker (std::string configPath);

  virtual EL::StatusCode initialize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(WlnuTupleMaker, 1);

protected:

  WlnuTuple * m_tuple = 0; //!

  int leptonIsoOption;

  TLorentzVector l1p4;

  virtual EL::StatusCode processEvent();

  EL::StatusCode fillJets ();

private:



};

#endif
