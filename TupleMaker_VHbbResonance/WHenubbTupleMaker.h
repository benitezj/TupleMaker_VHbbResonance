#ifndef WHENUBBTUPLEMAKER_H
#define WHENUBBTUPLEMAKER_H

#include "TupleMaker_VHbbResonance/WHlnubbTupleMaker.h"
#include "TupleMaker_VHbbResonance/WHenubbTuple.h"

class WHenubbTupleMaker : public WHlnubbTupleMaker
{
public:
  WHenubbTupleMaker ();//Default Do not use.
  WHenubbTupleMaker (std::string configPath);

  virtual EL::StatusCode initialize ();

  EL::StatusCode processEMETBB();  
  

  ClassDef(WHenubbTupleMaker, 1);

private:
  
  WHenubbTuple * m_tuple =0 ; //!

  virtual EL::StatusCode processEvent();
 
};

#endif
