#ifndef WHLNUBBTUPLEMAKER_H
#define WHLNUBBTUPLEMAKER_H

#include "TupleMaker_VHbbResonance/VHTupleMaker.h"
#include "TupleMaker_VHbbResonance/WHlnubbTuple.h"

class WHlnubbTupleMaker : public VHTupleMaker
{
public:
  WHlnubbTupleMaker ();//Default Do not use.
  WHlnubbTupleMaker (std::string configPath);

  virtual EL::StatusCode initialize ();
  
  // this is needed to distribute the algorithm to the workers
  ClassDef(WHlnubbTupleMaker, 1);

protected:
 
  WHlnubbTuple * m_tuple = 0 ; //!
  TLorentzVector lp4;
  

  virtual EL::StatusCode processEvent();

  bool jetOverlap(int jetindex); 
  bool fatJetOverlap(int jetindex); 
  bool trkJetOverlap(int jetindex); 

private:

};

#endif
