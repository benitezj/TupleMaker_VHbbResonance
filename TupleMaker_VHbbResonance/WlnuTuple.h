#ifndef WLNUTUPLE
#define WLNUTUPLE

#include "VTuple.h"

class WlnuTuple : public VTuple {
 public:
  WlnuTuple();
  ~WlnuTuple();


  int wlnu_lMet=0;//index of selected lepton-Met pair, should be the same as that of the lepton 

  virtual void DefineBranches(TTree * tr); 
  
 private:

};
#endif 
