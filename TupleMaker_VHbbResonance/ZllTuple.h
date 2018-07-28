#ifndef ZLLTUPLE
#define ZLLTUPLE

#include "VTuple.h"

class ZllTuple : public VTuple {
 public:
  ZllTuple();
  ~ZllTuple();

  int zll_ll=0;//index of ll pair

  int zll_l1=0;//leading leg
  int zll_l2=0;//subleading leg

  virtual void DefineBranches(TTree * tr); 

 private:
  
};
#endif 
