#ifndef ZHEMJTUPLE
#define ZHEMJTUPLE

#include "ZHllJTuple.h"

class ZHemJTuple : public ZHllJTuple {
 public:
  ZHemJTuple();
  ~ZHemJTuple();
  
  //This is the top level.
  //Any variables specific to e-mu here:

  virtual void DefineBranches(TTree * tr); 

 private:

};
#endif
