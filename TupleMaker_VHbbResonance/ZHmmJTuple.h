#ifndef ZHMMJTUPLE
#define ZHMMJTUPLE

#include "ZHllJTuple.h"

class ZHmmJTuple : public ZHllJTuple {
 public:
  ZHmmJTuple();
  ~ZHmmJTuple();
  
  //This is the top level.
  //Any variables specific to Z-->mu mu here:

  virtual void DefineBranches(TTree * tr); 

 private:

};
#endif
