#ifndef ZHEEJTUPLE
#define ZHEEJTUPLE

#include "ZHllJTuple.h"

class ZHeeJTuple : public ZHllJTuple {
 public:
  ZHeeJTuple();
  ~ZHeeJTuple();
  
  //This is the top level.
  //Any variables specific to Z-->ee  here:

  virtual void DefineBranches(TTree * tr); 

 private:

};
#endif
