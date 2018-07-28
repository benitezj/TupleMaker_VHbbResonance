#ifndef WHMNUBBTUPLE
#define WHMNUBBTUPLE

#include "WHlnubbTuple.h"

class WHmnubbTuple : public WHlnubbTuple {
 public:
  WHmnubbTuple();
  ~WHmnubbTuple();
  
  virtual void DefineBranches(TTree * tr); 

 private:

};
#endif
