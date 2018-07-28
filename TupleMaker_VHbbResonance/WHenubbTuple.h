#ifndef WHENUBBTUPLE
#define WHENUBBTUPLE

#include "WHlnubbTuple.h"

class WHenubbTuple : public WHlnubbTuple {
 public:
  WHenubbTuple();
  ~WHenubbTuple();
  
  virtual void DefineBranches(TTree * tr); 

 private:

};
#endif
