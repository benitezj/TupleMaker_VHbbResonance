#include "TupleMaker_VHbbResonance/WHmnubbTuple.h"

WHmnubbTuple::WHmnubbTuple():
  WHlnubbTuple()
{
}

WHmnubbTuple::~WHmnubbTuple() {;}

void WHmnubbTuple::DefineBranches(TTree * tr) {
  
  WHlnubbTuple::DefineBranches(tr);
  
}



