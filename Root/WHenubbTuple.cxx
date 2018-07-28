#include "TupleMaker_VHbbResonance/WHenubbTuple.h"

WHenubbTuple::WHenubbTuple():
  WHlnubbTuple()
{
}

WHenubbTuple::~WHenubbTuple() {;}

void WHenubbTuple::DefineBranches(TTree * tr) {
  
  WHlnubbTuple::DefineBranches(tr);
  
}



