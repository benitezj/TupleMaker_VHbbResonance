#include "TupleMaker_VHbbResonance/WlnuTuple.h"

WlnuTuple::WlnuTuple():
  VTuple()
{
}

WlnuTuple::~WlnuTuple() {;}

void WlnuTuple::DefineBranches(TTree * tr) {
  VTuple::DefineBranches(tr);

  tr->Branch("wlnu_lMet",&wlnu_lMet,"wlnu_lMet/I");

}

