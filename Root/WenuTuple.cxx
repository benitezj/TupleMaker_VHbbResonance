#include "TupleMaker_VHbbResonance/WenuTuple.h"

WenuTuple::WenuTuple():
  WlnuTuple()
{
}

WenuTuple::~WenuTuple() {;}

void WenuTuple::DefineBranches(TTree * tr) {
  WlnuTuple::DefineBranches(tr);
  
}


// template <class T> void WenuTuple::SetBranchAddressSafe(TTree * tr, const char * label, T * var) {
//   if( tr->GetBranch(label) ) {
//     tr->SetBranchAddress(label,var);
//   }
// }
// void WenuTuple::SetBranchAddresses(TTree * tr) {
//   SetBranchAddressSafe(tr,"eve_num",&eve_num);
// }


