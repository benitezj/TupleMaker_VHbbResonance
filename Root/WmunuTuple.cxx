#include "TupleMaker_VHbbResonance/WmunuTuple.h"

WmunuTuple::WmunuTuple():
  WlnuTuple()
{
}

WmunuTuple::~WmunuTuple() {;}

void WmunuTuple::DefineBranches(TTree * tr) {
  WlnuTuple::DefineBranches(tr);
}


// template <class T> void WmunuTuple::SetBranchAddressSafe(TTree * tr, const char * label, T * var) {
//   if( tr->GetBranch(label) ) {
//     tr->SetBranchAddress(label,var);
//   }
// }
// void WmunuTuple::SetBranchAddresses(TTree * tr) {
//   SetBranchAddressSafe(tr,"eve_num",&eve_num);
// }


