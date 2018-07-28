#include "TupleMaker_VHbbResonance/VTuple.h"

VTuple::VTuple():
  RecoTuple()
{
}

VTuple::~VTuple() {;}

void VTuple::DefineBranches(TTree * tr) {
  RecoTuple::DefineBranches(tr);

  tr->Branch("v_njet",&v_njet,"v_njet/I");
  tr->Branch("v_jet",&v_jet,"v_jet/I");

  tr->Branch("v_njetfw",&v_njetfw,"v_njetfw/I");
  tr->Branch("v_jetfw",&v_jetfw,"v_jetfw/I");

  tr->Branch("v_nfatjet",&v_nfatjet,"v_nfatjet/I");
  tr->Branch("v_fatjet",&v_fatjet,"v_fatjet/I");

  tr->Branch("v_ntrkjet",&v_ntrkjet,"v_ntrkjet/I");
  tr->Branch("v_trkjet",&v_trkjet,"v_trkjet/I");

}
