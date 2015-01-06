#include "ShowerDeconstruction/SDAlgorithm/interface/BModel.h"

#include <iostream>
#include <vector>
#include <fastjet/PseudoJet.hh>
#include "ShowerDeconstruction/SDAlgorithm/interface/Helper.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/Parameters.h"

using namespace Deconstruction;
using namespace std;

Deconstruction::BModel::BModel(Parameters &param, Flavour::id flavourb, Flavour::id flavourt)
  : Deconstruction::Model::Model(param), m_flavourb(flavourb), m_flavourt(flavourt) {
}

Deconstruction::BModel::~BModel() {
}

//double Deconstruction::BModel::weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sum) {
double Deconstruction::BModel::weight(const StoredJet &jets, fastjet::PseudoJet &sum) {
  // top is color connected partner von der ersten emission vom b
  // evtl sollte top bzw. bottom flavor mitgegeben werden - sind ja korrelliert.

  StoredJet topset(jets.store(), m_topset);

#ifdef DEBUGCODE
  printoutput(jets,"input in btopmodel", DEBUG);
  LOG(DEBUG) << "input flavor: " << m_flavourt << endl;
  LOG(DEBUG) << "top: " << endl;
  printoutput(topset, DEBUG);
#endif

  // check for theta cut if there is no grandmother:

  double bglusplitweight = 0;

  StoredJet empty(jets.store());

  // weight is the same for qquark and antiqquark...
  if ( (m_flavourt == Flavour::t) && (m_flavourb == Flavour::b) ) {
    bglusplitweight = make_splitting(jets,empty,topset,topset,Flavour::b,Flavour::t,Flavour::noflav,Flavour::t,Shower::b);  // b radiates right (top right col partner)
  } else if ( (m_flavourt == Flavour::t) && (m_flavourb == Flavour::bbar) ) {
    bglusplitweight = make_splitting(jets,topset,empty,topset,Flavour::b,Flavour::t,Flavour::t,Flavour::noflav,Shower::b); // antib radiates left (top left col.partner)
  } else
    LOG(ERROR) << "ERROR in btopmodel: either tflavor or bflavor wrong...";

#ifdef DEBUGCODE
  LOG(DEBUG) << "bglusplitweight: " << bglusplitweight << endl;
#endif

  return bglusplitweight;
}

//void Deconstruction::BModel::setTopset(const std::vector<fastjet::PseudoJet> &topset) {
void Deconstruction::BModel::setTopset(const StoredJet &topset) {
  m_topset = topset.getList();
}

double Deconstruction::BModel::hamiltonian(double pTsum) {
  return 1.0;
}

