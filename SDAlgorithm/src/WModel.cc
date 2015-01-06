#include "ShowerDeconstruction/SDAlgorithm/interface/WModel.h"

#include <vector>
#include <fastjet/PseudoJet.hh>
#include "ShowerDeconstruction/SDAlgorithm/interface/Helper.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/AnalysisParameters.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/Parameters.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/JetInfo.h"
#include <iostream>
#include <algorithm>

using namespace Deconstruction;
using namespace std;
using namespace fastjet;

Deconstruction::WModel::WModel(Parameters &param, Flavour::id flavour, bool includeAngleFactor)
  : Deconstruction::Model::Model(param), m_flavour(flavour), m_includeAngleFactor(includeAngleFactor) {
}

Deconstruction::WModel::~WModel() {
}

//void Deconstruction::WModel::setTopset(const std::vector<fastjet::PseudoJet> &topset) {
void Deconstruction::WModel::setTopset(const StoredJet &topset) {
  m_topset = topset.getList();
}

//double Deconstruction::WModel::weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sumJets) {
double Deconstruction::WModel::weight(const StoredJet &jets, fastjet::PseudoJet &sumJets) {

  double signalweight = 0;

  StoredJet topset(jets.store(), m_topset);

  PseudoJet Wjet = sumJets;
  PseudoJet tjet = topset.sum();


  // possible combinations for the set in [first, last] are given by permutations of the binary
  // number with p->mom().size() bits with any number of bits activated
  unsigned int iElements = (unsigned int) jets.size();
  // iElements cannot be greater than the maximum,
  // but this has already been tested before, when generating the powerset for ISR/SB
  // there is no point in testing it again, since these sets must be smaller or equal in size
  // than the subsets of the powerset of all microjets
  unsigned long long nSets = (unsigned long long) (0x1 << iElements); // 2^iElements

  StoredJet empty(jets.store());

  // the bits position in index i give the positions of the vector to be used
  // we are skipping the empty set combination here
  // that's why i starts at one and ends at nSets-1
  for (unsigned long long i = 1; i < nSets; ++i) {

    //std::vector<fastjet::PseudoJet> ubranch;
    //std::vector<fastjet::PseudoJet> dbranch;
    StoredJet ubranch(jets.store());
    StoredJet dbranch(jets.store());

    // check the bits set in i and include the elements corresponding to
    // that bit position in result
    for (unsigned int k = 0; k < iElements; ++k) {
      if ( (i & (0x1 << k)) != 0) {
        ubranch.push_back(jets.getIdx(k));
      } else {
        dbranch.push_back(jets.getIdx(k));
      }
    }

    // quark goes left, antiquark goes right
    // if top hypothesis -> W+ > dbar u
    // grandmother is set to be empty because the scale where new shower starts is set by hand

    double valueleft = 0;
    double valueright = 0;

    if (m_flavour == Flavour::Wp) {
      valueleft = make_splitting(ubranch, empty, dbranch, jets, Flavour::q, Flavour::noflav, Flavour::noflav, Flavour::qbar, Shower::W);
      valueright = make_splitting(dbranch, ubranch, empty, jets, Flavour::qbar, Flavour::noflav, Flavour::q, Flavour::noflav,Shower::W);
    } else if (m_flavour == Flavour::Wm) {
      valueleft = make_splitting(dbranch, empty, ubranch, jets, Flavour::q, Flavour::noflav, Flavour::noflav, Flavour::qbar, Shower::W);
      valueright = make_splitting(ubranch, dbranch, empty, jets, Flavour::qbar, Flavour::noflav,Flavour::q, Flavour::noflav, Shower::W);
    } else LOG(ERROR) << "ERROR: wrong flavor in Wmodel: " << m_flavour << endl;

    PseudoJet djet = dbranch.sum();

    double anglefacW = 1;
    if (m_includeAngleFactor) {
      anglefacW = 12.0*scalprod(djet,tjet)*square_p1minusp2(tjet,djet);
      anglefacW /= (square(param()[p_topmass]) - square(param()[p_wmass])) * (square(param()[p_topmass]) + 2*square(param()[p_wmass]));
    }

    double Hwsplit = 8*square(M_PI)*param()[p_wmass]*wwidth/std::atan(std::fabs(square(Wjet.m())-square(param()[p_wmass]))/(param()[p_wmass]*wwidth));
    Hwsplit /= square(square(Wjet.m()) - square(param()[p_wmass])) + square(param()[p_wmass])*square(wwidth);

    ///////// sudakov
    double expSW = 0;
    double Wsud = 0;

    if(std::fabs(square(Wjet.m()) - square(param()[p_wmass])) < param()[p_wmass]*delta_wmass) {
      Wsud=(std::log(std::atan(delta_wmass/wwidth)) - std::log(std::atan(std::fabs(square(Wjet.m())-square(param()[p_wmass]))/(param()[p_wmass]*wwidth))));
      expSW = exp(-Wsud);
    } else
      expSW = 1;

    double decayfactorW = Hwsplit *anglefacW* expSW;

    signalweight += decayfactorW*valueleft*valueright;

#ifdef DEBUGCODE
    LOG(DEBUG) << "top input in Wmodel" << endl;
    printoutput(topset, DEBUG);
    LOG(DEBUG) << "input in Wmodel" << endl;
    printoutput(jets, DEBUG);
    LOG(DEBUG) << "ubranch" << endl;
    printoutput(ubranch, DEBUG);
    LOG(DEBUG) << "dbranch" << endl;
    printoutput(dbranch, DEBUG);
    LOG(DEBUG) << "input flavor: " << m_flavour << endl;
    LOG(DEBUG) << "anglefacW: " << anglefacW << endl;
    LOG(DEBUG) << "Hwsplit: " << Hwsplit << endl;
    LOG(DEBUG) << "Wsud: " << Wsud << endl;
    LOG(DEBUG) << "decayfactorW: " << decayfactorW << endl;
    LOG(DEBUG) << "value left: " << valueleft << endl;

    LOG(DEBUG) << "value right: " << valueright << endl;
    LOG(DEBUG) << "W shower: decayfactorW*valueleft*valueright :  " << decayfactorW*valueleft*valueright << endl;
#endif

  }

#ifdef DEBUGCODE
  LOG(DEBUG) << "W showre total signalweight: " << signalweight << endl;
#endif

  return param()[p_br_wqq]*signalweight;
}

double Deconstruction::WModel::hamiltonian(double pTsum) {
  return 1.0;
}

