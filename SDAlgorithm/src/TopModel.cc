#include "ShowerDeconstruction/SDAlgorithm/interface/TopModel.h"

#include <vector>
#include <fastjet/PseudoJet.hh>
#include "ShowerDeconstruction/SDAlgorithm/interface/Helper.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/JetInfo.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/Parameters.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/AnalysisParameters.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/Message.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/BModel.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/WModel.h"

#include <iostream>

#include <vector>

using namespace Deconstruction;
using namespace fastjet;
using namespace std;

Deconstruction::TopModel::TopModel(Parameters &param, Flavour::id flavour)
  : Deconstruction::Model::Model(param), m_flavour(flavour) {
}

Deconstruction::TopModel::~TopModel() {
}

//double Deconstruction::TopModel::weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sumJets) {
double Deconstruction::TopModel::weight(const StoredJet &jets, fastjet::PseudoJet &sumJets) {
  double result = 0;

#ifdef DEBUGCODE
  LOG(DEBUG) << "In TopModel: " << endl;
  printoutput(jets,"input in TopModel", DEBUG);
#endif

  /////////////////////////////////////////////////////////////
  //////// Make Top decay:   ////////////////////////////////

  // checks if input is already smaller than top window
  // msp: include mass check after cross check with dave
  if ( (jets.size()<3) /*|| (offTop.m() < (topmass-delta_topmass)) */) {
#ifdef DEBUGCODE
    LOG(DEBUG) << "nr microjets < 3 or too light -> no top decay" << endl;
#endif
    return result;
  }

#ifdef DEBUGCODE
  LOG(DEBUG) << "passed lower topmass cut: " << sumJets.m() << endl;
#endif

  /////// split into left and right branches /////////////////////////

  /*
  std::vector<std::vector<fastjet::PseudoJet> > bbranches;

  int start=1;
  int end=jets.size()-1;   // msp: check if input.size()-1 or not..
  power_set_orig<fastjet::PseudoJet,std::vector<fastjet::PseudoJet>::iterator >(jets.begin(), jets.end(), bbranches, start, end);

  for(unsigned i=0; i<bbranches.size(); i++) {
    vector<PseudoJet> Wbranch(jets.size()-bbranches[i].size());
    vector<PseudoJet> bbranch(bbranches[i]);

    sort(jets.begin(), jets.end(),  lessThanIndex);
    sort(bbranch.begin(), bbranch.end(),  lessThanIndex);
    set_difference(jets.begin(), jets.end(), bbranch.begin(), bbranch.end(), Wbranch.begin(), lessThanIndex);
  */

  // possible combinations for the set in [first, last] are given by permutations of the binary
  // number with p->mom().size() bits with any number of bits activated
  unsigned int iElements = (unsigned int) jets.size();
  // iElements cannot be greater than the maximum,
  // but this has already been tested before, when generating the powerset for ISR/SB
  // there is no point in testing it again, since these sets must be smaller or equal in size
  // than the subsets of the powerset of all microjets
  unsigned long long nSets = (unsigned long long) (0x1 << iElements); // 2^iElements

  // the bits position in index i give the positions of the vector to be used
  // we are skipping the empty set combination here
  // that's why i starts at one and ends at nSets-1
  for (unsigned long long i = 1; i < nSets; ++i) {

    // reject combinations with less than two partons for the W-boson earlier
    if (numberOfSetBits(i) < 2)
      continue;

    //std::vector<fastjet::PseudoJet> Wbranch;
    //std::vector<fastjet::PseudoJet> bbranch;
    StoredJet Wbranch(jets.store());
    StoredJet bbranch(jets.store());

    // check the bits set in i and include the elements corresponding to
    // that bit position in result
    for (unsigned int k = 0; k < iElements; ++k) {
      if ( (i & (0x1 << k)) != 0) {
        Wbranch.push_back(jets.getIdx(k));
      } else {
        bbranch.push_back(jets.getIdx(k));
      }
    }

    PseudoJet topjet = sumJets;
    PseudoJet Wjet = Wbranch.sum();
    PseudoJet bjet = bbranch.sum();

    //msp: should be checked before branches are splitted
    if (std::fabs(square(topjet.m()) - square(param()[p_topmass])) >= param()[p_topmass]*delta_topmass)
      continue;

    if (std::fabs(square(Wjet.m()) - square(param()[p_wmass])) >= param()[p_wmass]*delta_wmass) continue;

    // already done previously
    //if (Wbranch.size() < 2)
    //  continue;

    // first decay factor for t -> bW

    double decayfactortop = 8.0*square(M_PI)*square(param()[p_topmass])*param()[p_topmass]*topwidth/std::atan(std::fabs(square(topjet.m()) - square(param()[p_topmass]))/(param()[p_topmass]*topwidth));

    decayfactortop /= (square(param()[p_topmass])-square(param()[p_wmass]))*(square(square(topjet.m())-square(param()[p_topmass]))+square(param()[p_topmass])*square(topwidth));

    double bweight = 0.0;
    double Wweight = 0.0;

#ifdef DEBUGCODE
    LOG(DEBUG) << "input b model: " << endl;
    printoutput(bbranch, DEBUG);

    LOG(DEBUG) << "input W model: " << endl;
    printoutput(Wbranch, DEBUG);
#endif

    if (m_flavour == Flavour::t) {
      BModel bm(param(), Flavour::b,Flavour::t);
      bm.setTopset(jets);
      bweight = bm.weight(bbranch, bjet); 
      WModel wm(param(), Flavour::Wp);
      wm.setTopset(jets);
      Wweight = wm.weight(Wbranch, Wjet);
    } else if (m_flavour == Flavour::tbar) {
      BModel bm(param(), Flavour::bbar,Flavour::tbar);
      bm.setTopset(jets);
      bweight = bm.weight(bbranch, bjet); 
      WModel wm(param(), Flavour::Wm);
      wm.setTopset(jets);
      Wweight = wm.weight(Wbranch, Wjet);
    } else {
      LOG(ERROR) << "TopModel: WRONG FLAVOR namely: " << m_flavour << endl;
    }

#ifdef DEBUGCODE
    LOG(DEBUG) << "tflavor: " << m_flavour << endl;
    LOG(DEBUG) << "decayfactortop: " << decayfactortop << endl;
    LOG(DEBUG) << "bweight: " << bweight << endl;
    LOG(DEBUG) << "W shower*decayfactorW: " << Wweight << endl;
    LOG(DEBUG) << "topweight: " << decayfactortop*Wweight*bweight << endl ;

    LOG(DEBUG) << "topweight: " << decayfactortop*Wweight*bweight << endl;
    LOG(DEBUG) << "topweight sum: " << result << endl << endl;
#endif

    result += decayfactortop*Wweight*bweight;
  }

#ifdef DEBUGCODE
  LOG(DEBUG) << "final signal weight in topmodel: " << result << endl << endl;
#endif

  return result;
}

double Deconstruction::TopModel::hamiltonian(double pTsum) {
  return 1.0;
}

