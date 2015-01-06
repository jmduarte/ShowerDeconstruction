#include "ShowerDeconstruction/SDAlgorithm/interface/WDecayModel.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/WModel.h"

#include <iostream>

#include <vector>
#include <fastjet/PseudoJet.hh>
#include "ShowerDeconstruction/SDAlgorithm/interface/Helper.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/Parameters.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/AnalysisParameters.h"

using namespace Deconstruction;
using namespace fastjet;
using namespace std;

using std::endl;

Deconstruction::WDecayModel::WDecayModel(Parameters &param, bool includeAngleFactor)
  : Deconstruction::Model::Model(param), m_includeAngleFactor(includeAngleFactor) {
}

Deconstruction::WDecayModel::~WDecayModel() {
}

//double Deconstruction::WDecayModel::weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sumJets) {
double Deconstruction::WDecayModel::weight(const StoredJet &jets, fastjet::PseudoJet &sumJets) {

  double result = 0;

  if (sumJets.perp() <= sumJets.m())
    return 0;

#ifdef DEBUGCODE
  LOG(DEBUG) << "input of WDecayModel: " << endl;
  printoutput(jets, DEBUG);
#endif

  if ( (jets.size()<2) ) {
#ifdef DEBUGCODE
    LOG(DEBUG) << "nr microjets < 2 or too light -> no W decay" << endl;
#endif
    return result;
  }

  PseudoJet Wjet = sumJets;
  if (std::fabs(square(Wjet.m()) - square(param()[p_wmass])) >= param()[p_wmass]*delta_wmass)
    return result;

  m_calc.clear();
  WModel wm(param(), Flavour::Wp, m_includeAngleFactor);
  wm.setTopset(jets); // if it is a top
  result = wm.weight(jets, Wjet);

  m_calc.clear();

#ifdef DEBUGCODE
  LOG(DEBUG) << "input of WDecayModel: " << endl;
  printoutput(jets, DEBUG);
  LOG(DEBUG) << "weight: " << result << endl;
#endif

  return result;

}

double Deconstruction::WDecayModel::hamiltonian(double pTsum) {
  double kH2 = square(pTsum);
  double resmass2 = square(param()[p_wmass]);
  double Hfj_sig = Cte::Npdf_signal*std::pow((Cte::pTmin2 + resmass2)/(kH2 + resmass2), Cte::Npdf_signal)/(kH2 + resmass2);
  return Hfj_sig;
}


