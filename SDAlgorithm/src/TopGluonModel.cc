#include "ShowerDeconstruction/SDAlgorithm/interface/TopGluonModel.h"

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

Deconstruction::TopGluonModel::TopGluonModel(Parameters &param)
  : Deconstruction::Model::Model(param) {
}

Deconstruction::TopGluonModel::~TopGluonModel() {
}

//double Deconstruction::TopGluonModel::weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sumJets) {
double Deconstruction::TopGluonModel::weight(const StoredJet &jets, fastjet::PseudoJet &sumJets) {

  double result = 0;

  if (sumJets.perp() <= sumJets.m())
    return 0;

#ifdef DEBUGCODE
  LOG(DEBUG) << "input of TopGluonModel: " << endl;
  printoutput(jets, DEBUG);
#endif

  //vector<PseudoJet> firstleftcol,firstrightcol,grandmother;
  //firstleftcol.clear();
  //firstrightcol.clear();
  //grandmother.clear();
  //result = make_splitting(jets, firstleftcol, firstrightcol, grandmother, Flavour::t, Flavour::noflav, Flavour::noflav, Flavour::noflav, Shower::t);

  result = start_splitting(jets, Flavour::t, Shower::t);
  m_calc.clear();

#ifdef DEBUGCODE
  LOG(DEBUG) << "input of TopGluonModel: " << endl;
  printoutput(jets, DEBUG);
  LOG(DEBUG) << "GluTopSplitweight: " << result << endl;
#endif

  return result;

}

double Deconstruction::TopGluonModel::hamiltonian(double pTsum) {
  double kH2 = square(pTsum);
  double resmass2 = square(param()[p_topmass]);
  double Hfj_sig = Cte::Npdf_signal*std::pow((Cte::pTmin2 + resmass2)/(kH2 + resmass2), Cte::Npdf_signal)/(kH2 + resmass2);
  return Hfj_sig;
}



