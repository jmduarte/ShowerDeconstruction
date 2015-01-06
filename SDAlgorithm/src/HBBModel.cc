#include "ShowerDeconstruction/SDAlgorithm/interface/HBBModel.h"

#include <vector>
#include <fastjet/PseudoJet.hh>
#include "ShowerDeconstruction/SDAlgorithm/interface/Helper.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/Parameters.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/AnalysisParameters.h"

using namespace Deconstruction;

Deconstruction::HBBModel::HBBModel(Parameters &param)
  : Deconstruction::Model::Model(param) {
}

Deconstruction::HBBModel::~HBBModel() {
}

double Deconstruction::HBBModel::weight(const StoredJet &jets, fastjet::PseudoJet &sum) {

  double signalweight = 0;

  // there must be 2 microjets for the Higgs decay
  if (jets.size() < 2) return signalweight;

  if (std::fabs(jets.sum().m() - param()[p_Higgsmass]) > param()[p_delta_Higgsmass]) {
    return signalweight;
  } else { //assign weight H*exp(-S):
    signalweight = 4.0*square(M_PI)/param()[p_Higgsmass]/param()[p_delta_Higgsmass];
    signalweight /= pow(pow(jets.sum().m(), 2) - pow(param()[p_Higgsmass], 2), 2) + pow(param()[p_Higgsmass], 2)*pow(param()[p_delta_Higgsmass], 2);
  }

  double w = 0;

  unsigned int iElements = (unsigned int) jets.size();
  unsigned long long nSets = (unsigned long long) (0x1 << iElements); // 2^iElements
  StoredJet empty(jets.store());
  for (unsigned long long i = 1; i < nSets; ++i) {
    StoredJet ubranch(jets.store());
    StoredJet dbranch(jets.store());
    for (unsigned int k = 0; k < iElements; ++k) {
      if ( (i & (0x1 << k)) != 0) {
        ubranch.push_back(jets.getIdx(k));
      } else {
        dbranch.push_back(jets.getIdx(k));
      }
    }

    double valueleft = 0;
    double valueright = 0;
    m_calc.clear();
    valueleft = make_splitting(ubranch, empty, dbranch, jets, Flavour::b, Flavour::noflav, Flavour::noflav, Flavour::noflav, Shower::b);
    m_calc.clear();
    valueright = make_splitting(dbranch, ubranch, empty, jets, Flavour::bbar, Flavour::noflav, Flavour::noflav, Flavour::noflav, Shower::b);

    w += valueleft*valueright;
  }

  return w*signalweight;
}

double Deconstruction::HBBModel::hamiltonian(double pTsum) {
  double kH2 = square(pTsum);
  double Higgsmass2 = square(param()[p_Higgsmass]);

  double Hfj_sig = Cte::Npdf_signal*std::pow((Cte::pTmin2 + Higgsmass2)/(kH2 + Higgsmass2),Cte::Npdf_signal)/(kH2 + Higgsmass2);

  return Hfj_sig;
}

