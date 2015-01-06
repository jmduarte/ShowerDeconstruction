#include "ShowerDeconstruction/SDAlgorithm/interface/BackgroundModel.h"

#include <vector>
#include <fastjet/PseudoJet.hh>
#include "ShowerDeconstruction/SDAlgorithm/interface/Helper.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/Parameters.h"

using namespace Deconstruction;

Deconstruction::BackgroundModel::BackgroundModel(Parameters &param, Flavour::id flav)
  : Deconstruction::Model::Model(param), m_flav(flav) {
}

Deconstruction::BackgroundModel::~BackgroundModel() {
}

//double Deconstruction::BackgroundModel::weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sum) {
double Deconstruction::BackgroundModel::weight(const StoredJet &jets, fastjet::PseudoJet &sum) {
  // Background model ---> calls make_splitting with gluon

  // the Background model consists simply of a FSRshower
  // keep in mind that the H for the generation of the 
  // fat jet is applied in Deconstruct::deconstruct()

  // check for theta cut if there is no grandmother:
  if (sum.perp() <= sum.m())
    return 0;

  double w = start_splitting(jets, m_flav, Shower::QCD);
  m_calc.clear();

  return w;
}

double Deconstruction::BackgroundModel::hamiltonian(double pTsum) {
  double kH2 = square(pTsum);
  double Hfj_sig = Cte::Npdf_signal*std::pow(Cte::pTmin2/kH2, Cte::Npdf_signal)/kH2;
  return Hfj_sig;
}

