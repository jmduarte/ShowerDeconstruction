#ifndef BACKGROUNDMODEL_H
#define BACKGROUNDMODEL_H

#include "ShowerDeconstruction/SDAlgorithm/interface/Model.h"
#include <fastjet/PseudoJet.hh>
#include <vector>

#include "ShowerDeconstruction/SDAlgorithm/interface/Parameters.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/StoredJet.h"

namespace Deconstruction {

  class BackgroundModel : public Model {
    public:
      BackgroundModel(Parameters &param, Flavour::id flav = Flavour::g);
      virtual ~BackgroundModel();

      //double weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sum);
      double weight(const StoredJet &jets, fastjet::PseudoJet &sum);
      double hamiltonian(double pTsum = 0);

    private:
      Flavour::id m_flav;
  };

}

#endif

