#ifndef HBBMODEL_H
#define HBBMODEL_H

#include "ShowerDeconstruction/SDAlgorithm/interface/Model.h"
#include <fastjet/PseudoJet.hh>
#include <vector>

#include "ShowerDeconstruction/SDAlgorithm/interface/Parameters.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/StoredJet.h"

namespace Deconstruction {

  class HBBModel : public Model {
    public:
      HBBModel(Parameters &param);
      virtual ~HBBModel();

      //double weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sum);
      double weight(const StoredJet &jets, fastjet::PseudoJet &sum);
      double hamiltonian(double pTsum = 0);

  };

}

#endif

