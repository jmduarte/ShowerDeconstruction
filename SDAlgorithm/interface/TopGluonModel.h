#ifndef TOPGLUONMODEL_H
#define TOPGLUONMODEL_H

#include "ShowerDeconstruction/SDAlgorithm/interface/Model.h"
#include <fastjet/PseudoJet.hh>
#include <vector>

#include "ShowerDeconstruction/SDAlgorithm/interface/Parameters.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/StoredJet.h"

namespace Deconstruction {

  class TopGluonModel : public Model {
    public:
      TopGluonModel(Parameters &param);
      virtual ~TopGluonModel();

      //double weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sum);
      double weight(const StoredJet &jets, fastjet::PseudoJet &sum);
      double hamiltonian(double aScale = 0);

  };

}

#endif

