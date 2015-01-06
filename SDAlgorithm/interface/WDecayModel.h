#ifndef WDECAYMODEL_H
#define WDECAYMODEL_H

#include "ShowerDeconstruction/SDAlgorithm/interface/Model.h"
#include <fastjet/PseudoJet.hh>
#include <vector>

#include "ShowerDeconstruction/SDAlgorithm/interface/Parameters.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/StoredJet.h"

namespace Deconstruction {

  class WDecayModel : public Model {
    public:
      WDecayModel(Parameters &param, bool includeAngleFactor = false);
      virtual ~WDecayModel();

      //double weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sum);
      double weight(const StoredJet &jets, fastjet::PseudoJet &sum);
      double hamiltonian(double aScale = 0);

    protected:
      bool m_includeAngleFactor;

  };

}

#endif

