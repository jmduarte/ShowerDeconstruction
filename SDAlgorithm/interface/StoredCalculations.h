#ifndef STOREDCALCULATIONS_H
#define STOREDCALCULATIONS_H

#include "ShowerDeconstruction/SDAlgorithm/interface/StoredJet.h"
#include <map>
#include <vector>
#include "ShowerDeconstruction/SDAlgorithm/interface/Helper.h"

namespace Deconstruction {

  #define StoredKey std::vector<long long>

  class StoredCalculations {
    public:
      StoredCalculations();
      virtual ~StoredCalculations();

      bool check(StoredKey, double &w);
      void store(StoredKey, double w);
      void clear();

    private:
      std::map<StoredKey, double> m_table;
  };

}

#endif

