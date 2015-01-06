#include "ShowerDeconstruction/SDAlgorithm/interface/StoredCalculations.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/StoredJet.h"

#include <algorithm>
#include <cmath>

using namespace Deconstruction;

Deconstruction::StoredCalculations::StoredCalculations() {
}
Deconstruction::StoredCalculations::~StoredCalculations() {
}

bool Deconstruction::StoredCalculations::check(StoredKey k, double &w) {
  std::map<StoredKey, double>::iterator it = m_table.find(k);
  if (it != m_table.end()) {
    w = it->second;
    return true;
  }
  return false;
}

void Deconstruction::StoredCalculations::store(StoredKey k, double w) {
  m_table.insert(std::pair<StoredKey, double>(k, w));
}

void Deconstruction::StoredCalculations::clear() {
  m_table.clear();
}

