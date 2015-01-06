#include "ShowerDeconstruction/SDAlgorithm/interface/Storage.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/Exception.h"

Deconstruction::Storage::Storage() {
}

Deconstruction::Storage::~Storage() {
}

std::vector<fastjet::PseudoJet> Deconstruction::Storage::get(const std::vector<int> &list) const {
  std::vector<fastjet::PseudoJet> out;
  for (int i = 0; i < list.size(); ++i) {
    if (list[i] >= m_input.size())
      throw NEW_EXCEPTION("Trying to access element out of bound in Storage.");

    out.push_back(m_input[list[i]]);
  }
  return out;
}

fastjet::PseudoJet Deconstruction::Storage::sum(const std::vector<int> &list) const {
  fastjet::PseudoJet out(0,0,0,0);
  out.set_user_index(0);
  if (list.size() == 1) return m_input[list[0]];
  for (int i = 0; i < list.size(); ++i) {
    if (list[i] >= m_input.size())
      throw NEW_EXCEPTION("Trying to access element out of bound in Storage.");
    out += m_input[list[i]];
    out.set_user_index(out.user_index() + m_input[list[i]].user_index());
  }
  return out;
}

int Deconstruction::Storage::size() const {
  return m_input.size();
}

void Deconstruction::Storage::set(const std::vector<fastjet::PseudoJet> &input) {
  m_input = input;
}

std::vector<fastjet::PseudoJet> &Deconstruction::Storage::get() {
  return m_input;
}

const std::vector<fastjet::PseudoJet> &Deconstruction::Storage::get() const {
  return m_input;
}

fastjet::PseudoJet &Deconstruction::Storage::operator[](int i) {
  return m_input[i];
}

const fastjet::PseudoJet &Deconstruction::Storage::operator[](int i) const {
  return m_input[i];
}


