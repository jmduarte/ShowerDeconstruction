#include "ShowerDeconstruction/SDAlgorithm/interface/JetInfo.h"

Deconstruction::JetInfo::JetInfo()
  : fastjet::PseudoJet::UserInfoBase(), m_i(0), m_user_index(0) {
}

Deconstruction::JetInfo::JetInfo(unsigned int i, int user_index)
  : fastjet::PseudoJet::UserInfoBase(), m_i(i), m_user_index(user_index) {
}

Deconstruction::JetInfo::~JetInfo() {
}

unsigned int &Deconstruction::JetInfo::i() {
  return m_i;
}

int &Deconstruction::JetInfo::user_index() {
  return m_user_index;
}

const unsigned int &Deconstruction::JetInfo::i() const {
  return m_i;
}

const int &Deconstruction::JetInfo::user_index() const {
  return m_user_index;
}

bool Deconstruction::lessThanIndex(const fastjet::PseudoJet &i, const fastjet::PseudoJet &j) {
  return (i.user_info<Deconstruction::JetInfo>().i() < j.user_info<Deconstruction::JetInfo>().i());
}

