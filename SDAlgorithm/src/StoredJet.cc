#include "ShowerDeconstruction/SDAlgorithm/interface/StoredJet.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/Exception.h"

Deconstruction::StoredJet::StoredJet()
  :m_store(0) {
}

Deconstruction::StoredJet::StoredJet(const Storage &store)
  :m_store(&store) {
}

Deconstruction::StoredJet::StoredJet(const Storage *store)
  :m_store(store) {
}

Deconstruction::StoredJet::StoredJet(const Storage &store, std::vector<int> list)
  :m_store(&store), m_list(list) {
}

Deconstruction::StoredJet::StoredJet(const StoredJet &s)
  :m_store(s.m_store), m_list(s.m_list) {
}

Deconstruction::StoredJet::~StoredJet() {
}

std::vector<int> Deconstruction::StoredJet::getList() const {
  return m_list;
}

const Deconstruction::Storage &Deconstruction::StoredJet::store() const {
  if (!m_store) throw NEW_EXCEPTION("No Storage available!");
  return *m_store;
}

int Deconstruction::StoredJet::getIdx(int i) const {
  return m_list[i];
}

const fastjet::PseudoJet &Deconstruction::StoredJet::operator[](int i) const {
  if (!m_store) throw NEW_EXCEPTION("No Storage available!");
  int j = m_list[i];
  return (*m_store)[j];
}

const fastjet::PseudoJet &Deconstruction::StoredJet::at(int i) const {
  if (!m_store) throw NEW_EXCEPTION("No Storage available!");
  int j = m_list.at(i);
  return (*m_store)[j];
}

int Deconstruction::StoredJet::size() const {
  return m_list.size();
}

fastjet::PseudoJet Deconstruction::StoredJet::sum() const {
  if (!m_store) throw NEW_EXCEPTION("No Storage available!");
  return m_store->sum(m_list);
}

void Deconstruction::StoredJet::push_back(int j) {
  m_list.push_back(j);
}

Deconstruction::StoredJet::operator std::vector<fastjet::PseudoJet>() const {
  if (!m_store) throw NEW_EXCEPTION("No Storage available!");
  return m_store->get(m_list);
}

