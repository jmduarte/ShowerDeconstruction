#include "ShowerDeconstruction/SDAlgorithm/interface/Parameters.h"

#include <string>
#include <fstream>
#include <sstream>

#include "ShowerDeconstruction/SDAlgorithm/interface/Exception.h"

#include <algorithm>

Deconstruction::Parameters::Parameters() {
}

Deconstruction::Parameters::Parameters(const std::string &input) {
  read(input);
}

Deconstruction::Parameters::~Parameters() {
}

int Deconstruction::Parameters::read(const std::string &input) {
  int nReadKeys = 0;

  std::ifstream infile(input.c_str());
  if (!infile.is_open()) {
    throw NEW_EXCEPTION("Failed to open [file] with parameters.").setParam("file", input);
  }
  while (infile.good()) {
    std::stringstream ss;
    infile.get(*ss.rdbuf(), '\n');
    if (infile.get() == EOF)
      break;

    std::string key;
    double value;
    ss >> key >> value;
    (*this)[key] = value;
  }
}

int Deconstruction::Parameters::insert(const std::string &key, const double value) {
  m_param.push_back(value);
  m_keys.push_back(key);
  return m_param.size() - 1;
}

double &Deconstruction::Parameters::operator[](const std::string &key) {
  std::vector<std::string>::iterator it = std::find(m_keys.begin(), m_keys.end(), key);
  if (it == m_keys.end()) {
    return (*this)[insert(key, 0)];
  }
  int index = (int) (it - m_keys.begin());
  return (*this)[index];
}

double &Deconstruction::Parameters::operator[](const int key) {
  return m_param[key];
}

