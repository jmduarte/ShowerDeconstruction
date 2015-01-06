#include "ShowerDeconstruction/SDAlgorithm/interface/Exception.h"

#include <string>
#include <map>
#include <iostream>

Deconstruction::Exception::Exception(const std::string &msg, int lineno, const std::string &sourceFile)
  : m_msg(msg), m_lineno(lineno), m_sourceFile(sourceFile) {
}

Deconstruction::Exception &Deconstruction::Exception::setParam(const std::string &key, const std::string &value) {
  m_param[key] = value;
  return *this;
}

std::string Deconstruction::Exception::what() const {
  std::string reason;
  reason += "Deconstruction::Exception in ";
  reason += m_sourceFile;
  reason += ", line ";
  reason += m_lineno;
  reason += ":\n";
  reason += fixedMsg();

  return reason;
}

std::string Deconstruction::Exception::fixedMsg() const {
  std::string msg = m_msg;

  for (std::map<std::string, std::string>::const_iterator it = m_param.begin(); it != m_param.end(); ++it) {
    std::string lookFor = "[";
    lookFor += it->first;
    lookFor += "]";
    findAndReplace(msg, lookFor, it->second);
  }
  return msg;
}

void Deconstruction::Exception::findAndReplace(std::string msg, const std::string &key, const std::string &value) const {
  size_t pos = 0;
  while ((pos = msg.find(key, pos)) != std::string::npos) {
     msg.replace(pos, msg.length(), value);
     pos += key.length();
  }
}

