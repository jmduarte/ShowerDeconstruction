#include "ShowerDeconstruction/SDAlgorithm/interface/Message.h"
#include <streambuf>

Deconstruction::nullbuf Deconstruction::MsgInterface::_nullbuf = Deconstruction::nullbuf();

Deconstruction::nullstream Deconstruction::MsgInterface::nullStream = Deconstruction::nullstream();
Deconstruction::Message Deconstruction::MsgInterface::msg = Deconstruction::Message();

Deconstruction::nullbuf::nullbuf()
  : std::streambuf() {
}

Deconstruction::nullbuf::nullbuf(const Deconstruction::nullbuf &c)
  : std::streambuf() {
}

Deconstruction::nullstream::nullstream()
: std::ostream(&Deconstruction::MsgInterface::_nullbuf) {
}

Deconstruction::nullstream::nullstream(const Deconstruction::nullstream &c)
: std::ostream(&Deconstruction::MsgInterface::_nullbuf) {
}

std::ostream &Deconstruction::nullstream::m(const Deconstruction::MsgLevel x) {
  return *this;
}

Deconstruction::Message::Message()
  : std::ostream(std::cout.rdbuf()) {
  m_level = INFO;
}

Deconstruction::Message::Message(const Deconstruction::Message &c)
  : std::ostream(std::cout.rdbuf()) {
  m_level = c.m_level;
}

Deconstruction::Message::~Message() {
}

void Deconstruction::Message::setLevel(const Deconstruction::MsgLevel level) {
  m_level = level;
}

const Deconstruction::MsgLevel &Deconstruction::Message::level() const {
  return m_level;
}

std::ostream &Deconstruction::Message::m(const Deconstruction::MsgLevel level) {
  if (level >= m_level)
    return *this;
  return Deconstruction::MsgInterface::nullStream;
}

