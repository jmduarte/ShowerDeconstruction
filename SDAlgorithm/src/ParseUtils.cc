#include "ShowerDeconstruction/SDAlgorithm/interface/ParseUtils.h"

#include <getopt.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <vector>
#include <algorithm>

void buildOptions(struct extendedOption *extOpt, struct option *&opt, std::string &shortOpts) {
  int numberOfOpts = 0;
  if (extOpt == 0) {
    opt = 0;
    return;
  }
  struct extendedOption *copyOfExtOpt = extOpt;
  while (copyOfExtOpt->name != 0) {
    ++copyOfExtOpt;
    ++numberOfOpts;
  }
  ++numberOfOpts; // for the zero-entry
  std::vector<int> check;
  opt = new option[numberOfOpts];
  std::stringstream ss;
  ss.str("");
  for (int i = 0; i < numberOfOpts; ++i) {
    if (std::find(check.begin(), check.end(), extOpt[i].val) != check.end()) {
      std::cout << "buildOptions: WARNING -- Repeated option identifier (identifier \'" << (char) extOpt[i].val << "\' [" << (int) extOpt[i].val << "])." << std::endl; 
    }

    if (extOpt[i].has_arg == no_argument) { // set a flag
    } else if (extOpt[i].has_arg == required_argument) { // have arguments -- pass them on shortOpts
      ss << (char) extOpt[i].val << ":";
      check.push_back(extOpt[i].val);
    }
    opt[i].name = extOpt[i].name;
    opt[i].has_arg = extOpt[i].has_arg;
    opt[i].flag = extOpt[i].flag;
    opt[i].val = extOpt[i].val;
  }
  shortOpts = ss.str();
}

void dumpHelp(const std::string &nameOfProgram, struct extendedOption *extOpt, const std::string &description) {
  std::cout << description << std::endl;
  std::cout << "Usage: " << nameOfProgram << " [options]" << std::endl;
  std::cout << std::endl;
  std::cout << "where [options] can be a combination of:" << std::endl;
  while (extOpt->name != 0) {
    if (extOpt->has_arg == no_argument) {
      std::cout << "  --" << std::left << std::setw(38) << extOpt->name << "    " << extOpt->description << std::endl;
    } else if (extOpt->has_arg == required_argument) {
      std::cout << "  --" << std::left << std::setw(25) << extOpt->name << "|-" << (char) extOpt->val << "   (arg)      " << extOpt->description << "  (default: ";
      if (extOpt->type == extendedOption::eOTFloat) {
        float *f = (float *) extOpt->pointerToValue;
        std::cout << *f << " [float]";
      } else if (extOpt->type == extendedOption::eOTString) {
        std::string *s = (std::string *) extOpt->pointerToValue;
        std::cout << "\"" << *s << "\" [string]";
      } else if (extOpt->type == extendedOption::eOTInt) {
        int *k = (int *) extOpt->pointerToValue;
        std::cout << *k << " [int]";
      } else {
        std::cout << "unknown";
      }
      std::cout << ")" << std::endl;
    }
    ++extOpt;
  }
  std::cout << std::endl;
  std::cout << "Choose wisely :)" << std::endl;
  std::cout << std::endl;
}

void dumpOptions(struct extendedOption *extOpt) {
  while (extOpt->name != 0) {
    std::cout << std::left << std::setw(25) << extOpt->name << " = ";
    if (extOpt->type == extendedOption::eOTFloat) {
      float *f = (float *) extOpt->pointerToValue;
      std::cout << *f << " [float]";
    } else if (extOpt->type == extendedOption::eOTString) {
      std::string *s = (std::string *) extOpt->pointerToValue;
      std::cout << "\"" << *s << "\" [string]";
    } else if (extOpt->type == extendedOption::eOTInt) {
      int *k = (int *) extOpt->pointerToValue;
      std::cout << *k << " [int]";
    } else {
      std::cout << "unknown";
    }
    std::cout << std::endl;
    ++extOpt;
  }
  std::cout << std::endl;
}

bool parseArguments(int argc, char **argv, struct extendedOption *extOpt) {
  struct option *long_options;
  std::string shortOpts = "";
  buildOptions(extOpt, long_options, shortOpts);

  int c;
  while (true) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long(argc, argv, shortOpts.c_str(), long_options, &option_index);
    if (c == -1) break; // end of parsing arguments

    // treat unrecognised option
    if (c == '?') return false;

    if (c == 0) { // flag setting option -- they are set automatically
    } else {
      struct extendedOption *eO = extOpt;
      while (eO->name != 0) { // look for its identifier in the array
        if (c == eO->val) { // if I found the identifier
          //std::cout << "Found ID " << (char) c << " for " << eO->name << " with contents \"" << optarg << "\", type == " << (int) eO->type << std::endl;
          // convert to the appropriate type and copy to the contents of the pointer storing it
          if (eO->type == extendedOption::eOTFloat) {
            std::stringstream ss;
            float *f = (float *) eO->pointerToValue;
            ss << optarg;
            ss >> *f;
          } else if (eO->type == extendedOption::eOTString) {
            std::stringstream ss;
            std::string s;
            ss << optarg;
            ss >> s;
            std::string *sp = (std::string *) eO->pointerToValue;
            *sp = s;
          } else if (eO->type == extendedOption::eOTInt) {
            std::stringstream ss;
            int *k = (int *) eO->pointerToValue;
            ss << optarg;
            ss >> *k;
          }
          break;
        } // if I have found the identifier
        ++eO;
      } // end of loop for finding the identifier of this option in the extended options array
    } // if it is NOT a flag setting option
  } // until they are all parsed

  delete [] long_options;

  return true;
}

