#include "ShowerDeconstruction/SDAlgorithm/interface/Helper.h"

#include <vector>
#include <cmath>
#include "ShowerDeconstruction/SDAlgorithm/interface/Exception.h"
#include <fastjet/PseudoJet.hh>

#include "ShowerDeconstruction/SDAlgorithm/interface/JetInfo.h"

#include <iostream>
#include <iomanip>
#include "stdlib.h"

#include <sstream>
#include <string>

#include "ShowerDeconstruction/SDAlgorithm/interface/StoredJet.h"

std::string d_to_string(unsigned int x) {
  std::stringstream ss;
  ss << x;
  return ss.str();
}

using namespace Deconstruction;
using std::cout;
using std::endl;
using std::max;
using std::setw;
using std::vector;

const double Deconstruction::Cte::pi          = 3.141592653589793;
const double Deconstruction::Cte::TR          = 0.5;
const double Deconstruction::Cte::CA          = 3.0;
const double Deconstruction::Cte::CF          = (square(Cte::CA) - 1.0)/2.0/Cte::CA;
const double Deconstruction::Cte::nf          = 5;
const double Deconstruction::Cte::b0          = (33.0 - 2.0*Cte::nf)/12.0/Cte::pi;

// default is 0.5
const double Deconstruction::Cte::cnp         = 1.0;
// default is 0.5
const double Deconstruction::Cte::kappa_np2   = 4.0;
// default is 1
const double Deconstruction::Cte::kappa_p2    = 4.0;
// default is 2
const double Deconstruction::Cte::nnp         = 1.5;
// default is 2
const double Deconstruction::Cte::cr          = 2.0;
// default is 4
const double Deconstruction::Cte::nr          = 1.0;

// minpT for top quark
const double Deconstruction::Cte::pTmin2      = square(400.0);
// new default is 1.4
const double Deconstruction::Cte::Npdf_backg  = 2.0;
const double Deconstruction::Cte::Npdf_signal = 2.0;

// sudakov top factor a
const double Deconstruction::Cte::alphattg    = 1.0;
// sudakov top factor b
const double Deconstruction::Cte::betattg     = 1.0;

const double Deconstruction::Cte::smalldouble = 1e-7;

void Deconstruction::dummy() {
  std::vector<fastjet::PseudoJet>::iterator it;
  std::vector<std::vector<fastjet::PseudoJet> > v;
  std::vector<fastjet::PseudoJet> x;
  power_set_orig<fastjet::PseudoJet,  std::vector<fastjet::PseudoJet>::iterator >(it, it, v, 0, 0);
  find_subset_orig(it, it, 1, x, v);
}

template <class T , class Iter>
void  Deconstruction::find_subset_orig( Iter first , Iter last , int n , vector<T>& foo, vector<vector<T> > & result ) {

  if (n == 0) {
    result.push_back(foo);        
  } else {
    for ( Iter iter = first ; iter != last ; ++iter ) {
      foo.push_back( *iter ) ;
      ++iter ;
      find_subset_orig( iter , last , n-1 , foo, result ) ;
      --iter ;
      foo.pop_back() ;
    }
  }
}

template <class T , class Iter>
void  Deconstruction::power_set_orig( Iter first , Iter last, vector < vector < T > > & result, int start, int end) {
  vector<T>  sets ;
  for ( int i = start ; i <= end ; ++i ) {
    sets.resize(0) ;
    find_subset_orig( first , last , i , sets, result ) ;
  }
}


void Deconstruction::copyJI(const fastjet::PseudoJet &from, fastjet::PseudoJet &to) {
  to.set_user_info(new JetInfo(from.user_info<JetInfo>().i(), from.user_info<JetInfo>().user_index()));
}

/// does the actual work for printing out a jet
void Deconstruction::printoutput(const std::vector<fastjet::PseudoJet> &leftbranch) {
  LOG(DEBUG) << " { ";
  for(unsigned ii=0; ii<leftbranch.size(); ii++) {
    LOG(DEBUG) << leftbranch[ii].user_info<JetInfo>().i() << " ";
  }
  LOG(DEBUG) << " } " << endl;
}

void Deconstruction::printoutput(const std::vector<fastjet::PseudoJet> &leftbranch, MsgLevel x) {
  LOG(x) << " { ";
  for(unsigned ii=0; ii<leftbranch.size(); ii++) {
    LOG(x) << leftbranch[ii].user_info<JetInfo>().i() << " ";
  }
  LOG(x) << " } " << endl;
}
void Deconstruction::printoutput(const std::vector<fastjet::PseudoJet> &leftbranch, const std::string &s) {
  LOG(DEBUG) << s<< endl;
  LOG(DEBUG) << " { ";
  for(unsigned ii=0; ii<leftbranch.size(); ii++) {
    LOG(DEBUG) << leftbranch[ii].user_info<JetInfo>().i() << " ";
  }
  LOG(DEBUG) << " } " << endl;
}

void Deconstruction::printoutput(const std::vector<fastjet::PseudoJet> &leftbranch, const std::string &s, MsgLevel x) {
  if (x == FORCE) {
    std::cout << s<< endl;
    std::cout << " { ";
    for(unsigned ii=0; ii<leftbranch.size(); ii++) {
      std::cout << leftbranch[ii].user_info<JetInfo>().i() << " ";
    }
    std::cout << " } " << endl;
    return;
  }
  LOG(x) << s<< endl;
  LOG(x) << " { ";
  for(unsigned ii=0; ii<leftbranch.size(); ii++) {
    LOG(x) << leftbranch[ii].user_info<JetInfo>().i() << " ";
  }
  LOG(x) << " } " << endl;
}

void Deconstruction::printoutput(const StoredJet &leftbranch) {
  std::vector<fastjet::PseudoJet> lj = (std::vector<fastjet::PseudoJet>) leftbranch;
  printoutput(lj);
}

void Deconstruction::printoutput(const StoredJet &leftbranch, MsgLevel x) {
  std::vector<fastjet::PseudoJet> lj = (std::vector<fastjet::PseudoJet>) leftbranch;
  printoutput(lj, x);
}
void Deconstruction::printoutput(const StoredJet &leftbranch, const std::string &s) {
  std::vector<fastjet::PseudoJet> lj = (std::vector<fastjet::PseudoJet>) leftbranch;
  printoutput(lj, s);
}

void Deconstruction::printoutput(const StoredJet &leftbranch, const std::string &s, MsgLevel x) {
  std::vector<fastjet::PseudoJet> lj = (std::vector<fastjet::PseudoJet>) leftbranch;
  printoutput(lj, s, x);
}

void Deconstruction::printjet (const fastjet::PseudoJet &jet) {
  LOG(DEBUG) << "E, px, py, pz = "
       << " " << setw(10) << jet.e()
       << " " << setw(10) << jet.px()
       << " " << setw(10) << jet.py()
       << " " << setw(10) << jet.pz() << endl;
}

void Deconstruction::printjet (const fastjet::PseudoJet &jet, const std::string &s, Deconstruction::MsgLevel x) {
  if (x == FORCE) {
    std::cout << s <<endl;
    std::cout << "M, Pt, y, phi = "
       << " " << setw(10) << jet.m()
       << " " << setw(10) << jet.perp()
       << " " << setw(10) << jet.rap()
       << " " << setw(10) << jet.phi() << endl;
    return;
  }
  LOG(x) << s <<endl;
  LOG(x) << "M, Pt, y, phi = "
       << " " << setw(10) << jet.m()
       << " " << setw(10) << jet.perp()
       << " " << setw(10) << jet.rap()
       << " " << setw(10) << jet.phi() << endl;
}

// This calculates the set of all combination of elements in v
// This algorithm takes O(n * 2^n)
template <class T>
void Deconstruction::powerset(const std::vector<T> &v,
                              std::vector<std::vector<T> > &result,
                              const unsigned int minElements,
                              const unsigned int maxElements) {
  // possible combinations for the set in [first, last] are given by permutations of the binary
  // number with v.size() bits with any number of bits activated
  // unsigned long long is 8 bytes long in most computers = 8*8 bits = 64 bits
  // this function works well as long as the maximum number of elements in the set is 64
  // (for more elements, a cleverer way can be thought of using more than one unsigned long long)
  unsigned int nElements = (unsigned int) v.size();
  if (maxElements > sizeof(unsigned long long)*8) {
    throw NEW_EXCEPTION("Using powerset function with maxElements = [maxElements] and number of elements \
                         in vector = [nElements]. The maximum size of sets to be used in the powerset \
                         function is [maxsize], given by the size in bits of the unsigned long long in \
                         this computer.")
                         .setParam("maxsize", d_to_string(sizeof(unsigned long long)*8))
                         .setParam("maxElements", d_to_string(maxElements))
                         .setParam("nElements", d_to_string(nElements));
  }
  unsigned int iElements = std::min(maxElements, nElements);
  unsigned long long nSets = (unsigned long long) (0x1 << iElements); // 2^iElements

  // the bits position in index i give the positions of the vector to be used
  for (unsigned long long i = 0; i <= nSets; ++i) {
    // calculate if the sum of bits in i is >= minElements
    if (minElements >= 1) {
      // ok, we really need to test this
      unsigned int sumBits = numberOfSetBits(i);
      if (sumBits < minElements)
        continue;
    }
    // } else { } // no minimum number of elements ... ignore the previous check

    // check the bits set in i and include the elements corresponding to
    // that bit position in result
    // creating a vector in result itself and then changing it saves copy time
    result.push_back(std::vector<T>());
    for (unsigned int k = 0; k < iElements; ++k) {
      if ( (i & (0x1 << k)) != 0)
        result.back().push_back(v[k]);
    }
  }
}

unsigned int Deconstruction::numberOfSetBits(unsigned long long i) {
  // use popcount algorithm to count number of bits set in i
  // some processors have an instruction that does that using a lookup table
  // if this is not available gcc has a very efficient implementation of this
  // This requires gcc >= 3.4
  return __builtin_popcount(i);
}

double Deconstruction::square(const double &x) {
  return x*x;
}

fastjet::PseudoJet Deconstruction::sum(const std::vector<fastjet::PseudoJet> &v) {
  fastjet::PseudoJet result(0,0,0,0);
  for (int i = 0; i < v.size(); ++i) {
    result += v[i];
  }
  return result;
}

double Deconstruction::square_p1minusp2(const fastjet::PseudoJet & p1, const fastjet::PseudoJet & p2) {
  return ((p1.e()-p2.e())*(p1.e()-p2.e()) - (p1.px()-p2.px())*(p1.px()-p2.px()) - (p1.py()-p2.py())*(p1.py()-p2.py()) - (p1.pz()-p2.pz())*(p1.pz()-p2.pz()));
}

double Deconstruction::scalprod(const fastjet::PseudoJet & p1, const fastjet::PseudoJet & p2) {
  return (p1.e()*p2.e() - p1.px()*p2.px() - p1.py()*p2.py() - p1.pz() * p2.pz());
}

double Deconstruction::alphas(const double scale2) {
  double mZ2= square(91.1876);
  double alphasMZ = 0.118;
  
  double als = alphasMZ/(1.0 + alphasMZ*Cte::b0*std::log(scale2/mZ2));
  return als;
}

double dot(const fastjet::PseudoJet &x, const fastjet::PseudoJet &y) {
  return x.e()*y.e() - x.px()*y.px() - x.py()*y.py() - x.pz()*y.pz();
}

