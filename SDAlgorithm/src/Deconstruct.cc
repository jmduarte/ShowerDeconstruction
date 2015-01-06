#include "ShowerDeconstruction/SDAlgorithm/interface/Deconstruct.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/Parameters.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/AnalysisParameters.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/Model.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/Exception.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/JetInfo.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/Helper.h"

#include <algorithm>

#include <fastjet/PseudoJet.hh>
#include "ShowerDeconstruction/SDAlgorithm/interface/StoredJet.h"

#include <map>

#include <string>

using namespace Deconstruction;
using namespace std;
using namespace fastjet;

Deconstruction::Deconstruct::Deconstruct(Parameters &param, Model &signal, Model &background, Model &isr)
  : m_param(param), m_signal(signal), m_background(background), m_isr(isr) {
  std::cout << std::endl;
  std::cout << "Shower Deconstruction" << std::endl;
  std::cout << "=======================" << std::endl;
  std::cout << "If you use this code, please cite the following papers:" << std::endl;
  std::cout << "D. E. Soper and M. Spannowsky, ``Finding physics signals with shower deconstruction,''  Phys. Rev. D 84 (2011) 074002 [arXiv:1102.3480 [hep-ph]]." << std::endl;
  std::cout << "D. E. Soper and M. Spannowsky, ``Finding top quarks with shower deconstruction,'' Phys. Rev. D 87 (2013) 5,  054012 [arXiv:1211.3140 [hep-ph]]." << std::endl;
  std::cout << "In case of questions, please contact Danilo.Enoque.Ferreira.De.Lima@cern.ch or Michael.Spannowsky@cern.ch" << std::endl;
  std::cout << std::endl;
}

Deconstruction::Deconstruct::~Deconstruct() {
}

Model &Deconstruction::Deconstruct::signal() {
  return m_signal;
}

Model &Deconstruction::Deconstruct::background() {
  return m_background;
}

Model &Deconstruction::Deconstruct::isr() {
  return m_isr;
}

Parameters &Deconstruction::Deconstruct::param() {
  return m_param;
}

double Deconstruction::Deconstruct::deconstruct(std::vector<fastjet::PseudoJet> &input,
                              double &wSignal, double &wBackground) {

  m_signalWeight.clear();

  Model::initlevel = -1;// to agree with dave's numbering scheme
  Model::binitlevel = -1;
  Model::tinitlevel = -1;
  Model::Winitlevel = -1;

  Model::level=Model::initlevel;
  Model::btopshower_level=Model::binitlevel;
  Model::tshower_level=Model::tinitlevel;
  Model::Wshower_level=Model::Winitlevel;

  // create an ordering for the input microjets
  // and calculate the fatjet momentum
  fastjet::PseudoJet fatjet(0,0,0,0);
  for (unsigned int i = 0; i < input.size(); ++i) {
    fatjet += input[i];
    // the new UserInfoBase type will be owned by PseudoJet, which deletes it when it goes out of scope
    input[i].set_user_info(new JetInfo(i, input[i].user_index()));
#ifdef DEBUGCODE
    LOG(DEBUG) << "mom " << std::endl;
    printjet(input[i]);
#endif
  }
  // no sort necessary so far ...
  //std::sort(input.begin(), input.end(), lessThanIndex);

  // according to eq. 3.1: Qsquare = pT_fatjet^2 + m_fatjet^2
  double Qsquare = square(fatjet.perp()) + square(fatjet.m());

  // store the input vector
  Storage store;
  store.set(input);

  wBackground = 0;
  wSignal = 0;

  // this embeds the powerset function to calculate ISR/SB combination
  // it will take much less time and memory to do this calculation inline, instead of creating
  // and storing all combinations (they are uncorrelated, so the following can be used)

  // This calculates the set of all combination of elements in v
  // This algorithm takes O(n * 2^n)

  // possible combinations for the set in [first, last] are given by permutations of the binary
  // number with v.size() bits with any number of bits activated
  // unsigned long long is 8 bytes long in most computers = 8*8 bits = 64 bits
  // this function works well as long as the maximum number of elements in the set is 64
  // (for more elements, a cleverer way can be thought of using more than one unsigned long long)

  unsigned int nElements = (unsigned int) input.size();
  unsigned int iElements = (unsigned int) param()[p_cut_n_sub];
  if (nElements < iElements)
    iElements = nElements;
  iElements = nElements;
  if (iElements > sizeof(unsigned long long)*8) {
    throw NEW_EXCEPTION("Number of microjets [maxElements] exceeds size of unsigned long long ([maxsize] bits)\
                         in this computer.")
                         .setParam("maxsize", d_to_string(sizeof(unsigned long long)*8))
                         .setParam("maxElements", d_to_string(iElements));
  }
  unsigned long long nSets = (unsigned long long) (0x1 << iElements); // 2^iElements

  // we generate the isrWeights and the signal weights first
  // then we sum the signal weight and check whether we should move on
  std::vector<double> isrWeight(nSets, 0);

  // the bits position in index i give the positions of the vector to be used
  for (unsigned long long i = 0; i < nSets; ++i) {
    // check the bits set in i and include the elements corresponding to
    // that bit position in result
    // creating a vector in result itself and then changing it saves copy time

    StoredJet isrJets(store);
    StoredJet hardJets(store);

    fastjet::PseudoJet isrSum(0,0,0,0);
    fastjet::PseudoJet hardSum(0,0,0,0);

    for (unsigned int k = 0; k < iElements; ++k) {
      if ( (i & (0x1 << k)) != 0) {
        isrJets.push_back(k);
        //copyJI(input[k], isrJets.back());
        isrSum += input[k];
      } else { // if a jet is not in the ISR set, it is in the hard proccess
        hardJets.push_back(k);
        //copyJI(input[k], hardJets.back());
        hardSum += input[k];
      }
    }

    // no sort needed so far ...
    //std::sort(isrJets.begin(), isrJets.end(), lessThanIndex);
    //std::sort(hardJets.begin(), hardJets.end(), lessThanIndex);
  
    // now go over this combination of isr/sb jets and calculate weights for ISR, signal and background
    // first for ISR:
    // fractionISR_of_hardInt_PT: value which determines the fraction of ISR_pT vs HardInt_pT
    if (square(isrSum.perp()) >= Qsquare*square(param()[p_fractionISR_of_hardInt_PT]) ) {
      continue; // demand that the ISR p_T is less than fractionISR_of_hardInt_PT times the virtuality
                // of the hard proccess
    }

#ifdef DEBUGCODE
    LOG(DEBUG) << "accepted ISR vs Signal branch" << endl;
    LOG(DEBUG) << "ISR_PT**2 : " << square(isrSum.perp()) << endl;
    LOG(DEBUG) << "Qsquare : " << Qsquare << endl;
    LOG(DEBUG) << "Qsquare*pow(fractionISR_of_hardInt_PT,2) : " << Qsquare*square(param()[p_fractionISR_of_hardInt_PT]) << endl;
    LOG(DEBUG) << square(isrSum.perp()) << " < " << Qsquare*square(param()[p_fractionISR_of_hardInt_PT]) << endl;

    LOG(DEBUG) << "ISRbranch" << endl; printoutput(isrJets);
    LOG(DEBUG) << "SBbranch" << endl; printoutput(hardJets);

    LOG(DEBUG) << "pT_ISR^2: " << square(isrSum.perp()) << "  <  " << Qsquare/4 << endl;
    LOG(DEBUG) << "input ISR: " << endl;
    printoutput(isrJets);
#endif

    double thisISRWeight = 1.0;
    if (param()[p_noISR] < 0.1) {
      m_isr.setQsquare(Qsquare);
      thisISRWeight = m_isr.weight(isrJets, hardSum)*m_isr.hamiltonian();
    }
    isrWeight[i] = thisISRWeight;

#ifdef DEBUGCODE
    LOG(DEBUG) << "signal hypothesis: " << endl;
    LOG(DEBUG) << "input hard interaction" << endl;
    printoutput(hardJets);
#endif

    m_signal.setQsquare(Qsquare);
    double thisSignalWeight = m_signal.weight(hardJets, hardSum);
    double thisSignalH = m_signal.hamiltonian(hardSum.perp());
    wSignal += thisSignalH*thisISRWeight*thisSignalWeight;


    m_signalWeight.insert(std::pair<double, std::vector<fastjet::PseudoJet> >(thisSignalH*thisISRWeight*thisSignalWeight, (std::vector<fastjet::PseudoJet>) (hardJets)));

#ifdef DEBUGCODE
    LOG(DEBUG) << "input Hard" << endl;
    printoutput(hardJets,DEBUG);
    LOG(DEBUG) << "isrweight: " << isrWeight[i] << endl;
    LOG(DEBUG) << "vH = topweight: " << thisSignalWeight << endl;
    LOG(DEBUG) << "topweight*isrweight: " << isrWeight[i]*thisSignalWeight << endl;
    LOG(DEBUG) << "Hamiltonian signal: " << thisSignalH << endl;
    LOG(DEBUG) << "Result of this history (topweight*Hamiltonina= h* vH): " << thisSignalH*thisSignalWeight << endl;
    LOG(DEBUG) << "Result of this history including ISR (topweight*Hamiltonina*isrweight): " << thisSignalH*thisSignalWeight*isrWeight[i] << endl;

    LOG(DEBUG) << "------------------------------------------------ " << endl;
    LOG(DEBUG) << endl;

    LOG(DEBUG) << endl;
#endif


  // don't bother calculating the background weight if the signal weight is negative
  // need to create two loops for this: is it worth it if the sample is signal-enriched?
  // commenting it out for now: it should be faster if the sample has mostly signal
  }

  if (!(wSignal > 0)) {
    wBackground = -1;
    return wSignal/wBackground;
  }

  // background weights
  for (unsigned long long i = 0; i < nSets; ++i) {

    // check the bits set in i and include the elements corresponding to
    // that bit position in result
    // creating a vector in result itself and then changing it saves copy time

    StoredJet hardJets(store);
    StoredJet isrJets(store);

    fastjet::PseudoJet isrSum(0,0,0,0);
    fastjet::PseudoJet hardSum(0,0,0,0);

    for (unsigned int k = 0; k < iElements; ++k) {
      if ( (i & (0x1 << k)) != 0) {
        isrJets.push_back(k);
        //copyJI(input[k], isrJets.back());
        isrSum += input[k];
      } else { // if a jet is not in the ISR set, it is in the hard proccess
        hardJets.push_back(k);
        //copyJI(input[k], hardJets.back());
        hardSum += input[k];
      }
    }

    // no sort needed so far ...
    //std::sort(isrJets.begin(), isrJets.end(), lessThanIndex);
    //std::sort(hardJets.begin(), hardJets.end(), lessThanIndex);
  
    // now go over this combination of isr/sb jets and calculate weights for ISR, signal and background
    // first for ISR:
    // fractionISR_of_hardInt_PT: value which determines the fraction of ISR_pT vs HardInt_pT
    if (square(isrSum.perp()) >= Qsquare*square(param()[p_fractionISR_of_hardInt_PT]) ) {
      continue; // demand that the ISR p_T is less than fractionISR_of_hardInt_PT times the virtuality
                // of the hard proccess
    }

#ifdef DEBUGCODE
    LOG(DEBUG) << "background hypothesis: " << endl;
    LOG(DEBUG) << "input hard interaction" << endl; printoutput(hardJets);
#endif


    double thisBackgroundWeight = m_background.weight(hardJets, hardSum);
    double thisBackgroundH = m_background.hamiltonian(hardSum.perp());
    wBackground += thisBackgroundH*isrWeight[i]*thisBackgroundWeight;

#ifdef DEBUGCODE
    LOG(DEBUG) << "input ISR" << endl; printoutput(isrJets);
    LOG(DEBUG) << "isrweight: " <<  isrWeight[i]  << endl;
    LOG(DEBUG) << "backgroundweight: " <<  thisBackgroundWeight << endl;
    LOG(DEBUG) << "Hamiltonian backg: " << thisBackgroundH << endl;
    LOG(DEBUG) << "background*Hamiltonina: " << thisBackgroundH*thisBackgroundWeight << endl;
    LOG(DEBUG) << "------------------------------------------------ " << endl;
    LOG(DEBUG) << endl;
    LOG(DEBUG) << "isrweight*background*Hamiltonian: " << thisBackgroundH*thisBackgroundWeight*isrWeight[i]  << endl;
#endif
  }

  double finalWeight = 0;
  if (wBackground > 0)
    finalWeight = wSignal/wBackground;

  return finalWeight;
}

const std::multimap<double, std::vector<fastjet::PseudoJet> > &Deconstruction::Deconstruct::signalWeight() const {
  return m_signalWeight;
}

