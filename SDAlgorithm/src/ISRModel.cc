#include "ShowerDeconstruction/SDAlgorithm/interface/ISRModel.h"

#include <vector>
#include <fastjet/PseudoJet.hh>
#include "ShowerDeconstruction/SDAlgorithm/interface/Helper.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/Parameters.h"

#include <memory>
#include <algorithm>

#include "ShowerDeconstruction/SDAlgorithm/interface/partition_deconstruction.h"

using namespace Deconstruction;
using namespace std;
using namespace fastjet;

Deconstruction::ISRModel::ISRModel(Parameters &param)
  : Deconstruction::Model::Model(param) {
}

Deconstruction::ISRModel::~ISRModel() {
}

//double Deconstruction::ISRModel::weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sumJets) {
double Deconstruction::ISRModel::weight(const StoredJet &jets, fastjet::PseudoJet &sumJets) {

  double Qsquare = m_Qsquare;

  double isr_sudakov = 1.0; // simplified version -- cancels anyway between signal and background!
  // 0) calculate total Sudakov, which tells how likely it is
  //    that vector<PseudoJet*> input ended up in the cone.
  //    This Sudakov cancels in S/B...

  if (jets.size() == 0)
    return isr_sudakov;

  vector<int> v = jets.getList();

  int iisize((int)jets.size());

  vector< vector< vector<int> > > partonsetmix;
  partonsetmix.clear();
  try
    {
      partition::iterator it(iisize);
      while(true)
        {
          // cout << it << " : " << it.subsets() << " : ";
          auto_ptr<vector<vector<int> > > part = it[v];
          partonsetmix.push_back(*part);

          ++it;
        }
    }
  catch (overflow_error&)
    {
      //return(0.0);
    }


  // now we have all combinations for the ISR which enter the FSR shower.

  // man hat: H*FSR_shower fuer jede moegliche ISR-Abstrahlung-- 
  //          der Sudakov ist nach unserer neuen definition global 
  //          gleich und cancelt zwischen Signal und Untergrund.
  //          ISR_sudakov damit = 1;

  double ISRweight = 0.0; // muss summiert werden
  for(unsigned ii=0; ii<partonsetmix.size(); ii++)
    {
      if(partonsetmix[ii].size()==0) continue; //faengt addition von ISRtmp=1 ab
      double ISRtmp(1.0); // muss multipliziert werden
      for(unsigned jj=0; jj<partonsetmix[ii].size(); jj++)
        {
          // all weights for each vector in partonsetmix is summed up

          // weight for each ISR radiation:
          // H*Sudakov_tot*FSRshower
          StoredJet jets_ij(jets.store(), partonsetmix[ii][jj]);
          PseudoJet tmp_J = jets_ij.sum();

          // theta cut in ISR hamiltonian:
          // has to jump out of the loop with ISRtmp=0.0
          // because if in {{1} {2}} the cut is failed by {2} only, the whole
          // shower history is rejected....
          if(tmp_J.perp() <= tmp_J.m() )
            {
              ISRtmp=0.0;
              break;
            }


          // ISR/UE hamiltonian, non-pert part und pert part
          double hamiltonian(8.0*M_PI*Cte::CA*alphas(square(tmp_J.perp())+Cte::kappa_p2)
                             /(square(tmp_J.perp()) + Cte::kappa_p2)
                             /pow(1.0+Cte::cr*tmp_J.perp()/sqrt(Qsquare),Cte::nr)
                             + 16.0*M_PI*Cte::cnp*pow(Cte::kappa_np2,Cte::nnp-1)/pow(square(tmp_J.perp()) + Cte::kappa_np2,Cte::nnp));

          //ISRtmp *= hamiltonian * start_splitting(partonsetmix[ii][jj], Flavour::g, Shower::QCD);
          ISRtmp *= hamiltonian * start_splitting(jets_ij, Flavour::g, Shower::QCD);
        }

      ISRweight += ISRtmp;
    }

  return ISRweight*isr_sudakov;

}

double Deconstruction::ISRModel::hamiltonian(double pTsum) {
  return 1.0;
}

