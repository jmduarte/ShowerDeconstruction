#include "ShowerDeconstruction/SDAlgorithm/interface/AnalysisParameters.h"

const Parameter listOfParameters[] = {
    CREATE_PARAMETER(R),
    CREATE_PARAMETER(ptmin_jet),
    CREATE_PARAMETER(n_of_jets),
    CREATE_PARAMETER(rapmax_jet),
    CREATE_PARAMETER(eta_cell),
    CREATE_PARAMETER(phi_cell),
    CREATE_PARAMETER(cell_ptcut),
    CREATE_PARAMETER(ptcut_microjets),
    CREATE_PARAMETER(masscut_microjets),
    CREATE_PARAMETER(lambda_mu_ext),
    CREATE_PARAMETER(lambda_theta_ext),
    CREATE_PARAMETER(cut_n_sub),
    CREATE_PARAMETER(fractionISR_of_hardInt_PT),
    CREATE_PARAMETER(noISR),
    CREATE_PARAMETER(topmass),
    CREATE_PARAMETER(topwindow),
    CREATE_PARAMETER(wmass),
    CREATE_PARAMETER(wwindow),
    CREATE_PARAMETER(Higgsmass),
    CREATE_PARAMETER(delta_Higgsmass),
    CREATE_PARAMETER(Rsmall),
    CREATE_PARAMETER(m_min_jet),
    CREATE_PARAMETER(br_wqq),
    CREATE_PARAMETER(useBtag),
    CREATE_PARAMETER(tagprob),
    CREATE_PARAMETER(fakeprob)
  };

AnalysisParameters::AnalysisParameters()
  : Deconstruction::Parameters() {
  init();
}

AnalysisParameters::AnalysisParameters(const std::string &input)
  : Deconstruction::Parameters() {
  init();
  read(input);
}

void AnalysisParameters::init() {
  for (int k = 0; k < sizeof(listOfParameters)/sizeof(Parameter); ++k) {
    insert(listOfParameters[k].m_key, 0);
  }
}

