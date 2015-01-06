#include "ShowerDeconstruction/SDAlgorithm/interface/Model.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/Helper.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/Parameters.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/AnalysisParameters.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/TopModel.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/Exception.h"

#include <fastjet/PseudoJet.hh>

#include <string>
#include <iostream>

#include <iomanip>

#include <cmath>

#include "ShowerDeconstruction/SDAlgorithm/interface/JetInfo.h"

#include "ShowerDeconstruction/SDAlgorithm/interface/StoredJet.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/StoredCalculations.h"

using namespace Deconstruction;

using std::endl;
using std::max;
using std::fabs;
using std::min;
using std::vector;
using namespace fastjet;

int Model::level = -1;
int Model::btopshower_level= -1;
int Model::tshower_level = -1;
int Model::Wshower_level = -1;

int Model::initlevel = -1;
int Model::binitlevel = -1;
int Model::tinitlevel = -1;
int Model::Winitlevel = -1;

StoredCalculations Model::m_calc = StoredCalculations();

Deconstruction::Model::Model(Parameters &p)
  : m_param(p) {

  nrbreitw = 2.0;

  wwidth = std::fabs(square(param()[p_wmass]+param()[p_wwindow]) - square(param()[p_wmass]))/nrbreitw/param()[p_wmass];
  delta_wmass = nrbreitw*wwidth;

  topwidth = std::fabs(square(param()[p_topmass]+param()[p_topwindow]) - square(param()[p_topmass]))/nrbreitw/param()[p_topmass];
  delta_topmass = nrbreitw*topwidth;

  lambda_mu = param()[p_lambda_mu_ext];
  topmass = param()[p_topmass];

  m_useBtag = param()[p_useBtag];
  m_tagprob = param()[p_tagprob];
  m_fakeprob = param()[p_fakeprob];

  dR2_outside = square(param()[p_R]);
}

void Deconstruction::Model::setQsquare(double q) {
  m_Qsquare = q;
}
 
Deconstruction::Model::~Model() {
}

Parameters &Deconstruction::Model::param() const {
  return m_param;
}

double Deconstruction::Model::start_splitting(const StoredJet &input, const Flavour::id flavour, const Shower::id shower) {
  m_calc.clear();

  StoredJet empty(input.store());
  return make_splitting(input, empty, empty, empty, flavour, Flavour::noflav, Flavour::noflav, Flavour::noflav, shower);
}

double Deconstruction::Model::make_splitting(const StoredJet &input, const StoredJet &leftcolpartner, const StoredJet &rightcolpartner, const StoredJet &grandmother, const Flavour::id flavor, const Flavour::id granflavor, const Flavour::id leftcolflav, const Flavour::id rightcolflav, const Shower::id shower) {

  // check if the calculation was already done
  std::vector<int> a = input.getList();
  std::sort(a.begin(), a.end());
  std::vector<int> b = leftcolpartner.getList();
  std::sort(b.begin(), b.end());
  std::vector<int> c = rightcolpartner.getList();
  std::sort(c.begin(), c.end());
  std::vector<int> d = grandmother.getList();
  std::sort(d.begin(), d.end());
  StoredKey sk;
  for (int i = 0; i < a.size(); ++i) sk.push_back(a[i]);
  sk.push_back(-1);
  for (int i = 0; i < b.size(); ++i) sk.push_back(b[i]);
  sk.push_back(-1);
  for (int i = 0; i < c.size(); ++i) sk.push_back(c[i]);
  sk.push_back(-1);
  for (int i = 0; i < d.size(); ++i) sk.push_back(d[i]);
  sk.push_back(-1);
  sk.push_back(flavor);
  sk.push_back(granflavor);
  sk.push_back(leftcolflav);
  sk.push_back(rightcolflav);
  sk.push_back(shower);
  double w = 0;
  if (m_calc.check(sk, w)) {
    return w;
  }

  ++Model::level;
  if (shower==Shower::W) ++Model::Wshower_level;
  if (shower==Shower::b) ++Model::btopshower_level;
  if (shower==Shower::t) ++Model::tshower_level;

#ifdef DEBUGCODE
  LOG(DEBUG) << "input make_splitting - level: " << Model::level<<  endl;
  if(shower==Shower::W) LOG(DEBUG) << "Wshower_level: " << Model::Wshower_level << endl;
  if(shower==Shower::b) LOG(DEBUG) << "Model::btopshower_level: " << Model::btopshower_level << endl;
  if(shower==Shower::t) LOG(DEBUG) << "tshower_level: " << Model::tshower_level << endl;

  printoutput(input,"input");
  printoutput(leftcolpartner,"leftcolpartner");
  printoutput(rightcolpartner,"rightcolpartner");
  printoutput(grandmother,"grandmother");
  LOG(DEBUG) << "flavor: " << flavor << endl;
  LOG(DEBUG) << "granflavor: " << granflavor << endl;
  LOG(DEBUG) << "leftcolflavor: " << leftcolflav << endl;
  LOG(DEBUG) << "rightcolflavor: " << rightcolflav << endl;
  LOG(DEBUG) << "shower: " << shower << endl;

#endif

  /// top has to decay into at least 3 particles /////////
  if ((flavor==Flavour::t || flavor==Flavour::tbar) && (input.size() < 3) ) {
    --Model::level;
    if (shower==Shower::W) --Model::Wshower_level;
    if (shower==Shower::b) --Model::btopshower_level;
    if (shower==Shower::t) --Model::tshower_level;
    m_calc.store(sk, 0);
    return 0;
  }

  dR2_outside = square(param()[p_R]);

  //////////////////////////////////////////////////////////////////////////////
  /////////// if last particle of one branch //////////////////////////////////
  if (input.size() < (unsigned) 2) {
    if (input.size() == 0) {
      --Model::level;
      if (shower==Shower::W) --Model::Wshower_level;
      if (shower==Shower::b) --Model::btopshower_level;
      if (shower==Shower::t) --Model::tshower_level;
      m_calc.store(sk, 0);
      return 0;
    }

    PseudoJet tmptot(input.sum());
    PseudoJet tmpgrandm(grandmother.sum());
    PseudoJet tmpmrightcol(rightcolpartner.sum());
    PseudoJet tmpmleftcol(leftcolpartner.sum());

    double weight_lastparton(0.0);

    if (flavor==Flavour::g)
      weight_lastparton = sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor, granflavor, leftcolflav,rightcolflav,shower);
    else if ((flavor == Flavour::q) || (flavor == Flavour::qbar) || (flavor == Flavour::b) || (flavor == Flavour::bbar))
      weight_lastparton = sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor,granflavor, leftcolflav, rightcolflav,shower);
    else
      LOG(ERROR) << "error in sudakov make_splitting: wrong flavor" << endl;

    // for radiation off b-quark after top decay we need to know 
    // the flavor of the color connected partner (here top)

    LOG(DEBUG) << "sudakov value last parton: " << weight_lastparton << endl;

    double btaggerweight(1.0);
    if (m_useBtag) {
      double isNotB = 0;
      if ( (flavor != Flavour::b) && (flavor != Flavour::bbar)) isNotB = 1.0;
      if (tmptot.user_index() == 0)   btaggerweight = 1;     // not tried to be b-tagged
      else if (tmptot.user_index() < 0) btaggerweight = 1 - (m_tagprob*(1-isNotB) + m_fakeprob*isNotB );
      else if (tmptot.user_index() > 0) btaggerweight = m_tagprob*(1-isNotB) + m_fakeprob*isNotB;
    }
    // construct here sudakov if input is last parton (microjet)
	  
    LOG(DEBUG) << "exp(-sudakov)*btag last parton: " << weight_lastparton*btaggerweight << endl;


    --Model::level;
    if(shower==Shower::W) --Model::Wshower_level;
    if(shower==Shower::b) --Model::btopshower_level;
    if(shower==Shower::t) --Model::tshower_level;
 
    m_calc.store(sk, weight_lastparton*btaggerweight);
    return(weight_lastparton*btaggerweight); // this has to be changed to the sudakov
  }
  
  // return value
  double valuetotal(0.0);
  
  /////////////////// the top quark is special: it decays always before hadronization 
  /////////////////// thus we do not allow the top to radiate down to the hadronization scale
  /////////////////// in every step it has to be checked if the top decays and the weight has
  /////////////////// to be added to the t -> t g splitting weight.

  double topdecayweight(0.0);
  double antitopdecayweight(0.0);

  if(flavor == Flavour::t && (input.size()>2) ) {
    PseudoJet tmptot(input.sum());
      
    if (std::fabs(square(tmptot.m())-square(topmass))<topmass*delta_topmass) {
      PseudoJet tmpgrandm(grandmother.sum());
      PseudoJet tmpmrightcol(rightcolpartner.sum());
      PseudoJet tmpmleftcol(leftcolpartner.sum());
	  
      topdecayweight = sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,Flavour::t, granflavor,leftcolflav,rightcolflav,shower);                             // multiply sudakov for top before decay
      double sudtop=topdecayweight;
      double topdecaymodel = 0;
      TopModel tm(param(), flavor);
      topdecaymodel=tm.weight(input, tmptot);
      topdecayweight *= topdecaymodel;

#ifdef DEBUGCODE
      LOG(DEBUG) << "input in top decay model: " << endl;
      printoutput(input);
      LOG(DEBUG) << "results: " << endl;
      LOG(DEBUG) << "exp(-St) top vor decay: " << sudtop << endl;
      LOG(DEBUG) << "Htopdecaymodel: " << topdecaymodel << endl;
      LOG(DEBUG) << "Htopdecaymodel*exp(-St): " << topdecayweight << endl;
#endif

      valuetotal += topdecayweight;
    }
  }


  if(flavor == Flavour::tbar && (input.size()>2)) {
    PseudoJet tmptot(input.sum());
      
    if (std::fabs(square(tmptot.m())-square(topmass))<topmass*delta_topmass) {
      PseudoJet tmpgrandm(grandmother.sum());
      PseudoJet tmpmrightcol(rightcolpartner.sum());
      PseudoJet tmpmleftcol(leftcolpartner.sum());
	  
      antitopdecayweight = sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,Flavour::tbar, granflavor,leftcolflav,rightcolflav,shower);

      // multiply sudakov for antitop before decay 
      
      TopModel tm(param(), flavor);
      antitopdecayweight *= tm.weight(input, tmptot);
      valuetotal += antitopdecayweight;
    }
  }
  ///////////// until here decays /////////////////////////////

  StoredJet empty(input.store());

  /////////// now start splittings ////////////////////////////
  /////// split into left and right branches /////////////////////////

  // possible combinations for the set in [first, last] are given by permutations of the binary
  // number with input.size() bits with any number of bits activated
  unsigned int iElements = (unsigned int) input.size();
  // iElements cannot be greater than the maximum,
  // but this has already been tested before, when generating the powerset for ISR/SB
  // there is no point in testing it again, since these sets must be smaller or equal in size
  // than the subsets of the powerset of all microjets
  unsigned long long nSets = (unsigned long long) (0x1 << iElements); // 2^iElements

  // the bits position in index i give the positions of the vector to be used
  // we are skipping the empty set combination here
  // that's why i starts at one and ends at nSets-1
  for (unsigned long long i = 1; i < nSets; ++i) {
    //std::vector<fastjet::PseudoJet> leftbranch;
    //std::vector<fastjet::PseudoJet> rightbranch;

    StoredJet leftbranch(input.store());
    StoredJet rightbranch(input.store());

    // check the bits set in i and include the elements corresponding to
    // that bit position in result
    for (unsigned int k = 0; k < iElements; ++k) {
      if ( (i & (0x1 << k)) != 0) {
        leftbranch.push_back(input.getIdx(k));
      } else {
        rightbranch.push_back(input.getIdx(k));
      }
    }

    PseudoJet tmptot(input.sum());
    PseudoJet tmpL(leftbranch.sum());
    PseudoJet tmpR(rightbranch.sum());
    PseudoJet tmpgrandm(grandmother.sum());
    PseudoJet tmpmrightcol(rightcolpartner.sum());
    PseudoJet tmpmleftcol(leftcolpartner.sum());
      
#ifdef DEBUGCODE
    LOG(DEBUG) << "split left right: " << endl;
    LOG(DEBUG) << "left: " << endl;
    printoutput(leftbranch);
    LOG(DEBUG) << "right: " << endl;
    printoutput(rightbranch);
#endif

    // next make_splittings have to be called with input (mother) as argument for grandmother...
    double valueleft(0.0);
    double valueright(0.0);
    double valuegbbleft(0.0);
    double valuegbbright(0.0);
    double valuegqqleft(0.0);
    double valuegqqright(0.0);
    
    double splitfactorgqq(0.0),splitfactorggg(0.0),splitfactorgbb(0.0);
    double splitfactorbbg(0.0),splitfactorqqg(0.0),splitfactorantibbg(0.0),splitfactorantiqqg(0.0);
    double splitfactorttg(0.0),splitfactorantittg(0.0);
    
    double checksudtot(0.0);
    double checksplittot(0.0);

    double gbb(0.0),ggg(0.0),gqq(0.0),qqg(0.0),antiqqg(0.0),bbg(0.0),antibbg(0.0);
    double gtotal(0.0), ttg(0.0),antittg(0.0);
     
    // IMPORTANT
    // if Flavour::g splits to quarks, the anti-quark goes left and the quark right
    // if Higgs splits to quarks, the anti-quark goes right and the quark left
    
    if(flavor == Flavour::g) // Flavour::g splits
    {
      //msp: ATTENTION! for massive b I might need to include new check for splitting (including massive b to calculate virtuality)

      //check if splitting allowed
	  
      if( ((tmpL.perp() != 0) && ((square(tmpL.m())/tmpL.perp()) >= lambda_mu*(square(tmptot.m())/tmptot.perp()*0.5))) ||
          ((tmpR.perp() != 0) && ((square(tmpR.m())/tmpR.perp()) >= lambda_mu*(square(tmptot.m())/tmptot.perp()*0.5))) )
        continue;

      // calc sudakov for this "propagator line"
      double sudtot(sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor, granflavor,leftcolflav,rightcolflav,shower));

      /// g-> g g splitting
#ifdef DEBUGCODE
      LOG(DEBUG) << "do g->g g splitting " << endl;
#endif

      splitfactorggg = HGluonGG(tmptot,tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor,granflavor,leftcolflav,rightcolflav,shower);

#ifdef DEBUGCODE
      LOG(DEBUG) << "for g->gg splitting: " << endl;
      LOG(DEBUG) << "Hgqq: " << splitfactorggg << endl;
      LOG(DEBUG) << "sudakov: " << sudtot << endl;
      LOG(DEBUG) << "Hggg*sudakov: " << sudtot*splitfactorggg << endl;
#endif

      valueleft = make_splitting(leftbranch,leftcolpartner, rightbranch,input,Flavour::g,flavor,leftcolflav,Flavour::g,shower);
      valueright = make_splitting(rightbranch,leftbranch, rightcolpartner,input,Flavour::g,flavor,Flavour::g,rightcolflav,shower);
	  
      ggg = valueleft*valueright*splitfactorggg;

      // g->q q splitting
#ifdef DEBUGCODE
      LOG(DEBUG) << "do g->q q splitting " << endl;	       
#endif

      // TO speed up the code HGluonQQ can be replaced by HGluonBB and simply multiplied by nf-1
      // HGluonQQ=HGluonBB*(nf-1)
      splitfactorgqq = (Cte::nf-1)*HGluonQQ(tmptot, tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor, granflavor, leftcolflav , rightcolflav, shower);

#ifdef DEBUGCODE	 
      LOG(DEBUG) << "for g->qq splitting: " << endl;
      LOG(DEBUG) << "Hgqq: " << splitfactorgqq/Cte::nf << endl;
      LOG(DEBUG) << "nf*Hgqq: " << splitfactorgqq << endl;
      LOG(DEBUG) << "sudakov: " << sudtot << endl;
      LOG(DEBUG) << "nf*Hgqq*sudakov: " << sudtot*splitfactorgqq << endl;
#endif

      valuegqqleft = make_splitting(leftbranch,leftcolpartner, empty,input,Flavour::qbar,flavor,leftcolflav,Flavour::noflav,shower);
      valuegqqright = make_splitting(rightbranch, empty, rightcolpartner,input,Flavour::q,flavor,Flavour::noflav,rightcolflav,shower);
	  
      gqq = valuegqqleft*valuegqqright*splitfactorgqq;
	
      // g->b b splitting
      LOG(DEBUG) << "do g->b b splitting " << endl;
      splitfactorgbb = HGluonBB(tmptot,tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor, granflavor, leftcolflav , rightcolflav, shower); 
      LOG(DEBUG) << "H g->bb : " << splitfactorgbb << endl;
      LOG(DEBUG) << "H exp(-S) : " << sudtot*splitfactorgbb << endl;

      valuegbbleft = make_splitting(leftbranch, leftcolpartner, empty, input,Flavour::bbar,flavor,leftcolflav,Flavour::noflav,shower);
      valuegbbright = make_splitting(rightbranch, empty, rightcolpartner,input,Flavour::b,flavor,Flavour::noflav,rightcolflav,shower);
	      	      
      gbb = valuegbbleft*valuegbbright*splitfactorgbb;

      gtotal = (gbb+gqq+ggg)*sudtot;

#ifdef DEBUGCODE
      LOG(DEBUG) << "flavor: " << flavor << endl;
      printoutput(input, "total branch");
      printoutput(leftbranch,"left branch");
      printoutput(rightbranch,"right branch");
      LOG(DEBUG) << "in Flavour::gsplitting: " << endl;
      LOG(DEBUG) << "ggg: valueleft - " << valueleft << endl;
      LOG(DEBUG) << "ggg: valueright - " << valueright << endl;
      LOG(DEBUG) << "sudakov factor for this propagator line: " << sudtot << endl;
      LOG(DEBUG) << "ggg: Hggg - " << splitfactorggg << endl;
      LOG(DEBUG) << "exp(-S) : " << sudtot << endl;
      LOG(DEBUG) << "Hggg*exp(-s): " << splitfactorggg*sudtot << endl;
      LOG(DEBUG) << "valueleft*valueright*H " << ggg << endl;

      LOG(DEBUG) << "gqq: valuegqqleft - " << valuegqqleft << endl;
      LOG(DEBUG) << "gqq: valuegqqright - " << valuegqqright << endl;
      LOG(DEBUG) << "gqq: Hgqq - " << splitfactorgqq << endl;
      LOG(DEBUG) << "exp(-S) : " << sudtot << endl;
      LOG(DEBUG) << "Hgqq*exp(-S) - " << splitfactorgqq*sudtot << endl;
      LOG(DEBUG) << "valueleft*valueright*H " << gqq << endl;

      LOG(DEBUG) << "Hgbb: " << splitfactorgbb << endl;
      LOG(DEBUG) << "Hgqq: " << splitfactorgqq << endl;
      LOG(DEBUG) << "Hggg: " << splitfactorggg << endl;
      LOG(DEBUG) << "Hggg*exp(-S): " << splitfactorggg*sudtot << endl;
      LOG(DEBUG) << "Hgbb*exp(-S): " << splitfactorgbb*sudtot << endl;
      LOG(DEBUG) << "Hgqq*exp(-S): " << splitfactorgqq*sudtot << endl;
      LOG(DEBUG) << "(Hggg+nf*Hgqq)*exp(-S): " << (splitfactorgqq+splitfactorgbb+splitfactorggg)*sudtot << endl;
      LOG(DEBUG) << "gbb : " << gbb << endl;
      LOG(DEBUG) << "gqq : " << gqq << endl;
      LOG(DEBUG) << "ggg : " << ggg << endl;
      LOG(DEBUG) << "sudakov: " << sudtot << endl;
      LOG(DEBUG) << "gbb*exp(-S) : " << gbb*sudtot << endl;
      LOG(DEBUG) << "gqq*exp(-S) : " << gqq*sudtot << endl;
      LOG(DEBUG) << "ggg*exp(-S) : " << ggg*sudtot << endl;
      LOG(DEBUG) << "(gqq+ggg)*exp(-S): " << (gqq+ggg)*sudtot << endl;
      LOG(DEBUG) << "for gtotal: " << gtotal << endl;
#endif

    } else if(flavor == Flavour::t) {
      // top can only radiate if for the decay 3 subjets remain
      if (input.size() < 3)
        continue;

      // definition of virtuality $\mu_J^2 = (p_A + p_B)^2 - m_J^2$
      double muJ2(square(tmptot.m())-square(topmass));
      double muh2(std::fabs(square(tmpL.m())-square(topmass))); //top

      //msp changed in last version check!! if(muh2 < topmass*topwidth) muh2 = topmass*topwidth;
      // muh2 = topmass*topwidth;
      double mus2(square(tmpR.m())); //Flavour::g
	  
      //check splitting with massive daughters
      //check if splitting allowed

      if ( ((tmpL.perp() != 0) && ((std::fabs(muh2)/tmpL.perp()) >= lambda_mu*(muJ2/tmptot.perp()*0.5))) ||
           ((tmpR.perp() != 0) && ((std::fabs(mus2)/tmpR.perp()) >= lambda_mu*(muJ2/tmptot.perp()*0.5))) ) 
	 continue;

      // new cut:
      /*
      if(granflavor == noflav) {
        if((2*lambda_mu*(muJ2/tmptot.perp())) > (2*(muJ2+square(topmass))/tmptot.perp())) continue;
      } else if(granflavor == Flavour::t) {
        if((2*lambda_mu*(muJ2/tmptot.perp())) > (square(tmpgrandm.m())-square(topmass))/tmpgrandm.perp()) continue;
      } else
        LOG(DEBUG) << "ERROR grandmother no Flavour::t flavor although jet is Flavour::t" << endl;
      */
      // hier muss man immer beide gewichte addieren: splitting t -> t + g und top decay

      double sudtot(sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor,granflavor, leftcolflav, rightcolflav,shower));

      valueleft = make_splitting(leftbranch, leftcolpartner, rightbranch,input,Flavour::t,flavor,Flavour::noflav,Flavour::g,shower);
      valueright = make_splitting(rightbranch,leftbranch, rightcolpartner,input,Flavour::g,flavor,Flavour::t,rightcolflav,shower);

      splitfactorttg = HQuark(tmptot,tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor,granflavor,leftcolflav,rightcolflav,shower);

      ttg = valueleft*valueright*splitfactorttg*sudtot;
    } else if(flavor == Flavour::tbar) {
      // antitop can only radiate if for the decay 3 subjets remain 
      //if (input.size() < 3)
      if (input.size() < 3)
        continue;

      // definition of virtuality $\mu_J^2 = (p_A + p_B)^2 - m_J^2$
      double muJ2(square(tmptot.m())-square(topmass));
      double muh2(std::fabs(square(tmpR.m())-square(topmass))); //top
      // if(muh2 < topmass*topwidth) muh2 = topmass*topwidth;
      // muh2 = topmass*topwidth;
      double mus2(square(tmpL.m())); //Flavour::g

      //check splitting with massive daughters
      //check if splitting allowed
      if (((tmpR.perp() != 0) && ((std::fabs(muh2)/tmpR.perp()) >= lambda_mu*(muJ2/tmptot.perp()*0.5))) ||
          ((tmpL.perp() != 0) && ((std::fabs(mus2)/tmpL.perp()) >= lambda_mu*(muJ2/tmptot.perp()*0.5))) )
        continue;

      double sudtot(sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor,granflavor, leftcolflav, rightcolflav,shower));

      valueleft = make_splitting(leftbranch, leftcolpartner, rightbranch,input,Flavour::g,flavor,leftcolflav,Flavour::tbar,shower);
      valueright = make_splitting(rightbranch,leftbranch, empty,input,Flavour::tbar,flavor,Flavour::g,Flavour::noflav,shower);
	  
      splitfactorantittg = HQuark(tmptot,tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor,granflavor,leftcolflav,rightcolflav,shower);

      antittg = valueleft*valueright*splitfactorantittg*sudtot;
    } else if(flavor == Flavour::bbar) // antib splits always antib (right) Flavour::g (left)
    {
      // definition of virtuality $\mu_J^2 = (p_A + p_B)^2 - m_J^2$
      // however take b massless 
      double muJ2(square(tmptot.m()));
      double muh2(square(tmpR.m())); //antib
      double mus2(square(tmpL.m())); //Flavour::g
	   
      if (((tmpR.perp() != 0) && ((std::fabs(muh2)/tmpR.perp()) >= lambda_mu*(muJ2/tmptot.perp()*0.5))) ||
          ((tmpL.perp() != 0) && ((std::fabs(mus2)/tmpL.perp()) >= lambda_mu*(muJ2/tmptot.perp()*0.5))) )
	      continue;

      double sudtot(sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor,granflavor, leftcolflav, rightcolflav,shower));

      // left splitting inherits leftcolpartner from mother
      // and the right branch is the right color partner
      valueleft = make_splitting(leftbranch, leftcolpartner, rightbranch,input,Flavour::g,flavor,leftcolflav,Flavour::bbar,shower);
      // right splitting inherits rightcolpartner from mother
      // and the left branch is the left color partner
      valueright = make_splitting(rightbranch, leftbranch, empty,input,Flavour::bbar,flavor,Flavour::g,Flavour::noflav,shower);

      splitfactorantibbg = HQuark(tmptot,tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor,granflavor,leftcolflav,rightcolflav,shower);

      antibbg = valueleft*valueright*splitfactorantibbg*sudtot;
    } else if(flavor == Flavour::b) // b splits always b (right) Flavour::g (left)
    {
      // definition of virtuality $\mu_J^2 = (p_A + p_B)^2 - m_J^2$
      double muJ2(square(tmptot.m()));
      double muh2(square(tmpL.m())); //b
      double mus2(square(tmpR.m())); //Flavour::g
	  
      //check splitting with massive daughters
      //check if splitting allowed
      if (((tmpL.perp() != 0) && ((std::fabs(muh2)/tmpL.perp()) >= lambda_mu*(muJ2/tmptot.perp()*0.5))) ||
          ((tmpR.perp() != 0) && ((std::fabs(mus2)/tmpR.perp()) >= lambda_mu*(muJ2/tmptot.perp()*0.5))) )
        continue;

      double sudtot(sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor,granflavor, leftcolflav, rightcolflav,shower));
	  
      valueleft = make_splitting(leftbranch, empty, rightbranch,input,Flavour::b,flavor,Flavour::noflav,Flavour::g,shower);
      valueright = make_splitting(rightbranch,leftbranch, rightcolpartner,input,Flavour::g,flavor,Flavour::b,rightcolflav,shower);

      splitfactorbbg = HQuark(tmptot,tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor,granflavor,leftcolflav,rightcolflav,shower);

      bbg = valueleft*valueright*splitfactorbbg*sudtot;
		
#ifdef DEBUGCODE
      LOG(DEBUG) << " flavor: " << flavor << endl;
      LOG(DEBUG) << " sudakov: " << sudtot << endl;
      LOG(DEBUG) << " leftbranch: " << endl;
      printoutput(leftbranch);
      LOG(DEBUG) << " rightbranch: " << endl;
      printoutput(rightbranch);
      LOG(DEBUG) << "splitfactorbbg: " << splitfactorbbg << endl;
      LOG(DEBUG) << "splitfactorbbg*sudtot: " << splitfactorbbg*sudtot << endl;
      LOG(DEBUG) << "bbg: " << bbg << endl;
#endif
    }  else if(flavor == Flavour::qbar) // antiq splits antiq (right) Flavour::g(left)
    {
      if (((tmpL.perp() != 0) && ((square(tmpL.m())/tmpL.perp()) >= lambda_mu*(square(tmptot.m())/tmptot.perp()*0.5))) ||
          ((tmpR.perp() != 0) && ((square(tmpR.m())/tmpR.perp()) >= lambda_mu*(square(tmptot.m())/tmptot.perp()*0.5))) )
        continue;

      double sudtot(sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor,granflavor, leftcolflav, rightcolflav,shower));

      splitfactorantiqqg = HQuark(tmptot,tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor,granflavor,leftcolflav,rightcolflav,shower);

      LOG(DEBUG) << "H_qbar->qbarg : " << splitfactorantiqqg << endl;
      LOG(DEBUG) << "exp(-S) : " << sudtot << endl;
      LOG(DEBUG) << "H exp(-S) : " << sudtot*splitfactorantiqqg << endl;

      // left splitting inherits leftcolpartner from mother
      // and the right branch is the right color partner
      valueleft = make_splitting(leftbranch, leftcolpartner, rightbranch,input,Flavour::g,flavor,leftcolflav,Flavour::qbar,shower);
      // right splitting inherits rightcolpartner from mother
      // and the left branch is the left color partner
      valueright = make_splitting(rightbranch, leftbranch, empty,input,Flavour::qbar,flavor,Flavour::g,Flavour::noflav,shower);

      antiqqg=valueleft*valueright*splitfactorantiqqg*sudtot;
    } else if(flavor == Flavour::q) // b splits (which is left branch)
    {

      if ( ((tmpL.perp() != 0) && ((square(tmpL.m())/tmpL.perp()) >= lambda_mu*(square(tmptot.m())/tmptot.perp()*0.5))) ||
           ((tmpR.perp() != 0) && ((square(tmpR.m())/tmpR.perp()) >= lambda_mu*(square(tmptot.m())/tmptot.perp()*0.5))) )
        continue;

      double sudtot(sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor,granflavor, leftcolflav, rightcolflav,shower));
      splitfactorqqg = HQuark(tmptot,tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor,granflavor,leftcolflav,rightcolflav,shower);

#ifdef DEBUGCODE
      LOG(DEBUG) << "H_qbar->qbarg : " << splitfactorqqg << endl;
      LOG(DEBUG) << "exp(-S) : " << sudtot << endl;
      LOG(DEBUG) << "H exp(-S) : " << sudtot*splitfactorqqg << endl;
#endif

      valueleft=make_splitting(leftbranch, empty, rightbranch,input,Flavour::q,flavor,Flavour::noflav,Flavour::g,shower);
      valueright=make_splitting(rightbranch,leftbranch, rightcolpartner,input,Flavour::g,flavor,Flavour::q,rightcolflav,shower);

      qqg=valueleft*valueright*splitfactorqqg*sudtot;

#ifdef DEBUGCODE
      LOG(DEBUG) << " flavor: " << flavor << endl;
      LOG(DEBUG) << " sudakov: " << sudtot << endl;
      LOG(DEBUG) << " leftbranch: " << endl;
      printoutput(leftbranch);
      LOG(DEBUG) << " rightbranch: " << endl;
      printoutput(rightbranch);
      LOG(DEBUG) << "splitfactorqqg: " << splitfactorqqg << endl;
      LOG(DEBUG) << "splitfactorqqg*sudtot: " << splitfactorqqg*sudtot << endl;
      LOG(DEBUG) << "qqg: " << qqg << endl;
#endif

    } else
      LOG(ERROR) << "in make_splittin: WRONG FLAVOR: " << flavor << " -- serious error" << endl;

    valuetotal += gtotal+bbg+antibbg+qqg+antiqqg+ttg+antittg;

#ifdef DEBUGCODE      
    LOG(DEBUG) << "in make_splitting value total: " << valuetotal << endl;
    LOG(DEBUG) << "flavor: " << flavor << endl;
    LOG(DEBUG) << "gtotal: " << gtotal << endl;
    LOG(DEBUG) << "bbg: " << bbg << endl;
    LOG(DEBUG) << "antibbg: " << antibbg << endl;
    LOG(DEBUG) << "qqg: " << qqg << endl;
    LOG(DEBUG) << "antiqqg: " << antiqqg << endl;
    LOG(DEBUG) << "ttg: " << ttg << endl;
    LOG(DEBUG) << "antittg: " << antittg << endl;
    LOG(DEBUG) << "topdecayweight: " << topdecayweight << endl;
    LOG(DEBUG) << "antitopdecayweight: " << antitopdecayweight << endl;
      
    //if((flavor == Flavour::t || flavor == Flavour::tbar) && (level == 0)) 
    LOG(DEBUG) << "level: " << level << endl;
    LOG(DEBUG) << "splitting combination nr. " << i << " of " << nSets<< endl;
    LOG(DEBUG) << "flavor: " << flavor << endl;
    LOG(DEBUG) << "left branch: " << endl;
    printoutput(leftbranch);
    LOG(DEBUG) << "right branch: " << endl;
    printoutput(rightbranch);
    LOG(DEBUG) << "checksudtot: " << checksudtot << endl;
    LOG(DEBUG) << "checksplittot: " << checksplittot << endl;
    LOG(DEBUG) << "topdecayweight: " << topdecayweight << endl;
    LOG(DEBUG) << "antitopdecayweight: " << antitopdecayweight << endl;
    LOG(DEBUG) << "ttg: " << ttg << endl;

    LOG(DEBUG) << "antittg: " << antittg << endl;
    LOG(DEBUG) << "gtotal: " << gtotal << endl;
    LOG(DEBUG) << "bbg: " << bbg << endl;
    LOG(DEBUG) << "antibbg: " << antibbg << endl;
    LOG(DEBUG) << "qqg: " << qqg << endl;
    LOG(DEBUG) << "antiqqg: " << antiqqg << endl;
    LOG(DEBUG) << "valuetotal: " << valuetotal << endl;

#endif

  } // end loop over left branches
  
  // if((flavor == Flavour::t || flavor == Flavour::tbar) && (level == 0))
#ifdef DEBUGCODE
  LOG(DEBUG) << "level: " << level << endl;
  printoutput(input);
  LOG(DEBUG) << "tshower: " << Shower::t << endl;
  LOG(DEBUG) << "topdecay: " << topdecayweight << endl;
  LOG(DEBUG) << "radiation off top: " << valuetotal - topdecayweight << endl;
  LOG(DEBUG) << "topdecay + radiation: " << valuetotal << endl << endl;
  LOG(DEBUG) << "____________" << endl;
  LOG(DEBUG) << "____________" << endl;
#endif

  --Model::level;
  if(shower==Shower::W) --Model::Wshower_level;
  if(shower==Shower::b) --Model::btopshower_level;
  if(shower==Shower::t) --Model::tshower_level;

  m_calc.store(sk, valuetotal);

  return valuetotal;
}

double Deconstruction::Model::HQuark(const fastjet::PseudoJet &tmptot,
                              const fastjet::PseudoJet & tmpleftcol,
                              const fastjet::PseudoJet & tmprightcol,
                              const fastjet::PseudoJet & tmpL,
                              const fastjet::PseudoJet & tmpR,
                              const fastjet::PseudoJet & tmpgrandm,
                              const Flavour::id flavour,
                              const Flavour::id granflavour,
                              const Flavour::id tmpleftcolflav,
                              const Flavour::id tmprightcolflav,
                              const Shower::id shower) {
  double Hamiltonian = 0.0;

  double anglefac = 1.0;  // should be one here because if there is no col. connected partner we take it 1;

  //note for fermions, left=hard, right=soft
  //     antifemrions, left=soft, right=hard

  double gfunc = 0.0;

  double mu2 = square(tmptot.m());

  if (flavour == Flavour::t) {
    if (shower != Shower::t)
      throw NEW_EXCEPTION("flavor = Flavour::t but not tshower, cant be...");
    if (tmprightcol.e() > Cte::smalldouble) {
      anglefac = ((tmpL.squared_distance(tmpR) + square(param()[p_topmass])/square(tmpL.perp()))*(tmpL.squared_distance(tmprightcol) + square(param()[p_topmass])/square(tmpL.perp()))
                  - square(param()[p_topmass])/square(tmpL.perp())*tmpR.squared_distance(tmprightcol));
      anglefac *= 1.0/((tmpL.squared_distance(tmpR) + square(param()[p_topmass])/square(tmpL.perp()))*(tmpR.squared_distance(tmprightcol) + tmpR.squared_distance(tmpL) + square(param()[p_topmass])/square(tmpL.perp())));
    } else {
      anglefac = tmpL.squared_distance(tmpR) / (tmpL.squared_distance(tmpR)+(square(param()[p_topmass])/square(tmpL.perp())));
    }

    mu2 -= square(param()[p_topmass]); // top virtuality is p^2-mt^2

    if (mu2 < 0.0)
      throw NEW_EXCEPTION("ERROR: MSP top in HQuark, mu2 cannot be smaller 0 ");

    Hamiltonian = 8.0*M_PI*Cte::CF*alphas(mu2)/mu2*tmptot.perp()/tmpR.perp()*(1.0+square(tmpL.perp()/tmptot.perp()))*anglefac;
  } else if (flavour == Flavour::tbar) {
    if (shower != Shower::t)
      throw NEW_EXCEPTION("flavour = Flavour::tbar but not tshower, cant be...");

    if (tmpleftcol.e() > Cte::smalldouble) {
      anglefac = ((tmpR.squared_distance(tmpL) + square(param()[p_topmass])/square(tmpR.perp()))*(tmpR.squared_distance(tmpleftcol) + square(param()[p_topmass])/square(tmpR.perp()))
                 - square(param()[p_topmass])/square(tmpR.perp())*tmpL.squared_distance(tmpleftcol));
      anglefac *= 1.0/((tmpR.squared_distance(tmpL) + square(param()[p_topmass])/square(tmpR.perp()))*(tmpL.squared_distance(tmpleftcol) + tmpL.squared_distance(tmpR) + square(param()[p_topmass])/square(tmpR.perp())));
    } else {
      anglefac = tmpL.squared_distance(tmpR) / (tmpL.squared_distance(tmpR)+(square(param()[p_topmass])/square(tmpR.perp())));
    }

    mu2 -= square(param()[p_topmass]); // top virtuality is p^2-mt^2
    if (mu2 < 0.0)
      throw NEW_EXCEPTION("ERROR: MSP antitop in HQuark, mu2 cannot be smaller 0 ");
    Hamiltonian = 8.0*M_PI*Cte::CF*alphas(mu2)/mu2*tmptot.perp()/tmpL.perp()*(1.0+square(tmpR.perp()/tmptot.perp()))*anglefac;

  } else if (flavour == Flavour::q || flavour == Flavour::b) //assume massless Flavour::bs, thus Flavour::q = Flavour::b for kinematic
  {

    if ((tmprightcolflav == Flavour::t) || (tmprightcolflav == Flavour::tbar)) {
      if (tmprightcol.e() < Cte::smalldouble)
        throw NEW_EXCEPTION("HQuark: color conn. right flavor is top but no energy?? cant be ... ");

      double numerator = (tmprightcol.squared_distance(tmpR)+square(param()[p_topmass])/square(tmprightcol.perp()))*
                         (tmprightcol.squared_distance(tmpL)+square(param()[p_topmass])/square(tmprightcol.perp()))
                         - (square(param()[p_topmass])/square(tmprightcol.perp()))*tmpR.squared_distance(tmpL);
      double denom = 0.0;
      if (shower == Shower::t) {
        denom = (tmpR.squared_distance(tmprightcol) + square(param()[p_topmass])/square(tmprightcol.perp()))*(tmpR.squared_distance(tmpL)+tmpR.squared_distance(tmprightcol) + square(param()[p_topmass])/square(tmprightcol.perp()));
      } else if (shower == Shower::b) {
        denom = square(tmprightcol.squared_distance(tmpR)+square(param()[p_topmass])/square(tmprightcol.perp()));
      } else
        throw NEW_EXCEPTION("ERROR in HQuark for quark: if top color connected has to be tshower or btopshower!! ");

      anglefac = numerator/denom;
    } else {
      if (tmprightcolflav == Flavour::noflav) {
        anglefac = 1.0;
      } else {
        anglefac = tmprightcol.squared_distance(tmpL)/(tmpR.squared_distance(tmpL) + tmpR.squared_distance(tmprightcol));
      }
    }

    Hamiltonian = 8.0*M_PI*Cte::CF*alphas(mu2)/mu2*tmptot.perp()/tmpR.perp()*(1.0+square(tmpL.perp()/tmptot.perp()))*anglefac;
  } else if (flavour == Flavour::qbar || flavour == Flavour::bbar) {
    //note for fermions, left=hard, right=soft
    //     antifemrions, left=soft, right=hard

    if ((tmpleftcolflav == Flavour::t) || (tmpleftcolflav == Flavour::tbar)) {
      if (tmpleftcol.e() < Cte::smalldouble)
        throw NEW_EXCEPTION("ERROR in HQuark: color conn. left flavor is top but no energy?? cant be ... ");

      double numerator = (tmpleftcol.squared_distance(tmpL)+square(param()[p_topmass])/square(tmpleftcol.perp()))*
                         (tmpleftcol.squared_distance(tmpR)+square(param()[p_topmass])/square(tmpleftcol.perp()))
                       - ((square(param()[p_topmass])/square(tmpleftcol.perp()))*tmpL.squared_distance(tmpR));

      double denom = 0;

      if (shower == Shower::t) {
        denom = ((tmpL.squared_distance(tmpleftcol) + square(param()[p_topmass])/square(tmpleftcol.perp()))*(tmpR.squared_distance(tmpL)+tmpL.squared_distance(tmpleftcol) + square(param()[p_topmass])/square(tmpleftcol.perp())));
      } else if (shower == Shower::b) {
        denom = square(tmpleftcol.squared_distance(tmpL)+square(param()[p_topmass])/square(tmpleftcol.perp()));
      } else
        throw NEW_EXCEPTION("ERROR in HQuark for antiquark: if antitop color connected has to be tshower or btopshower!! ");

      anglefac=numerator/denom;
    } else {
      if (tmpleftcolflav == Flavour::noflav) {
        anglefac = 1.0;
      } else {
        anglefac = tmpleftcol.squared_distance(tmpR)/(tmpR.squared_distance(tmpL) + tmpL.squared_distance(tmpleftcol));
      }
    }

    Hamiltonian = 8.0*M_PI*Cte::CF*alphas(mu2)/mu2*tmptot.perp()/tmpL.perp()*(1.0+square(tmpR.perp()/tmptot.perp()))*anglefac;
  }

  if (Hamiltonian < 0.0)
    throw NEW_EXCEPTION("Error - HQuark cannot be smaller 0.");

  return Hamiltonian;
}

double Deconstruction::Model::HGluonQQ(const fastjet::PseudoJet & tmptot,
                                const fastjet::PseudoJet & tmpmleftcol,
                                const fastjet::PseudoJet & tmpmrightcol,
                                const fastjet::PseudoJet & tmpL,
                                const fastjet::PseudoJet & tmpR,
                                const fastjet::PseudoJet & tmpgrandm,
                                const Flavour::id flavour,
                                const Flavour::id granflavour,
                                const Flavour::id tmpleftcolflav,
                                const Flavour::id tmprightcolflav,
                                const Shower::id shower) {
  double Hgqq = 8.0*M_PI*Cte::TR*alphas(square(tmptot.m()))*(square(tmpL.perp()) + square(tmpR.perp()))/square(tmptot.perp())/square(tmptot.m());

  return Hgqq;
}

double Deconstruction::Model::HGluonBB(const fastjet::PseudoJet & tmptot, 
				 const fastjet::PseudoJet & tmpmleftcol, 
				 const fastjet::PseudoJet & tmpmrightcol,
				 const fastjet::PseudoJet & tmpL,
				 const fastjet::PseudoJet & tmpR, 
				 const fastjet::PseudoJet & tmpgrandm,
				 const Flavour::id flavor,
				 const Flavour::id granflavor,
				 const Flavour::id tmpleftcolflav,
				 const Flavour::id tmprightcolflav,
				 const Shower::id shower)
{
  
  double mu2(square(tmptot.m())); // should be (q_b+q_bbar)^2 with q_b^2=q_bbar^2=mb^2 if bmass is considered
  //////// Hgbb //////////////////////////////
  double Hgbb(8.0*M_PI*Cte::TR*alphas(mu2)*(square(tmpL.perp()) + square(tmpR.perp()))/square(tmptot.perp())/mu2);

  if(Hgbb < 0.0)
    LOG(ERROR) << "MSP: Error - Hgbb cannot be smaller 0" << endl;

  return Hgbb;
}


double Deconstruction::Model::HGluonGG(const fastjet::PseudoJet & tmptot, 
				    const fastjet::PseudoJet & tmpmleftcol, 
				    const fastjet::PseudoJet & tmpmrightcol,
				    const fastjet::PseudoJet & tmpL,
				    const fastjet::PseudoJet & tmpR, 
				    const fastjet::PseudoJet & tmpgrandm,
				    const Flavour::id flavor,
				    const Flavour::id granflavor,
				    const Flavour::id tmpleftcolflav,
				    const Flavour::id tmprightcolflav,
				    const Shower::id shower)
{

  double Hggg(0.0);
  
  
  //////////////   Hggg  /////////////////////////////////

  double anglefac(1.0);

  fastjet::PseudoJet soft;
  fastjet::PseudoJet hard;
  fastjet::PseudoJet colcon;
  Flavour::id colconflav;

  if(tmpL.perp() < tmpR.perp())
    {
      soft = tmpL;
      hard = tmpR;
      colcon = tmpmleftcol;
      colconflav = tmpleftcolflav;
    }
  else
    {
      soft = tmpR;
      hard = tmpL;
      colcon = tmpmrightcol;
      colconflav = tmprightcolflav;
    }
    

  if(colconflav == Flavour::noflav)
    {
      anglefac = 1.0;
    }
  else
    {
    
      if(shower==Shower::t)
	{
	  if((colconflav == Flavour::t) || (colconflav == Flavour::tbar))
	    {
	      anglefac = ((soft.squared_distance(colcon)+square(param()[p_topmass])/square(colcon.perp()))*(hard.squared_distance(colcon) + square(param()[p_topmass])/square(colcon.perp())) - (square(param()[p_topmass])/square(colcon.perp()))*soft.squared_distance(hard));
	      
	      anglefac *= 1/((soft.squared_distance(colcon) + square(param()[p_topmass])/square(colcon.perp())) * (soft.squared_distance(hard) + soft.squared_distance(colcon) + square(param()[p_topmass])/square(colcon.perp())));	  
	    }
	  else
	    {
	      anglefac = hard.squared_distance(colcon)/(soft.squared_distance(hard)+soft.squared_distance(colcon));
	    }
	  
	}
      else if(shower == Shower::b)
	{
	  if((colconflav == Flavour::t) || (colconflav == Flavour::tbar))
	    {

	      anglefac = ((soft.squared_distance(colcon)+square(param()[p_topmass])/square(colcon.perp()))*(hard.squared_distance(colcon) + square(param()[p_topmass])/square(colcon.perp())) - (square(param()[p_topmass])/square(colcon.perp()))*soft.squared_distance(hard));
	      
	      anglefac *= 1/square(soft.squared_distance(colcon) + square(param()[p_topmass])/square(colcon.perp()));
	    }
	  else
	    {
	      anglefac = hard.squared_distance(colcon)/(soft.squared_distance(hard)+soft.squared_distance(colcon));
	    }
	}
      else
	{
	  if((colconflav == Flavour::t) || (colconflav == Flavour::tbar)) LOG(ERROR) << "colcon = top in wrong shower : " << shower << endl;
	  
	  anglefac = hard.squared_distance(colcon)/(soft.squared_distance(hard)+soft.squared_distance(colcon));
	}

    } // end if no color partner else
  
  Hggg = 8.0*M_PI*Cte::CA*alphas(square(tmptot.m()))*square(tmptot.perp())/(square(tmptot.m())*soft.perp()*hard.perp()) * square(1.0 - soft.perp()*hard.perp()/square(tmptot.perp())) * anglefac;

  if(Hggg < 0.0) LOG(ERROR) << "MSP: Error - Hggg cannot be small 0" << endl;
  
  return Hggg;
}  




double Deconstruction::Model::sudakovGluon(const fastjet::PseudoJet & tmptot, 
					const fastjet::PseudoJet & tmpmleftcol, 
					const fastjet::PseudoJet & tmpmrightcol,
					const fastjet::PseudoJet & tmpgrandm,
					const Flavour::id flavor,
					const Flavour::id granflavor,
					const Flavour::id tmpleftcolflav,
					const Flavour::id tmprightcolflav,
					const Shower::id shower)
{
  
  double mothermass2(square(tmptot.m()));  
  double kappaK(0.0);
  double mu2max(0.0);
  //following line checks for existence of grandmother
  if(tmpgrandm.e() < Cte::smalldouble)
    {
      kappaK = 2.0*square(tmptot.perp())/tmptot.perp();
    }
  else 
    {
      if(granflavor == Flavour::t || granflavor == Flavour::tbar)
	{
	  kappaK = (square(tmpgrandm.m())-square(param()[p_topmass]))/tmpgrandm.perp();
	}
      else
	{
	  kappaK = square(tmpgrandm.m())/tmpgrandm.perp();
	}
    }
  

  double Sggg(0.0), Sgqq(0.0);
  double deltaSggg(0.0);

  double colconAmass(0.0);
  double colconBmass(0.0);


  //////// Sggg //////////////////////////////////////////////////////////////////////
  double theta_kA(0.0);
  if(tmpmleftcol.e()<Cte::smalldouble)
    {
      theta_kA = sqrt(dR2_outside);
    }
  else
    {
      if(tmpleftcolflav == Flavour::t || tmpleftcolflav == Flavour::tbar) colconAmass = param()[p_topmass]; 
      theta_kA=sqrt(tmptot.squared_distance(tmpmleftcol) + square(colconAmass)/square(tmpmleftcol.perp()));
    }
  
  double theta_kB(0.0);
  if(tmpmrightcol.e()<Cte::smalldouble)
    {
      theta_kB = sqrt(dR2_outside);
    }
  else 
    {
      if(tmprightcolflav == Flavour::t || tmprightcolflav == Flavour::tbar) colconBmass = param()[p_topmass];
      theta_kB=sqrt(tmptot.squared_distance(tmpmrightcol) + square(colconBmass)/square(tmpmrightcol.perp()));
    }


  //// BparaL, BparaR, AparaL, AparaR

  double BparaL(0.0), BparaR(0.0), AparaL(0.0), AparaR(0.0);
  
  
  double fbargg((log(2.0)-11.0/12.0)*Cte::CA);
  double f0gg(Cte::CA);
  

  // sudakov for massive and massless color connected parnters the same except B-> B' and A->A'
  // the primed ones are for massive partners.


  //MSP: where do we need kappa?


     double SILmin(0.0);
     double SILmax(0.0);
     double SIRmin(0.0);
     double SIRmax(0.0);

     // left radiation:
     //mu2star= AparaL*square(tmptot.perp())*exp(BparaL+fbargg/f0gg);
     

     // right radiation:


  
  // if grandmass2 =0 ->error
 
  Sggg = Cte::CA/M_PI/square(Cte::b0)*(log(alphas(std::fabs(mothermass2))/alphas(tmptot.perp()*kappaK/2.0))*(1.0/alphas(theta_kA*theta_kB*square(tmptot.perp()))-11.0*Cte::b0/12.0) + 1.0/alphas(std::fabs(mothermass2)) - 1.0/alphas(tmptot.perp()*kappaK/2.0));
  
  
  if(Sggg<0) Sggg =0;
  
  if(shower == Shower::b)
    {
      if( tmpleftcolflav == Flavour::t || tmpleftcolflav == Flavour::tbar)
	{
	  deltaSggg = Cte::CA/2.0/M_PI/Cte::b0 *(log(square(theta_kA)/(square(param()[p_topmass])/square(tmpmleftcol.perp()))) + 1.0)*log(alphas(std::fabs(mothermass2))/alphas(tmptot.perp()*kappaK/2.0));
	}
      else if(tmprightcolflav == Flavour::t || tmprightcolflav == Flavour::tbar)
	{
	  deltaSggg = Cte::CA/2.0/M_PI/Cte::b0 *(log(square(theta_kB)/(square(param()[p_topmass])/square(tmpmrightcol.perp()))) + 1.0)*log(alphas(std::fabs(mothermass2))/alphas(tmptot.perp()*kappaK/2.0));
	}
      
      /*
      if(deltaSggg<0.0)
	{
	  LOG(DEBUG) << "ERROR: MSP deltaSggg<0 cannot be! "<< endl;
	}
      */

      if(deltaSggg < 0.0) deltaSggg = 0;

      //      if(deltaSggg < 0.0) LOG(DEBUG) << "MSP: ERROR deltaSggg < 0. Cannot be!" << endl;

      Sggg = Sggg + deltaSggg ;
    }

  ///////// Sgqq //////////////////////////////////////////////////////////////////////
  Sgqq = Cte::TR/3.0/M_PI/Cte::b0*log(alphas(std::fabs(mothermass2))/alphas(tmptot.perp()*kappaK/2.0));

  // flavor factor g-> all 5 quarks
  Sgqq = Cte::nf*Sgqq;  // Sgqq = Sgbb

  double Sg = Sggg + Sgqq;

  return(1.0);
  //  return(exp(-Sg));

}



double Deconstruction::Model::sudakovQuark(const fastjet::PseudoJet & tmptot, 
					const fastjet::PseudoJet & tmpmleftcol, 
					const fastjet::PseudoJet & tmpmrightcol,
					const fastjet::PseudoJet & tmpgrandm,
					const Flavour::id flavor,
					const Flavour::id granflavor,
					const Flavour::id tmpleftcolflav,
					const Flavour::id tmprightcolflav,
					const Shower::id shower)
{

  double kappaK(0);
  double virt2(square(tmptot.m()));
  if(flavor == Flavour::t || flavor == Flavour::tbar) virt2 -= square(param()[p_topmass]);

  //following line checks for existence of grandmother

  // use fix starting scale from top decay in W shower:
  if(shower == Shower::W)
    {
      if(Model::Wshower_level == Model::Winitlevel+1)
	{
	  if((tmpgrandm.e() < Cte::smalldouble)) LOG(ERROR) << "MSP: ERROR! Grandm shouldnt be empty in W shower" << endl;
	  kappaK = 2.0*square(param()[p_wmass])/tmpgrandm.perp(); // grandperp should be perp of W
	}
      else
	{
	  kappaK=square(tmpgrandm.m())/tmpgrandm.perp();
	}
    }


  // use fix starting scale from top decay in b shower:  
  if(shower == Shower::b)
    {
      if(Model::btopshower_level == Model::binitlevel+1)
	{
	  if((tmpgrandm.e() < Cte::smalldouble)) LOG(ERROR) << "MSP: ERROR! Grandm shouldnt be empty in btop shower" << endl;
	  kappaK = 2.0*(square(param()[p_topmass])-square(param()[p_wmass]))/tmpgrandm.perp(); //tmpgrand is top here      
	}
      else
	{
	  kappaK=square(tmpgrandm.m())/tmpgrandm.perp();
	}
    }

  if(shower==Shower::QCD) {      
    if(tmpgrandm.e() < Cte::smalldouble) kappaK=2.0*square(tmptot.m())/tmptot.perp();
    else kappaK=square(tmpgrandm.m())/tmpgrandm.perp();}


  // check if topshower at the beginning:
  if(shower==Shower::t)
    { 
      
      if(Model::tshower_level == Model::tinitlevel+1)
	{
	  kappaK = 2.0*(square(tmptot.perp())+square(param()[p_topmass]))/tmptot.perp();
	}
      else
	{
	  // MSP Following needs to be changed... minus param()[p_topmass]^2 is wrong!!
	  if(granflavor == Flavour::t || granflavor == Flavour::tbar)
	    {
	      kappaK = (square(tmpgrandm.m())-square(param()[p_topmass]))/tmpgrandm.perp();
	    }
	  else if(granflavor == Flavour::noflav)
	    { 
	      LOG(ERROR) << "in sudakov Quark: find no flav for GrandM ERROR" << endl;
	    }
	  else kappaK = square(tmpgrandm.m())/tmpgrandm.perp();
	}

    }
  
  double Sqqg(0.0);


  if(flavor == Flavour::t)
    {
      if(shower != Shower::t) LOG(ERROR) << "flavor Flavour::t but shower not tshower, cant be..." << endl;
      double drtotcol2(dR2_outside);
      if(tmpmrightcol.e()<Cte::smalldouble)
	{
	  drtotcol2 = dR2_outside;
	}
      else drtotcol2=tmptot.squared_distance(tmpmrightcol);
      
     
      // msp - check here the mass in the alphas term
      Sqqg = Cte::CF/M_PI/Cte::b0*log(alphas(std::fabs(virt2))/alphas(tmptot.perp()*kappaK/2.0))*log((Cte::alphattg*drtotcol2*square(tmptot.perp())+2.0*Cte::betattg*square(param()[p_topmass]))/square(param()[p_topmass]));
      
      
      if(Sqqg <0) Sqqg=0;

      double Sdecay(0.0);
      
      if(std::fabs(square(tmptot.m())-square(param()[p_topmass]))<param()[p_topmass]*delta_topmass)
	{
	  if(kappaK/2.0*tmptot.perp() < param()[p_topmass]*delta_topmass)
	    {
	      Sdecay = log(atan(kappaK/2.0*tmptot.perp()/param()[p_topmass]/topwidth));
	    }
	  else Sdecay = log(atan(delta_topmass/topwidth));
	  
	  Sdecay -= log(atan(std::fabs(square(tmptot.m())-square(param()[p_topmass]))/param()[p_topmass]/topwidth));
	}

            
      if(Sdecay < 0)
	{
	  LOG(ERROR) << "MSP: ERROR Sdecay for top cannot be < 0 " << endl;
	  LOG(DEBUG) << "splitting sudakov: " << Sqqg << endl;
	  LOG(DEBUG) << "kappaK: " << kappaK << endl;
	  LOG(DEBUG) << "topmass: " << param()[p_topmass] << endl;
	  LOG(DEBUG) << "topwidth: " << topwidth << endl;
	  LOG(DEBUG) << "delta_topmass: " << delta_topmass << endl;
	  LOG(DEBUG) << "jetmass: " << tmptot.m() << endl;
	  LOG(DEBUG) << "delta_top: " << delta_topmass << endl;
	  LOG(DEBUG) << "pT jet: " << tmptot.perp() << endl;
	  LOG(DEBUG) << "log(atan(kappaK/2.0*tmptot.perp()/param()[p_topmass]/delta_topmass)): " << log(atan(kappaK/2.0*tmptot.perp()/param()[p_topmass]/delta_topmass)) << endl;
	  LOG(DEBUG) << "log(atan(std::fabs(square(tmptot.m())-square(param()[p_topmass]))/param()[p_topmass]/topwidth)) " << log(atan(std::fabs(square(tmptot.m())-square(param()[p_topmass]))/param()[p_topmass]/topwidth)) << endl;
	}
	

    
      Sqqg += Sdecay;

    }
  else if(flavor == Flavour::tbar)
    {
      
      if(shower != Shower::t) LOG(ERROR) << "flavor Flavour::tbar but shower not tshower, cant be..." << endl;
      double drtotcol2(dR2_outside);
      if(tmpmleftcol.e()<Cte::smalldouble)
	{
	  drtotcol2 = dR2_outside;
	}
      else drtotcol2=tmptot.squared_distance(tmpmleftcol);

      // msp - check here the mass in the alphas term
      Sqqg = Cte::CF/M_PI/Cte::b0*log(alphas(std::fabs(virt2))/alphas(tmptot.perp()*kappaK/2.0))*log((drtotcol2*square(tmptot.perp())+2*square(param()[p_topmass]))/square(param()[p_topmass]));


      if(Sqqg <0) Sqqg=0;



      double Sdecay(0.0);
      
      if(std::fabs(square(tmptot.m())-square(param()[p_topmass]))<param()[p_topmass]*delta_topmass)
	{
	  if(kappaK/2.0*tmptot.perp() < param()[p_topmass]*delta_topmass)
	    {
	      Sdecay = log(atan(kappaK/2.0*tmptot.perp()/param()[p_topmass]/delta_topmass));
	    }
	  else Sdecay = log(atan(delta_topmass/topwidth));
	  
	  Sdecay -= log(atan(std::fabs(square(tmptot.m())-square(param()[p_topmass]))/param()[p_topmass]/topwidth));
	}
      
      if(Sdecay < 0.0) LOG(ERROR) << "MSP: ERROR Sdecay cannot be < 0 " << endl;
      
      Sqqg += Sdecay;


    }
  else if(flavor == Flavour::q || flavor == Flavour::b)
    {
      
      double drtotcol2(dR2_outside);
      if(tmpmrightcol.e()<Cte::smalldouble)
	{
	  drtotcol2 = dR2_outside;
	}
      else 
	{
	  if(tmprightcolflav == Flavour::t ||tmprightcolflav == Flavour::tbar)
	    {
	      drtotcol2=tmptot.squared_distance(tmpmrightcol) + square(param()[p_topmass])/square(tmpmrightcol.perp());
	    }
	  else
	    {
	      drtotcol2=tmptot.squared_distance(tmpmrightcol);
	    }
	}

      //LOG(DEBUG) << "drtotcol2 : " << drtotcol2 << endl;
      Sqqg = Cte::CF/M_PI/square(Cte::b0)*(log(alphas(std::fabs(virt2))/alphas(tmptot.perp()*kappaK/2.0))*
			       (1.0/alphas(drtotcol2*square(tmptot.perp())) - 3.0*Cte::b0/4.0) + 1.0/alphas(std::fabs(virt2))-1.0/alphas(tmptot.perp()*kappaK/2.0));

      if(Sqqg <0) Sqqg=0;

      if((shower == Shower::b) && (tmprightcolflav == Flavour::t ||tmprightcolflav == Flavour::tbar))
	{
	  double deltaSqqg = Cte::CF/M_PI/Cte::b0*(log((tmpmrightcol.squared_distance(tmptot)+square(param()[p_topmass])/square(tmpmrightcol.perp()))/(square(param()[p_topmass])/square(tmpmrightcol.perp()))) + 1.0);
	  deltaSqqg *= log(alphas(std::fabs(virt2))/(alphas(tmptot.perp()*kappaK/2.0)));

	  /*
	  if(deltaSqqg < 0.0)
	    {
	      LOG(DEBUG) << "MSP: ERROR deltaSqqg < 0: " << deltaSqqg << " Cannot be!" << endl;
	      LOG(DEBUG) << "flavor: "  << flavor  << endl;
	      LOG(DEBUG) << "granflavor: " << granflavor << endl;
	      LOG(DEBUG) << "tmpmleftcolflav: " << tmpleftcolflav << endl;
	      LOG(DEBUG) << "tmpmrightcolflav: " << tmprightcolflav << endl;
	      LOG(DEBUG) << "shower: " << shower << endl;
	      LOG(DEBUG) << "tmptot: " << endl;
	      printjet(tmptot);
	      LOG(DEBUG) << "tmpmleftcol: " << endl;
	      printjet(tmpmleftcol);
	      LOG(DEBUG) << "tmprightcol: " << endl;
	      printjet(tmpmrightcol);
	      LOG(DEBUG) << "kappaK : "<< kappaK << endl;
	    }
	  */
	  if(deltaSqqg < 0.0) deltaSqqg = 0.0;

	  Sqqg = Sqqg + deltaSqqg;
	}

    }
  else if(flavor == Flavour::qbar || flavor == Flavour::bbar)
    {

      double drtotcol2(dR2_outside);
      if(tmpmleftcol.e()<Cte::smalldouble)
	{
	  drtotcol2 = dR2_outside;
	}
      else 
	{
	  if(tmpleftcolflav == Flavour::t ||tmpleftcolflav == Flavour::tbar)
	    {
	      drtotcol2=tmptot.squared_distance(tmpmleftcol) + square(param()[p_topmass])/square(tmpmleftcol.perp());
	    }
	  else
	    {
	      drtotcol2=tmptot.squared_distance(tmpmleftcol);
	    }
	}

      //LOG(DEBUG) << "drtotcol2 : " << drtotcol2 << endl;
      Sqqg = Cte::CF/M_PI/square(Cte::b0)*(log(alphas(std::fabs(virt2))/alphas(tmptot.perp()*kappaK/2.0))*
			       (1.0/alphas(drtotcol2*square(tmptot.perp())) - 3.0*Cte::b0/4.0) + 1.0/alphas(std::fabs(virt2))-1.0/alphas(tmptot.perp()*kappaK/2.0));

      if(Sqqg <0) Sqqg=0;


      if((shower == Shower::b) && (tmpleftcolflav == Flavour::t ||tmpleftcolflav == Flavour::tbar))
	{
	  double deltaSqqg = Cte::CF/M_PI/Cte::b0*(log((tmpmleftcol.squared_distance(tmptot)+square(param()[p_topmass])/square(tmpmleftcol.perp()))/(square(param()[p_topmass])/square(tmpmleftcol.perp()))) + 1.0);
	  deltaSqqg *= log(alphas(square(tmptot.m()))/(alphas(tmptot.perp()*kappaK/2.0)));
	  
	  if(deltaSqqg< 0.0) deltaSqqg = 0.0;
	  Sqqg = Sqqg + deltaSqqg;
	}

    }



  //  return(exp(-Sqqg));
  return(1.0);

}




double Deconstruction::Model::sudakovTopEND(const fastjet::PseudoJet & tmptot, 
				      const fastjet::PseudoJet & tmpmleftcol, 
				      const fastjet::PseudoJet & tmpmrightcol,
				      const fastjet::PseudoJet & tmpgrandm,
				      const Flavour::id flavor,
				      const Flavour::id granflavor,
				      const Flavour::id tmpleftcolflav,
				      const Flavour::id tmprightcolflav,
				      const Shower::id shower)
{

  double grandmass2(square(tmpgrandm.m())-square(param()[p_topmass]));
  double grandperp(tmpgrandm.perp());
  //  double mothermass2(topmass*topwidth);
  double mothermass2(square(tmptot.m()));

  //following line checks for existence of grandmother
  if(tmpgrandm.e() < Cte::smalldouble)
    {
      grandmass2 = 2.0*(square(tmptot.perp())+square(param()[p_topmass]));
      grandperp = tmptot.perp();
    }

  double Sqqg(0.0);

  double drtotcol2(dR2_outside);
  if(flavor == Flavour::t)
    {
      if(tmpmrightcol.e()<Cte::smalldouble)
	{
	  drtotcol2 = dR2_outside;
	}
      else drtotcol2=tmptot.squared_distance(tmpmrightcol);
    }
  else if(flavor == Flavour::tbar)
    {
      if(tmpmleftcol.e()<Cte::smalldouble)
	{
	  drtotcol2 = dR2_outside;
	}
      else drtotcol2=tmptot.squared_distance(tmpmleftcol);
     
    }
  else LOG(ERROR) << "in sudakovTopEND: wrong flavor " << endl;

  Sqqg = Cte::CF/M_PI/Cte::b0*log(alphas(std::fabs(mothermass2))/alphas(tmptot.perp()*grandmass2/2.0/grandperp)) * log((drtotcol2*square(tmptot.perp())+2.0*square(param()[p_topmass]))/square(param()[p_topmass]));
  
  if(Sqqg <0) Sqqg=0;

  double Sdecay(0.0);
  
  if(std::fabs(square(tmptot.m())-square(param()[p_topmass]))<param()[p_topmass]*delta_topmass)
    {
      if(grandmass2/grandperp/2.0*tmptot.perp() < param()[p_topmass]*delta_topmass)
	{
	  Sdecay = log(atan(grandmass2/grandperp/2.0*tmptot.perp()/param()[p_topmass]/delta_topmass));
	}
      else Sdecay = log(atan(delta_topmass/topwidth));
      
      Sdecay -= log(atan(std::fabs(square(tmptot.m())-square(param()[p_topmass]))/param()[p_topmass]/topwidth));
    }

  if(Sdecay < 0.0) LOG(ERROR) << "MSP: Error - Sdecay cannot be smaller 0" << endl;

  //  return(exp(-Sqqg+Sdecay));
  return(1.0);
}




double Deconstruction::Model::sudakov(const fastjet::PseudoJet & tmptot, 
				   const fastjet::PseudoJet & tmpmleftcol, 
				   const fastjet::PseudoJet & tmpmrightcol,
				   const fastjet::PseudoJet & tmpgrandm,
				   const Flavour::id flavor,
				   const Flavour::id granflavor,
				   const Flavour::id tmpleftcolflav,
				   const Flavour::id tmprightcolflav,
				   const Shower::id shower)
{

  // calculate virtuality of J
  double mu(tmptot.m());
  if(flavor == Flavour::t || flavor == Flavour::tbar) {
    mu = std::sqrt(std::fabs(square(tmptot.m())-square(param()[p_topmass])));
  }

  // calculate mumax (maximum scale), depends if first splitting after decay,... (many cases)
  double mumax(0.0);
    
  if(shower == Shower::W)
    {
      if(Model::Wshower_level == Model::Winitlevel+1)
	{
	  if((tmpgrandm.e() < Cte::smalldouble)) LOG(ERROR) << "MSP: ERROR! Grandm shouldnt be empty in W shower" << endl;
	  mumax = sqrt(tmptot.perp()*square(param()[p_wmass])/tmpgrandm.perp()); // grandperp should be perp of W
	}
      else
	{
	  mumax=sqrt(tmptot.perp()*square(tmpgrandm.m())/(2.0*tmpgrandm.perp()));
	}
    }
  else if(shower == Shower::b)
    {
      if(Model::btopshower_level == Model::binitlevel+1)
	{
	  if((tmpgrandm.e() < Cte::smalldouble)) LOG(ERROR) << "MSP: ERROR! Grandm shouldnt be empty in btop shower" << endl;
	  mumax = sqrt(tmptot.perp()*(square(param()[p_topmass])-square(param()[p_wmass]))/tmpgrandm.perp()); //tmpgrand is top here      
	}
      else
	{
	  mumax=sqrt(tmptot.perp()*square(tmpgrandm.m())/(2.0*tmpgrandm.perp()));
	}
    }
  else if(shower==Shower::QCD) 
    {      
      if(tmpgrandm.e() < Cte::smalldouble)
	{
	  mumax=tmptot.perp();
	}
      else
	{
	  mumax=sqrt(tmptot.perp()/2.0*square(tmpgrandm.m())/tmpgrandm.perp());
	}
    }
  else if(shower==Shower::t)
    { 
      
      if(Model::tshower_level == Model::tinitlevel+1)
	{
	  mumax = sqrt(tmptot.perp()*(square(tmptot.perp())+square(param()[p_topmass]))/tmptot.perp());
	}
      else
	{
	  if(granflavor == Flavour::t || granflavor == Flavour::tbar)
	    {
	      mumax = sqrt(tmptot.perp()/2.0*(square(tmpgrandm.m())-square(param()[p_topmass]))/tmpgrandm.perp());
	    }
	  else if(granflavor == Flavour::noflav)
	    { 
	      LOG(ERROR) << "in sudakov Quark: find no flav for GrandM ERROR" << endl;
	    }
	  else mumax = sqrt(tmptot.perp()/2.0*square(tmpgrandm.m())/tmpgrandm.perp());
	}

    }
  else LOG(ERROR) << "in sudakov(): Wrong shower" << shower << endl;

  double f0(0.0);
  double fbar(0.0);
  
  double SIL(0.0);
  double SIR(0.0);
  double SILmumax(0.0);
  double SIRmumax(0.0);
  double Sgqq(0.0);

  double stotal(0.0);

  double stopdecay(0.0);
  double stoprad(0.0);

  if(flavor == Flavour::tbar)
    {
      // antitop can radiate to the left
      // antitop can decay
      double theta2kJL(0.0);
 
      if(tmpmleftcol.e()<Cte::smalldouble)
	{
	  theta2kJL = dR2_outside;
	}
      else theta2kJL=tmptot.squared_distance(tmpmleftcol);
      
      double rhotildeq(square(param()[p_topmass])/square(tmptot.perp()));
      double x0(2.0-1.0*theta2kJL/(theta2kJL+rhotildeq) + 0.1 *log(rhotildeq/(theta2kJL+rhotildeq)));
      double x1(-3.0-1.0*theta2kJL/(theta2kJL+rhotildeq) + 0.4 *log(rhotildeq/(theta2kJL+rhotildeq)));

      double Iphiz(2.0*Cte::CF*((theta2kJL+2.0*rhotildeq)/(theta2kJL+rhotildeq) * log(2.0 + theta2kJL/rhotildeq) - 1.0));

      double tscale0((theta2kJL*square(tmptot.perp())+square(param()[p_topmass]))*exp(x0));
      double tscale1((theta2kJL*square(tmptot.perp())+square(param()[p_topmass]))*exp(x1));

      if(tscale0 <= square(mu))
	{
	  //	  LOG(DEBUG) << "1" << endl;
	  SIL = 0.0;
	}
      else if((tscale1 < square(mu)) &&
	      (tscale0 >= square(mu)))
	{
	  //	  LOG(DEBUG) << "2" << endl;
	  SIL = 1.0/alphas(square(mu)) - 1.0/(alphas(tscale0));
	  SIL += (1/alphas(tscale0) + x0*Cte::b0)*log(alphas(square(mu))/(alphas(tscale0)));
	  SIL *= 1/2.0/M_PI/(x0-x1)/square(Cte::b0)*Iphiz;
	}
      else if(tscale1 >= square(mu))
	{
	  //	  LOG(DEBUG) << "3" << endl;
	  SIL = 1.0/alphas(tscale1) - 1.0/(alphas(tscale0));
	  SIL += (1/alphas(tscale0) + x0*Cte::b0)*log(alphas(tscale1)/(alphas(tscale0)));
	  SIL *= 1/2.0/M_PI/(x0-x1)/square(Cte::b0)*Iphiz;

	  SIL += 1/(2.0*M_PI*Cte::b0)*Iphiz*log(alphas(square(mu))/alphas(tscale1));
	}
      else LOG(ERROR) << "error in sudakov for antitop: Scales dont match" << endl;

      ////// now max part:

      if(tscale0 <= square(mumax))
	{
	  //	  LOG(DEBUG) << "1 max" << endl;
	  SILmumax = 0.0;
	}
      else if((tscale1 < square(mumax)) &&
	      (tscale0 >= square(mumax)))
	{
	  
	  //	  LOG(DEBUG) << "2 max" << endl;

	  SILmumax = 1.0/alphas(square(mumax)) - 1.0/(alphas(tscale0));
	  SILmumax += (1/alphas(tscale0) + x0*Cte::b0)*log(alphas(square(mumax))/(alphas(tscale0)));
	  SILmumax *= 1/2.0/M_PI/(x0-x1)/square(Cte::b0)*Iphiz;

	}
      else if(tscale1 >= square(mumax))
	{
	  //	  LOG(DEBUG) << "3 max" << endl;

	  SILmumax = 1.0/alphas(square(mumax)) - 1.0/(alphas(tscale0));
	  SILmumax += (1/alphas(tscale0) + x0*Cte::b0)*log(alphas(square(mumax))/(alphas(tscale0)));
	  SILmumax *= 1/2.0/M_PI/(x0-x1)/square(Cte::b0)*Iphiz;
	  SILmumax += 1/(2.0*M_PI*Cte::b0)*Iphiz*log(alphas(square(mumax))/alphas(tscale1));
	  
	}
      else LOG(ERROR) << "error in sudakov for antitop: Scales dont match" << endl;

      stoprad = SIL - SILmumax;
	  
      /// antitop decay
      double grandmass2(square(tmpgrandm.m())-square(param()[p_topmass]));
      double grandperp(tmpgrandm.perp());
      //  double mothermass2(topmass*topwidth);
      double mothermass2(square(tmptot.m()));
      

      
      double kappaK(0.0);
      
      if(granflavor == Flavour::t || granflavor == Flavour::tbar)
	{
	  kappaK = (square(tmpgrandm.m())-square(param()[p_topmass]))/tmpgrandm.perp();
	}
      else kappaK = square(tmpgrandm.m())/tmpgrandm.perp();


      if(fabs(square(tmptot.m())-square(param()[p_topmass]))<param()[p_topmass]*delta_topmass)
	{
	  if(kappaK/2.0*tmptot.perp() < param()[p_topmass]*delta_topmass)
	    {
	     stopdecay = log(atan(kappaK/2.0*tmptot.perp()/param()[p_topmass]/topwidth));
	    }
	  else stopdecay = log(atan(delta_topmass/topwidth));
	  
	  stopdecay -= std::log(std::atan(std::fabs(square(tmptot.m())-(param()[p_topmass]*param()[p_topmass]))/param()[p_topmass]/topwidth));
	}
      

      stotal = max(stoprad,0.0) + stopdecay;

    }
  else if(flavor == Flavour::t)
    {
       
      // top can radiate to the right
      // top can decay
      double theta2kJR(0.0);
 
      if(tmpmrightcol.e()<Cte::smalldouble)
	{
	  theta2kJR = dR2_outside;
	}
      else theta2kJR=tmptot.squared_distance(tmpmrightcol);
      
      
      double rhotildeq(square(param()[p_topmass])/square(tmptot.perp()));
      double x0(2.0-1.0*theta2kJR/(theta2kJR+rhotildeq) + 0.1 *log(rhotildeq/(theta2kJR+rhotildeq)));
      double x1(-3.0-1.0*theta2kJR/(theta2kJR+rhotildeq) + 0.4 *log(rhotildeq/(theta2kJR+rhotildeq)));

      double Iphiz(2.0*Cte::CF*((theta2kJR+2.0*rhotildeq)/(theta2kJR+rhotildeq) * log(2.0 + theta2kJR/rhotildeq) - 1.0));
  
      double tscale0((theta2kJR*square(tmptot.perp())+square(param()[p_topmass]))*exp(x0));
      double tscale1((theta2kJR*square(tmptot.perp())+square(param()[p_topmass]))*exp(x1));

      if( tscale0 <= square(mu))
	{
	  SIR=0.0;
	}      
      else if((tscale1 < square(mu)) &&
	      (tscale0 >= square(mu)) )
	{
	  SIR = 1.0/alphas(square(mu)) - 1.0/(alphas(tscale0));
	  SIR += (1/alphas(tscale0) + x0*Cte::b0)*log(alphas(square(mu))/(alphas(tscale0)));
	  SIR *= 1/2.0/M_PI/(x0-x1)/square(Cte::b0)*Iphiz;

	}
      else if(tscale1 >= square(mu))
	{
	  SIR = 1.0/alphas(tscale1) - 1.0/(alphas(tscale0));
	  SIR += (1/alphas(tscale0) + x0*Cte::b0)*log(alphas(tscale1)/(alphas(tscale0)));
	  SIR *= 1/2.0/M_PI/(x0-x1)/square(Cte::b0)*Iphiz;
	  
	  SIR += 1/(2.0*M_PI*Cte::b0)*Iphiz*log(alphas(square(mu))/alphas(tscale1));

	}
      else LOG(ERROR) << "ERROR in sudakovtop: scales do not match " << endl;
      

      if( tscale0 <= square(mumax))
	{
	  SIRmumax=0.0;
	}      
      else if((tscale1 < square(mumax)) &&
	      (tscale0 >= square(mumax)) )
	{
	  SIRmumax = 1.0/alphas(square(mumax)) - 1.0/(alphas(tscale0));
	  SIRmumax += (1/alphas(tscale0) + x0*Cte::b0)*log(alphas(square(mumax))/(alphas(tscale0)));
	  SIRmumax *= 1/2.0/M_PI/(x0-x1)/square(Cte::b0)*Iphiz;

	}
      else if(tscale1 >= square(mumax))
	{
	  SIRmumax = 1.0/alphas(tscale1) - 1.0/(alphas(tscale0));
	  SIRmumax += (1/alphas(tscale0) + x0*Cte::b0)*log(alphas(tscale1)/(alphas(tscale0)));
	  SIRmumax *= 1/2.0/M_PI/(x0-x1)/square(Cte::b0)*Iphiz;
	  
	  SIRmumax += 1/(2.0*M_PI*Cte::b0)*Iphiz*log(alphas(square(mumax))/alphas(tscale1));
	}
      else LOG(ERROR) << "ERROR in sudakovtop: scales do not match " << endl;

      stoprad = SIR - SIRmumax;

      /// top decay:

      double grandmass2(square(tmpgrandm.m())-square(param()[p_topmass]));
      double grandperp(tmpgrandm.perp());
      //  double mothermass2(topmass*topwidth);
      double mothermass2(square(tmptot.m()));
      
      //following line checks for existence of grandmother

      double kappaK(0.0);
      
      if (tmpgrandm.perp() != 0) {
        if(granflavor == Flavour::t || granflavor == Flavour::tbar)
  	  {
	    kappaK = (square(tmpgrandm.m())-square(param()[p_topmass]))/tmpgrandm.perp();
	  }
        else kappaK = square(tmpgrandm.m())/tmpgrandm.perp();
      } else {
        kappaK = 2.0*tmptot.perp()*param()[p_topmass]*delta_topmass+999; // to avoid NaN (DANILO)
      }

      if(fabs(square(tmptot.m())-square(param()[p_topmass]))<param()[p_topmass]*delta_topmass)
	{
	  if(kappaK/2.0*tmptot.perp() < param()[p_topmass]*delta_topmass) // this fails if tmpgrandm.perp() == 0 from above comment
	    {
	      stopdecay = log(atan(kappaK/2.0*tmptot.perp()/param()[p_topmass]/topwidth));
	    }
	  else stopdecay = std::log(std::atan(delta_topmass/topwidth));
	  
	  stopdecay -= std::log(std::atan(std::fabs(square(tmptot.m())-(param()[p_topmass]*param()[p_topmass]))/param()[p_topmass]/topwidth));
	}

      stotal = max(stoprad,0.0) + stopdecay;

    }
  else if(flavor == Flavour::q || flavor == Flavour::qbar || flavor == Flavour::b || flavor == Flavour::bbar || flavor == Flavour::g)
    {
      
      if(flavor == Flavour::g)
	{
	  f0 = Cte::CA;
	  fbar = (log(2.0)-11.0/12.0)*Cte::CA;
	}
      else if(flavor == Flavour::q || flavor == Flavour::qbar || flavor == Flavour::b || flavor == Flavour::bbar)
	{
	  f0 = 2.0*Cte::CF;
	  fbar = -3.0/2.0*Cte::CF;
	}
      else LOG(ERROR) << "in sudakov: Wrong flavor: " << flavor << endl;

      
      double rhoL(0.0);
      double theta2kJL(0.0);
      double rhoR(0.0);
      double theta2kJR(0.0);

      double AparaR(0.0);
      double AparaL(0.0);
      double BparaR(0.0);
      double BparaL(0.0);
      /// left color connected partner
      if(tmpleftcolflav == Flavour::t || tmpleftcolflav == Flavour::tbar) // massive
	{
	  rhoL = square(param()[p_topmass])/square(tmpmleftcol.perp());

	  if(tmpmleftcol.e()<Cte::smalldouble)
	    {
	      theta2kJL = dR2_outside;
	    }
	  else theta2kJL=tmptot.squared_distance(tmpmleftcol);

	  // q-t dipole (btopshower) or not
	  if(shower == Shower::b) // a q-t dipole
	    {
	      AparaL=square(rhoL+theta2kJL)/rhoL;
	      BparaL=-1;
	    }
	  else
	    {
	      AparaL=square(rhoL+theta2kJL)/(2*rhoL+theta2kJL);
	      BparaL=rhoL/(rhoL+theta2kJL)*log(rhoL/(2*rhoL+theta2kJL));
	    }
	}
      else ///// massless
	{
	  rhoL = 0.0;
	  
	  if(tmpmleftcol.e()<Cte::smalldouble)
	    {
	      theta2kJL = dR2_outside;
	    }
	  else theta2kJL=tmptot.squared_distance(tmpmleftcol);
	  
	  AparaL=theta2kJL;
	  BparaL=0.0;
	  
	}

      // right color connected partner
      if(tmprightcolflav == Flavour::t || tmprightcolflav == Flavour::tbar) // massive
	{
	  rhoR = square(param()[p_topmass])/square(tmpmrightcol.perp());

	  
	  if(tmpmrightcol.e()<Cte::smalldouble)
	    {
	      theta2kJR = dR2_outside;
	    }
	  else theta2kJR=tmptot.squared_distance(tmpmrightcol);

	  // q-t dipole (btopshower) or not
	  if(shower == Shower::b) // a q-t dipole
	    {
	      AparaR=square(rhoR+theta2kJR)/rhoR;
	      BparaR=-1;
	    }
	  else
	    {
	      AparaR=square(rhoR+theta2kJR)/(2*rhoR+theta2kJR);
	      BparaR=rhoR/(rhoR+theta2kJR)*log(rhoR/(2*rhoR+theta2kJR));
	    }
	}
      else  /// massless
	{
	  rhoR = 0.0;

	  if(tmpmrightcol.e()<Cte::smalldouble)
	    {
	      theta2kJR = dR2_outside;
	    }
	  else theta2kJR=tmptot.squared_distance(tmpmrightcol);
	  
	  AparaR=theta2kJR;
	  BparaR=0.0;
	}


      double mustarL(AparaL*square(tmptot.perp())*exp(BparaL+fbar/f0));
      double mustarR(AparaR*square(tmptot.perp())*exp(BparaR+fbar/f0));


      if((flavor == Flavour::g || flavor == Flavour::qbar || flavor == Flavour::bbar) && (mustarL > square(mu)))
	{	  
	  SIL = f0/2.0/M_PI/square(Cte::b0)/alphas(AparaL*square(tmptot.perp()))+(fbar+BparaL*f0)/2.0/M_PI/Cte::b0;
	  SIL *= log(alphas(square(mu))/alphas(mustarL));
	  SIL -= f0/2.0/M_PI/square(Cte::b0)*(1/alphas(mustarL)-1/alphas(square(mu))); 
	}
	 
      if((flavor == Flavour::g || flavor == Flavour::q || flavor == Flavour::b) && (mustarR > square(mu)))
	{
	  SIR = f0/2.0/M_PI/square(Cte::b0)/alphas(AparaR*square(tmptot.perp()))+(fbar+BparaR*f0)/2.0/M_PI/Cte::b0;
	  SIR *= log(alphas(square(mu))/alphas(mustarR));
	  SIR -= f0/2.0/M_PI/square(Cte::b0)*(1/alphas(mustarR)-1/alphas(square(mu)));
	}

      if((flavor == Flavour::g || flavor == Flavour::qbar || flavor == Flavour::bbar)  && (mustarL > square(mumax)))
	{
	  SILmumax = f0/2.0/M_PI/square(Cte::b0)/alphas(AparaL*square(tmptot.perp()))+(fbar+BparaL*f0)/2.0/M_PI/Cte::b0;
	  SILmumax *= log(alphas(square(mumax))/alphas(mustarL));
	  SILmumax -= f0/2.0/M_PI/square(Cte::b0)*(1/alphas(mustarL)-1/alphas(square(mumax))); 
	}
	 
      if((flavor == Flavour::g || flavor == Flavour::q || flavor == Flavour::b)  && (mustarR > square(mumax)))
	{
	  SIRmumax = f0/2.0/M_PI/square(Cte::b0)/alphas(AparaR*square(tmptot.perp()))+(fbar+BparaR*f0)/2.0/M_PI/Cte::b0;
	  SIRmumax *= log(alphas(square(mumax))/alphas(mustarR));
	  SIRmumax -= f0/2.0/M_PI/square(Cte::b0)*(1/alphas(mustarR)-1/alphas(square(mumax)));
	}

      if(flavor == Flavour::g)
	{
	  
	  Sgqq = Cte::TR/3.0/M_PI/Cte::b0*log(alphas(square(mu))/alphas(square(mumax)));
	  
	  // flavor factor g-> all 5 quarks
	  Sgqq = Cte::nf*Sgqq;  // Sgqq = Sgbb

	  stotal = max(SIL-SILmumax + SIR-SIRmumax,0.0) + Sgqq;
	}

      if(flavor == Flavour::q || flavor == Flavour::b) stotal = max(SIR-SIRmumax,0.0);

      if(flavor == Flavour::qbar || flavor == Flavour::bbar) stotal = max(SIL-SILmumax,0.0);
    }
  else 
    {
      LOG(ERROR) << "In sudakov: Wrong flavor!!" << endl;
    } 


  // Question to Dave:
  // Eq. 47 and Eq. 82 are both supposed to be good for massive color connected partner. The only difference is that Eq. 82 is valid in t-shower and Eq. 47 is valid in the b-shower??

  return(exp(-stotal));
}


double Deconstruction::Model::M2(std::vector<fastjet::PseudoJet> &transformed) {
  return 1.0;
}

