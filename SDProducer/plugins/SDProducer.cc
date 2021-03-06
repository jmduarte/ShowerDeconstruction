// -*- C++ -*-
//
// Package:    ShowerDeconstruction/SDProducer
// Class:      SDProducer
 
/**\class SDProducer SDProducer.cc ShowerDeconstruction/SDProducer/plugins/SDProducer.cc

 Description: Allow calling Shower Deconstruction (1211.3140) from CMSSW

*/
//
// Original Author:  Gregor Kasieczka
//         Created:  Mon, 05 Jan 2015 13:11:25 GMT
//
//


// system include files
#include <memory>
#include <limits>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"

#include "ShowerDeconstruction/SDAlgorithm/interface/AnalysisParameters.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/Exception.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/Message.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/TopGluonModel.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/BackgroundModel.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/ISRModel.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/Deconstruct.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/ParseUtils.h"


//
// class declaration
//

class SDProducer : public edm::EDProducer {
   public:
      explicit SDProducer(const edm::ParameterSet&);
      ~SDProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      // ----------member data ---------------------------
     double microjet_cone_;
     std::string fatjet_name_;  
     std::string input_card_;
};


//
// constructors and destructor
//
SDProducer::SDProducer(const edm::ParameterSet& iConfig) :
  microjet_cone_(iConfig.getParameter<double>("MicrojetCone")),
  fatjet_name_(iConfig.getParameter<std::string>("FatjetName")),
  input_card_(iConfig.getParameter<std::string>("InputCard"))
{  
  produces<edm::ValueMap<double> >("chi");  
  produces<edm::ValueMap<int> >("nmicrojets");  
}


SDProducer::~SDProducer()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
SDProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace Deconstruction;


   // Get Fatjets
   edm::Handle<reco::PFJetCollection> fatjets;
   iEvent.getByLabel(fatjet_name_, fatjets);

   // Store the output here (later transfer to value map)
   std::vector<float> values;
   std::vector<int> nmicrojets;
   values.reserve(fatjets->size());
   nmicrojets.reserve(fatjets->size());

   // Setup Shower Deconstruction 
   // Read input parameters
   AnalysisParameters param(input_card_);

   TopGluonModel *signal = 0;
   BackgroundModel *background = 0;
   ISRModel *isr = 0;
   Deconstruct *deconstruct = 0;
   
   signal = new TopGluonModel(param);
   background = new BackgroundModel(param);
   isr = new ISRModel(param);
   deconstruct = new Deconstruct(param, *signal, *background, *isr);

   // Loop over the fatjets
   for (std::vector<reco::PFJet>::const_iterator fj_it = fatjets->begin();
	fj_it != fatjets->end();
	++fj_it){
     
     // First get the constituents of the jet
     // these will be of type Constituents (which is a edm::Ptr<Candidate>)
     reco::Jet::Constituents tmp_constituents = fj_it->getJetConstituents();

     
     // Second: Convert them to fastjet PseudoJet objects
     std::vector<fastjet::PseudoJet> constituents;          
     for (reco::Jet::Constituents::const_iterator it = tmp_constituents.begin();
	  it != tmp_constituents.end();
	  ++it){     
       fastjet::PseudoJet pj((*it)->px(), 
			     (*it)->py(), 
			     (*it)->pz(), 
			     (*it)->energy());
       constituents.push_back(pj);
     } // end of constituent loop

     // Second and a half: If the microjet cone-size is smaller than zero
     // determine is a function of the fatjet-pT
     // Otherwise: use the passed value
     double microconesize;
     if  (microjet_cone_<0){
       // From Tobias:
       // 0..500   -> 0.3
       // 500..700 -> 0.2
       // 700..inf -> 0.
       if (fj_it->pt() < 500)
	 microconesize=0.3;
       else if (fj_it->pt() < 700)
	 microconesize=0.2;
       else
	 microconesize=0.1;             
     }
     else{
       microconesize = microjet_cone_;
     }
            
     // Third: recluster to microjets and count them
     fastjet::JetDefinition reclustering(fastjet::JetAlgorithm::kt_algorithm, microconesize);
     fastjet::ClusterSequence * cs_micro = new fastjet::ClusterSequence(constituents, reclustering);
     std::vector<fastjet::PseudoJet> microjets = fastjet::sorted_by_pt(cs_micro->inclusive_jets(10.));
     nmicrojets.push_back(microjets.size());
     
     //     std::cout << "SDProducer " << microjets.size() << std::endl;

     // Fourth: Make sure we have at most nine microjets
     if (microjets.size()>9) 
       microjets.erase(microjets.begin()+9,microjets.end()); 


     // Fifth: run shower deconstruction
     double Psignal = 0.0;
     double Pbackground = 0.0;
     try {
       double chi = deconstruct->deconstruct(microjets, Psignal, Pbackground);	       
       values.push_back(chi);
     } catch(Deconstruction::Exception &e) {       
       std::cout << "Exception while running SD: " << e.what() << std::endl;
       values.push_back(std::numeric_limits<double>::quiet_NaN());
     }
               
   } // end of fatjet loop

   // Create chi value map and add to event
   std::auto_ptr<edm::ValueMap<double> > out_chi(new edm::ValueMap<double>());
   edm::ValueMap<double>::Filler filler(*out_chi);
   filler.insert(fatjets, values.begin(), values.end());
   filler.fill();
   iEvent.put(out_chi, "chi");       

   // Create nmicrojets value map and add to event
   std::auto_ptr<edm::ValueMap<int> > out_nmj(new edm::ValueMap<int>());
   edm::ValueMap<int>::Filler filler_nmj(*out_nmj);
   filler_nmj.insert(fatjets, nmicrojets.begin(), nmicrojets.end());
   filler_nmj.fill();
   iEvent.put(out_nmj, "nmicrojets");       
}


// ------------ method called once each job just before starting event loop  ------------
void 
SDProducer::beginJob()
{
}


// ------------ method called once each job just after ending the event loop  ------------
void 
SDProducer::endJob() {
}

 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SDProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SDProducer);
