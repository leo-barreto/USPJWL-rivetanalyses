// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/JetShape.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "Rivet/Projections/SubtractedJewelEvent.hh"
#include "Rivet/Projections/SubtractedJewelFinalState.hh"

//Not sure if I must include these yet, probably not since Rivet already does it
#include "fstream"
#include "cmath"
#include "string.h"

namespace Rivet {

      

      class USPJWL_JET_MASS : public Analysis {
            public:

                  USPJWL_JET_MASS()
                        : Analysis("USPJWL_JET_MASS")
                  {   }

                  

                  //! Cells which make up the grid  
                  struct candidate{
                        //! these are the boundaries 
                        double etaMin;
                        double phiMin;
                        double etaMax;
                        double phiMax;
                        //! Pseudojets will have the objects inside 
                        Jets objects;
                        vector <int> jetID;
                        //! Four momentum to get useful informations like eta, phi, mass, pT etc... 
                        FourMomentum candMom;
                        FourMomentum bkgMom;
                        double sumnegpT;
                        double sumpT;
                  };
            

                  

                  

                  void init() {

                        
                        //! verbosity flag
                        verbose = false;

                        //! Jet Radius
                        _jetR=0.4;
                        //! Jet Kinematic cuts 
                        _pTCut=0.15;

                        //! Grid initialization  
                        _etaMax=0.9;
                        _phiMax=M_PI;
                        _delRMin=0.05;

                        //! cell boundaries 
                        _Nbounds_eta = (int)2*_etaMax/(_delRMin);
                        _Nbounds_phi = (int)2*_phiMax/(_delRMin);

                        for(int i = 0; i<=_Nbounds_eta; ++i){
                              double etax = -1*_etaMax + i*_delRMin;
                              _etabins.push_back(etax);
                        }

                        for(int i = 0; i<=_Nbounds_phi; ++i){
                              double phix = -1*_phiMax + i*_delRMin;
                              _phibins.push_back(phix);	
                        }

                        if(verbose)
                              std::cout<<"Grid we are using is ("<<_Nbounds_eta<<" x "<<_Nbounds_phi
                                    <<")   eta: [-"<<_etaMax<<", "<<_etaMax
                                    <<"]   phi: [-"<<_phiMax<<", "<<_phiMax<<"]"<<std::endl;      

                        // Initialise and register projections

                        //Final State Particles with pseudo-rapidity cuts
                        Cut cut(Cuts::abseta<0.9);
                        //Cut cut(Cuts::pt>0.300*GeV && Cuts::abseta<4.5);
                        
                        SubtractedJewelEvent sev(1.0);
                        SubtractedJewelFinalState fs(sev, cut);

                        //const FinalState fs(cut);
                        declare(fs, "FS");

                        //Final State Particles with pseudo-rapidity cuts
                        //ALICE <0.9
                        //ATLAS <2.something
                        Cut ChargedCut(Cuts::abseta<0.9);
                        const ChargedFinalState cfs(fs);
                        declare(cfs, "CFS");

                        //Aplying Fast-Jet algorithms
                        //Anti-kt Algorithm R=0.4
                        FastJets ak_04(fs, FastJets::ANTIKT, 0.4);
                        ak_04.useInvisibles();
                        declare(ak_04, "AntiKt_04");


                        
                        vector<double> mass_edges=linspace(200,0.0,100.0);
                        book(_hs_mass[0],"Jet_Mass_60_80",mass_edges);
                        book(_hs_mass[1],"Jet_Mass_80_100",mass_edges);
                        book(_hs_mass[2],"Jet_Mass_100_120",mass_edges);
                        book(_hs_mass[3],"Jet_Mass_120_140",mass_edges);
                        book(_hs_mass[4],"Jet_Mass_140_160",mass_edges);
                        book(_hs_mass[5],"Jet_Mass_160_180",mass_edges);
                        book(_hs_mass[6],"Jet_Mass_180_200",mass_edges);
                        book(_hs_mass[7],"Jet_Mass_200_220",mass_edges);
                        book(_hs_mass[8],"Jet_Mass_220_240",mass_edges);
                        book(_hs_mass[9],"Jet_Mass_240_260",mass_edges);
                        book(_hs_mass[10],"Jet_Mass_260_280",mass_edges);
                        book(_hs_mass[11],"Jet_Mass_280_300",mass_edges);
                        book(_hs_mass[12],"Jet_Mass_300",mass_edges);



                        vector<double> pt_edges=linspace(50,20.0,520.0);

                        book(_h_JetpT_NSub_04,"JetpT_NSub_04",pt_edges);
                       

                        
                  }

                  /// Perform the per-evt analysis
                  void analyze(const Event& evt){


                        //Here I create my array of jetsets to ensure I have
                        //fewer variable names, this array contains both
                        //jet w/o thermal subtraction and also jets
                        //with 3 different subtraction procedures,
                        //and also a CDFMidpoint.


                        PseudoJets jetAr[4];
                        

                        //! **************************************
                        //! Jet Collection from all particles in the evt
                        //! Used for w/o recoils and in vacuum 
                        if(verbose) std::cout<<"Jet Collection built without subtraction"<<std::endl;
                        Cut cuts = (Cuts::abseta < _etaMax) & (Cuts::pT > _pTCut*GeV);

                        const FastJets& AJets_04 = apply<FastJets>(evt, "AntiKt_04");
                        const Jets jets_noSub_04 = AJets_04.jetsByPt(cuts);


                        //! **************************************
                        //! BKG-Sub Jet Collection from all particles w/ recoils inclding the dummy and Scattering centers 
                        //Here begins analysis for pt of jets with R=0.4 whithout any subtraction procedure
                        //if(verbose) std::cout<<"Jet Collection w/ 4MomSub Subtraction"<<std::endl;
                        //PseudoJets jets_4MomSub_04 = do4MomSub(jets_noSub_04, pscat, doSubtraction);
                        //jetAr[1] = jets_4MomSub_04;  
                        for(const Jet& jet: jets_noSub_04) {
                              if(jet.abseta()<0.5 && jet.pt()>20.0){
                                    _h_JetpT_NSub_04->fill(jet.pt());
                              }
                        }





                        //Jet mass analysis here for jets with R=0.4
                        //using 4MomSub method
                        for(const Jet& jet: jets_noSub_04){
                              // The leading jet
                              //const PseudoJet& jet = jetAr[alg][0];
                              const double m   = jet.mass();
                              const double eta = jet.eta();
                              const double pt  = jet.pt();

                              if(m>=0 && abs(eta)<(_etaMax-_jetR)){
                                    if(60.0 <= pt && pt < 80.0){
                                          _hs_mass[0]->fill(m/GeV);
                                    }
                                    if(80.0 <= pt && pt < 100.0){
                                          _hs_mass[1]->fill(m/GeV);
                                    }
                                    if(100.0 <= pt && pt < 120.0){
                                          _hs_mass[2]->fill(m/GeV);
                                    }
                                    if(120.0 <= pt && pt < 140.0){
                                          _hs_mass[3]->fill(m/GeV);
                                    }
                                    if(140.0 <= pt && pt < 160.0){
                                          _hs_mass[4]->fill(m/GeV);
                                    }
                                    if(160.0 <= pt && pt < 180.0){
                                          _hs_mass[5]->fill(m/GeV);
                                    }
                                    if(180.0 <= pt && pt < 200.0){
                                          _hs_mass[6]->fill(m/GeV);
                                    }
                                    if(200.0 <= pt && pt < 220.0){
                                          _hs_mass[7]->fill(m/GeV);
                                    }
                                    if(220.0 <= pt && pt < 240.0){
                                          _hs_mass[8]->fill(m/GeV);
                                    }
                                    if(240.0 <= pt && pt < 260.0){
                                          _hs_mass[9]->fill(m/GeV);
                                    }
                                    if(260.0 <= pt && pt < 280.0){
                                          _hs_mass[10]->fill(m/GeV);
                                    }
                                    if(280.0 <= pt && pt < 300.0){
                                          _hs_mass[11]->fill(m/GeV);
                                    }
                                    if(300.0 <= pt){
                                          _hs_mass[12]->fill(m/GeV);
                                    }
                              }
                        }
                  }                  

                  void finalize(){
                        //std::cout << _h_NinPlane->sumW() << std::endl;
                        //std::cout << _h_Nout->sumW() << std::endl;
                  }

            protected:

                  double _pTCut;
                  double _jetR;
                  bool verbose;
                  vector<double> _etabins;
                  vector<double> _phibins;
                  int _Nbounds_eta;
                  int _Nbounds_phi;
                  double _delRMin;
                  double _etaMax;
                  double _phiMax;

            private:

                                    

                  vector<double> _ptedges; // This can't be a raw array if we want to initialise it non-painfully

                  vector<double> _massedges; // This can't be a raw array if we want to initialise it non-painfully

                 
                  

                  

                  Histo1DPtr _hs_mass[13];

                  Histo1DPtr _h_JetpT_NSub_04;

                  

                  
                  

                  
                  
      };

      // The hook for the plugin system
      DECLARE_RIVET_PLUGIN(USPJWL_JET_MASS);

}   
