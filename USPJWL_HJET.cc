// -*- C++ -*-
#include "Rivet/Analysis.hh"

#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/JetShape.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"


#include "Rivet/Projections/JetAlg.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"


#include "Rivet/Particle.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "Rivet/Tools/Logging.hh"

#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"

#include "Rivet/Projections/SubtractedJewelEvent.hh"
#include "Rivet/Projections/SubtractedJewelFinalState.hh"


#include "HepMC/PdfInfo.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/SimpleVector.h"
#include "HepMC/GenVertex.h"

#include "fstream"
#include "cmath"
#include "string.h"



namespace Rivet {

      

      class USPJWL_HJET : public Analysis {
            public:

                  USPJWL_HJET()
                        : Analysis("USPJWL_HJET")
                  {   }


            void init() {

                  


                  //output.open("output.dat");
                  //output << "pid \t phi \t eta \t pt \t counter_hadrons \t" << "\n";
                  //output << "===================================\n";

                  

                  //RJETS = getenv("RJETS"), RJETS_f = std::stof(RJETS);
                        

                  RJETS_f = 0.4;
                  Ntrig = 0;
                  etamax = 0.9;                  
                  etamax_jet = etamax - RJETS_f;


                  counter_hadrons = 0;                  
                  counter_hadrons8_9 = 0;
                  counter_hadrons6_7 = 0;
                  counter_hadrons12_50 = 0;
                  counter_hadrons1 = 0;
                  counter_hadrons_eta = 0;


                  std::cout << "\nR jet algorithm: " << RJETS_f << std::endl;



                  Cut cut(Cuts::abseta<0.9);


                  SubtractedJewelEvent sev(1.0);
                  SubtractedJewelFinalState fs(sev, cut);

                  //const FinalState fs(cut);
                  declare(fs, "FS");

                  const ChargedFinalState cfs(fs);
                  declare(cfs, "CFS");


                  // Apply FastJet
                  FastJets cfj(cfs, FastJets::ANTIKT, RJETS_f);                        
                  declare(cfj, "C_Jets");


                  vector<double> hjet_edges=linspace(100,0.0,100.0);
                  //vector<double> edges_20_50=linspace(100,19.5,50.5);
                  //vector<double> edges_8_9=linspace(100,7.5,9.5);
                  //vector<double> edges_1=linspace(100,0.5,100.5);

                  

                  
                  
                  // Book histograms
                  book(_hs_pTJet,"Njet_20_50",hjet_edges);
                  book(_hs_pTJet_12_50,"Njet_12_50",hjet_edges);
                  book(_hs_pTJet_6_7,"Njet_6_7",hjet_edges);
                  book(_hs_pTJet_8_9,"Njet_8_9",hjet_edges);
                  book(_hs_pTJet_1,"Njet_1",hjet_edges);
                  book(_hs_pTJet_eta,"Njet_eta",hjet_edges);

                  book(_hs_Ntrig,"hNtrig_20_50",hjet_edges);
                  book(_hs_Ntrig_12_50,"hNtrig_12_50",hjet_edges);
                  book(_hs_Ntrig_8_9,"hNtrig_8_9",hjet_edges);
                  book(_hs_Ntrig_6_7,"hNtrig_6_7",hjet_edges);
                  book(_hs_Ntrig_1,"hNtrig_1",hjet_edges);
                  book(_hs_Ntrig_eta,"hNtrig_eta",hjet_edges);

                  book(_hs_pTJet_all,"Njet_all_20_50",hjet_edges);
                  book(_hs_pTJet_all_12_50,"Njet_all_12_50",hjet_edges);
                  book(_hs_pTJet_all_8_9,"Njet_all_8_9",hjet_edges);
                  book(_hs_pTJet_all_6_7,"Njet_all_6_7",hjet_edges);
                  book(_hs_pTJet_all_1,"Njet_all_1",hjet_edges);
                  book(_hs_pTJet_all_eta,"Njet_all_eta",hjet_edges);

                  
                  





                  
                  
            }


            /// Perform the per-event analysis
            void analyze(const Event& evt) {

                  



                  //20 < pT,trig < 50 GeV/c, denoted by TT{20,50} 
                  Cut partcuts = Cuts::pT > 20 * GeV && Cuts::pT < 50 * GeV && Cuts::abseta < etamax;
                  //8 < pT,trig < 9 GeV/c, denoted by TT{8,9} 
                  Cut partcuts8_9 = Cuts::pT > 8 * GeV && Cuts::pT < 9 * GeV && Cuts::abseta < etamax;
                  Cut partcuts6_7 = Cuts::pT > 6 * GeV && Cuts::pT < 7 * GeV && Cuts::abseta < etamax;
                  Cut partcuts12_50 = Cuts::pT > 12 * GeV && Cuts::pT < 50 * GeV && Cuts::abseta < etamax;
                  Cut partcuts1 = Cuts::pT > 1 * GeV && Cuts::abseta < etamax;
                  Cut partcuts_eta = Cuts::abseta < etamax;


                  

                  

                  

                  //JETS
                  //Cut for jets by paper 20 < pT < 100 GeV/c for R=0.2 and R=0.4
                  Cut jetcuts = Cuts::pT >= 0.15 * GeV && Cuts::pT <= 100.0 * GeV && Cuts::abseta < etamax_jet;
                  const Jets alljets = apply<FastJets>(evt, "C_Jets").jetsByPt(jetcuts);



                  //PARTICLES
                  auto Particles  = evt.allParticles(partcuts);
                  auto Particles8_9  = evt.allParticles(partcuts8_9);
                  auto Particles6_7  = evt.allParticles(partcuts6_7);
                  auto Particles12_50  = evt.allParticles(partcuts12_50);
                  auto Particles1  = evt.allParticles(partcuts1);
                  auto Particles_eta  = evt.allParticles(partcuts_eta);

                  //output << "size Particles = " << Particles.size() << "\n";
                  //output << "size Particles8_9 = " << Particles8_9.size() << "\n"; 
                  //output << "size Particles8_15 = " << Particles8_15.size() << "\n"; 
                  //output << "size Particles_1 = " << Particles1.size() << "\n"; 
                  //output << "size Particles_eta = " << Particles_eta.size() << "\n";  




                  for (size_t i = 0; i < Particles.size(); i++) {

                    //double phi = Particles[i].phi(), eta = Particles[i].eta(), pt = Particles[i].pT();
                    double phi = Particles[i].phi(), pt = Particles[i].pT();
                    int pid = Particles[i].pid();

                    //Charged hadrons selection
                    if (PID::isHadron(pid) && PID::isCharged(pid)){
                        counter_hadrons+=1;
                        //Ntrig+=evt.weight();

                        //particles identification
                        //output << pid << "\t" << phi << "\t" << eta << "\t" << pt << "\t" << counter_hadrons<< "\t" << Ntrig << "\n";  
                        //output << "Particle (20-50) pT = " << pt << "\n";


                        _hs_Ntrig->fill(pt/GeV);
                        //JETS 
                        for(const Jet& j: alljets){
                              //Jets identification
                              //double phi_j = j.phi(), eta_j = j.eta(), pt_j = j.pT();
                              double phi_j = j.phi(), pt_j = j.pT();
                              

                              _hs_pTJet_all->fill(pt_j/GeV);
                              //selection condition for jets: difference in azimuthal angle >= pi - 0.6
                              if(deltaPhi(phi,phi_j) >= M_PI - 0.6){

                                    //output << "jet" << "\t" << phi_j << "\t" << eta_j << "\t" << pt_j << "\t" << counter_hadrons<< "\t" << "\n";

                                    //jet spectrum histogram 
                                    _hs_pTJet->fill(pt_j/GeV);
                              }
                        }
                        
                    }
                    
              
                  }

                  for (size_t i = 0; i < Particles8_9.size(); i++) {

                    //double phi = Particles8_9[i].phi(), eta = Particles8_9[i].eta(), pt = Particles8_9[i].pT();
                    int pid = Particles8_9[i].pid(), pt = Particles8_9[i].pT();
                    double phi = Particles8_9[i].phi();

                    //Charged hadrons selection
                    if (PID::isHadron(pid) && PID::isCharged(pid)){
                        counter_hadrons8_9+=1.;
                        //Ntrig8_9+=evt.weight();


                        _hs_Ntrig_8_9->fill(pt/GeV);
                        //particles identification
                        //output << pid << "\t" << phi << "\t" << eta << "\t" << pt << "\t" << counter_hadrons8_9 << "\t" << Ntrig8_9 << "\n";  

                        //JETS 
                        for(const Jet& j: alljets){
                              //Jets identification
                              //double phi_j = j.phi(), eta_j = j.eta(), pt_j = j.pT();
                              double phi_j = j.phi(), pt_j = j.pT();

                              _hs_pTJet_all_8_9->fill(pt_j/GeV);
                              //selection condition for jets: difference in azimuthal angle >= pi - 0.6
                              if(deltaPhi(phi,phi_j) >= M_PI - 0.6){

                                    //output << "jet" << "\t" << phi_j << "\t" << eta_j << "\t" << pt_j << "\t" << counter_hadrons8_9 << "\t" << "\n";

                                    //jet spectrum histogram 
                                    _hs_pTJet_8_9->fill(pt_j/GeV);
                              }
                        }
                        
                    }

                  for (size_t i = 0; i < Particles6_7.size(); i++) {

                    //double phi = Particles6_7[i].phi(), eta = Particles6_7[i].eta(), pt = Particles6_7[i].pT();
                    int pid = Particles6_7[i].pid(), pt = Particles6_7[i].pT();
                    double phi = Particles6_7[i].phi();

                    //Charged hadrons selection
                    if (PID::isHadron(pid) && PID::isCharged(pid)){
                        counter_hadrons6_7+=1;
                        //Ntrig+=evt.weight();

                        //particles identification
                        //output << pid << "\t" << phi << "\t" << eta << "\t" << pt << "\t" << counter_hadrons<< "\t" << Ntrig << "\n";  
                        //output << "Particle (6-7) pT = " << pt << "\n";


                        _hs_Ntrig_6_7->fill(pt/GeV);
                        //JETS 
                        for(const Jet& j: alljets){
                              //Jets identification
                              //double phi_j = j.phi(), eta_j = j.eta(), pt_j = j.pT();
                              double phi_j = j.phi(), pt_j = j.pT();

                              _hs_pTJet_all_6_7->fill(pt_j/GeV);
                              //selection condition for jets: difference in azimuthal angle >= pi - 0.6
                              if(deltaPhi(phi,phi_j) >= M_PI - 0.6){

                                    //output << "jet" << "\t" << phi_j << "\t" << eta_j << "\t" << pt_j << "\t" << counter_hadrons<< "\t" << "\n";

                                    //jet spectrum histogram 
                                    _hs_pTJet_6_7->fill(pt_j/GeV);
                              }
                        }
                        
                    }
                    
              
                  }  
                    
              
                  }
                  for (size_t i = 0; i < Particles1.size(); i++) {

                    //double phi = Particles1[i].phi(), eta = Particles1[i].eta(), pt = Particles1[i].pT();
                    int pid = Particles1[i].pid(), pt = Particles1[i].pT();
                    double phi = Particles1[i].phi();

                    //Charged hadrons selection
                    if (PID::isHadron(pid) && PID::isCharged(pid)){
                        counter_hadrons1+=1.;
                        //Ntrig1+=evt.weight();


                        _hs_Ntrig_1->fill(pt/GeV);
                        //particles identification
                        //output << pid << "\t" << phi << "\t" << eta << "\t" << pt << "\t" << counter_hadrons1 << "\t" << Ntrig1 << "\n";  

                        //JETS 
                        for(const Jet& j: alljets){
                              //Jets identification
                              //double phi_j = j.phi(), eta_j = j.eta(), pt_j = j.pT();
                              double phi_j = j.phi(), pt_j = j.pT();


                              _hs_pTJet_all_1->fill(pt_j/GeV);
                              //selection condition for jets: difference in azimuthal angle >= pi - 0.6
                              if(deltaPhi(phi,phi_j) >= M_PI - 0.6){

                                    //output << "jet" << "\t" << phi_j << "\t" << eta_j << "\t" << pt_j << "\t" << counter_hadrons1 << "\t" << "\n";

                                    //jet spectrum histogram 
                                    _hs_pTJet_1->fill(pt_j/GeV);
                              }
                        }

                    
                    }

              
                  }
                  for (size_t i = 0; i < Particles_eta.size(); i++) {

                    //double phi = Particles_eta[i].phi(), eta = Particles_eta[i].eta(), pt = Particles_eta[i].pT();
                    int pid = Particles_eta[i].pid(), pt = Particles_eta[i].pT();
                    double phi = Particles_eta[i].phi();

                    //Charged hadrons selection
                    if (PID::isHadron(pid) && PID::isCharged(pid)){
                        counter_hadrons_eta+=1.;
                        //Ntrig_eta+=evt.weight();


                        _hs_Ntrig_eta->fill(pt/GeV);
                        //particles identification
                        //output << pid << "\t" << phi << "\t" << eta << "\t" << pt << "\t" << counter_hadrons_eta << "\t" << Ntrig_eta << "\n";  

                        //JETS 
                        for(const Jet& j: alljets){
                              //Jets identification
                              //double phi_j = j.phi(), eta_j = j.eta(), pt_j = j.pT();
                              double phi_j = j.phi(), pt_j = j.pT();

                              _hs_pTJet_all_eta->fill(pt_j/GeV);
                              //selection condition for jets: difference in azimuthal angle >= pi - 0.6
                              if(deltaPhi(phi,phi_j) >= M_PI - 0.6){

                                    //output << "jet" << "\t" << phi_j << "\t" << eta_j << "\t" << pt_j << "\t" << counter_hadrons_eta << "\t" << "\n";

                                    //jet spectrum histogram 
                                    _hs_pTJet_eta->fill(pt_j/GeV);
                              }
                        }
                        //
                    }
                    
              
                  }
                  for (size_t i = 0; i < Particles12_50.size(); i++) {

                    //double phi = Particles12_50[i].phi(), eta = Particles12_50[i].eta(), pt = Particles12_50[i].pT();
                    int pid = Particles12_50[i].pid(), pt = Particles12_50[i].pT();
                    double phi = Particles12_50[i].phi();

                    //Charged hadrons selection
                    if (PID::isHadron(pid) && PID::isCharged(pid)){
                        counter_hadrons12_50+=1.;
                        //Ntrig_eta+=evt.weight();


                        _hs_Ntrig_12_50->fill(pt/GeV);
                        //particles identification
                        //output << pid << "\t" << phi << "\t" << eta << "\t" << pt << "\t" << counter_hadrons_eta << "\t" << Ntrig_eta << "\n";  

                        //JETS 
                        for(const Jet& j: alljets){
                              //Jets identification
                              //double phi_j = j.phi(), eta_j = j.eta(), pt_j = j.pT();
                              double phi_j = j.phi(), pt_j = j.pT();


                              _hs_pTJet_all_12_50->fill(pt_j/GeV);
                              //selection condition for jets: difference in azimuthal angle >= pi - 0.6
                              if(deltaPhi(phi,phi_j) >= M_PI - 0.6){

                                    //output << "jet" << "\t" << phi_j << "\t" << eta_j << "\t" << pt_j << "\t" << counter_hadrons_eta << "\t" << "\n";

                                    //jet spectrum histogram 
                                    _hs_pTJet_12_50->fill(pt_j/GeV);
                              }
                        }
                        //
                    }
                    
              
                  }


                  //hadron trigger spectrum histogram
                  //_hs_Ntrig->fill(counter_hadrons);    //N_trig spectrum
                  //_hs_Ntrig_8_9->fill(counter_hadrons8_9);    //N_trig spectrum
                  //_hs_Ntrig_1->fill(counter_hadrons1);
                  //_hs_Ntrig_eta->fill(counter_hadrons_eta);


                  
                   
                //output << "===================================\n";
            }


            void finalize() {
                  //output.close();



                  scale(_hs_pTJet,  1/(2*etamax_jet));      //Normalization and divide by dN/d_eta 
                  scale(_hs_pTJet_8_9,  1/(2*etamax_jet));
                  scale(_hs_pTJet_6_7,  1/(2*etamax_jet));
                  scale(_hs_pTJet_12_50,  1/(2*etamax_jet));
                  scale(_hs_pTJet_1,  1/(2*etamax_jet));
                  scale(_hs_pTJet_eta,  1/(2*etamax_jet));



                  
            }


            //constants
            double RJETS_f;
            double Ntrig;
            double Ntrig1;
            double Ntrig8_9;
            double Ntrig_eta;
            double etamax;
            double etamax_jet;
            double counter_hadrons;                  
            double counter_hadrons8_9;
            double counter_hadrons6_7;
            double counter_hadrons12_50;
            double counter_hadrons1;
            double counter_hadrons_eta;


            //std::ofstream output;

            Histo1DPtr _hs_Ntrig, _hs_pTJet, _hs_Ntrig_8_9, _hs_pTJet_8_9,
            _hs_Ntrig_1, _hs_pTJet_1, _hs_Ntrig_eta, _hs_pTJet_eta, _hs_pTJet_all,
            _hs_pTJet_all_8_9, _hs_pTJet_all_1, _hs_pTJet_all_eta, _hs_Ntrig_6_7, 
            _hs_pTJet_6_7, _hs_pTJet_all_6_7, _hs_Ntrig_12_50, _hs_pTJet_12_50,
            _hs_pTJet_all_12_50;


            

            //Scatter2DPtr _hs_hjet;
            

            };
                        

      // The hook for the plugin system
      DECLARE_RIVET_PLUGIN(USPJWL_HJET);


}
