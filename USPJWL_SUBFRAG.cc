// -*- C++ -*-

// This is a Rivet analysis for JEWEL  
// Subjet fragmentation based on ALICE arXiv:2204.10270 (hepdata: https://www.hepdata.net/record/ins2070434)
//
// --Leonardo Barreto, IF USP, 2022

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/PseudoJet.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/JetShape.hh"
#include "Rivet/Projections/SubtractedJewelEvent.hh"
#include "Rivet/Projections/SubtractedJewelFinalState.hh"
#include <string>

namespace Rivet {

  // Declaration of functions


  class USPJWL_SUBFRAG : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(USPJWL_SUBFRAG);

    void init() {

      // Subjet fragmentation based on arXiv:2204.10270 (hepdata: https://www.hepdata.net/record/ins2070434)

      // Grab variable jet R parameter from environment, default value of 0.4
      RJETS = getenv("RJETS") ? getenv("RJETS") : "0.4";
      RJETS_f = std::stof(RJETS);
      std::cout << "\nR chosen for jet algorithm: " << RJETS << std::endl;


      // abseta range on ALICE TPC:
      Cut cut(Cuts::abseta < 0.9);

      SubtractedJewelEvent sev(1.0);
      SubtractedJewelFinalState fs(sev, cut);
      declare(fs, "FS");
      ChargedFinalState cfs(fs);
      declare(cfs, "CFS");

      FastJets cfj(cfs, FastJets::ANTIKT, RJETS_f);
      declare(cfj, "ChargedJets");

      // Book histograms
      // Full: pp analysis, High: 80 < pT < 120 GeV, HighD: 100 < pT < 150 GeV, 
      // Custom: very detailed and full range 
      // for each r = [0.1, 0.2]

      book(zfull_1,"z_Full_r01", PTEDGES_FULL);
      book(zhigh_1,"z_High_r01", PTEDGES_HIGH);
      book(zhighd_1,"z_HighD_r01", PTEDGES_HIGHD);
      book(zcustom_1,"z_Custom_r01", 25, 0.50001, 1.00001);

      book(zfull_2,"z_Full_r02", PTEDGES_FULL);
      book(zhigh_2,"z_High_r02", PTEDGES_HIGH);
      book(zhighd_2,"z_HighD_r02", PTEDGES_HIGHD);
      book(zcustom_2,"z_Custom_r02", 25, 0.50001, 1.00001);
      
      // Counter for a better control on the inclusive and full range normalizations
      // First bin (0): 80 < pT < 120 GeV, second bin (1): 100 < pT < 150 GeV
      book(jetcount, "Number_Jets", 2, -0.5, 1.5);

    }


    /// Perform the per-event analysis
    void analyze(const Event& evt) {

      // Get jets of event
      double etamax = 0.9 - RJETS_f;
      Cut jetcuts = Cuts::pT > 80 * GeV && Cuts::pT < 150 * GeV 
                    && Cuts::abseta < etamax;

      const vector<double> rs = {0.1, 0.2};
      const Jets& jets = apply<FastJets>(evt, "ChargedJets").jetsByPt(jetcuts);
 
      for (const Jet& j : jets) {
        // Apply jet algorithm on jets constituents to calculate z_r
       
	Particles jetconsti = j.constituents();
        double jpt = j.pT();

        //std::cout << "\nJet pt = " << jpt << std::endl;

        for (double r : rs) {
 
          // Apply jet algorithm on jets constituents to calculate z_r
          // using the fastjet classes (arXiv:1111.6097)
          JetDefinition def_subjet(fastjet::kt_algorithm, r);
          ClusterSequence cs(jetconsti, def_subjet);

          // Grab leading subjet
          PseudoJets subjets = sorted_by_pt(cs.inclusive_jets());
          PseudoJet lead_subjet = subjets[0];

          double z_lead = lead_subjet.perp() / jpt;

          //std::cout << "z lead = " << z_lead << std::endl;

          vector<Histo1DPtr> histos;
          if (r == 0.1) { 
            histos = {zfull_1, zhigh_1, zhighd_1, zcustom_1}; 
          }
          else { 
            histos = {zfull_2, zhigh_2, zhighd_2, zcustom_2}; 
          } 

          // Select correct histogram
          if (jpt < 120 * GeV) { 
            histos[1] -> fill(z_lead);
            jetcount -> fill(0.);
          }
          
          if (jpt > 100 * GeV) {
            histos[2] -> fill(z_lead);
            jetcount -> fill(1.);
          }
          
          // Fill Custom for all pt
          histos[3] -> fill(z_lead);

          // Inclusive calculation
          for (auto subj : subjets) {
            double z = subj.perp() / jpt;
            //std::cout << "z = " << z << std::endl;
            histos[0] -> fill(z);
          }
        }
      }

    }


    void finalize() {
      // Scale only after yoda merge
    }


    /// @name Histograms
    // R_AA
    Histo1DPtr zfull_1, zhigh_1, zhighd_1, zcustom_1,
               zfull_2, zhigh_2, zhighd_2, zcustom_2,
               jetcount;


    double RJETS_f;
    std::string RJETS;

    std::vector<double> PTEDGES_FULL = {0., 0.02, 0.04, 0.1, 0.3, 0.6, 0.7, 
                                        0.77, 0.83, 0.89, 0.95, 1.00001};
    std::vector<double> PTEDGES_HIGH = {0.6, 0.7, 0.77, 0.83, 0.89, 0.95, 1.00001}; 
    std::vector<double> PTEDGES_HIGHD = {0.7, 0.75, 0.77, 0.8, 0.83, 0.86, 0.9, 
                                        0.92, 0.95, 0.98, 1.00001};

    

  };



  DECLARE_RIVET_PLUGIN(USPJWL_SUBFRAG);



}
