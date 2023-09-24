// -*- C++ -*-

// This is a Rivet analysis for JEWEL  
// Jet phi distribution based on ATLAS arXiv:2111.06606 (hepdata: https://www.hepdata.net/record/ins1967021)
// vn calculation is done post Rivet with auxiliary scripts
// Applications can be found in
// arXiv:2208.02061 and https://doi.org/10.11606/D.43.2021.tde-05112021-191914
// 
// --Leonardo Barreto, IF USP, 2021

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/PseudoJet.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "Rivet/Projections/SubtractedJewelEvent.hh"
#include "Rivet/Projections/SubtractedJewelFinalState.hh"
#include <string>

namespace Rivet {

  // Declaration of functions
  // Given a jetpT, returns in which interval it belongs to (0 = out of bounds)
  int pTRange(double jetpT);


  class USPJWL_PHIDIST : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(USPJWL_PHIDIST);


    int pTRange(double jetpT) {
      // Given a jetpT, returns in which interval it belongs to (0 = out of bounds)
      // ATLAS pT bins
      std::vector<double> PTBINS = {71., 79., 89., 100., 126., 158., 200., 251., 316., 398., 500., 650., 1000.};
      int pos = 0;
      for (double edge : PTBINS) {
        if (jetpT > edge) {
          pos++;
        }

        // Just need to check until fails
        else {
          return pos;
        }

      }

      // If pT is larger than every edge => out of bounds
      return 0;
    }


    void init() {

      // Jet anisotropies based on arXiv:2111.06606 (hepdata: https://www.hepdata.net/record/ins1967021)

      // Grab variable jet R parameter from environment, default value of 0.2
      RJETS = getenv("RJETS") ? getenv("RJETS") : "0.2";
      RJETS_f = std::stof(RJETS);
      std::cout << "\nR chosen for jet algorithm: " << RJETS << std::endl;

      Cut cut(Cuts::abseta < 3.2); 
  
      SubtractedJewelEvent sev(1.0); 
      SubtractedJewelFinalState fs(sev, cut); 
      declare(fs, "FS"); 

      // Apply FastJet
      FastJets fj(fs, FastJets::ANTIKT, RJETS_f);
      fj.useInvisibles();
      declare(fj, "Jets");


      // Book histograms, each for a pt bin
      book(_hist_1, "71_79_phi_R" + RJETS, 64, 0., 2 * M_PI);
      book(_hist_2, "79_89_phi_R" + RJETS, 64, 0., 2 * M_PI);
      book(_hist_3, "89_100_phi_R" + RJETS, 64, 0., 2 * M_PI);
      book(_hist_4, "100_126_phi_R" + RJETS, 64, 0., 2 * M_PI);
      book(_hist_5, "126_158_phi_R" + RJETS, 64, 0., 2 * M_PI);
      book(_hist_6, "158_200_phi_R" + RJETS, 64, 0., 2 * M_PI);
      book(_hist_7, "200_251_phi_R" + RJETS, 64, 0., 2 * M_PI);
      book(_hist_8, "251_316_phi_R" + RJETS, 64, 0., 2 * M_PI);
      book(_hist_9, "316_398_phi_R" + RJETS, 64, 0., 2 * M_PI);
      book(_hist_10, "398_500_phi_R" + RJETS, 64, 0., 2 * M_PI);
      book(_hist_11, "500_650_phi_R" + RJETS, 64, 0., 2 * M_PI);
      book(_hist_12, "650_1000_phi_R" + RJETS, 64, 0., 2 * M_PI);
    }


    /// Perform the per-event analysis
    void analyze(const Event& evt) {

      // Method definitions
      double etamax = 3.2 - RJETS_f;
      Cut jetcuts = Cuts::pT > 70 * GeV && Cuts::absrap < 1.2 && Cuts::abseta < etamax;
      const Jets jets = apply<FastJets>(evt, "Jets").jetsByPt(jetcuts);

      for (const Jet& j : jets) {
        // Jet properties
        double phi = j.phi(), pt = j.pT();


        // Fill the right histograms for each pT range
        switch (pTRange(pt)) {
          case 1:
            _hist_1 -> fill(phi);
            break;
          case 2:
            _hist_2 -> fill(phi);
            break;
          case 3:
            _hist_3 -> fill(phi);
            break;
          case 4:
            _hist_4 -> fill(phi);
            break;
          case 5:
            _hist_5 -> fill(phi);
            break;
          case 6:
            _hist_6 -> fill(phi);
            break;
          case 7:
            _hist_7 -> fill(phi);
            break;
          case 8:
            _hist_8 -> fill(phi);
            break;
          case 9:
            _hist_9 -> fill(phi);
            break;
          case 10:
            _hist_10 -> fill(phi);
            break;
          case 11:
            _hist_11 -> fill(phi);
            break;
          case 12:
            _hist_12 -> fill(phi);
            break;
          default:
             //std::cout << "pT out of bounds: " << pt << " GeV!" << std::endl;
            break;


        }

      }

    }


    void finalize() {
      // Scale only after yoda merge
    }


    /// @name Histograms
    // R_AA
    Histo1DPtr _hist_1, _hist_2, _hist_3, _hist_4, _hist_5, _hist_6, _hist_7,
    _hist_8, _hist_9, _hist_10, _hist_11, _hist_12;

    double RJETS_f;
    std::string RJETS;

  };



  DECLARE_RIVET_PLUGIN(USPJWL_PHIDIST);





}
