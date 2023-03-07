// -*- C++ -*-

// This is a Rivet analysis for JEWEL  
// Jet spectrum/RAA based on CMS arXiv:2102.13080 (hepdata: https://www.hepdata.net/record/ins1848440)
// and ALICE arXiv:1909.09718 (hepdata: https://www.hepdata.net/record/ins1755387)
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

  class EXTRA : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(EXTRA);


    void init() {

      // Grab variable jet R parameter from environment, default value of 0.4
      RJETS = getenv("RJETS") ? getenv("RJETS") : "0.4";
      RJETS_f = std::stof(RJETS);
      std::cout << "\nR chosen for jet algorithm: " << RJETS << std::endl;

      SubtractedJewelEvent sev(1.0);
      SubtractedJewelFinalState fs(sev, Cuts::abseta < 3.2);
      declare(fs, "FS");

      // Apply FastJet
      FastJets fj(fs, FastJets::ANTIKT, RJETS_f);
      fj.useInvisibles();
      declare(fj, "Jets");



      // Book histograms

      // For R_AA:
      _hist_jet = book(_hist_jet, "JetpT_R" + RJETS, PTEDGES);
      _hist_alice = book(_hist_alice, "ALICEpT_R" + RJETS, PTEDGES_ALICE);
      _hist_alice2 = book(_hist_alice2, "ALICEpT_nolead_R" + RJETS, PTEDGES_ALICE);
      _hist_cms = book(_hist_cms, "CMSpT_R" + RJETS, PTEDGES_CMS);
    }


    /// Perform the per-event analysis
    void analyze(const Event& evt) {

      // Method definitions
      Cut cutlead = Cuts::pT > (10 * RJETS_f + 3) * GeV;
      // Vector that will store size of array of particles that
      // passes cutlead
      std::vector<int> sizelead = {};
      double etamax = 3.2 - RJETS_f;
      double etaspace;
      if (RJETS_f <= 0.4) {
        etaspace = 0.7 - RJETS_f;
      }

      else {  // Only consider ALICE's eta space for reasonable R
        etaspace = 3.2 - RJETS_f;
      }


      // Get jets of event
      Cut jetcuts = Cuts::pT > 40 * GeV && Cuts::abseta < etamax;
      const Jets jets = apply<FastJets>(evt, "Jets").jetsByPt(jetcuts);

      // Need to loop through jets before substraction to access constituents
      for (const Jet& j : jets) {
        Particles plead = j.constituents(cutlead);
        sizelead.push_back(plead.size());
      }

      // CALCULATE JET PT FOR RAA
      // Counter to associate new jet (after subtraction) to
      // element in arraysizelead
      int counter_jets = 0;
      for (const Jet& j : jets) {
        // Jet properties
        double y = j.absrap(), pt = j.pT(), eta = j.abseta();

        if (y <= 1.2) {
          _hist_jet -> fill(pt);
        }

	if (eta <= 2) {
	  _hist_cms -> fill(pt);
        }

        if (eta <= etaspace) {
          _hist_alice2 -> fill(pt); // ALICE no lead method

          if (sizelead[counter_jets] > 0) {
            _hist_alice -> fill(pt);
          }
        }

        counter_jets++;
      }
    }


    void finalize() {
      // Scale only after yoda merge
    }


    /// @name Histograms
    // R_AA
    Histo1DPtr _hist_jet, _hist_alice, _hist_alice2, _hist_cms;

    double RJETS_f;
    std::string RJETS;
    std::vector<double> PTEDGES = {71., 79., 89., 100., 126., 158., 200., 251.,
                                   316., 398., 500., 650., 1000.};

    std::vector<double> PTEDGES_ALICE = {40, 50, 60, 70, 80, 100, 120, 140};

    std::vector<double> PTEDGES_CMS = {200, 250, 300, 400, 500, 1000};


  };



  DECLARE_RIVET_PLUGIN(EXTRA);

}
