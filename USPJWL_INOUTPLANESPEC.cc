// -*- C++ -*-

// This is a Rivet analysis for JEWEL  
// In- and out-of-plane charged jet spectrum based on ALICE arXiv:2307.14097 (hepdata: https://www.hepdata.net/record/ins2681682)
// It must receive the symmetry plane angles from hydro as environment variables (or assume 0 for all)
// Applications can be found in
// [Work in progress]
// 
// --Leonardo Barreto, Muenster, 2024

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/PseudoJet.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "Rivet/Projections/SubtractedJewelEvent.hh"
#include "Rivet/Projections/SubtractedJewelFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include <string>

namespace Rivet {
  
  // Converts the psi from [-pi, pi] to [0, 2pi]
  double planeConversion(double psi);

  // Returns the angular distance between phi1 and phi2
  double angDistance(double phi1, double phi2);

  // Returns the angle phi in the interval [0, 2pi]
  double angBoundary(double phi);

  // Calculates Dphi given phi and psi
  bool isInPlane(double phi, double psi, int n);

  class USPJWL_INOUTPLANESPEC : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(USPJWL_INOUTPLANESPEC);


    void init() {

      // Grab variable jet R parameter from environment, default value of 0.4
      RJETS = getenv("RJETS") ? getenv("RJETS") : "0.4";
      RJETS_f = std::stof(RJETS);
      std::cout << "\nR chosen for jet algorithm: " << RJETS << std::endl;
     
	    // Get soft symmetry planes, default value of 0
	    PSI2 = getenv("PSI2") ? std::stof(getenv("PSI2")) : 0.; 
	    PSI3 = getenv("PSI3") ? std::stof(getenv("PSI3")) : 0.; 
	    PSI4 = getenv("PSI4") ? std::stof(getenv("PSI4")) : 0.; 

	    PSI2 = planeConversion(PSI2);
      PSI3 = planeConversion(PSI3);
      PSI4 = planeConversion(PSI4);
      std::cout << "Psi angles conversion [-pi, pi] -> [0, 2pi]" << std::endl;
      std::cout << getenv("PSI2") << " -> " << PSI2 << std::endl;
      std::cout << getenv("PSI3") << " -> " << PSI3 << std::endl;
      std::cout << getenv("PSI4") << " -> " << PSI4 << std::endl;

      SubtractedJewelEvent sev(1.0);
      SubtractedJewelFinalState fs(sev, Cuts::abseta < 0.9);
	    ChargedFinalState cfs(fs);
      declare(cfs, "CFS");

      // Apply FastJet
      FastJets fj(cfs, FastJets::ANTIKT, RJETS_f);
      fj.useInvisibles();
      declare(fj, "Jets");

      // Book histograms
      _hist_inplane2 = book(_hist_inplane2, "InPlaneSpec_N2_R" + RJETS, PTEDGES);
      _hist_outplane2 = book(_hist_outplane2, "OutPlaneSpec_N2_R" + RJETS, PTEDGES);
      
      _hist_inplane3 = book(_hist_inplane3, "InPlaneSpec_N3_R" + RJETS, PTEDGES);
      _hist_outplane3 = book(_hist_outplane3, "OutPlaneSpec_N3_R" + RJETS, PTEDGES);
	  
      _hist_inplane4 = book(_hist_inplane4, "InPlaneSpec_N4_R" + RJETS, PTEDGES);
      _hist_outplane4 = book(_hist_outplane4, "OutPlaneSpec_N4_R" + RJETS, PTEDGES);
	  
	    _hist_allplane = book(_hist_allplane, "Spec_R" + RJETS, PTEDGES);
    }


    /// Perform the per-event analysis
    void analyze(const Event& evt) {

      // Method definitions
      Cut cutlead = Cuts::pT > 5 * GeV && Cuts::pT < 100 * GeV;
      double etamax = 0.9 - RJETS_f;

      // Get jets of event
      Cut jetcuts = Cuts::pT > 20 * GeV && Cuts::abseta < etamax;
      const Jets jets = apply<FastJets>(evt, "Jets").jetsByPt(jetcuts);

      for (const Jet& j : jets) {
        // Jet properties
        double pt = j.pT(), phi = j.phi();
		
		    // Check leading particle respects selection cuts
        Particles plead = j.constituents(cutlead);
        if (plead.size() == 0) continue;
	
		    // Fill histograms
		    // If not in-plane, check if it is in-plane considering out-of-plane angle
		    // i.e. if jet is out-of-plane, it is possible to be neither given
		    // out-of-plane angle is PSI_n + pi / n 
		    // ALICE definition
		    
		    // n = 2	
		    if (isInPlane(phi, PSI2, 2)) {
		    	_hist_inplane2 -> fill(pt);
		    } else if (isInPlane(phi, PSI2 + M_PI / 2, 2)) {
		    	_hist_outplane2 -> fill(pt);
		    }
		    
		    // n = 3	
		    if (isInPlane(phi, PSI3, 3)) {
		    	_hist_inplane3 -> fill(pt);
		    } else if (isInPlane(phi, PSI2 + M_PI / 3, 2)) {
		    	_hist_outplane3 -> fill(pt);
		    }
		    
		    // n = 4	
		    if (isInPlane(phi, PSI4, 4)) {
		    	_hist_inplane4 -> fill(pt);
		    } else if (isInPlane(phi, PSI2 + M_PI / 4, 2)) {
		    	_hist_outplane4 -> fill(pt);
		    }
	
		    _hist_allplane -> fill(pt);
      }
    }


    void finalize() {
      // Scale only after yoda merge
    }


    /// @name Histograms
    // R_AA
    Histo1DPtr _hist_inplane2, _hist_outplane2, _hist_inplane3, _hist_outplane3, _hist_inplane4, _hist_outplane4, _hist_allplane;

    double RJETS_f, PSI2, PSI3, PSI4;
    std::string RJETS;

    std::vector<double> PTEDGES = {20., 25., 35., 40., 50., 60., 80., 100., 120., 140., 200.};
  };



  DECLARE_RIVET_PLUGIN(USPJWL_INOUTPLANESPEC);

  
  // Auxiliary functions
  double planeConversion(double psi) {
    // Converts the psi to [0, 2pi]
	  return fmod(psi + 2 * M_PI, 2 * M_PI);
  }


  double angDistance(double phi1, double phi2) {
    // Returns the angular distance between phi1 and phi2
	  double phi1_conv = planeConversion(phi1), phi2_conv = planeConversion(phi2);
    double diff = std::abs(phi1_conv - phi2_conv);

	  return diff > M_PI ? 2 * M_PI - diff : diff;
  }


  bool isInPlane(double phi, double psi, int n) {
  	// Calculates if phi is in-plane, given harmonic n
  	// This is done by finding the minimum distance between phi
  	// and all symmentry angles of psi (mindist) and comparing to maximum
  	// distance of in-plane angle (maxinplanedist) 
    std::vector<double> dists = {};

  	// Usually, maxinpladist for n = 2 is pi / 4. 
  	// ALICE used pi / 6 for a better contrast between in- and out-of-plane yields
  	// Thus a 2/3 factor is added in the generalized formula
  	double maxinplanedist = (2. / 3.) * M_PI / (2 * n);
  
  	for (int i = 0; i < n - 1; i++) {
        dists.push_back(angDistance(phi, psi + 2 * M_PI * i / n));
  	}
  
  	double mindist = *std::min_element(std::begin(dists), std::end(dists));
  	return mindist < maxinplanedist;
  }

}
