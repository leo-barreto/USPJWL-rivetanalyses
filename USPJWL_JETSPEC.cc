// -*- C++ -*-

// This is a Rivet analysis for JEWEL  
// Jet spectrum/RAA based on ATLAS arxiv:1805.05635 (hepdata: https://www.hepdata.net/record/ins1673184)
// xJ based on ATLAS arXiv:2205.00682 (hepdata: missing?) 
// Applications can be found in
// arXiv:2208.02061 and https://doi.org/10.11606/D.43.2021.tde-05112021-191914
//
// --Leonardo Barreto, IF USP, 2021

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/PseudoJet.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/JetShape.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "Rivet/Projections/SubtractedJewelEvent.hh"
#include "Rivet/Projections/SubtractedJewelFinalState.hh"
#include <string>

namespace Rivet {

  // Declaration of functions

  // Given a abs rapidity, returns in which interval it belongs to (0 = out of bounds)
  int absrapRange(double jety);

  // Given a jetpT, returns in which interval it belongs to (0 = out of bounds)
  int pTRange(double jetpT);


  class RAA_ATLAS : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(RAA_ATLAS);

      // Necessary functions

      int absrapRange(double jety) {
        // Given a absrap, returns in which interval it belongs to (0 = out of bounds)

        // ATLAS absolute rapidity bin edges
        std::vector<double> ABSRAPEDGES = {0., 0.3, 0.8, 1.2, 1.6, 2.1, 2.8};

        int pos = 0;
        for (double edge : ABSRAPEDGES) {
          if (jety > edge) {
            pos++;
          }

          // Just need to check until fails
          else {
            return pos;
          }

        }

        // If absrap is larger than every edge => out of bounds
        return 0;
      }



      int pTRange(double jetpT) {
        // Given a jetpT, returns in which interval it belongs to (0 = out of bounds)

        // ATLAS pT bins + super lower testing
        std::vector<double> PTEDGES = {10., 30., 60., 90., 100., 112., 126., 141.,
                                       158., 178., 200., 224., 251., 282., 316.,
                                       398., 562., 630., 1000.};
        int pos = 0;
        for (double edge : PTEDGES) {
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


      // Jet spectrum/RAA based on ATLAS arxiv:1805.05635 (hepdata: https://www.hepdata.net/record/ins1673184)
      // xJ based on ATLAS arXiv:2205.00682 (hepdata: missing?)      

      // Grab variable jet R parameter from environment, default value of 0.4
      RJETS = getenv("RJETS") ? getenv("RJETS") : "0.4";
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



      // Book histograms

      // For R_AA:
      // Name convention: _hist_[rapidity range index], except for inclusive

      // absrap bins: 0–0.3, 0.3–0.8, 0.8–1.2, 1.2–1.6, 1.6–2.1, 2.1–2.8
      // inclusive: 0-2.1, 0-2.8
      book(_hist_1,"JetpT_0_0.3_R" + RJETS, PTEDGES);
      book(_hist_2,"JetpT_0.3_0.8_R" + RJETS, PTEDGES);
      book(_hist_3,"JetpT_0.8_1.2_R" + RJETS, PTEDGES);
      book(_hist_4,"JetpT_1.2_1.6_R" + RJETS, PTEDGES);
      book(_hist_5,"JetpT_1.6_2.1_R" + RJETS, PTEDGES);
      book(_hist_6,"JetpT_2.1_2.8_R" + RJETS, PTEDGES);
      book(_hist_7,"JetpT_0_2.1_R" + RJETS, PTEDGES);
      book(_hist_8,"JetpT_0_2.8_R" + RJETS, PTEDGES);
      book(_hist_9,"JetpT_0_1.2_R" + RJETS, PTEDGES);
      book(_hist_10,"JetpT_R" + RJETS, PTEDGES);

      // For x_J:
      // Name convention: _xj_[pT range index]

      // leading jet pt binning is: 158-178, 178-200, 200-224, 224-251, 251-282,
      // 282-316, 316-398, 398-562, 562-700, 700-1000
      // LOW PT EDGES = {10., 30., 60., 90., 120., 158.}
      book(_xj_1,"xJ_10_30_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_2,"xJ_30_60_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_3,"xJ_60_90_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_4,"xJ_90_100_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_5,"xJ_100_112_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_6,"xJ_112_126_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_7,"xJ_126_141_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_8,"xJ_141_158_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_9,"xJ_158_178_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_10,"xJ_178_200_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_11,"xJ_200_224_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_12,"xJ_224_251_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_13,"xJ_251_282_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_14,"xJ_282_316_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_15,"xJ_316_398_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_16,"xJ_398_562_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_17,"xJ_562_630_R" + RJETS, 20, 0.32, 1.0);
      book(_xj_18,"xJ_630_1000_R" + RJETS, 20, 0.32, 1.0);

      // For R_AA^Lead and R_AA^Sublead
      book(_lead,"JetpT1_R" + RJETS, PTEDGES_J);
      book(_sublead,"JetpT2_R" + RJETS, PTEDGES_J);
      book(_counter,"xJ_counter_R" + RJETS, 2., -0.5, 1.5);


    }


    /// Perform the per-event analysis
    void analyze(const Event& evt) {

      // Get jets of event
      double etamax = 3.2 - RJETS_f;
      Cut jetcuts = Cuts::pT > 20 * GeV && Cuts::abseta < etamax;

      const Jets& jets = apply<FastJets>(evt, "Jets").jetsByPt(jetcuts);
 
      // CALCULATE JET PT FOR RAA
      for (const Jet& j : jets) {
        // Jet properties
        double y = j.absrap(), pt = j.pT();

        // Fill the right histograms for each pT range
        switch (absrapRange(y)) {
          case 1:
            _hist_1 -> fill(pt);
            break;

          case 2:
            _hist_2 -> fill(pt);
            break;

          case 3:
            _hist_3 -> fill(pt);
            break;

          case 4:
            _hist_4 -> fill(pt);
            break;

          case 5:
            _hist_5 -> fill(pt);
            break;

          case 6:
            _hist_6 -> fill(pt);
            break;

          // If something weird happens, signalize and ignore from analysis
          default:
            //std::cout << "|y| out of bounds: " << y << std::endl;
            break;
        }

        // Fill inclusive histograms
        if (y <= 2.1) {
          _hist_7 -> fill(pt);
        }

        if (y <= 2.8) {
          _hist_8 -> fill(pt);
        }

        if (y <= 1.2) {
          _hist_9 -> fill(pt);
        }

        _hist_10 -> fill(pt);

      }


      // CALCULATE JET PT LEADING AND SUBLEADING FOR XJ
      // Apply new cuts (with selectors) before sorting -> no need anymore, as jets are Jets not PseudoJets (new subtraction method)
      
      Cut xjcuts = Cuts::abseta < 2.1 && Cuts::pT > 20 * GeV;
      vector<Jet> sorted_jets = select(jets, xjcuts);

      if (sorted_jets.size() >= 2) { // Need at least two jets
        double pTLead = sorted_jets[0].pt(), pTSubLead = sorted_jets[1].pt();
        double Dphi = deltaPhi(sorted_jets[0].phi(), sorted_jets[1].phi());


        // Two conditions must be satisfied: both |eta| < 2.1 and Dphi > 7pi / 8
        // We add pTSubLead > 20 GeV to eliminate weird events
        if (Dphi > 7 * M_PI / 8) {
          // Add to the counter if the event pass the criteria
          _counter -> fill(1.);
          _lead -> fill(pTLead);
          _sublead -> fill(pTSubLead);

          double xj = pTSubLead / pTLead;


          switch (pTRange(pTLead)) {
            case 1:
              _xj_1 -> fill(xj);
              break;

            case 2:
              _xj_2 -> fill(xj);
              break;

            case 3:
              _xj_3 -> fill(xj);
              break;

            case 4:
              _xj_4 -> fill(xj);
              break;

            case 5:
              _xj_5 -> fill(xj);
              break;

            case 6:
              _xj_6 -> fill(xj);
              break;

            case 7:
              _xj_7 -> fill(xj);
              break;

            case 8:
              _xj_8 -> fill(xj);
              break;

            case 9:
              _xj_9 -> fill(xj);
              break;

            case 10:
              _xj_10 -> fill(xj);
              break;

            case 11:
              _xj_11 -> fill(xj);
              break;

            case 12:
              _xj_12 -> fill(xj);
              break;

            case 13:
              _xj_13 -> fill(xj);
              break;

            case 14:
              _xj_14 -> fill(xj);
              break;

            case 15:
              _xj_15 -> fill(xj);
              break;

            case 16:
              _xj_16 -> fill(xj);
              break;

            case 17:
              _xj_17 -> fill(xj);
              break;

            case 18:
              _xj_18 -> fill(xj);
              break;

            default:
              break;
          }
        }

        else {
          _counter -> fill(0.);
        }
      }
    }


    void finalize() {
      // Scale only after yoda merge
    }


    /// @name Histograms
    // R_AA
    Histo1DPtr _hist_1, _hist_2, _hist_3, _hist_4, _hist_5, _hist_6, _hist_7,
               _hist_8, _hist_9, _hist_10;

    // x_J
    Histo1DPtr _xj_1, _xj_2, _xj_3, _xj_4, _xj_5, _xj_6, _xj_7, _xj_8, _xj_9,
               _xj_10, _xj_11, _xj_12, _xj_13, _xj_14, _xj_15, _xj_16, _xj_17,
               _xj_18;

    // J_AA
    Histo1DPtr _lead, _sublead, _counter;


    double RJETS_f;
    std::string RJETS;
    std::vector<double> PTEDGES = {30, 40, 50, 56, 63, 70, 79, 89, 100, 112,
                                   125, 141, 158, 177, 199, 223, 251, 281,
                                   316, 354, 398, 501, 630, 1000};

    std::vector<double> PTEDGES_J = {100, 112, 126, 141, 158, 178, 200, 224,
                                     251, 282, 316, 398, 562, 630, 1000};

    

  };



  DECLARE_RIVET_PLUGIN(RAA_ATLAS);



}
