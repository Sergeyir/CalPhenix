// $HEADER$
//------------------------------------------------------------------------------------------------
//                      SigmailzedResiduals functions declarations
//------------------------------------------------------------------------------------------------
// SigmalizedResiduals
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Basic code for calibrating sigmalized residuals (sdphi, sdz)
 **/
//------------------------------------------------------------------------------------------------

#ifndef SIGMALIZED_RESIDUALS_HPP
#define SIGMALIZED_RESIDUALS_HPP

#include <memory>
#include <algorithm>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"

#include "IOTools.hpp"
#include "MathTools.hpp"

#include "TCanvasPrinter.hpp"

#include "PBar.hpp"

#include "InputReader.hpp"

struct
{
   InputJSONReader inputJSONCal, inputJSONMain;

   std::string runName;
   
   // function for the first preliminary fit of means and sigmas; 
   // it is needed to evaluate the parameters limits ranges
   // the functions listed below are quite good at this first preliminary approximation
   const std::string meansFitPrelimFunc = "[0] - [1]*exp([2]*x) + [3]*exp([4]*x)";
   const std::string sigmasFitPrelimFunc = "[0] - [1]*exp([2]*x) + [3]*exp([4]*x)";
 
   const double minIntegralValue = 3e2; // minimum number of entries for 
                                        // the histogram to be approximated
                                        // if the requirement for this value is not met
                                        // warning will be printed
   
   // number of consequent fits of dphi and dz distributions for better approximation results
   const unsigned short fitNTries = 5;

   // useful object to employ for quick TLatex insertions
   TLatex texText;
} Par;

int main(int argc, char **argv);
void PerformFits(TH3F *hist, TGraphErrors& grMeans, TGraphErrors& grSigmas, 
                 const Json::Value& calibrationInput,
                 const Json::Value& detector, const Json::Value& variable,
                 const Json::Value& zDCBin, const Json::Value& particleType,
                 const Json::Value& centralityBin);

#endif /* SIGMALIZED_RESIDUALS_HPP */
