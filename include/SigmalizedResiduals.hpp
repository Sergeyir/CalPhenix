/** 
 *  @file   SigmalizedResiduals.hpp 
 *  @brief  Contains declarations of functions and variables that are used for estimation of values for calibration of sigmalized residuals dphi and dz
 *
 *  This file is a part of a project CalPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef SIGMALIZED_RESIDUALS_HPP
#define SIGMALIZED_RESIDUALS_HPP

#include <memory>
#include <thread>
#include <algorithm>
#include <filesystem>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMath.h"

#include "IOTools.hpp"
#include "MathTools.hpp"

#include "TF1Tools.hpp"
#include "TCanvasTools.hpp"

#include "PBar.hpp"

#include "InputYAMLReader.hpp"

/*! @namespace SigmalizedResiduals
 * @brief Contains all functions and containers for SigmalizedResiduals.cpp and CheckSigmalizedResiduals.cpp
 */
namespace SigmalizedResiduals
{
/*! @brief Calls PerformFitsForDifferentPT for the specified detector and variable for different centrality and zDC ranges
 * @param[in] detectorBin detector bin (i.e. element of array in "detectors_to_calibrate" field in input .yaml file) 
 * @param[in] variableBin variable bin (0 for dphi and 1 for dz) 
 */
   void PerformFitsForDifferentCentrAndZDC(const unsigned int detectorBin, 
                                           const unsigned int variableBin);
/*! @brief Performs approximations for the specified variable, zDC range, centrality, charge for different pT bins 
 * @param[in] hist histogram from which projections will be taken and approximated 
 * @param[in] grMeans graph in which means of approximations will be stored
 * @param[in] grSigmas graph in which sigmas of approximations will be stored
 * @param[in] detector container for the specified detector containing various data (name, approximation function, etc.; see "detectors_to_calibrate" field in input .yaml file)
 * @param[in] zDC container for the specified zDC range containing various data (minimum, maximum, etc.; see "zdc_bins" field in input .yaml file)
 * @param[in] charge charge of the particles that will be calibrated 
 * @param[in] centrality container for the specified centrality range containing various data (minimum, maximum, etc.; see "centrality_bins" field in input .yaml file)
 */
   void PerformFitsForDifferentPT(TH3F *hist, TGraphErrors& grMeans, TGraphErrors& grSigmas, 
                                  const YAML::Node& detector, 
                                  const unsigned int variableBin, const YAML::Node& zDC, 
                                  const int charge, const YAML::Node& centrality);
/*! @brief Returns yield of a signal of a distribution that can be characterised with FG+BG approximations
 * @param[in] hist histogram containing the distribution
 * @param[in] fitBG background approximation of the histogram 
 * @param[in] mean parameter of foreground (gaus) approximation 
 * @param[in] sigma parameter of foreground (gaus) approximation  
 */
   double GetYield(const TH1D* hist, const TF1& fitBG, const double mean, const double sigma);
   /// @brief Function for ProgressBar thread call
   void PBarCall();
   /// @brief If showProgress set to false in main this function will calculate how many calls have passed from .tmp files and set this value for Par::numberOfCalls
   void SetNumberOfCalls();
   /// Contents of input .yaml file for calibration
   InputYAMLReader inputYAMLCal;
   /// Contents of input .yaml file for run configuration
   InputYAMLReader inputYAMLMain;
   /// Name of run (e.g. Run14HeAu200 or Run7AuAu200)
   std::string runName;
   // Charges of particles to be analyzed independently
   const std::array<int, 2> particleCharges{1, -1};
   /// Names of variables to be calibrated
   std::array<std::string, 2> variableName{"dphi", "dz"};
   /// Names of variables to be calibrated in LaTex format
   std::array<std::string, 2> variableNameTex{"d#varphi", "dz_{DC}"};
   /// Useful objects to employ for quick TLatex insertions
   TLatex pTRangeTLatex, zDCRangeTLatex, chargeTLatex, centralityRangeTLatex;
   /// Input file (from taxi output)
   std::unique_ptr<TFile> inputFile;
   /// Output file
   std::unique_ptr<TFile> outputFile;
   /// Output directory
   std::string outputDir;
   /// Minimum pT of the whole pT range
   double pTMin;
   /// Maximum pT of the whole pT range
   double pTMax;
   /// pT ranges for ROOT TAxis
   std::vector<double> pTRanges;
   /// zDC ranges for ROOT TAxis
   std::vector<double> zDCRanges;
   /// Centrality ranges
   std::vector<double> centralityRanges;
   /// pProgress bar - shows progress (see ProgressBar)
   ProgressBar pBar{"FANCY1", "", PBarColor::BOLD_RED};
   /// Value that shows whether the computation part of this program is finished; the other part joins the threads and finishes the program
   bool isProcessFinished = false;
   /// Overall number of iterations that the program will make
   unsigned long numberOfIterations = 0;
   /// Number of calls (e.g. calls of a function in loop with numberOfIterations number of iterations)
   unsigned long numberOfCalls = 0;
   /// If true ProgressBar is printed; else the progress is written in tmp file (used if the program calls itself via shell)
   bool showProgress = true;
   /// Minimum number of entries for the histogram to be approximated. If this requirement for this value is not met warning will be printed but the program will not finish
   const double minIntegralValue = 3e2;  
   /// Number of consequent fits of dphi and dz distributions for better approximation results
   /// each consequent fit decreases the limits around value from previous fit for every parameter
   /// which makes bettter gradual gradient descent of approximation parameters since ROOT built in
   /// approximation algorithm has only limited resource to perform the gradient descent
   /// This value will be read and updated from .yaml calibration input file
   unsigned int fitNTries = 1;
   /// flag that tells the program whether dphi and dz distributions for all bins (pT, zDC, centrality, charge) should be drawn
   bool drawDValDistr = false;
   /// Mode in which the program was launched in; see main function description for more detail
   int programMode;
};
/*! @brief Main function
 *
 * Can be provided either 2 (mode 1) or 5 (mode 2) user passed input arguments (here we don't account for the name of executable as a first parameter when it is called).
 * 
 * When called in mode 1 (goes over all detectors specified in input file and values such as dphi and dz)
 * @param[in] argv[1] name of the .yaml input file or name of the directory containing .yaml input file 
 * @param[in] argv[2] number of threads the program will run on (if no value is passed this value is set to std::thread::hadrware_concurrency())
 *
 * When called in mode 2 (only for one detector and value - either dphi or dz)
 * @param[in] argv[1] name of the .yaml input file or name of the directory containing .yaml input file 
 * @param[in] argv[2] detector bin (i.e. element of array in "detectors_to_calibrate" field in input .yaml file) 
 * @param[in] argv[3] variable bin (0 for dphi and 1 for dz) 
 * @param[in] argv[4] number of threads the program will run on (if no value is passed this value is set to std::thread::hadrware_concurrency())
 * @param[in] argv[5] if it is doesn't equal to '0' the ProgressBar will show the progress (if no value is passed this value is set to true)
 */
int main(int argc, char **argv);

#endif /* SIGMALIZED_RESIDUALS_HPP */
