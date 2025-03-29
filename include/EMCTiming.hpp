/** 
 *  @file   EMCTTiming.hpp
 *  @brief  Contains declarations of functions and variables that are used for estimation for timing calibration parameters for EMCal
 *
 *  This file is a part of a project CalPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef EMC_TIMING_HPP
#define EMC_TIMING_HPP

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

/*! @namespace EMCTiming
 * @brief Contains all functions, variables, and containers for EMCTowerOffset.cpp
 */
namespace EMCTiming
{
   /*! @brief Performs t vs ADC fit for a single tower. If the data was not empty and fit was succesfully performed returns true; else fasle.
    *
    * @param[in] distr histogram containing t vs ADC distribution for the single tower
    * @param[in] fitFunc function that will be used for approximation 
    * @param[in] yTowerIndex y index of the tower 
    * @param[in] zTowerIndex y index of the tower 
    */
   bool PerformFitsForSingleTower(TH2D *distr, TF1& fitFunc, const std::string& sectorName,
                                  const int yTowerIndex, const int zTowerIndex);
   /*! @brief Calls PerformFitsForSingleTower for different towers in a given sector
    *
    * @param[in] sector EMCal sector
    */
   void ProcessSector(const int sectorBin);
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
   /// Output directory
   std::string outputDir;
   /// pProgress bar - shows progress (see ProgressBar)
   ProgressBar pBar{"FANCY1", "", PBarColor::BOLD_GREEN};
   /// Value that shows whether the computation part of this program is finished; the other part joins the threads and finishes the program
   bool isProcessFinished = false;
   /// Overall number of iterations that the program will make
   unsigned long numberOfIterations = 0;
   /// Number of calls (e.g. calls of a function in loop with numberOfIterations number of iterations)
   unsigned long numberOfCalls = 0;
   /// If true ProgressBar is printed; else the progress is written in tmp file (used if the program calls itself via shell)
   bool showProgress = true;
   /// Number of consequent fits of t vs ADC distribution
   /// each consequent fit decreases the limits around value from previous fit for every parameter
   /// which makes bettter gradual gradient descent of approximation parameters since ROOT built in
   /// approximation algorithm has only limited resource to perform the gradient descent
   /// This value will be read and updated from .yaml calibration input file
   unsigned int fitNTries = 5;
   /// minimum value of ADC for the fit
   double fitADCMin = 0.;
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

#endif /* EMC_TIMING_HPP */
