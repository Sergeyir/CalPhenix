/** 
 *  @file   EMCTRunByRunOffset.cpp
 *  @brief  Contains realistations of functions and variables that are used to check EMCal timing calibrations
 *
 *  This file is a part of a project CalPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef CHECK_EMC_TIMING_CPP
#define CHECK_EMC_TIMING_CPP

#include "../include/EMCTiming.hpp"

// This program works in 2 modes
// Mode2 analyzes track residual variable for the specific detector, variable, and charge
// Mode1 recursively calls the program outside the current process to analyze all configurations
int main(int argc, char **argv)
{
   using namespace EMCTiming;

   if (argc < 2 || argc > 5) 
   {
      std::string errMsg = "Expected 1-2 or 3-4 parameters while " + std::to_string(argc - 1) + 
                           " parameter(s) were provided \n";
      errMsg += "Usage: bin/EMCTRunByRunOffset inputFile numberOfThreads=" + 
                std::to_string(std::thread::hardware_concurrency()) + "*\n";
      errMsg += "Or**: bin/EMCTRunByRunOffset inputFile sectorBin numberOfThreads showProgress=true";
      errMsg += "*: default argument is the number of threads on the current machine \n";
      errMsg += "**: this mode processes only one sector \n";
      CppTools::PrintError(errMsg);
   }

   unsigned int numberOfThreads;

   // initializing ROOT parameters
   ROOT::EnableThreadSafety();
   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   // initializing this program parameters
   inputYAMLCal.OpenFile(argv[1], "emc_timing");
   inputYAMLCal.CheckStatus("emc_timing");

   runName = inputYAMLCal["run_name"].as<std::string>();

   CppTools::CheckInputFile("data/EMCTiming/" + runName + "/raw_sum.root");

   // opening input file with parameters of a run
   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   if (inputYAMLCal["sectors_to_calibrate"].size() == 0)
   {
      CppTools::PrintInfo("No sectors were specified for calibrations");
      CppTools::PrintInfo("Exiting the program");
      exit(1);
   }

   TDirectory::AddDirectory(kFALSE);

   const std::string inputDir = "data/EMCTiming/" + runName + "/";

   std::set<std::filesystem::path> inputDataFilePaths;

   for (const auto &file : std::filesystem::directory_iterator(inputDir))
   {
      const std::string fileName = static_cast<std::string>(file.path());
      if (fileName.substr(inputDir.size(), 3) == static_cast<std::string>("se-"))
      {
         inputDataFilePaths.insert(file);
      }
   }

   // iterating over set inputDataFiles to get run numbers sorted
   for (const auto &file : inputDataFilePaths)
   {
      const std::string fileName = static_cast<std::string>(file);
      runNumbers.emplace_back(std::stoi(fileName.substr(inputDir.size() + 3, 6)));
   }

   if (argc < 3) // Mode1
   {
      programMode = 1;
      if (argc > 2) numberOfThreads = std::stoi(argv[2]);
      else numberOfThreads = std::thread::hardware_concurrency();
      if (numberOfThreads == 0) CppTools::PrintError("Number of threads must be bigger than 0");

      system(("mkdir -p tmp/progress/EMCTRunByRunOffset/" + runName).c_str());
      system(("rm -rf tmp/progress/EMCTRunByRunOffset/" + runName + "/*").c_str());

      numberOfIterations = inputYAMLCal["sectors_to_calibrate"].size()*runNumbers.size();

      int subprocessNumberOfThreads = numberOfThreads/inputYAMLCal["sectors_to_calibrate"].size();
      if (subprocessNumberOfThreads < 1) subprocessNumberOfThreads = 1;

      auto SingleThreadCall = [&](const unsigned long sectorBin)
      {
         // man ROOT sucks (TF1::Fit is still not thread safe) so I have to call the same program 
         // recursively in shell outside of the current instance to implement multithreading
         system((static_cast<std::string>("./bin/EMCTRunByRunOffset ") + 
                 argv[1] + " " + std::to_string(sectorBin) + " " + 
                 std::to_string(subprocessNumberOfThreads) + " 0").c_str());;
      };

      std::vector<std::thread> thrCalls;
      std::thread pBarThr(PBarCall); 

      for (unsigned int sectorBin = 0; sectorBin < 
           inputYAMLCal["sectors_to_calibrate"].size(); sectorBin++)
      {
            if (thrCalls.size() >= numberOfThreads)
            {
               while (!thrCalls.empty())
               {
                  thrCalls.back().join();
                  thrCalls.pop_back();
               }
            }
            thrCalls.emplace_back(SingleThreadCall, sectorBin);
      }

      while (!thrCalls.empty())
      {
         thrCalls.back().join();
         thrCalls.pop_back();
      }

      isProcessFinished = true;
      pBarThr.join();
   }
   else // Mode2
   {
      programMode = 2;
      numberOfThreads = std::stoi(argv[3]);
      if (numberOfThreads <= 0) CppTools::PrintError("Number of threads must be bigger than 0");
 
      ROOT::EnableImplicitMT(numberOfThreads);

      if (argc > 4) showProgress = static_cast<bool>(std::stoi(argv[4]));

      outputDir = "output/EMCTCalibration/" + runName + "/";
      system(("mkdir -p " + outputDir + "CalibrationParameters").c_str());
      system(("mkdir -p " + outputDir + inputYAMLCal["sectors_to_calibrate"][std::stoi(argv[2])]
                                                    ["name"].as<std::string>()).c_str());

      fitNTries = inputYAMLCal["number_of_fit_tries"].as<unsigned int>();
      fitADCMin = inputYAMLCal["fit_adc_min"].as<double>();

      std::thread pBarThr(PBarCall); 

      numberOfIterations = runNumbers.size();
 
      ProcessSector(std::stoi(argv[2]));

      isProcessFinished = true;
      pBarThr.join();
   }
   return 0;
}

void EMCTiming::ProcessSector(const int sectorBin)
{
   const YAML::Node sector = inputYAMLCal["sectors_to_calibrate"][sectorBin];
   const std::string sectorName = sector["name"].as<std::string>();

   const std::string tPhotonFitFunc = 
      inputYAMLCal["t_photon_fit_func"].as<std::string>();
   const std::string tPhotonMeanVsADCFitFunc = 
      inputYAMLCal["tcorr_mean_vs_adc_fit_func"].as<std::string>();

   TFile outputFile((outputDir + sectorName + "/tcorr_fits.root").c_str(), "RECREATE");

   for (const int runNumber : runNumbers)
   {
      numberOfCalls++;

      TFile inputFile(("data/EMCTiming/" + runName + "/se-" + 
                       std::to_string(runNumber) + ".root").c_str());

      TH2D *tVsADC = static_cast<TH2D *>(inputFile.Get(("tcorr vs ADC: " + sectorName).c_str()));

      outputFile.mkdir(std::to_string(runNumber).c_str());
      outputFile.cd(std::to_string(runNumber).c_str());

      TGraphErrors meansTVsADC;
      TGraphErrors sigmasTVsADC;

      for (const auto& rangeADC : sector["adc_ranges"])
      {
         const double rangeADCMin = rangeADC["min"].as<double>();
         const double rangeADCMax = rangeADC["max"].as<double>();

         // used to mark the first bin in the adc range for 
         // cases when there is not enough statistics to merge bins
         int firstValidBinInRange = CppTools::Maximum(1, tVsADC->GetXaxis()->FindBin(rangeADCMin));

         for (int i = CppTools::Maximum(1, tVsADC->GetXaxis()->FindBin(rangeADCMin)); 
              i <= CppTools::Minimum(tVsADC->GetXaxis()->FindBin(rangeADCMax), 
                                    tVsADC->GetXaxis()->GetNbins()); i++)
         {
            if (tVsADC->Integral(firstValidBinInRange, i, tVsADC->GetYaxis()->FindBin(-10.), 
                                 tVsADC->GetYaxis()->FindBin(10.)) < 1000.) continue;

            // ADC value of the (current bin)/(merged bins)
            const double valADC = 
               CppTools::Average(tVsADC->GetXaxis()->GetBinCenter(i), 
                                 tVsADC->GetXaxis()->GetBinCenter(firstValidBinInRange));

            TH1D *tVsADCProj = 
               tVsADC->ProjectionY(("tcorr " + CppTools::DtoStr(valADC, 0)).c_str(), 
                                   firstValidBinInRange, i);

            // if statistics is sufficient first valid bin resets 
            // to the current bin on next iteration bin
            firstValidBinInRange = i + 1;

            TF1 tPhotonFit(("tcorr fit " + CppTools::DtoStr(valADC, 0)).c_str(), tPhotonFitFunc.c_str());

            tPhotonFit.SetRange(-10., 10.);
            tPhotonFit.SetParameters(tVsADCProj->GetMaximum(), 0, 0.5, 1., 1.);

            for (unsigned int j = 1; j <= fitNTries; j++)
            {
               if (j < fitNTries) tVsADCProj->Fit(&tPhotonFit, "RQMBN");
               else tVsADCProj->Fit(&tPhotonFit, "RQMB");

               const double parameterDeviationScale = 1. + 1./static_cast<double>(j*j);

               tPhotonFit.SetRange(tPhotonFit.GetParameter(1) - 
                                 fabs(tPhotonFit.GetParameter(2))*parameterDeviationScale,
                                 tPhotonFit.GetParameter(1) - 
                                 fabs(tPhotonFit.GetParameter(2))*parameterDeviationScale);

               for (int k = 0; k < tPhotonFit.GetNpar(); k++)
               {
                  if (k == 1)
                  {
                     tPhotonFit.SetParLimits(k, tPhotonFit.GetParameter(1) - 
                                           fabs(tPhotonFit.GetParameter(2))*
                                           (parameterDeviationScale - 1.),
                                           tPhotonFit.GetParameter(1) +
                                           fabs(tPhotonFit.GetParameter(2))*
                                           (parameterDeviationScale - 1.));
                  }
                  else
                  {
                     tPhotonFit.SetParLimits(k, tPhotonFit.GetParameter(k)/parameterDeviationScale,
                                           tPhotonFit.GetParameter(k)*parameterDeviationScale);
                  }
               }
            }

            tVsADCProj->Write();

            // skipping outliers
            if (fabs(tPhotonFit.GetParameter(1)) > 5. || 
                fabs(tPhotonFit.GetParameter(2)) > 3. || 
                fabs(tPhotonFit.GetParameter(2)) < 0.1) continue;

            meansTVsADC.AddPoint(valADC, tPhotonFit.GetParameter(1));
            sigmasTVsADC.AddPoint(valADC, fabs(tPhotonFit.GetParameter(2)));
         }
      }

      meansTVsADC.SetMarkerStyle(20);
      meansTVsADC.SetMarkerColor(kRed - 3);
      meansTVsADC.SetMarkerSize(0.5);

      sigmasTVsADC.SetMarkerStyle(20);
      sigmasTVsADC.SetMarkerColor(kAzure - 3);
      sigmasTVsADC.SetMarkerSize(0.5);

      TF1 tPhotonMeanVsADCFit("tcorr mean vs ADC fit", tPhotonMeanVsADCFitFunc.c_str(), 0., 10000);
      tPhotonMeanVsADCFit.SetLineWidth(3);
      tPhotonMeanVsADCFit.SetLineStyle(2);
      tPhotonMeanVsADCFit.SetLineColor(kBlack);

      if (meansTVsADC.GetN() == 0) // bad run that passed the first bad run check
      {
         continue;
      }
      else if (meansTVsADC.GetN() > 1) 
      {
         for (unsigned int i = 1; i < fitNTries; i++)
         {
            meansTVsADC.Fit(&tPhotonMeanVsADCFit, "RQMBN");

            for (int j = 0; j < tPhotonMeanVsADCFit.GetNpar(); j++)
            {
               tPhotonMeanVsADCFit.SetParLimits(j, tPhotonMeanVsADCFit.GetParameter(j)*
                                              (1. - 6./static_cast<double>(i*i*i)), 
                                              tPhotonMeanVsADCFit.GetParameter(j)*
                                              (1. + 4./static_cast<double>(i*i*i)));
            }
         }
      }
      else // root can't fit 1 point data
      {
         tPhotonMeanVsADCFit.SetParameter(0, meansTVsADC.GetPointY(0));
         for (int i = 1; i < tPhotonMeanVsADCFit.GetNpar(); i++)
         {
            tPhotonMeanVsADCFit.SetParameter(i, 0.);
         }
      }

      TCanvas parCanv("mean and sigma t parameters vs ADC", "", 600, 600);

      TH1 *frame = gPad->DrawFrame(meansTVsADC.GetPointX(0)/1.1, 
                                   CppTools::Minimum(TMath::MinElement(meansTVsADC.GetN(), 
                                                                       meansTVsADC.GetY()),
                                                     TMath::MinElement(sigmasTVsADC.GetN(), 
                                                                       sigmasTVsADC.GetY())) - 0.5,
                                   meansTVsADC.GetPointX(meansTVsADC.GetN() - 1)*1.1,
                                   CppTools::Maximum(TMath::MaxElement(meansTVsADC.GetN(), 
                                                                       meansTVsADC.GetY()),
                                                     TMath::MaxElement(sigmasTVsADC.GetN(), 
                                                                       sigmasTVsADC.GetY())) + 0.5);

      frame->GetXaxis()->SetTitle("ADC");

      tPhotonMeanVsADCFit.Draw("SAME");
      meansTVsADC.Draw("P");
      sigmasTVsADC.Draw("P");

      ROOTTools::PrintCanvas(&parCanv, "output/EMCTCalibration/" + runName + "/" + 
                             sectorName + "/tcorr_par_vs_adc_" + std::to_string(runNumber));

      
      frame->GetXaxis()->SetTitle("ADC");

      if (!showProgress)
      {
         std::ofstream progressFile("tmp/progress/EMCTRunByRunOffset/" + runName + 
                                    "/" + std::to_string(sectorBin));
         progressFile << numberOfCalls;
      }
   }

   outputFile.Close();
}

void EMCTiming::PBarCall()
{
   if (!showProgress) return;
   while (!isProcessFinished)
   {
      SetNumberOfCalls();
      pBar.Print(static_cast<double>(numberOfCalls)/
                 static_cast<double>(numberOfIterations));
      std::this_thread::sleep_for(std::chrono::milliseconds(200));
   }
   pBar.Print(1.);
};

void EMCTiming::SetNumberOfCalls()
{
   if (programMode != 1) return; // Only Mode1 passes
   numberOfCalls = 0;
   for (const auto &file : 
        std::filesystem::directory_iterator("tmp/progress/EMCTRunByRunOffset/" + runName))
   {
      std::string fileName = static_cast<std::string>(file.path());
      std::ifstream tmpFile(fileName.c_str());
      double currentFileNCalls = 0;
      if (tmpFile >> currentFileNCalls) numberOfCalls += currentFileNCalls;
   }
}

#endif /* CHECK_EMC_TIMING_CPP */
