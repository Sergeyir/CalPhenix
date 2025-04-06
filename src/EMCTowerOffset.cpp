/** 
 *  @file   EMCTowerOffset.cpp
 *  @brief  Contains realisations of functions and variables that are used for tower offset estimation from the sum of all runs
 *
 *  This file is a part of a project CalPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef EMC_TOWER_OFFSET_CPP
#define EMC_TOWER_OFFSET_CPP

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
      errMsg += "Usage: bin/EMCTowerOffset inputFile numberOfThreads=" + 
                std::to_string(std::thread::hardware_concurrency()) + "*\n";
      errMsg += "Or**: bin/EMCTowerOffset inputFile sectorBin numberOfThreads showProgress=true";
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

   if (argc < 3) // Mode1
   {
      programMode = 1;
      if (argc > 2) numberOfThreads = std::stoi(argv[2]);
      else numberOfThreads = std::thread::hardware_concurrency();
      if (numberOfThreads == 0) CppTools::PrintError("Number of threads must be bigger than 0");

      system(("mkdir -p tmp/EMCTowerOffset/" + runName).c_str());
      system(("rm -rf tmp/EMCTowerOffset/" + runName + "/*").c_str());

      numberOfIterations = 0;
      for (unsigned int sectorBin = 0; sectorBin < 
           inputYAMLCal["sectors_to_calibrate"].size(); sectorBin++)
      {
         numberOfIterations +=
            inputYAMLCal["sectors_to_calibrate"][sectorBin]["number_of_y_towers"].as<int>()*
            inputYAMLCal["sectors_to_calibrate"][sectorBin]["number_of_z_towers"].as<int>();
      }

      int subprocessNumberOfThreads = numberOfThreads/inputYAMLCal["sectors_to_calibrate"].size();
      if (subprocessNumberOfThreads < 1) subprocessNumberOfThreads = 1;

      auto SingleThreadCall = [&](const unsigned long sectorBin)
      {
         // man ROOT sucks (TF1::Fit is still not thread safe) so I have to call the same program 
         // recursively in shell outside of the current instance to implement multithreading
         system((static_cast<std::string>("./bin/EMCTowerOffset ") + 
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

      numberOfIterations = 
         inputYAMLCal["sectors_to_calibrate"][std::stoi(argv[2])]["number_of_y_towers"].as<int>()*
         inputYAMLCal["sectors_to_calibrate"][std::stoi(argv[2])]["number_of_z_towers"].as<int>();
 
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

   const int numberOfYTowers = sector["number_of_y_towers"].as<int>();
   const int numberOfZTowers = sector["number_of_z_towers"].as<int>();

   TFile inputFile(("data/EMCTiming/" + runName + "/raw_sum.root").c_str());

   parametersOutput.open(outputDir + "CalibrationParameters/tower_offset_" + sectorName + ".txt");

   parametersOutput << numberOfYTowers << " " << numberOfZTowers << std::endl;

   for (int i = 0; i < numberOfYTowers; i++)
   {
      TH3D *distrTVsADCVsZTower = static_cast<TH3D *>
         (inputFile.Get(("traw vs ADC vs iz: " + sectorName + ", iy" + std::to_string(i)).c_str()));

      if (distrTVsADCVsZTower->GetZaxis()->GetNbins() != numberOfZTowers)
      {
         std::string errMsg = "Mismatching number of z towers from input file";
         errMsg += "and histogram for sector " + sectorName;
         CppTools::PrintError(errMsg);
      }

      for (int j = 0; j < numberOfZTowers; j++)
      {
         numberOfCalls++;

         TF1 fitFunc("t vs ADC fit", inputYAMLCal["traw_vs_adc_fit_func"].as<std::string>().c_str());

         distrTVsADCVsZTower->GetZaxis()->SetRange(j + 1, j + 1); // j + 1 to get the bin

         if (PerformFitsForSingleTower(static_cast<TH2D *>(distrTVsADCVsZTower->Project3D("xy")), 
                                       fitFunc, sectorName, i, j))
         {
            parametersOutput << 1 << " ";
            for (int k = 0; k < fitFunc.GetNpar() - 1; k++)
            {
               parametersOutput << fitFunc.GetParameter(k) << " ";
            }
            parametersOutput << fitFunc.GetParameter(fitFunc.GetNpar() - 1);
            parametersOutput << std::endl;
         }
         else
         {
            parametersOutput << 0 << std::endl;
         }

         if (!showProgress)
         {
            std::ofstream progressFile("tmp/EMCTowerOffset/" + runName + 
                                       "/" + std::to_string(sectorBin));
            progressFile << numberOfCalls;
         }
      }
   }

   parametersOutput.close();
}

bool EMCTiming::PerformFitsForSingleTower(TH2D *distr, TF1& fitFunc, const std::string& sectorName,
                                          const int yTowerIndex, const int zTowerIndex)
{
   if (distr->Integral(distr->GetXaxis()->FindBin(fitADCMin), distr->GetXaxis()->GetNbins(),
                       1, distr->GetYaxis()->GetNbins()) < 1e-15) return false;

   // distribution of means of Y projections 
   TH1D meanDistr(("mean distribution of iy" + std::to_string(yTowerIndex) + 
                   " iz" + std::to_string(zTowerIndex)).c_str(), 
                  ("iy" + std::to_string(yTowerIndex) + 
                   " iz" + std::to_string(zTowerIndex)).c_str(), 
                  distr->GetXaxis()->GetNbins(), distr->GetXaxis()->GetBinLowEdge(1),
                  distr->GetXaxis()->GetBinUpEdge(distr->GetXaxis()->GetNbins()));

   double minT = 1e31;
   double maxT = -1e31;

   for (int i = distr->GetXaxis()->FindBin(fitADCMin); i <= distr->GetXaxis()->GetNbins(); i++)
   {
      TH1D *distrProj = distr->
         ProjectionY(((std::string) distr->GetName() + "_px_" + std::to_string(i)).c_str(), i, i);

      if (distrProj->Integral(1, distrProj->GetXaxis()->GetNbins()) < 1e-15) continue;

      meanDistr.SetBinContent(i, distrProj->GetMean());
      meanDistr.SetBinError(i, distrProj->GetMeanError());

      for (int j = 1; j <= distrProj->GetXaxis()->GetNbins(); j++)
      {
         if (distrProj->GetBinContent(j) < 1e-15) continue;

         if (minT < distrProj->GetXaxis()->GetBinLowEdge(j)) break;
         else minT = distrProj->GetXaxis()->GetBinLowEdge(j);
      }

      for (int j = distrProj->GetXaxis()->GetNbins(); j >= 1; j--)
      {
         if (distrProj->GetBinContent(j) < 1e-15) continue;

         if (maxT > distrProj->GetXaxis()->GetBinUpEdge(j)) break;
         else maxT = distrProj->GetXaxis()->GetBinUpEdge(j);
      }
   }

   fitFunc.SetRange(meanDistr.GetXaxis()->GetBinLowEdge(1), 
                    meanDistr.GetXaxis()->GetBinUpEdge(meanDistr.GetXaxis()->GetNbins()));

   meanDistr.SetMinimum(minT - 5.);
   meanDistr.SetMaximum(maxT + 5.);

   // ROOT can't fit 1 point data
   if (static_cast<int>(meanDistr.GetEntries()) < 2)
   {
      fitFunc.SetParameters(minT, 0., 0.);
   }
   else
   {
      fitFunc.SetParameters(minT, 50, -1.);

      distr->GetYaxis()->SetRange(distr->GetYaxis()->FindBin(minT - 5.), 
                                  distr->GetYaxis()->FindBin(maxT + 5.));

      for (unsigned int i = 1; i <= fitNTries; i++)
      {
         meanDistr.Fit(&fitFunc, "RQMBN");

         for (int j = 0; j < fitFunc.GetNpar(); j++)
         {
            fitFunc.SetParLimits(j, fitFunc.GetParameter(j)/(1. + 2./static_cast<double>(i*i)), 
                                 fitFunc.GetParameter(j)*(1. + 2./static_cast<double>(i*i)));
         }
      }
   }

   TCanvas meanCanv("mean distr", "",  1000, 500);
   meanCanv.Divide(2);

   meanCanv.cd(1);
   distr->DrawClone("COLZ");

   meanCanv.cd(2);
   meanDistr.DrawClone();
   fitFunc.DrawClone("SAME");

   ROOTTools::PrintCanvas(&meanCanv, outputDir + sectorName + "/mean_iy" + 
                                     std::to_string(yTowerIndex) + "_iz" + 
                                     std::to_string(zTowerIndex));
   return true;
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
        std::filesystem::directory_iterator("tmp/EMCTowerOffset/" + runName))
   {
      std::string fileName = static_cast<std::string>(file.path());
      std::ifstream tmpFile(fileName.c_str());
      double currentFileNCalls = 0;
      if (tmpFile >> currentFileNCalls) numberOfCalls += currentFileNCalls;
   }
}

#endif /* EMC_TOWER_OFFSET_CPP */
