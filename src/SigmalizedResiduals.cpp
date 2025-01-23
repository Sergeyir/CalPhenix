// $SOURCE$
//------------------------------------------------------------------------------------------------
//                      SigmailzedResiduals functions realisations
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

#ifndef SIGMALIZED_RESIDUALS_CPP
#define SIGMALIZED_RESIDUALS_CPP

#include "../include/SigmalizedResiduals.hpp"

int main(int argc, char **argv)
{
   if (argc < 2 || argc > 6 || (argc > 3 && argc < 5)) 
   {
      std::string errMsg = "Expected 1-2 or 5-6 parameters while " + std::to_string(argc - 1) + 
                           " parameter(s) were provided \n";
      errMsg += "Usage: bin/SigmalizedResiduals inputFile numberOfThreads=" + 
                std::to_string(std::thread::hardware_concurrency()) + "*\n";
      errMsg += "Or**: bin/SigmalizedResiduals inputFile detectorBin variableBin particleBin ";
      errMsg += "numberOfThreads=" + std::to_string(std::thread::hardware_concurrency()) + "*\n";
      errMsg += "*: default argument is the number of threads on the current machine \n";
      errMsg += "**: this method is not for user and is used for Multithreading"
      PrintError(errMsg);
   }

   unsigned int numberOfThreads;

   // initializing ROOT parameters
   ROOT::EnableImplicitMT(numberOfThreads);
   ROOT::EnableThreadSafety();
   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   // initializing this program parameters
   Par.inputJSONCal.OpenFile(argv[1], "sigmalized_residuals");
   Par.inputJSONCal.CheckStatus("sigmalized_residuals");
   
   Par.runName = Par.inputJSONCal["run_name"].asString();
   
   // opening input file with parameters of a run
   Par.inputJSONMain.OpenFile("input/" + Par.runName + "/main.json");
   Par.inputJSONMain.CheckStatus("main");

   const Json::Value calibrationInput = Par.inputJSONCal["SIGMALIZED_RESIDUALS"];
   
   if (calibrationInput["detectors_to_calibrate"].size() == 0)
   {
      PrintInfo("No detectors are specified for calibrations");
      PrintInfo("Exiting the program");
      exit(1);
   }
   
   if (argc < 4)
   {
      if (argc > 2) numberOfThreads = std::stoi(argv[2]);
      else numberOfThreads = std::thread::hardware_concurrency();
      if (numberOfThreads == 0) PrintError("Number of threads must be bigger than 0");

      Par.numberOfIterations = calibrationInput["detectors_to_calibrate"].size()*
                               calibrationInput["centrality_bins"].size()*
                               calibrationInput["variables_to_calibrate"].size()*
                               calibrationInput["particles_to_calibrate"].size()*
                               calibrationInput["zdc_bins"].size();

      for (const auto& detector : calibrationInput["detectors_to_calibrate"])
      {
         system(("mkdir -p output/SigmalizedResiduals/" + Par.runName + 
                 "/" + detector["name"].asString()).c_str());
      }
      
      bool isProcessFinished = false;
      
      auto PBarCall = [&]()
      {
         ProgressBar pBar("FANCY1", "", PBarColor::BOLD_GREEN);
         while (!isProcessFinished)
         {
            pBar.Print(static_cast<double>(Par.numberOfCalls)/
                           static_cast<double>(Par.numberOfIterations));
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
         }
         pBar.Print(1.);
      };
      
      auto SingleThreadCall = [&](const unsigned long detectorBin, 
                                  const unsigned long variableBin,
                                  const unsigned long particleBin)
      {
         // man ROOT sucks (TF1::Fit is still not thread safe) so I have to call this program 
         // recursively outside of the current instance to implement multithreading
         system(("./bin/SigmalizedResiduals " + argv[1] + 
                 std::to_string(detectorBin) + " " + std::to_string(variableBin) + " " +
                 std::to_string(particleBin) + " " + std::to_string(numberOfThreads)).c_str());
      }

      std::vector<std::thread> thrCalls;
      std::thread pBarThr(PBarCall); 

      for (unsigned long detectorBin = 0; detectorBin < 
           calibrationInput["detectors_to_calibrate"].size(); detectorBin++)
      {
         for (unsigned long variableBin = 0; variableBin < 
              calibrationInput["variables_to_calibrate"].size(); variableBin++)
         { 
            for (unsigned long particleBin = 0; particleBin < 
                 calibrationInput["particles_to_calibrate"].size(); particleBin++)
            {
               if (thrCalls.size() >= numberOfThreads)
               {
                  thrCalls.back().join();
                  thrCalls.pop_back();
               }
               thrCalls.emplace_back(SingleThreadCall, detectorBin, variableBin, particleBin);
               //PerformFitsForDifferentCentrAndZDC(detectorBin, variableBin, particleBin);
            }
         }
      }

      while (!thrCalls.empty())
      {
         thrCalls.back().join();
         thrCalls.pop_back();
      }
      
      isProcessFinished = true;
      pBarThr.join();
   }
   else
   {
      if (argc > 5) numberOfThreads = std::stoi(argv[6]);
      else numberOfThreads = std::thread::hardware_concurrency();
      if (numberOfThreads == 0) PrintError("Number of threads must be bigger than 0");

      Par.inputFile = std::unique_ptr<TFile>(TFile::Open(("data/Real/" + Par.runName + 
                                                          "/SingleTrack/sum.root").c_str(), "READ"));
      
      Par.texText.SetTextFont(52);
      Par.texText.SetTextSize(0.06);
      Par.texText.SetNDC();

      for (const Json::Value& pTBin: calibrationInput["pt_bins"])
      {
         Par.pTRanges.push_back(pTBin["min"].asDouble());
      }
      Par.pTRanges.push_back(calibrationInput["pt_bins"].back()["max"].asDouble());

      for (const Json::Value& centrality: calibrationInput["centrality_bins"])
      {
         Par.centralityRanges.push_back(centrality["min"].asDouble());
      }
      Par.centralityRanges.push_back(calibrationInput["centrality_bins"].back()["max"].asDouble());

      Par.outputDir = "output/SigmalizedResiduals/" + Par.runName + "/";
      Par.outputFile = std::unique_ptr<TFile>(TFile::Open((Par.outputDir + "all_fits.root").c_str(), 
                                                          "RECREATE"));
      
      Par.pTMin = calibrationInput["pt_bins"].front()["min"].asDouble();
      Par.pTMax = calibrationInput["pt_bins"].back()["max"].asDouble();
      
      PerformFitsForDifferentCentrAndZDC(std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]));
   }

   Par.outputFile->Close();
   
   return 0;
}
            
void PerformFitsForDifferentCentrAndZDC(const unsigned int detectorBin, 
                                        const unsigned int variableBin,
                                        const unsigned int particleBin)
{ 
   const Json::Value calibrationInput = Par.inputJSONCal["SIGMALIZED_RESIDUALS"];
   
   const Json::Value detector = calibrationInput["detectors_to_calibrate"][detectorBin];
   const Json::Value variable = calibrationInput["variables_to_calibrate"][variableBin];
   const Json::Value particle = calibrationInput["particles_to_calibrate"][particleBin];
   
   const std::string detectorName = detector["name"].asString();
   const std::string variableName = variable["name"].asString();
   const std::string particleName = particle["name"].asString();

   // TPad for current thread
   TVirtualPad *currentPad;
   
   for (const Json::Value& centrality : calibrationInput["centrality_bins"])
   {
      const std::string centralityRangeName = centrality["min"].asString() + "-" + 
                                              centrality["max"].asString() + "%"; 

      // unique name for this bin to prevent ROOT 
      // printing warnings about replacing same named objects
      const std::string thisBinUniqueName = 
         detectorName + variableName + particleName + centralityRangeName;
      
      // same as before but without percent; used for output file names
      const std::string centralityRangePathName = 
         "_c" + centrality["min"].asString() + "-" + centrality["max"].asString();
      
      // TDirectory pointer in output file for the current thread
      TDirectory *currentOutputTDir = 
         Par.outputFile->mkdir((detectorName + "/" + variableName + "/" + 
                               particleName + "/" + centralityRangePathName).c_str());
      
      std::vector<TGraphErrors> grVMeansVsPT, grVSigmasVsPT;
      std::vector<TF1> fVMeansVsPT, fVSigmasVsPT;
      
      // histograms with weights representing means and sigmas
      TH2D distrMeansVsZDCVsPT(("means" + thisBinUniqueName).c_str(), "#mu", 
                               calibrationInput["zdc_bins"].size(), 
                               calibrationInput["zdc_bins"].front()["min"].asDouble(), 
                               calibrationInput["zdc_bins"].back()["max"].asDouble(),
                               Par.pTRanges.size() - 1, &Par.pTRanges[0]);
      TH2D distrSigmasVsZDCVsPT(("sigmas" + thisBinUniqueName).c_str(), "#sigma", 
                                calibrationInput["zdc_bins"].size(), 
                                calibrationInput["zdc_bins"].front()["min"].asDouble(), 
                                calibrationInput["zdc_bins"].back()["max"].asDouble(),
                                Par.pTRanges.size() - 1, &Par.pTRanges[0]);
      // histograms with weights representing the difference between means and sigmas and
      // a fit parameter of means and sigmas
      TH2D distrMeansDiffVsZDCVsPT(("means diff" + thisBinUniqueName).c_str(), 
                                   "#cbar#mu - #mu_{fit}#cbar/#mu",
                                   calibrationInput["zdc_bins"].size(), 
                                   calibrationInput["zdc_bins"].front()["min"].asDouble(), 
                                   calibrationInput["zdc_bins"].back()["max"].asDouble(),
                                   Par.pTRanges.size() - 1, &Par.pTRanges[0]);
      TH2D distrSigmasDiffVsZDCVsPT(("sigmas diff" + thisBinUniqueName).c_str(), 
                                    "#cbar#sigma - #sigma_{fit}#cbar/#sigma",
                                    calibrationInput["zdc_bins"].size(), 
                                    calibrationInput["zdc_bins"].front()["min"].asDouble(), 
                                    calibrationInput["zdc_bins"].back()["max"].asDouble(),
                                    Par.pTRanges.size() - 1, &Par.pTRanges[0]);
      
      for (const Json::Value& zDC : calibrationInput["zdc_bins"])
      { 
         Par.numberOfCalls++;
          
         const std::string zDCRangeName = zDC["min"].asString() + "<zDC<" + 
                                          zDC["max"].asString();
         const std::string zDCRangePathName = "_zDC" + zDC["min"].asString() + 
                                              "-" + zDC["max"].asString();
         
         const double zDCMin = zDC["min"].asDouble();
         const double zDCMax = zDC["max"].asDouble();
         
         // name of histogram
         const std::string distrVariableName =  
            variableName + " vs pT vs centrality: " + detectorName + ", " + 
            particleName + ", " + zDCRangeName;
         
         TH3F *distrVariable = 
            static_cast<TH3F *>(Par.inputFile->Get(distrVariableName.c_str()));
         
         if (!distrVariable) 
         {
            PrintError("Histogram named \"" + distrVariableName + "\" does not exist in file " + 
                       Par.inputFile->GetName(), false);
            return;
         }
         
         std::string fitsOutputFileName = Par.outputDir + detectorName + "/" + variableName + 
                                          "_" + particleName + centralityRangePathName + 
                                          zDCRangePathName;

         grVMeansVsPT.emplace_back();
         grVSigmasVsPT.emplace_back();
         
         fVMeansVsPT.emplace_back((zDCRangeName + centralityRangeName + detectorName + 
                                   particleName + variableName).c_str(), 
                                  detector["means_fit"].asCString());
         fVSigmasVsPT.emplace_back((zDCRangeName + centralityRangeName + detectorName + 
                                    particleName + variableName).c_str(), 
                                   detector["sigmas_fit"].asCString());
                        
         PerformFitsForDifferentPT(distrVariable, grVMeansVsPT.back(), grVSigmasVsPT.back(), 
                                   currentOutputTDir, calibrationInput, detector, variable, 
                                   zDC, particle, centrality);
         
         fVMeansVsPT.back().SetRange(Par.pTMin/1.05, Par.pTMax*1.05);
         fVSigmasVsPT.back().SetRange(Par.pTMin/1.05, Par.pTMax*1.05);
         
         for (short j = 1; j <= Par.fitNTries; j++)
         {
            grVMeansVsPT.back().Fit(&fVMeansVsPT.back(), "RQN");
            grVSigmasVsPT.back().Fit(&fVSigmasVsPT.back(), "RQN");
            
            for (int k = 0; k < fVMeansVsPT.back().GetNpar(); k++)
            {
               fVMeansVsPT.back().SetParLimits(k, fVMeansVsPT.back().GetParameter(k)/
                                               (1. + 5./static_cast<double>(j*j)),
                                               fVMeansVsPT.back().GetParameter(k)*
                                               (1. + 5./static_cast<double>(j*j)));
            }
            for (int k = 0; k < fVSigmasVsPT.back().GetNpar(); k++)
            {
               fVSigmasVsPT.back().SetParLimits(k, fVSigmasVsPT.back().GetParameter(k)/
                                                (1. + 5./static_cast<double>(j*j)),
                                                fVSigmasVsPT.back().GetParameter(k)*
                                                (1. + 5./static_cast<double>(j*j)));
            }
         }

         // filling 2D histograms with weights as fit parameters means and sigmas
         for (int j = 0; j < grVMeansVsPT.back().GetN(); j++)
         {
            const double x = grVMeansVsPT.back().GetPointX(j);
            const int xBin = 
               distrMeansVsZDCVsPT.GetXaxis()->FindBin(Average(zDCMin, zDCMax));
            const int yBin = distrMeansVsZDCVsPT.GetYaxis()->FindBin(x);
            
            distrMeansVsZDCVsPT.
               SetBinContent(xBin, yBin, grVMeansVsPT.back().GetPointY(j));
            distrSigmasVsZDCVsPT.
               SetBinContent(xBin, yBin, grVSigmasVsPT.back().GetPointY(j));
            distrMeansDiffVsZDCVsPT.SetBinContent(xBin, yBin, 
                                         fabs((grVMeansVsPT.back().GetPointY(j) - 
                                               fVMeansVsPT.back().Eval(x))/
                                              grVMeansVsPT.back().GetPointY(j)));
            distrSigmasDiffVsZDCVsPT.SetBinContent(xBin, yBin, 
                                          fabs(grVSigmasVsPT.back().GetPointY(j) - 
                                               fVSigmasVsPT.back().Eval(x))/
                                          grVSigmasVsPT.back().GetPointY(j));
         }

         currentOutputTDir->WriteTObject(grVMeansVsPT.back().Clone(), 
                                         ("means: " + zDCRangeName).c_str());
         currentOutputTDir->WriteTObject(grVSigmasVsPT.back().Clone(), 
                                         ("sigmas: " + zDCRangeName).c_str());
         currentOutputTDir->WriteTObject(fVMeansVsPT.back().Clone(), 
                                         ("means fit: " + zDCRangeName).c_str());
         currentOutputTDir->WriteTObject(fVSigmasVsPT.back().Clone(), 
                                         ("sigmas fit: " + zDCRangeName).c_str());
      }
      
      double meanYMin = 1e31, meanYMax = -1e31;
      double sigmaYMin = 1e31, sigmaYMax = -1e31;
      
      for (unsigned long i = 0; i < calibrationInput["zdc_bins"].size(); i++)
      {
         meanYMin = Minimum(meanYMin, TMath::MinElement(grVMeansVsPT[i].GetN(), 
                                                        grVMeansVsPT[i].GetY()));
         meanYMax = Maximum(meanYMax, TMath::MaxElement(grVMeansVsPT[i].GetN(), 
                                                        grVMeansVsPT[i].GetY()));
         sigmaYMin = Minimum(sigmaYMin, TMath::MinElement(grVSigmasVsPT[i].GetN(), 
                                                          grVSigmasVsPT[i].GetY()));
         sigmaYMax = Maximum(sigmaYMax, TMath::MaxElement(grVSigmasVsPT[i].GetN(), 
                                                          grVSigmasVsPT[i].GetY()));

         const Json::Value zDC = calibrationInput["zdc_bins"][static_cast<int>(i)];
         
         grVMeansVsPT[i].SetMarkerStyle(zDC["marker_style"].asInt());
         grVMeansVsPT[i].SetMarkerSize(1.4);
         grVMeansVsPT[i].SetMarkerColorAlpha(zDC["color"].asInt(), 0.8);
         grVMeansVsPT[i].SetLineColorAlpha(zDC["color"].asInt(), 0.8);
         fVMeansVsPT[i].SetLineColorAlpha(zDC["color"].asInt(), 0.9);
         fVMeansVsPT[i].SetLineStyle(3);
         
         grVSigmasVsPT[i].SetMarkerStyle(zDC["marker_style"].asInt());
         grVSigmasVsPT[i].SetMarkerSize(1.4);
         grVSigmasVsPT[i].SetMarkerColorAlpha(zDC["color"].asInt(), 0.8);
         grVSigmasVsPT[i].SetLineColorAlpha(zDC["color"].asInt(), 0.8);
         fVSigmasVsPT[i].SetLineColorAlpha(zDC["color"].asInt(), 0.9);
         fVSigmasVsPT[i].SetLineStyle(2);
      }
      
      TCanvas canvValVsPTVsZDC("", "", 800, 800);
      currentPad = canvValVsPTVsZDC.cd();
      
      TLegend legend{0.15, 0.7, 0.88, 0.89};
      legend.SetNColumns(3);
      legend.SetLineColorAlpha(0, 0.);
      legend.SetFillColorAlpha(0, 0.);
      
      currentPad->SetLeftMargin(0.135);
      
      TH1F meansFrame(("means frame " + thisBinUniqueName).c_str(), "", 
                      10, Par.pTMin - 0.1, Par.pTMax*1.05);
      meansFrame.SetMinimum(meanYMin - (meanYMax - meanYMin)*0.05);
      meansFrame.SetMaximum(meanYMax + (meanYMax - meanYMin)*0.35);
      
      meansFrame.GetXaxis()->SetTitle("p_{T} [GeV/c]");
      meansFrame.GetYaxis()->
         SetTitle(("#mu_{" + variable["tex_name"].asString() + "}").c_str());
      meansFrame.GetXaxis()->SetTitleOffset(1.1);
      meansFrame.GetYaxis()->SetTitleOffset(2.0);
      currentPad->Add(meansFrame.Clone(), "AXIS");
      currentPad->Add(meansFrame.Clone(), "SAME AXIS X+ Y+");
      
      for (unsigned long i = 0; i < calibrationInput["zdc_bins"].size(); i++)
      {
         const std::string zDCRangeName = 
            calibrationInput["zdc_bins"][static_cast<int>(i)]["min"].asString() + 
            "<z_{DC}<" + 
            calibrationInput["zdc_bins"][static_cast<int>(i)]["max"].asString();
         
         legend.AddEntry(&grVMeansVsPT[i], zDCRangeName.c_str(), "P");
         currentPad->Add(grVMeansVsPT[i].Clone(), "SAME P");
         currentPad->Add(fVMeansVsPT[i].Clone(), "SAME");
      }

      currentPad->Add(legend.Clone());
      
      PrintCanvas(&canvValVsPTVsZDC, Par.outputDir + detectorName + "_means_" + 
                  variableName + "_" + particle["name_short"].asString() + 
                  centralityRangePathName);
      
      legend.Clear();
      canvValVsPTVsZDC.Clear();

      TH1F sigmasFrame(("sigmas frame " + thisBinUniqueName).c_str(), "", 
                       10, Par.pTMin - 0.1, Par.pTMax*1.05);
      sigmasFrame.SetMinimum(sigmaYMin/1.1);
      sigmasFrame.SetMaximum(sigmaYMax*1.4);
      
      sigmasFrame.GetXaxis()->SetTitle("p_{T} [GeV/c]");
      sigmasFrame.GetYaxis()->
         SetTitle(("#sigma_{" + variable["tex_name"].asString() + "}").c_str());
      currentPad->Add(sigmasFrame.Clone(), "AXIS");
      currentPad->Add(sigmasFrame.Clone(), "SAME AXIS X+ Y+");
      
      for (unsigned long i = 0; i < calibrationInput["zdc_bins"].size(); i++)
      { 
         const std::string zDCRangeName = 
            calibrationInput["zdc_bins"][static_cast<int>(i)]["min"].asString() + 
            "<z_{DC}<" + 
            calibrationInput["zdc_bins"][static_cast<int>(i)]["max"].asString();
         
         legend.AddEntry(&grVSigmasVsPT[i], zDCRangeName.c_str(), "P");
         currentPad->Add(grVSigmasVsPT[i].Clone(), "SAME P");
         currentPad->Add(fVSigmasVsPT[i].Clone(), "SAME");
      }

      currentPad->Add(legend.Clone());
      
      PrintCanvas(&canvValVsPTVsZDC, Par.outputDir + detectorName + "_sigmas_" + 
                  variableName + "_" + particle["name_short"].asString() + 
                  centralityRangePathName);
      
      TCanvas parCanv("", "", 800, 800);
      parCanv.Divide(2, 2);
      
      currentPad = parCanv.cd(1);
      currentPad->SetRightMargin(0.13);
      distrMeansVsZDCVsPT.GetXaxis()->SetTitle(variable["tex_name"].asCString());
      distrMeansVsZDCVsPT.GetYaxis()->SetTitle("p_{T}");
      currentPad->Add(distrMeansVsZDCVsPT.Clone(), "COLZ");

      currentPad = parCanv.cd(2);
      currentPad->SetRightMargin(0.13);
      distrSigmasVsZDCVsPT.GetXaxis()->SetTitle(variable["tex_name"].asCString());
      distrSigmasVsZDCVsPT.GetYaxis()->SetTitle("p_{T}");
      currentPad->Add(distrSigmasVsZDCVsPT.Clone(), "COLZ");

      currentPad = parCanv.cd(3);
      currentPad->SetLogz();
      currentPad->SetRightMargin(0.13);
      distrMeansDiffVsZDCVsPT.GetXaxis()->SetTitle(variable["tex_name"].asCString());
      distrMeansDiffVsZDCVsPT.GetYaxis()->SetTitle("p_{T}");
      currentPad->Add(distrMeansDiffVsZDCVsPT.Clone(), "COLZ");

      currentPad = parCanv.cd(4);
      currentPad->SetLogz();
      currentPad->SetRightMargin(0.13);
      distrSigmasDiffVsZDCVsPT.GetXaxis()->SetTitle(variable["tex_name"].asCString());
      distrSigmasDiffVsZDCVsPT.GetYaxis()->SetTitle("p_{T}");
      currentPad->Add(distrSigmasDiffVsZDCVsPT.Clone(), "COLZ");
      
      PrintCanvas(&parCanv, Par.outputDir + detectorName + "/fitPar_" + variableName + "_" + 
                  particle["name_short"].asString() + centralityRangePathName);

      currentOutputTDir->WriteTObject(distrMeansVsZDCVsPT.Clone(), "means: zDC vs pT");
      currentOutputTDir->WriteTObject(distrSigmasVsZDCVsPT.Clone(), "sigmas: zDC vs pT");
      currentOutputTDir->WriteTObject(distrMeansVsZDCVsPT.Clone(), "means: zDC vs pT");
      currentOutputTDir->WriteTObject(distrSigmasVsZDCVsPT.Clone(), "sigmas: zDC vs pT");
   }
}

void PerformFitsForDifferentPT(TH3F *hist, TGraphErrors &grMeans, TGraphErrors &grSigmas,
                               TDirectory *currentOutputTDir, const Json::Value& calibrationInput,
                               const Json::Value& detector, const Json::Value& variable,
                               const Json::Value& zDC, const Json::Value& particle,
                               const Json::Value& centrality)
{
   const double minBinX = hist->GetXaxis()->GetBinLowEdge(1);
   const double maxBinX = hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->GetNbins());
   const double binWidth = hist->GetXaxis()->GetBinWidth(1);

   const std::string centralityRangeName = centrality["min"].asString() + "-" + 
                                           centrality["max"].asString() + "%";
   const std::string centralityRangePathName = "_c" + centrality["min"].asString() + 
                                                 "-" + centrality["max"].asString();

   const std::string zDCRangeName = zDC["min"].asString() + "<zDC<" + zDC["max"].asString();
   const std::string zDCRangePathName = "_zDC" + zDC["min"].asString() + 
                                        "-" + zDC["max"].asString();

   // TPad for current thread
   TVirtualPad *currentPad;

   TCanvas canv("dval vs pT", "", 
                calibrationInput["pt_nbinsx"].asDouble()*400.,
                calibrationInput["pt_nbinsy"].asDouble()*400.);
   
   canv.Divide(calibrationInput["pt_nbinsx"].asInt(), calibrationInput["pt_nbinsy"].asInt());
 
   auto PerformFitsInRange = [&](const double pTMin, const double pTMax, 
                                 TF1 *meansFit, TF1 *sigmasFit,
                                 TF1 *limitMeansFitFunc = NULL, TF1 *limitSigmasFitFunc = NULL,
                                 bool limitParByFit = false)
   {
      // clearing the graphs after previous fits
      for (int i = grMeans.GetN() - 1; i >= 0; i--)
      {
         grMeans.RemovePoint(i);
         grSigmas.RemovePoint(i);
      }

      int iCanv = 1;
      
      for (const Json::Value& pTBin : calibrationInput["pt_bins"])
      {
         const double pT = Average(pTBin["min"].asDouble(), pTBin["max"].asDouble());
         if (pT < pTMin || pT > pTMax) continue;
         
         currentPad = canv.cd(iCanv);
         iCanv++;
      
         TH1D *distrVariableProj = hist->
            ProjectionX(((std::string) hist->GetName() + "_projX").c_str(), 
                        hist->GetYaxis()->FindBin(pTBin["min"].asDouble() + 1e-6), 
                        hist->GetYaxis()->FindBin(pTBin["max"].asDouble() - 1e-6),
                        hist->GetZaxis()->FindBin(centrality["min"].asDouble() + 1e-6),
                        hist->GetZaxis()->FindBin(centrality["max"].asDouble() - 1e-6));

         const std::string pTRangeName = DtoStr(pTBin["min"].asDouble(), 1) + "<pT<" + 
                                         DtoStr(pTBin["max"].asDouble(), 1);
         
         if (distrVariableProj->Integral(1, distrVariableProj->GetXaxis()->GetNbins()) < 
             Par.minIntegralValue) 
         {
            PrintInfo("Integral is insufficient for projection of " + variable["name"].asString() + 
                      ", " + detector["name"].asString() + ", " + particle["name"].asString() + 
                      " at " + zDCRangeName + ", " + centralityRangeName + ", " + pTRangeName);
            continue;
         }
         
         double minX = 0., maxX = -1.;
         for (int k = 1; k <= distrVariableProj->GetXaxis()->GetNbins(); k++)
         {
            if (distrVariableProj->GetBinContent(k) > 1e-7)
            {
               minX = distrVariableProj->GetXaxis()->GetBinLowEdge(k);
               break;
            }
         }
         
         for (int k = distrVariableProj->GetXaxis()->GetNbins(); k > 1; k--)
         {
            if (distrVariableProj->GetBinContent(k) > 1e-7)
            {
               maxX = distrVariableProj->GetXaxis()->GetBinUpEdge(k);
               break;
            }
         }
         
         if (minX > maxX) 
         {
            PrintWarning("Something wrong for projection of " + variable["name"].asString() + ", " + 
                         detector["name"].asString() + ", " + particle["name"].asString() + 
                         " at " + zDCRangeName + ", " + centralityRangeName + ", " + pTRangeName);
            continue;
         }

         TF1 fitFuncGaus("gaus", "gaus");
         TF1 fitFuncDVal("fitFunc", "gaus(0) + gaus(3)");
         TF1 fitFuncBG("bg", "gaus");
         fitFuncGaus.SetParameters(1., 0., binWidth*2.);
         fitFuncDVal.SetParameters(1., 0., binWidth*2., 1., 0., maxX/2.);
         
         if (!limitParByFit)
         {
            fitFuncGaus.SetParLimits(1, minBinX/5., maxBinX/5.);
            fitFuncGaus.SetParLimits(2, binWidth, maxBinX/5.);
            fitFuncDVal.SetParLimits(1, minX/10., maxX/10.);
            fitFuncDVal.SetParLimits(2, binWidth, Average(maxX, maxX, minX));
         }
         else
         {
            fitFuncGaus.SetParLimits(1, limitMeansFitFunc->Eval(pT)/1.5, 
                                     limitMeansFitFunc->Eval(pT)*1.5);
            fitFuncGaus.SetParLimits(2, limitSigmasFitFunc->Eval(pT)/1.5, 
                                     limitSigmasFitFunc->Eval(pT)*1.5);
            fitFuncDVal.SetParLimits(1, limitMeansFitFunc->Eval(pT)/1.5, 
                                     limitMeansFitFunc->Eval(pT)*1.5);
            fitFuncDVal.SetParLimits(2, limitSigmasFitFunc->Eval(pT)/1.5, 
                                     limitSigmasFitFunc->Eval(pT)*1.5);
         }
         
         fitFuncDVal.SetParLimits(4, minX*2., maxX*2.);
         fitFuncDVal.SetParLimits(5, maxX/3., maxX*3.);
         
         fitFuncDVal.SetLineColorAlpha(kRed+1, 0.6);
         fitFuncBG.SetLineColorAlpha(kGreen+1, 0.9);
         fitFuncBG.SetLineStyle(2);
         fitFuncGaus.SetLineColorAlpha(kAzure-3, 0.9);
         fitFuncGaus.SetLineStyle(2);

         distrVariableProj->GetXaxis()->SetTitle(variable["tex_name"].asCString());
         distrVariableProj->SetTitle("");
         distrVariableProj->SetTitleSize(0.06, "X");
         distrVariableProj->SetTitleSize(0.06, "Y");
         distrVariableProj->SetLabelSize(0.06, "X");
         distrVariableProj->SetLabelSize(0.06, "Y");
    
         distrVariableProj->GetXaxis()->SetRange(distrVariableProj->GetXaxis()->FindBin(minX+0.01),
                                                 distrVariableProj->GetXaxis()->FindBin(maxX-0.01));
         
         const double maxBinVal = 
            distrVariableProj->GetBinContent(distrVariableProj->GetMaximumBin());
         
         // scale limits
         fitFuncGaus.SetParLimits(0, maxBinVal/5., maxBinVal);
         fitFuncDVal.SetParLimits(0, maxBinVal/5., maxBinVal);
         fitFuncDVal.SetParLimits(3, maxBinVal/20., maxBinVal);
         
         fitFuncGaus.SetRange(minBinX/5., maxBinX/5.);
         
         distrVariableProj->Fit(&fitFuncGaus, "RQBN");

         fitFuncGaus.SetRange(minBinX, maxBinX);
         fitFuncBG.SetRange(minBinX, maxBinX);
         
         for (int k = 0; k < 3; k++)
         {
            fitFuncDVal.SetParameter(k, fitFuncGaus.GetParameter(k));
         }
         
         fitFuncDVal.SetRange(minBinX, maxBinX);
         distrVariableProj->Fit(&fitFuncDVal, "RQBN");

         for (unsigned short j = 1; j <= Par.fitNTries; j++)
         {
            for (int k = 0; k < fitFuncDVal.GetNpar(); k++)
            {
               if (k != 2) // sigmas cannot be negative 
               {
                  fitFuncDVal.SetParLimits(k, fitFuncDVal.GetParameter(k)*
                                           (1. - 5./static_cast<double>(j*j)),
                                           fitFuncDVal.GetParameter(k)*
                                           (1. + 5./static_cast<double>(j*j)));
               }
               else
               {
                  fitFuncDVal.SetParLimits(k, fitFuncDVal.GetParameter(k)/
                                           (1. + 5./static_cast<double>(j*j)),
                                           fitFuncDVal.GetParameter(k)*
                                           (1. + 5./static_cast<double>(j*j)));
               }
            }
            
            distrVariableProj->Fit(&fitFuncDVal, "RQBN");
         }

         for (int j = 0; j < 3; j++)
         {
            fitFuncGaus.SetParameter(j, fitFuncDVal.GetParameter(j));
            fitFuncBG.SetParameter(j, fitFuncDVal.GetParameter(j + 3));
         }
         
         distrVariableProj->SetMarkerStyle(20);
         distrVariableProj->SetMarkerSize(0.7);
         distrVariableProj->SetMarkerColorAlpha(kBlack, 0.8);
         distrVariableProj->SetLineColorAlpha(kBlack, 0.8);
         distrVariableProj->SetMaximum(maxBinVal*1.2);
         
         currentPad->SetLeftMargin(0.155);
         currentPad->SetBottomMargin(0.12);
         
         currentPad->Add(distrVariableProj->Clone(), "P");
         currentPad->Add(fitFuncDVal.Clone(), "SAME");
         currentPad->Add(fitFuncBG.Clone(), "SAME");
         currentPad->Add(fitFuncGaus.Clone(), "SAME");
         
         Par.texText.SetText(0.17, 0.85, pTRangeName.c_str());
         currentPad->Add(Par.texText.Clone());
         Par.texText.SetText(0.17, 0.79, zDCRangeName.c_str());
         currentPad->Add(Par.texText.Clone());
         Par.texText.SetText(0.17, 0.73, particle["name"].asCString());
         currentPad->Add(Par.texText.Clone());
         Par.texText.SetText(0.17, 0.66, centralityRangeName.c_str());
         currentPad->Add(Par.texText.Clone());
         
         if (fabs(fitFuncDVal.GetParameter(1)) < variable["abs_max_fit_val"].asDouble() && 
             fabs(fitFuncDVal.GetParameter(2)) < variable["abs_max_fit_val"].asDouble())
         {
            grMeans.AddPoint(pT, fitFuncDVal.GetParameter(1));
            grSigmas.AddPoint(pT, fabs(fitFuncDVal.GetParameter(2)));
            
            //grMeans.SetPointError(grMeans.GetN() - 1, 0, fitFuncDVal.GetParError(1));
            //grSigmas.SetPointError(grSigmas.GetN() - 1, 0, fitFuncDVal.GetParError(2));
         }
      }
      if (grMeans.GetN() == 0) 
      {
         PrintError("Graph is empty for " + variable["name"].asString() + ", " + 
                    detector["name"].asString() + ", " + particle["name"].asString() + 
                    " at " + zDCRangeName + ", " + centralityRangeName, false);
         return;
      }
      
      meansFit->SetRange(pTMin, pTMax);
      sigmasFit->SetRange(pTMin, pTMax);
      
      grMeans.Fit(meansFit, "RQNW");
      grSigmas.Fit(sigmasFit, "RQNW");
   };
   
   TF1 meansFit("meansFit", Par.meansFitPrelimFunc.c_str());
   TF1 sigmasFit("sigmasFit", Par.sigmasFitPrelimFunc.c_str());
   
   for (const Json::Value& pT : calibrationInput["consequetive_pt_fit_ranges"])
   {
      PerformFitsInRange(pT["min"].asDouble(), pT["max"].asDouble(), &meansFit, &sigmasFit);
   }

   TF1 *oldMeansFit = (TF1 *) meansFit.Clone();
   TF1 *oldSigmasFit = (TF1 *) sigmasFit.Clone();
   
   meansFit.SetCurrent((TF1 *) TF1("meansFit", detector["means_fit"].asCString()).Clone());
   sigmasFit.SetCurrent((TF1 *) TF1("sigmasFit", detector["sigmas_fit"].asCString()).Clone());

   PerformFitsInRange(calibrationInput["pt_bins"].front()["min"].asDouble()/1.05, 
                      calibrationInput["pt_bins"].back()["max"].asDouble()*1.05, 
                      &meansFit, &sigmasFit, oldMeansFit, oldSigmasFit, true);
   
   const std::string outputFileNameNoExt = 
      "output/SigmalizedResiduals/" + Par.runName + "/" + detector["name"].asString() + "/" + 
      variable["name"].asString() + "_" + particle["name_short"].asString() + 
      centralityRangePathName + zDCRangePathName;
   PrintCanvas(&canv, outputFileNameNoExt, false, true);
}

#endif /* SIGMALIZED_RESIDUALS_CPP */
