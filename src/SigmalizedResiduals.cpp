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

// This program works in 2 modes
// Mode2 analyzes track residual variable for the specific detector, variable, and charge
// Mode1 recursively calls the program outside the current process to analyze all configurations
int main(int argc, char **argv)
{
   if (argc < 2 || argc > 6 || (argc > 3 && argc < 4)) 
   {
      std::string errMsg = "Expected 1-2 or 5-7 parameters while " + std::to_string(argc - 1) + 
                           " parameter(s) were provided \n";
      errMsg += "Usage: bin/SigmalizedResiduals inputFile numberOfThreads=" + 
                std::to_string(std::thread::hardware_concurrency()) + "*\n";
      errMsg += "Or**: bin/SigmalizedResiduals inputFile detectorBin variableBin ";
      errMsg += "numberOfThreads=" + std::to_string(std::thread::hardware_concurrency()) + "* " +
                "showProgress=true";
      errMsg += "*: default argument is the number of threads on the current machine \n";
      errMsg += "**: this mode analyzes only one configuration";
      PrintError(errMsg);
   }

   unsigned int numberOfThreads;

   // initializing ROOT parameters
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

   if (Par.inputJSONCal["detectors_to_calibrate"].size() == 0)
   {
      PrintInfo("No detectors are specified for calibrations");
      PrintInfo("Exiting the program");
      exit(1);
   }
   
   if (argc < 4) // Mode1
   {
      Par.programMode = 1;
      if (argc > 2) numberOfThreads = std::stoi(argv[2]);
      else numberOfThreads = std::thread::hardware_concurrency();
      if (numberOfThreads == 0) PrintError("Number of threads must be bigger than 0");

      system("rm -rf tmp/SigmalizedResiduals/*");
      system(("mkdir -p tmp/SigmalizedResiduals/" + Par.runName).c_str());

      Par.numberOfIterations = Par.inputJSONCal["detectors_to_calibrate"].size()*
                               Par.inputJSONCal["centrality_bins"].size()*
                               Par.inputJSONCal["zdc_bins"].size()*4;

      for (const auto& detector : Par.inputJSONCal["detectors_to_calibrate"])
      {
         system(("mkdir -p output/SigmalizedResiduals/" + Par.runName + 
                 "/" + detector["name"].asString()).c_str());
      }
      
      auto SingleThreadCall = [&](const unsigned long detectorBin, 
                                  const unsigned long variableBin)
      {
         // man ROOT sucks (TF1::Fit is still not thread safe) so I have to call the same program 
         // recursively outside of the current instance to implement multithreading
         system((static_cast<std::string>("./bin/SigmalizedResiduals ") + 
                 argv[1] + " " + std::to_string(detectorBin) + " " + 
                 std::to_string(variableBin) + " 1 0").c_str());;
      };

      std::vector<std::thread> thrCalls;
      std::thread pBarThr(PBarCall); 

      for (unsigned int detectorBin = 0; detectorBin < 
           Par.inputJSONCal["detectors_to_calibrate"].size(); detectorBin++)
      {
         for (unsigned long variableBin = 0; variableBin < Par.variableName.size(); variableBin++)
         { 
               if (thrCalls.size() >= numberOfThreads)
               {
                  while (!thrCalls.empty())
                  {
                     thrCalls.back().join();
                     thrCalls.pop_back();
                  }
               }
               thrCalls.emplace_back(SingleThreadCall, detectorBin, variableBin);
            }
      }

      while (!thrCalls.empty())
      {
         thrCalls.back().join();
         thrCalls.pop_back();
      }
      
      Par.isProcessFinished = true;
      pBarThr.join();
   }
   else // Mode2
   {
      Par.programMode = 2;
      if (argc > 4) numberOfThreads = std::stoi(argv[4]);
      else numberOfThreads = std::thread::hardware_concurrency();
      if (numberOfThreads == 0) PrintError("Number of threads must be bigger than 0");
      
      ROOT::EnableImplicitMT(numberOfThreads);
      
      if (argc > 5) Par.showProgress = static_cast<bool>(std::stoi(argv[5]));
      
      Par.inputFile = std::unique_ptr<TFile>(TFile::Open(("data/Real/" + Par.runName + 
                                                          "/SingleTrack/sum.root").c_str(), "READ"));
      
      Par.texText.SetTextFont(52);
      Par.texText.SetTextSize(0.06);

      Par.numberOfIterations = 2.*Par.inputJSONCal["centrality_bins"].size()*
                               Par.inputJSONCal["zdc_bins"].size();
      Par.texText.SetNDC();

      for (const Json::Value& pTBin: Par.inputJSONCal["pt_bins"])
      {
         Par.pTRanges.push_back(pTBin["min"].asDouble());
      }
      Par.pTRanges.push_back(Par.inputJSONCal["pt_bins"].back()["max"].asDouble());

      for (const Json::Value& centrality: Par.inputJSONCal["centrality_bins"])
      {
         Par.centralityRanges.push_back(centrality["min"].asDouble());
      }
      Par.centralityRanges.push_back(Par.inputJSONCal["centrality_bins"].back()["max"].asDouble());

      Par.outputDir = "output/SigmalizedResiduals/" + Par.runName + "/";
      
      Par.pTMin = Par.inputJSONCal["pt_bins"].front()["min"].asDouble();
      Par.pTMax = Par.inputJSONCal["pt_bins"].back()["max"].asDouble();
      
      std::thread pBarThr(PBarCall); 
      
      PerformFitsForDifferentCentrAndZDC(std::stoi(argv[2]), std::stoi(argv[3]));
      
      Par.isProcessFinished = true;
      pBarThr.join();
   }
 
   return 0;
}
            
void PerformFitsForDifferentCentrAndZDC(const unsigned int detectorBin, 
                                        const unsigned int variableBin)
{   
   const Json::Value detector = Par.inputJSONCal["detectors_to_calibrate"][detectorBin];
   
   const std::string detectorName = detector["name"].asString();
   const std::string variableName = Par.variableName[variableBin];
   
   Par.outputFile = 
      std::unique_ptr<TFile>(TFile::Open((Par.outputDir + detectorName + "/all_fits_" + 
                                          variableName + ".root").c_str(), "RECREATE"));

   for (int charge : Par.particleCharges)
   {
      const std::string chargeName = ((charge > 0) ? "charge>0" : "charge<0");
      const std::string chargeNameShort = ((charge > 0) ? "pos" : "neg");

      // output file in which parameters of approximations will be written
      std::ofstream parametersOutputMeans(Par.outputDir + "par_means_" + detectorName + "_" + 
                                          variableName + "_" + chargeNameShort + ".txt");
      std::ofstream parametersOutputSigmas(Par.outputDir + "par_sigmas_" + detectorName + "_" + 
                                           variableName + "_" + chargeNameShort + ".txt");

      for (unsigned int centralityBin = 0; centralityBin < 
           Par.inputJSONCal["centrality_bins"].size(); centralityBin++)
      {
         const Json::Value centrality = Par.inputJSONCal["centrality_bins"][centralityBin];
         
         const std::string centralityRangeName = centrality["min"].asString() + "-" + 
                                                 centrality["max"].asString() + "%"; 
         
         // unique name for this bin to prevent ROOT 
         // printing warnings about replacing same named objects
         const std::string thisBinUniqueName = 
            detectorName + variableName + chargeName + centralityRangeName;
         
         // same as before but without percent; used for output file names
         const std::string centralityRangePathName = 
            "_c" + centrality["min"].asString() + "-" + centrality["max"].asString();
         
         // TDirectory pointer in output file for the current thread
         const std::string outputFileDirName = detectorName + "/" + variableName + "/" + 
                                               chargeName + "/" + centralityRangePathName;
         Par.outputFile->mkdir(outputFileDirName.c_str());
         Par.outputFile->cd(outputFileDirName.c_str());
         
         std::vector<TGraphErrors> grVMeansVsPT, grVSigmasVsPT;
         std::vector<TF1> fVMeansVsPT, fVSigmasVsPT;
         
         // histograms with weights representing means and sigmas
         TH2D distrMeansVsZDCVsPT(("means" + thisBinUniqueName).c_str(), "#mu", 
                                  Par.inputJSONCal["zdc_bins"].size(), 
                                  Par.inputJSONCal["zdc_bins"].front()["min"].asDouble(), 
                                  Par.inputJSONCal["zdc_bins"].back()["max"].asDouble(),
                                  Par.pTRanges.size() - 1, &Par.pTRanges[0]);
         TH2D distrSigmasVsZDCVsPT(("sigmas" + thisBinUniqueName).c_str(), "#sigma", 
                                   Par.inputJSONCal["zdc_bins"].size(), 
                                   Par.inputJSONCal["zdc_bins"].front()["min"].asDouble(), 
                                   Par.inputJSONCal["zdc_bins"].back()["max"].asDouble(),
                                   Par.pTRanges.size() - 1, &Par.pTRanges[0]);
         // histograms with weights representing the difference between means and sigmas and
         // a fit parameter of means and sigmas
         TH2D distrMeansDiffVsZDCVsPT(("means diff" + thisBinUniqueName).c_str(), 
                                      "#cbar#mu - #mu_{fit}#cbar/#mu",
                                      Par.inputJSONCal["zdc_bins"].size(), 
                                      Par.inputJSONCal["zdc_bins"].front()["min"].asDouble(), 
                                      Par.inputJSONCal["zdc_bins"].back()["max"].asDouble(),
                                      Par.pTRanges.size() - 1, &Par.pTRanges[0]);
         TH2D distrSigmasDiffVsZDCVsPT(("sigmas diff" + thisBinUniqueName).c_str(), 
                                       "#cbar#sigma - #sigma_{fit}#cbar/#sigma",
                                       Par.inputJSONCal["zdc_bins"].size(), 
                                       Par.inputJSONCal["zdc_bins"].front()["min"].asDouble(), 
                                       Par.inputJSONCal["zdc_bins"].back()["max"].asDouble(),
                                       Par.pTRanges.size() - 1, &Par.pTRanges[0]);
         
         for (const Json::Value& zDC : Par.inputJSONCal["zdc_bins"])
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
               chargeName + ", " + zDCRangeName;
            
            TH3F *distrVariable = 
               static_cast<TH3F *>(Par.inputFile->Get(distrVariableName.c_str()));
            
            if (!distrVariable) 
            {
               PrintError("Histogram named \"" + distrVariableName + "\" does not exist in file " + 
                          Par.inputFile->GetName());
            }
            
            std::string fitsOutputFileName = Par.outputDir + detectorName + "/" + 
                                             variableName + "_" + chargeName + 
                                             centralityRangePathName + zDCRangePathName;

            grVMeansVsPT.emplace_back();
            grVSigmasVsPT.emplace_back();
            
            // lambda expressions for TF1
            const std::string meansFitFunc = 
               "[](double *x, double *par) {return " + 
               detector["means_fit_func_" + variableName + "_" + chargeNameShort].asString() + ";}";
            const std::string sigmasFitFunc = 
               "[](double *x, double *par) {return " + 
               detector["sigmas_fit_func_" + variableName + "_" + chargeNameShort].asString() + ";}";
            
            fVMeansVsPT.
               emplace_back((zDCRangeName + centralityRangeName + 
                             detectorName + chargeName + variableName).c_str(), 
                            meansFitFunc.c_str(), 0., 1., GetNumberOfParameters(meansFitFunc));
            fVSigmasVsPT.
               emplace_back((zDCRangeName + centralityRangeName + 
                             detectorName + chargeName + variableName).c_str(), 
                            sigmasFitFunc.c_str(), 0., 1., GetNumberOfParameters(sigmasFitFunc));
                           
            PerformFitsForDifferentPT(distrVariable, grVMeansVsPT.back(), grVSigmasVsPT.back(), 
                                      detector, variableBin, zDC, charge, centrality);

            fVMeansVsPT.back().SetRange(Par.pTMin/1.05, Par.pTMax*1.05);
            fVSigmasVsPT.back().SetRange(Par.pTMin/1.05, Par.pTMax*1.05);

            // these TGraphsError will be fitted instead of actual ones; this is due to the fact
            // that some weights are really big for some points, usually somewhere in the middle.
            // This makes ROOT almost ignore other points and perform uncanny fits.
            // The idea behind these copies is that all points will have equal weight except for
            // those that have large uncertainty and they will have smaller weight. This will
            // allow for ROOT to not bother about some inconsistent points that "jumped out" of the 
            // distribution and not to overestimate very small weight points.
            TGraphErrors meansVsPTForFit(grVMeansVsPT.back());
            TGraphErrors sigmasVsPTForFit(grVSigmasVsPT.back());

            for (int i = 0; i < meansVsPTForFit.GetN(); i++)
            {
               // for means errors sigma measurements are used since uncertainty of mean would not
               // depend on it's position while value of mean can be 0 or it can be some much larger
               // number which makes the relative uncertainties inconsistent
               meansVsPTForFit.SetPointError(i, 0., 1. + meansVsPTForFit.GetErrorY(i)/
                                                    sigmasVsPTForFit.GetPointY(i));
               sigmasVsPTForFit.SetPointError(i, 0., 1. + sigmasVsPTForFit.GetErrorY(i)/
                                                     sigmasVsPTForFit.GetPointY(i));
            }
            
            for (unsigned int i = 1; i <= Par.fitNTries; i++)
            {
               meansVsPTForFit.Fit(&fVMeansVsPT.back(), "RQMBN");
               sigmasVsPTForFit.Fit(&fVSigmasVsPT.back(), "RQMBN");
               
               for (int j = 0; j < fVMeansVsPT.back().GetNpar(); j++)
               {
                  fVMeansVsPT.back().SetParLimits(j, fVMeansVsPT.back().GetParameter(j)*
                                                  (1. - 6./static_cast<double>(i*i*i)),
                                                  fVMeansVsPT.back().GetParameter(j)*
                                                  (1. + 4./static_cast<double>(i*i*i)));
               }
               for (int j = 0; j < fVSigmasVsPT.back().GetNpar(); j++)
               {
                  fVSigmasVsPT.back().SetParLimits(j, fVSigmasVsPT.back().GetParameter(j)*
                                                   (1. - 6./static_cast<double>(i*i*i)),
                                                   fVSigmasVsPT.back().GetParameter(j)*
                                                   (1. + 4./static_cast<double>(i*i*i)));
               }
            }
            
            /*
            grVMeansVsPT.back().Fit(&fVMeansVsPT.back(), "RQMBN");
            grVSigmasVsPT.back().Fit(&fVSigmasVsPT.back(), "RQMBN");
            */
            
            for (int i = 0; i < fVMeansVsPT.back().GetNpar(); i++)
            {
               parametersOutputMeans << fVMeansVsPT.back().GetParameter(i);
               if (i < fVMeansVsPT.back().GetNpar() - 1) parametersOutputMeans << " ";
            }
            parametersOutputMeans << std::endl;

            for (int i = 0; i < fVSigmasVsPT.back().GetNpar(); i++)
            {
               parametersOutputSigmas << fVSigmasVsPT.back().GetParameter(i);
               if (i < fVSigmasVsPT.back().GetNpar() - 1) parametersOutputSigmas << " ";
            }
            parametersOutputSigmas << std::endl;

            // filling 2D histograms with weights as fit parameters means and sigmas
            for (int i = 0; i < grVMeansVsPT.back().GetN(); i++)
            {
               const double x = grVMeansVsPT.back().GetPointX(i);
               const int xBin = 
                  distrMeansVsZDCVsPT.GetXaxis()->FindBin(Average(zDCMin, zDCMax));
               const int yBin = distrMeansVsZDCVsPT.GetYaxis()->FindBin(x);
               
               distrMeansVsZDCVsPT.
                  SetBinContent(xBin, yBin, grVMeansVsPT.back().GetPointY(i));
               distrSigmasVsZDCVsPT.
                  SetBinContent(xBin, yBin, grVSigmasVsPT.back().GetPointY(i));
               distrMeansDiffVsZDCVsPT.SetBinContent(xBin, yBin, 
                                            fabs((grVMeansVsPT.back().GetPointY(i) - 
                                                  fVMeansVsPT.back().Eval(x))/
                                                 grVMeansVsPT.back().GetPointY(i)));
               distrSigmasDiffVsZDCVsPT.SetBinContent(xBin, yBin, 
                                             fabs(grVSigmasVsPT.back().GetPointY(i) - 
                                                  fVSigmasVsPT.back().Eval(x))/
                                             grVSigmasVsPT.back().GetPointY(i));
            }

            grVMeansVsPT.back().Write(("means: " + zDCRangeName).c_str());
            grVSigmasVsPT.back().Write(("sigmas: " + zDCRangeName).c_str());
            fVMeansVsPT.back().Write(("means fit: " + zDCRangeName).c_str());
            fVSigmasVsPT.back().Write(("sigmas fit: " + zDCRangeName).c_str());
            
            if (!Par.showProgress)
            {
               std::ofstream progressFile("tmp/SigmalizedResiduals/" + Par.runName + "/" + 
                                     std::to_string(detectorBin) + std::to_string(variableBin));
               progressFile << Par.numberOfCalls;
            }
         }
         
         double meanYMin = 1e31, meanYMax = -1e31;
         double sigmaYMin = 1e31, sigmaYMax = -1e31;
         
         for (unsigned long i = 0; i < Par.inputJSONCal["zdc_bins"].size(); i++)
         {
            meanYMin = Minimum(meanYMin, TMath::MinElement(grVMeansVsPT[i].GetN(), 
                                                           grVMeansVsPT[i].GetY()));
            meanYMax = Maximum(meanYMax, TMath::MaxElement(grVMeansVsPT[i].GetN(), 
                                                           grVMeansVsPT[i].GetY()));
            sigmaYMin = Minimum(sigmaYMin, TMath::MinElement(grVSigmasVsPT[i].GetN(), 
                                                             grVSigmasVsPT[i].GetY()));
            sigmaYMax = Maximum(sigmaYMax, TMath::MaxElement(grVSigmasVsPT[i].GetN(), 
                                                             grVSigmasVsPT[i].GetY()));

            const Json::Value zDC = Par.inputJSONCal["zdc_bins"][static_cast<int>(i)];
            
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
         canvValVsPTVsZDC.cd();
         
         TLegend legend{0.15, 0.7, 0.88, 0.89};
         legend.SetNColumns(3);
         legend.SetLineColorAlpha(0, 0.);
         legend.SetFillColorAlpha(0, 0.);
         
         gPad->SetLeftMargin(0.135);
         
         TH1F meansFrame(("means frame " + thisBinUniqueName).c_str(), "", 
                         10, Par.pTMin - 0.1, Par.pTMax*1.05);
         meansFrame.SetMinimum(meanYMin - (meanYMax - meanYMin)*0.05);
         meansFrame.SetMaximum(meanYMax + (meanYMax - meanYMin)*0.35);
         
         meansFrame.GetXaxis()->SetTitle("p_{T} [GeV/c]");
         meansFrame.GetYaxis()->
            SetTitle(("#mu_{" + Par.variableNameTex[variableBin] + "}").c_str());
         meansFrame.GetXaxis()->SetTitleOffset(1.1);
         meansFrame.GetYaxis()->SetTitleOffset(2.0);
         gPad->Add(&meansFrame, "AXIS");
         gPad->Add(&meansFrame, "SAME AXIS X+ Y+");
         
         for (unsigned long i = 0; i < Par.inputJSONCal["zdc_bins"].size(); i++)
         {
            const std::string zDCRangeName = 
               Par.inputJSONCal["zdc_bins"][static_cast<int>(i)]["min"].asString() + 
               "<z_{DC}<" + 
               Par.inputJSONCal["zdc_bins"][static_cast<int>(i)]["max"].asString();
            
            legend.AddEntry(&grVMeansVsPT[i], zDCRangeName.c_str(), "P");
            gPad->Add(&grVMeansVsPT[i], "SAME P");
            gPad->Add(&fVMeansVsPT[i], "SAME");
         }

         gPad->Add(&legend);
         
         PrintCanvas(&canvValVsPTVsZDC, Par.outputDir + detectorName + "_means_" + 
                     variableName + "_" + chargeNameShort + 
                     centralityRangePathName);
         
         legend.Clear();
         canvValVsPTVsZDC.Clear();

         TH1F sigmasFrame(("sigmas frame " + thisBinUniqueName).c_str(), "", 
                          10, Par.pTMin - 0.1, Par.pTMax*1.05);
         sigmasFrame.SetMinimum(sigmaYMin/1.1);
         sigmasFrame.SetMaximum(sigmaYMax*1.4);
         
         sigmasFrame.GetXaxis()->SetTitle("p_{T} [GeV/c]");
         sigmasFrame.GetYaxis()->
            SetTitle(("#sigma_{" + Par.variableNameTex[variableBin] + "}").c_str());
         gPad->Add(&sigmasFrame, "AXIS");
         gPad->Add(&sigmasFrame, "SAME AXIS X+ Y+");
         
         for (unsigned long i = 0; i < Par.inputJSONCal["zdc_bins"].size(); i++)
         { 
            const std::string zDCRangeName = 
               Par.inputJSONCal["zdc_bins"][static_cast<int>(i)]["min"].asString() + 
               "<z_{DC}<" + 
               Par.inputJSONCal["zdc_bins"][static_cast<int>(i)]["max"].asString();
            
            legend.AddEntry(&grVSigmasVsPT[i], zDCRangeName.c_str(), "P");
            gPad->Add(&grVSigmasVsPT[i], "SAME P");
            gPad->Add(&fVSigmasVsPT[i], "SAME");
         }

         gPad->Add(&legend);
         
         PrintCanvas(&canvValVsPTVsZDC, Par.outputDir + detectorName + "_sigmas_" + 
                     variableName + "_" + chargeNameShort + 
                     centralityRangePathName);
         
         TCanvas canvPar("", "", 800, 800);

         canvPar.Divide(2, 2);
         
         canvPar.cd(1);
         gPad->SetRightMargin(0.13);
         distrMeansVsZDCVsPT.GetXaxis()->SetTitle(Par.variableNameTex[variableBin].c_str());
         distrMeansVsZDCVsPT.GetYaxis()->SetTitle("p_{T}");
         gPad->Add(&distrMeansVsZDCVsPT, "COLZ");
         
         canvPar.cd(2);
         gPad->SetRightMargin(0.13);
         distrSigmasVsZDCVsPT.GetXaxis()->SetTitle(Par.variableNameTex[variableBin].c_str());
         distrSigmasVsZDCVsPT.GetYaxis()->SetTitle("p_{T}");
         gPad->Add(&distrSigmasVsZDCVsPT, "COLZ");
         
         canvPar.cd(3);
         gPad->SetLogz();
         gPad->SetRightMargin(0.13);
         distrMeansDiffVsZDCVsPT.GetXaxis()->SetTitle(Par.variableNameTex[variableBin].c_str());
         distrMeansDiffVsZDCVsPT.GetYaxis()->SetTitle("p_{T}");
         gPad->Add(&distrMeansDiffVsZDCVsPT, "COLZ");

         canvPar.cd(4);
         gPad->SetLogz();
         gPad->SetRightMargin(0.13);
         distrSigmasDiffVsZDCVsPT.GetXaxis()->SetTitle(Par.variableNameTex[variableBin].c_str());
         distrSigmasDiffVsZDCVsPT.GetYaxis()->SetTitle("p_{T}");
         gPad->Add(&distrSigmasDiffVsZDCVsPT, "COLZ");
         
         PrintCanvas(&canvPar, Par.outputDir + detectorName + 
                     "/fitPar_" + variableName + "_" + 
                     chargeNameShort + centralityRangePathName);

         distrMeansVsZDCVsPT.Write("means: zDC vs pT");
         distrSigmasVsZDCVsPT.Write("sigmas: zDC vs pT");
         distrMeansVsZDCVsPT.Write("means: zDC vs pT");
         distrSigmasVsZDCVsPT.Write("sigmas: zDC vs pT");
      }
      
      parametersOutputMeans.close();
      parametersOutputSigmas.close();
   }
   
   Par.outputFile->Close();
}

void PerformFitsForDifferentPT(TH3F *hist, TGraphErrors &grMeans, TGraphErrors &grSigmas,
                               const Json::Value& detector, const unsigned int variableBin, 
                               const Json::Value& zDC, const int charge, 
                               const Json::Value& centrality)
{
   const double minBinX = hist->GetXaxis()->GetBinLowEdge(1);
   const double maxBinX = hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->GetNbins());
   const double binWidth = hist->GetXaxis()->GetBinWidth(1);

   const std::string variableName = Par.variableName[variableBin];
   const std::string chargeName = ((charge > 0) ? "charge>0" : "charge<0");
   const std::string chargeNameShort = ((charge > 0) ? "pos" : "neg");

   const std::string centralityRangeName = centrality["min"].asString() + "-" + 
                                           centrality["max"].asString() + "%";
   const std::string centralityRangePathName = "_c" + centrality["min"].asString() + 
                                                 "-" + centrality["max"].asString();

   const std::string zDCRangeName = zDC["min"].asString() + "<zDC<" + zDC["max"].asString();
   const std::string zDCRangePathName = "_zDC" + zDC["min"].asString() + 
                                        "-" + zDC["max"].asString();

   TCanvas canvDValVsPT("dval vs pT", "", 
                Par.inputJSONCal["pt_nbinsx"].asDouble()*400.,
                Par.inputJSONCal["pt_nbinsy"].asDouble()*400.);
   
   canvDValVsPT.Divide(Par.inputJSONCal["pt_nbinsx"].asInt(), 
                       Par.inputJSONCal["pt_nbinsy"].asInt());


   // functions for fits of dz and dphi distributions
   // vectors are needed for the object to not be deleted out of scope which 
   // would result in deleting it from canvas
   std::vector<TF1> fitFuncDVal, fitFuncGaus, fitFuncBG;
   // alternate fits in different region for estimation of uncertainty
 
   int iCanv = 1;
    
   for (const Json::Value& pTBin : Par.inputJSONCal["pt_bins"])
   {
      const double pT = Average(pTBin["min"].asDouble(), pTBin["max"].asDouble());
      if (pT < Par.pTMin || pT > Par.pTMax) continue;
       
      TH1D *distrVariableProj = hist->
         ProjectionX(((std::string) hist->GetName() + "_projX_" + std::to_string(pT)).c_str(), 
                     hist->GetYaxis()->FindBin(pTBin["min"].asDouble() + 1e-6), 
                     hist->GetYaxis()->FindBin(pTBin["max"].asDouble() - 1e-6),
                     hist->GetZaxis()->FindBin(centrality["min"].asDouble() + 1e-6),
                     hist->GetZaxis()->FindBin(centrality["max"].asDouble() - 1e-6));

      const std::string pTRangeName = DtoStr(pTBin["min"].asDouble(), 1) + "<pT<" + 
                                      DtoStr(pTBin["max"].asDouble(), 1);
      
      if (distrVariableProj->Integral(1, distrVariableProj->GetXaxis()->GetNbins()) < 
          Par.minIntegralValue) 
      {
         PrintInfo("Integral is insufficient for projection of " + variableName 
                   + ", " + detector["name"].asString() + ", " + chargeName + 
                   " at " + zDCRangeName + ", " + centralityRangeName + ", " + pTRangeName);
         continue;
      }
      
      double minX = 0., maxX = -1.;
      for (int i = 1; i <= distrVariableProj->GetXaxis()->GetNbins(); i++)
      {
         if (distrVariableProj->GetBinContent(i) > 1e-7)
         {
            minX = distrVariableProj->GetXaxis()->GetBinLowEdge(i);
            break;
         }
      }
      
      for (int i = distrVariableProj->GetXaxis()->GetNbins(); i > 1; i--)
      {
         if (distrVariableProj->GetBinContent(i) > 1e-7)
         {
            maxX = distrVariableProj->GetXaxis()->GetBinUpEdge(i);
            break;
         }
      }
      
      if (minX > maxX) 
      {
         PrintWarning("Something wrong for projection of " + variableName + 
                      ", " + detector["name"].asString() + ", " + chargeName + 
                      " at " + zDCRangeName + ", " + centralityRangeName + ", " + pTRangeName);
         continue;
      }

      const double maxBinVal = 
         distrVariableProj->GetBinContent(distrVariableProj->GetMaximumBin());

      // main fit; it will be drawn and it's parameters will be extracted for further analysis
      fitFuncDVal.emplace_back(("fitFuncDVal_" + std::to_string(pT)).c_str(), "gaus(0) + gaus(3)");
      fitFuncGaus.emplace_back(("fitGaus_" + std::to_string(pT)).c_str(), "gaus");
      fitFuncBG.emplace_back(("fitBg_" + std::to_string(pT)).c_str(), "gaus");
         
      // scale limits
      fitFuncGaus.back().SetParLimits(0, maxBinVal/2., maxBinVal);
      fitFuncDVal.back().SetParLimits(0, maxBinVal/2., maxBinVal);
      fitFuncDVal.back().SetParLimits(3, maxBinVal/20., maxBinVal);

      // set of alternative fit functions used for uncertainty estimation 
      // by varying ranges of approximation around mean by n*sigma of the main fit
      // for first vector approximation ranges are varied symmetrically within mean
      // for second vector approximation ranges are varied within right of mean; left range is 1sigma
      // for third vector approximation ranges are varied within left of mean; right range is 1sigma
      std::vector<TF1> fitFuncDValAlt, fitFuncDValAltRight, fitFuncDValAltLeft;
      
      for (unsigned long i = 0; i < 4; i++)
      {
         fitFuncDValAlt.emplace_back(("fitFuncDValAlt_" + std::to_string(i) + 
                                         "_" + std::to_string(pT)).c_str(), "gaus(0) + gaus(3)");
         fitFuncDValAlt.back().SetParLimits(0, maxBinVal/2., maxBinVal);
         fitFuncDValAlt.back().SetParLimits(3, maxBinVal/20., maxBinVal);
         fitFuncDValAltRight.emplace_back(("fitFuncDValAltRight_" + std::to_string(i) + 
                                         "_" + std::to_string(pT)).c_str(), "gaus(0) + gaus(3)");
         fitFuncDValAltRight.back().SetParLimits(0, maxBinVal/2., maxBinVal);
         fitFuncDValAltRight.back().SetParLimits(3, maxBinVal/20., maxBinVal);
         fitFuncDValAltLeft.emplace_back(("fitFuncDValAltLeft_" + std::to_string(i) + 
                                         "_" + std::to_string(pT)).c_str(), "gaus(0) + gaus(3)");
         fitFuncDValAltLeft.back().SetParLimits(0, maxBinVal/2., maxBinVal);
         fitFuncDValAltLeft.back().SetParLimits(3, maxBinVal/20., maxBinVal);
      }

      fitFuncGaus.back().SetParameters(1., 0., binWidth*2.);
      fitFuncDVal.back().SetParameters(1., 0., binWidth*2., 1., 0., maxX/2.);
      
      fitFuncGaus.back().SetParLimits(1, minBinX/5., maxBinX/5.);
      fitFuncGaus.back().SetParLimits(2, binWidth, maxBinX/5.);
      fitFuncDVal.back().SetParLimits(1, minX/10., maxX/10.);
      fitFuncDVal.back().SetParLimits(2, binWidth, Average(maxX, maxX, minX));
      
      fitFuncDVal.back().SetParLimits(4, minX*2., maxX*2.);
      fitFuncDVal.back().SetParLimits(5, maxX/3., maxX*3.);
      
      fitFuncDVal.back().SetLineColorAlpha(kRed+1, 0.6);
      fitFuncBG.back().SetLineColorAlpha(kGreen+1, 0.9);
      fitFuncBG.back().SetLineStyle(2);
      fitFuncGaus.back().SetLineColorAlpha(kAzure-3, 0.9);
      fitFuncGaus.back().SetLineStyle(2);

      distrVariableProj->GetXaxis()->SetTitle(Par.variableNameTex[variableBin].c_str());
      distrVariableProj->SetTitle("");
      distrVariableProj->SetTitleSize(0.06, "X");
      distrVariableProj->SetTitleSize(0.06, "Y");
      distrVariableProj->SetLabelSize(0.06, "X");
      distrVariableProj->SetLabelSize(0.06, "Y");
 
      distrVariableProj->GetXaxis()->SetRange(distrVariableProj->GetXaxis()->FindBin(minX+0.01),
                                              distrVariableProj->GetXaxis()->FindBin(maxX-0.01));

      //distrVariableProj->Sumw2();
       
      fitFuncGaus.back().SetRange(minBinX/5., maxBinX/5.);
      
      distrVariableProj->Fit(&fitFuncGaus.back(), "RQMBN");
 
      for (int i = 0; i < 3; i++)
      {
         fitFuncDVal.back().SetParameter(i, fitFuncGaus.back().GetParameter(i));
      }
      
      fitFuncDVal.back().SetRange(minBinX, maxBinX);
      distrVariableProj->Fit(&fitFuncDVal.back(), "RQMBN");

      // fit range
      int fitRangeXMinBin = 
         distrVariableProj->GetXaxis()->FindBin(fitFuncDVal.back().GetParameter(1) - 
                                                fitFuncDVal.back().GetParameter(2)*5.);
      int fitRangeXMaxBin = 
         distrVariableProj->GetXaxis()->FindBin(fitFuncDVal.back().GetParameter(1) + 
                                                fitFuncDVal.back().GetParameter(2)*5.);
      double fitRangeXMin = distrVariableProj->GetXaxis()->GetBinLowEdge(fitRangeXMinBin);
      double fitRangeXMax = distrVariableProj->GetXaxis()->GetBinUpEdge(fitRangeXMaxBin);

      for (unsigned short i = 1; i <= Par.fitNTries; i++)
      {
         fitFuncDVal.back().SetParLimits(0, fitFuncDVal.back().GetParameter(0)/
                                         (1. + 2./static_cast<double>(i*i*i)),
                                         fitFuncDVal.back().GetParameter(0)*
                                         (1. + 2./static_cast<double>(i*i*i)));
         fitFuncDVal.back().SetParLimits(1, fitFuncDVal.back().GetParameter(1)*
                                         (1. - 6./static_cast<double>(i*i*i)),
                                         fitFuncDVal.back().GetParameter(1)*
                                         (1. + 4./static_cast<double>(i*i*i)));
         fitFuncDVal.back().SetParLimits(2, fitFuncDVal.back().GetParameter(2)/
                                         (1. + 5./static_cast<double>(i*i*i)),
                                         fitFuncDVal.back().GetParameter(2)*
                                         (1. + 5./static_cast<double>(i*i*i)));
         fitFuncDVal.back().SetParLimits(3, fitFuncDVal.back().GetParameter(3)/
                                         (1. + 5./static_cast<double>(i*i)),
                                         fitFuncDVal.back().GetParameter(3)*
                                         (1. + 5./static_cast<double>(i*i)));
         fitFuncDVal.back().SetParLimits(4, fitFuncDVal.back().GetParameter(4)*
                                         (1. - 6./static_cast<double>(i*i)),
                                         fitFuncDVal.back().GetParameter(4)*
                                         (1. + 4./static_cast<double>(i*i)));
         fitFuncDVal.back().SetParLimits(5, fitFuncDVal.back().GetParameter(5)/
                                         (1. + 5./static_cast<double>(i*i)),
                                         fitFuncDVal.back().GetParameter(5)*
                                         (1. + 5./static_cast<double>(i*i)));
         
         fitRangeXMinBin = 
            distrVariableProj->GetXaxis()->FindBin(fitFuncDVal.back().GetParameter(1) - 
                                                   fitFuncDVal.back().GetParameter(2)*5.);
         fitRangeXMaxBin = 
            distrVariableProj->GetXaxis()->FindBin(fitFuncDVal.back().GetParameter(1) + 
                                                   fitFuncDVal.back().GetParameter(2)*5.);
         
         fitRangeXMin = distrVariableProj->GetXaxis()->GetBinLowEdge(fitRangeXMinBin);
         fitRangeXMax = distrVariableProj->GetXaxis()->GetBinUpEdge(fitRangeXMaxBin);
            
         fitFuncDVal.back().SetRange(fitRangeXMin, fitRangeXMax);
         distrVariableProj->Fit(&fitFuncDVal.back(), "RQMBNL");
      }

      fitFuncGaus.back().SetRange(fitRangeXMin, fitRangeXMax);
      fitFuncBG.back().SetRange(fitRangeXMin, fitRangeXMax);

      fitRangeXMinBin = 
         distrVariableProj->GetXaxis()->FindBin(fitFuncDVal.back().GetParameter(1) - 
                                                fitFuncDVal.back().GetParameter(2)*10.);
      fitRangeXMaxBin = 
         distrVariableProj->GetXaxis()->FindBin(fitFuncDVal.back().GetParameter(1) + 
                                                fitFuncDVal.back().GetParameter(2)*10.);
      distrVariableProj->GetXaxis()->SetRange(fitRangeXMinBin, fitRangeXMaxBin);
       
      for (int i = 0; i < 3; i++)
      {
         fitFuncGaus.back().SetParameter(i, fitFuncDVal.back().GetParameter(i));
         fitFuncBG.back().SetParameter(i, fitFuncDVal.back().GetParameter(i + 3));
      }
      
      distrVariableProj->SetMarkerStyle(20);
      distrVariableProj->SetMarkerSize(0.7);
      distrVariableProj->SetMarkerColorAlpha(kBlack, 0.8);
      distrVariableProj->SetLineColorAlpha(kBlack, 0.8);
      distrVariableProj->SetMaximum(maxBinVal*1.2);
      
      canvDValVsPT.cd(iCanv);
      
      gPad->SetLeftMargin(0.155);
      gPad->SetBottomMargin(0.128);
      
      gPad->Add(distrVariableProj, "P");
      //Print(fitFuncDVal.back().GetName(), fitFuncBG.back().GetName(), fitFuncGaus.back().GetName());
      fitFuncDVal.back().DrawClone("SAME");
      fitFuncBG.back().DrawClone("SAME");
      fitFuncGaus.back().DrawClone("SAME");
      
      Par.texText.SetText(0.17, 0.85, pTRangeName.c_str());
      gPad->Add(Par.texText.Clone());
      Par.texText.SetText(0.17, 0.79, zDCRangeName.c_str());
      gPad->Add(Par.texText.Clone());
      Par.texText.SetText(0.17, 0.73, chargeName.c_str());
      gPad->Add(Par.texText.Clone());
      Par.texText.SetText(0.17, 0.66, centralityRangeName.c_str());
      gPad->Add(Par.texText.Clone());

      iCanv++;

      for (unsigned long i = 0; i < fitFuncDValAlt.size(); i++)
      {
         fitFuncDValAlt[i].
            SetRange(fitFuncDVal.back().GetParameter(1) - fitFuncDVal.back().GetParameter(2)*
                     static_cast<double>(i + 1)*2., fitFuncDVal.back().GetParameter(1) + 
                     fitFuncDVal.back().GetParameter(2)*static_cast<double>(i + 1)*2.);
         fitFuncDValAltRight[i].
            SetRange(fitFuncDVal.back().GetParameter(1) - fitFuncDVal.back().GetParameter(2), 
                     fitFuncDVal.back().GetParameter(1) + 
                     fitFuncDVal.back().GetParameter(2)*static_cast<double>(i + 1)*2.);
         fitFuncDValAltLeft[i].
            SetRange(fitFuncDVal.back().GetParameter(1) - fitFuncDVal.back().GetParameter(2)*
                     static_cast<double>(i + 1)*2., fitFuncDVal.back().GetParameter(1) + 
                     fitFuncDVal.back().GetParameter(2));

         for (int j = 0; j < fitFuncDVal.back().GetNpar(); j++)
         {
            fitFuncDValAlt[i].SetParameter(j, fitFuncDVal.back().GetParameter(j)); 
            fitFuncDValAltRight[i].SetParameter(j, fitFuncDVal.back().GetParameter(j)); 
            fitFuncDValAltLeft[i].SetParameter(j, fitFuncDVal.back().GetParameter(j)); 
            
            if (j == 0 || j == 3)
            {
               fitFuncDValAlt[i].SetParLimits(j, fitFuncDVal.back().GetParameter(j)/1.2, 
                                              fitFuncDVal.back().GetParameter(j)*1.2); 
               fitFuncDValAltRight[i].SetParLimits(j, fitFuncDVal.back().GetParameter(j)/1.2, 
                                                   fitFuncDVal.back().GetParameter(j)*1.2); 
               fitFuncDValAltLeft[i].SetParLimits(j, fitFuncDVal.back().GetParameter(j)/1.2, 
                                                  fitFuncDVal.back().GetParameter(j)*1.2); 
            }
            else if (j == 2 || j == 4)
            {
               fitFuncDValAlt[i].SetParLimits(j, fitFuncDVal.back().GetParameter(j)/1.5, 
                                              fitFuncDVal.back().GetParameter(j)*1.5); 
               fitFuncDValAltRight[i].SetParLimits(j, fitFuncDVal.back().GetParameter(j)/1.5, 
                                                   fitFuncDVal.back().GetParameter(j)*1.5); 
               fitFuncDValAltLeft[i].SetParLimits(j, fitFuncDVal.back().GetParameter(j)/1.5, 
                                                   fitFuncDVal.back().GetParameter(j)*1.5); 
            }
         }
         distrVariableProj->Fit(&fitFuncDValAlt[i], "RQMBNL");
         distrVariableProj->Fit(&fitFuncDValAltRight[i], "RQMBNL");
         distrVariableProj->Fit(&fitFuncDValAltLeft[i], "RQMBNL");
      }
      
      if (fabs(fitFuncDVal.back().GetParameter(1)) < 
          detector["abs_max_fit_" + variableName].asDouble() && 
          fabs(fitFuncDVal.back().GetParameter(2)) < 
          detector["abs_max_fit_" + variableName].asDouble())
      {
         grMeans.AddPoint(pT, fitFuncDVal.back().GetParameter(1));
         grSigmas.AddPoint(pT, fabs(fitFuncDVal.back().GetParameter(2)));
         
         // relative uncertainties of means are inconsistent since absolute uncertainty does not
         // depend on the position of the mean (it can be close to 0 or much larger number) 
         // therefore for means difference of means divided by sigma is used as uncertainty
         // after CppTools is updated needs to be replaced
         const double meanError = 
            StandardError(fitFuncDVal.back().GetParameter(1),
                          fitFuncDValAlt[0].GetParameter(1),
                          fitFuncDValAlt[1].GetParameter(1),
                          fitFuncDValAlt[2].GetParameter(1),
                          fitFuncDValAlt[3].GetParameter(1), 
                          fitFuncDValAltRight[0].GetParameter(1),
                          fitFuncDValAltRight[1].GetParameter(1),
                          fitFuncDValAltRight[2].GetParameter(1),
                          fitFuncDValAltRight[3].GetParameter(1),
                          fitFuncDValAltLeft[0].GetParameter(1),
                          fitFuncDValAltLeft[1].GetParameter(1),
                          fitFuncDValAltLeft[2].GetParameter(1),
                          fitFuncDValAltLeft[3].GetParameter(1));
         
         const double sigmaError = 
            StandardError(fitFuncDVal.back().GetParameter(2),
                          fitFuncDValAlt[0].GetParameter(2),
                          fitFuncDValAlt[1].GetParameter(2),
                          fitFuncDValAlt[2].GetParameter(2),
                          fitFuncDValAlt[3].GetParameter(2), 
                          fitFuncDValAltRight[0].GetParameter(2),
                          fitFuncDValAltRight[1].GetParameter(2),
                          fitFuncDValAltRight[2].GetParameter(2),
                          fitFuncDValAltRight[3].GetParameter(2),
                          fitFuncDValAltLeft[0].GetParameter(2),
                          fitFuncDValAltLeft[1].GetParameter(2),
                          fitFuncDValAltLeft[2].GetParameter(2),
                          fitFuncDValAltLeft[3].GetParameter(2));
         
         grMeans.SetPointError(grMeans.GetN() - 1, 0, meanError);
         grSigmas.SetPointError(grSigmas.GetN() - 1, 0, sigmaError);
      }
   }
   if (grMeans.GetN() == 0) 
   {
      PrintError("Graph is empty for " + variableName + ", " + 
                 detector["name"].asString() + ", " + chargeName + 
                 " at " + zDCRangeName + ", " + centralityRangeName);
   }
   
   const std::string outputFileNameNoExt = 
      "output/SigmalizedResiduals/" + Par.runName + "/" + detector["name"].asString() + "/" + 
      variableName + "_" + chargeNameShort + 
      centralityRangePathName + zDCRangePathName;
   
   PrintCanvas(&canvDValVsPT, outputFileNameNoExt, false, true);
}

void PBarCall()
{
   if (!Par.showProgress) return;
   while (!Par.isProcessFinished)
   {
      SetNumberOfCalls();
      Par.pBar.Print(static_cast<double>(Par.numberOfCalls)/
                     static_cast<double>(Par.numberOfIterations));
      std::this_thread::sleep_for(std::chrono::milliseconds(200));
   }
   Par.pBar.Print(1.);
};

void SetNumberOfCalls()
{
   if (Par.programMode != 1) return; // Only Mode1 passes
   Par.numberOfCalls = 0;
   for (const auto &file : std::filesystem::directory_iterator("tmp/SigmalizedResiduals/" + 
                                                               Par.runName))
   {
      std::string fileName = static_cast<std::string>(file.path());
      std::ifstream tmpFile(fileName.c_str());
      double currentFileNCalls = 0;
      if (tmpFile >> currentFileNCalls) Par.numberOfCalls += currentFileNCalls;
   }
}

double GetNormRatio(const double ratio)
{
   if (ratio >= 1.) return ratio;
   return 1./ratio;
}

#endif /* SIGMALIZED_RESIDUALS_CPP */
