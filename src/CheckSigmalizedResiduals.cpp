/** 
 *  @file CheckSigmalizedResiduals.cpp 
 *  @brief Same as SigmalizedResiduals but used for estimation of means and sigmas of sdphi and sdz for check of correctness of calibration
 *
 *  This file is a part of a project CalPhenix (https://github.com/Sergeyir/CalPhenix).
 *
 *  @author Sergei Antsupov (antsupov0124@gmail.com)
 **/
#ifndef CHECK_SIGMALIZED_RESIDUALS_CPP
#define CHECK_SIGMALIZED_RESIDUALS_CPP

#include "../include/SigmalizedResiduals.hpp"


// This program works in 2 modes
// Mode2 analyzes track residual variable for the specific detector, variable, and charge
// Mode1 recursively calls the program outside the current process to analyze all configurations
int main(int argc, char **argv)
{
   using namespace SigmalizedResiduals;
   if (argc < 2 || argc > 6 || (argc > 3 && argc < 4)) 
   {
      std::string errMsg = "Expected 1-2 or 3-5 parameters while " + std::to_string(argc - 1) + 
                           " parameter(s) were provided \n";
      errMsg += "Usage: bin/SigmalizedResiduals inputFile numberOfThreads=" + 
                std::to_string(std::thread::hardware_concurrency()) + "*\n";
      errMsg += "Or**: bin/SigmalizedResiduals inputFile detectorBin variableBin ";
      errMsg += "numberOfThreads=" + std::to_string(std::thread::hardware_concurrency()) + "* " +
                "showProgress=true";
      errMsg += "*: default argument is the number of threads on the current machine \n";
      errMsg += "**: this mode analyzes only one configuration";
      CppTools::PrintError(errMsg);
   }

   unsigned int numberOfThreads;

   // initializing ROOT parameters
   ROOT::EnableThreadSafety();
   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   // initializing this program parameters
   inputYAMLCal.OpenFile(argv[1], "sigmalized_residuals");
   inputYAMLCal.CheckStatus("sigmalized_residuals");

   runName = inputYAMLCal["run_name"].as<std::string>();

   // opening input file with parameters of a run
   inputYAMLMain.OpenFile("input/" + runName + "/main.yaml");
   inputYAMLMain.CheckStatus("main");

   if (inputYAMLCal["detectors_to_calibrate"].size() == 0)
   {
      CppTools::PrintInfo("No detectors are specified for calibrations");
      CppTools::PrintInfo("Exiting the program");
      exit(1);
   }

   TDirectory::AddDirectory(kFALSE);

   drawDValDistr = inputYAMLCal["draw_dval_distr"].as<bool>();

   if (argc < 4) // Mode1
   {
      programMode = 1;
      if (argc > 2) numberOfThreads = std::stoi(argv[2]);
      else numberOfThreads = std::thread::hardware_concurrency();
      if (numberOfThreads == 0) CppTools::PrintError("Number of threads must be bigger than 0");

      system("rm -rf tmp/CheckSigmalizedResiduals/*");
      system(("mkdir -p tmp/CheckSigmalizedResiduals/" + runName).c_str());

      numberOfIterations = inputYAMLCal["detectors_to_calibrate"].size()*
                                inputYAMLCal["centrality_bins"].size()*
                                inputYAMLCal["zdc_bins"].size()*4;

      auto SingleThreadCall = [&](const unsigned long detectorBin, 
                                  const unsigned long variableBin)
      {
         // man ROOT sucks (TF1::Fit is still not thread safe) so I have to call the same program 
         // recursively in shell outside of the current instance to implement multithreading
         system((static_cast<std::string>("./bin/CheckSigmalizedResiduals ") + 
                 argv[1] + " " + std::to_string(detectorBin) + " " + 
                 std::to_string(variableBin) + " 1 0").c_str());;
      };

      std::vector<std::thread> thrCalls;
      std::thread pBarThr(PBarCall); 

      for (unsigned int detectorBin = 0; detectorBin < 
           inputYAMLCal["detectors_to_calibrate"].size(); detectorBin++)
      {
         for (unsigned long variableBin = 0; variableBin < variableName.size(); variableBin++)
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

      isProcessFinished = true;
      pBarThr.join();
   }
   else // Mode2
   {
      programMode = 2;
      if (argc > 4) numberOfThreads = std::stoi(argv[4]);
      else numberOfThreads = std::thread::hardware_concurrency();
      if (numberOfThreads == 0) CppTools::PrintError("Number of threads must be bigger than 0");

      ROOT::EnableImplicitMT(numberOfThreads);

      if (argc > 5) showProgress = static_cast<bool>(std::stoi(argv[5]));

      inputFile = 
         std::unique_ptr<TFile>(TFile::Open(("data/SigmalizedResiduals/" + runName + 
                                             "/sum.root").c_str(), "READ"));

      pTRangeTLatex.SetTextFont(52);
      pTRangeTLatex.SetTextSize(0.06);
      pTRangeTLatex.SetNDC();
      zDCRangeTLatex.SetTextFont(52);
      zDCRangeTLatex.SetTextSize(0.06);
      zDCRangeTLatex.SetNDC();
      chargeTLatex.SetTextFont(52);
      chargeTLatex.SetTextSize(0.06);
      chargeTLatex.SetNDC();
      centralityRangeTLatex.SetTextFont(52);
      centralityRangeTLatex.SetTextSize(0.06);
      centralityRangeTLatex.SetNDC();

      numberOfIterations = 2.*inputYAMLCal["centrality_bins"].size()*
                                inputYAMLCal["zdc_bins"].size();

      outputDir = "output/SigmalizedResiduals/" + runName + "/";
      system(("mkdir -p " + outputDir + "CalibrationParameters").c_str());

      pTMin = inputYAMLCal["pt_bins"][0]["min"].as<double>();
      pTMax = inputYAMLCal["pt_bins"][inputYAMLCal["pt_bins"].size() - 1]
                                    ["max"].as<double>();

      fitNTries = inputYAMLCal["number_of_fit_tries"].as<unsigned int>();

      std::thread pBarThr(PBarCall); 

      PerformFitsForDifferentCentrAndZDC(std::stoi(argv[2]), std::stoi(argv[3]));

      isProcessFinished = true;
      pBarThr.join();
   }
 
   return 0;
}

void SigmalizedResiduals::PerformFitsForDifferentCentrAndZDC(const unsigned int detectorBin, 
                                                             const unsigned int variableBin)
{
   const YAML::Node detector = inputYAMLCal["detectors_to_calibrate"][detectorBin];

   const std::string detectorName = detector["name"].as<std::string>();

   system(("mkdir -p " + outputDir + detectorName).c_str());

   outputFile = std::unique_ptr<TFile>
      (TFile::Open((outputDir + detectorName + "/all_fits_s" + 
                   variableName[variableBin] + ".root").c_str(), "RECREATE"));

   for (int charge : particleCharges)
   {
      const std::string chargeName = ((charge > 0) ? "charge>0" : "charge<0");
      const std::string chargeNameShort = ((charge > 0) ? "pos" : "neg");

      // output file in which shifts and scales to sigmalized value is written
      // this is needed to correct sigmalized means and sigmas distributions
      std::ofstream recalOutput(outputDir + "CalibrationParameters/recal_" + 
                                detectorName + "_s" + variableName[variableBin] + 
                                "_" + chargeNameShort + ".txt");

      recalOutput << "1 1" << std::endl;

      const std::string meansFitFunc = 
         detector["means_fit_func_s" + variableName[variableBin] + 
                  "_" + chargeNameShort].as<std::string>();
      const std::string sigmasFitFunc = 
         detector["sigmas_fit_func_s" + variableName[variableBin] + 
                  "_" + chargeNameShort].as<std::string>();

      for (unsigned int centralityBin = 0; centralityBin < 
           inputYAMLCal["centrality_bins"].size(); centralityBin++)
      {
         const YAML::Node centrality = inputYAMLCal["centrality_bins"][centralityBin];

         const std::string centralityRangeName = centrality["min"].as<std::string>() + "-" + 
                                                 centrality["max"].as<std::string>() + "%"; 

         // unique name for this bin to prevent ROOT 
         // printing warnings about replacing same named objects
         const std::string thisBinUniqueName = 
            detectorName + "s" + variableName[variableBin] + chargeName + centralityRangeName;

         // same as before but without percent; used for output file names
         const std::string centralityRangePathName = 
            "_c" + centrality["min"].as<std::string>() + "-" + centrality["max"].as<std::string>();

         outputFile->mkdir((chargeName + "/" + centralityRangePathName).c_str());
         outputFile->cd((chargeName + "/" + centralityRangePathName).c_str());

         std::vector<TGraphErrors> grVMeansVsPT, grVSigmasVsPT;

         for (const YAML::Node& zDC : inputYAMLCal["zdc_bins"])
         { 
            numberOfCalls++;
 
            const std::string zDCRangeName = zDC["min"].as<std::string>() + "<zDC<" + 
                                             zDC["max"].as<std::string>();
            const std::string zDCRangePathName = "_zDC" + zDC["min"].as<std::string>() + 
                                                 "-" + zDC["max"].as<std::string>();

            // name of histogram
            const std::string distrVariableName =  
               "s" + variableName[variableBin] + " vs pT vs centrality: " + 
               detectorName + ", " + chargeName + ", " + zDCRangeName;

            TH3F *distrVariable = 
               static_cast<TH3F *>(inputFile->Get(distrVariableName.c_str()));

            if (!distrVariable) 
            {
               CppTools::PrintError("Histogram named \"" + distrVariableName + 
                                    "\" does not exist in file " + inputFile->GetName());
            }

            std::string fitsOutputFileName = outputDir + detectorName + "/s" + 
                                             variableName[variableBin] + "_" + chargeName + 
                                             centralityRangePathName + zDCRangePathName;

            grVMeansVsPT.emplace_back();
            grVSigmasVsPT.emplace_back();

            PerformFitsForDifferentPT(distrVariable, grVMeansVsPT.back(), grVSigmasVsPT.back(), 
                                      detector, variableBin, zDC, charge, centrality);

            grVMeansVsPT.back().Write(("means: " + zDCRangeName).c_str());
            grVSigmasVsPT.back().Write(("sigmas: " + zDCRangeName).c_str());

            if (!showProgress)
            {
               std::ofstream progressFile("tmp/CheckSigmalizedResiduals/" + runName + "/" + 
                                     std::to_string(detectorBin) + std::to_string(variableBin));
               progressFile << numberOfCalls;
            }

            std::vector<double> grMeansVsPTWeights, grSigmasVsPTWeights;

            for (int i = 0; i < grVMeansVsPT.back().GetN(); i++)
            {
               // these errors will be used instead of actual ones; this is due to the fact
               // that some weights are really big for some points, usually somewhere in the middle.
               // which makes mean of the distribution to be equal to the values in these points
               grMeansVsPTWeights. push_back(1./pow(1. + grVMeansVsPT.back().GetErrorY(i)/
                                                    grVSigmasVsPT.back().GetPointY(i), 2));
               grSigmasVsPTWeights. push_back(1./pow(1. + grVMeansVsPT.back().GetErrorY(i)/
                                                     grVSigmasVsPT.back().GetPointY(i), 2));
            }

            const double meanShift = -TMath::Mean(grVMeansVsPT.back().GetN(), 
                                                  grVMeansVsPT.back().GetY(), 
                                                  &grMeansVsPTWeights[0]);
            const double sigmaScale = 1./TMath::Mean(grVSigmasVsPT.back().GetN(), 
                                                     grVSigmasVsPT.back().GetY(), 
                                                     &grSigmasVsPTWeights[0]);

            recalOutput << meanShift << " " << sigmaScale << std::endl;

            // correct the graphs for the shift to check the effect of the shift
            /*
            for (int i = 0; i < grVMeansVsPT.back().GetN(); i++)
            {
               grVMeansVsPT.back().SetPointY(i, grVMeansVsPT.back().GetPointY(i) + meanShift);
               grVSigmasVsPT.back().SetPointY(i, grVSigmasVsPT.back().GetPointY(i)*sigmaScale);
            }
            */
         }

         for (unsigned long i = 0; i < inputYAMLCal["zdc_bins"].size(); i++)
         {
            const YAML::Node zDC = inputYAMLCal["zdc_bins"][static_cast<int>(i)];

            grVMeansVsPT[i].SetMarkerStyle(zDC["marker_style"].as<short>());
            grVMeansVsPT[i].SetMarkerSize(1.4);
            grVMeansVsPT[i].SetMarkerColorAlpha(zDC["color"].as<short>(), 0.8);
            grVMeansVsPT[i].SetLineColorAlpha(zDC["color"].as<short>(), 0.8);

            grVSigmasVsPT[i].SetMarkerStyle(zDC["marker_style"].as<short>());
            grVSigmasVsPT[i].SetMarkerSize(1.4);
            grVSigmasVsPT[i].SetMarkerColorAlpha(zDC["color"].as<short>(), 0.8);
            grVSigmasVsPT[i].SetLineColorAlpha(zDC["color"].as<short>(), 0.8);
         }


         TCanvas canvValVsPTVsZDC("", "", 800, 800);
         canvValVsPTVsZDC.cd();
 
         TLegend legend{0.15, 0.7, 0.88, 0.89};
         legend.SetNColumns(3);
         legend.SetLineColorAlpha(0, 0.);
         legend.SetFillColorAlpha(0, 0.);

         gPad->SetLeftMargin(0.135);

         TH1F meansFrame(("means frame " + thisBinUniqueName).c_str(), "", 
                         10, pTMin - 0.1, pTMax*1.05);
         meansFrame.SetMinimum(-1.);
         meansFrame.SetMaximum(1.);

         meansFrame.GetXaxis()->SetTitle("p_{T} [GeV/c]");
         meansFrame.GetYaxis()->
            SetTitle(("#mu_{s" + variableNameTex[variableBin] + "}").c_str());
         meansFrame.GetXaxis()->SetTitleOffset(1.1);
         meansFrame.GetYaxis()->SetTitleOffset(2.0);
         gPad->Add(&meansFrame, "AXIS");
         gPad->Add(&meansFrame, "SAME AXIS X+ Y+");

         TLine expectedMeanValLine(pTMin - 0.1, 0., pTMax*1.05, 0.);
         expectedMeanValLine.SetLineColorAlpha(kGray+3, 0.5);
         expectedMeanValLine.SetLineWidth(3);
         expectedMeanValLine.SetLineStyle(2);
         gPad->Add(&expectedMeanValLine);

         for (unsigned long i = 0; i < inputYAMLCal["zdc_bins"].size(); i++)
         {
            const std::string zDCRangeName = 
               inputYAMLCal["zdc_bins"][static_cast<int>(i)]["min"].as<std::string>() + 
               "<z_{DC}<" + 
               inputYAMLCal["zdc_bins"][static_cast<int>(i)]["max"].as<std::string>();

            legend.AddEntry(&grVMeansVsPT[i], zDCRangeName.c_str(), "P");
            gPad->Add(&grVMeansVsPT[i], "SAME P");
         }

         gPad->Add(&legend);

         ROOTTools::PrintCanvas(&canvValVsPTVsZDC, outputDir + detectorName + "/means_s" + 
                                variableName[variableBin] + "_" + chargeNameShort + 
                                centralityRangePathName);

         legend.Clear();
         canvValVsPTVsZDC.Clear();

         TH1F sigmasFrame(("sigmas frame " + thisBinUniqueName).c_str(), "", 
                          10, pTMin - 0.1, pTMax*1.05);
         sigmasFrame.SetMinimum(0.);
         sigmasFrame.SetMaximum(2.);

         sigmasFrame.GetXaxis()->SetTitle("p_{T} [GeV/c]");
         sigmasFrame.GetYaxis()->
            SetTitle(("#sigma_{s" + variableNameTex[variableBin] + "}").c_str());
         gPad->Add(&sigmasFrame, "AXIS");
         gPad->Add(&sigmasFrame, "SAME AXIS X+ Y+");

         TLine expectedSigmaValLine(pTMin - 0.1, 1., pTMax*1.05, 1.);
         expectedSigmaValLine.SetLineColorAlpha(kGray+3, 0.5);
         expectedSigmaValLine.SetLineWidth(3);
         expectedSigmaValLine.SetLineStyle(2);
         gPad->Add(&expectedSigmaValLine);

         for (unsigned long i = 0; i < inputYAMLCal["zdc_bins"].size(); i++)
         { 
            const std::string zDCRangeName = 
               inputYAMLCal["zdc_bins"][static_cast<int>(i)]["min"].as<std::string>() + 
               "<z_{DC}<" + 
               inputYAMLCal["zdc_bins"][static_cast<int>(i)]["max"].as<std::string>();

            legend.AddEntry(&grVSigmasVsPT[i], zDCRangeName.c_str(), "P");
            gPad->Add(&grVSigmasVsPT[i], "SAME P");
         }

         gPad->Add(&legend);

         ROOTTools::PrintCanvas(&canvValVsPTVsZDC, outputDir + detectorName + "/sigmas_s" + 
                                variableName[variableBin] + "_" + chargeNameShort + 
                                centralityRangePathName);
      }

      recalOutput.close();
   }

   outputFile->Close();
}

void SigmalizedResiduals::PerformFitsForDifferentPT(TH3F *hist, TGraphErrors &grMeans, 
                                                    TGraphErrors &grSigmas, 
                                                    const YAML::Node& detector, 
                                                    const unsigned int variableBin, 
                                                    const YAML::Node& zDC, const int charge, 
                                                    const YAML::Node& centrality)
{
   const std::string chargeName = ((charge > 0) ? "charge>0" : "charge<0");
   const std::string chargeNameShort = ((charge > 0) ? "pos" : "neg");

   const std::string centralityRangeName = centrality["min"].as<std::string>() + "-" + 
                                           centrality["max"].as<std::string>() + "%";
   const std::string centralityRangePathName = "_c" + centrality["min"].as<std::string>() + 
                                                 "-" + centrality["max"].as<std::string>();

   const std::string zDCRangeName = zDC["min"].as<std::string>() + "<zDC<" + 
                                    zDC["max"].as<std::string>();
   const std::string zDCRangePathName = "_zDC" + zDC["min"].as<std::string>() + 
                                        "-" + zDC["max"].as<std::string>();

   TCanvas canvDValVsPT(("all fits, " + zDCRangeName).c_str(), "", 
                         inputYAMLCal["pt_nbinsx"].as<double>()*400.,
                         inputYAMLCal["pt_nbinsy"].as<double>()*400.);

   canvDValVsPT.Divide(inputYAMLCal["pt_nbinsx"].as<int>(), 
                       inputYAMLCal["pt_nbinsy"].as<int>());


   // functions for fits of dz and dphi distributions
   // vectors are needed for the object to not be deleted out of scope which 
   // would result in deleting it from canvas
   std::vector<TF1> fitFuncDVal, fitFuncGaus, fitFuncBG;
   // alternate fits in different region for estimation of uncertainty
 
   int iCanv = 1;
 
   for (const YAML::Node& pTBin : inputYAMLCal["pt_bins"])
   {
      const double pT = CppTools::Average(pTBin["min"].as<double>(), pTBin["max"].as<double>());
      if (pT < pTMin || pT > pTMax) continue;
 
      TH1D *distrVariableProj = hist->
         ProjectionX(((std::string) hist->GetName() + "_projX_" + std::to_string(pT)).c_str(), 
                     hist->GetYaxis()->FindBin(pTBin["min"].as<double>() + 1e-6), 
                     hist->GetYaxis()->FindBin(pTBin["max"].as<double>() - 1e-6),
                     hist->GetZaxis()->FindBin(centrality["min"].as<double>() + 1e-6),
                     hist->GetZaxis()->FindBin(centrality["max"].as<double>() - 1e-6));

      const std::string pTRangeName = CppTools::DtoStr(pTBin["min"].as<double>(), 1) + "<pT<" + 
                                      CppTools::DtoStr(pTBin["max"].as<double>(), 1);

      if (distrVariableProj->Integral(1, distrVariableProj->GetXaxis()->GetNbins()) < 
          minIntegralValue) 
      {
         CppTools::PrintInfo("Integral is insufficient for projection of s" + 
                             variableName[variableBin] + ", " + 
                             detector["name"].as<std::string>() + ", " + 
                             chargeName + " at " + zDCRangeName + ", " + 
                             centralityRangeName + ", " + pTRangeName);
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

      fitFuncGaus.back().SetParameters(1., 0., 1.);
      fitFuncDVal.back().SetParameters(1., 0., 1.);

      fitFuncGaus.back().SetParLimits(1, -0.5, 0.5);
      fitFuncGaus.back().SetParLimits(2, 0.5, 2.);
      fitFuncDVal.back().SetParLimits(1, -0.5, 0.5);
      fitFuncDVal.back().SetParLimits(2, 0.5, 2.);

      fitFuncDVal.back().SetParLimits(4, -10., 10.);
      fitFuncDVal.back().SetParLimits(5, 2., 50.);

      fitFuncDVal.back().SetLineColorAlpha(kRed+1, 0.6);
      fitFuncBG.back().SetLineColorAlpha(kGreen+1, 0.9);
      fitFuncBG.back().SetLineStyle(2);
      fitFuncGaus.back().SetLineColorAlpha(kAzure-3, 0.9);
      fitFuncGaus.back().SetLineStyle(2);

      distrVariableProj->GetXaxis()->SetTitle(("s" + variableNameTex[variableBin]).c_str());
      distrVariableProj->SetTitle("");
      distrVariableProj->SetTitleSize(0.06, "X");
      distrVariableProj->SetTitleSize(0.06, "Y");
      distrVariableProj->SetLabelSize(0.06, "X");
      distrVariableProj->SetLabelSize(0.06, "Y");
 
      fitFuncGaus.back().SetRange(-0.5, 0.5);

      distrVariableProj->Fit(&fitFuncGaus.back(), "RQMBN");
 
      for (int i = 0; i < 3; i++)
      {
         fitFuncDVal.back().SetParameter(i, fitFuncGaus.back().GetParameter(i));
      }

      fitFuncDVal.back().SetRange(-5., 5.);
      fitFuncGaus.back().SetRange(-5., 5.);
      fitFuncBG.back().SetRange(-5., 5.);

      distrVariableProj->Fit(&fitFuncDVal.back(), "RQMBN");

      for (unsigned short i = 1; i <= fitNTries; i++)
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

         distrVariableProj->Fit(&fitFuncDVal.back(), "RQMBNL");
      }

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
      fitFuncDVal.back().DrawClone("SAME");
      fitFuncBG.back().DrawClone("SAME");
      fitFuncGaus.back().DrawClone("SAME");

      pTRangeTLatex.SetText(0.17, 0.85, pTRangeName.c_str());
      zDCRangeTLatex.SetText(0.17, 0.79, zDCRangeName.c_str());
      chargeTLatex.SetText(0.17, 0.73, chargeName.c_str());
      centralityRangeTLatex.SetText(0.17, 0.66, centralityRangeName.c_str());
      gPad->Add(&pTRangeTLatex);
      gPad->Add(&zDCRangeTLatex);
      gPad->Add(&chargeTLatex);
      gPad->Add(&centralityRangeTLatex);

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

      if (fabs(fitFuncDVal.back().GetParameter(1)) > 1. ||
          fabs(fitFuncDVal.back().GetParameter(2) - 1.) > 1.) continue;

      grMeans.AddPoint(pT, fitFuncDVal.back().GetParameter(1));
      grSigmas.AddPoint(pT, fabs(fitFuncDVal.back().GetParameter(2)));

      // relative uncertainties of means are inconsistent since absolute uncertainty does not
      // depend on the position of the mean (it can be close to 0 or much larger number) 
      // therefore for means difference of means divided by sigma is used as uncertainty
      // after CppTools is updated needs to be replaced
      const double meanError = 
         CppTools::StandardError(fitFuncDValAlt[0].GetParameter(1),
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
                                 fitFuncDValAltLeft[3].GetParameter(1),
                                 fitFuncDVal.back().GetParameter(1));

      const double sigmaError = 
         CppTools::StandardError(fitFuncDValAlt[0].GetParameter(2),
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
                                 fitFuncDValAltLeft[3].GetParameter(2),
                                 fitFuncDVal.back().GetParameter(2));

      grMeans.SetPointError(grMeans.GetN() - 1, 0, meanError);
      grSigmas.SetPointError(grSigmas.GetN() - 1, 0, sigmaError);
   }

   if (grMeans.GetN() == 0) 
   {
      CppTools::PrintError("Graph is empty for s" + variableName[variableBin] + ", " + 
                           detector["name"].as<std::string>() + ", " + chargeName + 
                           " at " + zDCRangeName + ", " + centralityRangeName);
   }

   canvDValVsPT.Write();

   if (drawDValDistr)
   {
      ROOTTools::PrintCanvas(&canvDValVsPT, outputDir + detector["name"].as<std::string>() + "/s" + 
                                            variableName[variableBin] + "_" + chargeNameShort + 
                                            centralityRangePathName + zDCRangePathName, false);
   }
}

void SigmalizedResiduals::PBarCall()
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

void SigmalizedResiduals::SetNumberOfCalls()
{
   if (programMode != 1) return; // Only Mode1 passes
   numberOfCalls = 0;
   for (const auto &file : 
        std::filesystem::directory_iterator("tmp/CheckSigmalizedResiduals/" + runName))
   {
      std::string fileName = static_cast<std::string>(file.path());
      std::ifstream tmpFile(fileName.c_str());
      double currentFileNCalls = 0;
      if (tmpFile >> currentFileNCalls) numberOfCalls += currentFileNCalls;
   }
}

#endif /* CHECK_SIGMALIZED_RESIDUALS_CPP */
