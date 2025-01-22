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
   if (argc != 2)
   {
      PrintError("Expected 1 parameters while " + std::to_string(argc - 1) + " were provided");
   }
   
   // opening input file with parameters for calibrations
   Par.inputJSONCal.OpenFile(argv[1], "calibration");
   Par.inputJSONCal.CheckStatus("calibration");
   
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

   ROOT::EnableImplicitMT();
   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   
   const double pTMin = calibrationInput["pt_bins"].front()["min"].asDouble();
   const double pTMax = calibrationInput["pt_bins"].back()["max"].asDouble();
   
   ProgressBar pBar("FANCY1");
   pBar.SetColor(PBarColor::BOLD_GREEN);
   
   const unsigned long numberOfIterations = calibrationInput["detectors_to_calibrate"].size()*
                                            calibrationInput["centrality_bins"].size()*
                                            calibrationInput["variables_to_calibrate"].size()*
                                            calibrationInput["particles_to_calibrate"].size()*
                                            calibrationInput["zdc_bins"].size();
   
   unsigned long numberOfCalls = 0.;

   const std::string inputFileName = "data/Real/" + Par.runName + "/SingleTrack/sum.root";

   for (const auto& detector : calibrationInput["detectors_to_calibrate"])
   {
      const std::string outputDir = "output/ResidualCal/" + Par.runName + 
                                    "/" + detector["name"].asString() + "/";
      system(("mkdir -p " + outputDir).c_str());
   }

   TFile inputFile(inputFileName.c_str(), "READ");

   Par.texText.SetTextFont(52);
   Par.texText.SetTextSize(0.06);

   // ranges for TAxis
   std::vector<double> pTRanges, centralityRanges;
   
   for (const Json::Value& pTBin: calibrationInput["pt_bins"])
   {
      pTRanges.push_back(pTBin["min"].asDouble());
   }
   pTRanges.push_back(calibrationInput["pt_bins"].back()["max"].asDouble());

   for (const Json::Value& centralityBin: calibrationInput["centrality_bins"])
   {
      centralityRanges.push_back(centralityBin["min"].asDouble());
   }
   centralityRanges.push_back(calibrationInput["centrality_bins"].back()["max"].asDouble());

   for (const Json::Value& detector : calibrationInput["detectors_to_calibrate"])
   {
      const std::string outputDir = "output/ResidualCal/" + Par.runName + "/";
      const std::string detectorName = detector["name"].asString();
      
      TFile outputFile((outputDir + detectorName + "/all_fits.root").c_str(), "RECREATE");

      for (const Json::Value& variable : calibrationInput["variables_to_calibrate"])
      { 
         for (const Json::Value& particleType : calibrationInput["particles_to_calibrate"])
         {
            TCanvas canvAllCentr("dval means vs centrality and pT", "", 
                                 calibrationInput["pt_nbinsx"].asDouble()*400.,
                                 calibrationInput["pt_nbinsy"].asDouble()*400.);
            canvAllCentr.Divide(calibrationInput["pt_nbinsx"].asInt(), 
                                calibrationInput["pt_nbinsy"].asInt());
            
            std::vector<TH2D> distrVMeansVsCentralityVsPT, distrVSigmasVsCentralityVsPT;
            
            for (const Json::Value& zDCBin : calibrationInput["zdc_bins"])
            {
               const std::string zDCRangeName = zDCBin["min"].asString() + "<zDC<" + 
                                                zDCBin["max"].asString();
               distrVMeansVsCentralityVsPT.
                  emplace_back((variable["name"].asString() + " means: centrality vs pT " + 
                                zDCRangeName).c_str(), "",
                               centralityRanges.size() - 1, &centralityRanges[0],
                               pTRanges.size() - 1, &pTRanges[0]);
               distrVSigmasVsCentralityVsPT.
                  emplace_back((variable["name"].asString() + "sigmas: centrality vs pT" + 
                                zDCRangeName).c_str(), "",
                               centralityRanges.size() - 1, &centralityRanges[0],
                               pTRanges.size() - 1, &pTRanges[0]);
            }
            
            unsigned int iCentralityBin = 0;
            for (const Json::Value& centralityBin : calibrationInput["centrality_bins"])
            {
               const std::string centralityRangeName = centralityBin["min"].asString() + "-" + 
                                                       centralityBin["max"].asString() + "%"; 
               // same as before but without percent; used for output file names
               const std::string centralityRangePathName = 
                  "_c" + centralityBin["min"].asString() + "-" + centralityBin["max"].asString();

               std::vector<TGraphErrors> grVMeansVsPT, grVSigmasVsPT;
               std::vector<TF1> fVMeansVsPT, fVSigmasVsPT;

               // histograms with weights representing means and sigmas
               TH2D distrMeansVsZDCVsPT("means", "#mu", 
                                        calibrationInput["zdc_bins"].size(), 
                                        calibrationInput["zdc_bins"].front()["min"].asDouble(), 
                                        calibrationInput["zdc_bins"].back()["max"].asDouble(),
                                        pTRanges.size() - 1, &pTRanges[0]);
               TH2D distrSigmasVsZDCVsPT("sigmas", "#sigma", 
                                         calibrationInput["zdc_bins"].size(), 
                                         calibrationInput["zdc_bins"].front()["min"].asDouble(), 
                                         calibrationInput["zdc_bins"].back()["max"].asDouble(),
                                         pTRanges.size() - 1, &pTRanges[0]);
               // histograms with weights representing the difference between means and sigmas and
               // a fit parameter of means and sigmas
               TH2D distrMeansDiffVsZDCVsPT("means diff", "#cbar#mu - #mu_{fit}#cbar/#mu",
                                            calibrationInput["zdc_bins"].size(), 
                                            calibrationInput["zdc_bins"].front()["min"].asDouble(), 
                                            calibrationInput["zdc_bins"].back()["max"].asDouble(),
                                            pTRanges.size() - 1, &pTRanges[0]);
               TH2D distrSigmasDiffVsZDCVsPT("sigmas diff", "#cbar#sigma - #sigma_{fit}#cbar/#sigma",
                                             calibrationInput["zdc_bins"].size(), 
                                             calibrationInput["zdc_bins"].front()["min"].asDouble(), 
                                             calibrationInput["zdc_bins"].back()["max"].asDouble(),
                                             pTRanges.size() - 1, &pTRanges[0]);
         
               unsigned int iZDCBin = 0;
               for (const Json::Value& zDCBin : calibrationInput["zdc_bins"])
               { 
                  pBar.Print(static_cast<double>(numberOfCalls)/
                             static_cast<double>(numberOfIterations));
                  numberOfCalls++;
                   
                  const std::string zDCRangeName = zDCBin["min"].asString() + "<zDC<" + 
                                                   zDCBin["max"].asString();
                  const std::string zDCRangePathName = "_zDC" + zDCBin["min"].asString() + 
                                                       "-" + zDCBin["max"].asString();
                  
                  const double zDCMin = zDCBin["min"].asDouble();
                  const double zDCMax = zDCBin["max"].asDouble();
                  
                  // name of histogram
                  const std::string distrVariableName =  
                     variable["name"].asString() + " vs pT vs centrality: " + 
                     detector["name"].asString() + ", " + particleType["name"].asString() + ", " + 
                     zDCRangeName;
                  
                  TH3F *distrVariable = 
                     static_cast<TH3F *>(inputFile.Get(distrVariableName.c_str()));
                  
                  if (!distrVariable) PrintError("Histogram named \"" + distrVariableName + 
                                                 "\" does not exist in file " + inputFileName);
                  
                  std::string fitsOutputFileName = outputDir + detector["name"].asString() + "/" + 
                                                   variable["name"].asString() + "_" +
                                                   particleType["name_short"].asString() +
                                                   centralityRangePathName + zDCRangePathName;

                  grVMeansVsPT.emplace_back();
                  grVSigmasVsPT.emplace_back();
                  
                  fVMeansVsPT.emplace_back((zDCRangeName + centralityRangeName + 
                                        detector["name"].asString() + 
                                        particleType["name"].asString() + 
                                        variable["name"].asString()).c_str(), 
                                       detector["means_fit"].asCString());
                  fVSigmasVsPT.emplace_back((zDCRangeName + centralityRangeName + 
                                         detector["name"].asString() + 
                                         particleType["name"].asString() +
                                         variable["name"].asString()).c_str(), 
                                        detector["sigmas_fit"].asCString());
                                 
                  PerformFits(distrVariable, grVMeansVsPT.back(), grVSigmasVsPT.back(), 
                              calibrationInput, detector, variable, 
                              zDCBin, particleType, centralityBin);
                  
                  fVMeansVsPT.back().SetRange(pTMin/1.05, pTMax*1.05);
                  fVSigmasVsPT.back().SetRange(pTMin/1.05, pTMax*1.05);
                  
                  for (short j = 1; j <= Par.fitNTries; j++)
                  {
                     grVMeansVsPT.back().Fit(&fVMeansVsPT.back(), "RQMN");
                     grVSigmasVsPT.back().Fit(&fVSigmasVsPT.back(), "RQMN");
                     
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

                  grVMeansVsPT.back().Clone()->
                     Write(("Data means: " + variable["name"].asString() + 
                            ", " + particleType["name"].asString() + ", " + 
                            centralityRangeName + ", " + zDCRangeName).c_str());
                  
                  grVSigmasVsPT.back().Clone()->
                     Write(("Data sigmas: " + variable["name"].asString() + 
                            ", " + particleType["name"].asString() + ", " + 
                            centralityRangeName + ", " + zDCRangeName).c_str());
                  
                  fVMeansVsPT.back().Clone()->
                     Write(("Fit means: " + variable["name"].asString() + 
                            ", " + particleType["name"].asString() + ", " + 
                            centralityRangeName + ", " + zDCRangeName).c_str());
                  
                  fVSigmasVsPT.back().Clone()->
                     Write(("Fit sigmas: " + variable["name"].asString() + 
                            ", " + particleType["name"].asString() + ", " + 
                            centralityRangeName + ", " + zDCRangeName).c_str());
                  iZDCBin++;
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

                  const Json::Value zDCBin = calibrationInput["zdc_bins"][static_cast<int>(i)];
                  
                  grVMeansVsPT[i].SetMarkerStyle(zDCBin["marker_style"].asInt());
                  grVMeansVsPT[i].SetMarkerSize(1.4);
                  grVMeansVsPT[i].SetMarkerColorAlpha(zDCBin["color"].asInt(), 0.8);
                  grVMeansVsPT[i].SetLineColorAlpha(zDCBin["color"].asInt(), 0.8);
                  fVMeansVsPT[i].SetLineColorAlpha(zDCBin["color"].asInt(), 0.9);
                  fVMeansVsPT[i].SetLineStyle(3);
                  
                  grVSigmasVsPT[i].SetMarkerStyle(zDCBin["marker_style"].asInt());
                  grVSigmasVsPT[i].SetMarkerSize(1.4);
                  grVSigmasVsPT[i].SetMarkerColorAlpha(zDCBin["color"].asInt(), 0.8);
                  grVSigmasVsPT[i].SetLineColorAlpha(zDCBin["color"].asInt(), 0.8);
                  fVSigmasVsPT[i].SetLineColorAlpha(zDCBin["color"].asInt(), 0.9);
                  fVSigmasVsPT[i].SetLineStyle(2);

                  for (int j = 0; j < grVMeansVsPT[i].GetN(); j++)
                  {
                     distrVMeansVsCentralityVsPT[i].
                        SetBinContent(iCentralityBin + 1, 
                                      distrVMeansVsCentralityVsPT[i].GetYaxis()->
                                      FindBin(grVMeansVsPT[i].GetPointX(j)), 
                                      grVMeansVsPT[i].GetPointY(j));
                     distrVSigmasVsCentralityVsPT[i].
                        SetBinContent(iCentralityBin + 1, 
                                      distrVSigmasVsCentralityVsPT[i].GetYaxis()->
                                      FindBin(grVSigmasVsPT[i].GetPointX(j)), 
                                      grVSigmasVsPT[i].GetPointY(j));
                  }
               }
               
               TCanvas canv("", "", 800, 800);
               canv.cd();
               
               TLegend legend{0.15, 0.7, 0.88, 0.89};
               legend.SetNColumns(3);
               legend.SetLineColorAlpha(0, 0.);
               legend.SetFillColorAlpha(0, 0.);
               
               gPad->SetLeftMargin(0.135);
               
               TH1 *meansFrame = 
                  gPad->DrawFrame(pTMin - 0.1, meanYMin - (meanYMax - meanYMin)*0.05, 
                                  pTMax*1.05, meanYMax + (meanYMax - meanYMin)*0.35);
               
               meansFrame->GetXaxis()->SetTitle("p_{T} [GeV/c]");
               meansFrame->GetYaxis()->
                  SetTitle(("#mu_{" + variable["tex_name"].asString() + "}").c_str());
               meansFrame->GetXaxis()->SetTitleOffset(1.1);
               meansFrame->GetYaxis()->SetTitleOffset(2.0);
               meansFrame->Draw("SAME AXIS X+ Y+");
               
               for (unsigned long i = 0; i < calibrationInput["zdc_bins"].size(); i++)
               {
                  const std::string zDCRangeName = 
                     calibrationInput["zdc_bins"][static_cast<int>(i)]["min"].asString() + 
                     "<z_{DC}<" + 
                     calibrationInput["zdc_bins"][static_cast<int>(i)]["max"].asString();
                  
                  legend.AddEntry(&grVMeansVsPT[i], zDCRangeName.c_str(), "P");
                  grVMeansVsPT[i].Clone()->Draw("SAME P");
                  fVMeansVsPT[i].Clone()->Draw("SAME");
               }

               legend.DrawClone();
               PrintCanvas(&canv, outputDir + detector["name"].asString() + "_means_" + 
                           variable["name"].asString() + "_" + 
                           particleType["name_short"].asString() +
                           centralityRangePathName);

               canv.Clone()->Write((variable["name"].asString() + " means: " + 
                                    particleType["name"].asString() + ", " + 
                                    centralityRangeName).c_str());
               
               legend.Clear();
               canv.Clear();

               TH1 *sigmasFrame = 
                  gPad->DrawFrame(pTMin - 0.1, sigmaYMin/1.1, 
                                  pTMax*1.05, sigmaYMax*1.4);
               sigmasFrame->GetXaxis()->SetTitle("p_{T} [GeV/c]");
               sigmasFrame->GetYaxis()->
                  SetTitle(("#sigma_{" + variable["name"].asString() + "}").c_str());
               sigmasFrame->Draw("SAME AXIS X+ Y+");
               
               for (unsigned long i = 0; i < calibrationInput["zdc_bins"].size(); i++)
               { 
                  const std::string zDCRangeName = 
                     calibrationInput["zdc_bins"][static_cast<int>(i)]["min"].asString() + 
                     "<z_{DC}<" + 
                     calibrationInput["zdc_bins"][static_cast<int>(i)]["max"].asString();
                  
                  legend.AddEntry(&grVSigmasVsPT[i], zDCRangeName.c_str(), "P");
                  grVSigmasVsPT[i].Clone()->Draw("SAME P");
                  fVSigmasVsPT[i].Clone()->Draw("SAME");
               }

               legend.DrawClone();
               PrintCanvas(&canv, outputDir + detector["name"].asString() + "_sigmas_" + 
                           variable["name"].asString() + "_" + 
                           particleType["name_short"].asString() +
                           centralityRangePathName);

               canv.Clone()->Write((variable["name"].asString() + " sigmas: " + 
                                    particleType["name"].asString() + ", " + 
                                    centralityRangeName).c_str());

               TCanvas parCanv("", "", 800, 800);
               parCanv.Divide(2, 2);
               
               parCanv.cd(1);
               gPad->SetRightMargin(0.13);
               distrMeansVsZDCVsPT.GetXaxis()->SetTitle(variable["tex_name"].asCString());
               distrMeansVsZDCVsPT.GetYaxis()->SetTitle("p_{T}");
               distrMeansVsZDCVsPT.Draw("COLZ");

               parCanv.cd(2);
               gPad->SetRightMargin(0.13);
               distrSigmasVsZDCVsPT.GetXaxis()->SetTitle(variable["tex_name"].asCString());
               distrSigmasVsZDCVsPT.GetYaxis()->SetTitle("p_{T}");
               distrSigmasVsZDCVsPT.Draw("COLZ");

               parCanv.cd(3);
               gPad->SetLogz();
               gPad->SetRightMargin(0.13);
               distrMeansDiffVsZDCVsPT.GetXaxis()->SetTitle(variable["tex_name"].asCString());
               distrMeansDiffVsZDCVsPT.GetYaxis()->SetTitle("p_{T}");
               distrMeansDiffVsZDCVsPT.Draw("COLZ");

               parCanv.cd(4);
               gPad->SetLogz();
               gPad->SetRightMargin(0.13);
               distrSigmasDiffVsZDCVsPT.GetXaxis()->SetTitle(variable["tex_name"].asCString());
               distrSigmasDiffVsZDCVsPT.GetYaxis()->SetTitle("p_{T}");
               distrSigmasDiffVsZDCVsPT.Draw("COLZ");
               
               PrintCanvas(&parCanv, outputDir + detectorName + 
                           "/fitPar_" + variable["name"].asString() + "_" + 
                           particleType["name_short"].asString() + 
                           centralityRangePathName);
               
               iCentralityBin++;
            }

            for (unsigned long i = 0; i < distrVMeansVsCentralityVsPT.size(); i++)
            {
               distrVMeansVsCentralityVsPT[i].Write();
               distrVSigmasVsCentralityVsPT[i].Write();
            }
         }
      }
      outputFile.Close();
      //break;
   }
   pBar.Print(1.);
}

void PerformFits(TH3F *hist, TGraphErrors &grMeans, TGraphErrors &grSigmas,
                 const Json::Value& calibrationInput,
                 const Json::Value& detector, const Json::Value& variable,
                 const Json::Value& zDCBin, const Json::Value& particleType,
                 const Json::Value& centralityBin)
{ 
   const double minBinX = hist->GetXaxis()->GetBinLowEdge(1);
   const double maxBinX = hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->GetNbins());
   const double binWidth = hist->GetXaxis()->GetBinWidth(1);

   const std::string centralityRangeName = centralityBin["min"].asString() + "-" + 
                                           centralityBin["max"].asString() + "%";
   const std::string centralityRangePathName = "_c" + centralityBin["min"].asString() + 
                                                 "-" + centralityBin["max"].asString();

   const std::string zDCRangeName = zDCBin["min"].asString() + "<zDC<" + zDCBin["max"].asString();
   const std::string zDCRangePathName = "_zDC" + zDCBin["min"].asString() + 
                                        "-" + zDCBin["max"].asString();

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
         
         canv.cd(iCanv);
         iCanv++;
      
         TH1D *distrVariableProj = hist->
            ProjectionX(((std::string) hist->GetName() + "_projX").c_str(), 
                        hist->GetYaxis()->FindBin(pTBin["min"].asDouble() + 1e-6), 
                        hist->GetYaxis()->FindBin(pTBin["max"].asDouble() - 1e-6),
                        hist->GetZaxis()->FindBin(centralityBin["min"].asDouble() + 1e-6),
                        hist->GetZaxis()->FindBin(centralityBin["max"].asDouble() - 1e-6));

         const std::string pTRangeName = DtoStr(pTBin["min"].asDouble(), 1) + "<pT<" + 
                                         DtoStr(pTBin["max"].asDouble(), 1);
         
         if (distrVariableProj->Integral(1, distrVariableProj->GetXaxis()->GetNbins()) < 
             Par.minIntegralValue) 
         {
            PrintInfo("Integral is insufficient for projection of " + variable["name"].asString() + 
                      ", " + detector["name"].asString() + ", " + particleType["name"].asString() + 
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
                         detector["name"].asString() + ", " + particleType["name"].asString() + 
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
         
         distrVariableProj->Fit(&fitFuncGaus, "RQMBN");

         fitFuncGaus.SetRange(minBinX, maxBinX);
         fitFuncBG.SetRange(minBinX, maxBinX);
         
         for (int k = 0; k < 3; k++)
         {
            fitFuncDVal.SetParameter(k, fitFuncGaus.GetParameter(k));
         }
         
         fitFuncDVal.SetRange(minBinX, maxBinX);
         distrVariableProj->Fit(&fitFuncDVal, "RQMBN");

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
            
            distrVariableProj->Fit(&fitFuncDVal, "RQMBN");
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
         
         gPad->SetLeftMargin(0.155);
         gPad->SetBottomMargin(0.115);
         
         distrVariableProj->Clone()->Draw("P");
         fitFuncDVal.Clone()->Draw("SAME");
         fitFuncBG.Clone()->Draw("SAME");
         fitFuncGaus.Clone()->Draw("SAME");
         
         Par.texText.DrawLatexNDC(0.17, 0.85, pTRangeName.c_str());
         Par.texText.DrawLatexNDC(0.17, 0.79, zDCRangeName.c_str());
         Par.texText.DrawLatexNDC(0.17, 0.73, particleType["name"].asCString());
         Par.texText.DrawLatexNDC(0.17, 0.66, centralityRangeName.c_str());
         
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
                    detector["name"].asString() + ", " + particleType["name"].asString() + 
                    " at " + zDCRangeName + ", " + centralityRangeName);
      }
      
      meansFit->SetRange(pTMin, pTMax);
      sigmasFit->SetRange(pTMin, pTMax);
      
      grMeans.Fit(meansFit, "RQMNW");
      grSigmas.Fit(sigmasFit, "RQMNW");
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
      "output/ResidualCal/" + Par.runName + "/" + detector["name"].asString() + "/" + 
      variable["name"].asString() + "_" + particleType["name_short"].asString() + 
      centralityRangePathName + zDCRangePathName;
   PrintCanvas(&canv, outputFileNameNoExt, false, true);
}

#endif /* SIGMALIZED_RESIDUALS_CPP */
