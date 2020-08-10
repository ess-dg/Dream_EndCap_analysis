//   This fitting program was written by Oleg Ivanov
//   modified and updated by Irina Stefanescu
//      How to use:
//   1. Load the program in root by typing .L FitPeaks.C - the code will be compiled and all the functions loaded to memory
//   2. Call the function of choice, e.g., FitSinglePeaks(...)

#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"

Double_t gBgConstant, gBgSlope, gContent, gMean, gContent_1, gMean_1, gContent_2, 
        gMean_2, gSigma, gSigma_1, gSigma_2, gBinW, gChi2pNDF;

Double_t gExp1, gExp2;

//********************

Double_t gaus(Double_t *x, Double_t *par)
{
/*
  par[0]   gauss width
  par[1]   gauss0 constant
  par[2]   gauss0 mean
*/
   static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
   Double_t arg;
   if (par[0] == 0) par[0]=1;                 //  force widths /= 0
   arg = (x[0] - par[2])/(sqrt2*par[0]);
   Double_t fitval = gBinW/(sqrt2pi*par[0]) * par[1] * exp(-arg*arg);
   return fitval;
}
//_____________________________________________________

Double_t gaus_lbg(Double_t *x, Double_t *par)
{
/*
  par[0]   background constant
  par[1]   background slope
  par[2]   gauss width
  par[3]   gauss0 constant
  par[4]   gauss0 mean
*/
   static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
   Double_t arg;
   if (par[2] == 0) par[2]=1;                 //  force widths /= 0
   arg = (x[0] - par[4])/(sqrt2*par[2]);
   Double_t fitval = par[0] + x[0]*par[1]
              + gBinW/(sqrt2pi*par[2]) * par[3] * exp(-arg*arg);
   return fitval;
}
//_____________________________________________________

Double_t db_gaus_lbg_diff(Double_t *x, Double_t *par)
{
/*
  par[0]   background constant
  par[1]   background slope
  par[2]   gauss0 width
  par[3]   gauss0 constant
  par[4]   gauss0 mean
  par[5]   gauss1 width
  par[6]   gauss1 constant
  par[7]   gauss1 mean
*/
   static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
   Double_t arg_1, arg_2;
   if (par[2] == 0) par[2]=1;                 //  force widths /= 0
   arg_1 = (x[0] - par[4])/(sqrt2*par[2]);
   arg_2 = (x[0] - par[7])/(sqrt2*par[5]);
   Double_t fitval = par[0] + x[0]*par[1]
                   + gBinW/(sqrt2pi*par[2]) * par[3] * exp(-arg_1*arg_1)
                   + gBinW/(sqrt2pi*par[5]) * par[6] * exp(-arg_2*arg_2);
   return fitval;
}

//_____________________________________________________

Double_t lorentz(Double_t *x, Double_t *par)
{
/*
  par[0]   background constant
  par[1]   background slope
  par[2]   lorentz area
  par[3]   lorentz FWHM
  par[4]   lorentz mean
*/
   if (par[3] == 0) par[3]=1;                 //  force widths /= 0
   
// function is area normalised 

   Double_t fitval = par[0] + x[0]*par[1] + (0.5*par[2]*par[3]/TMath::Pi()) /
                     TMath::Max( 1.e-10,(x[0]-par[4])*(x[0]-par[4])
                     + .25*par[3]*par[3])*gBinW;
   
   return fitval;
}
//_____________________________________________________

Double_t gaus_lorentz_lbg(Double_t *x, Double_t *par)
{
/*
  par[0]   background constant
  par[1]   background slope
  par[2]   gauss width
  par[3]   gauss0 area
  par[4]   gauss0 mean
  par[5]   lorentz width
  par[6]   lorentz mean
  par[7]   lorentz height
  par[8]   gaus/lorentz mixing
*/
   static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
   Double_t arg,arg1,arg2;
   if (par[2] == 0) par[2]=1;                 //  force gauss widths /= 0
   if (par[5] == 0) par[5]=1;                //  force lorentz width /= 0
   Double_t m = 0.2;               
   
   arg = (x[0] - par[4])/(sqrt2*par[2]);  // gauss
   
   arg1 = 16*(x[0] - par[6])*(x[0] - par[6]);  // lorentz
   arg2 = par[5] * par[5];   //lorentz

// sum of 3 functions: linear bkg +  Gauss normalised to area + lorentz normalized to area

   Double_t fitval = par[0] + x[0]*par[1]
              + m*gBinW/(sqrt2pi*par[2]) * par[3] * exp(-arg*arg)
              + (1-m)*arg2/(arg2 + arg1)/TMath::Pi()/par[5];
   return fitval;
}
//_____________________________________________________


Double_t db_gaus_lbg(Double_t *x, Double_t *par)
{
/*
  par[0]   background constant
  par[1]   background slope
  par[2]   gauss width
  par[3]   gauss0 constant
  par[4]   gauss0 mean
  par[5]   gauss1 constant
  par[6]   gauss1 mean
*/
   static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
   Double_t arg_1, arg_2;
   if (par[2] == 0) par[2]=1;                 //  force widths /= 0
   arg_1 = (x[0] - par[4])/(sqrt2*par[2]);
   arg_2 = (x[0] - par[6])/(sqrt2*par[2]);
   Double_t fitval = par[0] + x[0]*par[1]
                   + gBinW/(sqrt2pi*par[2]) * par[3] * exp(-arg_1*arg_1)
                   + gBinW/(sqrt2pi*par[2]) * par[5] * exp(-arg_2*arg_2);
   return fitval;
}
//_____________________________________________________

// Written by Irina for use in Heimdal NoCath
// inspired from the
// https://doxygen.mantidproject.org/nightly/d5/d56/BackToBackExponential_8cpp_source.html


Double_t db_exp_lbg(Double_t *x, Double_t *par)
{
/*
  par[0]   background constant
  par[1]   background slope
  par[2]   integrated intensity of the peak (area)
  par[3]   exp1 constant of the rising part of neutron pulse  (a)
  par[4]   exp2 constant of decaying part of neutron pulse    (b)
  par[5]   peak position
  par[6]   standard deviation of gaussian part of peakshape function   (s)
*/
   static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
   
   Double_t arg,arg1,arg2,arg3,arg4;

//   arg = par[3]*par[4]/2/(par[3] + par[4]);   //a*b/[2*(a+b)], normalisation factor
   
   if (arg ==0) arg = 1;
   
   arg1 = par[3]/2 *(par[3]*par[6]*par[6] + 2*(x[0]-par[5]));
   arg2 = par[4]/2 *(par[4]*par[6]*par[6] - 2*(x[0]-par[5]));
   arg3 = par[3]*par[6]*par[6] + (x[0]-par[5])/TMath::Sqrt(2*par[6]*par[6]);
   arg4 = par[4]*par[6]*par[6] - (x[0]-par[5])/TMath::Sqrt(2*par[6]*par[6]);
   
//  Double_t val =  exp(arg1 + gsl_sf_log_erfc(arg3))
//                + exp(arg2 + gsl_sf_log_erfc(arg4));
                        
  Double_t val =  exp(arg1 + log(ROOT::Math::erfc(arg3)))
                + exp(arg2 + log(ROOT::Math::erfc(arg4)));

   Double_t fitval = par[0] + x[0]*par[1] + gBinW * par[2] * val * arg;  
                      
   return fitval;
}
//_____________________________________________________

void FitSinglePeak(TH1F *hist, Double_t gLowX, Double_t gUpX)
{
 if(gLowX < gUpX)
  { 
// *** Creating the function of the form 'gaus_lbg' defined above ***
   TF1 fitfunc("gauss_linbg",gaus_lbg, 0, 1, 5);
// *** Obtaining and specifying the start values for the fit ***
   gContent = hist->Integral(hist->FindBin(gLowX),hist->FindBin(gUpX)); 
   gMean    = 0.5 * ( gLowX + gUpX);  
   gSigma   = 0.3 * ( gUpX  - gLowX); 
   gBinW = hist->GetBinWidth(1);
//   printf("__________________\n_The Start Values_\n  Bin Width: %d\n Mean Value: %d\n    Content: %d\n      Sigma: %d\n__________________\n",gBinW,gMean,gContent,gSigma);
   fitfunc.SetParameters(0, 0, gSigma, gContent, gMean); 
   fitfunc.SetRange(gLowX, gUpX);
  
   fitfunc.SetParName(0,"BgConstant");
   fitfunc.SetParName(1,"BgSlope   ");
   fitfunc.SetParName(2,"Sigma     ");
   fitfunc.SetParName(3,"Content   ");
   fitfunc.SetParName(4,"Mean      ");

// *** Fitting: 'R' means within the range specified above ***
   hist->Fit("gauss_linbg", "R", "SAME");

   gBgConstant = fitfunc.GetParameter(0);
   gBgSlope    = fitfunc.GetParameter(1);
   gSigma      = fitfunc.GetParameter(2);
   gContent    = fitfunc.GetParameter(3);
   gMean       = fitfunc.GetParameter(4);
   gChi2pNDF   = fitfunc.GetChisquare() / fitfunc.GetNDF();

   printf("      Chi Square: %f\n",fitfunc.GetChisquare());
   printf("            FWHM: %f +- %f\n",2*gSigma*sqrt(2*log(2)),2*sqrt(2*log(2))*fitfunc.GetParError(2));
  } else cout << "Couldn't fit! Error: The Lower Limit is larger than the Upper Limit!" << endl;
}
//_____________________________________________________


void FitSinglePeak(TH1D *hist, Double_t gLowX, Double_t gUpX)
{
 if(gLowX < gUpX)
  { 
// *** Creating the function of the form 'gaus_lbg' defined above ***
   TF1 fitfunc("gauss_linbg",gaus_lbg, 0, 1, 5);
// *** Obtaining and specifying the start values for the fit ***
   gContent = hist->Integral(hist->FindBin(gLowX),hist->FindBin(gUpX)); 
   gMean    = 0.5 * ( gLowX + gUpX);  
   gSigma   = 0.3 * ( gUpX  - gLowX); 
   gBinW = hist->GetBinWidth(1);
//   printf("__________________\n_The Start Values_\n  Bin Width: %d\n Mean Value: %d\n    Content: %d\n      Sigma: %d\n__________________\n",gBinW,gMean,gContent,gSigma);
   fitfunc.SetParameters(0, 0, gSigma, gContent, gMean); 
   fitfunc.SetRange(gLowX, gUpX);
  
   fitfunc.SetParName(0,"BgConstant");
   fitfunc.SetParName(1,"BgSlope   ");
   fitfunc.SetParName(2,"Sigma     ");
   fitfunc.SetParName(3,"Content   ");
   fitfunc.SetParName(4,"Mean      ");

// *** Fitting: 'R' means within the range specified above ***
   hist->Fit("gauss_linbg", "R", "SAME");

   gBgConstant = fitfunc.GetParameter(0);
   gBgSlope    = fitfunc.GetParameter(1);
   gSigma      = fitfunc.GetParameter(2);
   gContent    = fitfunc.GetParameter(3);
   gMean       = fitfunc.GetParameter(4);
   gChi2pNDF   = fitfunc.GetChisquare() / fitfunc.GetNDF();
   
   gStyle->SetOptFit(1);

   printf("      Chi Square: %f\n",fitfunc.GetChisquare());
   printf("            FWHM: %f +- %f\n",2*gSigma*sqrt(2*log(2)),2*sqrt(2*log(2))*fitfunc.GetParError(2));
  } else cout << "Couldn't fit! Error: The Lower Limit is larger than the Upper Limit!" << endl;
}

void FitDoublePeak(TH1D *hist, Double_t gLowX_1, Double_t gUpX_1, Double_t gLowX_2, Double_t gUpX_2)
{
 if(gLowX_1 < gUpX_1)
  {
 if(gLowX_2 < gUpX_2)
  {
// *** Creating the function of the form '2_gaus_lbg' defined above ***
   TF1 fitfunc("db_gauss_linbg_diff",db_gaus_lbg_diff, 0, 1, 8);
// *** Obtaining and specifying the start values for the fit ***
   gBinW      = hist->GetBinWidth(1);
   gContent_1 = gBinW*(hist->Integral(hist->FindBin(gLowX_1),hist->FindBin(gUpX_1)));
   gContent_2 = gBinW*(hist->Integral(hist->FindBin(gLowX_2),hist->FindBin(gUpX_2)));
// *** Searching for maximum Y value through the bins specified by limits
   int i, i_1, i_2;
   Double_t V, V_max;
   
   i_1 = int(hist->FindBin(gLowX_1));
   i_2 = int(hist->FindBin(gUpX_1));
   V_max = hist->GetBinContent(i_1);
   for(i = i_1; i <= i_2; i++)
    {
     V = hist->GetBinContent(i);
     if(V > V_max)
      {
       V_max = V;
       gMean_1 = double(i);
      }
    }
   gMean_1 = gBinW*gMean_1; 

   i_1 = int(hist->FindBin(gLowX_2));
   i_2 = int(hist->FindBin(gUpX_2));
   V_max = hist->GetBinContent(i_1);
   for(i = i_1; i <= i_2; i++)
    {
     V = hist->GetBinContent(i);
     if(V > V_max)
      {
       V_max = V;
       gMean_2 = double(i);
      }
    }
   gMean_2 = gBinW*gMean_2; 

   gSigma_1 = 0.3 * (gUpX_1 - gLowX_1);
   gSigma_2 = 0.3 * (gUpX_2 - gLowX_2);
   
   printf("__________________\n_Peak 1: The Start Values_\n  Bin Width: %f\n Mean Value: %f\n    Content: %f\n      Sigma: %f\n__________________\n",gBinW,gMean_1,gContent_1,gSigma);
   printf("__________________\n_Peak 2: The Start Values_\n  Bin Width: %f\n Mean Value: %f\n    Content: %f\n      Sigma: %f\n__________________\n",gBinW,gMean_2,gContent_2,gSigma);
   fitfunc.SetParameters(0, 0, gSigma_1, gContent_1, gMean_1, gSigma_2, gContent_2, gMean_2);
   fitfunc.SetRange(gLowX_1, gUpX_2);
  
   fitfunc.SetParName(0,"BgConstant");
   fitfunc.SetParName(1,"BgSlope   ");
   fitfunc.SetParName(2,"Sigma 1   ");
   fitfunc.SetParName(3,"Content 1 ");
   fitfunc.SetParName(4,"Mean 1    ");
   fitfunc.SetParName(5,"Sigma 2   ");
   fitfunc.SetParName(6,"Content 2 ");
   fitfunc.SetParName(7,"Mean 2    ");

// *** Fitting: 'R' means within the range specified above ***
   hist->Fit("db_gauss_linbg_diff", "R", "SAME");

   gBgConstant = fitfunc.GetParameter(0);
   gBgSlope    = fitfunc.GetParameter(1);
   gSigma_1    = fitfunc.GetParameter(2);
   gContent_1  = fitfunc.GetParameter(3);
   gMean_1     = fitfunc.GetParameter(4);
   gSigma_2    = fitfunc.GetParameter(5);
   gContent_2  = fitfunc.GetParameter(6);
   gMean_2     = fitfunc.GetParameter(7);
   gChi2pNDF   = fitfunc.GetChisquare() / fitfunc.GetNDF();

   printf("      Chi Square: %f\n",fitfunc.GetChisquare());
   printf("            FWHM 1: %f +- %f\n",2*gSigma_1*sqrt(2*log(2)),2*sqrt(2*log(2))*fitfunc.GetParError(2));
   printf("            FWHM 2: %f +- %f\n",2*gSigma_2*sqrt(2*log(2)),2*sqrt(2*log(2))*fitfunc.GetParError(5));
  } else cout << "Couldn't fit! Error: Peak 2: The Lower Limit is larger than the Upper Limit!" << endl;
  } else cout << "Couldn't fit! Error: Peak 1: The Lower Limit is larger than the Upper Limit!" << endl;
}

void FitDoublePeak(TH1F *hist, Double_t gLowX_1, Double_t gUpX_1, Double_t gLowX_2, Double_t gUpX_2)
{
 if(gLowX_1 < gUpX_1)
  {
 if(gLowX_2 < gUpX_2)
  {
// *** Creating the function of the form '2_gaus_lbg' defined above ***
   TF1 fitfunc("db_gauss_linbg_diff",db_gaus_lbg_diff, 0, 1, 8);
// *** Obtaining and specifying the start values for the fit ***
   gBinW      = hist->GetBinWidth(1);
   gContent_1 = gBinW*(hist->Integral(hist->FindBin(gLowX_1),hist->FindBin(gUpX_1)));
   gContent_2 = gBinW*(hist->Integral(hist->FindBin(gLowX_2),hist->FindBin(gUpX_2)));
// *** Searching for maximum Y value through the bins specified by limits
   int i, i_1, i_2;
   Double_t V, V_max;
   
   i_1 = int(hist->FindBin(gLowX_1));
   i_2 = int(hist->FindBin(gUpX_1));
   V_max = hist->GetBinContent(i_1);
   for(i = i_1; i <= i_2; i++)
    {
     V = hist->GetBinContent(i);
     if(V > V_max)
      {
       V_max = V;
       gMean_1 = double(i);
      }
    }
   gMean_1 = gBinW*gMean_1; 

   i_1 = int(hist->FindBin(gLowX_2));
   i_2 = int(hist->FindBin(gUpX_2));
   V_max = hist->GetBinContent(i_1);
   for(i = i_1; i <= i_2; i++)
    {
     V = hist->GetBinContent(i);
     if(V > V_max)
      {
       V_max = V;
       gMean_2 = double(i);
      }
    }
   gMean_2 = gBinW*gMean_2; 

   gSigma_1 = 0.3 * (gUpX_1 - gLowX_1);
   gSigma_2 = 0.3 * (gUpX_2 - gLowX_2);
   
   printf("__________________\n_Peak 1: The Start Values_\n  Bin Width: %f\n Mean Value: %f\n    Content: %f\n      Sigma: %f\n__________________\n",gBinW,gMean_1,gContent_1,gSigma);
   printf("__________________\n_Peak 2: The Start Values_\n  Bin Width: %f\n Mean Value: %f\n    Content: %f\n      Sigma: %f\n__________________\n",gBinW,gMean_2,gContent_2,gSigma);
   fitfunc.SetParameters(0, 0, gSigma_1, gContent_1, gMean_1, gSigma_2, gContent_2, gMean_2);
   fitfunc.SetRange(gLowX_1, gUpX_2);
  
   fitfunc.SetParName(0,"BgConstant");
   fitfunc.SetParName(1,"BgSlope   ");
   fitfunc.SetParName(2,"Sigma 1   ");
   fitfunc.SetParName(3,"Content 1 ");
   fitfunc.SetParName(4,"Mean 1    ");
   fitfunc.SetParName(5,"Sigma 2   ");
   fitfunc.SetParName(6,"Content 2 ");
   fitfunc.SetParName(7,"Mean 2    ");

// *** Fitting: 'R' means within the range specified above ***
   hist->Fit("db_gauss_linbg_diff", "R", "SAME");

   gBgConstant = fitfunc.GetParameter(0);
   gBgSlope    = fitfunc.GetParameter(1);
   gSigma_1    = fitfunc.GetParameter(2);
   gContent_1  = fitfunc.GetParameter(3);
   gMean_1     = fitfunc.GetParameter(4);
   gSigma_2    = fitfunc.GetParameter(5);
   gContent_2  = fitfunc.GetParameter(6);
   gMean_2     = fitfunc.GetParameter(7);
   gChi2pNDF   = fitfunc.GetChisquare() / fitfunc.GetNDF();

   printf("      Chi Square: %f\n",fitfunc.GetChisquare());
   printf("            FWHM 1: %f +- %f\n",2*gSigma_1*sqrt(2*log(2)),2*sqrt(2*log(2))*fitfunc.GetParError(2));
   printf("            FWHM 2: %f +- %f\n",2*gSigma_2*sqrt(2*log(2)),2*sqrt(2*log(2))*fitfunc.GetParError(5));
  } else cout << "Couldn't fit! Error: Peak 2: The Lower Limit is larger than the Upper Limit!" << endl;
  } else cout << "Couldn't fit! Error: Peak 1: The Lower Limit is larger than the Upper Limit!" << endl;
}



void DeconvoluteDoublePeak(TH1D *hist, Double_t gLowX_1, Double_t gUpX_1, Double_t gLowX_2, Double_t gUpX_2)
{
 if(gLowX_1 < gUpX_1)
  {
 if(gLowX_2 < gUpX_2)
  {
// *** Creating the function of the form '2_gaus_lbg' defined above ***
   TF1 fitfunc("db_gauss_linbg",db_gaus_lbg, 0, 1, 7);
// *** Obtaining and specifying the start values for the fit ***
   gBinW      = hist->GetBinWidth(1);
   gContent_1 = gBinW*(hist->Integral(hist->FindBin(gLowX_1),hist->FindBin(gUpX_1)));
   gContent_2 = gBinW*(hist->Integral(hist->FindBin(gLowX_2),hist->FindBin(gUpX_2)));
// *** Searching for maximum Y value through the bins specified by limits
   int i, i_1, i_2;
   Double_t V, V_max;
   
   i_1 = int(hist->FindBin(gLowX_1));
   i_2 = int(hist->FindBin(gUpX_1));
   V_max = hist->GetBinContent(i_1);
   for(i = i_1; i <= i_2; i++)
    {
     V = hist->GetBinContent(i);
     if(V > V_max)
      {
       V_max = V;
       gMean_1 = double(i);
      }
    }
   gMean_1 = gBinW*gMean_1; 

   i_1 = int(hist->FindBin(gLowX_2));
   i_2 = int(hist->FindBin(gUpX_2));
   V_max = hist->GetBinContent(i_1);
   for(i = i_1; i <= i_2; i++)
    {
     V = hist->GetBinContent(i);
     if(V > V_max)
      {
       V_max = V;
       gMean_2 = double(i);
      }
    }
   gMean_2 = gBinW*gMean_2; 

   gSigma     = 0.5 * (0.3 * (gUpX_1 - gLowX_1) + 0.3 * (gUpX_2 - gLowX_2));
   printf("__________________\n_Peak 1: The Start Values_\n  Bin Width: %f\n Mean Value: %f\n    Content: %f\n      Sigma: %f\n__________________\n",gBinW,gMean_1,gContent_1,gSigma);
   printf("__________________\n_Peak 2: The Start Values_\n  Bin Width: %f\n Mean Value: %f\n    Content: %f\n      Sigma: %f\n__________________\n",gBinW,gMean_2,gContent_2,gSigma);
   fitfunc.SetParameters(0, 0, gSigma, gContent_1, gMean_1, gContent_2, gMean_2);
   fitfunc.SetRange(gLowX_1, gUpX_2);
  
   fitfunc.SetParName(0,"BgConstant");
   fitfunc.SetParName(1,"BgSlope   ");
   fitfunc.SetParName(2,"Sigma     ");
   fitfunc.SetParName(3,"Content 1 ");
   fitfunc.SetParName(4,"Mean 1    ");
   fitfunc.SetParName(5,"Content 2 ");
   fitfunc.SetParName(6,"Mean 2    ");

// *** Fitting: 'R' means within the range specified above ***
   hist->Fit("db_gauss_linbg", "R", "SAME");

   gBgConstant = fitfunc.GetParameter(0);
   gBgSlope    = fitfunc.GetParameter(1);
   gSigma      = fitfunc.GetParameter(2);
   gContent_1  = fitfunc.GetParameter(3);
   gMean_1     = fitfunc.GetParameter(4);
   gContent_2  = fitfunc.GetParameter(5);
   gMean_2     = fitfunc.GetParameter(6);
   gChi2pNDF   = fitfunc.GetChisquare() / fitfunc.GetNDF();

   printf("      Chi Square: %f\n",fitfunc.GetChisquare());
   printf("            FWHM: %f +- %f\n",2*gSigma*sqrt(2*log(2)),2*sqrt(2*log(2))*fitfunc.GetParError(2));
  } else cout << "Couldn't fit! Error: Peak 2: The Lower Limit is larger than the Upper Limit!" << endl;
  } else cout << "Couldn't fit! Error: Peak 1: The Lower Limit is larger than the Upper Limit!" << endl;
}
//_____________________________________________________

//
// from here down functions written by Irina
//

Double_t powgen(Double_t *x, Double_t *par)
{
/*
  par[0]   area G
  par[1]   FWHM
  par[2]   mean G
  par[3]   area L
  par[4]   mean L
  par[5]   eta

 */
   
   static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);

  

  // 1. Calculate Gauss part  
    
   Double_t argG = (x[0] - par[2])/(sqrt2*par[0]);
   
   Double_t fitvalG = gBinW/(sqrt2pi*par[0]) * par[1] * exp(-argG*argG);
   
  // 2. Calculate Lorentz part

   Double_t fitvalL = gBinW*(0.5*par[3]*par[1]/TMath::Pi()) /
                     TMath::Max( 1.e-10,(x[0]-par[4])*(x[0]-par[4])
                     + .25*par[1]*par[1]);

   Double_t fitval = par[5] * fitvalG + (1 - par[5]) * fitvalL;  
                      
   return fitval;
}
//_____________________________________________________

Double_t exp_gaus(Double_t *x, Double_t *par)
{
/*
  par[0]   area G
  par[1]   sigma G
  par[2]   mean G
  par[3]   exponential constant
  par[4]   Bkg constant
  par[5]   Bkg slope


 */
   
   static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
   
   Double_t y = x[0] - par[2];
   Double_t sigma2 = par[1]*par[1];
   Double_t lambda = -par[3];

    
   double u = 0.5 * lambda * (y - lambda * sigma2 * 0.5);  // exponential
   double yy = (y - lambda * sigma2)/(sqrt2 * par[1]);    // arg err function

   Double_t fitval = par[4] + x[0]*par[5]  
                      + gBinW * par[0] * TMath::Exp(u) * ( 1 + ROOT::Math::erf(yy));  
                      
   return fitval;
}
//_____________________________________________________

void GetFitHelp(void)
{
/* printf("=============================================================================================\n");
 printf("   The following functions are available:\n");
 printf("---------------------------------------------------------------------------------------------\n");
 printf("          FindSinglePeak(TH1F *HistogramName, Double_t LeftLimit,        Double_t Energy)\n");
 printf("           FitSinglePeak(TH1F *HistogramName, Double_t LeftLimit,        Double_t RightLimit)\n");
 printf("           FitDoublePeak(TH1F *HistogramName, Double_t LeftLimit_Peak_1, Double_t RightLimit_Peak_1,\n");
 printf("                                              Double_t LeftLimit_Peak_2, Double_t RightLimit_Peak_2)\n");
 printf("   DeconvoluteDoublePeak(TH1F *HistogramName, Double_t LeftLimit_Peak_1, Double_t RightLimit_Peak_1,\n");
 printf("                                              Double_t LeftLimit_Peak_2, Double_t RightLimit_Peak_2)\n");
 printf("---------------------------------------------------------------------------------------------\n");
 printf("---------------------------------------------------------------------------------------------\n");
 printf(" * Information\n");
 printf("---------------------------------------------------------------------------------------------\n");
 printf(" *     TH1F : histograms with one float per channel. Maximum precision 7 digits.\n");
 printf(" * Double_t : the same as 'double' - Float 8 bytes\n");
 printf(" *  Float_t : the same as 'float'  - Float 4 bytes\n");
 printf("---------------------------------------------------------------------------------------------\n");
 printf("   Have a Nice Work!\n");
 printf("=============================================================================================\n");*/
}

void FindSinglePeak(TH1D *hist, Double_t gEnergy)
{
 Double_t gBinW = hist->GetBinWidth(1);
 Double_t gLowX = gEnergy - gBinW*10.;
 Double_t gUpX  = gEnergy + gBinW*10.;
 hist->GetXaxis()->SetRange(gLowX,gUpX);
 FitSinglePeak(hist, gLowX, gUpX);
}

//_____________________________________________________


void FindSinglePeak(TH1F *hist, Double_t gEnergy)
{
 Double_t gBinW = hist->GetBinWidth(1);
 Double_t gLowX = gEnergy - gBinW*10.;
 Double_t gUpX  = gEnergy + gBinW*10.;
 hist->GetXaxis()->SetRange(gLowX,gUpX);
 FitSinglePeak(hist, gLowX, gUpX);
}

//_____________________________________________________


void FitPowgen(TH1F *hist, Double_t gLowX, Double_t gUpX)
{

Double_t gScaling, gEta, gFWHM1, gFWHM2, gSigma, gContent1, gContent2;

 if(gLowX < gUpX)
  { 
// *** Creating the function of the form 'gaus_lbg' defined above ***
   TF1 fitfunc("powgen",powgen, 0, 1, 6);
   
// *** Obtaining and specifying the start values for the fit ***
   gContent = hist->Integral(hist->FindBin(gLowX),hist->FindBin(gUpX)); 
   gMean    = gLowX + 0.5 * ( gLowX + gUpX);  
   gSigma   = 0.3 * ( gUpX  - gLowX); 
   gBinW = hist->GetBinWidth(1);
   
// *** Searching for maximum Y value through the bins specified by limits
   int i, i_1, i_2;
   Double_t V, V_max;
   
   i_1 = int(hist->FindBin(gLowX));
   i_2 = int(hist->FindBin(gUpX));
      
   V_max = hist->GetBinContent(i_1);
   for(i = i_1; i <= i_2; i++)
    {
     V = hist->GetBinContent(i);
     if(V > V_max)
      {
       V_max = V;
       gMean = double(i);
      }
    }
    
    printf("max peak %f\n",V_max);
    printf("mean peak %f\n",gMean*gBinW);

//   printf("__________________\n_The Start Values_\n  Bin Width: %d\n Mean Value: %d\n    Content: %d\n      Sigma: %d\n__________________\n",gBinW,gMean,gContent,gSigma);

//   fitfunc.SetParameters(1., 0.1, gSigma, gMean, gContent/2, gContent/2); 
   fitfunc.SetParameters(gContent/2, gSigma, gMean, gContent/2, gMean, gEta); 
   fitfunc.SetRange(gLowX, gUpX);
  
   fitfunc.SetParName(0,"Area G              ");
   fitfunc.SetParName(1,"FWHM                ");
   fitfunc.SetParName(2,"Mean G                ");
   fitfunc.SetParName(3,"Area L              ");
   fitfunc.SetParName(4,"Mean L               ");
   fitfunc.SetParName(5,"Eta                 ");


// *** Fitting: 'R' means within the range specified above ***
   hist->Fit("powgen", "R", "SAME");

   gContent1        = fitfunc.GetParameter(0);
   gFWHM1           = fitfunc.GetParameter(1);
   gMean           = fitfunc.GetParameter(2);
   gContent2        = fitfunc.GetParameter(3);   
   gFWHM2           = fitfunc.GetParameter(4);
//   gEta            = fitfunc.GetParameter(5);
   
   gStyle->SetOptFit(1);

   printf("      Chi Square: %f\n",fitfunc.GetChisquare());
  } else cout << "Couldn't fit! Error: The Lower Limit is larger than the Upper Limit!" << endl;
}
//_____________________________________________________

void FitHeimdal(TH1F *hist, Double_t gLowX, Double_t gUpX)
{

Double_t gExp, gMean, gSigma, gContent;

 if(gLowX < gUpX)
  { 
// *** Creating the function of the form 'gaus_lbg' defined above ***
   TF1 fitfunc("exp_gaus",exp_gaus, 0, 1, 6);
   
// *** Obtaining and specifying the start values for the fit ***
   gContent = hist->Integral(hist->FindBin(gLowX),hist->FindBin(gUpX)); 
   gMean    = gLowX + 0.5 * ( gLowX + gUpX);  
   gSigma   = 0.3 * ( gUpX  - gLowX); 
   gBinW = hist->GetBinWidth(1);
   
// *** Searching for maximum Y value through the bins specified by limits
   int i, i_1, i_2;
   Double_t V, V_max;
   
   i_1 = int(hist->FindBin(gLowX));
   i_2 = int(hist->FindBin(gUpX));
      
   V_max = hist->GetBinContent(i_1);
   for(i = i_1; i <= i_2; i++)
    {
     V = hist->GetBinContent(i);
     if(V > V_max)
      {
       V_max = V;
       gMean = double(i);
      }
    }
    
    printf("max peak %f\n",V_max);
    printf("mean peak %f\n",gMean*gBinW);
    printf("bins per channel  %f\n",gBinW);

//   printf("__________________\n_The Start Values_\n  Bin Width: %d\n Mean Value: %d\n    Content: %d\n      Sigma: %d\n__________________\n",gBinW,gMean,gContent,gSigma);

//   fitfunc.SetParameters(1., 0.1, gSigma, gMean, gContent/2, gContent/2); 
   fitfunc.SetParameters(gContent, gSigma, gMean*gBinW, 1, 0 , 0); 
   fitfunc.SetRange(gLowX, gUpX);
  
   fitfunc.SetParName(4,"BgConstant                ");
   fitfunc.SetParName(5,"BgSlope                   ");
   fitfunc.SetParName(0,"Area G                    ");
   fitfunc.SetParName(1,"Sigma                     ");
   fitfunc.SetParName(2,"Mean                      ");
   fitfunc.SetParName(3,"Exp Constant              ");


// *** Fitting: 'R' means within the range specified above ***
   hist->Fit("exp_gaus", "R", "SAME");

   gContent         = fitfunc.GetParameter(0);
   gSigma           = fitfunc.GetParameter(1);
   gMean            = fitfunc.GetParameter(2);
   gExp             = fitfunc.GetParameter(3);   
   gBgConstant      = fitfunc.GetParameter(4);
   gBgSlope         = fitfunc.GetParameter(5);
   gChi2pNDF        = fitfunc.GetChisquare() / fitfunc.GetNDF();
   
   gStyle->SetOptFit(1);

   printf("      Chi Square: %f\n",fitfunc.GetChisquare());
   printf("      Reduced Chi Square: %f\n",gChi2pNDF);
  } else cout << "Couldn't fit! Error: The Lower Limit is larger than the Upper Limit!" << endl;
}
//_____________________________________________________

void FitLorentzPeak(TH1F *hist, Double_t gLowX, Double_t gUpX)
{

Double_t gArea, gFWHM;

 if(gLowX < gUpX)
  { 
// *** Creating the function of the form 'gaus_lbg' defined above ***
   TF1 fitfunc("lorentz",lorentz, 0, 1, 5);
// *** Obtaining and specifying the start values for the fit ***
   gContent = hist->Integral(hist->FindBin(gLowX),hist->FindBin(gUpX)); 
   gBinW = hist->GetBinWidth(1);

// *** Searching for maximum Y value through the bins specified by limits
   int i, i_1, i_2;
   Double_t V, V_max, gHeight;
   
   i_1 = int(hist->FindBin(gLowX));
   i_2 = int(hist->FindBin(gUpX));
   
   printf("Bins for the peak:  %d   %d\n",i_1,i_2);
   
   V_max = hist->GetBinContent(i_1);
   for(i = i_1; i <= i_2; i++)
    {
     V = hist->GetBinContent(i);
     if(V > V_max)
      {
       V_max = V;
       gMean = double(i);
      }
    }
    
   gMean = gBinW*gMean; 
   gFWHM   = 0.3 * ( gUpX  - gLowX); 
   gArea = gContent;

   printf("gArea = %f\n",gArea);
   printf("gMean = %f\n",gMean);
   printf("gFWHM = %f\n",gFWHM);

   fitfunc.SetParameters(0, 0, gArea, gFWHM, gMean); 
   fitfunc.SetRange(gLowX, gUpX);
  
   fitfunc.SetParName(0,"BgConstant");
   fitfunc.SetParName(1,"BgSlope   ");
   fitfunc.SetParName(2,"Area      ");
   fitfunc.SetParName(3,"Sigma      ");
   fitfunc.SetParName(4,"Mean      ");

// *** Fitting: 'R' means within the range specified above ***
   hist->Fit("lorentz", "R", "SAME");

   gBgConstant  = fitfunc.GetParameter(0);
   gBgSlope     = fitfunc.GetParameter(1);
   gArea        = fitfunc.GetParameter(2);
   gFWHM        = fitfunc.GetParameter(3)/2.35;
   gMean        = fitfunc.GetParameter(4);
   gChi2pNDF    = fitfunc.GetChisquare() / fitfunc.GetNDF();

   gStyle->SetOptFit(1);

   printf("      Chi Square: %f\n",fitfunc.GetChisquare());
   printf("      Reduced Chi Square: %f\n",gChi2pNDF);
  } else cout << "Couldn't fit! Error: The Lower Limit is larger than the Upper Limit!" << endl;
}

//_____________________________________________________


// next function created by Irina to fit diffraction peaks with a PseudoVoigt function.
// Description: (taken from the Mantid website)
//
//The Pseudo-Voigt function is an approximation for the Voigt function, 
//which is a convolution of Gaussian and Lorentzian function. 
//It is often used as a peak profile in powder diffraction for cases where neither a 
// pure Gaussian or Lorentzian function appropriately describe a peak.
//
// Instead of convoluting those two functions, the Pseudo-Voigt function is defined as 
// the sum of a Gaussian peak G(x) and a Lorentzian peak L(x), weighted by a fourth 
// parameter \eta (values between 0 and 1) which shifts the profile more towards pure 
// Gaussian or pure Lorentzian when approaching 1 or 0 respectively:
// PV(x) = \eta G(x) + (1 - \eta)L(x)
//
// Both functions share three parameters: Height (height of the peak at the maximum), 
// PeakCentre (position of the maximum) and FWHM (full width at half maximum of the peak).



void FitPseudoVoigt(TH1F *hist, Double_t gLowX_1, Double_t gUpX_1, Double_t gLowX_2, Double_t gUpX_2)
{

//Double_t gY1,gY2,gLowY,gUpY,gSigma,gMean,gContent,lContent,lSigma,lMean,Mixing,gBinW;
Double_t gSigma,gMean,gHeight,lHeight,lSigma,lMean,Mixing,gBinW,gContent,lContent;

 if(gLowX_1 < gUpX_1)
  {
 if(gLowX_2 < gUpX_2)
  {
// *** Creating the function of the form 'gaus_lbg' defined above ***
   TF1 fitfunc("gauss_lorentz_lbg",gaus_lorentz_lbg, 0, 1, 8);
   
// *** Obtaining and specifying the start values for the fit ***

   gBinW = hist->GetBinWidth(1);

//   gY1 = hist->GetBinContent(hist->FindBin(gX1));   // the maximum in the gaus peak
//   gY2 = hist->GetBinContent(hist->FindBin(gX2));   // the maximum in the lorentz peak
//   gLowY = hist->GetBinContent(hist->FindBin(gLowX));   // the basis, left, bkg left 
//   gUpY = hist->GetBinContent(hist->FindBin(gUpX));   // the basis, right, bkg right 

   printf(" \n");
   printf("Guess parameters:\n");

   gContent = (hist->Integral(hist->FindBin(gLowX_1),hist->FindBin(gUpX_1)));
   lContent = (hist->Integral(hist->FindBin(gLowX_2),hist->FindBin(gUpX_2)));

   printf("gContent = %f\n",gContent);
   printf("lContent = %f\n",lContent);

   // *** Searching for maximum Y value through the bins specified by limits
   int i, i_1, i_2;
   Double_t V, V_max;
   
   i_1 = int(hist->FindBin(gLowX_1));
   i_2 = int(hist->FindBin(gUpX_1));
   
   printf("Bins for the Gauss part of the peak:  %d   %d\n",i_1,i_2);
   
   V_max = hist->GetBinContent(i_1);
   for(i = i_1; i <= i_2; i++)
    {
     V = hist->GetBinContent(i);
     if(V > V_max)
      {
       V_max = V;
       gMean = double(i);
      }
    }
   gMean = gBinW*gMean; 
   gHeight = V_max;

   printf("gHeight = %f\n",gHeight);

   i_1 = int(hist->FindBin(gLowX_2));
   i_2 = int(hist->FindBin(gUpX_2));
   printf("Bins for the Lorentz part of the peak:  %d   %d\n",i_1,i_2);

   V_max = hist->GetBinContent(i_1);
   for(i = i_1; i <= i_2; i++)
    {
     V = hist->GetBinContent(i);
     if(V > V_max)
      {
       V_max = V;
       lMean = double(i);
      } else {lMean = double(i_1);}
    }
   lMean = gBinW*lMean; 
   lHeight = V_max;
   printf("lHeight = %f\n",lHeight);

   gSigma = 0.3 * (gUpX_1 - gLowX_1)/gBinW;
   lSigma = 0.3 * (gUpX_2 - gLowX_2)/gBinW;

//   Mixing   = 1;
   
   printf("gMean = %f\n",gMean);
   printf("lMean = %f\n",lMean);

   printf("gSigma = %f\n",gSigma);
   printf("lSigma = %f\n",lSigma);
   
//   printf("Mixing = %e\n",Mixing);
   printf("  \n");
    
   fitfunc.SetParameters(0.0,0.0, gSigma, gHeight, gMean, lMean, lSigma, lHeight); 
   fitfunc.SetRange(gLowX_1, gUpX_2);
  
   fitfunc.SetParName(0,"Bkg            ");
   fitfunc.SetParName(1,"BkgSlope       ");
   fitfunc.SetParName(2,"GaussWidth     ");
   fitfunc.SetParName(3,"GaussContent   ");
   fitfunc.SetParName(4,"GaussMean      ");
   fitfunc.SetParName(5,"LorentzWidth   ");
   fitfunc.SetParName(6,"LorentzMean    ");
   fitfunc.SetParName(7,"LorentzContent ");
//   fitfunc.SetParName(8,"Mixing         ");

// *** Fitting: 'R' means within the range specified above ***
   hist->Fit("gauss_lorentz_lbg", "R", "SAME");

   gBgConstant = fitfunc.GetParameter(0);
   gBgSlope    = fitfunc.GetParameter(1);
   gSigma      = fitfunc.GetParameter(2);
   gContent    = fitfunc.GetParameter(3);
   gMean       = fitfunc.GetParameter(4);
   lMean       = fitfunc.GetParameter(5);
   lSigma      = fitfunc.GetParameter(6);
   lContent    = fitfunc.GetParameter(7);
//   Mixing      = fitfunc.GetParameter(8);
   
   gChi2pNDF   = fitfunc.GetChisquare() / fitfunc.GetNDF();
   
   gStyle->SetOptFit(1);

   printf("      Chi Square: %f\n",fitfunc.GetChisquare());
   printf("            FWHM: %f +- %f\n",2*gSigma*sqrt(2*log(2)),2*sqrt(2*log(2))*fitfunc.GetParError(2));
  } else cout << "Couldn't fit! Error: The Lower Limit2 is larger than the Upper Limit2!" << endl;
  } else cout << "Couldn't fit! Error: The Lower Limit1 is larger than the Upper Limit1!" << endl;
}

//_____________________________________________________


// next function created by Irina to fit diffraction peaks with two back2back exponentials

void FitBk2BkExp(TH1F *hist, Double_t gLowX, Double_t gUpX)
{

Double_t gY,gContent,gSigma,gMean,gA,gB;

 if(gLowX < gUpX)
  { 
// *** Creating the function of the form 'gaus_lbg' defined above ***
   TF1 fitfunc("db_exp_lbg",db_exp_lbg, 0, 1, 7);
// *** Obtaining and specifying the start values for the fit ***
   gBinW = hist->GetBinWidth(1);


   gContent = (hist->Integral(hist->FindBin(gLowX),hist->FindBin(gUpX)));

   gSigma   = 0.3 * ( gUpX  - gLowX); 
   gMean   = gLowX + 0.5 * ( gUpX  - gLowX); 
   gA = 1.0;
   gB = 0.05;
   
   printf("gBinW = %f\n",gBinW);
   printf("gContent = %f\n",gContent);
   printf("gSigma = %f\n",gSigma);
   printf("gMean = %f\n",gMean);

   fitfunc.SetParameters(0.0, 0.0, gContent, gA, gB, gMean, gSigma); 
   fitfunc.SetRange(gLowX, gUpX);
  
   fitfunc.SetParName(0,"BgConstant        ");
   fitfunc.SetParName(1,"BgSlope           ");
   fitfunc.SetParName(2,"Area              ");
   fitfunc.SetParName(3,"Exp1 constant     ");
   fitfunc.SetParName(4,"Exp2 constant     ");
   fitfunc.SetParName(5,"Mean   ");
   fitfunc.SetParName(6,"Sigma   ");

// *** Fitting: 'R' means within the range specified above ***
   hist->Fit("db_exp_lbg", "R", "SAME");

   gBgConstant   = fitfunc.GetParameter(0);
   gBgSlope      = fitfunc.GetParameter(1);
   gContent      = fitfunc.GetParameter(2);
   gA            = fitfunc.GetParameter(3);
   gB            = fitfunc.GetParameter(4);
   gMean         = fitfunc.GetParameter(5);
   gSigma        = fitfunc.GetParameter(6);
   gChi2pNDF     = fitfunc.GetChisquare() / fitfunc.GetNDF();
   
   gStyle->SetOptFit(1);

   printf("      Chi Square: %f\n",fitfunc.GetChisquare());
  } else cout << "Couldn't fit! Error: The Lower Limit is larger than the Upper Limit!" << endl;
}
//_____________________________________________________



