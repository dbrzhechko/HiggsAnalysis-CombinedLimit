
#include "RooFit.h"

#include "Riostream.h"
#include <TMath.h>
#include <cassert>
#include <cmath>
#include <math.h>

#include "../interface/RooDijetFisherAlt7Pdf.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"
#include "Math/GSLIntegrator.h"

using namespace std;
using namespace RooFit;

ClassImp(RooDijetFisherAlt7Pdf)
//---------------------------------------------------------------------------
RooDijetFisherAlt7Pdf::RooDijetFisherAlt7Pdf(const char *name, const char *title,
				   RooAbsReal& _th1x,  
				   RooAbsReal& _p1, 
                   RooAbsReal& _p2, 
				   RooAbsReal& _p3, 
                   RooAbsReal& _p4, 
                   RooAbsReal& _p5, 
                   RooAbsReal& _p6, 
//                   RooAbsReal& _p7, 
//                   RooAbsReal& _p8, 
//                   RooAbsReal& _p9, 
                   RooAbsReal& _sqrts) : RooAbsPdf(name, title), 
//TH3* _Hnominal) : RooAbsPdf(name, title), 
  th1x("th1x", "th1x Observable", this, _th1x),
  p1("p1", "p1", this, _p1),
  p2("p2", "p2", this, _p2),
  p3("p3", "p3", this, _p3),
  p4("p4", "p4", this, _p4),
  p5("p5", "p5", this, _p5),
  p6("p6", "p6", this, _p6),
//  p7("p7", "p7", this, _p7),
//  p8("p8", "p8", this, _p8),
//  p9("p9", "p9", this, _p9),
  sqrts("sqrts", "sqrts", this, _sqrts),
  xBins(0),
  xMax(0),
  xMin(0),
  relTol(1E-12),
  absTol(1E-12)
{
  memset(&xArray, 0, sizeof(xArray));
}
//---------------------------------------------------------------------------
RooDijetFisherAlt7Pdf::RooDijetFisherAlt7Pdf(const RooDijetFisherAlt7Pdf& other, const char* name) :
   RooAbsPdf(other, name), 
   th1x("th1x", this, other.th1x),  
   p1("p1", this, other.p1),
   p2("p2", this, other.p2),
   p3("p3", this, other.p3),
   p4("p4", this, other.p4),
   p5("p5", this, other.p5),
   p6("p6", this, other.p6),
//   p7("p7", this, other.p7),
//   p8("p8", this, other.p8),
//   p9("p9", this, other.p9),
   sqrts("sqrts", this, other.sqrts),
   xBins(other.xBins),
   xMax(other.xMax),
   xMin(other.xMin),
   relTol(other.relTol),
   absTol(other.absTol)
{
  //memset(&xArray, 0, sizeof(xArray));
  for (Int_t i=0; i<xBins+1; i++){
    xArray[i] = other.xArray[i];
  }
}
//---------------------------------------------------------------------------
void RooDijetFisherAlt7Pdf::setTH1Binning(TH1* _Hnominal){
  xBins = _Hnominal->GetXaxis()->GetNbins();
  xMin = _Hnominal->GetXaxis()->GetBinLowEdge(1);
  xMax = _Hnominal->GetXaxis()->GetBinUpEdge(xBins);
  memset(&xArray, 0, sizeof(xArray));
  for (Int_t i=0; i<xBins+1; i++){
    xArray[i] =  _Hnominal->GetXaxis()->GetBinLowEdge(i+1);
  }
}
//---------------------------------------------------------------------------
void RooDijetFisherAlt7Pdf::setRelTol(double _relTol){
  relTol = _relTol;
}
//---------------------------------------------------------------------------
void RooDijetFisherAlt7Pdf::setAbsTol(double _absTol){
  absTol = _absTol;
}
//---------------------------------------------------------------------------
Double_t RooDijetFisherAlt7Pdf::evaluate() const
{
  Double_t integral = 0.0;
  

  Int_t iBin = (Int_t) th1x;
  if(iBin < 0 || iBin >= xBins) {
    //cout << "in bin " << iBin << " which is outside of range" << endl;
    return 0.0;
  }

  
  Double_t xLow = xArray[iBin];
  Double_t xHigh = xArray[iBin+1];
    
  // define the function to be integrated numerically
  DijetFisherAlt7Function func;
  double params[7];
  params[0] = sqrts;    
  params[1] = p1;
  params[2] = p2;       
  params[3] = p3;
  params[4] = p4;       
  params[5] = p5;
  params[6] = p6;
//  params[7] = p7;
//  params[8] = p8;
//  params[9] = p9;
  func.SetParameters(params);

  ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,absTol,relTol);
  ig.SetFunction(func,false);

  
  integral = ig.Integral(xLow,xHigh);
  //Double_t total_integral = ig.Integral(xMin,xMax);

  if (integral>0.0) {
    return integral;
  } else return 0;

}

// //---------------------------------------------------------------------------
Int_t RooDijetFisherAlt7Pdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const{
  if (matchArgs(allVars, analVars, th1x)) return 1;
  return 0;
}

// //---------------------------------------------------------------------------
Double_t RooDijetFisherAlt7Pdf::analyticalIntegral(Int_t code, const char* rangeName) const{

   Double_t th1xMin = th1x.min(rangeName); Double_t th1xMax = th1x.max(rangeName);
   Int_t iBinMin = (Int_t) th1xMin; Int_t iBinMax = (Int_t) th1xMax;


   Double_t integral = 0.0;
      
   //cout <<  "iBinMin = " << iBinMin << ",iBinMax = " << iBinMax << endl;

   
   // define the function to be integrated numerically  
   DijetFisherAlt7Function func;
   double params[8];
   params[0] = sqrts;    
   params[1] = p1;
   params[2] = p2;       
   params[3] = p3;
   params[4] = p4;       
   params[5] = p5;
   params[6] = p6;
//   params[7] = p7;
//   params[8] = p8;
//   params[9] = p9;
   func.SetParameters(params);

   
  ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,absTol,relTol);
  ig.SetFunction(func,false);
    

   if (code==1 && iBinMin<=0 && iBinMax>=xBins){
     integral = ig.Integral(xMin,xMax);
     
   }
   else if(code==1) { 
     for (Int_t iBin=iBinMin; iBin<iBinMax; iBin++){
       
       if(iBin < 0 || iBin >= xBins) {
	 integral += 0.0;
       }
       else{	 
	 Double_t xLow = xArray[iBin];
	 Double_t xHigh = xArray[iBin+1];    
	 integral += ig.Integral(xLow,xHigh);
       }
     }
   } else {
     cout << "WARNING IN RooDijetFisherAlt7Pdf: integration code is not correct" << endl;
     cout << "                           what are you integrating on?" << endl;
     return 1.0;
   }

   if (integral>0.0) {
     
     return integral;
   } else return 1.0;
}
