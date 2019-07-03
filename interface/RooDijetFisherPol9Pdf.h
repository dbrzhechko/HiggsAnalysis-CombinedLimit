
#ifndef HiggsAnalysis_CombinedLimit_RooDijetFisherPol9Pdf_h
#define HiggsAnalysis_CombinedLimit_RooDijetFisherPol9Pdf_h
//---------------------------------------------------------------------------
#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "RooRealProxy.h"
//---------------------------------------------------------------------------
class RooRealVar;
class RooAbsReal;

#include "Riostream.h"
#include "TMath.h"
#include <TH1.h>
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

//---------------------------------------------------------------------------
class RooDijetFisherPol9Pdf : public RooAbsPdf
{
public:
   RooDijetFisherPol9Pdf() {} ;
   RooDijetFisherPol9Pdf(const char *name, const char *title,
        RooAbsReal& _th1x, 
        RooAbsReal& _p1,
        RooAbsReal& _p2, 
        RooAbsReal& _p3, 
        RooAbsReal& _p4, 
        RooAbsReal& _p5,
        RooAbsReal& _p6,
        RooAbsReal& _p7,
        RooAbsReal& _p8,
//        RooAbsReal& _p9,
        RooAbsReal& _sqrts);
   RooDijetFisherPol9Pdf(const RooDijetFisherPol9Pdf& other,
      const char* name = 0);
   void setTH1Binning(TH1* _Hnominal);
   void setAbsTol(double _absTol);
   void setRelTol(double _relTol);
   virtual TObject* clone(const char* newname) const { return new RooDijetFisherPol9Pdf(*this,newname); }
   inline virtual ~RooDijetFisherPol9Pdf() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:   

   RooRealProxy th1x;        // dependent variable
   RooRealProxy p1;       // p1
   RooRealProxy p2;        // p2
   RooRealProxy p3;        // p3
   RooRealProxy p4;        // p4
   RooRealProxy p5;        // p5
   RooRealProxy p6;        // p6
   RooRealProxy p7;        // p7
   RooRealProxy p8;        // p8
//   RooRealProxy p9;        // p9
   RooRealProxy sqrts;        // sqrts
   Int_t xBins;        // X bins
   Double_t xArray[2000]; // xArray[xBins+1]
   Double_t xMax;        // X max
   Double_t xMin;        // X min
   Double_t relTol;      //relative tolerance for numerical integration
   Double_t absTol;      //absolute tolerance for numerical integration

   Double_t evaluate() const;
private:
   ClassDef(RooDijetFisherPol9Pdf,1) // RazorDijetFisherPol9Pdf function
    
};
//---------------------------------------------------------------------------
#endif

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
 
class DijetFisherPol9Function: public ROOT::Math::IParametricFunctionOneDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(double x,const double* p) const
   {
     double pdf = TMath::Power(1 + p[2] * x/p[0] + p[3] * TMath::Power(x/p[0],2) + p[4] * TMath::Power(x/p[0],3) + p[5] * TMath::Power(x/p[0],4) + p[6] * TMath::Power(x/p[0],5) + p[7] * TMath::Power(x/p[0],6) + p[8] * TMath::Power(x/p[0],7),-p[1]);
     return pdf;
   }
   double DoEval(double x) const
   {
     double pdf = TMath::Power(1 + pars[2] * x/pars[0] + pars[3] * TMath::Power(x/pars[0],2) + pars[4] * TMath::Power(x/pars[0],3) + pars[5] * TMath::Power(x/pars[0],4) + pars[6] * TMath::Power(x/pars[0],5) + pars[7] * TMath::Power(x/pars[0],6) + pars[8] * TMath::Power(x/pars[0],7),-pars[1]);
     return pdf;
   }
 
   ROOT::Math::IBaseFunctionOneDim* Clone() const
   {
      return new DijetFisherPol9Function();
   }
 
   const double* Parameters() const
   {
      return pars;
   }
 
   void SetParameters(const double* p)
   {
      pars = p;
   }
 
   unsigned int NPar() const
   {
      return 9;
   }
};
