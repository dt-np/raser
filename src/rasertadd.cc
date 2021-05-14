// 
// Author: SHI Xin <shixin@ihep.ac.cn> 
// Created [2021-03-31 Sun 15:05] 
// Based on KDetSim : http://kdetsim.org 
// 


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TH2D.h> 
#include <TFile.h> 
#include <TCanvas.h> 
#include <TStyle.h> 
#include <Riostream.h>
#include <TTree.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TApplication.h> 
#include <TMath.h> 
#include <TLine.h> 
#include <TVector3.h>
#include "TH1.h"
#include "TRandom3.h"
#include "TGaxis.h"
// KDetSim
using namespace std;
// nrutil
// double *dvector(long nl, long nh);
#define NR_END 1 
#define FREE_ARG char *
#define EPS 1.0e-14
#define MAXPOINT 10001



static double e_0 = 1.60217733e-19;          //As
static double perm0 = 8.854187817e-12; // F/m
static double Kboltz = 8.617385e-5; // eV/K

void nrerror(char error_text[])
	/* Numerical Recipes standard error handler */
{
	fprintf(stderr, "Numerical Recipes run-time error...\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr, "...now exiting to system...\n");
	exit(1);
}


double *dvector(long nl, long nh)
	/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v = (double *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double)));
	if (!v)
		nrerror("allocation failure in dvector()");
	return v - nl + NR_END;
}

void free_dvector(double *v, long nl, long nh)
	/* free a double vector allocated with dvector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v = (float *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(float)));
	if (!v)
		nrerror("allocation failure in vector()");
	return v - nl + NR_END;
}
void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}
/***************************************************** 
atimes : How to multiply the vector with the matrices 
Modified in this form by GK 4.10.2012
 ******************************************************/

// #define C1(x,y,z) y3[n]=x+y+z;


double *b,*y6,*y2,*y3,*y4,*y5,*y7,*y8; 


/**********************************************************
snrm: Calculates the norm of vector 
unsigned long n - dimension of the vector
double sx[]     - components
itol            - <=3 real norm , otherwise max component.
 ************************************************************/

double snrm(unsigned long n, double sx[], int itol)
{
	unsigned long i,isamax;
	double ans;

	if (itol <= 3) {
		ans = 0.0;
		for (i=1;i<=n;i++) ans += sx[i]*sx[i];
		return sqrt(ans);
	} else {
		isamax=1;
		for (i=1;i<=n;i++) {
			if (fabs(sx[i]) > fabs(sx[isamax])) isamax=i;
		}
		return fabs(sx[isamax]);
	}
}


/***************************************************
  Solve a very simple system of equations 
  with the diagonal elements to be the only ones 
 ****************************************************/
void asolve(unsigned long n, double b[], double x[], int itrnsp)
{
	//Solve a very simple system of n equations
	//with the diagonal elements to be the only ones!
	unsigned long i;

	for(i=1;i<=n;i++) x[i]=(y3[i] != 0.0 ? b[i]/y3[i] : b[i]);
}

void atimes(unsigned long n,int dim[], double x[],double r[],int itrnsp)
{
	// This is a function used to multuply the vectors with matrices!
	// Used for calculation of the electric and ramo field
	int nx,ny,nz;
	int i,j,k,q=0;
	double C,L,D,O,R,U,I;

	nx=dim[0]; ny=dim[1]; nz=dim[2];

	for(k=1; k<=nz; k++)
		for(j=1; j<=ny; j++)          /*mnozenje po stolpcu*/ 
			for(i=1; i<=nx; i++)           /*mnozenje po vrstici*/ 
			{ 
				q++;
				C=y3[q]*x[q];
				if(q-1>1)       L=y2[q]*x[q-1]; else L=0;
				if(q-nx>1)      D=y6[q]*x[q-nx]; else D=0;
				if((q-nx*ny)>1) I=y7[q]*x[q-ny*nx]; else I=0;

				if(q+1<=n)       R=y4[q]*x[q+1]; else R=0;
				if(q+nx<=n)      U=y5[q]*x[q+nx]; else U=0;
				if((q+nx*ny)<=n) O=y8[q]*x[q+ny*nx]; else O=0;

				r[q]=C+L+D+O+R+U+I;
			}
	if(n!=q) printf("\n Error in matrix solving!");
	return;
}



void linbcg(unsigned long n, int dim[], double b[], double x[], int itol, double tol,
		int itmax, int *iter, double *err)
{
	// The main function for electric field calcualtion
	void asolve(unsigned long n, double b[], double x[], int itrnsp);
	void atimes(unsigned long n, int dim[],double x[], double r[], int itrnsp);
	double snrm(unsigned long n, double sx[], int itol);
	double *dvector(long, long);
	void free_dvector(double *, long, long);
	void nrerror(char error_text[]);
	unsigned long j;
	double ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
	double *p,*pp,*r,*rr,*z,*zz;

	p=dvector(1,n);
	pp=dvector(1,n);
	r=dvector(1,n);
	rr=dvector(1,n);
	z=dvector(1,n);
	zz=dvector(1,n);

	*iter=0;
	atimes(n,dim,x,r,0);
	for (j=1;j<=n;j++) {
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
	atimes(n,dim,r,rr,0); // minimal residual invariant
	znrm=1.0;
	if (itol == 1) bnrm=snrm(n,b,itol);
	else if (itol == 2) {
		asolve(n,b,z,0);
		bnrm=snrm(n,z,itol);
	}
	else if (itol == 3 || itol == 4) {
		asolve(n,b,z,0);
		bnrm=snrm(n,z,itol);
		asolve(n,r,z,0);
		znrm=snrm(n,z,itol);
	} else nrerror("illegal itol in linbcg");
	asolve(n,r,z,0);
	while (*iter <= itmax) {
		++(*iter);
		zm1nrm=znrm;
		asolve(n,rr,zz,1);
		for (bknum=0.0,j=1;j<=n;j++) bknum += z[j]*rr[j];
		if (*iter == 1) {
			for (j=1;j<=n;j++) {
				p[j]=z[j];
				pp[j]=zz[j];
			}
		}
		else {
			bk=bknum/bkden;
			for (j=1;j<=n;j++) {
				p[j]=bk*p[j]+z[j];
				pp[j]=bk*pp[j]+zz[j];
			}
		}
		bkden=bknum;
		atimes(n,dim,p,z,0);
		for (akden=0.0,j=1;j<=n;j++) akden += z[j]*pp[j];
		ak=bknum/akden;
		atimes(n,dim,pp,zz,1);
		for (j=1;j<=n;j++) {
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
			rr[j] -= ak*zz[j];
		}
		asolve(n,r,z,0);
		if (itol == 1 || itol == 2) {
			znrm=1.0;
			*err=snrm(n,r,itol)/bnrm;
		} else if (itol == 3 || itol == 4) {
			znrm=snrm(n,z,itol);
			if (fabs(zm1nrm-znrm) > EPS*znrm) {
				dxnrm=fabs(ak)*snrm(n,p,itol);
				*err=znrm/fabs(zm1nrm-znrm)*dxnrm;
			} else {
				*err=znrm/bnrm;
				continue;
			}
			xnrm=snrm(n,x,itol);
			if (*err <= 0.5*xnrm) *err /= xnrm;
			else {
				*err=znrm/bnrm;
				continue;
			}
		}
		if(*iter==itmax) {printf("\n Number of iterations exceeded the max. value \n"); 
			printf("iter=%4d err=%12.6f\n",*iter,*err);}

		if (*err <= tol) break;
	}

	free_dvector(p,1,n);
	free_dvector(pp,1,n);
	free_dvector(r,1,n);
	free_dvector(rr,1,n);
	free_dvector(z,1,n);
	free_dvector(zz,1,n);
}

Double_t KAlpha(Double_t E, Double_t T, Short_t Charg, Int_t which)
{
	// Function calculates impact ionization coefficientf
	// for a given E [V/um].
	// Short_t Charg;  ---> Charg=1; holes
	//                 ---> Charg=-1; electrons
	// Int_t which;    ---> 0 -> silicon
	//                 ---> 10 -> diamond Trew parametrization
	//                 ---> 11 -> diamond Watanabe parametrization
	//                 ---> 12 -> diamond Hiraiwa parametrization

	Double_t alp, A, B, a = TMath::Sqrt(10), b = TMath::Sqrt(10);
	switch (which)
	{
	case 0: // R.J. R.J. McIntyre, IEEE Trans. Electron Dev. 46 (8) (1999).
		if (Charg > 0)
			alp = 1.3e-3 * TMath::Exp(-13.2 * (2e7 / (E * 1e6) - 1));
		else
			alp = 2.3e-1 * TMath::Exp(-6.78 * (2e7 / (E * 1e6) - 1));
		break;
	case 1: // Van Oversraeten model Solid-State Electronics 1970
		if (Charg > 0)
			if (E < 40)
				alp = 158.2 * TMath::Exp(-203.6 / E);
			else
				alp = 67.1 * TMath::Exp(-169.6 / E);
		else
			alp = 70.3 * TMath::Exp(-123.1 / E);
		break;
	case 2: // Temperature dependence of the impact ionization - Massey
		// D.J Massey et al., Temperature Dependence of Impact Ionization in Submicrometer Silicon Devices,
		// IEEE Transactions on Electron Devices, vol. 53, no. 9, September 2006.
		if (Charg > 0)
			alp = 113 * TMath::Exp(-(171 + 0.109 * T) / E);
		else
			alp = 44.3 * TMath::Exp(-(96.6 + 0.0499 * T) / E);
		break;
	case 3: // Chynoweth original 1958 silicon
		break;
		//

	case 10:
		// Trew parametrization
		A = 1.935e4;
		B = 7.749e2;
		alp = A * TMath::Exp(-B / E);
		break;
	case 11:
		// Watanabe parametrization
		if (Charg > 0)
		{
			A = 19.3;
			B = 4.41e2;
		}
		else
		{
			A = 46.2;
			B = 7.59e2;
		}
		alp = A * TMath::Exp(-B / E);
		break;
	case 12:
		// Hiraiwa parametrization
		if (Charg > 0)
		{
			A = 19.3 / a;
			B = 4.41e2 * b;
		}
		else
		{
			A = 46.2 / a;
			B = 7.59e2 * b;
		}
		alp = A * TMath::Exp(-B / E);
		break;
	}
	return alp;
}
// KMaterial

class KMaterial
{
private:
public:
	static Int_t Mat;			   // Material index
	static Float_t Temperature;	   // Temperature
	static Int_t Mobility;		   // mobility model for each material
	static Int_t ImpactIonization; // impact ionization model

	//////////////////////////////////////////////////////

	KMaterial() { Mat = 1; } // MobMod=1;}
	~KMaterial(){};
	static Float_t Perm(Int_t = 1);
	static Float_t Loss_energy(Int_t = 1);
	static Int_t MobMod();
	// ClassDef(KMaterial, 1)
};

// ClassImp(KMaterial)
Int_t KMaterial::Mat = 1;
Float_t KMaterial::Temperature = 293;
Int_t KMaterial::Mobility = 1;
Int_t KMaterial::ImpactIonization = 0;

Float_t KMaterial::Perm(Int_t Material)
{
	Float_t perm;
	switch (Material)
	{
	case 0:
		perm = 9.76;
		break; //silicon
	case 1:
		perm = 9.76;
		break; //poly silicon
	// case 0:
	// 	perm = 11.7;
	// 	break; //silicon
	// case 1:
	// 	perm = 11.7;
	// 	break; //poly silicon
	case 10:
		perm = 5.7;
		break; //diamond
	case 20:
		perm = 1;
		break; //air
	case 100:
		perm = 1;
		break; //aluminium
	default:
		perm = 1;
		break;
	}
	return perm;
}
Float_t KMaterial::Loss_energy(Int_t Material)
{
	Float_t loss_energy;
	switch (Material)
	{
	case 0:
		loss_energy = 3.6;
		break; //silicon
	case 1:
		loss_energy = 8.4;
		break; // silicon carbide
	return loss_energy;
}
}

Int_t KMaterial::MobMod()
{
	Int_t ret;
	switch (Mat)
	{
	case 0:
		ret = Mobility;
		break; //if(Mobility==1) ret=1; else ret=0; break; //silicon
	case 1:
		ret = 8;
		break; //poly silicon
	case 2:
		ret = 9;
		break; //silicon oxide 2.648
	case 10:
		ret = 10;
		break; //diamond
	}
	return ret;
}

// KStruct
class KStruct 

{
	public:
		Int_t PCharge;
		Int_t Steps;
		Int_t DStrip;
		Float_t Xlenght;
		Float_t Ylenght;
		Float_t Zlenght;
		Float_t TTime;
		Float_t TCharge;
		Float_t Xtrack[MAXPOINT];
		Float_t Ytrack[MAXPOINT];
		Float_t Ztrack[MAXPOINT];
		Float_t Charge[MAXPOINT];
		Float_t Time[MAXPOINT];
		Float_t Efield[MAXPOINT];
		Float_t MulCar[MAXPOINT];


		KStruct();
		~KStruct(){};
		void GetCH(TH1F *, Int_t = 0, Float_t = 1, Float_t = -1);		 //,Int_t=200, Float_t=100e-9); Changed from 2.23 -> 3.0
		Float_t GetCHMult(TH1F *, Int_t = 0, Float_t = 1, Float_t = -1); //,Int_t=200, Float_t=100e-9); Changed from 2.23 -> 3.0
		void Clear();  
};

KStruct::KStruct()
{
	Clear();
}


void KStruct::Clear()
{
	Int_t i = 0;
	for (i = 0; i < MAXPOINT; i++) {
		Xtrack[i] = 0;
		Ytrack[i] = 0;
		Ztrack[i] = 0;
		Charge[i] = 0;
		Time[i] = 0;
		Efield[i] = 0;
		MulCar[i] = 0;
	}
	Ylenght = 0;
	Xlenght = 0;
	TTime = 0;
	Steps = 0;
	DStrip = 0;
	TCharge = 0;
}

void KStruct::GetCH(TH1F *histo, Int_t Update, Float_t Mult, Float_t tau)
{
	//
	//
	//   Float_t tau;  trapping time [s]

	Int_t i;
	TH1F *his = new TH1F("ch", "Charge Histogram", histo->GetNbinsX(), histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax());
	Double_t *ch = new Double_t[Steps + 1];
	Axis_t *ti = new Axis_t[Steps + 1]; // Changed when migrating from 2.23 to 2.25 or higher
	for (i = 1; i < Steps + 1; i++)
	{
		ch[i] = Charge[i]*e_0;
		ti[i] = Time[i];
	}								   // Changed when migrating from 2.23 to 2.25
	his->FillN(Steps, &ti[1], &ch[1]); // Changed when migrating from 2.23 to 2.25
	for (i = 1; i < his->GetNbinsX(); i++)
		his->SetBinContent(i, his->GetBinContent(i) /((histo->GetXaxis()->GetXmax() - histo->GetXaxis()->GetXmin()) / histo->GetNbinsX()));
	// Trappping is included if tau>0   // added 20.1.2010
	if (tau > 0)
	{
		for (i = 1; i < his->GetNbinsX(); i++)
			his->SetBinContent(i, his->GetBinContent(i) * TMath::Exp(-(his->GetBinCenter(i) - Time[0]) / tau));
	}
	////////////////////////////////////////////////////////

	if (Update)
	{
		his->Scale(Mult);
		histo->Add(his);
	}
	else
	{
		if (histo)
			delete histo;

		histo = (TH1F *)histo->Clone("CHGet");
	}

	delete his;
	delete[] ti; // Changed when migrating from 2.23 to 2.25
	delete[] ch; // Changed when migrating from 2.23 to 2.25
}

Float_t KStruct::GetCHMult(TH1F *histo, Int_t Update, Float_t Mult, Float_t tau)
{
	//
	//
	//   Float_t tau;  trapping time [s]
	Double_t dx, dy, dif, dift, sum = 1, summ = 1., mulf, traf;
	Int_t i;
	TH1F *his = new TH1F("ch", "Charge Histogram", histo->GetNbinsX(), histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax());
	Double_t *ch = new Double_t[Steps + 1];
	Double_t *Multi = new Double_t[Steps + 1];
	Axis_t *ti = new Axis_t[Steps + 1]; // Changed when migrating from 2.23 to 2.25 or higher

	for (i = 1; i < Steps + 1; i++)
	{
		if (PCharge < 0)
			dif = KAlpha(0.5 * (Efield[i + 1] + Efield[i]), KMaterial::Temperature, -1, KMaterial::ImpactIonization);
		else
			dif = KAlpha(0.5 * (Efield[i + 1] + Efield[i]), KMaterial::Temperature, 1, KMaterial::ImpactIonization);

		dift = Time[i + 1] - Time[i];
		dx = TMath::Sqrt(TMath::Power((Xtrack[i + 1] - Xtrack[i]), 2) + TMath::Power((Ytrack[i + 1] - Ytrack[i]), 2));

		//---- Calculation of multiplication and trapping factors in given step //
		mulf = (1 + dif * dx);
		traf = (1 - dift / tau);
		MulCar[i] = (mulf - 1) * sum * traf;
		summ *= mulf;
		sum *= mulf; // multiplication 8.3.2011
		if (tau > 0)
			sum *= traf; // trapping 8.3.2011
		//---- END ---- //

		ch[i] = Charge[i] * sum;
		ti[i] = Time[i];
		// Trappping is included if tau>0   // added 20.1.2010
		////////////////////////////////////////////////////////

		// printf("%d :: X=%4.1f , Y=%4.1f :: E=%4.2e ::  Time=%4.1e ; Charge=%4.1e ; dif=%5.2e ; MultT=%5.4e Mult=%5.4f hole=%5.3e\n",i,Xtrack[i],Ytrack[i],Efield[i],ti[i],ch[i],dif,sum,summ,MulCar[i]);
	}

	// Changed ti time when migrating from 2.23 to 2.25
	his->FillN(Steps, &ti[1], &ch[1]); // Changed when migrating from 2.23 to 2.25

	if (Update)
	{
		his->Scale(Mult);
		histo->Add(his);
	}
	else
	{
		if (histo)
			delete histo;

		histo = (TH1F *)histo->Clone("CHMult");
	}

	delete his;
	delete[] ti;	// Changed when migrating from 2.23 to 2.25
	delete[] ch;	// Changed when migrating from 2.23 to 2.25
	delete[] Multi; // 8.3.2011
	return summ;
}

// KGeometry 

#include "TH3I.h"
#include "TH2F.h"

// #include "TMath.h"
#include "math.h"

TH2F *KHisProject(void *hisIn,Int_t axis,Int_t Bin1)
{ 
	//Projects any quantity maped to geometry in different views
	Int_t i;
	TH2F *his2D;
	TH3F *his=(TH3F *) hisIn;

	Int_t Nx=his->GetNbinsX();
	Int_t Ny=his->GetNbinsY();
	Int_t Nz=his->GetNbinsZ();

	Double_t *Xbins=new Double_t [Nx+1];
	Double_t *Ybins=new Double_t [Ny+1];;
	Double_t *Zbins=new Double_t [Nz+1];;

	//  printf("%d %d %d\n", Nx,Ny,Nz);
	for(i=0;i<=Nx;i++) 
	{
		Xbins[i]=his->GetXaxis()->GetBinLowEdge(i+1);
		   //printf("x:: %d %f\n",i,Xbins[i]);
	}
	for(i=0;i<=Ny;i++) 
	{
		Ybins[i]=his->GetYaxis()->GetBinLowEdge(i+1);
		//   printf("x:: %d %f\n",i,Ybins[i]);
	}

	for(i=0;i<=Nz;i++) 
	{
		Zbins[i]=his->GetZaxis()->GetBinLowEdge(i+1);
		//    printf("x:: %d %f\n",i,Zbins[i]);
	}

	switch(axis)
	{
		case 1:
			his2D=new TH2F("YZ plane","YZ",Ny,Ybins,Nz,Zbins);
			for(int i=1;i<=Ny;i++)
				for(int j=1;j<=Nz;j++)
					his2D->SetBinContent(i,j,his->GetBinContent(Bin1,i,j));
			his2D->GetXaxis()->SetTitle("y [#mum]");
			his2D->GetYaxis()->SetTitle("z [#mum]");
			break;
		case 2:
			his2D=new TH2F("XZ plane","XZ",Nx,Xbins,Nz,Zbins);
			for(int i=1;i<=Nx;i++)
				for(int j=1;j<=Nz;j++)
					his2D->SetBinContent(i,j,his->GetBinContent(i,Bin1,j));
			his2D->GetXaxis()->SetTitle("x [#mum]");
			his2D->GetYaxis()->SetTitle("z [#mum]");
			break;
		case 3:
			his2D=new TH2F("XY plane","XY",Nx,Xbins,Ny,Ybins);
			for(int i=1;i<=Nx;i++)
				for(int j=1;j<=Ny;j++)
					his2D->SetBinContent(i,j,his->GetBinContent(i,j,Bin1));
			his2D->GetXaxis()->SetTitle("x [#mum]");
			his2D->GetYaxis()->SetTitle("y [#mum]");
			break;
	}
	return his2D;
}

class KGeometry
{
	private:
	public:

		TH3I    *EG;        //electrode geometry
		TH3I    *DM;        //detector material
		Int_t nx;           //x-divisions
		Int_t ny;           //y-divisions
		Int_t nz;           //z-divisions 

		// Constructors of the class
		KGeometry();
		~KGeometry();
		KGeometry(TH3I *x){GetGrid(x,0); };
		KGeometry(TH3I *x, TH3I *y){GetGrid(x,0); GetGrid(x,1);};
		void GetGrid(TH3I *,Short_t =0); 
		void ElRectangle(Float_t *Pos, Float_t *Size, Int_t Wei, Int_t Mat);
		void ElCylinder(Float_t *Pos,Float_t R, Float_t L,Int_t O, Int_t Wei, Int_t Mat); 
		Int_t SetBoundaryConditions();
		TH3F *MapToGeometry(Double_t *, Double_t =1);
		TH3F *GetGeom();
		Float_t GetLowEdge(Int_t);
		Float_t GetUpEdge(Int_t);
		Double_t GetStepSize(Int_t, Int_t);
		//   ClassDef(KGeometry,1) 
};



// ClassImp(KGeometry)
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KGeometry                                                            //
//                                                                      //
// Class for description of detector geometry                           //
// The class defines the geometry of the detector.                      //
// It is based upon two TH3S histograms                                 // 
// It contains some functions defined to design the electrodes          //  
//                                                                      //
//////////////////////////////////////////////////////////////////////////

KGeometry::KGeometry()
{
	EG=NULL;
	DM=NULL;
	nx=1;
	ny=1;
	nz=1;
}

KGeometry::~KGeometry()
{
	if(EG!=NULL) delete EG;
	if(DM!=NULL) delete DM;
}

void KGeometry::GetGrid(TH3I *x, Short_t which)
{
	// EG grid
	//bit 1 = 1  -> 1st electrode
	//bit 2 = 2  -> 2nd electrode

	switch(which)
	{
		case 0: 
			if(EG!=NULL) delete EG;
			EG=new TH3I(); x->Copy(*EG);
			nx=EG->GetNbinsX();
			ny=EG->GetNbinsY();
			nz=EG->GetNbinsZ();
			break;
		case 1:
			if(DM!=NULL) delete DM;
			DM=new TH3I(); x->Copy(*DM); 
			break;
	}

	if(DM!=NULL) if(DM->GetNbinsX()!=nx) printf("Warning: dimenssions mismatch - X !\n");
	if(DM!=NULL) if(DM->GetNbinsY()!=ny) printf("Warning: dimenssions mismatch - Y !\n");
	if(DM!=NULL) if(DM->GetNbinsZ()!=nz) printf("Warning: dimenssions mismatch - Z !\n");

}

void KGeometry::ElRectangle(Float_t *Pos, Float_t *Size, Int_t Wei, Int_t Mat)
{
	// Sets Up 
	// Int_t i,j,k,q;
	Int_t i,j,k;

	Int_t xpl,ypl,zpl,xpr,ypr,zpr;
	//Set up left edge of the box
	xpl=EG->GetXaxis()->FindBin(Pos[0]-Size[0]);
	ypl=EG->GetYaxis()->FindBin(Pos[1]-Size[1]);
	zpl=EG->GetZaxis()->FindBin(Pos[2]-Size[2]);
	//Set up rigth edge of the box	   
	xpr=EG->GetXaxis()->FindBin(Pos[0]+Size[0]);
	ypr=EG->GetYaxis()->FindBin(Pos[1]+Size[1]);
	zpr=EG->GetZaxis()->FindBin(Pos[2]+Size[2]);
	//printf("(%d %d),(%d %d),(%d %d)\n",xpl,xpr,ypl,ypr,zpl,zpr);
	// Fill the geometry histogram
	for(k=zpl;k<=zpr;k++)	   
		for(j=ypl;j<=ypr;j++)
			for(i=xpl;i<=xpr;i++)
			{
				//	  if(x0+i>1 && x0+i<=nx && y0+i>1 && y0+i<=ny)
				//	  printf("%d %d %d %d\n",i,j,k,Wei);
				if(EG!=NULL) EG->SetBinContent(i,j,k,Wei);
				if(DM!=NULL) DM->SetBinContent(i,j,k,Mat);
			}

}


void KGeometry::ElCylinder(Float_t *Pos,Float_t R, Float_t L,Int_t O, Int_t Wei, Int_t Mat)
{
	// Cylindrical electrode 
	// Float_t *Pos;  - postion of the cone center 
	Float_t Dist,D,x,y,z,Bu,Bd;
	//    Int_t i,j,k,q;
	Int_t i,j,k;
	for(k=1;k<=nz;k++)
		for(j=1;j<=ny;j++)
			for(i=1;i<=nx;i++)
			{
				x=EG->GetXaxis()->GetBinCenter(i);
				y=EG->GetYaxis()->GetBinCenter(j);
				z=EG->GetZaxis()->GetBinCenter(k);

				switch(O)
				{
					case 3:   D=TMath::Sqrt(TMath::Power(x-Pos[0],2)+TMath::Power(y-Pos[1],2))-R; break;
					case 2:   D=TMath::Sqrt(TMath::Power(x-Pos[0],2)+TMath::Power(z-Pos[2],2))-R; break;
					case 1:   D=TMath::Sqrt(TMath::Power(y-Pos[1],2)+TMath::Power(z-Pos[2],2))-R; break;
				}

				if(D<=0)
				{
					switch(O)
					{
						case 3:   Dist=EG->GetZaxis()->GetBinCenter(k); 
							  Bu=Pos[2]+L; Bd=Pos[2]-L; break;
						case 2:   Dist=EG->GetYaxis()->GetBinCenter(j); 
							  Bu=Pos[1]+L; Bd=Pos[1]-L; break;
						case 1:   Dist=EG->GetXaxis()->GetBinCenter(i); 
							  Bu=Pos[0]+L; Bd=Pos[0]-L; break;
					}


					if(Dist<=Bu && Dist>=Bd)
					{
						if(EG!=NULL) EG->SetBinContent(i,j,k,Wei); 
						if(DM!=NULL) DM->SetBinContent(i,j,k,Mat);
					}

				}


			}
}

Int_t KGeometry::SetBoundaryConditions()
{
	Int_t i,j,k,val,cval,nval;
	if(EG==NULL) {printf("Please set the geometry first ! \n"); return -1;}
	nx=EG->GetNbinsX();
	ny=EG->GetNbinsY();
	nz=EG->GetNbinsZ();

	//Bit 1 = 1 -> GND - 0 V bias
	//Bit 2 = 2 -> Voltage (usual bias volage)
	//Bit 15-32 = 32768 -> Additional Votlages 

	//Bits determining boundary conditions:
	//bit 2 = 4  -> down val
	//bit 3 = 8  -> up val 
	//bit 4 = 16 -> left val
	//bit 5 = 32 -> rigth val
	//bit 6 = 64  -> down der
	//bit 7 = 128  -> up der 
	//bit 8 = 256 -> left der
	//bit 9 = 512 -> rigth der
	//the 3D section is separated
	//bit 10= 1024 -> out val
	//bit 11= 2048 -> in val
	//bit 12= 4096 -> out der
	//bit 13= 8192 -> in der
	//bit 14= 16384-> read out node
	for(k=1;k<=nz;k++)
	{
		for(j=1;j<=ny;j++)
		{
			for(i=1;i<=nx;i++)
			{
				cval=EG->GetBinContent(i,j,k);

				if(!(cval&1 || cval&2 || cval>=32768))
				{

					nval=0;
					if(i+1<=nx) 
					{
						val=EG->GetBinContent(i+1,j,k);
						if(val&1 || val&2 || val>=32768) nval|=32;
					}
					else nval|=512;

					if(i-1>0)   
					{
						val=EG->GetBinContent(i-1,j,k);
						if(val&1 || val&2 || val>=32768) nval|=16;
					}
					else nval|=256;

					if(j+1<=ny) 
					{

						val=EG->GetBinContent(i,j+1,k);
						if(val&1 || val&2 || val>=32768) nval|=8;
					}
					else nval|=128;

					if(j-1>0)   
					{

						val=EG->GetBinContent(i,j-1,k); 
						if(val&1 || val&2 || val>=32768)  nval|=4; 
					}
					else nval|=64;


					if(k+1<=nz)
					{
						val=EG->GetBinContent(i,j,k+1);
						if(val&1 || val&2 || val>=32768) nval|=2048;
					}
					else if(nz!=1) nval|=8192;

					if(k-1>0)   
					{
						val=EG->GetBinContent(i,j,k-1);
						if(val&1 || val&2 || val>=32768) nval|=1024;
					}
					else  if(nz!=1) nval|=4096;


					EG->SetBinContent(i,j,k,nval);
					//      	if(k==nz) printf("%d %d %d :: %d \n",i,j,k,nval);
				}

			}
		}
	}

	return 0;
}

TH3F *KGeometry::MapToGeometry(Double_t *x,Double_t Scale)
{
	TH3F *fhis=new TH3F();
	EG->Copy(*fhis);
	fhis->Reset();

	// Map the array of values: E, U, W ... to the geometry.
	int i,j,k,n;
	//   Double_t xb=EG->GetXaxis()->GetBinUpEdge(nx);
	//   Double_t yb=EG->GetYaxis()->GetBinUpEdge(ny);
	//   Double_t zb=EG->GetZaxis()->GetBinUpEdge(nz);
	//   Double_t bx=EG->GetXaxis()->GetBinLowEdge(1);
	//   Double_t by=EG->GetYaxis()->GetBinLowEdge(1);
	//   Double_t bz=EG->GetZaxis()->GetBinLowEdge(1);

	//   TH3F *fhis=new TH3F("Pot3d","Pot3d",nx,bx,xb,ny,by,yb,nz,bz,zb);

	for (k=1;k<=nz;k++)
		for (j=1;j<=ny;j++)
			for(i=1;i<=nx;i++)
			{
				n=(k-1)*nx*ny+(j-1)*nx+i; 
				fhis->SetBinContent(i,j,k,x[n]*Scale);
			}
	return fhis;
}

TH3F *KGeometry::GetGeom()
{
	// Map the array of values: E, U, W ... to the geometry.
	//
	int i,j,k,n,bin,col;
	TH3F *dhis=new TH3F();
	EG->Copy(*dhis);

	for (k=1;k<=nz;k++)
		for (j=1;j<=ny;j++)
			for(i=1;i<=nx;i++)
			{
				bin=dhis->GetBinContent(i,j,k);
				col=0;
				if(bin>=32768) col=1; 
				if(bin==1 || bin==2) col=2;
				if(bin==16385) col=3;		 
				dhis->SetBinContent(i,j,k,col);
			}
	//  dhis->Draw("glbox");
	return dhis;
}


Float_t KGeometry::GetUpEdge(Int_t dir)
{
	Float_t ret=0;
	switch(dir)
	{
		case 0: ret=EG->GetXaxis()->GetBinUpEdge(nx); break;
		case 1: ret=EG->GetYaxis()->GetBinUpEdge(ny); break;
		case 2: ret=EG->GetZaxis()->GetBinUpEdge(nz); break;
		default: printf("Index out of scope!\n"); ret=0; break;
	}
	return ret;
}

Float_t KGeometry::GetLowEdge(Int_t dir)
{
	Float_t ret=0;
	switch(dir)
	{
		case 0: ret=EG->GetXaxis()->GetBinLowEdge(1); break;
		case 1: ret=EG->GetYaxis()->GetBinLowEdge(1); break;
		case 2: ret=EG->GetZaxis()->GetBinLowEdge(1); break;
		default: printf("Index out of scope!\n"); ret=0; break;
	}
	return ret;
}
Double_t KGeometry::GetStepSize(Int_t dir, Int_t i)
{
  Double_t Lo,Hi;
  Double_t ret;
  switch(dir)
    {
    case 0:
      Hi=fabs(EG->GetXaxis()->GetBinCenter(i+1)-EG->GetXaxis()->GetBinCenter(i));
      Lo=fabs(EG->GetXaxis()->GetBinCenter(i)-EG->GetXaxis()->GetBinCenter(i-1));
      break;
    case 1:
       Hi=fabs(EG->GetYaxis()->GetBinCenter(i+1)-EG->GetYaxis()->GetBinCenter(i));
       Lo=fabs(EG->GetYaxis()->GetBinCenter(i)-EG->GetYaxis()->GetBinCenter(i-1));
      break;
    case 2:
       Hi=fabs(EG->GetZaxis()->GetBinCenter(i+1)-EG->GetZaxis()->GetBinCenter(i));
       Lo=fabs(EG->GetZaxis()->GetBinCenter(i)-EG->GetZaxis()->GetBinCenter(i-1));
      break;
    default:
      Hi=-1; Lo=-1;
      break;
    }
  ret=0.5*Hi+0.5*Lo;
  return ret;
}




// KField 

class KField
{
	private:
		Int_t Method;   // Method to calculate the intermediate points
		Int_t dim;
	public:
		TH3F *U;
		TH3F *Ex;
		TH3F *Ey;
		TH3F *Ez;
		TH3F *E;

		KField() {U=NULL; Ex=NULL; Ey=NULL; Ez=NULL;};
		~KField();
		Int_t CalField();
		static Float_t GetFieldPoint(Float_t *, Float_t *);
		TVector3 *CalFieldXYZ(Float_t x, Float_t y, Float_t z);
		void  CalFieldXYZ(Float_t x, Float_t y, Float_t z, Float_t *E);  
		void  CalFieldXYZ(Float_t , Float_t , Float_t , TVector3 *); 
		Float_t CalPotXYZ(Float_t x, Float_t y, Float_t z);
		Double_t DriftVelocity(Float_t E,Float_t Charg, Float_t T, Double_t Neff, Int_t which);
		Double_t Mobility(Float_t E,Float_t T,Float_t Charg,Double_t Neff, Int_t which);
		Float_t KInterpolate2D(TH3F *, Float_t ,Float_t, Int_t=3, Int_t=1);

		//   ClassDef(KField,1) 
};

KField::~KField()
{
	if(U!=NULL)  delete U; 
	if(Ex!=NULL) delete Ex; 
	if(Ey!=NULL) delete Ey; 
	if(Ez!=NULL) delete Ez;
}

Float_t KField::KInterpolate2D(TH3F *his, Float_t x, Float_t y, Int_t dir, Int_t bin)
{
	Int_t EX1,EX2,EY1,EY2;
	Float_t t,u,ret;
	TAxis *ax1,*ax2;
	Float_t v11,v21,v22,v12;

	ret=0;

	switch(dir)
	{
		case 3: ax1=his->GetXaxis(); ax2=his->GetYaxis(); break;
		case 2: ax1=his->GetXaxis(); ax2=his->GetZaxis(); break;
		case 1: ax1=his->GetYaxis(); ax2=his->GetZaxis(); break;
	}

	EX1=ax1->FindBin(x); 
	if(ax1->GetBinCenter(EX1)<=x) {EX2=EX1+1;} else {EX2=EX1; EX1--;}
	EY1=ax2->FindBin(y);
	if(ax2->GetBinCenter(EY1)<=y) {EY2=EY1+1;} else {EY2=EY1; EY1--;}

	if(EY2>ax2->GetNbins()) {u=0; EY2=ax2->GetNbins();} else
		if(EY1<1) {u=0; EY1=1;} else
			u=(y-ax2->GetBinCenter(EY1))/(ax2->GetBinCenter(EY2)-ax2->GetBinCenter(EY1));		

	if(EX2>ax1->GetNbins()) {t=0; EX2=ax1->GetNbins();} else
		if(EX1<1) {t=0; EX1=1;} else  
			t=(x-ax1->GetBinCenter(EX1))/(ax1->GetBinCenter(EX2)-ax1->GetBinCenter(EX1));


	//     printf("Points are:: %d %d %d %d [dir=%d, bin=%d] (t=%f u=%f)\n",EX1,EX2,EY1,EY2, dir, bin, t,u);

	switch(dir)
	{
		case 3: 
			v11=his->GetBinContent(EX1,EY1,bin); v21=his->GetBinContent(EX2,EY1,bin);
			v22=his->GetBinContent(EX2,EY2,bin); v12=his->GetBinContent(EX1,EY2,bin);
			break;
		case 2: 
			v11=his->GetBinContent(EX1,bin,EY1); v21=his->GetBinContent(EX2,bin,EY1);
			v22=his->GetBinContent(EX2,bin,EY2); v12=his->GetBinContent(EX1,bin,EY2);
			break;
		case 1: 
			v11=his->GetBinContent(bin,EX1,EY1); v21=his->GetBinContent(bin,EX2,EY1);
			v22=his->GetBinContent(bin,EX2,EY2); v12=his->GetBinContent(bin,EX1,EY2);
			break;
	}

	ret=(1-t)*(1-u)*v11;
	ret+=t*(1-u)*v21;
	ret+=t*u*v22;
	ret+=(1-t)*u*v12;
	return ret;
}


Int_t KField::CalField()
{
	Float_t X[3],Y[3],EE;
	Int_t q,i,j,k;

	if(U==NULL) 
	{printf("Can not calculate field - no potential array!"); return -1;};

	Int_t nx=U->GetNbinsX();
	Int_t ny=U->GetNbinsY();
	Int_t nz=U->GetNbinsZ();

	if(nz==1) {printf("2D field!\n"); dim=2;} else dim=3;



	Ex=new TH3F(); U->Copy(*Ex); Ex->Reset();

	Ey=new TH3F(); U->Copy(*Ey); Ey->Reset();

	Ez=new TH3F(); U->Copy(*Ez); Ez->Reset();

	E=new TH3F(); U->Copy(*E);  E->Reset();


	for(k=1;k<=nz;k++)
		for(j=1;j<=ny;j++)
			for(i=1;i<=nx;i++)
			{

				// Get X field
				if(i==1 || i==nx) Ex->SetBinContent(i,j,k,0); else 
				{
					for(q=0;q<=2;q++) 
					{
						X[q]=U->GetXaxis()->GetBinCenter(i+q-1);
						Y[q]=U->GetBinContent(i+q-1,j,k);
					}

					Ex->SetBinContent(i,j,k,GetFieldPoint(X,Y));
				}


				// Get Y field
				if(j==1 || j==ny) Ey->SetBinContent(i,j,k,0); else 
				{
					for(q=0;q<=2;q++) 
					{
						X[q]=U->GetYaxis()->GetBinCenter(j+q-1);
						
						Y[q]=U->GetBinContent(i,j+q-1,k);

					}
					//std::cout << "i:" << i << "j:" << j << "k:" << k << std::endl;
					//std::cout << "X[0]:" << X[0] << "X[1]:" << X[1] << "X[2]:" << X[2] << std::endl;
					//std::cout << "Y[0]:" << Y[0] << "Y[1]:" << Y[1] << "Y[2]:" << Y[2] << std::endl;

					Ey->SetBinContent(i,j,k,GetFieldPoint(X,Y));
					// std::cout << GetFieldPoint(X, Y) << std::endl;
				}
				// Get Z field
				if(k==1 || k==nz) Ez->SetBinContent(i,j,k,0); else 
				{
					for(q=0;q<=2;q++) 
					{
						X[q]=U->GetZaxis()->GetBinCenter(k+q-1);
						Y[q]=U->GetBinContent(i,j,k+q-1);
					}
					Ez->SetBinContent(i,j,k,GetFieldPoint(X,Y));
				}

				EE=TMath::Sqrt(TMath::Power(Ex->GetBinContent(i,j,k),2)+TMath::Power(Ey->GetBinContent(i,j,k),2)+TMath::Power(Ez->GetBinContent(i,j,k),2));
				//std::cout << "i:" << i << "j:" << j << "k:" << k << std::endl;
				//std::cout << "EX:"<<Ex->GetBinContent(i,j,k)<<",Ey:"<<Ey->GetBinContent(i,j,k)<<",EE:" <<EE<< std::endl;
				E->SetBinContent(i,j,k,EE);

			}

	return 0;
}

Float_t KField::GetFieldPoint(Float_t *X, Float_t *Y)
{
	//electric potential to field
	Float_t a,b,k12,k23;
	k12=(Y[0]-Y[1])/(X[0]-X[1]);
	k23=(Y[1]-Y[2])/(X[1]-X[2]);
	a=(k23-k12)/(X[2]-X[0]);
	b=k12-a*(X[0]+X[1]);
	//std::cout<< "b:" <<b<<"field:"<<2*a*X[1]+b<<std::endl;
	return 2*a*X[1]+b;
}


void  KField::CalFieldXYZ(Float_t x, Float_t y, Float_t z, Float_t *E)
{
	if(dim==2)
	{
		E[1]=KInterpolate2D(Ex,x,y);
		E[2]=KInterpolate2D(Ey,x,y);
		E[3]=0;
	}
	else
	{
		E[1]=Ex->Interpolate(x,y,z);
		E[2]=Ey->Interpolate(x,y,z);
		E[3]=Ez->Interpolate(x,y,z);
	}

	E[0]=TMath::Sqrt(E[1]*E[1]+E[2]*E[2]+E[3]*E[3]);
}

TVector3 *KField::CalFieldXYZ(Float_t x, Float_t y, Float_t z)
{
	Float_t E[4];
	CalFieldXYZ(x,y,z,E); 
	TVector3 *vec=new TVector3(E[1],E[2],E[3]);
	return vec;
}

void KField::CalFieldXYZ(Float_t x, Float_t y, Float_t z, TVector3 *vec)
{
	Float_t E[4];
	CalFieldXYZ(x,y,z,E); 
	vec->SetXYZ(E[1],E[2],E[3]);
}

Double_t KField::DriftVelocity(Float_t E,Float_t Charg, Float_t T, Double_t Neff, Int_t which)
{
	E*=1e4;//um->cm
	return(Mobility(E,T,Charg,Neff,which)*(Double_t)E);
}

Float_t KField::CalPotXYZ(Float_t x, Float_t y, Float_t z)
{
	Float_t ret=0;
	Int_t nx,ny,nz,bx,by,bz;
	if(dim==2) ret=KInterpolate2D(U,x,y);
	else 
	{
		nz=U->GetZaxis()->GetNbins(); 
		ny=U->GetYaxis()->GetNbins(); 
		nx=U->GetXaxis()->GetNbins(); 

		if( z >= U->GetZaxis()->GetBinCenter(nz) ) { bz=U->GetZaxis()->FindBin(nz); ret=KInterpolate2D(U,x,y,3,nz); } else
			if( y >= U->GetYaxis()->GetBinCenter(ny) ) { by=U->GetYaxis()->FindBin(ny); ret=KInterpolate2D(U,x,z,2,ny); } else
				if( x >= U->GetXaxis()->GetBinCenter(nx) ) { bx=U->GetXaxis()->FindBin(nx); ret=KInterpolate2D(U,y,z,1,nx); } else 
					if( z <= U->GetZaxis()->GetBinCenter(1)  ) { ret=KInterpolate2D(U,x,y,3,1);                                 } else
						if( y <= U->GetYaxis()->GetBinCenter(1)  ) { ret=KInterpolate2D(U,x,z,2,1);                                 } else
							if( x <= U->GetXaxis()->GetBinCenter(1)  ) { ret=KInterpolate2D(U,y,z,1,1);                                 } else 
								ret=U->Interpolate(x,y,z);
	}
	return ret;
}

Double_t KField::Mobility(Float_t E,Float_t T,Float_t Charg,Double_t Neff, Int_t which)
{
	Double_t lfm=0,hfm=0;
	Double_t vsatn,vsatp,vsat;
	Double_t betap,betan;
	Double_t alpha;
	//std::cout<<"mobility choose:"<<which<<std::endl;
	switch(which)
	{
		case 0:
			//printf("%e ",par[0]);
			if(Charg>0)
			{
				lfm=8.54e5*TMath::Power(T,-1.075)*TMath::Exp(1-T/124.);
				vsatp=1.445e7*TMath::Exp(-T/435.9);
				betap=2.49*TMath::Exp(-T/270.3);
				hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatp,1/betap),betap);
			}
			else
			{
				lfm=2.712e8*TMath::Power(T,-2.133);
				vsatn=1.586e7*TMath::Exp(-T/723.6);
				betan=-8.262e-8*TMath::Power(T,3)+6.817e-5*TMath::Power(T,2)-1.847e-2*T+2.429;
				hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatn,1/betan),betan);
			}
			break;
		//silicon carbide
		case 1:
			if (Charg > 0)
			{
				// unit cm
				alpha = 0.34;
				Double_t ulp = 124 * TMath::Power(T / 300, -2);
				Double_t uminp = 15.9;
				Double_t Crefp = 1.76e19;
				betap = 1.213 * TMath::Power(T / 300, 0.17);
				vsatp = 2e7 * TMath::Power(T / 300, 0.52);
				lfm = uminp + ulp/ (1 + TMath::Power(Neff / Crefp, alpha));
				hfm = lfm / (TMath::Power(1 + TMath::Power(lfm * E / vsatp, betap), 1 / betap));
			}
			else
			{
				alpha = 0.61;
				Double_t ulp = 947 * TMath::Power(T / 300, -2);
				Double_t uminp = 40.0 * TMath::Power(T / 300, -0.5);
				Double_t Crefp = 1.94e19;
				betap = 1 * TMath::Power(T / 300, 0.66);
				vsatp = 2e7 * TMath::Power(T / 300, 0.87);
				lfm = ulp/ (1 + TMath::Power(Neff / Crefp, alpha));
				hfm = lfm / (TMath::Power(1 + TMath::Power(lfm * E / vsatp, betap), 1/betap));
			}
			break;
		//silicon
		// case 1:       
		// 	alpha=0.72*TMath::Power(T/300,0.065);
		// 	if(Charg>0)
		// 	{
		// 		Double_t ulp=460*TMath::Power(T/300,-2.18);
		// 		Double_t uminp=45*TMath::Power(T/300,-0.45);      
		// 		Double_t Crefp=2.23e17*TMath::Power(T/300,3.2);
		// 		betap=1;
		// 		vsatp=9.05e6*TMath::Sqrt(TMath::TanH(312/T));
		// 		lfm=uminp+(ulp-uminp)/(1+TMath::Power(Neff/Crefp,alpha));
		// 		hfm=2*lfm/(1+TMath::Power(1+TMath::Power(2*lfm*E/vsatp,betap),1/betap));
		// 	}
		// 	else
		// 	{
		// 		Double_t uln=1430*TMath::Power(T/300,-2); 
		// 		Double_t uminn=80*TMath::Power(T/300,-0.45);
		// 		Double_t Crefn=1.12e17*TMath::Power(T/300,3.2);
		// 		betan=2;      
		// 		vsatn=1.45e7*TMath::Sqrt(TMath::TanH(155/T));
		// 		lfm=uminn+(uln-uminn)/(1+TMath::Power(Neff/Crefn,alpha));
		// 		hfm=2*lfm/(1+TMath::Power(1+TMath::Power(2*lfm*E/vsatn,betan),1/betan));
		// 	}
		// 	break;
 
		case 2:   // WF2
			if(Charg>0)
			{
				lfm=480;
				vsatp=9.5e6;
				betap=1;
				hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatp,1/betap),betap);
			}
			else
			{
				lfm=1350;
				vsatn=1.1e7;
				betan=0.5;
				hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatn,1/betan),betan);
			}
			break; 
		case 3:  // Klanner Scharf
			Double_t bb,cc,E0;

			if(Charg>0)
			{
				E0=2970*TMath::Power(T/300,5.63);
				bb=9.57e-8*TMath::Power(T/300,-0.155);
				cc=-3.24e-13;
				lfm=457*TMath::Power(T/300,-2.80);
				if(E>E0) hfm=1./(1/lfm+bb*(E-E0)+cc*TMath::Power(E-E0,2)); else hfm=lfm;
			}
			else
			{
				E0=2970*TMath::Power(T/300,5.63);
				lfm=1430*TMath::Power(T/300,-1.99);
				vsatn=1.05e7*TMath::Power(T/300,-0.302); //corrected by suggestions of Scharf
				if(E>E0) hfm=1./(1/lfm+1/vsatn*(E-E0)); else hfm=lfm;
			}
			break;
		case 4:   //Jacoboni
			if (Charg>0)
			{
				lfm = 474 * TMath::Power(T/300., -2.619);
				vsatp = 0.940e7  * TMath::Power(T/300., -0.226);
				betap = 1.181 * TMath::Power(T/300., 0.633 ); // <100> orientation
				hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatp,betap),1/betap);
			} 
			else 
			{
				lfm = 1440*  TMath::Power(T/300., -2.260);
				vsatn = 1.054e7  *  TMath::Power(T/300., -0.602);
				betan = 0.992 *  TMath::Power(T/300., 0.572); // <100> orientation
				hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatn,betan),1/betan);
			}
			break;
		case 5:   //Jacoboni
			if (Charg>0)
			{
				lfm = 474 * TMath::Power(T/300., -2.619);
				vsatp = 0.940e7  * TMath::Power(T/300., -0.226);
				betap = 1.181 * TMath::Power(T/300., 0.633 ); // <100> orientation
				hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatp,betap),1/betap);
			} 
			else 
			{
				lfm = 1440*  TMath::Power(T/300., -2.260);
				vsatn = 1.054e7  *  TMath::Power(T/300., -0.602);
				betan = 0.992 *  TMath::Power(T/300., 0.572); // <100> orientation
				hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatn,betan),1/betan);
			}
			break;


		case 9:
			if(Charg>0) {hfm=0;} else {hfm=0;}
			break;
		case 10: //Diamond parametrization
			if(Charg>0) {lfm=2064; vsat=14.1e6;} else {lfm=1714; vsat=9.6e6;};
			hfm=lfm/(1+(lfm*E)/vsat);
			break;
	}
	return hfm; 
}




// KDetector 

#include "TRandom.h"
#include "TF3.h"

class KDetector : public KGeometry, public KMaterial { 


	private:
		Double_t Deps;
		TRandom *ran;               //random number generator
		Double_t CalErr;               //Error of the solver
		Int_t MaxIter;              //Maximum number of iterations in eq solver
		Short_t Debug;              //Print information of drift calculation etc.

	public:
		Float_t Voltage;  //Voltage
		Float_t Voltage2; //Voltage2 
		TArrayF Voltages; //Array of voltages

		// Definition of space charge 
		TF3     *NeffF;     //effective dopping concentration function
		TH3F    *NeffH;     //effective dopping concentration histogram

		// Weigthing, electric and magnetic field
		KField *Ramo;       // ramo field 
		KField *Real;       // electric field
		Float_t B[3];      // magnetic field

		// Trapping and variables used for multiplication studies
		Float_t taue;      // effective trapping time constants 
		Float_t tauh;      // effective trapping time constants 
		TF3 *TauE;         // Function of TauE(x,y,z);
		TF3 *TauH;         // Function of TauH(x,y,z);

		Int_t BreakDown;     // if break down occurs it goes to 1 otherwise is 0
		Float_t MTresh;      // treshold for taking multiplication into account
		Float_t BDTresh;     // hole multiplication - break down treshold 
		Float_t DiffOffField;// electric field where diffusion is switched off [V/um] !!
		// Drift parameters

		Float_t enp[3];      //entry point for the charge drift
		Float_t exp[3];      //exit point for the cahrge drift
		Int_t diff;          // Diffusion simulation (yes=1, no=0)
		Int_t average;       // Average (over how many events)
		Float_t SStep;       // Simulation step size;
		Float_t MaxDriftLen; // Maximum drift lenght before stopping the drift

		// Output histograms
		TH1F *pos;           // contribution of the holes to the total drift current
		TH1F *neg;           // contribution of the electrons  to the total drift current
		TH1F *sum;	       // total drift current
		TH1F *pairs;		// e-h pairs each micros
		// Constructors and destructor
		KDetector();
		~KDetector();
		//_______________________________________________________________________________
		void SetDriftHisto(Float_t x,Int_t=200);
		// start declaration followed by solving Poisson's equation.
		void CalField(Int_t);         
		void Declaration(Int_t);                 // declaration of boundary conditions
		Double_t V(int ,int);                    // defining voltage
		Double_t kappa(int ,int , int , int);    // defining space charge

		void ShowMipIR(Int_t, Int_t=14, Int_t=1);
		void ShowUserIonization(Int_t, Float_t *, Float_t *, Float_t *, Float_t *, Int_t=14, Int_t=1);
		void MipIR(Int_t = 20, Float_t = 0);
		void Drift(Double_t, Double_t, Double_t, Float_t, KStruct *, Double_t = 0);
  		void SetEntryPoint(Float_t x, Float_t y, Float_t z) {enp[0]=x; enp[1]=y; enp[2]=z;};
  		void SetExitPoint(Float_t x, Float_t y, Float_t z) {exp[0]=x; exp[1]=y; exp[2]=z;};


		//Configuration functions 
		// ClassDef(KDetector,1) 
};

#include "TPolyLine3D.h"
// #include "KDetector.h"
#include "TFile.h"

#define ABS(x) x>0?x:-x
#define PREDZNAK(x) x>0?1:-1

#define C1(x,y,z) y3[n]=x+y+z;

#define L1(x) y2[n]=x;
#define R1(x) y4[n]=x;
#define U1(x) y5[n]=x;
#define D1(x) y6[n]=x;
#define I1(x) y7[n]=x;
#define O1(x) y8[n]=x;

#define L2(x) y2[n]=2*x;
#define R2(x) y4[n]=2*x;
#define U2(x) y5[n]=2*x;
#define D2(x) y6[n]=2*x;
#define I2(x) y7[n]=2*x;
#define O2(x) y8[n]=2*x;

#define C0 y3[n]=1.;
#define U0 y5[n]=0.;
#define D0 y6[n]=0.;
#define R0 y4[n]=0.;
#define L0 y2[n]=0.;
#define I0 y7[n]=0.;
#define O0 y8[n]=0.;

#define  PI  3.1415927
#define EPS 1.0e-14
#define STEP_DET_CUR 25e-9

// ClassImp(KDetector)


// double **a,*b,*y6,*y2,*y3,*y4,*y5,*y7,*y8;


void KDetector::SetDriftHisto(Float_t x, Int_t numbins)
{
	if(x<0 || x>10000e-9) 
		printf("Selected time range out of scope !\n"); 
	else
	{
		if(pos!=NULL) delete pos;
		pos  = new TH1F("charge+","Positive Charge",numbins,0,x);
		if(neg!=NULL) delete neg;
		neg  = new TH1F("charge-","Negative Charge",numbins,0,x); 	
		if(sum!=NULL) delete sum;
		sum  = new TH1F("charge","Total Charge",numbins,0,x); 
		if(pairs!=NULL) delete pairs;
		pairs = new TH1F("pairs", "pairs", numbins/5, 0, 150);
		sum->SetXTitle("t [s]");  neg->SetXTitle("t [s]"); pos->SetXTitle("t [s]"); pairs->SetXTitle("e-h pairs/um");
		sum->SetYTitle("current [A]");neg->SetYTitle("current [A]");pos->SetYTitle("current [A]"); pairs->SetXTitle("events");
		sum->GetYaxis()->SetTitleOffset(1.4); neg->GetYaxis()->SetTitleOffset(1.4); pos->GetYaxis()->SetTitleOffset(1.4);
		pairs->GetYaxis()->SetTitleOffset(1.4);
		pos->SetLineColor(2);
		neg->SetLineColor(4);
		sum->SetLineColor(1);
		pairs->SetLineColor(6);
	}
}



KDetector::KDetector()
{
	//////////////////////////////////////////////////////////////////////
	// Author: Gregor Kramberger                                        //
	// Default constructor for KDetector class                           //
	// Default = no space charge -> default space charge is function     //
	//////////////////////////////////////////////////////////////////////

	NeffF=new TF3("Profile","x[0]*x[1]*x[2]*0+[0]",0,1000,0,1000,0,1000);
	NeffF->SetParameter(0,0);
	NeffH=NULL; //Neff histogram se to NULL

	//set magnetic field to zero
	for(Int_t i=0;i<3;i++) B[i]=0;

	// setting up default random generator - diffusion
	ran=new TRandom3(0);  

	// Calculation parameters
	CalErr=1e-6;
	MaxIter=2000;

	// histograms for storing the drift
	pos=NULL; neg=NULL; sum=NULL; pairs=NULL;
	SetDriftHisto(25e-9);

	// setting up general variables
	// multiplication 
	TauE=NULL;
	TauH=NULL;
	taue=-1;   //no electron trapping 
	tauh=-1;   //no hole trapping 
	MTresh=-1; //no multiplication 
	BDTresh=-1;
	DiffOffField=8;  // critical field for diffusion to be switched off

	// drift

	Deps=1e-5;       //precision of tracking
	MaxDriftLen=1e9; // maximum driftlenght in [um]

	//MobMod=1;  //Mobility parametrization
	average=1; //average over waveforms
	diff=0;    //diffusion
	SStep=1;   //Step size of simulation
	Temperature=263; //temperature
	BreakDown=0; // no breakdown
	Debug=0;     // bo printing of debug information
	Voltage2=0;
	Ramo=new KField();
	Real=new KField();

}


KDetector::~KDetector()
{
	//Destructor of the detector class

	if(NeffF!=NULL) delete NeffF;
	if(NeffH!=NULL) delete NeffH;
	if(ran!=NULL) delete ran;
	if(pos!=NULL) delete pos;
	if(pairs!=NULL) delete pairs;
	if(neg!=NULL) delete neg;
	if(sum!=NULL) delete sum;

}

void KDetector::CalField(Int_t what)
{
	double err;
	int iteracije,i;
	int num=nx*ny*nz,dim[3];
	Double_t *x;
	//booking memory
	b=dvector(1,num); y6=dvector(1,num); y2=dvector(1,num); 
	y3=dvector(1,num); y4=dvector(1,num); y5=dvector(1,num);
	y7=dvector(1,num); y8=dvector(1,num);
	x=dvector(1,num);    
	// Setting up the boundary conditions
	printf("Setting up matrix ... \n");
	Declaration(what);
	printf("Solving matrix ...\n");
	// matrix solving
	for(i=1;i<=num;i++) x[i]=1.;
	dim[0]=nx; dim[1]=ny; dim[2]=nz;
	linbcg(num,dim,b,x,1,CalErr,MaxIter,&iteracije,&err);
	// Calculating the field
	if(!what)
	{
		Real->U=MapToGeometry(x);
		Real->CalField();
	}
	else
	{
		// Scale back to 1 from 1000 when storing the potential
		Ramo->U=MapToGeometry(x,1e-4);

		Ramo->CalField();
	}
	// Freeing the memory
	free_dvector(x, 1,num);   free_dvector(b, 1,num);   free_dvector(y2, 1,num); 
	free_dvector(y3, 1,num);  free_dvector(y4, 1,num);  free_dvector(y5, 1,num); 
	free_dvector(y6, 1,num);  free_dvector(y7, 1,num);  free_dvector(y8, 1,num);
}

void KDetector::Declaration(Int_t dowhat)
{
	// New declaration for a general detector class 
	Int_t i,j,k,val;
	Double_t Rd,Ld,Dd,Ud,Od,Id;
	Double_t PRd,PLd,PDd,PUd,POd,PId;

	Double_t Xr=0,Yr=0,Zr=0;
	Double_t Xc=0,Yc=0,Zc=0;
	Double_t Xl=0,Yl=0,Zl=0;
	Int_t ii,jj,kk;
	Double_t fac;

	long n=0;
	Int_t num=nx*ny*nz;

	for (k=1;k<=nz;k++)
		for (j=1;j<=ny;j++)
			for(i=1;i<=nx;i++)
			{
				n=(k-1)*nx*ny+(j-1)*nx+i; //Get index of the matrix element
				if(j-1<1) jj=1;  else jj=j-1;
				if(i-1<1) ii=1;  else ii=i-1; 
				if(k-1<1) kk=1;  else kk=k-1; 

				/////////// DEFINE STEPS IN X //////////////////////////////////////
				Rd=fabs(EG->GetXaxis()->GetBinCenter(i+1)-EG->GetXaxis()->GetBinCenter(i));
				Ld=fabs(EG->GetXaxis()->GetBinCenter(i)-EG->GetXaxis()->GetBinCenter(i-1));
				if(i+1>nx) Rd=Ld; if(i-1<1) Ld=Rd;

				////////// DEFINE PEMITIVITY IN X - normal surface ////////////////////////////
				PRd=Perm(DM->GetBinContent(i,j,k))+Perm(DM->GetBinContent(i,jj,k));

				if(nz!=1) 
				{
					PRd+=Perm(DM->GetBinContent(i,j,kk))+Perm(DM->GetBinContent(i,jj,kk));
					PRd/=4;
				} else PRd/=2;

				PLd=Perm(DM->GetBinContent(ii,j,k))+Perm(DM->GetBinContent(ii,jj,k));
				if(nz!=1) 
				{
					PLd+=Perm(DM->GetBinContent(ii,j,kk))+Perm(DM->GetBinContent(ii,jj,kk));
					PLd/=4;
				} else PLd/=2;

				/////////// DEFINE STEPS IN Y //////////////////////////////////////

				Ud=fabs(EG->GetYaxis()->GetBinCenter(j+1)-EG->GetYaxis()->GetBinCenter(j));
				Dd=fabs(EG->GetYaxis()->GetBinCenter(j)-EG->GetYaxis()->GetBinCenter(j-1));
				if(j+1>ny) Ud=Dd; if(j-1<1) Dd=Ud;

				////////// DEFINE PEMITIVITY IN Y ////////////////////////////
				PUd=Perm(DM->GetBinContent(i,j,k))  +Perm(DM->GetBinContent(ii,j,k));
				if(nz!=1) 
				{
					PUd+=Perm(DM->GetBinContent(i,j,kk))+Perm(DM->GetBinContent(ii,j,kk));
					PUd/=4;
				} else PUd/=2;

				PDd=Perm(DM->GetBinContent(i,jj,k))+Perm(DM->GetBinContent(ii,jj,k));
				if(nz!=1) 
				{
					PDd+=Perm(DM->GetBinContent(i,jj,kk))+Perm(DM->GetBinContent(ii,jj,kk));
					PDd/=4;
				} else PDd/=2;

				/////////// DEFINE STEPS IN Z //////////////////////////////////////

				Od=fabs(EG->GetZaxis()->GetBinCenter(k+1)-EG->GetZaxis()->GetBinCenter(k));
				Id=fabs(EG->GetZaxis()->GetBinCenter(k)-EG->GetZaxis()->GetBinCenter(k-1));
				if(k+1>nz) Od=Id; 
				if(k-1<1) Id=Od;

				//////////DEFINE PEMITIVITY IN Z ////////////////////////////
				if(nz!=1)
				{
					POd=Perm(DM->GetBinContent(i,jj,k))+Perm(DM->GetBinContent(i,j,k))+
						Perm(DM->GetBinContent(ii,j,k))+Perm(DM->GetBinContent(ii,jj,k));

					PId=Perm(DM->GetBinContent(i,jj,kk))+Perm(DM->GetBinContent(i,j,kk))+
						Perm(DM->GetBinContent(ii,j,kk))+Perm(DM->GetBinContent(ii,jj,kk));

					POd/=4;
					PId/=4;

				}
				///////////

				if(dowhat==1) {PRd=1; PLd=1; PUd=1; PDd=1; POd=1; PId=1;}

				Xr=PRd/(0.5*Rd*(Rd+Ld));
				Xl=PLd/(0.5*Ld*(Rd+Ld));  Xc=-(Xr+Xl);
				Yr=PUd/(0.5*Ud*(Ud+Dd));
				Yl=PDd/(0.5*Dd*(Ud+Dd));  Yc=-(Yr+Yl);

				if(nz!=1)
				{
					Zr=POd/(0.5*Od*(Od+Id));
					Zl=PId/(0.5*Id*(Od+Id));  Zc=-(Zr+Zl);
				} 

				//////////

				b[n]=0.;
				val=EG->GetBinContent(i,j,k);

				if(nz==1) { C1(Xc,Yc,0) I0 O0 } 
				else      { C1(Xc,Yc,Zc) I1(Zl) O1(Zr) }

				R1(Xr) U1(Yr) L1(Xl) D1(Yl) 


					if(val&4)            {D0 b[n]-=V(EG->GetBinContent(i,j-1,k),dowhat)*Yl;}
				if(val&8)            {U0 b[n]-=V(EG->GetBinContent(i,j+1,k),dowhat)*Yr;}
				if(val&16)           {L0 b[n]-=V(EG->GetBinContent(i-1,j,k),dowhat)*Xl;}
				if(val&32)           {R0 b[n]-=V(EG->GetBinContent(i+1,j,k),dowhat)*Xr;}
				if(val&1024)         {I0 b[n]-=V(EG->GetBinContent(i,j,k-1),dowhat)*Zl;}
				if(val&2048)         {O0 b[n]-=V(EG->GetBinContent(i,j,k+1),dowhat)*Zr;}


				if(val&64)           {U2(Yr) D0 if(val&8)    {U0 b[n]-=V(EG->GetBinContent(i,j+1,k),dowhat)*Yr;}}
				if(val&128)          {D2(Yl) U0 if(val&4)    {D0 b[n]-=V(EG->GetBinContent(i,j-1,k),dowhat)*Yl;}}
				if(val&256)          {R2(Xr) L0 if(val&32)   {R0 b[n]-=V(EG->GetBinContent(i+1,j,k),dowhat)*Xr;}}
				if(val&512)          {L2(Xl) R0 if(val&16)   {L0 b[n]-=V(EG->GetBinContent(i-1,j,k),dowhat)*Xl;}}
				if(val&4096)         {O2(Zr) I0 if(val&2048) {O0 b[n]-=V(EG->GetBinContent(i,j,k+1),dowhat)*Zr;}}
				if(val&8192)         {I2(Zl) O0 if(val&1024) {I0 b[n]-=V(EG->GetBinContent(i,j,k-1),dowhat)*Zl;}}


				b[n]-=kappa(i,j,k,dowhat);
				if(val&1 || val&2 || val>=32768)   {U0 D0 L0 R0 C0 O0 I0 b[n]=V(val,dowhat);}   
				//if(j<=2 && i<=2) printf("stevilki: i=%d, j=%d, k=%d X=(%f %f ::%f %f), Y(%f %f :: %f %f), Z(%f %f :: %f %f) y[2,3,4,5,6,7,8]=%f %f %f %f %f %f %f :: b[n]=%f :: %d\n",i,j,k,Xr,Xl,Ld,Rd,Yr,Yl,Dd,Ud,Zr,Zl,Id,Od,y2[n],y3[n],y4[n],y5[n],y6[n],y7[n],y8[n],b[n],Mat);
				//      if(k==nz && (j==2 || j==ny-1)) printf("stevilki: i=%d, j=%d, k=%d X=(%f %f ::%f %f), Y(%f %f :: %f %f), Z(%f %f :: %f %f) y[2,3,4,5,6,7,8]=%f %f %f %f %f %f %f :: b[n]=%f\n",i,j,k,Xr,Xl,Ld,Rd,Yr,Yl,Dd,Ud,Zr,Zl,Id,Od,y2[n],y3[n],y4[n],y5[n],y6[n],y7[n],y8[n],b[n]);
			}

}

double KDetector::V(int val, int dowhat)
{
	double voltage;
	int k=0;
	if(dowhat==0) 
	{
		if(val&1) voltage=0;  
		if(val&2) voltage=Voltage;
		if(val & 32768) 
			//if bit 15 is on - many voltages
			voltage=Voltages[val>>16];
	}
	else
	{
		// numerical calculation converges faster if 1000 is used instead of 1
		// therefore the potential is scaled after calculation to 1
		if(val&16384) voltage=10000; else voltage=0;
	}
	return voltage;
}

Double_t KDetector::kappa(int i,int j, int k,  int dowhat )
{
	//Sets the effective space charge values for given point in the mesh!
	Double_t x,y,z,ret;

	//  if(NeffF!=NULL && NeffH!=NULL) printf("Warning:: Histogram values will be taken for Neff!\n");

	//Position in space
	x=EG->GetXaxis()->GetBinCenter(i);
	y=EG->GetYaxis()->GetBinCenter(j);
	z=EG->GetZaxis()->GetBinCenter(k);

	if(DM!=NULL) KMaterial::Mat=DM->GetBinContent(i,j,k); else KMaterial::Mat=0;

	if (dowhat==0) 
	{
		if(NeffF!=NULL)  // Neff=v enotah [um-3]
			//            ret=(NeffF->Eval(x,y,z)*1e6*e_0)/(KMaterial::Perm()*perm0); /*printf("i=%d,j=%d,y=%e\n",i,j,y);*/
			ret=(NeffF->Eval(x,y,z)*1e6*e_0)/(perm0); /*printf("i=%d,j=%d,y=%e\n",i,j,y);*/

		if(NeffH!=NULL)
			//   ret=(NeffH->GetBinContent(i,j,k)*1e6*e_0)/(KMaterial::Perm()*perm0);
			ret=(NeffH->GetBinContent(i,j,k)*1e6*e_0)/(perm0);

	}

	//if (dowhat==0) if(j>nc) y=(Step*Step*1e-12)/(Si_mue*Ro*perm*perm0); else y=-(Step*Step*1e-12)/(Si_mue*Ro*perm*perm0);
	else 
		ret=0.;
	return ret;
}

void KDetector::ShowMipIR(Int_t div, Int_t color,Int_t how)
{
	Float_t *x=new Float_t [div];
	Float_t *y=new Float_t [div];
	Float_t *z=new Float_t [div];
	Float_t *Q=new Float_t [div];

	for(Int_t i=0;i<div;i++) 
	{
		x[i]=((exp[0]-enp[0])/div)*i+enp[0]+(exp[0]-enp[0])/(2*div);
		y[i]=((exp[1]-enp[1])/div)*i+enp[1]+(exp[1]-enp[1])/(2*div);
		z[i]=((exp[2]-enp[2])/div)*i+enp[2]+(exp[2]-enp[2])/(2*div);
	}
	ShowUserIonization(div, x,y,z,Q,color,how);
}

void KDetector::ShowUserIonization(Int_t div, Float_t *x, Float_t *y, Float_t *z, Float_t *Q, Int_t color,Int_t how)
{
	// The simulation of the drift for the minimum ionizing particles. 
	// A track is devided into Int_ div buckets. Each bucket is drifted in the field. The
	// induced currents for each carrier is calculated as the sum  all buckets. 

	TGraph *gr;
	TPolyLine3D *gr3D;
	Float_t sp[3];
	int i,j;
	KStruct seg;  //???
	TLine *line;
	TH3F *hh;

	// Draw histograms 

	if(EG!=NULL) 
	{ 
		
		if(nz==1) KHisProject(EG,3,how)->Draw("COL");
		else
		{
			//TH3F *hh=GetGeom();
			
			hh=GetGeom();
			hh->SetFillColor(color);
			hh->Draw("iso");
		}
	} 

	// Draw drift paths

	for(i=0;i<div;i++) 
	{

		// hole drift
		if(Debug)   printf("Entry Point: %f %f %f \n",x[i],y[i],z[i]);
		Drift(x[i],y[i],z[i],1,&seg);

		if(nz==1)
			gr=new TGraph(seg.Steps,&seg.Xtrack[1],&seg.Ytrack[1]); 
		else
			gr3D=new TPolyLine3D(seg.Steps,&seg.Xtrack[1],&seg.Ytrack[1],&seg.Ztrack[1]); 	

		if(nz==1)
		{
			gr->SetLineColor(2);
			gr->SetLineStyle(1);
			gr->Draw("L");
		}
		else
		{
			gr3D->SetLineStyle(1);  
			gr3D->SetLineColor(2); 
			gr3D->Draw("SAME"); 
		}

		// electron drift

		Drift(x[i],y[i],z[i],-1,&seg);

		if(nz==1)
			gr=new TGraph(seg.Steps,&seg.Xtrack[1],&seg.Ytrack[1]); 
		else
		{
			gr3D=new TPolyLine3D(seg.Steps,&seg.Xtrack[1],&seg.Ytrack[1],&seg.Ztrack[1]); 	
		}

		if(nz==1)
		{
			gr->SetLineColor(4);
			gr->SetLineStyle(3); 
			gr->Draw("L");
		}
		else
		{
			gr3D->SetLineStyle(1);  
			gr3D->SetLineColor(4); 
			gr3D->Draw("SAME"); 
		}

	}
}

void KDetector::MipIR(Int_t div, Float_t lambda)
{
	// The simulation of the drift for the laser induced carriers - attenuation of the light.
	// A track is devided into Int_ div buckets. Each bucket is drifted in the field. The
	// induced currents for each carrier is calculated as the sum  all buckets.
	//	Int_t MobMod; mobility model
	//	Float_t B; magnetic field

	Double_t data[4];
	Float_t sp[3],sp_2[3];
	Float_t drift_d; //drift distance each bucket
	Float_t n_pairs=0,total_pairs=0;
	Double_t loss_energy=0;
	Float_t Io, len = 0, lent = 0, lenlam, scalef = 0, pscalef = 0, mule, mulh;
	int i, j, k, e,mm;
	KStruct seg, segmul;
	// silicon
		
	if (Perm(0)== Float_t(9.76) ) loss_energy=Loss_energy(1);
	if (Perm(0) == Float_t(11.7) ) loss_energy = Loss_energy(0);
		//std::cout << "loss_energy:" << loss_energy << std::endl;
	TH1F *histop = new TH1F("ch-", "charge+", pos->GetNbinsX(), pos->GetXaxis()->GetXmin(), pos->GetXaxis()->GetXmax()); //2.23 -> 3.0
	TH1F *histon = new TH1F("ch+", "charge-", neg->GetNbinsX(), neg->GetXaxis()->GetXmin(), neg->GetXaxis()->GetXmax());
	sum->Reset();
	pos->Reset();
	neg->Reset();
	pairs->Reset();


	/////////////////////////////////////////////////////////////////
	// If there is no attenuation coefficient we consider that as mip
	/////////////////////////////////////////////////////////////////
	TFile Myf("Geant_Vin.root");
	TH1F *EDist;
	EDist = (TH1F *)Myf.Get("Silicon_Vin");

	if (lambda != 0)
	{
		for (k = 0; k < 3; k++)
			lent += TMath::Power(exp[k] - enp[k], 2);
		lent = TMath::Sqrt(lent);
		lenlam = lent / lambda;
		Io = div / (lambda * (1 - TMath::Exp(-lenlam)));

			//     if(Debug)
			//  printf("Calculated length=%f um, lenDIVlamda=%f, I0=%f \n", lent,lenlam,Io );
	}
	/////////////////////////////////////////////////////////////////

	for (i = 0; i < div-1; i++)
	{
		for (j = 0; j < 3; j++)
		{
			sp[j] = ((exp[j] - enp[j]) / div) * i + enp[j] + (exp[j] - enp[j]) / (2 * div);
		}
		for (j = 0; j < 3; j++)
		{
			sp_2[j] = ((exp[j] - enp[j]) / div) * (i+1) + enp[j] + (exp[j] - enp[j]) / (2 * div);
		}
		//     printf("#i=%d div=%d, pointx=%f, pointy=%f pointz=%f \n",i,div,sp[0],sp[1],sp[2]);
		drift_d = TMath::Sqrt(TMath::Power(sp_2[0]-sp[0], 2) + TMath::Power(sp_2[1]-sp[1], 2)+ TMath::Power(sp_2[2]-sp[2], 2));
		gRandom = new TRandom3(0);
		float ran_pairs = EDist->GetRandom();
		n_pairs =drift_d*ran_pairs*1e6/(5*loss_energy)*1.582;
		total_pairs+=n_pairs;
		//1.582 is change silicon loss energy to silicon carbide 
		// if (i<10) std::cout<<"drift_d:"<<drift_d<<std::endl;
		// if (i<10) std::cout<<"loss_energy:"<<loss_energy<<std::endl;
		// if (i<10) std::cout<<"n_pairs:"<<n_pairs<<std::endl;
			//if (mm==0) std::cout<<"n_pairs:"<<n_pairs<<std::endl;
		for (j = 0; j < average; j++)
		{

			Drift(sp[0], sp[1], sp[2], 1, &seg);
			seg.GetCH(histop, 1, 1, tauh);
			Drift(sp[0], sp[1], sp[2], -1, &seg);

			if (MTresh > 1)
			{
				mule = seg.GetCHMult(histon, 1, 1, taue); // performing multiplication
															//      if(Debug)  printf(":: Mstep = %f ::",mule);
			}
			else
				seg.GetCH(histon, 1, 1, taue);

			if (mule > MTresh && MTresh > 1) //if the multiplication is large enough then do the hole tracking
			{
				for (e = 1; e < seg.Steps + 1; e++)
				{
					Drift(seg.Xtrack[e], seg.Ytrack[e], seg.Ztrack[e], 1, &segmul, seg.Time[e]);
					mulh = segmul.GetCHMult(histop, 1, seg.MulCar[e], tauh);
					if (mulh > BDTresh && Debug)
					{
						printf("HOLE MULTIPLICATION - BREAKDOWN\n");
						BreakDown = 1;
					}
				}
			}
		}
		histop->Scale(n_pairs);
		histon->Scale(n_pairs);

		///////////////////////////////////////////////////////////////////
		// If there is no attenuation coefficient we consider that as mip//
		// Changes made by GK 3.3.2020                                 ////
		///////////////////////////////////////////////////////////////////
		if (lambda != 0)
		{
			len = 0;
			for (k = 0; k < 3; k++)
				len += TMath::Power(sp[k] - enp[k], 2);
			len = TMath::Sqrt(len);
			scalef = Io * TMath::Exp(-len / lambda) * SStep;
			//printf("step=%d ::  len=%f um, scalef=%f pscalef=%f\n",i,len,scalef,pscalef);
			histon->Scale(scalef);
			histop->Scale(scalef);
			pscalef += scalef;
			histop->Scale(lent / div);
			histon->Scale(lent / div);
		}
		/////////////////////////////////////////////////////////////////////////////////
		pairs->Fill(ran_pairs*1e6/(5*loss_energy)*1.582);
		pos->Add(histop);
		neg->Add(histon);
		histop->Reset();
		histon->Reset();
		
	}

	/// total laudau distribution
	float d = TMath::Sqrt(TMath::Power((exp[0] - enp[0]), 2) + TMath::Power((exp[1] - enp[1]), 2) + TMath::Power((exp[2] - enp[2]), 2));
	float mpv = (0.027 * log(d) + 0.126)*1.582; //keV/micron ==> log base "e"
	float FWHM = 0.31 * pow(d, -0.19);	   // this is the FWHM keV/micron, 4 times the sigma parameter in root.
	float DA = FWHM / 4.;
	TRandom3 Lan;
	TDatime *time; // current time
	time = new TDatime();
	Lan.SetSeed(time->TDatime::GetTime());
	float LanConst = 0;
	LanConst = Lan.Landau(mpv, DA);

	if (LanConst > 5. * mpv)
		LanConst = Lan.Landau(mpv, DA);
	LanConst *=1000/loss_energy*d;

	pos->Scale(1/total_pairs*LanConst);
	neg->Scale(1/total_pairs*LanConst);
	
	sum->Add(neg);
	sum->Add(pos);

	delete histop;
	delete histon;
}

void KDetector::Drift(Double_t sx, Double_t sy, Double_t sz, Float_t charg, KStruct *seg, Double_t t0)
{
	//Drift simulation for a point charge (Float_t charg;)
	//starting from ( sx,sy, sz)
	//KStruct *seg is the structure where the  the drift paths, drift times and induced cahrges are stored

	Double_t Stime=0;                       // Step time
	Double_t difx=0,dify=0,difz=0;          // diffusion steps in all directions
	Double_t sigma;                         // sigma of the diffusion step  
	Float_t *xcv,*ycv,*zcv,*time,*charge;   // drift variables to be stored in KStruct
	Double_t cx=0,cy=0,cz=0;                // current position of the charge bucket
	Double_t vel=0;                         // drift velocity
	Double_t vth2=0;                        // thermal velocity sqared
	Double_t tfc;                           // trapped charge fraction   
	Double_t t=0;                           // drift time
	Double_t sumc=0;                        // total induced charge
	Double_t ncx=0,ncy=0,ncz=0;             // next position of the charge bucket   
	Double_t deltacx,deltacy,deltacz;       // drift step due to drift

	Int_t st=0;                             // current step
	Int_t ishit=0;                          // local counter of the step
	TVector3 *EE;                           // Set up electric field vector
	TVector3 *EEN;                          // Set up electric field vector
	TVector3 FF;                            // Combined drift field
	Float_t pathlen=0;                      // pathlength
	Float_t WPot;                           // current ramo potential

	// Inclusion of Magnetic field 28.8.2001 - revised 15.10.2012
	TVector3 BB(B);                           // Create a magnetic field vector
	Float_t muhe=1650;                        // parametrization of the magnetic field 
	Float_t muhh=310;	                  // is based on simply this mobility parameters

	// Start time in the  absolute domain (used for delayed charge generation in multiplication 

	t=t0;

	// Intitialize KStruct class and its members 

	xcv=seg->Xtrack; 
	ycv=seg->Ytrack;  
	zcv=seg->Ztrack; 
	time=seg->Time; 
	charge=seg->Charge; 
	seg->Clear();
	seg->PCharge=(Int_t) charg;

	// start drift

	cx=sx; cy=sy; cz=sz;                   // set current coordinates
	xcv[st]=cx; ycv[st]=cy; zcv[st]=cz;    // put the first point in the KStruct 
	time[st]=t;  charge[st]=0; 

	EE=Real->CalFieldXYZ(cx,cy,cz);         // Get the electric field vector 
	EEN=Real->CalFieldXYZ(cx,cy,cz);       // Get the electric field vector for next step - here the default is the same 12.9.2018
	seg->Efield[st]=EE->Mag();             // Store the magnitude of E field






	//printf("%d %f : (%f %f %f)\n",st,charg,cx,cy,cz);

	do 
	{
		//    printf("Calculate field\n");
		st++;
		if(charg>0)
			FF=(*EE)+muhh*EE->Cross(BB); else 
				FF=(*EE)-muhe*EE->Cross(BB); 
		// "-muhe" stands for the fact that at the same 
		// field the drift direction has changed due to different charge

		//printf("Field : %f %f %f (%f %f %f)  ---- ",FF[0],FF[1],FF[2],(*EE)[0],(*EE)[1],(*EE)[2]);
		if(FF.Mag()!=0)
		{
			deltacy=-SStep*charg*FF.y()/FF.Mag();   // deltay 
			deltacx=-SStep*charg*FF.x()/FF.Mag();   // deltax
			deltacz=-SStep*charg*FF.z()/FF.Mag();   // deltaz
		}
		else { deltacx=0; deltacy=0; deltacz=0; }

		//    printf("Calculate velocity \n");

		if(DM!=NULL)
			KMaterial::Mat=DM->GetBinContent(DM->FindBin(cx,cy,cz)); 
		else KMaterial::Mat=0;
						 //      EEN=Real->CalFieldXYZ(cx+deltacx,cy+deltacy,cz+deltacz); // get field & velocity at new location //12.9.2018
		Real->CalFieldXYZ(cx + deltacx, cy + deltacy, cz + deltacz, EEN); // get field & velocity at new location
		vel=Real->DriftVelocity( (EEN->Mag()+EE->Mag())/2,charg,Temperature,TMath::Abs(NeffF->Eval(cx,cy,cz)),MobMod());

		//printf("Calculate vel: %e EEN = %e ::: ",vel, EEN->Mag());
		if(vel==0) {
			deltacx=0; deltacy=0; deltacz=0; 
			difx=0; dify=0; difz=0;
			ishit=9;
		} 
		else 
			if(diff && !(DiffOffField<EEN->Mag() && MTresh>1))  
				// is diffusion ON - if yes then include it
				// if multiplication is ON the diffusion must be switched 
				// off when the field gets large enough
			{
				Stime=SStep*1e-4/vel; // calcualte step time  
				sigma=TMath::Sqrt(2*Kboltz*Real->Mobility(EE->Mag(),Temperature,charg,TMath::Abs(NeffF->Eval(cx,cy,cz)),MobMod())*Temperature*Stime); 
				dify=ran->Gaus(0,sigma)*1e4; 
				difx=ran->Gaus(0,sigma)*1e4;
				if(nz!=1) difz=ran->Gaus(0,sigma)*1e4; else difz=0;
			} else {difx=0; dify=0; difz=0;}

		if((cx+deltacx+difx)>=GetUpEdge(0)) ncx=GetUpEdge(0); else
			if((cx+deltacx+difx)<GetLowEdge(0)) ncx=GetLowEdge(0); else
				ncx=cx+(deltacx+difx);;

		if((cy+deltacy+dify)>=GetUpEdge(1)) ncy=GetUpEdge(1); else
			if((cy+deltacy+dify)<GetLowEdge(1)) ncy=GetLowEdge(1); else
				ncy=cy+(deltacy+dify);

		if((cz+deltacz+difz)>=GetUpEdge(2)) ncz=GetUpEdge(2); else
			if((cz+deltacz+difz)<GetLowEdge(2)) ncz=GetLowEdge(2); else
				ncz=cz+(deltacz+difz);

		if(Debug) printf("%d %f E=%e (%e %e %e): x:%f->%f y:%f->%f z:%f->%f (%f %f %f)(%f %f %f) : Mat=%d :: ",st,charg,EEN->Mag(),EEN->x(),EEN->y(),EEN->z(),cx,ncx,cy,ncy,cz,ncz,deltacx,deltacy,deltacz,dify,dify,difz,KMaterial::Mat);

		charge[st]=charg*(Ramo->CalPotXYZ(ncx,ncy,ncz)-Ramo->CalPotXYZ(cx,cy,cz));
		cx=ncx; cy=ncy; cz=ncz;

		//////////////////// calculate strict trapping, e.g. depending on position ///////////////////////
		if(TauE!=NULL && TauH!=NULL)
		{
			if(charg<0) {
				// vth2=3*Kboltz*Temperature*Clight*Clight/(511e3*EmeC(KMaterial::Mat))*1e4*0;
				tfc=1e4*(TauE->Eval((ncx+cx)/2,(ncy+cy)/2,(ncz+cz)/2)*TMath::Sqrt(vel*vel+vth2));	              
			}
			else
			{
				//   vth2=3*Kboltz*Temperature*Clight*Clight/(511e3*EmhC(KMaterial::Mat))*1e4*0;
				tfc=1e4*(TauH->Eval((ncx+cx)/2,(ncy+cy)/2,(ncz+cz)/2)*TMath::Sqrt(vel*vel+vth2));
			}	    
			if(ran->Rndm()>TMath::Exp(-SStep/tfc)) ishit=12;
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////

		sumc+=charge[st];


		if(vel!=0) 
		{
			t=t+SStep*1e-4/vel; //else t+=Stime;
			pathlen+=SStep;
		}

		xcv[st]=cx;
		ycv[st]=cy;
		zcv[st]=cz;
		time[st]=t;

		//    EE=Real->CalFieldXYZ(cx,cy,cz);     //12.9.2018
		Real->CalFieldXYZ(cx,cy,cz,EE);     //12.9.2018
		seg->Efield[st]=EE->Mag();

		// Checking for termination of the drift //// 

		WPot=Ramo->CalPotXYZ(cx,cy,cz);
		if(WPot>(1-Deps)) ishit=1;
		// if(TMath::Abs(WPot)<Deps) ishit=2;
		if(cx<= GetLowEdge(0)) ishit=3;
		if(cx>=  GetUpEdge(0)) ishit=4;
		if(cy<= GetLowEdge(1)) ishit=5;
		if(cy>= GetUpEdge(1))  ishit=6;
		if(cz<= GetLowEdge(2)) ishit=7;
		if(cz>= GetUpEdge(2))  ishit=8;
		if(pathlen>MaxDriftLen) ishit=11;
		if(st>=MAXPOINT-1) ishit=20;   

		if(Debug) printf("(t=%e, vel=%e, velth=%e) [Ch=%f ChInt=%f TFC=%f] Ishit=%d \n",t,vel,TMath::Sqrt(vth2),charge[st],sumc,tfc,ishit);

	} while (!ishit); // Do until the end of drift


	(*seg).Xlenght=pathlen; (*seg).Ylenght=pathlen; 
	(*seg).TTime=t; (*seg).TCharge=sumc; (*seg).Steps=st;
	delete EE; delete EEN;

	return;
}





// K3D

#include "Rtypes.h"
#include "TArray.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TGraph.h"
#include <stdio.h>
#include <stdlib.h>
#include "TMinuit.h"



class K3D : public KDetector
{
	private:
	public:
		Int_t Col;
		Float_t CellZ;
		Float_t CellX;
		Float_t CellY;

		Float_t *PosD; //[Col]
		Float_t *PosX; //[Col]
		Float_t *PosY; //[Col]
		Float_t *PosR; //[Col]
		Short_t *PosW; //[Col]
		Short_t *PosM; //[Col]

		K3D(Int_t, Float_t = 100, Float_t = 100, Float_t = 105);
		~K3D();
		void SetUpVolume(Float_t, Float_t);
		void SetUpColumn(Int_t, Float_t, Float_t, Float_t, Float_t, Short_t, Short_t);
		void SetUpElectrodes(Int_t = 0);

		// ClassDef(K3D, 1)
};

// ClassImp(K3D)
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// K3D                                                                  //
//                                                                      //
// Description of the 3D detector.                                      //
// The geometry of the                                                  //
// is defined by                                                        //
// dasdasd                                                              //
//////////////////////////////////////////////////////////////////////////

K3D::~K3D()
{
	delete PosD;
	delete PosX;
	delete PosY;
	delete PosR;
	delete PosW;
}

K3D::K3D(Int_t x1, Float_t x, Float_t y, Float_t z)
{
	Col = x1;
	PosD = new Float_t[Col];
	PosX = new Float_t[Col];
	PosY = new Float_t[Col];
	PosR = new Float_t[Col];
	PosW = new Short_t[Col];
	PosM = new Short_t[Col];

	for (Int_t i = 0; i < Col; i++) {
		PosD[i] = 0;
		PosX[i] = 0;
		PosY[i] = 0;
		PosR[i] = 0;
		PosW[i] = 0;
		PosM[i] = 0;
	}

	CellZ = z;
	CellX = x;
	CellY = y;
}

void K3D::SetUpVolume(Float_t St1, Float_t St2)
{
	nx = (int)(CellX / St1);
	ny = (int)(CellY / St1);
	nz = (int)(CellZ / St2);

	EG = new TH3I("EG", "EG", nx, 0, CellX, ny, 0, CellY, nz, 0, CellZ);
	EG->GetXaxis()->SetTitle("x [#mum]");
	EG->GetYaxis()->SetTitle("y [#mum]");
	EG->GetZaxis()->SetTitle("z [#mum]");

	GetGrid(EG, 1);
	//printf("nz：%i\n", nz);
}

void K3D::SetUpColumn(Int_t n, Float_t posX, Float_t posY, Float_t R, Float_t Depth, Short_t Wei, Short_t Mat)
{
	PosD[n] = Depth;
	PosR[n] = R;
	PosX[n] = posX;
	PosY[n] = posY;
	PosW[n] = Wei;
	PosM[n] = Mat;
}

void K3D::SetUpElectrodes(Int_t back)
{
	Float_t Pos[3], L = 0;
	Int_t i, j, k;
	if (back) {
		for (k = 1; k <= nz; k++)
			for (j = 1; j <= ny; j++)
				for (i = 1; i <= nx; i++)
					if (k == 1)
						EG->SetBinContent(i, j, k, 2);
					else
						EG->SetBinContent(i, j, k, 0);
	}

	for (Int_t i = 0; i < Col; i++) {
		Pos[0] = PosX[i];
		Pos[1] = PosY[i];
		Pos[2] = PosD[i] >= 0 ? Pos[2] = PosD[i] / 2. : (2 * CellZ + PosD[i]) / 2.;
		L = TMath::Abs(PosD[i] / 2.);
		ElCylinder(Pos, PosR[i], L, 3, PosW[i], PosM[i]);
	}
}


// KPad
#define JMAX 40
class KPad : public KDetector
{
private:
//Runge Kutta method for solving the field
  void           rk4(float *,float *,int,float,float,float*); 
  Float_t        rtbis(float, float, float);
  Float_t        PoEqSolve(Float_t);
  void           Derivs(float x,float *,float *);
  TArrayF PhyPot;       //electric potential
  TArrayF PhyField;     //electric field 
public:
   TF1     *Neff;   // effective dopping concentration 
   Float_t CellY;   // thickness of the diode
   Float_t CellX;   // width of the diode

   KPad(Float_t=50,Float_t=301);
  ~KPad(); 
   	void SetUpVolume(Float_t,Int_t =0); 
	void SetUpElectrodes();
	void GetRamoField();
	void GetField();
	TGraph   *DrawPad(char*);
};

KPad::KPad(Float_t x,Float_t y)
{
  	//Constructor of the class KPad:
  	//		Float_t x ; width of the diode in um. For the simulation it is not neccessary to put real dimensions
  	//		Float_t y ; thickness in um
 	Neff=NULL;
  	CellX=x;  CellY=y;
}
 
KPad::~KPad()
{
}

void KPad::Derivs(float x, float y[], float dydx[])
{
	//Double_t permMat=(material==0)?perm:permDi;
	Double_t permMat = Perm(KMaterial::Mat);
	Double_t permV = perm0 * 1e-6;
	;
	dydx[1] = y[2];
	dydx[2] = -Neff->Eval(x) * e_0 / (permMat * permV);
	// if(x>150) dydx[2]*=permMat;
}

Float_t KPad::PoEqSolve(Float_t der)
{
	//TArrayF PhyField(ny+2);
	//TArrayF PhyPot(ny+2);
	//Float_t Step=1;

	Int_t i;
	Float_t h, x = 0;

	Float_t y[3], dydx[3], yout[3];

	y[1] = Voltage;
	y[2] = der;
	PhyField[0] = y[2];
	PhyPot[0] = y[1];
	
	Derivs(x, y, dydx);
	
	for (i = 1; i <= ny; i++)
	{
		h = GetStepSize(1, i);
		rk4(y, dydx, 2, x, h, yout);
		//printf("%f %f %f\n",x+h,yout[1],yout[2]);
		y[1] = yout[1]; //electric potential
		y[2] = yout[2];	// electric filed
		PhyField[i] = y[2];
		PhyPot[i] = y[1];
		//std::cout << "x:" << x << ",y[x]:" << y[1] << "," << y[2] << std::endl;
		x = x + h;
		Derivs(x, y, dydx);
		
	}
	//     printf("y[1]=%f\n",xp1);
	return y[1];
}

void KPad::SetUpVolume(Float_t St1, Int_t Mat)
{
  	nx=(Int_t) (CellX/St1);
  	ny=(Int_t) (CellY/St1);
  	nz=1;
  	//Set the boundary condition matrix //
  	EG=new TH3I("EG","EG",nx,0,CellX,ny,0,CellY,1,0,1);
  	EG->GetXaxis()->SetTitle("x [#mum]");
  	EG->GetYaxis()->SetTitle("y [#mum]");

  	DM=new TH3I("DM","DM",nx,0,CellX,ny,0,CellY,1,0,1);
  	DM->GetXaxis()->SetTitle("x [#mum]");
  	DM->GetYaxis()->SetTitle("y [#mum]");
  	for(int i=0;i<nx;i++) 
     	for(int j=0;j<ny;j++)
        	DM->SetBinContent(i,j,1,Mat);
}

void KPad::SetUpElectrodes()
{
 
 	for(int i=1;i<=nx;i++){ EG->SetBinContent(i,1,1,2);  EG->SetBinContent(i,ny,1,16385);} 
 	// KMaterial::Mat=10;
  	//Default track
 	// enp[0]=CellX/2;  exp[0]=enp[0];
 	// enp[1]=1;        exp[1]=CellY;

 	if(Neff!=NULL )
 	{
		GetField();
		GetRamoField();
 	}
 	else 
   		printf("Please define space charge function Neff before field calculation\n");
 
}

Float_t KPad::rtbis(float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j;
	float dx, f, fmid, xmid, rtb;

	f = PoEqSolve(x1);
	fmid = PoEqSolve(x2);
	if (f * fmid >= 0.0)
		nrerror("Root must be bracketed for bisection in rtbis");
	rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
	for (j = 1; j <= JMAX; j++)
	{
		//xmid find a good electric field intial value
		fmid = PoEqSolve(xmid = rtb + (dx *= 0.5));
		// fmid is the electric field at the readout electrodes
		if (fmid <= 0.0)
			rtb = xmid;
		if (fabs(dx) < xacc || fmid == 0.0)
			return rtb;
	}
	nrerror("Too many bisections in rtbis");
	return 0.0;
}
void KPad::rk4(float y[], float dydx[], int n, float x, float h, float yout[])
{
	float *vector(long, long);
	void free_vector(float *, long, long);
	int i;
	float xh, hh, h6, *dym, *dyt, *yt;

	dym = vector(1, n);
	dyt = vector(1, n);
	yt = vector(1, n);
	hh = h * 0.5;
	h6 = h / 6.0;
	xh = x + hh;
	for (i = 1; i <= n; i++)
		yt[i] = y[i] + hh * dydx[i];
	Derivs(xh, yt, dyt);
	for (i = 1; i <= n; i++)
		yt[i] = y[i] + hh * dyt[i];
	Derivs(xh, yt, dym);
	for (i = 1; i <= n; i++)
	{
		yt[i] = y[i] + h * dym[i];
		dym[i] += dyt[i];
	}
	Derivs(x + h, yt, dyt);
	for (i = 1; i <= n; i++)
		yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
	free_vector(yt, 1, n);
	free_vector(dyt, 1, n);
	free_vector(dym, 1, n);
}

void KPad::GetRamoField()
{
	Int_t i, j;
	Double_t *x = new Double_t[nx * ny + 1];
	for (j = 1; j <= ny; j++)
	{
		for (i = 1; i <= nx; i++)
		{
			x[i + (j - 1) * nx] = (Double_t)(j - 1) / (Double_t)(ny - 1);
			
		}
	}
	Ramo->U = MapToGeometry(x);
	Ramo->CalField();
	delete[] x;
}

void KPad::GetField()
{
	Int_t i, j;
	//Float_t Step=1;
	Float_t aa;
	PhyPot = TArrayF(ny + 2);
	PhyField = TArrayF(ny + 2);
	Double_t *x = new Double_t[nx * ny + 1];
	//TArrayD PhyPot2D(nx*ny+1);
	//TArrayI StripPosition=TArrayI(2); StripPosition[0]=1; StripPosition[1]=nx;

	aa = rtbis(-100, 100, 0.000001);
	//if(!Invert && GetDepletionVoltage()>Voltage) aa=rtbis(1,300,0.000001); else aa=rtbis(-10,10,0.000001);

	for (j = 1; j <= ny; j++)
		for (i = 1; i <= nx; i++)
		{
			x[i + (j - 1) * nx] = (Double_t)PhyPot[j - 1];
		}
	//Real=EField(PhyPot2D,nx,ny);
	//Real->CalField(Step,1);

	// for(i=0;i<nx*ny+1;i++)     {printf("%f ",PhyPot2D[i]); if(i%nx==0) printf("\n");}

	Real->U = MapToGeometry(x);
	Real->CalField();
	delete[] x;
}

TGraph *KPad::DrawPad(char *option)
{
  //Draws potential "p" or  electric field "f"
	Char_t Opt[10];
	TGraph *gr;
	TArrayF xx=TArrayF(ny+1);
 	for(Int_t i=0;i<=ny;i++) xx[i]=(Float_t)i*GetStepSize(1,i);

	if(strchr(option,'p')!=NULL) gr=new TGraph(ny+1,xx.GetArray(),PhyPot.GetArray());
	if(strchr(option,'f')!=NULL) gr=new TGraph(ny+1,xx.GetArray(),PhyField.GetArray());
	if(strchr(option,'s')!=NULL) sprintf(Opt,"CP"); else  sprintf(Opt,"ACP"); 
	//if(! strcmp(option,"p")) gr=new TGraph(ny,xx.GetArray(),PhyPot.GetArray());
	//if(! strcmp(option,"f")) gr=new TGraph(ny,xx.GetArray(),PhyField.GetArray());
	gr->Draw(Opt);
	gr->GetHistogram()->Draw();
	gr->Draw(Opt);
	return gr;
}


void set_root_style(int stat=1110, int grid=0){
	gROOT->Reset();

	gStyle->SetTitleFillColor(0) ; 
	gStyle->SetTitleBorderSize(0); 

	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetCanvasDefX(0); 
	gStyle->SetCanvasDefY(0); 
	gStyle->SetFrameBorderMode(0); 
	gStyle->SetFrameBorderSize(1); 
	gStyle->SetFrameFillColor(0); 
	gStyle->SetFrameFillStyle(0); 
	gStyle->SetFrameLineColor(1); 
	gStyle->SetFrameLineStyle(1); 
	gStyle->SetFrameLineWidth(1); 

	// gStyle->SetPadTopMargin(PadTopMargin);  
	gStyle->SetPadLeftMargin(0.10);  
	gStyle->SetPadRightMargin(0.05);  

	gStyle->SetLabelSize(0.03, "XYZ");  
	gStyle->SetTitleSize(0.04, "XYZ");  
	gStyle->SetTitleOffset(1.2, "Y");  

	gStyle->SetPadBorderMode(0);  
	gStyle->SetPadColor(0);  
	gStyle->SetPadTickX(1); 
	gStyle->SetPadTickY(1); 
	gStyle->SetPadGridX(grid); 
	gStyle->SetPadGridY(grid); 

	gStyle->SetOptStat(stat); 
	gStyle->SetStatColor(0); 
	gStyle->SetStatBorderSize(1); 
}


TGraph * get_graph_from_log(TString inputFile, TString& err_msg) {
	std::ifstream in;
	in.open(inputFile); 
	Float_t x, y; 
	Float_t factor_I(1.0);
	if (inputFile.Contains("_uA_") ) {
		factor_I = 1e-6;
		std::cout << "Using current factor: " << factor_I << " for " << inputFile << std::endl;  
	} 
	std::vector<Float_t> voltages; 
	std::vector<Float_t> currents; 

	std::string line;
	int nlines = 0;

	// float i_v150, i_v100;
	double i_v150, i_v100;
	bool pass1(false), pass2(false); 

	while (getline(in, line)) {
		std::istringstream iss(line);
		if ( line.find("#") == 0 ) continue; 
		if (!(iss >> x >> y )) break; 
		if (!in.good()) break;
		voltages.push_back(fabs(x));
		currents.push_back(fabs(y)*factor_I);
		// Pick up values like: -100.043	
		if ( fabs(fabs(x)-150) < 1) i_v150 = fabs(y); 
		if ( fabs(fabs(x)-100) < 1) i_v100 = fabs(y); 
		nlines ++; 
	}

	if (nlines < 1) {
		std::cerr << "No valid data found in : " << inputFile << std::endl;
		return NULL; 
	}

	if ( i_v150 < 2E-6) pass1 = true;
	if ( i_v150/i_v100 < 2 ) pass2 = true;

	if (!pass1) err_msg = Form("I(150V) >= 2uA (%.1e)",i_v150) ;
	if (!pass2) err_msg += Form("I(150V)/I(100V) >= 2 (%.1f)", i_v150/i_v100) ;

	in.close();
	TGraph *gr = new TGraph(nlines, &voltages[0], &currents[0]);
	return gr; 
}

TCanvas* drawIV(std::vector<TString> inputFiles){
	set_root_style();

	TCanvas *c = new TCanvas("c", "IV scan", 800, 800);
	c->SetGrid();

	TMultiGraph *mg = new TMultiGraph();
	TLegend *leg = new TLegend(0.2, 0.6, 0.5, 0.8);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->SetNColumns(1);
	leg->SetTextSize(0.02);
	leg->SetTextSizePixels(25);

	for (std::vector<int>:: size_type i = 0; i != inputFiles.size(); i++) {
		TString err_msg = "";  
		TGraph *gr = get_graph_from_log(inputFiles[i], err_msg);
		if (gr == NULL) continue; 
		gr->SetMarkerStyle(20+i);
		gr->SetMarkerSize(0.9);
		int color = i+1;
		if (color >= 5) color ++; // bypass the yellow  
		if (color >= 10) color = color % 10 + 1 ; // reuse the first 9 colors
		gr->SetMarkerColor(color);
		leg->AddEntry(gr, Form("%s %s", inputFiles[i].Data(),
					err_msg.Data()), "p"); 
		mg->Add(gr); 
	}

	mg->Draw("APL"); 
	mg->GetXaxis()->SetTitle("Bias Voltage [V]");
	mg->GetYaxis()->SetTitle("Leakage Current [A]");
	mg->GetYaxis()->SetRangeUser(1e-10, 1e-4); 
	leg->Draw(); 

	c->SetLogy();
	c->Update(); 
	return c;
}
#define EXPORT

class EXPORT KElec
{
private:
	Int_t Method;
	Double_t Cp;
	Double_t Rp;
	Double_t Crc;
	Double_t R1rc;
	Double_t R2rc;
	Double_t Ccr;
	Double_t R1cr;
	Double_t R2cr;
	Double_t PeakTime;
	Double_t IntTime;

public:
	KElec(Double_t = 5e-12, Double_t = 50, Double_t = 25e-9, Double_t = 1e-12, Double_t = 200, Double_t = 200, Double_t = 1e-12, Double_t = 200, Double_t = 200, Double_t = 25e-9, Int_t = 0);
	virtual ~KElec();
	Double_t Trapez(TH1F *, Int_t, Double_t);
	Double_t Simpson(TH1F *, Int_t, Double_t);
	void preamp(TH1F *his) { preamp(Cp, Rp, his, IntTime, Method); };
	void preamp(Double_t, Double_t, TH1F *, Double_t = -1111, Int_t = 0);
	void RCshape(Double_t, Double_t, Double_t, TH1F *, Int_t = 0);
	void CRshape(Double_t, Double_t, Double_t, TH1F *, Int_t = 0);
	Double_t CSAamp(TH1F *his, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
};



KElec::KElec(Double_t x1, Double_t x2, Double_t x3, Double_t x4, Double_t x5, Double_t x6, Double_t x7, Double_t x8, Double_t x9, Double_t x10, Int_t met)
{
	Cp = x1;
	Rp = x2;
	IntTime = x3;
	Crc = x4;
	R1rc = x5;
	R2rc = x6;
	Ccr = x7;
	R1cr = x8;
	R2cr = x9;
	PeakTime = x10;
	Method = met;
}
KElec::~KElec()
{
	//Clear();
}

Double_t KElec::CSAamp(TH1F *histo, Double_t TRise, Double_t TFall, Double_t C_detector, Double_t CSATransImp, Double_t CSA_Noise, Double_t CSAVth)
{
	//test graph//
	TH1F *whisto = new TH1F();
	histo->Copy(*whisto);
	//test graph //

	int IMaxSh=histo->GetNbinsX();
	int Step=1;
	int DStep = Step;
	float Qdif_Shaper = 0;
	double *itot = new double[IMaxSh];
	double *PreAmp_Q = new double[IMaxSh];
	double *ShaperOut_Q = new double[IMaxSh];
	double *ShaperOut_V = new double[IMaxSh];
	double TauRise = TRise / 2.2*1e-9;
	double TauFall = TFall / 2.2*1e-9;
	double sh_max = 0.0;
	double sh_min = 0.0;
	int NMax_sh = 0;
	int Nmin_sh = 0;
	double thre_time=0.0;

	double qtot=0;
	int SC_CSAOn=1;



	if(TauRise == TauFall) TauRise *= 0.9;


	double TIMEUNIT=(histo->GetXaxis()->GetXmax() - histo->GetXaxis()->GetXmin()) / histo->GetNbinsX();
	

	for (int k=0;k<IMaxSh;k++)
 	{
		PreAmp_Q[k]=0.0;
		itot[k]=0.0;
		ShaperOut_Q[k]=0.0;
		ShaperOut_V[k]=0.0;
	}

	//get the induced charge anc current
	for (int i = 0; i < IMaxSh - Step; i += Step)
	{
		if (i > 0)
		{
			PreAmp_Q[i] = 0;
			itot[i]=histo->GetBinContent(i);
			PreAmp_Q[i] = itot[i] * TIMEUNIT;
			qtot += itot[i] * TIMEUNIT;
			//histo->SetBinContent(i,PreAmp_Q[i]);   // induced charge each step
		}
		else if (i == 0)
		{
			PreAmp_Q[i] = 0; // final scale to zero
			//histo->SetBinContent(i, PreAmp_Q[i]);
		}
		// Iout_temp reproduces correctly Speiler pag 127//
	}

	// get the CSA charge.   transfer induced charge to CAS charge
	for (int i = 0; i < IMaxSh - Step; i += Step)
	{
		Qdif_Shaper = PreAmp_Q[i];

		if (Qdif_Shaper !=0)
			for (int ll = 0; ll<IMaxSh-i;ll+=Step)  // valid only up to IMaxSh 
			{
				ShaperOut_Q[i+ll] += TauFall/(TauFall+TauRise)*Qdif_Shaper*		    
				(TMath::Exp(-ll*TIMEUNIT/TauFall)-TMath::Exp(-ll*TIMEUNIT/TauRise));					 			
			}
		if (fabs(ShaperOut_Q[i]) > fabs(sh_max))
		{
			sh_max = ShaperOut_Q[i];
		}
		
	}


	//TOFFEE_gain  change charge to amplitude ??
	double Ci = 500 * 70*1E-15;
	double Qfrac = 1. / (1. + C_detector * 1E-12 / Ci);
	for (int i = 0; i < IMaxSh; i += Step)
	{
		if (SC_CSAOn)
			ShaperOut_V[i] = ShaperOut_Q[i] * CSATransImp * 1e15 * qtot / sh_max; // [mV/fQ *Q/Q]
		else
			ShaperOut_V[i] = ShaperOut_Q[i] * CSATransImp * 1e15 * qtot * Qfrac / sh_max; // [mV/fQ *Q/Q]											  //cout << ShaperOut_V[i] << endl;
		histo->SetBinContent(i, fabs(ShaperOut_V[i]));
	}

	sh_max = 0.;

	for (int i=Step;i<IMaxSh-2*Step;i+=Step)
	{	
		if (ShaperOut_V[i]> sh_max) 
		{
			sh_max = ShaperOut_V[i];		    
			NMax_sh = i;
		}
	
		if ( ShaperOut_V[i] < sh_min) 
		{
			sh_min = ShaperOut_V[i];		    
			Nmin_sh = i;
		}
	}
	// // CSA noise

	int NCSA_der0 = 0;
	double CSAx1 = 0;
	double CSAx2 = 0;
	double Jitter = 0;
	bool FVTh = true;
	float DT = (TIMEUNIT*1e9*2.*DStep);
	double STime = 0;

	NCSA_der0 = (fabs(sh_max) > fabs(sh_min)) ? NMax_sh : Nmin_sh;

	TRandom3 r;
	r.SetSeed(0);
	CSAx1 = r.Gaus(CSA_Noise,2.36);
	
	//histo->Reset();

	for (int i = 2 * DStep; i < IMaxSh - 2 * DStep; i++)
	{
		if ((fabs(ShaperOut_V[i]) + CSAx1) > fabs(CSAVth) && FVTh && i < NCSA_der0 - 20)
		{
			Jitter = 0;
			FVTh = false;
			STime = (double)(i)*TIMEUNIT * 1e9;
			float dVdt = fabs(ShaperOut_V[i + DStep] - ShaperOut_V[i - DStep]) / DT; //mV/nSec											 //mV /nSec
			float Jit = 0;
			if (dVdt != 0)
			{
				Jitter = CSA_Noise / dVdt; // in ns
				Jit = gRandom->Gaus(0, Jitter);
			}
			thre_time=STime+Jit;
		}
	}
	//record the waveform above the threshold
	// if (!FVTh)
	// {
	// 	for (int i = 2 * DStep; i < IMaxSh - 2 * DStep; i++)
	// 	{
	// 		CSAx2 = r.Gaus(CSA_Noise, 2.36);
	// 		//histo->SetBinContent(i,fabs(ShaperOut_V[i]) + CSAx2);
	// 	}
	// }

	delete whisto;
	return thre_time;
	}

void KElec::preamp(Double_t C, Double_t R, TH1F *histo, Double_t cut, Int_t method)
{
	//static double e_0 = 1.60217733e-19;
	Double_t suma = 0;
	Double_t tau = R * C;
	Double_t t;
	Int_t i;
	TH1F *whisto = new TH1F();
	histo->Copy(*whisto);

	if (cut == -1111)
		cut = histo->GetNbinsX() * histo->GetBinWidth(1);

	for (i = 1; i < histo->GetNbinsX() - 1; i++)
	{
		t = whisto->GetBinCenter(i);
		if (t <= cut)
		{
			if (method == 0)
				suma += Trapez(whisto, i, tau);
			if (method == 1)
				suma += Simpson(whisto, i, tau);
		}
		histo->SetBinContent(i, (Float_t)(1 / C * suma * TMath::Exp(-t / tau)));
		// printf("%d,t=%e,suma=%e,tau=%e,tapez=%e\n",i,t,suma,tau,Trapez(histo,i,tau));
	}
	delete whisto;
}

Double_t KElec::Trapez(TH1F *histo, Int_t i, Double_t tau)
{
	Double_t f1 = 0, f2 = 0, h = 0, t1 = 0, t2 = 0;
	if (i == 1)
		return (histo->GetBinContent(i) * histo->GetBinWidth(i) / 2);
	else
	{
		t1 = histo->GetBinCenter(i - 1);
		f1 = histo->GetBinContent(i - 1) * TMath::Exp(t1 / tau);
		t2 = histo->GetBinCenter(i);
		f2 = histo->GetBinContent(i) * TMath::Exp(t2 / tau);
		h = histo->GetBinWidth(i); //printf("tutut %d,t1=%e,t2=%e,f1=%e,f2=%e\n",i,t1,t2,f1,f2);
		return (h * 0.5 * (f1 + f2));
	}
}

Double_t KElec::Simpson(TH1F *histo, Int_t i, Double_t tau)
{
	Double_t f1 = 0, f2 = 0, f3 = 0, h = 0, t1 = 0, t2 = 0, t3 = 0;
	if (i < 3 || i % 2 == 0)
		return (Trapez(histo, i, tau));
	else
	{
		h = histo->GetBinWidth(i);
		t1 = histo->GetBinCenter(i - 2);
		f1 = histo->GetBinContent(i - 2) * TMath::Exp(t1 / tau);
		t2 = histo->GetBinCenter(i - 1);
		f2 = histo->GetBinContent(i - 1) * TMath::Exp(t2 / tau);
		t3 = histo->GetBinCenter(i);
		f3 = histo->GetBinContent(i) * TMath::Exp(t3 / tau);
		return (h * (0.33333333333 * f1 + 1.333333333 * f2 + 0.3333333333333 * f3) - Trapez(histo, i - 1, tau));
	}
}

void KElec::RCshape(Double_t C, Double_t R1, Double_t R2, TH1F *histo, Int_t method)
{
	Double_t suma = 0;
	Double_t tau = (R1 * R2) * C / (R1 + R2);
	Double_t tau1 = R1 * C;
	Double_t t;
	Int_t i;
	TH1F *whisto = new TH1F();
	histo->Copy(*whisto);

	for (i = 1; i < histo->GetNbinsX() - 1; i++)
	{
		t = whisto->GetBinCenter(i);
		if (method == 0)
			suma += Trapez(whisto, i, tau);
		if (method == 1)
			suma += Simpson(whisto, i, tau);
		histo->SetBinContent(i, (Float_t)(1 / tau1 * suma * TMath::Exp(-t / tau)));
	}
	delete whisto;
}

void KElec::CRshape(Double_t C, Double_t R1, Double_t R2, TH1F *histo, Int_t method)
{
	Double_t suma = 0;
	Double_t tau = (R1 * R2) * C / (R1 + R2);
	Double_t tau1 = R1 * C;
	Double_t t, f;
	Int_t i;
	TH1F *whisto = new TH1F();
	histo->Copy(*whisto);

	for (i = 1; i < histo->GetNbinsX() - 1; i++)
	{
		t = whisto->GetBinCenter(i);
		f = whisto->GetBinContent(i);
		if (method == 0)
			suma += Trapez(whisto, i, tau);
		if (method == 1)
			suma += Simpson(whisto, i, tau);
		histo->SetBinContent(i, f - (1 / tau - 1 / tau1) * suma * TMath::Exp(-t / tau));
	}
	delete whisto;
}

#ifndef __CINT__ 

void print_usage(){
	printf("NAME\n\traser - REdiation SEmiconductoR\n");
	printf("\nSYNOPSIS\n\traser arg ...\n");
	printf("\nOPTIONS\n");
	printf("\t%-5s  %-40s\n", "-h", "Print this message");
	printf("\t%-5s  %-40s\n", "3D", "3D KDetSim");
	printf("\t%-5s  %-40s\n", "2D", "2D KDetSim");
	printf("\nAUTHOR\n\tXin Shi <shixin@ihep.ac.cn>\n");
}


int main(int argc, char** argv) {

	if (argc < 2)
	{
		print_usage();
		return -1;
	}

	TApplication theApp("App", 0, 0);
	theApp.SetReturnFromRun(true);
	gStyle->SetCanvasPreferGL(kTRUE);

	for (int i = 0; i < argc; i++)
	{
		if (!strcmp(argv[i], "-h"))
		{
			print_usage();
			break;
		}

		if (!strcmp(argv[i], "3D"))
		{
			// define a 3D detector with 5 electrodes
			// x=100 , y is 50 and thickness 120
			K3D *det = new K3D(7, 80, 80, 300);
			det->Voltage = 50;
			// define the drift mesh size and simulation mesh size in microns
			det->SetUpVolume(0.1, 1);
			// define  columns #, postions, weigthing factor 2=0 , material Al=1
			det->SetUpColumn(0, 40, 15, 4, 280, 2, 1);
			det->SetUpColumn(1, 40, 65, 4, 280, 2, 1);
			det->SetUpColumn(2, 61.65, 27.5, 4, 280, 2, 1);
			det->SetUpColumn(3, 61.65, 52.5, 4, 280, 2, 1);
			det->SetUpColumn(4, 18.35, 27.5, 4, 280, 2, 1);
			det->SetUpColumn(5, 18.35, 52.5, 4, 280, 2, 1);
			det->SetUpColumn(6, 40, 40, 4, -280, 16385, 1);
			det->Temperature = 300;
			
			// Float_t Pos[3] = {80, 80, 1};
			// Float_t Size[3] = {80, 80, 2};
			// det->ElRectangle(Pos, Size, 0, 20); ///how to use?

			det->SetUpElectrodes();
			det->SetBoundaryConditions();
			//define the space charge
			TF3 *f2 = new TF3("f2", "x[0]*x[1]*x[2]*0+[0]", 0, 3000, 0, 3000, 0, 3000);
			f2->SetParameter(0, -2);
			det->NeffF = f2;

			// calculate weigting field
			// calculate electric field
			det->CalField(0);
			det->CalField(1);
			// set entry points of the track
			det->enp[0] = 40;
			det->enp[1] = 27.5;
			det->enp[2] = 250;
			det->exp[0] = 40;
			det->exp[1] = 27.5;
			det->exp[2] = 50;

			// switch on the diffusion
			det->diff = 1;
			// Show mip track
			TCanvas c1;
			c1.cd();
			det->ShowMipIR(30);
			TCanvas c3;
			c3.cd();
			// int number_pairs_si =76;	//silicon each um
			// int number_pairs = 51;		//silicon carbide
			det->MipIR(1520);     // hole-electron paris of silicon
			//det->MipIR(1020);			 // hole-electron paris of silicon carbide
			det->sum->Draw("HIST");		 //total current
			// det->neg->Draw("HIST same"); //electrons current
			// det->pos->Draw("HIST same"); // hole current
			c3.SaveAs("txt/Si_induced_current.C");
			//theApp.Run();
		} // End "3D"

		if (!strcmp(argv[i], "2D"))
		{
			TF1 *neff = new TF1("neff", "[0]+x[0]*0", 0, 1000);
			neff->SetParameter(0, 10);
			KPad det(100,100);
			det.Neff = neff;
			det.SetDriftHisto(10e-9, 5000);
			det.Voltage = -500;
			det.SetUpVolume(1);
			det.SetEntryPoint(50, 0, 0.5); //set entry point of the track
			det.SetExitPoint(50, 100, 0.5);
			det.SetUpElectrodes();
			det.SStep = 0.1; // set the drift step of e-h pairs
			det.Temperature = 300; // set the operation temperature
			
 			//set exit point of the track
			// switch on the diffusion
			det.diff = 1;

			//--------------------------------------basic information---------------------------------------//

			 TGraph *ElField;						  // electric field
			// TGraph *ElPotential;					  // electric potential
			 TCanvas c2("Plots", "Plots", 1400, 1000); //open canvas
			c2.Divide(2, 3);						  //divide canvas

			// // // // // // electic field
			 c2.cd(1);
			 ElField = det.DrawPad("f");
			 ElField->SetTitle("Electric field");
			 ElField->GetXaxis()->SetTitle("depth [#mum] (0 is voltage applied electrode)");
			 ElField->GetYaxis()->SetTitle("E [V/#mum]");

			// // // // // // electric potential  or weighting potential ???
			// c2.cd(2);
			// ElPotential=det.DrawPad("p");
			// ElPotential->SetTitle("Electric Potential");
			// ElPotential->GetXaxis()->SetTitle("depth [#mum] (0 is voltage applied electrode)");
			// ElPotential->GetYaxis()->SetTitle("U [V]");

			// // // // // // doping distribution
			 c2.cd(3);
			 TF1 *neffGc;
			 neffGc = neff->DrawCopy();
			 neffGc->SetRange(0, 21);
			 neffGc->SetTitle("Doping profile (full detector)");
			 neffGc->GetXaxis()->SetTitle("depth [#mum] (0 is voltage applied electrode)");
			 neffGc->GetYaxis()->SetTitle("N_{eff} [10^{12} cm^{-3}]");
			 neffGc->GetYaxis()->SetTitleOffset(1.5);

			// // // // // // induced current
			//c2.cd(4);
			// gStyle->SetOptStat(kFALSE);
			// det.MipIR(1000);
			// TH1F *ele_current = (TH1F *)det.sum->Clone();
			// det.sum->SetLineColor(3);
			// det.neg->SetLineColor(4);
			// det.pos->SetLineColor(2);
			// det.sum->Draw("HIST");		//plot total induced current
			// det.neg->Draw("SAME HIST"); //plot electrons' induced current
			// det.pos->Draw("SAME HIST"); //plot holes' induced current
			// c2.Update();
			// KElec tct;
			// double test_out=tct.CSAamp(ele_current, 1, 1, 21.7, 0.66, 4, 5);
			// cout<<test_out<<endl;
			//tct.preamp(ele_current);
			// tct.CRshape(40e-12, 50, 10000, ele_current);
			// // tct.RCshape(40e-12, 50, 10000, ele_current);
			// // tct.RCshape(40e-12, 50, 10000, ele_current);
			// // tct.RCshape(40e-12, 50, 10000, ele_current);
			// float rightmax = 1.1 * ele_current->GetMinimum();
			// float scale = gPad->GetUymin() / rightmax;
			// ele_current->SetLineColor(6);
			// ele_current->Scale(scale);
			//ele_current->Draw("HIST");
			// auto axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
			// 					   gPad->GetUxmax(), gPad->GetUymax(), rightmax, 0, 510, "+L");
			// axis->SetLineColor(6);
			// axis->SetTextColor(6);
			// axis->SetTextSize(0.02);
			// axis->SetLabelColor(6);
			// axis->SetLabelSize(0.02);
			// axis->SetTitle("ampl(mV)");
			// axis->Draw();

			// auto legend = new TLegend(0.6, 0.3, 0.9, 0.6);
			// legend->AddEntry(det.sum, "sum", "l");
			// legend->AddEntry(det.neg, "electron", "l");
			// legend->AddEntry(det.pos, "hole", "l");
			// legend->AddEntry(ele_current, "current after electric", "l");
			// legend->SetBorderSize(0);
			// legend->Draw();

			// // // // // // charge drift
			c2.cd(5);
			det.ShowMipIR(150);

			// // // // // // e-h pairs
			// c2.cd(6);
			// det.MipIR(5100);
			// det.pairs->Draw("HIST");

			//------------------------------end------------------------------------------------------


			//----------------------------current before electric and after electric---------------------------------------//
			
			// TCanvas c3("Plots", "Plots", 1000, 1000);
			// det.MipIR(50);
			// TH1F *Icurrent = (TH1F *)det.sum->Clone();
			// Icurrent->SetLineColor(2);
			// Icurrent->SetTitle("induced current");
			// Icurrent->Draw("HIST");
			// c3.Update();
			// KElec tct;
			// tct.preamp(det.sum);
			// tct.CRshape(40e-12, 50, 10000, det.sum);
			// // tct.RCshape(40e-12, 50, 10000, det.sum);
			// // tct.RCshape(40e-12, 50, 10000, det.sum);
			// // tct.RCshape(40e-12, 50, 10000, det.sum);
			// det.sum->Scale(1000);
			// float rightmax = 1.1 * det.sum->GetMinimum();
			// float scale = gPad->GetUymin() / rightmax;
			// det.sum->SetLineColor(4);
			// det.sum->Scale(scale);
			// det.sum->Draw("HIST SAME");
			// auto axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
			// 					   gPad->GetUxmax(), gPad->GetUymax(), rightmax, 0, 510, "+L");
			// axis->SetLineColor(4);
			// axis->SetTextColor(4);
			// axis->SetTextSize(0.02);
			// axis->SetLabelColor(4);
			// axis->SetLabelSize(0.02);
			// axis->SetTitle("ampl(mV)");
			// axis->Draw();

			//--------------------------------------end--------------------------------//


			// // //--------------------timing scan------------------------------
			//ofstream outfile;
			//outfile.open("time_resolution_2021_4_21.txt");
			
			//TH1F *charge_total = new TH1F("t_charge", "t_charge", 200, -2, 0);
			//TH1F *CSA_time_resolution = new TH1F("CSA_time_resolution", "CSA_time_resolution", 2000, 0, 50);
			//Int_t i=0;
			/*do
			{
		 	TH1::AddDirectory(kFALSE);

			// // // induced current
			TCanvas c4("Plots", "Plots", 1000, 1000);
			c4.cd();
			det.MipIR(1000);
			TH1F *Icurrent = (TH1F *)det.sum->Clone();
			Icurrent->Draw("HIST");
			Double_t charge_t;
			charge_t = Icurrent->Integral() * ((Icurrent->GetXaxis()->GetXmax() - Icurrent->GetXaxis()->GetXmin()) / Icurrent->GetNbinsX()) * 1e15;
			charge_total->Fill(charge_t);
			std::cout << "charge:" << charge_t << std::endl;

			std::string output;
			output = "out/sic_2021_4_21_ic/sic_events_" + std::to_string(i) + ".C";
			const char *out = output.c_str();
			c4.SaveAs(out);

			c4.Update();
			// //current after electric
			TCanvas c5("Plots", "Plots", 1000, 1000);
			c5.cd();
			KElec tct;
			Double_t CSA_thre_time = tct.CSAamp(det.sum, 0.3, 0.6, 75, 10, 0.66, 3);
			// // taurise=1ns taufall=1ns  detector capacitance=75 pF  CSATransImp=10
			// // CSA noise=0.66  CSA_threshold=3

			std::cout << "CSA_thre_time:" << CSA_thre_time << std::endl;
			CSA_time_resolution->Fill(CSA_thre_time);
			//outfile<<CSA_thre_time<<endl;
			// //tct.preamp(det.sum);
			// // tct.CRshape(40e-12, 50, 10000, det.sum);
			// // tct.RCshape(40e-12, 50, 10000, det.sum);
			// // tct.RCshape(40e-12, 50, 10000, det.sum);
			// // tct.RCshape(40e-12, 50, 10000, det.sum);
			det.sum->Draw("HIST");

			std::string output_1;
			output_1 = "out/sic_2021_4_21_ec/sic_events_" + std::to_string(i)+".C";
			const char *out_1 = output_1.c_str();
			c5.SaveAs(out_1);
			std::cout<<output_1<<std::endl;
			i++;

			} while (i<3000);
			*/
			///////////////e-h pairs
			//TCanvas c6("Plots", "Plots", 1000, 1000);
			//c6.cd();
			//charge_total->Draw("HIST");
			//c6.SaveAs("eh_pairs.pdf");

			//TCanvas c7("Plots", "Plots", 1000, 1000);
			//c7.cd();
			//CSA_time_resolution->Draw("HIST");
			//c7.SaveAs("time_resolution_2021_4_21.pdf");
			//outfile.close();

			// --------------------end---------------

		theApp.Run();
		} // End "2D"
		if(!strcmp(argv[i], "tracs"))
		{
			std::string carrierFile;
			std::vector<std::string> carrierThread_fileNames;
			TString transferFun;
			TH1D *i_rc;
			TH1D *i_conv;
			//	TH1D *i_ramo;
		

			int counted_numLines = 0;
			int number_of_lines = 0;
			int nns, count2, count;
			std::string line;
			std::string temp;


			if(argc==1)
			{
				num_threads = std::thread::hardware_concurrency(); // No. of threads = No. of cores
			}

			if(argc==2)
			{
				num_threads = atoi(argv[1]);
				if (num_threads == 0)
				{
					num_threads = 1;
				}
			}

			if(argc==3)
			{
				num_threads = atoi(argv[1]);
				fnm = argv[2];
			}

			//file creation
			for (int i = 0; i < num_threads; ++i)
			{
				temp = std::to_string(i);
				carrierThread_fileNames.push_back("file" + temp);
			}

			utilities::parse_config_file(fnm, carrierFile);
			std::ifstream in(carrierFile);
	
			while (in)
			{
				for (int i = 0 ; i < num_threads ; i++)
				{
					std::ofstream outfile;
					const char * c = carrierThread_fileNames[i].c_str();
					outfile.open(c, std::ofstream::out | std::ofstream::app);
					if (in) std::getline(in, line);
					outfile << line << std::endl;
					outfile.close();
				}
		
			}
			in.close();

			spread_into_threads();
			double timeSteps = (int) std::floor(max_time / dTime);
			int total_sizeZ = vector_voltValues.size() * vector_zValues.size();
			int total_sizeY = vector_voltValues.size() * vector_yValues.size();


			if (scanType == "edge")
			{
				if (vector_yValues.size() > 1)
				{
					std::cout << "This execution is wrongly configured with edge-TCT; Check parameters." << std::endl;
					std::quick_exit(1);
				}

				i_rc_array.resize(total_sizeZ);
				i_ramo_vector.resize(total_sizeZ);
				i_conv_vector.resize(total_sizeZ);
				vItotals.resize(total_sizeZ);
				for (int i = 0; i < total_sizeZ ; i++)
				{
					vItotals[i].resize(timeSteps);
				}

			}


			if (scanType == "top" || scanType == "bottom")
			{
				if (vector_zValues.size() > 1)
				{
					std::cout << "This execution is wrongly configured with top/bottom-TCT; Check parameters." << std::endl;
					std::quick_exit(1);
				}

				i_rc_array.resize(total_sizeY);
				i_ramo_vector.resize(total_sizeY);
				i_conv_vector.resize(total_sizeY);
				vItotals.resize(total_sizeY);
				for (int i = 0; i < total_sizeY ; i++)
				{
					vItotals[i].resize(timeSteps);
				}
			}

			for (int i = 0 ; i < vItotals.size() ; i++)
			{
				for (int j = 0 ; j < vItotals[i].size() ; j++)
				{
					vItotals[i][j] = 0;
				}
			}


			TRACSsim.resize(num_threads);
			t.resize(num_threads);
			for (int i = 0; i < num_threads; ++i)
			{
				t[i] = std::thread(call_from_thread, i, std::ref(carrierThread_fileNames[i]), vector_zValues, vector_yValues, vector_voltValues);
			}
			for (int i = 0; i < num_threads; ++i) 
			{
				t[i].join();
			}



			for (int i = 0 ; i < vItotals.size(); i++)
			{
				for (int j = 0; j < num_threads; j++)
				{
					vItotals[i] = vItotals[i] + TRACSsim[j]->vSemiItotals[i];// + temp_s;
				}
			}


			//Current to rc array of TH1D -> root file
			if (global_TF != "NO_TF")
			{
				for (int i = 0 ; i < i_ramo_vector.size(); i++ )
				{
		  			std::cout << "Yes" << std::endl;
					transferFun = global_TF;
					TString htit, hname;
					TString htit2, hname2;

					htit.Form("ramo_%d%d", 0, count2);
					hname.Form("Ramo_current_%d_%d", 0, count2);
					//i_ramo = new TH1D(htit,hname, timeSteps, 0.0, max_time);

					TH1D *i_ramo = new TH1D(htit.Data(), hname.Data(),  timeSteps, 0.0, max_time);

					htit2.Form("ramo_conv%d%d", 0, count2);
					hname2.Form("Ramo_current_%d_%d", 0, count2);

					for (int k = 0 ; k < timeSteps; k++ )
					{
						i_ramo->SetBinContent(k+1, vItotals[i][k]);
					}
			
					i_conv = H1DConvolution(i_ramo , capacitance*1.e12, count, transferFun);

/*
			// ******************************** D I R T Y   FIX : Shift histogram back !!! DELETE ME !!! *********************************************************************
			// ******************************** D I R T Y   FIX : Shift histogram back !!! DELETE ME !!! *********************************************************************
			// ******************************** D I R T Y   FIX : Shift histogram back !!! DELETE ME !!! *********************************************************************
			Int_t Nbins = TMath::Nint( max_time/dTime ) ;
			TH1D *hitf = new TH1D("hitf","hitf", Nbins ,0.,max_time) ;
			Int_t istart = i_conv->FindBin( -20e-9) , iend=i_conv->FindBin( 20.e-9) ;
			Int_t cont =  hitf->FindBin( 0.102141e-09 );
			for (Int_t k=istart ; k < iend ; k++ ) 
			{
				//std::cout<< i_conv->GetBinContent( k ) <<std::endl;
				hitf->SetBinContent( cont , i_conv->GetBinContent( k+1 ) );
				cont++;
			}
			i_conv_vector[i] = hitf;

			// ****************************************************************************************************************************************************************
*/

					i_conv_vector[i] = i_conv;
					vItotals[i].resize(i_conv_vector[i]->GetNbinsX());
					i_ramo = nullptr;
					i_conv = nullptr;
					count2++;
					count++;
				}

				for (int i = 0 ; i < i_conv_vector.size() ; i++)
				{
					for (int j = 0; j < i_conv_vector[i]->GetNbinsX(); j++  )
					{
						vItotals[i][j] = i_conv_vector[i]->GetBinContent(j+1);
					}
				}
			}
			if (global_TF == "NO_TF")
			{
				for (int i = 0 ; i < i_rc_array.size(); i++ )
				{
					TString htit, hname;
					htit.Form("ramo_rc%d%d", 0, count2);
					hname.Form("Ramo_current_%d_%d", 0, count2);
					i_rc = new TH1D(htit,hname, timeSteps, 0.0, max_time);
					for (int k = 0 ; k < timeSteps; k++ )
					{
						i_rc->SetBinContent(k+1, vItotals[i][k]);
					}
					i_rc_array[i] = i_rc;
					vItotals[i].resize(i_rc_array[i]->GetNbinsX());
					i_rc = nullptr;
					count2++;
				}
				for (int i = 0 ; i < i_rc_array.size() ; i++)
				{
					for (int j = 0; j < i_rc_array[i]->GetNbinsX(); j++  )
					{
						vItotals[i][j] = i_rc_array[i]->GetBinContent(j+1);
					}
				}

			}

			//write output to single file!
			TRACSsim[0]->write_to_file(0);
			for (int i = 0 ; i < num_threads ; i++)
			{
				const char * c = carrierThread_fileNames[i].c_str();
				remove(c);
			}


			for (uint i = 0; i < TRACSsim.size(); i++)
			{
				delete TRACSsim[i];
			}

			std::quick_exit(1);
		}//end tracs
	}

}
 
#endif

// // Random noise simulation // //

// TRandom3 r;
// double CSAx1 = 0;
// //
// TH1F *charge_total = new TH1F("sim_noise", "sim_noise", 1000, -13, 12);
// for (int i = 0; i < 300000; i += 1)
// {

// 	r.SetSeed(0);

// 	CSAx1 = r.Gaus(0.66, 2.36);

// 	charge_total->Fill(CSAx1);
// 	//
// }
// charge_total->Draw();

//endl












///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
















// // code from here based on TRACS  //  //

// Head files appear in Kdetsim-based part may be repeated for now 2021.04.28  

#include <thread>//C++ standard thread class
#include <boost/asio.hpp>//do not have this file
#include <stdio.h>//C++/C standard IO

//#include <TRACSInterface.h>
/*------------------------------------TRACSInterface.h start----------------------------------------*/
#ifndef TRACSINTERFACE_H
#define TRACSINTERFACE_H

#include <stdlib.h>
#include <iostream>
#include <ctime>

#include <iterator>
#include <limits>  // std::numeric_limits
#include <cmath>
#include <functional>
#include <vector>
#include <array>

#include <TFile.h>
#include <TF1.h>
#include <TH1D.h> // 1 Dimesional ROOT histogram
#include <TTree.h>


//#include <TMeas.h>
/*------------------------------------TMeas.h start------------------------------------------*/


#include "TObject.h"
#include "TDatime.h"

#define TMEAS_H
#ifndef TMEAS_H
#define TMEAS_H

#include "TObject.h"
#include "TDatime.h"


//CONSTANTS      ----------------------
#define NXCHAR 512

using namespace std ;       //Evita usar std::cin;

class TMeas: public TObject
{
public:
     TMeas() ;
     ~TMeas() ;
     UShort_t   NV ;	          //! (DNS) Number of voltage scans	   
     UShort_t   Nx ;	          //! (DNS) Number of different x positions 		   
     UShort_t   Ny ;	          //! (DNS) Number of different y positions 			    
     UShort_t   Nz ;	          //! (DNS) Number of different z positions 		    
     
     Double_t   At ;              //! (DNS) Time step [ns]
     Double_t   Ax ;              //  X step [mum]
     Double_t   Ay ;	          //  Y step [mum]		   
     Double_t   Az ;	          //  Z step [mum]		   

     TDatime    utc   ;           //  Time in UTC 
     UInt_t     event ;           //  DATA-START-event iterative
     
     Double_t   Itot ;            //  Current from HV supply [A]
     Double_t   x    ;            //  X [mm]
     Double_t   y    ;            //  Y [mm]
     Double_t   z    ;            //  Z [mm]

     Int_t   ix    ;              // Slice in X 
     Int_t   iy    ;	          // Slice in Y 
     Int_t   iz    ;	          // Slice in Z 

     Double_t   Temp ;            // Temperature
     Double_t   Vbias ;           //  Bias voltage [V]
     
     //I-DLTS specific
     Double_t pWidth ;            //I-DLTS pulse width
     Double_t pAmpl ;             //I-DLTS pulse amplitude
     

     char       comment[NXCHAR] ; //! (DNS) File comment
     
     Int_t      Nt       ;	  //Number of bins in scope		    
     Double_t   *volt	 ;        //[Nt]
     Double_t   *time	 ;        //[Nt]
     //One possible change: Qt is not neccessary to be declared as a pointer.
     Double_t * Qt	 ;            //[Nt]

     
     UShort_t   Setup    ;          //! 0=OldTCT , 1=eTCT (default), 2=TCT+
                                    //! 3=Hamburg, 4=IFCA , 5=TRACS, 6=IDLTS
             
     UInt_t     Ntevent ;          //! Total number of DATA-START

    ClassDef(TMeas,1)  //Edge-TCT data class
} ;


#endif
/*------------------------------------TMeas.h end------------------------------------------*/

//#include <TWaveform.h>
/*------------------------------------TWaveform.h start------------------------------------------*/

#include <TH1.h>
#define TWAVEFORM_H
#ifndef TWAVEFORM_H
#define TWAVEFORM_H

//#include <TObject.h>
#include <TH1.h>


#define RTIME 0.8

//#include "TMeas.h"

using namespace std ;       //Evita usar std::cin;


class TWaveform : public TObject
{		
  public:
	  	  
       TWaveform() ;
       TWaveform( TMeas *em ) ;
       TWaveform( int nel , double *tin , double *vin , double Bias  ) ;
       ~TWaveform() ;
	  
	  /*------------------- DATA MEMBERS ---------------------*/

       Double_t *volt   ;			        //!
       Double_t *time  ;			            //!
	  double   Vbias   ;			        //!

	  /* Histograms */
	  TH1D *hvt    ;                        //! Histo of Voltage vs time
	  TH1D *hbl    ;		        //! Only bline
	  TH1D *hpbl   ;			//! Projected bline (vplot of bline)
	  TH1D *hvtbl  ;			//! Voltage vs time, baseline corrected
	  
	  double LPower ;                       //Laser power calculated as: Ipd [A]/0.77 [A/W], 
	                                        //where Ipd=Photocurrent is simply Vmax [V]/50 [Ohm]
	  double LNph ;                         //Number of photons injected in photodiode: Int(Vpd*dt)/(50*0.77*Ephoton)
	  
	  int      GetNbins()	 ;  //Ner of time points      //!
	  int      GetPolarity() ;			      //!
	  Double_t GetAbsVmax()  ;			      //!
	  double   GetVmax()	 ;			      //!
	  double   GetVmin()	 ;			      //!
	  double   GettVmax()	 ;			      //!
	  double   GettVmin()	 ;			      //!
	  double   BlineGetMean();			      //!
	  double   BlineGetRMS() ;			      //!
	  double   GetTleft( )   ;			      //!
	  double   GetTright( )  ;			      //!
	  double   GetTrms()	 ;			      //!
	  double   GetRiseTime( double fraction=0.9 ) ;	      //!
	  Double_t GetCharge( double tinf , double tsup ) ;   
	  Double_t RGetCharge( double tinf , double tsup ) ;   
          void     CalcRunningCharge( TMeas *em ) ;           //!
	  
	  

  private:						     
          int    Nbins     ;				    //!
	  int    Polarity  ;				    //!
	  
	  Double_t Vmax      ;
	  Double_t Vmin      ;
	  int      iVmax     ;
	  int      iVmin     ;
	  Double_t tVmax     ;
	  Double_t tVmin     ;
	   
	  double tleft     ;
	  double tright    ;
	  double trms      ;
	  int    itleft    ;
	  int    itright   ;
	  
	  double BlineMean ; 
	  double BlineRMS  ;
	  
	  Double_t Q50       ;   //From [tleft,tright+10 ns]
	  Double_t Qtot      ;   //From [time[0],time[Nt]]
	  
	  double RiseTime  ;
	  
	  
	  void   CalcVmaxmin()  ;						       //!
	  void   CalcPolarity() ;						       //!
	  void   CalcBline()  ; 						       //!
	  double CalcRiseTime( double fraction=RTIME )  ;			       //!
	  void   CreateHistos() ;						       //!		       //!
   
  /*protected:  //classes that inherit from TWaveform can access these methods*/
    
	  void   CalcSignalTimeL( double fraction , double &TimeL , int &iTimeL )  ;   //!
	  void   FitSignalTimeL( double &TimeL , int &iTimeL )  ;   //!
	  void   CalcSignalTimeR( double fraction , double &TimeR , int &iTimeR )  ;   //!
          void   CalcQTimeR(  double fraction , double &TimeR , int &iTimeR ) ;        //!
	  void   FitSignalTimeR( double &TimeR , int &iTimeR )  ;   //!
	  void   CalcSignalTimeLR( )  ;			                               //!
	  
  ClassDef(TWaveform,3) ; //ROOT RTTI
	
};

#endif

/*------------------------------------TWaveform.h end------------------------------------------*/
//#include <TMeasHeader.h>
/*------------------------------------TMeasHeader.h start------------------------------------------*/

#include <TString.h>
#include <TVectorD.h>
#define TMEASHEADER_H
#ifndef TMEASHEADER_H
#define TMEASHEADER_H

//#include <TObject.h>
#include <TString.h>
//#include <TDatime.h>
#include <TVectorD.h>

using namespace std ;       //Evita usar std::cin;

class TMeasHeader: public TObject
{

   public:
     TMeasHeader() ;
     TMeasHeader( Int_t NV ) ;
     ~TMeasHeader();
     
     UShort_t   NV ;	          //Number of voltage scans	   
     UShort_t   Nx ;	          //Number of different x positions 		   
     UShort_t   Ny ;	          //Number of different y positions 			    
     UShort_t   Nz ;	          //Number of different z positions 		    
     Short_t    Polarity ;	  //Polarity of the FILE! 		    
     
     Double_t   At ;              //Time step [ns]
     Double_t   Ax ;              //X step [mum]
     Double_t   Ay ;              //Y step [mum]                 
     Double_t   Az ;              //Z step [mum]                 
     Double_t   Nt ;              //Number of steps in time
     
     TString    comment ;         //File comment
     Double_t   Temp ;            //Temperature
     Double_t   Hum ;             //Humidity
     Double_t   Lambda ;          //Lambda
     Double_t   Power ;           //Pulse Power
     Double_t   Gain ;            //Amp. Gain
     Double_t   Freq ;            //Frequency
     Int_t      Nav ;             //Ner of averages
     Double_t   t0 ;              //Scope t0
     Double_t   Position ;        //Position of the sensor in mm (mostly for TCT standard)
  
     Double_t   Cend ;            //Capacitance used for RC deconvolution of the file
     UShort_t   Setup ;           //0=OldTCT , 1=eTCT (default), 2=TCT+
                                  //3=Hamburg, 4=IFCA , 5=TRACS, 6=TPA
     
     Int_t      Illum ;           //Illumination 1 = top , 0=side, -1 = bottom

     TDatime date0 ;              //Start date
     TDatime datef ;              //End date
     
     Double_t   iann ;            //annealing step
     Double_t   tann ;            //Specific annealing time at TempAnn
     Double_t   Etann ;           //Accumulated annealing time at 60C for this sample
     Double_t   TempAnn ;         //Annealing temperature
     Double_t   Fluence ;         //Fluence

     TVectorD   vVbias ;          //Vector with bias voltages
     TVectorD   vIleak ;          //Vector with leakage current
     TVectorD   vTemp ;           //Vector with temperature readings
     
     //Methods
     UShort_t GetNV() { return NV ;}
     UShort_t GetNx() { return Nx ;}
     UShort_t GetNy() { return Ny ;}
     UShort_t GetNz() { return Nz ;}
     UShort_t GetNav(){ return Nav ;}

     Double_t GetAt() { return At ;}
     Double_t GetAx() { return Ax ;}
     Double_t GetAy() { return Ay ;}
     Double_t GetAz() { return Az ;}
     Double_t GetGain()    { return Gain ;}
     Double_t Getiann()    { return iann ;   }  
     Double_t Gettann()    { return tann ;   }  
     Double_t GetEtann()   { return Etann ;  } 
     Double_t GetTempAnn() { return TempAnn ;} 
     Double_t GetFluence()  { return Fluence ;} 
     Double_t GetFrequency()    { return Freq;}
     Double_t GetFilePolarity() { return Polarity ;}

     Double_t GetTemperature()    {return Temp ;} 	    
     TString  GetComment() { return comment; } 

    ClassDef(TMeasHeader,3)  //Edge-TCT data header class
} ;

#endif
/*------------------------------------TMeasHeader.h end------------------------------------------*/

#include <TMath.h>
#include <TPad.h>

//#include <SMSDetector.h>
/*------------------------------------SMSDetector.h start------------------------------------------*/
#ifndef SMSDETECTOR_H
#define SMSDETECTOR_H


//#include <cmath> 
//#include <limits>  // std::numeric_limits

//#include <array>
#include <string>

#include <dolfin.h>

#include "Poisson.h"// This code conforms with the UFC specification version 1.5.0 and was automatically generated by FFC version 1.5.0.
#include "Gradient.h"// This code conforms with the UFC specification version 1.5.0 and was automatically generated by FFC version 1.5.0.

//#include <SMSDSubDomains.h>
/*------------------------------------SMSDSubDomains.h start------------------------------------------*/
#ifndef SMSSUBDOMAINS_H
#define SMSSUBDOMAINS_H

//#include <dolfin.h>

using namespace dolfin;

class CentralStripBoundary: public SubDomain
{
  private:
    double _pitch;      // strip pitch
    double _width;      // strip width
    int _nns;           // number of neighbor strips
  public:
    CentralStripBoundary(double pitch, double width, int nns );
    bool inside(const Array<double>& x, bool on_boundary) const;
};

class NeighbourStripBoundary: public SubDomain
{
  private:
    double _pitch;      // strip pitch
    double _width;      // strip width
    int _nns;           // number of neighbor strips
  public:
    NeighbourStripBoundary(double pitch, double width, int nns );
    bool inside(const Array<double>& x, bool on_boundary) const;
};


class BackPlaneBoundary: public SubDomain
{
  private:
    double _x_min;      // min x value
    double _x_max;      // max x value
    double _depth;      // detector depth
  public:
    BackPlaneBoundary(double x_min, double x_max, double depth);
    bool inside(const Array<double>& x, bool on_boundary) const;
};

class PeriodicLateralBoundary: public SubDomain
{
  private:
    double _x_min;      // min x value
    double _x_max;      // max x value
    double _depth;      // detector depth
  public:
    PeriodicLateralBoundary(double x_min, double x_max, double depth);
    bool inside(const Array<double>& x, bool on_boundary) const;
    void map(const Array<double>& x, Array<double>& y) const;
};


/********************************/

class CentralStripBoundaryWP: public SubDomain
{
  private:
    double _pitch;      // strip pitch
    double _width;      // strip width
    int _nns;           // number of neighbor strips
    double _depletion_width;
  public:
    CentralStripBoundaryWP(double pitch, double width, int nns, double depletion_width);
    bool inside(const Array<double>& x, bool on_boundary) const;
};

class NeighbourStripBoundaryWP: public SubDomain
{
  private:
    double _pitch;      // strip pitch
    double _width;      // strip width
    int _nns;           // number of neighbor strips
    double _depletion_width;

  public:
    NeighbourStripBoundaryWP(double pitch, double width, int nns, double depletion_width );
    bool inside(const Array<double>& x, bool on_boundary) const;
};

class BackPlaneBoundaryWP: public SubDomain
{
  private:
    double _x_min;      // min x value
    double _x_max;      // max x value
    double _depletion_width;      // detector depth
  public:
    BackPlaneBoundaryWP(double x_min, double x_max, double depletion_width);
    bool inside(const Array<double>& x, bool on_boundary) const;
};

/********************************/

#endif // SMSSUBDOMAINS_H
/*------------------------------------SMSDSubDomains.h end------------------------------------------*/

using namespace dolfin;

class SMSDetector
{
private:
	// detector characteristics
	double _pitch; // in microns
	double _width; // in microns
	double _depth; // in microns
	double _tempK; // temperature in Kelvin
	double _trapping_time; // radiation damage effect
	double _fluence; // irradiation fluence (in neq)
	int _nns; // numbre of neighbouring strips
	char _bulk_type; // p or n
	char _implant_type; // n or p
	std::string _neff_type;
	std::vector<double> _neff_param; // Neff parametrization

    // For avalanche region
    std::string _avalanche_flag;
    std::array<std::array<double, 3>, 2> _doping_param;
    
	// some useful derived variables
	double _x_min; // in microns
	double _x_max; // in microns
	double _y_min; // in microns
	double _y_max; // in microns

	double _vdep; // depletion voltage
	double _v_bias;// bias
	double _v_backplane;
	double _v_strips;
	bool _depleted;

	// poisson term for solving electric field
	double _f_poisson;

	//Variables used for diffusion
	int _diffusion; //If diffusion is ON by user. 0 for NO, 1 for YES
	double _depletion_width;
	double _dt;

	// Meshing parameters
	int _n_cells_x;
	int _n_cells_y;

	// meshes (one for each could be used)
	RectangleMesh _mesh; // mesh for both weighing and drifing potential

	// mesh subdomains
	PeriodicLateralBoundary _periodic_boundary;
	CentralStripBoundary _central_strip;
	NeighbourStripBoundary _neighbour_strips;
	BackPlaneBoundary _backplane;

	//New subdomains to solve no depleted detectors
	//CentralStripBoundaryWP _central_stripUnderDep;
	//NeighbourStripBoundaryWP _neighbour_stripsUnderDep;
	//BackPlaneBoundaryWP _backplaneUnderDep;

	// Poisson PDE Function Space
	Poisson::FunctionSpace _V_p;
	Poisson::BilinearForm _a_p;
	Poisson::LinearForm _L_p;

	// Gradient PDE Function Space
	Gradient::FunctionSpace _V_g;
	Gradient::BilinearForm _a_g;
	Gradient::LinearForm _L_g;

	// potentials
	Function _w_u;  // function to store the weighting potential
	Function _d_u;  // function to store the drifting potential

	// fields
	Function _w_f_grad; // function to store the weighting field (vectorial)
	Function _d_f_grad; // function to store the drifting field (vectorial)

public:
	// default constructor and destructor
	SMSDetector(double pitch, double width, double depth, int nns, char bulk_type, char implant_type,
                std::string avalanche_flag, std::array<std::array<double, 3>, 2> doping_param,
                int n_cells_x = 100, int n_cells_y = 100, double tempK = 253., double trapping = 9e300,
			double fluence = 0.0, std::vector<double> neff_param = {0}, std::string neff_type = "Trilinear", int diffusion = 0, double dt = 300);
	~SMSDetector();
	// set methods
	//void calculate_mod();
	void set_voltages(double v_bias, double v_depletion);
	void set_pitch(double pitch);
	void set_width(double width);
	void set_depth(double depth);
	void set_nns(int nns);
	void set_bulk_type(char bulk_type);
	void set_implant_type(char implant_type);
	void set_n_cells_x(int n_cells_x);
	void set_n_cells_y(int n_cells_y);
	void set_derived(); // Properly sets all derived quantities
	void set_temperature(double temperature);
	void set_trapping_time(double trapping_tau);
	void set_fluence(double fluencia);
	void setFitParameters(std::vector<double> fitParameters);
	void set_neff_type(std::string newApproach);
	// solve potentials
	void solve_w_u();
	void solve_d_u();
	void solve_w_f_grad();
	void solve_d_f_grad();

	// get methods
	Function * get_w_u();
	Function * get_d_u();
	Function * get_w_f_grad();
	Function * get_d_f_grad();
	RectangleMesh * get_mesh();
	double get_x_min();
	double get_x_max();
	double get_y_min();
	double get_y_max();
	int get_n_cells_x();
	int get_n_cells_y();
	double get_temperature();
	double get_trapping_time();
	double get_fluence();
	double get_depth();
	double get_pitch();
	double get_width();
	int get_nns();
	char get_bulk_type();
	char get_implant_type();
	double get_vbias();
	double get_vdep();
	int diffusionON();
	double get_neff();
	double get_depletionWidth();
	double calculate_depletionWidth();
	double get_dt();



	// some other methods
	bool is_out(const std::array< double,2> &x);
	bool is_out_dep(const std::array< double,2> &x);

};

#endif // SMSDETECTOR_H
/*-------------------------------------SMSDetector.h end-----------------------------------------*/
//#include <Source.h>
/*------------------------------------Source.h start------------------------------------------*/
#ifndef SOURCE_H // SOURCE_H 
#define SOURCE_H
//#include <cmath>
//#include <array>

//#include <TH1D.h>
//#include <TFile.h>

//#include <dolfin.h>

using namespace dolfin;


/*
 * SOURCE TERM
 *
 * This class contains the source term for the poisson equation
 * Space charge distribution is parametrized here.
 *
 *
 */
 class Source:public Expression
 {
 public:


	 //Declaration and default values. They will taken from steering file.
	// Concentration in transition and extremal points
	double y0 = -25.; // Neff(z0)
	double y1 = 0.02; // Neff(z1)
	double y2 = y1+y1*10; // Neff(z2)
	double y3 = 33; // Neff(z3)

	// Define diferent zones in depth
	double z0 = 0.;
	double z1 = 120.;
	double z2 = 220.;
	double z3 = 300.;
	std::string NeffApproach = "Triconstant";

    
	void eval(Array<double>& values, const Array<double>& x) const;

    std::string get_NeffApproach() const
    {
        return NeffApproach;
    }        
	void set_NeffApproach(std::string Neff_type)
	{
		NeffApproach = Neff_type;
	}

    void set_avalanche_doping_param( const std::array<std::array<double, 3>, 2>& param )
    {
        _peak_height[0]    = param[0][0];
        _peak_position[0] = param[0][1];
        _gauss_sigma[0]  = param[0][2];

        _peak_height[1]    = param[1][0];
        _peak_position[1] = param[1][1];
        _gauss_sigma[1]  = param[1][2];
    }
    void set_bulk_doping_param( const double& param )
    {
        _f_poisson = param;
    }
    
	void set_y0(double newValue)
	{
		y0 = newValue;
	}

	void set_y1(double newValue)
	{
		y1 = newValue;
	}

	void set_y2(double newValue)
	{
		y2 = newValue;
	}

	void set_y3(double newValue)
	{
		y3 = newValue;
	}

	void set_z0(double newValue)
	{
		z0 = newValue;
	}

	void set_z1(double newValue)
	{
		z1 = newValue;
	}

	void set_z2(double newValue)
	{
		z2 = newValue;
	}

	void set_z3(double newValue)
	{
		z3 = newValue;
	}

    void save_Neff_dist(double ymin, double ymax);
        
private:

    // For the effective doping parameters
    std::array<double, 2> _peak_height;
    std::array<double, 2> _peak_position;
    std::array<double, 2> _gauss_sigma;
    double _f_poisson;

    static const double _elementary_charge;
    static const double _vacuum_permittivity;
    static const double _relative_permittivity_silicon; 
    
 };

 #endif // SOURCE_H 
/*-------------------------------------Source.h end-----------------------------------------*/
//#include <Utilities.h>
/*------------------------------------Utilities.h start------------------------------------------*/
#ifndef UTILITIES_H
#define UTILITIES_H


//#include <dolfin.h>

//#include "Poisson.h"
//#include "Gradient.h"

#include <TH2D.h>
//#include <TH1D.h>
//#include <TString.h>
#include <fstream>
#include <iomanip>
#include <utility>
//#include "SMSDetector.h"
//#include <cmath> 
#include <valarray>
//#include <TMath.h>

//#include <Threading.h>
/*------------------------------------Threading.h start--------------------------------------------*/

#ifndef SRC_THREADING_H_
#define SRC_THREADING_H_

//#include <vector>
//#include <TRACSFit.h>

void call_from_thread(int, std::string&, const std::vector<double>& z, const std::vector<double>& y, const std::vector<double>& volt);
void call_from_thread_FitPar(int, std::string&, const std::vector<double>& z, const std::vector<double>& y, const std::vector<double>& volt, const std::vector<Double_t>& par);
void call_from_thread_FitNorm(int, std::string&, const std::vector<double>& z, const std::vector<double>& y, const std::vector<double>& volt, const std::vector<Double_t>& par);

#endif /* SRC_THREADING_H_ */

/*-------------------------------------Threading.h end-----------------------------------------------*/

//#include <TRACSInterface.h>
//#include <Global.h>
/*------------------------------------Global.h start------------------------------------------*/
#ifndef GLOBAL_H
#define GLOBAL_H

//#include <vector>
//#include <valarray>
#include <mutex>
#include <atomic>

//#include <TH1D.h> // 1 Dimensional ROOT histogram


//extern std::vector<std::vector <TH1D*> >  i_ramo_array, i_conv_array;//, i_rc_array;
extern std::vector<TH1D*> i_rc_array, i_ramo_vector, i_conv_vector;
extern int num_threads;
extern std::mutex mtx;
extern std::string fnm;
extern std::string global_TF;
extern std::valarray<std::valarray <double> > vItotals;
//extern std::valarray<std::valarray <double> > i_temp;
extern std::ofstream fileDiffDrift;
//extern std::atomic<double> numberDs;
//extern std::atomic<int> tempNumberDs;
//extern bool printa;

#endif // GLOBAL_H
/*-------------------------------------Global.h end-----------------------------------------*/
//#include "qcustomplot.h"
//#include <QFile>

//#include <array>

using namespace dolfin;

namespace utilities
{
	TH2D export_to_histogram(Function &func, TString hist_name, TString hist_title, int n_bins_x , double x_min, double x_max, int n_bins_y, double y_min, double y_max);

	TH2D export_mod_to_histogram(Function &func, TString hist_name, TString hist_title, int n_bins_x , double x_min, double x_max, int n_bins_y, double y_min, double y_max);

	TH1D export_1D_histogram(Function &func, TString hist_name, TString hist_title, int n_bins_x , double x_min, double x_max, int n_bins_y , double y_min, double y_max);
        
	void write_to_file_row(std::string filename, TH1D *hconv, double temp, double yShift, double height, double voltage);

	void write_to_hetct_header(std::string filename, SMSDetector detector, double C, double dt, std::vector<double> y_shifts, std::vector<double> z_shifts, double landa,
			std::string type, std::string carriers_file, std::vector<double> voltages);

	void write_to_hetct_header(std::string filename, SMSDetector * detector, double C, double dt, std::vector<double> y_shifts, std::vector<double> z_shifts, double landa,
			std::string type, std::string carriers_file, std::vector<double> voltages);

	std::string vector_to_string(std::vector<double> input_list);

    
	void parse_config_file(std::string fileName, double &depth, double &width, double &pitch, int &nns, double &temp, double &trapping, double &fluence,
                           int &nThreads, int &n_cells_x, int &n_cells_y, char &bulk_type, char &implant_type, std::string& skip_flag,
                           std::string& set_avalanche_flag,  std::array<double, 2>& doping_peakheight, std::array<double, 2>& doping_peakpos,
                           std::array<double, 2>& doping_gauss_sigma, double& max_multiplication_factor,
                           int &waveLength, std::string &scanType, double &C, double &dt, double &max_time,
                           double &v_init, double &deltaV, double &v_max, double &v_depletion, double &zInit, double &zMax, double &deltaZ, double &yInit, double &yMax, double &deltaY,
                           std::vector<double> &neff_param, std::string &neffType, double &tolerance, double &chiFinal, int &diffusion, double &fitNorm, std::string& simulation_polarity_flag);
    
	void parse_config_file(std::string fileName, std::string &scanType, double &v_init, double &deltaV, double &v_max, double &v_depletion,
							double &zInit, double &zMax, double &deltaZ, double &yInit, double &yMax, double &deltaY, double &dt, double &max_time, double &capacitance, std::string &transferFun);
	void parse_config_file(std::string fileName, std::string &carrierFile);


	void valarray2Hist(TH1D *hist, std::valarray<double> &valar);

}

#endif // UTILITIES_H

/*-------------------------------------Utilities.h end-----------------------------------------*/
//#include <Carrier.h>
/*------------------------------------Carrier.h start------------------------------------------*/
#ifndef CARRIER_H
#define CARRIER_H

//#include  <valarray>
//#include  <mutex>
//#include <iostream>
//#include <fstream>

#ifndef Q_MOC_RUN  // See: https://bugreports.qt-project.org/browse/QTBUG-22829
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>/////////////////////////////////////////////??????//////////////////////////////////////////
#endif

#include <TRandom3.h>

//#include <CarrierTransport.h>
/*------------------------------------CarrierTransport.h start-----------------------------------------*/
#ifndef CARRIERTRANSPORT_H
#define CARRIERTRANSPORT_H


//#include <dolfin.h>
//#include <TRandom3.h>
//#include <CarrierMobility.h>
/*-------------------------------------CarrierMobility.h start-----------------------------------------*/
#ifndef CARRIERMOBILITY_H
#define CARRIERMOBILITY_H

//#include <cmath>

/*
 ****************JACOBONI MOBILITY**************
 *
 * This class is able to calculate the mobility 
 * of a charge carrier given its type, the 
 * detector thickness and the electric field at
 * de desired point.
 *
 *
 */

class JacoboniMobility
{
  private:
    double _T; // Temperature of the detector
    double _mu0; // Mobility
    double _vsat;
    double _beta;

  public:
    JacoboniMobility(char carrier_type, double T);
		JacoboniMobility();
    ~JacoboniMobility();
    double obtain_mobility(double e_field_mod);
    double obtain_mu0();

};

#endif // CARRIERMOBILITY_H
/*-------------------------------------CarrierMobility.h end-----------------------------------------*/
//#include <Constants.h>
/*-------------------------------------Constants.h start-----------------------------------------*/
#ifndef CONSTANTS_H
#define CONSTANTS_H

//#include <iostream>
//#include "TString.h"
//#include "TMath.h"
/** Constants definition used for diffusion.
 */
//#define ECH		1.602177e-19	// elementary charge in C
//#define kB	    1.38065e-23	    // Boltzmann constant (J/K)
//#define EPS     11.9            // EPS silicon relative permitivity
//#define EPS0    8.85e-14       // EPS0 F/cm  vacuum permitivity

//Mobility electrons	≤1400 cm2 V-1s-1
//Mobility holes	≤450 cm2 V-1s-1


	const double ECH=1.602177e-19;
	const double kB=1.38065e-23;
	const double EPS_=11.9;
	//in micrometers
	const double EPS0=8.85e-18;

	//Lenght in micrometers
	const double cm=10000.;
	const double m =100*cm;
	const double um=1.;

	//Time in secs
	const double s = 1.;
	const double us= 1.e-6;
	const double ns= 1.e-9;

	const double coefficientElectron = 36*cm*cm/s;
	const double coefficientHole = 12*cm*cm/s;
	const double keV =1.;
	const double MeV=1000.;

	const double V=1;
	const double uA=1.;
	const double mA=1000.;
	const double A=1.e6;

	const double T=V*s/m/m;
	const double pi = TMath::Pi();

	#endif

/*-------------------------------------Constants.h end-----------------------------------------*/
using namespace dolfin;

class DriftTransport
{
private:
	JacoboniMobility _mu;
	Function * _d_f_grad;
	int _sign;
	int _diffusion;
	double _dt;
	double _temp;
	TRandom3 gRandom;


public:
	DriftTransport(char carrier_type, Function * d_f_grad, double givenT = 253., int difussion = 0, double dt = 300);
	DriftTransport();
	~DriftTransport();
	void operator() ( const std::array< double,2> &x , std::array< double,2> &dxdt , const double /* t */ );



};

#endif // CARRIERTRANSPORT_H
/*-------------------------------------CarrierTransport.h end-----------------------------------------*/

//#include <SMSDetector.h>
//#include <Constants.h>

//#include <Global.h>

//#include <vector>
#include <memory>

using namespace boost::numeric::odeint;
//using namespace dolfin;

/************************************CARRIER***********************************
 *
 * Used for instansciating every carrier read from the file.
 * Each carrier is treated independently, calculating its position, diffusion and drifting, then, its contribution to the global induced current is taken into account.
 *
 * Carriers have different properties like carrier_type (e or h), trapping time, etc...
 *
 */

class Carrier
{
	//It is not very clear why, but this class does not accept more variables.
	//The size is limited somewhere by someone, maybe the carrier vector list is touching the maximum space of the program's stack...
private:

	char _carrier_type;
	double _q; // charge
	double _dy;
	double _dx; //to calculate diffusion step
	double _gen_time; // instant of generation of the carrier
	std::array< double,2> _x; // carrier position array
	//std::array< double,2> _vectordx;
	std::array< double,2> _temp;
	std::array< double,2> _e_field; // electric field at the carrier position
	std::array< double,2> _w_field; // weighting field at the carrier positions
	double _e_field_mod;
	int _sign; // sign to describe if carrier moves in e field direction or opposite
	mutable std::mutex safeRead;
	SMSDetector * _detector;
	double _myTemp; // Temperature of the detector
	DriftTransport _drift;
	JacoboniMobility _mu;
	double _trapping_time;
	double diffDistance;
	bool _crossed;
	TRandom3 gRandom;

public:



	Carrier( char carrier_type, double q, double x_init, double y_init, SMSDetector * detector, double gen_time);
	Carrier(Carrier&& other); // Move declaration
	Carrier& operator = (Carrier&& other); // Move assignment
	Carrier(const Carrier& other); // Copy declaration
	Carrier& operator = (const Carrier& other); // Copy Assignment
	~Carrier();

	char get_carrier_type();

    inline double get_gen_time() const
    {
        return _gen_time;
    }
    
	//		double get_gen_time();
	//    std::array< double,2> get_e_field;
	//    std::array< double,2> get_w_field;
	//		double get_e_field_mod;
	//    int get_sign;
	//    SMSDetector *get_detector;
	//    double get_myTemp;
	//    DriftTransport get_drift;
	//    JacoboniMobility get_mu;
	//    double get_trapping_time;
	//
	std::array< double,2> get_x();
	double get_q();
	double get_diffDistance();
	bool crossed();

	//double get_diffH(){return diffH;};

	void calculateDiffusionStep(double dt);

    //double calculateAlpha(double efield); 
    std::valarray<double> calculateAlpha(double efield);    

	//void calculateDiffusionH(double dt);
	std::valarray<double> simulate_drift( double dt, double max_time);
	std::valarray<double> simulate_drift(double dt, double max_time, double x_init, double y_init, const std::string &scanType,
        std::vector<std::unique_ptr<Carrier>>& gen_carrier_list );    
};

#endif // CARRIER_H

/*-------------------------------------Carrier.h end-----------------------------------------*/
//#include <CarrierCollection.h>
/*------------------------------------CarrierCollection.h start------------------------------------------*/
#ifndef CARRIER_COLLECTION_H
#define CARRIER_COLLECTION_H

//#include "Carrier.h"
//#include <CarrierMobility.h>

//#include <string>
#include <sstream>
//#include <fstream>
//#include <vector>

//#include <QString>
//#include <TH2D.h>
//#include <TString.h>
//#include <TMath.h>
//#include <TRandom3.h>

#include <deque>
//#include <memory>

/*
 ***********************************CARRIER COLLECTION***********************************
 *
 *Carrier collection get carriers information from a file and store it in a vector of carriers.
 *This vector of carriers will be used by each of the threads in each of the steps of the TCT scan.
 *
 */

class CarrierCollection
{
private:
	std::deque<std::vector<std::unique_ptr<Carrier>>> _queue_carrier_list;
	SMSDetector * _detector;
	TRandom3 gRandom;
public:
	CarrierCollection(SMSDetector * detector);
	~CarrierCollection();
	double beamy = 0. , beamz = 0.; //Mean position of the injected carriers in detector plane (y,z)
	void add_carriers_from_file(const std::string& filename,const std::string& scanType, double depth);
	void simulate_drift( double dt, double max_time, double shift_x, double shift_y,  std::valarray<double> &curr_elec,     
			std::valarray<double> &curr_hole, int &totalCrosses, const std::string &scantype);

	void simulate_drift( double dt, double max_time, double shift_x, double shift_y,
                         std::valarray<double> &curr_elec, std::valarray<double> &curr_hole,
                         std::valarray<double> &curr_gen_elec, std::valarray<double> &curr_gen_hole, double max_mul_factor, 
                         int &totalCrosses, const std::string &scantype, std::string skip_event_loop);

	TH2D get_e_dist_histogram(int n_bins_x, int n_bins_y, TString hist_name = "e_dist", TString hist_title ="e_dist");
	TH2D get_e_dist_histogram(int n_bins_x, int n_bins_y, double shift_x, double shift_y, TString hist_name = "e_dist", TString hist_title ="e_dist");

	void record_carrier_gen_time(double max_time, int n_time_slice);
};


#endif // CARRIER_COLLECTION_H



/*-------------------------------------CarrierCollection.h end-----------------------------------------*/
//#include <Global.h>


using std::vector;

extern TH1D *H1DConvolution( TH1D *htct, Double_t Cend=0. , int tid=0, TString transferFun = "NO_TF") ;

class TRACSInterface
{

	//friend TTree * GetTree( ) ;
private:

	// Declaring external convolution function
	double pitch;
	double width;
	double depth;
	double temp;
	double trapping;
	double fluence;
	double C;
	double dt;
	double max_time;
	double vInit; //added v
	double deltaV;
	double vMax;
	double v_strips;
	double v_backplane;
	//double v_depletion;
	double deltaZ;
	double zInit;
	double zMax;
	double yInit;
	double yMax; //((2*nns)+1)*pitch,
	double deltaY; //added ^
	double vBias;
	double vDepletion;
	double zPos;
	double yPos;
	double tolerance;
	double chiFinal;
	double depletion_width;
	int diffusion;
	double fitNorm;
	//double gen_time;

    // For avalanche region
    std::string _set_avalanche_flag;
    std::array<double, 2>      _doping_peakheight;
    std::array<double, 2>      _doping_peakpos;    
    std::array<double, 2>      _doping_gauss_sigma;
    double      _max_multiplication_factor;
    
    std::array<std::array<double, 3>, 2> _doping_param;    // For two Gaussian sets

    std::string _skip_event_loop;

    std::string _simulation_polarity_inverse_flag;

    int _after_fitting_process;   // indicating whether the status is before or after fitting procedure. 
        
	int total_crosses;
	bool underDep;

	int nThreads;
	int nns;
	int n_cells_y;
	int n_cells_x;
	int n_tSteps;
	int waveLength; //added v
	int n_vSteps;
	int n_zSteps, n_zSteps1, n_zSteps2, n_zSteps_per_thread, n_zSteps_iter, n_balance;
	int n_ySteps;
	//int num_threads;


	uint n_par0;
	uint n_par1;
	uint n_par2;
	std::vector<uint> params = {0, 0, 0};
	int tcount;
	int count1, count2, count3;

	char bulk_type;
	char implant_type;

	std::vector<double> neff_param = {0};
	std::valarray<double> i_total;
	std::valarray<double> i_elec;
	std::valarray<double> i_hole;
    
	std::valarray<double> i_gen_elec;        
	std::valarray<double> i_gen_hole;
    
	std::valarray<double> i_shaped;
	std::valarray<double> i_temp;



	//double arrayTemp[300][300];

	vector<double> yVector;
	vector<double> zVector;
	vector<double> voltVector;

	std::vector<double>  z_shifts;
	vector<vector <double> >  z_shifts_array;

	//double z_shifts_array[10][10];
	std::vector<double>  z_shifts1, z_shifts2;
	std::vector<double>  y_shifts; // laser shift in X axis to center laser focus over read-out strip
	std::vector<double>  voltages;

	//vector of i_total
	//std::vector<double> vI_totals;
	//std::vector<double> valItotals;

	std::string carrierFile;
	std::string neffType;
	std::string scanType;

	//file naming
	std::string trap, start;
	// Convert relevant simulation numbers to string for fileNaming
	std::string dtime;
	std::string neighbrg_strips;
	std::string stepV;
	std::string stepZ;
	std::string stepY;
	std::string capacitance;
	//std::string z_step  = std::to_string((int) std::floor(deltaZ));
	std::string voltage;

	// filename for data analysis
	std::string hetct_conv_filename;
	std::string hetct_noconv_filename;
	std::string hetct_rc_filename;

	//TH1D i_ramo;
	TH1D *i_ramo;
	TH1D *i_rc;
	TH1D *i_conv;

	vector<vector <TH1D*> > i_rc_vector;


	//TH1D *hnoconv , *hconv;
	// Pointer to detector and carrier collection
	SMSDetector * detector;
	//SMSDetector * pDetector;
	CarrierCollection * carrierCollection;

	//vector of i_total
	//std::valarray<std::valarray <double> > vItotals;

	//Time variables
	UShort_t year, month, day, hour, min, sec;

public:
	std::valarray<std::valarray <double> > vSemiItotals;

	// Constructor
	TRACSInterface(std::string filename, const std::string& cFile, const std::vector<double>& z, const std::vector<double>& y, const std::vector<double>& v); // Reads values, initializes detector

	// Destructor
	~TRACSInterface();

	// Getters
	//TH1D GetItRamo();
	TH1D *GetItRamo();
	//TH1D *GetItRc();
	void GetItRc();
	TH1D *GetItConv();

	void GetItRc(std::valarray<double>& i_out, const std::valarray<double>& i_in);

	//Tree functions
	//friend TTree * GetTree(); //Returns the pointer to the TRACS simulated tree
	//friend void DumpToTree( TMeas em , TTree *tree );

	// Simulations
	void simulate_ramo_current();
	void calculate_fields();

	//Calculate time
	UShort_t GetYear();
	UShort_t GetMonth();
	UShort_t GetDay();
	UShort_t GetHour();
	UShort_t GetMinute();
	UShort_t GetSecond();
	double get_vDep();
	double get_depth();
	double get_fitNorm();
	double get_capacitance();
	double get_fluence();
	double get_vBias();
	double get_dt();
	//double get_genTime();
	std::vector<double> get_NeffParam();
	std::string get_neff_type();
	int GetnSteps();
	double GetTolerance();
	double GetchiFinal();
	int GettotalCrosses();
	double get_Itemp(int, int);



	//Loops
	void loop_on(int tid = 0); //MULTITHREADING
	void loop_on_topBottom();

	// Setters
	void set_Fit_Norm(std::vector<double> vectorFitTri);
	void set_FitParam(std::vector<double> newParam);
	void set_trappingTime(double newTrapTime);
	void set_zPos(double newZPos);
	void set_yPos(double newYPos);
	void set_vBias(double newVBias);
	void set_tcount(int tid = 0);
	void write_header(int tid = 0);
	void resize_array();
	void write_to_file(int tid = 0);
	void fields_hist_to_file(int, int);
	void currents_hist_to_file(int tid, int vpos, int nscan);
	void set_neffType(std::string newParametrization);
	void set_carrierFile(std::string newCarrFile);
	void set_vItotals(double);
	void resetAll();

	inline void set_fit_status(){ _after_fitting_process = 1; };
        
	//ROOT related
	//void DumpToTree( TMeas *em , TTree *tree ) ;
	void GetTree( TTree * tree ) ;
};

#endif // TRACSINTERFACE_H
/*-------------------------------------TRACSInterface.h end-----------------------------------------*/



//#include "../include/Threading.h"
//#include <Utilities.h>
//#include <TH1D.h>

//for H1DConvolution.C
#include <cstdarg>
#include "TKey.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAttLine.h"
#include "TLegend.h"
#include "THStack.h"



////////////////////////boarder between class (head file) and function code (.cpp)////////////////////////////////////////
//global var 

//Main variables of TRACS to store the induced current during the whole execution
std::valarray<std::valarray <double> > vItotals;
vector<TH1D*> i_rc_array, i_ramo_vector, i_conv_vector;
std::string global_TF;

//Define here the steering file you want to use. Store it in myApp folder.
std::string fnm="MyConfigTRACS";
//For mutex areas
std::mutex mtx;
int num_threads;
std::ofstream fileDiffDrift;


std::vector<TRACSInterface*> TRACSsim;
std::vector<std::thread> t;
std::vector<double> vector_zValues, vector_voltValues, vector_yValues;
std::string scanType;
double dTime;
double max_time;
double capacitance;
void spread_into_threads();
//global var end


//function needed list
void utilities::parse_config_file(std::string fileName, double &depth, double &width, double &pitch, int &nns, double &temp, double &trapping, double &fluence,
                                  int &nThreads, int &n_cells_x, int &n_cells_y, char &bulk_type, char &implant_type, std::string& skip_flag,
                                  std::string& set_avalanche_flag, std::array<double, 2>& doping_peakheight, std::array<double, 2>& doping_peakpos,  // for avalanche region
                                  std::array<double, 2>& doping_gauss_sigma, double& max_multiplication_factor,                                                             // for avalanche region
                                  int &waveLength, std::string &scanType, double &C, double &dt, double &max_time,
                                  double &v_init, double &deltaV, double &v_max, double &v_depletion, double &zInit, double &zMax, double &deltaZ, double &yInit, double &yMax, double &deltaY,
                                  std::vector<double> &neff_param, std::string &neffType, double &tolerance, double &chiFinal, int &diffusion, double &fitNorm, std::string& simulation_polarity_flag)
{
	// Creat map to hold all values as strings 
	std::map< std::string, std::string> valuesMap;
	std::string id, eq, val;
	std::stringstream converter;
	std::string tempString;

	std::ifstream configFile(fileName, std::ios_base::in);

	if (configFile.is_open())
	{

		std::string line;
		char comment = '#';
		char empty = '\0';
		char tab = '\t';
		while(std::getline(configFile, line))
		{
			char start = line[0];
			if (start == comment || start == empty || start == tab) continue;  // skip comments
			std::istringstream isstream(line);
			isstream >> id >> eq >> val;
			if (eq != "=") 
			{
				std::cout << "Error ecountered while reading '" << id << std::endl;
				break;
			}
			else
			{
				// Store value on map as map[variabel] = value
				valuesMap[id] = val;
			}
		}

		tempString = std::string("Depth");
		converter << valuesMap[tempString];
		converter >> depth;
		converter.str("");
		converter.clear();
		tempString = std::string("");

		tempString = std::string("Width");
		converter << valuesMap[tempString];
		converter >> width;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("Pitch");
		converter << valuesMap[tempString];
		converter >> pitch;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("NeighbouringStrips");
		converter << valuesMap[tempString];
		converter >> nns;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("Temperature");
		converter << valuesMap[tempString];
		converter >> temp;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("TrappingTime");
		converter << valuesMap[tempString];
		converter >> trapping;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("Fluence");
		converter << valuesMap[tempString];
		converter >> fluence;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("NumberOfThreads");
		converter << valuesMap[tempString];
		converter >> nThreads;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("CellsX");
		converter << valuesMap[tempString];
		converter >> n_cells_x;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("CellsY");
		converter << valuesMap[tempString];
		converter >> n_cells_y;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("Bulk");
		converter << valuesMap[tempString];
		converter >> bulk_type;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("Implant");
		converter << valuesMap[tempString];
		converter >> implant_type;
		converter.clear();
		converter.str("");
		tempString = std::string("");

        
		tempString = std::string("SkipEventLoop");
        skip_flag = valuesMap[tempString];
		tempString = std::string("");        

        // For avalanche region
		tempString = std::string("SetAvalanche");
		converter << valuesMap[tempString];
		converter >> set_avalanche_flag;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("DopingProfile_PeakHeight");
		converter << valuesMap[tempString];
		converter >> doping_peakheight[0];
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("DopingProfile_PeakPosition");
		converter << valuesMap[tempString];
		converter >> doping_peakpos[0];
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("DopingProfile_GaussSigma");
		converter << valuesMap[tempString];
		converter >> doping_gauss_sigma[0];
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("DopingProfile_PeakHeight2");
		converter << valuesMap[tempString];
		converter >> doping_peakheight[1];
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("DopingProfile_PeakPosition2");
		converter << valuesMap[tempString];
		converter >> doping_peakpos[1];
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("DopingProfile_GaussSigma2");
		converter << valuesMap[tempString];
		converter >> doping_gauss_sigma[1];
		converter.clear();
		converter.str("");
		tempString = std::string("");        

		tempString = std::string("MaxMultiplicationRatio");
		converter << valuesMap[tempString];
		converter >> max_multiplication_factor;
		converter.clear();
		converter.str("");
		tempString = std::string("");        
        
        //
        
		tempString = std::string("Lambda");
		converter << valuesMap[tempString];
		converter >> waveLength;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("Capacitance");
		converter << valuesMap[tempString];
		converter >> C;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("TimeStep");
		converter << valuesMap[tempString];
		converter >> dt;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("TotalTime");
		converter << valuesMap[tempString];
		converter >> max_time;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("InitialVoltage");
		converter << valuesMap[tempString];
		converter >> v_init;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("VoltageStep");
		converter << valuesMap[tempString];
		converter >> deltaV;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("MaxVoltage");
		converter << valuesMap[tempString];
		converter >> v_max;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("DepletionVoltage");
		converter << valuesMap[tempString];
		converter >> v_depletion;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("StepInZ");
		converter << valuesMap[tempString];
		converter >> deltaZ;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("InitialZ");
		converter << valuesMap[tempString];
		converter >> zInit;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("MaximumZ");
		converter << valuesMap[tempString];
		converter >> zMax;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("StepInY");
		converter << valuesMap[tempString];
		converter >> deltaY;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("InitialY");
		converter << valuesMap[tempString];
		converter >> yInit;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("MaximumY");
		converter << valuesMap[tempString];
		converter >> yMax;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("ScanType");
		scanType = valuesMap[tempString];
		tempString = std::string("");

		tempString = std::string("NeffParametrization");
		neffType = valuesMap[tempString];
		tempString = std::string("");

		tempString = std::string("y0");
		converter << valuesMap[tempString];
		converter >> neff_param[0];
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("y1");
		converter << valuesMap[tempString];
		converter >> neff_param[1];
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("y2");
		converter << valuesMap[tempString];
		converter >> neff_param[2];
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("y3");
		converter << valuesMap[tempString];
		converter >> neff_param[3];
		converter.clear();
		converter.str("");
		tempString = std::string("");


		tempString = std::string("z0");
		converter << valuesMap[tempString];
		converter >> neff_param[4];
		converter.clear();
		converter.str("");
		tempString = std::string("");
		//z0 = 0.0
		//neff_param[4] = 0.0;

		tempString = std::string("z1");
		converter << valuesMap[tempString];
		converter >> neff_param[5];
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("z2");
		converter << valuesMap[tempString];
		converter >> neff_param[6];
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("z3");
		converter << valuesMap[tempString];
		converter >> neff_param[7];
		converter.clear();
		converter.str("");
		tempString = std::string("");
		//z3 = depth
		//neff_param[7] = depth;

		tempString = std::string("tolerance");
		converter << valuesMap[tempString];
		converter >> tolerance;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("chiFinal");
		converter << valuesMap[tempString];
		converter >> chiFinal;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("diffusion");
		converter << valuesMap[tempString];
		converter >> diffusion;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("normalizator");
		converter << valuesMap[tempString];
		converter >> fitNorm;
		converter.clear();
		converter.str("");
		tempString = std::string("");

        tempString = std::string("SimulationPolarityInversion");
        simulation_polarity_flag = valuesMap[tempString];
		tempString = std::string("");
        
		configFile.close();
	}
	else
	{
		std::cout << "Error opening the file. Does the file " << fileName << " exist?" << std::endl;
	}
}

void utilities::parse_config_file(std::string fileName, std::string &scanType, double &v_init, double &deltaV, double &v_max, double &v_depletion,
		double &zInit, double &zMax, double &deltaZ, double &yInit, double &yMax, double &deltaY, double &dt, double &max_time, double &capacitance, std::string &transferFun)
{
	// Creat map to hold all values as strings
	std::map< std::string, std::string> valuesMap;
	std::string id, eq, val;
	std::stringstream converter;
	std::string tempString;

	std::ifstream configFile(fileName, std::ios_base::in);

	if (configFile.is_open())
	{

		std::string line;
		char comment = '#';
		char empty = '\0';
		char tab = '\t';
		while(std::getline(configFile, line))
		{
			char start = line[0];
			if (start == comment || start == empty || start == tab) continue;  // skip comments
			std::istringstream isstream(line);
			isstream >> id >> eq >> val;
			if (eq != "=")
			{
				std::cout << "Error ecountered while reading '" << id << std::endl;
				break;
			}
			else
			{
				// Store value on map as map[variabel] = value
				valuesMap[id] = val;
			}
		}


		tempString = std::string("InitialVoltage");
		converter << valuesMap[tempString];
		converter >> v_init;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("VoltageStep");
		converter << valuesMap[tempString];
		converter >> deltaV;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("MaxVoltage");
		converter << valuesMap[tempString];
		converter >> v_max;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("DepletionVoltage");
		converter << valuesMap[tempString];
		converter >> v_depletion;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("StepInZ");
		converter << valuesMap[tempString];
		converter >> deltaZ;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("InitialZ");
		converter << valuesMap[tempString];
		converter >> zInit;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("MaximumZ");
		converter << valuesMap[tempString];
		converter >> zMax;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("StepInY");
		converter << valuesMap[tempString];
		converter >> deltaY;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("InitialY");
		converter << valuesMap[tempString];
		converter >> yInit;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("MaximumY");
		converter << valuesMap[tempString];
		converter >> yMax;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("ScanType");
		scanType = valuesMap[tempString];
		tempString = std::string("");

		tempString = std::string("TimeStep");
		converter << valuesMap[tempString];
		converter >> dt;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("TotalTime");
		converter << valuesMap[tempString];
		converter >> max_time;
		converter.clear();
		converter.str("");
		tempString = std::string("");

		tempString = std::string("Capacitance");
		converter << valuesMap[tempString];
		converter >> capacitance;
		converter.clear();
		converter.str("");
		tempString = std::string("");

	    tempString = std::string("TransferFunction");
		transferFun = valuesMap[tempString];
		tempString = std::string("");



		configFile.close();
	}
	else
	{
		std::cout << "Error opening the file. Does the file " << fileName << " exist?" << std::endl;
	}
}
void utilities::parse_config_file(std::string fileName, std::string &carrierFile){

	// Creat map to hold all values as strings
	std::map< std::string, std::string> valuesMap;
	std::string id, eq, val;
	std::stringstream converter;
	std::string tempString;

	std::ifstream configFile(fileName, std::ios_base::in);

	if (configFile.is_open())
	{

		std::string line;
		char comment = '#';
		char empty = '\0';
		char tab = '\t';
		while(std::getline(configFile, line))
		{
			char start = line[0];
			if (start == comment || start == empty || start == tab) continue;  // skip comments
			std::istringstream isstream(line);
			isstream >> id >> eq >> val;
			if (eq != "=")
			{
				std::cout << "Error ecountered while reading '" << id << std::endl;
				break;
			}
			else
			{
				// Store value on map as map[variabel] = value
				valuesMap[id] = val;
			}
		}


		tempString = std::string("CarrierFile");
		carrierFile = valuesMap[tempString];
		tempString = std::string("");

		configFile.close();
	}
	else
	{
		std::cout << "Error opening the file. Does the file " << fileName << " exist?" << std::endl;
	}

}


void call_from_thread(int tid, std::string& carrier_name_FileThr, const std::vector<double>& zVector, const std::vector<double>& yVector, const std::vector<double>& voltVector)
{
	// every thread instantiates a new TRACSInterface object

	TRACSsim[tid] = new TRACSInterface(fnm, carrier_name_FileThr, zVector, yVector, voltVector);///////////////////////////////////////<<<<<<<<<<<working on 2
	TRACSsim[tid]->set_tcount(tid);
	if(tid==0)TRACSsim[tid]->write_header(tid);
	TRACSsim[tid]->loop_on(tid);
}

std::mutex mtx2;
TRACSInterface::TRACSInterface(std::string filename, const std::string& carrFile, const std::vector<double>& zVector, const std::vector<double>& yVector,
                               const std::vector<double>& voltVector):zVector(zVector), yVector(yVector), voltVector(voltVector), carrierFile(carrFile), _after_fitting_process(0)
{
	neff_param = std::vector<double>(8,0);
	total_crosses = 0;

    utilities::parse_config_file(filename, depth, width,  pitch, nns, temp, trapping, fluence, nThreads, n_cells_x, n_cells_y, bulk_type, implant_type, _skip_event_loop,
                                 _set_avalanche_flag, _doping_peakheight, _doping_peakpos, _doping_gauss_sigma, _max_multiplication_factor,   // new parameters for avalanche regions
                                 waveLength, scanType, C, dt, max_time, vInit, deltaV, vMax, vDepletion, zInit, zMax, deltaZ, yInit, yMax, deltaY, neff_param, neffType,
                                 tolerance, chiFinal, diffusion, fitNorm, _simulation_polarity_inverse_flag);
    
    // Pack into a structure .  In future, it can be a Structure or Class.
    _doping_param[0][0] = _doping_peakheight[0];
    _doping_param[0][1] = _doping_peakpos[0];
    _doping_param[0][2] = _doping_gauss_sigma[0];

    _doping_param[1][0] = _doping_peakheight[1];
    _doping_param[1][1] = _doping_peakpos[1];
    _doping_param[1][2] = _doping_gauss_sigma[1];

    
	if (fluence == 0) // if no fluence -> no trapping
	{
		//Trapping configuration
		trapping = std::numeric_limits<double>::max();
		trap = "NOtrapping";
		start = "NOirrad";
	}
	else{
		//trapping configuration
		trap = std::to_string((int) std::floor(1.e9*trapping));
		start = "irrad";

	}

	//Following variables are used for writing data in the output file.
	//Convert relevant simulation numbers to string for fileNaming
	dtime = std::to_string((int) std::floor(dt*1.e12));
	neighbrg_strips = std::to_string(nns);
	stepV = std::to_string((int) std::floor(deltaV));
	stepZ = std::to_string((int) std::floor(deltaZ));
	stepY = std::to_string((int) std::floor(deltaY));
	capacitance = std::to_string((int) std::floor(C*1.e12));
	voltage = std::to_string((int) std::floor(vInit));


	//Dolfin instruction por mesh boundary extrapolation
	parameters["allow_extrapolation"] = true;

	std::unique_lock<std::mutex> guard(mtx2);
	detector = new SMSDetector(pitch, width, depth, nns, bulk_type, implant_type, _set_avalanche_flag, _doping_param,
                               n_cells_x, n_cells_y, temp, trapping, fluence, neff_param, neffType, diffusion, dt);
	guard.unlock();
	
	carrierCollection = new CarrierCollection(detector);
	carrierCollection->add_carriers_from_file(carrierFile, scanType, depth);

	n_tSteps = (int) std::floor(max_time / dt);
	//currents vectors used to store temporal and final values.
	i_elec.resize((size_t) n_tSteps);
	i_hole.resize ((size_t) n_tSteps);

	i_gen_elec.resize((size_t) n_tSteps);
	i_gen_hole.resize ((size_t) n_tSteps);
    
	i_total.resize((size_t) n_tSteps);
	i_shaped.resize((size_t) n_tSteps);
	i_temp.resize((size_t) n_tSteps);

	for (int k = 1 ; k < n_tSteps; k++ ){
		i_temp[k] = 0 ;
	}


	if (scanType == "edge"){
		vSemiItotals.resize(zVector.size() * voltVector.size());

		for (int i = 0; i < (zVector.size() * voltVector.size()) ; i++){

			vSemiItotals[i].resize(n_tSteps);

		}
	}
	if (scanType == "top" || scanType == "bottom"){
		vSemiItotals.resize(yVector.size() * voltVector.size());

		for (int i = 0; i < (yVector.size() * voltVector.size()) ; i++){

			vSemiItotals[i].resize(n_tSteps);

		}
	}


	for (int i = 0 ; i < vSemiItotals.size() ; i++){
		for (int j = 0 ; j < vSemiItotals[i].size() ; j++)
			vSemiItotals[i][j] = 0;
	}


	vBias = vInit;
	set_tcount(0);

	i_ramo  = NULL;
	i_rc    = NULL;
	i_conv  = NULL;

}

/**
 *
 * @param pitch
 * @param width
 * @param depth
 * @param nns
 * @param bulk_type
 * @param implant_type
 * @param n_cells_x
 * @param n_cells_y
 * @param tempK
 * @param trapping
 * @param fluence
 * @param neff_param
 * @param neff_type
 * @param diffusion
 * @param dt
 */
SMSDetector::SMSDetector(double pitch, double width, double depth, int nns, char bulk_type, char implant_type, std::string avalanche_flag, std::array<std::array<double, 3>, 2> doping_param,
                         int n_cells_x, int n_cells_y, double tempK, double trapping,
                         double fluence, std::vector<double> neff_param, std::string neff_type, int diffusion, double dt) :

				_pitch(pitch), //Distance between implants
				_width(width), //Size of the implant
				_depth(depth), //Vertical size of the pad (typically 300microns)
				_tempK(tempK), // Temperature of the detector
				_trapping_time(trapping), // Trapping constant simulates radiation-induced defects (traps)
				_fluence(fluence), // Irradiation fluence in neutron equivalent (neq)
				_nns(nns), // Number of neighbouring strips
				_bulk_type(bulk_type), //Dopant type of the silicon (p/n)
				_implant_type(implant_type), //Dopant type of the implant, normally opposite of the bulk (n/p)
				_neff_type(neff_type), // Select aproach to parametrize Neff (irradiation only)
				_neff_param(neff_param), // Parametrized description of the Space Charge distribution
				_x_min(0.0), // Starting horizontal position for carrier generation (hereafter CG)
				_x_max(_pitch * (2*_nns+1)), // Endingvertical positio for CG
				_y_min(0.0), // Starting vertical position for CG (microns)
				_y_max(_depth), // Ending vertical position for CG (microns)
				_vdep(0),
				_v_bias(0),
				_v_backplane(0),
				_v_strips(0),
				_f_poisson(0),
				_diffusion(diffusion),
				_depletion_width(0),
				_depleted(false),
				_dt(dt),
                _avalanche_flag(avalanche_flag),
                _doping_param(doping_param),
				// Mesh properties
				_n_cells_x(n_cells_x),
				_n_cells_y(n_cells_y),
#if DOLFIN_VERSION_MINOR>=6
				_mesh(Point(_x_min,_y_min),Point(_x_max,_y_max), _n_cells_x, _n_cells_y),///////////////////////////////<<<<<<<<<working on 3
#else
				_mesh(_x_min,_y_min,_x_max,_y_max, _n_cells_x, _n_cells_y),
#endif
				_periodic_boundary(_x_min, _x_max, _depth),
				_central_strip(_pitch, _width, _nns),
				_neighbour_strips(_pitch, _width, _nns),
				_backplane(_x_min, _x_max, _depth),

				// Functions & variables to solve the PDE
				_V_p(_mesh, _periodic_boundary),
				_a_p(_V_p, _V_p),
				_L_p(_V_p),
				_V_g(_mesh),
				_a_g(_V_g, _V_g),
				_L_g(_V_g),
				_w_u(_V_p),
				_d_u(_V_p),
				_w_f_grad(_V_g), // Weighting field
				_d_f_grad(_V_g)  // Drifting field
{
}


/*
 * Getter method for the temperature of the diode
 */
double  SMSDetector::get_temperature()
{
	return _tempK;
}

/*
 * Getter for the drift field
 */
Function * SMSDetector::get_d_f_grad()
{
	return &_d_f_grad;
}

int SMSDetector::diffusionON()
{
	return _diffusion;
}

double SMSDetector::get_dt(){

	return _dt;
}

/*
 * Getter for the trapping factor
 */
double SMSDetector::get_trapping_time()
{
	return _trapping_time;
}

/*
 * Getter for the fluence
 */
double SMSDetector::get_fluence()
{
	return _fluence;
}


/*
 * Getter for the bulk type
 */
char SMSDetector::get_bulk_type()
{
	return _bulk_type;
}


/*
 * Getter for the implant type
 */
char SMSDetector::get_implant_type()
{
	return _implant_type;
}


/*
 * Getter for the number of neighbouring strips
 */
int SMSDetector::get_nns()
{
	return _nns;
}


/*
 * Getter for the pitch
 */
double SMSDetector::get_pitch()
{
	return _pitch;
}


/*
 * Getter for the width
 */
double SMSDetector::get_width()
{
	return _width;
}


/*
 * Getter for the pitch
 */
double SMSDetector::get_depth()
{
	return _depth;
}


/*
 * Getter for the depletion voltage
 */
double SMSDetector::get_vdep()
{
	return _vdep;
}



CarrierCollection::CarrierCollection(SMSDetector * detector) :
						_detector(detector)
{

}

/*
 * Adds carriers from the file. The name of the file must be set in the steering file, in the field CarrierFile, more instructions can be found there.
 * The calculation of average positions is neccesary for fitting purposes, DumpToTree function.
 */
/**
 *
 * @param filename
 */
void CarrierCollection::add_carriers_from_file(const std::string& filename, const std::string& scanType, double depth)
{
	double extra_y;
	std::vector<std::unique_ptr<Carrier>>_carrier_list_sngl;
    
	// get char representation and make ifstream
	//char * char_fn = filename.toLocal8Bit().data();
	std::ifstream infile(filename);
	bool once = true;
	char carrier_type;
	double q, x_init, y_init, gen_time;

	// process line by line
	std::string line;

    
	//Preprocessing for fitting bottom scan_type: Fixing mismatch between detector depth and y_init.
	//Outside the loop for performance purpose
	if (scanType == "bottom"){
		std::getline(infile, line);
		std::istringstream iss(line);
		if (!(iss >> carrier_type >> q >> x_init >> y_init >> gen_time)) {
			//std::cout << "Error while reading file" << std::endl;
		}
		//Extra in micrometers to shift y_init position
		extra_y = depth - y_init;

		//Calculate average beam position
		beamy += x_init;
		beamz += y_init + extra_y;

        //std::cout << "Charge is registered at x=" << x_init << ", y = " << y_init + extra_y << ", gen_time = " << gen_time << std::endl;
        std::unique_ptr<Carrier> ptr( new Carrier(carrier_type, q, x_init, y_init + extra_y, _detector, gen_time) );
        _carrier_list_sngl.push_back( std::move(ptr) );
        
		if ( _carrier_list_sngl.size()!=0 ) {
			beamy = beamy / _carrier_list_sngl.size();
			beamz = beamz / _carrier_list_sngl.size();
		}
	}
	while (std::getline(infile, line))
	{

		std::istringstream iss(line);
		if (!(iss >> carrier_type >> q >> x_init >> y_init >> gen_time)) { 
			//std::cout << "Error while reading file" << std::endl;
			break;
		} 
		extra_y = depth;
		//Calculate average beam position
		beamy += x_init;
		beamz += y_init + extra_y;

        std::unique_ptr<Carrier> ptr( new Carrier(carrier_type, q, x_init, y_init + extra_y, _detector, gen_time) );
        _carrier_list_sngl.push_back( std::move(ptr) );        
	}
	if ( _carrier_list_sngl.size()!=0 ) {
		beamy = beamy / _carrier_list_sngl.size();
		beamz = beamz / _carrier_list_sngl.size();
	}

    // Push the default carrier list to the queue
    _queue_carrier_list.push_back( std::move(_carrier_list_sngl) );
    
}



/*
 * Sets a number (current thread).
 * Used to index the different output files 
 *
 */
/**
 *
 * @param tid
 */
void TRACSInterface::set_tcount(int tid)
{
	tcount = tid;
}


/*
 *
 * Write to file header. The input int is used to label files (multithreading)!
 *
 *
 *
 */
/**
 *
 * @param tid
 */
void TRACSInterface::write_header(int tid)
{
	// Convert Z to milimeters
	std::vector<double> aux_zsh(zVector.size());
	aux_zsh = zVector;
	std::transform(aux_zsh.begin(), aux_zsh.end(), aux_zsh.begin(), std::bind1st(std::multiplies<double>(),(1./1000.)));

	// Convert Z to milimeters
	std::vector<double> aux_ysh(yVector.size());
	aux_ysh = yVector;
	std::transform(aux_ysh.begin(), aux_ysh.end(), aux_ysh.begin(), std::bind1st(std::multiplies<double>(),(1./1000.)));
	hetct_rc_filename = start+"_dt"+dtime+"ps_"+capacitance+"pF_t"+trap+"ns_dz"+stepZ+"um_dy"+stepY+"dV"+stepV+"V_"+neighbrg_strips+"nns_"+scanType+"_"+std::to_string(tcount)+"_rc.hetct";
	// write header for data analysis
	utilities::write_to_hetct_header(hetct_rc_filename, detector, C, dt, aux_ysh, aux_zsh, waveLength, scanType, carrierFile, voltVector);

}


/*
 * Constructor for Carrier.cpp that sets and stores the values given in their respective places.
 *
 *
 *
 */
/**
 * @param carrier_type
 * @param q
 * @param x_init
 * @param y_init
 * @param detector
 * @param gen_time
 */
Carrier::Carrier( char carrier_type, double q,  double x_init, double y_init , SMSDetector * detector, double gen_time = 1.e-9):

		_carrier_type(carrier_type), // Charge carrier(CC)  type. Typically  electron/positron
		_q(q), //Charge in electron units. Always positive.
		_dy(0.),
		_dx(0.),
		_gen_time(gen_time), // Instant of CC generation
		_e_field_mod(0.),
		_detector(detector), // Detector type and characteristics
		//	_electricField(_detector->get_d_f_grad(),
		//	_weightingField(_detector->get_w_f_grad(),
		_myTemp(_detector->get_temperature()), // Temperature of the diode
		_drift(_carrier_type, detector->get_d_f_grad(), _myTemp, _detector->diffusionON(), _detector->get_dt()), // Carrier Transport object
		_mu(_carrier_type, _myTemp),// Mobility of the CC
		_trapping_time(_detector->get_trapping_time()),
		diffDistance(0.),
		_crossed(false)

{

	_x[0] = x_init; // Starting horizontal position
	_x[1] = y_init; // Starting vertical position

	if (_carrier_type == 'e')
	{ // If electron-like
		_sign = -1; // Negative charge
	}
	else
	{ // it's hole-like
		_sign = 1; // Positive charge
	}
}


//SLIGHTLY MODIFIED TO WORK WITH SMSDetector *detector
// function to write results to file (in rows)
// overloaded (now from TH1D)
void utilities::write_to_hetct_header(std::string filename, SMSDetector * detector, double C, double dt,std::vector<double> y_shifts, std::vector<double> z_shifts, double landa, std::string type, std::string carriers_file, std::vector<double> voltages)
{
	// Initialize stream for outputting to file
	std::ofstream header;  

	// Derived quantities
	int nZ = z_shifts.size();
	int nY = y_shifts.size();
	int nV = voltages.size();
	int nScans = nZ * nY * nV;
	double deltaZ = 0;
	double deltaV = 0;	
	double deltaY = 0;	
	if (nZ > 1)	deltaZ = std::abs((z_shifts.back()-z_shifts.front())/(nZ-1));
	if (nV > 1) deltaV = std::abs((voltages.back()-voltages.front())/(nV-1));
	if (nY > 1) deltaY = std::abs((y_shifts.back()-y_shifts.front())/(nY-1));
	double temp = detector->get_temperature() - 273;

	// Convert relevant quantities to string for outputting
	std::string s_date = "2015-06-15_19-25";
	std::string s_zVector = utilities::vector_to_string(z_shifts);
	std::string s_vVector = utilities::vector_to_string(voltages);
	std::string s_yVector = utilities::vector_to_string(y_shifts);

	// Open file
	header.open (filename);  

	// Check the file was open
	if (header.is_open())
	{ 
		// Store header string 
		header << "================\n";
		header << "SSD simulation file\n";
		header << "version: 1.0\n";
		header << "================\n";
		header << "scanType: " << "TRACS\n";
		header << "startTime: " << s_date << "\n" ;
		header << "landaLaser: " << landa << "\n";
		header << "illumDirect: " << type << "\n";
		header << "ampGain: 0\n";
		header << "fluence: " << detector->get_fluence() << "\n";
		header << "annealing: 0\n";
		header << "numOfScans: " << nScans << "\n";
		header << "temperatureMin: " << temp  << "\n";
		header << "temperatureMax: " << temp  << "\n";
		header << "deltaTemperature: " << 0. << "\n";
		header << "nTemperature: 1\n";
		header << "temperatureVector: " << temp  << "\n";
		header << "voltageMin: " << voltages.front() << "\n";
		header << "voltageMax: " << voltages.back() << "\n";
		header << "deltaVoltage: " << deltaV << "\n";
		header << "nVoltage: " << nV << "\n";
		header << "voltageVector: " << s_vVector << "\n";
		header << "pulsePowerIntensityMin: 60.000\n";
		header << "pulsePowerIntensityMax: 60.000\n";
		header << "deltaPulsePowerIntensity: 0\n";
		header << "nPulsePowerIntensity: 1\n";
		header << "pulsePowerIntensityVector: 60.000 \n";
		header << "pulseWidthMin: 0.000\n";
		header << "pulseWidthMax: 0.000\n";
		header << "deltaPulseWidth: 0\n";
		header << "nPulseWidth: 1\n";
		header << "pulseWidthVector: 0.000 \n";
		header << "xMin: 0.000\n";
		header << "xMax: 0.000\n";
		header << "deltaX: 0\n";
		header << "nX: 1\n";
		// X o Y tengo que aclararme porque no lo entiendo!!!
		header << "yMin: " << y_shifts.front() << "\n";
		header << "yMax: " << y_shifts.back() << "\n";
		header << "deltaY: " << deltaY << "\n";
		header << "nY: " << nY << "\n";
		header << "yVector: " << s_yVector << "\n";
		header << "zMin: " << z_shifts.front() << "\n";
		header << "zMax: " << z_shifts.back() << "\n";
		header << "deltaZ: " << deltaZ << "\n";
		header << "nZ: " << nZ << "\n";
		header << "zVector: " << s_zVector << "\n";
		header << "At: " << std::fixed << std::setprecision(15) << dt << "\n";
		header << "Capacitance[F]: " <<  std::fixed << std::setprecision(15) << C << "\n";
		header << "Bulk: " << detector->get_bulk_type() << "\n";
		header << "Implant: " << detector->get_implant_type() << "\n";
		header << "NStrips: " << detector->get_nns() << "\n";
		header << "Pitch: " << std::setprecision(0) << detector->get_pitch() << "\n";
		header << "Width: " << std::setprecision(0) << detector->get_width() << "\n";
		header << "Depth: " << std::setprecision(0) << detector->get_depth() << "\n";
		header << "Vdep: " << detector->get_vdep() << "\n";
		header << "Carriers File: " << carriers_file << "\n";
		header << "================\n";
		header <<	 "Nt T[C] Vset[V] x[mm] y[mm] z[mm] I(t)[A]\n";
		header <<	 "================\n";

		header.close();
	}
	else // Error output
	{
		std::cout << "File could not be created"<<std::endl;
	}
}


// Utility to convert large arrays/vecto into strings
std::string utilities::vector_to_string(std::vector<double> input_list)
{
	// Use stringstream to convert doubles to string
	unsigned long int points = input_list.size();
	std::stringstream converter; 

	for (unsigned int i = 0; i < points; i++) 
	{
		converter << input_list[i] << " ";
	}

	return converter.str();
}


/**
 *
 * @param tid
 */
void TRACSInterface::loop_on(int tid)
{

	int index_total = 0;
	int i,j;
	i = 0; 
	j = 0;

	if (scanType == "edge"){
		//Voltage scan
		for (int index_volt = 0; index_volt < voltVector.size() ; index_volt++){

			detector->set_voltages(voltVector[index_volt], vDepletion);
			std::unique_lock<std::mutex> guard(mtx2);
			detector->solve_w_u();
			detector->solve_d_u();
			detector->solve_w_f_grad();
			guard.unlock();
			detector->solve_d_f_grad();
			detector->get_mesh()->bounding_box_tree();

			for (int index_zscan = 0; index_zscan < zVector.size(); index_zscan++){

				//simulate_ramo_current();
                carrierCollection->simulate_drift( dt, max_time, yInit, zVector[index_zscan], i_elec, i_hole, i_gen_elec, i_gen_hole, _max_multiplication_factor, total_crosses, scanType,
                                                   _skip_event_loop); 
                
				//i_total = i_elec + i_hole;
                i_total = i_elec + i_hole + i_gen_elec + i_gen_hole;

                if( _after_fitting_process )
                {
                    // Reflect the scaling factor from the fitting. 
                    i_total *= fitNorm;
                    i_elec *= fitNorm;
                    i_hole *= fitNorm;
                    i_gen_elec *= fitNorm;
                    i_gen_hole *= fitNorm;
                    
                    // Temporally, reflect the polarity setting for the simulation.
                    if( _simulation_polarity_inverse_flag == "yes" )
                    {
                        i_total *= -1.0;
                        i_elec *= -1.0;
                        i_hole *= -1.0;
                        i_gen_elec *= -1.0;
                        i_gen_hole *= -1.0;                    
                    }
                }
                
				GetItRc();
				vSemiItotals[index_total] = i_shaped;

                // Save current distribution 
                currents_hist_to_file(tid, index_volt, index_zscan);
                
				index_total++;
				i_total = 0 ; i_elec = 0; i_hole = 0; i_shaped = 0; i_temp = 0;
                i_gen_elec = 0; i_gen_hole = 0;
			}

			if (tid == 0) fields_hist_to_file(tid, index_volt);
		}


	}

	if (scanType == "top" || scanType == "bottom"){
		//Voltage scan
		for (int index_volt = 0; index_volt < voltVector.size() ; index_volt++){

			detector->set_voltages(voltVector[index_volt], vDepletion);
			std::unique_lock<std::mutex> guard(mtx2);
			detector->solve_w_u();
			detector->solve_d_u();
			detector->solve_w_f_grad();
			guard.unlock();
			detector->solve_d_f_grad();
			detector->get_mesh()->bounding_box_tree();

			for (int index_yscan = 0; index_yscan < yVector.size(); index_yscan++){

				carrierCollection->simulate_drift( dt, max_time, yVector[index_yscan], zInit, i_elec, i_hole, i_gen_elec, i_gen_hole, _max_multiplication_factor, total_crosses, scanType,
                    _skip_event_loop);
                
				//i_total = i_elec + i_hole;
                i_total = i_elec + i_hole + i_gen_elec + i_gen_hole;

                if( _after_fitting_process )
                {
                    // Reflect the scaling factor from the fitting. Default value of "fitNorm" is 1.
                    i_total *= fitNorm;
                    i_elec *= fitNorm;
                    i_hole *= fitNorm;
                    i_gen_elec *= fitNorm;
                    i_gen_hole *= fitNorm;
                    
                    // Temporally, reflect the polarity setting for the simulation.
                    if( _simulation_polarity_inverse_flag == "yes" )
                    {
                        i_total *= -1.0;
                        i_elec *= -1.0;
                        i_hole *= -1.0;
                        i_gen_elec *= -1.0;
                        i_gen_hole *= -1.0;                    
                    }
                }
                
				GetItRc();
				vSemiItotals[index_total] = i_shaped;

                // Save current distribution 
                currents_hist_to_file(tid, index_volt, index_yscan);
                                
				index_total++;
				i_total = 0 ; i_elec = 0; i_hole = 0; i_shaped = 0; i_temp = 0;
                i_gen_elec = 0; i_gen_hole = 0;
			}

			if (tid == 0)fields_hist_to_file(tid, index_volt);
		}

	}

}

/*
 **
 ** Set right polarity for detector depending on the type (p-on-n/n-on-p)
 */
/**
 *
 * @param v_bias
 * @param v_depletion
 */
void SMSDetector::set_voltages(double v_bias, double v_depletion)
{
	_vdep = v_depletion; // Store depletion voltage
	_v_bias = v_bias;
	_depleted = (_v_bias >= _vdep) ? true : false;
	_v_strips = (_implant_type == 'n') ? v_bias : 0.0;
	_v_backplane = (_implant_type == 'p') ? v_bias : 0.0;

	// neff defined in F/microns
	//original only taking depth of the detector (case not depleted is ignored): _f_poisson = ((_bulk_type== 'p') ? +1.0 : -1.0)*(-2.0*v_depletion)/(_depth*_depth);

	if (_fluence == 0)
	{

		if (_depleted)
		{//Detector is depleted
			_depletion_width = _depth;
			_f_poisson = ((_bulk_type== 'p') ? +1.0 : -1.0)*(-2.0*_vdep)/(_depth*_depth);
			std::cout << "fp depleted: " << _f_poisson << std::endl;
			//std::cout << "Depletion Width: " << _depletion_width << std::endl;

		}

		else
		{//Detector is not depleted
			_depletion_width = _depth * sqrt((abs(_v_strips-_v_backplane))/_vdep);
			_f_poisson = ((_bulk_type== 'p') ? +1.0 : -1.0)*(-2.0*_v_bias)/(_depletion_width * _depletion_width);
			std::cout << "fp NO depleted: " << _f_poisson << std::endl;
			std::cout << "Depletion Width: " << _depletion_width << std::endl;

		}
	}
	//if fluence > 0, parameters are taken in the constructor inicialization list, coming from the config file
	if (_fluence > 0)
	{
		_depletion_width = _neff_param[7] - _neff_param[4];
	}

}

/*
 * Method for solving the weighting potential using Laplace triangles 
 */
void SMSDetector::solve_w_u()
//Weighting Potential
{
	// Solving Laplace equation f = 0
	std::vector<const DirichletBC*> bcs;
	Constant f(0.0);
	_L_p.f = f;

	// Set BC values
	Constant central_strip_V(1.0);
	Constant neighbour_strip_V(0.0);
	Constant backplane_V(0.0);

	BackPlaneBoundaryWP backplaneUnderDep(_x_min, _x_max, _depletion_width);

	DirichletBC central_strip_BC(_V_p, central_strip_V, _central_strip);
	DirichletBC neighbour_strip_BC(_V_p, neighbour_strip_V, _neighbour_strips);
	DirichletBC backplane_BC(_V_p, backplane_V, backplaneUnderDep);

	bcs.push_back(&central_strip_BC);
	bcs.push_back(&neighbour_strip_BC);
	bcs.push_back(&backplane_BC);

	solve(_a_p == _L_p , _w_u, bcs);

}

/*
 * Method for solving the d_u (drifting potential) using Poisson's equation
 */
void SMSDetector::solve_d_u()
{
	Constant fpois(_f_poisson);

	std::vector<const DirichletBC*> bcs;
	Source f;


	if (_fluence == 0) //NO irrad but YES depleted. fpoisson, charge distribution, is a constant during the whole detector
	{
        if ( _avalanche_flag == "yes" )
        {
            std::cout << "Avalanche layer is inserted " << std::endl;

            f.set_NeffApproach("AvalancheMode");
            f.set_avalanche_doping_param( _doping_param );
            f.set_bulk_doping_param( _f_poisson );

            // Save Neff distribution into a ROOT file. Now, it is only 1D histo.
            f.save_Neff_dist(_y_min, _y_max);
            
            _L_p.f = f;
        }
        else 
        {
            std::cout << "Avalanche layer is not set " << std::endl;
            
            // Setting for "save_Neff_dist"  function . last parameter (Gaussian sigma) should be >0 
            std::array<std::array<double, 3>, 2> doping_param_null{
                {  {0.0, 0.0, 0.001},
                    {0.0, 0.0, 0.001} }
            };
            
            f.set_avalanche_doping_param( doping_param_null );
            f.set_bulk_doping_param( _f_poisson );            

            // Save Neff distribution into a ROOT file. Now, it is only 1D histo.
            f.save_Neff_dist(_y_min, _y_max);            
            
            _L_p.f = fpois;   // original in TRACS V1.0
        }
        //std::cout << "NeffApproch = " << f.get_NeffApproach() << std::endl;
       
		_trapping_time = std::numeric_limits<double>::max();
	}

	else
		//If YES Irrad OR NO depleted, charge distribution is not a constant. Parameters from steering file of from above if non-depleted.
	{

		f.set_NeffApproach(_neff_type);
		f.set_y0(_neff_param[0]);
		f.set_y1(_neff_param[1]);
		f.set_y2(_neff_param[2]);
		f.set_y3(_neff_param[3]);
		f.set_z0(_neff_param[4]);
		f.set_z1(_neff_param[5]);
		f.set_z2(_neff_param[6]);
		f.set_z3(_neff_param[7]);
		_L_p.f = f;
	}


	// Set BC values
	Constant central_strip_V(_v_strips);
	Constant neighbour_strip_V(_v_strips);
	Constant backplane_V(_v_backplane);


	BackPlaneBoundaryWP backplaneUnderDep(_x_min, _x_max, _depletion_width);

	DirichletBC central_strip_BC(_V_p, central_strip_V, _central_strip);
	DirichletBC neighbour_strip_BC(_V_p, neighbour_strip_V, _neighbour_strips);
	DirichletBC backplane_BC(_V_p, backplane_V, backplaneUnderDep);

	bcs.push_back(&central_strip_BC);
	bcs.push_back(&neighbour_strip_BC);
	bcs.push_back(&backplane_BC);

	solve(_a_p == _L_p , _d_u, bcs);

}


// Static member
const double Source::_elementary_charge = 1.60217662e-19;    // [ C ]
const double Source::_vacuum_permittivity = 8.85418782e-12;  // [ F/m ]
const double Source::_relative_permittivity_silicon = 11.9;  // no unit.
void Source::save_Neff_dist(double zmin, double zmax)
{
    constexpr int nstep = 500;
    
    TFile *fout = new TFile("Neff_dist.root","RECREATE");
    TH1D hist = TH1D("neff_dist", "Effective Doping Profile", nstep, zmin, zmax);
    
    for( int i=0; i<nstep; i++)
    {            
        double tmp_z = zmin + (zmax-zmin)/static_cast<double>(nstep) * i;
        
        auto gauss_term1 = _peak_height[0] * std::exp(- std::pow((_peak_position[0] - tmp_z), 2.0)/(2*std::pow(_gauss_sigma[0], 2.0))) ;  // [ /cm^3 ]
        auto gauss_term2 = _peak_height[1] * std::exp(- std::pow((_peak_position[1] - tmp_z), 2.0)/(2*std::pow(_gauss_sigma[1], 2.0))) ;  // [ /cm^3 ]        

        auto permittivity = _vacuum_permittivity * _relative_permittivity_silicon * 1e-2;  // [ F/cm ]        
        auto base_term = ( _f_poisson * 1e4*1e4 ) * permittivity / _elementary_charge ;  // [ /cm^3 ]

        hist.SetBinContent(i+1, std::abs(base_term)+std::abs(gauss_term1)+std::abs(gauss_term2) );
    }

    hist.Write();    
    fout->Close();
}

/*
 * Method that calculates the weighting field inside the detector
 */
void SMSDetector::solve_w_f_grad()
{

	_L_g.u = _w_u;
	solve(_a_g == _L_g, _w_f_grad);
	// Change sign E = - grad(u)
	_w_f_grad = _w_f_grad * (-1.0);
}

/*
 * Method for calculating the Electric field
 */
void SMSDetector::solve_d_f_grad()
{
	_L_g.u = _d_u;
	solve(_a_g == _L_g, _d_f_grad);
	// Change sign E = - grad(u)
	_d_f_grad = _d_f_grad * (-1.0);

}

/*
 * Getter for the mesh of both fields
 */
RectangleMesh * SMSDetector::get_mesh()
{
	return &_mesh;
}

/*
 * From the carrier list taken from the file of carriers, this method runs through the vector of carriers and calls the simulate.drift for each one of them, depending if the carrier
 * is a hole or an electron. Trapping effects are included at the end directly in the valarry of currents.
 * Info related to diffusion is displayed using this method as well. Whenever a carrier crosses to depleted region, the acumulative variable totalCross grows.
 */
/**
 *
 * @param dt
 * @param max_time
 * @param shift_x
 * @param shift_y
 * @param curr_elec
 * @param curr_hole
 * @param totalCrosses
 *
 * Modified version !
 */
void CarrierCollection::simulate_drift( double dt, double max_time, double shift_x /*yPos*/, double shift_y /*zPos*/,
                                        std::valarray<double>&curr_elec, std::valarray<double> &curr_hole,
                                        std::valarray<double>&curr_gen_elec, std::valarray<double> &curr_gen_hole, double max_mul_factor,                               
                                        int &totalCrosses, const std::string &scanType, std::string skip_event_loop)
{
    double x_init, y_init;
    int totalCross = 0;
    bool control = true;

    //auto itr_queue = _queue_carrier_list.begin();
    //std::vector<std::unique_ptr<Carrier>> gen_carrier_list;    // container for carriers which will be generated by the impact ionization effect
    bool continue_loop = true;

    // Skip the following loop, if the corresponding flag is set
    if( skip_event_loop == "yes")
    {
        std::cout << std::endl;
        std::cout << "Event loop is going to be skipped ! " << std::endl;
        std::cout << "If you want to run normally, please change the setting (\"SkipEventLoop\") in config file " << std::endl;        
        std::cout << std::endl;
        
        continue_loop = false;
    }
    int loop = 0;
    
    //auto initial_n_carrier = (*itr_queue).size();
    auto generated_n_carrier = 0;
    
    while( continue_loop )
    {        
        std::vector<std::unique_ptr<Carrier>> gen_carrier_list;    // container for carriers which will be generated by the impact ionization effect

        // Declare the iterator and increment it here to avoid the invalidation of the iterator.
        // Probably, it is fine, as long as std::deque with "puch_back" insertion method, thus it is kinds of safety reason.
        auto itr_queue = _queue_carrier_list.begin();
        auto initial_n_carrier = (*itr_queue).size();
        for( int i=0; i<loop; i++)
        {
            ++itr_queue;
            //std::cout << "Loop=" << loop << " : Size of std::vector<std::unique_ptr<Carrier> =  " << (*itr_queue).size() << std::endl;            
            //std::cout << "Loop=" << loop << " : Size of std::deque<std::vector<std::unique_ptr<Carrier>>> =  " << _queue_carrier_list.size() << std::endl;
        }
        
        for (auto& carrier : *itr_queue )   // modified            
        {               
            char carrier_type = carrier->get_carrier_type();  
            // simulate drift and add to proper valarray
            if (carrier_type == 'e')
            {                
                // get and shift carrier position
                std::array<double,2> x = carrier->get_x();     
                x_init = x[0] + shift_x;
                y_init = x[1] + shift_y;

                //std::cout << "x[0]=" << x[0] << std::endl;
                
                //x[0] represents the X position read from the carriers file
                //shift_x represents the shift applied to X read from the steering file, namely, where the laser points to.
                //x[1] represents the Y position read from the carriers file. Y is seen as Z.
                //shift_y represents the Z (y) position read from the steering file. Since the program (usually) does edge-TCT, Z can be defined in one or more steps.                

                if( loop == 0 ) // Current made by initial carriers
                {
                    curr_elec += carrier->simulate_drift( dt , max_time, x_init, y_init, scanType,  gen_carrier_list); 
                }
                else                //  Current made by secondary carriers
                {
                    curr_gen_elec += carrier->simulate_drift( dt , max_time, x_init, y_init, scanType,  gen_carrier_list); 
                }
            }            
            else if (carrier_type =='h')
            {
                // get and shift carrier position
                std::array<double,2> x = carrier->get_x();      
                double x_init = x[0] + shift_x;
                double y_init = x[1] + shift_y ;               

                if( loop == 0 ) // Current made by initial carriers
                {
                    curr_hole += carrier->simulate_drift( dt , max_time, x_init, y_init, scanType,  gen_carrier_list);   
                }
                else
                {
                    curr_gen_hole += carrier->simulate_drift( dt , max_time, x_init, y_init, scanType,  gen_carrier_list);    
                }
            }
            //Let's see how many carriers from the carrier list for this step in Z have crossed to the depleted region.
            //A flag is switched to true on carrier.simulate_drift when a carrier filfull the requirements. See carrier.simulate_drift
            if (carrier->crossed()){    
                totalCross += 1;
            }
            
        } // end of for( auto carrier : *itr_queue )

        
        // push back the newly generated carrier list to the Queue
        if( !gen_carrier_list.empty() )
        {
            generated_n_carrier += gen_carrier_list.size();
            std::cout << "Initial number of carriers = " << initial_n_carrier << std::endl;
            std::cout << "Number of carriers generated in this loop = " << gen_carrier_list.size() << std::endl;            
            std::cout << "Total number of generated carriers = " << generated_n_carrier << std::endl;
            std::cout << "loop number = " << loop << std::endl;
            std::cout << std::endl;
            
            if( generated_n_carrier > initial_n_carrier * max_mul_factor )  // If number of carriers exceeds a thredhols, stop the loop. (to avoid endless loop)
            {
                std::cout << "Number of generated charges exceeds the maxmimum limit set in the configuration file !" << std::endl;
                std::cout << "Therefore, finishing the iteration of carrier drift loop. The results obtained by now is successfully recorded. " << std::endl;
                continue_loop = false;
            }
            else
            {
                continue_loop = true;
                _queue_carrier_list.push_back( std::move(gen_carrier_list) );                
                loop++;
            }
        }
        else
            continue_loop = false;
            
    }   // end of while loop 

    
	//fileDiffDrift.close();
	//std::cout << "Number of carriers crossed to DR in last Z step with Height " << shift_y << ": " << totalCross << std::endl;
	totalCrosses += totalCross;

	double trapping_time = _detector->get_trapping_time();

	for (double i = 0.; i < curr_hole.size(); i ++)
	{
		double elapsedT = i*dt;
		curr_elec[i] *= exp(-elapsedT/trapping_time);
		curr_hole[i] *= exp(-elapsedT/trapping_time);

        curr_gen_elec[i] *= exp(-elapsedT/trapping_time);
		curr_gen_hole[i] *= exp(-elapsedT/trapping_time);
	}

    // Record number of carriers vs gen_time
    record_carrier_gen_time( max_time, curr_hole.size() );
}

/*
 * Getter for the type of the CC (electro / hole)
 */
/**
 *
 * @return
 */
char Carrier::get_carrier_type()
{
	return _carrier_type; // electron or hole
}

/*
 * Getter for the position of the CC
 */

std::array< double,2> Carrier::get_x()
{
	return _x;
}

/*
 ******************** CARRIER DRIFT SIMULATION METHOD**************************
 *
 * Simulates how the CC drifts inside the detector in the 
 * desired number of steps
 *
 */
/**
 *
 * @param dt
 * @param max_time
 * @param x_init
 * @param y_init
 * @return
 */
std::valarray<double> Carrier::simulate_drift(double dt, double max_time, double x_init/*y_pos*/, double y_init /*z_pos*/, const std::string &scanType,
                                              std::vector<std::unique_ptr<Carrier>>& gen_carrier_list)    
//First Approach
{
    bool debug_flag = false;
    
    // For confirmation
    if(debug_flag)std::cout << "_gen_time = " << _gen_time << std::endl;
    
    // Newly Added for impact ionization effect
    std::vector<std::unique_ptr<Carrier>> _carrier_list_secondary;
    
	_x[0] = x_init;
	_x[1] = y_init;
	//std::ofstream fileDiffDrift;

	bool regularCarrier = true;
	double tDiff=0.;
	double tDep = 0;

	// get number of steps from time
	int max_steps = (int) std::floor(max_time / dt);
	std::valarray<double>  i_n(max_steps); // valarray to save intensity
	runge_kutta4<std::array< double,2>> stepper;
	// wrapper for the arrays using dolphin array class
	Array<double> wrap_x(2, _x.data());
	Array<double> wrap_e_field(2, _e_field.data());
	Array<double> wrap_w_field(2, _w_field.data());
        
	//Carrier is in NO depleted area
	if ( (_x[1] > _detector->get_depletionWidth()) && (_x[1] <= _detector->get_y_max()) && (_x[0] >= _detector->get_x_min()) && (_x[0] <= _detector->get_x_max())
			&& (_detector->diffusionON()) && (_x[1] >= _detector->get_y_min())){
        
		regularCarrier = false;
		//Four times the trapping time represent almost 100% of the signal.
		while( (tDiff < (4 * _trapping_time)) &&  (tDiff<(max_time))  ){
			_e_field_mod = 0;
			calculateDiffusionStep(dt); //Carrier movement due to diffusion
			tDiff += dt;
			if (_x[1] > _detector->get_y_max()) _x[1] = _detector->get_y_max();
			if (_x[1] < _detector->get_y_min()) _x[1] = _detector->get_y_min();
			if (_x[0] > _detector->get_x_max()) _x[0] = _detector->get_x_max();
			if (_x[0] < _detector->get_x_min()) _x[0] = _detector->get_x_min();

			if ((_x[1] < _detector->get_depletionWidth())){ //&& (_x[1] < _detector->get_y_max()) && (_x[0] > _detector->get_x_min()) && (_x[0] < _detector->get_x_max()) && (_x[1] > _detector->get_y_min())){
				regularCarrier = true;
				//To count carriers passing to the depleted region
				_crossed = true;

				break;
			}

		}

	}
	//End NO depleted area
	if ((regularCarrier)  && (_x[1] <= _detector->get_depletionWidth())){

		int it0 = ( _detector->diffusionON() ) ? TMath::Nint( (_gen_time + tDiff)/dt ) : TMath::Nint( _gen_time/dt ) ;

         // indicator of generating carriers.  i.e. ngen==2 shows that one e/h pair can be generated.    
        double ngen = 1.0; 
        
		for ( int i = it0 ; i < max_steps; i++)
		{
            
            if(debug_flag)std::cout << "_x[0] = " << _x[0] << " ,  _x[1] = " << _x[1] << ",  time = " << i*dt << std::endl;
            
			if (_x[1] <= _detector->get_depletionWidth()){

                // 1. Carrier position before drift
                auto pre_x0 = _x[0];
                auto pre_x1 = _x[1];
                
				//	safeRead.lock();
				_detector->get_d_f_grad()->eval(wrap_e_field, wrap_x);
				_detector->get_w_f_grad()->eval(wrap_w_field, wrap_x);
				//	safeRead.unlock();
				_e_field_mod = sqrt(_e_field[0]*_e_field[0] + _e_field[1]*_e_field[1]);
				stepper.do_step(_drift, _x, tDep, dt); //Carrier movement due to drift
				i_n[i] = _q *_sign* _mu.obtain_mobility(_e_field_mod) * (_e_field[0]*_w_field[0] + _e_field[1]*_w_field[1]);

                // 2. Carrier position diff between before and after its drift
                auto diff_x0 = _x[0] - pre_x0;
                auto diff_x1 = _x[1] - pre_x1;

                std::valarray<double> alpha_tmp = calculateAlpha( _e_field[1] * 1e+4 );    // efield , [V/um] --> [V/cm]
                double local_alpha;
                if(  _carrier_type == 'e' )
                    local_alpha = alpha_tmp[0];
                else
                    local_alpha = alpha_tmp[1];
                
                ngen *= std::exp( local_alpha * std::fabs(diff_x1) * 1e-4 );     // distance,  [um] --> [cm]

                if(debug_flag) std::cout << "time = " << i*dt << ",  ngen = " << ngen << "  _e_field[1] = " << _e_field[1] << ", local_alpha = " << local_alpha << std::endl;
                
                while( ngen > 2.0 )
                {
                    if(debug_flag)std::cout << "ngen = " << ngen << "_x[0] = " << _x[0] << " ,  _x[1] = " << _x[1]  << " :  time = " << i*dt << std::endl;

                    if( _carrier_type == 'e' )    // comment. Bellow,  the member variable "_sign" can be used to make the code simple. 
                    {
                        std::unique_ptr<Carrier> ptr_e( new Carrier('e', _q, _x[0], _x[1], _detector, i*dt) );  
                        gen_carrier_list.push_back( std::move(ptr_e) );

                        std::unique_ptr<Carrier> ptr_h( new Carrier('h', std::fabs(_q), _x[0], _x[1], _detector, i*dt) );  // Hole & electron pair should be created                    
                        gen_carrier_list.push_back( std::move(ptr_h) );
                    }
                    if( _carrier_type == 'h' )
                    {
                        std::unique_ptr<Carrier> ptr_e( new Carrier('e', -std::fabs(_q), _x[0], _x[1], _detector, i*dt) );  
                        gen_carrier_list.push_back( std::move(ptr_e) );

                        std::unique_ptr<Carrier> ptr_h( new Carrier('h', _q, _x[0], _x[1], _detector, i*dt) );  // Hole & electron pair should be created                    
                        gen_carrier_list.push_back( std::move(ptr_h) );
                    }                    

                    ngen -= 1.0;    // indicator decrement one. 
                }
			}
			if  (_detector->diffusionON()){
                
				calculateDiffusionStep(dt); //Carrier movement due to diffusion                

				if (_x[1] > _detector->get_depletionWidth() && !(_carrier_type == 'h' && scanType == "top")){
					_crossed = false;							//Holes are not followed anymore in top scan when reach no depleted area
					_e_field_mod = 0;
					//Come back to depleted area
					if (_x[1] > _detector->get_y_max()) _x[1] = _detector->get_y_max();
					if (_x[1] < _detector->get_y_min()) _x[1] = _detector->get_y_min();
					if (_x[0] > _detector->get_x_max()) _x[0] = _detector->get_x_max();
					if (_x[0] < _detector->get_x_min()) _x[0] = _detector->get_x_min();

				}
			}

			if (_detector->is_out_dep(_x)) // If CC outside detector
			{
                if(debug_flag)std::cout << "Carrier is out of detector !" << std::endl;
				break; // Finish (CC gone out)
			}

			tDep+=dt;//Could be inside first IF.TODO
		}
        
		return i_n;    
	}
	return i_n=0.;

}

double SMSDetector::get_depletionWidth()
{

	//_depletion_width = _depth * sqrt((_v_strips-_v_backplane)/_vdep);
	return _depletion_width;
}

/*
 * Getter for the maximum Y value.
 */
double SMSDetector::get_y_max()
{
	return _y_max;
}

/*
 * Getter for the minimum X value
 */
double  SMSDetector::get_x_min()
{
	return _x_min;
}

/*
 * Getter for the maximum X value
 */
double  SMSDetector::get_x_max()
{
	return _x_max;
}

/*
 * Getter for the minimum Y value.
 */
double  SMSDetector::get_y_min()
{
	return _y_min;
}


/*Diffusion method: When calling obtain_mobility method, the _mu object has been instanciated knowing whether it is a H or an E.
 * Thus, both obtain_mobility are the same method but with different variables.
 * The result of the diffusion is added to both coordenates Z(Y) and X modifying the position of the carrier.
 *
 */
/**
 *
 * @param dt
 */
void Carrier::calculateDiffusionStep(double dt){

	TRandom3 Rand(0);

	diffDistance = pow(2*_mu.obtain_mobility(_e_field_mod)*kB*_myTemp/(ECH)*dt,0.5);
	_dx = diffDistance * Rand.Gaus(0,1);
	_dy = diffDistance * Rand.Gaus(0,1);

	_x[0] += _dx;
	_x[1] += _dy;

	/*Becker Thesis, pag. 50. Lutz book, pag 18. Ejercicio de Lutz, pag. 36.
	diffDistance = pow(2*(1440*cm*cm/(V*s))*kB*_myTemp/(ECH)*dt,0.5);

	Other approach could be following MCTSi, Tim Janssen
	diff = sqrt(2*coefficientElectron*dt) * sin(2*pi*gRandom.Uniform()) * sqrt(-2*log(gRandom.Uniform()));*/


}


/*
 * Getter for the weighting field
 */
Function * SMSDetector::get_w_f_grad()
{
	return &_w_f_grad;
}

/*
 * OBTAIN MOBILITY - Getter method
 *
 * This method provides the value of the mobility using all the 
 * varibles initialized in the invokation.
 *
 */
/**
 *
 * @param e_field_mod
 * @return
 */
double JacoboniMobility::obtain_mobility(double e_field_mod)
{
  return _mu0/std::pow(1.0+std::pow(_mu0*e_field_mod/_vsat,_beta), 1.0/_beta); // mum**2/ Vs
}


std::valarray<double> Carrier::calculateAlpha(double efield)    
{
    std::valarray<double> a0 = { 4.3383           ,  2.376 };
    std::valarray<double> a1 = { -2.42e-12      ,  0.01033 };
    std::valarray<double> a2 = { 4.1233           ,  1.0 };

    std::valarray<double> b0 = { 0.235             ,  0.17714};
    std::valarray<double> b1 = { 0.0                 ,  -0.002178 };

    std::valarray<double> c0 = { 1.6831e+4     , 0.0 };
    std::valarray<double> c1 = { 4.3796            , 0.00947 };
    std::valarray<double> c2 = { 1.0                  , 2.4924 };
    std::valarray<double> c3 = { 0.13005          , 0.0 };

    std::valarray<double> d0 = { 1.233735e+6 ,  1.4043e+6 };
    std::valarray<double> d1 = { 1.2039e+3     ,  2.9744e+3 };
    std::valarray<double> d2 = { 0.56703          ,  1.4829 };
    
    //double T = 300.0; // 273K + 25K ~ 300K
    double T = _myTemp;
    
    std::valarray<double> aT = a0 + a1*std::pow(T, a2);
    std::valarray<double> bT = b0*std::exp(b1*T);
    std::valarray<double> cT = c0 + c1*std::pow(T, c2) + c3*T*T;
    std::valarray<double> dT = d0 + d1*T + d2*T*T;
    
    std::valarray<double> alpha = efield/( aT + bT*std::exp( dT/( efield + cT) ) );
 
    return alpha;
}

bool SMSDetector::is_out_dep(const std::array< double,2> &x)
{
	bool out = true;
	if ( (x[0] > _x_min) && (x[0] < _x_max) && (x[1] > _y_min) && (x[1] < _depth) )
	{
		out = false;
	}
	return out;
}

bool Carrier::crossed(){

	return _crossed;
}

/**
 * Record number of carriers
 **/
void CarrierCollection::record_carrier_gen_time(double max_time, int n_time_slice)
{
    TFile *f = new TFile("ncarrier.root", "RECREATE");

    TH1D *h_e_gen_time = new TH1D("e_gentime", "electron gen_time", n_time_slice, 0.0, max_time);
    TH1D *h_h_gen_time = new TH1D("h_gentime", "hole gen_time", n_time_slice, 0.0, max_time);    
   
    for( auto itr = _queue_carrier_list.begin();  itr != _queue_carrier_list.end(); itr++)
    {
        for (auto& carrier : *itr ) 
        {
            char carrier_type = carrier->get_carrier_type();
            double gen_time = carrier->get_gen_time();

            if( carrier_type == 'e') 
                h_e_gen_time->Fill( gen_time );
            else
                h_h_gen_time->Fill( gen_time );
        }
    }

    h_e_gen_time->Write();
    h_h_gen_time->Write();

    f->Close();  
}


/*
 * Convert i_total to TH1D after simulating simple RC circuit. ROOT based method.
 */
void TRACSInterface::GetItRc()
//TH1D * TRACSInterface::GetItRc()
{

	double RC = 50.*C; // Ohms*Farad
	double alfa = dt/(RC+dt);


	for (int j = 1; j <n_tSteps; j++)
	{
		i_shaped[j]=i_shaped[j-1]+alfa*(i_total[j]-i_shaped[j-1]);

	}
	count2++;


}

void TRACSInterface::currents_hist_to_file(int tid, int vPos, int nscan)
{
    TString file_name;
    if( ! _after_fitting_process )
    {
        file_name.Form("current%dV_scan%d.root", (int)( voltVector[vPos] ) , nscan);
    }
    else
    {
        file_name.Form("current%dV_scan%d_wFit.root", (int)( voltVector[vPos] ) , nscan);
    }
    TFile *fout = new TFile(file_name, "RECREATE");
    
    TH1D *h_i_total = new TH1D("i_total", "total current", n_tSteps, 0.0, max_time);
    TH1D *h_i_init_elec = new TH1D("i_init_elec", "current induced by initial electron", n_tSteps, 0.0, max_time);
    TH1D *h_i_init_hole = new TH1D("i_init_hole", "current induced by initial hole", n_tSteps, 0.0, max_time);
    TH1D *h_i_gen_elec = new TH1D("i_gen_elec", "current induced by secondary electron", n_tSteps, 0.0, max_time);
    TH1D *h_i_gen_hole = new TH1D("i_gen_hole", "current induced by secondary hole", n_tSteps, 0.0, max_time);     
    
    std::valarray<double> i_total_shaped;
    std::valarray<double> i_init_elec_shaped;
    std::valarray<double> i_init_hole_shaped;
    std::valarray<double> i_gen_elec_shaped;
    std::valarray<double> i_gen_hole_shaped;                    
    
    i_total_shaped.resize( static_cast<size_t>(n_tSteps) );
    i_init_elec_shaped.resize( static_cast<size_t>(n_tSteps) );
    i_init_hole_shaped.resize( static_cast<size_t>(n_tSteps) );
    i_gen_elec_shaped.resize( static_cast<size_t>(n_tSteps) );
    i_gen_hole_shaped.resize( static_cast<size_t>(n_tSteps) );                    

    // Using defalt RC shaping. May need to consider in future, which signal transformation is suitable.
    GetItRc( i_total_shaped, i_total );
    GetItRc( i_init_elec_shaped, i_elec );
    GetItRc( i_init_hole_shaped, i_hole );
    GetItRc( i_gen_elec_shaped, i_gen_elec );
    GetItRc( i_gen_hole_shaped, i_gen_hole );                                                            
    
    for( int bin_index=0; bin_index < n_tSteps ;  bin_index++)
    {
        h_i_total->SetBinContent( bin_index+1, i_total_shaped[bin_index] );
        h_i_init_elec->SetBinContent( bin_index+1, i_init_elec_shaped[bin_index] );
        h_i_init_hole->SetBinContent( bin_index+1, i_init_hole_shaped[bin_index] );
        h_i_gen_elec->SetBinContent( bin_index+1, i_gen_elec_shaped[bin_index] );
        h_i_gen_hole->SetBinContent( bin_index+1, i_gen_hole_shaped[bin_index] );                                                 
    }
    
    h_i_total->Write();
    h_i_init_elec->Write();
    h_i_init_hole->Write();
    h_i_gen_elec->Write();
    h_i_gen_hole->Write();
    
    fout->Close();                     
}

/*
 * Convert i_total to TH1D after simulating simple RC circuit. ROOT based method.
 * Not only for specific current, i_total, but for general usage.
 */
void TRACSInterface::GetItRc(std::valarray<double>& i_out,  const std::valarray<double>& i_in)
{

	double RC = 50.*C; // Ohms*Farad
	double alfa = dt/(RC+dt);

    if( i_out.size() != n_tSteps || i_in.size() != n_tSteps )
    {
        for(int i=0; i<50; i++){ std::cout << "" << std::endl; }
    }
    
	for (int j = 1; j <n_tSteps; j++)
	{
		i_out[j]=i_out[j-1]+alfa*(i_in[j]-i_out[j-1]);
	}
	//count2++;
}

void TRACSInterface::fields_hist_to_file(int tid, int vPos)
{

	/*Exporting 2D Histograms to file*/
	// get plot and set new data
	TString file_name;
	file_name.Form("wf%dV", TMath::Nint(voltVector[vPos]));//"2Dhistos"+voltages[vPos]+"V"+".root";
	TFile *fout = new TFile(file_name,"RECREATE");

	TString wpm;
	wpm.Form("WP_%d_V", TMath::Nint(voltVector[vPos]) ) ;
	TString wfm;
	wfm.Form("WF_%d_V", TMath::Nint(voltVector[vPos]) ) ;
	TString efm;
	efm.Form("EF_%d_V", TMath::Nint(voltVector[vPos]) ) ;

	TH2D h_w_u      = utilities::export_to_histogram    ( *detector->get_w_u(), "h_w_u", wpm /*"Weighting potential"*/, detector->get_n_cells_x(), detector->get_x_min(), detector->get_x_max(), detector->get_n_cells_y(), detector->get_y_min(), detector->get_y_max());
	TH2D h_w_f_grad = utilities::export_mod_to_histogram( *detector->get_w_f_grad(), "h_w_f_grad", wfm /*"Weighting field"*/, detector->get_n_cells_x(), detector->get_x_min(), detector->get_x_max(), detector->get_n_cells_y(), detector->get_y_min(), detector->get_y_max());
	TH2D h_d_f_grad = utilities::export_mod_to_histogram( *detector->get_d_f_grad(), "h_d_f_grad", efm /*"Electric field"*/, detector->get_n_cells_x(), detector->get_x_min(), detector->get_x_max(), detector->get_n_cells_y(), detector->get_y_min(), detector->get_y_max());

    TH1D h_d_f_grad_y = utilities::export_1D_histogram( *detector->get_d_f_grad(), "h_d_f_grad_Y", efm /*"Electric field"*/, detector->get_n_cells_x(), detector->get_x_min(), detector->get_x_max(), detector->get_n_cells_y(), detector->get_y_min(), detector->get_y_max() );
    
	h_w_u.Write();
	h_w_f_grad.Write();
	h_d_f_grad.Write();

    h_d_f_grad_y.Write();
    
	fout->Close();
}

TH2D utilities::export_to_histogram(Function &func, TString hist_name, TString hist_title, int n_bins_x , double x_min, double x_max, int n_bins_y, double y_min, double y_max)
{
	TH2D hist = TH2D(hist_name,hist_title, n_bins_x , x_min, x_max, n_bins_y, y_min, y_max);
	double step_x = (x_max -x_min)/n_bins_x;
	double step_y = (y_max -y_min)/n_bins_y;
	for (int i=1; i<=n_bins_x; i++) {
		for (int j=1; j<=n_bins_y; j++) {
			double x_value = (i - 0.5)*step_x;
			double y_value = (j - 0.5)*step_y;
			hist.Fill(x_max-x_value,y_max-y_value, func(x_value, y_value));
		}
	}
	return hist;
}

TH2D utilities::export_mod_to_histogram(Function &func, TString hist_name, TString hist_title, int n_bins_x , double x_min, double x_max, int n_bins_y, double y_min, double y_max)
{

	double x_value;
	double y_value;
	Function * e_f_grad = &func;
	double e_field_mod;
	double e_f_x;
	double e_f_y;


	double step_x = (x_max -x_min)/n_bins_x;
	double step_y = (y_max -y_min)/n_bins_y;
	TH2D hist = TH2D(hist_name,hist_title, n_bins_x , x_min-0.5*step_x, x_max+0.5*step_x, n_bins_y, y_min-0.5*step_y, y_max+0.5*step_y);

	double xval , yval ;
	for (int x=1; x<=n_bins_x; x++) {
		xval=x_min+(x-1)*step_x;
		for (int y=1; y<=n_bins_y; y++) {
			yval=y_min+(y-1)*step_y;
			e_f_x = ((*e_f_grad)[0])(xval , yval);
			e_f_y = ((*e_f_grad)[1])(xval , yval);
			e_field_mod = sqrt(e_f_x*e_f_x +  e_f_y *e_f_y);

			hist.SetBinContent( n_bins_x-x+1,n_bins_y-y+1, e_field_mod);
		}
	}
	return hist;
}

TH1D utilities::export_1D_histogram(Function &func, TString hist_name, TString hist_title, int n_bins_x , double x_min, double x_max, int n_bins_y , double y_min, double y_max)
{

	double x_value;
	double y_value;
	Function * e_f_grad = &func;
	double e_field_mod;
	double e_f_x;
	double e_f_y;


	double step_x = (x_max -x_min)/n_bins_x;
	double step_y = (y_max -y_min)/n_bins_y;
	TH1D hist = TH1D(hist_name,hist_title, n_bins_y , y_min-0.5*step_y, y_max+0.5*step_y);

	double xval , yval ;
	for (int x=1; x<=n_bins_x; x++) {
		xval=x_min+(x-1)*step_x;
		for (int y=1; y<=n_bins_y; y++) {
			yval=y_min+(y-1)*step_y;
			e_f_x = ((*e_f_grad)[0])(xval , yval);
			e_f_y = ((*e_f_grad)[1])(xval , yval);
			e_field_mod = sqrt(e_f_x*e_f_x +  e_f_y *e_f_y);

            // if it is center of X-axis
            if( x == static_cast<int>( n_bins_x/2.0 ) )
            {
                hist.SetBinContent( n_bins_y-y+1, e_f_y);
            }
		}
	}
	return hist;
}


#define EXE 1      //1 for shared lib inside root. Note that still a ".o" can be produced in this mode
                   //by compiling the source, as needed for plotvd.C              

		   //2 for executable outside root
#define CPF 2.0		   

#define CTRLPLOT 0      //1 if intermediate control plots wanted


//std::mutex mtx_conv;           // mutex for critical section


TH1D *H1DConvolution( TH1D *htf , TH1D *htct , Double_t Cend , int tid) { 
      
   //------------>Both input histograms should have the same bin width<-----------------
   //Here you can apply an extra LPFiltering
   //if (Cend!=0) htct = LPFilter( htct , Cend );
   int count=0;
   //Convolute (commutative)
   //C(t) = Int[ tct(x) transferfunction(t-x) dx ]
   Double_t bw = htct->GetBinCenter(2) - htct->GetBinCenter(1);
   //Double_t bw = htct->GetXaxis()->GetBinCenter(2) - htct->GetXaxis()->GetBinCenter(1);
   //TF1 *f1 = new TF1("f1","abs(sin(x)/x)*sqrt(x)",0,10);
   //   float r = f1->GetRandom();
   TString tftit, tfname;
   tftit.Form("hConv_%d_%d", tid, count);
   tfname.Form("conv_%d_%d", tid, count);
   TH1D *hConv = new TH1D(tftit,tfname,2*htct->GetNbinsX(),-htct->GetNbinsX()*bw,htct->GetNbinsX()*bw);
   
   std::cout << "test1" << std::endl;
   //The convoluted response to the TCT signal is going to be another histogram sized similar to htct
   //Int_t Ntf = htf->GetNbinsX() , Ntct = htct->GetNbinsX();
   Int_t Ntf = htf->GetNbinsX();//20210328 first bug fix: uncomment it
   TH1D *nhtf = (TH1D *) htf->Clone();
   //int Ntf;//20210328 fisrt bug fix: comment it
   std::cout << "test2" << std::endl;
   Int_t Ntct = htct->GetNbinsX();
   
   //Create the reverse histogram of the transfer function
   TH1D *hinv = (TH1D *) htf->Clone(); hinv->Reset();
   for (Int_t j=1; j<= Ntf ; j++) hinv->SetBinContent( j , htf->GetBinContent(Ntf-j+1) );  
  // std::basic_string<char> t= std::to_string(tid);
   TH1D *hg = (TH1D *) htct->Clone(); 
   
   hinv->Draw();
   #if CTRLPLOT == 1
   gPad->Print( "hgcontrol.pdf[" ) ; gPad->Print( "hgcontrol.pdf" );
   #endif
   
    for ( Int_t i=1 ; i<=2*Ntct ; i++ )
    { 
          
     //Create a shifted histogram version of the inverse
     hg->Reset();
     //for (Int_t j=TMath::Nint(-0.5*Ntf); j<= TMath::Nint(0.5*Ntf) ; j++) hg->SetBinContent( i-j , hinv->GetBinContent( j+TMath::Nint(0.5*Ntf) ) );  
     if ( i<=Ntf ) 
    {
       //Histogram is shifting in from the left
       for (Int_t j=1; j<=i ; j++) hg->SetBinContent( j , hinv->GetBinContent( Ntf-i+j ) );
    } else
   {
       //Histogram is shifting out. Leaving from the right
       Int_t cont=1 ;
       for (Int_t j=i-Ntf+1; j<=2*Ntf ; j++) 
      {
         hg->SetBinContent( j , hinv->GetBinContent( cont ) );
	 cont++;
       }
  
    }
     
     
     //Multiply f(tau)*g(t-tau)
     hg->Multiply( htct );
     
     //Double_t fxg = hg->Integral("width");
     Double_t fxg = hg->Integral();
     
     hConv->SetBinContent(i,fxg);
     
     #if CTRLPLOT==1
       THStack *hst=new THStack("hst","conv");
       hConv->SetLineColor(2);hConv->SetLineWidth(2);
       hst->Add(hg) ; hst->Add(hConv);
       hst->Draw("nostack");
       if   (i==2*Ntct) { 
         gPad->Print( "hgcontrol.pdf" )  ; gPad->Print( "hgcontrol.pdf]" ) ; 
       } else if (i%10==0)  gPad->Print( "hgcontrol.pdf" )  ;
     #endif

   }

   return hConv;  

}


TH1D *H1DConvolution( TH1D *htct , Double_t Cend, int tid, TString transferFun) 
{
   
   //mtx_conv.lock();

   TFile  *ftf = new TFile(transferFun);
   TH1D   *htf = (TH1D *) ftf->Get("htf");
   TH1D *hConv = H1DConvolution(htf,htct,Cend,tid);
   hConv->SetDirectory(0);
   ftf->Close();
   //mtx_conv.unlock();
   return hConv;
   
}

/*
 * Writing to a single file
 *
 *
 */
/**
 *
 * @param tid
 */
void TRACSInterface::write_to_file(int tid)
{

	int index_conv = 0;
	if (scanType == "edge"){

		for (int i = 0 ; i < voltVector.size() ;  i++){
			for (int j = 0 ; j < zVector.size() ; j++){
				if (global_TF != "NO_TF")
					utilities::write_to_file_row(hetct_rc_filename, i_conv_vector[index_conv], detector->get_temperature(), yInit, zVector[j], voltVector[i]);
				else utilities::write_to_file_row(hetct_rc_filename, i_rc_array[index_conv], detector->get_temperature(), yInit, zVector[j], voltVector[i]);
				index_conv++;
			}

		}

	}



	if (scanType == "top" || scanType == "bottom"){

		for (int i = 0 ; i < voltVector.size() ;  i++){
			for (int j = 0 ; j < yVector.size() ; j++){
				if (global_TF != "NO_TF")
					utilities::write_to_file_row(hetct_rc_filename, i_conv_vector[index_conv], detector->get_temperature(), yVector[j], zInit, voltVector[i]);
				else utilities::write_to_file_row(hetct_rc_filename, i_rc_array[index_conv], detector->get_temperature(), yVector[j], zInit, voltVector[i]);
				index_conv++;
				//std::cout << i_rc_array[i][j] << std::endl;
			}

		}

	}




}



// function to write results to file (in rows)
// overloaded (now from TH1D)
void utilities::write_to_file_row(std::string filename, TH1D *hconv, double temp, double yShift, double height, double voltage)
{
	unsigned long int steps = hconv->GetNbinsX();
	height = height/1000.;
	yShift = yShift/1000.;

	std::ofstream out; // open file
	out.open(filename, std::ios_base::app);
	if (out.is_open())
	{
		out << steps << " ";
		out << temp-273. << " ";
		out << voltage << " ";
		out << "0 " << " " << yShift << " " << height << " ";


		//Shift simulated data to gen_time by padding 0's before gen_time
		//for (unsigned int i = 0; i < TMath::Nint(TRACSsim[0]->get_genTime()/TRACSsim[0]->get_dt()); i++ ) out << std::fixed << std::setprecision(9) << 0. << " ";

		// Scan on times
		for (unsigned int i = 1; i <= steps; i++ )
		{
			out << std::fixed << std::setprecision(9) << hconv->GetBinContent(i) << " ";
		}

		out << std::endl;
		out.close();
	}
	else // Error output
	{
		std::cout << "File could not be read/created"<<std::endl;
	}
}

/**
 *
 * @param x_min
 * @param x_max
 * @param depth
 */
PeriodicLateralBoundary::PeriodicLateralBoundary(double x_min, double x_max, double depth)
{
  _x_min = x_min;
  _x_max = x_max;
  _depth = depth;
}

/**
 *
 * @param pitch
 * @param width
 * @param nns
 */
CentralStripBoundary::CentralStripBoundary(double pitch, double width, int nns )
{
  _pitch = pitch;
  _width = width;
  _nns = nns;
}

/**
 *
 * @param pitch
 * @param width
 * @param nns
 */
NeighbourStripBoundary::NeighbourStripBoundary(double pitch, double width, int nns )
{
  _pitch = pitch;
  _width = width;
  _nns = nns;
}


/**
 *
 * @param x_min
 * @param x_max
 * @param depth
 */
BackPlaneBoundary::BackPlaneBoundary(double x_min, double x_max, double depth)
{
  _x_min = x_min;
  _x_max = x_max;
  _depth = depth;
}


/**
 *
 * @param carrier_type
 * @param d_f_grad
 * @param givenT
 * @param diffusion
 * @param dt
 */
DriftTransport::DriftTransport(char carrier_type, Function * d_f_grad, double givenT, int diffusion, double dt) :
_mu(carrier_type, givenT),
_diffusion(diffusion),
_dt(dt),
_temp(givenT)
{
	_d_f_grad = d_f_grad;
	if (carrier_type == 'e') {
		_sign = -1;
	}
	else {
		_sign = 1;
	}
}

JacoboniMobility::JacoboniMobility( char carrier_type, double T)
{
  _T = T;
  if (carrier_type == 'e') // Electrons
  {
    _mu0 = 1440.e8 * std::pow(_T/300., -2.260); //um**2/ Vs
    _vsat = 1.054e11  * std::pow(_T/300., -0.602); //um**2/ Vs
    _beta = 0.992 * std::pow(_T/300., 0.572); // <100> orientation

	//Extended Canali model. Following Sentaurus manual. Pag 217-224.
    //_mu0  = 1417.e8 * std::pow(_T/300., -2.5); //um**2/ Vs
    //_vsat = 1.07e11  * std::pow(_T/300., -0.87); //um**2/ Vs
    //_beta = 1.109 * std::pow(_T/300., 0.66); // <100> orientation
  }

  //Extended Canali model. Following Sentaurus manual. Pag 217-224.
  else if (carrier_type == 'h') // Holes
  {
    _mu0 = 474.e8 * std::pow(_T/300., -2.619); //um**2/ Vs
    _vsat = 0.940e11  * std::pow(_T/300., -0.226); //um**2/ Vs
    _beta = 1.181 * std::pow(_T/300., 0.633 ); // <100> orientation

	  //Extended Canali model. Following Sentaurus manual. Pag 217-224.
    //_mu0  = 470.e8 * std::pow(_T/300., -2.2); //um**2/ Vs
    //_vsat = 0.837e11  * std::pow(_T/300., -0.52); //um**2/ Vs
    //_beta = 1.213 * std::pow(_T/300., 0.17 ); // <100> orientation
  }
}

JacoboniMobility::~JacoboniMobility()
{

}

DriftTransport::~DriftTransport()
{

}

BackPlaneBoundaryWP::BackPlaneBoundaryWP(double x_min, double x_max, double depletion)
{
  _x_min = x_min;
  _x_max = x_max;
  _depletion_width = depletion;
}

int  SMSDetector::get_n_cells_x()
{
	return _n_cells_x;
}
int  SMSDetector::get_n_cells_y()
{
	return _n_cells_y;
}

Function * SMSDetector::get_w_u()
{
	return &_w_u;
}

// Destructor
TRACSInterface::~TRACSInterface()
{
	delete i_ramo;
	delete i_rc;
	delete i_conv;
	delete carrierCollection;
	delete detector;
}

/*
 ********************** DESTRUCTOR OF THE CLASS CARRIER	**************************
 */
Carrier::~Carrier()
{

}


/**
 *
 * @param x
 * @param dxdt
 * @param
 */
void DriftTransport::operator() ( const std::array<double,2>  &x , std::array<double,2>  &dxdt , const double /* t */ )
{
	Array<double> e_field((std::size_t) 2); // temp wrap for e. field
	double e_field_mod;
	//TRandom3 Rand(0);
	Point eval_point(x[0],x[1],0.0);
	(*_d_f_grad)(e_field, eval_point);
	e_field_mod = sqrt(e_field[0]*e_field[0] + e_field[1]*e_field[1]);
	dxdt[0] = _sign*_mu.obtain_mobility(e_field_mod) * e_field[0];
	dxdt[1] = _sign*_mu.obtain_mobility(e_field_mod) * e_field[1];
}

/*
 ********************** DESTRUCTOR OF THE CLASS CARRIER	COLLECTION **************************
 */
CarrierCollection::~CarrierCollection()
{

}


SMSDetector::~SMSDetector()
{

}

/**
 *
 * @param x
 * @param on_boundary
 * @return
 */
// Left boundary is "target domain"
bool PeriodicLateralBoundary::inside(const Array<double>& x, bool on_boundary) const
{
  return (std::abs(x[0]) < DOLFIN_EPS); // left lateral domain
}

// Map right boundary to left boundary
/**
 *
 * @param x
 * @param y
 */
void PeriodicLateralBoundary::map(const Array<double>& x, Array<double>& y) const
{
  y[0] = x[0] - _x_max; //translate x coordinate
  y[1] = x[1];  // leave y equal
}

/**
 *
 * @param x
 * @param on_boundary
 * @return
 */
bool NeighbourStripBoundary::inside(const Array<double>& x, bool on_boundary) const
{
  bool is_inside = false;
  if ((x[1] < DOLFIN_EPS ) && on_boundary) // y = 0 condition
  {
    int t_nns = 1 + 2*_nns; // total number of strips
    for ( int count = 0; count < t_nns; count++) // for loop for each strip
    {
      if (count != _nns) // not central strip
      {
        double x_translated = x[0] - _pitch*count; // reference system change
        double l_lim = (_pitch - _width) / 2.0;
        double r_lim = l_lim + _width;
        if ((x_translated > l_lim*(1-DOLFIN_EPS)) && (x_translated < r_lim*(1+DOLFIN_EPS))) // check if strip
        {
          is_inside = true;
        }
      }
    }
  }
  return is_inside;
}

bool BackPlaneBoundary::inside(const Array<double>& x, bool on_boundary) const
{
  bool is_inside = false;
  if ((x[1] > (_depth - DOLFIN_EPS*_depth) ) && on_boundary) // y = depth condition
  {
    if ((x[0] > _x_min - DOLFIN_EPS) && (x[0] < _x_max + DOLFIN_EPS )) // within boundaries
    {
      is_inside = true;
    }
  }
  return is_inside;
}


bool BackPlaneBoundaryWP::inside(const Array<double>& x, bool on_boundary) const
{
  bool is_inside = false;
  if (x[1] > ((_depletion_width) - DOLFIN_EPS* (_depletion_width)) )
  {
      is_inside = true;
 }
 return is_inside;
}


void Source::eval(Array<double>& values, const Array<double>& x) const
{
    if (NeffApproach == "Triconstant") 
    {
        /*
         * 3 ZONE constant space distribution
         *
         * We define here a Neff distribution consisting in 3 different zones
         * each zone is defined as a constant within the given region.
         * It uses all but the last parameter (y3 = Neff(z3)). It takes zX 
         * values as boundaries of the zones and the three first yX as the 
         * value of Neff in each region
         *
         * Even though a function like this is generally not continuous, we 
         * add the hyperbolic tangent bridges to ensure not only continuity 
         * but also derivability.
         *
         */
        double neff_1 = y0;
        double neff_2 = y1;
        double neff_3 = y2;
        
        // For continuity and smoothness purposes
        double bridge_1 = tanh(1000*(x[1]-z0)) - tanh(1000*(x[1]-z1));
        double bridge_2 = tanh(1000*(x[1]-z1)) - tanh(1000*(x[1]-z2));
        double bridge_3 = tanh(1000*(x[1]-z2)) - tanh(1000*(x[1]-z3));
        
        double neff = 0.5*((neff_1*bridge_1)+(neff_2*bridge_2)+(neff_3*bridge_3));
        values[0] = neff*0.00152132;
    }
    else if (NeffApproach == "Linear") 
    {
        /*
         * 1 ZONE approximatin
         *
         * First aproximation to the after-irradiation space charge distribution
         * Consists on a simple straight line defined by the points (z0, y0) and 
         * (z3, y3) and neglects the rest of the values.
         *
         */
        
        double neff = ((y0-y3)/(z0-z3))*(x[1]-z0) + y0;
        values[0] = neff*0.00152132;
    }
    else if ( NeffApproach == "AvalancheMode" )
    {
        /*
         * Gaussian approximation. 
         * Pedestal part is "_f_poisson"  parameter
         */
        
        auto permittivity = _vacuum_permittivity * _relative_permittivity_silicon;  // [ F/m ]
        
        auto poisson_term1 =  ((std::signbit(_f_poisson)== false) ? +1.0 : -1.0) * ( _elementary_charge * _peak_height[0] * 1e6 / permittivity );  // [ V/m/m ]
        auto poisson_term_unit_in_microm1 = poisson_term1 * 1e-12;  // [ V/um/um ]

        auto poisson_term2 =  ((std::signbit(_f_poisson)== false) ? +1.0 : -1.0) * ( _elementary_charge * _peak_height[1] * 1e6 / permittivity );  // [ V/m/m ]
        auto poisson_term_unit_in_microm2 = poisson_term2 * 1e-12;  // [ V/um/um ]
                
        auto gauss_term1 = poisson_term_unit_in_microm1 * std::exp(- std::pow((_peak_position[0] - x[1]), 2.0)/(2*std::pow(_gauss_sigma[0], 2.0))) ;
        auto gauss_term2 = poisson_term_unit_in_microm2 * std::exp(- std::pow((_peak_position[1] - x[1]), 2.0)/(2*std::pow(_gauss_sigma[1], 2.0))) ;
        
        
        auto base_term = _f_poisson ;  // [ V/um/um ]
        
        values[0] = base_term + gauss_term1 + gauss_term2;
    }
    else 
    {
        /*
         * 3 ZONE space distribution
         *
         * It consists in 3 different straight lines corresponding to 3 different
         * charge distributions. It uses all 8 parameters to compute the Neff.
         *
         * Continuity is assumed as straight lines have common points, continuity 
         * is ensured by the hyperbolic tangent bridges
         */
        double neff_1 = ((y0-y1)/(z0-z1))*(x[1]-z0) + y0;
        double neff_2 = ((y1-y2)/(z1-z2))*(x[1]-z1) + y1;
        double neff_3 = ((y2-y3)/(z2-z3))*(x[1]-z2) + y2;
        
        // For continuity and smoothness purposes
        double bridge_1 = tanh(1000*(x[1]-z0)) - tanh(1000*(x[1]-z1));
        double bridge_2 = tanh(1000*(x[1]-z1)) - tanh(1000*(x[1]-z2));
        double bridge_3 = tanh(1000*(x[1]-z2)) - tanh(1000*(x[1]-z3));
        
        double neff = 0.5*((neff_1*bridge_1)+(neff_2*bridge_2)+(neff_3*bridge_3));
        values[0] = neff*0.00152132;
        
    }
    // Fix units from the PdC version
    
    
}

/**
 *
 * @param x
 * @param on_boundary
 * @return
 */
bool CentralStripBoundary::inside(const Array<double>& x, bool on_boundary) const
{
  bool is_inside = false;
  if ((x[1] < DOLFIN_EPS ) && on_boundary) // y = 0 condition
  {
    double x_translated = x[0] - _pitch*_nns; // reference system change
    double l_lim = (_pitch - _width) / 2.0;
    double r_lim = l_lim + _width;
    if ((x_translated > l_lim*(1-DOLFIN_EPS)) && (x_translated < r_lim*(1+DOLFIN_EPS))) // check if strip
    {
      is_inside = true;
    }
  }
  return is_inside;
}





























void spread_into_threads()
{

	double vInit, deltaV, vMax, v_depletion, zInit, zMax, deltaZ, yInit, yMax, deltaY;
	int zSteps, ySteps, vSteps, seeker;
	bool lastElement = false;


	uint index = 1;

	//Reading file to get values that will determine scan vectors
	utilities::parse_config_file(fnm, scanType, vInit, deltaV, vMax, v_depletion, zInit, zMax, deltaZ, yInit,
			yMax, deltaY, dTime, max_time, capacitance, global_TF);

	//Reserve to 300 since is a common measurement for a detector Z dimension.
	//Reserve does not implies anything, better for performance.
	//For sure there is at least one value of Z coordinates in an edge-TCT scan
	vector_zValues.reserve(300);
	vector_zValues.push_back(zInit);
	seeker = zInit + deltaZ;

	//Filling Z scan vector
	while (seeker < zMax)
	{
		vector_zValues.push_back(seeker);
		++index;
		seeker = seeker + deltaZ;
		if (seeker > zMax)
		{
			vector_zValues.push_back(zMax);
			lastElement = true;
		}

	}
	if (seeker >= zMax && !lastElement && zMax !=zInit) vector_zValues.push_back(zMax);
	lastElement = false;

	//Reserve to 100 for differente voltages.
	//There will be at least one bias voltage
	vector_voltValues.reserve(300);
	vector_voltValues.push_back(vInit);
	seeker = vInit + deltaV;

	//Filling voltage scan vector
	while (seeker < vMax)
	{
		vector_voltValues.push_back(seeker);
		++index;
		seeker = seeker + deltaV;
		if (seeker > vMax)
		{
			vector_voltValues.push_back(vMax);
			lastElement = true;
		}
	}
	if (seeker >= vMax && !lastElement && vMax != vInit) vector_voltValues.push_back(vMax);
	lastElement = false;

	//Reserve to 400 since is a common measurement for a detector 4 dimension.
	//For sure there is at least one value of Z coordinates in an edge-TCT scan
	vector_yValues.reserve(400);
	vector_yValues.push_back(yInit);
	seeker = yInit + deltaY;

	//Filling Y scan vector
	while (seeker < yMax)
	{
		vector_yValues.push_back(seeker);
		++index;
		seeker = seeker + deltaY;
		if (seeker > yMax)
		{
			vector_yValues.push_back(yMax);
			lastElement = true;
		}
	}
	if (seeker >= yMax && !lastElement && yMax != yInit) vector_yValues.push_back(yMax);
	lastElement = false;

	//Checking for a right user configuration
	if ((vector_yValues.size() > 1) && (vector_zValues.size() > 1))
	{
		std::cout << "TRACS does not allow this configuration, please, choose zScan OR yScan" << std::endl;
		std::quick_exit(1);
	}

}





