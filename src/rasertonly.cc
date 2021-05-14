






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
  return (std::abs(x[0]) < DOLFIN_EPS ); // left lateral domain
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
  if ((x[1] > (_depth - DOLFIN_EPS *_depth) ) && on_boundary) // y = depth condition
  {
    if ((x[0] > _x_min - DOLFIN_EPS ) && (x[0] < _x_max + DOLFIN_EPS )) // within boundaries
    {
      is_inside = true;
    }
  }
  return is_inside;
}


bool BackPlaneBoundaryWP::inside(const Array<double>& x, bool on_boundary) const
{
  bool is_inside = false;
  if (x[1] > ((_depletion_width) - DOLFIN_EPS * (_depletion_width)) )
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
























int main( int argc, char *argv[]) 
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
		t[i] = std::thread(call_from_thread, i, std::ref(carrierThread_fileNames[i]), vector_zValues, vector_yValues, vector_voltValues);///////////////////////////////////////<<<<<<<<<<<working on1
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






