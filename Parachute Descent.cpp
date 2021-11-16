/******************************************************************************

Welcome to GDB Online.
GDB online is an online compiler and debugger tool for C, C++, Python, Java, PHP, Ruby, Perl,
C#, VB, Swift, Pascal, Fortran, Haskell, Objective-C, Assembly, HTML, CSS, JS, SQLite, Prolog.
Code, Compile, Run and Debug online from anywhere in world.

*******************************************************************************/
#include <iostream>
#include <cmath>
#include <string>

using namespace std;

// { Input data
double const pi = 3.1415926535;  // Pi
double xmu = 28.966 / 8315;      // mu/R  (  P = (rho/mu)RT  )
double Cx;              // drag cefficient
double rho;             // air dencity  (kg/m^3)
double Sp;              // area of parachute opening (m^2)
double Pg;              // overload  (N)
double Pg_;             // overload if parchute were completely inflated
double mg;              // overload  (kg)
double mg_;             // overload if parchute were completely inflated (kg)
double acc;             // acceleration
double P;               // air pressure at the current altitude     (Pa)
double T;               // air temperature at the current altitude  (K)
double g = 9.80665;     // free fall acceleration
double acc0 = g;        // initial/previous acceleration
double N_gore = 8;      // the number of parachute gores
double m_box =  1.65;   //1.108;  //0.710;    // 0.5; 1.3;     // payload  (kg)
double Cx_box = 1.4;    // Cx of payload box
double a_box = 0.22;    // payload box width
double b_box = 0.22;    // payload box length
int    k_box = 1;       // 1 - account for payload box drag, 0 - do not account for payload box drag
double m_balloon = 0.0; // 0.05; 1.2; // balloon mass  (kg) 
double Cx_all = 1.1;    // Cx of completely inflated parachute
double Cx_rif = 0.6;    // Cx of riffled parachute:  0.94 at 15 m/s ; 0.5 at max speed 
double Sc_S = 0.984385; //0.984385; //0.9821943; // dragging part of Sp  
double R  = 0.388239;        //0.388239; //0.485394;  // 0.1 ; 0.3;  radius of the inflated parachute  (m)
double Ro = 0.05;     // limiting radius of parachute rifling  (m)
double t0 = 0;          // start time    (s)
double t  = 0;          // current time  (s)
double t_par  = -0.1;   // the moment of unrifling the parachute   (s)
double H_open = 1600;   // the altitude of unrifling the parachute (m);
double v0 = 0; //-5.424;          // initial velocity  (m/s)       
double v  = v0;         // current velocity  (m/s)
double dt = 0.001;      // numerical integration step (s)
double dt_print = 0.5;  // cout step (s)
double h0 = 33000;       // initial altitude (m);
double h0_loc = 1680;   // altitude of launch location (m);
double n_par = 1;       // number of parachutes
// Input data }
 
// { initialization
double m = m_box + m_balloon;
double Rc = R * sqrt(1 - Sc_S);
double Rx = (2*pi*R/N_gore-pi*Ro*sin(pi/N_gore))/(2+pi*sin(pi/N_gore));
double S_par = n_par*pi * R * R * Sc_S;
double S_rif = pi * (Rx + Ro)*(Rx + Ro) + N_gore * (pi / 2) * pow(((Rx + Ro) * sin(pi / N_gore)), 2) - pi * Rc * Rc;
double dv = 0;
double dh = 0;
double h  = h0;
int n, i, Kg=0, Pg_flag = 0, result = 0;
// initialization }

// { standard atmosphere
int Hatm[] = {
  0,
  500,
  1000,
  1500,
  2000,
  2500,
  3000,
  4000,
  5000,
  6000,
  7000,
  8000,
  9000,
  10000,
  11000,
  12000,
  14000,
  16000,
  18000,
  20000,
  24000,
  28000,
  32000,
  36000
};
int Patm[] = {
  101330,
  95464,
  89877,
  84559,
  79499,
  74690,
  70123,
  61661,
  54052,
  47217,
  41106,
  35653,
  30801,
  26500,
  22700,
  19399,
  14170,
  10353,
  7565,
  5529,
  2971,
  1616,
  889,
  499
};
double Tatm[] = {
  288.2,
  284.9,
  281.7,
  278.4,
  275.2,
  271.9,
  268.7,
  262.2,
  255.7,
  249.2,
  242.7,
  236.2,
  292.7,
  223.3,
  216.8,
  216.7,
  216.7,
  216.7,
  216.7,
  216.7,
  220.6,
  224.5,
  228.5,
  239.3
};
double Beta[] = {
  0.000119267,
  0.000120614,
  0.000121985,
  0.00012341,
  0.000124796,
  0.000126191,
  0.000128599,
  0.000131705,
  0.000135193,
  0.0001386,
  0.000142321,
  0.000146286,
  0.000150402,
  0.00015478,
  0.000157143,
  0.000157047,
  0.000156925,
  0.000156872,
  0.000156763,
  0.000155277,
  0.000152236,
  0.000149403,
  0.000144373,
  0
};

// standard atmosphere }

int main(){
  while (h >= h0_loc){
    if (h0 <= H_open){  
        if (t <= t_par){
           Cx = Cx_rif;
           Sp = S_rif;
        } 
        else{
        Cx = Cx_all;
        Sp = S_par;
        //Cx = Cx_rif;
        //Sp = S_rif;
        m = m_box;
        Kg = 1;
        }
    }
    else{
        if (h >= H_open){
           Cx = Cx_rif;
           Sp = S_rif;
        } 
        else{
        Cx = Cx_all;
        Sp = S_par;
        //Cx = Cx_rif;
        //Sp = S_rif;
        m = m_box;
        Kg = 1;
        }
    }
    
    for (n = 0; n <= 24; ++n){
      if (h > Hatm[n]){
        i = n;
      }
    }

    T = Tatm[i] + (Tatm[i + 1] - Tatm[i]) * (h - Hatm[i]) / (Hatm[i + 1] - Hatm[i]);
    P = Patm[i] * exp(-Beta[i] * (h - Hatm[i]));
    rho = P * xmu / T;
    
//    mg_= Cx_all * rho * v * v * S_par / 2 / g/;
    mg_= Cx_all * rho * v * v / 4 / g;  // kg on parachute's 1 m^2

    //if(Kg == 0 && h < H_open){
    if(Kg == 1 && Pg_flag == 0){
        //m = m_box;
        Pg = Cx * rho * v * v * Sp / 2 / m / g;
        mg = Pg * m;
        Pg_flag = 1;
        cout << "Pg---" <<endl;
        //Kg = 1;
    }
    
    acc0 = acc;
    acc = (g - (Cx * Sp + k_box * Cx_box * a_box * b_box)*(rho * v * v) / 2 / m);
    //dv = (acc + acc0) * dt/2;
    dv = acc * dt;
    v0 = v;
    v = v + dv;
    dh = ((v + v0) / 2) * dt;
    h = h - dh;
    
    t = t + dt;
    if ((t - t0) >= dt_print){
    t0 = t;
    cout << t << "\t" << h << "\t" << "  " <<  v << "\t" << acc <<" \t" << P <<" \t" << mg_<<endl;
    }
    
  }
  
    cout << endl << "Pg = " << Pg << " g" << "\t" << "mg = " << mg << " kg" << "   " << Pg_flag;
    cout << endl << "Sp = " << Sp << "\t" << " Srif =" << S_rif << " m^2";

  return 0;

    
}





















