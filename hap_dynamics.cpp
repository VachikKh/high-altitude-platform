/******************************************************************************

This code simulates weather balloon ascent using the data of Standard Atmosphere. 

*******************************************************************************/
#include <iostream>
#include <cmath>
#include <string>

using namespace std;

// Constants

    double const pi   = 3.1415926535;    // Pi
    double const R_id = 8314.462;        // gas constant                                   (J/K/mol*1000)
    double const pix  = pi/180;          // coverts degree to radian
    double const Ksph = 4.0*pi/3;        // sphere volume coefficient 

    double const g    = 9.80665;         // free fall acceleration at release location     (m/s^2)
    double const k_Sp = 0.739668778;     // coefficient of pumpkin shape area
    double const k_Vp = 1.21852421611856;// coefficient of pumpkin shape volume
    double const k_p_rate = (5/3)/100000;// converts pumpimg rate unit from l/min to m^3/s

// Balloon and ballonet data
    double ro_f = 920;                   // density of the balloon film                    (kg/m^3)
    double d_f  = 0.0000254;             // thickness of the balloon film                  (m)
    double R_b  = 6.122831782671;        // radius of fully the inflated balloon pumpkin   (m)
    double R_b_x;                        // current drag radius of the partly inflated balloon     (m)
    double R_bl = 2;                     // ballonet pumpkin shape radius                  (m)
    double So_b;                         // 1/2 of the balloon pumpkin meridian length     (m)
    double So_bl;                        // 1/2 of the ballonet pumpkin meridian length    (m)
    double S_b  = pi*R_b*pi*R_b;         // balloon surface area                           (m^2)
    double S_bl = pi*R_bl*pi*R_bl;       // ballonet surface area                          (m^2)
    double Sb_drag;                      // balloon drag area ( = pi*R_b_x^2 )             (m^2)
    double V_b;                          // balloon volume                                 (m^3)
    double V_b_max;                      // volume of fully inflated balloon pumpkin       (m^3)
    double V_bl;                         // ballonet volume                                (m^3)
    double V_bl_max;                     // volume of fully inflated ballonet pumpkin      (m^3)
    double m_b;                          // balloon mass                                   (kg)
    double m_bl;                         // ballonet mass                                  (kg)
//    double Sb0;                          // initial drag area of the balloon               (m^2)
//    double Sb;                           // current drag area of the balloon               (m^2)
    double rb;                           // current radius of the partly inflated balloon's top part         (m) 
    double Cd_z_up;                      // drag cefficient for ascend
    double Cd_z_down;                    // drag cefficient for descent
    double Cd_z;                         // drag cefficient along z_0
    double Cd_z_max = 0.47;              // max drag cefficient along z_0
    double Cd_xy;                        // drag cefficient for horizontal motion
    double eta;                          // semi-vertex angle of the partly inflated balloon at its bottom  (rad) 
    double eta0;                         // lower limit of the search range for eta        (rad)    
    double eta1;                         // uper  limit of the search range for eta        (rad) 
    double eta_max = 31.2*pix;           // max eta, beyond wich Cd_z does not change      (deg*pix -> rad)
    double eta_error;                    // accuracy of defining eta                       (%)
    double V_b_eq;                       // equation for eta
    double P_bl;                         // air pressure in the ballonet                   (Pa)
    double dP;                           // dP = P_gas - P_atm ; balloon  diff pressure    (Pa)
    double dP_bl;                        // dP = P_bl  - P_gas ; ballonet diff pressure    (Pa)    
    double dP_max    =  100;             // balloon  max diff pressure = P_b-P_atm         (Pa)
    double dP_max_bl = 7500;             // ballonet max diff pressure = P_bl-P_gas        (Pa)
    double m_bl_air   = 0.0;             // current air mass in the ballonet               (kg)   0.9
    double m_bl_air_1;                   // ballonet air mass at which P_gas = Patm        (kg)
    double m_bl_air_2;                   // ballonet air mass at which P_bl = P_gas_max    (kg)
    double V_bl_air_1;                   // ballonet volume at m_bl_air = m_bl_air_1       (m^3)
    double T_bl_air;                     // air temperature in the ballonet                (K)
    double K_b_bl;                       // auxiliary coefficient for defining V_bl
    int diff_Pmode = 0;                  // differential pressure mode; values: 1/2/3/4
    bool diff_1;                         // flag for diff pressure mode 1  
    bool diff_2;                         // flag for diff pressure mode 2  
    bool diff_3;                         // flag for diff pressure mode 3
    bool diff_4;                         // flag for diff pressure mode 4
  
// Covering net data (net with diamond-shaped cells)
    double F_tendon_b =  4;              // breaking strength of the balloon  tendon       (kg)
    double F_tendon_bl= 20;              // breaking strength of the ballonet tendon       (kg)
    double k_sf       =  2;              // safety factor
    double L_h;                          // horizontal size of a mesh                      (m)
    double L_v;                          // verticlal  size of a mesh                      (m)
    double n_h;                          // number of cells along the equator
    double n_v;                          // number of cells along the meridian    
    double Ad         = 0.38;            // aspect ratio of a diamond-shaped mesh (L_h/L_v)
    double F_mesh;                       // force on the mesh film caused by diff pressure (kg)
    double n_meshes;                     // number of net mehses
    double m_tendons;                    // total mass of tendons of the net               (kg)

// Payload data
    double m_p  = 8.3;                   // payload                                        (kg)
    double S_p  = 0.16;                  // payload drag area                              (m^2)
    double Cd_p = 0.25;                  // payload drag cefficient
    
// Pump data
    int N_pump          = 6;             // number of pumps                                
    int N_pumps_on      = 3;             // number of pumps running at a time
    double m_pump       = 0.544;         // weight of one pump                             (kg)
    double P_diff_pump  = 8000;          // maximum differential pressure of the pump      (Pa)
    double m_pump_all;                   // weight of all pumps                            (kg)
    double pump_rate_V  = 600;           // air volume pumping rate of one pump            (l/min)
    double pumping_rate_V;               // overall air volume pumping rate                (m^3/s)
    double pumping_rate_m;               // overall air mass pumping rate                  (kg/s)
    double k_pump;                       // pump gain - due to more powerfull BLDC motor
    double S_pump       = pi*0.01*0.01;  // pump outlet cross section area                 (m^2)
    int pump_on_off     = 0;             // 0 - off, 1 - on
    
// LTA gas data
    double mu_gas;                       // molar mass of the LTA gas                      (gram)
    int k_He_H2 = 1;                     // LTA gas: 1 - helium , 0 - hydrogen
    double xmu_gas;                      // R/mu ratio for the LTA gas 
                                         // ( P = (ro_gas/mu)*R*T = (R/mu)*ro_gas*T )
    double m_gas_0;                      // initial LTA gass mass                          (kg)    // = 0.569;
    double m_gas;                        // working LTA gass mass                          (kg)    //  = 3.414515;
    double dT_gas  = 5.0;                // additional temperature of the LTA gas due to greenhouse effect
    double V_gas;                        // volume of the LTA gas                          (m^3)
    double P_gas;                        // pressure of the LTA gas                        (Pa)
    double T_gas;                        // temperature of the LTA gas                     (Pa)
    double P_gas_max;                    // max pressure of the LTA gas at V_bl = V_bl_max (Pa)

// Atmospheric data and variables
    double mu_air  = 28.966;             // air molar mass                                 (gram)
    double xmu_air = R_id / mu_air ;     // mu/R ratio for air    (  P = (rho/mu)RT  )
    double P_atm;                        // air pressure at the current altitude           (Pa)
    double T_atm;                        // air temperature at the current altitude        (K)
    double rho_atm;                      // air density at the current altitude            (kg/m^3)
    double P_atm_hmax;                   // air pressure at maximum altitude               (Pa)
    double T_atm_hmax;                   // air temperature at at maximum altitude         (K)
    double Vw_xy0  = 0;                  // initial horizontal velocity of wind            (m/s)
//  double Vw_z0   = 0;                  // initial vertival velocity of wind              (m/s)
    double Vw_xy;                        // current horizontal velocity of wind            (m/s)  
//  double Vw_z;                         // current vertical velocity of wind              (m/s)
    double dA_dh   = 0.4*pix;            // wind azimuth change per km of altitude         (deg/km*pix -> rad/km)
    double A0      = 50*pix;             // initial wind azimuth at h0 (ref: North vector) (deg*pix -> rad)
    double A;                            // current wind azimuth at h0 (ref: North vector) (rad)
    double dVw_dh  = 0.77;               // change of Vw_xy per km of altitude increase    (m/s/km)        


// Overload    
    double Pg0;                          // previuos overload                              (g)
    double Pg;                           // current overload                               (g)
    double Pg_max;                       // maximum overload                               (g)
    double mg;                           // current overload force                         (kg)
    double mg_max;                       // maximum overload force                         (kg)

// Declaration and initialization of variables 
    double t0   = 0;                     // previous time                                  (s)
    double t    = 0;                     // current time                                   (s)
    double ad_x;                         // air drag acceleration along x                  (m/s^2)
    double ad_y;                         // air drag acceleration along y                  (m/s^2)
    double ad_z;                         // air drag acceleration along z                  (m/s^2)
    double vx_0 = 0;                     // previous velocity along x                      (m/s)
    double vy_0 = 0;                     // previous velocity along y                      (m/s) 
    double vz_0 = 0;                     // previous velocity along z                      (m/s) 
    double vx   = 0;                     // current velocity along x                       (m/s)
    double vy   = 0;                     // current velocity along y                       (m/s) 
    double vz   = 0;                     // current velocity along z                       (m/s) 
    double v_0  = 0;                     // previous velocity                              (m/s)
    double v    = 0;                     // current velocity                               (m/s)
    double h0   = 15000; //1120;  24000           // initial altitude                 ( < 36000 )   (m)
    double h_max= 25000;                 // the max altitude to stop ascent  ( < 36000 )   (km)
    double h_min= 15000;                 // the min altitude to stop descent ( < 36000 )   (km)
    double x_0  = 0;                     // previous x coordinate (East)                   (m/s)
    double y_0  = 0;                     // previous y coordinate (North)                  (m/s) 
    double z_0  = h0;                    // initial/previous z coordinate (altitude)       (m/s) 
    double x    = 0;                     // current x coordinate  ( -> East  )             (m/s)
    double y    = 0;                     // current y coordinate  ( -> North )             (m/s) 
    double z    = h0;                    // current z coordinate  (altitude)               (m/s) 
    double r    = 0;                     // current distance                               (m)
    double a_x;                          // acceleration along x                           (m/s^2)
    double a_y;                          // acceleration along y                           (m/s^2)
    double a_z;                          // acceleration along z                           (m/s^2)
    double mu_;                          // gas/air molar mass ratio
    double dt   = 0.001;                // numerical integration step                     (s)
    double dt_print = 1;                // cout step                                      (s)
    double eps  = 1;                     // error % for solving the equation for eta       (%)
    double Fa;                           // buoyant force                                  (N)
    double Fg;                           // gravitational force                            (N)
    double Fd_z;                         // air drag force for vertical motion             (N)

    double m_net;                        // overall net mass                   (kg)
    double m;                            // overall mass                       (kg)
    double dvx      = 0;                 // x increnment                       (m)
    double dvy      = 0;                 // y increnment                       (m)
    double dvz      = 0;                 // z increnment (altitude)            (m)
    double dv       = 0;                 // velocity increnment                (m/s)
    double dx       = 0;                 // x increnment                       (m)
    double dy       = 0;                 // y increnment                       (m)
    double dz       = 0;                 // z increnment (altitude)            (m)
    int n;                               // number of atmospհeric layers in the standard model 
    int i;                               // index of the current atmospհeric layer
    int i_max;                           // index of the atmospհeric layer for the max altitude
    int j = 1;
    
    double u[5]    = {0,0,1./2,1./2,1};
    double k_vy[5] = {0};
    double k_vx[5] = {0};
    double k_vz[5] = {0};
    double k_y[5]  = {0};
    double k_x[5]  = {0};
    double k_z[5]  = {0};
    double Vx_temp[5] = {0};
    double Vy_temp[5] = {0};
    double t_temp[5]  = {0};
    double x_temp[5]  = {0};
    double y_prev     =  0;
    double y_temp[5]  = {0};
    double m_temp[5]  = {0};
    double ms_temp[5] = {0};
    double W_temp[5]  = {0};
    double F_temp[5]  = {0};
    double k_Imp[5]   = {0};

    
// { Standard Atmosphere

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
/* double Aw[] = {
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
double Vw[] = {
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
*/

// Standard Atmosphere }

// Input data }

int main() {

    Vw_xy = Vw_xy0;                    
    A     = A0;        

    So_b  = k_Sp * pow(pi,0.5) * R_b; 
    So_bl = k_Sp * pow(pi,0.5) * R_bl; 
    V_b_max  = k_Vp * pow(So_b,3);
    V_bl_max = k_Vp * pow(So_bl,3); 
    m_b  = ro_f * S_b  * d_f;
    m_bl = ro_f * S_bl * d_f;
    //   rb  = 4 * So_b / (pi + 2 / sin(eta)); 
    m_pump_all      = N_pump * m_pump;
   // pumping_rate_V  = k_pump * N_pumps_on * pump_rate_V;
    
    m_net    = m_p + m_b + m_bl + m_pump_all; // overall net mass                   (kg)

    if (k_He_H2 = 1){
        mu_gas = 4;
        }
    else{
        mu_gas = 2;
        }
        
    mu_  = mu_gas/mu_air;
    xmu_gas = R_id / mu_gas;
        
    for (n = 0; n <= 24; ++n){
      if (h_max > Hatm[n]){
        i_max = n;
      }
    }
    
    T_atm_hmax = Tatm[i_max] + (h_max - Hatm[i_max]) * (Tatm[i_max + 1] - Tatm[i_max]) / (Hatm[i_max + 1] - Hatm[i_max]);    
    P_atm_hmax = Patm[i_max] * exp(-Beta[i_max] * (h_max - Hatm[i_max]));
    /*  cout << "i_max=" << i_max <<endl;
        cout << "T_atm_hmax=" << T_atm_hmax <<endl;
        cout << "P_atm_hmax=" << P_atm_hmax <<endl;
        cout << "Hatm[i_max]=" << Hatm[i_max] <<endl;
        cout << "Hatm[i_max+1]=" << Hatm[i_max+1] <<endl;
        cout << "Tatm[i_max]=" << Tatm[i_max] <<endl;
        cout << "z=" << z <<endl;
    */
    
    while (z <= h_max){
        for (n = 0; n <= 24; ++n){
          if (z > Hatm[n]){
            i = n;
          }
        }
    
        T_atm = Tatm[i] + (z - Hatm[i]) * (Tatm[i + 1] - Tatm[i]) / (Hatm[i + 1] - Hatm[i]);
        P_atm = Patm[i] * exp(-Beta[i] * (z - Hatm[i]));
        rho_atm = P_atm / xmu_air / T_atm;
        pumping_rate_m = rho_atm * pumping_rate_V;
        T_gas = T_atm + dT_gas;
        
        if (z < h_min){
            m_gas = m_gas_0;
            }
        else{
            m_gas = V_b_max * P_atm_hmax / xmu_gas / (T_atm_hmax+dT_gas);
            }
        
        P_gas_max = m_gas * xmu_gas * T_gas / (V_b_max - V_bl_max);
        
        T_bl_air   = T_atm;
        // V_bl_air_1 = V_b_max - m_gas * xmu_gas * T_gas / P_atm;
        // m_bl_air_1 = P_atm * V_bl_air_1 / xmu_air / T_bl_air;
        
        m_bl_air_1 = P_atm * (mu_air/R_id/T_atm) * (V_b_max-(m_gas/mu_gas)*R_id*(T_atm+dT_gas)/P_atm);
        // m_bl_air_2 = IF(P_atm<P_gas_max,m_gas*(mu_air/mu_gas)*V_bl_max/(V_b_max-V_bl_max)*(T_atm+dT)/T_atm,P_atm*V_bl_max*mu_air/R_id/T_atm)
        if(P_atm < P_gas_max) {
            m_bl_air_2 = (m_gas / mu_) * (V_bl_max / (V_b_max - V_bl_max)) * T_gas / T_bl_air;
        }
        else {
            m_bl_air_2 = P_atm * V_bl_max / xmu_air / T_bl_air;
        }
        //cout << "m_bl_air_1=" << m_bl_air_1 <<endl;
        //cout << "m_bl_air_2=" << m_bl_air_2 <<endl;
        
        if (z > h_max) {
            cout <<" m_bl_air = ";   
            cin >> m_bl_air;   
            }
        
        K_b_bl = (m_bl_air/m_gas)*(mu_gas/mu_air)*(T_atm/T_gas);
        
        diff_1 = false;
        diff_2 = false;
        diff_3 = false;
        diff_4 = false;
       
       //  IF(OR(AND(P_atm<P_gas_max,m_bl_air<mI_bl_air),AND(P_atm>=P_gas_max,m_bl_air<=mII_bl_air),AND(P_atm<P_gas_max,m_bl_air=mII_bl_air)),1,0)
       
        if (((P_atm < P_gas_max) && (m_bl_air <= m_bl_air_1)) || ((P_atm >= P_gas_max) && (m_bl_air <= m_bl_air_2)) || ((P_atm < P_gas_max) && (m_bl_air == m_bl_air_2))) {  
            diff_Pmode = 1;
            diff_1 = true;
        } 
            
        //  IF(AND(P_atm>=P_gas_max,m_bl_air>mII_bl_air),1,0)
            
        if (((P_atm >= P_gas_max) && (m_bl_air > m_bl_air_2))) {  
            diff_Pmode = 2;
            diff_2 = true;
        } 
            
        //  IF(AND(P_atm<P_gas_max,(m_bl_air-mI_bl_air)>0,(m_bl_air-mII_bl_air)<=0),1,0)
            
        if (((P_atm < P_gas_max) && ((m_bl_air - m_bl_air_1) > 0)) && ((m_bl_air-m_bl_air_2) <= 0)) {  
            diff_Pmode = 3;
            diff_3 = true;
        }  
            
        //  IF(AND(P_atm<P_gas_max,(m_bl_air-mI_bl_air)>0,(m_bl_air-mII_bl_air)<=0),1,0)
            
        if ((P_atm < P_gas_max) && (m_bl_air > m_bl_air_2)) {  
            diff_Pmode = 4;
            diff_4 = true;
        }  
            
        if (diff_1 + diff_2 + diff_3 + diff_4 != 1) {
          cout << "j=" << j <<endl;
            cout << "P_atm=" << P_atm <<endl;
            cout << "T_atm=" << T_atm <<endl;
            cout << "mu_=" << mu_ <<endl;
            cout << "m_gas=" << m_gas <<endl;
            cout << "V_b_max=" << V_b_max <<endl;
            cout << "P_gas_max=" << P_gas_max <<endl;
            cout << "m_bl_air_1=" << m_bl_air_1 <<endl;
            cout << "m_bl_air_2=" << m_bl_air_2 <<endl;
            cout << "diff_Pmode=" << diff_Pmode <<endl;
            cout << "diff_1=" << diff_1 <<endl;
            cout << "diff_2=" << diff_2 <<endl;
            cout << "diff_3=" << diff_3 <<endl;
            cout << "diff_4=" << diff_4 <<endl; 
            cout << "Mode failure";
            exit(0);
        }  
            //cout << "diff_Pmode" << diff_Pmode <<endl;
            //j++;
        
            
        if (diff_1) {
            V_b   = ((m_gas / mu_gas) * (R_id * T_gas / P_atm) + V_bl);
            V_bl  = (m_bl_air / mu_air) * (R_id *T_atm / P_atm);
            P_gas = P_atm;
            P_bl  = P_atm;
        }
            
        if (diff_2) {
            V_b   = ((m_gas / mu_gas) * (R_id * T_gas / P_atm) + V_bl_max);
            V_bl = V_bl_max;
            P_gas = P_atm;
            P_bl  = (m_bl_air / mu_air) * (R_id * T_atm / V_bl_max);
        }
        
        if (diff_3) {
            V_b   = V_b_max;
            V_bl  = V_b_max * K_b_bl / (K_b_bl + 1);
            P_gas = (m_gas / mu_gas) * R_id * T_gas / (V_b_max - V_bl);
            if (V_bl == 0) P_bl = P_atm;
            else P_bl = (m_bl_air / mu_air) * (R_id * T_atm / V_bl);
        }
            
            
        if (diff_4) {
            V_b   = V_b_max;
            V_bl  = V_bl_max;
            P_gas = P_gas_max;
            P_bl  = (m_bl_air / mu_air) * R_id * T_atm / V_bl_max;
        }
            
        dP = P_gas - P_atm;
        /*
        if (dP > max_dP) {
            cout << "Balloon burst! Q -> @";
            exit(0);
        }
        */
        
        dP_bl = P_bl  - P_gas;
        /*
        if (dP_bl > max_dP_bl) {
            cout << "Ballonet burst! o -> @";
            exit(0);
        }
           */
        // Solve the following equation to get eta
        // V_b = (pi/3)*(2+1/tan(eta)*(4/(pi+2/sin(eta)))^3*So^3
        eta0 =  1*pix;
        eta1 = 32*pix;
        eta_error = (abs((eta1-eta0)/eta0))*100;
        while (eta_error >= eps) {
               eta = (eta0 + eta1) / 2;
               V_b_eq = ((pi/3)*(2+1/tan(eta0))*pow(4*So_b/(pi+2/sin(eta0)),3)-V_b)*((pi/3)*(2+1/tan(eta))*pow(4*So_b/(pi+2/sin(eta)),3)-V_b);
            //   cout << "V_b_eq=" << V_b_eq <<endl;
            //   cout << "eta=" << eta <<endl;
               if (V_b_eq < 0) eta1 = eta;
               else eta0 = eta;
               eta_error = (abs((eta1-eta0)/eta0))*100;
        }
              
        R_b_x = 4 * So_b / (pi + 2 / sin(eta)); 
        if (R_b_x > R_b) R_b_x = R_b;
        Sb_drag  = pi * R_b_x * R_b_x;
              
        if (vz < 0) {
            if (eta <= eta_max) {
                Cd_z = 1 - pow(cos(eta),4); 
            }
            else Cd_z = Cd_z_max;    
        }
        else Cd_z = Cd_z_max;
        
        /*
           while ((b-a)>=eps)
              {
                   c=(a+b)/2;
                   if (f(a)*f(c)<0)
                   b=c;
                   else a=c;
        */
    
        /*   Rb  = Rb0*pow(((P0+547.7)/(P+162)), (1.0/3));
            Sb  = pi * Rb * Rb;
            Vb  = Ksph * Rb * Rb * Rb;
            Nz  = rho * Vb - (m_balloon + m_gas);
            Vas = sqrt( 2 * (rho * Vb - m) * g / Cx / rho / Sb);
        */
        
        // Pg0 = Cx * rho * v * v * Sp / 2 / m / g;
        // ax = Cx * (rho * v * v/ 2) * Sp / m;
        m    = m_net + m_gas + m_bl_air;
        Fa   = rho_atm * V_b * g;
        Fg   = m * g;
        Fd_z = ((vz > 0) - (vz < 0))* Cd_z * (rho_atm * vz * vz / 2) * Sb_drag;
        a_z  = (Fa - Fg - Fd_z) / m;
        dvz  = a_z * dt;
        vz_0 = vz;
        vz   = vz + dvz;
        dz   = ((vz + vz_0) / 2) * dt;
        z    = z + dz;
        
        Pg   = a_z / g;
        mg   = Pg * m;
        
        /*
        if(Pg > Pg_max) {
            Pg_max = Pg;
            mg_max = Pg * m;
        }
            
        if(Kg == 1 && Pg_flag == 0){
            Pg_flag = 1;
            cout << "Parachute inflated" <<endl;
            cout << t << "\t" << h << "\t" << v << "\t" << acc << "\t" << P << "\t" << Cx << "\t" << Kg << "\t" << Pg <<endl;
        }  
        */
        
        t = t + dt;
        if ((t - t0) >= dt_print){
            t0 = t;
            cout << t << "\t" << z << "\t" << vz << "\t" << a_z << "\t" << Pg << "\t" << P_atm << "\t" << T_atm << "\t" << rho_atm << "\t" << m_bl_air << "\t" << diff_Pmode <<"\t" << dP << "\t" << dP_bl << "\t" << m_gas << "\t" << eta/pix <<  "\t" << Cd_z <<  "\t" << m << endl;
            //cout << Rb << "\t" << Sb << "\t" << Vb << "\t" << rho << "\t" << m*g << "\t" << Cx * (rho * v * v / 2) * Sb <<endl;
            //cout << h << "\t" << v0 << "\t" << v << "\t" << acc << "\t" << dh << "\t" << m <<endl;
        }
   }
  
   // cout << endl << "Sb0 = " << Sb0 << " m^2" << "\t" << "Vb0 = " << Vb0 << " m^3" << "   " << Pg_flag;
   // cout << endl << "Nz = " << Nz << " m^2" << "\t" << "Spar = " << S_par << "\t" << " Srif =" << S_rif << " m^2\t" << " dh =" << dh << " m";

  return 0;

    
}






















