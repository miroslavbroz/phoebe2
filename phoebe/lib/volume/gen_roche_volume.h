#if !defined(__volume_h)
#define __volume_h

#include <iostream>
#include <cmath>
#include <limits>
#include <fstream>


/*
  Calculating volume and area of generalized Roche lobes implicitely 
  defined by
   
    Omega(x,y,z) = Omega_0
  
  where Omega is generalized Kopal potential

    Omega = 
    1/rho 
    + q [(delta^2 + rho^2 - 2 rho lambda delta)^(-1/2) - rho lambda/delta^2]
    + 1/2 F^2(1 + q) rho^2 (1 - nu^2)

  with position in spherical coordinates is given as

  x = rho lambda      lambda = sin(theta) cos(phi)
  y = rho  mu         mu = sin(theta) sin(phi)
  z = rho nu          nu = cos(theta)

  Author: Martin Horvat,  April 2016 
*/


namespace gen_roche {
  
  /*
    Computing area of the surface and the volume of the Roche lobes
    intersecting x-axis
    
      {x0,0,0}  {x1,0,0}
    
    The range on x-axis is [x0, x1].
  
    Input:
      x_bounds[2] = {x0,x1} 
      Omega0 
      q - mass ratio M2/M1
      F - synchronicity parameter
      delta - separation between the two objects
      m - number of steps in x - direction
      choice -
        1  - area
        2  - volume
        3  - both
      polish - if true than after each RK step 
                we perform reprojection onto surface
                 
    Using: Integrating surface in cylindric geometry
      a. Gauss-Lagrange integration in phi direction
      b. RK4 in x direction
    
    Precision:
      Relative precision should be better than 1e-4.
      
    Stability:
      It works for overcontact and well detached cases, but could be problematic
      d(polar radius)/dz is too small.
       
    Output:
      av[2] = {area, volume}
  
    Ref: 
      * https://en.wikipedia.org/wiki/Gaussian_quadrature
      * https://en.wikipedia.org/wiki/Gauss–Kronrod_quadrature_formula
      * http://mathworld.wolfram.com/LobattoQuadrature.html <-- this would be better  
  */
  
  template<class T> 
  void area_volume(
    T av[2],
    T xrange[2],
    const T & Omega0,
    const T & q,
    const T & F = 1,
    const T & delta = 1,
    const int & m = 1 << 14,
    const unsigned choice = 3,
    const bool polish = false)
  {
    
    #if 1
    const int n = 10;

    // cos()^2 of Gauss-Lagrange nodes x_i in [0, Pi/2] 
    const T c[n]={0.999580064408504868392251519574,0.988810441236620370626267163745,0.93792975088501989705468221942,0.814698074016316269208833703789,0.615862835349225523334990208675,0.384137164650774476665009791325,0.185301925983683730791166296211,0.06207024911498010294531778058,0.011189558763379629373732836255,0.0004199355914951316077484804259};

    const T d[n]={0.0004197592455941272416985899940,0.011064352538060503513192603217,0.058217533289784414692374340457,0.150965122210421127009962838622,0.236575803384838256502140293303,0.236575803384838256502140293303,0.150965122210421127009962838622,0.058217533289784414692374340457,0.011064352538060503513192603217,0.0004197592455941272416985899940};
        
    // Gauss-Lagrange weights * Pi
    const T weights[n]={0.209454205485130325204687978,0.46951526056054717731928526,0.68828010698191943974759705,0.845926347240509469003431984,0.92841673332168682718764111324,0.92841673332168682718764111324,0.845926347240509469003431984,0.68828010698191943974759705,0.46951526056054717731928526,0.209454205485130325204687978};
   
    #else
    const int n = 15;
    
    // cos()^2 of Gauss-Lagrange nodes x_i in [0, Pi/2] 
    const T c[n]={0.999911065396171632736007253518,0.997574886997681325232904530099,0.985854212895255771578898079594,0.953879938591329487419660590834,0.890692147886632646280933330332,0.790164039217058712341619458721,0.655400162984313217687403793532,0.5,0.344599837015686782312596206468,0.209835960782941287658380541279,0.109307852113367353719066669668,0.046120061408670512580339409166,0.014145787104744228421101920406,0.0024251130023186747670954699007,0.0000889346038283672639927464819};
    
    const T d[n]={0.0000889266944646091553555375073,0.0024192318292446596704491832554,0.013945683811931480320682058875,0.04399300134433097347511477941,0.097359645579729565862303734566,0.165804830345241213606518615509,0.225850789344448888056399859737,0.25,0.225850789344448888056399859737,0.165804830345241213606518615509,0.097359645579729565862303734566,0.04399300134433097347511477941,0.013945683811931480320682058875,0.0024192318292446596704491832554,0.0000889266944646091553555375073};
        
    // Gauss-Lagrange weights * Pi
    const T weights[n]={0.0966141591290711189794453,0.22106145785079100848207,0.336650619784076362391549,0.4384742164293535129419684,0.5223501155128774575517226,0.584842030033819626900854873,0.6233908965445645786498187508,0.636418316610479145130430187624,0.6233908965445645786498187508,0.584842030033819626900854873,0.5223501155128774575517226,0.4384742164293535129419684,0.336650619784076362391549,0.22106145785079100848207,0.0966141591290711189794453};
    #endif
    
    const int dim = n + 2;
       
    T d2 = delta*delta,
      d3 = d2*delta,
      a = d3*F*F,
      b = a*(1 + q),
      omega = delta*Omega0,
      
      t = xrange[0]/delta, 
      dt = (xrange[1] - xrange[0])/(m*delta),
      y1[dim], y[dim], k[4][dim], w[n];
    
    
    //
    // What is calculated
    //
    av[0] = av[1] = 0;
    
    bool 
      b_area = (choice & 1u) == 1u,
      b_volume = (choice & 2u) == 2u;
    
    if (!b_area && !b_volume) return;
       
    //
    // Integration over the surface with RK4
    //  Def: 
    //    y = R^2, F = dR/dt, G = dR/dphi, 
    //    x = delta t, c = cos(phi)^2
    //  Eq:
    //    y_i(t)' = F(t, y_i(t), c_i) i = 0, .., n -1
    //    A'(t) = delta^2 sqrt(R + G^2/(4R) + F^2/4) df dt
    //    V'(t) = delta^3 1/2 R df dx
   
    for (int i = 0; i < n; ++i) w[i] = dt*weights[i];
      
    // init point
    for (int i = 0; i < dim; ++i) y[i] = 0;
    
    for (int j = 0; j < m; ++j){
      
      {
        T t1, f, g, f1, f2, s, s1, s2; // auxiliary variables
        
        //
        // RK iteration
        //
        
        // 1. step
        k[0][n] = k[0][n+1] = 0;
        s1 = t*t, s2 = (t-1)*(t-1);
        for (int i = 0; i < n; ++i) {
          s = y[i];
          
          // f1 = (R + t^2)^(-3/2)
          f1 = s + s1;
          f1 = 1/(f1*std::sqrt(f1));
          
          // f1 = (R + (t-1)^2)^(-3/2)
          f2 = s + s2;
          f2 = 1/(f2*std::sqrt(f2));
          
          // k0 = -dx F_t/F_R
          g = 1/(b*c[i] - f1 - q*f2);
          f = (q*(1 + (t - 1)*f2) + t*(f1 - b))*g; // = 1/2 dR/dt, x = delta*t
          k[0][i] = dt*2*f;
          
          if (b_area) {    
            g *= b; // = -(dR/dc)/R   Note: dR/dphi = -dR/dc*2*sqrt(c(1-c)) 
            k[0][n] += w[i]*std::sqrt(s*(1 + g*g*d[i]) + f*f); // = dA
          }
         
          if (b_volume)  k[0][n+1] += w[i]*s;  // = dV/2
        }
        
       
        // prepare: y1 = y + k0/2
        for (int i = 0; i < n; ++i) y1[i] = y[i] + k[0][i]/2;
        
        // 2. step
        t1 = t + dt/2;
        
        k[1][n] = k[1][n+1] = 0; 
        s1 = t1*t1, s2 = (t1-1)*(t1-1);
        for (int i = 0; i < n; ++i){
          s = y1[i];
          
          // f1 = (R + t1^2)^(-3/2)
          f1 = s + s1;
          f1 = 1/(f1*std::sqrt(f1));
          
          // f2 = (R + (t1-1)^2)^(-3/2)
          f2 = s + s2;
          f2 = 1/(f2*std::sqrt(f2));
          
          // k0 = -dx F_t/F_R
          g = 1/(b*c[i] - f1 - q*f2);
          f = (q*(1 + (t1 - 1)*f2) + t1*(f1 - b))*g; // = 1/2 dR/dt, x = delta*t
          k[1][i] = dt*2*f;
        
          if (b_area) {
            g *= b; // = -(dR/dc)/R   Note: dR/dphi = -dR/dc*2*sqrt(c(1-c)) 
            k[1][n] += w[i]*std::sqrt(s*(1 + g*g*d[i]) + f*f);  // = dA
          }
          
          if (b_volume) k[1][n+1] += w[i]*s; // = dV/2
          
        }
        
        // prepare: y1 = y + k1/2
        for (int i = 0; i < n; ++i) y1[i] = y[i] + k[1][i]/2;
        
        // 3. step
        k[2][n] = k[2][n+1] = 0;
        for (int i = 0; i < n; ++i) {
          s = y1[i];
          
          // f1 = (R + t1^2)^(-3/2)
          f1 = s + s1;
          f1 = 1/(f1*std::sqrt(f1));
          
          // f2 = (R + (t1-1)^2)^(-3/2)
          f2 = s + s2;
          f2 = 1/(f2*std::sqrt(f2));
          
          // k0 = -dx F_t/F_R
          T g = 1/(b*c[i] - f1 - q*f2);       
          T f = (q*(1 + (t1 - 1)*f2) + t1*(f1 - b))*g; // = 1/2 dR/dt, x = delta*t
          k[2][i] = dt*2*f;
                  
          if (b_area) {        
            g *= b; // = -(dR/dc)/R   Note: dR/dphi = -dR/dc*2*sqrt(c(1-c))
            k[2][n] += w[i]*std::sqrt(s*(1 + g*g*d[i]) + f*f); // = dA
          }
          if (b_volume) k[2][n+1] += w[i]*s; // = dV/2
        }
        
      
        // y1 = y + k2
        for (int i = 0; i < n; ++i) y1[i] = y[i] + k[2][i];
        
        // 4. step
        t1 = t + dt;
        k[3][n] = k[3][n+1] = 0;
        s1 = t1*t1, s2 = (t1-1)*(t1-1);
        for (int i = 0; i < n; ++i){
          s = y1[i];
          
          // f1 = (R + t1^2)^(-3/2)
          f1 = s + s1;
          f1 = 1/(f1*std::sqrt(f1));
          
          // f2 = (R + (t1-1)^2)^(-3/2)
          f2 = s + s2;
          f2 = 1/(f2*std::sqrt(f2));
          
          // k0 = -dx F_t/F_R
          g = 1/(b*c[i] - f1 - q*f2);
          f = (q*(1 + (t1 - 1)*f2) + t1*(f1 - b))*g; // = 1/2 dR/dt, x = delta*t
          k[3][i] = dt*2*f;
          
          if (b_area) {       
            g *= b; // = (dR/dc)/R   Note: dR/dphi = dR/dc*2*sqrt(c(1-c))
            k[3][n] += w[i]*std::sqrt(s*(1 + g*g*d[i]) + f*f); // = dA
          }
          
          if (b_volume) k[3][n+1] += w[i]*s; // = dV/2
        }
      }
            
      for (int i = 0; i < dim; ++i)
        y[i] += (k[0][i] + 2*(k[1][i] + k[2][i]) + k[3][i])/6;  
      
      t += dt;
      
      //
      // Polishing solutions with Newton-Rapson iteration
      //
      if (polish){
        
        const int it_max = 10;
        const T eps = 10*std::numeric_limits<T>::epsilon();
        const T min = 10*std::numeric_limits<T>::min();
        
        int it;
        
        T s1 = t*t, s2 = (t - 1)*(t - 1),
          s, f, df, ds, f1, f2, g1, g2;
          
        for (int i = 0; i < n; ++i) {
          
          s = y[i];

          it = 0;
          do {
            
            // g1 = (R + (t-1)^2)^(-1/2), f1 = (R + t^2)^(-3/2)
            f1 = 1/(s + s1);
            g1 = std::sqrt(f1);
            f1 *= g1;
            
            // g2 = (R + (t-1)^2)^(-1/2),  f2 = (R + (t-1)^2)^(-3/2)
            f2 = 1/(s + s2);
            g2 = std::sqrt(f2);
            f2 *= g2;
             
            df = -(b*c[i] - f1 - q*f2)/2;                     // = dF/dR
            f = omega - q*(g2 - t) - g1 - b*(c[i]*s + s1)/2;  // =F
            
            s -= (ds = f/df);
          
          } while (std::abs(f) > eps && 
                  std::abs(ds) > eps*std::abs(s) + min && 
                  ++it < it_max);
          
          if (!(it < it_max)) 
            std::cerr << "Volume: Polishing did not succeed\n";
        
          y[i] = s;
        }
      }
    }   
    
    if (b_area) av[0] = d2*y[n];
    if (b_volume) av[1] = d3*y[n+1]/2;
  }

  /*
    Computing area of the surface and the volume of the Roche lobes
    intersecting x-axis
    
      {x0,0,0}  {x1,0,0}
    
    The range on x-axis is [x0, x1].
  
    Input:
      x_bounds[2] = {x0,x1} 
      Omega0 
      q - mass ratio M2/M1
      F - synchronicity parameter
      delta - separation between the two objects
      m - number of steps in x - direction
      polish - if true than after each RK step 
                we perform reprojection onto surface
                
    Using: Integrating surface in cylindric geometry
      a. Gauss-Lagrange integration in phi direction
      b. RK4 in x direction to get the surface points
      c. Romberg integration scheme to compute the surface and volume
    
    Precision:
      Relative precision should be better than 1e-4.
      
    Stability:
      It works for overcontact and well detached cases, but could be problematic
      d(polar radius)/dz is too small.
       
    Output:
      av[2] = {area, volume}
  
    Ref: 
      * https://en.wikipedia.org/wiki/Gaussian_quadrature
      * https://en.wikipedia.org/wiki/Gauss–Kronrod_quadrature_formula
      * http://mathworld.wolfram.com/LobattoQuadrature.html <-- this would be better  
  */


template<class T> 
  void area_volume_romberg(
    T av[2],
    T xrange[2],
    const T & Omega0,
    const T & q,
    const T & F = 1,
    const T & delta = 1,
    const int & m = 1<<12,
    const bool & polish = false)
  {
   
    #if 1
    const int n = 10;

    // cos()^2 of Gauss-Lagrange nodes x_i in [0, Pi/2] 
    const T c[n]={0.999580064408504868392251519574,0.988810441236620370626267163745,0.93792975088501989705468221942,0.814698074016316269208833703789,0.615862835349225523334990208675,0.384137164650774476665009791325,0.185301925983683730791166296211,0.06207024911498010294531778058,0.011189558763379629373732836255,0.0004199355914951316077484804259};

    const T d[n]={0.0004197592455941272416985899940,0.011064352538060503513192603217,0.058217533289784414692374340457,0.150965122210421127009962838622,0.236575803384838256502140293303,0.236575803384838256502140293303,0.150965122210421127009962838622,0.058217533289784414692374340457,0.011064352538060503513192603217,0.0004197592455941272416985899940};
    
        
    // Gauss-Lagrange weights * Pi
    const T weights[n]={0.209454205485130325204687978,0.46951526056054717731928526,0.68828010698191943974759705,0.845926347240509469003431984,0.92841673332168682718764111324,0.92841673332168682718764111324,0.845926347240509469003431984,0.68828010698191943974759705,0.46951526056054717731928526,0.209454205485130325204687978};
   
    #else
    const int n = 15;
    
    // cos()^2 of Gauss-Lagrange nodes x_i in [0, Pi/2] 
    const T c[n]={0.999911065396171632736007253518,0.997574886997681325232904530099,0.985854212895255771578898079594,0.953879938591329487419660590834,0.890692147886632646280933330332,0.790164039217058712341619458721,0.655400162984313217687403793532,0.5,0.344599837015686782312596206468,0.209835960782941287658380541279,0.109307852113367353719066669668,0.046120061408670512580339409166,0.014145787104744228421101920406,0.0024251130023186747670954699007,0.0000889346038283672639927464819};
    
    const T d[n]={0.0000889266944646091553555375073,0.0024192318292446596704491832554,0.013945683811931480320682058875,0.04399300134433097347511477941,0.097359645579729565862303734566,0.165804830345241213606518615509,0.225850789344448888056399859737,0.25,0.225850789344448888056399859737,0.165804830345241213606518615509,0.097359645579729565862303734566,0.04399300134433097347511477941,0.013945683811931480320682058875,0.0024192318292446596704491832554,0.0000889266944646091553555375073};
        
    // Gauss-Lagrange weights * Pi
    const T weights[n]={0.0966141591290711189794453,0.22106145785079100848207,0.336650619784076362391549,0.4384742164293535129419684,0.5223501155128774575517226,0.584842030033819626900854873,0.6233908965445645786498187508,0.636418316610479145130430187624,0.6233908965445645786498187508,0.584842030033819626900854873,0.5223501155128774575517226,0.4384742164293535129419684,0.336650619784076362391549,0.22106145785079100848207,0.0966141591290711189794453};
    #endif
  
        
    //
    // Integration over the surface with RK4
    //  Def: 
    //    y = R^2, F = dR/dt, G = dR/dphi, 
    //    x = delta t, c = cos(phi)^2
    //  Eq:
    //    y_i(t)' = F(t, y_i(t), c_i) i = 0, .., n -1
    //    A'(t) = delta^2 sqrt(R + G^2/(4R) + F^2/4) df dt
    //    V'(t) = delta^3 1/2 R df dx
    //
    
    T d2 = delta*delta,
      d3 = d2*delta,
      a = d3*F*F,
      b = a*(1 + q),
      omega = delta*Omega0,
      
      t = xrange[0]/delta, 
      dt = (xrange[1] - xrange[0])/(m*delta),
      
      *V = new T [m + 1],
      *A = new T [m + 1],
      
      y1[n], y[n], k[4][n];

    
    // init point
    for (int i = 0; i < n; ++i) y[i] = 0;
           
    for (int j = 0; j < m; ++j){
      
      //
      // RK iteration
      //
      {
        
        T t1, f, g, f1, f2, s, s1, s2; // auxiliary variables
        
        // 1. step
        A[j] = V[j] = 0;
        
        s1 = t*t;
        s2 = (t-1)*(t-1);
        for (int i = 0; i < n; ++i) {
          s = y[i];
          
          // f1 = (R + t^2)^(-3/2)
          f1 = s + s1;
          f1 = 1/(f1*std::sqrt(f1));
          
          // f1 = (R + (t-1)^2)^(-3/2)
          f2 = s + s2;
          f2 = 1/(f2*std::sqrt(f2));
          
          // k0 = -dx F_t/F_R
          g = 1/(b*c[i] - f1 - q*f2);
          f = (q*(1 + (t - 1)*f2) + t*(f1 - b))*g; // = 1/2 dR/dt, x = delta*t
          k[0][i] = dt*2*f;
                
          g *= b; // = -(dR/dc)/R   Note: dR/dphi = -dR/dc*2*sqrt(c(1-c)) 
          A[j] += weights[i]*std::sqrt(s*(1 + g*g*d[i]) + f*f);
          V[j] += weights[i]*s;
        }
        
        if (j == m - 1) break; 
       
        // prepare: y1 = y + k0/2
        for (int i = 0; i < n; ++i) y1[i] = y[i] + k[0][i]/2;
        
        // 2. step
        t1 = t + dt/2;
        s1 = t1*t1;
        s2 = (t1-1)*(t1-1);
        for (int i = 0; i < n; ++i){
          s = y1[i];
          
          // f1 = (R + t1^2)^(-3/2)
          f1 = s + s1;
          f1 = 1/(f1*std::sqrt(f1));
          
          // f2 = (R + (t1-1)^2)^(-3/2)
          f2 = s + s2;
          f2 = 1/(f2*std::sqrt(f2));
          
          // k0 = -dx F_t/F_R
          k[1][i] = dt*2*(q*(1 + (t1 - 1)*f2) + t1*(f1 - b))/(b*c[i] - f1 - q*f2);
        }
        
        // prepare: y1 = y + k1/2
        for (int i = 0; i < n; ++i) y1[i] = y[i] + k[1][i]/2;
        
        // 3. step
        for (int i = 0; i < n; ++i) {
          s = y1[i];
          
          // f1 = (R + t1^2)^(-3/2)
          f1 = s + s1;
          f1 = 1/(f1*std::sqrt(f1));
          
          // f2 = (R + (t1-1)^2)^(-3/2)
          f2 = s + s2;
          f2 = 1/(f2*std::sqrt(f2));
          
          // k0 = -dx F_t/F_R
          k[2][i] = dt*2*(q*(1 + (t1 - 1)*f2) + t1*(f1 - b))/(b*c[i] - f1 - q*f2);
        }
        
      
        // y1 = y + k2
        for (int i = 0; i < n; ++i) y1[i] = y[i] + k[2][i];
        
        // 4. step
        t1 = t + dt;
        s1 = t1*t1;
        s2 = (t1-1)*(t1-1);
        for (int i = 0; i < n; ++i){
          s = y1[i];
          
          // f1 = (R + t1^2)^(-3/2)
          f1 = s + s1;
          f1 = 1/(f1*std::sqrt(f1));
          
          // f2 = (R + (t1-1)^2)^(-3/2)
          f2 = s + s2;
          f2 = 1/(f2*std::sqrt(f2));
          
          // k0 = -dx F_t/F_R
          k[3][i] = dt*2*(q*(1 + (t1 - 1)*f2) + t1*(f1 - b))/(b*c[i] - f1 - q*f2);
        }
              
        for (int i = 0; i < n; ++i)
          y[i] += (k[0][i] + 2*(k[1][i] + k[2][i]) + k[3][i])/6;  
      }
      
      t += dt;
      
      //
      // Polishing solutions with Newton-Rapson iteration
      //
      if (polish){
        
        const int it_max = 10;
        const T eps = 10*std::numeric_limits<T>::epsilon();
        const T min = 10*std::numeric_limits<T>::min();
        
        int it;
        
        T s1 = t*t, s2 = (t - 1)*(t - 1),
          s, f, df, ds, f1, f2, g1, g2;
          
        for (int i = 0; i < n; ++i) {
          
          s = y[i];

          it = 0;
          do {
            
            // g1 = (R + (t-1)^2)^(-1/2), f1 = (R + t^2)^(-3/2)
            f1 = 1/(s + s1);
            g1 = std::sqrt(f1);
            f1 *= g1;
            
            // g2 = (R + (t-1)^2)^(-1/2),  f2 = (R + (t-1)^2)^(-3/2)
            f2 = 1/(s + s2);
            g2 = std::sqrt(f2);
            f2 *= g2;
             
            df = -(b*c[i] - f1 - q*f2)/2;                     // = dF/dR
            f = omega - q*(g2 - t) - g1 - b*(c[i]*s + s1)/2;  // =F
            
            s -= (ds = f/df);
          
          } while (std::abs(f) > eps && 
                  std::abs(ds) > eps*std::abs(s) + min && 
                  ++it < it_max);
          
          if (!(it < it_max)) 
            std::cerr << "Volume: Polishing did not succeed\n";
        
          y[i] = s;
        }
      }
    }   
    
    A[m] = V[m] = 0;
        
    //
    // Romberg integration scheme -- repeated Richarson extrapolation
    //
    // Ref: 
    //  https://en.wikipedia.org/wiki/Romberg's_method
    //  https://en.wikipedia.org/wiki/Richardson_extrapolation
    //  https://ece.uwaterloo.ca/~dwharder/NumericalAnalysis/13Integration/romberg/complete.html
    
    // m = 2^K, determine K, initial step di in summation
    int K = 0, di = m;
    
    while ((m >> K) != 1) ++K;

    T *tA = new T[K+1], 
      *tV = new T[K+1],
      s1, s2, f1, f2;
    
    tA[0] = tV[0] = 0;
    
    for (int k = 1; k <= K; ++k) {
        
      //
      // compozite Trapezoidal rule
      //
            
      s1 = tA[0];     // temporary store integrals of previous division
      s2 = tV[0];  
      
      // calculate new integrals summing over only odd sites
      f1 = f2 = 0;
      for (int i = (di >> 1); i < m; i += di) {
        f1 += A[i]; 
        f2 += V[i];
      }
      
      // new integrals
      tA[0] = (s1 + di*f1)/2;
      tV[0] = (s2 + di*f2)/2;
      
      //
      // Richardson extrapolation
      //

      for (int i = 1, l = 4; i <= k; ++i, l <<= 2){
        
        f1 = (l*tA[i-1] - s1)/(l - 1);
        f2 = (l*tV[i-1] - s2)/(l - 1);
        
        // temporary store old integrals
        s1 = tA[i];
        s2 = tV[i];
        
        // store new values
        tA[i] = f1;
        tV[i] = f2;
      }

      // decrease step size
      di >>= 1;
    }
   
    av[0] = d2*dt*tA[K];
    av[1] = d3*dt*tV[K]/2;
   
    delete [] tA;
    delete [] tV;
    delete [] A;
    delete [] V;
  }


 
}


#endif // #if !defined(__volume_h)
