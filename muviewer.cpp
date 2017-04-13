#include "math.h"
#include "stdlib.h"
#include "stdio.h"
extern "C"
{
  void setv(double v1[3], double x, double y, double z)
  {
    v1[0] = x;
    v1[1] = y;
    v1[2] = z;
  }

  double* subv(double newv[3], double v1[3], double v2[3])
  {
    newv[0] = v1[0] - v2[0];
    newv[1] = v1[1] - v2[1];
    newv[2] = v1[2] - v2[2];
    return newv;
  }

  double dot(double v1[3], double v2[3])
  {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
  }

  double* cross(double newv[3], double u[3], double v[3])
  {
    newv[0] = (u[1]*v[2] - u[2]*v[1]);
    newv[1] = (u[2]*v[0] - u[0]*v[3]);
    newv[2] = (u[0]*v[1] - u[1]*v[0]);
    return newv;
  }

  double magnitude(double v1[3])
  {
    return sqrt(dot(v1,v1));
  }

  double angler(double u[3], double v[3])
  {
    double cResult[3];
    cross(cResult, u, v);
    return atan2(magnitude(cResult), dot(u, v));
  }

  double iso(double* regions, int size)
  {
    double mean;
    for(int i=0; i<size; i++)
      mean += regions[i] / size;
    double variance;
    for(int i=0; i<size; i++)
      variance += (regions[i] - mean)*(regions[i] - mean);
    return variance;
  }

  double* findBestPMT(double* charge, double* px, double* py, double* pz, 
      double cx, double cy, double cz, int npmt)
  {
    //printf("pmt count: %i\n", npmt);
    //for(int i=0; i<3; i++)
    //{
    //  printf("c++: charge %f, px %f, py %f, pz %f\n", charge[i], px[i], py[i], pz[i]);
    //}
    double pout[3], pin[3], pi[3], v1[3], v2[3];
    double* fOfq = (double*)malloc(npmt * sizeof(double));
    double costheta = 1/1.33;
    double rCut = 10000;
    setv(pout, cx, cy, cz);
    //printf("cx: %f, cy: %f, cz: %f, pout: (%f, %f, %f)\n",cx, cy, cz, pout[0], pout[1], pout[2]);
    //// test subv, dot, and magnitude
    //setv(pin, 100, 100, 100);
    //printf("mag(100,100,100): %f\n", magnitude(pin));
    //subv(v1, pout, pin);
    //printf("subv: (%f, %f, %f)", v1[0], v1[1], v1[2]);
    for(int testpmt=0; testpmt<npmt; testpmt++)
    {
      double qval = 0;
      //double isoval[10];// Split into 10 wedges (in theta)
      //for(int i=0; i<10; i++)
      //  isoval[i] = 0;
      setv(pin, px[testpmt], py[testpmt], pz[testpmt]);
      subv(v1,pout,pin);
      //subv(v1,pin,pout);
      double pVector[3];
      double sVector[3];
      int firstPVector = 1;
      if((magnitude(v1) > 0) && magnitude(pin) < 9000 )
      {
        for(int opmt=0; opmt<npmt; opmt++)
        {
          setv(pi, px[opmt], py[opmt], pz[opmt]);
          subv(v2, pi, pin);
          //subv(v2, pin, pi);
          if( (magnitude(v2) > 0) && magnitude(pi) < 9000 )
          {
            double comp = dot(v1,v2)/magnitude(v1)/magnitude(v2);
            //printf("comp: %f\n", comp);
            if( comp > costheta )
            {
              // If this is the first PMT to contribute, define it as pVector
              //if( firstPVector )
              //{
              //  firstPVector = 0;
              //  cross(pVector, v2, v1);
              //}
              //cross(sVector, v2, v1);
              //double angle = angler(pVector, sVector);
              //int wedge = int((angle + M_PI) / (2*M_PI) * 10);
              //isoval[wedge] += charge[opmt];
              if(charge[opmt] > 0)
                qval += 1;
              else
                qval -= 4;
              //qval += 1;
            }
          }
        }
      }
      //double isotropy = iso(isoval,10);
      //fOfq[testpmt] = isotropy;
      fOfq[testpmt] = qval;
      //printf("%f\n", qval);
    }
    //for(int i=0; i<10; i++)
    //  printf("%f\n", fOfq[i]);

    return fOfq;
  }

  double* findBestSpot(double* charge, double* px, double* py, double* pz,
      double cx, double cy, double cz, int npmt)
  {
    // This is similar to findBestPMT but instead of returning a pmt, it returns
    // the exact position on a sphere.
    // Detector is in mm, with 8000mm radius
    double costheta = 1/1.33;
    double precision = 10; // mm
    double start_precision = precision * pow(2,4);
    double lprec = start_precision;
    double radius = 8500; // mm
    double pi = M_PI;
    double ct_min = -1;
    double ct_max = 1;
    double phi_min = -pi/2;
    double phi_max = pi/2;
    double pout[3], pin[3], ptest[3], v1[3], v2[3];
    double best_ct = ct_min;
    double best_phi = phi_min;
    setv(pout, cx, cy, cz);
    while( lprec >= precision )
    {
      // Calculate the enclosed charge from point on grid
      double precPhi = (phi_max - phi_min) / 10;
      double precCt = (ct_max - ct_min) / 10;
      best_ct = ct_min;
      best_phi = phi_min;
      double best_qval = 0;
      // Find the best point onthe coordinate map
      for(double ct = ct_min; ct < ct_max; ct += precCt)
        for(double phi = phi_min; phi < phi_max; phi += precPhi)
        {
          // Sum charge at this location
          double xpos = radius * sin( acos(ct) ) * cos(phi);
          double ypos = radius * sin( acos(ct) ) * sin(phi);
          double zpos = radius * ct;
          setv(pin, xpos, ypos, zpos);
          subv(v1, pout, pin);
          double qval = 0;
          for(int opmt=0; opmt<npmt; opmt++)
          {
            setv(ptest, px[opmt], py[opmt], pz[opmt]);
            subv(v2, ptest, pin);
            if( magnitude(v2) > 0 && magnitude(ptest) < 9000 )
            {
              double comp = dot(v1,v2)/magnitude(v1)/magnitude(v2);
              if( comp > costheta )
                //qval += charge[opmt];
                if( charge[opmt] > 0 )
                  qval += 1;
                else
                  qval -= 4;
            }
          }
          if( qval > best_qval )
          {
            best_qval = qval;
            best_ct = ct;
            best_phi = phi;
          }
        }
      // Update the map to be smaller and more precise
      ct_min = best_ct   - 4*precCt;
      ct_max = best_ct   + 4*precCt;
      phi_min = best_phi - 4*precPhi;
      phi_max = best_phi + 4*precPhi;
      lprec = lprec / 2;
    }
    printf("best phi: %f, best ctheta: %f\n",best_phi,best_ct);
    double* returnVector = (double*)malloc(3*sizeof(double));
    returnVector[0] = radius * sin( acos(best_ct) ) * cos(best_phi);
    returnVector[1] = radius * sin( acos(best_ct) ) * sin(best_phi);
    returnVector[2] = radius * best_ct;

    return returnVector;
  }

  double* findNormBestSpot(double* charge, double* px, double* py, double* pz,
      double cx, double cy, double cz, int npmt)
  {
    // This is similar to findBestPMT but instead of returning a pmt, it returns
    // the exact position on a sphere.
    // Detector is in mm, with 8000mm radius
    double costheta = 1/1.33;
    double precision = 10; // mm
    double start_precision = precision * pow(2,4);
    double lprec = start_precision;
    double radius = 8500; // mm
    double pi = M_PI;
    double ct_min = -1;
    double ct_max = 1;
    double phi_min = -pi/2;
    double phi_max = pi/2;
    double pout[3], pin[3], ptest[3], v1[3], v2[3];
    double best_ct = ct_min;
    double best_phi = phi_min;
    setv(pout, cx, cy, cz);
    while( lprec >= precision )
    {
      // Calculate the enclosed charge from point on grid
      double precPhi = (phi_max - phi_min) / 10;
      double precCt = (ct_max - ct_min) / 10;
      best_ct = ct_min;
      best_phi = phi_min;
      double best_qval = 0;
      // Find the best point onthe coordinate map
      for(double ct = ct_min; ct < ct_max; ct += precCt)
        for(double phi = phi_min; phi < phi_max; phi += precPhi)
        {
          // Sum charge at this location
          double xpos = radius * sin( acos(ct) ) * cos(phi);
          double ypos = radius * sin( acos(ct) ) * sin(phi);
          double zpos = radius * ct;
          setv(pin, xpos, ypos, zpos);
          subv(v1, pout, pin);
          double qval = 0;
          double area = 0;
          for(int opmt=0; opmt<npmt; opmt++)
          {
            setv(ptest, px[opmt], py[opmt], pz[opmt]);
            subv(v2, ptest, pin);
            if( magnitude(v2) > 0 )
            {
              double comp = dot(v1,v2)/magnitude(v1)/magnitude(v2);
              if( comp > costheta )
              {
                area += 1;
                qval += charge[opmt];
              }
            }
          }
          if( area > 0 )
            qval /= area;
          if( qval > best_qval )
          {
            best_qval = qval;
            best_ct = ct;
            best_phi = phi;
          }
        }
      // Update the map to be smaller and more precise
      ct_min = best_ct   - 2*precCt;
      ct_max = best_ct   + 2*precCt;
      phi_min = best_phi - 2*precPhi;
      phi_max = best_phi + 2*precPhi;
      lprec = lprec / 2;
    }
    printf("best phi: %f, best ctheta: %f\n",best_phi,best_ct);
    double* returnVector = (double*)malloc(3*sizeof(double));
    returnVector[0] = radius * sin( acos(best_ct) ) * cos(best_phi);
    returnVector[1] = radius * sin( acos(best_ct) ) * sin(best_phi);
    returnVector[2] = radius * best_ct;

    return returnVector;
  }
}
