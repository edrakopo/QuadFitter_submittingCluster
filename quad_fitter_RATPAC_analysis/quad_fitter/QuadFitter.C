///////////  Reads a RAT-PAC output file and does a quick analysis ////////

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLeaf.h>
#include <Rtypes.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TH2.h>
#include <TH3.h>
#include <TPad.h>
#include <TVector3.h>
#include <TString.h>
#include <TPRegexp.h>
#include <TGraph.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TClonesArray.h>

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>

// #include "rootstart.h"

// Header file for the classes stored in the TTree if any.
#include <RAT/DS/Root.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/Calib.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/PMTInfo.hh>
#include <RAT/DS/RunStore.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DSReader.hh>
#include <RAT/TrackNav.hh>
#include <RAT/TrackCursor.hh>
#include <RAT/TrackNode.hh>
#include <RAT/DB.hh>

using namespace std;
using namespace TMath;

#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TH2.h"
#include "TGraph.h"
#include <TCanvas.h>

double  R_SPHERE=2000.; //sphere diameter [cm]
double  R_LENGTH=2000.; //length of detector [cm]
double N_REF=1.34; //average index of refraction
double C_VAC=29.9792458; //speed of light in vacuum [cm/ns]
int NSeedsTarget=50; //number of quadruplets
double TSIGMA=0.5; //total time spread (including detector TTS chromatic dispersions)
int fNDigits=0;
int fThisDigit=0;
int fLastEntry=0;
int fLastEntry_once=0;
int fCounter=0;
int fMinTime=0;
double fBaseFOM=100.0; //Figure of merit. Borrowed from WCSim: the higher it is the better
double meanTime=0.;
double seedTime=0.;

std::vector<double> fDigitX;
std::vector<double> fDigitY;
std::vector<double> fDigitZ;
std::vector<double> fDigitT;
std::vector<double> fDigitQ;
std::vector<double> fDigitPE;
std::vector<double> fDigitW;
std::vector<double> fDigitV;
std::vector<double> fDelta;

std::map<int, double> INDEX;
TRandom RND;
// store seed vertex calculated from quaruplets
std::vector<vector<double>> vSeed;
std::vector<vector<double>> vSeeds;
std::vector<double> vSeedVtxX;
std::vector<double> vSeedVtxY;
std::vector<double> vSeedVtxZ;
std::vector<double> vSeedVtxTime;
std::vector<int> vSeedDigitList;
std::vector<double> fLastEntry_v;

// This function solves system of 4 equations with 4 unknowns to find the seed vertex for each quadruple
int FindVertex(Double_t x0, Double_t y0, Double_t z0, Double_t t0, Double_t x1, Double_t y1, Double_t z1, Double_t t1, Double_t x2, Double_t y2, Double_t z2, Double_t t2, Double_t x3, Double_t y3, Double_t z3, Double_t t3, Double_t& vxm, Double_t& vym, Double_t& vzm, Double_t& vtm, Double_t& vxp, Double_t& vyp, Double_t& vzp, Double_t& vtp)
{
    vxm = -99999.9;
    vym = -99999.9;
    vzm = -99999.9;
    vtm = -99999.9;

    vxp = -99999.9;
    vyp = -99999.9;
    vzp = -99999.9;
    vtp = -99999.9;

    // speed of light in water
    // =======================
    Double_t c = C_VAC/N_REF;

    // causality checks
    // ================
    if( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0) >= c*c*(t1-t0)*(t1-t0)
        && (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) >= c*c*(t2-t1)*(t2-t1)
        && (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) + (z3-z2)*(z3-z2) >= c*c*(t3-t2)*(t3-t2)
        && (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0) >= c*c*(t2-t0)*(t2-t0)
        && (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) + (z3-z1)*(z3-z1) >= c*c*(t3-t1)*(t3-t1)
        && (x3-x0)*(x3-x0) + (y3-y0)*(y3-y0) + (z3-z0)*(z3-z0) >= c*c*(t3-t0)*(t3-t0) ){

        // [Note: for causality, require that |x_{i}-x_{j}| >= c*|t_{i}-t_{j}|
        //        for each pair of points]
        Double_t dx1 = x1-x0;  Double_t dy1 = y1-y0;  Double_t dz1 = z1-z0;  Double_t dt1 = c*(t1-t0);
        Double_t dx2 = x2-x0;  Double_t dy2 = y2-y0;  Double_t dz2 = z2-z0;  Double_t dt2 = c*(t2-t0);
        Double_t dx3 = x3-x0;  Double_t dy3 = y3-y0;  Double_t dz3 = z3-z0;  Double_t dt3 = c*(t3-t0);

        Double_t epsilon = 1.0e-7;

        // check that points don't all lie in a plane
        if( !( fabs(dx1)<epsilon && fabs(dx2)<epsilon && fabs(dx3)<epsilon )
            && !( fabs(dy1)<epsilon && fabs(dy2)<epsilon && fabs(dy3)<epsilon )
            && !( fabs(dz1)<epsilon && fabs(dz2)<epsilon && fabs(dz3)<epsilon )
            && !( fabs(dx1)<epsilon && fabs(dy1)<epsilon && fabs(dz1)<epsilon )
            && !( fabs(dx2)<epsilon && fabs(dy2)<epsilon && fabs(dz2)<epsilon )
            && !( fabs(dx3)<epsilon && fabs(dy3)<epsilon && fabs(dz3)<epsilon ) ){

            // [Note: this is a problem for detectors with flat faces!]

            Double_t Mdata[9] = { dx1+epsilon, dy1, dz1,
                                  dx2, dy2+epsilon, dz2,
                                  dx3, dy3, dz3+epsilon };

            Double_t Qdata[3] = { 0.5*( dx1*dx1 + dy1*dy1 + dz1*dz1 - dt1*dt1 ),
                                  0.5*( dx2*dx2 + dy2*dy2 + dz2*dz2 - dt2*dt2 ),
                                  0.5*( dx3*dx3 + dy3*dy3 + dz3*dz3 - dt3*dt3 ) };

            Double_t Tdata[3] = { dt1,
                                  dt2,
                                  dt3 };

            TMatrixD M(3,3,Mdata);
            TMatrixD Q(3,1,Qdata);
            TMatrixD T(3,1,Tdata);

            if( M.Determinant() != 0.0 ){

                TMatrixD A(3,1);
                TMatrixD B(3,1);

                M.Invert();
                A.Mult(M,T);
                B.Mult(M,Q);

                Double_t ax = A(0,0);
                Double_t ay = A(1,0);
                Double_t az = A(2,0);

                Double_t bx = B(0,0);
                Double_t by = B(1,0);
                Double_t bz = B(2,0);

                Double_t ab = ax*bx + ay*by + az*bz;
                Double_t a2 = ax*ax + ay*ay + az*az;
                Double_t b2 = bx*bx + by*by + bz*bz;

                Double_t qa = a2-1.0;
                Double_t qb = 2.0*ab;
                Double_t qc = b2;

                // check for solutions
                if( qb*qb-4.0*qa*qc>0.0 ){

                    // The common vertex is given by a quadratic equation, which has two solutions.
                    // Typically, one solution corresponds to photons travelling forwards in time,
                    // and the other solution corresponds to photons travelling backwards in time.
                    // However, sometimes there appear to be two valid solutions.

                    Double_t ctm = ( -qb - sqrt(qb*qb-4.0*qa*qc) ) / ( 2.0*qa );
                    Double_t ctp = ( -qb + sqrt(qb*qb-4.0*qa*qc) ) / ( 2.0*qa );

                    Double_t tm = t0 + ctm/c;
                    Double_t xm = x0 + ctm*ax + bx;
                    Double_t ym = y0 + ctm*ay + by;
                    Double_t zm = z0 + ctm*az + bz;

                    if( tm<t0 && tm<t1
                        && tm<t2 && tm<t3 ){
                        vxm = xm;
                        vym = ym;
                        vzm = zm;
                        vtm = tm;
                    }

                    Double_t tp = t0 + ctp/c;
                    Double_t xp = x0 + ctp*ax + bx;
                    Double_t yp = y0 + ctp*ay + by;
                    Double_t zp = z0 + ctp*az + bz;

                    if( tp<t0 && tp<t1
                        && tp<t2 && tp<t3 ){
                        vxp = xp;
                        vyp = yp;
                        vzp = zp;
                        vtp = tp;
                    }

                }
                else
                {
                   //std::cout << "qb*qb-4.0*qa*qc<0.0." << std::endl;
                }
            }
            else
            {
                //std::cout << "M.Determinant() == 0.0. " << std::endl;
            }

        }
        else
        {
           // std::cout << "The 4 points lie on a plane. " << std::endl;
        }
    }
    else
    {
        // std::cout << "Causality check for this 4-group of digits has not been passed. " << std::endl;
    }


    return 0;
}

Int_t ChooseNextDigitOnce(Double_t& xpos, Double_t& ypos, Double_t& zpos, Double_t& time, std::vector<double> fDigitX,std::vector<double> fDigitY,std::vector<double> fDigitZ,std::vector<double> fDigitT, std::vector<int> vSeedDigitList){
	
    xpos=0.; ypos=0.; zpos=0.; time=0.;
    // ROOT random number generator
    Double_t r = RND.Rndm();
    // pseudo-random number generator
    Int_t numEntries = (Int_t) vSeedDigitList.size();

    fLastEntry = (Int_t)(r*numEntries);
    //check if the element already exists in fLastEntry vector
    if( (std::count(fLastEntry_v.begin(), fLastEntry_v.end(), fLastEntry)) ==0){      
      fLastEntry_v.push_back(fLastEntry);
      //cout<<"fLastEntry: "<<fLastEntry<<endl;
      fLastEntry_once=fLastEntry;
    }else{ //cout<<"entry already exists!!! "<<fLastEntry<<endl; 
       return -1;
     }
    
    xpos = fDigitX[fLastEntry_once];
    ypos = fDigitY[fLastEntry_once];
    zpos = fDigitZ[fLastEntry_once];
    time = fDigitT[fLastEntry_once];
    //cout<<"return-- fLastEntry_once: "<<fLastEntry_once<<endl;
    //std::cout<<"xpos = "<<xpos<<" ypos = "<<ypos<<" zpos = "<<zpos<<"   time = "<<time<<std::endl;
    //return fLastEntry;          
    return fLastEntry_once;
}

Int_t ChooseNextDigit(Double_t& xpos, Double_t& ypos, Double_t& zpos, Double_t& time, std::vector<double> fDigitX,std::vector<double> fDigitY,std::vector<double> fDigitZ,std::vector<double> fDigitT, std::vector<int> vSeedDigitList){

    xpos=0.; ypos=0.; zpos=0.; time=0.;
    // ROOT random number generator
    Double_t r = RND.Rndm();
    //cout<<"random number is: "<<r<<endl;

    //  std::cout<<"I'm inside ChooseNextDitig"<<std::endl
    // pseudo-random number generato
    Int_t numEntries = (Int_t) vSeedDigitList.size();
    //cout<<"In ChooseNextDigit numEntries: "<<numEntries<<endl;
/*  cout<<"fCounter = "<<fCounter<<"   fLastEntry = "<<fLastEntry<<endl;
  fCounter++
  if( fCounter>=fNDigits ) fCounter = 0;
  fThisDigit = vSeedDigitList.at(fLastEntry);

  Double_t t0 = 0.5 + fDigitT[fCounter] - fMinTime;
  Double_t q0 = 0.5 + fDigitQ[fCounter];

  Double_t t1 = 0.5 + fDigitT[fThisDigit] - fMinTime;
  Double_t q1 = 0.5 + fDigitQ[fThisDigit];

  Double_t tq = 100.0*(t0*q0+t1*q1);
  Double_t r = tq - TMath::Floor(tq);
*/
//  r = gRandom->Uniform(); // Christoph Aberle, August 14, 2013: use of a proper RN generator since I saw that quadruplets were duplicated with the pseudo-random number generator used in the lines above

    fLastEntry = (Int_t)(r*numEntries);

    //std::cout<<"In ChooseNextDigit- fLastEntry = "<<fLastEntry<<"   r = "<<r<<"   numEntries = "<<numEntries<<std::endl;

    // return the new digit
    //fThisDigit = vSeedDigitList.at((unsigned long) fLastEntry);
    //fThisDigit = vSeedDigitList[fLastEntry];
    //  std::cout<<"fThisDigit = "<<fThisDigit<<std::endl;
    /*xpos = fDigitX[fThisDigit];
    ypos = fDigitY[fThisDigit];
    zpos = fDigitZ[fThisDigit];
    time = fDigitT[fThisDigit];
    std::cout<<"xpos = "<<xpos<<"   ypos = "<<ypos<<"   zpos = "<<zpos<<"   time = "<<time<<std::endl;
    return fThisDigit;*/
    xpos = fDigitX[fLastEntry];
    ypos = fDigitY[fLastEntry];
    zpos = fDigitZ[fLastEntry];
    time = fDigitT[fLastEntry];
    //std::cout<<"xpos = "<<xpos<<"   ypos = "<<ypos<<"   zpos = "<<zpos<<"   time = "<<time<<std::endl;
    return fLastEntry;
}


// Make a quadruple
Int_t ChooseNextQuadruple(Double_t& x0, Double_t& y0, Double_t& z0, Double_t& t0, Double_t& x1, Double_t& y1, Double_t& z1, Double_t& t1, Double_t& x2, Double_t& y2, Double_t& z2, Double_t& t2, Double_t& x3, Double_t& y3, Double_t& z3, Double_t& t3, std::vector<double> fDigitX,std::vector<double> fDigitY,std::vector<double> fDigitZ,std::vector<double> fDigitT,std::vector<int> vSeedDigitList)
{
    int code=0; // 0 -if OK, 1 -if failed to chose 4 different digits
    int digit0=0;
    int digit1=-1;
    int digit2=-2;
    int digit3=-3;
    int counter1=0;
    int counter2=0;
    int counter3=0;
    digit0 = ChooseNextDigit(x0,y0,z0,t0,fDigitX,fDigitY,fDigitZ,fDigitT,vSeedDigitList);
    //digit0 = ChooseNextDigitOnce(x0,y0,z0,t0,fDigitX,fDigitY,fDigitZ,fDigitT,vSeedDigitList);

///  digit1 = ChooseNextDigit(x1,y1,z1,t1);
///  digit2 = ChooseNextDigit(x2,y2,z2,t2);
///  digit3 = ChooseNextDigit(x3,y3,z3,t3);

//cout<<"I'm inside ChooseNextQuadruple: "<<vSeedDigitList.size()<<endl;
//check that selected digits for the quadruple are not identical
//if they are after 100 attempts then be it, the quadruple will not
//be used later
    while( (digit1<0 || digit1==digit0) && counter1<100)
    {
        digit1 = ChooseNextDigit(x1,y1,z1,t1,fDigitX,fDigitY,fDigitZ,fDigitT,vSeedDigitList);
        //digit1 = ChooseNextDigitOnce(x1,y1,z1,t1,fDigitX,fDigitY,fDigitZ,fDigitT,vSeedDigitList);
        counter1++;
    }

    while( (digit2<0 || digit2==digit0 || digit2==digit1) && counter2<100)
    {
        digit2 = ChooseNextDigit(x2,y2,z2,t2,fDigitX,fDigitY,fDigitZ,fDigitT,vSeedDigitList);
        //digit2 = ChooseNextDigitOnce(x2,y2,z2,t2,fDigitX,fDigitY,fDigitZ,fDigitT,vSeedDigitList);
        counter2++;
    }

    while( (digit3<0 || digit3==digit0 || digit3==digit1 || digit3==digit2) && counter3<100)
    {
        digit3 = ChooseNextDigit(x3,y3,z3,t3,fDigitX,fDigitY,fDigitZ,fDigitT,vSeedDigitList);
        //digit3 = ChooseNextDigitOnce(x3,y3,z3,t3,fDigitX,fDigitY,fDigitZ,fDigitT,vSeedDigitList);
        counter3++;
    }

    if(counter1>=100 || counter2>=100 || counter3>=100) code=1;

    return code;
}

std::vector<double> GetSeedVtx(vector<vector<double>> vSeeds, int wvert){

    if(wvert==0) return vSeeds[0]; //vSeedVtxX;
    if(wvert==1) return vSeeds[1]; //vSeedVtxY;
    if(wvert==2) return vSeeds[2]; //vSeedVtxZ;
    else return vSeeds[3]; //vSeedVtxTime;
}

// Calculate NSeedsTarget vertices
//Int_t vector<vector<int> > vec;
std::vector<vector<double>> CalcVertexSeeds(std::vector<double> fDigitX,std::vector<double> fDigitY,std::vector<double> fDigitZ,std::vector<double> fDigitT,std::vector<double> fDigitQ,std::vector<double> fDigitW,std::vector<double> fDigitV,std::vector<int> vSeedDigitList)
{
    // reset list of seeds
    // ===================
    vSeed.clear();
    vSeedVtxX.clear();
    vSeedVtxY.clear();
    vSeedVtxZ.clear();
    vSeedVtxTime.clear();

    Double_t x0 = 0.0;
    Double_t y0 = 0.0;
    Double_t z0 = 0.0;
    Double_t t0 = 0.0;

    Double_t x1 = 0.0;
    Double_t y1 = 0.0;
    Double_t z1 = 0.0;
    Double_t t1 = 0.0;

    Double_t x2 = 0.0;
    Double_t y2 = 0.0;
    Double_t z2 = 0.0;
    Double_t t2 = 0.0;

    Double_t x3 = 0.0;
    Double_t y3 = 0.0;
    Double_t z3 = 0.0;
    Double_t t3 = 0.0;

    double fVtxX1, fVtxY1, fVtxZ1, fVtxTime1;
    double fVtxX2, fVtxY2, fVtxZ2, fVtxTime2;

    int counter=0; int counter0=0;
    //std::cout<<"I'm inside CalcVertexSeeds"<<std::endl;
    //std::cout<<"vSeedVtxX.size(): "<<vSeedVtxX.size()<<" NSeedsTarget "<<NSeedsTarget<<" counter "<<counter<<" 100*NSeedsTarget "<<100*NSeedsTarget<<std::endl;
    while( vSeedVtxX.size()<NSeedsTarget && counter<100*NSeedsTarget && counter0<400*NSeedsTarget)
    {
        //    std::cout<<"counter = "<<std::endl;
        ChooseNextQuadruple(x0,y0,z0,t0,
                            x1,y1,z1,t1,
                            x2,y2,z2,t2,
                            x3,y3,z3,t3,
			    fDigitX,fDigitY,fDigitZ,fDigitT,vSeedDigitList);
        
        /*std::cout<<"counter = "<<counter<<std::endl;
        std::cout << "   digit0: (x,y,z,t)=(" << x0 << "," << y0 << "," << z0 << "," << t0 << ") " << std::endl;
        std::cout << "   digit1: (x,y,z,t)=(" << x1 << "," << y1 << "," << z1 << "," << t1 << ") " << std::endl;
        std::cout << "   digit2: (x,y,z,t)=(" << x2 << "," << y2 << "," << z2 << "," << t2 << ") " << std::endl;
        std::cout << "   digit3: (x,y,z,t)=(" << x3 << "," << y3 << "," << z3 << "," << t3 << ") " << std::endl;
        */

        FindVertex(x0,y0,z0,t0,
                   x1,y1,z1,t1,
                   x2,y2,z2,t2,
                   x3,y3,z3,t3,
                   fVtxX1,fVtxY1,fVtxZ1,fVtxTime1,
                   fVtxX2,fVtxY2,fVtxZ2,fVtxTime2);
        
         /*std::cout << "   result1: (x,y,z,t)=(" << fVtxX1 << "," << fVtxY1 << "," << fVtxZ1 << "," << fVtxTime1 << ") " << std::endl
                   << "   result2: (x,y,z,t)=(" << fVtxX2 << "," << fVtxY2 << "," << fVtxZ2 << "," << fVtxTime2 << ") " << std::endl;
         std::cout<<"-------------"<<std::endl;
         */

        if(fVtxX1==-99999.9 && fVtxX2==-99999.9) //no solutions try anouther quadruple
        {
            //counter++;
            counter0++;
            continue;
        }


        bool inside_det=0;
        inside_det= (sqrt(fVtxX1*fVtxX1+fVtxY1*fVtxY1)<R_SPHERE && abs(fVtxZ1)<R_LENGTH);
        //std::cout<<"Solution1: inside_det = "<<inside_det<<"  X = "<<fVtxX1<<"  Y = "<<fVtxY1<<"  Z = "<<fVtxZ1<<"  T = "<<fVtxTime1<<std::endl;
//    if(!inside_det)
//    {
//      std::cout<<"Solution outside detector"<<std::endl;
//      continue;
//    }

        // add first digit
        if( inside_det ){
            vSeedVtxX.push_back(fVtxX1);
            vSeedVtxY.push_back(fVtxY1);
            vSeedVtxZ.push_back(fVtxZ1);
            vSeedVtxTime.push_back(fVtxTime1);
            //std::cout << "New vertex seed 1: x= " << fVtxX1 << ", " << fVtxY1 << ", " << fVtxZ1 << ", " <<fVtxTime1 << std::endl;
            //std::cout<<"Solution1: inside_det = "<<inside_det<<"  X = "<<fVtxX1<<"  Y = "<<fVtxY1<<"  Z = "<<fVtxZ1<<"  T = "<<fVtxTime1<<std::endl;
            //cout<<"fVtxZ1<abs(R_LENGTH) = "<<fVtxZ1<<"<"<<abs(R_LENGTH)<<endl;
            //cout<<"R= "<<sqrt(fVtxX1*fVtxX1+fVtxY1*fVtxY1)<<endl;
        }


        inside_det= (sqrt(fVtxX2*fVtxX2+fVtxY2*fVtxY2)<R_SPHERE && abs(fVtxZ2)<R_LENGTH);
        //std::cout<<"Solution2: inside_det = "<<inside_det<<"  X = "<<fVtxX2<<"  Y = "<<fVtxY2<<"  Z = "<<fVtxZ2<<"  T = "<<fVtxTime2<<std::endl;
//    if(!inside_det)
//    {
//      std::cout<<"Solution outside detector"<<std::endl;
//      continue;
//    }

        // add second digit
        if( inside_det ){
            vSeedVtxX.push_back(fVtxX2);
            vSeedVtxY.push_back(fVtxY2);
            vSeedVtxZ.push_back(fVtxZ2);
            vSeedVtxTime.push_back(fVtxTime2);
            //std::cout << "New vertex seed 2: x= " << fVtxX2 << ", " << fVtxY2 << ", " << fVtxZ2 << ", " <<fVtxTime2 << std::endl;
            //std::cout<<"Solution2: inside_det = "<<inside_det<<"  X = "<<fVtxX2<<"  Y = "<<fVtxY2<<"  Z = "<<fVtxZ2<<"  T = "<<fVtxTime2<<std::endl;
            //cout<<"fVtxZ2<abs(R_LENGTH) = "<<fVtxZ2<<"<"<<abs(R_LENGTH)<<endl;
            //cout<<"R= "<<sqrt(fVtxX2*fVtxX2+fVtxY2*fVtxY2)<<endl;
        }
        counter++;

    }

     //std::cout << "The number of calculated seeds in vSeedVtxX, vSeedVtxY, vSeedVtxZ, vSeedVtxTime for this event is = " << vSeedVtxX.size() << std::endl;
     vSeed.push_back(vSeedVtxX);
     vSeed.push_back(vSeedVtxY);
     vSeed.push_back(vSeedVtxZ);
     vSeed.push_back(vSeedVtxTime);
     //cout<<"vSeed[0].size(): "<<vSeed[0].size()<<" vSeed[3].size(): "<<vSeed[3].size()<<endl;

    //return 0;
    return vSeed;
}

double fomtotal;
double weightedRecoX;
double weightedRecoY;
double weightedRecoZ;
std::vector<double> fomarray;

// Calculates FOM by looking at time residuals
// make sure fDelta has been filled before calling this function
Int_t TimePropertiesLnL(double & vtx_time, double & fom, std::vector<double> fDigitX, std::vector<double> fDelta)
{
    double A = 1.0 / ( 2.0*TSIGMA*TMath::Sqrt(0.5*TMath::Pi()) );
    double Preal;
    fom=0.;
    double chi2=0.;
    double ndof=0.0;

    //std::cout<<"fDelta.size() = "<<fDelta.size()<<std::endl;
    for( Int_t idigit=0; idigit<fDigitX.size(); idigit++ )
    {
      double delta = fDelta[idigit] - vtx_time;
      //std::cout << delta << std::endl;
      Preal = A*exp(-(delta*delta)/(2.0*TSIGMA*TSIGMA));
      // P = (1.0-Pnoise)*Preal + Pnoise;
      // chi2 += -2.0*log(P);
      //std::cout<<"A = "<<A<<"   delta = "<<delta<<"   TSIGMA = "<<TSIGMA<<"   Preal = "<<Preal<<std::endl;
      //std::cout<<"delta = "<<delta<<"   Preal = "<<Preal<<"   fDelta[idigit] = "<<fDelta[idigit]<<"   vtx_time = "<<vtx_time<<std::endl;
      // chi2 += -2.0*log(Preal);
      // chi2 += 2.0*(-Preal);
      
      chi2 += Preal;
        ndof += 1.0;
    }
     //std::cout<<"chi2 = "<<chi2<<"   ndof = "<<ndof<<std::endl;
    // if( ndof>0.0 ){
    // fom = fBaseFOM - 1.0*chi2/ndof;
    // fom = fBaseFOM - 100.0*chi2/ndof;
    fom = chi2;
    fomtotal +=fom;
     //std::cout << "chi2/fom=" << fom << std::endl;
    // fom = chi2/ndof;
    // }

       return 0;
}

Int_t SelectBestSeed(std::vector<double> fDigitX,std::vector<double> fDigitY,std::vector<double> fDigitZ,std::vector<double> fDigitT,std::vector<double> vSeedVtxX,std::vector<double> vSeedVtxY,std::vector<double> vSeedVtxZ,std::vector<double> vSeedVtxTime)
{
    Int_t bestSeed = -1;
    Double_t bestFOM = -1.0;
    double dRmin=200.;
    int seedRmin =0;

    fomtotal = 0.;
    weightedRecoX = 0.;
    weightedRecoY = 0.;
    weightedRecoZ = 0.;
    fomarray.clear();

    for(int i=0;i!=vSeedVtxX.size();++i)
    {
     // loop over digits
     double Swx=0.;
     double Sw=0.;     
     fDelta.clear();
        for( Int_t idigit=0; idigit<fDigitX.size(); idigit++ )
        {
            Double_t dx = fDigitX[idigit]-vSeedVtxX[i];
            Double_t dy = fDigitY[idigit]-vSeedVtxY[i];
            Double_t dz = fDigitZ[idigit]-vSeedVtxZ[i];
            Double_t ds = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
            //need to check what is the proper variable out of the two listed below
            Double_t time0 = fDigitT[idigit] - 0; //this is what was done in WCSim for the JINST paper
            Double_t time1 = fDigitT[idigit] - vSeedVtxTime[i];
      
            double dRseed=sqrt((vSeedVtxX[i]*vSeedVtxX[i])+(vSeedVtxY[i]*vSeedVtxY[i])+(vSeedVtxZ[i]*vSeedVtxZ[i]));//valid only if mcx,y,z=0.
            if(dRseed<dRmin){
               dRmin=dRseed;
               seedRmin=i;
            } 
         /////   
            double fPointResidual0 = time0 - ds/(C_VAC/N_REF);
            fDelta.push_back(fPointResidual0);
            double weight = 1.0/(TSIGMA*TSIGMA);
            Swx += time1*weight; // here is some room for upgrade id TSIGMA is not always the same
            Sw += weight;
        ////comment till here
        }
       //////
        meanTime=Swx/Sw;
        seedTime=vSeedVtxTime[i];

        double vtx_time=0.;
        double fom=0.;
        TimePropertiesLnL(vtx_time,fom,fDigitX,fDelta); 
        //cout<<"fom "<<fom<<endl;
        fomarray.push_back(fom);
         if( fom>bestFOM ){
            bestSeed = i;
            bestFOM = fom;
        }/////comment till here
     }
     //std::cout<<"fomtotal: " << fomtotal << std::endl;
     /////////////
     for(int iweight=0;iweight!=vSeedVtxX.size();++iweight){
        double seedweight= fomarray[iweight]/fomtotal;
        weightedRecoX += vSeedVtxX[iweight]*seedweight;
        weightedRecoY += vSeedVtxY[iweight]*seedweight;
        weightedRecoZ += vSeedVtxZ[iweight]*seedweight;
        // std::cout <<  seedweight << "  " << weightedRecoX << "  " << vSeedVtxX[iweight] << std::endl;
     }////comment till here
     //std::cout << "--------------------:" << weightedRecoX << "  " << weightedRecoY << "  " << weightedRecoZ << std::endl;
     cout<<"minimum dR is: "<<dRmin<<endl;

   return bestSeed;
   //return seedRmin; //returning the seed resulting in the minimum dR
}
//---------------------------------------------------------
//analyse the output file from ratpac 
void QuadFitter(const char* filename_ratpac) {
     
      //create the output .csv file
      ofstream csvfile;
      csvfile.open("seeds.csv");
      
      TH1D* plot = new TH1D("plot", "plot", 300, 0., 600.);
      TH1D* plot_VSeedVtxX= new TH1D("plot_VSeedVtxX", "plot_VSeedVtxX", 1000, -1000., 1000.); //300, 0., 600.);
      TH1D* plot_VSeedVtxY= new TH1D("plot_VSeedVtxY", "plot_VSeedVtxY", 1000, -1000., 1000.); //300, 0., 600.);
      TH1D* plot_VSeedVtxZ= new TH1D("plot_VSeedVtxZ", "plot_VSeedVtxZ", 1000, -1000., 1000.); //300, 0., 600.);
      TH1D* plot_VSeedVtxdr= new TH1D("plot_VSeedVtxdr", "plot_VSeedVtxdr", 1000, 0., 2000.);
      TH1D* plot_VSeedVtxdR= new TH1D("plot_VSeedVtxdR", "plot_VSeedVtxdR", 1000, 0., 2000.);
      TH1D* plot_fdigitX= new TH1D("plot_fdigitX", "plot_fdigitX", 1600, -800., 800.);
      TH1D* plot_fdigitY= new TH1D("plot_fdigitY", "plot_fdigitY", 1600, -800., 800.);
      TH1D* plot_fdigitZ= new TH1D("plot_fdigitZ", "plot_fdigitZ", 1600, -800., 800.);

      //initialise variables
      std::vector<double> fDigitX;
      std::vector<double> fDigitY;
      std::vector<double> fDigitZ;
      std::vector<double> fDigitT;
      std::vector<double> fDigitQ;
      std::vector<double> fDigitPE;
      std::vector<double> fDigitW;
      std::vector<double> fDigitV;
      std::vector<double> fDelta; // time residual 
      std::vector<double> VSeedVtxX;
      std::vector<double> VSeedVtxY;
      std::vector<double> VSeedVtxZ;
      std::vector<double> VSeedVtxTime;
      std::vector<int> vSeedDigitList;

      //WChRecoLite* mRec = WChRecoLite::Instance();

      const double N_REF=1.34; //average index of refraction
      const double C_VAC=29.9792458; //speed of light in vacuum [cm/ns]

      double recoVtxX;
      double recoVtxY;
      double recoVtxZ;
      double recoVtxTime;	
      // ------------------------------------------------------------------------------------------- //	
      // Need to seperate the Inner-Detector tubes from the Outer-Detector tubes
      static const int innerPMTcode = 1;
      static const int vetoPMTcode = 2;
      //
      TFile *f;
      TTree *rat_tree, *run_tree, *data;
      Int_t n_events;
      TTree *run_summary;

      RAT::DS::Root *ds = new RAT::DS::Root();
      RAT::DS::Run *run = new RAT::DS::Run();
      RAT::DS::EV *ev;
      RAT::DS::PMTInfo *pmtinfo;
      RAT::DS::MC *mc;
      RAT::DS::MCParticle *prim;
      RAT::DS::PMT *pmt;
   
      int event, sub_event, n, count;
      int tot_inner, tot_veto, id;
      int hit, nhit, veto_count;
      int inpmt;
      int cables[5000], veto_cables[5000];
      int cables_win[500], veto_cables_win[5000];
      float times[5000], veto_times[5000], offsetT;
      float charges[5000], veto_charges[5000];
      float hitpos_x[5000], hitpos_y[5000], hitpos_z[5000];
      Double_t mc_x = 0., mc_y = 0., mc_z = 0., mc_tim = 0., mc_u = 0., mc_v = 0., mc_w = 0.;
      Int_t mcid = 0, gtid = 0, inner_hit = 0;
      Int_t timestamp_ns = 0, timestamp_s = 0, code = 0;
      Double_t totPE = 0.;
      Double_t timestamp=0., mc_energy = 0.; 

      //Open input file     
       f = new TFile(filename_ratpac); //argv[1]);

       rat_tree = (TTree *)f->Get("T");
       rat_tree->SetBranchAddress("ds", &ds);
       run_tree = (TTree *)f->Get("runT");
       if (rat_tree == 0x0 || run_tree == 0x0) {
          printf("can't find trees T and runT in this file\n");
       return -1;
       }
       run_tree->SetBranchAddress("run", &run);
       if (run_tree->GetEntries() != 1) {
          printf("More than one run! abort...\n");
       return -1;
       }
     
        run_tree->GetEntry(0);

     // loop over PMTs and find positions and location of PMT support
     pmtinfo = run->GetPMTInfo();
     n = pmtinfo->GetPMTCount(); //number of PMTs in the detector
     tot_inner = 0;
     tot_veto = 0;

    // Determines the number of inner and veto pmts
    for (hit = count = 0; hit < n; hit++) {
       if (pmtinfo->GetType(hit) == innerPMTcode)
           ++tot_inner;
       else if (pmtinfo->GetType(hit) == vetoPMTcode)
           ++tot_veto;
      else
          printf("PMT does not have valid identifier: %d \n",
          pmtinfo->GetType(hit));
    }
    if (n != (tot_inner + tot_veto))
        printf("Mis-match in total PMT numbers: %d, %d \n", n,
        tot_inner + tot_veto);
   inpmt = tot_inner;
   printf("In total there are  %d PMTs in WATCHMAN\n", n);
   printf("%d PMTs are inner PMTs\n", inpmt);

   n_events = rat_tree->GetEntries();
   // loop over all events
   for (event = 0; event < n_events; event++) {
   //for (event = 0; event < 10; event++) { //loop over a fixed number of events
        int fThisDigit=0;
        fDigitX.clear();
        fDigitY.clear();
        fDigitZ.clear();
        fDigitT.clear();
        fDigitW.clear();
        fDigitV.clear();
        fDigitQ.clear();
        fDigitPE.clear();
        vSeedDigitList.clear();
        vSeeds.clear();
        fLastEntry_v.clear();
      
    if (event % 100 == 0)
       printf("Evaluating event %d of %d (%d sub events)\n", event, n_events, ds->GetEVCount());
       rat_tree->GetEntry(event);

    // loop over all subevents
    for (sub_event = 0; sub_event < ds->GetEVCount(); sub_event++) {
      gtid += 1;
      mcid = event;

      ev = ds->GetEV(sub_event);
      totPE = ev->GetTotalCharge(); //total charge per event

      TVector3 temp;
      //MC particle information:
      mc = ds->GetMC();
      prim = mc->GetMCParticle(0);
      mc_energy = prim->GetKE();
      temp = prim->GetPosition(); //interaction vertex position
      mc_x = temp.X();
      mc_y = temp.Y();
      mc_z = temp.Z();
      mc_tim = prim->GetTime();                   
      //mc_tim = prim->GetTime() - ev->GetDeltaT();  
      //timestamp_s = ev->GetUTC().GetSec();
      //timestamp_ns = ev->GetUTC().GetNanoSec();
      //timestamp = 1e6*mc->GetUTC().GetSec()+1e-3*mc->GetUTC().GetNanoSec() + 1e-3*ev->GetDeltaT();

      nhit = ev->GetPMTCount(); //number of PMTs with hits per event
      cout<<"event: "<<event<<" sub_event: "<<sub_event<<" number of hits"<<nhit<<" totPE "<<totPE<<endl; 
      // loop over all PMT hits for each subevent
      for (hit = count = veto_count = 0; hit < nhit; hit++) {
        pmt = ev->GetPMT(hit);
        id = pmt->GetID(); //PMT id
        // only use information from the inner pmts
        if (pmtinfo->GetType(id) == innerPMTcode) {
          cables[count] = pmt->GetID() + 1;
          times[count] = pmt->GetTime(); //+ offsetT;
          charges[count] = pmt->GetCharge(); //charge of each PMT

          TVector3 hitpos = pmtinfo->GetPosition(ev->GetPMT(hit)->GetID());
          hitpos_x[count] = hitpos.X()*0.1; //converting to cm
          hitpos_y[count] = hitpos.Y()*0.1; //converting to cm
          hitpos_z[count] = hitpos.Z()*0.1; //converting to cm
          cout<<"Hit: "<<hit<<"|"<<count<<" PMTid "<<cables[count]<<" hitpos "<<hitpos_x[count]<<","<<hitpos_y[count]<<","<<hitpos_z[count]<<" time: "<<times[count]<<" charge: "<<charges[count]<<endl;	  
          
          fDigitX.push_back(hitpos_x[count]); 
          fDigitY.push_back(hitpos_y[count]);
          fDigitZ.push_back(hitpos_z[count]);
          fDigitT.push_back(times[count]);// - min_PE_time; 
          fDigitQ.push_back(charges[count]);
          fDigitW.push_back(500);//photon_wavelength_v[iphot]
          fDigitV.push_back(C_VAC/N_REF);
          vSeedDigitList.push_back(fThisDigit);
          fThisDigit++;

           count++;
        } /*else if (pmtinfo->GetType(id) == vetoPMTcode) {
          veto_cables[veto_count] = pmt->GetID() + 1;
          veto_times[veto_count] = pmt->GetTime() + offsetT;
          veto_charges[veto_count] = pmt->GetCharge();
          veto_count++;
        }*/ else
          printf("Unidentified PMT type: (%d,%d) \n", count,
                 pmtinfo->GetType(id));
      }  // end of loop over all PMT hits
      //veto_hit = veto_count;
      inner_hit = count;
      nhit = count;
      cout<<"inner_hits= "<<inner_hit<<endl;
      cout<<"fThisDigit= "<<fThisDigit<<" count= "<<count<<endl;
      std::cout << "vSeedDigitList:" << vSeedDigitList.size() << std::endl;
      if(vSeedDigitList.size()==0) return -1;

      //calculate seed vertices forposition and timing and store them in vSeeds
       vSeeds = CalcVertexSeeds(fDigitX,fDigitY,fDigitZ,fDigitT,fDigitQ,fDigitW,fDigitV,vSeedDigitList);

      //Get the vertices of seeds:
      VSeedVtxX = GetSeedVtx(vSeeds,0);
      VSeedVtxY = GetSeedVtx(vSeeds,1);
      VSeedVtxZ = GetSeedVtx(vSeeds,2);
      VSeedVtxTime = GetSeedVtx(vSeeds,3);
      
     //now find the best vertex
     int best_seed = SelectBestSeed(fDigitX,fDigitY,fDigitZ,fDigitT,VSeedVtxX,VSeedVtxY,VSeedVtxZ,VSeedVtxTime);
     std::cout<<"best_seed = "<<best_seed<<"  vSeedVtxX.size() = "<<vSeedVtxX.size()<<std::endl;

     //check fdigits_plots
     for(int i=0; i<fDigitX.size(); i++){
        plot_fdigitX->Fill(fDigitX[i]);
     }
     for(int i=0; i<fDigitY.size(); i++){
        plot_fdigitY->Fill(fDigitY[i]);
     }
     for(int i=0; i<fDigitZ.size(); i++){
        plot_fdigitZ->Fill(fDigitZ[i]);
     }

     //filling the csvfile:
     csvfile<<inner_hit<<",";
     csvfile<<mc_x<<",";
     csvfile<<mc_y<<",";
     csvfile<<mc_z<<",";
     csvfile<<VSeedVtxX.size()<<",";
     for(int i=0; i<vSeedVtxX.size(); i++)
     {  csvfile<<VSeedVtxX[i]<<",";      
        plot_VSeedVtxX->Fill((vSeedVtxX[i]));
        plot_VSeedVtxdr->Fill(sqrt((vSeedVtxX[i]-mc_x)*(vSeedVtxX[i]-mc_x) + (vSeedVtxY[i]-mc_y)*(vSeedVtxY[i]-mc_y)));
        plot_VSeedVtxdR->Fill(sqrt((vSeedVtxX[i]-mc_x)*(vSeedVtxX[i]-mc_x) + (vSeedVtxY[i]-mc_y)*(vSeedVtxY[i]-mc_y) + (vSeedVtxZ[i]-mc_z)*(vSeedVtxZ[i]-mc_z)));

     }
     for(int i=0; i<vSeedVtxX.size(); i++)
     { csvfile<<VSeedVtxY[i]<<",";       
       plot_VSeedVtxY->Fill((vSeedVtxY[i]));
     }
     for(int i=0; i<vSeedVtxX.size(); i++)
     { csvfile<<VSeedVtxZ[i]<<",";       
       plot_VSeedVtxZ->Fill((vSeedVtxZ[i]));
     }
     for(int i=0; i<vSeedVtxX.size(); i++)
     { csvfile<<VSeedVtxTime[i]<<",";    }
     csvfile<<'\n';
     //if(best_seed == -1) return -1;
 
     //the reconstructed vertex with this best seed is:
     recoVtxX = VSeedVtxX[best_seed];
     recoVtxY = VSeedVtxY[best_seed];
     recoVtxZ = VSeedVtxZ[best_seed];
     recoVtxTime = VSeedVtxTime[best_seed];
     cout<<"mc_x= "<<mc_x<<" mc_y= "<<mc_y<<" mc_z= "<<mc_z<<" mc_time= "<<mc_tim<<endl;
     cout<<"recoVtxX= "<<recoVtxX<<" recoVtxY= "<<recoVtxY<<" recoVtxZ= "<<recoVtxZ<<" recoVtxTime= "<<recoVtxTime<<endl; 
     cout<<"dR= "<<sqrt((recoVtxX-mc_x)*(recoVtxX-mc_x) + (recoVtxY-mc_y)*(recoVtxY-mc_y) + (recoVtxZ-mc_z)*(recoVtxZ-mc_z))<<endl;
     plot->Fill(sqrt((recoVtxX-mc_x)*(recoVtxX-mc_x) + (recoVtxY-mc_y)*(recoVtxY-mc_y) + (recoVtxZ-mc_z)*(recoVtxZ-mc_z)));
   }
  }//closes the event loop              

  csvfile.close();
 /* new TCanvas();
  plot_fdigitX->Draw();
  new TCanvas();
  plot_fdigitY->Draw();
  new TCanvas();
  plot_fdigitZ->Draw();
*/
  new TCanvas();
  plot->Draw();
  new TCanvas();
  plot_VSeedVtxX->Draw();
  new TCanvas();
  plot_VSeedVtxY->Draw();
  new TCanvas();
  plot_VSeedVtxZ->Draw();
  new TCanvas();
  plot_VSeedVtxdr->Draw();
  new TCanvas();
  plot_VSeedVtxdR->Draw();

}
