/*
** @author Krista McCord, 2016 - @Revised Shahram Talei 2017
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"
#include <gsl/gsl_sf_bessel.h>


void rotation_angles(int galaxy)
{

   float Jmag, Jmag2;
   float thetay;
   float thetaz;
   float Jx, Jy, Jz, Jxp, Jyp, Jzp, Jxpp, Jypp, Jzpp;
   
#ifdef Extrap
   if(snap_oldgal_count > 1)
   {
      Jx = AllGal[galaxy].NewJ[0];
      Jy = AllGal[galaxy].NewJ[1];
      Jz = AllGal[galaxy].NewJ[2];
   }
   else
   {
   
      Jx = AllGal[galaxy].J[0];
      Jy = AllGal[galaxy].J[1];
      Jz = AllGal[galaxy].J[2];
   }
#else

   Jx = AllGal[galaxy].J[0];
   Jy = AllGal[galaxy].J[1];
   Jz = AllGal[galaxy].J[2];
#endif
   

   Jmag = 0.0;
   Jmag2 = 0.0;
   

//Calculate Rotation angles about z and y axis to have J parallel to Z-axis
   


//Rotate about z axis based on Jx and Jy

   Jmag = sqrt(Jx*Jx + Jy*Jy + Jz*Jz);
   
   thetaz = fabs(atanf(Jy / Jx));

   if((Jx < 0.0 && Jy < 0.0) || (Jx >= 0.0 && Jy >= 0.0))
   {
      Jxp = Jx * cosf(thetaz) - Jy * sinf(thetaz);
      Jyp = Jx * sinf(thetaz) + Jy * cosf(thetaz);
      Jzp = Jz;

   }
   else if ((Jx >= 0.0 && Jy < 0.0) || (Jx < 0.0 && Jy >= 0.0))
   {
      Jxp = Jx * cosf(thetaz) + Jy * sinf(thetaz);
      Jyp = Jy * cosf(thetaz) - Jx * sinf(thetaz);
      Jzp = Jz;

   }


//Rotate about y axis based on new Jx, Jy, Jz   

   Jmag2 = sqrt(Jxp * Jxp + Jyp * Jyp + Jzp * Jzp);
 
//Dot product of J vector to z axis (0,0,1) to find angle of rotation about y axis

   thetay = acosf(Jz / Jmag2);

   if(Jxp >= 0.0)
   {
     Jxpp = Jxp * cosf(thetay) - Jzp * sinf(thetay);
     Jypp = Jyp;
     Jzpp = Jxp * sinf(thetay) + Jzp * cosf(thetay);

   }
   else if (Jxp < 0.0)
   {
     Jxpp = Jxp * cosf(thetay) + Jzp * sinf(thetay);
     Jypp = Jyp;
     Jzpp = Jzp * cosf(thetay) - Jxp * sinf(thetay);

   }
  
   AllGal[galaxy].angley = thetay;
   AllGal[galaxy].anglez = thetaz;

return;
}


void RotateR(int Pindex, int gal)
{
   
   float a;
   float rxp, ryp, rzp;
   float thetaz, thetay;
   float rxpp, rypp, rzpp;
   float rx, ry, rz, Jx, Jy, Jz;
 
   thetaz = AllGal[gal].anglez;
   thetay = AllGal[gal].angley;
   a = -1.0;
   rx = P[Pindex].Pos[0];
   ry = P[Pindex].Pos[1];
   rz = P[Pindex].Pos[2];

#ifdef Extrap
   if(snap_oldgal_count > 1)
   {
      Jx = AllGal[gal].NewJ[0];
      Jy = AllGal[gal].NewJ[1];
      Jz = AllGal[gal].NewJ[2];
   }
   else
   {

      Jx = AllGal[gal].J[0];
      Jy = AllGal[gal].J[1];
      Jz = AllGal[gal].J[2];
   }
#else

   Jx = AllGal[gal].J[0];
   Jy = AllGal[gal].J[1];
   Jz = AllGal[gal].J[2];

#endif
 if((Jx < 0.0 && Jy < 0.0) || (Jx >= 0.0 && Jy >= 0.0))
   {
      rxp = rx * cosf(thetaz) - ry * sinf(thetaz);
      ryp = rx * sinf(thetaz) + ry * cosf(thetaz);
      rzp = rz;

   }
   else if ((Jx >= 0.0 && Jy < 0.0) || (Jx < 0.0 && Jy >= 0.0))
   {
      rxp = rx * cosf(thetaz) + ry * sinf(thetaz);
      ryp = ry * cosf(thetaz) - rx * sinf(thetaz);
      rzp = rz;

   }


  if(Jx >= 0.0)
   {
     rxpp = rxp * cosf(thetay) - rzp * sinf(thetay);
     rypp = ryp;
     rzpp = rxp * sinf(thetay) + rzp * cosf(thetay);

   }
   else if (Jx < 0.0)
   {
     rxpp = rxp * cosf(thetay) + rzp * sinf(thetay);
     rypp = ryp;
     rzpp = rzp * cosf(thetay) - rxp * sinf(thetay);

   }

   NewP.Px = rxpp;
   NewP.Py = rypp;
   NewP.Pz = rzpp;
   
  
return;
}

void Final_Disk_Forces(int target, int gal)
{
   double DiskRad, DiskRadx, DiskRady;
   double DiskVerz, diskR;
   double rx,ry,rz;
   double Px, Py, Pz;
   double Softhalo;   
   Px = NewP.Px;
   Py = NewP.Py;
   Pz = NewP.Pz;

   Softhalo = All.SofteningHalo / All.Time;

//#ifdef Extrap
   rx = -1.0 *( Px - AllGal[gal].CM_Pos[0]);
   ry = -1.0*(Py - AllGal[gal].CM_Pos[1]);
   rz = Pz - AllGal[gal].CM_Pos[2];
//#else

//   rx =  -1.0 *( Px - AllGal[gal].Pos[0]);
//   ry = -1.0*(Py - AllGal[gal].Pos[1]);
//   rz = Pz - AllGal[gal].Pos[2];

//#endif
  
   diskR = sqrt(rx*rx + ry*ry + Softhalo * Softhalo);        //All.SofteningHalo * All.SofteningHalo);

   DiskRad = Calculate_Disk_RadF(target, gal);
	//printf("DO211\n");
   DiskVerz = Calculate_Disk_VerF(target, gal);
	// test output
   //printf("DO213\n");
   DiskRadx = (rx / diskR) * DiskRad;
   DiskRady = (ry/ diskR) * DiskRad;
  
   Rotate_Forces(DiskVerz, DiskRadx, DiskRady, gal);

return;
}


double Calculate_Disk_RadF(int target, int gal)
{

    double diskmass, diskR, diskR2, y, surfaceden,zo;
    double ycheck, FDiskRad;
    double rx, ry, rz;
    double diskscaler;
    float Px, Py, Pz;
    double xmin;
    double Softhalo;
	double M1,M2,M3,a1,a2,a3,b,x_MN;
    double x3,x4;


    Softhalo = All.SofteningHalo / All.Time;

#ifdef Extrap
    if(snap_oldgal_count > 1)
    {
       diskmass = (AllGal[gal].NewStellarMass - AllGal[gal].NewBulgeMass) + AllGal[gal].NewColdGas;
       diskscaler =  AllGal[gal].NewDiskScaleRadius / All.Time;
    }
    else
    {
       diskmass = (AllGal[gal].StellarMass - AllGal[gal].BulgeMass) + AllGal[gal].ColdGas;
       diskscaler =  AllGal[gal].DiskScaleRadius / All.Time;
    }
#else

    diskmass = (AllGal[gal].StellarMass - AllGal[gal].BulgeMass) + AllGal[gal].ColdGas;
    diskscaler =  AllGal[gal].DiskScaleRadius / All.Time;

#endif


    xmin = DBL_MIN; 
    diskR = 0;
    diskR2 = 0;
    y = 0;
    surfaceden = 0;
    ycheck = 0;
    FDiskRad = 0;
    Px = NewP.Px;
    Py = NewP.Py;
    Pz = NewP.Pz;

 
//#ifdef Extrap
    rx = -1.0 *( Px - AllGal[gal].CM_Pos[0]);
    ry = -1.0*(Py - AllGal[gal].CM_Pos[1]);
    rz = Pz - AllGal[gal].CM_Pos[2];
//#else

  //  rx = -1.0 *( Px - AllGal[gal].Pos[0]);
   // ry = -1.0*(Py - AllGal[gal].Pos[1]);
 //   rz = Pz - AllGal[gal].Pos[2];


//#endif

    diskR = sqrt(rx*rx + ry*ry + Softhalo * Softhalo);         //All.SofteningHalo * All.SofteningHalo);
    diskR2 = sqrt(rx*rx + ry*ry);

    surfaceden = diskmass / (M_PI * diskscaler*diskscaler);

    y = diskR / (2 * diskscaler);
     zo =diskscaler/9.4;
	b=fabs(zo);
        x_MN=b/diskscaler;
        x3=pow(x_MN,3.0);
        x4=x3*x_MN;
        //this potential is taken from Smith et al 2015
        M1=(0.0036*x4-0.0330*x3+0.1117*x_MN*x_MN-0.1335*x_MN+0.1749)*diskmass;
        M2=(-0.0131*x4+0.1090*x3-0.3035*x_MN*x_MN+0.2921*x_MN-5.7976)*diskmass;
        M3=(-0.0048*x4+0.0454*x3-0.1425*x_MN*x_MN+0.1012*x_MN+6.7120)*diskmass;
        a1=(-0.0158*x4+0.0993*x3-0.2070*x_MN*x_MN-0.7089*x_MN+0.6445)*diskscaler;
        a2=(-0.0319*x4+0.1514*x3-0.1279*x_MN*x_MN-0.9325*x_MN+2.6836)*diskscaler;
        a3=(-0.0326*x4+0.1816*x3-0.2943*x_MN*x_MN-0.6329*x_MN+2.3193)*diskscaler;
        //                                                        //
        //



    if(diskR <= 200.0*diskscaler)
    {

        ycheck = fabs(y);

        if(ycheck  > xmin)
        {
            //FDiskRad = (4 * M_PI  * surfaceden * (diskscaler/diskR) * y * y * (gsl_sf_bessel_I0_scaled(y) * gsl_sf_bessel_K0_scaled(y) - gsl_sf_bessel_I1_scaled(y) * gsl_sf_bessel_K1_scaled(y))) / All.HubbleParam;
            FDiskRad=Miyamoto_Nagai_Rforce(diskR,rz,M1,a1,b)+Miyamoto_Nagai_Rforce(diskR,rz,M2,a2,b)+Miyamoto_Nagai_Rforce(diskR,rz,M3,a3,b);

        }
       
        else
        {

            FDiskRad = 0;
           
        }
          
   

     }

   
return FDiskRad;
}



double Calculate_Disk_VerF(int target, int gal)
{
    double diskmass, diskR, diskR2, y, surfaceden;
    double ycheck;
    double rx, ry, rz;
    double dirz;
    double diskscaler;
    double Px, Py, Pz;
    double z1, ztot, zsoft, zo;
    double diskz, FDiskVer, Fvcorr, Fvdamp;
    double Softhalo;
    double M1,M2,M3,a1,a2,a3,b,x_MN;
    double x3,x4;
    Softhalo = All.SofteningHalo / All.Time;    


#ifdef Extrap
    if(snap_oldgal_count > 1)
    {
       diskmass = (AllGal[gal].NewStellarMass - AllGal[gal].NewBulgeMass) + AllGal[gal].NewColdGas;
       diskscaler =  AllGal[gal].NewDiskScaleRadius / All.Time;
    }
    else
    {
       diskmass = (AllGal[gal].StellarMass - AllGal[gal].BulgeMass) + AllGal[gal].ColdGas;
       diskscaler =  AllGal[gal].DiskScaleRadius / All.Time;
    }
#else

    diskmass = (AllGal[gal].StellarMass - AllGal[gal].BulgeMass) + AllGal[gal].ColdGas;
    diskscaler =  AllGal[gal].DiskScaleRadius / All.Time;

#endif

    dirz = -1.0;
    Fvcorr = 1.0;
    Fvdamp = 1.0;
    diskR = 0;
    diskR2 = 0;
    y = 0;
    surfaceden = 0;
    ycheck = 0;
    FDiskVer = 0;
    Px = NewP.Px;
    Py = NewP.Py;
    Pz = NewP.Pz;

//#ifdef Extrap

    rx = -1.0 *( Px - AllGal[gal].CM_Pos[0]);
    ry = -1.0*(Py - AllGal[gal].CM_Pos[1]);
    rz = Pz - AllGal[gal].CM_Pos[2];
//#else

 //   rx = -1.0 *( Px - AllGal[gal].Pos[0]);
 //   ry = -1.0*(Py - AllGal[gal].Pos[1]);
 //   rz = Pz - AllGal[gal].Pos[2];

//#endif

    diskR = sqrt(rx*rx + ry*ry + Softhalo * Softhalo);           //All.SofteningHalo * All.SofteningHalo);
    diskR2 = sqrt(rx*rx + ry*ry);

    surfaceden = diskmass / (M_PI * diskscaler*diskscaler);

    diskz = sqrt(rx*rx + ry*ry + rz*rz);

   
       
     ztot = fabs(rz);  //z1 - AllGal[gal].Pos[2];
       

   // z1 = sqrt(fabs(diskz*diskz - diskR2*diskR2));
    
  //  ztot = z1 - AllGal[gal].Pos[2];
         
    zsoft = sqrt(ztot*ztot + Softhalo * Softhalo);      //All.SofteningHalo*All.SofteningHalo);

    //zo = 0.7 / All.Time;  //scale height of disk .2kpc
zo =diskscaler/9.4;
   x_MN=b/diskscaler;
        x3=pow(x_MN,3.0);
        x4=x3*x_MN;
        //this potential is taken from Smith et al 2015
         M1=(0.0036*x4-0.0330*x3+0.1117*x_MN*x_MN-0.1335*x_MN+0.1749)*diskmass;
         M2=(-0.0131*x4+0.1090*x3-0.3035*x_MN*x_MN+0.2921*x_MN-5.7976)*diskmass;
         M3=(-0.0048*x4+0.0454*x3-0.1425*x_MN*x_MN+0.1012*x_MN+6.7120)*diskmass;
         a1=(-0.0158*x4+0.0993*x3-0.2070*x_MN*x_MN-0.7089*x_MN+0.6445)*diskscaler;
         a2=(-0.0319*x4+0.1514*x3-0.1279*x_MN*x_MN-0.9325*x_MN+2.6836)*diskscaler;
         a3=(-0.0326*x4+0.1816*x3-0.2943*x_MN*x_MN-0.6329*x_MN+2.3193)*diskscaler;
    //FDiskVer = fabs((2 * M_PI  * surfaceden * exp(-diskR / diskscaler) * ztot * (1 - exp(-zsoft / zo))) /  (All.HubbleParam * zsoft));
	FDiskVer=fabs(Miyamoto_Nagai_zforce(diskR,rz,M1,a1,b)+Miyamoto_Nagai_zforce(diskR,rz,M2,a2,b)+Miyamoto_Nagai_zforce(diskR,rz,M3,a3,b));	

     /***Vertical force correction for particles within certain regimes where this particle force equation assumption fails
         These particles end up with too large of a verticle force when the cylindrical radius becomes small*******/

       //Fvcorr = 1 / pow((1 + pow(fabs(rz) / (1.7*diskscaler),1.65)), 1.212);

       /***The vertical force is too high for a small percentage of particles that come close to a particular regime.  A damping
          value is applied to slow down these particles.  This effect seems to be a characteristic of this force equation***/

       //Fvdamp =1;// exp(-0.049729 / (diskR*diskR));                      //-1.0*pow(0.223,2.0) / pow(diskR,2.0) );



    if(rz < 0.0)
    {
        dirz = 1.0;
    }
    else
    {
        dirz = -1.0;
    }

    //FDiskVer = FDiskVer * dirz * Fvcorr * Fvdamp;
      FDiskVer = FDiskVer * dirz;// * Fvcorr;
	//FDiskVer=ForceDiskVer* dirz;
return FDiskVer;
}

void Rotate_Forces(double Diskz, double Diskx, double Disky, int gal)
{

   double Fxpp, Fypp, Fzpp;
   double Fxp, Fyp, Fzp;
   double thetay, thetaz;
   double Jx, Jy, Jz;

#ifdef Extrap
   if(snap_oldgal_count > 1)
   {
      Jx = AllGal[gal].NewJ[0];
      Jy = AllGal[gal].NewJ[1];
      Jz = AllGal[gal].NewJ[2];
   }
   else
   {

      Jx = AllGal[gal].J[0];
      Jy = AllGal[gal].J[1];
      Jz = AllGal[gal].J[2];
   }
#else

  Jx = AllGal[gal].J[0];
  Jy = AllGal[gal].J[1];
  Jz = AllGal[gal].J[2];

#endif

   thetaz = AllGal[gal].anglez;
   thetay = AllGal[gal].angley;

   if(Jx >= 0.0)
   {
     Fxp = Diskx * cosf(thetay) + Diskz * sinf(thetay);
     Fyp = Disky;
     Fzp = Diskz * cosf(thetay) - Diskx * sinf(thetay);

   }
   else if (Jx < 0.0)
   {
     Fxp = Diskx * cosf(thetay) - Diskz * sinf(thetay);
     Fyp = Disky;
     Fzp = Diskz * cosf(thetay) + Diskx * sinf(thetay);

   }


   if((Jx < 0.0 && Jy < 0.0) || (Jx >= 0.0 && Jy >= 0.0))
   {

      Fxpp = Fxp * cosf(thetaz) + Fyp * sinf(thetaz);
      Fypp = Fyp * cosf(thetaz) - Fxp * sinf(thetaz);
      Fzpp = Fzp;
      

   }
   else if ((Jx >= 0.0 && Jy < 0.0) || (Jx < 0.0 && Jy >= 0.0))
   {

      Fxpp = Fxp * cosf(thetaz) - Fyp * sinf(thetaz);
      Fypp = Fxp * sinf(thetaz) + Fyp * cosf(thetaz);
      Fzpp = Fzp;

   }

 
   DiskF.Fx = Fxpp;
   DiskF.Fy = Fypp;
   DiskF.Fz = Fzpp;

return;
}





