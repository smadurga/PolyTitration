// Module with subroutines to perform a certain rotation (in degrees) in each axis given the cartesian coordinates.

#include "lib.h"
#include "global.h"

/******************************************************/


void gencoordp(int nrot, double ang_Rot)
{
double ab_un[1][3], bc_un[1][3], mod_ab, mod_bc, cross_ab3[1][3], mod_cross, cross_2[1][3], rot_M[3][3], pre_rot[1][3], disp[1][3];
int i;
double alea;
double phi, theta, dang=2, dang2=0.1;

// adapted from python subroutine 
// https://github.com/rbhatta8/protein-design/blob/master/nerf.py


// Provisional Phi of the rotated monomer

Phi_p=Phi[nrot]+ang_Rot;


coordp[1][0]=0;
coordp[1][1]=0;
coordp[1][2]=0;

if (nrot == 1){

alea=ran2(&idum);

if (alea < 0.5) { 

Theta_p=Angles[0]+dang;



}
if (alea > 0.5) { 

Theta_p=Angles[0]-dang;



}
if (alea == 0.5) { Theta_p=Angles[0];}

alea=ran2(&idum);

if (alea < 0.5) {
 Phi_p=acos(Phi[1]+dang2);
if ( Phi[1]+dang2 > 1 ) { Phi_p=acos(Phi[1]-dang2);}
if (Phi[1]+dang2 == 1)  { Phi_p=0;} 

}
if (alea > 0.5) { 

Phi_p=acos(Phi[1]-dang2);
if ( Phi[1]-dang2 < -1 ) { Phi_p=acos(Phi[1]+dang2);}
if (Phi[1]-dang2 == -1)  { Phi_p=PI;}


}



if (alea == 0.5) { Phi_p=acos(Phi[1]);}

theta=Theta_p/360*2*PI;
phi=Phi_p;

coordp[2][0]=Bounds[1]*sin(phi)*cos(theta);
coordp[2][1]=Bounds[1]*sin(phi)*sin(theta);
coordp[2][2]=Bounds[1]*cos(phi);

}
else {

coordp[2][0]=coord[2][0];
coordp[2][1]=coord[2][1];
coordp[2][2]=coord[2][2];

}

if (nrot == 2){

alea=ran2(&idum);

if (alea < 0.5) {
 Phi_p=acos(Phi[2]+dang2);
if ( Phi[2]+dang2 > 1 ) { Phi_p=acos(Phi[2]-dang2);}
if (Phi[2]+dang2 == 1)  { Phi_p=0;}

}
if (alea > 0.5) {

Phi_p=acos(Phi[2]-dang2);
if ( Phi[2]-dang2 < -1 ) { Phi_p=acos(Phi[2]+dang2);}
if (Phi[2]-dang2 == -1)  { Phi_p=PI;}




}


theta=Angles[1]/360*2*PI;
phi=Phi_p;//360*2*PI; 

//alea=ran2(&idum);

//if (alea <= 0.5) {phi=Phi_p;} else { phi=2*PI - Phi_p;}



// Unit vector between atom i-2 and i-3


ab_un[0][0]=0;
ab_un[0][1]=0;
ab_un[0][2]=1;

if (coordp[2][2] == Bounds[1]) { ab_un[0][1]=0.044710; ab_un[0][2]=0.999; }
if (coordp[2][2] == -Bounds[1]) { ab_un[0][1]=0.044710; ab_un[0][2]=0.999;}



// Unit vector between atom i-1 and i-2

mod_bc=pow(pow(coordp[2][0]-coordp[1][0],2)+pow(coordp[2][1]-coordp[1][1],2)+pow(coordp[2][2]-coordp[1][2],2),1./2.);

bc_un[0][0]=(coordp[2][0]-coordp[1][0])/mod_bc;
bc_un[0][1]=(coordp[2][1]-coordp[1][1])/mod_bc;
bc_un[0][2]=(coordp[2][2]-coordp[1][2])/mod_bc;

// Cross product of vector ab and vector bc

cross_ab3[0][0]=ab_un[0][1]*bc_un[0][2]-ab_un[0][2]*bc_un[0][1];
cross_ab3[0][1]=ab_un[0][2]*bc_un[0][0]-ab_un[0][0]*bc_un[0][2];
cross_ab3[0][2]=ab_un[0][0]*bc_un[0][1]-ab_un[0][1]*bc_un[0][0];

// Normalize the resulting vector

mod_cross=pow(pow(cross_ab3[0][0],2)+pow(cross_ab3[0][1],2)+pow(cross_ab3[0][2],2),1./2.);

cross_ab3[0][0]=cross_ab3[0][0]/mod_cross;
cross_ab3[0][1]=cross_ab3[0][1]/mod_cross;
cross_ab3[0][2]=cross_ab3[0][2]/mod_cross;


// Cross product between this vector and vector bc

cross_2[0][0]=cross_ab3[0][1]*bc_un[0][2]-cross_ab3[0][2]*bc_un[0][1];
cross_2[0][1]=cross_ab3[0][2]*bc_un[0][0]-cross_ab3[0][0]*bc_un[0][2];
cross_2[0][2]=cross_ab3[0][0]*bc_un[0][1]-cross_ab3[0][1]*bc_un[0][0];

// create rotation matrix [vbc; cross_2; cross_ab3]

rot_M[0][0]=bc_un[0][0]; rot_M[0][1]=cross_2[0][0]; rot_M[0][2]=cross_ab3[0][0];
rot_M[1][0]=bc_un[0][1]; rot_M[1][1]=cross_2[0][1]; rot_M[1][2]=cross_ab3[0][1];
rot_M[2][0]=bc_un[0][2]; rot_M[2][1]=cross_2[0][2]; rot_M[2][2]=cross_ab3[0][2];


// calculate coord pre rotation row

pre_rot[0][0]= -Bounds[2]*cos(theta);
pre_rot[0][1]= Bounds[2]*sin(theta)*cos(phi);
pre_rot[0][2]= Bounds[2]*sin(theta)*sin(phi);

// displacement vector rot_M*pre_rot

disp[0][0]= rot_M[0][0]*pre_rot[0][0]+rot_M[0][1]*pre_rot[0][1]+rot_M[0][2]*pre_rot[0][2];
disp[0][1]= rot_M[1][0]*pre_rot[0][0]+rot_M[1][1]*pre_rot[0][1]+rot_M[1][2]*pre_rot[0][2];
disp[0][2]= rot_M[2][0]*pre_rot[0][0]+rot_M[2][1]*pre_rot[0][1]+rot_M[2][2]*pre_rot[0][2];

// Final coordinates of new atom

coordp[3][0]=coordp[2][0]+disp[0][0];
coordp[3][1]=coordp[2][1]+disp[0][1];
coordp[3][2]=coordp[2][2]+disp[0][2];



}

else{


theta=Angles[1]/360*2*PI;
phi=acos(Phi[2]); //360*2*PI; 

// Unit vector between atom i-2 and i-3


ab_un[0][0]=0;
ab_un[0][1]=0;
ab_un[0][2]=1;

if (coordp[2][2] == Bounds[1]) { ab_un[0][1]=0.044710; ab_un[0][2]=0.999;}
if (coordp[2][2] == -Bounds[1]) { ab_un[0][1]=0.044710; ab_un[0][2]=0.999;}


// Unit vector between atom i-1 and i-2

mod_bc=pow(pow(coordp[2][0]-coordp[1][0],2)+pow(coordp[2][1]-coordp[1][1],2)+pow(coordp[2][2]-coordp[1][2],2),1./2.);

bc_un[0][0]=(coordp[2][0]-coordp[1][0])/mod_bc;
bc_un[0][1]=(coordp[2][1]-coordp[1][1])/mod_bc;
bc_un[0][2]=(coordp[2][2]-coordp[1][2])/mod_bc;

// Cross product of vector ab and vector bc

cross_ab3[0][0]=ab_un[0][1]*bc_un[0][2]-ab_un[0][2]*bc_un[0][1];
cross_ab3[0][1]=ab_un[0][2]*bc_un[0][0]-ab_un[0][0]*bc_un[0][2];
cross_ab3[0][2]=ab_un[0][0]*bc_un[0][1]-ab_un[0][1]*bc_un[0][0];

// Normalize the resulting vector

mod_cross=pow(pow(cross_ab3[0][0],2)+pow(cross_ab3[0][1],2)+pow(cross_ab3[0][2],2),1./2.);

cross_ab3[0][0]=cross_ab3[0][0]/mod_cross;
cross_ab3[0][1]=cross_ab3[0][1]/mod_cross;
cross_ab3[0][2]=cross_ab3[0][2]/mod_cross;


// Cross product between this vector and vector bc

cross_2[0][0]=cross_ab3[0][1]*bc_un[0][2]-cross_ab3[0][2]*bc_un[0][1];
cross_2[0][1]=cross_ab3[0][2]*bc_un[0][0]-cross_ab3[0][0]*bc_un[0][2];
cross_2[0][2]=cross_ab3[0][0]*bc_un[0][1]-cross_ab3[0][1]*bc_un[0][0];

// create rotation matrix [vbc; cross_2; cross_ab3]

rot_M[0][0]=bc_un[0][0]; rot_M[0][1]=cross_2[0][0]; rot_M[0][2]=cross_ab3[0][0];
rot_M[1][0]=bc_un[0][1]; rot_M[1][1]=cross_2[0][1]; rot_M[1][2]=cross_ab3[0][1];
rot_M[2][0]=bc_un[0][2]; rot_M[2][1]=cross_2[0][2]; rot_M[2][2]=cross_ab3[0][2];


// calculate coord pre rotation row

pre_rot[0][0]= -Bounds[2]*cos(theta);
pre_rot[0][1]= Bounds[2]*sin(theta)*cos(phi);
pre_rot[0][2]= Bounds[2]*sin(theta)*sin(phi);

// displacement vector rot_M*pre_rot

disp[0][0]= rot_M[0][0]*pre_rot[0][0]+rot_M[0][1]*pre_rot[0][1]+rot_M[0][2]*pre_rot[0][2];
disp[0][1]= rot_M[1][0]*pre_rot[0][0]+rot_M[1][1]*pre_rot[0][1]+rot_M[1][2]*pre_rot[0][2];
disp[0][2]= rot_M[2][0]*pre_rot[0][0]+rot_M[2][1]*pre_rot[0][1]+rot_M[2][2]*pre_rot[0][2];

// Final coordinates of new atom

coordp[3][0]=coordp[2][0]+disp[0][0];
coordp[3][1]=coordp[2][1]+disp[0][1];
coordp[3][2]=coordp[2][2]+disp[0][2];




}



// Transformation to radians

for (i=4; i <= Natom; i++){

theta=Angles[i-2]/360*2*PI;

if ( i-1 != nrot){

  phi=Phi[i-1]/360*2*PI; // Phi[i] stores the internal angle between i+1,i,i-1 and i-2

}else{

  phi=Phi_p/360*2*PI;

}


// Unit vector between atom i-2 and i-3

mod_ab=pow(pow(coordp[i-2][0]-coordp[i-3][0],2)+pow(coordp[i-2][1]-coordp[i-3][1],2)+pow(coordp[i-2][2]-coordp[i-3][2],2),1./2.);

ab_un[0][0]=(coordp[i-2][0]-coordp[i-3][0])/mod_ab;
ab_un[0][1]=(coordp[i-2][1]-coordp[i-3][1])/mod_ab;
ab_un[0][2]=(coordp[i-2][2]-coordp[i-3][2])/mod_ab;

// Unit vector between atom i-1 and i-2

mod_bc=pow(pow(coordp[i-1][0]-coordp[i-2][0],2)+pow(coordp[i-1][1]-coordp[i-2][1],2)+pow(coordp[i-1][2]-coordp[i-2][2],2),1./2.);

bc_un[0][0]=(coordp[i-1][0]-coordp[i-2][0])/mod_bc;
bc_un[0][1]=(coordp[i-1][1]-coordp[i-2][1])/mod_bc;
bc_un[0][2]=(coordp[i-1][2]-coordp[i-2][2])/mod_bc;

// Cross product of vector ab and vector bc

cross_ab3[0][0]=ab_un[0][1]*bc_un[0][2]-ab_un[0][2]*bc_un[0][1];
cross_ab3[0][1]=ab_un[0][2]*bc_un[0][0]-ab_un[0][0]*bc_un[0][2];
cross_ab3[0][2]=ab_un[0][0]*bc_un[0][1]-ab_un[0][1]*bc_un[0][0];

// Normalize the resulting vector

mod_cross=pow(pow(cross_ab3[0][0],2)+pow(cross_ab3[0][1],2)+pow(cross_ab3[0][2],2),1./2.);

cross_ab3[0][0]=cross_ab3[0][0]/mod_cross;
cross_ab3[0][1]=cross_ab3[0][1]/mod_cross;
cross_ab3[0][2]=cross_ab3[0][2]/mod_cross;


// Cross product between this vector and vector bc

cross_2[0][0]=cross_ab3[0][1]*bc_un[0][2]-cross_ab3[0][2]*bc_un[0][1];
cross_2[0][1]=cross_ab3[0][2]*bc_un[0][0]-cross_ab3[0][0]*bc_un[0][2];
cross_2[0][2]=cross_ab3[0][0]*bc_un[0][1]-cross_ab3[0][1]*bc_un[0][0];

// create rotation matrix [vbc; cross_2; cross_ab3]

rot_M[0][0]=bc_un[0][0]; rot_M[0][1]=cross_2[0][0]; rot_M[0][2]=cross_ab3[0][0];
rot_M[1][0]=bc_un[0][1]; rot_M[1][1]=cross_2[0][1]; rot_M[1][2]=cross_ab3[0][1];
rot_M[2][0]=bc_un[0][2]; rot_M[2][1]=cross_2[0][2]; rot_M[2][2]=cross_ab3[0][2];


// calculate coord pre rotation row

pre_rot[0][0]= -Bounds[i-1]*cos(theta);
pre_rot[0][1]= Bounds[i-1]*sin(theta)*cos(phi);
pre_rot[0][2]= Bounds[i-1]*sin(theta)*sin(phi);

// displacement vector rot_M*pre_rot

disp[0][0]= rot_M[0][0]*pre_rot[0][0]+rot_M[0][1]*pre_rot[0][1]+rot_M[0][2]*pre_rot[0][2];
disp[0][1]= rot_M[1][0]*pre_rot[0][0]+rot_M[1][1]*pre_rot[0][1]+rot_M[1][2]*pre_rot[0][2];
disp[0][2]= rot_M[2][0]*pre_rot[0][0]+rot_M[2][1]*pre_rot[0][1]+rot_M[2][2]*pre_rot[0][2];

// Final coordinates of new atom

coordp[i][0]=coordp[i-1][0]+disp[0][0];
coordp[i][1]=coordp[i-1][1]+disp[0][1];
coordp[i][2]=coordp[i-1][2]+disp[0][2];

}

}

void gencoordpstr(int nstr)
{
double ab_un[1][3], bc_un[1][3], mod_ab, mod_bc, cross_ab3[1][3], mod_cross, cross_2[1][3], rot_M[3][3], pre_rot[1][3], disp[1][3];
int i,signe;
double alea, theta, phi;
double  dl=0.01;

// adapted from python subroutine 
// https://github.com/rbhatta8/protein-design/blob/master/nerf.py


// Provisional length of the stretched bond

alea=ran2(&idum);
if (alea <= 0.5) {signe=1;} else {signe=-1;}

p_length= Bounds[nstr]+dl*signe;

coordp[1][0]=0;
coordp[1][1]=0;
coordp[1][2]=0;

if (nstr == 1){


      coordp[2][0]=p_length*sin(acos(Phi[1]))*cos(Angles[0]/360*2*PI);
      coordp[2][1]=p_length*sin(acos(Phi[1]))*sin(Angles[0]/360*2*PI);
      coordp[2][2]=p_length*cos(acos(Phi[1]));
}
else{
      coordp[2][0]=Bounds[1]*sin(acos(Phi[1]))*cos(Angles[0]/360*2*PI);
      coordp[2][1]=Bounds[1]*sin(acos(Phi[1]))*sin(Angles[0]/360*2*PI);
      coordp[2][2]=Bounds[1]*cos(acos(Phi[1]));

}

if (nstr == 2){


theta=Angles[1]/360*2*PI;
phi=acos(Phi[2]);

ab_un[0][0]=0;
ab_un[0][1]=0;
ab_un[0][2]=1;

if (coordp[2][2] == p_length) { ab_un[0][1]=0.044710; ab_un[0][2]=0.999;}
if (coordp[2][2] == -p_length) { ab_un[0][1]=0.044710; ab_un[0][2]=0.999;}


// Unit vector between atom i-1 and i-2

mod_bc=pow(pow(coordp[2][0]-coordp[1][0],2)+pow(coordp[2][1]-coordp[1][1],2)+pow(coordp[2][2]-coordp[1][2],2),1./2.);

bc_un[0][0]=(coordp[2][0]-coordp[1][0])/mod_bc;
bc_un[0][1]=(coordp[2][1]-coordp[1][1])/mod_bc;
bc_un[0][2]=(coordp[2][2]-coordp[1][2])/mod_bc;

// Cross product of vector ab and vector bc

cross_ab3[0][0]=ab_un[0][1]*bc_un[0][2]-ab_un[0][2]*bc_un[0][1];
cross_ab3[0][1]=ab_un[0][2]*bc_un[0][0]-ab_un[0][0]*bc_un[0][2];
cross_ab3[0][2]=ab_un[0][0]*bc_un[0][1]-ab_un[0][1]*bc_un[0][0];

// Normalize the resulting vector

mod_cross=pow(pow(cross_ab3[0][0],2)+pow(cross_ab3[0][1],2)+pow(cross_ab3[0][2],2),1./2.);

cross_ab3[0][0]=cross_ab3[0][0]/mod_cross;
cross_ab3[0][1]=cross_ab3[0][1]/mod_cross;
cross_ab3[0][2]=cross_ab3[0][2]/mod_cross;


// Cross product between this vector and vector bc

cross_2[0][0]=cross_ab3[0][1]*bc_un[0][2]-cross_ab3[0][2]*bc_un[0][1];
cross_2[0][1]=cross_ab3[0][2]*bc_un[0][0]-cross_ab3[0][0]*bc_un[0][2];
cross_2[0][2]=cross_ab3[0][0]*bc_un[0][1]-cross_ab3[0][1]*bc_un[0][0];

// create rotation matrix [vbc; cross_2; cross_ab3]

rot_M[0][0]=bc_un[0][0]; rot_M[0][1]=cross_2[0][0]; rot_M[0][2]=cross_ab3[0][0];
rot_M[1][0]=bc_un[0][1]; rot_M[1][1]=cross_2[0][1]; rot_M[1][2]=cross_ab3[0][1];
rot_M[2][0]=bc_un[0][2]; rot_M[2][1]=cross_2[0][2]; rot_M[2][2]=cross_ab3[0][2];


// calculate coord pre rotation row

pre_rot[0][0]= -p_length*cos(theta);
pre_rot[0][1]= p_length*sin(theta)*cos(phi);
pre_rot[0][2]= p_length*sin(theta)*sin(phi);

// displacement vector rot_M*pre_rot

disp[0][0]= rot_M[0][0]*pre_rot[0][0]+rot_M[0][1]*pre_rot[0][1]+rot_M[0][2]*pre_rot[0][2];
disp[0][1]= rot_M[1][0]*pre_rot[0][0]+rot_M[1][1]*pre_rot[0][1]+rot_M[1][2]*pre_rot[0][2];
disp[0][2]= rot_M[2][0]*pre_rot[0][0]+rot_M[2][1]*pre_rot[0][1]+rot_M[2][2]*pre_rot[0][2];

// Final coordinates of new atom

coordp[3][0]=coordp[2][0]+disp[0][0];
coordp[3][1]=coordp[2][1]+disp[0][1];
coordp[3][2]=coordp[2][2]+disp[0][2];




}
else{

phi=acos(Phi[2]);
theta=Angles[1]/360*2*PI;

// Unit vector between atom i-2 and i-3

ab_un[0][0]=0;
ab_un[0][1]=0;
ab_un[0][2]=1;

if (coordp[2][2] == Bounds[1]) { ab_un[0][1]=0.044710; ab_un[0][2]=0.999; }
if (coordp[2][2] == -Bounds[1]) { ab_un[0][1]=0.044710; ab_un[0][2]=0.999;}


// Unit vector between atom i-1 and i-2

mod_bc=pow(pow(coordp[2][0]-coordp[1][0],2)+pow(coordp[2][1]-coordp[1][1],2)+pow(coordp[2][2]-coordp[1][2],2),1./2.);

bc_un[0][0]=(coordp[2][0]-coordp[1][0])/mod_bc;
bc_un[0][1]=(coordp[2][1]-coordp[1][1])/mod_bc;
bc_un[0][2]=(coordp[2][2]-coordp[1][2])/mod_bc;

// Cross product of vector ab and vector bc

cross_ab3[0][0]=ab_un[0][1]*bc_un[0][2]-ab_un[0][2]*bc_un[0][1];
cross_ab3[0][1]=ab_un[0][2]*bc_un[0][0]-ab_un[0][0]*bc_un[0][2];
cross_ab3[0][2]=ab_un[0][0]*bc_un[0][1]-ab_un[0][1]*bc_un[0][0];

// Normalize the resulting vector

mod_cross=pow(pow(cross_ab3[0][0],2)+pow(cross_ab3[0][1],2)+pow(cross_ab3[0][2],2),1./2.);

cross_ab3[0][0]=cross_ab3[0][0]/mod_cross;
cross_ab3[0][1]=cross_ab3[0][1]/mod_cross;
cross_ab3[0][2]=cross_ab3[0][2]/mod_cross;


// Cross product between this vector and vector bc

cross_2[0][0]=cross_ab3[0][1]*bc_un[0][2]-cross_ab3[0][2]*bc_un[0][1];
cross_2[0][1]=cross_ab3[0][2]*bc_un[0][0]-cross_ab3[0][0]*bc_un[0][2];
cross_2[0][2]=cross_ab3[0][0]*bc_un[0][1]-cross_ab3[0][1]*bc_un[0][0];

// create rotation matrix [vbc; cross_2; cross_ab3]

rot_M[0][0]=bc_un[0][0]; rot_M[0][1]=cross_2[0][0]; rot_M[0][2]=cross_ab3[0][0];
rot_M[1][0]=bc_un[0][1]; rot_M[1][1]=cross_2[0][1]; rot_M[1][2]=cross_ab3[0][1];
rot_M[2][0]=bc_un[0][2]; rot_M[2][1]=cross_2[0][2]; rot_M[2][2]=cross_ab3[0][2];


// calculate coord pre rotation row

pre_rot[0][0]= -Bounds[2]*cos(theta);
pre_rot[0][1]= Bounds[2]*sin(theta)*cos(phi);
pre_rot[0][2]= Bounds[2]*sin(theta)*sin(phi);

// displacement vector rot_M*pre_rot

disp[0][0]= rot_M[0][0]*pre_rot[0][0]+rot_M[0][1]*pre_rot[0][1]+rot_M[0][2]*pre_rot[0][2];
disp[0][1]= rot_M[1][0]*pre_rot[0][0]+rot_M[1][1]*pre_rot[0][1]+rot_M[1][2]*pre_rot[0][2];
disp[0][2]= rot_M[2][0]*pre_rot[0][0]+rot_M[2][1]*pre_rot[0][1]+rot_M[2][2]*pre_rot[0][2];

// Final coordinates of new atom

coordp[3][0]=coordp[2][0]+disp[0][0];
coordp[3][1]=coordp[2][1]+disp[0][1];
coordp[3][2]=coordp[2][2]+disp[0][2];

}

for (i=4; i <= Natom; i++){
phi=Phi[i-1]/360*2*PI;

theta=Angles[i-2]/360*2*PI;

if (nstr == i-1){


// Unit vector between atom i-2 and i-3

mod_ab=pow(pow(coordp[i-2][0]-coordp[i-3][0],2)+pow(coordp[i-2][1]-coordp[i-3][1],2)+pow(coordp[i-2][2]-coordp[i-3][2],2),1./2.);

ab_un[0][0]=(coordp[i-2][0]-coordp[i-3][0])/mod_ab;
ab_un[0][1]=(coordp[i-2][1]-coordp[i-3][1])/mod_ab;
ab_un[0][2]=(coordp[i-2][2]-coordp[i-3][2])/mod_ab;

// Unit vector between atom i-1 and i-2

mod_bc=pow(pow(coordp[i-1][0]-coordp[i-2][0],2)+pow(coordp[i-1][1]-coordp[i-2][1],2)+pow(coordp[i-1][2]-coordp[i-2][2],2),1./2.);

bc_un[0][0]=(coordp[i-1][0]-coordp[i-2][0])/mod_bc;
bc_un[0][1]=(coordp[i-1][1]-coordp[i-2][1])/mod_bc;
bc_un[0][2]=(coordp[i-1][2]-coordp[i-2][2])/mod_bc;

// Cross product of vector ab and vector bc

cross_ab3[0][0]=ab_un[0][1]*bc_un[0][2]-ab_un[0][2]*bc_un[0][1];
cross_ab3[0][1]=ab_un[0][2]*bc_un[0][0]-ab_un[0][0]*bc_un[0][2];
cross_ab3[0][2]=ab_un[0][0]*bc_un[0][1]-ab_un[0][1]*bc_un[0][0];

// Normalize the resulting vector

mod_cross=pow(pow(cross_ab3[0][0],2)+pow(cross_ab3[0][1],2)+pow(cross_ab3[0][2],2),1./2.);

cross_ab3[0][0]=cross_ab3[0][0]/mod_cross;
cross_ab3[0][1]=cross_ab3[0][1]/mod_cross;
cross_ab3[0][2]=cross_ab3[0][2]/mod_cross;


// Cross product between this vector and vector bc

cross_2[0][0]=cross_ab3[0][1]*bc_un[0][2]-cross_ab3[0][2]*bc_un[0][1];
cross_2[0][1]=cross_ab3[0][2]*bc_un[0][0]-cross_ab3[0][0]*bc_un[0][2];
cross_2[0][2]=cross_ab3[0][0]*bc_un[0][1]-cross_ab3[0][1]*bc_un[0][0];

// create rotation matrix [vbc; cross_2; cross_ab3]

rot_M[0][0]=bc_un[0][0]; rot_M[0][1]=cross_2[0][0]; rot_M[0][2]=cross_ab3[0][0];
rot_M[1][0]=bc_un[0][1]; rot_M[1][1]=cross_2[0][1]; rot_M[1][2]=cross_ab3[0][1];
rot_M[2][0]=bc_un[0][2]; rot_M[2][1]=cross_2[0][2]; rot_M[2][2]=cross_ab3[0][2];


// calculate coord pre rotation row

pre_rot[0][0]= -p_length*cos(theta);
pre_rot[0][1]= p_length*sin(theta)*cos(phi);
pre_rot[0][2]= p_length*sin(theta)*sin(phi);

// displacement vector rot_M*pre_rot

disp[0][0]= rot_M[0][0]*pre_rot[0][0]+rot_M[0][1]*pre_rot[0][1]+rot_M[0][2]*pre_rot[0][2];
disp[0][1]= rot_M[1][0]*pre_rot[0][0]+rot_M[1][1]*pre_rot[0][1]+rot_M[1][2]*pre_rot[0][2];
disp[0][2]= rot_M[2][0]*pre_rot[0][0]+rot_M[2][1]*pre_rot[0][1]+rot_M[2][2]*pre_rot[0][2];

// Final coordinates of new atom

coordp[i][0]=coordp[i-1][0]+disp[0][0];
coordp[i][1]=coordp[i-1][1]+disp[0][1];
coordp[i][2]=coordp[i-1][2]+disp[0][2];



}



else{


// Unit vector between atom i-2 and i-3

mod_ab=pow(pow(coordp[i-2][0]-coordp[i-3][0],2)+pow(coordp[i-2][1]-coordp[i-3][1],2)+pow(coordp[i-2][2]-coordp[i-3][2],2),1./2.);

ab_un[0][0]=(coordp[i-2][0]-coordp[i-3][0])/mod_ab;
ab_un[0][1]=(coordp[i-2][1]-coordp[i-3][1])/mod_ab;
ab_un[0][2]=(coordp[i-2][2]-coordp[i-3][2])/mod_ab;

// Unit vector between atom i-1 and i-2

mod_bc=pow(pow(coordp[i-1][0]-coordp[i-2][0],2)+pow(coordp[i-1][1]-coordp[i-2][1],2)+pow(coordp[i-1][2]-coordp[i-2][2],2),1./2.);

bc_un[0][0]=(coordp[i-1][0]-coordp[i-2][0])/mod_bc;
bc_un[0][1]=(coordp[i-1][1]-coordp[i-2][1])/mod_bc;
bc_un[0][2]=(coordp[i-1][2]-coordp[i-2][2])/mod_bc;

// Cross product of vector ab and vector bc

cross_ab3[0][0]=ab_un[0][1]*bc_un[0][2]-ab_un[0][2]*bc_un[0][1];
cross_ab3[0][1]=ab_un[0][2]*bc_un[0][0]-ab_un[0][0]*bc_un[0][2];
cross_ab3[0][2]=ab_un[0][0]*bc_un[0][1]-ab_un[0][1]*bc_un[0][0];

// Normalize the resulting vector

mod_cross=pow(pow(cross_ab3[0][0],2)+pow(cross_ab3[0][1],2)+pow(cross_ab3[0][2],2),1./2.);

cross_ab3[0][0]=cross_ab3[0][0]/mod_cross;
cross_ab3[0][1]=cross_ab3[0][1]/mod_cross;
cross_ab3[0][2]=cross_ab3[0][2]/mod_cross;


// Cross product between this vector and vector bc

cross_2[0][0]=cross_ab3[0][1]*bc_un[0][2]-cross_ab3[0][2]*bc_un[0][1];
cross_2[0][1]=cross_ab3[0][2]*bc_un[0][0]-cross_ab3[0][0]*bc_un[0][2];
cross_2[0][2]=cross_ab3[0][0]*bc_un[0][1]-cross_ab3[0][1]*bc_un[0][0];

// create rotation matrix [vbc; cross_2; cross_ab3]

rot_M[0][0]=bc_un[0][0]; rot_M[0][1]=cross_2[0][0]; rot_M[0][2]=cross_ab3[0][0];
rot_M[1][0]=bc_un[0][1]; rot_M[1][1]=cross_2[0][1]; rot_M[1][2]=cross_ab3[0][1];
rot_M[2][0]=bc_un[0][2]; rot_M[2][1]=cross_2[0][2]; rot_M[2][2]=cross_ab3[0][2];


// calculate coord pre rotation row

pre_rot[0][0]= -Bounds[i-1]*cos(theta);
pre_rot[0][1]= Bounds[i-1]*sin(theta)*cos(phi);
pre_rot[0][2]= Bounds[i-1]*sin(theta)*sin(phi);

// displacement vector rot_M*pre_rot

disp[0][0]= rot_M[0][0]*pre_rot[0][0]+rot_M[0][1]*pre_rot[0][1]+rot_M[0][2]*pre_rot[0][2];
disp[0][1]= rot_M[1][0]*pre_rot[0][0]+rot_M[1][1]*pre_rot[0][1]+rot_M[1][2]*pre_rot[0][2];
disp[0][2]= rot_M[2][0]*pre_rot[0][0]+rot_M[2][1]*pre_rot[0][1]+rot_M[2][2]*pre_rot[0][2];

// Final coordinates of new atom

coordp[i][0]=coordp[i-1][0]+disp[0][0];
coordp[i][1]=coordp[i-1][1]+disp[0][1];
coordp[i][2]=coordp[i-1][2]+disp[0][2];



}
}
}

void gencoordpbnd(int nbnd)
{
double ab_un[1][3], bc_un[1][3], mod_ab, mod_bc, cross_ab3[1][3], mod_cross, cross_2[1][3], rot_M[3][3], pre_rot[1][3], disp[1][3];
int i,signe;
double alea, theta, phi;
double   dang=0.5;

// adapted from python subroutine 
// https://github.com/rbhatta8/protein-design/blob/master/nerf.py


// Provisional angle of the bended angle


alea=ran2(&idum);
if (alea <= 0.5) {signe=1;} else {signe=-1;}

p_angle= Angles[nbnd]+dang*signe;

coordp[1][0]=0;
coordp[1][1]=0;
coordp[1][2]=0;

coordp[2][0]=Bounds[1]*sin(acos(Phi[1]))*cos(Angles[0]/360*2*PI);
coordp[2][1]=Bounds[1]*sin(acos(Phi[1]))*sin(Angles[0]/360*2*PI);
coordp[2][2]=Bounds[1]*cos(acos(Phi[1]));





if (nbnd == 1) {
  

theta=p_angle/360*2*PI;   
phi=acos(Phi[2]);

ab_un[0][0]=0;
ab_un[0][1]=0;
ab_un[0][2]=1;

if (coordp[2][2] == Bounds[1]) { ab_un[0][1]=0.044710; ab_un[0][2]=0.999;}
if (coordp[2][2] == -Bounds[1]) { ab_un[0][1]=0.044710; ab_un[0][2]=0.999;}


// Unit vector between atom i-1 and i-2

mod_bc=pow(pow(coordp[2][0]-coordp[1][0],2)+pow(coordp[2][1]-coordp[1][1],2)+pow(coordp[2][2]-coordp[1][2],2),1./2.);

bc_un[0][0]=(coordp[2][0]-coordp[1][0])/mod_bc;
bc_un[0][1]=(coordp[2][1]-coordp[1][1])/mod_bc;
bc_un[0][2]=(coordp[2][2]-coordp[1][2])/mod_bc;

// Cross product of vector ab and vector bc

cross_ab3[0][0]=ab_un[0][1]*bc_un[0][2]-ab_un[0][2]*bc_un[0][1];
cross_ab3[0][1]=ab_un[0][2]*bc_un[0][0]-ab_un[0][0]*bc_un[0][2];
cross_ab3[0][2]=ab_un[0][0]*bc_un[0][1]-ab_un[0][1]*bc_un[0][0];

// Normalize the resulting vector

mod_cross=pow(pow(cross_ab3[0][0],2)+pow(cross_ab3[0][1],2)+pow(cross_ab3[0][2],2),1./2.);

cross_ab3[0][0]=cross_ab3[0][0]/mod_cross;
cross_ab3[0][1]=cross_ab3[0][1]/mod_cross;
cross_ab3[0][2]=cross_ab3[0][2]/mod_cross;


// Cross product between this vector and vector bc

cross_2[0][0]=cross_ab3[0][1]*bc_un[0][2]-cross_ab3[0][2]*bc_un[0][1];
cross_2[0][1]=cross_ab3[0][2]*bc_un[0][0]-cross_ab3[0][0]*bc_un[0][2];
cross_2[0][2]=cross_ab3[0][0]*bc_un[0][1]-cross_ab3[0][1]*bc_un[0][0];

// create rotation matrix [vbc; cross_2; cross_ab3]

rot_M[0][0]=bc_un[0][0]; rot_M[0][1]=cross_2[0][0]; rot_M[0][2]=cross_ab3[0][0];
rot_M[1][0]=bc_un[0][1]; rot_M[1][1]=cross_2[0][1]; rot_M[1][2]=cross_ab3[0][1];
rot_M[2][0]=bc_un[0][2]; rot_M[2][1]=cross_2[0][2]; rot_M[2][2]=cross_ab3[0][2];


// calculate coord pre rotation row

pre_rot[0][0]= -Bounds[2]*cos(theta);
pre_rot[0][1]= Bounds[2]*sin(theta)*cos(phi);
pre_rot[0][2]= Bounds[2]*sin(theta)*sin(phi);

// displacement vector rot_M*pre_rot

disp[0][0]= rot_M[0][0]*pre_rot[0][0]+rot_M[0][1]*pre_rot[0][1]+rot_M[0][2]*pre_rot[0][2];
disp[0][1]= rot_M[1][0]*pre_rot[0][0]+rot_M[1][1]*pre_rot[0][1]+rot_M[1][2]*pre_rot[0][2];
disp[0][2]= rot_M[2][0]*pre_rot[0][0]+rot_M[2][1]*pre_rot[0][1]+rot_M[2][2]*pre_rot[0][2];

// Final coordinates of new atom

coordp[3][0]=coordp[2][0]+disp[0][0];
coordp[3][1]=coordp[2][1]+disp[0][1];
coordp[3][2]=coordp[2][2]+disp[0][2];




}
else{

phi=acos(Phi[2]);
theta=Angles[1]/360*2*PI;


// Unit vector between atom i-2 and i-3


ab_un[0][0]=0;
ab_un[0][1]=0;
ab_un[0][2]=1;

if (coordp[2][2] == Bounds[1]) { ab_un[0][1]=0.044710; ab_un[0][2]=0.999;}
if (coordp[2][2] == -Bounds[1]) { ab_un[0][1]=0.044710; ab_un[0][2]=0.999;}


// Unit vector between atom i-1 and i-2

mod_bc=pow(pow(coordp[2][0]-coordp[1][0],2)+pow(coordp[2][1]-coordp[1][1],2)+pow(coordp[2][2]-coordp[1][2],2),1./2.);

bc_un[0][0]=(coordp[2][0]-coordp[1][0])/mod_bc;
bc_un[0][1]=(coordp[2][1]-coordp[1][1])/mod_bc;
bc_un[0][2]=(coordp[2][2]-coordp[1][2])/mod_bc;

// Cross product of vector ab and vector bc

cross_ab3[0][0]=ab_un[0][1]*bc_un[0][2]-ab_un[0][2]*bc_un[0][1];
cross_ab3[0][1]=ab_un[0][2]*bc_un[0][0]-ab_un[0][0]*bc_un[0][2];
cross_ab3[0][2]=ab_un[0][0]*bc_un[0][1]-ab_un[0][1]*bc_un[0][0];

// Normalize the resulting vector

mod_cross=pow(pow(cross_ab3[0][0],2)+pow(cross_ab3[0][1],2)+pow(cross_ab3[0][2],2),1./2.);

cross_ab3[0][0]=cross_ab3[0][0]/mod_cross;
cross_ab3[0][1]=cross_ab3[0][1]/mod_cross;
cross_ab3[0][2]=cross_ab3[0][2]/mod_cross;


// Cross product between this vector and vector bc

cross_2[0][0]=cross_ab3[0][1]*bc_un[0][2]-cross_ab3[0][2]*bc_un[0][1];
cross_2[0][1]=cross_ab3[0][2]*bc_un[0][0]-cross_ab3[0][0]*bc_un[0][2];
cross_2[0][2]=cross_ab3[0][0]*bc_un[0][1]-cross_ab3[0][1]*bc_un[0][0];

// create rotation matrix [vbc; cross_2; cross_ab3]

rot_M[0][0]=bc_un[0][0]; rot_M[0][1]=cross_2[0][0]; rot_M[0][2]=cross_ab3[0][0];
rot_M[1][0]=bc_un[0][1]; rot_M[1][1]=cross_2[0][1]; rot_M[1][2]=cross_ab3[0][1];
rot_M[2][0]=bc_un[0][2]; rot_M[2][1]=cross_2[0][2]; rot_M[2][2]=cross_ab3[0][2];


// calculate coord pre rotation row

pre_rot[0][0]= -Bounds[2]*cos(theta);
pre_rot[0][1]= Bounds[2]*sin(theta)*cos(phi);
pre_rot[0][2]= Bounds[2]*sin(theta)*sin(phi);

// displacement vector rot_M*pre_rot

disp[0][0]= rot_M[0][0]*pre_rot[0][0]+rot_M[0][1]*pre_rot[0][1]+rot_M[0][2]*pre_rot[0][2];
disp[0][1]= rot_M[1][0]*pre_rot[0][0]+rot_M[1][1]*pre_rot[0][1]+rot_M[1][2]*pre_rot[0][2];
disp[0][2]= rot_M[2][0]*pre_rot[0][0]+rot_M[2][1]*pre_rot[0][1]+rot_M[2][2]*pre_rot[0][2];

// Final coordinates of new atom

coordp[3][0]=coordp[2][0]+disp[0][0];
coordp[3][1]=coordp[2][1]+disp[0][1];
coordp[3][2]=coordp[2][2]+disp[0][2];

}

for (i=4; i <= Natom; i++){
phi=Phi[i-1]/360*2*PI;

if (nbnd == i-2) {
theta=p_angle/360*2*PI;

// Unit vector between atom i-2 and i-3

mod_ab=pow(pow(coordp[i-2][0]-coordp[i-3][0],2)+pow(coordp[i-2][1]-coordp[i-3][1],2)+pow(coordp[i-2][2]-coordp[i-3][2],2),1./2.);

ab_un[0][0]=(coordp[i-2][0]-coordp[i-3][0])/mod_ab;
ab_un[0][1]=(coordp[i-2][1]-coordp[i-3][1])/mod_ab;
ab_un[0][2]=(coordp[i-2][2]-coordp[i-3][2])/mod_ab;

// Unit vector between atom i-1 and i-2

mod_bc=pow(pow(coordp[i-1][0]-coordp[i-2][0],2)+pow(coordp[i-1][1]-coordp[i-2][1],2)+pow(coordp[i-1][2]-coordp[i-2][2],2),1./2.);

bc_un[0][0]=(coordp[i-1][0]-coordp[i-2][0])/mod_bc;
bc_un[0][1]=(coordp[i-1][1]-coordp[i-2][1])/mod_bc;
bc_un[0][2]=(coordp[i-1][2]-coordp[i-2][2])/mod_bc;

// Cross product of vector ab and vector bc

cross_ab3[0][0]=ab_un[0][1]*bc_un[0][2]-ab_un[0][2]*bc_un[0][1];
cross_ab3[0][1]=ab_un[0][2]*bc_un[0][0]-ab_un[0][0]*bc_un[0][2];
cross_ab3[0][2]=ab_un[0][0]*bc_un[0][1]-ab_un[0][1]*bc_un[0][0];

// Normalize the resulting vector

mod_cross=pow(pow(cross_ab3[0][0],2)+pow(cross_ab3[0][1],2)+pow(cross_ab3[0][2],2),1./2.);

cross_ab3[0][0]=cross_ab3[0][0]/mod_cross;
cross_ab3[0][1]=cross_ab3[0][1]/mod_cross;
cross_ab3[0][2]=cross_ab3[0][2]/mod_cross;


// Cross product between this vector and vector bc

cross_2[0][0]=cross_ab3[0][1]*bc_un[0][2]-cross_ab3[0][2]*bc_un[0][1];
cross_2[0][1]=cross_ab3[0][2]*bc_un[0][0]-cross_ab3[0][0]*bc_un[0][2];
cross_2[0][2]=cross_ab3[0][0]*bc_un[0][1]-cross_ab3[0][1]*bc_un[0][0];

// create rotation matrix [vbc; cross_2; cross_ab3]

rot_M[0][0]=bc_un[0][0]; rot_M[0][1]=cross_2[0][0]; rot_M[0][2]=cross_ab3[0][0];
rot_M[1][0]=bc_un[0][1]; rot_M[1][1]=cross_2[0][1]; rot_M[1][2]=cross_ab3[0][1];
rot_M[2][0]=bc_un[0][2]; rot_M[2][1]=cross_2[0][2]; rot_M[2][2]=cross_ab3[0][2];


// calculate coord pre rotation row

pre_rot[0][0]= -Bounds[i-1]*cos(theta);
pre_rot[0][1]= Bounds[i-1]*sin(theta)*cos(phi);
pre_rot[0][2]= Bounds[i-1]*sin(theta)*sin(phi);

// displacement vector rot_M*pre_rot

disp[0][0]= rot_M[0][0]*pre_rot[0][0]+rot_M[0][1]*pre_rot[0][1]+rot_M[0][2]*pre_rot[0][2];
disp[0][1]= rot_M[1][0]*pre_rot[0][0]+rot_M[1][1]*pre_rot[0][1]+rot_M[1][2]*pre_rot[0][2];
disp[0][2]= rot_M[2][0]*pre_rot[0][0]+rot_M[2][1]*pre_rot[0][1]+rot_M[2][2]*pre_rot[0][2];

// Final coordinates of new atom

coordp[i][0]=coordp[i-1][0]+disp[0][0];
coordp[i][1]=coordp[i-1][1]+disp[0][1];
coordp[i][2]=coordp[i-1][2]+disp[0][2];



}



else{

theta=Angles[i-2]/360*2*PI;

// Unit vector between atom i-2 and i-3

mod_ab=pow(pow(coordp[i-2][0]-coordp[i-3][0],2)+pow(coordp[i-2][1]-coordp[i-3][1],2)+pow(coordp[i-2][2]-coordp[i-3][2],2),1./2.);

ab_un[0][0]=(coordp[i-2][0]-coordp[i-3][0])/mod_ab;
ab_un[0][1]=(coordp[i-2][1]-coordp[i-3][1])/mod_ab;
ab_un[0][2]=(coordp[i-2][2]-coordp[i-3][2])/mod_ab;

// Unit vector between atom i-1 and i-2

mod_bc=pow(pow(coordp[i-1][0]-coordp[i-2][0],2)+pow(coordp[i-1][1]-coordp[i-2][1],2)+pow(coordp[i-1][2]-coordp[i-2][2],2),1./2.);

bc_un[0][0]=(coordp[i-1][0]-coordp[i-2][0])/mod_bc;
bc_un[0][1]=(coordp[i-1][1]-coordp[i-2][1])/mod_bc;
bc_un[0][2]=(coordp[i-1][2]-coordp[i-2][2])/mod_bc;

// Cross product of vector ab and vector bc

cross_ab3[0][0]=ab_un[0][1]*bc_un[0][2]-ab_un[0][2]*bc_un[0][1];
cross_ab3[0][1]=ab_un[0][2]*bc_un[0][0]-ab_un[0][0]*bc_un[0][2];
cross_ab3[0][2]=ab_un[0][0]*bc_un[0][1]-ab_un[0][1]*bc_un[0][0];

// Normalize the resulting vector

mod_cross=pow(pow(cross_ab3[0][0],2)+pow(cross_ab3[0][1],2)+pow(cross_ab3[0][2],2),1./2.);

cross_ab3[0][0]=cross_ab3[0][0]/mod_cross;
cross_ab3[0][1]=cross_ab3[0][1]/mod_cross;
cross_ab3[0][2]=cross_ab3[0][2]/mod_cross;


// Cross product between this vector and vector bc

cross_2[0][0]=cross_ab3[0][1]*bc_un[0][2]-cross_ab3[0][2]*bc_un[0][1];
cross_2[0][1]=cross_ab3[0][2]*bc_un[0][0]-cross_ab3[0][0]*bc_un[0][2];
cross_2[0][2]=cross_ab3[0][0]*bc_un[0][1]-cross_ab3[0][1]*bc_un[0][0];

// create rotation matrix [vbc; cross_2; cross_ab3]

rot_M[0][0]=bc_un[0][0]; rot_M[0][1]=cross_2[0][0]; rot_M[0][2]=cross_ab3[0][0];
rot_M[1][0]=bc_un[0][1]; rot_M[1][1]=cross_2[0][1]; rot_M[1][2]=cross_ab3[0][1];
rot_M[2][0]=bc_un[0][2]; rot_M[2][1]=cross_2[0][2]; rot_M[2][2]=cross_ab3[0][2];


// calculate coord pre rotation row

pre_rot[0][0]= -Bounds[i-1]*cos(theta);
pre_rot[0][1]= Bounds[i-1]*sin(theta)*cos(phi);
pre_rot[0][2]= Bounds[i-1]*sin(theta)*sin(phi);

// displacement vector rot_M*pre_rot

disp[0][0]= rot_M[0][0]*pre_rot[0][0]+rot_M[0][1]*pre_rot[0][1]+rot_M[0][2]*pre_rot[0][2];
disp[0][1]= rot_M[1][0]*pre_rot[0][0]+rot_M[1][1]*pre_rot[0][1]+rot_M[1][2]*pre_rot[0][2];
disp[0][2]= rot_M[2][0]*pre_rot[0][0]+rot_M[2][1]*pre_rot[0][1]+rot_M[2][2]*pre_rot[0][2];

// Final coordinates of new atom

coordp[i][0]=coordp[i-1][0]+disp[0][0];
coordp[i][1]=coordp[i-1][1]+disp[0][1];
coordp[i][2]=coordp[i-1][2]+disp[0][2];



}
}
}
