// Module with all the subroutines needed to inicialize the MC_simulation

#include "lib.h"
#include "global.h"

/*******************************************************/
void fes_polymer_trans(int Natom)
{
double ab_un[1][3], bc_un[1][3], mod_ab, mod_bc, cross_ab3[1][3], mod_cross, cross_2[1][3], rot_M[3][3], pre_rot[1][3], disp[1][3];
int i, nat;
double theta, phi;
long idum_aux; 
MPI_Request request;
time_t rawtime;
struct tm * timeinfo;

//Subroutine to generate the initial conformation of the polymer. The conformation is settled to have all the chain bounds in trans conformation. 


/* First atom */

conformacioatoms[1]=TRANS;
coord[1][0]=0;
coord[1][1]=0;
coord[1][2]=0;


/* Second atom */

conformacioatoms[2]=TRANS;
Phi[1]=0.0;
Angles[0]=0;
coord[2][0]=Bounds[1];
coord[2][1]=0;
coord[2][2]=0;

/* Third atom */

conformacioatoms[3]=TRANS;
Phi[2]=0.0;
coord[3][0]=coord[2][0]+Bounds[2]*cos((180-Angles[1])*PI/180.);
coord[3][1]=coord[2][1]+Bounds[2]*sin((180-Angles[1])*PI/180.);
coord[3][2]=0;

/* The remaining atoms */

for (i=4; i <= Natom; i++){

   conformacioatoms[i]=GAUCHEP;
   Phi[i-1]=120.0;
   theta=Angles[i-2]/360*2*PI;
   phi=Phi[i-1]/360*2*PI;
// Unit vector between atom i-2 and i-3

mod_ab=pow(pow(coord[i-2][0]-coord[i-3][0],2)+pow(coord[i-2][1]-coord[i-3][1],2)+pow(coord[i-2][2]-coord[i-3][2],2),1./2.);

ab_un[0][0]=(coord[i-2][0]-coord[i-3][0])/mod_ab;
ab_un[0][1]=(coord[i-2][1]-coord[i-3][1])/mod_ab;
ab_un[0][2]=(coord[i-2][2]-coord[i-3][2])/mod_ab;

// Unit vector between atom i-3 and i-2

mod_bc=pow(pow(coord[i-1][0]-coord[i-2][0],2)+pow(coord[i-1][1]-coord[i-2][1],2)+pow(coord[i-1][2]-coord[i-2][2],2),1./2.);

bc_un[0][0]=(coord[i-1][0]-coord[i-2][0])/mod_bc;
bc_un[0][1]=(coord[i-1][1]-coord[i-2][1])/mod_bc;
bc_un[0][2]=(coord[i-1][2]-coord[i-2][2])/mod_bc;

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

coord[i][0]=coord[i-1][0]+disp[0][0];
coord[i][1]=coord[i-1][1]+disp[0][1];
coord[i][2]=coord[i-1][2]+disp[0][2];


 }

/*  All the charges are set to the reference state. This is the state with lower energy */

for(nat=1;nat<= Natom; nat++){
   
        carrega[nat]=L_charge[nat];

        Ebond_prot[nat]=1; 
}

// If the random seed is set to 0, the program reads the system hour

if(idum_ref==0){

time ( &rawtime ); idum=-abs(rawtime);


idum=idum+nrank*100;
  cout << "cpu= " << nrank << " idum= " << idum << endl;

}
else{idum=idum_ref;}



}

void init_Emat(){
// Subroutine tu inicialice the conformational energy matrix



// Energy Matrix 1

Econf1[1][1]=-Rkcal*temp*log(Conf1[1][1]); Econf1[1][2]=-Rkcal*temp*log(Conf1[1][2]); Econf1[1][3]=-Rkcal*temp*log(Conf1[1][3]);
Econf1[2][1]=-Rkcal*temp*log(Conf1[2][1]); Econf1[2][2]=-Rkcal*temp*log(Conf1[2][2]); Econf1[2][3]=-Rkcal*temp*log(Conf1[2][3]);
Econf1[3][1]=-Rkcal*temp*log(Conf1[3][1]); Econf1[3][2]=-Rkcal*temp*log(Conf1[3][2]); Econf1[3][3]=-Rkcal*temp*log(Conf1[3][3]);


// Energy Matrix 2

Econf2[1][1]=-Rkcal*temp*log(Conf2[1][1]); Econf2[1][2]=-Rkcal*temp*log(Conf2[1][2]); Econf2[1][3]=-Rkcal*temp*log(Conf2[1][3]);
Econf2[2][1]=-Rkcal*temp*log(Conf2[2][1]); Econf2[2][2]=-Rkcal*temp*log(Conf2[2][2]); Econf2[2][3]=-Rkcal*temp*log(Conf2[2][3]);
Econf2[3][1]=-Rkcal*temp*log(Conf2[3][1]); Econf2[3][2]=-Rkcal*temp*log(Conf2[3][2]); Econf2[3][3]=-Rkcal*temp*log(Conf2[3][3]);



// Energy Matrix 3

Econf3[1][1]=-Rkcal*temp*log(Conf3[1][1]); Econf3[1][2]=-Rkcal*temp*log(Conf3[1][2]); Econf3[1][3]=-Rkcal*temp*log(Conf3[1][3]);
Econf3[2][1]=-Rkcal*temp*log(Conf3[2][1]); Econf3[2][2]=-Rkcal*temp*log(Conf3[2][2]); Econf3[2][3]=-Rkcal*temp*log(Conf3[2][3]);
Econf3[3][1]=-Rkcal*temp*log(Conf3[3][1]); Econf3[3][2]=-Rkcal*temp*log(Conf3[3][2]); Econf3[3][3]=-Rkcal*temp*log(Conf3[3][3]);



// Energy Matrix 4

Econf4[1][1]=-Rkcal*temp*log(Conf4[1][1]); Econf4[1][2]=-Rkcal*temp*log(Conf4[1][2]); Econf4[1][3]=-Rkcal*temp*log(Conf4[1][3]);
Econf4[2][1]=-Rkcal*temp*log(Conf4[2][1]); Econf4[2][2]=-Rkcal*temp*log(Conf4[2][2]); Econf4[2][3]=-Rkcal*temp*log(Conf4[2][3]);
Econf4[3][1]=-Rkcal*temp*log(Conf4[3][1]); Econf4[3][2]=-Rkcal*temp*log(Conf4[3][2]); Econf4[3][3]=-Rkcal*temp*log(Conf4[3][3]);



// Energy Matrix 5

Econf5[1][1]=-Rkcal*temp*log(Conf5[1][1]); Econf5[1][2]=-Rkcal*temp*log(Conf5[1][2]); Econf5[1][3]=-Rkcal*temp*log(Conf5[1][3]);
Econf5[2][1]=-Rkcal*temp*log(Conf5[2][1]); Econf5[2][2]=-Rkcal*temp*log(Conf5[2][2]); Econf5[2][3]=-Rkcal*temp*log(Conf5[2][3]);
Econf5[3][1]=-Rkcal*temp*log(Conf5[3][1]); Econf5[3][2]=-Rkcal*temp*log(Conf5[3][2]); Econf5[3][3]=-Rkcal*temp*log(Conf5[3][3]);



// Energy Matrix 6

Econf6[1][1]=-Rkcal*temp*log(Conf6[1][1]); Econf6[1][2]=-Rkcal*temp*log(Conf6[1][2]); Econf6[1][3]=-Rkcal*temp*log(Conf6[1][3]);
Econf6[2][1]=-Rkcal*temp*log(Conf6[2][1]); Econf6[2][2]=-Rkcal*temp*log(Conf6[2][2]); Econf6[2][3]=-Rkcal*temp*log(Conf6[2][3]);
Econf6[3][1]=-Rkcal*temp*log(Conf6[3][1]); Econf6[3][2]=-Rkcal*temp*log(Conf6[3][2]); Econf6[3][3]=-Rkcal*temp*log(Conf6[3][3]);


// Energy Matrix 7

Econf7[1][1]=-Rkcal*temp*log(Conf7[1][1]); Econf7[1][2]=-Rkcal*temp*log(Conf7[1][2]); Econf7[1][3]=-Rkcal*temp*log(Conf7[1][3]);
Econf7[2][1]=-Rkcal*temp*log(Conf7[2][1]); Econf7[2][2]=-Rkcal*temp*log(Conf7[2][2]); Econf7[2][3]=-Rkcal*temp*log(Conf7[2][3]);
Econf7[3][1]=-Rkcal*temp*log(Conf7[3][1]); Econf7[3][2]=-Rkcal*temp*log(Conf7[3][2]); Econf7[3][3]=-Rkcal*temp*log(Conf7[3][3]);

// Energy Matrix 8

Econf8[1][1]=-Rkcal*temp*log(Conf8[1][1]); Econf8[1][2]=-Rkcal*temp*log(Conf8[1][2]); Econf8[1][3]=-Rkcal*temp*log(Conf8[1][3]);
Econf8[2][1]=-Rkcal*temp*log(Conf8[2][1]); Econf8[2][2]=-Rkcal*temp*log(Conf8[2][2]); Econf8[2][3]=-Rkcal*temp*log(Conf8[2][3]);
Econf8[3][1]=-Rkcal*temp*log(Conf8[3][1]); Econf8[3][2]=-Rkcal*temp*log(Conf8[3][2]); Econf8[3][3]=-Rkcal*temp*log(Conf8[3][3]);

// Energy Matrix 9

Econf9[1][1]=-Rkcal*temp*log(Conf9[1][1]); Econf9[1][2]=-Rkcal*temp*log(Conf9[1][2]); Econf9[1][3]=-Rkcal*temp*log(Conf9[1][3]);
Econf9[2][1]=-Rkcal*temp*log(Conf9[2][1]); Econf9[2][2]=-Rkcal*temp*log(Conf9[2][2]); Econf9[2][3]=-Rkcal*temp*log(Conf9[2][3]);
Econf9[3][1]=-Rkcal*temp*log(Conf9[3][1]); Econf9[3][2]=-Rkcal*temp*log(Conf9[3][2]); Econf9[3][3]=-Rkcal*temp*log(Conf9[3][3]);

// Energy Matrix 10

Econf10[1][1]=-Rkcal*temp*log(Conf10[1][1]); Econf10[1][2]=-Rkcal*temp*log(Conf10[1][2]); Econf10[1][3]=-Rkcal*temp*log(Conf10[1][3]);
Econf10[2][1]=-Rkcal*temp*log(Conf10[2][1]); Econf10[2][2]=-Rkcal*temp*log(Conf10[2][2]); Econf10[2][3]=-Rkcal*temp*log(Conf10[2][3]);
Econf10[3][1]=-Rkcal*temp*log(Conf10[3][1]); Econf10[3][2]=-Rkcal*temp*log(Conf10[3][2]); Econf10[3][3]=-Rkcal*temp*log(Conf10[3][3]);





}



void init_energy(){
int nat,nbonds,i,j;
double dx,dy,dz,ener_bound,  prom_angle,ener_prot_aux;
double Aener,Cener,dist2,dist6,dist12,Aenertot,Cenertot;
double enerJoules,enJ,dist;
int nat1,nat2;
char C[2]="C", N[2]="N";

// Subroutine to calculate the initial energy of the system



init_Emat(); // set the conformational energy matrix

/* *** Initial energy calculation *** */


// The first and last bound do not have any conformational energy

Ebond[1]=0; 
Ebond[Natom]=0; 
Ebond[Natom-1]=0;
ener_prot=0;
ener_conf=0;
EDebye=0;
W_mech=0;
// Conformational initial energy


for(nat=Ntranex;nat<=Natom-Ntranex;nat++)
 {get_Ebound(nat,&ener_bound);
   Ebond[nat]=ener_bound;
   ener_conf= ener_conf+ener_bound;
 }



// Protonation initial energy

  ener_prot_aux=get_Eprot(confMC);
  ener_prot=ener_prot+ener_prot_aux;

// r2 matrix calculation

for(i=1;i< Natom;i++)
 for(j=i+1;j<= Natom;j++)
   {
      dx=coord[i][0]-coord[j][0];
      dy=coord[i][1]-coord[j][1];
      dz=coord[i][2]-coord[j][2];
      matriu_dist2[i][j]=dx*dx+dy*dy+dz*dz;
   }


// energy of the long range interaction (Debye)

if(Debye==0)
{
        // Ener=q*q/(4*pi*epsi_0) * (exp(-kappa*r) / r

        // Kappa calculation

        kappa=sqrt((2*NA*q_elect*q_elect*I)/(epsi0*DIEBULK*kBoltzmann*temp));

        /* kappa units: 10**(-17/2) A-2 */

        kappa=kappa*pow(10,-17./2);

        /* kappa units: A-1 */

        // Debye prefactor calculation
 
        prefactDebye=1e10*q_elect*q_elect/(4*PI*epsi0*DIEBULK);  /* Units: J*A  (1 particule) */

// 12DCE: 331.842*qso*qh : kcal/mol
//      R= 0.00198589
//      RT=R*T


        enerJoules=0;

        for(nat1=1;nat1<=Natom;nat1++)
         {
          for(nat2=nat1+1;nat2<=Natom;nat2++)
            {
            if( abs(nat1-nat2)>1  )
               {
                if(carrega[nat1] != 0 && carrega[nat2] != 0)
                  {
                   dist2=matriu_dist2[nat1][nat2];
                   dist=sqrt(dist2);

                   enJ=(carrega[nat1]*carrega[nat2]*prefactDebye/dist)*(exp(-kappa*dist));
                   matriu_EDebye[nat1][nat2]=enJ;
                   enerJoules=enerJoules+enJ; 
                  }
               }
             }
          }

EDebye=EDebye+enerJoules/JoulesCal/1000*NA;

}


if ( F != 0 ) {

W_mech= get_Wmech();

}


// TOTAL ENERGY OF THE SYSTEM

etotkcal=ener_conf+ener_prot+EDebye+W_mech;




      movsi=0;
      movno=0;
    for(nat=1;nat<=Natom;nat++)
     {
      vectormovsi[nat]=0;
      vectormovno[nat]=0;
     }


/* Accumulators inicializations */
acc_r2_100=0; acc_s2_100=0; acc_qt_100=0; acc_Lx_100=0, acc_cov_100=0, acc_N_CC_trans_100=0, acc_N_NC_trans_100=0, acc_N_trans_100=0;
acc_r2_90=0; acc_s2_90=0; acc_qt_90=0;  acc_Lx_90=0, acc_cov_90=0;
acc_r2_80=0; acc_s2_80=0; acc_qt_80=0;   acc_Lx_80=0, acc_cov_80=0;
acc_r2_50=0; acc_s2_50=0; acc_qt_50=0;    acc_Lx_50=0, acc_cov_50=0;
n_acc_100=0; n_acc_90=0; n_acc_80=0; n_acc_50=0;
EDebye_acc=0;
EincprotDeb_acc=0;

nbounds=Natom-1;

// Aproximation of the  r2 and s2 (Flory)

prom_bound=0;
prom_angle=0;

for(i=1; i <= Natom-1 ; i++){prom_bound=prom_bound+Bounds[i];}

prom_bound=prom_bound/(Natom-1);

for(i=1; i <= Natom-2 ; i++){prom_angle=prom_angle+Angles[i];}

prom_angle=prom_angle/(Natom-2);



r2endToend=r2_FreelyRotatingChain(nbounds,prom_bound,prom_angle);
s2_free=s2_FreelyRotatingChain(nbounds,prom_bound,prom_angle);


// Equilibrium state of the bonds

for(i=1; i <= Natom-1 ; i++){Bounds_eq[i]=Bounds[i];}

// Equilibrium state of the angles

for(i=1; i <= Natom-2 ; i++){Angles_eq[i]=Angles[i];}

// Inicialization of the stretching control variable

control_str=0;

// Inicialization of the bending control variable

control_bnd=0;

// Computation of the total number of protonated sites 

N_prot=0;
N_sites=0;

for(i=1; i <= Natom; i++){
if (pKa[i] != -100)
{N_sites=N_sites+1;}}


// Inicialization of % of acceptance

N_test_str=0;
N_test_rot=0;
N_accepted_str=0;
N_accepted_rot=0;

// Inicialization of trans counters

N_CC=0;
N_NC=0;

for(i=2; i <= Natom-2 ; i++){

if ( chain[i] == C[0] && chain[i+1] == C[0]){N_CC++;}
if ( chain[i] == N[0] && chain[i+1] == C[0]){N_NC++;}
if ( chain[i] == C[0] && chain[i+1] == N[0]){N_NC++;}



}


}
