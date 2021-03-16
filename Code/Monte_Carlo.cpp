// Module with all the subroutines concerning the MonteCarlo algorithm

#include "lib.h"
#include "global.h"

/*******************************************************/


void generaconf(void){

float alea;
double angle_Rot_Conformacio;
int signe;
double Probrot=0.1;

if (nconf < 100 ) {

nmog=1;

}
else{
alea=ran2(&idum);

if (alea < Probrot) {
alea=ran2(&idum);
if (alea < 0.5) {nmog=1; } else {nmog =2;}
if(nmog == 1) {  N_test_rot=N_test_rot+1; }

}
else{
alea=ran2(&idum);
nmog=alea*(Natom-3)+3; /* Bound rotated */

if(nmog == Natom) { nmog= Natom-1; }



alea=ran2(&idum); /* orientation of the rotation */

if(alea<0.5) { signe=1; } else {signe=-1; }
conformacionmogp=conformacioatoms[nmog]+signe;
if(conformacionmogp>GAUCHEN) {conformacionmogp=TRANS;}
if(conformacionmogp<TRANS) {conformacionmogp=GAUCHEN;}

angle_Rot_Conformacio=-signe*120.;   /* criterio: Mas cercano-Sentido horario: t->g+->g- */
}
}
gencoordp(nmog,angle_Rot_Conformacio);

}

/**********************************************/
int acceptar(void)
{
float alea;
double difE,bf,MComega;



if(moviment!=OK)
 {  return(NO); }   /* NO per moviment=SOBREPOSAT */

     difE=etotkcalp-etotkcal;
     bf=-difE/(Rkcal*temp);

     if(bf<-709)
      {
return(NO);}   /* REBUTJAT PER INCREMENT MASSA GRAN DE L'ENERGIA */
     else
      {
      MComega=exp(bf);

      alea=ran2(&idum);

      if(MComega>=alea)
        { return(SI); }
      else
        { return(NO); }
      }

}

/**********************************************/
void actualitza(void)
{
int nat,nat1,nat2;
char C[2]="C", N[2]="N";

/* Actualize structure */

if(nmog > 2) {conformacioatoms[nmog]=conformacionmogp;}
Phi[nmog]=Phi_p;

if (nmog == 1) {Angles[0]=Theta_p; Phi[nmog]=cos(Phi_p); N_accepted_rot=N_accepted_rot+1; }
if (nmog == 2) {Phi[nmog]=cos(Phi_p);}



for(nat=1;nat<= Natom;nat++)
  {coord[nat][0]=coordp[nat][0];
   coord[nat][1]=coordp[nat][1];
   coord[nat][2]=coordp[nat][2]; }

for(nat1=1;nat1<= nmog-1;nat1++)
 {
  for(nat2=nmog+1;nat2<=Natom;nat2++)
    {
    matriu_dist2[nat1][nat2]=matriu_dist2p[nat1][nat2];
    }
 }

dist_extrems=dist_extremsp;

/* Actualize conformational energy */
etotkcal=etotkcalp;
Ebond[nmog]=ener_bond_nmog;
Ebond[nmog+1]=ener_bond_nmog1;
ener_conf=ener_conf+ener_confp;
ener_prot=ener_protp;

/* Actualize  Debye interaction Energy */
if(Debye== 0)
{
for(nat1=1;nat1<=nmog-1;nat1++)
 {
  for(nat2=nmog+2;nat2<=Natom;nat2++)
    {  // Done even if the atoms are uncharged. 
        matriu_EDebye[nat1][nat2]=matriu_EDebyep[nat1][nat2];
     }
  }
EDebye=EDebye+EDebyep-EDebyec;
}


/* Actualize Mechanical work  */

if ( F != 0) {

W_mech=W_mechp;

}


/* Actualize Bond from stretching     */

if ( control_str != 0) {
Bounds[nbond] = p_length;
control_str=0;

}

/* Actualize angle from bending     */

if ( control_bnd != 0) {
Angles[nangle] = p_angle;
control_bnd=0; 

}


}





/**************************************************************/
void MC_SGC_Aminas(void)
{
float alea;
int procesSGC,nat;
double SGComega,difE,bf,RT;

RT=Rkcal*temp;

/* If the atom is non-ionizable, another atom changes its ionization state */




while(pKa[nmog] == NonProt){
alea=ran2(&idum);
nmog= (Natom+1)*alea+1;


if(nmog==Natom+1) { nmog=Natom; }


}







ener_protp=get_Eprot(confSGC);
etotkcalp=etotkcal-ener_prot+ener_protp;


/* Debye long range interaction*/


if(Debye == 0)
 {
 EDebyec=get_EDebye_SGC(nmog);
 EDebyep=get_EDebye_SGCp(nmog);

 etotkcalp=etotkcalp+EDebyep-EDebyec;
 }

//////////  Is it ionization change acepted?  //////////


       difE=etotkcalp-etotkcal;

      if( carrega[nmog] == L_charge[nmog])
      {
      bf=-difE/RT+log(10)*(pKa[nmog]-pH);
      }
      else
      {
      bf=-difE/RT-log(10)*(pKa[nmog]-pH);
      }




     if(bf<-709)

     { } /* The ionization change is rejected because the energy increase is too big */

     else
     {

      SGComega=exp(bf);

      alea=ran2(&idum);


      if(SGComega>=alea)
        { /* Actualitza Energia */
             etotkcal=etotkcalp;
             ener_prot=ener_protp;
             if(carrega[nmog]==L_charge[nmog])
                { carrega[nmog]=H_charge[nmog];
                  N_prot=N_prot+1;}
              else
                { carrega[nmog]=L_charge[nmog];
                       N_prot=N_prot-1; }

                 /* Debye energy actualization */
                 for(nat=1;nat<=Natom;nat++)
                   {
                      if(nat<nmog)
                      { matriu_EDebye[nat][nmog]=matriu_EDebyep[nat][nmog]; }
                      else
                      { matriu_EDebye[nmog][nat]=matriu_EDebyep[nmog][nat]; }
                    }
                EDebye=EDebye+EDebyep-EDebyec;
                EincprotDeb=fabs(EDebyep-EDebyec);
        }
      else
        { } /* The ionization change is rejected */

  }




}

void generabondangle(void){

float alea;
double angle_Rot_Conformacio, Prob_str=0.5;
int signe;


alea=ran2(&idum);

if (alea < Prob_str) {

alea=ran2(&idum);

nbond=alea*(Natom-1)+1; /* Bond streteched */

if(nbond == Natom) { nbond= Natom-1; }

gencoordpstr(nbond);

control_str=1; // to indicate that we have done a stretching 


}
else

{

alea=ran2(&idum);

nangle=alea*(Natom-2)+1; /* angle bended */

if(nangle == Natom-1) { nangle = Natom-2; }

gencoordpbnd(nangle);

control_bnd=1; // to indicate that we have done a bending 


}








}


void actualitza_strbnd(void)
{
int nat,nat1,nat2;
/* Actualize structure */

for(nat=1;nat<= Natom;nat++)
  {coord[nat][0]=coordp[nat][0];
   coord[nat][1]=coordp[nat][1];
   coord[nat][2]=coordp[nat][2]; }

for(nat1=1;nat1<= nmog-1;nat1++)
 {
  for(nat2=nmog+2;nat2<=Natom;nat2++)
    {
    matriu_dist2[nat1][nat2]=matriu_dist2p[nat1][nat2];
    }
 }

dist_extrems=dist_extremsp;

/* Actualize conformational energy */
etotkcal=etotkcalp;
Ebond[nmog]=ener_bond_nmog;
Ebond[nmog+1]=ener_bond_nmog1;
ener_conf=ener_conf+ener_confp;
ener_prot=ener_protp;

/* Actualize  Debye interaction Energy */
if(Debye== 0)
{
for(nat1=1;nat1<=nmog-1;nat1++)
 {
  for(nat2=nmog+2;nat2<=Natom;nat2++)
    {  // Done even if the atoms are uncharged. 
        matriu_EDebye[nat1][nat2]=matriu_EDebyep[nat1][nat2];
     }
  }
EDebye=EDebye+EDebyep-EDebyec;

}


/* Actualize Mechanical work  */

if ( F != 0) {

W_mech=W_mechp;

}


/* Actualize Bond from stretching     */

if ( control_str != 0) {
Bounds[nbond] = p_length;
control_str=0;

}

/* Actualize angle from bending     */

if ( control_bnd != 0) {
Angles[nangle] = p_angle;
control_bnd=0;

}








}





