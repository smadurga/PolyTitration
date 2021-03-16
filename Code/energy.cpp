// Module with all the subroutines needed to calculate the energy of the system

#include "lib.h"
#include "global.h"

                                                                                                               
void get_Ebound(int nbound, double *EnerBound){
int Conf_Prev, Conf_bound;


// The first bound always is in trans

if(nbound==1)
       {
        *EnerBound=0.;
       }
      else
       {

          // Conformation of the bound of interest

          Conf_Prev=conformacioatoms[nbound-1];
          Conf_bound=conformacioatoms[nbound];

          // to calculate the energy asociated to a rotated bound 
          if(nmog==nbound-1){ Conf_Prev=conformacionmogp;}
          if(nmog==nbound){ Conf_bound=conformacionmogp;}



          // Energy transfer matrix of the bound

          if(1 == Mat_ind[nbound-1]){*EnerBound = Econf1[Conf_Prev][Conf_bound];}
          if(2 == Mat_ind[nbound-1]){*EnerBound = Econf2[Conf_Prev][Conf_bound];}
          if(3 == Mat_ind[nbound-1]){*EnerBound = Econf3[Conf_Prev][Conf_bound];}
          if(4 == Mat_ind[nbound-1]){*EnerBound = Econf4[Conf_Prev][Conf_bound];}
          if(5 == Mat_ind[nbound-1]){*EnerBound = Econf5[Conf_Prev][Conf_bound];}
          if(6 == Mat_ind[nbound-1]){*EnerBound = Econf6[Conf_Prev][Conf_bound];}
          if(7 == Mat_ind[nbound-1]){*EnerBound = Econf7[Conf_Prev][Conf_bound];}
          if(8 == Mat_ind[nbound-1]){*EnerBound = Econf8[Conf_Prev][Conf_bound];}
          if(9 == Mat_ind[nbound-1]){*EnerBound = Econf9[Conf_Prev][Conf_bound];}
          if(10 == Mat_ind[nbound-1]){*EnerBound = Econf10[Conf_Prev][Conf_bound];}
            }

}

double get_Eprot(int tip)
{
int nat;
double eprot;
int carrega_aux[MAXATOMS], conf;

// Subroutine to calculate the energy of a given protonation state of the polymeric chain
// tip accounts for the previous change of the system (if tip == confSGC, it has changed the protonation state, if tip == confMC it only has been a conformational change
// carrega_aux account for the new charge state if the nmog atom changes its charge


// If the previous change of state is a change in the charge, calculates carrega_aux


if(tip==confSGC)
 { for(nat=1;nat<=Natom;nat++)
   {
    if(nat!=nmog)
     { carrega_aux[nat]=carrega[nat]; }
     else
     { if(carrega[nat] == H_charge[nat]) 
        {carrega_aux[nat]=L_charge[nat];} else {carrega_aux[nat]=H_charge[nat];}
     }
   }
 }


// If only has been a conformational change, it asignates carrega_aux=carrega


if(tip==confMC)
 { for(nat=1;nat<=Natom;nat++)
     { carrega_aux[nat]=carrega[nat]; }
 }



// Calculate the energical contribution of charge interaction of this state


eprot=0;

for(nat=1;nat<=Natom-N_neigh;nat=nat+N_neigh)
{


 if(carrega_aux[nat] != 0 && carrega_aux[nat+N_neigh] != 0)
   {
   if(tip==confMC && nmog==nat+2){conf=conformacionmogp;}else{conf=conformacioatoms[nat+2]; }
   if(conf == TRANS){eprot=eprot+(-Rkcal*temp*log(u_t[nat]));}
   if(conf == GAUCHEP){eprot=eprot+(-Rkcal*temp*log(u_gp[nat]));}
   if(conf == GAUCHEN){eprot=eprot+(-Rkcal*temp*log(u_gn[nat]));}
   }
}


return(eprot);
}

void calculener(void)
{
int nat1,nat2,i;
double dx,dy,dz,dist2,dist;



/* Calcul distancies */

moviment=OK;
for(nat1=1;nat1<=nmog-1;nat1++)
 {
  for(nat2=nmog+1;nat2<=Natom;nat2++)
    {
    dx=coordp[nat2][0]-coordp[nat1][0];
    dy=coordp[nat2][1]-coordp[nat1][1];
    dz=coordp[nat2][2]-coordp[nat1][2];
    dist2=dx*dx+dy*dy+dz*dz;

    matriu_dist2p[nat1][nat2]=dist2;

    dist=sqrt(dist2);
    if(abs(nat2-nat1)+1>=INTERRADIS_AT_AT)
     {



     if(dist<radiatom[nat1]+radiatom[nat2])
       {

 
moviment=SOBREPOSAT; }  /* Calcul per tots encara que 1 SOBREPOSAT!!! */
     }

    }
 }



/* Distance calculation for the chain ends */
    dx=coordp[Natom][0]-coordp[1][0];
    dy=coordp[Natom][1]-coordp[1][1];
    dz=coordp[Natom][2]-coordp[1][2];
    dist2=dx*dx+dy*dy+dz*dz;
    dist_extremsp=sqrt(dist2);

/* energy calculation */


   /* 1- Conformational contribution */
 get_Ebound(nmog,&ener_bond_nmog);
 get_Ebound(nmog+1,&ener_bond_nmog1);

ener_confp=+ener_bond_nmog+ener_bond_nmog1-Ebond[nmog]-Ebond[nmog+1];



 etotkcalp=etotkcal+ener_bond_nmog+ener_bond_nmog1-Ebond[nmog]-Ebond[nmog+1];


   /* 2- Protonation contribution */
 ener_protp=get_Eprot(confMC);

 etotkcalp=etotkcalp+ener_protp-ener_prot;



/* 3- long range interaction contribution*/
if(Debye==0)
 {
 EDebyec=get_EDebye(nmog,nmog+1);
 EDebyep=get_EDebyep(nmog,nmog+1);

 etotkcalp=etotkcalp+EDebyep-EDebyec;

 }

/* 4- Mechanical work contribution */

if ( F != 0)
{

W_mech= get_Wmech();
W_mechp= get_Wmechp();

etotkcalp=etotkcalp+W_mechp-W_mech;

}

/* 5- Stretching contribution */

if ( control_str != 0)
{

E_str= get_Estr();
E_strp= get_Estrp();

etotkcalp=etotkcalp+E_strp-E_str;
}

/* 6- Bending contribution */

if ( control_bnd != 0)
{

E_bdn= get_Ebdn();
E_bdnp= get_Ebdnp();

etotkcalp=etotkcalp+E_bdnp-E_bdn;

}




}


/**********************************************/
// Long range Debye energy betwen atoms 1 to atomfix1 and atomfix2 to Natom
// Usually: get_EDebye(nmog,nmog+1)
double get_EDebye(int atomfix1,int atomfix2)
{
int nat1,nat2;
double dist2,dist;
double enerJoules,EDebyekcal;


enerJoules=0;

for(nat1=1;nat1<=atomfix1-1;nat1++)
 {
  for(nat2=atomfix2+1;nat2<=Natom;nat2++)
    {
    if((abs(nat1-nat2)+1>= INTERDH_AT_AT))
      {
      if(carrega[nat1] != 0 && carrega[nat2] != 0)
        { enerJoules+=matriu_EDebye[nat1][nat2]; }

       }
     }
  }


EDebyekcal=enerJoules/JoulesCal/1000*NA;

return(EDebyekcal);

}



/**********************************************/
    
// Provisional Debye interaction energy between atom 1 to atomfix1 and atomfix2 to Natom
// Usually: get_EDebyep(nmog,nmog+1)
double get_EDebyep(int atomfix1,int atomfix2)
{
int nat1,nat2;
double dist2,dist;
double enJ,enerJoules,EDebyekcal;



enerJoules=0;

for(nat1=1;nat1<=atomfix1-1;nat1++)
 {
  for(nat2=atomfix2+1;nat2<=Natom;nat2++)
    {
    if((abs(nat2-nat1)+1>= INTERDH_AT_AT))
      {
      if(carrega[nat1] != 0 && carrega[nat2] != 0)
        {
         dist2=matriu_dist2p[nat1][nat2];
         dist=sqrt(dist2);

         enJ=(prefactDebye/dist)*(exp(-kappa*dist));
         matriu_EDebyep[nat1][nat2]=enJ;
         enerJoules+=enJ;
        }
       }
     }
  }


EDebyekcal=enerJoules/JoulesCal/1000*NA;



return(EDebyekcal);


}


/**********************************************/
// Calculation of the Debye energy change when atomfix1 changes its ionization state
// Usually: get_EDebye_SGC(nmog)  for SGC
double get_EDebye_SGC(int atomfix1)
{
int nat,nat1;
double dist2,dist;
double enJ,enerJoules,EDebyekcal;
int carrega2[MAXATOMS];

enerJoules=0;

for(nat1=1;nat1<=Natom;nat1++)
 {
    if((abs(atomfix1-nat1)+1>= INTERDH_AT_AT))
      {
      if(carrega[nat1] != 0  && carrega[atomfix1] != 0 )    /* The ionization state is permently changed, is not provisional */
        {

         if(nat1<atomfix1)
           { enJ=matriu_EDebye[nat1][atomfix1]; }
          else
           { enJ=matriu_EDebye[atomfix1][nat1]; }

         enerJoules+=enJ;
        }
       }
  }

EDebyekcal=enerJoules/JoulesCal/1000*NA;

return(EDebyekcal);

}

/**********************************************/
// Calculation of the provisional Debye energy change when atomfix1 changes its ionization state
// Usually: get_EDebye_SGCp(nmog)  for SGC
double get_EDebye_SGCp(int atomfix1)
{
int nat,nat1;
double dist2,dist;
double enJ,enerJoules,EDebyekcal;
int carrega2[MAXATOMS];

// Provisional ionization state (carrega2)
  for(nat=1;nat<= Natom;nat++)
   {
    if(nat!=nmog)
     { carrega2[nat]=carrega[nat]; }
     else
     { if(carrega[nat]== H_charge[nat])
        {carrega2[nat]=L_charge[nat];} else {carrega2[nat]=H_charge[nat];}
     }
   }

enerJoules=0;

for(nat1=1;nat1<=Natom;nat1++)
 {
    if((abs(atomfix1-nat1)+1>= INTERDH_AT_AT))
      {
      if(carrega2[nat1] != 0  && carrega2[atomfix1] != 0)
        {
         if(nat1<atomfix1)
           { dist2=matriu_dist2[nat1][atomfix1]; }
          else
           { dist2=matriu_dist2[atomfix1][nat1]; }

         dist=sqrt(dist2);

         enJ=(prefactDebye/dist)*(exp(-kappa*dist));
         matriu_EDebyep[nat1][atomfix1]=enJ;  /* If nat1 < atomfix1 */
         matriu_EDebyep[atomfix1][nat1]=enJ;  /* If atomfix1 < nat1 */
         enerJoules+=enJ;
        }
       }
  }


EDebyekcal=enerJoules/JoulesCal/1000*NA;


return(EDebyekcal);

}


/**********************************************/
// 
/**********************************************/

double get_Wmech()
{
double EForcekcal,W,dist2,dist;


dist=coord[Natom][0]-coord[1][0];


W=-F*dist;  /* Fext en [pN] ;  dist en [angs] */
               /* canvi pN a N-> 1E-12 */
               /* canvi angs a m -> 1E-10 */
               /* canvi part a mol -> NA */
               /* canvi Jouls a calories -> JoulesCal */
               /* canvi cal a kcal -> *1E-3 */

EForcekcal=W*1E-25*NA/JoulesCal;
return(EForcekcal);

}

double get_Wmechp()
{
double EForcekcal,W,dist2,dist;


dist=coordp[Natom][0]-coordp[1][0];


W=-F*dist;  /* Fext en [pN] ;  dist en [angs] */
               /* canvi pN a N-> 1E-12 */
               /* canvi angs a m -> 1E-10 */
               /* canvi part a mol -> NA */
               /* canvi Jouls a calories -> JoulesCal */
               /* canvi cal a kcal -> *1E-3 */

EForcekcal=W*1E-25*NA/JoulesCal;

return(EForcekcal);

}

double get_Estr()
{
double Estrkcal;
int nat;

Estrkcal=0;

for(nat=1;nat<=Natom-1;nat++)
 {

 Estrkcal=Estrkcal+kstr*pow(Bounds[nat]-Bounds_eq[nat],2)/2.;

}
return(Estrkcal);

}


double get_Estrp()
{
double Estrkcalp;
int nat;

Estrkcalp=0;

for(nat=1;nat<=Natom-1;nat++)
 {
 if ( nat == nbond){
 Estrkcalp=Estrkcalp+kstr*pow(p_length-Bounds_eq[nat],2)/2.;

 }
 else{
 Estrkcalp=Estrkcalp+kstr*pow(Bounds[nat]-Bounds_eq[nat],2)/2.;

 }
}
return(Estrkcalp);

}

double get_Ebdn()
{
double Ebdnkcal;
int nat;

Ebdnkcal=0;

for(nat=1;nat<=Natom-2;nat++)
 {

 Ebdnkcal=Ebdnkcal+kbdn*pow(Angles[nat]-Angles_eq[nat],2)/2.;

}
return(Ebdnkcal);

}


double get_Ebdnp()
{
double Ebdnkcalp;
int nat;

Ebdnkcalp=0;

for(nat=1;nat<=Natom-2;nat++)
 {
 if ( nat == nangle){
 Ebdnkcalp=Ebdnkcalp+kbdn*pow(p_angle-Angles_eq[nat],2)/2.;

 }
 else{
 Ebdnkcalp=Ebdnkcalp+kbdn*pow(Angles[nat]-Angles_eq[nat],2)/2.;

 }
}
return(Ebdnkcalp);

}


