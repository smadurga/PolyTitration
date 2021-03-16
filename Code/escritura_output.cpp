// Module with all the subroutines concerning the output writing

#include "lib.h"
#include "global.h"

/*******************************************************/


void initial_output()
{

// Output with all relevant input data

  int n,i;
  ofstream posIni("../Outputs/initial_data.txt");

  

// Global simulation parameters

  posIni << "----------------------------------------------------------------" << endl;
  posIni << "***********************Input data check*************************" << endl;
  posIni << "----------------------------------------------------------------" << endl;
  posIni << "" << endl;
  posIni << "Number of atoms: " << Natom << endl;
  posIni << "System temperature: " << temp << endl;
  posIni << "Number of Configurations: " << conffin << endl;
  posIni << "Random seed: " << idum << endl;
  posIni << "System pH: " << pH << endl;
  posIni << "System ionic strengh: " << I << endl;
  posIni << "Chain:                ";
  for(i=1; i <= Natom; i++){if( i < Natom){posIni << chain[i]<< "-";}else{posIni << chain[i]<< endl;}}

  posIni << "Pka list:             ";
  for(i=1; i <= Natom; i++){posIni << pKa[i]<< " ";}
  posIni << endl;


  posIni << "Atomic Radius:        ";
  for(i=1; i <= Natom; i++){posIni << radiatom[i]<< " ";}
  posIni << endl;



  posIni << "Higher charge values: ";
  for(i=1; i <= Natom; i++){posIni << H_charge[i]<< " ";}
  posIni << endl;
  posIni << "Lower charge values:  ";
  for(i=1; i <= Natom; i++){posIni << L_charge[i]<< " ";}
  posIni << endl;



  posIni << "Transfer Matrix list:  ";
  for(i=1; i <= Natom-1; i++){posIni << Mat_ind[i]<< " ";}
  posIni << endl;


  posIni << "Are the long range interaction (Debye potential) included? "; 

  if(Debye == 0){posIni<< "Yes" << endl;}else{posIni<< "No" << endl;}

// Initial polymer propierties

  posIni << "" << endl;
  posIni << "--------------Initial properties of the polymer---------------" << endl;
  posIni << "" << endl;

  posIni << "Conformational energy: " << ener_conf << endl;
  posIni << " Protonation energy: " << ener_prot << endl;
  posIni << " Debye long range interaction energy: " << EDebye << endl;
  posIni << " Mechanical Work: " << W_mech << endl;
  posIni << " Total energy: " << etotkcal;
  posIni << "" << endl;
  posIni << "Flory aproximated value of r^2: " << r2endToend <<  endl;
  posIni << "Flory aproximated value of s^2: " << s2_free <<  endl;



// Initial structure of the polymer in xyz format



  posIni << "" << endl;
  posIni << "------Initial structure of the polymer in xyz format-----------" << endl;
  posIni << "" << endl;

  posIni << Natom << endl;
  posIni << "POLIAMINA $ TRANS  $ NUMATOMS " << Natom << endl;

  for(n=1; n <= Natom; n++){ posIni << chain[n] << " " << coord[n][0] << " " << coord[n][1] << " " << coord[n][2] << endl;}


// Transfer Matrix list


  posIni << "" << endl;
  posIni << "----------------------Transfer Matrix-------------------------" << endl;
  posIni << "" << endl;


  posIni << "Matrix 1: " ;
  if ( Nmat >= 1){posIni << labelA1[1] << " " << labelA2[1] << " " << labelA3[1] << " " << labelA4[1]<< endl;}else{posIni << labelA1[1] << " " << labelA2[1] << " " << labelA3[1] << endl;}
  posIni << "" << endl;
  posIni << Conf1[1][1] << " " << Conf1[1][2]<< " " << Conf1[1][3] << endl;
  posIni << Conf1[2][1] << " " << Conf1[2][2]<< " " << Conf1[2][3] << endl;
  posIni << Conf1[3][1] << " " << Conf1[3][2]<< " " << Conf1[3][3] << endl;
  posIni << "" << endl;
  posIni << "u_trans: " <<u_taux[1] << " u_gosh+: " <<u_gpaux[1] <<  " u_gosh-: " <<u_gnaux[1] << endl;
  posIni << "" << endl;
  

  posIni << "Matrix 2: ";
  if ( Nmat >= 2){posIni << labelA1[2] << " " << labelA2[2] << " " << labelA3[2] << " " << labelA4[2]<< endl;}else{posIni << labelA1[2] << " " << labelA2[2] << " " << labelA3[2] << endl;}
  posIni << "" << endl;
  posIni << Conf2[1][1] << " " << Conf2[1][2]<< " " << Conf2[1][3] << endl;
  posIni << Conf2[2][1] << " " << Conf2[2][2]<< " " << Conf2[2][3] << endl;
  posIni << Conf2[3][1] << " " << Conf2[3][2]<< " " << Conf2[3][3] << endl;
  posIni << "" << endl;
  posIni << "u_trans: " <<u_taux[2] << " u_gosh+: " <<u_gpaux[2] <<  " u_gosh-: " <<u_gnaux[2] << endl;
  posIni << "" << endl;  

  posIni << "Matrix 3: " ;
  if ( Nmat >= 3){posIni << labelA1[3] << " " << labelA2[3] << " " << labelA3[3] << " " << labelA4[3]<< endl;}else{posIni << labelA1[3] << " " << labelA2[3] << " " << labelA3[3] << endl;}
  posIni << "" << endl;
  posIni << Conf3[1][1] << " " << Conf3[1][2]<< " " << Conf3[1][3] << endl;
  posIni << Conf3[2][1] << " " << Conf3[2][2]<< " " << Conf3[2][3] << endl;
  posIni << Conf3[3][1] << " " << Conf3[3][2]<< " " << Conf3[3][3] << endl;
  posIni << "" << endl;
  posIni << "u_trans: " <<u_taux[3] << " u_gosh+: " <<u_gpaux[3] <<  " u_gosh-: " <<u_gnaux[3] << endl;
  posIni << "" << endl;

  posIni << "Matrix 4: ";
  if ( Nmat >= 4){posIni << labelA1[4] << " " << labelA2[4] << " " << labelA3[4] << " " << labelA4[4]<< endl;}else{posIni << labelA1[4] << " " << labelA2[4] << " " << labelA3[4] << endl;}
  posIni << "" << endl;
  posIni << Conf4[1][1] << " " << Conf4[1][2]<< " " << Conf4[1][3] << endl;
  posIni << Conf4[2][1] << " " << Conf4[2][2]<< " " << Conf4[2][3] << endl;
  posIni << Conf4[3][1] << " " << Conf4[3][2]<< " " << Conf4[3][3] << endl;
  posIni << "" << endl;
  posIni << "u_trans: " <<u_taux[4] << " u_gosh+: " <<u_gpaux[4] <<  " u_gosh-: " <<u_gnaux[4] << endl;
  posIni << "" << endl;

  posIni << "Matrix 5: ";
  if ( Nmat >= 5){posIni << labelA1[5] << " " << labelA2[5] << " " << labelA3[5] << " " << labelA4[5]<< endl;}else{posIni << labelA1[5] << " " << labelA2[5] << " " << labelA3[5] << endl;}
  posIni << "" << endl;
  posIni << Conf5[1][1] << " " << Conf5[1][2]<< " " << Conf5[1][3] << endl;
  posIni << Conf5[2][1] << " " << Conf5[2][2]<< " " << Conf5[2][3] << endl;
  posIni << Conf5[3][1] << " " << Conf5[3][2]<< " " << Conf5[3][3] << endl;
  posIni << "" << endl;
  posIni << "u_trans: " <<u_taux[5] << " u_gosh+: " <<u_gpaux[5] <<  " u_gosh-: " <<u_gnaux[5] << endl;
  posIni << "" << endl;

  posIni << "Matrix 6: " ;
  if ( Nmat >= 6){posIni << labelA1[6] << " " << labelA2[6] << " " << labelA3[6] << " " << labelA4[6]<< endl;}else{posIni << labelA1[6] << " " << labelA2[6] << " " << labelA3[6] << endl;}
  posIni << "" << endl;
  posIni << Conf6[1][1] << " " << Conf6[1][2]<< " " << Conf6[1][3] << endl;
  posIni << Conf6[2][1] << " " << Conf6[2][2]<< " " << Conf6[2][3] << endl;
  posIni << Conf6[3][1] << " " << Conf6[3][2]<< " " << Conf6[3][3] << endl;
  posIni << "" << endl;
  posIni << "u_trans: " <<u_taux[6] << " u_gosh+: " <<u_gpaux[6] <<  " u_gosh-: " <<u_gnaux[6] << endl;
  posIni << "" << endl;

  posIni << "Matrix 7: " ;
  if ( Nmat >= 7){posIni << labelA1[7] << " " << labelA2[7] << " " << labelA3[7] << " " << labelA4[7]<< endl;}else{posIni << labelA1[7] << " " << labelA2[7] << " " << labelA3[7] << endl;}
  posIni << "" << endl;
  posIni << Conf7[1][1] << " " << Conf7[1][2]<< " " << Conf7[1][3] << endl;
  posIni << Conf7[2][1] << " " << Conf7[2][2]<< " " << Conf7[2][3] << endl;
  posIni << Conf7[3][1] << " " << Conf7[3][2]<< " " << Conf7[3][3] << endl;
  posIni << "" << endl;
  posIni << "u_trans: " <<u_taux[7] << " u_gosh+: " <<u_gpaux[7] <<  " u_gosh-: " <<u_gnaux[7] << endl;
  posIni << "" << endl;

  posIni << "Matrix 8: ";
  if ( Nmat >= 8){posIni << labelA1[8] << " " << labelA2[8] << " " << labelA3[8] << " " << labelA4[8]<< endl;}else{posIni << labelA1[8] << " " << labelA2[8] << " " << labelA3[8] << endl;}
  posIni << "" << endl;
  posIni << Conf8[1][1] << " " << Conf8[1][2]<< " " << Conf8[1][3] << endl;
  posIni << Conf8[2][1] << " " << Conf8[2][2]<< " " << Conf8[2][3] << endl;
  posIni << Conf8[3][1] << " " << Conf8[3][2]<< " " << Conf8[3][3] << endl;
  posIni << "" << endl;
  posIni << "u_trans: " <<u_taux[8] << " u_gosh+: " <<u_gpaux[8] <<  " u_gosh-: " <<u_gnaux[8] << endl;
  posIni << "" << endl;

  posIni << "Matrix 9: " ;
  if ( Nmat >= 9){posIni << labelA1[9] << " " << labelA2[9] << " " << labelA3[9] << " " << labelA4[9]<< endl;}else{posIni << labelA1[9] << " " << labelA2[9] << " " << labelA3[9] << endl;}
  posIni << "" << endl;
  posIni << Conf9[1][1] << " " << Conf9[1][2]<< " " << Conf9[1][3] << endl;
  posIni << Conf9[2][1] << " " << Conf9[2][2]<< " " << Conf9[2][3] << endl;
  posIni << Conf9[3][1] << " " << Conf9[3][2]<< " " << Conf9[3][3] << endl;
  posIni << "" << endl;
  posIni << "u_trans: " <<u_taux[9] << " u_gosh+: " <<u_gpaux[9] <<  " u_gosh-: " <<u_gnaux[9] << endl;
  posIni << "" << endl;

  posIni << "Matrix 10: " ;
  if ( Nmat >= 10){posIni << labelA1[10] << " " << labelA2[10] << " " << labelA3[10] << " " << labelA4[10]<< endl;}else{posIni << labelA1[10] << " " << labelA2[10] << " " << labelA3[10] << endl;}
  posIni << "" << endl;
  posIni << Conf10[1][1] << " " << Conf10[1][2]<< " " << Conf10[1][3] << endl;
  posIni << Conf10[2][1] << " " << Conf10[2][2]<< " " << Conf10[2][3] << endl;
  posIni << Conf10[3][1] << " " << Conf10[3][2]<< " " << Conf10[3][3] << endl;
  posIni << "" << endl;
  posIni << "u_trans: " <<u_taux[10] << " u_gosh+: " <<u_gpaux[10] <<  " u_gosh-: " <<u_gnaux[10] << endl;


}

void escriuredades()
{
  int n,nat1,nat2,qt,i;
  double s2;
  char C[2]="C", N[2]="N";

if (nrank == 0 && sample == Nini) {
if(nconf%intervalOUT==0) {
ofstream posEne("../Outputs/energy_output.txt", std::ios_base::app);
posEne << "NCONF: " << nconf << ", ETOT(kcal/mol)= " << etotkcal <<",  Econf=  "<< ener_conf << ",   Eprot= " << ener_prot << ",   EDebye=  " << EDebye << ",   W=  "<< W_mech  << ",   Estr=  "<< E_str  <<  ",   Ebdn=  "<< E_bdn  << endl;
posEne.close();
 }
}


if(nconf%intervalAVG==0 || nconf%intervalOUT==0)
{

     /* s2 calculation*/
     s2=0;
     for(nat1=1;nat1< Natom;nat1++)
       { for(nat2=nat1+1;nat2<=Natom;nat2++)
         { s2+=matriu_dist2[nat1][nat2]; }
       }
      s2=s2/(Natom*Natom);

      s2_global=s2;

      /*  Global charge calculation */
     qt=0;
     for(nat1=1;nat1<= Natom;nat1++)
      {qt+=carrega[nat1];}


      /* Lx calculation */

     Lx=(coord[Natom][0]-coord[1][0])/(prom_bound*(Natom-1));

     
     /* number of bonds in trans computation */

     N_CC_trans=0;
     N_NC_trans=0;
     N_trans=0;
     for (i=3; i <= Natom-1 ; i++){ 
     if (/*conformacioatoms[i] == TRANS*/ int(Phi[i]) % 180 == 0) {
     N_trans++;
     if ( chain[i-1] == C[0] && chain[i] == C[0]){N_CC_trans++;}
     if ( chain[i-1] == N[0] && chain[i] == C[0]){N_NC_trans++;}
     if ( chain[i-1] == C[0] && chain[i] == N[0]){N_NC_trans++;}
     }
 

     }
}

if (nrank == 0 && sample == Nini) {
if(nconf%intervalOUT==0)
{

  cout << "I'm in Configuration = " << nconf << endl;
  ofstream posEns("../Outputs/ensemble.xyz", std::ios_base::app);

    posEns << Natom << endl;
    posEns << "NCONF  " << nconf  << "   ";

    for(n=1;n<= Natom;n++)
     { if(conformacioatoms[n]==TRANS) 
         { if(carrega[n]== L_charge[n]) { posEns << "t";}else{ posEns << "T";}}
       if(conformacioatoms[n]==GAUCHEP) 
         { if(carrega[n]== L_charge[n]) { posEns << "g";}else{ posEns << "G";}}
       if(conformacioatoms[n]==GAUCHEN) 
         { if(carrega[n]== L_charge[n]) { posEns << "j";}else{ posEns << "J";}}
     }



    posEns << " r=  " << dist_extrems << " r2= " << dist_extrems*dist_extrems << " s2= " << s2 << " Q= " << qt;

    posEns << " ETOT(kcal/mol)= " << etotkcal << " Eprot= " << ener_prot <<   " EDebye= " << EDebye << " W= "<< W_mech  << endl;

    for(n=1;n<= Natom;n++)
     { posEns << chain[n] << " " << coord[n][0] << " " << coord[n][1] << " " << coord[n][2] << endl;}


  posEns.close();
}
}

/* Averages calculation */

if(nconf%intervalAVG==0)
{
if (dist_extrems == dist_extrems && s2 == s2){
acc_r2_100+=dist_extrems*dist_extrems;
acc_s2_100+=s2;
acc_qt_100+=qt;
acc_Lx_100+=Lx;
acc_cov_100+=N_prot/N_sites;
acc_N_CC_trans_100+=double(N_CC_trans)/N_CC;
acc_N_NC_trans_100+=double(N_NC_trans)/N_NC;
acc_N_trans_100+=double(N_trans)/(Natom-3);
EDebye_acc+=EDebye;
EincprotDeb_acc+=EincprotDeb;
n_acc_100++;


if(100*nconf/conffin>10)
{ acc_r2_90+=dist_extrems*dist_extrems;
  acc_s2_90+=s2;
  acc_qt_90+=qt;
  acc_Lx_90+=Lx;
  acc_cov_90+=N_prot/N_sites;
  n_acc_90++; }

if(100*nconf/conffin>20)
{ acc_r2_80+=dist_extrems*dist_extrems;
  acc_s2_80+=s2;
  acc_qt_80+=qt;
  acc_Lx_80+=Lx;
  acc_cov_80+=N_prot/N_sites;
  n_acc_80++; }

if(100*nconf/conffin>50)
{ acc_r2_50+=dist_extrems*dist_extrems;
  acc_s2_50+=s2;
  acc_qt_50+=qt;
  acc_Lx_50+=Lx;
  acc_cov_50+=N_prot/N_sites;
  n_acc_50++; }


}


}

if (nrank == 0 && sample == Nini) {
if(nconf%intervalAVG==0)
{
ofstream posCon("../Outputs/covergence_output.txt", std::ios_base::app);
posCon << "NCONF: " << nconf << " r2: " << acc_r2_100/n_acc_100 <<" s2:  "<< acc_s2_100/n_acc_100 << " Lx: " << acc_Lx_100/n_acc_100 << ",   cov:  " << acc_cov_100/n_acc_100 <<  endl;
posCon.close();


}
}

}

/**********************************************/
// Subroutine to write the final output
/**********************************************/



extern void write_final_data(void)
{
int i,nat;
double L_pers100, L_pers90, L_pers80, L_pers50;
double L_kuhn100, L_kuhn90, L_kuhn80, L_kuhn50;


L_pers100= r2_av100/(2*nbounds*prom_bound)+prom_bound/2.;
L_pers90= r2_av90/(2*nbounds*prom_bound)+prom_bound/2.;
L_pers80= r2_av80/(2*nbounds*prom_bound)+prom_bound/2.;
L_pers50= r2_av50/(2*nbounds*prom_bound)+prom_bound/2.;

L_kuhn100= r2_av100/(nbounds*prom_bound);
L_kuhn90= r2_av90/(nbounds*prom_bound);
L_kuhn80= r2_av80/(nbounds*prom_bound);
L_kuhn50= r2_av50/(nbounds*prom_bound);



ofstream final_output("../Outputs/final_data.txt");


 final_output << "-------------------" << endl;
 final_output << "N_Bond MOV_SI/NO TOTAL" << endl;
    for(nat=1;nat<= Natom-1;nat++)
     {
     final_output << nat << " " << " " << vectormovsi[nat] << " " << vectormovno[nat] << " " << vectormovsi[nat]+vectormovno[nat] << endl;
     }

 final_output << "-------------------" << endl;
 final_output << "Configuracions MC acceptades:  "<<  100.*(float)movsi/((float)(movsi+movno))<< endl;


 final_output << "Configuracions analitzades per AVG 100,90,80,50: " << n_acc_100 << " " << n_acc_90 << " " << n_acc_80 << " " << n_acc_50 << endl;


  final_output << "AVG100%: r2= " << r2_av100 << " s2= " << s2_av100 << " q= " << q_av100 <<  " Cn= " << r2_av100/(pow(prom_bound,2)*(Natom-1))  << "  Lx = " << Lx_av100   <<" L_pers= " << L_pers100 << " L_kuhn= " << L_kuhn100 <<  " Cov= " << cov_av100 << " log(Kc)= " << log10(cov_av100/((1-cov_av100)*pow(10,-pH))) << endl;
  final_output << "AVG90%: r2= " << r2_av90 << " s2= " << s2_av90 << " q= " << q_av90 << " Cn= " << r2_av90/(pow(prom_bound,2)*(Natom-1))<< "  Lx = " << Lx_av90  <<" L_pers= " << L_pers90 << " L_kuhn= " << L_kuhn90  << " Cov= " << cov_av90 << " log(Kc)= " << log10(cov_av90/((1-cov_av90)*pow(10,-pH))) << endl;
  final_output << "AVG80%: r2= " << r2_av80 << " s2= " << s2_av80 << " q= " << q_av80 << " Cn= " << r2_av80/(pow(prom_bound,2)*(Natom-1))<< "  Lx = " << Lx_av80 <<" L_pers= " << L_pers80 << " L_kuhn= " << L_kuhn80  << " Cov= " << cov_av80 << " log(Kc)= " << log10(cov_av80/((1-cov_av80)*pow(10,-pH))) <<  endl;
  final_output << "AVG50%: r2= " << r2_av50 << " s2= " << s2_av50 << " q= " << q_av50 << " Cn= " << r2_av50/(pow(prom_bound,2)*(Natom-1))<< "  Lx = " << Lx_av50   <<" L_pers= " << L_pers50 << " L_kuhn= " << L_kuhn50  << " Cov= " << cov_av50 << " log(Kc)= " << log10(cov_av50/((1-cov_av50)*pow(10,-pH))) <<  endl;

  final_output << "% Rotation acceptance: " << double(N_accepted_rot)/N_test_rot*100 << endl;
  final_output << "% stretching and bending acceptance: " << double(N_accepted_str)/N_test_str*100 << endl;
  final_output << "% C-C in trans: " << perCC << endl;
  final_output << "% C-N in trans: " << perNC << endl;
  final_output << "% bonds in trans: " << pertrans << endl;
  final_output << "Total average EDebye (kcal/mol): " <<  EDebye_av100 << endl;
  final_output << "Average EDebye increase in protonation/desprotonation: " << EincprotDeb_av100 << endl;

}

void compute_averages(void)
{


// r2 averages

r2_av100=r2_av100+acc_r2_100/n_acc_100;
r2_av90=r2_av90+acc_r2_90/n_acc_90;
r2_av80=r2_av80+acc_r2_80/n_acc_80;
r2_av50=r2_av50+acc_r2_50/n_acc_50;


// s2 averages

s2_av100=s2_av100+acc_s2_100/n_acc_100;
s2_av90=s2_av90+acc_s2_90/n_acc_90;
s2_av80=s2_av80+acc_s2_80/n_acc_80;
s2_av50=s2_av50+acc_s2_50/n_acc_50;


// cov averages

cov_av100=cov_av100+acc_cov_100/n_acc_100;
cov_av90=cov_av90+acc_cov_90/n_acc_90;
cov_av80=cov_av80+acc_cov_80/n_acc_80;
cov_av50=cov_av50+acc_cov_50/n_acc_50;

// q averages

q_av100=q_av100+acc_qt_100/n_acc_100;
q_av90=q_av90+acc_qt_90/n_acc_90;
q_av80=q_av80+acc_qt_80/n_acc_80;
q_av50=q_av50+acc_qt_50/n_acc_50;

// Lx averages

Lx_av100=Lx_av100+acc_Lx_100/n_acc_100;
Lx_av90=Lx_av90+acc_Lx_90/n_acc_90;
Lx_av80=Lx_av80+acc_Lx_80/n_acc_80;
Lx_av50=Lx_av50+acc_Lx_50/n_acc_50;

// EDebye average

EDebye_av100=EDebye_av100+EDebye_acc/n_acc_100;
EincprotDeb_av100=EincprotDeb_av100+EincprotDeb_acc/n_acc_100;

// % of trans average

perCC=perCC+acc_N_CC_trans_100/n_acc_100*100;
perNC=perNC+acc_N_NC_trans_100/n_acc_100*100;
pertrans=pertrans+acc_N_trans_100/n_acc_100*100;


}

void send_info(void)
{
double gr2_av100,gr2_av90,gr2_av80,gr2_av50, gs2_av100,gs2_av90,gs2_av80,gs2_av50, gcov_av100,gcov_av90,gcov_av80,gcov_av50, gq_av100,gq_av90,gq_av80,gq_av50, gLx_av100,gLx_av90,gLx_av80,gLx_av50,gEDebye_av100, gEincprotDeb_av100, gperCC, gperNC, gpertrans;


// r2 averages
r2_av100=r2_av100/(Nfin-Nini);
r2_av90=r2_av90/(Nfin-Nini);
r2_av80=r2_av80/(Nfin-Nini);
r2_av50=r2_av50/(Nfin-Nini);

// s2 averages

s2_av100=s2_av100/(Nfin-Nini);
s2_av90=s2_av90/(Nfin-Nini);
s2_av80=s2_av80/(Nfin-Nini);
s2_av50=s2_av50/(Nfin-Nini);

// cov averages

cov_av100=cov_av100/(Nfin-Nini);
cov_av90=cov_av90/(Nfin-Nini);
cov_av80=cov_av80/(Nfin-Nini);
cov_av50=cov_av50/(Nfin-Nini);


// q averages

q_av100=q_av100/(Nfin-Nini);
q_av90=q_av90/(Nfin-Nini);
q_av80=q_av80/(Nfin-Nini);
q_av50=q_av50/(Nfin-Nini);


// Lx averages
Lx_av100=Lx_av100/(Nfin-Nini);
Lx_av90=Lx_av90/(Nfin-Nini);
Lx_av80=Lx_av80/(Nfin-Nini);
Lx_av50=Lx_av50/(Nfin-Nini);

// EDebye average

EDebye_av100=EDebye_av100/(Nfin-Nini);
EincprotDeb_av100=EincprotDeb_av100/(Nfin-Nini);

// % of trans average

perCC=perCC/(Nfin-Nini);
perNC=perNC/(Nfin-Nini);
pertrans=pertrans/(Nfin-Nini);



// sum all the cpu information and send to master
MPI_Reduce(&r2_av100, &gr2_av100,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&r2_av90, &gr2_av90,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&r2_av80, &gr2_av80,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&r2_av50, &gr2_av50,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

MPI_Reduce(&s2_av100, &gs2_av100,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&s2_av90, &gs2_av90,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&s2_av80, &gs2_av80,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&s2_av50, &gs2_av50,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

MPI_Reduce(&cov_av100, &gcov_av100,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&cov_av90, &gcov_av90,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&cov_av80, &gcov_av80,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&cov_av50, &gcov_av50,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

MPI_Reduce(&q_av100, &gq_av100,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&q_av90, &gq_av90,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&q_av80, &gq_av80,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&q_av50, &gq_av50,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

MPI_Reduce(&Lx_av100, &gLx_av100,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&Lx_av90, &gLx_av90,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&Lx_av80, &gLx_av80,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&Lx_av50, &gLx_av50,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

MPI_Reduce(&EDebye_av100, &gEDebye_av100,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&EincprotDeb_av100, &gEincprotDeb_av100,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

MPI_Reduce(&perCC, &gperCC,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&perNC, &gperNC,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
MPI_Reduce(&pertrans, &gpertrans,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);





// master average over all cpus

if (nrank == 0){

r2_av100=gr2_av100/totalcpu;
r2_av90=gr2_av90/totalcpu;
r2_av80=gr2_av80/totalcpu;
r2_av50=gr2_av50/totalcpu;


s2_av100=gs2_av100/totalcpu;
s2_av90=gs2_av90/totalcpu;
s2_av80=gs2_av80/totalcpu;
s2_av50=gs2_av50/totalcpu;


cov_av100=gcov_av100/totalcpu;
cov_av90=gcov_av90/totalcpu;
cov_av80=gcov_av80/totalcpu;
cov_av50=gcov_av50/totalcpu;

q_av100=gq_av100/totalcpu;
q_av90=gq_av90/totalcpu;
q_av80=gq_av80/totalcpu;
q_av50=gq_av50/totalcpu;

Lx_av100=gLx_av100/totalcpu;
Lx_av90=gLx_av90/totalcpu;
Lx_av80=gLx_av80/totalcpu;
Lx_av50=gLx_av50/totalcpu;

EDebye_av100=gEDebye_av100/totalcpu;
EincprotDeb_av100=gEincprotDeb_av100/totalcpu;

perCC=gperCC/totalcpu;
perNC=gperNC/totalcpu;
pertrans=gpertrans/totalcpu;


}

}


