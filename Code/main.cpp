#include "lib.h"
#include "global.h"

int Natom, carrega[MAXATOMS], conformacioatoms[MAXATOMS], Mat_ind[MAXATOMS], H_charge[MAXATOMS], L_charge[MAXATOMS], Debye, Nmat, Ntranex, conformacionmogp, Ref_State[MAXATOMS], nbounds, vectormovsi[MAXATOMS],vectormovno[MAXATOMS],movsi,movno, n_acc_100,n_acc_90,n_acc_80,n_acc_50, moviment, N_neigh, nbond, nangle, control_str, control_bnd, N_test_rot, N_test_str, N_accepted_rot, N_accepted_str, N_CC, N_NC, N_CC_trans, N_NC_trans, N_trans, Nsim, sample;
char  chain[MAXATOMS], labelA1[MAXATOMS], labelA2[MAXATOMS], labelA3[MAXATOMS], labelA4[MAXATOMS];
double Bounds[MAXATOMS], Angles[MAXATOMS], pKa[MAXATOMS], coord[MAXATOMS][3], Ebond[MAXATOMS], Ebond_prot[MAXATOMS], u_t[MAXATOMS], u_gp[MAXATOMS], u_gn[MAXATOMS], matriu_dist2[MAXATOMS][MAXATOMS], acc_r2_100,acc_s2_100,acc_qt_100, acc_r2_90,acc_s2_90,acc_qt_90,acc_r2_80,acc_s2_80,acc_qt_80,acc_r2_50,acc_s2_50,acc_qt_50, u_taux[MAXTRANSFER+2], u_gpaux[MAXTRANSFER+2], u_gnaux[MAXTRANSFER+2], coordp[MAXATOMS][3], radiatom[MAXATOMS], kappa, prefactDebye, Phi[MAXATOMS], acc_Lx_100, acc_Lx_90, acc_Lx_80, acc_Lx_50, ProbConf, Bounds_eq[MAXATOMS], Angles_eq[MAXATOMS], acc_cov_100, acc_cov_90, acc_cov_80, acc_cov_50, acc_N_CC_trans_100, acc_N_NC_trans_100, acc_N_trans_100;
long nconf, idum, conffin, intervalAVG, intervalOUT, nmog, idum_ref; 
double temp, pH, DIEBULK, I,  probTest, etotkcal, etotkcalp,  ener_prot, ener_conf,  matriu_EDebye[MAXATOMS][MAXATOMS], r2endToend, s2_free, ener_bond_nmog, ener_bond_nmog1, ener_protp, matriu_EDebyep[MAXATOMS][MAXATOMS], dist_extrems,dist_extremsp,  Phi_p, F, W_mech, W_mechp, EDebye, EDebyep, ener_confp, Lx, prom_bound, Theta_p, Ly, Lz, kstr, kbdn, E_str, E_strp, p_length, p_angle, E_bdn, E_bdnp, N_prot, N_sites, r2_av100, s2_av100, cov_av100, perCC,  perNC, pertrans, r2_av90, s2_av90, cov_av90, r2_av80, s2_av80, cov_av80, r2_av50, s2_av50, cov_av50,  Lx_av100, q_av100, Lx_av90, q_av90, Lx_av80, q_av80, Lx_av50, q_av50, EDebyec, EDebye_acc, EDebye_av100, EincprotDeb, EincprotDeb_acc, EincprotDeb_av100;
double Conf1[3][3], Conf2[3][3], Conf3[3][3], Conf4[3][3], Conf5[3][3], Conf6[3][3], Conf7[3][3], Conf8[3][3], Conf9[3][3], Conf10[3][3], Econf1[3][3], Econf2[3][3],Econf3[3][3],Econf4[3][3],Econf5[3][3],Econf6[3][3],Econf7[3][3],Econf8[3][3],Econf9[3][3],Econf10[3][3];
double matriu_dist2p[MAXATOMS][MAXATOMS], s2_global;

// Open MPI variables

int totalcpu,nrank;
int MASTER=0;
int ndiv,res,icpu,Nini,Nfin;
MPI_Request request;



int main(int narg, char **arg){
int i, nat, n;
float alea;


ofstream posEne("../Outputs/energy_output.txt");
ofstream posEns("../Outputs/ensemble.xyz");
ofstream posCon("../Outputs/covergence_output.txt");
 MPI_Init(&narg,&arg);  // Inicializa MPI
 MPI_Comm_size(MPI_COMM_WORLD,&totalcpu);
 MPI_Comm_rank(MPI_COMM_WORLD,&nrank);
 if(nrank==MASTER) {
    printf("Total cpus = %d, proc ID = %d\n",totalcpu,nrank);
    }


lectura_input();

// Work division between cores. If it is not exact, last core takes extra work

ndiv=Nsim/totalcpu;
if (nrank != totalcpu-1){
Nini=nrank*ndiv;
Nfin=(nrank+1)*ndiv;


}
else
{
Nini=nrank*ndiv;
Nfin=Nsim;


}

// Inicialize all global averages to 0

r2_av100=0; s2_av100=0; cov_av100=0; perCC=0; perNC=0; pertrans=0; 
r2_av90=0; s2_av90=0; cov_av90=0;
r2_av80=0; s2_av80=0; cov_av80=0;
r2_av50=0; s2_av50=0; cov_av50=0;
EDebye_av100=0; EincprotDeb_av100=0;


for (sample=Nini; sample<Nfin; sample++) {

fes_polymer_trans(Natom);


init_energy();


 if (nrank == 0 && sample == Nini) {initial_output();}

 if  ( nrank == MASTER && sample == Nini) {cout << "Begining Monte Carlo simulation" << endl;}

for(nconf=1;nconf<=conffin;nconf++)
    {
     
    alea=ran2(&idum);
    if(alea<ProbConf){
       if(alea>probTest)
         {
          generaconf();
          calculener();
          if(acceptar()==SI)
           {vectormovsi[nmog]++; movsi++;  actualitza();  }
          else
           {vectormovno[nmog]++; movno++;  }
         }
         else
         {
         MC_SGC_Aminas();
         }
       }
     else{

     /* bond and angle stretching goes here  */

        generabondangle();
        calculener();
        N_test_str=N_test_str+1;
        if(acceptar()==SI) { actualitza_strbnd(); N_accepted_str=N_accepted_str+1;}



     }
     escriuredades();  /* always done */
      
} 

compute_averages();


}


send_info();

if (nrank == 0) {write_final_data();} // Final output

MPI_Finalize();

}
