// Module with all the global variables

#ifndef global_h
#define global_h

// Constants

#define MAXATOMS 400
#define MAXTORS 400
#define PI 3.14159265358979323846
#define SI  1
#define NO -1
#define OK  1
#define SOBREPOSAT -1
#define TRANS    1
#define GAUCHEP  2
#define GAUCHEN  3
#define ROT_A  1
#define ROT_B  2
#define ROT_C  3
#define PROTONAR    1
#define DESPROTONAR 2
#define confMC   1
#define confSGC  2
#define Rkcal 0.001987
#define NA 6.02214179E23    // mol-1
#define epsi0 8.854187817E-12   // A2s4kg-1 = C2N-1m-2 = F m-1
#define kBoltzmann 1.3806504E-23 // JK-1
#define q_elect 1.602176487E-19  // C
#define JoulesCal   4.184
#define NonProt -100 // Reference value to avoid atom ionization
#define MAXTRANSFER 10 // Maximum number of transfer matrix
#define INTERRADIS_AT_AT 5 // Minima distancia entre atomos permitida al realizar un movimiento 
#define INTERDH_AT_AT 2 // Minima distancia entre atomos para calcular su interacion de Debye

//Global variables

extern int Natom; // Total number of atoms of the polymer
extern char  chain[MAXATOMS]; // labels of each atom
extern double Bounds[MAXATOMS]; // bound lenght list
extern double Angles[MAXATOMS]; // bound angle list
extern double Bounds_eq[MAXATOMS]; // bound lenght of equilibrium list (initial)
extern double Angles_eq[MAXATOMS]; // bound angle of equilibrium list (initial
extern double Phi[MAXATOMS]; // internal angle list
extern double pKa[MAXATOMS]; // pKa list
extern long idum; // seed of the random number generator
extern long conffin; // Number of configurations
extern long intervalAVG; // number of steps until calculating the averages
extern long intervalOUT; // number of steps until writing the outputs
extern double temp;  // Temperature of the system (K)
extern double pH;    // pH of the media
extern double DIEBULK; // Relative dielectrical constant
extern double I;      // Ionic strengh of the media
extern double probTest; // Ionization probability
extern int conformacioatoms[MAXATOMS]; // Indicate the conformational state of the atom (TRANS GAUCHE+ GAUCHE-)
extern double coord[MAXATOMS][3]; // Atomic Cartesian coordinates
extern double coordp[MAXATOMS][3]; // Provisional atomic Cartesian coordinates
extern int carrega[MAXATOMS]; // Vector with the atomic charge
extern double Ebond[MAXATOMS]; // Vector with the bound energy due to conformation
extern double Ebond_prot[MAXATOMS]; // Vector with the bound energy due to protonation
extern double u_t[MAXATOMS]; // vector with all the interactions between protonated sites in trans 
extern double u_gp[MAXATOMS]; // vector with all the interactions between protonated sites in gosh+
extern double u_gn[MAXATOMS]; // vector with all the interactions between protonated sites in gosh-
extern int Mat_ind[MAXATOMS]; // vector with the rotation matrix asociated to each bond
extern int Debye; //  variable concerning if long range debye interactions are inclouded ( 0 = Yes, else = No) 
extern double etotkcal; // total energy of the system
extern double etotkcalp; // provisional total energy of the system
extern double EDebye; // energy of the system due to Debye long range interaction
extern double EDebyec; // change in Debye energy if nmog is rotated
extern double EDebyep; // provisional change in  energy of the system due to Debye long range interaction
extern double EincprotDeb; // change in the Debye energy due to a protonation
extern int Nmat; // number of final chain transfer matrix
extern char labelA1[MAXATOMS], labelA2[MAXATOMS],labelA3[MAXATOMS], labelA4[MAXATOMS]; // labels of the transfer matrix atoms
extern int Ntranex; // number of the first bound with transfer matrix
extern long nmog; //  number of the bound that have been rotated
extern int conformacionmogp; // new conformation of the rotated bound
extern double ener_prot; // energy due to protonation
extern double ener_conf; // energy due to the conformation
extern double ener_confp; // provisional energy increase due to the new conformation
extern double ener_protp; // provisional energy due to protonation
extern double W_mech; // work due to the mechanical force
extern double W_mechp; // provisional work due to the mechanical force
extern int H_charge[MAXATOMS]; // Higher protonation state of the atoms of the polymeric chain
extern int L_charge[MAXATOMS]; // Lower protonation state of the atoms of the polymeric chain
extern double matriu_dist2[MAXATOMS][MAXATOMS];  //matrix with the r2 distance between all the chain atoms
extern double matriu_dist2p[MAXATOMS][MAXATOMS]; //matrix with the provisional r2 distance between all the chain atoms
extern double matriu_EDebye[MAXATOMS][MAXATOMS]; // matrix with the long range interaction energy between two atoms
extern double matriu_EDebyep[MAXATOMS][MAXATOMS]; // matrix with the provisional long range interaction energy between two atoms
extern int nbounds; // Number of bounds of the chain 
extern int vectormovsi[MAXATOMS]; // Number of MC changes acepted for each atom
extern int vectormovno[MAXATOMS]; // Number of MC changes rejected for each atom
extern int movsi; // Total number of MC changes accepted
extern int movno; // Total number of MC changes rejected
extern double r2endToend; // end to end distance
extern double s2_free; // s2 calculated using the freely rotating chain model
extern double u_taux[MAXTRANSFER+2]; // auxiliar vector with all the interactions between protonated sites in trans 
extern double u_gpaux[MAXTRANSFER+2]; // auxiliar vector with all the interactions between protonated sites in gosh+
extern double u_gnaux[MAXTRANSFER+2]; // auxiliar vector with all the interactions between protonated sites in gosh-
extern long nconf; // Actual conformation of the simulation
extern int moviment; //integer containing if the conformational change is accepted or not
extern double radiatom[MAXATOMS]; // row containing the atomic radius of all the chain atom
extern double ener_bond_nmog; // energy change in the monomer nmog due to its conformational change
extern double ener_bond_nmog1; // energy change in the monomer nmog+1 due to the conformational change of nmog
extern double prefactDebye; // Debye prefactor 
extern double kappa; // Kappa constant
extern double dist_extrems; // end-to-end distance
extern double dist_extremsp; // provisional end-to-end distance
extern double s2_global; // global variable with the average s2 calculation
extern int N_neigh;
extern double Phi_p; // provisional value nmog internal angle
extern double Theta_p; // for the first bound, provisional value of its angle with x axis
extern double F; // Mechanical Force applied en x axis in pN
extern double Lx; // Normalized elongation in the x axis Lx/(Natom-1)
extern double prom_bound; // average bond lenght
extern double ProbConf; // Probability of conformational state change 
extern int nbond; // Number of the bond stretched
extern int nangle; // Number of the angle stretched
extern int control_str; // control variable of stretching
extern double kstr; //harmonic constant of bond strething in kcal/mol A-2
extern double kbdn; // harmonic constant of bond bending in kcal/mol deg-2
extern double E_str; // Energy due to bond stretching
extern double E_strp; // Provisional energy due to bond stretching
extern double p_length; // Provisional lenght of nbond
extern double p_angle; // Provisional angle of nangle
extern double E_bdn; // Energy due to bond bending
extern double E_bdnp; // Provisional energy due to bond bending
extern double N_prot; // Number of protonated sites
extern double N_sites; // Total number of sites
extern int control_bnd; // control variable of bending
extern int N_test_str; // number of times stretching and bending have been tested
extern int N_accepted_str; // number of times stretching and bending have been accepted
extern int N_test_rot; // number of times a rotation have been tested
extern int N_accepted_rot; // number of times a rotation  have been accepted
extern int N_CC; // total number of C-C bonds
extern int N_NC; // total number of C-N bonds
extern int N_CC_trans; // total number of C-C bonds in trans 
extern int N_NC_trans; // total number of N-C bonds in trans
extern int N_trans; // total bonds in trans
extern int Nsim; // Number of simulations performed and averaged
extern long idum_ref; // initial idum, never changed
extern int sample; // number of sample performed

// r2, s2, Lx, covering  and  qt accumulators

extern double acc_r2_100,acc_s2_100,  acc_qt_100, acc_Lx_100, acc_cov_100, acc_N_CC_trans_100, acc_N_NC_trans_100, acc_N_trans_100;
extern double acc_r2_90,acc_s2_90,  acc_qt_90, acc_Lx_90, acc_cov_90;
extern double acc_r2_80,acc_s2_80,  acc_qt_80, acc_Lx_80, acc_cov_80;
extern double acc_r2_50,acc_s2_50,  acc_qt_50, acc_Lx_50, acc_cov_50;
extern int n_acc_100,n_acc_90,n_acc_80,n_acc_50;
extern double EDebye_acc, EincprotDeb_acc;


// Transfer  matrix accounting for conformational energy

extern double Conf1[3][3];
extern double Conf2[3][3];
extern double Conf3[3][3];
extern double Conf4[3][3];
extern double Conf5[3][3];
extern double Conf6[3][3];
extern double Conf7[3][3];
extern double Conf8[3][3];
extern double Conf9[3][3];
extern double Conf10[3][3];

// Energy matrix accounting for conformational energy

extern double Econf1[3][3], Econf2[3][3],Econf3[3][3],Econf4[3][3],Econf5[3][3],Econf6[3][3],Econf7[3][3],Econf8[3][3],Econf9[3][3],Econf10[3][3];

// MPI variables

extern int totalcpu,nrank;
extern int MASTER;
extern int ndiv,res,icpu,Nini,Nfin;

// Global averages

extern double r2_av100, s2_av100, cov_av100, Lx_av100, q_av100, perCC,  perNC, pertrans;
extern double r2_av90, s2_av90, cov_av90, Lx_av90, q_av90;
extern double r2_av80, s2_av80, cov_av80, Lx_av80, q_av80;
extern double r2_av50, s2_av50, cov_av50, Lx_av50, q_av50;
extern double EDebye_av100;
extern double EincprotDeb_av100;


//Subroutines


extern void lectura_input();
extern void rotar_Rx(double xi,double yi,double zi,double angrota,double *xf,double *yf,double *zf);
extern void rotar_Ry(double xi,double yi,double zi,double angrota,double *xf,double *yf,double *zf);
extern void rotar_Rz(double xi,double yi,double zi,double angrota,double *xf,double *yf,double *zf);
extern void initial_output();
extern void fes_polymer_trans(int Natom);
extern void init_Emat();
extern void init_energy();
extern void get_Ebound(int nbound, double *Enerbound);
extern double get_Eprot(int tip); 
extern double r2_FreelyRotatingChain(int nbounds,double distancebond,double angle3sites);
extern double s2_FreelyRotatingChain(int nbounds,double distancebond,double angle3sites);
extern void generaconf(void);
extern void calculener(void);
extern int acceptar(void);
extern void actualitza(void);
extern void MC_SGC_Aminas(void);
extern void escriuredades(void);
extern void gencoordp(int nrot, double ang_Rot);
extern float ran2(long *idum);
extern double get_EDebye(int atomfix1,int atomfix2);
extern double get_EDebyep(int atomfix1,int atomfix2);
extern double get_EDebye_SGC(int atomfix1);
extern double get_EDebye_SGCp(int atomfix1);
extern double get_Wmech(void);
extern void write_final_data(void);
extern double get_Wmechp(void);
extern void generabondangle(void);
extern double get_Estr(void);
extern double get_Estrp(void);
extern double get_Ebdn(void);
extern double get_Ebdnp(void);
extern void gencoordpstr(int nstr);
extern void gencoordpbnd(int nbnd);
extern void actualitza_strbnd(void);
extern void compute_averages(void);
extern void send_info(void);


#endif

