#include "lib.h"
#include "global.h"


void lectura_input(){
char cadena[2000];
int period, i, j, k, cont, Nex,  Nend, assignment;
double Blenght[MAXATOMS],Bangle[MAXATOMS], Radius[MAXATOMS];
char fin[20] = "E";
ifstream file_input ("../Inputs/Global_input.data");

if (file_input)
   { }  // FILE OK
else
   { cout << "Fatal error: Please check ""Global_input.data"" input file"; exit(1);  } // CORRUPTED FILE

// Lecture of the header of the Input

for (i=1; i<= 4; i++){file_input.getline(cadena,500);}


//Number of atoms in the monomeric chain

file_input >> cadena;
file_input >>  Natom;
file_input.getline(cadena,500);

//Period

file_input >> cadena;
file_input >>  period;
file_input.getline(cadena,500);

//Chain structure

file_input >> cadena;

if (period == 0){

   for (i=1; i <= Natom; i++){file_input >> chain[i];}

}
else{
   for (i=1; i <= period; i++){file_input >> chain[i];}
   i=0;
   for (j=period+1; j <= Natom; j++) {
         i++;  
         if(i>period){i=1;}
         chain[j]=chain[i];      
   }

}

file_input.getline(cadena,500);

// LECTURE OF THE PKA OF THE POLYMERIC CHAIN (FOR ATOM)




// Lecture of the pKa of the initial atoms of the chain. If the chain is not periodic, the lecture is skipped.

if ( period != 0) {

     file_input >> cadena;
     file_input >>  Nex;
     file_input.getline(cadena,500);

}
else{file_input.getline(cadena,500);}


if ( period != 0) {
     
     file_input >> cadena;
     for (i=1; i <= Nex; i++){file_input >>  pKa[i];}
     file_input.getline(cadena,500);
}
else{file_input.getline(cadena,500);}

// Asignation of the pKa of the final atoms of the chain. If the chain is not periodic, the asignation is skipped.



if ( period != 0) {
  for (i=Natom-Nex+1; i <= Natom; i++){
    for( j=1 ; j <= Nex ; j++){
      if( chain[i] == chain[j]  && chain[i-1] == chain[j+1] && chain[i-2] == chain[j+2]  ){pKa[i]=pKa[j];}
    } 
  }
}



// Lecture of the pKa of the polymeric chain. If the chain is periodic, only the pKa of the periodic unit will be readen. the pKa are expected to continue the periodic serie where the initial pKa specifications finished.


// If the pKa of an atom is -100, it is considered non-ionizable 

if ( period != 0) {
          
     file_input >> cadena;
     for (i=1; i <= period; i++){file_input >>  pKa[Nex+i];}
     file_input.getline(cadena,500);
     cont=0;     
     for (i=Nex+period+1; i <= Natom-Nex; i++){
         cont++;
         pKa[i]=pKa[Nex+cont];
         if (cont == period){cont=0;}
     }

}
else{
     file_input >> cadena;
     for (i=1; i <= Natom; i++){file_input >>  pKa[i];}
     file_input.getline(cadena,500);
}






//CHARGE LECTURE

file_input >> cadena;




if (period == 0){

   for (i=1; i <= Natom; i++){file_input >> H_charge[i];}

}
else{



   for (i=1; i <= period; i++){file_input >> H_charge[i];}
   i=0;


   for (j=period+1; j <= Natom; j++) {
         i++;
         if(i>period){i=1;}
         H_charge[j]=H_charge[i];

   }

}


file_input.getline(cadena,500);

file_input >> cadena;
if (period == 0){

   for (i=1; i <= Natom; i++){file_input >> L_charge[i];}

}
else{
   for (i=1; i <= period; i++){file_input >> L_charge[i];}
   i=0;
   for (j=period+1; j <= Natom; j++) {
         i++;
         if(i>period){i=1;}
         L_charge[j]=L_charge[i];
   }

}

file_input.getline(cadena,500);


// Lecture of the number of monomers between two chargeable residues

file_input >> cadena;
file_input >> N_neigh ;
file_input.getline(cadena,500);



// Lecture of the seed of the random number generator

file_input >> cadena;
file_input >>  idum_ref;
file_input.getline(cadena,500);


// Lecture of the number of configurations

file_input >> cadena;
file_input >>  conffin;
file_input.getline(cadena,500);

// Lecture of the number of simulations

file_input >> cadena;
file_input >>  Nsim;
file_input.getline(cadena,500);

// Lecture of the number of steps until calculating the averages

file_input >> cadena;
file_input >>  intervalAVG;
file_input.getline(cadena,500);

// Lecture of the number of steps until writing the outputs

file_input >> cadena;
file_input >>  intervalOUT;
file_input.getline(cadena,500);

// Lecture of the Temperature

file_input >> cadena;
file_input >> temp;
file_input.getline(cadena,500);


// Lecture of the pH

file_input >> cadena;
file_input >>  pH;
file_input.getline(cadena,500);

// Lecture of the relative dielectrical constant

file_input >> cadena;
file_input >>  DIEBULK;
file_input.getline(cadena,500);

// Lecture of the ionic strengh

file_input >> cadena;
file_input >>  I;
file_input.getline(cadena,500);

// Lecture of ProbConf

file_input >> cadena;
file_input >> ProbConf; 
file_input.getline(cadena,500);


// Lecture of Probtest

file_input >> cadena;
file_input >>  probTest;
file_input.getline(cadena,500);

// Lecture of the long range interaction posibility

file_input >> cadena;
file_input >>  Debye;
file_input.getline(cadena,500);


// Lecture of the mechanical Force applied in x axis in pN

file_input >> cadena;
file_input >>  F;
file_input.getline(cadena,500);

// Lecture of the harmonic constant of bond strething in kcal/mol A-2

file_input >> cadena;
file_input >>  kstr;
file_input.getline(cadena,500);

// Lecture of the harmonic constant of bond bending in kcal/mol deg-2

file_input >> cadena;
file_input >>  kbdn;
file_input.getline(cadena,500);



//  GLOBAL LECTURE END//

file_input.close();



// LECTURE OF THE BOUND LENGHTS

ifstream bound_input ("../Inputs/Bound_lengh.data");

if (bound_input)
   { }  // FILE OK
else
   { cout << "Fatal error: Please check ""Bound_lengh.data"" input file"; exit(1);  } // CORRUPTED FILE

// Lecture of the header of the Input




for (i=1; i<= 5; i++){bound_input.getline(cadena,500);}

cont=0;

while( !bound_input.eof()){

cont++;

bound_input >> labelA1[cont];
bound_input >>  labelA2[cont];
bound_input >>  Blenght[cont];


}

// Asignation of each bound lenght in a vector (Bounds)



for (i=1; i <= Natom-1 ; i++ ){

  for (j=1; j <= cont ; j++){


        if (labelA1[j] == chain[i] &&  labelA2[j] == chain[i+1]){Bounds[i]=Blenght[j];}
        if (labelA1[j] == chain[i+1] &&  labelA2[j] == chain[i]){Bounds[i]=Blenght[j];}

  }


}


//  BOUND LENGHT LECTURE END//

bound_input.close();

// Computation of the average bond lenght

for(i=1; i <= Natom-1 ; i++){prom_bound=prom_bound+Bounds[i];}

prom_bound=prom_bound/(Natom-1);



// LECTURE OF THE ATOMIC RADIUS

ifstream radius_input ("../Inputs/Radius_input.data");

if (radius_input)
   { }  // FILE OK
else
   { cout << "Fatal error: Please check ""Radius_input.data"" input file"; exit(1);  } // CORRUPTED FILE

// Lecture of the header of the Input

for (i=1; i<= 5; i++){radius_input.getline(cadena,500);}

cont=0;

while( !radius_input.eof()){

cont++;

radius_input >> labelA1[cont];
radius_input >>  Radius[cont];

}

// Asignation of each atomic radius in a vector (radiatom)


for (i=1; i <= Natom ; i++ ){

  for (j=1; j <= cont ; j++){
        if (labelA1[j] == chain[i] ){radiatom[i]=Radius[j];}
  }

}


//  ATOMIC RADIUS LECTURE END//

radius_input.close();




// LECTURE OF THE BOUND ANGLE

ifstream angle_input ("../Inputs/Bound_angle.data");

if (angle_input)
   { }  // FILE OK
else
   { cout << "Fatal error: Please check ""Bound_angle.data"" input file"; exit(1);  } // CORRUPTED FILE

// Lecture of the header of the Input

for (i=1; i<= 5; i++){angle_input.getline(cadena,500);}

cont=0;

while( !angle_input.eof()){

cont++;

angle_input >> labelA1[cont];
angle_input >>  labelA2[cont];
angle_input >>  labelA3[cont];
angle_input >>  Bangle[cont];


}

// Asignation of each bound angle  in a vector (Bangle)



for (i=1; i <= Natom-2 ; i++ ){

  for (j=1; j <= cont ; j++){


        if (labelA1[j] == chain[i] &&  labelA2[j] == chain[i+1] && labelA3[j] == chain[i+2]){Angles[i]=Bangle[j];}
        if (labelA1[j] == chain[i+2] &&  labelA2[j] == chain[i+1] && labelA3[j] == chain[i] ){Angles[i]=Bangle[j];}
  }

}



//  BOUND ANGLE LECTURE END//

angle_input.close();



// TRANSFER MATRIX LECTURE //

ifstream trans_input ("../Inputs/Transfer_Matrix.data");

for(i=1 ; i <= 8 ; i++){trans_input.getline(cadena,500);}


// Number of final chain transfer matrix

trans_input >> cadena;
trans_input >>  Nmat;
trans_input.getline(cadena,500);

// Number of atoms with final chain transfer matrix

trans_input >> cadena;
trans_input >>  Nend;
trans_input.getline(cadena,500);

// Number of atoms with final chain transfer matrix

trans_input >> cadena;
trans_input >>  Ntranex;
trans_input.getline(cadena,500);
trans_input.getline(cadena,500);


// First Matrix


if ( Nmat >= 1){
trans_input >> labelA1[1];
trans_input >>  labelA2[1];
trans_input >>  labelA3[1];
trans_input >>  labelA4[1];
trans_input.getline(cadena,500);


}
else{
trans_input >> labelA1[1];
trans_input >>  labelA2[1];
trans_input >>  labelA3[1];
trans_input.getline(cadena,500);
}


trans_input >> Conf1[1][1];
trans_input >> Conf1[1][2];
trans_input >> Conf1[1][3];
trans_input.getline(cadena,500);


trans_input >> Conf1[2][1];
trans_input >> Conf1[2][2];
trans_input >> Conf1[2][3];
trans_input.getline(cadena,500);

trans_input >> Conf1[3][1];
trans_input >> Conf1[3][2];
trans_input >> Conf1[3][3];
trans_input.getline(cadena,500);


trans_input >> u_taux[1];
trans_input >> u_gpaux[1];
trans_input >> u_gnaux[1];
trans_input.getline(cadena,500);



trans_input.getline(cadena,500);

// Second Matrix



if ( Nmat >= 2){
trans_input >> labelA1[2];
trans_input >>  labelA2[2];
trans_input >>  labelA3[2];
trans_input >>  labelA4[2];
trans_input.getline(cadena,500);


}
else{

trans_input >> labelA1[2];
trans_input >>  labelA2[2];
trans_input >>  labelA3[2];
trans_input.getline(cadena,500);
}

trans_input >> Conf2[1][1];
trans_input >> Conf2[1][2];
trans_input >> Conf2[1][3];
trans_input.getline(cadena,500);

trans_input >> Conf2[2][1];
trans_input >> Conf2[2][2];
trans_input >> Conf2[2][3];
trans_input.getline(cadena,500);

trans_input >> Conf2[3][1];
trans_input >> Conf2[3][2];
trans_input >> Conf2[3][3];
trans_input.getline(cadena,500);

trans_input >> u_taux[2];
trans_input >> u_gpaux[2];
trans_input >> u_gnaux[2];
trans_input.getline(cadena,500);

trans_input.getline(cadena,500);

// Third Matrix

if ( Nmat >= 3){
trans_input >> labelA1[3];
trans_input >>  labelA2[3];
trans_input >>  labelA3[3];
trans_input >>  labelA4[3];
trans_input.getline(cadena,500);


}
else{
trans_input >> labelA1[3];
trans_input >>  labelA2[3];
trans_input >>  labelA3[3];
trans_input.getline(cadena,500);
}

trans_input >> Conf3[1][1];
trans_input >> Conf3[1][2];
trans_input >> Conf3[1][3];
trans_input.getline(cadena,500);

trans_input >> Conf3[2][1];
trans_input >> Conf3[2][2];
trans_input >> Conf3[2][3];
trans_input.getline(cadena,500);

trans_input >> Conf3[3][1];
trans_input >> Conf3[3][2];
trans_input >> Conf3[3][3];
trans_input.getline(cadena,500);

trans_input >> u_taux[3];
trans_input >> u_gpaux[3];
trans_input >> u_gnaux[3];
trans_input.getline(cadena,500);

trans_input.getline(cadena,500);

// Fourth Matrix

if ( Nmat >= 4){
trans_input >> labelA1[4];
trans_input >>  labelA2[4];
trans_input >>  labelA3[4];
trans_input >>  labelA4[4];
trans_input.getline(cadena,500);


}
else{

trans_input >> labelA1[4];
trans_input >>  labelA2[4];
trans_input >>  labelA3[4];
trans_input.getline(cadena,500);
}
trans_input >> Conf4[1][1];
trans_input >> Conf4[1][2];
trans_input >> Conf4[1][3];
trans_input.getline(cadena,500);

trans_input >> Conf4[2][1];
trans_input >> Conf4[2][2];
trans_input >> Conf4[2][3];
trans_input.getline(cadena,500);

trans_input >> Conf4[3][1];
trans_input >> Conf4[3][2];
trans_input >> Conf4[3][3];
trans_input.getline(cadena,500);

trans_input >> u_taux[4];
trans_input >> u_gpaux[4];
trans_input >> u_gnaux[4];
trans_input.getline(cadena,500);

trans_input.getline(cadena,500);

// Fifth Matrix

if ( Nmat >= 5){
trans_input >> labelA1[5];
trans_input >>  labelA2[5];
trans_input >>  labelA3[5];
trans_input >>  labelA4[5];
trans_input.getline(cadena,500);


}
else{

trans_input >> labelA1[5];
trans_input >>  labelA2[5];
trans_input >>  labelA3[5];
trans_input.getline(cadena,500);
}

trans_input >> Conf5[1][1];
trans_input >> Conf5[1][2];
trans_input >> Conf5[1][3];
trans_input.getline(cadena,500);

trans_input >> Conf5[2][1];
trans_input >> Conf5[2][2];
trans_input >> Conf5[2][3];
trans_input.getline(cadena,500);

trans_input >> Conf5[3][1];
trans_input >> Conf5[3][2];
trans_input >> Conf5[3][3];
trans_input.getline(cadena,500);

trans_input >> u_taux[5];
trans_input >> u_gpaux[5];
trans_input >> u_gnaux[5];
trans_input.getline(cadena,500);

trans_input.getline(cadena,500);

// Sixth Matrix

if ( Nmat >= 6){
trans_input >> labelA1[6];
trans_input >>  labelA2[6];
trans_input >>  labelA3[6];
trans_input >>  labelA4[6];
trans_input.getline(cadena,500);


}
else{

trans_input >> labelA1[6];
trans_input >>  labelA2[6];
trans_input >>  labelA3[6];
trans_input.getline(cadena,500);
}

trans_input >> Conf6[1][1];
trans_input >> Conf6[1][2];
trans_input >> Conf6[1][3];
trans_input.getline(cadena,500);

trans_input >> Conf6[2][1];
trans_input >> Conf6[2][2];
trans_input >> Conf6[2][3];
trans_input.getline(cadena,500);

trans_input >> Conf6[3][1];
trans_input >> Conf6[3][2];
trans_input >> Conf6[3][3];
trans_input.getline(cadena,500);

trans_input >> u_taux[6];
trans_input >> u_gpaux[6];
trans_input >> u_gnaux[6];
trans_input.getline(cadena,500);

trans_input.getline(cadena,500);

// Seventh Matrix

if ( Nmat >= 7){
trans_input >> labelA1[7];
trans_input >>  labelA2[7];
trans_input >>  labelA3[7];
trans_input >>  labelA4[7];
trans_input.getline(cadena,500);


}
else{

trans_input >> labelA1[7];
trans_input >>  labelA2[7];
trans_input >>  labelA3[7];
trans_input.getline(cadena,500);
}

trans_input >> Conf7[1][1];
trans_input >> Conf7[1][2];
trans_input >> Conf7[1][3];
trans_input.getline(cadena,500);

trans_input >> Conf7[2][1];
trans_input >> Conf7[2][2];
trans_input >> Conf7[2][3];
trans_input.getline(cadena,500);

trans_input >> Conf7[3][1];
trans_input >> Conf7[3][2];
trans_input >> Conf7[3][3];
trans_input.getline(cadena,500);

trans_input >> u_taux[7];
trans_input >> u_gpaux[7];
trans_input >> u_gnaux[7];
trans_input.getline(cadena,500);

trans_input.getline(cadena,500);

// Eight Matrix

if ( Nmat >= 8){
trans_input >> labelA1[8];
trans_input >>  labelA2[8];
trans_input >>  labelA3[8];
trans_input >>  labelA4[8];
trans_input.getline(cadena,500);


}
else{

trans_input >> labelA1[8];
trans_input >>  labelA2[8];
trans_input >>  labelA3[8];
trans_input.getline(cadena,500);
}

trans_input >> Conf8[1][1];
trans_input >> Conf8[1][2];
trans_input >> Conf8[1][3];
trans_input.getline(cadena,500);

trans_input >> Conf8[2][1];
trans_input >> Conf8[2][2];
trans_input >> Conf8[2][3];
trans_input.getline(cadena,500);

trans_input >> Conf8[3][1];
trans_input >> Conf8[3][2];
trans_input >> Conf8[3][3];
trans_input.getline(cadena,500);

trans_input >> u_taux[8];
trans_input >> u_gpaux[8];
trans_input >> u_gnaux[8];
trans_input.getline(cadena,500);

trans_input.getline(cadena,500);

// Nineth Matrix

if ( Nmat >= 9){
trans_input >> labelA1[9];
trans_input >>  labelA2[9];
trans_input >>  labelA3[9];
trans_input >>  labelA4[9];
trans_input.getline(cadena,500);


}
else{

trans_input >> labelA1[9];
trans_input >>  labelA2[9];
trans_input >>  labelA3[9];
trans_input.getline(cadena,500);
}

trans_input >> Conf9[1][1];
trans_input >> Conf9[1][2];
trans_input >> Conf9[1][3];
trans_input.getline(cadena,500);

trans_input >> Conf9[2][1];
trans_input >> Conf9[2][2];
trans_input >> Conf9[2][3];
trans_input.getline(cadena,500);

trans_input >> Conf9[3][1];
trans_input >> Conf9[3][2];
trans_input >> Conf9[3][3];
trans_input.getline(cadena,500);

trans_input >> u_taux[9];
trans_input >> u_gpaux[9];
trans_input >> u_gnaux[9];
trans_input.getline(cadena,500);

trans_input.getline(cadena,500);


// Tenth Matrix

if ( Nmat >= 10){
trans_input >> labelA1[10];
trans_input >>  labelA2[10];
trans_input >>  labelA3[10];
trans_input >>  labelA4[10];
trans_input.getline(cadena,500);


}
else{

trans_input >> labelA1[10];
trans_input >>  labelA2[10];
trans_input >>  labelA3[10];
trans_input.getline(cadena,500);
}

trans_input >> Conf10[1][1];
trans_input >> Conf10[1][2];
trans_input >> Conf10[1][3];
trans_input.getline(cadena,500);

trans_input >> Conf10[2][1];
trans_input >> Conf10[2][2];
trans_input >> Conf10[2][3];
trans_input.getline(cadena,500);

trans_input >> Conf10[3][1];
trans_input >> Conf10[3][2];
trans_input >> Conf10[3][3];
trans_input.getline(cadena,500);



trans_input >> u_taux[10];
trans_input >> u_gpaux[10];
trans_input >> u_gnaux[10];


trans_input.getline(cadena,500);



// Transfer Matrix lecture end

trans_input.close();






cont=0;




for (i=2; i <= Natom-2 ; i++ ){

  assignment = 1;
  j=0;

  while ( assignment !=  0){

          j++;
          if (labelA1[j] == fin[0] && labelA2[j] == chain[i] &&  labelA3[j] == chain[i+1] && labelA4[j] == chain[i+2] && assignment != 0 && cont < Nend-1){Mat_ind[i]= j; assignment = 0; cont++; }
          if (labelA1[j] == fin[0] && labelA2[j] == chain[i] &&  labelA3[j] == chain[i+1] &&  labelA4[j] == chain[i+2] && assignment !=  0 && i > Natom - Nend - 1){Mat_ind[i]= j; assignment = 0;}
          if(labelA1[j] == chain[i] && labelA2[j] == chain[i+1] &&  labelA3[j] == chain[i+2] && assignment != 0) {Mat_ind[i]= j; assignment = 0;}
          if (j > MAXTRANSFER){cout << "WARNING TRANSFER MATRIX " << i << " NOT ASIGNED" << endl; assignment = 0;}
          }
  }


// First and last bound have a different assigment


for ( j=1; j <= MAXTRANSFER ; j++){

if (labelA1[j] == fin[0] && labelA2[j] == chain[1] &&  labelA3[j] == chain[2] && labelA4[j] == chain[3] ){Mat_ind[1]= j;}
if (labelA1[j] == fin[0] && labelA2[j] == chain[Natom] &&  labelA3[j] == chain[Natom-1] && labelA4[j] == chain[Natom-2] ){Mat_ind[Natom-1]= j;}



}





for (i=1; i <= Natom-1 ; i++ ){


u_t[i] = u_taux[Mat_ind[i]];
u_gp[i] = u_gnaux[Mat_ind[i]];
u_gn[i]=u_gnaux[Mat_ind[i]];


}


}
