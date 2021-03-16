// Module with subroutines to perform2a certain rotation (in degrees) in each axis given the cartesian coordinates.

#include "lib.h"
#include "global.h"

/******************************************************/

// Eq.20, pag. 16, Statistical Mechanics of Chain Molecules, Paul J. Flory, 1989

double r2_FreelyRotatingChain(int nbonds,double distancebond,double angle3sites)
{
double alfa,den,val;

alfa=cos((180.-angle3sites)*PI/180.);
den=(1-alfa);
val=(1.+alfa)/den-(2*alfa/(double)nbonds)*(1-pow(alfa,nbonds))/(den*den);
val=val*nbonds*distancebond*distancebond;

return(val);
}

/**************************************************************************/

// Eq.21, pag. 17, Statistical Mechanics of Chain Molecules, Paul J. Flory, 1989

double s2_FreelyRotatingChain(int nbonds,double distancebond,double angle3sites)
{
double alfa,l2,n,t1,t2,sum,val;
int j;
alfa=cos((180.-angle3sites)*PI/180.);
l2=distancebond*distancebond;
n=nbonds;

t1=l2*n*(n+2)*(1+alfa)/6/(n+1)/(1-alfa);
t2=2*l2*alfa/(n+1)/(n+1)/(1-alfa)/(1-alfa);
sum=0;
for(j=1;j<=n;j++)
 {sum+= (double)j-(alfa-pow(alfa,j+1))/(1-alfa);}
val=t1-t2*sum;
return(val);
}

