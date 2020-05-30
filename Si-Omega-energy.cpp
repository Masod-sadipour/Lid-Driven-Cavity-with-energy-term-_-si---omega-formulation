 
// Solving Ψ-Ω Formulation in a Lid-Driven Cavity with energy equatoion
// April 2018.

// necessary libraries

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Variables and discritization parameters

int main() {
int imax=20, jmax=20 ,nmax=500, i, j, n,coef=500;
float ut=1.,	l=1.,	h=1.,	Re=1,	RI=0.1,	Pr=1.,	t=1.3, nu=0.001, Gr=Re*Re*Pr , utt=Re*nu/l;
float re=(ut*l)/nu, dx=l/(imax-1),	dy=h/(jmax-1),  betta=dx/dy;
float si[imax+1][jmax+1][nmax+1], T[imax+2][jmax+2][nmax+1], omega[imax+1][jmax+1][nmax+1];
float u[imax+1][jmax+1], v[imax+1][jmax+1];
float k,z;
float cxp, cxo, cxm, cyp, cyo, cym;
float g, g1;
float p1[imax+1][jmax+1];
float errorsi=0., erroromega=0., errorT=0.;

// initial conditions for "si", "omega", and "Temperature" 

for(i=1;i<=imax;i++){
for(j=1;j<=jmax;j++){
si[i][j][1]=0.1;
T[i][j][1]=0.5;
omega[i][j][1]=(-j/jmax)*((2*ut)/dy);
}
}


for(n=1;n<=nmax;n++){
// boundary condition on y bond

for (j=1;j<=jmax;j++){
si[1][j][n+1]=0;
si[imax][j][n+1]=0;
u[1][j]=0;
u[imax][j]=0;
T[1][j][n+1]=0;        //Tempreture B.C. in the left side of the cavity T=0
T[imax][j][n+1]=1;    // Tempreture B.C. in the right side of the cavity T=1
v[1][j]=0;
v[imax][j]=0;
}
// boundary condition on x boundaries

for(i=1;i<=imax;i++){
u[i][1]=0;			// down
u[i][jmax]=ut;                 // up
v[i][1]=0;
v[i][jmax]=0;
si[i][1][n+1]=0;
si[i][jmax][n+1]=0;
}

for(i=2;i<=imax-1;i++){
omega[i][1][n+1]=(-2*si[i][2][n])/(dy*dy);
omega[i][jmax][n+1]=-((2*si[i][jmax-1][n]+2*ut*dy)/(dy*dy));
}
// boundary condition on y boundaries

for (j=2;j<=jmax-1;j++){
omega[imax][j][n+1]=(-2*si[imax-1][j][n])/(dx*dx);
omega[1][j][n+1]=(-2*si[2][j][n])/(dx*dx);
}

for (j=2;j<=jmax-1;j++){
for(i=2;i<=imax-1;i++){
u[i][j]=(si[i][j+1][n]-si[i][j-1][n])/(2*dy);
v[i][j]=-((si[i+1][j][n]-si[i-1][j][n])/(2*dx));
k=abs((u[i][j]*dx)/nu);
z=abs((v[i][j]*dy)/nu);
}
}

// ***************************************   loop   ***************************************

for (j=2;j<=jmax-1;j++){
for(i=2;i<=imax-1;i++){
if (abs((u[i][j]*utt*dx)/nu)<=2){
cxp=0.5;
cxo=0;
cxm=-0.5;
}
else if (abs((u[i][j]*utt*dx)/nu)>2){
if (u[i][j]>0){
cxp=0;
cxo=1;
cxm=-1;
}
else if (u[i][j]<0){
cxp=1;
cxo=-1;
cxm=0;
}
}

if (abs((v[i][j]*utt*dy)/nu)<=2){
cyp=0.5;
cyo=0;
cym=-0.5;
}
else if (abs((v[i][j]*utt*dy)/nu)>2){
if (v[i][j]>0){
cyp=0;
cyo=1;
cym=-1;
}
else if (v[i][j]<0){
cyp=1;
cyo=-1;
cym=0;
}
}
g=((2*(1+(betta*betta)))+((u[i][j]*dx*cxo)*Re)+((betta*v[i][j]*dx*cyo)*Re));
omega[i][j][n+1]=omega[i][j][n]+((t/(g))*(((1-((cxp*dx*u[i][j])*Re))*omega[i+1][j][n])+((1-((u[i][j]*dx*cxm)*Re))*omega[i-1][j][n+1])+(((betta*betta)-((betta*v[i][j]*dx*cyp)*Re))*omega[i][j+1][n])+(((betta*betta)-((betta*dx*v[i][j]*cym)*Re))*omega[i][j-1][n+1])-(g*omega[i][j][n])))+(RI)*Re*(dx*dx/4)*(T[i+1][j][n]-T[i-1][j][n])/(2*dx);
}
}
for (j=2;j<=jmax-1;j++){
for(i=2;i<=imax-1;i++){
si[i][j][n+1]=si[i][j][n]+(t/(2*(1+(betta*betta))))*(si[i+1][j][n]+si[i-1][j][n+1]+(betta*betta)*si[i][j+1][n]+(betta*betta)*si[i][j-1][n+1]-2*(1+(betta*betta))*si[i][j][n]+(dx*dx)*(omega[i][j][n+1]));
}
}

for(i=2;i<=imax-1;i++){
	for (j=1;j<=jmax;j++){
		T[i][0][n+1]=T[i][2][n];
		T[i][jmax+1][n]=T[i][jmax-1][n];
	}
}
for(i=2;i<=imax-1;i++){
	for (j=1;j<=jmax;j++){
		g1=((2*(1+(betta*betta)))+((u[i][j]*dx*cxo)*Re*Pr)+((betta*v[i][j]*dx*cyo)*Re*Pr));
		T[i][j][n+1]=T[i][j][n]+((t/(g1))*(((1-((cxp*dx*u[i][j])*Re*Pr))*T[i+1][j][n])+((1-((u[i][j]*dx*cxm)*Re*Pr))*T[i-1][j][n+1])+(((betta*betta)-((betta*v[i][j]*dx*cyp)*Re*Pr))*T[i][j+1][n])+(((betta*betta)-((betta*dx*v[i][j]*cym)*Re*Pr))*T[i][j-1][n+1])-(g1*T[i][j][n])));
	}
}

for(i=1;i<=imax;i++){
for(j=1;j<=jmax;j++){
	errorsi=errorsi+fabs(si[i][j][n+1]-si[i][j][n]);
	erroromega=erroromega+fabs(omega[i][j][n+1]-omega[i][j][n]);
	errorT=errorT+fabs(T[i][j][n+1]-T[i][j][n]);
}}
if(errorsi<=10e-6 && errorT<=10e-6 && erroromega<=10e-6 ){
	break;
}

}

// ***************************************   end of the loop   ***************************************


//****************************************      pressure    ******************************************


	p1[1][1]=1;
	for(i=1;i<=imax-1;i++){
		p1[i+1][1]=p1[i][1]+dx*(-u[i][4]+4*u[i][3]-5*u[i][2]+2*u[i][1])/(dy*dy*Re);
		}
	for(j=1;j<=jmax-1;j++){
	
		p1[1][j+1]=p1[1][j]+dx*(-u[4][j]+4*u[3][j]-5*u[2][j]+2*u[1][j])/(dy*dy*Re);}
		
	for(j=1;j<=jmax-1;j++){
	
		p1[imax][j+1]=p1[imax][j]+dx*(-u[imax-3][j]+4*u[imax-2][j]-5*u[imax-1][j]+2*u[imax][j])/(dy*dy*Re);}

	for(i=1;i<=imax-1;i++){
	
		p1[i+1][jmax]=p1[i][jmax-1]+dx*(-u[i][jmax-3]+4*u[i][jmax-2]-5*u[i][jmax-1]+2*u[i][jmax])/(dy*dy*Re);}

	for(i=2;i<=imax-1;i++){
		for(j=2;j<=jmax-1;j++){
		
			p1[i][j]=p1[i-1][j]+(v[i+1][j]-2*v[i][j]+v[i-1][j])/(Re*dx*betta)+(v[i][j-1]-2*v[i][j]+v[i][j+1])/(Re*dy)-u[i][j]*(v[i+1][j]-v[i][j-1])/(2*betta)-v[i][j]*(v[i][j+1]-v[i][j-1])/2;}}


// ******************************************** tecplot  **********************************************
float x[imax+1][jmax+1], y[imax+1][jmax+1];

for(j=1;j<=jmax;j++){
for(i=1;i<=imax;i++){
x[i][j]=(j-1)*dx;
y[i][j]=(i-1)*dy;
}
}

FILE *fid;
fid = fopen("UVT.plt", "w+t");

fprintf(fid, "TITILE=" "\n");
fprintf(fid, "VARIABLES=X,Y,U,V,T\n");
fprintf(fid, "ZONE\n");
fprintf(fid, "I= %d J= %d\n",imax,jmax);
fprintf(fid, "F=POINT\n");

	for (i=1;i<=imax;i++){
		for (j=1;j<=jmax;j++){
		fprintf(fid, "%f\t%f\t%f\t%f\t%f\t\n ",y[i][j],x[i][j],u[i][j],v[i][j],T[i][j][n-2]);
	}
}
fclose(fid);

	return 0;
}
