#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>  /* clock_t, clock, CLOCKS_PER_SEC */
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include "Eigen/IterativeLinearSolvers"
#include <vector>
#include "ioNDA.h"
#include "defNDA.h"
using namespace Eigen;
using namespace std;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

/*
  2D Transport Solution With Step Characteristics and Non-linear Diffusion Acceleration
  
  Programer : Luke Cornejo
  Version   : 1.0
  Date      : 2-19-14
  
                         Changes
  ******************************************************

  
  To complile
  Windows
c++ -I .\eigen SC_NDA1.0.cpp ioNDA.cpp defNDA.cpp -o SC_NDA1.0.exe -std=c++0x
  Linux
c++ -I .\eigen SC_NDA1.0.cpp ioNDA.cpp defNDA.cpp -o SC_NDA1.0.exe
  
  hpc
c++ -I .\eigen -O3 SC_NDA1.0.cpp ioNDA.cpp defNDA.cpp -o SC_NDA1.0b.exe
-O
-O2
-O3 better
-Os
-O3 -ffast-math


icpc -I .\eigen -O3 SC_NDA1.0.cpp ioNDA.cpp defNDA.cpp -o SC_NDA1.0a.exe
-O1
-O2
-O3 best
-xO
-fast

  Boundary type options
  1: incoming according to side
  2: incoming according to angle
  3: reflective on all sides
  4: reflective on Left and Bottom, incoming on Right and Top according to side
  5: reflective on Right and Top, incoming on Left and Bottom according to side
  
  BC input order
  according to side : Left, Bottom, Right, Top
  according to angle: quad 1, quad 2, quad 3, quad 4
  
  Quadratures
  S4, S6, S8, S12, S16
  Q20, Q20
  

*/

void Iterations();    // perform source iterations on problem
void initialAngleSweep();   // find initial high-order solution using constant angular flux
void angleSweep();          // determine boundary conditions and how to sequence quadrant solutions
void quad1();               // sweep through angles and cells in first quadrant
void quad2();               // sweep through angles and cells in second quadrant
void quad3();               // sweep through angles and cells in third quadrant
void quad4();               // sweep through angles and cells in fourth quadrant
void cellSolution(double, double, double, double, double, double, double, double, double, double&, double&, double& );
void NDAsolution();          // NDA solution function
void output(char[]);          // output code
// data output functions

const double pi=3.141592654;
int nx, ny, sn, n_iterations;
int maxiter_lo=1000;
bool o_angular=false; // option variables
double epsilon_si=1e-5, epsilon_lo=1e-10; // default tolerances
int N=8;                     // quadrature
double *mu, *eta, *xi, *w;   // quadrature
double *x, *y, *hx, *hy, *xe, *ye; // grid arrays
double **sigmaT, **sigmaS, **sigmaF, **nuF, **s_ext; // material arrays
int **material;
double **psi , **psi_x , **psi_y; // angular flux
double **phi , **j_x , **j_y; // NDA scalar flux and current solution
double *phiB_L, *phiB_R, *phiB_B, *phiB_T; // edge scalar flux from NDA
double **phiT, **phi_xT, **phi_yT, **j_xT, **j_yT; // scalar flux and current from transport
double **phiL; // scalar flux from last iteration
int    kbc;
double *bcL, *bcR, *bcB, *bcT;
double ***psiL, ***psiR, ***psiB, ***psiT;
double *FL, *FR, *FB, *FT;
double **D, **D_x, **D_y; // Diffusion coefficients
double **D_xT, **D_yT, **D_xP, **D_xN, **D_yP, **D_yN; // Tilde, Positive and Negative
double max_res;
// iteration data
vector<double> rho_si (1,0.5), err_lo, dt_ho, dt_lo, dt_pc;
vector<int> num_lo;
ofstream outfile;
ofstream datfile;
ofstream temfile;
clock_t t;
clock_t t_lo;
clock_t t_ho;
clock_t t_pc;


//======================================================================================//
//++ Read Input ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void input(char fn[]) //
{
	char inpf[25]={""}, outf[25]={""};
	const int nl=250;
	char namel[nl], inl[nl];
	int *regnum, *matnum;
	double *xn, *yn, *xp, *yp, *st, *ss, *sf, *nf, *se, *xg, *yg;
	double *xbc_e, *ybc_e, *bc_l, *bc_r, *bc_b, *bc_t, centre; // bc temp
	int i, j, k, p, nmat, nreg, ngx=0, ngy=0,  *nxt, *nyt, nbc, xbc_n, ybc_n;
	
	// input file name
	strcat (inpf,fn); strcat (inpf,".inp"); // input file name
	cout << inpf << endl;
	
	// open file +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	ifstream infile;
	infile.open (inpf); // open input file
	
	// read in data
	infile.getline(namel, nl); // read in name of the input
	
	// read x grid data
	infile.getline(inl, nl); // x grid edges
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if ( ngx==0 ) {
				ngx=iparse(inl,i,nl); // find number of mesh regions
				nxt=new int[ngx];
				xg=new double[ngx+1];
				p=0;
			}
			else {
				xg[p]=dparse(inl,i,nl); p++;
			}
		}
		if ( inl[i]==';' ) break;
	}
	if ( p!=ngx+1 )	cout<<"xgrid input error\n";
	
	infile.getline(inl, nl); // # cells in x grid
	p=0;
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			nxt[p]=iparse(inl,i,nl); p++;
		}
		if ( inl[i]==';' )	break;
	}
	if ( p!=ngx ) cout<<"xgrid input error\n";
	
	// read y grid data
	infile.getline(inl, nl); // read in grid data
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if ( ngy==0 ) {
				ngy=iparse(inl,i,nl);
				nyt=new int[ngy];
				yg=new double[ngy+1];
				p=0;
			}
			else {
				yg[p]=dparse(inl,i,nl); p++; // y grid zones
			}
		}
		if ( inl[i]==';' )	break;
	}
	if ( p!=ngy+1 )	cout<<"ygrid input error\n";
	
	infile.getline(inl, nl); // read in grid data
	p=0;
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			nyt[p]=iparse(inl,i,nl); p++; // # cells in y grid zone
		}
		if ( inl[i]==';' ) break;
	}
	if ( p!=ngy ) cout<<"y grid input error\n";
	
	
	// read in boundary conditions
	// Get Type of BC
	infile.getline(inl, nl);
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) kbc=iparse(inl,i,nl); // kind of BC
		if ( inl[i]==';' ) break;
	}
	
	p=0; // Get Bottom and Top BC
	infile.getline(inl, nl);
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if ( p==0 ) {
				nbc=iparse(inl,i,nl); // kind of BC
				xbc_n=nbc;
				xbc_e=new double[nbc+1];
				bc_b=new double[nbc];
				bc_t=new double[nbc];
			}
			else if ( p<nbc ) 	xbc_e[p]=dparse(inl,i,nl); // Left   BC or Quadrant 1
			else if ( p<2*nbc ) bc_b[p-nbc]  =dparse(inl,i,nl); // Bottom BC or Quadrant 2
			else                bc_t[p-2*nbc]=dparse(inl,i,nl); // Right  BC or Quadrant 3
			p++;
		}
		if ( inl[i]==';' ) break;
	}
	
	p=0; // Get Left and Right BC
	infile.getline(inl, nl);
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if ( p==0 ) {
				nbc=iparse(inl,i,nl); // kind of BC
				ybc_n=nbc;
				ybc_e=new double[nbc+1];
				bc_l=new double[nbc];
				bc_r=new double[nbc];
			}
			else if ( p<nbc )   ybc_e[p]=dparse(inl,i,nl); // Left   BC or Quadrant 1
			else if ( p<2*nbc ) bc_l[p-nbc]  =dparse(inl,i,nl); // Bottom BC or Quadrant 2
			else                bc_r[p-2*nbc]=dparse(inl,i,nl); // Right  BC or Quadrant 3
			p++;
		}
		if ( inl[i]==';' ) break;
	}
	
	
	infile.getline(inl, nl); // space
	
	infile.getline(inl, nl);  // number of materials
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			nmat=iparse(inl,i,nl); break;
		}
	}
	
	// cross sections
	matnum=new int[nmat];
	st=new double[nmat];
	ss=new double[nmat];
	sf=new double[nmat];
	nf=new double[nmat];
	se=new double[nmat];
	
	for (k=0; k<nmat; k++) {
		p=0;
		infile.getline(inl, nl);
		for (i=0; i<nl; i++) {
			if ( isdigit(inl[i]) ) {
				if ( p==0 )  matnum[k]=iparse(inl,i,nl); // Material #
				else if ( p==1 ) st[k]=dparse(inl,i,nl); // sigmaT total cross section
				else if ( p==2 ) ss[k]=dparse(inl,i,nl); // sigmaS scattering cross section
				else if ( p==3 ) sf[k]=dparse(inl,i,nl); // sigmaF fission cross section
				else if ( p==4 ) nf[k]=dparse(inl,i,nl); // nuF
				else if ( p==5 ) se[k]=4.0*pi*dparse(inl,i,nl); // external source
				else cout<<" Material Input Error \n";
				p++;
			}
			if ( inl[i]==';' ) break;
		}
	}
	
	infile.getline(inl, nl);  // number of material regions
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			nreg=iparse(inl,i,nl); break;
		}
	}
	
	// material edges
	regnum=new int[nreg];
	xn=new double[nreg];  xp=new double[nreg];
	yn=new double[nreg];  yp=new double[nreg];
	
	for (k=0; k<nreg; k++) {
		p=0;
		infile.getline(inl, nl);
		for (i=0; i<nl; i++) {
			if ( isdigit(inl[i]) ) {
				if ( p==0 )  regnum[k]=iparse(inl,i,nl); // Material Number in Region
				else if ( p==1 ) xn[k]=dparse(inl,i,nl); // Left   Boundary of Region
				else if ( p==2 ) yn[k]=dparse(inl,i,nl); // Bottom Boundary of Region
				else if ( p==3 ) xp[k]=dparse(inl,i,nl); // Right  Boundary of Region
				else if ( p==4 ) yp[k]=dparse(inl,i,nl); // Top    Boundary of Region
				else cout<<" Region Input Error \n";
				p++;
			}
			if ( inl[i]==';' ) break;
		}
	}
	p=0;
	infile.getline(inl, nl); // Read additional data
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if ( p==0 )      epsilon_si=dparse(inl,i,nl); // get convergence criteria
			else if ( p==1 ) epsilon_lo=dparse(inl,i,nl); // get low-order convergence criteria 
			p++;
		}
		if ( inl[i]==';' ) break;
	}
	
	if ( epsilon_si<1e-12 ) {
		cout<<">> Error in SI convergence criteria input!\n"; epsilon_si=1e-5;
	}
	if ( epsilon_lo<1e-14 or epsilon_lo>epsilon_si) {
		cout<<">> Error in low-order convergence criteria input!\n"; epsilon_lo=1e-10;
	}
	p=0;
	infile.getline(inl, nl); // Read additional data
	for (i=0; i<nl; i++) {
		if ( isdigit(inl[i]) ) {
			if ( p==0 ) N=iparse(inl,i,nl);
			else if ( p==1) {
				j=iparse(inl,i,nl);
				if ( j==1 ) o_angular=true; 
			}
		}
		if ( inl[i]==';' ) break;
	}
	if ( N!=4 and N!=6 and N!=8 and N!=12 and N!=16 and N!=20 and N!=36 ) {
		cout<<">> Error in Quadrature input!\n"; N=36;
	}
	
	infile.close(); // close input file ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	
	// total grid
	nx=0; for (i=0; i<ngx; i++) nx+=nxt[i]; // find total number of cells in x grid
	ny=0; for (i=0; i<ngy; i++) ny+=nyt[i]; // find total number of cells in y grid
	
	x=new (nothrow) double[nx+1]; hx=new (nothrow) double[nx]; xe=new double[nx+1];
	y=new (nothrow) double[ny+1]; hy=new (nothrow) double[ny]; ye=new double[ny+1];
	
	x[0]=xg[0]; p=0;
	for (i=0; i<ngx; i++) {
		for (j=p; j<p+nxt[i]; j++) {
			hx[j]=(xg[i+1]-xg[i])/nxt[i];
			x[j+1]=x[j]+hx[j];
		}
		p+=nxt[i];
	}
	y[0]=yg[0]; p=0;
	for (i=0; i<ngy; i++) {
		for (j=p; j<p+nyt[i]; j++) {
		hy[j]=(yg[i+1]-yg[i])/nyt[i];
		y[j+1]=y[j]+hy[j];
		}
		p+=nyt[i];
	}
	xe[0]=0.5*hx[0]; xe[nx]=0.5*hx[nx-1];
	ye[0]=0.5*hy[0]; ye[ny]=0.5*hy[ny-1];
	for (i=1; i<nx; i++) xe[i]=0.5*(hx[i-1]+hx[i]);
	for (j=1; j<ny; j++) ye[j]=0.5*(hy[j-1]+hy[j]);
	
	// BC
	bcB=new double[nx]; // Initialize Bottom BC
	bcT=new double[nx]; // Initialize Top BC
	bcL=new double[ny]; // Initialize Left BC
	bcR=new double[ny]; // Initialize Right BC
	xbc_e[0]=x[0]; xbc_e[xbc_n]=x[nx];
	ybc_e[0]=y[0]; ybc_e[ybc_n]=y[ny];
	
	for (i=0; i<nx; i++) {
		centre=(x[i]+x[i+1])/2;
		for (p=0; p<xbc_n; p++) {
			if ( xbc_e[p]<centre and xbc_e[p+1]>centre ) {
				bcB[i]=bc_b[p];
				bcT[i]=bc_t[p];
			}
		}
	}
	for (j=0; j<ny; j++) {
		centre=(y[j]+y[j+1])/2;
		for (p=0; p<ybc_n; p++) {
			if ( ybc_e[p]<centre and ybc_e[p+1]>centre ) {
				bcL[j]=bc_l[p];
				bcR[j]=bc_r[p];
			}
		}
	}
	
	
	// cross section data
	sigmaT=new double*[nx];
	sigmaS=new double*[nx];
	sigmaF=new double*[nx];
	nuF   =new double*[nx];
	s_ext =new double*[nx];
	material =new int*[nx];
	for (i=0; i<nx; i++) {
		sigmaT[i]=new double[ny];
		sigmaS[i]=new double[ny];
		sigmaF[i]=new double[ny];
		nuF[i]   =new double[ny];
		s_ext[i] =new double[ny];
		material[i] =new int[ny];
		for (j=0; j<ny; j++) {
			k=matRegion(matnum, regnum, i, j, nmat, nreg, xn, yn, xp, yp);
			sigmaT[i][j]=st[k];
			sigmaS[i][j]=ss[k];
			sigmaF[i][j]=sf[k];
			nuF[i][j]   =nf[k];
			s_ext[i][j] =se[k];
			material[i][j]=matnum[k];
			if ( sigmaT[i][j]==0.0 ) sigmaT[i][j]=1e-23;
		}
	}
	
	strcat (outf,fn);
	strcat (outf,".out");
	cout<<outf<<endl;
	outfile.open(outf); // open output file. closed in output function
	
	outfile<<"Output of File : "<<outf<<endl;
	outfile<<"2D Transport by Step Characteristics with NDA\n";
	outfile<<"Version: 1.0\n";
	outfile<<"Programer : Luke Cornejo\n";
	outfile<<"Case Name : "<<namel<<endl;
	// current date/time based on current system
	time_t now = time(0);
	// convert now to string form
	char* dt = ctime(&now);
	
	outfile<<"Program Ran On: "<<dt<<endl;
	
	outfile<<"+-----------------------------------------------+\n";
	outfile<<"Iteration Convergence Criteria: "<<epsilon_si<<endl;
	outfile<<"Low-order Tolerance: "<<epsilon_lo<<endl;
	outfile<<"+-----------------------------------------------+\n";
	
	outfile<<"\n -- Material Properties --\n";
	outfile<<" Material |  Total   |Scattering|  Fission |   nuF    | External \n";
	outfile<<"    #     |   XS     |    XS    |    XS    |          |  Source  \n";
	outfile.precision(6);
	for (k=0; k<nmat; k++) 
		outfile<<setw(10)<<matnum[k]<<"|"<<setw(10)<<st[k]<<"|"<<setw(10)<<ss[k]<<"|"<<setw(10)<<sf[k]<<"|"<<setw(10)<<nf[k]<<"|"<<setw(10)<<se[k]<<endl;
	
	outfile<<"\n -- Material Map -- \n";
	
	for (j=ny-1; j>=0; j--) {
		for (i=0; i<nx; i++) outfile<<setw(3)<<material[i][j];
		outfile<<endl;
	}
	
	
}
//======================================================================================//

//======================================================================================//
//++ Main program ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
int main (int argc, char* argv[])
{
	char fn[50], outf[50]={""};
	int i, j, m;
	double nu1, c;
	strcpy (fn,argv[1]);
	if ( argc!=2 ) cout << "usage: " << argv[0] << " input file name";
	else cout << fn << endl;
	
	input(fn); // get input data
	
	quadSet(); // find quadrature 
	
	initialize(); // initialize memory space for solution
	
	strcat (outf,fn);
	strcat (outf,".temp.csv");
	cout<<outf<<endl;
	temfile.open(outf); // open temporary file
	
	cout<<"++++++++++++++++++++++\n";
	t = clock();    // start timer
	// ------------------------------------------
	Iterations(); // Call Iterations
	// ------------------------------------------
	t = clock() -t; // stop timer
	cout<<"++++++++++++++++++++++\n";
	temfile.close(); // close temporary file
	
	output(fn);     // write out solutions
	
	return 0;
}
//======================================================================================//

//======================================================================================//
//++ Calculate the Diffusion Coefficients +++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void calcDiffusion()
{
	int i, j;
	// Calculate D^tilde and F factors
	
	// D_x
	for (i=1; i<nx; i++) {
		for (j=0; j<ny; j++) D_xT[i][j]=2*(j_xT[i][j]+D_x[i][j]*(phiT[i][j]-phiT[i-1][j])/xe[i])/(phiT[i][j]+phiT[i-1][j]); // X Grid
	}
	for ( j=0; j<ny; j++ ) { // Left and Right
		D_xT[0][j] =2*(j_xT[0][j] +D_x[0][j] *(phiT[0][j]   -phi_xT[0][j]) /xe[0]) /(phiT[0][j]   +phi_xT[0][j]); // Left Side
		FL[j]=j_xT[0][j]/phi_xT[0][j]; // Left F
		D_xT[nx][j]=2*(j_xT[nx][j]+D_x[nx][j]*(phi_xT[nx][j]-phiT[nx-1][j])/xe[nx])/(phi_xT[nx][j]+phiT[nx-1][j]); // Right Side
		FR[j]=j_xT[nx][j]/phi_xT[nx][j]; // Right F
	}
	// D_y
	for (i=0; i<nx; i++) {
		for (j=1; j<ny; j++) D_yT[i][j]=2*(j_yT[i][j]+D_y[i][j]*(phiT[i][j]-phiT[i][j-1])/ye[j])/(phiT[i][j]+phiT[i][j-1]); // Y Grid
	}
	for ( i=0; i<nx; i++ ) { // Bottom and Top
		D_yT[i][0] =2*(j_yT[i][0] +D_y[i][0] *(phiT[i][0]   -phi_yT[i][0]) /ye[0]) /(phiT[i][0]   +phi_yT[i][0]); // Bottom Side
		FB[i]=j_yT[i][0]/phi_yT[i][0]; // Bottom F
		D_yT[i][ny]=2*(j_yT[i][ny]+D_y[i][ny]*(phi_yT[i][ny]-phiT[i][ny-1])/ye[ny])/(phi_yT[i][ny]+phiT[i][ny-1]); // Top
		FT[i]=j_yT[i][ny]/phi_yT[i][ny]; // Top F
	}
	// D_xP
	for (i=0; i<nx+1; i++) {
		for (j=0; j<ny; j++) D_xP[i][j]=D_x[i][j]-0.5*D_xT[i][j]*xe[i]; // X Grid
	}
	// D_yP
	for (i=0; i<nx; i++) {
		for (j=0; j<ny+1; j++) D_yP[i][j]=D_y[i][j]-0.5*D_yT[i][j]*ye[j]; // Y Grid
	}
	// D_xN
	for (i=0; i<nx+1; i++) {
		for (j=0; j<ny; j++) D_xN[i][j]=D_x[i][j]+0.5*D_xT[i][j]*xe[i]; // X Grid
	}
	// D_yN
	for (i=0; i<nx; i++) {
		for (j=0; j<ny+1; j++) D_yN[i][j]=D_y[i][j]+0.5*D_yT[i][j]*ye[j]; // Y Grid
	}
	
	// In case of reflective BC set edge values to zero
	switch (kbc) {
		case 3: // BC type 3 Reflective on all sides 
			for (i=0; i<nx; i++) {
				FB[i]=0;	FT[i]=0; // Bottom and Top
			}
			for (j=0; j<ny; j++) {
				FL[j]=0;	FR[j]=0; // Left and Right
			}
			break;
		case 4: // BC type 4 Reflective on Bottom and Left
			for (i=0; i<nx; i++) FB[i]=0;
			for (j=0; j<ny; j++) FL[j]=0;
			break;
		case 5: // BC type 5 Reflective on Top and Right
			for (i=0; i<nx; i++) FT[i]=0;
			for (j=0; j<ny; j++) FR[j]=0;
			break;
		default:
			break;
	}
}
//======================================================================================//

//======================================================================================//
//++ Iterate to converge on solution +++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void Iterations()
{
	int i, j, ii, jj, m, s;
	double norm_si, norm_siL, phi_sum, res, psi_const=1.0;
	
	norm_si=4.0*pi;
	norm_siL=norm_si;
	s=0;
	// begin iterations
	while ( norm_si>epsilon_si*(1/rho_si[s]-1) or s==0 ) { //========================================================================
		cout<<"Iteration # "<<s<<endl;
		// set previous iteration to phiLast
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) {
				phiL[i][j]=phiT[i][j];
			}
		}
		
		norm_siL=norm_si; // set previous norm to normLast
		
		temfile<<"Iteration # "<<s<<endl;
		t_ho=clock(); // start high-order timer
		if ( s==0 ) initialAngleSweep(); // calculate initial transport solution from constant angular flux <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		else               angleSweep();        // perform sweep through angles <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		t_ho=clock()-t_ho; // stop high-order timer
		dt_ho.push_back(((double)t_ho)/CLOCKS_PER_SEC); // add high-order solution time to vector
		
		// Calculate D^tilde and F factors
		calcDiffusion(); // Calculate diffusion factors and Boundary conditions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		
		// Calculate Transport Residuals
		max_res=0;
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) {
				res=abs((j_xT[i+1][j]-j_xT[i][j])*hy[j]+(j_yT[i][j+1]-j_yT[i][j])*hx[i]+
					(sigmaT[i][j]*phiT[i][j]-(sigmaS[i][j]+nuF[i][j]*sigmaF[i][j])*phi[i][j]-s_ext[i][j])*hx[i]*hy[j]);
				if ( res>max_res ) {
					max_res=res; ii=i; jj=j;
				}
			}
		}
		t_lo=clock(); // start low-order timer
		NDAsolution(); // call NDA function <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		t_lo=clock()-t_lo; // stop low-order timer
		dt_lo.push_back(((double)t_lo)/CLOCKS_PER_SEC); // add high-order solution time to vector
		cout<<"Iteration Completed in "<<dt_ho[s-1]+dt_lo[s-1]<<" sec\n";
		
		// Find new norm
		norm_si=0;
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) {
				if ( abs(phiT[i][j]-phiL[i][j])>norm_si ) norm_si=abs(phiT[i][j]-phiL[i][j]);
			}
		}
		rho_si.push_back(norm_si/norm_siL);
		s++;
	}
	n_iterations=s;
	
	// for periodic case normalize solution
	if ( kbc==3 ) {
		phi_sum=0.0;
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) phi_sum+=phi[i][j]*hx[i]*hy[j];
		}
		phi_sum/=4*pi; // normalize to integrate to 4pi
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) {
				phi[i][j] /=phi_sum;
				phiT[i][j]/=phi_sum;
			}
		}
		for (i=0; i<nx+1; i++) {
			for (j=0; j<ny; j++) {
				j_x[i][j]   /=phi_sum;
				phi_xT[i][j]/=phi_sum;
				j_xT[i][j]  /=phi_sum;
			}
		}
		for (i=0; i<nx; i++) {
			for (j=0; j<ny+1; j++) {
				j_y[i][j]   /=phi_sum;
				phi_yT[i][j]/=phi_sum;
				j_yT[i][j]  /=phi_sum;
			}
		}
	}
}
//======================================================================================//

//======================================================================================//
//++ Initial Transport Solution Found ++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void initialAngleSweep()
{
	int i, j, m;
	double psi_const=1.0;
	
	// Find Initial Flux and Current from constant angular flux
	for ( i=0; i<nx; i++ ) {
		for ( j=0; j<ny; j++ ) {
			phiT[i][j]=0.0; phi[i][j]=0.0;
			for ( m=0; m<sn; m++ ) phiT[i][j]+=psi_const*w[m];
			phiT[i][j]*=4.0;
		}
	}
	for ( i=0; i<nx+1; i++ ) { // X Grid
		for ( j=0; j<ny; j++ ) {
			phi_xT[i][j]=0.0; j_xT[i][j]=0.0;
			for ( m=0; m<sn; m++ ) {
				j_xT[i][j]+=(psi_const-psi_const)*w[m]*mu[m];
				phi_xT[i][j]+=psi_const*w[m];
			}
			phi_xT[i][j]*=4.0; j_xT[i][j]*=2.0;
		}
	}
	for ( i=0; i<nx; i++ ) { // Y Grid
		for ( j=0; j<ny+1; j++ ) {
			phi_yT[i][j]=0.0; j_yT[i][j]=0.0;
			for ( m=0; m<sn; m++ ) {
				j_yT[i][j]+=(psi_const-psi_const)*w[m]*mu[m];
				phi_yT[i][j]+=psi_const*w[m];
			}
			phi_yT[i][j]*=4.0; j_yT[i][j]*=2.0;
		}
	}
}

//======================================================================================//
//++ Determine how to sweep thought cells ++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void angleSweep()
{
	int i, j, m;
	cout<<"Transport Sweep Started : ";
	// Zero out Eddington Factors and Currents
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) phiT[i][j]=0.0;
	}
	for (i=0; i<nx+1; i++) {
		for (j=0; j<ny; j++) {
			phi_xT[i][j]=0.0; j_xT[i][j]=0.0;
		}
	}
	for (i=0; i<nx; i++) {
		for (j=0; j<ny+1; j++) {
			phi_yT[i][j]=0.0; j_yT[i][j]=0.0;
		}
	}
	
	switch ( kbc ) {
	case 1: // incoming boundary conditions on face only 111111111111111111111111111111111111111111111111111111111111111111
		for (m=0; m<sn; m++) {
			for (i=0; i<nx; i++) { // top and bottom boundary conditions
				psiB[i][m][1]=bcB[i]; // Quad 1 Bottom boundary condition
				psiB[i][m][2]=bcB[i]; // Quad 2 Bottom boundary condition
				psiT[i][m][3]=bcT[i]; // Quad 3 Top boundary condition
				psiT[i][m][4]=bcT[i]; // Quad 4 Top boundary condition
			}
			for (j=0; j<ny; j++) { // left and right boundary conditions
				psiL[j][m][1]=bcL[j]; // Quad 1 Left boundary condition
				psiL[j][m][4]=bcL[j]; // Quad 4 Left boundary condition
				psiR[j][m][2]=bcR[j]; // Quad 2 Right boundary condition
				psiR[j][m][3]=bcR[j]; // Quad 3 Right boundary condition
			}
		}
		
		// Start solution sweep///////////////////////////////////////
		quad1(); // solution sweep through quadrant 1
		quad2(); // solution sweep through quadrant 2
		quad3(); // solution sweep through quadrant 3
		quad4(); // solution sweep through quadrant 4
		break;
	case 2: // incoming boundary conditions in Quadrants 22222222222222222222222222222222222222222222222222222222222222222222
		for (m=0; m<sn; m++) {
			for (i=0; i<nx; i++) { // top and bottom boundary conditions
				psiB[i][m][1]=bcL[i]; // Quad 1 Bottom boundary condition
				psiB[i][m][2]=bcB[i]; // Quad 2 Bottom boundary condition
				psiT[i][m][3]=bcR[i]; // Quad 3 Top    boundary condition
				psiT[i][m][4]=bcT[i]; // Quad 4 Top    boundary condition
			}
			for (j=0; j<ny; j++) { // left and right boundary conditions
				psiL[j][m][1]=bcL[j]; // Quad 1 Left  boundary condition
				psiL[j][m][4]=bcT[j]; // Quad 4 Left  boundary condition
				psiR[j][m][2]=bcB[j]; // Quad 2 Right boundary condition
				psiR[j][m][3]=bcR[j]; // Quad 3 Right boundary condition
			}
		}
		
		// Start solution sweep///////////////////////////////////////
		quad1(); // solution sweep through quadrant 1
		quad2(); // solution sweep through quadrant 2
		quad3(); // solution sweep through quadrant 3
		quad4(); // solution sweep through quadrant 4
		break;
	case 3: // All Reflective BC 333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333

		// Start solution sweep
		for (m=0; m<sn; m++) { // quad 1
			for (i=0; i<nx; i++) {
				psiB[i][m][1]=psiB[i][m][4]; // Quad 1 Bottom Reflective BC
				psiB[i][m][2]=psiB[i][m][3]; // Quad 2 Bottom Reflective BC
				psiT[i][m][3]=psiT[i][m][2]; // Quad 3 Top Reflective BC
				psiT[i][m][4]=psiT[i][m][1]; // Quad 4 Top Reflective BC
			}
			for (j=0; j<ny; j++) {
				psiL[j][m][1]=psiL[j][m][2]; // Quad 1 Left Reflective BC
				psiR[j][m][2]=psiR[j][m][1]; // Quad 2 Right Reflective BC
				psiR[j][m][3]=psiR[j][m][4]; // Quad 3 Right Reflective BC
				psiL[j][m][4]=psiL[j][m][3]; // Quad 4 Left Reflective BC
			}
		}
		
		quad1(); // solution sweep through quadrant 1
		quad2(); // solution sweep through quadrant 2
		quad3(); // solution sweep through quadrant 3
		quad4(); // solution sweep through quadrant 4
			
		break;
	case 4: //Reflective BC on bottom and left, face BC 444444444444444444444444444444444444444444444444444444444444444444444444
		// Start solution sweep ///////////////////////////////////////////////////////////
		for (m=0; m<sn; m++) { // quad 3
			for (i=0; i<nx; i++) psiT[i][m][3]=bcT[i]; // Quad 3 Top BC
			for (j=0; j<ny; j++) psiR[j][m][3]=bcR[j]; // Quad 3 Right BC
		}
		quad3(); // solution sweep through quadrant 3
		
		for (m=0; m<sn; m++) { // quad 4
			for (i=0; i<nx; i++) psiT[i][m][4]=bcT[i];           // Quad 4 Top BC
			for (j=0; j<ny; j++) psiL[j][m][4]=psiL[j][m][3]; // Quad 4 Left Reflective BC
		}
		quad4(); // solution sweep through quadrant 4
		
		for (m=0; m<sn; m++) { // quad 2
			for (i=0; i<nx; i++) psiB[i][m][2]=psiB[i][m][3]; // Quad 2 Bottom Reflective BC
			for (j=0; j<ny; j++) psiR[j][m][2]=bcR[j];           // Quad 2 Right BC
		}
		quad2(); // solution sweep through quadrant 2
		
		for (m=0; m<sn; m++) { // quad 1
			for (i=0; i<nx; i++) psiB[i][m][1]=psiB[i][m][4]; // Quad 1 Bottom Reflective BC
			for (j=0; j<ny; j++) psiL[j][m][1]=psiL[j][m][2]; // Quad 1 Left Reflective BC
		}
		quad1(); // solution sweep through quadrant 1
		break;
	case 5: // Reflective BC on top and right, face BC 55555555555555555555555555555555555555555555555555555555555555555555555555
		// Start solution sweep ////////////////////////////////////////////////////////
		for (m=0; m<sn; m++) {
			for (i=0; i<nx; i++) psiB[i][m][1]=bcB[i]; // Quad 1 Bottom boundary condition
			for (j=0; j<ny; j++) psiL[j][m][1]=bcL[j]; // Quad 1 Left boundary condition
		}
		quad1(); // solution sweep through quadrant 1
		
		for (m=0; m<sn; m++) {
			for (i=0; i<nx; i++) psiB[i][m][2]=bcB[i];           // Quad 2 Bottom BC
			for (j=0; j<ny; j++) psiR[j][m][2]=psiR[j][m][1]; // Quad 2 Right reflective BC
		}
		quad2(); // solution sweep through quadrant 2
		
		for (m=0; m<sn; m++) { // quad 4
			for (i=0; i<nx; i++) psiT[i][m][4]=psiT[i][m][1]; // Quad 4 Top reflective BC
			for (j=0; j<ny; j++) psiL[j][m][4]=bcL[j];           // Quad 4 Left BC
		}
		quad4(); // solution sweep through quadrant 4
		
		for (m=0; m<sn; m++) { // quad 3
			for (i=0; i<nx; i++) psiT[i][m][3]=psiT[i][m][2]; // Quad 3 Top reflective BC
			for (j=0; j<ny; j++) psiR[j][m][3]=psiR[j][m][4]; // Quad 3 Right reflective BC
		}
		quad3(); // solution sweep through quadrant 3
		break;
	default:
		cout<<"bad boundary conditions: incorrect boundary type"<<endl;
		break;
	
	}
	cout<<"Completed \n";
}
//======================================================================================//

//======================================================================================//
//++ sweep through angles and cells in each angular quadrant +++++++++++++++++++++++++++//
//======================================================================================//
void quad1() // solution in quadrant 1
{
	int i, j, m, outw=16;
	double psiA, SA;
	double omega_x, omega_y;
	
	for (m=0; m<sn; m++) { // first quadrant
		omega_x=mu[m]; omega_y=eta[m];
		
		for (i=0; i<nx; i++) psi_y[i][0]=psiB[i][m][1]; // Bottom In BC
		for (j=0; j<ny; j++) psi_x[0][j]=psiL[j][m][1]; // Left In BC
		
		for (j=0; j<ny; j++) { // bottom to top
			for (i=0; i<nx; i++) { // left to right
				//psiInL=psi_x[i][j]; // incoming angular flux on the left
				//psiInB=psi_y[i][j]; // incoming angular flux on the bottom
				SA=((sigmaS[i][j] + nuF[i][j]*sigmaF[i][j])*phi[i][j]+s_ext[i][j])/(4*pi); // source in the cell
				cellSolution( psi_y[i][j], psi_x[i][j], SA, sigmaT[i][j], mu[m], eta[m], xi[m], hx[i], hy[j], psi_y[i][j+1], psi_x[i+1][j], psiA );
				//psi_x[i+1][j]=psiOutR; // outgoing angular flux on the right
				//psi_y[i][j+1]=psiOutT; // outgoing angular flux on the top
				psi[i][j]=psiA;        // cell average angular flux
				// Calculate cell centre values
				phiT[i][j]+=psiA*w[m];
			}
		}
		// Calculate Cell Edge values
		for (i=0; i<nx+1; i++) {
			for (j=0; j<ny; j++) {
				phi_xT[i][j]+=psi_x[i][j]*w[m];
				j_xT[i][j]  +=omega_x*psi_x[i][j]*w[m];
			}
		}
		for (i=0; i<nx; i++) {
			for (j=0; j<ny+1; j++) {
				phi_yT[i][j]+=psi_y[i][j]*w[m];
				j_yT[i][j]  +=omega_y*psi_y[i][j]*w[m];
			}
		}
		
		for (i=0; i<nx; i++) psiT[i][m][1]=psi_y[i][ny]; // Top Out BC
		for (j=0; j<ny; j++) psiR[j][m][1]=psi_x[nx][j]; // Right Out BC
		
		// option to print out angular flux
		if ( o_angular ) {
			temfile<<" Quadrant 1 Direction # ,"<<m<<",\n";
			temfile<<"  Qmega_x ,"<<print_csv(omega_x)<<"  Omega_y ,"<<print_csv(omega_y)<<"  Weight ,"<<print_csv(w[m])<<endl;
			temfile<<" -- Cell Averaged Angular Flux -- \n";
			write_cell_average_dat(psi, outw, temfile); // call function to write out cell average scalar flux
			temfile<<" -- X Vertical Cell Edge Angular Flux -- \n";
			write_cell_edge_x_dat(psi_x, outw, temfile); // call function to write out cell edge scalar flux on x grid
			temfile<<" -- Y Horizontal Cell Edge Angular Flux -- \n";
			write_cell_edge_y_dat(psi_y, outw, temfile); // call function to write out cell edge scalar flux on y grid
		}
	}
}
//======================================================================================//
void quad2() // solution in quadrant 2
{
	int i, j, m, outw=16;
	double psiA, SA;
	double omega_x, omega_y;
	
	for (m=0; m<sn; m++) { // second quadrant
		omega_x=-mu[m]; omega_y=eta[m];
		
		for (i=0; i<nx; i++) psi_y[i][0]=psiB[i][m][2]; // Bottom In BC
		for (j=0; j<ny; j++) psi_x[nx][j]=psiR[j][m][2]; // Right In BC
		
		for (j=0; j<ny; j++) { // bottom to top
			for (i=nx-1; i>=0; i--) { // right to left
				//psiInL=psi_x[i+1][j]; // incoming angular flux on the left
				//psiInB=psi_y[i][j];   // incoming angular flux on the bottom
				SA=((sigmaS[i][j] + nuF[i][j]*sigmaF[i][j])*phi[i][j]+s_ext[i][j])/(4*pi); // source in the cell
				cellSolution( psi_y[i][j], psi_x[i+1][j], SA, sigmaT[i][j], mu[m], eta[m], xi[m], hx[i], hy[j], psi_y[i][j+1], psi_x[i][j], psiA );
				//psi_x[i][j]=psiOutR;   // outgoing angular flux on the right
				//psi_y[i][j+1]=psiOutT; // outgoing angular flux on the top
				psi[i][j]=psiA;        // cell average angular flux
				// Calculate cell centre values
				phiT[i][j]+=psiA*w[m];
			}
		}
		// Calculate Cell Edge values
		for (i=0; i<nx+1; i++) {
			for (j=0; j<ny; j++) {
				phi_xT[i][j]+=psi_x[i][j]*w[m];
				j_xT[i][j]  +=omega_x*psi_x[i][j]*w[m];
			}
		}
		for (i=0; i<nx; i++) {
			for (j=0; j<ny+1; j++) {
				phi_yT[i][j]+=psi_y[i][j]*w[m];
				j_yT[i][j]  +=omega_y*psi_y[i][j]*w[m];
			}
		}
		
		// Out BC
		for (i=0; i<nx; i++) psiT[i][m][2]=psi_y[i][ny]; // Top Out BC
		for (j=0; j<ny; j++) psiL[j][m][2]=psi_x[0][j]; // Left Out BC
		
		// option to print out angular flux
		if ( o_angular ) {
			temfile<<" Quadrant 2 Direction # ,"<<m<<",\n";
			temfile<<"  Qmega_x ,"<<print_csv(omega_x)<<"  Omega_y ,"<<print_csv(omega_y)<<"  Weight ,"<<print_csv(w[m])<<endl;
			temfile<<" -- Cell Averaged Angular Flux -- \n";
			write_cell_average_dat(psi, outw, temfile); // call function to write out cell average scalar flux
			temfile<<" -- X Vertical Cell Edge Angular Flux -- \n";
			write_cell_edge_x_dat(psi_x, outw, temfile); // call function to write out cell edge scalar flux on x grid
			temfile<<" -- Y Horizontal Cell Edge Angular Flux -- \n";
			write_cell_edge_y_dat(psi_y, outw, temfile); // call function to write out cell edge scalar flux on y grid
		}
		
	}
}
//======================================================================================//
void quad3() // solution in quadrant 3
{
	int i, j, m, outw=16;
	double psiA, SA;
	double omega_x, omega_y;
	
	for (m=0; m<sn; m++) { // third quadrant
		omega_x=-mu[m]; omega_y=-eta[m];
		
		for (i=0; i<nx; i++) psi_y[i][ny]=psiT[i][m][3]; // Top In BC
		for (j=0; j<ny; j++) psi_x[nx][j]=psiR[j][m][3]; // Right In BC
		
		for (j=ny-1; j>=0; j--) { // top to bottom
			for (i=nx-1; i>=0; i--) { // right to left
				//psiInL=psi_x[i+1][j]; // incoming angular flux on the left
				//psiInB=psi_y[i][j+1]; // incoming angular flux on the bottom
				SA=((sigmaS[i][j] + nuF[i][j]*sigmaF[i][j])*phi[i][j]+s_ext[i][j])/(4*pi); // source in the cell
				cellSolution( psi_y[i][j+1], psi_x[i+1][j], SA, sigmaT[i][j], mu[m], eta[m], xi[m], hx[i], hy[j], psi_y[i][j], psi_x[i][j], psiA );
				//psi_x[i][j]=psiOutR; // outgoing angular flux on the right
				//psi_y[i][j]=psiOutT; // outgoing angular flux on the top
				psi[i][j]=psiA;      // cell average angular flux
				// Calculate cell centre values
				phiT[i][j]+=psiA*w[m];
			}
		}
		// Calculate Cell Edge values
		for (i=0; i<nx+1; i++) {
			for (j=0; j<ny; j++) {
				phi_xT[i][j]+=psi_x[i][j]*w[m];
				j_xT[i][j]  +=omega_x*psi_x[i][j]*w[m];
			}
		}
		for (i=0; i<nx; i++) {
			for (j=0; j<ny+1; j++) {
				phi_yT[i][j]+=psi_y[i][j]*w[m];
				j_yT[i][j]  +=omega_y*psi_y[i][j]*w[m];
			}
		}
		
		for (i=0; i<nx; i++) psiB[i][m][3]=psi_y[i][0]; // Bottom Out BC
		for (j=0; j<ny; j++) psiL[j][m][3]=psi_x[0][j]; // Left Out BC
		
		// option to print out angular flux
		if ( o_angular ) {
			temfile<<" Quadrant 3 Direction # ,"<<m<<",\n";
			temfile<<"  Qmega_x ,"<<print_csv(omega_x)<<"  Omega_y ,"<<print_csv(omega_y)<<"  Weight ,"<<print_csv(w[m])<<endl;
			temfile<<" -- Cell Averaged Angular Flux -- \n";
			write_cell_average_dat(psi, outw, temfile); // call function to write out cell average scalar flux
			temfile<<" -- X Vertical Cell Edge Angular Flux -- \n";
			write_cell_edge_x_dat(psi_x, outw, temfile); // call function to write out cell edge scalar flux on x grid
			temfile<<" -- Y Horizontal Cell Edge Angular Flux -- \n";
			write_cell_edge_y_dat(psi_y, outw, temfile); // call function to write out cell edge scalar flux on y grid
		}
	}
}
//======================================================================================//
void quad4() // solution in quadrant 4
{
	int i, j, m, outw=16;
	double psiA, SA;
	double omega_x, omega_y;
	
	for (m=0; m<sn; m++) { // fourth quadrant
		omega_x=mu[m]; omega_y=-eta[m];
		
		for (i=0; i<nx; i++) psi_y[i][ny]=psiT[i][m][4]; // Top In BC
		for (j=0; j<ny; j++) psi_x[0][j]=psiL[j][m][4]; // Left In BC
		
		for (j=ny-1; j>=0; j--) { // top to bottom
			for (i=0; i<nx; i++) { // left to right
				//psiInL=psi_x[i][j];   // incoming angular flux on the left
				//psiInB=psi_y[i][j+1]; // incoming angular flux on the bottom
				SA=((sigmaS[i][j] + nuF[i][j]*sigmaF[i][j])*phi[i][j]+s_ext[i][j])/(4*pi); // source in the cell
				cellSolution( psi_y[i][j+1], psi_x[i][j], SA, sigmaT[i][j], mu[m], eta[m], xi[m], hx[i], hy[j], psi_y[i][j], psi_x[i+1][j], psiA );
				//psi_x[i+1][j]=psiOutR; // outgoing angular flux on the right
				//psi_y[i][j]=psiOutT;   // outgoing angular flux on the top
				psi[i][j]=psiA;        // cell average angular flux
				// Calculate cell centre values
				phiT[i][j]+=psiA*w[m];
			}
		}
		// Calculate Cell Edge values
		for (i=0; i<nx+1; i++) {
			for (j=0; j<ny; j++) {
				phi_xT[i][j]+=psi_x[i][j]*w[m];
				j_xT[i][j]  +=omega_x*psi_x[i][j]*w[m];
			}
		}
		for (i=0; i<nx; i++) {
			for (j=0; j<ny+1; j++) {
				phi_yT[i][j]+=psi_y[i][j]*w[m];
				j_yT[i][j]  +=omega_y*psi_y[i][j]*w[m];
			}
		}
		
		for (i=0; i<nx; i++) psiB[i][m][4]=psi_y[i][0]; // Bottom Out BC
		for (j=0; j<ny; j++) psiR[j][m][4]=psi_x[nx][j]; // Right Out BC
		
		// option to print out angular flux
		if ( o_angular ) {
			temfile<<" Quadrant 4 Direction # ,"<<m<<",\n";
			temfile<<"  Qmega_x ,"<<print_csv(omega_x)<<"  Omega_y ,"<<print_csv(omega_y)<<"  Weight ,"<<print_csv(w[m])<<endl;
			temfile<<" -- Cell Averaged Angular Flux -- \n";
			write_cell_average_dat(psi, outw, temfile); // call function to write out cell average scalar flux
			temfile<<" -- X Vertical Cell Edge Angular Flux -- \n";
			write_cell_edge_x_dat(psi_x, outw, temfile); // call function to write out cell edge scalar flux on x grid
			temfile<<" -- Y Horizontal Cell Edge Angular Flux -- \n";
			write_cell_edge_y_dat(psi_y, outw, temfile); // call function to write out cell edge scalar flux on y grid
		}
	}
}
//======================================================================================//

//======================================================================================//
//++ function to solve transport in a single general cell ++++++++++++++++++++++++++++++//
//======================================================================================//
void cellSolution(double psiInB, double psiInL, double SA, double sigma, double mut, double etat, double xit, 
double LT, double LR, double& psiOutT, double& psiOutR, double& psiA  )
{
	double epsilon, exp_epsilon, epsilon_2, mup, muc, mu, du;
	double psiOut1, psiOut2, psiOut3, psiA1, psiA2, psiA3;
	double A1, A2, Lout1, Lout2, Lout3;
	
	mup=sqrt(1-xit*xit);             // mu'
	muc=LT/sqrt(LT*LT+LR*LR);        // mu of cell
	mu=mut/sqrt(mut*mut+etat*etat); // projection onto x-y plane
	
	if ( abs(mu-muc)<1e-15 ) { // ray passes through both corners of cell
		du=sqrt(LT*LT+LR*LR);
		if ( sigma<1e-10 ) {
			// triangle A
			psiOutT=psiInL+SA*du/mup/2.0; // find out going angular flux
			psiA1  =psiInL+SA*du/mup/3; // find cell angular flux
			
			// triangle C
			psiOutR=psiInB+SA*du/mup/2.0; // find out going angular flux
			psiA3  =psiInB+SA*du/mup/3; // find cell angular flux
			
			psiA=0.5*(psiInL+psiInB)+SA*du/mup/3.0;
		}
		else {
			epsilon=sigma*du/mup; // optical thickness
			exp_epsilon=exp(-epsilon); // exponent of epsilon
			epsilon_2=epsilon*epsilon; // square of epsilon
			// triangle A
			psiOutT=(psiInL*(1-exp_epsilon)+SA*(epsilon+exp_epsilon-1)/sigma)/epsilon; // find out going angular flux
			psiA1=2*(psiInL*(epsilon+exp_epsilon-1)+SA*(1+0.5*epsilon_2-exp_epsilon-epsilon)/sigma)/epsilon_2; // find cell angular flux
			
			// triangle C
			psiOutR=(psiInB*(1-exp_epsilon)+SA*(epsilon+exp_epsilon-1)/sigma)/epsilon; // find out going angular flux
			psiA3=2*(psiInB*(epsilon+exp_epsilon-1)+SA*(1+0.5*epsilon_2-exp_epsilon-epsilon)/sigma)/epsilon_2; // find cell angular flux
			
			psiA=((psiInL+psiInB)*(epsilon+exp_epsilon-1.0)+2.0*SA*(1.0+0.5*epsilon_2-exp_epsilon-epsilon)/sigma)/epsilon_2;
		}
		
		// total
		//psiOutT=psiOut1;
		//psiOutR=psiOut3;
		//psiA=0.5*(psiA1+psiA3);
		
		//cout<<"balance 1"<<endl;
		//cout<<"T1 balance "<<SA*du/mup-2*(psiOutT-psiInL)-epsilon*psiA1<<endl; // triangle 1 balance
		//cout<<"T3 balance "<<SA*du/mup-2*(psiOutR-psiInB)-epsilon*psiA3<<endl; // triangle 3 balance
		//cout<<"cell angular balance "<<SA*LT*LR-mut*LR*(psiOutR-psiInL)-etat*LT*(psiOutT-psiInB)-LT*LR*sigma*psiA<<endl;
		//cout<<"cell ang bal "<<SA*LT*LR-LR*(psiOutR-psiInL)-LT*(psiOutT-psiInB)-LT*LR*sigma*psiA<<endl;
		//cout<<SA*du-mup*(psiOutR-psiInL)-mup*(psiOutT-psiInB)-sigma*du*psiA<<endl;
	}
	else if ( mu<muc ) { // ray splits the top and bottom of the cell
		Lout1=mu*LR/sqrt(1.0-mu*mu);
		du=Lout1/mu;
		A1=Lout1*LR/2.0; // Triangle 1 Area
		Lout2=LT-Lout1;
		A2=LT*LR-2*A1; // Parallelogram 2 Area
		if ( sigma<1e-10 ) {
			// triangle A
			psiOut1=psiInL+SA*du/mup/2; // find out going angular flux
			psiA1  =psiInL+SA*du/mup/3; // find cell angular flux
			
			// parallelogram B
			psiOut2=psiInB+SA*du/mup;   // find out going angular flux
			psiA2  =psiInB+SA*du/mup/2; // find cell angular flux
			
			// triangle C
			psiOutR=psiA2;              // find out going angular flux
			psiA3  =psiInB+SA*du/mup/3; // find cell angular flux
		}
		else {
			epsilon=sigma*du/mup; // optical thickness
			exp_epsilon=exp(-epsilon); // exponent of epsilon
			epsilon_2=epsilon*epsilon; // square of epsilon
			// triangle A
			psiOut1=(psiInL*(1-exp_epsilon)+SA*(epsilon+exp_epsilon-1)/sigma)/epsilon; // find out going angular flux
			psiA1=2*(psiInL*(epsilon+exp_epsilon-1)+SA*(1+0.5*epsilon_2-exp_epsilon-epsilon)/sigma)/epsilon_2; // find cell angular flux
			
			// parallelogram B
			psiOut2=psiInB*exp_epsilon+SA*(1-exp_epsilon)/sigma;                       // find out going angular flux
			psiA2  =(psiInB*(1-exp_epsilon)+SA*(epsilon+exp_epsilon-1)/sigma)/epsilon; // find cell angular flux
			
			// triangle C
			psiOutR=psiA2;                                                                                       // find out going angular flux
			psiA3  =2*(psiInB*(epsilon+exp_epsilon-1)+SA*(1+0.5*epsilon_2-exp_epsilon-epsilon)/sigma)/epsilon_2; // find cell angular flux
			
		}
		
		// total
		psiOutT=(Lout1*psiOut1+Lout2*psiOut2)/LT;
		psiA=(A1*psiA1+A2*psiA2+A1*psiA3)/(LT*LR);
		
		//cout<<"balance 2"<<endl;
		//cout<<"T1 balance "<<SA*du/mup-2*(psiOut1-psiInL)-epsilon*psiA1<<endl; // triangle 1 balance
		//cout<<"P2 balance "<<SA*du/mup-(psiOut2-psiInB)-epsilon*psiA2<<endl; // parallelogram 2 balance
		//cout<<"T3 balance "<<SA*du/mup-2*(psiOutR-psiInB)-epsilon*psiA3<<endl; // triangle 3 balance
		//cout<<"cell angular balance "<<SA*LT*LR-mut*LR*(psiOutR-psiInL)-etat*LT*(psiOutT-psiInB)-LT*LR*sigma*psiA<<endl;
	}
	else { // ray splits the right and left side of the cell
		du=LT/mu;
		A1=sqrt(du*du-LT*LT)*LT/2.0; // Triangle 1 Area
		Lout3=sqrt(du*du-LT*LT); // Triangle 3 Length
		Lout2=LR-Lout3;
		A2=LT*LR-2*A1; // Parallelogram 2 Area
		
		if ( sigma<1e-10 ) {
			// triangle A
			psiOutT=psiInL+SA*du/mup/2; // find out going angular flux
			psiA1  =psiInL+SA*du/mup/3; // find cell angular flux
			
			// parallelogram B
			psiOut2=psiInL+SA*du/mup; // find out going angular flux
			psiA2  =psiOutT;          // find cell angular flux
			
			// triangle C
			psiOut3=psiInB+SA*du/mup/2; // find out going angular flux
			psiA3  =psiInB+SA*du/mup/2; // find cell angular flux
			
		}
		else {
			epsilon=sigma*du/mup; // optical thickness
			exp_epsilon=exp(-epsilon); // exponent of epsilon
			epsilon_2=epsilon*epsilon; // square of epsilon
			// triangle A
			psiOutT=(psiInL*(1-exp_epsilon)+SA*(epsilon+exp_epsilon-1)/sigma)/epsilon; // find out going angular flux
			psiA1=2*(psiInL*(epsilon+exp_epsilon-1)+SA*(1+0.5*epsilon_2-exp_epsilon-epsilon)/sigma)/epsilon_2; // find cell angular flux
			
			// parallelogram B
			psiOut2=psiInL*exp_epsilon+SA*(1-exp_epsilon)/sigma; // find out going angular flux
			psiA2=psiOutT;                                       // find cell angular flux
			
			// triangle C
			psiOut3=(psiInB*(1-exp_epsilon)+SA*(epsilon+exp_epsilon-1)/sigma)/epsilon; // find out going angular flux
			psiA3=2*(psiInB*(epsilon+exp_epsilon-1)+SA*(1+0.5*epsilon_2-exp_epsilon-epsilon)/sigma)/epsilon_2; // find cell angular flux
		}
		
		// total
		psiOutR=(Lout3*psiOut3+Lout2*psiOut2)/LR;
		psiA=(A1*psiA1+A2*psiA2+A1*psiA3)/(LT*LR);
		
		//cout<<"balance 3"<<endl;
		//cout<<"T1 balance "<<SA*du/mup-2*(psiOutT-psiInL)-epsilon*psiA1<<endl; // triangle 1 balance
		//cout<<"P2 balance "<<SA*du/mup-(psiOut2-psiInL)-epsilon*psiA2<<endl; // parallelogram 2 balance
		//cout<<"T3 balance "<<SA*du/mup-2*(psiOut3-psiInB)-epsilon*psiA3<<endl; // triangle 3 balance
		//cout<<"cell angular balance "<<SA*LT*LR-mut*LR*(psiOutR-psiInL)-etat*LT*(psiOutT-psiInB)-LT*LR*sigma*psiA<<endl;
	}
	//cout<<"cell angular balance "<<SA*LT*LR-mut*LR*(psiOutR-psiInL)-etat*LT*(psiOutT-psiInB)-LT*LR*sigma*psiA<<endl;
}
//======================================================================================//

//======================================================================================//
//++ function to solve NDA problem ++++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void NDAsolution()
{
	int i, j, p, c, N_ukn, N_c, N_xf, N_yf;
	double res=0;
	
	cout<<"NDA Solution Started : ";
	temfile<<"NDA Solution\n";
	
	N_c=nx*ny;
	N_xf=2*nx;
	N_yf=2*ny;
	N_ukn=nx*ny+2*nx+2*ny;
	
	VectorXd x(N_ukn);
	VectorXd b(N_ukn);
    SpMat A(N_ukn,N_ukn);
	
	p=0; // initialize solution vector guess to transport solution
	for (j=0; j<ny; j++) {
		for (i=0; i<nx; i++) {
			x[p]=phiL[i][j]; p++; // cell centre flux
		}
	}
	for (j=0; j<ny; j++) {
		x[p]=phiB_L[j]; p++; // Left   BC Flux
	}
	for (i=0; i<nx; i++) {
		x[p]=phiB_B[i]; p++; // Bottom BC Flux
	}
	for (j=0; j<ny; j++) {
		x[p]=phiB_R[j]; p++; // Right BC Flux
	}
	for (i=0; i<nx; i++) {
		x[p]=phiB_T[i]; p++; // Top   BC Flux
	}
	
	// Assign matrix A and vector b
	p=0; // ++++++++ Central Cell ++++++++
	for ( j=0; j<ny; j++ ) {
		for ( i=0; i<nx; i++ ) {
			A.insert(p,p)=hy[j]*D_xN[i+1][j]/xe[i+1] + hy[j]*D_xP[i][j]/xe[i] + hx[i]*D_yN[i][j+1]/ye[j+1] + hx[i]*D_yP[i][j]/ye[j] +
			hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]);
			b[p]=hx[i]*hy[j]*s_ext[i][j]; p++;
		}
	}
	// ++++++++ Periphery Cells ++++++++
	p=0;
	if ( nx==1 and ny==1 ) { // Single Cell
		// Balance Equations
		A.insert(0,1)=-hy[0]*D_xN[0][0]/xe[0];   A.insert(0,2)=-hx[0]*D_yN[0][0]/ye[0];
		A.insert(0,3)=-hy[0]*D_xP[1][0]/xe[1];   A.insert(0,4)=-hx[0]*D_yP[0][1]/ye[1];
	}
	else if ( nx==1 ) { // One Cell Wide
		// Bottom Cell
		A.insert(p,N_c)      =-hy[0]*D_xN[0][0]/xe[0];   A.insert(p,N_c+ny)=-hx[0]*D_yN[0][0]/ye[0];
		A.insert(p,N_c+ny+nx)=-hy[0]*D_xP[1][0]/xe[1];   A.insert(p,1)     =-hx[0]*D_yP[0][1]/ye[1];
		p++;
		for ( j=1; j<ny-1; j++ ) { // Middle Cells
			A.insert(p,N_c+j)      =-hy[j]*D_xN[0][j]/xe[0];   A.insert(p,j-1)=-hx[0]*D_yN[0][j]  /ye[j];
			A.insert(p,N_c+ny+nx+j)=-hy[j]*D_xP[1][j]/xe[1];   A.insert(p,j+1)=-hx[0]*D_yP[0][j+1]/ye[j+1];
			p++;
		}
		// Top Cell
		A.insert(p,N_c+ny-1)     =-hy[ny-1]*D_xN[0][ny-1]/xe[0];   A.insert(p,ny-2)           =-hx[0]*D_yN[0][ny-1]/ye[ny-1];
		A.insert(p,N_c+2*ny+nx-1)=-hy[ny-1]*D_xP[1][ny-1]/xe[1];   A.insert(p,N_c+2*ny+2*nx-1)=-hx[0]*D_yP[0][ny]  /ye[ny];
		p++;		
	}
	else if ( ny==1 ) { // One Cell Tall
		// Left Cell
		A.insert(p,N_c)=-hy[0]*D_xN[0][0]/xe[0];   A.insert(p,N_c+ny)     =-hx[0]*D_yN[0][0]/ye[0];
		A.insert(p,1)  =-hy[0]*D_xP[1][0]/xe[1];   A.insert(p,N_c+2*ny+nx)=-hx[0]*D_yP[0][1]/ye[1];
		p++;
		for ( i=1; i<nx-1; i++ ) { // Middle Cell
			A.insert(p,i-1)=-hy[0]*D_xN[i][0]  /xe[i];    A.insert(p,N_c+ny+i)     =-hx[i]*D_yN[i][0]/ye[0];
			A.insert(p,i+1)=-hy[0]*D_xP[i+1][0]/xe[i+1];  A.insert(p,N_c+2*ny+nx+i)=-hx[i]*D_yP[i][1]/ye[1];
			p++;
		}
		// Right Cell
		A.insert(p,nx-2)       =-hy[0]*D_xN[nx-1][0]/xe[nx-1];  A.insert(p,N_c+ny+nx-1)    =-hx[nx-1]*D_yN[nx-1][0]/ye[0];
		A.insert(p,N_c+2*ny+nx-1)=-hy[0]*D_xP[nx][0]  /xe[nx];    A.insert(p,N_c+2*ny+2*nx-1)=-hx[nx-1]*D_yP[nx-1][1]/ye[1];
		p++;
		
	}
	else {
		// Bottom Left Corner
		A.insert(p,N_c)=-hy[0]*D_xN[0][0]/xe[0];   A.insert(p,N_c+ny)=-hx[0]*D_yN[0][0]/ye[0];
		A.insert(p,1)  =-hy[0]*D_xP[1][0]/xe[1];   A.insert(p,nx)    =-hx[0]*D_yP[0][1]/ye[1];
		p++;
		for ( i=1; i<nx-1; i++ ) { // Bottom
			A.insert(p,p-1)=-hy[0]*D_xN[i][0]  /xe[i];    A.insert(p,N_c+ny+i)=-hx[i]*D_yN[i][0]/ye[0];
			A.insert(p,p+1)=-hy[0]*D_xP[i+1][0]/xe[i+1];  A.insert(p,nx+i)    =-hx[i]*D_yP[i][1]/ye[1];
			p++;
		}
		// Bottom Right Corner
		A.insert(p,nx-2)     =-hy[0]*D_xN[nx-1][0]/xe[nx-1]; A.insert(p,N_c+ny+nx-1)=-hx[nx-1]*D_yN[nx-1][0]/ye[0];
		A.insert(p,N_c+ny+nx)=-hy[0]*D_xP[nx][0]  /xe[nx];   A.insert(p,2*nx-1)     =-hx[nx-1]*D_yP[nx-1][1]/ye[1];
		p++;
		for ( j=1; j<ny-1; j++ ) { // Middle
			i=0; // Left Side
			A.insert(p,N_c+j)=-hy[j]*D_xN[i][j]  /xe[i];    A.insert(p,p-nx)=-hx[i]*D_yN[i][j]  /ye[j];
			A.insert(p,p+1)  =-hy[j]*D_xP[i+1][j]/xe[i+1];  A.insert(p,p+nx)=-hx[i]*D_yP[i][j+1]/ye[j+1];
			p++;
			for ( i=1; i<nx-1; i++ ) { // Centre Cells
				A.insert(p,p-1)=-hy[j]*D_xN[i][j]  /xe[i];    A.insert(p,p-nx)=-hx[i]*D_yN[i][j]  /ye[j];
				A.insert(p,p+1)=-hy[j]*D_xP[i+1][j]/xe[i+1];  A.insert(p,p+nx)=-hx[i]*D_yP[i][j+1]/ye[j+1];
				p++;
			}
			i=nx-1; // Right Side
			A.insert(p,p-1)        =-hy[j]*D_xN[i][j]  /xe[i];    A.insert(p,p-nx)=-hx[i]*D_yN[i][j]  /ye[j];
			A.insert(p,N_c+ny+nx+j)=-hy[j]*D_xP[i+1][j]/xe[i+1];  A.insert(p,p+nx)=-hx[i]*D_yP[i][j+1]/ye[j+1];
			p++;
		}
		j=ny-1; // Top Left Corner
		A.insert(p,N_c+j)=-hy[j]*D_xN[0][j]/xe[0];   A.insert(p,p-nx)       =-hx[0]*D_yN[0][j]  /ye[j];
		A.insert(p,p+1)  =-hy[j]*D_xP[1][j]/xe[1];   A.insert(p,N_c+2*ny+nx)=-hx[0]*D_yP[0][j+1]/ye[j+1];
		p++;
		for ( i=1; i<nx-1; i++ ) { // Top
			A.insert(p,p-1)=-hy[j]*D_xN[i][j]  /xe[i];    A.insert(p,p-nx)         =-hx[i]*D_yN[i][j]/ye[j];
			A.insert(p,p+1)=-hy[j]*D_xP[i+1][j]/xe[i+1];  A.insert(p,N_c+2*ny+nx+i)=-hx[i]*D_yP[i][j+1]/ye[j+1];
			p++;
		}
		i=nx-1; // Top Right Corner
		A.insert(p,p-1)          =-hy[j]*D_xN[i][j]  /xe[i];    A.insert(p,p-nx)         =-hx[i]*D_yN[i][j]/ye[j];
		A.insert(p,N_c+2*ny+nx-1)=-hy[j]*D_xP[i+1][j]/xe[i+1];  A.insert(p,N_c+2*ny+nx+i)=-hx[i]*D_yP[i][j+1]/ye[j+1];
		p++;
	}
	// ++++++++ Boundary Conditions ++++++++
	p=N_c;
	// Left BC
	for ( j=0; j<ny; j++ ) {
		A.insert(p,j*nx)= D_xP[0][j]/xe[0];   A.insert(p,N_c+j)=FL[j]-D_xN[0][j]/xe[0];   b[p]=0; // Left   BC
		p++;
	}
	// Bottom BC
	for ( i=0; i<nx; i++ ) {
		A.insert(p,i)= D_yP[i][0]/ye[0];   A.insert(p,N_c+ny+i)=FB[i]-D_yN[i][0]/ye[0];   b[p]=0; // Bottom BC
		p++;
	}
	// Right BC
	for ( j=0; j<ny; j++ ) {
		A.insert(p,(j+1)*nx-1)=-D_xN[nx][j]/xe[nx];   A.insert(p,N_c+ny+nx+j)=FR[j]+D_xP[nx][j]/xe[nx];   b[p]=0; // Right  BC
		p++;
	}
	// Top BC
	for ( i=0; i<nx; i++ ) {
		A.insert(p,(ny-1)*nx+i)=-D_yN[i][ny]/ye[ny];   A.insert(p,N_c+2*ny+nx+i)=FT[i]+D_yP[i][ny]/ye[ny];   b[p]=0; // Top    BC
		p++;
	}
	
	A.prune(1e-17, 10);
	//cout<<"b "<<b<<endl;
	//cout<<A<<endl;
	// Solve Ax=b iteratively with BiCGSTAB
	t_pc=clock(); // start low-order timer
	
	BiCGSTAB<SpMat,IncompleteLUT<double> > solver;
	solver.preconditioner().setFillfactor(11);
	solver.setTolerance(epsilon_lo);     // set convergence criteria 
	solver.setMaxIterations(maxiter_lo); // set the max number of lo iterations
	solver.compute(A);
	
	t_pc=clock()-t_pc; // stop low-order timer
	dt_pc.push_back(((double)t_pc)/CLOCKS_PER_SEC); // add high-order solution time to vector
	
	x = solver.solveWithGuess(b,x);
	err_lo.push_back(solver.error());      // error in lo solution
	num_lo.push_back(solver.iterations()); // number of lo iterations
	
	// set solution back to problem values
	p=0;
	for (j=0; j<ny; j++) {
		for (i=0; i<nx; i++) {
			phi[i][j]=x[p]; // Cell Centre Flux
			p++;
		}
	}
	for (j=0; j<ny; j++) {
		phiB_L[j]=x[p]; // Left   BC Flux
		p++;
	}
	for (i=0; i<nx; i++) {
		phiB_B[i]=x[p]; // Bottom BC Flux
		p++;
	}
	for (j=0; j<ny; j++) {
		phiB_R[j]=x[p]; // Right BC Flux
		p++;
	}
	for (i=0; i<nx; i++) {
		phiB_T[i]=x[p]; // Top   BC Flux
		p++;
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Calculate current ///////////////////////////////////////////////////////////////////////////////////////
	// Boundary currents ///////////////////////////////////////////////////////////////////////////////////////
	for (j=0; j<ny; j++) j_x[0][j]=FL[j]*phiB_L[j]; // Left boundary current J_x
	for (i=0; i<nx; i++) j_y[i][0]=FB[i]*phiB_B[i]; // Bottom boundary current J_y
	
	// Inner currents J_x J_y
	for (j=0; j<ny; j++) {
		for (i=1; i<nx; i++) j_x[i][j]=-(D_xP[i][j]*phi[i][j]-D_xN[i][j]*phi[i-1][j])/xe[i];
	}
	for (j=1; j<ny; j++) {
		for (i=0; i<nx; i++) j_y[i][j]=-(D_yP[i][j]*phi[i][j]-D_yN[i][j]*phi[i][j-1])/ye[j];
	}
	
	// Boundary currents
	for (j=0; j<ny; j++) j_x[nx][j]=FR[j]*phiB_R[j]; // Right boundary current
	for (i=0; i<nx; i++) j_y[i][ny]=FT[i]*phiB_T[i]; // Top boundary current
	cout<<"Completed \n";
}
//======================================================================================//

//======================================================================================//
//++ function to output data +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void output(char fn[])
{
	char outf[25]={""};
	int i, j, m, outw=16;
	double phi_sum=0.0, psi_max=0.0, res;
	double conv, conv_p, conv_px, conv_py, conv_jx, conv_jy;
	// residual data
	int ii=0, i_mbal=0, i_mx=0, i_my=0, i_mb=0, i_mt=0;
	int jj=0, j_mbal=0, j_mx=0, j_my=0, j_ml=0, j_mr=0;
	double res_mbal=0, res_mx=0, res_my=0, res_ml=0, res_mr=0, res_mb=0, res_mt=0;
	double res_bal=0, res_x=0, res_y=0, res_lbc=0, res_rbc=0, res_bbc=0, res_tbc=0;
	int i_bal=0, i_xx=0, i_yy=0, i_bbc=0, i_tbc=0;
	int j_bal=0, j_xx=0, j_yy=0, j_lbc=0, j_rbc=0;

	cout<<"Begin Output\n";
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// Matrix Residuals //////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// left BC residuals
	for (j=0; j<ny; j++) {
		// Left BC residuals
		res=abs(D_xP[0][j]*phi[0][j]/xe[0]+(FL[j]-D_xN[0][j]/xe[0])*phiB_L[j]);
		if ( res>res_ml ) { res_ml=res; j_ml=j; }
		// Right BC residuals
		res=abs(-D_xN[nx][j]*phi[nx-1][j]/xe[nx]+(FR[j]+D_xP[nx][j]/xe[nx])*phiB_R[j]);
		if ( res>res_mr ) { res_mr=res; j_mr=j; }
	}
	for (i=0; i<nx; i++) {
		// Bottom BC residuals
		res=abs(D_yP[i][0]*phi[i][0]/ye[0]+(FB[i]-D_yN[i][0]/ye[0])*phiB_B[i]);
		if ( res>res_mb ) { res_mb=res; i_mb=i; }
		// Top BC residuals
		res=abs(-D_yN[i][ny]*phi[i][ny-1]/ye[ny]+(FT[i]+D_yP[i][ny]/ye[ny])*phiB_T[i]);
		if ( res>res_mt ) { res_mt=res; i_mt=i; }
	}
	// balance residuals
	if ( nx==1 and ny==1 ) { // Single Cell
		i=0; j=0; // Balance Equations
		res=abs(-hy[j]*D_xP[i+1][j]*phiB_R[j]/xe[i+1]-hx[i]*D_yP[i][j+1]*phiB_T[i]/ye[j+1]+
		(hy[j]*D_xN[i+1][j]/xe[i+1]+hy[j]*D_xP[i][j]/xe[i]+hx[i]*D_yN[i][j+1]/ye[j+1]+hx[i]*D_yP[i][j]/ye[j]+
		hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]))*phi[i][j]
		-hy[j]*D_xN[i][j]*phiB_L[j]/xe[i]-hx[i]*D_yN[i][j]*phiB_B[i]/ye[j]-hx[i]*hy[j]*s_ext[i][j]);
		if ( res>res_mbal ) { res_mbal=res; i_mbal=i; j_mbal=j; }
	}
	else if ( nx==1 ) { // One Cell Wide
		i=0; j=0; // Bottom Cell
		res=abs(-hy[j]*D_xP[i+1][j]*phiB_R[j]/xe[i+1]-hx[i]*D_yP[i][j+1]*phi[i][j+1]/ye[j+1]+
		(hy[j]*D_xN[i+1][j]/xe[i+1]+hy[j]*D_xP[i][j]/xe[i]+hx[i]*D_yN[i][j+1]/ye[j+1]+hx[i]*D_yP[i][j]/ye[j]+
		hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]))*phi[i][j]
		-hy[j]*D_xN[i][j]*phiB_L[j]/xe[i]-hx[i]*D_yN[i][j]*phiB_B[i]/ye[j]-hx[i]*hy[j]*s_ext[i][j]);
		if ( res>res_mbal ) { res_mbal=res; i_mbal=i; j_mbal=j; }
		for ( j=1; j<ny-1; j++ ) { // Middle Cells
			res=abs(-hy[j]*D_xP[i+1][j]*phiB_R[j]/xe[i+1]-hx[i]*D_yP[i][j+1]*phi[i][j+1]/ye[j+1]+
			(hy[j]*D_xN[i+1][j]/xe[i+1]+hy[j]*D_xP[i][j]/xe[i]+hx[i]*D_yN[i][j+1]/ye[j+1]+hx[i]*D_yP[i][j]/ye[j]+
			hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]))*phi[i][j]
			-hy[j]*D_xN[i][j]*phiB_L[j]/xe[i]-hx[i]*D_yN[i][j]*phi[i][j-1]/ye[j]-hx[i]*hy[j]*s_ext[i][j]);
			if ( res>res_mbal ) { res_mbal=res; i_mbal=i; j_mbal=j; }
		}
		j=ny-1; // Top Cell
		res=abs(-hy[j]*D_xP[i+1][j]*phiB_R[j]/xe[i+1]-hx[i]*D_yP[i][j+1]*phiB_T[i]/ye[j+1]+
		(hy[j]*D_xN[i+1][j]/xe[i+1]+hy[j]*D_xP[i][j]/xe[i]+hx[i]*D_yN[i][j+1]/ye[j+1]+hx[i]*D_yP[i][j]/ye[j]+
		hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]))*phi[i][j]
		-hy[j]*D_xN[i][j]*phiB_L[j]/xe[i]-hx[i]*D_yN[i][j]*phi[i][j-1]/ye[j]-hx[i]*hy[j]*s_ext[i][j]);
		if ( res>res_mbal ) { res_mbal=res; i_mbal=i; j_mbal=j; }	
	}
	else if ( ny==1 ) { // One Cell Tall
		j=0; i=0; // Left Cell
		res=abs(-hy[j]*D_xP[i+1][j]*phi[i+1][j]/xe[i+1]-hx[i]*D_yP[i][j+1]*phiB_T[i]/ye[j+1]+
		(hy[j]*D_xN[i+1][j]/xe[i+1]+hy[j]*D_xP[i][j]/xe[i]+hx[i]*D_yN[i][j+1]/ye[j+1]+hx[i]*D_yP[i][j]/ye[j]+
		hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]))*phi[i][j]
		-hy[j]*D_xN[i][j]*phiB_L[j]/xe[i]-hx[i]*D_yN[i][j]*phiB_B[i]/ye[j]-hx[i]*hy[j]*s_ext[i][j]);
		if ( res>res_mbal ) { res_mbal=res; i_mbal=i; j_mbal=j; }
		for ( i=1; i<nx-1; i++ ) { // Middle Cell
			res=abs(-hy[j]*D_xP[i+1][j]*phi[i+1][j]/xe[i+1]-hx[i]*D_yP[i][j+1]*phiB_T[i]/ye[j+1]+
			(hy[j]*D_xN[i+1][j]/xe[i+1]+hy[j]*D_xP[i][j]/xe[i]+hx[i]*D_yN[i][j+1]/ye[j+1]+hx[i]*D_yP[i][j]/ye[j]+
			hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]))*phi[i][j]
			-hy[j]*D_xN[i][j]*phi[i-1][j]/xe[i]-hx[i]*D_yN[i][j]*phiB_B[i]/ye[j]-hx[i]*hy[j]*s_ext[i][j]);
			if ( res>res_mbal ) { res_mbal=res; i_mbal=i; j_mbal=j; }
		}
		i=nx-1; // Right Cell 
		res=abs(-hy[j]*D_xP[i+1][j]*phiB_R[j]/xe[i+1]-hx[i]*D_yP[i][j+1]*phiB_T[i]/ye[j+1]+
		(hy[j]*D_xN[i+1][j]/xe[i+1]+hy[j]*D_xP[i][j]/xe[i]+hx[i]*D_yN[i][j+1]/ye[j+1]+hx[i]*D_yP[i][j]/ye[j]+
		hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]))*phi[i][j]
		-hy[j]*D_xN[i][j]*phi[i-1][j]/xe[i]-hx[i]*D_yN[i][j]*phiB_B[i]/ye[j]-hx[i]*hy[j]*s_ext[i][j]);
		if ( res>res_mbal ) { res_mbal=res; i_mbal=i; j_mbal=j; }
		
	}
	else {
		i=0; j=0; // Bottom Left Corner
		res=abs(-hy[j]*D_xP[i+1][j]*phi[i+1][j]/xe[i+1]-hx[i]*D_yP[i][j+1]*phi[i][j+1]/ye[j+1]+
		(hy[j]*D_xN[i+1][j]/xe[i+1]+hy[j]*D_xP[i][j]/xe[i]+hx[i]*D_yN[i][j+1]/ye[j+1]+hx[i]*D_yP[i][j]/ye[j]+
		hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]))*phi[i][j]
		-hy[j]*D_xN[i][j]*phiB_L[j]/xe[i]-hx[i]*D_yN[i][j]*phiB_B[i]/ye[j]-hx[i]*hy[j]*s_ext[i][j]);
		if ( res>res_mbal ) { res_mbal=res; i_mbal=i; j_mbal=j; }
		for ( i=1; i<nx-1; i++ ) { // Bottom
			res=abs(-hy[j]*D_xP[i+1][j]*phi[i+1][j]/xe[i+1]-hx[i]*D_yP[i][j+1]*phi[i][j+1]/ye[j+1]+
			(hy[j]*D_xN[i+1][j]/xe[i+1]+hy[j]*D_xP[i][j]/xe[i]+hx[i]*D_yN[i][j+1]/ye[j+1]+hx[i]*D_yP[i][j]/ye[j]+
			hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]))*phi[i][j]
			-hy[j]*D_xN[i][j]*phi[i-1][j]/xe[i]-hx[i]*D_yN[i][j]*phiB_B[i]/ye[j]-hx[i]*hy[j]*s_ext[i][j]);
			if ( res>res_mbal ) { res_mbal=res; i_mbal=i; j_mbal=j; }
		}
		i=nx-1; // Bottom Right Corner
		res=abs(-hy[j]*D_xP[i+1][j]*phiB_R[j]/xe[i+1]-hx[i]*D_yP[i][j+1]*phi[i][j+1]/ye[j+1]+
		(hy[j]*D_xN[i+1][j]/xe[i+1]+hy[j]*D_xP[i][j]/xe[i]+hx[i]*D_yN[i][j+1]/ye[j+1]+hx[i]*D_yP[i][j]/ye[j]+
		hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]))*phi[i][j]
		-hy[j]*D_xN[i][j]*phi[i-1][j]/xe[i]-hx[i]*D_yN[i][j]*phiB_B[i]/ye[j]-hx[i]*hy[j]*s_ext[i][j]);
		if ( res>res_mbal ) { res_mbal=res; i_mbal=i; j_mbal=j; }
		for ( j=1; j<ny-1; j++ ) { // Middle
			i=0; // Left Side
			res=abs(-hy[j]*D_xP[i+1][j]*phi[i+1][j]/xe[i+1]-hx[i]*D_yP[i][j+1]*phi[i][j+1]/ye[j+1]+
			(hy[j]*D_xN[i+1][j]/xe[i+1]+hy[j]*D_xP[i][j]/xe[i]+hx[i]*D_yN[i][j+1]/ye[j+1]+hx[i]*D_yP[i][j]/ye[j]+
			hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]))*phi[i][j]
			-hy[j]*D_xN[i][j]*phiB_L[j]/xe[i]-hx[i]*D_yN[i][j]*phi[i][j-1]/ye[j]-hx[i]*hy[j]*s_ext[i][j]);
			if ( res>res_mbal ) { res_mbal=res; i_mbal=i; j_mbal=j; }
			for ( i=1; i<nx-1; i++ ) { // Centre Cells
				res=abs(-hy[j]*D_xP[i+1][j]*phi[i+1][j]/xe[i+1]-hx[i]*D_yP[i][j+1]*phi[i][j+1]/ye[j+1]+
				(hy[j]*D_xN[i+1][j]/xe[i+1]+hy[j]*D_xP[i][j]/xe[i]+hx[i]*D_yN[i][j+1]/ye[j+1]+hx[i]*D_yP[i][j]/ye[j]+
				hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]))*phi[i][j]
				-hy[j]*D_xN[i][j]*phi[i-1][j]/xe[i]-hx[i]*D_yN[i][j]*phi[i][j-1]/ye[j]-hx[i]*hy[j]*s_ext[i][j]);
				if ( res>res_mbal ) { res_mbal=res; i_mbal=i; j_mbal=j; }
			}
			i=nx-1; // Right Side
			res=abs(-hy[j]*D_xP[i+1][j]*phiB_R[j]/xe[i+1]-hx[i]*D_yP[i][j+1]*phi[i][j+1]/ye[j+1]+
			(hy[j]*D_xN[i+1][j]/xe[i+1]+hy[j]*D_xP[i][j]/xe[i]+hx[i]*D_yN[i][j+1]/ye[j+1]+hx[i]*D_yP[i][j]/ye[j]+
			hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]))*phi[i][j]
			-hy[j]*D_xN[i][j]*phi[i-1][j]/xe[i]-hx[i]*D_yN[i][j]*phi[i][j-1]/ye[j]-hx[i]*hy[j]*s_ext[i][j]);
			if ( res>res_mbal ) { res_mbal=res; i_mbal=i; j_mbal=j; }
		}
		j=ny-1; i=0; // Top Left Corner
		res=abs(-hy[j]*D_xP[i+1][j]*phi[i+1][j]/xe[i+1]-hx[i]*D_yP[i][j+1]*phiB_T[i]/ye[j+1]+
		(hy[j]*D_xN[i+1][j]/xe[i+1]+hy[j]*D_xP[i][j]/xe[i]+hx[i]*D_yN[i][j+1]/ye[j+1]+hx[i]*D_yP[i][j]/ye[j]+
		hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]))*phi[i][j]
		-hy[j]*D_xN[i][j]*phiB_L[j]/xe[i]-hx[i]*D_yN[i][j]*phi[i][j-1]/ye[j]-hx[i]*hy[j]*s_ext[i][j]);
		if ( res>res_mbal ) { res_mbal=res; i_mbal=i; j_mbal=j; }
		for ( i=1; i<nx-1; i++ ) { // Top
			res=abs(-hy[j]*D_xP[i+1][j]*phi[i+1][j]/xe[i+1]-hx[i]*D_yP[i][j+1]*phiB_T[i]/ye[j+1]+
			(hy[j]*D_xN[i+1][j]/xe[i+1]+hy[j]*D_xP[i][j]/xe[i]+hx[i]*D_yN[i][j+1]/ye[j+1]+hx[i]*D_yP[i][j]/ye[j]+
			hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]))*phi[i][j]
			-hy[j]*D_xN[i][j]*phi[i-1][j]/xe[i]-hx[i]*D_yN[i][j]*phi[i][j-1]/ye[j]-hx[i]*hy[j]*s_ext[i][j]);
			if ( res>res_mbal ) { res_mbal=res; i_mbal=i; j_mbal=j; }
		}
		i=nx-1; // Top Right Corner
		res=abs(-hy[j]*D_xP[i+1][j]*phiB_R[j]/xe[i+1]-hx[i]*D_yP[i][j+1]*phiB_T[i]/ye[j+1]+
		(hy[j]*D_xN[i+1][j]/xe[i+1]+hy[j]*D_xP[i][j]/xe[i]+hx[i]*D_yN[i][j+1]/ye[j+1]+hx[i]*D_yP[i][j]/ye[j]+
		hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j]))*phi[i][j]
		-hy[j]*D_xN[i][j]*phi[i-1][j]/xe[i]-hx[i]*D_yN[i][j]*phi[i][j-1]/ye[j]-hx[i]*hy[j]*s_ext[i][j]);
		if ( res>res_mbal ) { res_mbal=res; i_mbal=i; j_mbal=j; }
	}
	
	// calculate residuals ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Equation residuals ///////////////////////////////////////////////////////////////////////////////////
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) {
			// Balance in Cell
			res=abs(hy[j]*(j_x[i+1][j]-j_x[i][j])+hx[i]*(j_y[i][j+1]-j_y[i][j])+
			hx[i]*hy[j]*(sigmaT[i][j]-sigmaS[i][j]-nuF[i][j]*sigmaF[i][j])*phi[i][j]-hx[i]*hy[j]*s_ext[i][j]);
			if ( res>res_bal ) { res_bal=res; i_bal=i; j_bal=j; }
		}
	}
	for (i=1; i<nx-1; i++) {
		for (j=0; j<ny; j++) {
			// First Moment X Grid
			res=abs(j_x[i][j]+(D_xP[i][j]*phi[i][j]-D_xN[i][j]*phi[i-1][j])/xe[i]); if ( res>res_x ) { res_x=res; i_xx=i; j_xx=j; }
		}
	}
	for (i=0; i<nx; i++) {
		for (j=1; j<ny-1; j++) {
			// First Moment Y Grid
			res=abs(j_y[i][j]+(D_yP[i][j]*phi[i][j]-D_yN[i][j]*phi[i][j-1])/ye[j]); if ( res>res_y ) { res_y=res; i_yy=i; j_yy=j; }
		}
	}
	for (j=0; j<ny; j++) { // Left and Right
		res=abs(FL[j]*phiB_L[j]-j_x[0][j]);  if ( res>res_lbc ) { res_lbc=res; j_lbc=j; } // Left BC residual
		res=abs(FR[j]*phiB_R[j]-j_x[nx][j]); if ( res>res_rbc ) { res_rbc=res; j_rbc=j; } // Right BC residual
	}
	for (i=0; i<nx; i++) { // Bottom and Top
		res=abs(FB[i]*phiB_B[i]-j_y[i][0]);  if ( res>res_bbc ) { res_bbc=res; i_bbc=i; } // Bottom BC residual
		res=abs(FT[i]*phiB_T[i]-j_y[i][ny]); if ( res>res_tbc ) { res_tbc=res; i_tbc=i; } // Top BC residual
	}
	
	cout.precision(8);
	// find area averaged flux
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) phi_sum+=phiT[i][j]*hx[i]*hy[j];
	}
	cout<<"Area Averaged Flux :"<<print_out(phi_sum)<<"-------"<<endl;
	
	cout<<"output\n";
	// file output ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// boundary conditions
	outfile<<"\n -- Boundary Conditions -- \n";
	outfile<<"Type of BC: "<<kbc;   // Write type of BC
	
	// output grid
	outfile<<"\n -- Solution Grid -- \n";
	outfile<<" X grid \n"<<" Index     Cell Edge       Width Avg     Cell Centre      Cell Width   ";
	switch ( kbc ) {
	case 1:
		outfile<<"  Bottom BC In  "; // Write Bottom BC
		outfile<<"  Top    BC In  "; // Write Top BC
		break;
	case 2:
		outfile<<"  Quad 2 BC In  "; // Write Bottom BC
		outfile<<"  Quad 4 BC In  "; // Write Top BC
		break;
	case 3:
		outfile<<" Bottom BC REFL "; // Write Bottom BC
		outfile<<" Top    BC REFL "; // Write Top BC
		break;
	case 4:
		outfile<<" Bottom BC REFL "; // Write Bottom BC
		outfile<<"  Top    BC In  "; // Write Top BC
		break;
	case 5:
		outfile<<"  Bottom BC In  "; // Write Bottom BC
		outfile<<" Top    BC REFL "; // Write Top BC
		break;
	default:
		cout<<">>Invalid BC Type !!!!!\n";
		break;
	}
	outfile<<endl;
	for (i=0; i<nx; i++) outfile<<setw(6)<<i+1<<print_out(x[i])<<print_out(xe[i])<<print_out((x[i]+x[i+1])/2)<<print_out(hx[i])
	<<print_out(bcB[i])<<print_out(bcT[i])<<endl; // write x grid
	outfile<<setw(6)<<nx+1<<print_out(x[nx])<<print_out(xe[nx])<<endl;
	outfile<<" Y grid \n"<<" Index     Cell Edge       Width Avg     Cell Centre      Cell Width   ";
	switch ( kbc ) {
	case 1:
		outfile<<"  Left   BC In  "; // Write Left BC
		outfile<<"  Right  BC In  "; // Write Right BC
		break;
	case 2:
		outfile<<"  Quad 1 BC In  "; // Write Left BC
		outfile<<"  Quad 3 BC In  "; // Write Right BC
		break;
	case 3:
		outfile<<" Left   BC REFL "; // Write Left BC
		outfile<<" Right  BC REFL "; // Write Right BC
		break;
	case 4:
		outfile<<" Left   BC REFL "; // Write Left BC
		outfile<<"  Right  BC In  "; // Write Right BC
		break;
	case 5:
		outfile<<"  Left   BC In  "; // Write Left BC
		outfile<<" Right  BC REFL "; // Write Right BC
		break;
	default:
		cout<<">>Invalid BC Type !!!!!\n";
		break;
	}
	outfile<<endl;
	for (j=0; j<ny; j++) outfile<<setw(6)<<j+1<<print_out(y[j])<<print_out(ye[j])<<print_out((y[j]+y[j+1])/2)<<print_out(hy[j])
	<<print_out(bcL[j])<<print_out(bcR[j])<<endl; // write y grid
	outfile<<setw(6)<<ny+1<<print_out(y[ny])<<print_out(ye[ny])<<endl;
	
	// output quadrature
	outfile<<"\n -- Quadrature -- \n";
	if ( N==20 or N==36 ) outfile<<" Octant-Range Q";
	else outfile<<"Level Symmeteric S";
	outfile<<N<<" Normalized to Integrate to 4*pi \n";
	outfile<<"   m        mu              eta              xi            weight   \n";
	for (m=0; m<sn; m++) outfile<<setw(5)<<m+1<<print_out(mu[m])<<print_out(eta[m])<<print_out(xi[m])<<print_out(w[m])<<endl; // quadrature
	outfile<<endl;
	
	outfile<<endl;
	outfile<<" -------------- \n";
	outfile<<" -- Solution -- \n";
	outfile<<" -------------- \n";
	outfile<<"Program run time : "<<((float)t)/CLOCKS_PER_SEC<<" seconds"<<endl;
	outfile<<"\nNumber of iterations "<<n_iterations<<endl;
	
	outfile<<"\n -- Iteration Data -- \n";
	outfile<<"  Iter.   Convergence       # LO         LO Solution     High-order       Low-order    Preconditioner\n";
	outfile<<"    #        Rate        Iterations      Matrix Error   Sol. Time [s]   Sol. Time [s]  Sol. Time [s]\n";
	for (i=0; i<n_iterations; i++) outfile<<setw(4)<<i<<setw(2)<<" "<<print_out(rho_si[i])<<
		setw(outw/2)<<num_lo[i]+1<<setw(outw/2-1)<<" "<<print_out(err_lo[i])<<
		print_out(dt_ho[i])<<print_out(dt_lo[i])<<print_out(dt_pc[i])<<endl;
	
	
	outfile<<"\n -- Residuals -- \n";
	outfile<<"High-order Residual \n";
	outfile<<"Balance Residual:"<<print_out(max_res)<<" at "<<ii<<" , "<<jj<<endl;
	
	outfile<<"Low-order Residuals \n";
	outfile<<"Residual of Equations Solved In Matrix \n";
	outfile<<"Cell Balance Residual:"<<print_out(res_mbal)<<" at "<<i_mbal<<" , "<<j_mbal<<endl; 
	outfile<<"Left    BC   Residual:"<<print_out(res_ml)<<" at "<<j_ml<<endl;
	outfile<<"Right   BC   Residual:"<<print_out(res_mr)<<" at "<<j_mr<<endl;
	outfile<<"Bottom  BC   Residual:"<<print_out(res_mb)<<" at "<<i_mb<<endl;
	outfile<<"Top     BC   Residual:"<<print_out(res_mt)<<" at "<<i_mt<<endl;
	
	outfile<<"Residual of General Equations \n";
	outfile<<"Cell Balance Residual:"<<print_out(res_bal)<<" at "<<i_bal<<" , "<<j_bal<<endl; 
	outfile<<"X       Grid Residual:"<<print_out(res_x)<<" at "<<i_xx<<" , "<<j_xx<<endl;
	outfile<<"Y       Grid Residual:"<<print_out(res_y)<<" at "<<i_yy<<" , "<<j_yy<<endl;
	outfile<<"Left    BC   Residual:"<<print_out(res_lbc)<<" at "<<j_lbc<<endl;
	outfile<<"Right   BC   Residual:"<<print_out(res_rbc)<<" at "<<j_rbc<<endl;
	outfile<<"Bottom  BC   Residual:"<<print_out(res_bbc)<<" at "<<i_bbc<<endl;
	outfile<<"Top     BC   Residual:"<<print_out(res_tbc)<<" at "<<i_tbc<<endl;
	
	// Find difference between Transport solution and NDA solution
	conv_p=0.0;
	for ( i=0; i<nx; i++ ) {
		for ( j=0; j<ny; j++ ) {
			conv=abs(phi[i][j]-phiT[i][j]); if ( conv>conv_p ) conv_p=conv;
		}
	}
	conv_px=0.0;
	for ( j=0; j<ny; j++ ) {
		conv=abs(phiB_L[j]-phi_xT[0][j]);  if ( conv>conv_px ) conv_px=conv;
		conv=abs(phiB_R[j]-phi_xT[nx][j]); if ( conv>conv_px ) conv_px=conv;
	}
	conv_py=0.0;
	for ( i=0; i<nx; i++ ) {
		conv=abs(phiB_B[i]-phi_yT[i][0]);  if ( conv>conv_py ) conv_py=conv;
		conv=abs(phiB_T[i]-phi_yT[i][ny]); if ( conv>conv_py ) conv_py=conv;
	}
	conv_jx=0.0;
	for ( i=0; i<nx+1; i++ ) {
		for ( j=0; j<ny; j++ ) {
			conv=abs(j_x[i][j]-j_xT[i][j]); if ( conv>conv_jx ) conv_jx=conv;
		}
	}
	conv_jy=0.0;
	for ( i=0; i<nx; i++ ) {
		for ( j=0; j<ny+1; j++ ) {
			conv=abs(j_y[i][j]-j_yT[i][j]); if ( conv>conv_jy ) conv_jy=conv;
		}
	}
	
	outfile<<"\n -- Difference Between Transport and NDA Solution -- \n";
	outfile<<"Cell      Averaged Flux : "<<print_out(conv_p)<<endl;
	outfile<<"L and R   Boundary Flux : "<<print_out(conv_px)<<endl;
	outfile<<"B and T   Boundary Flux : "<<print_out(conv_py)<<endl;
	outfile<<"X Grid Face Avg Current : "<<print_out(conv_jx)<<endl;
	outfile<<"Y Grid Face Avg Current : "<<print_out(conv_jy)<<endl;
	
	// output flux
	outfile<<"\n ------------------- ";
	outfile<<"\n -- NDA Solution -- "; // Write NDA solution
	outfile<<"\n ------------------- \n";
	outfile<<"\n -- Cell Averaged Scalar Flux -- \n";
	write_cell_average_out(phi, outw, outfile); // call function to write out cell average scalar flux
	
	outfile<<"\n -- Left and Right Boundary Flux -- \n";
	outfile<<"  index    sol. grid     Left   Flux     Right  Flux  \n";
	for ( j=ny-1; j>=0; j-- ) outfile<<setw(6)<<j+1<<print_out((y[j]+y[j+1])/2)<<print_out(phiB_L[j])<<print_out(phiB_R[j])<<endl;
	
	outfile<<"\n -- Bottom and Top Boundary Flux -- \n";
	outfile<<"  index    sol. grid     Bottom  Flux     Top   Flux  \n";
	for ( i=nx-1; i>=0; i-- ) outfile<<setw(6)<<i+1<<print_out((x[i]+x[i+1])/2)<<print_out(phiB_B[i])<<print_out(phiB_T[i])<<endl;
	
	outfile<<"\n -- X Face Average Normal Current J_x -- \n";
	write_cell_edge_x_out(j_x, outw, outfile); // call function to write out cell edge current on x grid
	
	outfile<<"\n -- Y Face Average Normal Current J_y -- \n";
	write_cell_edge_y_out(j_y, outw, outfile); // call function to write out cell edge scalar flux on y grid
	//////////////////////////////////////////////////////////////////////////////////////////////////
	outfile<<"\n ------------------------- ";
	outfile<<"\n -- High Order Solution -- "; // Write Step Characteristics solution
	outfile<<"\n ------------------------- \n";
	outfile<<"\n -- Cell Averaged Scalar Flux -- \n";
	write_cell_average_out(phiT, outw, outfile); // call function to write out cell average scalar flux
	
	outfile<<"\n -- X Grid Face Averaged Scalar Flux -- \n";
	write_cell_edge_x_out(phi_xT, outw, outfile); // call function to write out cell edge scalar flux on x grid
	
	outfile<<"\n -- Y Grid Face Averaged Scalar Flux -- \n";
	write_cell_edge_y_out(phi_yT, outw, outfile); // call function to write out cell edge scalar flux on y grid
	
	outfile<<"\n -- X Face Average Normal Current J_x -- \n";
	write_cell_edge_x_out(j_xT, outw, outfile); // call function to write out cell edge current on x grid
	
	outfile<<"\n -- Y Face Average Normal Current J_y -- \n";
	write_cell_edge_y_out(j_yT, outw, outfile); // call function to write out cell edge scalar flux on y grid
	
	outfile<<" ---------------------------- \n";
	outfile<<" -- Diffusion Coefficients -- \n";
	outfile<<" ---------------------------- \n";
	
	outfile<<" -- Cell Average Diffusion Coefficient -- \n";
	write_cell_average_out(D, outw, outfile);
	
	outfile<<" -- X Grid Cell-edge Diffusion Coefficients -- \n";
	write_cell_edge_x_out(D_x, outw, outfile);
	
	outfile<<" -- Y Grid Cell-edge Diffusion Coefficients -- \n";
	write_cell_edge_y_out(D_y, outw, outfile);
	
	outfile<<"\n ----------------------- ";
	outfile<<"\n -- Consistency Terms -- "; // Write Step Characteristics solution
	outfile<<"\n ----------------------- \n";
	
	outfile<<"\n -- X Grid Consistency Term D_xTilde -- \n";
	write_cell_edge_x_out(D_xT, outw, outfile);
	
	outfile<<"\n -- Y Grid Consistency Term D_yTilde -- \n";
	write_cell_edge_y_out(D_yT, outw, outfile);
	
	outfile<<"\n -- X Grid Positive Diffusion Coefficient D^+_x -- \n";
	write_cell_edge_x_out(D_xP, outw, outfile);
	
	outfile<<"\n -- Y Grid Positive Diffusion Coefficient D^+_y -- \n";
	write_cell_edge_y_out(D_yP, outw, outfile);
	
	outfile<<"\n -- X Grid Negative Diffusion Coefficient D^-_x -- \n";
	write_cell_edge_x_out(D_xN, outw, outfile);
	
	outfile<<"\n -- Y Grid Negative Diffusion Coefficient D^-_y -- \n";
	write_cell_edge_y_out(D_yN, outw, outfile);
	
	// Boundary conditions
	outfile<<"\n -- Bottom and Top Boundary Factors -- \n";
	outfile<<" index "<<"     x centre    "<<"   F  Bottom    "<<"     F  Top     \n";
	for (i=0; i<nx; i++) outfile<<setw(6)<<i<<" "<<print_out((x[i]+x[i+1])/2)<<print_out(FB[i])<<print_out(FT[i])<<endl;
	
	outfile<<"\n -- Left and Right Boundary Factors -- \n";
	outfile<<" index "<<"     y centre    "<<"    F  Left     "<<"    F  Right    \n";
	for (j=0; j<ny; j++) outfile<<setw(6)<<i<<" "<<print_out((y[j]+y[j+1])/2)<<print_out(FL[j])<<print_out(FR[j])<<endl;
	
	outfile.close(); // close output file opened in input file +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	
	strcat (outf,fn); strcat (outf,".csv"); // name output file
	cout<<outf<<endl;
	datfile.open(outf); // open output data file ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// output flux
	datfile<<"# of iterations ,"<<n_iterations<<endl;
	datfile<<" -- Iteration Data -- \n";
	datfile<<"  Iter.,   Convergence,       # LO    ,   LO Solution ,  High-order   ,   Low-order   , Preconditioner\n";
	datfile<<"    #  ,      Rate    ,    Iterations ,   Matrix Error, Sol. Time [s] ,  Sol. Time [s], Sol. Time [s]\n";
	for (i=0; i<n_iterations; i++) datfile<<setw(4)<<i<<setw(2)<<" "<<","<<print_csv(rho_si[i])
		<<setw(outw/2)<<num_lo[i]+1<<setw(outw/2-1)<<" "<<","<<print_csv(err_lo[i])<<
		print_csv(dt_ho[i])<<print_csv(dt_lo[i])<<print_csv(dt_pc[i])<<endl;
	
	datfile<<" # of x cells , # of y cells ,\n";
	datfile<<nx<<" , "<<ny<<", \n";
	
	
	// output grid
	datfile<<" -- Solution Grid -- \n";
	datfile<<" X grid \n"<<" Index , Cell Edge , Cell Centre ,\n";
	for (i=0; i<nx; i++) datfile<<setw(6)<<i+1<<","<<print_csv(x[i])<<print_csv((x[i]+x[i+1])/2)<<endl; // write x grid
	datfile<<setw(6)<<nx+1<<","<<print_csv(x[nx])<<endl;
	
	datfile<<" Y grid \n"<<" Index , Cell Edge , Cell Centre ,\n";
	for (j=0; j<ny; j++) datfile<<setw(6)<<j+1<<","<<print_csv(y[j])<<print_csv((y[j]+y[j+1])/2)<<endl; // write y grid
	datfile<<setw(6)<<ny+1<<","<<print_csv(y[ny])<<endl;
	
	
	datfile<<" ------------------- \n";
	datfile<<" -- NDA Solution -- \n"; // Write NDA solution
	datfile<<" ------------------- \n";
	datfile<<" -- Cell Averaged Scalar Flux -- \n";
	write_cell_average_dat(phi, outw, datfile); // call function to write out cell average scalar flux
	
	datfile<<" -- Left and Right Boundary Flux -- \n";
	datfile<<" index,    sol. grid  ,  Left   Flux  ,  Right  Flux  ,\n";
	for ( j=ny-1; j>=0; j-- ) datfile<<setw(6)<<j+1<<","<<print_csv((y[j]+y[j+1])/2)<<print_csv(phiB_L[j])<<print_csv(phiB_R[j])<<endl;
	
	datfile<<" -- Bottom and Top Boundary Flux -- \n";
	datfile<<" index,    sol. grid  ,  Bottom  Flux ,   Top   Flux  ,\n";
	for ( i=nx-1; i>=0; i-- ) datfile<<setw(6)<<i+1<<","<<print_csv((x[i]+x[i+1])/2)<<print_csv(phiB_B[i])<<print_csv(phiB_T[i])<<endl;
	
	datfile<<" -- X Face Average Normal Current J_x -- \n";
	write_cell_edge_x_dat(j_x, outw, datfile); // call function to write out cell edge current on x grid
	
	datfile<<" -- Y Face Average Normal Current J_y -- \n";
	write_cell_edge_y_dat(j_y, outw, datfile); // call function to write out cell edge current on y grid
	
	datfile<<" ------------------------- \n";
	datfile<<" -- High Order Solution -- \n"; // Write Step Characteristics solution
	datfile<<" ------------------------- \n";
	datfile<<" -- Cell Averaged Scalar Flux -- \n";
	write_cell_average_dat(phiT, outw, datfile); // call function to write out cell average scalar flux
	
	datfile<<" -- X Vertical Cell Edge Scalar Flux -- \n";
	write_cell_edge_x_dat(phi_xT, outw, datfile); // call function to write out cell edge scalar flux on x grid
	
	datfile<<" -- Y Horizontal Cell Edge Scalar Flux -- \n";
	write_cell_edge_y_dat(phi_yT, outw, datfile); // call function to write out cell edge scalar flux on y grid
	
	datfile<<" -- X Face Average Normal Current J_x -- \n";
	write_cell_edge_x_dat(j_xT, outw, datfile); // call function to write out cell edge current on x grid
	
	datfile<<" -- Y Face Average Normal Current J_y -- \n";
	write_cell_edge_y_dat(j_yT, outw, datfile); // call function to write out cell edge current on y grid
	
	// Write D terms
	datfile<<" ---------------------------- \n";
	datfile<<" -- Diffusion Coefficients -- \n";
	datfile<<" ---------------------------- \n";
	
	datfile<<" -- Cell Average Diffusion Coefficient -- \n";
	write_cell_average_dat(D, outw, datfile);
	
	datfile<<" -- X Grid Cell-edge Diffusion Coefficients -- \n";
	write_cell_edge_x_dat(D_x, outw, datfile);
	
	datfile<<" -- Y Grid Cell-edge Diffusion Coefficients -- \n";
	write_cell_edge_y_dat(D_y, outw, datfile);
	
	datfile<<" ----------------------- \n";
	datfile<<" -- Consistency Terms -- \n"; // 
	datfile<<" ----------------------- \n";
	
	datfile<<" -- X Grid Consistency Term D_xTilde -- \n";
	write_cell_edge_x_dat(D_xT, outw, datfile);
	
	datfile<<" -- Y Grid Consistency Term D_yTilde -- \n";
	write_cell_edge_y_dat(D_yT, outw, datfile);
	
	datfile<<" -- X Grid Positive Diffusion Coefficient D^+_x -- \n";
	write_cell_edge_x_dat(D_xP, outw, datfile);
	
	datfile<<" -- Y Grid Positive Diffusion Coefficient D^+_y -- \n";
	write_cell_edge_y_dat(D_yP, outw, datfile);
	
	datfile<<" -- X Grid Negative Diffusion Coefficient D^-_x -- \n";
	write_cell_edge_x_dat(D_xN, outw, datfile);
	
	datfile<<" -- Y Grid Negative Diffusion Coefficient D^-_y -- \n";
	write_cell_edge_y_dat(D_yN, outw, datfile);
	
	// Boundary conditions
	datfile<<" -- Bottom and Top Boundary Factors -- \n";
	datfile<<" index,"<<"    x centre   ,"<<"    F Bottom   ,"<<"     F Top     ,\n";
	for (i=0; i<nx; i++) datfile<<setw(6)<<i<<","<<print_csv((x[i]+x[i+1])/2)<<print_csv(FB[i])<<print_csv(FT[i])<<endl;
	
	datfile<<" -- Left and Right Boundary Factors -- \n";
	datfile<<" index,"<<"    y centre   ,"<<"     F Left    ,"<<"    F Right    ,\n";
	for (j=0; j<ny; j++) datfile<<setw(6)<<i<<","<<print_csv((y[j]+y[j+1])/2)<<print_csv(FL[j])<<print_csv(FR[j])<<endl;
	
	datfile.close(); // close output data file ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
}
//======================================================================================//
