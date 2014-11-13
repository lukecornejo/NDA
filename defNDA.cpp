#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cmath>
#include <stdio.h>

using namespace std;
// file that contains definition functions for 2d_trans

extern const double pi=3.141592654;
extern int N;                // quadrature
extern double *mu, *eta, *xi, *w;
extern int nx, ny, sn;
extern double * x, * y, * hx, * hy, *xe, *ye;
extern double **sigmaT, **sigmaS, **sigmaF, **nuF, **s_ext;
extern int **material;
extern double **psi , **psi_x, **psi_y; // angular flux
extern double **phi , **j_x , **j_y; // NDA scalar flux and current solution
extern double **phiT, **phi_xT, **phi_yT, **j_xT, **j_yT; // scalar flux and current from transport
extern double *phiB_L, *phiB_R, *phiB_B, *phiB_T; // edge scalar flux from NDA
extern double **phiL; // scalar flux from last iteration
extern double ***psiL, ***psiR, ***psiB, ***psiT;
extern double *FL, *FR, *FB, *FT;
extern double **D, **D_x, **D_y; // Diffusion coefficients
extern double **D_xT, **D_yT, **D_xP, **D_xN, **D_yP, **D_yN; // Tilde, Positive and Negative

//======================================================================================//
//++ function to find the quadrature +++++++++++++++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void quadSet()
{
	double mup[N/2], c, wt[8], test=0.0;
	int i, j, m, nw, nt;
	int wi[N*(N+2)/8];
	double *wp[N*(N+2)/8];
	
	nt=N*(N+2)/8;
	
	if ( N<20 ) {
		switch ( N ) {
		case 2:
			nw=1;
			wi[0]=1;
			break;
		case 4: // S4 quadrature
			nw=1;
			wi[0]=1; wi[1]=1; 
				 wi[2]=1;
			mup[0]=0.3500212;
			wt[0]=0.333333333;
			break;
		case 6: // S6 quadrature
			nw=2;
			wi[0]=1; wi[1]=2; wi[2]=1; 
				 wi[3]=2; wi[4]=2; 
					wi[5]=1;
			mup[0]=0.2666355;
			wt[0]=0.1761263; wt[1]=0.1572071;
			break;
		case 8: // S8 quadrature
			nw=3;
			wi[0]=1; wi[1]=2; wi[2]=2; wi[3]=1; 
				  wi[4]=2; wi[5]=3; wi[6]=2; 
					  wi[7]=2; wi[8]=2; 
						  wi[9]=1;
			mup[0]=0.2182179;
			wt[0]=0.1209877; wt[1]=0.0907407; wt[2]=0.0925926;
			break;
		case 10:
			nw=4;
			wi[0]=1; wi[1]=2; wi[2]=3; wi[3]=2; wi[4]=1;
				  wi[5]=2; wi[6]=4; wi[7]=4; wi[8]=2; 
					  wi[9]=3; wi[10]=4; wi[11]=3; 
						  wi[12]=2; wi[13]=2; 
							  wi[14]=1;
			break;
		case 12: // S12 quadrature
			nw=5;
			wi[0]=1; wi[1]=2; wi[2]=3; wi[3]=3; wi[4]=2; wi[5]=1;
				 wi[6]=2; wi[7]=4; wi[8]=5; wi[9]=4; wi[10]=2;
					wi[11]=3; wi[12]=5; wi[13]=5; wi[14]=3;
						 wi[15]=3; wi[16]=4; wi[17]=3;
							 wi[18]=2; wi[19]=2;
								 wi[20]=1;
			mup[0]=0.1672126;
			wt[0]=0.0707626; wt[1]=0.0558811; wt[2]=0.0373377; wt[3]=0.0502819; wt[4]=0.0258513;
			break;
		case 14:
			nw=7;
			wi[0]=1; wi[1]=2; wi[2]=3; wi[3]=4; wi[4]=3; wi[5]=2; wi[6]=1;
			   wi[7]=2; wi[8]=5; wi[9]=6; wi[10]=6; wi[11]=5; wi[12]=2;
				   wi[13]=3; wi[14]=6; wi[15]=7; wi[16]=6; wi[17]=3;
					  wi[18]=4; wi[19]=6; wi[20]=6; wi[21]=4;
						   wi[22]=3; wi[23]=5; wi[24]=3;
							   wi[25]=2; wi[26]=2;
								   wi[27]=1;
			break;
		case 16: // S16 quadrature
			nw=8;
			wi[0]=1; wi[1]=2; wi[2]=3; wi[3]=4; wi[4]=4; wi[5]=3; wi[6]=2; wi[7]=1;
			  wi[8]=2; wi[9]=5; wi[10]=6; wi[11]=7; wi[12]=6; wi[13]=5; wi[14]=2;
				  wi[15]=3; wi[16]=6; wi[17]=8; wi[18]=8; wi[19]=6; wi[20]=3;
					 wi[21]=4; wi[22]=7; wi[23]=8; wi[24]=7; wi[25]=4;
						 wi[26]=4; wi[27]=6; wi[28]=6; wi[29]=4;
							 wi[30]=3; wi[31]=5; wi[32]=3;
								 wi[33]=2; wi[34]=2;
									 wi[35]=1;
			mup[0]=0.1389568;
			wt[0]=0.0489872; wt[1]=0.0413296; wt[2]=0.0212326; wt[3]=0.0256207; wt[4]=0.0360486; wt[5]=0.0144589; wt[6]=0.0344958; wt[7]=0.0085179;
			break;
		
		default:
			cout<<"invalid quadrature\n";
			break;
		}
		
		// find weights in quadrent
		for (i=0; i<nt; i++) wp[i]=wt+wi[i]-1; // pointers to angle weights
		
		c=2.0*(1-3*pow(mup[0],2))/(N-2);
		for (i=1; i<N/2; i++) mup[i]=sqrt(pow(mup[i-1],2)+c); // find cosines
		
		// quaderature
		sn=N*(N+2)/8; // total number of number combinations
		mu=new double[sn];
		eta=new double[sn];
		xi=new double[sn];
		w=new double[sn];
		
		// first quadrant
		m=0;
		for (i=N/2; i>0; i--) { // first quadrant
			for (j=0; j<i; j++) {
				mu[m] = mup[i-j-1]; // cosine on x  >
				eta[m]= mup[j];     // cosine on y  <
				xi[m] = mup[N/2-i]; // cosine on z  <
				w[m]=*wp[m]*pi;
				m++;
			}
		}
		
	}
	else {
		sn=N;
		mu=new double[sn];
		eta=new double[sn];
		xi=new double[sn];
		w=new double[sn];
		
		switch ( N ) {
		case 36: // Q36=q461214 quadrature
			
			w[ 0]=8.454511187252e-03; mu[ 0]=9.717784813336e-01; eta[ 0]=1.096881837272e-02; xi[ 0]=2.356401244281e-01;
			w[ 1]=1.913728513580e-02; mu[ 1]=9.701698603928e-01; eta[ 1]=5.695764868253e-02; xi[ 1]=2.356401244312e-01;
			w[ 2]=2.863542971348e-02; mu[ 2]=9.622473153642e-01; eta[ 2]=1.362124657777e-01; xi[ 2]=2.356401244295e-01;
			w[ 3]=3.648716160597e-02; mu[ 3]=9.410672772109e-01; eta[ 3]=2.426233944222e-01; xi[ 3]=2.356401244311e-01;
			w[ 4]=4.244873302980e-02; mu[ 4]=8.997294996538e-01; eta[ 4]=3.673697853806e-01; xi[ 4]=2.356401244316e-01;
			w[ 5]=4.642823955812e-02; mu[ 5]=8.335743322378e-01; eta[ 5]=4.996274255819e-01; xi[ 5]=2.356401244286e-01;
			w[ 6]=4.841339013884e-02; mu[ 6]=7.417637460141e-01; eta[ 6]=6.279014865859e-01; xi[ 6]=2.356401244320e-01;
			w[ 7]=4.841339013884e-02; mu[ 7]=6.279014865859e-01; eta[ 7]=7.417637460141e-01; xi[ 7]=2.356401244320e-01;
			w[ 8]=4.642823955812e-02; mu[ 8]=4.996274255819e-01; eta[ 8]=8.335743322378e-01; xi[ 8]=2.356401244286e-01;
			w[ 9]=4.244873302980e-02; mu[ 9]=3.673697853806e-01; eta[ 9]=8.997294996538e-01; xi[ 9]=2.356401244316e-01;
			w[10]=3.648716160597e-02; mu[10]=2.426233944222e-01; eta[10]=9.410672772109e-01; xi[10]=2.356401244311e-01;
			w[11]=2.863542971348e-02; mu[11]=1.362124657777e-01; eta[11]=9.622473153642e-01; xi[11]=2.356401244295e-01;
			w[12]=1.913728513580e-02; mu[12]=5.695764868253e-02; eta[12]=9.701698603928e-01; xi[12]=2.356401244312e-01;
			w[13]=8.454511187252e-03; mu[13]=1.096881837272e-02; eta[13]=9.717784813336e-01; xi[13]=2.356401244281e-01;
			w[14]=8.352354145856e-03; mu[14]=7.656319455497e-01; eta[14]=1.160393058611e-02; xi[14]=6.431742164832e-01;
			w[15]=1.873220073879e-02; mu[15]=7.633693960835e-01; eta[15]=5.995074957044e-02; xi[15]=6.431742164834e-01;
			w[16]=2.759429759588e-02; mu[16]=7.524467626583e-01; eta[16]=1.419535016004e-01; xi[16]=6.431742164829e-01;
			w[17]=3.442681426024e-02; mu[17]=7.241384940891e-01; eta[17]=2.488983098158e-01; xi[17]=6.431742164835e-01;
			w[18]=3.901232700510e-02; mu[18]=6.711819639118e-01; eta[18]=3.685670882907e-01; xi[18]=6.431742164829e-01;
			w[19]=4.130171453748e-02; mu[19]=5.909368760506e-01; eta[19]=4.869502395267e-01; xi[19]=6.431742164829e-01;
			w[20]=4.130171453748e-02; mu[20]=4.869502395267e-01; eta[20]=5.909368760506e-01; xi[20]=6.431742164829e-01;
			w[21]=3.901232700510e-02; mu[21]=3.685670882907e-01; eta[21]=6.711819639118e-01; xi[21]=6.431742164829e-01;
			w[22]=3.442681426024e-02; mu[22]=2.488983098158e-01; eta[22]=7.241384940891e-01; xi[22]=6.431742164835e-01;
			w[23]=2.759429759588e-02; mu[23]=1.419535016004e-01; eta[23]=7.524467626583e-01; xi[23]=6.431742164829e-01;
			w[24]=1.873220073879e-02; mu[24]=5.995074957044e-02; eta[24]=7.633693960835e-01; xi[24]=6.431742164834e-01;
			w[25]=8.352354145856e-03; mu[25]=1.160393058611e-02; eta[25]=7.656319455497e-01; xi[25]=6.431742164832e-01;
			w[26]=1.460888798152e-02; mu[26]=4.445439440056e-01; eta[26]=2.447911451942e-02; xi[26]=8.954225007226e-01;
			w[27]=2.995376809966e-02; mu[27]=4.288508824476e-01; eta[27]=1.196054590036e-01; xi[27]=8.954225007227e-01;
			w[28]=3.798783310581e-02; mu[28]=3.670788892962e-01; eta[28]=2.519357740235e-01; xi[28]=8.954225007226e-01;
			w[29]=3.798783310581e-02; mu[29]=2.519357740235e-01; eta[29]=3.670788892962e-01; xi[29]=8.954225007226e-01;
			w[30]=2.995376809966e-02; mu[30]=1.196054590036e-01; eta[30]=4.288508824476e-01; xi[30]=8.954225007227e-01;
			w[31]=1.460888798152e-02; mu[31]=2.447911451942e-02; eta[31]=4.445439440056e-01; xi[31]=8.954225007226e-01;
			w[32]=6.404244616724e-03; mu[32]=1.483114568272e-01; eta[32]=1.670387759191e-02; xi[32]=9.887996218887e-01;
			w[33]=1.162080754372e-02; mu[33]=1.293388490485e-01; eta[33]=7.447663982495e-02; xi[33]=9.887996218887e-01;
			w[34]=1.162080754372e-02; mu[34]=7.447663982495e-02; eta[34]=1.293388490485e-01; xi[34]=9.887996218887e-01;
			w[35]=6.404244616724e-03; mu[35]=1.670387759191e-02; eta[35]=1.483114568272e-01; xi[35]=9.887996218887e-01;
			
			break;
		case 20: // Q20=q2468 quadrature
			
			w[ 0]=2.419260514149E-02; mu[ 0]=9.713274064903E-01; eta[ 0]=3.157215799340E-02; xi[ 0]=2.356401244281E-01;
			w[ 1]=5.213067212540E-02; mu[ 1]=9.586898685237E-01; eta[ 1]=1.593344524838E-01; xi[ 1]=2.356401244307E-01;
			w[ 2]=7.185542471164E-02; mu[ 2]=9.028558915298E-01; eta[ 2]=3.596178122512E-01; xi[ 2]=2.356401244304E-01;
			w[ 3]=8.182604839076E-02; mu[ 3]=7.770210099715E-01; eta[ 3]=5.837054752370E-01; xi[ 3]=2.356401244296E-01;
			w[ 4]=8.182604839076E-02; mu[ 4]=5.837054752370E-01; eta[ 4]=7.770210099715E-01; xi[ 4]=2.356401244296E-01;
			w[ 5]=7.185542471164E-02; mu[ 5]=3.596178122512E-01; eta[ 5]=9.028558915298E-01; xi[ 5]=2.356401244304E-01;
			w[ 6]=5.213067212540E-02; mu[ 6]=1.593344524838E-01; eta[ 6]=9.586898685237E-01; xi[ 6]=2.356401244307E-01;
			w[ 7]=2.419260514149E-02; mu[ 7]=3.157215799340E-02; eta[ 7]=9.713274064903E-01; xi[ 7]=2.356401244281E-01;
			w[ 8]=2.998205782366E-02; mu[ 8]=7.645615896150E-01; eta[ 8]=4.210110375297E-02; xi[ 8]=6.431742164827E-01;
			w[ 9]=6.147460425028E-02; mu[ 9]=7.375714298063E-01; eta[ 9]=2.057068622698E-01; xi[ 9]=6.431742164831E-01;
			w[10]=7.796304620960E-02; mu[10]=6.313311043797E-01; eta[10]=4.332989313333E-01; xi[10]=6.431742164827E-01;
			w[11]=7.796304620960E-02; mu[11]=4.332989313333E-01; eta[11]=6.313311043797E-01; xi[11]=6.431742164827E-01;
			w[12]=6.147460425028E-02; mu[12]=2.057068622698E-01; eta[12]=7.375714298063E-01; xi[12]=6.431742164831E-01;
			w[13]=2.998205782366E-02; mu[13]=4.210110375297E-02; eta[13]=7.645615896150E-01; xi[13]=6.431742164827E-01;
			w[14]=2.932993043666E-02; mu[14]=4.424202396002E-01; eta[14]=4.982847370367E-02; xi[14]=8.954225007227E-01;
			w[15]=5.322055875020E-02; mu[15]=3.858240341629E-01; eta[15]=2.221674140412E-01; xi[15]=8.954225007227E-01;
			w[16]=5.322055875020E-02; mu[16]=2.221674140412E-01; eta[16]=3.858240341629E-01; xi[16]=8.954225007227E-01;
			w[17]=2.932993043666E-02; mu[17]=4.982847370367E-02; eta[17]=4.424202396002E-01; xi[17]=8.954225007227E-01;
			w[18]=1.802505216045E-02; mu[18]=1.409476441875E-01; eta[18]=4.908227124734E-02; xi[18]=9.887996218887E-01;
			w[19]=1.802505216045E-02; mu[19]=4.908227124734E-02; eta[19]=1.409476441875E-01; xi[19]=9.887996218887E-01;
			
			break;
		default:
			cout<<"invalid quadrature\n";
			break;
		}
		
		// first quadrant
		for (m=0; m<N; m++) { // first quadrant
			w[m]=w[m]*pi;
		}
		
	}
	
	for (m=0; m<sn; m++) test+=w[m];
	cout<<"weight test "<<test<<endl;
	// normalize weight
	for (m=0; m<sn; m++) w[m]*=pi/test;
	
}
//======================================================================================//

//======================================================================================//
//++ function to initialize problem memory space +++++++++++++++++++++++++++++++++++++++//
//======================================================================================//
void initialize()
{
	int i, j, m;
	// initialize memory space
	
	FL=new double[ny];
	FR=new double[ny];
	FB=new double[nx];
	FT=new double[nx];
	
	phiB_L=new double[ny];
	phiB_R=new double[ny];
	phiB_B=new double[nx];
	phiB_T=new double[nx];
	
	// initialize boundary conditions
	psiB=new double**[nx]; 
	psiT=new double**[nx];
	psiL=new double**[ny];
	psiR=new double**[ny];
	for (i=0; i<nx; i++) {
		psiB[i]=new double*[sn];
		psiT[i]=new double*[sn];
		for (m=0; m<sn; m++) {
			psiB[i][m]=new double[4];
			psiT[i][m]=new double[4];
		}
	}
	for (j=0; j<ny; j++) {
		psiL[j]=new double*[sn];
		psiR[j]=new double*[sn];
		for (m=0; m<sn; m++) {
			psiL[j][m]=new double[4];
			psiR[j][m]=new double[4];
		}
	}
	
	
	// cell average flux
	psi =new double*[nx];
	phiT=new double*[nx];
	phiL=new double*[nx];
	phi =new double*[nx];
	D   =new double*[nx];
	for (i=0; i<nx; i++) {
		psi[i] =new double[ny];
		phiT[i]=new double[ny];
		phiL[i]=new double[ny];
		phi[i] =new double[ny];
		D[i]   =new double[ny];
	}
	// cell edge values on x grid
	psi_x =new double*[nx+1];
	phi_xT=new double*[nx+1];
	j_x   =new double*[nx+1];
	j_xT  =new double*[nx+1];
	D_x   =new double*[nx+1];
	D_xT  =new double*[nx+1];
	D_xP  =new double*[nx+1];
	D_xN  =new double*[nx+1];
	for (i=0; i<nx+1; i++) {
		psi_x[i] =new double[ny];
		phi_xT[i]=new double[ny];
		j_x[i]   =new double[ny];
		j_xT[i]  =new double[ny];
		D_x[i]   =new double[ny];
		D_xT[i]  =new double[ny];
		D_xP[i]  =new double[ny];
		D_xN[i]  =new double[ny];
	}
	// cell edge flux on x grid
	psi_y =new double*[nx];
	phi_yT=new double*[nx];
	j_y   =new double*[nx];
	j_yT  =new double*[nx];
	D_y   =new double*[nx];
	D_yT  =new double*[nx];
	D_yP  =new double*[nx];
	D_yN  =new double*[nx];
	for (i=0; i<nx; i++) {
		psi_y[i] =new double[ny+1];
		phi_yT[i]=new double[ny+1];
		j_y[i]   =new double[ny+1];
		j_yT[i]  =new double[ny+1];
		D_y[i]   =new double[ny+1];
		D_yT[i]  =new double[ny+1];
		D_yP[i]  =new double[ny+1];
		D_yN[i]  =new double[ny+1];
	}
	// find values for D's
	for ( i=0; i<nx; i++ ) {
		for ( j=0; j<ny; j++ ) D[i][j]=1/(3*sigmaT[i][j]); // Diffusion Coefficient. Isotropic scattering so sigma_tr = sigma_t
	}
	for ( i=1; i<nx; i++ ) {
		for ( j=0; j<ny; j++ ) D_x[i][j]=2*xe[i]*D[i-1][j]*D[i][j]/(D[i-1][j]*hx[i]+D[i][j]*hx[i-1]); // X Grid D
	}
	for ( j=0; j<ny; j++ ) { // Top and Bottom
		D_x[0][j] =2*xe[0] *D[0][j]   /hx[0];
		D_x[nx][j]=2*xe[nx]*D[nx-1][j]/hx[nx-1];
	}
	for ( i=0; i<nx; i++ ) {
		for ( j=1; j<ny; j++ ) D_y[i][j]=2*ye[j]*D[i][j-1]*D[i][j]/(D[i][j-1]*hy[j]+D[i][j]*hy[j-1]); // Y Grid D
	}
	for ( i=0; i<nx; i++ ) { // Left and Right
		D_y[i][0] =2*ye[0] *D[i][0]   /hy[0];
		D_y[i][ny]=2*ye[ny]*D[i][ny-1]/hy[ny-1];
	}
	
}

//======================================================================================//

