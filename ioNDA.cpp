#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cstring>
#include <cstdio>

using namespace std;

extern int nx, ny;
extern double *x, *y;
extern ofstream outfile;
extern ofstream datfile;

// Function to parse Integer input
int iparse(char str[], int& i, int nl) {
	int p, t, num=0;
	while ( isdigit(str[i+1]) ) {
		num+=str[i]-'0';
		num*=10;
		i++;
	}
	num+=str[i]-'0';
	return num;
}

// Function to parse Decimal input
double dparse(char str[], int& i, int nl) {
	double num=0.0;
	int p, t, power;
	while ( isdigit(str[i+1]) ) {
		t=str[i]-'0';
		num+=t;
		num*=10.0;
		i++;
	}
	t=str[i]-'0';
	num+=t;
	if ( str[i+1]=='.' ) {
		i++; p=1;
		while ( isdigit(str[i+1]) ) {
			i++;
			t=str[i]-'0';
			num+=t/pow(10.0,double(p));
			p++;
		}
	}
	if ( str[i+1]=='e' ) {
		i++;
		if ( str[i+1]=='-' ) {
			i+=2;
			if ( isdigit(str[i]) ) power=iparse(str,i,nl);
			num/=pow(10.0,double(power));
		}
		else if ( str[i+1]=='+' ) {
			i+=2;
			if ( isdigit(str[i]) ) power=iparse(str,i,nl);
			num*=pow(10.0,double(power));
		}
		else if ( isdigit(str[i+1]) ) {
			i++;
			power=iparse(str,i,nl);
			num*=pow(10.0,double(power));
		}
	}
	return num;
}

// Function to find the material region
int matRegion(int *matnum, int *regnum, int i, int j, int nmat, int nreg, double xn[], double yn[], double xp[], double yp[]) {
	int k, m, reg, ret;
	double xc, yc;
	
	xc=(x[i+1]+x[i])/2.0;
	yc=(y[j+1]+y[j])/2.0;
	for (k=0; k<nreg; k++) {
		if ( xc>=xn[k] && xc<=xp[k] && yc>=yn[k] && yc<=yp[k] ) {
			reg=k;
			break;
		}
	}
	for (m=0; m<nmat; m++) {
		if ( regnum[reg]==matnum[m] ) {
			ret=m;
			break;
		}
	}
	return ret;
}

string print_out(double value)
{
	string dest;
	char buffer[50];
	
    // First print out using scientific notation with 0 mantissa digits
    sprintf(buffer, "%.8e", value);
	/*
    // Find the exponent and skip the "e" and the sign
    char *exponent = strchr(buffer, 'e') + 2;
	
    // If we have an exponent starting with 0, drop it
    if (exponent != NULL && exponent[0] == '0')
    {
        exponent[0] = exponent[1];
        exponent[1] = exponent[2];
		exponent[2] = '\0';
    }
	*/
	dest=buffer;
	if ( value>=0.0 ) dest=" " + dest;
	dest=" " + dest;
	return dest;
}

string print_csv(double value)
{
	string dest;
	char buffer[50];
	
    // First print out using scientific notation with 0 mantissa digits
    sprintf(buffer, "%.8e", value);
	/*
    // Find the exponent and skip the "e" and the sign
    char *exponent = strchr(buffer, 'e') + 2;
	
    // If we have an exponent starting with 0, drop it
    if (exponent != NULL && exponent[0] == '0')
    {
        exponent[0] = exponent[1];
        exponent[1] = exponent[2];
		exponent[2] = '\0';
    }
	*/
	dest=buffer;
	if ( value>=0.0 ) dest=" " + dest;
	dest=dest+",";
	return dest;
}

//======================================================================================//
//++ function to write cell average values in viewer friendly format +++++++++++++++++++//
//======================================================================================//
void write_cell_average_out(double **outp, int outw, ofstream& file) {
	int i, j, m;
	// cell average values
	for (m=0; m<int(nx/10); m++) {
		file<<setw(6)<<"index"<<setw(outw)<<" ";
		for (i=m*10; i<(m+1)*10; i++) file<<setw(outw)<<i+1;
		file<<"\n"<<setw(6)<<" "<<setw(outw)<<"sol. grid";
		for (i=m*10; i<(m+1)*10; i++) file<<print_out((x[i]+x[i+1])/2);
		file<<endl;
		for (j=ny-1; j>=0; j--) {
			file<<setw(6)<<j+1<<print_out((y[j]+y[j+1])/2);
			for (i=m*10; i<(m+1)*10; i++) file<<print_out(outp[i][j]);
			file<<endl;
		}
	}
	file<<setw(6)<<"index"<<setw(outw)<<" ";
	for (i=m*10; i<nx; i++) file<<setw(outw)<<i+1;
	file<<"\n"<<setw(6)<<" "<<setw(outw)<<"sol. grid";
	for (i=m*10; i<nx; i++) file<<print_out((x[i]+x[i+1])/2);
	file<<endl;
	for (j=ny-1; j>=0; j--) {
		file<<setw(6)<<j+1<<print_out((y[j]+y[j+1])/2);
		for (i=m*10; i<nx; i++) file<<print_out(outp[i][j]);
		file<<endl;
	}
}
//======================================================================================//
//++ function to write cell edge x values in viewer friendly format ++++++++++++++++++++//
//======================================================================================//
void write_cell_edge_x_out(double **outp, int outw, ofstream& file) {
	int i, j, m;
	// cell edge values on x grid
	for (m=0; m<int(nx/10); m++) {
		file<<setw(6)<<"index"<<setw(outw)<<" ";
		for (i=m*10; i<(m+1)*10; i++) file<<setw(outw)<<i+1;
		file<<"\n"<<setw(6)<<" "<<setw(outw)<<"sol. grid";
		for (i=m*10; i<(m+1)*10; i++) file<<print_out(x[i]);
		file<<endl;
		for (j=ny-1; j>=0; j--) {
			file<<setw(6)<<j+1<<print_out((y[j]+y[j+1])/2);
			for (i=m*10; i<(m+1)*10; i++) file<<print_out(outp[i][j]);
			file<<endl;
		}
	}
	file<<setw(6)<<"index"<<setw(outw)<<" ";
	for (i=m*10; i<nx+1; i++) file<<setw(outw)<<i+1;
	file<<"\n"<<setw(6)<<" "<<setw(outw)<<"sol. grid";
	for (i=m*10; i<nx+1; i++) file<<print_out(x[i]);
	file<<endl;
	for (j=ny-1; j>=0; j--) {
		file<<setw(6)<<j+1<<print_out((y[j]+y[j+1])/2);
		for (i=m*10; i<nx+1; i++) file<<print_out(outp[i][j]);
		file<<endl;
	}
}
//======================================================================================//
//++ function to write cell edge y values in viewer friendly format ++++++++++++++++++++//
//======================================================================================//
void write_cell_edge_y_out(double **outp, int outw, ofstream& file) {
	int i, j, m;
	// cell edge values on y grid
	for (m=0; m<int(nx/10); m++) {
		file<<setw(6)<<"index"<<setw(outw)<<" ";
		for (i=m*10; i<(m+1)*10; i++) file<<setw(outw)<<i+1;
		file<<"\n"<<setw(6)<<" "<<setw(outw)<<"sol. grid";
		for (i=m*10; i<(m+1)*10; i++) file<<print_out((x[i]+x[i+1])/2);
		file<<endl;
		for (j=ny; j>=0; j--) {
			file<<setw(6)<<j+1<<print_out(y[j]);
			for (i=m*10; i<(m+1)*10; i++) file<<print_out(outp[i][j]);
			file<<endl;
		}
	}
	file<<setw(6)<<"index"<<setw(outw)<<" ";
	for (i=m*10; i<nx; i++) file<<setw(outw)<<i+1;
	file<<"\n"<<setw(6)<<" "<<setw(outw)<<"sol. grid";
	for (i=m*10; i<nx; i++) file<<print_out((x[i]+x[i+1])/2);
	file<<endl;
	for (j=ny; j>=0; j--) {
		file<<setw(6)<<j+1<<print_out(y[j]);
		for (i=m*10; i<nx; i++) file<<print_out(outp[i][j]);
		file<<endl;
	}
}
//======================================================================================//
//++ function to write cell average values in data format ++++++++++++++++++++++++++++++//
//======================================================================================//
void write_cell_average_dat(double **outp, int outw, ofstream& file) {
	int i, j;
	// cell average values
	file<<setw(6)<<"index"<<","<<setw(outw)<<",";
	for (i=0; i<nx; i++) file<<setw(outw-1)<<i+1<<",";
	file<<"\n"<<setw(6)<<" "<<","<<setw(outw)<<" sol. grid,";
	for (i=0; i<nx; i++) file<<print_csv((x[i]+x[i+1])/2);
	file<<endl;
	for (j=0; j<ny; j++) {
		file<<setw(6)<<j+1<<","<<print_csv((y[j]+y[j+1])/2);
		for (i=0; i<nx; i++) file<<print_csv(outp[i][j]);
		file<<endl;
	}
}
//======================================================================================//
//++ function to write cell edge x values in data format +++++++++++++++++++++++++++++++//
//======================================================================================//
void write_cell_edge_x_dat(double **outp, int outw, ofstream& file) {
	int i, j;
	// cell edge values on x grid
	file<<setw(6)<<"index"<<","<<setw(outw)<<",";
	for (i=0; i<nx+1; i++) file<<setw(outw-1)<<i+1<<",";
	file<<"\n"<<setw(6)<<" "<<","<<setw(outw)<<" sol. grid,";
	for (i=0; i<nx+1; i++) file<<print_csv(x[i]);
	file<<endl;
	for (j=0; j<ny; j++) {
		file<<setw(6)<<j+1<<","<<print_csv((y[j]+y[j+1])/2);
		for (i=0; i<nx+1; i++) file<<print_csv(outp[i][j]);
		file<<endl;
	}
}
//======================================================================================//
//++ function to write cell edge y values in data format +++++++++++++++++++++++++++++++//
//======================================================================================//
void write_cell_edge_y_dat(double **outp, int outw, ofstream& file)
{
	int i, j;
	
	// cell edge values on y grid
	file<<setw(6)<<"index"<<","<<setw(outw)<<",";
	for (i=0; i<nx; i++) file<<setw(outw-1)<<i+1<<",";
	file<<"\n"<<setw(6)<<" "<<","<<setw(outw)<<" sol. grid,";
	for (i=0; i<nx; i++) file<<print_csv((x[i]+x[i+1])/2);
	file<<endl;
	for (j=0; j<ny+1; j++) {
		file<<setw(6)<<j+1<<","<<print_csv(y[j]);
		for (i=0; i<nx; i++) file<<print_csv(outp[i][j]);
		file<<endl;
	}
}

//======================================================================================//

