// ioNDA1.1.h
#include <fstream>

double dparse(char[], int&, int);

int iparse(char[], int&, int);

int matRegion(int*, int*, int, int, int, int, double [], double [], double [], double []);

std::string print_out(double);

std::string print_csv(double);

void write_cell_average_out(double**, int, std::ofstream&); 
void write_cell_edge_x_out(double**, int, std::ofstream&);
void write_cell_edge_y_out(double**, int, std::ofstream&);
void write_cell_average_dat(double**, int, std::ofstream&);
void write_cell_edge_x_dat(double**, int, std::ofstream&);
void write_cell_edge_y_dat(double**, int, std::ofstream&);

