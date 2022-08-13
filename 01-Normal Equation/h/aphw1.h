#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
using Matrix = std::vector<std::vector<double>>;
Matrix multiply(Matrix& a, Matrix& b);
Matrix transpose(Matrix& a);
double det(Matrix& A);
Matrix inv(Matrix& a);
void show(Matrix& a); 

Matrix getData( const char* filename);
Matrix getX(Matrix& a);
Matrix gety(Matrix& data);
Matrix solve(char* filename);


