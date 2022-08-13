#include<iostream>
#include<iomanip>
#include"aphw1.h"

using Matrix = std::vector<std::vector<double>>;

Matrix multiply(Matrix& a, Matrix& b)
{
    //Getting shapes of matrices
    double r1{a.size()};
    double r2{b.size()};
    double c1{a[0].size()};
    double c2{b[0].size()};

    //Defining the product matrix 
    Matrix ans (r1 , std::vector<double>(c2,0));

    //Check if first matrix's columns is equal to second matrix's rows
    if(r2 != c1) 
        std::cout<<"First matrix columns and second matrix rows are not equal";
    else
    {    
        for (int i = 0; i < r1; i++)
            for (int j = 0; j < c2; j++)
                for (int k = 0; k < c1; k++)
                    ans[i][j] += a[i][k] * b[k][j];
    }

    return ans;
}

void show(Matrix& m1)
{
    double r{m1.size()};
    double c{m1[0].size()};

    for(int i=0 ; i<r; i++)
    {
        for (int j=0 ; j<c ; j++)
            std::cout<<std::setw(25)<<std::setprecision(18)<<m1[i][j];
        std::cout<<std::endl;
    }
}

Matrix transpose(Matrix& a)
{
    double r{a.size()};
    double c{a[0].size()};

    Matrix out (c, std::vector<double>(r , 0));

    for(int i=0 ; i<c; i++)
        for (int j=0 ; j<r ; j++)
            out[i][j] = a[j][i];

    return out;
}

Matrix getCofactor( Matrix A,Matrix x ,int p, int q, int n) 
{ 
    int i{0}, j{0}; 

    //Looping for each element of the matrix 
    for (int row = 0; row < n; row++) 
    { 
        for (int col = 0; col < n; col++) 
        { 
            //Copying into temporary matrix only those element which are not in given row and column 
            if (row != p && col != q) 
            { 
                x[i][j++] = A[row][col]; 

                //Row is filled, so increase row index and reset col index 
                if (j == n - 1) 
                { 
                    j = 0; 
                    i++; 
                } 
            } 
        } 
    }
    return x; 
}

double det(Matrix& A)
{
    //Getting the dimension
    int n = static_cast<int> (A.size());

    //Define the answer variable
    double ans{0};

    //Check if the matrix is one-element
    if(n == 1) 
        return A[0][0];

    //Define a matrix for cofactors
    Matrix x (n-1 , std::vector<double>(n-1 , 0));

    int sign{1};  //Variable for saving the sign

    //Iterate on first row elements
    for (int k = 0; k < n; k++)
    {
        //Getting cofactors
        x = getCofactor (A, x ,0, k, n);

        ans += sign * A[0][k] * det(x);
        sign = -sign;
    }

    return ans;
    
}


Matrix adjoint(Matrix& m)
{
    //Get the dimension
    double n{m.size()};

    //Define the adjoint matrix
    Matrix adj(n , std::vector<double>(n , 0));
    
    if (n == 1)
    {
        adj[0][0] = 1;
        return adj;
    }
 
    //Defining sign and temp for cofactors
    int sign = 1;
    Matrix temp(n-1 , std::vector<double>(n-1 ,0));

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            // Get cofactor of A[i][j]
            temp = getCofactor(m, temp, i, j, n);
 
            // sign of adj[j][i] positive if sum of row and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;
 
            // Interchanging rows and columns to get the transpose of the cofactor matrix
            adj[j][i] = (sign)*(det(temp));
        }
    }
    return adj;
}


Matrix inv(Matrix& m)
{
    //Find the dimension
    double d {m.size()};

    //Calculate determinant
    double det_value{det(m)};

    //Define inverse matrix
    Matrix ans(d, std::vector<double>(d , 0));
    Matrix adj(d, std::vector<double>(d , 0));
    //Find adjoint matrix
    adj = adjoint(m);

    //Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i=0; i<d; i++)
        for (int j=0; j<d; j++)
            ans[i][j] = adj[i][j]/det_value;

    return ans;
}


Matrix getData( const char* filename)
{
    Matrix output(10 , std::vector<double>(5 , 1));
    std::ifstream ifile{filename};

    long double *x{ new long double[10]} , *y { new long double[10]};
    char* c{new char [10]};
    
    if (ifile.fail())
        std::cout << "not found"<< std::endl;
    for(size_t i{}; i<10; i++)
    {
        ifile >> x[i]>> c[i]>> y[i];
        //std::cout << x[i]<< "\t" << y[i] << std::endl;
    }
    for(size_t i{0}; i<10 ; i++)
    {
        size_t j{1};
        output[i][j++] =x[i];
        output[i][j++] = sin(x[i]) ;
        output[i][j++] = sqrt(x[i]);
        output[i][j] =y[i];
    } 
    delete[] x;
    delete[] y;
    delete[] c;
    ifile.close();

    return output;
}


Matrix getX(Matrix& data)
{
    Matrix output (10 ,std::vector<double>(4 ,0));
    for (size_t i{}  ; i< 10 ; i++ )
        for (size_t j{}; j < 5; j++)
        {
            output [i][j] = data[i][j];
        }
    return output;
}
Matrix gety(Matrix& data)
{
    Matrix output (10 ,std::vector<double>(1 ,0));
    for (size_t i{}  ; i< 10 ; i++ )
        for (size_t j{}; j <1; j++)
        {
            output [i][j] = data[i][j+4];
        }
    return output;   
}

Matrix solve(char* filename)
{
    Matrix W (4 , std::vector<double>(10,0));
    Matrix data (10 , std::vector<double>(5,0));
    Matrix X (10 , std::vector<double>(4,0));
    Matrix y (10 , std::vector<double>(1,0));
    Matrix temp1 (4 , std::vector<double>(10,0));
    Matrix temp2 (5 , std::vector<double>(5,0));
    data = getData(filename);
    X = getX(data) ;
    y= gety(data);
    temp1 = transpose(X);
    temp2 = multiply(temp1 , X);
    temp2 = inv(temp2);
    W = multiply(temp2 , temp1);
    W = multiply(W , y);
    
    return W;
}

