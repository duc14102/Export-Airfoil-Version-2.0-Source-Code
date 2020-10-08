#pragma once
#include "matrix.h"

// 'vectorMT' inherits from 'matrix' in a public way
class vectorMT : public matrix
{

public:
    //////////////////////// Constructors ////////////////////////////////////

    // default constructor
    vectorMT();
    // Constructor that takes 1 argument (size of the vectorMT)
    vectorMT(int no_of_elements);

    ////////////////// Binary Operators //////////////////////////////////////

    friend vectorMT operator+(const vectorMT& A, const vectorMT& B);
    friend vectorMT operator-(const vectorMT& A, const vectorMT& B);

    friend vectorMT operator*(const double& p, const vectorMT& A);
    friend vectorMT operator*(const vectorMT& A, const double& p);

    friend vectorMT operator/(const vectorMT& A, const double& p);

    ////////////////////// Unary operators //////////////////////////////////

    friend vectorMT operator+(const vectorMT& A);
    friend vectorMT operator-(const vectorMT& A);

    /////////////////////// Other operators /////////////////////////////////

    // Overloads (), so x(i) returns the ith entry a la MATLAB
    double& operator()(int i);

    // Overloads the assignment operator, '=' for vectorMT RHS
    vectorMT& operator=(const vectorMT& v);

    ////////////////////// Functions that are friends //////////////////////

    friend vectorMT mat2vec(matrix A);

    // Default call is norm(v) and returns 2−norm
    friend double norm(vectorMT v, int p);//p = 2

    friend vectorMT GMRES(matrix A, matrix b, double tol);//tol = 1e-6

    // friend vectorMT GMRESout(matrix A, matrix b, matrix x0, double tol = 1e-6);//tol = 1e-6

    //friend vectorMT resize(vectorMT v, int m);
    //friend vectorMT LeastSquares(matrix A, matrix b);
};
vectorMT LeastSquares(matrix A, matrix b);