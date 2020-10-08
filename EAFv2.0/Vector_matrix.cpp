#include "vector_matrix.h"

//////////////////////// Constructors /////////////////////////////////////

// in the class, 'vectorMT' there is a constructor of the same name
vectorMT::vectorMT()

// runs the default matrix constructor
    : matrix()
{
}

vectorMT::vectorMT(int no_of_elements)

    : matrix(no_of_elements, 1)
{
}

/////////////////////// Other operators //////////////////////////////////

// Overloads (), so x(i) returns the ith entry a la MATLAB
double& vectorMT::operator()(int i)
{ // Can call or assign values.
    if (i < 1)
    {
        std::cout << "Error: Your index may be too small \n\n";
    }

    else if (i > rows)
    {
        std::cout << "Error: Your index may be too large \n\n";
    }
    return xx[i - 1][0];
}

// Operator returns a matrix equal to the RHS
vectorMT& vectorMT::operator=(const vectorMT& v)
{
    // Destruct previous entries
    for (int i = 0; i < rows; i++)
    {
        delete[] xx[i];
    }
    delete[] xx;

    // Assign new dimensions to be equal to that of the RHS
    rows = v.rows;
    columns = v.columns;

    // Allocate the memory as in the constructor
    xx = new double* [v.rows];

    for (int i = 0; i < v.rows; i++)
    {
        xx[i] = new double[v.columns];
    }

    // Copy the values across from the RHS
    for (int i = 0; i < v.rows; i++)
    {
        for (int j = 0; j < v.columns; j++)
        {
            // Set entries to be the same as the RHS matrix
            xx[i][j] = v.xx[i][j];
        }
    }
    return *this;
}
////////////////// Binary Operators ///////////////////////////////////////

vectorMT operator+(const vectorMT& A, const vectorMT& B)
{
    int m, n;
    m = rows(A);
    n = rows(B);

    if (m != n)
    {
        std::cout << "Error: Matrices of different dimensions.";
        std::cout << " Returned first argument";
        return A;
    }
    else
    {
        vectorMT v(m);
        for (int i = 0; i < m; i++)
        {
            v(i + 1) = A.xx[i][0] + B.xx[i][0];
        }
        return v;
    }
}

vectorMT operator-(const vectorMT& A, const vectorMT& B)
{
    int m, n;
    m = rows(A);
    n = rows(B);

    if (m != n)
    {
        std::cout << "Error: Matrices of different dimensions.";
        std::cout << " Returned first argument";
        return A;
    }
    else
    {
        vectorMT v(m);
        for (int i = 0; i < m; i++)
        {
            v(i + 1) = A.xx[i][0] - B.xx[i][0];
        }
        return v;
    }
}

//
vectorMT operator*(const double& p, const vectorMT& A)
{
    int m = rows(A);
    vectorMT v(m);
    for (int i = 0; i < m; i++)
    {
        v(i + 1) = p * A.xx[i][0];
    }
    return v;
}
//
vectorMT operator*(const vectorMT& A, const double& p)
{
    int m = rows(A);
    vectorMT v(m);
    for (int i = 0; i < m; i++)
    {
        v(i + 1) = p * A.xx[i][0];
    }
    return v;
}
//
vectorMT operator/(const vectorMT& A, const double& p)
{
    int m = rows(A);
    vectorMT v(m);
    for (int i = 0; i < m; i++)
    {
        v(i + 1) = A.xx[i][0] / p;
    }
    return v;
}
//

//////////////////////// Unary operators /////////////////////////////////

vectorMT operator+(const vectorMT& A)
{
    int m = rows(A);
    vectorMT v(m);
    for (int i = 0; i < m; i++)
    {
        v(i + 1) = A.xx[i][0];
    }
    return v;
}
//
vectorMT operator-(const vectorMT& A)
{
    int m = rows(A);
    vectorMT v(m);
    for (int i = 0; i < m; i++)
    {
        v(i + 1) = -A.xx[i][0];
    }
    return v;
}

////////////////////// Functions that are friends ////////////////////////

// Function that returns the first column of a matrix, A as a vectorMT

vectorMT mat2vec(matrix A)
{
    // create vectorMT with same no of rows as A and 1 column
    vectorMT v(rows(A));

    for (int i = 1; i <= rows(A); i++)
    {
        v(i) = A(i, 1); // copy only first column
    }

    return v;
}

double norm(vectorMT v, int p = 2)
{
    // p = 2;
    // define variables and initialise sum
    double temp, value, sum = 0.0;

    for (int i = 1; i <= rows(v); i++)
    {
        // floating point absolute value
        temp = fabs(v(i));
        sum += pow(temp, p);
    }

    value = pow(sum, 1.0 / ((double)(p)));

    return value;
}
// GMRES
vectorMT GMRES(matrix A, matrix b, double tol)
{
    //double tol=1e-6;
    vectorMT x0(rows(b));
    // determine initial residual, r0 in vectorMT form
    vectorMT r0 = mat2vec(b - A * x0);

    //std::cout << "initial residual vectorMT, r0 = b-A*x0 : \n\n" << r0;

    // need this in least square part later
    double normr0 = norm(r0, 2);

    // initialise to enter while loop
    double residual = 1.0;

    // intialise vectorMT v
    vectorMT v = r0 / normr0;

    //std::cout << "initial vectorMT v = r0 / || ro || : \n\n" << v;

    // Arnoldi/GMRES step index
    int k = 1;

    // Declare Givens rotation matrix, initialise Jtotal;
    matrix J, Jtotal(2, 2);
    // Jtotal = eye(2);
    for (int i = 1; i <= 2; i++)
    {
        Jtotal(i, i) = 1;
    }

    // intialise H, declare tempMat, V, w
    matrix H(1, 1), Htemp, HH, bb(1, 1), c, cc;
    matrix tempMat, V, Vold, hNewCol;
    vectorMT w, vj(rows(v));

    bb(1, 1) = normr0;

    // initialise matrix V (matrix of orthogonal basis vectorMTs)
    V = v;

    while (residual > tol)
    {
        //std::cout<< " \n\n";

        // update Vold (used for checking Arnoldi later)
        Vold = V;

        H = resize(H, k + 1, k);

        // Arnoldi steps (using Gram-Schmidt process)
        w = mat2vec(A * v);

        //std::cout<< "(k = " << k <<") : vectorMT w=Av : \n\n" << w;

        for (int j = 1; j <= k; j++)
        {
            for (int i = 1; i <= rows(V); i++)
            {

                // set the vectorMT vj to be jth column of V
                vj(i) = V(i, j);
            }

            tempMat = (~vj) * w;

            // these two lines calculate the inner product
            H(j, k) = tempMat(1, 1);
            //std::cout<< "H("<<j<<","<<k<<")= "<<H(j,k)<<"\n\n";

            w = w - H(j, k) * vj;
            //std::cout<< "Gramm-Schmidt update of vectorMT w: \n\n" << w;
        }

        H(k + 1, k) = norm(w, 2);
        //std::cout<< "H(" << k+1 << "," << k << ")= " << H(k+1,k) << "\n\n";

        v = w / H(k + 1, k);
        //std::cout<< "(k = " << k <<") :new vectorMT v: \n\n" << v;

        // add one more column to matrix V
        V = resize(V, rows(V), k + 1);

        for (int i = 1; i <= rows(V); i++)
        {
            // copy entries of v to new column of V
            V(i, k + 1) = v(i);
        }

        //std::cout<< "(k = " << k << ") :latest matrix, V: \n\n" << V;

        //std::cout << "(k = " << k <<") :latest matrix, H: \n\n" << H;

        //std::cout << "check: AV[k] = V[k+1]H: \n\n" << A*Vold << V*H;

        /////////////////////////////// Least squares step //////////////////////

        if (k == 1)
        {
            // First pass through, Htemp=H
            Htemp = H;
        }
        else
        {
            // for subsequent passes, Htemp=Jtotal*H
            Jtotal = resize(Jtotal, k + 1, k + 1);
            Jtotal(k + 1, k + 1) = 1;
            Htemp = Jtotal * H;
        }

        // Form next Givens rotation matrix
        // J = eye(k - 1);
        J(k - 1, k - 1);
        for (int i = 1; i <= k - 1; i++)
        {
            J(i, i) = 1;
        }

        J = resize(J, k + 1, k + 1);

        J(k, k) = Htemp(k, k) / pow(pow(Htemp(k, k), 2) + pow(Htemp(k + 1, k), 2), 0.5);
        J(k, k + 1) = Htemp(k + 1, k) / pow(pow(Htemp(k, k), 2) + pow(Htemp(k + 1, k), 2), 0.5);
        J(k + 1, k) = -Htemp(k + 1, k) / pow(pow(Htemp(k, k), 2) + pow(Htemp(k + 1, k), 2), 0.5);
        J(k + 1, k + 1) = Htemp(k, k) / pow(pow(Htemp(k, k), 2) + pow(Htemp(k + 1, k), 2), 0.5);

        //std::cout<< "J: \n\n" << J;

        // combine together with previous Givens rotations

        Jtotal = J * Jtotal;

        //std::cout<< "Check orthogonality of Jtotal \n\n" << ~Jtotal*Jtotal;

        HH = Jtotal * H;

        for (int i = 1; i <= k + 1; i++)
        {
            for (int j = 1; j <= k; j++)
            {
                // set all 'small' values to zero
                if (fabs(HH(i, j)) < 1e-15)
                {
                    HH(i, j) = 0;
                }
            }
        }

        //std::cout<< "Check Jtotal*H is upper triangular: \n\n" << HH;

        bb = resize(bb, k + 1, 1);

        //std::cout<< "bb: \n\n" << bb;

        c = Jtotal * bb;

        //std::cout<< "c=J*bb: \n\n" << c;

        residual = fabs(c(k + 1, 1));

        //std::cout<< k << "th residual: \n\n" << residual << "\n\n";

        k++;
    }

    // std::cout << "GMRES iteration converged in " << k - 1 << " steps\n\n";

    // Extract upper triangular square matrix
    HH = resize(HH, rows(HH) - 1, columns(HH));

    //std::cout<< "HH: \n\n" << HH;

    cc = resize(c, rows(HH), 1);

    //std::cout<< "cc: \n\n" << cc;

    matrix yy = cc / HH; // solve linear system

    vectorMT y = mat2vec(yy);

    //std::cout<< "y: \n\n" << y;

    // chop the newest column off of matrix V
    V = resize(V, rows(V), columns(V) - 1);

    vectorMT x = mat2vec(x0 + V * y);

    return x;
}
vectorMT LeastSquares(matrix A, matrix b)
{
    //double tol=1e-6;
    matrix x;
    x = inverse(~A * A) * (~A * b);
    vectorMT out = mat2vec(x);
    return out;
}