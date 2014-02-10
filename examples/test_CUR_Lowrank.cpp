/**
 \author        Sivaram Ambikasaran
 \date          February 1st, 2014
 \brief         Example file for CUR Low-rank routine.
 */

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "Matrix_Generator.hpp"
#include "CUR_Lowrank.hpp"

using std::cout;
using std::endl;

using namespace Eigen;

#ifdef EXACTLOWRANK
class Test_Matrix : public Matrix_Generator {
public:
        MatrixXd U;
        MatrixXd V;

        Test_Matrix(int M, int N, int p) {
                Matrix_Generator(M, N);
                U       =       MatrixXd::Random(M,p);
                V       =       MatrixXd::Random(p,N);
        }

        double get_Matrix_Entry(int i, int j) {
                return U.row(i)*V.col(j);
        }
};

#elif KERNEL
class Kernel_Matrix : public Matrix_Generator {
public:
        MatrixXd R;

        Kernel_Matrix(int M, int N) {
                Matrix_Generator(M, N);
                VectorXd x      =       VectorXd::Random(M);
                VectorXd y      =       4.0*VectorXd::Ones(N) + VectorXd::Random(N);

                std::sort(x.data(), x.data()+x.size());
                std::sort(y.data(), y.data()+y.size());

                R               =       MatrixXd(M,N);
                for (int j=0; j<M; ++j) {
                        for (int k=0; k<N; ++k) {
                                R(j,k)  = fabs(x(j)-y(k));
                        }
                }
        }

        double get_Matrix_Entry(int i, int j) {
                return log(R(i,j));
        }
};
#endif

#ifdef EQUISPACED
void get_Equi_Spaced_Index(const int N, const int n, vector<int>& index) {
        int spacing     =       N/n;
        for (int i=0; i<n; ++i) {
                index.push_back(i*spacing);
        }
}

#elif CHEBSPACED
void get_Cheb_Spaced_Index(const int N, const int n, vector<int>& index) {
        int temp;
        temp    =       N*(0.5+0.5*cos(M_PI/(2.0*n)));
        index.push_back(temp);
        for (int i=1; i<n; ++i) {
                temp    =       N*(0.5+0.5*cos((2.0*i+1.0)*M_PI/(2.0*n)));
                if (temp!=index.back()) {
                        index.push_back(temp);
                }
        }
}
#endif

int main() {
        srand(time(NULL));

        int M   =       4000;
        int N   =       4000;
        int n   =       50;

        vector<int> row_Index;
        vector<int> col_Index;

        #ifdef EQUISPACED
        get_Equi_Spaced_Index(M, n, row_Index);
        get_Equi_Spaced_Index(N, n, col_Index);

        #elif CHEBSPACED
        get_Cheb_Spaced_Index(M, n, row_Index);
        get_Cheb_Spaced_Index(N, n, col_Index);
        #endif

        double tolerance        =       1e-11;

        #ifdef EXACTLOWRANK
        Test_Matrix* matrix             =       new Test_Matrix(M, N, n);
        CUR_Lowrank<Test_Matrix>* LR    =       new CUR_Lowrank<Test_Matrix>(matrix, M, N, row_Index, col_Index, tolerance);

        #elif KERNEL
        Kernel_Matrix* matrix           =       new Kernel_Matrix(M, N);
        CUR_Lowrank<Kernel_Matrix>* LR  =       new CUR_Lowrank<Kernel_Matrix>(matrix, M, N, row_Index, col_Index, tolerance);

        #endif

        cout << endl << "Rank obtaining using CUR decomposition is: " << LR->rank << endl;

        MatrixXd A(M,N);
        for (int j=0; j<M; ++j) {
                for (int k=0; k<N; ++k) {
                        A(j,k) = matrix->get_Matrix_Entry(j, k);
                }
        }

        MatrixXd diff   =       A-LR->U*LR->S.asDiagonal()*LR->V.transpose();

        cout << endl << "True rank is: " << A.fullPivLu().rank() << endl;


        cout << endl << "Maximum error coefficient wise is: " << diff.cwiseAbs().maxCoeff() << endl;
        cout << endl << "Relative error based on 2-norm is: " << diff.norm()/A.norm() << endl;
}