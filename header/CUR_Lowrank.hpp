/**
 \class         CUR_Lowrank
 \author        Sivaram Ambikasaran
 \date          February 1st, 2014
 \brief         Obtains the low-rank decomposition by CUR method. Needs a set of row and column indices.
 */

#ifndef __CUR_Lowrank_hpp__
#define __CUR_Lowrank_hpp__

#include <Eigen/Dense>
#include <vector>
#include "Matrix_Generator.hpp"

using std::vector;
using std::cout;
using std::endl;
using namespace Eigen;

template <typename MatrixType>

class CUR_Lowrank {
public:
        MatrixType* matrix;     //      The matrix.
        int nRows;              //      Number of rows.
        int nCols;              //      Number of columns.
        MatrixXd U;             //      Unitary matrix, which forms the column basis.
        MatrixXd V;             //      Unitary matrix, which forms the row basis.
        VectorXd S;             //      The matrix in the middle.
        int rank;               //      Rank of the approximation.
        double tolerance;       //      Tolerance needed for inverting the sub-matrix.

        MatrixXd columns;       //      Desired columns.
        MatrixXd rows;          //      Desired rows.
        MatrixXd Ainv;          //      The sub-matrix, which is the intersection of columns and rows.

        /*!
         Obtains the pseudo-inverse of the matrix.
         */
        void pseudo_Inverse(MatrixXd A, double tolerance, MatrixXd& pinv) {
                JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
                VectorXd Sigma  =       svd.singularValues();
                int n           =       Sigma.size();
                VectorXd Sigma_Inverse  =       VectorXd::Zero(n);
                int k           =       0;
                while (k<n && Sigma(k) > tolerance) {
                        Sigma_Inverse(k)=       1.0/Sigma(k);
                        ++k;
                }
                pinv            =       svd.matrixV()*Sigma_Inverse.asDiagonal()*svd.matrixU().transpose();
        }

        /*!
         Obtains the QR decomposition as A=QR.
         */
        void get_QR(MatrixXd A, MatrixXd& Q, MatrixXd& R) {
                HouseholderQR<MatrixXd> qr(A);
                Q  =   qr.householderQ()*(MatrixXd::Identity(A.rows(),A.cols()));
                R  =   qr.matrixQR().block(0,0,A.cols(),A.cols()).triangularView<Upper>();
        }

        void get_SVD(MatrixXd A, double tolerance, MatrixXd& U, VectorXd& S, MatrixXd& V, int& rank) {
                JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
                rank    =       0;
                VectorXd Sigma  =       svd.singularValues();
                while (rank < Sigma.size() && Sigma(rank)>tolerance) {
                        ++rank;
                }
                U       =       svd.matrixU().block(0,0,A.rows(),rank);
                V       =       svd.matrixV().block(0,0,A.cols(),rank);
                S       =       Sigma.segment(0,rank);
        }

        /*!
         Constructs a low-rank representation from a given set of rows and columns as 'U S V^T', where 'U', 'V' are unitray matrices and 'S' is a vector containing the singular values.
         */
        CUR_Lowrank(MatrixType* matrix, int nRows, int nCols, vector<int> row_Index, vector<int> col_Index, double tolerance) {
                this->nRows             =       nRows;
                this->nCols             =       nCols;
                this->matrix            =       matrix;

                ///     Desired columns of the matrix.
                matrix->get_Matrix_Cols(0, nRows, col_Index, columns);


                ///     Desired rows of the matrix.
                matrix->get_Matrix_Rows(0, nCols, row_Index, rows);


                ///     The sub-matrix, which is the intersection of columns and rows.
                MatrixXd A;
                matrix->get_Matrix(row_Index, col_Index, A);


                //      A       =       UA*SA*VA'
                MatrixXd UA, VA;
                VectorXd SA;
                get_SVD(A, tolerance, UA, SA, VA, rank);

                //      SAinverse       =       SA^{-1}
                VectorXd SAinverse(rank);
                for (int i=0; i<rank; ++i) {
                        SAinverse(i)    =       1.0/SA(i);
                }

                //      Ainv    =       VA*SAinverse*UA'
                Ainv    =       VA*SAinverse.asDiagonal()*UA.transpose();
                
                ///     QR of the columns, i.e., columns = U*R1.
                MatrixXd R1;
                get_QR(columns, U, R1);

                ///     QR of the rows, i.e., rows.transpose() = V*R2.
                MatrixXd R2;
                get_QR(rows.transpose(), V, R2);

                //      A       =       R1*Ainv*R2'.
                A       =       R1*Ainv*R2.transpose();

                //      A       =       UA*S*VA'.
                get_SVD(A, tolerance, UA, S, VA, rank);

                //      U       =       U*UA.
                U       =       U*UA;

                //      V       =       V*VA.
                V       =       V*VA;
        }
};

#endif /*(__CUR_Lowrank_hpp__)*/