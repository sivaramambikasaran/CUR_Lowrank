/**
 \class         Matrix_Generator
 \author        Sivaram Ambikasaran
 \date          February 1st, 2014
 \brief         Matrix generator, which obtains matrix entries, rows, columns, sub-matrices.
 */

#ifndef __Matrix_Generator_hpp__
#define __Matrix_Generator_hpp__

#include <Eigen/Dense>
#include <vector>
#include <cstdlib>
#include <cassert>

using std::vector;
using namespace Eigen;

class Matrix_Generator {
public:
        int nRows;
        int nCols;

        Matrix_Generator() {};
        /*!
         Constructor assigning the number of rows and columns of the matrix.
         */
        Matrix_Generator(const int nRows, const int nCols) {
                this->nRows     =       nRows;
                this->nCols     =       nCols;
        }

        /*!
         Returns the matrix element at row 'i' and column 'j'.
         */
        virtual double get_Matrix_Entry(const int i, const int j) {
                assert(i>=0 && "Row index must be non-negative");
                assert(j>=0 && "Column index must be non-negative");
                
                assert(i<nRows && "Row index must be less than number of rows");
                assert(j<nCols && "Column index must be less than number of columns");
                return double(rand())/RAND_MAX;
        };

        /*!
         Allows access to the sub-matrix, i.e., returns a 'n_Rows' by 'n_Cols' sub-matrix starting at the index (start_Row, start_Column) of the matrix.
         */
	void get_Matrix(const int start_Row, const int start_Col, const int n_Rows, const int n_Cols, MatrixXd& A) {
		A   =   MatrixXd(n_Rows,n_Cols);
		for (int i=0; i<n_Rows; ++i) {
			for (int j=0; j<n_Cols; ++j) {
				A(i,j)  =   get_Matrix_Entry(start_Row+i, start_Col+j);
			}
		}
	};

        /*!
         Allows access to a staggered sub-matrix, with the elements having the desired row-index and column-index.
         */
	void get_Matrix(const vector<int> row_Index, const vector<int> col_Index, MatrixXd& A) {
                int n_Rows      =       row_Index.size();
                int n_Cols      =       col_Index.size();
		A               =       MatrixXd(n_Rows,n_Cols);
		for (int i=0; i<n_Rows; ++i) {
			for (int j=0; j<n_Cols; ++j) {
				A(i,j)  =   get_Matrix_Entry(row_Index[i], col_Index[j]);
			}
		}
	};

        /*!
         Allows access to the rows of the matrix, i.e., returns the 'row_Index'th row with the column starting at 'start_Col' and ending at 'start_Col+n_Cols'.
         */
	void get_Matrix_Row(const int start_Col, const int n_Cols, const int row_Index, VectorXd& v) {
		v   =   VectorXd(n_Cols);
		for (int j=0; j<n_Cols; ++j) {
			v(j)    =   get_Matrix_Entry(row_Index, start_Col+j);
		}
	}

        /*!
         Allows access to the rows of the matrix, i.e., returns the 'row_Index'th row with the column starting at 'start_Col' and ending at 'start_Col+n_Cols'.
         */
	void get_Matrix_Rows(const int start_Col, const int n_Cols, const vector<int> row_Index, MatrixXd& A) {
		A   =   MatrixXd(row_Index.size(), n_Cols);
                for (int i=0; i<row_Index.size(); ++i) {
                        for (int j=0; j<n_Cols; ++j) {
                                A(i,j)  =   get_Matrix_Entry(row_Index[i], start_Col+j);
                        }
                }
	}

        /*!
         Allows access to the columns of the matrix, i.e., returns the 'col_Index'th column with the row starting at 'start_Row' and ending at 'start_Row+n_Rows'.
         */
	void get_Matrix_Col(const int start_Row, const int n_Rows, const int col_Index, VectorXd& v) {
		v   =   VectorXd(n_Rows);
		for (int j=0; j<n_Rows; ++j) {
			v(j)    =   get_Matrix_Entry(start_Row+j, col_Index);
		}
	}
        
        /*!
         Allows access to the columns of the matrix, i.e., returns the 'col_Index'th column with the row starting at 'start_Row' and ending at 'start_Row+n_Rows'.
         */
	void get_Matrix_Cols(const int start_Row, const int n_Rows, const vector<int> col_Index, MatrixXd& A) {
		A   =   MatrixXd(n_Rows, col_Index.size());
                for (int i=0; i<n_Rows; ++i) {
                        for (int j=0; j<col_Index.size(); ++j) {
                                A(i,j)  =   get_Matrix_Entry(start_Row+i, col_Index[j]);
                        }
                }
	}
};

#endif/*(__Matrix_Generator_hpp__)*/