#ifndef CSR_MATRIX_HPP
#define CSR_MATRIX_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <map>
#include <fstream>
#include <string>
#include <numeric> // For std::inner_product
#include <random>  // For random number generation  
#include <unordered_set>
#include <ctime>
#include <cstdlib>
#include <algorithm>


class Matrix
{
    public:
        Matrix(const size_t &rows_, const size_t &cols_):
            n_rows(rows_),
            n_cols(cols_),
            num_non_zeros(0)
        {
            non_zero_cells.resize(rows_);
        }

        double &at(const size_t &row, const size_t &col);
        
        void count_non_zeros();
        const size_t non_zeros(){ return num_non_zeros; }
        const size_t rows(){ return n_rows; }
        const size_t cols(){ return n_cols; }
        const std::vector< std::map<size_t, double> > &data()
        { return non_zero_cells; }

        void print()
        {
            for (size_t i = 0; i < n_rows; ++i)
            {
                for (size_t j = 0; j < n_cols; ++j)
                {
                    std::cout << at(i, j) << "\t\t";
                }
                std::cout <<  std::endl;
            }
            std::cout << std::endl;
        }
        

    private:
        std::vector< std::map<size_t, double> > non_zero_cells;

        const size_t n_rows;
        const size_t n_cols;

        size_t num_non_zeros;
};

class CSRMatrix
{
    public:
        CSRMatrix(const size_t &rows_, const size_t &cols_, const size_t &nnz):
        n_rows(rows_), n_cols(cols_)
        {
            row_ptrs = new std::pair<size_t, double>*[rows_];
            non_zero_elements = new std::pair<size_t, double>[nnz];
            end = &non_zero_elements[nnz];
            std::vector<bool> xx(rows_,true);
            mask_amg.push_back(xx);        }

        CSRMatrix(Matrix&A):
        n_rows(A.rows()), n_cols(A.cols())
        {
            row_ptrs = new std::pair<size_t, double>*[A.rows()];
            non_zero_elements = new std::pair<size_t, double>[A.non_zeros()];
            end = &non_zero_elements[A.non_zeros()];
            std::vector<bool> xx(A.rows(),true);
            mask_amg.push_back(xx);
        }

        ~CSRMatrix()
        {
            delete[] row_ptrs;
            delete[] non_zero_elements;
        }



        void copy_from(Matrix &A);   
        const double coeff(const size_t &row, const size_t &col);
        const std::vector< std::pair<size_t, double> > nonZerosInRow(const size_t &row);
        const size_t rows(){ return mask_amg[mask_amg.size() - 1].size(); }
        const size_t cols(){ return mask_amg[mask_amg.size() - 1].size(); }

        void print()
        {
            for (size_t i = 0; i < n_rows; ++i)
            {
                for (size_t j = 0; j < n_cols; ++j)
                {
                    std::cout << coeff(i, j) << " ";
                }
                std::cout <<  std::endl;
            }
            std::cout << std::endl;
        }


    private:
        size_t n_rows;
        size_t n_cols;

        std::vector<std::vector<bool>> mask_amg;

        std::pair<size_t, double> **row_ptrs;
        std::pair<size_t, double> *non_zero_elements;

        std::pair<size_t, double> *end;

};


#endif