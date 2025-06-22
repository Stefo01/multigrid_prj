#ifndef UTILITIES_HPP
#define UTILITIES_HPP

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
#include <memory>
#include <functional>
#include "CSRMatrix.hpp"


template<class Vector>
class SmootherClass{
    protected:
        std::function<void(std::vector<double> &)> apply_iteration_to_vec;

    public:

        inline friend std::vector<double>& operator*(std::vector<double> &x_k, SmootherClass &B)
        {
            B.apply_iteration_to_vec(x_k);
            return x_k;
        } // x^(k+1) = x^(k) * B

};

template<class Vector>
class Gauss_Seidel_iteration : public SmootherClass<Vector>{
    private:    
        CSRMatrix &m_A;
        Vector &b; // Ax = b


        void apply_iteration_to_vec_no_mask(std::vector<double> &sol)
        {
            for(size_t i = 0; i < m_A.rows(); i++)
            {
                double sum = 0;
                for (const auto &element : m_A.nonZerosInRow(i))
                {
                    if (element.first != i)
                    {
                        sum += element.second * sol[element.first];
                    }
                }
                sol[i] = (this->b[i] - sum) / m_A.coeff(i,i);
            }
        }

        void apply_iteration_to_vec_with_mask(std::vector<double> &sol)
        {
            for(size_t i = 0; i < m_A.rows(); i++)
            {
                size_t masked_i = m_A.component_mask.at(i);
                double sum = 0;
                for (const auto &element : m_A.nonZerosInRow(i))
                {
                    if (element.first != i)
                    {
                        sum += element.second * sol.at(m_A.component_mask.at(element.first));
                    }
                }
                sol[masked_i] = (this->b[masked_i] - sum) / m_A.coeff(i,i);
            }
        }

        

    public:

        Gauss_Seidel_iteration(CSRMatrix &A, Vector &f) : m_A(A), b(f) 
        {
            if (0 == A.component_mask.size())
                this->apply_iteration_to_vec = [this](std::vector<double> &sol)
                {
                    this->apply_iteration_to_vec_no_mask(sol);
                };
            
            else
                this->apply_iteration_to_vec = [this](std::vector<double> &sol)
                {
                    this->apply_iteration_to_vec_with_mask(sol);
                };

        }
        
        
};

const double boundary_function(const double &x, const double &y);

const double forcing_term(const double &x, const double &y);

const double alpha(const double &/*x*/, const double &/*y*/);

int getRandomInit(int max);

template <typename T>
void printVector(std::vector<T> result);


#endif