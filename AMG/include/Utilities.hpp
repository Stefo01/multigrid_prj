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
#include "CSRMatrix.hpp"

namespace Utilities
{
    typedef struct
    {
        double x;
        double y;
        bool is_on_boundary;
        size_t set_index;
    } Point;






    class TriangularMesh
    {
        public:
            TriangularMesh()
            {}

            void import_from_msh(const std::string &mesh_file_name);
            void export_to_vtu(const std::vector<double> &sol);

            const size_t n_nodes() const { return nodes.size(); }
            const size_t n_elements() const { return element_indexes.size(); }
            const size_t n_b_nodes() const { return num_boundary_nodes; }

            const std::vector< std::array<int, 3> > &get_elements_indexes() const { return element_indexes; }
            const std::vector<Point> &get_nodes(){ return nodes; };

        private:
            std::vector<Point> nodes;
            std::vector< std::array<int, 3> > element_indexes;
            size_t num_boundary_nodes;
    };


    

    template<class Vector>
    class SmootherClass{
        public:

            virtual void apply_iteration_to_vec(std::vector<double> &sol) = 0;

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
        public:

            Gauss_Seidel_iteration(CSRMatrix &A, Vector &f) : m_A(A), b(f) {}
            
            void apply_iteration_to_vec(std::vector<double> &sol) override{
                for(size_t i = 0; i < m_A.rows(); i++){
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
    };

    const double boundary_function(const double &x, const double &y)
    {
        /*
        const double k = 1;
        double r = sqrt(x * x + y * y);
        return - k * (cos(k * r) / r - k * sin(k * r));
        */
        //return x + y;
        return  std::sin(5 * std::sqrt(x * x + y * y));
        //return exp(x) * exp(-2 * y);
        //return 0.0;
    }

    const double forcing_term(const double &x, const double &y)
    {
        //return -5 * exp(x) * exp(-2 * y);
        //return 0.0;
        return -5 * ((std::cos(5 * std::sqrt(x * x + y * y)) / std::sqrt(x * x + y * y)) - 
            (5 * std::sin(5 * std::sqrt(x * x + y * y))));
    }

    const double alpha(const double &/*x*/, const double &/*y*/)
    {
        return 1.0;
    }


    int getRandomInit(int max){
        std::random_device rd;

        // Initialize a random number generator with the random device
        std::mt19937 gen(rd());
        // Create a uniform distribution within the specified range
        std::uniform_int_distribution<> distr(0, max);

        // Generate a random number within the range
        return distr(gen);
    }

    template <typename T>
    void printVector(std::vector<T> result){
        std::cout << "Vector elements: ";
        for (int i : result) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }

    
};

#endif