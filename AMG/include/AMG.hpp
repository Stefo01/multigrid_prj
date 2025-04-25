#ifndef AMG_HPP
#define AMG_HPP

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
#include "Utilities.hpp"

#define EPSILON 0.25


// matrix A is the correspondence of the connections of a graph, aij represent the edje between node i and j. So, if the connection exist, 
// it will be represented by a value different from 0, otherwise the two nodes are disconnected.
// So the "primi vicini" are represented by all values in a specific row differents from zero
// all spd matrix are all square graph and the AMG method correspond to geometric multigrid


class AMG
{
    public:
        AMG(CSRMatrix &A, std::vector<double> soll, size_t number_of_levels_, std::vector<double> rhs_): number_of_levels(number_of_levels_), rhs(rhs_)
        {
            size_t num_rows = A.rows();
            general_mask_for_rhs = std::vector<bool>(num_rows, true);

            // Initialize mask_nodes: vector of size 'number_of_levels', each with 'num_rows' set to false
            mask_nodes.resize(number_of_levels, std::vector<bool>(num_rows, false));

            // Initialize TOTStrongConnections: number_of_levels x num_rows x num_rows filled with false
            tot_strong_connections.resize(number_of_levels);
            
            for (size_t level = 0; level < number_of_levels; ++level) {
                tot_strong_connections[level].resize(num_rows, std::vector<bool>(num_rows, false));
            }
            levels_matrix.push_back(A); // Initialize the first level with the input matrix
            solution.push_back(soll);
        }

        ~AMG()
        {
        }


        void value_strong_connections(const size_t elementI, std::vector<bool> &Ret, int level);
        double evaluate_node(std::vector<std::vector<double>> allNodes, std::vector<double> V, size_t elementI );
        int apply_AMG();
        int apply_smoother_operator(int level, int iter_number);
        

    private:
    
        int apply_restriction_operator(int level);
        int apply_prolungation_operator(int level);
        
        size_t number_of_levels;

        std::vector<std::vector<bool>>              mask_nodes;                 // for each level, we'll save the vector of choosen Course/Fine nodes
        std::vector<std::vector<std::vector<bool>>> tot_strong_connections;     // for each level, we'll save the matrix of strong connections
        std::vector<CSRMatrix>                      levels_matrix;              // for each level, we'll save also the solution matrix                
        std::vector<std::vector<double>>            x_levels;                   // for each level, we'll save the solution vector
        std::vector<bool>                           general_mask_for_rhs;
        std::vector<double>                         rhs;
};

#endif