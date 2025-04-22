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

#define EPSILON 0.25


// matrix A is the correspondence of the connections of a graph, aij represent the edje between node i and j. So, if the connection exist, 
// it will be represented by a value different from 0, otherwise the two nodes are disconnected.
// So the "primi vicini" are represented by all values in a specific row differents from zero
// all spd matrix are all square graph and the AMG method correspond to geometric multigrid


class AMG
{
    public:
        AMG(CSRMatrix &A, size_t number_of_levels_): number_of_levels(number_of_levels_)
        {
            size_t num_rows = A.rows();

            // Initialize fine_course_nodes: vector of size 'number_of_levels', each with 'num_rows' set to false
            fine_course_nodes.resize(number_of_levels, std::vector<bool>(num_rows, false));

            // Initialize TOTStrongConnections: number_of_levels x num_rows x num_rows filled with false
            tot_strong_connections.resize(number_of_levels);
            
            for (size_t level = 0; level < number_of_levels; ++level) {
                tot_strong_connections[level].resize(num_rows, std::vector<bool>(num_rows, false));
            }
            levels_matrix.push_back(A); // Initialize the first level with the input matrix
        }

        ~AMG()
        {
        }


        void value_strong_connections(const size_t elementI, std::vector<int> &Ret, int level);
        double evaluate_node(std::vector<std::vector<double>> allNodes, std::vector<double> V, size_t elementI );
        int apply_AMG();
        int apply_restriction_operator(int level);
        int apply_prolungation_operator(int level);
        int apply_smoother_operator(int level);

    private:
        size_t number_of_levels;

        std::vector<std::vector<bool>>              fine_course_nodes;          // for each level, we'll save the vector of choosen Course/Fine nodes
        std::vector<std::vector<std::vector<bool>>> tot_strong_connections;     // for each level, we'll save the matrix of strong connections
        std::vector<CSRMatrix>                      levels_matrix;              // for each level, we'll save also the solution matrix
};

#endif