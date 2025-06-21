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
#include <memory>
#include "CSRMatrix.hpp"
#include "Utilities.hpp"

constexpr double EPSILON = 0.3;


// matrix A is the correspondence of the connections of a graph, aij represent the edje between node i and j. So, if the connection exist, 
// it will be represented by a value different from 0, otherwise the two nodes are disconnected.
// So the "primi vicini" are represented by all values in a specific row differents from zero
// all spd matrix are all square graph and the AMG method correspond to geometric multigrid


class AMG
{
    public:
        AMG(Matrix &A, std::vector<double> &soll, size_t number_of_levels_, std::vector<double> &rhs_): number_of_levels(number_of_levels_)
        {
            rhs.push_back(rhs_);
            //levels_matrix.push_back(&A); // Initialize the first level with the input matrix
            std::unique_ptr<CSRMatrix> tmp = std::make_unique<CSRMatrix>(A);
            tmp->copy_from(A);
            levels_matrix.push_back(std::move(tmp));
            x_levels.push_back(soll);
        }

        ~AMG()
        {

        }


        bool value_strong_connections(const size_t elementI, std::vector<bool> &Ret, int level, int &tot_strong_connections);
        double evaluate_node(std::vector<std::vector<double>> allNodes, std::vector<double> V, size_t elementI );
        int apply_AMG();
        int apply_smoother_operator(int level, int iter_number);
        void print_strong_connections(int level);
        void print_CSRmatrix(int level);
        void print_x_levels(int level);
        void print_mask_nodes(int level);
        bool is_connected(int i, int j, int level);
        double compute_residual(int level);

        std::vector<std::vector<bool>> get_strong_connections(int level) {
            return tot_strong_connections[level];
        }
        
        std::vector<double> get_x_levels(int level) {
            return x_levels[level];
        }

        std::vector<double> get_solution() {
            return x_levels[0];
        }

    private:
        
        int apply_restriction_operator(int level);
        int apply_prolungation_operator(int level);
        double compute_weight(int i, int j, int level);
        double compute_weight_real(int i, int j, int level);
        
        size_t number_of_levels;

        std::vector<std::vector<bool>>              mask_nodes;                 // for each level, we'll save the vector of choosen Course/Fine nodes
        std::vector<std::vector<std::vector<bool>>> tot_strong_connections;     // for each level, we'll save the matrix of strong connections
        std::vector<std::unique_ptr<CSRMatrix>>     levels_matrix;              // for each level, we'll save also the solution matrix                
        std::vector<std::vector<double>>            x_levels;                   // for each level, we'll save the solution vector
        std::vector<std::vector<double>>            rhs;
        std::vector<std::vector<std::vector<double>>> P_matrices;
};



class RestrictionOperator
{
    public:
        RestrictionOperator()
        {

        }

        ~RestrictionOperator()
        {

        }

        std::vector<size_t> strong_connections_in_row(CSRMatrix &current_matrix, size_t &row)
        {
            double max_value = 0.0;
            const auto &current_row = current_matrix.nonZerosInRow(row);
            for (const auto &pair : current_row)
            {
                if (row == pair.first)
                    continue;
                
                if (max_value < - pair.second)
                    max_value = - pair.second;
            }

            std::vector<size_t> strong_connections;

            for (const auto &pair : current_row)
            {
                if (- pair.second >= EPSILON * max_value)
                    strong_connections.push_back(pair.first);
            }

            return strong_connections;
        }

        void select_strong_connections
        (
            CSRMatrix &current_matrix, 
            std::vector<std::vector<size_t>> &strong_connections,
            std::vector<unsigned char> &coarse_mask
        )
        {
            for (size_t i = 0; i < current_matrix.rows(); ++i)
            {
                strong_connections.at(i) = 
                    strong_connections_in_row(current_matrix, i);

                coarse_mask.at(i) = static_cast<unsigned char>(strong_connections.at(i).size());

            }
        }


        void select_coarse_nodes(CSRMatrix &current_matrix, std::vector<unsigned char> &coarse_mask)
        {
            const size_t n_nodes = current_matrix.rows();
            std::vector<std::vector<size_t>> strong_connections(coarse_mask.size());

            select_strong_connections(current_matrix, strong_connections, coarse_mask);

            size_t index = getRandomInit(n_nodes);
            //bool finished = false;

            while(coarse_mask.at(index) & 0x0F)     // I do this in order to use only a variable
            {
                coarse_mask.at(index) = 0;

                for (const auto &connection : strong_connections.at(index))
                {
                    if (coarse_mask.at(connection) & 0x0F)
                    {
                        coarse_mask.at(connection) |= 0xF0;     // Set state as visited
                        coarse_mask.at(connection) &= 0xF0;     // Set the number of strong con. to 0

                        for (const auto &second_connection : strong_connections.at(connection))
                        {
                            if (coarse_mask.at(second_connection) & 0x0F)
                            {
                                coarse_mask.at(second_connection) += 1;
                            }
                        }
                    }
                }

                signed char max = -1;
                for (size_t i = 0; i < n_nodes; ++i)
                {
                    if ((signed char)(coarse_mask.at(i) & 0x0F))
                    {
                        max = (signed char)(coarse_mask.at(i) & 0x0F);
                        index = i;
                    }
                }
            }

        }
        

    private:
        std::vector<unsigned char> coarse_nodes;

        
};

#endif