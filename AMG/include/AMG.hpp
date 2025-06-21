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

constexpr double EPSILON = 0.2;


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
                
                if (max_value < std::abs(pair.second))
                    max_value =  std::abs(pair.second);
            }

            std::vector<size_t> strong_connections;

            for (const auto &pair : current_row)
            {
                if (row == pair.first)
                    continue;

                if (std::abs(pair.second) >= EPSILON * max_value)
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


        size_t select_coarse_nodes(CSRMatrix &current_matrix, std::vector<unsigned char> &coarse_mask)
        {
            const size_t n_nodes = current_matrix.rows();
            std::vector<std::vector<size_t>> strong_connections(coarse_mask.size());

            select_strong_connections(current_matrix, strong_connections, coarse_mask);

            size_t index = getRandomInit(n_nodes);
            //bool finished = false;

            size_t counter_fine = 0;

            while(coarse_mask.at(index) & 0x3F)     // I do this in order to use only a variable
            {
                coarse_mask.at(index) = 0;

                for (const auto &connection : strong_connections.at(index))
                {
                    if (coarse_mask.at(connection) & 0x3F)
                    {
                        coarse_mask.at(connection) |= 0xC0;     // Set state as fine
                        coarse_mask.at(connection) &= 0xC0;     // Set the number of strong con. to 0
                        counter_fine++;

                        for (const auto &second_connection : strong_connections.at(connection))
                        {
                            if (coarse_mask.at(second_connection) & 0x3F)
                            {
                                coarse_mask.at(second_connection) += 2;
                            }
                        }
                    }
                }

                signed char max = -1;
                for (size_t i = 0; i < n_nodes; ++i)
                {
                    if ((signed char)(coarse_mask.at(i) & 0x3F))
                    {
                        max = (signed char)(coarse_mask.at(i) & 0x3F);
                        index = i;
                    }
                }
            }

            //return the number of coarse nodes
            return coarse_mask.size() - counter_fine;

        }


        void build_component_mask
        (
            const std::vector<unsigned char> &coarse_mask, 
            std::vector<size_t> &component_mask,
            const size_t &num_coarse_nodes
        )
        {
            component_mask.resize(num_coarse_nodes);

            size_t index = 0;
            
            for (size_t i = 0; i < coarse_mask.size(); ++i)
            {
                // Check if the ith node is fine
                if (coarse_mask.at(i) & 0xC0)
                    continue;

                component_mask.at(index) = i; 
            }
        }

        void build_prolongation_matrix
        (
            CSRMatrix &current_matrix,
            std::unique_ptr<CSRMatrix> &P,
            std::vector<size_t> &component_mask,
            std::vector<unsigned char> &coarse_mask
        )
        {
            Matrix temporary_p_matrix(current_matrix.rows(), component_mask.size());
            
            for (size_t i = 0; i < current_matrix.rows(); ++i)
            {
                if (!(coarse_mask.at(i) & 0xC0))
                {
                    temporary_p_matrix.data().at(i)[i] = 1.0;
                    continue;
                }

                double alpha_num = 0.0;
                for (const auto &val : current_matrix.nonZerosInRow(i))
                {
                    if (val.first == i)
                        continue;

                    alpha_num += val.second;
                }

                double alpha_denum = 0.0;
                const auto strong_connections = strong_connections_in_row(current_matrix, i); 
                for (const auto &index : strong_connections)
                {
                    if (coarse_mask.at(index) & 0xC0)   // If the node is fine
                        continue;
                    
                    alpha_denum += current_matrix.coeff(i, index);
                }

                double alpha = alpha_num / alpha_denum;
                
                //double aii = current_matrix.coeff(i, i);
                double sum = 0.0;   //normalization constant

                for (const auto &index : strong_connections)
                {
                    if (coarse_mask.at(index) & 0xC0)
                        continue;


                    //temporary_p_matrix.data().at(i)[index] = 
                    //    alpha * current_matrix.coeff(i, index) / aii;
                    sum += alpha * current_matrix.coeff(i, index);
                }

                for (const auto &index : strong_connections)
                {
                    if (coarse_mask.at(index) & 0xC0)
                        continue;


                    temporary_p_matrix.data().at(i)[index] = 
                        alpha * current_matrix.coeff(i, index) / (-sum);  
                }

            }

            temporary_p_matrix.count_non_zeros();
            
            P = std::make_unique<CSRMatrix>(temporary_p_matrix);
            P->copy_from(temporary_p_matrix);
        }


        void build_coarse_matrix
        (
            CSRMatrix &current_matrix,
            CSRMatrix &P,
            std::unique_ptr<CSRMatrix> &coarse_matrix
        )
        {
            //Matrix PtA_temp(P.cols(), P.rows());
            //// Since A is symmetric cols of A are equal to cols of A
//
            //for (size_t j = 0; j < current_matrix.rows(); ++j)
            //{
            //    const auto &Aj_column = current_matrix.nonZerosInRow(j);
            //    for (size_t i = 0; i < P.cols(); ++i)
            //    {
            //        double sum = 0.0;
            //        
            //        for (const auto &non_zero_entry : Aj_column)
            //        {
            //            sum += non_zero_entry.second * P.coeff(non_zero_entry.first, i);
            //        }
//
            //        if (1e-10 < std::abs(sum))       // Chck if it's zero
            //            PtA_temp.data().at(i)[j] = sum;
            //    }
            //}

            //PtA_temp.count_non_zeros();

            //CSRMatrix PtA(PtA_temp);
            //PtA.copy_from(PtA_temp);

            //Matrix Ac(P.cols(), P.cols());
//
            //for (size_t i = 0; i < Ac.rows(); ++i)
            //{
            //    const auto &PtAi_row = PtA.nonZerosInRow(i);
//
            //    for (size_t j = 0; j < Ac.cols(); ++j)
            //    {
            //        double sum = 0.0;
//
            //        for (const auto &non_zero_entry : PtAi_row)
            //        {
            //            sum += non_zero_entry.second * P.coeff(non_zero_entry.first, j);
            //        }
            //        if (1e-10 < std::abs(sum))       // Chck if it's zero
            //            Ac.data().at(i)[j] = sum;
            //    }
            //}
            //Ac.count_non_zeros();

            //coarse_matrix = std::make_unique<CSRMatrix>(Ac);
            //coarse_matrix->copy_from(Ac);

        }
        

    private:
        std::vector<unsigned char> coarse_nodes;

        
};

#endif