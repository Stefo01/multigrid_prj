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
        void initialization();

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
        std::vector<std::unique_ptr<CSRMatrix>>     P_matrices;
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
            const size_t &num_coarse_nodes,
            std::map<size_t, size_t> &reversed_component_mask
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
                reversed_component_mask[i] = index;

                ++index;
            }
        }


        //void change_component_mask_level


        void build_prolongation_matrix
        (
            CSRMatrix &current_matrix,
            std::unique_ptr<CSRMatrix> &P,
            size_t &coarse_size,
            std::vector<unsigned char> &coarse_mask,
            std::map<size_t, size_t> &reversed_component_mask
        )
        {
            Matrix temporary_p_matrix(current_matrix.rows(), coarse_size);
            
            for (size_t i = 0; i < current_matrix.rows(); ++i)
            {
                if (!(coarse_mask.at(i) & 0xC0))
                {
                    temporary_p_matrix.data().at(i)[reversed_component_mask[i]] = 1.0;
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

                    //std::cout << "Index : " << index << " -> " 
                    //    << reversed_component_mask[index] << std::endl;
                    temporary_p_matrix.data().at(i)[reversed_component_mask[index]] = 
                        alpha * current_matrix.coeff(i, index) / (sum);  
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
            Matrix PtA_temp(P.cols(), P.rows());
            // Since A is symmetric cols of A are equal to cols of A

            // j from 0 to N
            #pragma omp parallel for default(none) shared(PtA_temp, current_matrix, P)
            for (size_t j = 0; j < current_matrix.rows(); ++j)
            {
                const auto &Aj_column = current_matrix.nonZerosInRow(j);
                //std::cout << "Aj nnz = " << Aj_column.size() << std::endl;
                // i from 0 to Nc
                for (size_t i = 0; i < P.cols(); ++i)
                {
                    double sum = 0.0;
                    
                    for (const auto &non_zero_entry : Aj_column)
                    {
                        sum += non_zero_entry.second * P.coeff(non_zero_entry.first, i);
                    }
                    //std::cout << "sum = " << sum << std::endl;
                    //if (1e-10 < std::abs(sum))       // Chck if it's zero
                        PtA_temp.data().at(i)[j] = sum;
                }
            }

            PtA_temp.count_non_zeros();
            //std::cout << "PtA non zeros : " << PtA_temp.non_zeros() << std::endl;

            //std::cout << "Start matrix compression" << std::endl;
            CSRMatrix PtA(PtA_temp);
            PtA.copy_from(PtA_temp);
            //std::cout << "Matrix compressed" << std::endl;

            Matrix Ac(P.cols(), P.cols());
            
            // i from 0 to Nc
            #pragma omp parallel for default(none) shared(Ac, PtA, P)
            for (size_t i = 0; i < Ac.rows(); ++i)
            {
                const auto &PtAi_row = PtA.nonZerosInRow(i);

                // j from 0 to Nc
                for (size_t j = 0; j < Ac.cols(); ++j)
                {
                    double sum = 0.0;

                    for (const auto &non_zero_entry : PtAi_row)
                    {
                        sum += non_zero_entry.second * P.coeff(non_zero_entry.first, j);
                    }
                    //if (1e-10 < std::abs(sum))       // Chck if it's zero
                        Ac.data().at(i)[j] = sum;
                }
            }
            Ac.count_non_zeros();

            //std::cout << "Start matrix compression" << std::endl;
            coarse_matrix = std::make_unique<CSRMatrix>(Ac);
            coarse_matrix->copy_from(Ac);
            //std::cout << "Matrix compressed" << std::endl;
        }




        void build_coarse_rhs
        (
            const std::vector<double> &fine_rhs,
            CSRMatrix &P,
            std::vector<double> &coarse_rhs,
            std::vector<size_t> &component_mask
        )
        {
            //coarse_rhs = std::make_unique<std::vector<double>>(P.cols(), 0);
            
            for (size_t j = 0; j < P.rows(); ++j)
            {
                const auto &Pj_column = P.nonZerosInRow(j);

                for (const auto &non_zero_entry : Pj_column)
                {
                    coarse_rhs.at(component_mask.at(non_zero_entry.first)) +=
                        non_zero_entry.second * fine_rhs.at(j);
                }
            }

        }

        double compute_residual
        (
            CSRMatrix &A,
            const std::vector<double> &x,
            const std::vector<double> &rhs,
            std::vector<double> &residual
        )
        {
            double cumulative_squared_residual = 0.0;
            for (size_t i = 0; i < A.rows(); ++i)
            {
                double Ax_i = 0.0;
                for (const auto &val : A.nonZerosInRow(i))
                {
                    Ax_i += val.second * x.at(val.first);
                }
                Ax_i = rhs.at(i) - Ax_i;
                residual.at(i) = Ax_i;
                cumulative_squared_residual += Ax_i * Ax_i;
            }
            return std::sqrt(cumulative_squared_residual);
        }

        double compute_residual_with_mask
        (
            CSRMatrix &A,
            const std::vector<double> &x,
            const std::vector<double> &rhs,
            std::vector<double> &residual
        )
        {
            double cumulative_squared_residual = 0.0;
            for (size_t i = 0; i < A.rows(); ++i)
            {
                size_t masked_i = A.component_mask.at(i);
                double Ax_i = 0.0;
                for (const auto &val : A.nonZerosInRow(i))
                {
                    Ax_i += val.second * x.at(A.component_mask.at(val.first));
                }
                Ax_i = rhs.at(masked_i) - Ax_i;
                residual.at(masked_i) = Ax_i;
                cumulative_squared_residual += Ax_i * Ax_i;
            }
            return std::sqrt(cumulative_squared_residual);
        }
        

    private:
        std::vector<unsigned char> coarse_nodes;

        
};

#endif

