#include "AMG.hpp"

#include <iostream> // For std::cout, std::cerr, std::endl
#include <vector>   // For std::vector
#include <string>   // For std::string (used in printMatrix for names)
#include <iomanip>  // For std::fixed, std::setprecision
#include <unordered_map>


// Function to print a matrix (for std::vector<std::vector<double>>)
// This function takes a matrix and a name string for printing purposes.
void printMatrix(const std::vector<std::vector<double>>& matrix, const std::string& name) {
    std::cout << name << " (" << matrix.size() << "x";
    if (!matrix.empty()) {
        std::cout << matrix[0].size();
    } else {
        std::cout << "0";
    }
    std::cout << "):" << std::endl;

    if (matrix.empty()) {
        std::cout << "Matrix is empty." << std::endl;
        return;
    }

    // Set precision for double output
    std::cout.precision(4); // Adjust precision as needed (e.g., 2 for two decimal places)
    std::cout << std::fixed; // Use fixed-point notation for doubles

    for (const auto& row : matrix) {
        for (double val : row) {
            std::cout << val << "\t";
        }
        std::cout << std::endl;
    }
}



bool AMG::is_connected(int i, int j, int level) {
    auto neighbors = levels_matrix[level]->nonZerosInRow(i);
    for (const auto &p : neighbors) {
        if (p.first == j)
            return true;
    }
    return false;
}

// Like input, it gets the level that you want to restrict the matrix. So, the future one
int AMG::apply_restriction_operator(int level){
    if (P_matrices.size() < level) {
        std::cerr << "Level " << level << " does not exist!" << std::endl;
        return -1;
    } else if (level == 0){
        std::cerr << "Level 0 is the coarsest level, restriction operator cannot be applied!" << std::endl;
        return -1;
    }

    int nn = x_levels[level - 1].size();

    std::vector<double> new_sol (P_matrices.at(level - 1)->cols());

    for(int i = 0; i < nn; i++){
        for(const auto &m_neighbor : P_matrices.at(level - 1)->nonZerosInRow(i))
        {
            size_t m = m_neighbor.first;
            new_sol[m] += P_matrices.at(level - 1)->coeff(i, m) * x_levels[level - 1][i];
        }
    }

    x_levels.push_back(new_sol);

    return 0;
}

void AMG::initialization(){

    RestrictionOperator R;

    for (int level = 1; level < number_of_levels; level++)
    {
        std::vector<unsigned char> coarse_mask(levels_matrix.at(level - 1)->rows(), 0);
        size_t num_coarse = R.select_coarse_nodes(*levels_matrix.at(level - 1), coarse_mask);

        std::cout << "There are " << num_coarse << " coarse nodes at level " 
            << level << std::endl;
        
        std::vector<size_t> component_mask;
        std::map<size_t, size_t> reversed_component_mask;
        R.build_component_mask(coarse_mask, component_mask, num_coarse, reversed_component_mask);

        std::unique_ptr<CSRMatrix> P;
        R.build_prolongation_matrix
        (
            *levels_matrix.at(level - 1), P, 
            num_coarse, coarse_mask, 
            reversed_component_mask
        );

        std::vector<double> rhs_new (P->cols());

        for(int i = 0; i < P->rows(); i++){
            for(const auto &m_neighbor : P->nonZerosInRow(i))
            {
                size_t m = m_neighbor.first;
                rhs_new[m] += P->coeff(i, m) * rhs[level - 1][i];
            }
        }
        rhs.push_back(rhs_new);

        std::cout << "P size : " << (P->rows()) << " x " << (P->cols()) << std::endl;

        std::unique_ptr<CSRMatrix> Ac;
        R.build_coarse_matrix(*levels_matrix.at(level - 1), *P, Ac);

        P_matrices.push_back(std::move(P));
        levels_matrix.push_back(std::move(Ac));
        
    }
}


double AMG::compute_weight(int i, int j, int level) {
    // Interpolazione classica: peso = a_ij / somma(|a_ik|) su tutti i coarse k collegati a i
    double a_ij = levels_matrix[level]->coeff(i, j);
    double sum = 0.0;

    for (const auto& neighbor : levels_matrix[level]->nonZerosInRow(i)) {
        size_t k = neighbor.first;
        if (mask_nodes[level][k] == 0 && k != i) { // solo coarse nodes
            sum += std::abs(levels_matrix[level]->coeff(i, k));
        }
    }
    //std::cout << "a_ij: " << a_ij << ", sum: " << sum << std::endl;

    if (std::abs(sum) < 1e-12) return 0.0; // evita divisione per zero
    return -a_ij / sum;
}


// int sign_func(double i){
//     if (i == 0){
//         return 0;
//     } else if (i < 0){
//         return -1;
//     } else{
//         return 1;
//     }
// }

// double AMG::compute_weight_real(int i, int j, int level) {
//     double weight = 0.0; // final result
//     double sum_weak = 0.0;
//     double sum_strong = 0.0;
//     double a_ik, a_kj, a_km, a_kk, denominator;
//     double a_sign_km,a_sign_kj,sum_a_sign_km; //these to manage sign 
//     //std::cout << "Computing weight for nodes (" << i << ", " << j << ") at level " << level << std::endl;
//     // start computation
//     double a_ij = levels_matrix[level]->coeff(i, j);
//     // 1) Compute den: a_ij + sum of weak connections
//     for (const auto &neighbor : levels_matrix[level]->nonZerosInRow(i)) {
        
//         size_t k = neighbor.first;
//         if (tot_strong_connections[level][i][k] == 0 && k != i && mask_nodes[level][k]) { // Weak connections
//             a_ik = levels_matrix[level]->coeff(i, k);
//             sum_weak += a_ik;
//             //std::cout << "Weak connection: a_ik = " << a_ik << " at node (" << i << "," << k << ")" << "sum of weak"<<sum_weak << std::endl;
//         }
//     }
//     denominator = a_ij + sum_weak;
//     //just to check den value 
//     if (std::abs(denominator) < 1e-12) {
//         //std::cerr << "Warning: Denominator very small at node (" << i << "," << j << ")" << std::endl;
//         //::cout << "Warning: Denominator very small at node (" << i << "," << j << "), setting to a small value." << std::endl;
//         denominator = (denominator <= 0) ? 1e-12 : -1e-12; // mantieni segno corretto
//     }

//     // 2) a_ij:
//     sum_strong += a_ij;

//       // 3) strong connections:
//     //for (const auto &neighbor : levels_matrix[level]->nonZerosInRow(i))
//     //{
//     //    size_t = k = neighbor.first;
//     //    double a_ik = neighbor.second;
//     //}
//     for (const auto &neighbor : levels_matrix[level]->nonZerosInRow(i)) {
//         size_t k = neighbor.first;
//         a_kk = sign_func(levels_matrix[level]->coeff(k, k));
        
//         if (tot_strong_connections[level][i][k] && k != i && mask_nodes[level][k]) { // Strong connections
//             a_ik = levels_matrix[level]->coeff(i, k);
//             a_kj = levels_matrix[level]->coeff(k, j);
//             a_kk = levels_matrix[level]->coeff(k, k);

//             // Somma a_sign_km on all course nodes
//             sum_a_sign_km = 0.0;
//             for (const auto &m_neighbor : levels_matrix[level]->nonZerosInRow(k)) {
//                 size_t m = m_neighbor.first;
//                 if (mask_nodes[level][m] == 0) { // m is course node
//                     sum_a_sign_km += levels_matrix[level]->coeff(k, m);//a_sign_km;
//                 }
//             }

//             if (std::abs(sum_a_sign_km) > 1e-12) {
//                 //sum_strong += a_ik * (a_sign_kj / sum_a_sign_km);
//                 sum_strong += a_ik * (a_kj / sum_a_sign_km);
//             }
//         }
//     }

//     // 4) Computing final weight
//     //std::cout << "sum_strong: " << sum_strong << ", denominator: " << denominator << std::endl;
//     weight =  -(sum_strong) / denominator;
//     return weight;
// }              

int AMG::apply_prolungation_operator(int level){
    int nn = x_levels[level].size(); //we start to interpolate from level n-1

    for(int i = 0; i < nn; i++){
        for(const auto &m_neighbor : P_matrices[level]->nonZerosInRow(i))
        {
            size_t m = m_neighbor.first;
            x_levels[level][i] += P_matrices[level]->coeff(i, m) * x_levels[level + 1][m];
        }
    }

    x_levels.pop_back();

    return 0;
}



int AMG::apply_smoother_operator(int level, int iter_number){
    if (level < 0 || level >= levels_matrix.size()) {
        std::cerr << "Invalid level: " << level << std::endl;
        return -1;
    }
    if (iter_number <= 0) {
        std::cerr << "Invalid number of iterations: " << iter_number << std::endl;
        return -1;
    }
    
    Gauss_Seidel_iteration< std::vector<double> > GS(*levels_matrix[level], rhs[level]);

    for (int i = 0; i < iter_number; ++i)
    {
        GS.apply_iteration_to_vec(x_levels[level]);
    }
    
    return 0;
}

double AMG::compute_residual(int level){
    int nn = x_levels[level].size();
    std::vector<double> residual(nn, 0.0);
    for (size_t i = 0; i < nn; ++i) {
        double Ax_i = 0.0;
        for (const auto &neighbor : levels_matrix.at(level)->nonZerosInRow(i)) {
            Ax_i += neighbor.second * x_levels[level][neighbor.first];
        }
        residual[i] = rhs[level][i] - Ax_i;
    }

    // Optional: Print or return residual norm
    double norm = 0.0;
    for (double r_i : residual) {
        norm += r_i * r_i;
    }
    norm = std::sqrt(norm);
    std::cout << "Residual norm: " << norm << std::endl;
    return norm;
}

int AMG::apply_AMG(){

    // TODO : implement the AMG algorithm. This class is the main class of the algorithm
    int i;
    initialization();
    std::cout << "Initialization done"<< std::endl;
    for (i = 0; i < number_of_levels - 1; ++i)
    {
        std::cout << "Applying AMG on level " << i << std::endl;
        std::cout << "PRE-SMOOTHING" << std::endl;
        apply_smoother_operator(i, 10);
        std::cout << "COARSENING" << std::endl;
        apply_restriction_operator(i+1);  // from 0 to 1
        //print_strong_connections(i);
        //print_x_levels(1); 
        //print_mask_nodes(i);
    }
    std::cout << "solution on course grid" << std::endl;
    apply_smoother_operator(i, 200);

    std::cout << "PROLUNGATION AND POST-SMOOTHING" << std::endl;
    for(i--; i >= 0; --i){
        std::cout << "PROLONGATION ON LEVEL " << i << std::endl;
        apply_prolungation_operator(i);
        std::cout << "POST-SMOOTHING level: " << i << std::endl;
        apply_smoother_operator(i, 10);
    }
    compute_residual(0);

    std::cout << "AMG applied successfully!" << std::endl;
    return 0;
}

void AMG::print_strong_connections(int level){
    if (level < 0 || level > tot_strong_connections.size()) {
        std::cerr << "Invalid level: " << level << std::endl;
        return;
    }
    for (int i = 0; i < tot_strong_connections[level].size(); ++i)
    {
        std::cout << "Node " << i << ": ";
        for (int j = 0; j < tot_strong_connections[level][i].size(); ++j)
        {
            std::cout << tot_strong_connections[level][i][j] << " ";
        }
        std::cout << std::endl;
    }
}



void AMG::print_CSRmatrix(int level){
    if (level < 0 || level >= levels_matrix.size()) {
        std::cerr << "Invalid level: " << level << std::endl;
        return;
    }
    levels_matrix[level]->print();
}


void AMG::print_x_levels(int level){
    if (level < 0 || level >= x_levels.size()) {
        std::cerr << "Invalid level: " << level << std::endl;
        return;
    }

    std::cout << "Solution at level " << level << ": ";
    for (int i = 0; i < x_levels[level].size(); ++i)
    {
        std::cout << x_levels[level][i] << " ";
    }
    std::cout << std::endl;
}

void AMG::print_mask_nodes(int level) {
    if (level < 0 || level >= mask_nodes.size()) {
        std::cerr << "Invalid level: " << level << std::endl;
        return;
    }
    std::cout << "Coarse nodes (C=coarse, F=fine) at level " << level << ": ";
    for (size_t i = 0; i < mask_nodes[level].size(); ++i) {
        if (mask_nodes[level][i] == 0)
            std::cout << "C ";
        else
            std::cout << "F ";
    }
    std::cout << std::endl;
}
