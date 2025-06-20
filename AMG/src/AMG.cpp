#include "AMG.hpp"

#include <iostream> // For std::cout, std::cerr, std::endl
#include <vector>   // For std::vector
#include <string>   // For std::string (used in printMatrix for names)
#include <iomanip>  // For std::fixed, std::setprecision
#include <unordered_map>

std::vector<double> mat_vec_multiply(const std::vector<std::vector<double>>& A,
                                     const std::vector<double>& x,
                                     std::vector<double>& Res) {
    size_t m = A.size();
    size_t n = x.size();
    std::vector<double> b(m, 0.0);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            b[i] += A[i][j] * x[j];
            Res[i] += A[i][j] * x[j];
        }
    }

    return b;
}

std::vector<double> mat_vec_multiply_real(const std::vector<std::vector<double>>& A,
                                     const std::vector<double>& x) {
    size_t m = A.size();
    size_t n = x.size();
    std::vector<double> b(m, 0.0);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            b[i] += A[i][j] * x[j];
        }
    }

    return b;
}

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

// Function to transpose a matrix
// Takes a constant reference to the original matrix and returns a new transposed matrix.
std::vector<std::vector<double>> transposeMatrix(const std::vector<std::vector<double>>& matrix) {
    if (matrix.empty()) {
        return {}; // Return an empty matrix if the input is empty
    }

    int rows = matrix.size();       // Number of rows in the original matrix
    int cols = matrix[0].size();    // Number of columns in the original matrix

    // Create a new matrix with swapped dimensions (cols x rows) for the transpose
    std::vector<std::vector<double>> transposed_matrix(cols, std::vector<double>(rows));

    // Populate the transposed matrix
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            transposed_matrix[j][i] = matrix[i][j];
        }
    }
    return transposed_matrix;
}

// Function to perform matrix multiplication for matrices of doubles
// C = A * B
// A: R1 x C1
// B: R2 x C2
// C: R1 x C2 (where C1 must be equal to R2 for multiplication to be valid)
// Returns the resulting product matrix.
std::vector<std::vector<double>> multiplyMatrices(const std::vector<std::vector<double>>& A,
                                                  const std::vector<std::vector<double>>& B) {
    int R1 = A.size();
    if (R1 == 0) {
        std::cerr << "Error in multiplyMatrices: First matrix (A) has 0 rows." << std::endl;
        return {}; // Return empty for invalid input
    }
    int C1 = A[0].size();
    if (C1 == 0) {
        std::cerr << "Error in multiplyMatrices: First matrix (A) has 0 columns." << std::endl;
        return {}; // Return empty for invalid input
    }

    int R2 = B.size();
    if (R2 == 0) {
        std::cerr << "Error in multiplyMatrices: Second matrix (B) has 0 rows." << std::endl;
        return {}; // Return empty for invalid input
    }
    int C2 = B[0].size();
    if (C2 == 0) {
        std::cerr << "Error in multiplyMatrices: Second matrix (B) has 0 columns." << std::endl;
        return {}; // Return empty for invalid input
    }

    // Check if multiplication is possible (number of columns in A must equal number of rows in B)
    if (C1 != R2) {
        std::cerr << "Error in multiplyMatrices: The number of columns in the first matrix (" << C1
                  << ") must be equal to the number of rows in the second matrix (" << R2
                  << ") for multiplication." << std::endl;
        return {}; // Return an empty matrix to indicate an error
    }

    // Initialize result matrix C with dimensions R1 x C2, filled with zeros
    std::vector<std::vector<double>> C(R1, std::vector<double>(C2, 0.0));

    // Perform the matrix multiplication using the standard algorithm
    // C[i][j] = sum(A[i][k] * B[k][j]) for k from 0 to C1-1 (or R2-1)
    for (int i = 0; i < R1; ++i) { // Iterate over rows of matrix A (and result C)
        for (int j = 0; j < C2; ++j) { // Iterate over columns of matrix B (and result C)
            for (int k = 0; k < C1; ++k) { // Iterate over columns of A / rows of B for the sum
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C; // Return the computed product matrix
}

bool AMG::value_strong_connections(const size_t elementI, std::vector<bool> &Ret, int level, int &tot_strong_connections){

    double max = 0.0;
    bool valid = false;
    std::vector<std::pair<size_t, double> > NonZR = levels_matrix[level]->nonZerosInRow(elementI);

    for(std::pair<size_t, double> el : NonZR)
    {  
        if(el.first != elementI && std::abs(el.second) > max)
            max = std::abs(el.second);
    }

    for(std::pair<size_t, double> el : NonZR)
    {
        if( el.first != elementI && std::abs(el.second) >= EPSILON*max)
        {
            
            Ret[el.first] = 1;
            tot_strong_connections+=1;
            if (!valid)
                valid = true;
        }
    }
    return valid;
}

double AMG::evaluate_node(std::vector<std::vector<double>> allNodes, std::vector<double> V, size_t elementI ){

    if (allNodes[elementI].size() != V.size()) {
        std::cerr << "Vectors must be of the same length!" << std::endl;
        return -1;
    }

    return std::inner_product(V.begin(), V.end(), allNodes[elementI].begin(), 0);
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
    if (levels_matrix.size() < level) {
        std::cerr << "Level " << level << " does not exist!" << std::endl;
        return -1;
    } else if (level == 0){
        std::cerr << "Level 0 is the coarsest level, restriction operator cannot be applied!" << std::endl;
        return -1;
    }

    int nn = levels_matrix[level-1]->rows();

    std::vector<std::vector<bool>> strong_connections_temp(nn, std::vector<bool>(nn, 0));
    std::vector<int> all_nodes_not_yet_selected(nn);

    bool continue_AMG = false;
    bool temp_AMG;
    for(int i = 0; i< nn; i++){
        temp_AMG = value_strong_connections(i, strong_connections_temp[i], level-1, all_nodes_not_yet_selected[i]); 
        if(!continue_AMG)
            continue_AMG = temp_AMG;
    }

    if (!continue_AMG)
        return 1;

    tot_strong_connections.push_back(strong_connections_temp);
    std::cout<<"Strong connection matrix setted for level "<< level <<std::endl;


    std::vector<bool> mask_nodes_temp(nn, false);
    bool GoOn = true;
    size_t index = getRandomInit(nn);

    while(GoOn)
    {
        GoOn = false;
        //std::cout<< "Node selected "  << index << "\n";
        // step 1:
        all_nodes_not_yet_selected[index] = 0;
        for(size_t i=0; i< nn; ++i){
            if(strong_connections_temp[index][i]  && (i != index) && all_nodes_not_yet_selected[i]){
                mask_nodes_temp[i] = true;
                all_nodes_not_yet_selected[i] = 0;
                for(size_t j=0; j< nn; ++j){
                    if((strong_connections_temp[i][j] != 0) && (all_nodes_not_yet_selected[j] != 0))
                    {
                        all_nodes_not_yet_selected[j] += 2;  
                    }
                }
            }
        }

        // auto iii = std::max_element(all_nodes_not_yet_selected.begin(), all_nodes_not_yet_selected.end());
        int maxx = -1;
        for(size_t j=0; j< nn; ++j){
            if(all_nodes_not_yet_selected[j] > maxx){
                maxx = all_nodes_not_yet_selected[j];
                index = j;
            }
        }


        if (all_nodes_not_yet_selected[index]) {
            GoOn = true;
        } else {
            break;
        }
    }

    mask_nodes.push_back(mask_nodes_temp);
    //std::cout<<"Mask node setted for level "<< level <<std::endl;

    // print_mask_nodes(level-1);


    
    // step 2:
    int conter_new_mat = 0;
    std::vector<double> sol_temp;
    std::vector<double> rhs_temp;

    for(int i = 0; i < mask_nodes_temp.size(); i++){
        if (mask_nodes_temp[i] == 0){
            // sol_temp.push_back(x_levels[level-1][i]);
            // rhs_temp.push_back(rhs[level-1][i]);
            conter_new_mat++;
        }
    }

    std::vector<std::vector<double>> P_mat (nn, std::vector<double>(conter_new_mat)); // n x m


    std::vector<int> coarse_nodes;
    std::unordered_map<int, int> coarse_index_map; // fine_index → coarse_col
    int coarse_cnt = 0;
    for (int i = 0; i < nn; ++i) {
        if (mask_nodes_temp[i] == 0) {
            coarse_nodes.push_back(i);
            coarse_index_map[i] = coarse_cnt++;
        }
    }

    for (int i = 0; i < nn; ++i) {
        if (mask_nodes_temp[i] == 0) {
            int col = coarse_index_map[i];
            P_mat[i][col] = 1.0;
        } else {
            for (int cj : coarse_nodes) {
                if (is_connected(i, cj, level - 1)) {
                    int col = coarse_index_map[cj];
                    P_mat[i][col] = compute_weight_real(i, cj, level - 1);
                }
            }
        }
    }
    std::cout<< "Matrix setted for level "<< level << std::endl;
    P_matrices.push_back(P_mat);

    // std::cout<< "Matrix P is " << std::endl;  
    // for(int i = 0; i < nn; i++){
       
    //     // for(int j = 0; j < conter_new_mat; j++){
    //     //     P_mat[i][j] = compute_weight_real(i,j,level - 1);
    //     //     std::cout<< P_mat[i][j] << " ";
    //     // }
    //     // // if (mask_nodes_temp[i] == 0){
    //     // //     P_mat[i][rows_cnt] = 1.0;
    //     // //     rows_cnt++;
    //     // //     for(int j = 0; j < conter_new_mat; ++j){
    //     // //         std::cout<< P_mat[i][j] << " ";
    //     // //     }
    //     // // } else {
    //         for(int j = 0; j < conter_new_mat; ++j){
    //             //P_mat[i][j] = compute_weight_real(i,j,level - 1);
    //             std::cout<< P_mat[i][j] << " ";
    //         }
    //     std::cout<<std::endl;
    // }
    // std::cout<<std::endl;

    

    std::vector<std::vector<double>> current_levels_matrix_vec (nn, std::vector<double>(nn));

    // std::cout<<std::endl;
    // std::cout<< "Matrix was "<<std::endl;
    for (int i = 0; i < nn; i++){
        for (int j = 0; j < nn; j++){
            current_levels_matrix_vec[i][j] = levels_matrix[level-1]->coeff(i, j);
            //std::cout<< levels_matrix[level-1]->coeff(i, j) << " ";
        }
        //std::cout<<std::endl;
    }
    //std::cout<<std::endl;


    std::vector<std::vector<double>> P_mat_T = transposeMatrix(P_mat);
    std::cout<< "Trasposta done" << std::endl;
    std::vector<std::vector<double>> result_matrix = multiplyMatrices(P_mat_T, current_levels_matrix_vec);
    std::cout<< "Matrix mult done" << std::endl;
    if (!result_matrix.empty()) {
        //printMatrix(result_matrix, "Result of (P_mat_T * levels_matrix[level-1])");
    } else {
        std::cout << "Multiplication resulted in an empty matrix (likely an error during calculation)." << std::endl;
    }
    std::vector<std::vector<double>> final_result_matrix = multiplyMatrices(result_matrix, P_mat);
    if (!final_result_matrix.empty()) {
        //printMatrix(final_result_matrix, "Result of ( res_before* P_mat)");

    } else {
        std::cout << "Multiplication resulted in an empty matrix (likely an error during calculation)." << std::endl;
    }




    
    int i_rows_temp = 0;
    int j_cols_temp = 0;

    Matrix new_matrix(conter_new_mat, conter_new_mat);

    
    for (int i = 0; i<nn; i++)
    {
        if (mask_nodes_temp[i] == 0)
        {
            for (int j = 0; j < nn; j++){
                if (mask_nodes_temp[j] == 0){
                    // if (i != j){
                    //     double weight_ij = compute_weight_real(i,j,level - 1);
                    //     new_matrix.at(i_rows_temp, j_cols_temp) = weight_ij;
                    // } else {
                        new_matrix.at(i_rows_temp, j_cols_temp) = final_result_matrix[i_rows_temp][j_cols_temp];
                    //}
                    j_cols_temp++;
                }
            }
            j_cols_temp = 0;
            i_rows_temp++;
        }
    }

    new_matrix.count_non_zeros();

    std::unique_ptr<CSRMatrix> A = std::make_unique<CSRMatrix> (new_matrix);
    A->copy_from(new_matrix);
    levels_matrix.push_back(std::move(A));


    sol_temp = mat_vec_multiply_real(P_mat_T, x_levels[level - 1]);
    rhs_temp = mat_vec_multiply_real(P_mat_T, rhs[level - 1]);

    x_levels.push_back(sol_temp);
    rhs.push_back(rhs_temp);

    return 0;
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


int sign_func(double i){
    if (i == 0){
        return 0;
    } else if (i < 0){
        return -1;
    } else{
        return 1;
    }
}

double AMG::compute_weight_real(int i, int j, int level) {
    double weight = 0.0; // final result
    double sum_weak = 0.0;
    double sum_strong = 0.0;
    double a_ik, a_kj, a_km, a_kk, denominator;
    double a_sign_km,a_sign_kj,sum_a_sign_km; //these to manage sign 
    //std::cout << "Computing weight for nodes (" << i << ", " << j << ") at level " << level << std::endl;
    // start computation
    double a_ij = levels_matrix[level]->coeff(i, j);
    // 1) Compute den: a_ij + sum of weak connections
    for (const auto &neighbor : levels_matrix[level]->nonZerosInRow(i)) {
        
        size_t k = neighbor.first;
        if (tot_strong_connections[level][i][k] == 0 && k != i && mask_nodes[level][k]) { // Weak connections
            a_ik = levels_matrix[level]->coeff(i, k);
            sum_weak += a_ik;
            //std::cout << "Weak connection: a_ik = " << a_ik << " at node (" << i << "," << k << ")" << "sum of weak"<<sum_weak << std::endl;
        }
    }
    denominator = a_ij + sum_weak;
    //just to check den value 
    if (std::abs(denominator) < 1e-12) {
        //std::cerr << "Warning: Denominator very small at node (" << i << "," << j << ")" << std::endl;
        //::cout << "Warning: Denominator very small at node (" << i << "," << j << "), setting to a small value." << std::endl;
        denominator = (denominator <= 0) ? 1e-12 : -1e-12; // mantieni segno corretto
    }

    // 2) a_ij:
    sum_strong += a_ij;

      // 3) strong connections:
    //for (const auto &neighbor : levels_matrix[level]->nonZerosInRow(i))
    //{
    //    size_t = k = neighbor.first;
    //    double a_ik = neighbor.second;
    //}
    for (const auto &neighbor : levels_matrix[level]->nonZerosInRow(i)) {
        size_t k = neighbor.first;
        a_kk = sign_func(levels_matrix[level]->coeff(k, k));
        
        if (tot_strong_connections[level][i][k] && k != i && mask_nodes[level][k]) { // Strong connections
            a_ik = levels_matrix[level]->coeff(i, k);
            a_kj = levels_matrix[level]->coeff(k, j);
            a_kk = levels_matrix[level]->coeff(k, k);

            // Somma a_sign_km on all course nodes
            sum_a_sign_km = 0.0;
            for (const auto &m_neighbor : levels_matrix[level]->nonZerosInRow(k)) {
                size_t m = m_neighbor.first;
                if (mask_nodes[level][m] == 0) { // m is course node
                    sum_a_sign_km += levels_matrix[level]->coeff(k, m);//a_sign_km;
                }
            }

            if (std::abs(sum_a_sign_km) > 1e-12) {
                //sum_strong += a_ik * (a_sign_kj / sum_a_sign_km);
                sum_strong += a_ik * (a_kj / sum_a_sign_km);
            }
        }
    }

    // 4) Computing final weight
    //std::cout << "sum_strong: " << sum_strong << ", denominator: " << denominator << std::endl;
    weight =  -(sum_strong) / denominator;
    return weight;
}              

int AMG::apply_prolungation_operator(int level){
    int nn = x_levels[level].size(); //we start to interpolate from level n-1
    std :: cout<<nn<<std::endl;
    int temp = 0;
    int temp1 = 0;
    double weight_ij; 

    mat_vec_multiply(P_matrices[level], x_levels[level + 1], x_levels[level]);

    // for(int i=0; i < nn; i++ ){
    //     //std::cout<<"for cicle "<< i << " " << x_levels[level+1][temp] << " " << mask_nodes[level][i] << std::endl;
    //     if(mask_nodes[level][i] == 0){
    //         x_levels[level][i] += x_levels[level+1][temp] ; //just to copy the value
    //         temp += 1;
    //     } else{
    //         //x_levels[level][i] = 0.0;
    //         double sum = 0.0;
    //         temp1 = 0;
    //         for (const auto &neighbor : levels_matrix[level]->nonZerosInRow(i)) {
    //             size_t j = neighbor.first;
    //             if(mask_nodes[level][j] == 0){
    //                 //std::cout<<"for cicle "<< j << " " << x_levels[level+1][temp1] << " " << mask_nodes[level][j] << std::endl;
    //                 weight_ij = compute_weight_real(i,j,level);
    //                 //std::cout<<"weight_ij: "<< weight_ij <<std::endl;
    //                 //std::cout<<"after weight"<<std::endl;
    //                 sum += weight_ij * x_levels[level+1][temp1];
    //                 //x_levels[level][i] += weight_ij * x_levels[level+1][temp1];
    //                 temp1++;
    //             }
    //         }
    //         x_levels[level][i] += sum; 
    //     }
    // }

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
        for (const auto &neighbor : levels_matrix[level]->nonZerosInRow(i)) {
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
    for (int aaa = 0; aaa < 3; aaa++){
        for (i = 0; i < number_of_levels; ++i)
    {
        std::cout << "Applying AMG on level " << i << std::endl;
        std::cout << "PRE-SMOOTHING" << std::endl;
        apply_smoother_operator(i, 10);
        std::cout << "COARSENING" << std::endl;
        int res = apply_restriction_operator(i+1);  // from 0 to 1
        //print_strong_connections(i);
        //print_x_levels(1); 
        print_mask_nodes(i);
        if (res == 1)
            break;
    }
    std::cout << "i is "<< i <<std::endl;
    std::cout << "solution on course grid" << std::endl;
    apply_smoother_operator(i, 50);


    std::cout << "PROLUNGATION AND POST-SMOOTHING" << std::endl;
    for(i--; i >= 0; --i){
        std::cout << "PROLONGATION ON LEVEL " << i << std::endl;
        apply_prolungation_operator(i);
        std::cout << "POST-SMOOTHING level: " << i << std::endl;
        apply_smoother_operator(i, 10);
    }
    compute_residual(0);
    }

    //apply_smoother_operator(0, 60);
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

/*
void AMG::print_mask_nodes(int level){
    if (level < 0 || level >= mask_nodes.size()) {
        std::cerr << "Invalid level: " << level << std::endl;
        return;
    }

    std::cout << "Mask_node at Level " << level << ": ";
    for (int i = 0; i < mask_nodes[level].size(); ++i)
    {
        std::cout << mask_nodes[level][i] << " ";
        
    }
    std::cout << std::endl;
}
    */

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