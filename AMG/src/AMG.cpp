#include "AMG.hpp"

void AMG::value_strong_connections(const size_t elementI, std::vector<bool> &Ret, int level){

    double max = 0.0;
    std::vector<std::pair<size_t, double> > NonZR = levels_matrix[level].nonZerosInRow(elementI);

    for(std::pair<size_t, double> el : NonZR)
    {  
        if(el.first != elementI && std::abs(el.second) > max)
            max = el.second;
    }

    for(std::pair<size_t, double> el : NonZR)
    {
        if( el.first != elementI && std::abs(el.second) >= EPSILON*max)
        {
            Ret[el.first] = 1;
        }
        
    }
}

double AMG::evaluate_node(std::vector<std::vector<double>> allNodes, std::vector<double> V, size_t elementI ){

    if (allNodes[elementI].size() != V.size()) {
        std::cerr << "Vectors must be of the same length!" << std::endl;
        return -1;
    }

    return std::inner_product(V.begin(), V.end(), allNodes[elementI].begin(), 0);
}

int AMG::apply_restriction_operator(int level){
    if (levels_matrix.size() <= level) {
        std::cerr << "Level " << level << " does not exist!" << std::endl;
        return -1;
    } else if (level == 0){
        std::cerr << "Level 0 is the coarsest level, restriction operator cannot be applied!" << std::endl;
        return -1;
    }

    int nn = levels_matrix[level].rows();
    
    for(int i = 0; i< nn; i++){
        value_strong_connections(i, tot_strong_connections[level][i], level); 
    }
    std::cout<<"Strong connection matrix setted for level "<< level <<std::endl;

    std::vector<int> all_nodes_not_yet_selected(nn);

    for(int i = 0; i< nn; i++){
        all_nodes_not_yet_selected[i] = std::accumulate(tot_strong_connections[level][i].begin(), tot_strong_connections[level][i].end(), 0) -1;
    }

    bool GoOn = true;
    size_t index = Utilities::getRandomInit(nn);
    while(GoOn)
    {
        GoOn = false;

        // step 1:
        all_nodes_not_yet_selected[index] = 0;
        tot_strong_connections[level][index][index] = 0;
        for(size_t i=0; i< nn; ++i){
            if(tot_strong_connections[level][index][i]  && (i != index) && all_nodes_not_yet_selected[i]){
                mask_nodes[level][i] = 1;
                all_nodes_not_yet_selected[i] = 0;
                for(size_t j=0; j< nn; ++j){
                    if((tot_strong_connections[level][i][j] != 0) && (all_nodes_not_yet_selected[j] != 0))
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
            std::cout << "Vector is empty!" << std::endl;
            break;
        }
    }

    // step 2:
    int conter_new_mat = 0;
    std::vector<double> sol_temp;

    for(int i = 0; i < mask_nodes[level].size(); i++){
        if (mask_nodes[level][i] == 0){
            sol_temp.push_back(solution[level-1][i]);
            conter_new_mat++;
        }
    }

    Matrix new_matrix(conter_new_mat, conter_new_mat);
    

    int i_rows_temp = 0;
    int j_cols_temp = 0;
    for (int i = 0; i<nn; i++)
    {
        if (mask_nodes[level][i] == 0)
        {
            for (int j = 0; j < nn; j++){
                if (mask_nodes[level][j] == 0){
                    new_matrix.at(i_rows_temp, j_cols_temp) = levels_matrix[level].coeff(i, j);
                    j_cols_temp++;
                }
            }
            i_rows_temp++;
        }
    }

    conter_new_mat = 0;
    for(int i= 0; i < general_mask_for_rhs.size(); i++){
        if (general_mask_for_rhs[i] == 1){
            if (mask_nodes[level][conter_new_mat] == 1) // fine node
                general_mask_for_rhs[i] = 0;
            conter_new_mat++;
        }
    }

    new_matrix.count_non_zeros();

    CSRMatrix A(new_matrix);
    A.copy_from(new_matrix);
    levels_matrix.push_back(A);
    solution.push_back(sol_temp);

    return 0;
}

int AMG::apply_prolungation_operator(int level){

    return 0;
}

int AMG::apply_smoother_operator(int level, int iter_number){

    // Possubility: choose between different solvers
    std::vector<double> rhs_temp;
    for (int i= 0; i < general_mask_for_rhs.size(); i++){
        if (general_mask_for_rhs[i] == 1){
            rhs_temp.push_back(rhs[i]);
        }
    }
    Utilities::Gauss_Seidel_iteration< std::vector<double> > GS(levels_matrix[level], rhs);

    for (int i = 0; i < iter_number; ++i)
    {
        solution[level] * GS;
    }

    return 0;
}

int AMG::apply_AMG(){

    // TODO : implement the AMG algorithm. This class is the main class of the algorithm

    return 0;
}