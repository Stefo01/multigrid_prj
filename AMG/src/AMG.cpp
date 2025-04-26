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

// Like input, it gets the level that you want to restrict the matrix. So, the future one
int AMG::apply_restriction_operator(int level){
    if (levels_matrix.size() < level) {
        std::cerr << "Level " << level << " does not exist!" << std::endl;
        return -1;
    } else if (level == 0){
        std::cerr << "Level 0 is the coarsest level, restriction operator cannot be applied!" << std::endl;
        return -1;
    }

    int nn = levels_matrix[level-1].rows();

    std::vector<std::vector<bool>> strong_connections_temp(nn, std::vector<bool>(nn, 0));
    
    for(int i = 0; i< nn; i++){
        value_strong_connections(i, strong_connections_temp[i], level-1); 
    }
    tot_strong_connections.push_back(strong_connections_temp);
    std::cout<<"Strong connection matrix setted for level "<< level <<std::endl;

    std::vector<int> all_nodes_not_yet_selected(nn);

    for(int i = 0; i< nn; i++){
        for (int j = 0; j< nn; j++){
            if (strong_connections_temp[i][j]){
                all_nodes_not_yet_selected[i] += 1;
            }
        }
    }

    std::vector<bool> mask_nodes_temp(nn, false);

    bool GoOn = true;
    size_t index = getRandomInit(nn);

    while(GoOn)
    {
        GoOn = false;

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
            std::cout << "Vector is empty!" << std::endl;
            break;
        }
    }

    mask_nodes.push_back(mask_nodes_temp);
    std::cout<<"Mask node setted for level "<< level <<std::endl;

    // print_mask_nodes(level-1);

    
    // step 2:
    int conter_new_mat = 0;
    std::vector<double> sol_temp;
    std::vector<double> rhs_temp;

    for(int i = 0; i < mask_nodes_temp.size(); i++){
        if (mask_nodes_temp[i] == 0){
            sol_temp.push_back(x_levels[level-1][i]);
            rhs_temp.push_back(rhs[level-1][i]);
            conter_new_mat++;
        }
    }

    Matrix new_matrix(conter_new_mat, conter_new_mat);

    int i_rows_temp = 0;
    int j_cols_temp = 0;
    for (int i = 0; i<nn; i++)
    {
        if (mask_nodes_temp[i] == 0)
        {
            for (int j = 0; j < nn; j++){
                if (mask_nodes_temp[j] == 0){
                    new_matrix.at(i_rows_temp, j_cols_temp) = levels_matrix[level-1].coeff(i, j);
                    j_cols_temp++;
                }
            }
            j_cols_temp = 0;
            i_rows_temp++;
        }
    }

    new_matrix.count_non_zeros();

    CSRMatrix A(new_matrix);
    A.copy_from(new_matrix);
    levels_matrix.push_back(A);
    x_levels.push_back(sol_temp);
    rhs.push_back(rhs_temp);

    return 0;
}

int AMG::compute_weight(int i, int j, int level){
    double a_ij = levels_matrix[level].coeff(i,j);
    double a_ii = levels_matrix[level].coeff(i,i);
    double weight = 0;
    double a_ik = 0;
    double a_kj = 0;
    //course node strong connected are useles for the interpolation
    
    //fine node strong connected with i
    int ns = tot_strong_connections[level][i].size();

    // Sum for strong connections
    for (int k = 0; k < ns; k++) {
        if (tot_strong_connections[level][i][k] == 1) {
            // Strong connection found, accumulate contribution to the weight
            a_ik = levels_matrix[level].coeff(i, k);
            a_kj = levels_matrix[level].coeff(k, j);
            weight += (a_ik * a_kj) / a_ii;
        }
    }
    
    // Sum for weak connections
    for (int k = 0; k < ns; k++) {
        if (tot_strong_connections[level][i][k] == 0) {
            // Weak connection found, accumulate contribution to the weight
            a_ik = levels_matrix[level].coeff(i, k);
            a_kj = levels_matrix[level].coeff(k, j);
            weight += (a_ik * a_kj) / (a_ii + a_kj);  // Adjusting the weak connections contribution
        }
    }
    return weight;
}

int AMG::apply_prolungation_operator(int level){
    int nn = x_levels[level+1].size(); //we start to interpolate from level n-1
    double weight_ij; 
    for(int i=0; i < nn; i++ ){
        if(mask_nodes[level][i] == 0){
            x_levels[level][i] += x_levels[level+1][i] ; //just to copy the value
        } else{
            //x_levels[level][i] = 0.0;
            for (int j = 0; j < nn; j++) { 
                weight_ij = compute_weight(i,j,level);
                x_levels[level][i] += weight_ij * x_levels[level+1][j];
            }
        }
    }

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

    Gauss_Seidel_iteration< std::vector<double> > GS(levels_matrix[level], rhs[level]);

    for (int i = 0; i < iter_number; ++i)
    {
        GS.apply_iteration_to_vec(x_levels[level]);
    }
    
    return 0;
}

int AMG::apply_AMG(){

    // TODO : implement the AMG algorithm. This class is the main class of the algorithm
    apply_restriction_operator(1);  // from 0 to 1
    std::cout << "Here we go 1!" << std::endl;

    // print_strong_connections(0);
    // print_x_levels(1);  
    // print_mask_nodes(0);

    apply_smoother_operator(1, 10); // Does not work yet
    std::cout << "Here we go 2!" << std::endl;
    apply_prolungation_operator(1);
    std::cout << "Here we go 3!" << std::endl;
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
    levels_matrix[level].print();
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