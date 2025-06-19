#include "AMG.hpp"

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
            sol_temp.push_back(x_levels[level-1][i]);
            rhs_temp.push_back(rhs[level-1][i]);
            conter_new_mat++;
        }
    }

    std::vector<std::vector<double>> P_mat (nn, std::vector<double>(conter_new_mat)); // n x m

    
    int rows_cnt = 0;

    std::cout<< "Matrix P is " << std::endl;
    for(int i = 0; i < nn; ++i){
        if (mask_nodes_temp[i] == 0){
            P_mat[i][rows_cnt] = 1.0;
            rows_cnt++;
            for(int j = 0; j < conter_new_mat; ++j){
                std::cout<< P_mat[i][j] << " ";
            }
        } else {
            for(int j = 0; j < conter_new_mat; ++j){
                P_mat[i][j] = compute_weight_real(i,j,level - 1);
                std::cout<< P_mat[i][j] << " ";
            }
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;


    
    int i_rows_temp = 0;
    int j_cols_temp = 0;

    Matrix new_matrix(conter_new_mat, conter_new_mat);

    
    for (int i = 0; i<nn; i++)
    {
        if (mask_nodes_temp[i] == 0)
        {
            for (int j = 0; j < nn; j++){
                if (mask_nodes_temp[j] == 0){
                    if (i != j){
                        double weight_ij = compute_weight_real(i,j,level - 1);
                        new_matrix.at(i_rows_temp, j_cols_temp) = weight_ij;
                    } else {
                        new_matrix.at(i_rows_temp, j_cols_temp) = levels_matrix[level-1]->coeff(i, j);
                    }
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


double AMG::compute_weight_real(int i, int j, int level) {
    double weight = 0.0; // final result
    double sum_weak = 0.0;
    double sum_strong = 0.0;
    double a_ik, a_kj, a_km, a_kk, denominator;
    double a_sign_km,a_sign_kj,sum_a_sign_km; //these to manage sign 
    
    // start computation
    double a_ij = levels_matrix[level]->coeff(i, j);
    // 1) Compute den: a_ij + sum of weak connections
    for (const auto &neighbor : levels_matrix[level]->nonZerosInRow(i)) {
        size_t k = neighbor.first;
        if (tot_strong_connections[level][i][k] == 0 && k != i) { // Weak connections
            a_ik = levels_matrix[level]->coeff(i, k);
            sum_weak += a_ik;
        }
    }

    denominator = a_ij + sum_weak;
    //just to check den value 
    if (std::abs(denominator) < 1e-12) {
        //std::cerr << "Warning: Denominator very small at node (" << i << "," << j << ")" << std::endl;
        denominator = (denominator >= 0) ? 1e-12 : -1e-12; // mantieni segno corretto
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
        if (tot_strong_connections[level][i][k] && k != i) { // Strong connections
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
    weight = - (sum_strong) / denominator;
    return weight;
}


int AMG::apply_prolungation_operator(int level){
    int nn = x_levels[level].size(); //we start to interpolate from level n-1
    std :: cout<<nn<<std::endl;
    int temp = 0;
    int temp1 = 0;
    double weight_ij; 
    for(int i=0; i < nn; i++ ){
        //std::cout<<"for cicle "<< i << " " << x_levels[level+1][temp] << " " << mask_nodes[level][i] << std::endl;
        if(mask_nodes[level][i] == 0){
            x_levels[level][i] += x_levels[level+1][temp] ; //just to copy the value
            temp += 1;
        } else{
            //x_levels[level][i] = 0.0;
            double sum = 0.0;
            temp1 = 0;
            for (const auto &neighbor : levels_matrix[level]->nonZerosInRow(i)) {
                size_t j = neighbor.first;
                if(mask_nodes[level][j] == 0){
                    //std::cout<<"for cicle "<< j << " " << x_levels[level+1][temp1] << " " << mask_nodes[level][j] << std::endl;
                    weight_ij = compute_weight(i,j,level);
                    //std::cout<<"weight_ij: "<< weight_ij <<std::endl;
                    //std::cout<<"after weight"<<std::endl;
                    sum += weight_ij * x_levels[level+1][temp1];
                    //x_levels[level][i] += weight_ij * x_levels[level+1][temp1];
                    temp1++;
                }
            }
            x_levels[level][i] += sum; 
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
    
    Gauss_Seidel_iteration< std::vector<double> > GS(*levels_matrix[level], rhs[level]);

    for (int i = 0; i < iter_number; ++i)
    {
        GS.apply_iteration_to_vec(x_levels[level]);
    }
    
    return 0;
}

int AMG::apply_AMG(){

    // TODO : implement the AMG algorithm. This class is the main class of the algorithm
    int i;
    for (i = 0; i < number_of_levels; ++i)
    {
        std::cout << "Applying AMG on level " << i << std::endl;
        std::cout << "PRE-SMOOTHING" << std::endl;
        apply_smoother_operator(i, 10);
        std::cout << "COARSENING" << std::endl;
        int res = apply_restriction_operator(i+1);  // from 0 to 1
        //print_strong_connections(i);
        // print_x_levels(1);  
        //print_mask_nodes(i);
        if (res == 1)
            break;
    }
    std::cout << "i is "<< i <<std::endl;
    std::cout << "solution on course grid" << std::endl;
    apply_smoother_operator(i, 200);


    std::cout << "PROLUNGATION AND POST-SMOOTHING" << std::endl;
    for(i--; i >= 0; --i){
        std::cout << "PROLONGATION ON LEVEL " << i << std::endl;
        apply_prolungation_operator(i);
        std::cout << "POST-SMOOTHING level: " << i << std::endl;
        apply_smoother_operator(i, 10);
    }
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