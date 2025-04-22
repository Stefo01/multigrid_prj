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
                fine_course_nodes[level][i] = 1;
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

    // TODO: create the coarse matrix
    // levels_matrix.push_back(coarse_matrix);

    return 0;
}

int AMG::apply_prolungation_operator(int level){

    return 0;
}

int AMG::apply_smoother_operator(int level){

    return 0;
}

int AMG::apply_AMG(){

    // TODO : implement the AMG algorithm. This class is the main class of the algorithm

    return 0;
}