#include "AMG.hpp"

void AMG::value_strong_connections(const size_t elementI, std::vector<int> &Ret, int level){

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

    return 0;
}

int AMG::apply_prolungation_operator(int level){

    return 0;
}

int AMG::apply_smoother_operator(int level){

    return 0;
}

int AMG::apply_AMG(){


    return 0;
}