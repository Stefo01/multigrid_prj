#include "Utilities.hpp"

const double boundary_function(const double &x, const double &y)
{
    /*
    const double k = 1;
    double r = sqrt(x * x + y * y);
    return - k * (cos(k * r) / r - k * sin(k * r));
    */
    //return x + y;
    return  std::sin(5 * std::sqrt(x * x + y * y));
    //return exp(x) * exp(-2 * y);
    //return 0.0;
}

const double forcing_term(const double &x, const double &y)
{
    //return -5 * exp(x) * exp(-2 * y);
    //return 0.0;
    return -5 * ((std::cos(5 * std::sqrt(x * x + y * y)) / std::sqrt(x * x + y * y)) - 
        (5 * std::sin(5 * std::sqrt(x * x + y * y))));
}

const double alpha(const double &/*x*/, const double &/*y*/)
{
    return 1.0;
}


int getRandomInit(int max){
    std::random_device rd;

    // Initialize a random number generator with the random device
    std::mt19937 gen(rd());
    // Create a uniform distribution within the specified range
    std::uniform_int_distribution<> distr(0, max);

    // Generate a random number within the range
    return distr(gen);
}

template <typename T>
void printVector(std::vector<T> result){
    std::cout << "Vector elements: ";
    for (int i : result) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}


