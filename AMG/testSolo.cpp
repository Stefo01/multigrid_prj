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




typedef struct
{
    double x;
    double y;
    bool is_on_boundary;
    size_t set_index;
} Point;



class TriangularMesh
{
    public:
        TriangularMesh()
        {}

        void import_from_msh(const std::string &mesh_file_name);
        void export_to_vtu(const std::vector<double> &sol);

        const size_t n_nodes() const { return nodes.size(); }
        const size_t n_elements() const { return element_indexes.size(); }
        const size_t n_b_nodes() const { return num_boundary_nodes; }

        const std::vector< std::array<int, 3> > &get_elements_indexes() const { return element_indexes; }
        const std::vector<Point> &get_nodes(){ return nodes; };

    private:
        std::vector<Point> nodes;
        std::vector< std::array<int, 3> > element_indexes;
        size_t num_boundary_nodes;
};

class Matrix
{
    public:
        Matrix(const size_t &rows_, const size_t &cols_):
            n_rows(rows_),
            n_cols(cols_),
            num_non_zeros(0)
        {
            non_zero_cells.resize(rows_);
        }

        double &at(const size_t &row, const size_t &col);
        
        void count_non_zeros();
        const size_t non_zeros(){ return num_non_zeros; }
        const size_t rows(){ return n_rows; }
        const size_t cols(){ return n_cols; }
        const std::vector< std::map<size_t, double> > &data()
        { return non_zero_cells; }

        void print()
        {
            for (size_t i = 0; i < n_rows; ++i)
            {
                for (size_t j = 0; j < n_cols; ++j)
                {
                    std::cout << at(i, j) << "\t\t";
                }
                std::cout <<  std::endl;
            }
            std::cout << std::endl;
        }
        

    private:
        std::vector< std::map<size_t, double> > non_zero_cells;

        const size_t n_rows;
        const size_t n_cols;

        size_t num_non_zeros;
};

class CSRMatrix
{
    public:
        CSRMatrix(const size_t &rows_, const size_t &cols_, const size_t &nnz):
        n_rows(rows_), n_cols(cols_)
        {
            row_ptrs = new std::pair<size_t, double>*[rows_];
            non_zero_elements = new std::pair<size_t, double>[nnz];
            end = &non_zero_elements[nnz];
            std::vector<bool> xx(rows_,true);
            mask_amg.push_back(xx);
        }

        CSRMatrix(Matrix&A):
        n_rows(A.rows()), n_cols(A.cols())
        {
            row_ptrs = new std::pair<size_t, double>*[A.rows()];
            non_zero_elements = new std::pair<size_t, double>[A.non_zeros()];
            end = &non_zero_elements[A.non_zeros()];
            std::vector<bool> xx(A.rows(),true);
            mask_amg.push_back(xx);
        }

        ~CSRMatrix()
        {
            delete[] row_ptrs;
            delete[] non_zero_elements;
        }



        void copy_from(Matrix &A);   
        const double coeff(const size_t &row, const size_t &col);
        const std::vector< std::pair<size_t, double> > nonZerosInRow(const size_t &row);
        const size_t rows(){ return mask_amg[mask_amg.size() - 1].size(); }
        const size_t cols(){ return mask_amg[mask_amg.size() - 1].size(); }

        void print()
        {
            for (size_t i = 0; i < n_rows; ++i)
            {
                for (size_t j = 0; j < n_cols; ++j)
                {
                    std::cout << coeff(i, j) << " ";
                }
                std::cout <<  std::endl;
            }
            std::cout << std::endl;
        }


    private:
        size_t n_rows;
        size_t n_cols;

        std::vector<std::vector<bool>> mask_amg;

        std::pair<size_t, double> **row_ptrs;
        std::pair<size_t, double> *non_zero_elements;

        std::pair<size_t, double> *end;

};

std::ostream &operator<<(std::ostream &os, const std::vector<double> &vec)
{
    os << "[ ";
    for (const auto &p : vec)
    {
        os << p << " ";
    }
    os << "]";
    return os;
}

template<class Vector>
class SmootherClass{
    public:

        virtual void apply_iteration_to_vec(std::vector<double> &sol) = 0;

        inline friend std::vector<double>& operator*(std::vector<double> &x_k, SmootherClass &B)
        {
            B.apply_iteration_to_vec(x_k);
            return x_k;
        } // x^(k+1) = x^(k) * B

};

template<class Vector>
class Gauss_Seidel_iteration : public SmootherClass<Vector>{
    private:    
        CSRMatrix &m_A;
        Vector &b; // Ax = b
    public:

        Gauss_Seidel_iteration(CSRMatrix &A, Vector &f) : m_A(A), b(f) {}
        
        void apply_iteration_to_vec(std::vector<double> &sol) override{
            for(size_t i = 0; i < m_A.rows(); i++){
                double sum = 0;
                for (const auto &element : m_A.nonZerosInRow(i))
                {
                    if (element.first != i)
                    {
                        sum += element.second * sol[element.first];
                    }
                }
                sol[i] = (this->b[i] - sum) / m_A.coeff(i,i);
            }
        }
};

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

/// @brief 
/// @tparam T type in input
/// @param R is the input CSRmatrix
/// @param elementI is the node index to be evaluate 
/// @param multi is equal to 1 when I want only the strong connection, all other are set to 0. If multi equal to 2, all elements are equal to 1, exept strong connections, to 2
/// @return a vector with position i-th equal to 1 if strongly connected, 0 otherwise

template <typename T>
void valueStrongConnection(CSRMatrix &R, const size_t elementI, int multi, std::vector<T> &Ret){
    // compute max
    T max = 0;
    std::vector<std::pair<size_t, T> > NonZR = R.nonZerosInRow(elementI);


    for(std::pair<size_t, T> el : NonZR)
    {  
        if(el.first != elementI && std::abs(el.second) > max)
            max = el.second;
    }
    
    // compute strong connections
    
    //to do:
    double epsilon = 0.25;


    for(std::pair<size_t, T> el : NonZR)
    {
        if(multi== 2){
            if(Ret[el.first]==1 && el.first != elementI && std::abs(el.second) >= epsilon*max)
            {
                Ret[el.first] = multi;
            }else if(el.first == elementI){
                Ret[el.first] = 0;
            }
        }else{
            if( el.first != elementI && std::abs(el.second) >= epsilon*max)
            {
                Ret[el.first] = 1;
            }
        }
        
    }
}

template <typename T>
void valueStrongConnectionV2(CSRMatrix &R, const size_t elementI, std::vector<int> &Ret){
    // compute max
    T max = 0;
    std::vector<std::pair<size_t, T> > NonZR = R.nonZerosInRow(elementI);

    for(std::pair<size_t, T> el : NonZR)
    {  
        if(el.first != elementI && std::abs(el.second) > max)
            max = el.second;
    }
    double epsilon = 0.25;

    for(std::pair<size_t, T> el : NonZR)
    {
        if( el.first != elementI && std::abs(el.second) >= epsilon*max)
        {
            Ret[el.first] = 1;
        }
        
    }
}

template <typename T>
T evaluateNodeV2(std::vector<std::vector<T>> allNodes, std::vector<T> V, size_t elementI ){

    if (allNodes[elementI].size() != V.size()) {
        std::cerr << "Vectors must be of the same length!" << std::endl;
        return -1;
    }

    return std::inner_product(V.begin(), V.end(), allNodes[elementI].begin(), 0);
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
std::vector<int> AMGV2(CSRMatrix &A){
    std::vector<int> R(A.rows(), 1); // numero righe = numero nodi
    size_t index = getRandomInit(R.size());
    bool GoOn;
    std::vector<std::vector<int>> TOTStrongConnections (A.rows(), std::vector<int>(A.rows(), 0));
    for(int i = 0; i< R.size(); i++){
        valueStrongConnectionV2<T>(A,i, TOTStrongConnections[i]); 
    }
    std::cout<<"Strong connection matrix setted\n"<<std::endl;
    std::vector<size_t> Indextemp;
    do{
        if(GoOn == false){
            for(size_t i= 0; i < Indextemp.size(); i++)
            {
                R[Indextemp[i]] = 0;
            }
            break;
        }
        GoOn = false;
        R[index] = 0;
        for(size_t i= 0; i < R.size(); i++)
        {
            if(R[i] == 1 && TOTStrongConnections[index][i] == 1)
            {
                R[i] = 2;
                //Indextemp.push_back(i);
                //
            }
        }

        int maxTemp = 0;
        for(size_t i= 0; i < Indextemp.size(); i++)
        {
            int rettt = evaluateNodeV2<int>(TOTStrongConnections, R, Indextemp[i]);
            if(rettt%2 == 0)
            {
                Indextemp.erase(Indextemp.begin() + i);
            }else if (rettt > maxTemp){
                maxTemp = rettt;
                index = i;
                GoOn = true;
            }
        }
    }while(true);
    return R;
    // A.size() is the row number
}


template <typename T>
std::vector<bool> AMGV3(CSRMatrix &A){
    // initialization of return vector R
    int TotNumn = A.rows();
    std::vector<int> R(TotNumn); // numero righe = numero nodi
    for (int i = 0; i < TotNumn; ++i)
        R[i] = i;

    // randomize it
    std::srand(static_cast<unsigned>(std::time(0)));
    std::random_shuffle(R.begin(), R.end());

    // build the nodes graph
    std::vector<std::vector<int>> TOTStrongConnections (TotNumn, std::vector<int>(TotNumn, 0));
    for(int i = 0; i< TotNumn; i++){
        valueStrongConnectionV2<T>(A,i, TOTStrongConnections[i]); 
    }
    std::cout<<"Strong connection graph build\n"<<std::endl;

    // compute the maximal independent set (MIS)

    std::unordered_set<int> selected;
    std::vector<bool> isSelected(TotNumn, false);

    for (int node : R) {
        // Check if any neighbor is already selected
        bool canSelect = true;
        for (int neighbor = 0; neighbor < TotNumn; ++neighbor) {
            if (TOTStrongConnections[node][neighbor] && isSelected[neighbor]) {
                canSelect = false;
                break;
            }
        }

        // If no neighbors are selected, add this node to the selected set
        if (canSelect) {
            selected.insert(node);
            isSelected[node] = true;

            // Stop if we have selected approximately half the nodes
            if (selected.size() >= TotNumn / 2) break;
        }
    }
    return isSelected;
    // A.size() is the row number
}



template <typename T>
std::vector<bool> AMGV4(CSRMatrix &A){   // non funziona u cazz
    // initialization of return vector R
    int nn = A.rows();
    std::vector<bool> R(nn, 0); // numero righe = numero nodi
    size_t index = getRandomInit(nn);
    bool GoOn = true;
    std::vector<std::vector<int>> TOTStrongConnections (nn, std::vector<int>(nn, 0));
    for(int i = 0; i< nn; i++){
        valueStrongConnectionV2<T>(A,i, TOTStrongConnections[i]); 
    }
    std::cout<<"Strong connection matrix setted\n"<<std::endl;

    std::cout<< "\n\nOur strong connection matrix:\n";

    for (int i = 0; i < nn; i++){
        for (int j = 0; j < nn; j++){
            std::cout << " " << TOTStrongConnections[i][j] << " ";
        }
        std::cout<< "\n";
    }
    std::cout << "\n\n";

    std::vector<int> all_nodes_not_yet_selected(nn);

    for(int i = 0; i< nn; i++){
        all_nodes_not_yet_selected[i] = std::accumulate(TOTStrongConnections[i].begin(), TOTStrongConnections[i].end(), 0) -1;
    }


    // init with index the random index
    while(GoOn)
    {
        GoOn = false;

        // step 1:
        all_nodes_not_yet_selected[index] = 0;
        TOTStrongConnections[index][index] = 0;
        for(size_t i=0; i< nn; ++i){
            if(TOTStrongConnections[index][i]  && (i != index) && all_nodes_not_yet_selected[i]){
                R[i] = 1;
                all_nodes_not_yet_selected[i] = 0;
                for(size_t j=0; j< nn; ++j){
                    if((TOTStrongConnections[i][j] != 0) && (all_nodes_not_yet_selected[j] != 0))
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


    int count_dim_col = 0;
    for (int i = 0; i< nn; i++){
        if(R[i] == 0)
            count_dim_col++;

    }

    std::vector<std::vector<double>> I_k (nn, std::vector<double>(count_dim_col, 0.0));

    int ccc = 0;
    for (int i = 0; i < nn; i++){
        if(R[i] == 0){
            I_k[i][ccc] = 1;
            ccc++;
        } else {
            int tempss = 0;
            for (int j = 0; j < nn; j++){
                
                if ( (TOTStrongConnections[i][j] == 1) && (R[j] == 0) ){
                    tempss++;
                    std::cout << "diooo\n";

                }
            }
            std::cout<< tempss;

            int tempss2 = 0;
            for (int j = 0; j < nn; j++) {
                if (TOTStrongConnections[i][j] == 1 ){
                    tempss2++;
                    if (R[j] == 0)
                        if(tempss == 1)
                            I_k[i][tempss2] = 0.8;
                        else
                            I_k[i][tempss2] = 1 / tempss;

                }
            }
        }
    }

    std::cout<< "\n\nI_k matrix:\n";

    for (int i = 0; i < nn; i++){
        for (int j = 0; j < count_dim_col; j++){
            std::cout << " " << I_k[i][j] << " ";
        }
        std::cout<< "\n";
    }
    std::cout << "\n\n";
    

    std::cout << "Gesù " << count_dim_col << std::endl;


    return R;
}


/// @brief this function evaluate the importance measure  of node i-th
/// @tparam T 
/// @param NodeToEval correspond to row of matrix A, is the node i-th
/// @param allNodes[i] is set to 2 if is in F and 1 if it is still undecided. 0 if in course 
/// @return 
template <typename T>
T evaluateNode(const std::vector<T> allNodes, std::vector<T> V ){

    if (allNodes.size() != V.size()) {
        std::cerr << "Vectors must be of the same length!" << std::endl;
        return -1;
    }

    return std::inner_product(V.begin(), V.end(), allNodes.begin(), 0);
}



template <typename T>
void printVector(std::vector<T> result){
    std::cout << "Vector elements: ";
    for (int i : result) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}

// matrix A is the correspondence of the connections of a graph, aij represent the edje between node i and j. So, if the connection exist, 
// it will be represented by a value different from 0, otherwise the two nodes are disconnected.
// So the "primi vicini" are represented by all values in a specific row differents from zero
// all spd matrix are all square graph and the AMG method correspond to geometric multigrid

template <typename T>
std::vector<T> AMG(CSRMatrix &A){
    std::vector<T> R(A.rows(), 1); // numero righe = numero nodi
    size_t index = getRandomInit(A.rows() - 1 );
    bool GoOn;
    std::vector<std::vector<T>> TOTStrongConnections (A.rows(), std::vector<T>(A.rows(), 0));
    for(int i = 0; i< R.size(); i++){
        valueStrongConnection<T>(A,i,1, TOTStrongConnections[i]); 
    }
    std::cout<<"Finish first part\n\n"<<std::endl;
    do{
        GoOn = false;
        std::vector<T> Tem(A.cols(), 1);
        valueStrongConnection<T>(A,index,2, R); 
        //printVector(R);

        int tempMax = 0;
        
        for(int i = 0; i< R.size(); i++)
        {
            if(R[i] == 1)
            {
                int tt = evaluateNode<T>(R, TOTStrongConnections[i]);

                if(tt > tempMax)
                {
                    GoOn = true;
                    tempMax = tt;
                    index = i;
                }
            }
        }
    }while(GoOn);
    return R;
    // A.size() is the row number
}

int main()
{    
    TriangularMesh mesh;
    mesh.import_from_msh("mesh/mesh2.msh");
    //mesh.export_to_vtu();
    std::cout << "Mesh imported! There are " << mesh.n_nodes() << " nodes and "
        << mesh.n_elements() << " elements." << std::endl;


    // Initializing the matrix
    std::cout << "Initializing the matrix" << std::endl;
    Matrix A_temp(mesh.n_nodes() - mesh.n_b_nodes(), mesh.n_nodes() - mesh.n_b_nodes());
    //Matrix B_temp(mesh.n_nodes() - mesh.n_b_nodes(), mesh.n_nodes());
    std::vector<double> rhs(mesh.n_nodes() - mesh.n_b_nodes());

    for (const auto &element : mesh.get_elements_indexes())
    {
        bool element_on_boundary = false;
        double element_area = fabs(
            ((mesh.get_nodes()[element[1]].x * mesh.get_nodes()[element[2]].y) - 
                (mesh.get_nodes()[element[2]].x * mesh.get_nodes()[element[1]].y)) +
            ((mesh.get_nodes()[element[0]].y * mesh.get_nodes()[element[2]].x) - 
                (mesh.get_nodes()[element[0]].x * mesh.get_nodes()[element[2]].y)) +
            ((mesh.get_nodes()[element[0]].x * mesh.get_nodes()[element[1]].y) - 
                (mesh.get_nodes()[element[0]].y * mesh.get_nodes()[element[1]].x))
        );

        double alpha_integral = element_area *
            (alpha(mesh.get_nodes()[element[0]].x, mesh.get_nodes()[element[0]].y) + 
            alpha(mesh.get_nodes()[element[1]].x, mesh.get_nodes()[element[1]].y) + 
            alpha(mesh.get_nodes()[element[2]].x, mesh.get_nodes()[element[2]].y)) / 3;
            
        std::array< std::array<double, 2>, 3> gradients;
        
        for (size_t i = 0; i < 3; ++i)
        {
            size_t j, k;

            bool node_on_boundary = false;
            if (mesh.get_nodes()[element[i]].is_on_boundary)
            {
                element_on_boundary = true;
                node_on_boundary = true;
            }

            switch (i)
            {
            case 0:
                j = 1;
                k = 2;
                break;
            case 1:
                j = 2;
                k = 0;
                break;     
            default:
                j = 0;
                k = 1;
                break;
            }

            std::array<double, 3> vj; 
            vj[0] = mesh.get_nodes()[element[j]].x - mesh.get_nodes()[element[i]].x;
            vj[1] = mesh.get_nodes()[element[j]].y - mesh.get_nodes()[element[i]].y;
            vj[2] = - 1.0;

            std::array<double, 3> vk; 
            vk[0] = mesh.get_nodes()[element[k]].x - mesh.get_nodes()[element[i]].x;
            vk[1] = mesh.get_nodes()[element[k]].y - mesh.get_nodes()[element[i]].y;
            vk[2] = - 1.0;

            // normal vector
            std::array<double, 3> n_vec;
            n_vec[0] = (vj[1] * vk[2]) - (vk[1] * vj[2]);
            n_vec[1] = (vj[2] * vk[0]) - (vk[2] * vj[0]);
            n_vec[2] = (vj[0] * vk[1]) - (vk[0] * vj[1]);

            gradients[i][0] = - n_vec[0] / n_vec[2];
            gradients[i][1] = - n_vec[1] / n_vec[2];

        }
        for (size_t i = 0; i < 3; ++i)
        {
            if (!mesh.get_nodes()[element[i]].is_on_boundary)
            {

                for (size_t j = 0; j < 3; ++j)
                {
                    if (!mesh.get_nodes()[element[j]].is_on_boundary)
                    {
                        // still not implemented the quadrature rule 
                        A_temp.at(mesh.get_nodes()[element[i]].set_index, mesh.get_nodes()[element[j]].set_index) += 
                            alpha_integral *                              // integral of alpha on the element
                            (gradients[i][0] * gradients[j][0] + 
                            gradients[i][1] * gradients[j][1]) / 3;
                                              
                    }
                }

                rhs.at(mesh.get_nodes()[element[i]].set_index) += 
                    forcing_term(mesh.get_nodes()[element[i]].x, mesh.get_nodes()[element[i]].y) * 
                    element_area / 3;               // volume of the corresponding tetrahedron

            }

        }
        // Dirichlet boundary condition
        if (element_on_boundary)
        {
            for (size_t i = 0; i < 3; ++i)
            {
                if (!mesh.get_nodes()[element[i]].is_on_boundary)
                {
                    // F - Bg
                    for (size_t j = 0; j < 3; ++j)
                    {
                        if(mesh.get_nodes()[element[j]].is_on_boundary)
                        {

                            rhs.at(mesh.get_nodes()[element[i]].set_index) -=
                                boundary_function(mesh.get_nodes()[element[j]].x,
                                    mesh.get_nodes()[element[j]].y) * alpha_integral *
                                (gradients[i][0] * gradients[j][0] +
                                gradients[i][1] * gradients[j][1]) / 3;
                        }
                    }
                }
            }
        }

    }
    std::cout << "Matrix created succesfully!" << std::endl;
    //A_temp.print();
    std::cout << "Counting non zero elements..." << std::endl;
    A_temp.count_non_zeros();
    //B_temp.count_non_zeros();
    std::cout << "There are " << A_temp.non_zeros() << " non zero elements." << std::endl;
    
    std::cout << "Compressing the matrix..." << std::endl;
    CSRMatrix A(A_temp);
    //CSRMatrix B(B_temp);
    A.copy_from(A_temp);
    std::vector<bool> result = AMGV4<double>(A); // Specify template type
    std::cout << "Result: ";
    for (const auto &val : result) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    //B.copy_from(B_temp);
    std::cout << "Matrix compressed successfully!"<< std::endl;

    Gauss_Seidel_iteration< std::vector<double> > GS(A, rhs);
    std::vector<double> sol(mesh.n_nodes() - mesh.n_b_nodes());

    //std::cout << sol << std::endl;
    std::cout << " " << std::endl;
    //A.print();
    //B.print();

    for (int i = 0; i < 5000; ++i)
    {
        sol * GS;
    }

    mesh.export_to_vtu(sol);

    //std::cout << sol << std::endl;

    return 0;
}


double &Matrix::at(const size_t &row, const size_t &col)
{
    if (row >= n_rows || col >= n_cols)
    {
        throw std::invalid_argument("Index out of range");
    }

    std::map<size_t, double> &row_map = non_zero_cells.at(row);
    return row_map[col];
}

void Matrix::count_non_zeros(){
    for (const auto &row : non_zero_cells)
    {
        for(const auto &val : row)
        {
            if (std::get<1>(val) != 0)
                ++num_non_zeros;
        }
    }
}


void CSRMatrix::copy_from(Matrix &A)
{
    if (A.rows() != n_rows || A.cols() != n_cols)
        throw std::invalid_argument("Input matrix doesn't match the size!");
    size_t index = 0;
    size_t row_index = 0;
    for (const auto &row : A.data())
    {
        row_ptrs[row_index] = &non_zero_elements[index];
        for (const auto &val : row)
        {
            if(std::get<1>(val) == 0)
                continue;
            
            non_zero_elements[index] = val;
            ++index;          
        }
        ++row_index;
    }
}

const double CSRMatrix::coeff(const size_t &row, const size_t &col)
{
    if (row >= n_rows || col >= n_cols)
        throw std::invalid_argument("Input out of range!!");
    

    const std::pair<size_t, double> *row_start = row_ptrs[row];
    const std::pair<size_t, double> *row_end = ((row < n_rows - 1) ? row_ptrs[row + 1] : end);

    for (auto *element = row_start; element != row_end; ++element)
    {
        if (element->first == col)
            return element->second;          
    }
    return 0.0;
    
}

const std::vector< std::pair<size_t, double> > CSRMatrix::nonZerosInRow(const size_t &row)
{
    if (row >= n_rows)
        throw std::invalid_argument("Input out of range");
    std::pair<size_t, double> *row_start = row_ptrs[row];
    std::pair<size_t, double> *row_end =
        ((row < n_rows - 1) ? row_ptrs[row + 1] : end);

    std::vector< std::pair<size_t,  double> > result(row_start, row_end);
    return result;
}


void TriangularMesh::import_from_msh(const std::string &mesh_file_name)
{
    std::ifstream mesh_file(mesh_file_name);
    std::string current_word;
    bool flag = false;
    size_t num;
    size_t min_tag;
    size_t max_tag;

    while(mesh_file >> current_word)
    {   
        // read nodes according to gmsh documentation
        if (current_word == "$Nodes")
        {
            std::string temp;

            // numEntityBlocks(size_t)
            mesh_file >> temp;

            // numNodes(size_t)
            mesh_file >> temp;
            num = stoi(temp);

            // minNodeTag(size_t)
            mesh_file >> temp;
            min_tag = stoi(temp);

            // maxNodeTag(size_t)
            mesh_file >> temp;
            max_tag = stoi(temp);
            
            size_t current_tag = min_tag;
            while (current_tag < max_tag)
            {
                // entityDim(int)
                mesh_file >> temp;

                // entityTag(int)
                mesh_file >> temp;

                // parametric(int; 0 or 1)
                mesh_file >> temp;

                // numNodesInBlock(size_t)
                mesh_file >> temp;
                size_t nodes_in_block = stoi(temp);
                current_tag += nodes_in_block;
                for (int i = 0; i < nodes_in_block; ++i)
                {
                    // nodeTag(size_t)
                    mesh_file >> temp;
                    //std::cout << "Node tag = " << temp << std::endl;
                }
                for (int i = 0; i < nodes_in_block; ++i)
                {
                    Point p;
                    // x(double)
                    mesh_file >> temp;
                    //std::cout << "\tx = " << temp << std::endl;
                    p.x = stod(temp);

                    // y(double)
                    mesh_file >> temp;
                    //std::cout << "\ty = " << temp << std::endl;
                    p.y = stod(temp);

                    // z(double)
                    mesh_file >> temp;
                    //std::cout << "\tz = " << temp << std::endl;

                    p.is_on_boundary = false;
                    nodes.push_back(p);
                }

            }

        }
        if (current_word == "$EndNodes")
            continue;
        
        if (current_word == "$Elements" )
        {
            std::string temp;

            // numEntityBlocks(size_t)
            mesh_file >> temp;

            // numElements(size_t)
            mesh_file >> temp;
            num = stoi(temp);

            // minElementTag(size_t)
            mesh_file >> temp;
            min_tag = stoi(temp);

            // maxElementTag(size_t)
            mesh_file >> temp;
            max_tag = stoi(temp);

            size_t current_tag = min_tag;
            while (current_tag < max_tag)
            {
                // entityDim(int)
                mesh_file >> temp;

                // entityTag(int)
                mesh_file >> temp;

                // elementType(int)
                mesh_file >> temp;
                int element_type = stoi(temp);

                // numElementsInBlock(size_t)
                mesh_file >> temp;
                size_t elements_in_block = stoi(temp);

                current_tag += elements_in_block;
                for (int i = 0; i < elements_in_block; ++i)
                {
                    std::array<int, 3> triangle;
                    // elementTag(size_t)
                    mesh_file >> temp;
                    //std::cout << "Element tag = " << temp << std::endl;
                    
                    if (1 == element_type)
                    {
                        // These are boundary elements
                        mesh_file >> temp;
                        nodes.at(stoi(temp) - 1).is_on_boundary = true;

                        mesh_file >> temp;
                        nodes.at(stoi(temp) - 1).is_on_boundary = true;
                        continue;
                    }
                    else if(2 == element_type)
                    {
                        // nodeTag(size_t)
                        mesh_file >> temp;
                        triangle[0] = stoi(temp) - 1;
                        
                        // nodeTag(size_t)
                        mesh_file >> temp;
                        //std::cout << temp << ", ";
                        triangle[1] = stoi(temp) - 1;

                        // nodeTag(size_t)
                        mesh_file >> temp;
                        //std::cout << temp << std::endl;
                        triangle[2] = stoi(temp) - 1;
                        
                        
                        element_indexes.push_back(triangle);
                    }
                    else
                        throw std::invalid_argument("This element type is not implemented");

                    
                }
            }
        }

        if (current_word == "$EndElements")
            continue;
        
    }
    mesh_file.close();

    size_t counter_on_boundary = 0;
    size_t counter_not_on_boundary = 0;
    num_boundary_nodes = 0;
    for (auto &node : nodes)
    {
        if (node.is_on_boundary)
        {
            node.set_index = counter_on_boundary;
            ++counter_on_boundary;
            ++num_boundary_nodes;
        }else
        {
            node.set_index = counter_not_on_boundary;
            ++counter_not_on_boundary;
        }
    }
}


void TriangularMesh::export_to_vtu(const std::vector<double> &sol)
// Output file
{
    std::ofstream outfile("output.vtu");

    outfile << "<?xml version=\"1.0\"?>" << std::endl;
    outfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << 
        std::endl;
    outfile << "\t<UnstructuredGrid>" << std::endl;
    outfile << "\t\t<Piece NumberOfPoints=\"" << nodes.size() << "\" NumberOfCells=\"" << 
        element_indexes.size() << "\">" << std::endl;

    outfile << "\t\t\t<Points>" << std::endl;
    outfile << "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << 
        std::endl;
    outfile << "\t\t\t\t\t";
    // Insert points here
    for (const auto &point : nodes)
    {
        outfile << point.x << " " << point.y << " " << 0 << std::endl;
    }
    outfile << std::endl;
    outfile << "\t\t\t\t</DataArray>" << std::endl;
    outfile << "\t\t\t</Points>" << std::endl;
    

    outfile << "\t\t\t<Cells>" << std::endl;
    outfile << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    // Insert cells here
    
    for (const auto &triang : element_indexes)
    {
        outfile << triang[0] << " " << triang[1] << " " << triang[2] << std::endl;
    }
    outfile << std::endl;
    
    outfile << "\t\t\t\t</DataArray>" << std::endl;
    outfile << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    // Insert offsets here
    
    for (int i = 0; i < element_indexes.size(); ++i)
    {
        outfile << 3 * (i + 1) << std::endl;
    }
    outfile << std::endl;
    
    outfile << "\t\t\t\t</DataArray>" << std::endl;
    outfile << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
    // Insert types here
    
    for (int i = 0; i < element_indexes.size(); ++i)
    {
        outfile << 5 << std::endl;
    }
    outfile << std::endl;
    
    outfile << "\t\t\t\t</DataArray>" << std::endl;
    outfile << "\t\t\t</Cells>" << std::endl;


    outfile << "\t\t\t<PointData>" << std::endl;
    outfile << "<DataArray type=\"Float32\" Name=\"u\" NumberOfComponents=\"1\" format=\"ascii\">" <<
        std::endl;
    // Insert values here
    for (const auto &node : nodes)
    {   
        double value;
        if (node.is_on_boundary)
            value = boundary_function(node.x, node.y);
        else
            value = sol.at(node.set_index);
        outfile << value << std::endl;
    }
    outfile << std::endl;

    outfile << "</DataArray>" << std::endl;
    outfile << "\t\t\t</PointData>" << std::endl;


    outfile << "\t\t\t<CellData>" << std::endl;

    outfile << "\t\t\t</CellData>" << std::endl;


    outfile << "\t\t</Piece>" << std::endl;
    outfile << "\t</UnstructuredGrid>" << std::endl;
    outfile << "</VTKFile>" << std::endl;


    outfile << "" << std::endl;



    outfile.close();
}