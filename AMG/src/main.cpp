#include "CSRMatrix.hpp"
#include "Utilities.hpp"
#include "AMG.hpp"

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


int main()
{    
    Utilities::TriangularMesh mesh;
    mesh.import_from_msh("../mesh/mesh2.msh");
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
            (Utilities::alpha(mesh.get_nodes()[element[0]].x, mesh.get_nodes()[element[0]].y) + 
            Utilities::alpha(mesh.get_nodes()[element[1]].x, mesh.get_nodes()[element[1]].y) + 
            Utilities::alpha(mesh.get_nodes()[element[2]].x, mesh.get_nodes()[element[2]].y)) / 3;
            
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
                    Utilities::forcing_term(mesh.get_nodes()[element[i]].x, mesh.get_nodes()[element[i]].y) * 
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
                                Utilities::boundary_function(mesh.get_nodes()[element[j]].x,
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
    //std::vector<bool> result = AMGV4<double>(A); // Specify template type
    // std::cout << "Result: ";
    // for (const auto &val : result) {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;

    //B.copy_from(B_temp);
    std::cout << "Matrix compressed successfully!"<< std::endl;

    Utilities::Gauss_Seidel_iteration< std::vector<double> > GS(A, rhs);
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