#include <iostream>
#include "CSRMatrix.hpp"
#include "Utilities.hpp"




int main(int argc, char **argv)
{
    LinearFE fe;
    //QuadraticFE fe;
    //ThirdOrderFE fe;

    TriangularMesh mesh(fe);
    mesh.import_from_msh("../mesh/mesh2.msh");
    //mesh.export_to_vtu();
    std::cout << "Mesh imported! There are " << mesh.n_nodes() << " dofs and "
        << mesh.n_elements() << " elements." << std::endl;


    // We need to define the finite element space
    // could be a good idea to use a sort of overlay on top of the mesh
    // the number of elements doesn't change
    // the number of dofs changes
    //size_t n_dofs = mesh.n_nodes() + 

    // Initializing the matrix


    std::cout << "Initializing the matrix" << std::endl;
    Matrix A_temp(mesh.n_nodes() - mesh.n_b_nodes(), mesh.n_nodes() - mesh.n_b_nodes());
    //Matrix B_temp(mesh.n_nodes() - mesh.n_b_nodes(), mesh.n_nodes());
    std::vector<double> rhs(mesh.n_nodes() - mesh.n_b_nodes());



    
    for (const auto &element : mesh.element_iterators())
    {
        std::vector<Point> local_dofs(fe.get_ndofs());
        for (int i = 0; i < local_dofs.size(); ++i)
        {
            local_dofs.at(i) = mesh.get_nodes()[element.at(i)];
        }

        fe.set_dofs(local_dofs);

        bool &element_on_boundary = fe.is_on_boundary();
        double element_area = fe.get_area();

             
        std::vector<Point> &quadrature_points = fe.get_quadrature_points();
        std::vector<double> &quadrature_weights = fe.get_quadrature_weights();

        std::vector<std::array<double, 2>> &gradients = fe.get_gradients();

        for (size_t i = 0; i < local_dofs.size(); ++i)
        {
            if (!local_dofs.at(i).is_on_boundary)
            {
                //std::cout << "[ ";
                for (size_t j = 0; j < local_dofs.size(); ++j)
                {
                    //std::cout << fe.get_basis_function(i)(local_dofs.at(j)) << " ";
                    if (!local_dofs.at(j).is_on_boundary)
                    {
                        for (size_t q = 0; q < quadrature_points.size(); ++q)
                        {
                            A_temp.at(local_dofs.at(i).set_index, local_dofs.at(j).set_index) += 
                                alpha(quadrature_points.at(q).x, quadrature_points.at(q).y) *                              // integral of alpha on the element
                                (gradients.at(i)[0] * gradients.at(j)[0] + 
                                gradients.at(i)[1] * gradients.at(j)[1]) * quadrature_weights.at(q);
                        }
                                              
                    }
                    

                }
                //std::cout << "]" << std::endl;

                for (size_t q = 0; q < quadrature_points.size(); ++q)
                {
                    rhs.at(local_dofs.at(i).set_index) += 
                        forcing_term(local_dofs.at(i).x, local_dofs.at(i).y) *
                        fe.get_basis_function(i)(quadrature_points.at(q)) * 
                        quadrature_weights.at(q);
                }
                

            }

        }
        // Dirichlet boundary condition
        if (element_on_boundary)
        {
            for (size_t i = 0; i < local_dofs.size(); ++i)
            {
                if (!local_dofs.at(i).is_on_boundary)
                {
                    // F - Bg
                    for (size_t j = 0; j < local_dofs.size(); ++j)
                    {
                        if(local_dofs.at(j).is_on_boundary)
                        {
                            for (size_t q = 0; q < quadrature_points.size(); ++q)
                            {
                                rhs.at(local_dofs.at(i).set_index) -=
                                    boundary_function(local_dofs.at(j).x,
                                        local_dofs.at(j).y) * 
                                    alpha(quadrature_points.at(q).x, quadrature_points.at(q).y) *
                                    (gradients[i][0] * gradients[j][0] +
                                    gradients[i][1] * gradients[j][1]) * quadrature_weights.at(q);
                            }
                            
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

    std::cout << "Matrix compressed successfully!"<< std::endl;

    std::vector<double> sol(mesh.n_nodes() - mesh.n_b_nodes());  

    Gauss_Seidel_iteration< std::vector<double> > GS(A, rhs);

    //std::cout << sol << std::endl;
    std::cout << " " << std::endl;
    //A.print();
    //B.print();

    for (int i = 0; i < 5000; ++i)
    {
        sol * GS;
    }

    mesh.export_to_vtu(sol);


    return 0;
}


