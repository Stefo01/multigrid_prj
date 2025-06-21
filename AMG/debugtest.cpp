#include <iostream>
#include "CSRMatrix.hpp"
#include "Utilities.hpp"
#include "AMG.hpp"
#include "FEM.hpp"



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
    //CSRMatrix A(A_temp);
    //CSRMatrix B(B_temp);
    //A.copy_from(A_temp);

    std::cout << "Matrix compressed successfully!"<< std::endl;

    std::vector<double> sol(mesh.n_nodes() - mesh.n_b_nodes());  

    //Gauss_Seidel_iteration< std::vector<double> > GS(A, rhs);

    std::cout << " " << std::endl;

    //for (int i = 0; i < 5000; ++i)
    //{
    //    sol * GS;
    //}




    std::vector<std::unique_ptr<CSRMatrix>> system_matrices;
    std::vector<std::unique_ptr<CSRMatrix>> prolongation_matrices;
    int amg_levels = 4;
    std::unique_ptr<CSRMatrix> tmp = std::make_unique<CSRMatrix>(A_temp);
    tmp->copy_from(A_temp);

    system_matrices.push_back(std::move(tmp));
    RestrictionOperator R;

    for (int level = 1; level < amg_levels; ++level)
    {
        std::vector<unsigned char> coarse_mask(system_matrices.at(level - 1)->rows(), 0);
        size_t num_coarse = R.select_coarse_nodes(*system_matrices.at(level - 1), coarse_mask);

        std::cout << "There are " << num_coarse << " coarse nodes at level " 
            << level << std::endl;
        
        std::vector<size_t> component_mask;
        std::map<size_t, size_t> reversed_component_mask;
        R.build_component_mask(coarse_mask, component_mask, num_coarse, reversed_component_mask);

        std::unique_ptr<CSRMatrix> P;
        R.build_prolongation_matrix
        (
            *system_matrices.at(level - 1), P, 
            num_coarse, coarse_mask, 
            reversed_component_mask
        );

        std::cout << "P size : " << (P->rows()) << " x " << (P->cols()) << std::endl;

        std::unique_ptr<CSRMatrix> Ac;
        R.build_coarse_matrix(*system_matrices.at(level - 1), *P, Ac);

        prolongation_matrices.push_back(std::move(P));
        system_matrices.push_back(std::move(Ac));
    }

    mesh.export_to_vtu(sol);


    //for (size_t i = 0; i < Ac->rows(); ++i)
    //{
    //    std::cout << "[ ";
    //    for (const auto &val : Ac->nonZerosInRow(i))
    //    {
    //        std::cout << val.first << " : " << val.second << " ";
    //    }
    //    std::cout << "]" << std::endl;
    //}


    return 0;
}


