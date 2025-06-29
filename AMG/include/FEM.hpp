#ifndef FEM_HPP
#define FEM_HPP

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
#include <memory>
#include <functional>
#include "CSRMatrix.hpp"
#include "Utilities.hpp"


struct Point
{
    double x;
    double y;
    size_t set_index;
    bool is_on_boundary;
    bool is_virtual;


    Point operator+(const Point &other)
    {
        Point output;
        output.x = x + other.x;
        output.y = y + other.y;
        output.is_on_boundary = is_on_boundary & other.is_on_boundary;
        output.is_virtual = true;
        return output;
    }

    Point operator-(const Point &other)
    {
        Point output;
        output.x = x - other.x;
        output.y = y - other.y;
        output.is_on_boundary = is_on_boundary & other.is_on_boundary;
        output.is_virtual = true;
        return output;
    }

    Point operator/(const double &val)
    {
        Point output;
        output.x = x / val;
        output.y = y / val;
        output.is_virtual = true;
        return output;
    }

    Point operator*(const double &val)
    {
        Point output;
        output.x = x * val;
        output.y = y * val;
        output.is_virtual =  true;
        return output;
    }

};






class FiniteElement
{
    public:
        FiniteElement(const int &r_) : order(r_)
        {
            n_dofs_per_element = (order + 1) * (order + 2) / 2;
        }

        ~FiniteElement()
        {

        }
        const int get_order() { return order; }
        const int get_ndofs() { return n_dofs_per_element; }

        void set_vertices(std::array<Point, 3> vertices)
        {
            triangle_vertices = vertices;
        }

        virtual void set_dofs(const std::vector<Point> &dofs_) {}
        

        const double get_area()
        {
            return element_area;
        }

        const std::vector<Point> &get_dofs()
        {
            return *dofs;
        }

        virtual std::vector<Point> get_virtual_dofs() = 0;

        std::vector<double> &get_quadrature_weights() { return quadrature_weights; }


        std::vector<std::array<double, 2>> &get_gradients() { return gradients; }

        std::array<Point, 3> &get_triang()
        {
            return triangle_vertices;
        }

        bool &is_on_boundary() { return element_on_boundary; }

        std::vector<Point> &get_quadrature_points() { return quadrature_points; }

        std::function<double(Point&)> &get_basis_function(const int &local_dof_index)
        {
            return basis_functions.at(local_dof_index);
        }
        

        void print()
        {
            for (const auto &p : triangle_vertices)
                std::cout << "x = " << p.x << "; y = " << p.y << "; b = " << p.is_on_boundary << "; v = " 
                    << p.is_virtual << std::endl;
        }
    
    protected:
        const int order;
        int n_dofs_per_element;
        std::array<Point, 3> triangle_vertices;
        const std::vector<Point> *dofs;
        double element_area;
        bool element_on_boundary;
        std::vector<std::array<double, 2>> gradients;
        std::vector<Point> quadrature_points;
        std::vector<double> quadrature_weights;
        std::vector<std::function<double(Point&)>> basis_functions;

};

class LinearFE : public FiniteElement
{
    public:
        LinearFE() : FiniteElement(1)
        {
            gradients.resize(3);
            quadrature_points.resize(3);
            quadrature_weights.resize(3);
            basis_functions.resize(3);
        }
        
        ~LinearFE()
        {

        }

        std::vector<Point> get_virtual_dofs() 
        {
            return {};
        }

        void set_dofs(const std::vector<Point> &dofs_)
        {
            dofs = &dofs_;

            element_area = fabs
            (
                ((dofs->at(1).x * dofs->at(2).y) - 
                    (dofs->at(2).x * dofs->at(1).y)) +
                ((dofs->at(0).y * dofs->at(2).x) - 
                    (dofs->at(0).x * dofs->at(2).y)) +
                ((dofs->at(0).x * dofs->at(1).y) - 
                    (dofs->at(0).y * dofs->at(1).x))
            );

            // Compute the normal vectors and gradients
            element_on_boundary = false;

            for (size_t i = 0; i < 3; ++i)
            {
                size_t j, k;

                bool node_on_boundary = false;
                if (dofs->at(i).is_on_boundary)
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
                vj[0] = dofs->at(j).x - dofs->at(i).x;
                vj[1] = dofs->at(j).y - dofs->at(i).y;
                vj[2] = - 1.0;

                std::array<double, 3> vk; 
                vk[0] = dofs->at(k).x - dofs->at(i).x;
                vk[1] = dofs->at(k).y - dofs->at(i).y;
                vk[2] = - 1.0;

                // normal vector
                std::array<double, 3> n_vec;
                n_vec[0] = (vj[1] * vk[2]) - (vk[1] * vj[2]);
                n_vec[1] = (vj[2] * vk[0]) - (vk[2] * vj[0]);
                n_vec[2] = (vj[0] * vk[1]) - (vk[0] * vj[1]);

                gradients.at(i)[0] = - n_vec[0] / n_vec[2];
                gradients.at(i)[1] = - n_vec[1] / n_vec[2];

                // Set quadrature points and weigths
                quadrature_points.at(i) = dofs->at(i);
                quadrature_weights.at(i) = element_area / 3.0;

                //double D = - n_vec[0] * dofs->at(i).x - n_vec[1] * dofs->at(i).y - n_vec[2] * 1.0;
                
                
                double cx = n_vec[0] / n_vec[2];
                double cy = n_vec[1] / n_vec[2];
                double c0 = 1 - cx * dofs->at(i).x - cy * dofs->at(i).y;

                basis_functions.at(i) = 
                    [c0, cx, cy](Point &p)
                    {
                        return 2 - (c0 + cx * p.x + cy * p.y); 
                    };
                
                

            }

        }
};

class QuadraticFE : public FiniteElement
{
    public:
        QuadraticFE() : FiniteElement(2)
        {

        }

        ~QuadraticFE()
        {

        }

        std::vector<Point> get_virtual_dofs() 
        { 
            Point p01 = 
                ((triangle_vertices[1] + triangle_vertices[0]) / 2.0);
            
            Point p12 =
                ((triangle_vertices[2] + triangle_vertices[1]) / 2.0);
            
            Point p20 =
                ((triangle_vertices[0] + triangle_vertices[2]) / 2.0);
            return {p01, p20, p12};
        }
};

class ThirdOrderFE : public FiniteElement
{
    public:
        ThirdOrderFE() : FiniteElement(3)
        {

        }

        ~ThirdOrderFE()
        {

        }

        std::vector<Point> get_virtual_dofs() 
        { 
            Point p01_1 = triangle_vertices[0] +
                ((triangle_vertices[1] - triangle_vertices[0]) / 3.0);
            
            Point p01_2 = triangle_vertices[0] +
                ((triangle_vertices[1] - triangle_vertices[0]) * 2.0 / 3.0);
            
            Point p12_1 = triangle_vertices[1] + 
                ((triangle_vertices[2] - triangle_vertices[1]) / 3.0);
            
            Point p12_2 = triangle_vertices[1] + 
                ((triangle_vertices[2] - triangle_vertices[1]) * 2.0 / 3.0);
            
            Point p02_1 = triangle_vertices[0] +
                ((triangle_vertices[2] - triangle_vertices[0]) / 3.0);
            
            Point p02_2 = triangle_vertices[0] +
                ((triangle_vertices[2] - triangle_vertices[0]) * 2.0 / 3.0);
            
            Point p_internal = triangle_vertices[0] +
                ((triangle_vertices[1] - triangle_vertices[0]) * 2.0 / 3.0) +
                ((triangle_vertices[2] - triangle_vertices[1]) / 3.0);

            return {p01_1, p01_2, p02_1, p02_2, p12_1, p12_2, p_internal};
        }
};

class TriangularMesh
{
    public:
        TriangularMesh(FiniteElement &fe_) : fe(fe_)
        {}

        void import_from_msh(const std::string &mesh_file_name);
        void import_from_msh2(const std::string &mesh_file_name);
        void export_to_vtu(const std::vector<double> &sol);

        const size_t n_nodes() const { return nodes.size(); }
        const size_t n_elements() const { return num_elements; }
        const size_t n_b_nodes() const { return num_boundary_nodes; }

        const std::vector< std::array<int, 3> > &get_elements_indexes() const { return element_indexes; }
        const std::vector<Point> &get_nodes(){ return nodes; };

        class ElementIterator
        {
            private:
                int* current;
                unsigned short num_dofs;

            public:
                using iterator_category = std::forward_iterator_tag;
                using value_type        = std::vector<int>;
                using difference_type   = std::ptrdiff_t;
                using pointer           = value_type*;
                using reference         = value_type;

                ElementIterator(int* ptr, unsigned short n_dofs_) : current(ptr), num_dofs(n_dofs_)
                {

                }

                value_type operator*() const 
                {
                    return value_type(current, current + num_dofs);
                }

                ElementIterator &operator++()
                {
                    current += num_dofs;
                    return *this;
                }

                ElementIterator operator++(int)
                {
                    ElementIterator tmp = *this;
                    ++(*this);
                    return tmp;
                }

                friend bool operator==(const ElementIterator &A, const ElementIterator &B)
                {
                    return A.current == B.current;
                }

                friend bool operator!=(const ElementIterator &A, const ElementIterator &B)
                {
                    return A.current != B.current;
                }

        };

        // Range wrapper
        class ElementIterable
        {
            private:
                int* data;
                unsigned short num_dofs;
                int num_elements;
            
            public:
                ElementIterable(int* data_, unsigned short num_dofs_, int num_elements_)
                : data(data_), num_dofs(num_dofs_), num_elements(num_elements_)
                {

                }

                ElementIterator begin() const
                {
                    return ElementIterator(data, num_dofs);
                }

                ElementIterator end() const
                {
                    return ElementIterator(data + num_dofs * num_elements, num_dofs);
                }
        };

        ElementIterable element_iterators()
        {
            return ElementIterable(test_indices.get(), fe.get_ndofs(), n_elements());
        }

    private:
        std::vector<Point> nodes;
        std::vector<Point> dofs;
        std::vector< std::array<int, 3> > element_indexes;
        std::unique_ptr<int[]> test_indices;
        size_t num_boundary_nodes;
        FiniteElement &fe;
        int num_elements;
};


#endif