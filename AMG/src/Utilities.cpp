#include "Utilities.hpp"

void TriangularMesh::import_from_msh(const std::string &mesh_file_name)
{
    std::ifstream mesh_file(mesh_file_name);
    std::string current_word;
    bool flag = false;
    size_t num;
    size_t min_tag;
    size_t max_tag;

    std::map<int, std::vector<std::pair<int, int>>> visited_pairs;  // this is for creating additional dofs

    unsigned short virtual_dofs_per_edge = fe.get_order() - 1;
    unsigned short virtual_dofs_internal = fe.get_ndofs() - 3 * (virtual_dofs_per_edge + 1);

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

            nodes.resize(num);          // We allocate the number of nodes needed

            // minNodeTag(size_t)
            mesh_file >> temp;
            min_tag = stoi(temp);

            // maxNodeTag(size_t)
            mesh_file >> temp;
            max_tag = stoi(temp);
            
            size_t current_tag = min_tag;
            size_t counter = 0;
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
                    p.is_virtual = false;
                    nodes.at(counter + i) = p;  // We do this in order to add the virtual nodes later
                    //nodes.push_back(p);
                }
                counter += nodes_in_block;
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


                if (2 == element_type)
                {
                    test_indices = std::make_unique<int[]>(elements_in_block * fe.get_ndofs());
                    num_elements = elements_in_block;
                }

                current_tag += elements_in_block;
                for (int i = 0; i < elements_in_block; ++i)
                {
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
                        int *triang = &test_indices.get()[i * fe.get_ndofs()];
                        int vertex1, vertex2, vertex3;
                        // nodeTag(size_t)
                        mesh_file >> temp;
                        vertex1 = stoi(temp) - 1;    

                        // nodeTag(size_t)
                        mesh_file >> temp;
                        //std::cout << temp << ", ";
                        vertex2 = stoi(temp) - 1;  

                        // nodeTag(size_t)
                        mesh_file >> temp;
                        //std::cout << temp << std::endl;
                        vertex3 = stoi(temp) - 1;  
                        

                        // Order the vertices in topological order
                        if (vertex1 > vertex3)
                            std::swap(vertex1, vertex3);
                        if (vertex1 > vertex2)
                            std::swap(vertex1, vertex2);
                        if (vertex2 > vertex3)
                            std::swap(vertex2, vertex3);

                        // Store the vertex indics in the structure
                        triang[0] = vertex1;
                        triang[1] = vertex2;
                        triang[2] = vertex3;

                        // add additional nodes if fem order > 1
                        if (1 < fe.get_order())
                        {
                            fe.set_vertices(
                                {nodes.at(vertex1), 
                                nodes.at(vertex2),
                                nodes.at(vertex3)}
                            );

                            std::vector<Point> virtual_dofs = fe.get_virtual_dofs();

                            // Check if the dofs were already created
                            unsigned short current_idx_local = 0;
                            for (int j = 0; j < 2; ++j)
                            {                             
                                for (int k = j + 1; k < 3; ++k)
                                {
                                    int index1 = triang[j];
                                    int index2 = triang[k];
                                    
                                    
                                    if (visited_pairs.contains(index1))
                                    {
                                        std::vector<std::pair<int, int>> &tmp = visited_pairs[index1];
                                        
                                        auto it = std::find_if
                                        (
                                            tmp.begin(), tmp.end(), 
                                            [index2](const std::pair<int, int> &p)
                                            {
                                                return p.first == index2;
                                            }
                                        );

                                        if (it != tmp.end())
                                        {
                                            // Pair already visited
                                            int current_idx_global = it->second;
                                            
                                            for (int l = 0; l < virtual_dofs_per_edge; ++l)
                                            {
                                                triang[current_idx_local + 3] = current_idx_global + l;
                                                current_idx_local++;
                                            }
                                        }
                                        else
                                        {
                                            int current_idx_global = nodes.size();
                                            tmp.push_back({index2, current_idx_global});

                                            for (int l = 0; l < virtual_dofs_per_edge; ++l)
                                            {
                                                nodes.push_back(virtual_dofs.at(current_idx_local));
                                                triang[current_idx_local + 3] = current_idx_global + l;
                                                current_idx_local++;
                                            }
                                            
                                        }
                                    }
                                    else
                                    {
                                        int current_idx_global = nodes.size();
                                        visited_pairs[index1] = {{index2, current_idx_global}};

                                        for (int l = 0; l < virtual_dofs_per_edge; ++l)
                                        {
                                            nodes.push_back(virtual_dofs.at(current_idx_local));
                                            triang[current_idx_local + 3] = current_idx_global + l;
                                            current_idx_local++;
                                        }

                                    }

                                }
                                
                            }

                            // Add the rest of the dofs that are not in common to other elements
                            for (int j = 3 * (virtual_dofs_per_edge); j < 3 * (virtual_dofs_per_edge) + virtual_dofs_internal; ++j)
                            {                              
                                int current_idx_global = nodes.size();
                                nodes.push_back(virtual_dofs.at(j));
                                triang[j + 3] = current_idx_global;
                            }

                        }

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
    std::vector<bool> visited(false, nodes.size());

    // Optimize here the dof distribution
    bool tst = true;
    for (size_t i = 0; i < element_indexes.size(); ++i)
    {
        int *triang = &test_indices.get()[i * fe.get_ndofs()];
        tst &= (element_indexes.at(i)[0] == triang[0]);
        tst &= (element_indexes.at(i)[1] == triang[1]);
        tst &= (element_indexes.at(i)[2] == triang[2]);
    }
    std::cout << "Test " << (tst ? "PASSED" : "FAILED") << std::endl;
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
        num_elements << "\">" << std::endl;

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
    
    for (const auto &triang : element_iterators())
    {
        outfile << triang[0] << " " << triang[1] << " " << triang[2] << std::endl;
    }
    outfile << std::endl;
    
    outfile << "\t\t\t\t</DataArray>" << std::endl;
    outfile << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    // Insert offsets here
    
    for (int i = 0; i < num_elements; ++i)
    {
        outfile << 3 * (i + 1) << std::endl;
    }
    outfile << std::endl;
    
    outfile << "\t\t\t\t</DataArray>" << std::endl;
    outfile << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
    // Insert types here
    
    for (int i = 0; i < num_elements; ++i)
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


