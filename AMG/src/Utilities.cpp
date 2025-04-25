#include "Utilities.hpp"

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