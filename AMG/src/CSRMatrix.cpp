#include "../include/CSRMatrix.hpp"

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