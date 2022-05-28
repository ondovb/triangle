#include <iostream>
#include <string.h>
#include <string>
#include <fstream>
#include <ctime>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ofstream;
using std::to_string;

int ** tri;
int ** exp;
int rows;

void fill(int row, int col)
{
    if ( row >= 0 && row < rows && col >= 0 && col < row * 2 + 1 )
    {
        exp[row][col] = 1;
    }
}

void expandNaive(int r)
{
    for ( int i = 0; i < rows; i++ )
    {
        int cols = i * 2 + 1;
        
        for ( int j = 0; j < cols; j++ )
        {
            if ( tri[i][j] == 1 )
            {
                for ( int k = 0; k < r + 1; k++ )
                {
                    for ( int l = 0; l < 2 * (r * 2 - k) + 1; l++ )
                    {
                        if ( j % 2 )
                        {
                            fill(i + k, j + l - 2 * (r - k));
                            
                            if ( k < r )
                            {
                                fill(i - k - 1, j + l - 2 * (r) - 1);
                            }
                        }
                        else
                        {
                            if ( k < r )
                            {
                                fill(i + k + 1, j + l - 2 * (r - k) + 1);
                            }
                            
                            fill(i - k, j + l - 2 * (r));
                        }
                    }
                }
            }
        }
    }
}

int abs(int a, int b) {return a > b ? a - b : b - a;}

class Node
{
    public:
    Node() : color(0), children(0) {}
    
    int color;
    
    Node ** children;

    void buildUp(int ** img, int size, int row = 0, int col = 0, int depth = 0)
    {
        for ( int i = 0; i < depth; i++ )
        {
            //cout << '+';
        }
    
        //cout << "buildUp\t" << size << '\t' << row << '\t' << col << endl;
        color = img[row][col];
    
        for ( int i = 0; i < size; i++ )
        {
            for ( int j = 0; j < 2 * (i + 1) - 1; j++ )
            {
                if ( img[i + row][j + col] != color )
                {
                    // split
                
                    children = new Node * [4];
                
                    for ( int k = 0; k < 4; k++ )
                    {
                        children[k] = new Node();
                    }
                
                    children[0]->buildUp(img, size / 2, row, col, depth + 1);
                    children[1]->buildUp(img, size / 2, row + size / 2, col, depth + 1);
                    children[2]->buildUp(img, size / 2, row + size / 2, col + size, depth + 1);
                
                    children[3]->buildDown(img, size / 2, row + size - 1, col + size - 1, depth + 1);
                
                    return;
                }
            }
        }
    }

    void buildDown(int ** img, int size, int row = 0, int col = 0, int depth = 0)
    {
        for ( int i = 0; i < depth; i++ )
        {
            //cout << '+';
        }
    
        //cout << "buildDown\t" << size << '\t' << row << '\t' << col << endl;
        color = img[row][col];
    
        for ( int i = 0; i < size; i++ )
        {
            for ( int j = 0; j < 2 * (i + 1) - 1; j++ )
            {
                if ( img[row - i][col - i * 2 + j] != color )
                {
                    // split
                
                    children = new Node * [4];
                
                    for ( int k = 0; k < 4; k++ )
                    {
                        children[k] = new Node();
                    }
                
                    children[0]->buildDown(img, size / 2, row, col, depth + 1);
                    children[1]->buildDown(img, size / 2, row - size / 2, col - size, depth + 1);
                    children[2]->buildDown(img, size / 2, row - size / 2, col, depth + 1);
                
                    children[3]->buildUp(img, size / 2, row - size + 1, col - size + 1, depth + 1);
                
                    return;
                }
            }
        }
    }

    void expandUp(int radius, int size, int row = 0, int col = 0, int depth = 0)
    {
        if ( (size <= radius + 1 && children != 0) || (children == 0 && color == 1) )
        {
            // all black
            
            for ( int i = 0; i < size; i++ )
            {
                for ( int j = 0; j < 2 * (i + 1) - 1; j++ )
                {
                    fill(i + row, j + col);
                }
            }
            
            int radiusSW = 0;
            int radiusN = 0;
            int radiusSE = 0;
            int radiusS = 0;
            int radiusNE = 0;
            int radiusNW = 0;
            
            int dists[size];
            
            if ( children == 0 )
            {
                radiusSW = radius;
                radiusN = radius;
                radiusSE = radius;
                radiusS = radius;
                radiusNE = radius;
                radiusNW = radius;
            }
            else
            {
                radiusN = radius - distNup(size, 0, size);
                radiusSW = radius - distSWup(size, 0, size);
                radiusSE = radius - distSEup(size, 0, size);
                radiusS = radius - distNdown(size, 0, size);
                radiusNE = radius - distSEdown(size, 0, size);
                radiusNW = radius - distSWdown(size, 0, size);
            }
            
            // SSE
            //
            if ( children )
            {
                memset(&dists, 0, size * sizeof(int));
                frontierSSUp(dists, size, size, 1);
            }
            else
            {
                for ( int i = 0; i < size; i++ )
                {
                    dists[i] = size;
                }
            }
            expandFrontierDiagBase(dists, size, radius, row, col, 1, 1);
            
            // SSW
            //
            if ( children )
            {
                memset(&dists, 0, size * sizeof(int));
                frontierSSUp(dists, size, size, -1);
            }
            expandFrontierDiagBase(dists, size, radius, row, col + size * 2, 1, -1);
            
            // NNW
            //
            if ( children )
            {
                memset(&dists, 0, size * sizeof(int));
                frontierSSDown(dists, size, size, -1);
            }
            expandFrontierDiagTip(dists, size, radius, row, col, -1, -1);
            
            // NNE
            //
            if ( children )
            {
                memset(&dists, 0, size * sizeof(int));
                frontierSSDown(dists, size, size, 1);
            }
            expandFrontierDiagTip(dists, size, radius, row, col, -1, 1);
            
            // W
            //
            if ( children )
            {
                memset(&dists, 0, size * sizeof(int));
                frontierHorzUp(dists, size, size, -1);
            }
            else
            {
                for ( int i = 0; i < size; i++ )
                {
                    dists[i] = size * 2 - i;
                }
            }
            expandFrontierHorzWUp(dists, size, radius, row, col);
            
            // E
            //
            if ( children )
            {
                memset(&dists, 0, size * sizeof(int));
                frontierHorzUp(dists, size, size, 1);
            }
            expandFrontierHorzEUp(dists, size, radius, row, col);
            
            for ( int i = 0; i < radius * 2; i++ )
            {
                for ( int j = 0; j < 2 * i + 1; j++ )
                {
                    if ( i < radiusN )
                    {
                        fill(row - i - 1, col - j - 1);
                    }
                    
                    if ( i < radiusSE )
                    {
                       fill(row + size + radiusSE - i - 1, col + size * 2 + radiusSE * 2 - j - 1);
                    }
                    
                    if ( i < radiusSW )
                    {
                       fill(row + size + radiusSW - i - 1, col - j - 1);
                    }
                }
            }
            
            // S
            //
            for ( int i = 0; i < (radiusS < size ? radiusS : size); i++ )
            {
                for ( int j = 0; j < (size - i) * 2 - 1; j++ )
                {
                    fill(row + size + i, col + size * 2 - j - 1);
                }
            }
            //
            if ( radiusS > size )
            {
                for ( int i = 0; i < radiusS - size; i++ )
                {
                    for ( int j = 0; j < i * 2 + 1; j++ )
                    {
                        fill(row + size * 2 + i, col + size * 2 + j);
                    }
                }
            }
            
            // NW
            //
            for ( int i = 0; i < (radiusNW < size ? radiusNW : size) * 2; i++ )
            {
                for ( int j = 0; j < size - (i + 1) / 2; j++ )
                {
                    fill(row + j, col - i - 1);
                }
            }
            //
            if ( radiusNW > size )
            {
                for ( int i = 0; i < radiusNW - size; i++ )
                {
                    for ( int j = 0; j < i * 2 + 1; j++ )
                    {
                        fill(row - radiusNW + size + i, col - radiusNW * 2 + j);
                    }
                }
            }
            
            // NE
            //
            for ( int i = 0; i < (radiusNE < size ? radiusNE : size) * 2; i++ )
            {
                for ( int j = 0; j < size - (i + 1) / 2; j++ )
                {
                    fill(row + j, col + i + j * 2 + 1);
                }
            }
            //
            if ( radiusNE > size )
            {
                for ( int i = 0; i < radiusNE - size; i++ )
                {
                    for ( int j = 0; j < i * 2 + 1; j++ )
                    {
                        fill(row - radiusNE + size + i, col + size * 2 + j);
                    }
                }
            }
            
            return;
        }
        
        if ( children != 0 )
        {
            children[0]->expandUp(radius, size / 2, row, col, depth + 1);
            children[1]->expandUp(radius, size / 2, row + size / 2, col, depth + 1);
            children[2]->expandUp(radius, size / 2, row + size / 2, col + size, depth + 1);
            
            children[3]->expandDown(radius, size / 2, row + size - 1, col + size - 1, depth + 1);
        }
    }
    
    void expandDown(int radius, int size, int row = 0, int col = 0, int depth = 0)
    {
        if ( (size <= radius + 1 && children != 0) || (children == 0 && color == 1) )
        {
            // all black
            
            for ( int i = 0; i < size; i++ )
            {
                for ( int j = 0; j < 2 * i + 1; j++ )
                {
                    fill(row - i, col - i * 2 + j);
                }
            }
            
            int radiusNW = 0;
            int radiusS = 0;
            int radiusNE = 0;
            int radiusSW = 0;
            int radiusN = 0;
            int radiusSE = 0;
            
            int dists[size];
            
            if ( children == 0 )
            {
                radiusNW = radius;
                radiusS = radius;
                radiusNE = radius;
                radiusSW = radius;
                radiusN = radius;
                radiusSE = radius;
            }
            else
            {
                radiusS = radius - distNup(size, 0, size);
                radiusNW = radius - distSWup(size, 0, size);
                radiusNE = radius - distSEup(size, 0, size);
                radiusN = radius - distNdown(size, 0, size);
                radiusSW = radius - distSWdown(size, 0, size);
                radiusSE = radius - distSEdown(size, 0, size);
            }
            
            // NNE
            //
            if ( children )
            {
                memset(&dists, 0, size * sizeof(int));
                frontierSSUp(dists, size, size, 1);
            }
            else
            {
                for ( int i = 0; i < size; i++ )
                {
                    dists[i] = size;
                }
            }
            expandFrontierDiagBase(dists, size, radius, row, col - size * 2, -1, 1);
            
            // NNW
            //
            if ( children )
            {
                memset(&dists, 0, size * sizeof(int));
                frontierSSUp(dists, size, size, -1);
            }
            expandFrontierDiagBase(dists, size, radius, row, col, -1, -1);
            
            // SSE
            //
            if ( children )
            {
                memset(&dists, 0, size * sizeof(int));
                frontierSSDown(dists, size, size, 1);
            }
            expandFrontierDiagTip(dists, size, radius, row, col, 1, 1);
            
            // SSW
            //
            if ( children )
            {
                memset(&dists, 0, size * sizeof(int));
                frontierSSDown(dists, size, size, -1);
            }
            expandFrontierDiagTip(dists, size, radius, row, col, 1, -1);
            
            // W
            //
            if ( children )
            {
                memset(&dists, 0, size * sizeof(int));
                frontierHorzDown(dists, size, size, -1, 1);
            }
            else
            {
                for ( int i = 0; i < size; i++ )
                {
                    dists[i] = size * 2 - i;
                }
            }
            expandFrontierHorzWDown(dists, size, radius, row, col);
            
            // E
            //
            if ( children )
            {
                memset(&dists, 0, size * sizeof(int));
                frontierHorzDown(dists, size, size, 1, 1);
            }
            expandFrontierHorzEDown(dists, size, radius, row, col);
            
            for ( int i = 0; i < size * 2; i++ )
            {
                for ( int j = 0; j < 2 * i + 1; j++ )
                {
                    if ( i < radiusS )
                    {
                        fill(row + i + 1, col + j + 1);
                    }
                    
                    if ( i < radiusNW )
                    {
                        fill(row - size - radiusNW + i + 1, col - (size + radiusNW) * 2 + j + 1);
                    }
                    
                    if ( i < radiusNE )
                    {
                        fill(row - size - radiusNE + i + 1, col + j + 1);
                    }
                }
            }
            
            // N
            //
            for ( int i = 0; i < (radiusN < size ? radiusN : size); i++ )
            {
                for ( int j = 0; j < (size - i) * 2 - 1; j++ )
                {
                    fill(row - size - i, col - size * 2 + j + 1);
                }
            }
            //
            if ( radiusN > size )
            {
                for ( int i = 0; i < radiusN - size; i++ )
                {
                    for ( int j = 0; j < i * 2 + 1; j++ )
                    {
                        fill(row - size * 2 - i, col - size * 2 - j);
                    }
                }
            }
            
            // SW
            //
            for ( int i = 0; i < (radiusSW < size ? radiusSW : size) * 2; i++ )
            {
                for ( int j = 0; j < size - (i + 1) / 2; j++ )
                {
                    fill(row - j, col - i - j * 2 - 1);
                }
            }
            //
            if ( radiusSW > size )
            {
                for ( int i = 0; i < radiusSW - size; i++ )
                {
                    for ( int j = 0; j < i * 2 + 1; j++ )
                    {
                        fill(row + radiusSW - size - i, col - size * 2 - j);
                    }
                }
            }
            
            // SE
            //
            for ( int i = 0; i < (radiusSE < size ? radiusSE : size) * 2; i++ )
            {
                for ( int j = 0; j < size - (i + 1) / 2; j++ )
                {
                    fill(row - j, col + i + 1);
                }
            }
            //
            if ( radiusSE > size )
            {
                for ( int i = 0; i < radiusSE - size; i++ )
                {
                    for ( int j = 0; j < i * 2 + 1; j++ )
                    {
                        fill(row + radiusSE - size - i, col + radiusSE * 2 - j);
                    }
                }
            }
            
            return;
        }
        
        if ( children != 0 )
        {
            children[0]->expandDown(radius, size / 2, row, col, depth + 1);
            children[1]->expandDown(radius, size / 2, row - size / 2, col - size, depth + 1);
            children[2]->expandDown(radius, size / 2, row - size / 2, col, depth + 1);
            
            children[3]->expandUp(radius, size / 2, row - size + 1, col - size + 1, depth + 1);
        }
    }
    
    void expandFrontierDiagBase(int * dists, int size, int radius, int row, int col, int dirrow, int dircol)
    {
        for ( int i = 0; i < radius; i++ )
        {
            int width = radius - i < size ? dists[radius - i - 1] == 0 ? 0 : (dists[radius - i - 1] + radius - i - 1) * 2 + 1 : size * 2;
            int min = i < size ? (size - i) * 2 - 1 : 0;
            int offset = (dirrow + dircol) * i + dircol;
            
            for ( int j = min; j < width && j < size * 2; j++ )
            {
                fill(row + dirrow * (size + i), col + offset + dircol * j);
            }
        }
    }
    
    void expandFrontierDiagTip(int * dists, int size, int radius, int row, int col, int dirrow, int dircol)
    {
        for ( int i = 0; i < radius; i++ )
        {
            int index = radius - i - 1;
            
            if ( index > size - 1 )
            {
                index = size - 1;
            }
            int width = dists[index] == 0 ? 0 : (dists[index] - size + radius - i - 1) * 2 + 1;
            //cout << i << '\t' << width << endl;
            int offset = (dirrow + dircol) * i + dircol + dirrow;
            
            for ( int j = 0; j < width && j < size * 2; j++ )
            {
                fill(row + dirrow * (1 + i), col + offset + dircol * j);
            }
        }
    }
    
    void expandFrontierHorzWUp(int * dists, int size, int radius, int row, int col)
    {
        for ( int i = 0; i < size; i++ )
        {
            int width = (dists[i] + radius - size * 2) * 2 + 1;
            
            for ( int j = (size - i) * 2; j < width; j += 2 )
            {
                fill(row + i, col - j);
                
                if ( i == 0 && j > (size - i) * 2 )
                {
                    fill(row, col - j + 1);
                }
                
                if ( i < size - 1 )
                {
                    fill(row + i + 1, col - j + 1);
                }
            }
        }
    }
    
    void expandFrontierHorzEUp(int * dists, int size, int radius, int row, int col)
    {
        for ( int i = 0; i < size; i++ )
        {
            int width = (dists[i] + radius - size * 3 + i) * 2 + 1;
            
            for ( int j = 0; j < width; j += 2 )
            {
                int colj = col + size * 2 + j;
                
                fill(row + i, colj);
                
                if ( i == 0 && j > 0 )
                {
                    fill(row, colj - 1);
                }
                
                if ( i < size - 1 )
                {
                    fill(row + i + 1, colj + 1);
                }
            }
        }
    }
    
    void expandFrontierHorzWDown(int * dists, int size, int radius, int row, int col)
    {
        for ( int i = 0; i < size; i++ )
        {
            int width = (dists[i] + radius - size * 2) * 2;
            
            for ( int j = 0; j < width; j += 2 )
            {
                int rowi = row - size + i + 1;
                int colj = col - size * 2 - j;
                
                fill(rowi, colj);
                
                if ( i == size - 1 && j > 0 )
                {
                    fill(rowi, colj + 1);
                }
                
                if ( i > 0 )
                {
                    fill(rowi - 1, colj - 1);
                }
            }
        }
    }
    
    void expandFrontierHorzEDown(int * dists, int size, int radius, int row, int col)
    {
        for ( int i = 0; i < size; i++ )
        {
            int width = (dists[i] + radius - size * 2) * 2;
            
            for ( int j = 0; j < width; j += 2 )
            {
                int rowi = row - size + i + 1;
                int colj = col + i * 2 + j + 2;
                
                fill(rowi, colj);
                
                if ( i == size - 1 && j > 0 )
                {
                    fill(rowi, colj - 1);
                }
                
                if ( i > 0 )
                {
                    fill(rowi - 1, colj - 1);
                }
            }
        }
    }
    
    void frontierSSUp(int * dists, int ndists, int size, int dir, int offrow = 0, int offcol = -1)
    {
        if ( offcol == -1 )
        {
            offcol = size;
        }
        
        if ( children != 0 )
        {
            children[dir == 1 ? 2 : 1]->frontierSSUp(dists, ndists, size / 2, dir, offrow, offcol);
            
            size /= 2;
            offcol -= size;
            
            if ( dists[offrow] < offcol )
            {
                children[dir == 1 ? 1 : 2]->frontierSSUp(dists, ndists, size, dir, offrow, offcol);
            }
            
            if ( dists[offrow] < offcol )
            {
                children[3]->frontierSSDown(dists, ndists, size, dir, offrow, offcol);
            }
            
            offrow += size;
            
            if ( dists[offrow] < offcol )
            {
                children[0]->frontierSSUp(dists, ndists, size, dir, offrow, offcol);
            }
        }
        else if ( color == 1 )
        {
            for ( int i = offrow; i < ndists && dists[i] < offcol; i++ )
            {
                dists[i] = offcol;
            }
        }
    }
    
    void frontierSSDown(int * dists, int ndists, int size, int dir, int offrow = 0, int offcol = -1)
    {
        if ( offcol == -1 )
        {
            offcol = size;
        }
        
        if ( children != 0 )
        {
            children[0]->frontierSSDown(dists, ndists, size / 2, dir, offrow, offcol);
            
            size /= 2;
            offrow += size;
            
            if ( dists[offrow] < offcol )
            {
                children[3]->frontierSSUp(dists, ndists, size, dir, offrow, offcol);
            }
            
            if ( dists[offrow] < offcol )
            {
                children[dir == 1 ? 2 : 1]->frontierSSDown(dists, ndists, size, dir, offrow, offcol);
            }
            
            offcol -= size;
            
            if ( dists[offrow] < offcol )
            {
                children[dir == 1 ? 1 : 2]->frontierSSDown(dists, ndists, size, dir, offrow, offcol);
            }
        }
        else if ( color == 1 )
        {
            for ( int i = offrow; i < ndists && dists[i] < offcol; i++ )
            {
                dists[i] = offcol;
            }
        }
    }
    
    void frontierNNUp(int * dists, int ndists, int size, int dir, int offrow = 0, int offcol = -1)
    {
        if ( offcol == -1 )
        {
            offcol = size;
        }
        
        if ( children != 0 )
        {
            children[dir == 1 ? 2 : 1]->frontierNNUp(dists, ndists, size / 2, dir, offrow, offcol);
            
            size /= 2;
            offcol -= size;
            
            if ( dists[offrow] < offcol )
            {
                children[dir == 1 ? 1 : 2]->frontierNNUp(dists, ndists, size, dir, offrow, offcol);
            }
            
            if ( dists[offrow] < offcol )
            {
                children[3]->frontierNNDown(dists, ndists, size, dir, offrow, offcol);
            }
            
            offrow += size;
            
            if ( dists[offrow] < offcol )
            {
                children[0]->frontierNNUp(dists, ndists, size, dir, offrow, offcol);
            }
        }
        else if ( color == 1 )
        {
            for ( int i = offrow; i < ndists && dists[i] < offcol; i++ )
            {
                dists[i] = offcol;
            }
        }
    }
    
    void frontierNNDown(int * dists, int ndists, int size, int dir, int offrow = 0, int offcol = 0)
    {
        if ( children != 0 )
        {
            children[0]->frontierNNDown(dists, ndists, size / 2, offrow, offcol);
            
            size /= 2;
            offrow += size;
            
            if ( dists[offrow] < offcol )
            {
                children[3]->frontierNNUp(dists, ndists, size, dir, offrow, offcol);
            }
            
            if ( dists[offrow] < offcol )
            {
                children[dir == 1 ? 2 : 1]->frontierNNDown(dists, ndists, size, dir, offrow, offcol);
            }
            
            offcol -= size;
            
            if ( dists[offrow] < offcol )
            {
                children[dir == 1 ? 1 : 2]->frontierNNDown(dists, ndists, size, dir, offrow, offcol);
            }
        }
        else if ( color == 1 )
        {
            for ( int i = offrow; i < ndists && dists[i] < offcol; i++ )
            {
                dists[i] = offcol;
            }
        }
    }
    
    void frontierHorzUp(int * dists, int ndists, int size, int dir, int offrow = -1, int offcol = 0)
    {
        if ( offrow == -1 )
        {
            offrow = size;
            offcol = size;
        }
        
        if ( children != 0 )
        {
            children[dir == 1 ? 2 : 1]->frontierHorzUp(dists, ndists, size / 2, dir, offrow, offcol);
            
            size /= 2;
            offrow -= size;
            
            if ( dists[offrow - 1] < offcol + ndists )
            {
                children[dir == 0]->frontierHorzUp(dists, ndists, size, dir, offrow, offcol);
            }
            
            if ( dists[offrow - 1] < offcol + ndists )
            {
                children[3]->frontierHorzDown(dists, ndists, size, dir, offrow, offcol);
            }
            
            offrow += size;
            offcol -= size;
            
            if ( dists[offrow - 1] < offcol + ndists )
            {
                children[dir == 1 ? 1 : 2]->frontierHorzUp(dists, ndists, size, dir, offrow, offcol);
            }
        }
        else if ( color == 1 )
        {
            for ( int i = offrow - 1; i >= 0 && dists[i] < offcol + ndists; i-- )
            {
                dists[i] = offcol + ndists;
            }
            for ( int i = offrow; i < ndists && dists[i] < offcol + ndists - offrow + i - 1; i++ )
            {
                dists[i] = offcol + ndists + offrow - i - 1;
            }
        }
    }
    
    void frontierHorzDown(int * dists, int ndists, int size, int dir, int offrow = 0, int offcol = -1)
    {
        if ( offcol == -1 )
        {
            offcol = size;
        }
        
        if ( children != 0 )
        {
            children[dir == 1 ? 2 : 1]->frontierHorzDown(dists, ndists, size / 2, dir, offrow, offcol);
            
            size /= 2;
            offrow += size;
            offcol -= size;
            
            if ( dists[offrow - 1] < offcol + ndists )
            {
                children[0]->frontierHorzDown(dists, ndists, size, dir, offrow, offcol);
            }
            
            if ( dists[offrow - 1] < offcol + ndists )
            {
                children[3]->frontierHorzUp(dists, ndists, size, dir, offrow, offcol);
            }
            
            offrow -= size;
            
            if ( dists[offrow - 1] < offcol + ndists )
            {
                children[dir == 1 ? 1 : 2]->frontierHorzDown(dists, ndists, size, dir, offrow, offcol);
            }
        }
        else if ( color == 1 )
        {
            for ( int i = offrow - 1; i >= 0 && dists[i] < offcol + ndists; i-- )
            {
                dists[i] = offcol + ndists;
            }
            
            for ( int i = offrow; i < ndists && dists[i] < offcol + ndists + offrow - i; i++ )
            {
                dists[i] = offcol + ndists + offrow - i - 1;
            }
        }
    }
    
    int distNup(int size, int offset, int dist)
    {
        if ( children != 0 )
        {
            dist = children[0]->distNup(size / 2, offset, dist);
            
            if ( offset + size / 2 < dist )
            {
                dist = children[1]->distNup(size / 2, offset + size / 2, dist);
                dist = children[2]->distNup(size / 2, offset + size / 2, dist);
                dist = children[3]->distNdown(size / 2, offset + size / 2, dist);
            }
        }
        else
        {
            if ( color == 1 ) dist = offset;
        }
        
        return dist;
    }
    
    int distNdown(int size, int offset, int dist)
    {
        if ( children != 0 )
        {
            dist = children[1]->distNdown(size / 2, offset, dist);
            dist = children[2]->distNdown(size / 2, offset, dist);
            
            dist = children[3]->distNup(size / 2, offset, dist);
            
            if ( offset + size / 2 < dist )
            {
                dist = children[0]->distNdown(size / 2, offset + size / 2, dist);
            }
        }
        else
        {
            if ( color == 1 ) dist = offset;
        }
        
        return dist;
    }
    
    int distSWup(int size, int offset, int dist)
    {
        if ( children != 0 )
        {
            dist = children[1]->distSWup(size / 2, offset, dist);
            
            if ( offset + size / 2 < dist )
            {
                dist = children[0]->distSWup(size / 2, offset + size / 2, dist);
                dist = children[2]->distSWup(size / 2, offset + size / 2, dist);
                dist = children[3]->distSWdown(size / 2, offset + size / 2, dist);
            }
        }
        else
        {
            if ( color == 1 ) dist = offset;
        }
        
        return dist;
    }
    
    int distSWdown(int size, int offset, int dist)
    {
        if ( children != 0 )
        {
            dist = children[0]->distSWdown(size / 2, offset, dist);
            dist = children[1]->distSWdown(size / 2, offset, dist);
            
            dist = children[3]->distSWup(size / 2, offset, dist);
            
            if ( offset + size / 2 < dist )
            {
                dist = children[2]->distSWdown(size / 2, offset + size / 2, dist);
            }
        }
        else
        {
            if ( color == 1 ) dist = offset;
        }
        
        return dist;
    }
    
    int distSEup(int size, int offset, int dist)
    {
        if ( children != 0 )
        {
            dist = children[2]->distSEup(size / 2, offset, dist);
            
            if ( offset + size / 2 < dist )
            {
                dist = children[0]->distSEup(size / 2, offset + size / 2, dist);
                dist = children[1]->distSEup(size / 2, offset + size / 2, dist);
                dist = children[3]->distSEdown(size / 2, offset + size / 2, dist);
            }
        }
        else
        {
            if ( color == 1 ) dist = offset;
        }
        
        return dist;
    }
    
    int distSEdown(int size, int offset, int dist)
    {
        if ( children != 0 )
        {
            dist = children[0]->distSEdown(size / 2, offset, dist);
            dist = children[2]->distSEdown(size / 2, offset, dist);
            
            dist = children[3]->distSEup(size / 2, offset, dist);
            
            if ( offset + size / 2 < dist )
            {
                dist = children[1]->distSEdown(size / 2, offset + size / 2, dist);
            }
        }
        else
        {
            if ( color == 1 ) dist = offset;
        }
        
        return dist;
    }
    
    void drawUp(string & s, string col, int size, int r = 0, int c = 0, int depth = 0)
    {
        if ( children != 0 )
        {
            size /= 2;
            
            children[0]->drawUp(s, col, size, r, c, depth + 1);
            children[1]->drawUp(s, col, size, r + size, c, depth + 1);
            children[2]->drawUp(s, col, size, r + size, c + size * 2, depth + 1);
            
            children[3]->drawDown(s, col, size, r + size * 2, c + size * 2, depth + 1);
        }
        else
        {
            for ( int i = 0; i < depth; i++ )
            {
                s += "  ";
            }
            
            float scale = 1.2;
            float xoff = rows * scale;
            float x = xoff + (float)c * scale - (float)r * scale;
            float y = 50 + (float)r * 1.732 * scale;
            float w = (float)size * scale;
            float h = (float)size * 1.732 * scale;
            
            s += "<polygon points=\"" + to_string(x) + "," + to_string(y) + " " +
                to_string(x - w) + "," + to_string(y + h) + " " +
                to_string(x + w) + "," + to_string(y + h) +
                "\" fill-opacity=\"0.4\" style=\"stroke:black;stroke-width:" + to_string(1./depth) + ";fill:" + (color == 1 ? col : "white") +
                ";\" />\n";
        }
    }
    
    void drawDown(string & s, string col, int size, int r = 0, int c = 0, int depth = 0)
    {
        if ( children != 0 )
        {
            size /= 2;
            
            children[0]->drawDown(s, col, size, r, c, depth + 1);
            children[1]->drawDown(s, col, size, r - size, c - size * 2, depth + 1);
            children[2]->drawDown(s, col, size, r - size, c, depth + 1);
            
            children[3]->drawUp(s, col, size, r - size * 2, c - size * 2, depth + 1);
        }
        else
        {
            for ( int i = 0; i < depth; i++ )
            {
                s += "  ";
            }
            
            float scale = 1.2;
            float xoff = rows * scale;
            float x = xoff + (float)c * scale - (float)r * scale;
            float y = 50 + (float)r * 1.732 * scale;
            float w = (float)size * scale;
            float h = (float)size * 1.732 * scale;
            
            s += "<polygon points=\"" + to_string(x) + "," + to_string(y) + " " +
                to_string(x + w) + "," + to_string(y - h) + " " +
                to_string(x - w) + "," + to_string(y - h) +
                "\" fill-opacity=\"0.4\" style=\"stroke:black;stroke-width:" + to_string(1./depth) + ";fill:" + (color == 1 ? col : "white") +
                ";\" />\n";
        }
    }
};

int ** loadBmp(const char * file, int & width, int & height)
{
    #define HEADER_SIZE 54;

    std::ifstream bmp(file, std::ios::binary);

    char header[54];
    bmp.read(header, 54);

    uint32_t fileSize = *reinterpret_cast<uint32_t *>(&header[2]);
    uint32_t dataOffset = *reinterpret_cast<uint32_t *>(&header[10]);
    width = *reinterpret_cast<uint32_t *>(&header[18]);
    height = *reinterpret_cast<uint32_t *>(&header[22]);
    uint16_t depth = *reinterpret_cast<uint16_t *>(&header[28]);

    if ( width < 0 ) width = -width;
    if ( height < 0 ) height = -height;
    
    char img[fileSize];
    bmp.read(img, dataOffset - 54);
    
    int dataSize = ((width * 3 + 3) & (~3)) * height;
    bmp.read(img, dataSize);

    char temp = 0;

    int x = 0;
    int y = 0;
    
    int ** map = new int * [height];
    map[0] = new int[width];
    
    for (int i = dataSize - 6; i >= 1; i -= 3)
    {
        temp = img[i];
        img[i] = img[i+2];
        img[i+2] = temp;

        int r = int(img[i] & 0xff);
        int g = int(img[i+1] & 0xff);
        int b = int(img[i+2] & 0xff);
        
        map[y][width - x - 1] = r < 255 || g < 255 || b < 255;
        
        x++;
        
        if ( x == width )
        {
            x = 0;
            y++;
            map[y] = new int[width];
            i -= 3;
        }
    }

    for ( int i = 0; i < height; i++ )
    {
        for ( int j = 0; j < width; j++ )
        {
            //cout << map[i][j];
        }
        
        //cout << endl;
    }
    
    return map;
}

void triangulate(int ** map, int w, int h, int & rows)
{
    int horz = w + (float)h / 0.866;
    int vert = h + (float)w;
    
    int min = horz < vert ? horz : vert;
    
    rows = 1;
    
    while ( rows < min )
    {
        rows *= 2;
    }
    
    //cout << w << '\t' << h << '\t' << min << '\t' << rows << endl;
    
    tri = new int * [rows];
    exp = new int * [rows];
    
    for ( int i = 0; i < rows; i++ )
    {
        tri[i] = new int[2 * i + 1];
        exp[i] = new int[2 * i + 1];
        memset(tri[i], 0, (2 * i + 1) * sizeof(int));
        memset(exp[i], 0, (2 * i + 1) * sizeof(int));
    }
     
     int cols = rows * 2 - 1;
     float factor = (float)min / rows / 0.866;
     
     for ( int i = h; i < rows; i++ )
     {
        for ( int j = 0; j < cols - i; j++)
        {
            //cout << " ";
        }
        
        for ( int j = 0; j < 2 * (i + 1) - 1; j++ )
        {
            int x = (float)w / 2 - (float)i / 1.732 * factor + (float)j / 1.732 * factor;
            int y = (i - (rows - h) * factor) * factor;
            
            if
            (
                x < 0 ||
                x >= w ||
                y < 0 ||
                y >= h
            ) {continue;};
            
            tri[i][j] = map[y][x];
            
            //cout << map[y][x];// << ' ';
        }
        
        //cout << endl;
     }
}

int main(int argc, const char ** argv)
{
    if ( argc == 1 )
    {
        cerr << "Usage:" << endl;
        cerr << "   triangle <img.bmp>           Run time trials" << endl;
        cerr << "   triangle <img.bmp> <radius>  Run once for radius" << endl; 
        exit(0);
    }
    
    int rmin;
    int rmax;
    int iters;
    
    if ( argc == 2 )
    {
        rmin = 1;
        rmax = 33;
        iters = 10000;
        cout << "Radius\tB2 (s)\tB1 (s)" << endl;
    }
    else
    {
        rmin = atoi(argv[2]);
        rmax = rmin + 1;
        iters = 1;
    }
    
    for ( int r = rmin; r < rmax; r++ )// *= 2 )
    {
        int wbmp;
        int hbmp;

        int ** image = loadBmp(argv[1], wbmp, hbmp);
        triangulate(image, wbmp, hbmp, rows);

        Node * tree = new Node();

        tree->buildUp(tri, rows);
    
        std::time_t time = std::time(nullptr);
        for ( int i = 0; i < iters; i++ )
        {
            tree->expandUp(r, rows);
        }
        
        if ( iters > 1 )
        {
            time_t timenew = std::time(nullptr);
            for ( int i = 0; i < iters; i++ )
            {
                expandNaive(r);
            }
            cout << r << '\t' << timenew - time << '\t' << std::time(nullptr) - timenew << endl;
        }
        
        Node * tree2 = new Node();

        tree2->buildUp(exp, rows);

        string head = "<svg version=\"1.1\" baseProfile=\"full\" width=\"640\" height=\"640\" xmlns=\"http://www.w3.org/2000/svg\">\n";
        string foot = "</svg>";
        string s;
        string t;
        tree->drawUp(s, "blue", rows);
        tree2->drawUp(t, "red", rows);
    
        ofstream out;
        out.open("out-" + to_string(r) + "-orig.svg");
        out << head << s << foot;
        out.close();
    
        out.open("out-" + to_string(r) + "-grow.svg");
        out << head << t << foot;
        out.close();
    
        out.open("out-" + to_string(r) + "-comb.svg");
        out << head << t << s << foot;
        out.close();
    }
}
