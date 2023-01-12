/*
 *  ParticleContainer.h
 *
 *  Created on: 5.12.2022
 *      Author: borisov
 */

#pragma once

#include <vector>
#include <array>
#include "ArrayUtils.h"

/**
 * @brief namespace that includes utils for the particle container
 */
namespace PContainer
{
    // Predeclarations
    std::vector<int> getSurroundingZ(int index, std::array<int, 3> &numCells);

    std::vector<int> getHaloX(int index, std::array<int, 3> &numCells);

    std::vector<int> getHaloY(int index, std::array<int, 3> &numCells);

    std::vector<int> getHaloZ(int index, std::array<int, 3> &numCells);

    std::array<int, 3> convert1DTo3D(int index, std::array<int, 3> &numCells);

    int convert3DTo1D(std::array<int, 3> index, std::array<int, 3> &numCells);

    std::vector<std::array<double, 3>> getMirroringOffsets(int currentCell, std::vector<int> &boundaries, std::array<int, 3> &numCells, std::array<double, 3> &cellSize);

    std::array<int, 3> applyBoundaries(std::array<int, 3> &index3D, std::array<int, 3> &numCells, std::vector<int> boundaries);

    int getMirroredHaloIndex(int currentCellIndex3D, int numCell);

    /**
     * @brief calculates the indices of the neighbours for a given cell, considering Newton's 3rd Law
     *
     * @param index index of the cell
     * @param numCells dimension of the cell array
     * @return indices of the neighbours
     */
    inline std::vector<int> getNeighboursNewton(int index, std::array<int, 3> &numCells)
    {
        std::vector<int> result;

        // transform index to 3d, so its easier to find the index of the cell right behind
        std::array<int, 3> index3D = convert1DTo3D(index, numCells);

        // all neighbours in the front layer
        std::vector<int> frontLayerAll = getSurroundingZ(index, numCells);

        // filter out smaller indices (according to newtons 3rd law)
        for (int i : frontLayerAll)
        {
            if (i > index)
                result.push_back(i);
        }

        // all neighbours in the back layer
        // z out of bounds
        if (index3D[2] + 1 >= numCells[2])
            return result;

        // z + 1
        index3D[2] = index3D[2] + 1;

        // get the 1 dimensional index of the cell right behind
        int backLayerIndex = convert3DTo1D(index3D, numCells);
        // add it right away to the result
        result.push_back(backLayerIndex);

        // add all neighbours of the cell right behind to the result
        for (int i : getSurroundingZ(backLayerIndex, numCells))
        {
            result.push_back(i);
        }
        return result;
    }


    /**
     * @brief calculates indices of neighbouring halo cells at a periodic boundary
     * @param index 1D index of cell
     * @param numCells number of cells in each dimension
     * @param pb positions of periodic boundaries
     * @return 1D indices of halo cells neighbouring given cell at a periodic boundary
    */
    inline std::vector<int> getPeriodicHaloNeighbours(int index, std::array<int,3> &numCells, std::vector<int> &pb) {
        std::array<int, 3> index3D = convert1DTo3D(index, numCells);
        std::set<int> result;

        //cell at front periodic boundary 
        if (index3D[2] == 1 && std::find(pb.begin(), pb.end(), 4) != pb.end()) {
            int tmp_index = convert3DTo1D({index3D[0], index3D[1], 0}, numCells);
            std::vector<int> halo = getHaloZ(tmp_index, numCells);
            std::copy(halo.begin(), halo.end(), std::inserter(result, result.end()));
        }
        // cell at back periodic boundary
        if (index3D[2] == numCells[2] - 2 && std::find(pb.begin(), pb.end(), 5) != pb.end())
        {
            int tmp_index = convert3DTo1D({index3D[0], index3D[1], index3D[2] + 1}, numCells);
            std::vector<int> halo = getHaloZ(tmp_index, numCells);
            std::copy(halo.begin(), halo.end(), std::inserter(result, result.end()));
        }
        // cell at left periodic boundary
        if (index3D[0] == 1 && std::find(pb.begin(), pb.end(), 0) != pb.end())
        {
            int tmp_index = convert3DTo1D({0, index3D[1], index3D[2]}, numCells);
            std::vector<int> halo = getHaloX(tmp_index, numCells);
            std::copy(halo.begin(), halo.end(), std::inserter(result, result.end()));
        }
        // cell at right periodic boundary
        if (index3D[0] == numCells[0] - 2 && std::find(pb.begin(), pb.end(), 1) != pb.end())
        {
            int tmp_index = convert3DTo1D({index3D[0] + 1, index3D[1], index3D[2]}, numCells);
            std::vector<int> halo = getHaloX(tmp_index, numCells);
            std::copy(halo.begin(), halo.end(), std::inserter(result, result.end()));
        }
        // cell at bottom periodic boundary
        if (index3D[1] == 1 && std::find(pb.begin(), pb.end(), 2) != pb.end())
        {
            int tmp_index = convert3DTo1D({index3D[0], 0, index3D[2]}, numCells);
            std::vector<int> halo = getHaloY(tmp_index, numCells);
            std::copy(halo.begin(), halo.end(), std::inserter(result, result.end()));
        }
        // cell at top periodic boundary
        if (index3D[1] == numCells[1] - 2 && std::find(pb.begin(), pb.end(), 3) != pb.end())
        {
            int tmp_index = convert3DTo1D({index3D[0], index3D[1] + 1, index3D[2]}, numCells);
            std::vector<int> halo = getHaloY(tmp_index, numCells);
            std::copy(halo.begin(), halo.end(), std::inserter(result, result.end()));
        }

        std::vector<int> res(result.begin(), result.end());
        return res;
    }


    /**
     * @brief calculates indices of neighbouring halo cells of given cell
     * @param index 1D index of cell
     * @param numCells number of cells in each dimension
     * @return 1D indices of halo cells neighbouring given cell
     */
    inline std::vector<int> getHaloNeighbours(int index, std::array<int, 3> &numCells)
    {
        std::set<int> result;
        std::array<int, 3> index3D = convert1DTo3D(index, numCells);

        // cell at front boundary
        if (index3D[2] == 1)
        {
            int tmp_index = convert3DTo1D({index3D[0], index3D[1], 0}, numCells);
            std::vector<int> halo = getHaloZ(tmp_index, numCells);
            std::copy(halo.begin(), halo.end(), std::inserter(result, result.end()));
        }
        // cell at back boundary
        if (index3D[2] == numCells[2] - 2)
        {
            int tmp_index = convert3DTo1D({index3D[0], index3D[1], index3D[2] + 1}, numCells);
            std::vector<int> halo = getHaloZ(tmp_index, numCells);
            std::copy(halo.begin(), halo.end(), std::inserter(result, result.end()));
        }
        // cell at left boundary
        if (index3D[0] == 1)
        {
            int tmp_index = convert3DTo1D({0, index3D[1], index3D[2]}, numCells);
            std::vector<int> halo = getHaloX(tmp_index, numCells);
            std::copy(halo.begin(), halo.end(), std::inserter(result, result.end()));
        }
        // cell at right boundary
        if (index3D[0] == numCells[0] - 2)
        {
            int tmp_index = convert3DTo1D({index3D[0] + 1, index3D[1], index3D[2]}, numCells);
            std::vector<int> halo = getHaloX(tmp_index, numCells);
            std::copy(halo.begin(), halo.end(), std::inserter(result, result.end()));
        }
        // cell at bottom boundary
        if (index3D[1] == 1)
        {
            int tmp_index = convert3DTo1D({index3D[0], 0, index3D[2]}, numCells);
            std::vector<int> halo = getHaloY(tmp_index, numCells);
            std::copy(halo.begin(), halo.end(), std::inserter(result, result.end()));
        }
        // cell at top boundary
        if (index3D[1] == numCells[1] - 2)
        {
            int tmp_index = convert3DTo1D({index3D[0], index3D[1] + 1, index3D[2]}, numCells);
            std::vector<int> halo = getHaloY(tmp_index, numCells);
            std::copy(halo.begin(), halo.end(), std::inserter(result, result.end()));
        }

        std::vector<int> res(result.begin(), result.end());
        return res;
    }

    /**
     * @brief calculates the indices of the neighbours for a given cell, considering Newton's 3rd Law
     *
     * @param index index of the cell
     * @param numCells dimension of the cell array
     * @return indices of the neighbours
     */
    inline std::vector<int> getDomainNeighboursNewton(int index, std::array<int, 3> &numCells)
    {
        std::vector<int> result;

        // transform index to 3d, so its easier to find the index of the cell right behind
        std::array<int, 3> index3D = convert1DTo3D(index, numCells);

        // all neighbours in the front layer
        std::vector<int> frontLayerAll = getSurroundingZ(index, numCells);

        // filter out smaller indices (according to newtons 3rd law)
        for (int i : frontLayerAll)
        {
            if (i > index)
                result.push_back(i);
        }

        // all neighbours in the back layer
        // z out of bounds (pointing to halo cell)
        if (index3D[2] + 1 >= numCells[2] - 1)
            return result;

        // z + 1
        index3D[2] = index3D[2] + 1;

        // get the 1 dimensional index of the cell right behind
        int backLayerIndex = convert3DTo1D(index3D, numCells);
        // add it right away to the result
        result.push_back(backLayerIndex);

        // add all neighbours of the cell right behind to the result
        for (int i : getSurroundingZ(backLayerIndex, numCells))
        {
            result.push_back(i);
        }
        return result;
    }

    /**
     * @brief calculates the indices of all surrounding neighbours for a given cell
     *
     * @param index index of the cell
     * @param numCells dimension of the cell array
     * @return indices of the neighbours
     */
    inline std::vector<int> getSurroundingZ(int index, std::array<int, 3> &numCells)
    {
        std::vector<int> result;

        // transform index to 3 dim index for easier use
        // but we will only use x and y in this method
        std::array<int, 3> index3D = convert1DTo3D(index, numCells);

        // get all surrounding cells
        // from bottom left to top right
        std::array<int, 3> tmpIndex3D;
        for (int y = index3D[1] - 1; y <= index3D[1] + 1; y++)
        {
            // y out of bounds
            if (y < 1 || y > numCells[1] - 2)
                continue;

            for (int x = index3D[0] - 1; x <= index3D[0] + 1; x++)
            {
                // x out of bounds
                if (x < 1 || x > numCells[0] - 2)
                    continue;

                tmpIndex3D[0] = x;
                tmpIndex3D[1] = y;
                tmpIndex3D[2] = index3D[2];

                int n = convert3DTo1D(tmpIndex3D, numCells);
                // index is not neighbour of itself
                if (n != index)
                    result.push_back(n);
            }
        }
        return result;
    }

    /**
     * @brief calculates indices of neighboring halo cells in x-plane
     * @param index 1D index of cell
     * @param numCells number of cells in each dimension
     * @return 1D indices of neighbouring halo cells in x-plane
     */
    inline std::vector<int> getHaloX(int index, std::array<int, 3> &numCells)
    {
        std::vector<int> result;
        std::array<int, 3> index3D = convert1DTo3D(index, numCells);

        std::array<int, 3> tmpIndex3D;
        for (int y = index3D[1] - 1; y <= index3D[1] + 1; y++)
        {
            for (int z = index3D[2] - 1; z <= index3D[2] + 1; z++)
            {
                tmpIndex3D[0] = index3D[0];
                tmpIndex3D[1] = y;
                tmpIndex3D[2] = z;
                int n = convert3DTo1D(tmpIndex3D, numCells);
                result.push_back(n);
            }
        }
        return result;
    }

    /**
     * @brief calculates indices of neighboring halo cells in y-plane
     * @param index 1D index of cell
     * @param numCells number of cells in each dimension
     * @return 1D indices of neighbouring halo cells in y-plane
     */
    inline std::vector<int> getHaloY(int index, std::array<int, 3> &numCells)
    {
        std::vector<int> result;
        std::array<int, 3> index3D = convert1DTo3D(index, numCells);

        std::array<int, 3> tmpIndex3D;
        for (int x = index3D[0] - 1; x <= index3D[0] + 1; x++)
        {
            for (int z = index3D[2] - 1; z <= index3D[2] + 1; z++)
            {
                tmpIndex3D[0] = x;
                tmpIndex3D[1] = index3D[1];
                tmpIndex3D[2] = z;
                int n = convert3DTo1D(tmpIndex3D, numCells);
                result.push_back(n);
            }
        }
        return result;
    }

    /**
     * @brief calculates indices of neighboring halo cells in z-plane
     * @param index 1D index of cell
     * @param numCells number of cells in each dimension
     * @return 1D indices of neighbouring halo cells in z-plane
     */
    inline std::vector<int> getHaloZ(int index, std::array<int, 3> &numCells)
    {
        std::vector<int> result;
        std::array<int, 3> index3D = convert1DTo3D(index, numCells);

        std::array<int, 3> tmpIndex3D;
        for (int y = index3D[1] - 1; y <= index3D[1] + 1; y++)
        {
            for (int x = index3D[0] - 1; x <= index3D[0] + 1; x++)
            {
                tmpIndex3D[0] = x;
                tmpIndex3D[1] = y;
                tmpIndex3D[2] = index3D[2];
                int n = convert3DTo1D(tmpIndex3D, numCells);
                result.push_back(n);
            }
        }
        return result;
    }

    /*
     * Helper Functions
     */

    /**
     * @brief converts a 1 dimensional index to a 3 dimensional index
     *
     * @param index 1 dimensional index
     * @param numCells dimension of the array
     * @return 3 dimensional index
     */
    inline std::array<int, 3> convert1DTo3D(int index, std::array<int, 3> &numCells)
    {
        std::array<int, 3> result;

        // x
        result[0] = index % numCells[0];
        // y
        result[1] = (index / numCells[0]) % numCells[1];
        // z
        result[2] = index / (numCells[0] * numCells[1]);

        return result;
    }

    /**
     * @brief converts a 3 dimensional index to a 1 dimensional index
     *
     * @param index 3 dimensional index
     * @param numCells dimension fo the array
     * @return 1 dimensional index
     */
    inline int convert3DTo1D(std::array<int, 3> index, std::array<int, 3> &numCells)
    {
        return index[0] + (index[1] * numCells[0]) + (index[2] * numCells[0] * numCells[1]);
    }

    /**
     * @brief calculates which border of the cell the boundary crossed
     *
     * @param prev_cell the cell where the particle was
     * @param new_cell the cell where the particle currently is
     * @param numCells the number of cells for the linked cell
     * @returns a vector with all boundaries that are crossed
     */
    inline std::vector<int> crossedBoundary(int prev_cell, int new_cell, std::array<int, 3> &numCells)
    {
        std::array<int, 3> prev_cell_3D = convert1DTo3D(prev_cell, numCells);
        std::array<int, 3> new_cell_3D = convert1DTo3D(new_cell, numCells);

        std::vector<int> boundaries = {};

        int boundary = 0;
        for (int i = 0; i < 3; i++)
        {
            if (new_cell_3D[i] - prev_cell_3D[i] == -1)
            {
                boundaries.push_back(boundary);
                boundary++;
            }
            else
                boundary++;
            if (new_cell_3D[i] - prev_cell_3D[i] == 1)
            {
                boundaries.push_back(boundary);
                boundary++;
            }
            else
                boundary++;
        }
        if (boundaries.size() == 0)
            throw std::invalid_argument("Cells not adjacent, cannot compute crossed boundary between " + std::to_string(prev_cell) + " and " + std::to_string(new_cell));

        return boundaries;
    }

    /**
     * @brief calculates all the mirroring offsets for the periodic boundary for a specific cell
     *
     * @param currentCell the cell from which particles are mirrored
     * @param boundaries the boundaries of the cell
     * @param numCells the number of cells for the linked cell
     * @param cellSize the size off the cells in each dimension
     * @returns a vector of arrays with mirroring positions in all dimensions
     */
    inline std::vector<std::array<double, 3>> getMirroringOffsets(int currentCell, std::vector<int> &boundaries, std::array<int, 3> &numCells, std::array<double, 3> &cellSize)
    {
        std::vector<std::array<int, 3>> cellsToMirror;
        auto index3D = convert1DTo3D(currentCell, numCells);

        // hardcoded as it is more efficient and there no more than 3 dimensions
        if (boundaries.size() == 1)
        {
            cellsToMirror.emplace_back(applyBoundaries(index3D, numCells, boundaries));
        }
        else if (boundaries.size() == 2)
        {
            cellsToMirror.emplace_back(applyBoundaries(index3D, numCells, {boundaries[0]}));
            cellsToMirror.emplace_back(applyBoundaries(index3D, numCells, {boundaries[1]}));
            cellsToMirror.emplace_back(applyBoundaries(index3D, numCells, boundaries));
        }
        else if (boundaries.size() == 3)
        {
            cellsToMirror.emplace_back(applyBoundaries(index3D, numCells, {boundaries[0]}));
            cellsToMirror.emplace_back(applyBoundaries(index3D, numCells, {boundaries[1]}));
            cellsToMirror.emplace_back(applyBoundaries(index3D, numCells, {boundaries[2]}));
            cellsToMirror.emplace_back(applyBoundaries(index3D, numCells, {boundaries[0], boundaries[1]}));
            cellsToMirror.emplace_back(applyBoundaries(index3D, numCells, {boundaries[1], boundaries[2]}));
            cellsToMirror.emplace_back(applyBoundaries(index3D, numCells, {boundaries[0], boundaries[2]}));
            cellsToMirror.emplace_back(applyBoundaries(index3D, numCells, boundaries));
        }

        std::vector<std::array<double, 3>> mirroringOffsets;

        for (auto &cellToMirror : cellsToMirror)
        {
            std::array<double, 3> mirroringOffset;
            mirroringOffset[0] = cellSize[0] * (cellToMirror[0] - index3D[0]);
            mirroringOffset[1] = cellSize[1] * (cellToMirror[1] - index3D[1]);
            mirroringOffset[2] = cellSize[2] * (cellToMirror[2] - index3D[2]);

            mirroringOffsets.emplace_back(mirroringOffset);
        }

        return mirroringOffsets;
    }

    /**
     * @brief applies all mirroring boundaries to the right dimension
     *
     * @param index3D the 3D index of the cell
     * @param numCells the number of cells in each dimension
     * @param boundaries the boundaries of the current cell
     * @returns an array with cell setoff in each direction
     */
    inline std::array<int, 3> applyBoundaries(std::array<int, 3> &index3D, std::array<int, 3> &numCells, std::vector<int> boundaries)
    {
        std::array<int, 3> newIndex;
        newIndex[0] = index3D[0];
        newIndex[1] = index3D[1];
        newIndex[2] = index3D[2];

        for (auto b : boundaries)
        {
            if (b < 2)
            {
                newIndex[0] = getMirroredHaloIndex(index3D[0], numCells[0]);
            }
            else if (b < 4)
            {
                newIndex[1] = getMirroredHaloIndex(index3D[1], numCells[1]);
            }
            else if (b < 6)
            {
                newIndex[2] = getMirroredHaloIndex(index3D[2], numCells[2]);
            }
        }

        return newIndex;
    }

    /**
     * @brief calculates where  in the halo the cell should be mirrored to
     *
     * @param currentCellIndex3D the cell index of the current cell in 3D
     * @param numCell the number of cells in each dimension
     * @returns the index of the mirrored cell
     */
    inline int getMirroredHaloIndex(int currentCellIndex3D, int numCell)
    {
        if (currentCellIndex3D == numCell - 2)
            return 0;
        else
            return numCell - 1;
    }
}