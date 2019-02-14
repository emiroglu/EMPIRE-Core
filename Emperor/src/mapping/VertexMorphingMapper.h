/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Munich
 *
 *  All rights reserved.
 *
 *  This file is part of EMPIRE.
 *
 *  EMPIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EMPIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EMPIRE.  If not, see http://www.gnu.org/licenses/.
 */
/******************************************************************************//**
 * \file VertexMorphingMapper.h
 * The header file of class MortarMapper.
 * \date 14/02/2019
 * \author Altug Emiroglu
 *********************************************************************************/
#ifndef VERTEXMORPHINGMAPPER_H_
#define VERTEXMORPHINGMAPPER_H_

#include "MathLibrary.h"
#include <vector>
#include <set>
#include <map>
#include "AbstractMapper.h"

class ANNkd_tree;

namespace flann {
template<typename Distance> class Index;
template<class T> struct L2;
template<typename T> class Matrix;
}

enum EMPIRE_VMM_FilterType {
    EMPIRE_VMM_HatFilter,
    EMPIRE_VMM_GaussianFilter
};

namespace EMPIRE {

class FEMesh;

/***********************************************************************************************
 * \brief This class performs mortar mapping
 * ***********/
class VertexMorphingMapper : public AbstractMapper {

private:
    /// Name of the mapper
    std::string name;

    /// Mesh
    FEMesh *mesh;

    /// Filter function name
    EMPIRE_VMM_FilterType filterType;

    /// Search radius;
    double filterRadius;

    /// nearest neighbors searching tree of FLAN libraray
    flann::Index<flann::L2<double> > *FLANNkd_tree;
    /// nodes constructing the searching tree
    flann::Matrix<double> *FLANNNodes;

    /// nearest neighbors searching tree of ANN libraray
    ANNkd_tree *nodesTree;
    /// nodes constructing the searching tree
    double **ANNNodes;

    /// directElemTable means the entries is not the node number, but the position in nodeCoors
    std::vector<int> **directElemTable;
    /// given a node, all the elements containing it are got
    std::vector<int> **nodeToElemTable;

    /// New sparse matrix.
    MathLibrary::SparseMatrix<double> *C_BB;
    /// New sparse matrix.
    MathLibrary::SparseMatrix<double> *C_BA;

    /// number of Gauss points used for computing triangle element mass matrix
    static const int numGPsMassMatrixTri;
    /// number of Gauss points used for computing quad element mass matrix
    static const int numGPsMassMatrixQuad;
    /// number of Gauss points used for computing shape function (tri) products on a clip
    static const int numGPsOnClipTri;
    /// number of Gauss points used for computing shape function (quad) products on a clip
    static const int numGPsOnClipQuad;

    /// pardiso variable
    void *pt[64]; // this is related to internal memory management, see PARDISO manual
    /// pardiso variable
    int iparm[64];
    /// pardiso variable
    int mtype;
    /// pardiso variable
    int maxfct;
    /// pardiso variable
    int mnum;
    /// pardiso variable
    int msglvl;
    /// pardiso variable
    int neq;
    /// pardiso variable
    int nrhs;

public:

    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name The name of the mapper
     * \param[in] _mesh Mesh (only FEM mesh for now)
     * \author Altug Emiroglu
     ***********/
    VertexMorphingMapper(std::string _name, EMPIRE::FEMesh *_mesh, EMPIRE_VMM_FilterType _filterType, double _filterRadius);
    /***********************************************************************************************
     * \brief Destructor
     * \author Altug Emiroglu
     ***********/
    virtual ~VertexMorphingMapper();

    /***********************************************************************************************
     * \brief Build Coupling Matrices
     * \author Altug Emiroglu
     ***********/
    void buildCouplingMatrices();

    /***********************************************************************************************
     * \brief Do consistent mapping on fields (e.g. displacements or tractions) --- C_BB * masterField = C_BA * slaveField
     * \param[in] slaveField the field of the slave side (e.g. x-displacements on all structure nodes)
     * \param[out] masterField the field of the master side (e.g. x-displacements on all fluid nodes)
     * \author Altug Emiroglu
     ***********/
    void consistentMapping(const double *slaveField, double *masterField) {}
    void conservativeMapping(const double *slaveField, double *masterField) {}
    void computeErrorsConsistentMapping(const double *slaveField, const double *masterField) {}

    /// defines number of threads used for MKL routines
    static int mklSetNumThreads;
    /// defines number of threads used for mapper routines
    static int mapperSetNumThreads;

private:

    /***********************************************************************************************
     * \brief Initialize all tables that help referring from an element to its nodes or vice versa
     * \author Altug Emiroglu
     ***********/
    void initTables();
    /***********************************************************************************************
     * \brief Deallocate the memory of all tables
     * \author Altug Emiroglu
     ***********/
    void deleteTables();
    /***********************************************************************************************
     * \brief Initialize the ANN nearest neighbor searching tree
     * \author Altug Emiroglu
     ***********/
    void initANNTree();
    /***********************************************************************************************
     * \brief Deallocate the memory of the searching tree
     * \author Altug Emiroglu
     ***********/
    void deleteANNTree();

    /***********************************************************************************************
     * \brief Find the overlapping candidates of the master element with the searching radius.
     * \param[in] controlNode the master node
     * \param[out] infNodeIdxs the influenced node indices
     * \param[out] infElemIdxs the influenced element indices
     * \author Altug Emiroglu
     ***********/
    void findCandidates(double* controlNode, std::set<int> *infNodeIdxs, std::set<int> *infElemIdxs);



};

} /* namespace EMPIRE */
#endif /* MORTARMAPPER_H_ */
