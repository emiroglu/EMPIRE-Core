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

namespace EMPIRE {

class FEMesh;
/***********************************************************************************************
 * \brief This class performs mortar mapping
 * ***********/
class VertexMorphingMapper : public AbstractMapper {
private:

    /**
     * \brief Class AbstractFilterFunction for the derived filter functions.
     *        The important is the overloading of FilterFunction::computeFunction in the derived classes
     **/
    class AbstractFilterFunction {
    public:

        double filterRadius;

    public:

        /***********************************************************************************************
         * \brief Constructor
         * \param[in] _filterRadius The filter radius
         * \author Altug Emiroglu
         ***********/
        AbstractFilterFunction(double _filterRadius): filterRadius(_filterRadius) {}
        /***********************************************************************************************
         * \brief Destructor
         * \author Altug Emiroglu
         ***********/
        virtual ~AbstractFilterFunction() {}
        /***********************************************************************************************
         * \brief Computes the function value given the support center and the global cartesian coordinates
         *        Only this function is overloaded in the derived classes
         * \param[in] _supportCenter the support center of the filter function
         * \param[in] _globalCoor the global coordinates to compute the filter function value
         * \return the function value
         * \author Altug Emiroglu
         ***********/
        virtual double computeFunction(double* _supportCenter, double* _globalCoor) = 0;

    };

    /**
     * \brief Class HatFilterFunction
     *        Computes the filter function value using a hat filter function
     **/
    class HatFilterFunction: public AbstractFilterFunction {

    public:

        /***********************************************************************************************
         * \brief Constructor
         * \param[in] _filterRadius The filter radius
         * \author Altug Emiroglu
         ***********/
        HatFilterFunction(double _filterRadius): AbstractFilterFunction(_filterRadius) {}
        /***********************************************************************************************
         * \brief Destructor
         * \author Altug Emiroglu
         ***********/
        virtual ~HatFilterFunction() {}
        /***********************************************************************************************
         * \brief Computes the function value given the support center and the global cartesian coordinates
         * \param[in] _supportCenter the support center of the filter function
         * \param[in] _globalCoor the global coordinates to compute the filter function value
         * \return the function value
         * \author Altug Emiroglu
         ***********/
        double computeFunction(double* _supportCenter, double* _globalCoor);

    };

    /**
     * \brief Class GaussianFilterFunction
     *        Computes the filter function value using a Gaussian filter function
     **/
    class GaussianFilterFunction: public AbstractFilterFunction {

    public:

        /***********************************************************************************************
         * \brief Constructor
         * \param[in] _filterRadius The filter radius
         * \author Altug Emiroglu
         ***********/
        GaussianFilterFunction(double _filterRadius): AbstractFilterFunction(_filterRadius) {}
        /***********************************************************************************************
         * \brief Destructor
         * \author Altug Emiroglu
         ***********/
        virtual ~GaussianFilterFunction(){}
        /***********************************************************************************************
         * \brief Computes the function value given the support center and the global cartesian coordinates
         * \param[in] _supportCenter the support center of the filter function
         * \param[in] _globalCoor the global coordinates to compute the filter function value
         * \return the function value
         * \author Altug Emiroglu
         ***********/
        double computeFunction(double* _supportCenter, double* _globalCoor);

    };

    /**
     * \brief Class FilterFunctionProduct computes the filter and shape function products on an element
     *        This class is mainly a modification of class MortarMapper::ShapeFunctionProduct
     **/
    // Integrand for the computation of the matrices
    class FilterFunctionProduct: public EMPIRE::MathLibrary::IntegrandFunction {
    private:

        /// control node global coordinates
        double* controlNode;
        /// element nodes
        double* elem;
        /// number of Gauss points
        int numGaussPoints;
        /// coordinates of all gauss points
        double * gaussPoints;
        /// shape function ID of the slave element
        int funcID;
        /// different shape function products on all Gauss points
        double **functionProducts;

    public:

        /***********************************************************************************************
         * \brief Constructor This class is derived from
         * \param[in] _elem Triangular element to perform integration
         * \param[in] _controlNode The control node global cartesian coordinates
         * \author Altug Emiroglu
         ***********/
        FilterFunctionProduct(double* _elem, double* _controlNode);
        /***********************************************************************************************
         * \brief Destructor
         * \author Altug Emiroglu
         ***********/
        virtual ~FilterFunctionProduct();
        /***********************************************************************************************
         * \brief Computes the product of the filter function and the shape functions and stores them in functionProducts
         * \param[in] _filterFunction The filter function object to compute the products with
         * \author Altug Emiroglu
         ***********/
        void computeFunctionProducts(AbstractFilterFunction* _filterFunction);
        /***********************************************************************************************
         * \brief set all Gauss points
         * \param[in] _gaussPoints x,y,z coordinates of all Gauss points
         * \param[in] _numGaussPoints x,y,z coordinates of all Gauss points
         * \author Altug Emiroglu
         ***********/
        void setGaussPoints(const double* _gaussPoints, int _numGaussPoints);
        /***********************************************************************************************
         * \brief Set the function ID
         * \param[in] _funcID function ID related to the respective node of the element
         * \author Altug Emiroglu
         ***********/
        void setFunctionID(int _funcID);
        /***********************************************************************************************
         * \brief Retrieve the filter and shape function product on the Gauss point. In fact, the shape function
         *        products have been computed by computeFunctionProducts().
         * \param[in] gaussPoint x,y,z coordinates of the Gauss point
         * \return the product value
         * \author Altug Emiroglu
         ***********/
        double operator()(double* gaussPoint);

    };

private:
    /// Name of the mapper
    std::string name;

    /// Mesh
    FEMesh *meshA;
    FEMesh *meshB;

    /// Filter function pointer
    AbstractFilterFunction* filterFunction;

    /// Search radius;
    double filterRadius;

    /// mortar type of computation through C_BB calculation
    bool mortar;

    /// nearest neighbors searching tree of FLAN libraray
    flann::Index<flann::L2<double> > *FLANNkd_tree;
    /// nodes constructing the searching tree
    flann::Matrix<double> *FLANNSlaveNodes;

    /// nearest neighbors searching tree of ANN libraray
    ANNkd_tree *slaveNodesTree;
    /// nodes constructing the searching tree
    double **ANNNodes;

    /// directElemTable means the entries is not the node number, but the position in nodeCoors
    std::vector<int> **masterDirectElemTable;
    std::vector<int> **slaveDirectElemTable;

    /// given a node, all the elements containing it are got
    std::vector<int> **masterNodeToElemTable;
    std::vector<int> **slaveNodeToElemTable;

    /// New sparse matrix.
    MathLibrary::SparseMatrix<double> *C_BB;
    /// New sparse matrix.
    MathLibrary::SparseMatrix<double> *C_BA;

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

    /// number of Gauss points used for computing shape function (tri) products on a clip
    static const int numGPsOnTri;
    /// number of Gauss points used for computing triangle element mass matrix
    static const int numGPsMassMatrixTri;

    friend class TestVertexMorphingMapper;

public:

    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name The name of the mapper
     * \param[in] _meshA Mesh (only FEM mesh for now)
     * \param[in] _meshB Mesh (only FEM mesh for now)
     * \param[in] _filterType choose from one of the EMPIRE_VMM_FilterType types
     * \param[in] _filterRadius effective radius of the filter function
     * \param[in] _mortar through C_BB computation
     * \author Altug Emiroglu
     ***********/
    VertexMorphingMapper(std::string _name, AbstractMesh *_meshA, AbstractMesh *_meshB, EMPIRE_VMM_FilterType _filterType, double _filterRadius, bool _mortar=false);
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
    void consistentMapping(const double *slaveField, double *masterField);
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
     * \param[in] _controlNode the control field node
     * \param[out] _infNodeIdxs the influenced node indices
     * \param[out] _infElemIdxs the influenced element indices
     * \author Altug Emiroglu
     ***********/
    void findCandidates(double* _controlNode, std::set<int>* _infNodeIdxs, std::set<int>* _infElemIdxs);
    /***********************************************************************************************
     * \brief Performs integration of the filter and the shape function product on a full element
     * \param[in] _controlNode the control field node
     * \param[in] _elemIdx the influenced element index
     * \param[out] _contributions the result of the integration
     * \author Altug Emiroglu
     ***********/
    void doFullIntegration(double* _controlNode, int _elemIdx, double* _contributions, double& _intFilter);
    /***********************************************************************************************
     * \brief Performs integration of the filter and the shape function product on a clipped element
     * \param[in] _controlNode the control field node
     * \param[in] _elemIdx the influenced element index
     * \param[in] _elemNodeInside the list of if the nodes are inside the filter radius
     * \param[out] _contributions the result of the integration
     * \author Altug Emiroglu
     ***********/
    void doPartialIntegration(double* _controlNode, int _elemIdx, std::vector<bool>* _elemNodeInside, double* _contributions, double& _intFilter);
    /***********************************************************************************************
     * \brief Assembles the contributions to C_BA
     * \param[in] _controlNodeIdx the control field node index
     * \param[in] _elemIdx the influenced element index
     * \param[in] _contributions the result of the integration
     * \author Altug Emiroglu
     ***********/
    void assemble_C_BA(int _controlNodeIdx, double* _contributions, double intFilter);
    /***********************************************************************************************
     * \brief Clamps line parameter xsi to the limits 0 <= xsi <= 1
     * \param[in/out] _xsi line parameter
     * \author Altug Emiroglu
     ***********/
    void clampXsi(double &_xsi);


private:

    /// Tolerance for cleaning a triangle before integrating
    static double EPS_XSI;

};



} /* namespace EMPIRE */
#endif /* VERTEXMORPHINGMAPPER_H_ */
