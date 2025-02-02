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
     * \brief Class Polygon
     *        given a set of counter clockwise ordered points of a polygon
     * 	      make a triangulation always starting from the first point
     **/
    class Polygon {
    private:

        /// Tolerance for adding new points to the polygon
        static double EPS_PolygonPoint;

    public:

        /// Counter-clockwise ordered polygon points
        std::vector<double*> points;

        /// Generated triangles from the given polygon
        std::vector<double*> triangles;

        /***********************************************************************************************
         * \brief Constructor
         * \author Altug Emiroglu
         ***********/
        Polygon(){}
        /***********************************************************************************************
         * \brief Destructor
         * \author Altug Emiroglu
         ***********/
        virtual ~Polygon();
        /***********************************************************************************************
         * \brief addPoint Adds the given point by copy
         * \param[in] _point The point to be added to the polygon
         * \author Altug Emiroglu
         ***********/
        void addPoint(double* _point);
        /***********************************************************************************************
         * \brief triangulate Triangulates the formed polygon.
     *        Each triangle starts from the first point
         * \author Altug Emiroglu
         ***********/
        void triangulate();
        /***********************************************************************************************
         * \brief printPolygon Prints the polygon on the terminal
         * \author Altug Emiroglu
         ***********/
        void printPolygon();
        
    };

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
        /// number of nodes
        int numNodes;
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
        FilterFunctionProduct(int _numNodes, double* _elem, double* _controlNode);
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
         * \param[in] _numGaussPoints number of Gauss points
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

    /// nearest neighbors searching tree of FLAN libraray
    flann::Index<flann::L2<double> > *FLANNkd_slaveTree;
    /// nodes constructing the searching tree
    flann::Matrix<double> *FLANNSlaveNodes;

    /// nearest neighbors searching tree of ANN libraray
    ANNkd_tree *slaveNodesTree;
    /// nodes constructing the searching tree
    double **ANNSlaveNodes;

    /// directElemTable means the entries is not the node number, but the position in nodeCoors
    std::vector<int> **slaveDirectElemTable;

    /// given a node, all the elements containing it are listed
    std::vector<int> **slaveNodeToElemTable;
    
    // The value of the default integration of the filter function (used for adjusting the filter function)
    double* masterFilterFunctionIntegrationOnSlave;

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
    int numGPsOnTri;
    /// number of Gauss points used for computing shape function (quad) products
    int numGPsOnQuad;

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
    VertexMorphingMapper(std::string _name, AbstractMesh *_meshA, AbstractMesh *_meshB);
    /***********************************************************************************************
     * \brief Destructor
     * \author Altug Emiroglu
     ***********/
    virtual ~VertexMorphingMapper();

    void setParameters(EMPIRE_VMM_FilterType _filterType = EMPIRE_VMM_HatFilter, double _filterRadius = 0.0);

    void initialize();

    /***********************************************************************************************
     * \brief Initialize and build Coupling Matrices
     * \author Altug Emiroglu
     ***********/
    void buildCouplingMatrices();

    /***********************************************************************************************
     * \brief Build matrix C_BA
     * \author Altug Emiroglu
     ***********/
    void build_C_BA();

    /***********************************************************************************************
     * \brief Find slave elements in the influence radius of master node
     * \param[in] _masterNodeIdx the index of the master node
     * \param[out] _slaveElemIdxs the indices of the influenced slave elements
     * \author Altug Emiroglu
     ***********/
    void findSlaveElements(const int _masterNodeIdx, std::vector<int>& _slaveElemIdxs);

    /***********************************************************************************************
     * \brief Find slave nodes in the influence radius of master node.
     * The first function extracts the coordinates and feeds to the second
     * \param[in] _masterNodeIdx the index of the master node
     * \param[out] _slaveNodeIdxs the indices of the influenced slave nodes
     * \author Altug Emiroglu
     ***********/
    void findSlaveNodes(const int _masterNodeIdx, std::vector<int>& _slaveNodeIdxs);
    void findSlaveNodes(double* _masterNodeCoords, std::vector<int>& _slaveNodeIdxs);
    
    
    /***********************************************************************************************
     * \brief Do consistent mapping on fields masterField = C_BA * slaveField
     * \param[in] slaveField the field of the slave side
     * \param[out] masterField the field of the master side
     * \author Altug Emiroglu
     ***********/
    void consistentMapping(const double *slaveField, double *masterField);
    /***********************************************************************************************
     * \brief Do conservative mapping on integral fields slaveField = C_BA^T * masterField
     * \param[in] masterField the field of the master side
     * \param[out] slaveField the field of the slave side
     * \author Altug Emiroglu
     * ***********/
    void conservativeMapping(const double *masterField, double *slaveField);

    /***********************************************************************************************
     * \brief Computes consistent mapping errors (Not implemented)
     * \author Altug Emiroglu
     ***********/
    void computeErrorsConsistentMapping(const double *slaveField, const double *masterField);

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
     * \brief Initialize the ANN nearest neighbor searching trees
     * \author Altug Emiroglu
     ***********/
    void initANNTree();
    /***********************************************************************************************
     * \brief Deallocate the memory of the searching trees
     * \author Altug Emiroglu
     ***********/
    void deleteANNTree();
    /***********************************************************************************************
     * \brief Find the influenced slave elements by the master node (C_BA)
     * \param[in] _masterNodeIdx the master node index
     * \author Altug Emiroglu
     ***********/
    // void findSlaveElemInfluencingMasterNodes(const int _masterNodeIdx);
    /***********************************************************************************************
     * \brief Performs clipping of a given element with the filter radius of the master node
     * \param[in] _masterNodeIdx the master node index
     * \param[in] _numNodes number of nodes
     * \param[in] _elem the global cartesian coordinates of the element vertices
     * \param[out] _polygon the resulting polygon
     * \author Altug Emiroglu
     ***********/
    void clipSlaveElementWithFilterRadius(const int _masterNodeIdx, const int _numNodes, double* _elem, VertexMorphingMapper::Polygon* _polygon);

    /***********************************************************************************************
     * \brief Performs integration of the filter and the shape function product on a full element
     * \param[in] _slaveElemIdx the influenced slave element index
     * \param[in] _masterNodeIdx the master node index
     * \author Altug Emiroglu
     ***********/
    void doFullIntegration(const int _slaveElemIdx, const int _masterNodeIdx);
    /***********************************************************************************************
     * \brief Performs integration of the filter and the shape function product on a clipped element
     * \param[in] _slaveElemIdx the influenced slave element index
     * \param[in] _masterNodeIdx
     * \author Altug Emiroglu
     ***********/
    void doIntegration_C_BA(const int _slaveElemIdx, const int _masterNodeIdx);

    /***********************************************************************************************
     * \brief Adjusts the filter function value by manipulating C_BA such that the unit integration property is satisfied
     * \author Altug Emiroglu
     ***********/
    void adjustFilterFunctions_BA();

    /***********************************************************************************************
     * \brief Finds clipping between P0 and Pn
     * \param[in] _masterNodeIdx the master node index
     * \param[in] P0 global cartesian coordinates of the first point in line
     * \param[in] Pn global cartesian coordinates of the second point in line
     * \param[in/out] _xsi line parameter
     * \return if the parameter is in acceptable range -EPS_XSI < xsi < 1+EPS_XSI
     * \author Altug Emiroglu
     ***********/
    bool findClippingWithMaster(const int _masterNodeIdx, double* _P0, double* _Pn, std::vector<double>& _xsi);

    /***********************************************************************************************
     * \brief Clamps line parameter xsi to the limits 0 <= xsi <= 1
     * \param[in/out] _xsi line parameter
     * \return if the parameter is in acceptable range -EPS_XSI < xsi < 1+EPS_XSI
     * \author Altug Emiroglu
     ***********/
    bool clampXsi(double &_xsi);
    /***********************************************************************************************
     * \brief Writes the polygons related to the master nodes
     * \author Altug Emiroglu
     ***********/
    void writeCartesianPolygons(std::string _fileName, std::map<int, std::vector<std::vector<double*> > >& _polygons);


private:

    /// Tolerance for clamping xsi to the bounds [0,1]
    static double EPS_XSI;

};



} /* namespace EMPIRE */
#endif /* VERTEXMORPHINGMAPPER_H_ */
