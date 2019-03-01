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
#include <string>
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"

#include "MappingFilter.h"
#include "DataField.h"
#include "FEMesh.h"
#include "ConnectionIO.h"
#include "ConnectionIOSetup.h"
#include "MapperAdapter.h"
#include "Message.h"
#include <math.h>
#include <iostream>

using namespace std;

namespace EMPIRE {
/********//**
 * \brief Test vertex morphing mapper. Also some test cases used to cause problems can be put here.
 ***********/
class TestVertexMorphingMapper: public CppUnit::TestFixture {
private:

public:
    void setUp() {

    }
    void tearDown() {
    }
    /***********************************************************************************************
     * \brief compareC_BBandC_BA
     ***********/
    void compareC_BBandC_BA() {

    }
    /***********************************************************************************************
     * \brief rowsum of C_BA_dual is negative
     ***********/
    void testComputeMappingOnTrias() {
      
     /* 
          Mesh A and Mesh B are the same
        1___2___6
        |4 /|\ 3|
        | /8|7\ |
       3|/__4__\|8
        |\ 5|6 /|
        |1\ | /2|
        |__\|/__|
        5   7   9
     */
      
        FEMesh *meshA;
        FEMesh *meshB;
        DataField *a1;
        DataField *b1;
        static const double EPS = 1E-8;
        int numNodesPerElem = 3;
        int numNodesA = 9;
        int numElemsA = 8;
        meshA = new FEMesh("", numNodesA, numElemsA, false);
        for (int i = 0; i < numElemsA; i++)
            meshA->numNodesPerElem[i] = numNodesPerElem;
        meshA->initElems();

        for (int i = 0; i < numNodesA; i++)
            meshA->nodeIDs[i] = i + 1;

        meshA->nodes[0 * 3 + 0] = 0.0;
        meshA->nodes[0 * 3 + 1] = 1.0;
        meshA->nodes[0 * 3 + 2] = 0.0;

        meshA->nodes[1 * 3 + 0] = 0.5;
        meshA->nodes[1 * 3 + 1] = 1.0;
        meshA->nodes[1 * 3 + 2] = 0.0;

        meshA->nodes[2 * 3 + 0] = 0.0;
        meshA->nodes[2 * 3 + 1] = 0.5;
        meshA->nodes[2 * 3 + 2] = 0.0;

        meshA->nodes[3 * 3 + 0] = 0.5;
        meshA->nodes[3 * 3 + 1] = 0.5;
        meshA->nodes[3 * 3 + 2] = 0.0;
	
	meshA->nodes[4 * 3 + 0] = 0.0;
        meshA->nodes[4 * 3 + 1] = 0.0;
        meshA->nodes[4 * 3 + 2] = 0.0;
	
	meshA->nodes[5 * 3 + 0] = 1.0;
        meshA->nodes[5 * 3 + 1] = 1.0;
        meshA->nodes[5 * 3 + 2] = 0.0;
	
	meshA->nodes[6 * 3 + 0] = 0.5;
        meshA->nodes[6 * 3 + 1] = 0.0;
        meshA->nodes[6 * 3 + 2] = 0.0;
	
	meshA->nodes[7 * 3 + 0] = 1.0;
        meshA->nodes[7 * 3 + 1] = 0.5;
        meshA->nodes[7 * 3 + 2] = 0.0;
        
        meshA->nodes[8 * 3 + 0] = 1.0;
        meshA->nodes[8 * 3 + 1] = 0.0;
        meshA->nodes[8 * 3 + 2] = 0.0;

        meshA->elems[0] = 5;
        meshA->elems[1] = 7;
        meshA->elems[2] = 3;
	
        meshA->elems[3] = 9;
	meshA->elems[4] = 8;
	meshA->elems[5] = 7;
	
	meshA->elems[6] = 6;
	meshA->elems[7] = 2;
	meshA->elems[8] = 8;
	
	meshA->elems[9] = 1;
	meshA->elems[10] = 3;
	meshA->elems[11] = 2;
	
	meshA->elems[12] = 3;
	meshA->elems[13] = 7;
	meshA->elems[14] = 4;
	
	meshA->elems[15] = 4;
	meshA->elems[16] = 7;
	meshA->elems[17] = 8;
	
	meshA->elems[18] = 4;
	meshA->elems[19] = 8;
	meshA->elems[20] = 2;
	
	meshA->elems[21] = 4;
	meshA->elems[22] = 2;
	meshA->elems[23] = 3;
	
        int numNodesB = 9;
        int numElemsB = 8;
        meshB = new FEMesh("", numNodesB, numElemsB, false);
        for (int i = 0; i < numElemsB; i++)
            meshB->numNodesPerElem[i] = numNodesPerElem;
        meshB->initElems();

        for (int i = 0; i < numNodesB; i++)
            meshB->nodeIDs[i] = i + 1;

        meshB->nodes[0 * 3 + 0] = 0.0;
        meshB->nodes[0 * 3 + 1] = 1.0;
        meshB->nodes[0 * 3 + 2] = 0.0;

        meshB->nodes[1 * 3 + 0] = 0.5;
        meshB->nodes[1 * 3 + 1] = 1.0;
        meshB->nodes[1 * 3 + 2] = 0.0;

        meshB->nodes[2 * 3 + 0] = 0.0;
        meshB->nodes[2 * 3 + 1] = 0.5;
        meshB->nodes[2 * 3 + 2] = 0.0;

        meshB->nodes[3 * 3 + 0] = 0.5;
        meshB->nodes[3 * 3 + 1] = 0.5;
        meshB->nodes[3 * 3 + 2] = 0.0;
	
	meshB->nodes[4 * 3 + 0] = 0.0;
        meshB->nodes[4 * 3 + 1] = 0.0;
        meshB->nodes[4 * 3 + 2] = 0.0;
	
	meshB->nodes[5 * 3 + 0] = 1.0;
        meshB->nodes[5 * 3 + 1] = 1.0;
        meshB->nodes[5 * 3 + 2] = 0.0;
	
	meshB->nodes[6 * 3 + 0] = 0.5;
        meshB->nodes[6 * 3 + 1] = 0.0;
        meshB->nodes[6 * 3 + 2] = 0.0;
	
	meshB->nodes[7 * 3 + 0] = 1.0;
        meshB->nodes[7 * 3 + 1] = 0.5;
        meshB->nodes[7 * 3 + 2] = 0.0;
        
        meshB->nodes[8 * 3 + 0] = 1.0;
        meshB->nodes[8 * 3 + 1] = 0.0;
        meshB->nodes[8 * 3 + 2] = 0.0;

        meshB->elems[0] = 5;
        meshB->elems[1] = 7;
        meshB->elems[2] = 3;
	
        meshB->elems[3] = 9;
	meshB->elems[4] = 8;
	meshB->elems[5] = 7;
	
	meshB->elems[6] = 6;
	meshB->elems[7] = 2;
	meshB->elems[8] = 8;
	
	meshB->elems[9] = 1;
	meshB->elems[10] = 3;
	meshB->elems[11] = 2;
	
	meshB->elems[12] = 3;
	meshB->elems[13] = 7;
	meshB->elems[14] = 4;
	
	meshB->elems[15] = 4;
	meshB->elems[16] = 7;
	meshB->elems[17] = 8;
	
	meshB->elems[18] = 4;
	meshB->elems[19] = 8;
	meshB->elems[20] = 2;
	
	meshB->elems[21] = 4;
	meshB->elems[22] = 2;
	meshB->elems[23] = 3;
	
	EMPIRE_VMM_FilterType filterType = EMPIRE_VMM_GaussianFilter;
	double filterRadius = 0.5;
	bool mortar = false;
	
	MapperAdapter *mapper = new MapperAdapter("testVM", meshA, meshB);
        mapper->initVertexMorphingMapper(filterType, filterRadius, mortar);
	
	a1 = new DataField("a1", EMPIRE_DataField_atNode, meshA->numNodes, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);
        b1 = new DataField("b1", EMPIRE_DataField_atNode, meshB->numNodes, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);
        for (int i = 0; i < meshA->numNodes; i++)
            a1->data[i] = 1.0;
	
// 	a1->data[3]=1.0;

	
	AbstractFilter *filterConsistent = new MappingFilter(mapper);
        ConnectionIOSetup::setupIOForFilter(filterConsistent, meshA, a1, meshB, b1);
        filterConsistent->filtering();	
	
// 	cout << "Results" << endl;
// 	for (int i = 0; i < meshB->numNodes; i++)
// 	  cout << b1->data[i] << endl;
	
        delete mapper;
        delete meshA;
        delete meshB;
	delete a1;
	delete b1;
	delete filterConsistent;
	
// 	exit(0);

      
    }
    
    void testComputeMappingOnQuads() {
      
     /* 
          Mesh A and Mesh B are the same
        
        1___2___6
        |   |   |
        | 1 | 2 |
       3|___4__ |8
        |   |   |
        | 3 | 4 |
        |__ | __|
        5   7   9
        
     */
      
        FEMesh *meshA;
        FEMesh *meshB;
        DataField *a1;
        DataField *b1;
        static const double EPS = 1E-8;
        int numNodesPerElem = 4;
        int numNodesA = 9;
        int numElemsA = 4;
        meshA = new FEMesh("", numNodesA, numElemsA, false);
        for (int i = 0; i < numElemsA; i++)
            meshA->numNodesPerElem[i] = numNodesPerElem;
        meshA->initElems();

        for (int i = 0; i < numNodesA; i++)
            meshA->nodeIDs[i] = i + 1;

        meshA->nodes[0 * 3 + 0] = 0.0;
        meshA->nodes[0 * 3 + 1] = 1.0;
        meshA->nodes[0 * 3 + 2] = 0.0;

        meshA->nodes[1 * 3 + 0] = 0.5;
        meshA->nodes[1 * 3 + 1] = 1.0;
        meshA->nodes[1 * 3 + 2] = 0.0;

        meshA->nodes[2 * 3 + 0] = 0.0;
        meshA->nodes[2 * 3 + 1] = 0.5;
        meshA->nodes[2 * 3 + 2] = 0.0;

        meshA->nodes[3 * 3 + 0] = 0.5;
        meshA->nodes[3 * 3 + 1] = 0.5;
        meshA->nodes[3 * 3 + 2] = 0.0;
	
	meshA->nodes[4 * 3 + 0] = 0.0;
        meshA->nodes[4 * 3 + 1] = 0.0;
        meshA->nodes[4 * 3 + 2] = 0.0;
	
	meshA->nodes[5 * 3 + 0] = 1.0;
        meshA->nodes[5 * 3 + 1] = 1.0;
        meshA->nodes[5 * 3 + 2] = 0.0;
	
	meshA->nodes[6 * 3 + 0] = 0.5;
        meshA->nodes[6 * 3 + 1] = 0.0;
        meshA->nodes[6 * 3 + 2] = 0.0;
	
	meshA->nodes[7 * 3 + 0] = 1.0;
        meshA->nodes[7 * 3 + 1] = 0.5;
        meshA->nodes[7 * 3 + 2] = 0.0;
        
        meshA->nodes[8 * 3 + 0] = 1.0;
        meshA->nodes[8 * 3 + 1] = 0.0;
        meshA->nodes[8 * 3 + 2] = 0.0;

        meshA->elems[0] = 1;
        meshA->elems[1] = 3;
        meshA->elems[2] = 4;
	meshA->elems[3] = 2;
	
	meshA->elems[4] = 2;
	meshA->elems[5] = 4;	
	meshA->elems[6] = 8;
	meshA->elems[7] = 6;
	
	meshA->elems[8] = 3;	
	meshA->elems[9] = 5;
	meshA->elems[10] = 7;
	meshA->elems[11] = 4;
	
	meshA->elems[12] = 4;
	meshA->elems[13] = 7;
	meshA->elems[14] = 9;	
	meshA->elems[15] = 8;
	
        int numNodesB = 9;
        int numElemsB = 4;
        meshB = new FEMesh("", numNodesB, numElemsB, false);
        for (int i = 0; i < numElemsB; i++)
            meshB->numNodesPerElem[i] = numNodesPerElem;
        meshB->initElems();

        for (int i = 0; i < numNodesB; i++)
            meshB->nodeIDs[i] = i + 1;

        meshB->nodes[0 * 3 + 0] = 0.0;
        meshB->nodes[0 * 3 + 1] = 1.0;
        meshB->nodes[0 * 3 + 2] = 0.0;
            
        meshB->nodes[1 * 3 + 0] = 0.5;
        meshB->nodes[1 * 3 + 1] = 1.0;
        meshB->nodes[1 * 3 + 2] = 0.0;
            
        meshB->nodes[2 * 3 + 0] = 0.0;
        meshB->nodes[2 * 3 + 1] = 0.5;
        meshB->nodes[2 * 3 + 2] = 0.0;
            
        meshB->nodes[3 * 3 + 0] = 0.5;
        meshB->nodes[3 * 3 + 1] = 0.5;
        meshB->nodes[3 * 3 + 2] = 0.0;
	    
	meshB->nodes[4 * 3 + 0] = 0.0;
        meshB->nodes[4 * 3 + 1] = 0.0;
        meshB->nodes[4 * 3 + 2] = 0.0;
	    
	meshB->nodes[5 * 3 + 0] = 1.0;
        meshB->nodes[5 * 3 + 1] = 1.0;
        meshB->nodes[5 * 3 + 2] = 0.0;
	    
	meshB->nodes[6 * 3 + 0] = 0.5;
        meshB->nodes[6 * 3 + 1] = 0.0;
        meshB->nodes[6 * 3 + 2] = 0.0;
	    
	meshB->nodes[7 * 3 + 0] = 1.0;
        meshB->nodes[7 * 3 + 1] = 0.5;
        meshB->nodes[7 * 3 + 2] = 0.0;
            
        meshB->nodes[8 * 3 + 0] = 1.0;
        meshB->nodes[8 * 3 + 1] = 0.0;
        meshB->nodes[8 * 3 + 2] = 0.0;
            
        meshB->elems[0] = 1;
        meshB->elems[1] = 3;
        meshB->elems[2] = 4;
	meshB->elems[3] = 2;
	    
	meshB->elems[4] = 2;
	meshB->elems[5] = 4;
	meshB->elems[6] = 8;
	meshB->elems[7] = 6;
	    
	meshB->elems[8] = 3;
	meshB->elems[9] = 5;
	meshB->elems[10] = 7;
	meshB->elems[11] = 4;
	    
	meshB->elems[12] = 4;
	meshB->elems[13] = 7;
	meshB->elems[14] = 9;
	meshB->elems[15] = 8;
	
	EMPIRE_VMM_FilterType filterType = EMPIRE_VMM_HatFilter;
	double filterRadius = 0.51;
	bool mortar = false;
	
	MapperAdapter *mapper = new MapperAdapter("testVM", meshA, meshB);
        mapper->initVertexMorphingMapper(filterType, filterRadius, mortar);
	
	a1 = new DataField("a1", EMPIRE_DataField_atNode, meshA->numNodes, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);
        b1 = new DataField("b1", EMPIRE_DataField_atNode, meshB->numNodes, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);
	for (int i = 0; i < meshA->numNodes; i++)
	  a1->data[i] = 1.0;
	
// 	a1->data[3]=1.0;
	
	AbstractFilter *filterConsistent = new MappingFilter(mapper);
        ConnectionIOSetup::setupIOForFilter(filterConsistent, meshA, a1, meshB, b1);
        filterConsistent->filtering();	
	
// 	for (int i = 0; i < meshB->numNodes; i++)
// 	  cout << b1->data[i] << endl;
	
// 	exit(0);
	
        delete mapper;
        delete meshA;
        delete meshB;
	delete a1;
	delete b1;
	delete filterConsistent;
      
    }
    
    CPPUNIT_TEST_SUITE( TestVertexMorphingMapper );
    CPPUNIT_TEST( testComputeMappingOnTrias);
    CPPUNIT_TEST( testComputeMappingOnQuads);
        //CPPUNIT_TEST( testMemoryLeakOfConstructor); // test memory leak, comment it except when checking memory leak
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestVertexMorphingMapper);
