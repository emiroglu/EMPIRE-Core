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

    FEMesh* meshA;
    FEMesh* meshAcopy;
    FEMesh* meshB;
    FEMesh* meshBcopy;

public:
    void setUp() {

        /*

         Mesh A                  Mesh B

        1___2___6               1___2___6
        |4 /|\ 3|               |   |   |
        | /8|7\ |               | 1 | 2 |
       3|/__4__\|8             3|___4__ |8
        |\ 5|6 /|               |   |   |
        |1\ | /2|               | 3 | 4 |
        |__\|/__|               |__ | __|
        5   7   9               5   7   9

      */

        // Mesh A
        {
            int numNodesA = 9;
            int numElemsA = 8;
            int numNodesPerElemA = 3;
            meshA = new FEMesh("meshA", numNodesA, numElemsA, false);
            for (int i = 0; i < numElemsA; i++)
                meshA->numNodesPerElem[i] = numNodesPerElemA;
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
        }

        // Mesh A copy
        {
            meshAcopy = new FEMesh("meshAcopy", meshA->numNodes, meshA->numElems, false);
            for (int i = 0; i < meshA->numElems; i++)
                meshAcopy->numNodesPerElem[i] = meshA->numNodesPerElem[i];
            meshAcopy->initElems();

            for (int i = 0; i < meshA->numNodes; i++)
                meshAcopy->nodeIDs[i] = meshA->nodeIDs[i];
            for (int i = 0; i < meshA->numNodes * 3; i++)
                meshAcopy->nodes[i] = meshA->nodes[i];
            for (int i = 0; i < meshA->numElems * 3; i++)
                meshAcopy->elems[i] = meshA->elems[i];
        }

        {
            int numNodesB = 9;
            int numElemsB = 4;
            int numNodesPerElemB = 4;
            meshB = new FEMesh("meshB", numNodesB, numElemsB, false);
            for (int i = 0; i < numElemsB; i++)
                meshB->numNodesPerElem[i] = numNodesPerElemB;
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
        }

        // Mesh B copy
        {
            meshBcopy = new FEMesh("meshBcopy", meshB->numNodes, meshB->numElems, false);
            for (int i = 0; i < meshB->numElems; i++)
                meshBcopy->numNodesPerElem[i] = meshB->numNodesPerElem[i];
            meshBcopy->initElems();

            for (int i = 0; i < meshB->numNodes; i++)
                meshBcopy->nodeIDs[i] = meshB->nodeIDs[i];
            for (int i = 0; i < meshB->numNodes * 3; i++)
                meshBcopy->nodes[i] = meshB->nodes[i];
            for (int i = 0; i < meshB->numElems * 4; i++)
                meshBcopy->elems[i] = meshB->elems[i];
        }

    }
    void tearDown() {

        delete meshA;
        delete meshAcopy;
        delete meshB;
        delete meshBcopy;

    }
    /***********************************************************************************************
     * \brief compareC_BBandC_BA
     ***********/
    void compareC_BBandC_BA() {

    }
    /***********************************************************************************************
     * \brief meshA to meshA
     ***********/
    void testHatTriasToTrias() {

        double filterRadius = 0.45;
        MapperAdapter* mapper = new MapperAdapter("testTriasToTrias", meshA, meshAcopy);
        mapper->initVertexMorphingMapper(EMPIRE_VMM_HatFilter, filterRadius);

        DataField *a1;
        DataField *a2;
        a1 = new DataField("a1", EMPIRE_DataField_atNode, meshA->numNodes, EMPIRE_DataField_scalar,
                           EMPIRE_DataField_field);
        a1->data[3] = 1.0;

        a2 = new DataField("b1", EMPIRE_DataField_atNode, meshAcopy->numNodes, EMPIRE_DataField_scalar,
                           EMPIRE_DataField_field);

        AbstractFilter *filterConsistent = new MappingFilter(mapper);
        ConnectionIOSetup::setupIOForFilter(filterConsistent, meshA, a1, meshAcopy, a2);
        filterConsistent->filtering();

        delete mapper;
        delete a1;
        delete a2;
        delete filterConsistent;

    }
    
    /***********************************************************************************************
     * \brief meshB to meshB
     ***********/
    void testHatQuadsToQuads() {

        double filterRadius = 0.45;
        MapperAdapter* mapper = new MapperAdapter("testTriasToTrias", meshB, meshBcopy);
        mapper->initVertexMorphingMapper(EMPIRE_VMM_HatFilter, filterRadius);

        DataField *b1;
        DataField *b2;
        b1 = new DataField("a1", EMPIRE_DataField_atNode, meshB->numNodes, EMPIRE_DataField_scalar,
                           EMPIRE_DataField_field);
        b1->data[3] = 1.0;

        b2 = new DataField("b1", EMPIRE_DataField_atNode, meshBcopy->numNodes, EMPIRE_DataField_scalar,
                           EMPIRE_DataField_field);

        AbstractFilter *filterConsistent = new MappingFilter(mapper);
        ConnectionIOSetup::setupIOForFilter(filterConsistent, meshB, b1, meshBcopy, b1);
        filterConsistent->filtering();

        delete mapper;
        delete b1;
        delete b2;
        delete filterConsistent;

    }
    
    /***********************************************************************************************
     * \brief meshA to meshB
     ***********/
    void testHatTriasToQuads() {

        double filterRadius = 0.45;
        MapperAdapter *mapper = new MapperAdapter("testTriasToQuads", meshA, meshB);
        mapper->initVertexMorphingMapper(EMPIRE_VMM_HatFilter, filterRadius);

        DataField *a1;
        DataField *b1;
        a1 = new DataField("a1", EMPIRE_DataField_atNode, meshA->numNodes, EMPIRE_DataField_scalar,
                           EMPIRE_DataField_field);
        a1->data[3] = 1.0;
        b1 = new DataField("b1", EMPIRE_DataField_atNode, meshB->numNodes, EMPIRE_DataField_scalar,
                           EMPIRE_DataField_field);

        AbstractFilter *filterConsistent = new MappingFilter(mapper);
        ConnectionIOSetup::setupIOForFilter(filterConsistent, meshA, a1, meshB, b1);
        filterConsistent->filtering();

        delete mapper;
        delete a1;
        delete b1;
        delete filterConsistent;

    }
    
    /***********************************************************************************************
     * \brief meshB to meshA
     ***********/
    void testHatQuadsToTrias() {

        double filterRadius = 0.45;
        MapperAdapter *mapper = new MapperAdapter("testTriasToQuads", meshB, meshA);
        mapper->initVertexMorphingMapper(EMPIRE_VMM_HatFilter, filterRadius);

        DataField *a1;
        DataField *b1;
        a1 = new DataField("a1", EMPIRE_DataField_atNode, meshA->numNodes, EMPIRE_DataField_scalar,
                           EMPIRE_DataField_field);

        b1 = new DataField("b1", EMPIRE_DataField_atNode, meshB->numNodes, EMPIRE_DataField_scalar,
                           EMPIRE_DataField_field);
        b1->data[3] = 1.0;

        AbstractFilter *filterConsistent = new MappingFilter(mapper);
        ConnectionIOSetup::setupIOForFilter(filterConsistent, meshB, b1, meshA, a1);
        filterConsistent->filtering();

        delete mapper;
        delete a1;
        delete b1;
        delete filterConsistent;

    }
    
    /***********************************************************************************************
     * \brief check consistency
     ***********/
    void testConsistency() {

        // Initiate fields
        DataField *a1;
        DataField *a1copy;
        DataField *b1;
        DataField *b1copy;
        a1 = new DataField("a1", EMPIRE_DataField_atNode, meshA->numNodes, EMPIRE_DataField_scalar,
                           EMPIRE_DataField_field);
        a1copy = new DataField("a1copy", EMPIRE_DataField_atNode, meshAcopy->numNodes, EMPIRE_DataField_scalar,
                               EMPIRE_DataField_field);
        b1 = new DataField("b1", EMPIRE_DataField_atNode, meshB->numNodes, EMPIRE_DataField_scalar,
                           EMPIRE_DataField_field);
        b1copy = new DataField("b1copy", EMPIRE_DataField_atNode, meshBcopy->numNodes, EMPIRE_DataField_scalar,
                               EMPIRE_DataField_field);

        // Mapper properties
        double filterRadius = 0.45;

        // Mappers
        { // Hat A - A
            MapperAdapter *mapper = new MapperAdapter("", meshA, meshAcopy);
            mapper->initVertexMorphingMapper(EMPIRE_VMM_HatFilter, filterRadius);

            for (int i = 0; i < a1->numLocations; i++)
                a1->data[i] = 1.0;
            for (int i = 0; i < a1copy->numLocations; i++)
                a1copy->data[i] = 0.0;

            AbstractFilter *filterConsistent = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConsistent, meshA, a1, meshAcopy, a1copy);
            filterConsistent->filtering();

            const double EPS = 1e-10;
            for (int i = 0; i < a1copy->numLocations; i++)
                CPPUNIT_ASSERT(fabs(a1copy->data[i] - 1.0) < EPS);

            delete mapper;
            delete filterConsistent;
        }


        { // Gaussian A - A
            MapperAdapter *mapper = new MapperAdapter("", meshA, meshAcopy);
            mapper->initVertexMorphingMapper(EMPIRE_VMM_GaussianFilter, filterRadius);

            for (int i = 0; i < a1->numLocations; i++)
                a1->data[i] = 1.0;
            for (int i = 0; i < a1copy->numLocations; i++)
                a1copy->data[i] = 0.0;

            AbstractFilter *filterConsistent = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConsistent, meshA, a1, meshAcopy, a1copy);
            filterConsistent->filtering();

            const double EPS = 1e-10;
            for (int i = 0; i < a1copy->numLocations; i++)
                CPPUNIT_ASSERT(fabs(a1copy->data[i] - 1.0) < EPS);

            delete mapper;
            delete filterConsistent;

        }

        { // Hat B - B
            MapperAdapter *mapper = new MapperAdapter("", meshB, meshBcopy);
            mapper->initVertexMorphingMapper(EMPIRE_VMM_HatFilter, filterRadius);

            for (int i = 0; i < b1->numLocations; i++)
                b1->data[i] = 1.0;
            for (int i = 0; i < b1copy->numLocations; i++)
                b1copy->data[i] = 0.0;

            AbstractFilter *filterConsistent = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConsistent, meshB, b1, meshBcopy, b1copy);
            filterConsistent->filtering();

            const double EPS = 1e-10;
            for (int i = 0; i < b1copy->numLocations; i++)
                CPPUNIT_ASSERT(fabs(b1copy->data[i] - 1.0) < EPS);

            delete mapper;
            delete filterConsistent;

        }

        { // Gaussian B - B
            MapperAdapter *mapper = new MapperAdapter("", meshB, meshBcopy);
            mapper->initVertexMorphingMapper(EMPIRE_VMM_GaussianFilter, filterRadius);

            for (int i = 0; i < b1->numLocations; i++)
                b1->data[i] = 1.0;
            for (int i = 0; i < b1copy->numLocations; i++)
                b1copy->data[i] = 0.0;

            AbstractFilter *filterConsistent = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConsistent, meshB, b1, meshBcopy, b1copy);
            filterConsistent->filtering();

            const double EPS = 1e-10;
            for (int i = 0; i < b1copy->numLocations; i++)
                CPPUNIT_ASSERT(fabs(b1copy->data[i] - 1.0) < EPS);

            delete mapper;
            delete filterConsistent;

        }

        delete a1;
        delete a1copy;
        delete b1;
        delete b1copy;

    }
    
    void testMemoryLeakOfConstructor(){

        double filterRadius = 0.45;
        for (int i = 0; i < 100000; i++){

            MapperAdapter *mapper = new MapperAdapter("", meshA, meshAcopy);
            mapper->initVertexMorphingMapper(EMPIRE_VMM_GaussianFilter, filterRadius);
            delete mapper;

        }

    }
    
    CPPUNIT_TEST_SUITE( TestVertexMorphingMapper );
    CPPUNIT_TEST( testHatTriasToTrias);
    CPPUNIT_TEST( testHatQuadsToQuads);
    CPPUNIT_TEST( testHatTriasToQuads);
    CPPUNIT_TEST( testHatQuadsToTrias);
    CPPUNIT_TEST( testConsistency);
    
    //     CPPUNIT_TEST( testMemoryLeakOfConstructor); // test memory leak, comment it except when checking memory leak
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestVertexMorphingMapper);
