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
/***********************************************************************************************//**
 * \file MetaDataStructures.h
 * This file holds the the data structure of setting (metadata)
 * \author Tianyang Wang
 * \date 6/5/2012
 **************************************************************************************************/

#ifndef METADATASTRUCTURES_H_
#define METADATASTRUCTURES_H_

#include "EMPEROR_Enum.h"

namespace EMPIRE {
/*
 * The structs defined here should correspond to the counterpart in the XML inputDataField file.
 * Please mark the difference: sometimes a struct owns another struct, sometimes it refers to another struct.
 * There is such a difference also in the XML inputDataField file.
 */
struct structClientCode {
    struct structMesh {
        struct structDataField {
            std::string name;
            EMPIRE_DataField_location location;
            EMPIRE_DataField_dimension dimension;
            EMPIRE_DataField_typeOfQuantity typeOfQuantity;
        };
        std::string name;

        EMPIRE_Mesh_type type;
        std::string clientNameToCopyFrom;
        std::string meshNameToCopyFrom;
        bool sendMeshToClient;
        bool triangulateAll;
        std::vector<structDataField> dataFields;
    };
    struct structSignal {
        std::string name;
        int size3D[3];
    };
    std::string name;
    bool isRestart; // Flag mentioning if the code is restarted.
    std::vector<structMesh> meshes;
    std::vector<std::string> initialDataFields;
    std::vector<structSignal> signals;
};

struct structMeshRef {
    std::string clientCodeName;
    std::string meshName;
};

struct structDataFieldRef {
    std::string clientCodeName;
    std::string meshName;
    std::string dataFieldName;
};

struct structSignalRef {
    std::string clientCodeName;
    std::string signalName;
};

struct structConnectionIO {
    EMPIRE_ConnectionIO_Type type;
    structSignalRef signalRef;
    structDataFieldRef dataFieldRef;
};

struct structDataOutput {
    std::string name;
    int interval;
    std::vector<structConnectionIO> connectionIOs;
};

struct structResidual {
    struct structComponent {
        int coefficient;
        std::string timeToUpdate;
        structConnectionIO connectionIO;
    };
    int index;
    std::vector<structComponent> components;
};

struct structResidualRef {
    std::string couplingAlgorithmName;
    int index;
};

struct structMapper {
    struct structMortarMapper {
        bool oppositeSurfaceNormal;
        bool dual;
        bool enforceConsistency;
    };
    struct structIGAMortarMapper {
        struct propConsistency {
            bool enforceConsistency;
            double tolConsistency;
        } propConsistency;
        struct propProjection {
            double maxProjectionDistance;
            int noInitialGuess;
            double maxProjectionDistanceOnDifferentPatches;
        } propProjection;
        struct propNonlinearScheme {
            int noIterations;
            double tolProjection;
        } propNewtonRaphson, propNewtonRaphsonBoundary, propBisection;
        struct propIntegration {
            bool isAutomaticNoGPTriangle;
            int noGPTriangle;
            bool isAutomaticNoGPQuadrilateral;
            int noGPQuadrilateral;
        } propIntegration;
        struct propWeakCurveDirichletConditions {
                bool isWeakCurveDirichletConditions;
                bool isAutomaticPenaltyParameters;
                bool isPrimPrescribed;
                bool isSecBendingPrescribed;
                bool isSecTwistingPrescribed;
                double alphaPrim;
                double alphaSecBending;
                double alphaSecTwisting;
        } propWeakCurveDirichletConditions;
        struct propWeakSurfaceDirichletConditions {
            bool isWeakSurfaceDirichletConditions;
            bool isPrimPrescribed;
            bool isSecBendingPrescribed;
            bool isAutomaticPenaltyParameters;
            double alphaPrim;
        } propWeakSurfaceDirichletConditions;
        struct propWeakPatchContinuityConditions {
                bool isWeakPatchContinuityConditions;
                bool isAutomaticPenaltyParameters;
                bool isPrimCoupled;
                bool isSecBendingCoupled;
                bool isSecTwistingCoupled;
                double alphaPrim;
                double alphaSecBending;
                double alphaSecTwisting;
        } propWeakPatchContinuityConditions;
        struct propStrongCurveDirichletConditions {
            bool isStrongCurveDirichletConditions;
        } propStrongCurveDirichletConditions;
        struct propErrorComputation {
            bool isErrorComputation;
            bool isDomainError;
            bool isCurveError;
            bool isInterfaceError;
        } propErrorComputation;
    };
    struct structIGABarycentricMapper {
        struct propProjection {
            double maxProjectionDistance;
            int noInitialGuess;
            double maxProjectionDistanceOnDifferentPatches;
        } propProjection;
        struct propNonlinearScheme {
            int noIterations;
            double tolProjection;
        } propNewtonRaphson;
    };
    struct structCurveSurfaceMapper {
        EMPIRE_CurveSurfaceMapper_type type;
    };
    std::string name;
    structMeshRef meshRefA;
    structMeshRef meshRefB;
    EMPIRE_Mapper_type type;
    structMortarMapper mortarMapper;
    structIGAMortarMapper IGAMortarMapper;
    structIGABarycentricMapper IGABarycentricMapper;
    structCurveSurfaceMapper curveSurfaceMapper;
};

struct structCouplingAlgorithm {
    struct structAitken {
        double initialRelaxationFactor;
    };
    struct structConstantRelaxation {
        double relaxationFactor;
    };
    struct structOutput {
        int index;
        structConnectionIO connectionIO;
    };
    struct structInterfaceJacobian {
        unsigned int indexRow;
        unsigned int indexColumn;
        bool isConstant;
        double value;
        bool isAutoDiff;
        structConnectionIO functionInput;
        structConnectionIO functionOutput;
        bool isSignal;
        structConnectionIO interfaceJacobianSignal;
        double coefficient;
    };

    // Properties of GMRES algorithm
    // Aditya Ghantasala.
    // TODO :: Change accordingly.
    struct structGMRES {
    	 int maxOuterItter;
    	 int maxInnerItter;
    	 double residualTolerance;
         std::vector<structConnectionIO> inputs;
         std::vector<structConnectionIO> outputs;
/*       unsigned int indexRow;
         unsigned int indexColumn;
         bool isConstant;
         double value;
         bool isAutoDiff;
         structConnectionIO functionInput;
         structConnectionIO functionOutput;
         bool isSignal;
         structConnectionIO interfaceJacobianSignal;
         double coefficient;*/
     };


    std::string name;
    EMPIRE_CouplingAlgorithm_type type;
    std::vector<structOutput> outputs;
    std::vector<structInterfaceJacobian> interfaceJacobians;
    std::vector<structResidual> residuals;
    structAitken aitken;
    structConstantRelaxation constantRelaxation;
    structGMRES gmres;
};

struct structFilter {
    struct structMappingFilter {
        std::string mapperName;
    };
    struct structScalingFilter {
        double factor;
    };
    struct structWeakCouplingFilter {
        double beta;
    };
    struct structCopyFilter {
        int signalOffset;
    };
    struct structSetFilter {
        std::vector<double> value;
    };
    struct structDataFieldIntegrationFilter {
        structMeshRef meshRef;
    };
    struct structAdditionFilter {
        double a, b;
    };
    EMPIRE_DataFieldFilter_type type;
    structMappingFilter mappingFilter;
    structScalingFilter scalingFilter;
    structWeakCouplingFilter weakCouplingFilter;
    structSetFilter setFilter;
    structCopyFilter copyFilter;
    structDataFieldIntegrationFilter dataFieldIntegrationFilter;
    structAdditionFilter additionFilter;

    std::vector<structConnectionIO> inputs;
    std::vector<structConnectionIO> outputs;
};

struct structExtrapolator {
    std::string name;
    EMPIRE_Extrapolator_type type;
    std::vector<structConnectionIO> connectionIOs;
};

struct structConnection {
    std::string name;
    std::vector<structFilter> filterSequence;
    std::vector<structConnectionIO> inputs;
    std::vector<structConnectionIO> outputs;
};

struct structCouplingLogic {
    struct structConnectionRef {
        std::string connectionName;
    };
    struct structTimeStepLoop {
        int numTimeSteps;
        std::pair<bool, std::string> extrapolatorRef; // bool: has the ref or not
        std::vector<std::string> dataOutputRefs;
    };
    struct structIterativeCouplingLoop {
        struct structConvergenceChecker {
            struct structCheckResidual {
                double relativeTolerance;
                double absoluteTolerance;
                bool isAbsolute;
                structResidualRef residualRef;
            };
            double maxNumOfIterations;
            std::vector<structCheckResidual> checkResiduals;
        };
        structConvergenceChecker convergenceChecker;
        std::vector<std::string> convergenceObservers;
        std::vector<std::string> couplingAlgorithmRefs;
        std::vector<std::string> dataOutputRefs;
    };
    struct structOptimizationLoop {
        int maxNumOfIterations;
        std::vector<std::string> convergenceSignalReceivers;
        std::string convergenceSignalSender;
        std::vector<std::string> dataOutputRefs;
    };
    EMPIRE_CouplingLogic_type type;
    std::vector<structCouplingLogic> sequence;
    structConnectionRef connectionRef;
    structTimeStepLoop timeStepLoop;
    structIterativeCouplingLoop iterativeCouplingLoop;
    structOptimizationLoop optimizationLoop;
};

} /* namespace EMPIRE */
#endif /* METADATASTRUCTURES_H_ */
