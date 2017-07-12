/*  Copyright &copy; 2014, TU Muenchen, Chair of Structural Analysis,
 *  Fabien Pean, Andreas Apostolatos, Chenshen Wu, Munich
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
 * \file IGAPatchSurfaceTrimming.h
 * This file holds the class IGAPatchSurfaceTrimming.h
 * \date 28/5/2013
 **************************************************************************************************/
 
 #ifndef IGAPatchSurfaceTrimming_H_
 #define IGAPatchSurfaceTrimming_H_
 
 // Inclusion of user defined libraries
 #include "IGAPatchCurve.h"
 #include "IGAControlPoint.h"
 #include <vector>
 #include <utility>
 #include <assert.h>

 
 namespace EMPIRE {
     class IGAPatchSurface;
     class IGAControlPoint;
     class Message;
     class IGAPatchSurfaceTrimmingLoop;
     
     /********//**
     * \brief class IGAPatchSurfaceTrimming is a container/proxy to the set of loops trimming the parent Patch
     ***********/
     class IGAPatchSurfaceTrimming {
     private:

     protected:
    	 std::vector<IGAPatchSurfaceTrimmingLoop*> loops;

    	 std::vector<int> outter;
         
         /// The constructor and the destructor and the copy constructor
     public:
         /***********************************************************************************************
          * \brief  Default constructor
          * \author Fabien Pean
          ***********/
         IGAPatchSurfaceTrimming();
         
         /***********************************************************************************************
          * \brief Destructor
          * \author Fabien Pean
          ***********/
         ~IGAPatchSurfaceTrimming();

         /***********************************************************************************************
          * \brief Setup information about the loop soon to be received
          * \param[in] inner 0 for outter and 1 for inner
          * \param[in] numCurves Number of curves to be received for this loop 
          * \author Fabien Pean
          ***********/
         void addTrimLoop(int _inner, int _numCurves);
         
         /***********************************************************************************************
          * \brief Add a Nurbs curve for the current loop and its attached information
          * \param[in] direction The direction of the curve if is following standard or not
          * \param[in] ID The id of the curve
          * \param[in] _pDegree The polynomial degree of the IGA 1D curve in the u-direction
          * \param[in] _uNoKnots The number of knots for the knot vector in the u-direction
          * \param[in] _uKnotVector The underlying knot vector of the IGA 1D curve in the u-direction
          * \param[in] _uNoControlPoints The number of the Control Points for the 1D NURBS patch in the u-direction
          * \param[in] _controlPointNet The set of the Control Points related to the 1D NURBS patch
          * \author Fabien Pean
          ***********/
         void addTrimCurve(int _direction,int _IDBasis, int _pDegree, int _uNoKnots, double* _uKnotVector,
                           int _uNoControlPoints, double* _controlPointNet);
         
         /***********************************************************************************************
          * \brief Create a linear approximation for every loop using position computed at Greville abscissae
          * 	   And remove non-unique points or aligned points from the set
          * \author Fabien Pean
          ***********/
         void linearizeLoops();
         
         /// Get and set functions
     public:
         /***********************************************************************************************
          * \brief Boolean indicating if there are trimming information
          * \author Fabien Pean
          ***********/
         inline bool isTrimmed() const {
            return (loops.size()>0);
         }
         /***********************************************************************************************
          * \brief Get the outter loop
          * \author Fabien Pean
          ***********/
         inline const std::vector<const IGAPatchSurfaceTrimmingLoop*> getOutterLoop() {
        	 std::vector<const IGAPatchSurfaceTrimmingLoop*> out;
        	 for(std::vector<int>::const_iterator it = outter.begin(); it != outter.end(); it++)
        		 out.push_back(loops[*it]);
             return out;
         }
         /***********************************************************************************************
          * \brief Get the outter loop index
          * \author Fabien Pean
          ***********/
         inline const std::vector<int>& getOutterLoopIndex() const {
            return outter;
         }
         inline int getOutterLoopIndex(int i) const {
            return outter.at(i);
         }
         /***********************************************************************************************
          * \brief Get a specific loop
          * \author Fabien Pean
          ***********/
         inline const IGAPatchSurfaceTrimmingLoop& getFirstLoop() const {
             return *(loops.front());
         }
         inline const IGAPatchSurfaceTrimmingLoop& getLoop(int i) const {
             return *(loops.at(i));
         }
         inline IGAPatchSurfaceTrimmingLoop& getLoop(int i) {
             return *(loops.at(i));
         }
         inline const IGAPatchSurfaceTrimmingLoop& getLastLoop() const {
             return *(loops.back());
         }
         /***********************************************************************************************
          * \brief Get vector of loops
          * \author Fabien Pean
          ***********/
         inline const std::vector<IGAPatchSurfaceTrimmingLoop*>& getLoops() const {
			 return loops;
         }
         /***********************************************************************************************
          * \brief Get the number of loops.
          * \author Fabien Pean
          ***********/
         inline int getNumOfLoops() const {
        	 return loops.size();
         }
     };

     /********//**
     * \brief class IGAPatchSurfaceTrimmingLoop holds all the curves for one trimming
     ***********/
	 class IGAPatchSurfaceTrimmingLoop {
    	 /// The parent class can access it freely
		 friend class IGAPatchSurfaceTrimming;
	 public:
         /***********************************************************************************************
          * \brief Connstructor, reserve data storage for n curves
          * \param[in] _numCurves The number of curves in the loop
          * \author Fabien Pean
          ***********/
		 IGAPatchSurfaceTrimmingLoop(int _numCurves);
         /***********************************************************************************************
          * \brief Destructor
          * \author Fabien Pean
          ***********/
		 ~IGAPatchSurfaceTrimmingLoop();
	 private:
         /// The basis functions of the curves
         std::vector<IGAPatchCurve*> IGACurves;
         /// Direction = if curve is in the correct orientation C1(1)=C2(0) or not C1(1)=C2(1)
         ///	WARNING : this has not been tested, so code may break if it is used
         std::vector<bool> direction;
         // List of points making up the linearized version of the trimming loop
         std::vector<double> polylines;

         /// Linearizing related functions
     public:
         /***********************************************************************************************
          * \brief Create a linear approximation for this loop
          * \author Fabien Pean
          ***********/
         void linearize();
         /***********************************************************************************************
          * \brief Create a linear approximation for this loop using 1 node per control point at Greville position
          * 	This is a coarse linearization whose deviation is limited by using the Greville abscissae
          * \author Fabien Pean
          ***********/
         void linearizeUsingGreville();
         /***********************************************************************************************
          * \brief Create a linear approximation for this loop using (#CPs times degree of the curve) nodes
          * 	This allows a scalable linearization, quite refined.
          * \author Fabien Pean
          ***********/
         void linearizeUsingNCPxP();

         void linearizeCombined();


         /// get functions
	 public:
         /***********************************************************************************************
          * \brief Get the underlying IsoGeometric curve of index i of the loop
          * \author Fabien Pean
          ***********/
         inline const IGAPatchCurve& operator[](int i) const {
             return (const IGAPatchCurve&)*IGACurves.at(i);
         }
         inline IGAPatchCurve& operator[](int i) {
             return *IGACurves.at(i);
         }
         /***********************************************************************************************
          * \brief Get the underlying IsoGeometric curve of index i of the loop
          * \author Fabien Pean
          ***********/
         inline const IGAPatchCurve& getIGACurve(int i) const {
             return (const IGAPatchCurve&)*IGACurves.at(i);
         }
         inline IGAPatchCurve& getIGACurve(int i) {
             return (IGAPatchCurve&)*IGACurves.at(i);
         }
         /***********************************************************************************************
          * \brief Get the underlying IsoGeometric curve of the loop
          * \author Fabien Pean
          ***********/
         inline const std::vector<IGAPatchCurve*>& getIGACurves() const {
             return IGACurves;
         }
         /***********************************************************************************************
          * \brief Get the number of curves in this loop
          * \author Fabien Pean
          ***********/
         inline int getNoCurves() const  {
             return IGACurves.size();
         }
         /***********************************************************************************************
          * \brief 	Get the direction for curve i
          * \author Fabien Pean
          ***********/
         inline bool getDirection(int i) const {
             return direction.at(i);
         }
         /***********************************************************************************************
          * \brief Get the number of the Control Points of the curve i
          * \author Fabien Pean
          ***********/
         inline int getNoControlPoints(int i) const {
             return IGACurves.at(i)->getNoControlPoints();
         }
         /***********************************************************************************************
          * \brief Get the Control Points of the curve i
          * \author Andreas Apostolatos
          ***********/
         inline const std::vector<IGAControlPoint>& getControlPointNet(int i) const {
             return IGACurves.at(i)->getControlPointNet();
         }
         /***********************************************************************************************
          * \brief Get the linearized version of the curve
          * \author Fabien Pean
          ***********/
         inline const std::vector<double>& getPolylines() const {
             return polylines;
         }
         inline const double* getPolylines(int* size) const {
        	 if(size!=NULL)
        		 *size=polylines.size();
        	 return &polylines[0];
         }

         /// The tolerance in the curve parameter space for the linearization of the trimming loops
         static double TOL_LINEARIZATION;
	 };
     
     /***********************************************************************************************
      * \brief Allows for nice debug output
      * \author Fabien Pean
      ***********/
     Message &operator<<(Message &message, const IGAPatchSurfaceTrimming &trim);
     Message &operator<<(Message &message, const IGAPatchSurfaceTrimmingLoop &trim);

     
 }/* namespace EMPIRE */
 
 #endif /* IGAPatchSurfaceTrimming_H_ */
 
