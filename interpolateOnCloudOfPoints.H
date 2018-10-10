/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
-------------------------------------------------------------------------------

Author: Ajit Kumar, Shiv Nadar University, ajitkumar22@gmail.com
Date:   May 29, 2018
*/


#ifndef functionObjects_interpolateOnCloudOfPoints_H
#define functionObjects_interpolateOnCloudOfPoints_H

#include "fvMeshFunctionObject.H"

#include <fstream>
#include <iostream>
#include "pointIOField.H"
#include "dictionary.H"
 
#include "labelListIOList.H"
#include "scalar.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class interpolateOnCloudOfPoints Declaration
\*---------------------------------------------------------------------------*/

class interpolateOnCloudOfPoints
:
    public fvMeshFunctionObject
{


        wordList patchNames_;

	    pointIOField pointcloud_;

        pointIOField current_pointcloud_; 
        
        labelListIOList connectivity_matrix_;       
  
    
    // Private member functions

        //- Disallow default bitwise copy construct
        interpolateOnCloudOfPoints(const interpolateOnCloudOfPoints&);

        //- Disallow default bitwise assignment
        void operator=(const interpolateOnCloudOfPoints&);
		  



private:
       scalarField area_;   


public:

    //- Runtime type information
    TypeName("interpolateOnCloudOfPoints");


    // Constructors

        //- Construct from Time and dictionary
        interpolateOnCloudOfPoints
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
             
        );


    //- Destructor
    virtual ~interpolateOnCloudOfPoints();


    // Member Functions

        //- Read the cell-centre rate data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write the cell-centre fields
        virtual bool write();

        //- Computes weight at the grid point location for given cloud point coordinate

        scalar weight(vector,vector);

        //- Computes neighbour cells of a given Lagrangian point(probe point)

		Foam::labelHashSet findNeighbourCells(const vector probePoint) const;


         

  
         

   


        

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
