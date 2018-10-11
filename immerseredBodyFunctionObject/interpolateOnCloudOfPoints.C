/*---------------------------------------------------------------------------*\
    =========                 |
    \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
    \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
    -------------------------------------------------------------------------------
    License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

    \*---------------------------------------------------------------------------*/

#include "interpolateOnCloudOfPoints.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "surfaceFields.H"
#include "scalar.H"
#include "tensor2D.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
defineTypeNameAndDebug(interpolateOnCloudOfPoints, 0);
addToRunTimeSelectionTable(functionObject, interpolateOnCloudOfPoints, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::interpolateOnCloudOfPoints::interpolateOnCloudOfPoints
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
     
)
    :

     
    fvMeshFunctionObject(name, runTime, dict),
    patchNames_(dict.lookup("patches")),

    pointcloud_
    (
        IOobject
        (
            "cloud_of_points",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE                    //NO_WRITE
        )
    ),

    current_pointcloud_
    (
        IOobject
        (
            "currentcloudofpoints",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE                    //NO_WRITE
        )
    ),
 
    connectivity_matrix_
    (
        IOobject
        (
            "connectivity",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE                    //NO_WRITE
        )
    ),
    area_(connectivity_matrix_.size(),0.0) 	 
{
    read(dict);

    forAll(area_,i)
    {
        vector  p0=pointcloud_[connectivity_matrix_[i][0]];
        vector  p1=pointcloud_[connectivity_matrix_[i][1]];
        vector  p2=pointcloud_[connectivity_matrix_[i][2]];

        area_[i]= 0.5*mag( ((p1-p0)^(p2-p0)) );
 
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::interpolateOnCloudOfPoints::~interpolateOnCloudOfPoints()
{}
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::interpolateOnCloudOfPoints::read(const dictionary& dict)
{
    return fvMeshFunctionObject::read(dict);
}


Foam::scalar Foam::functionObjects::interpolateOnCloudOfPoints::weight(Foam::vector xcloud, Foam::vector xgrid)
{

    scalar hx=0.01;  //TODO read from mesh
    scalar hy=0.01;
    scalar phi1;
    scalar phi2;

    ///computing scaled distance between grid point and cloud point 
    scalar r1 = (1/hx)*(xgrid[0]-xcloud[0]);
    scalar r2 = (1/hy)*(xgrid[1]-xcloud[1]);

    r1=fabs(r1);
    r2=fabs(r2);
  

    ///computing phi1

    if (r1<=1.0)
    {
        phi1=(1.0/8.0)*(3.0-2.0*r1+sqrt(1.0+4.0*r1-4.0*r1*r1));
    } 
    else if(r1>=1.0 && r1<=2.0)
    {
        phi1=(1.0/8.0)*(5.0-2.0*r1-sqrt(-7.0+12.0*r1-4.0*r1*r1));
    } 
    else
    {
        phi1=0.0;
    } 


    ///computing phi2

    if (r2<=1.0)
    {
        phi2=(1.0/8.0)*(3.0-2.0*r2+sqrt(1.0+4.0*r2-4.0*r2*r2));
    } 
    else if(r2>=1.0 && r2<=2.0)
    {
        phi2=(1.0/8.0)*(5.0-2.0*r2-sqrt(-7.0+12.0*r2-4.0*r2*r2));
    } 
    else
    {
        phi2=0.0;
    } 


    scalar w=(phi1/hx)*(phi2/hy);

    return w;

}




//**********************************************************************************************************************************//




bool Foam::functionObjects::interpolateOnCloudOfPoints::execute()
{

    scalar hx=0.01;  //TODO read from mesh
    scalar hy=0.01;
    scalar delta=.001;  //TODO read time step
    //FOR LOOP RUNNING OVER ALL CLOUD POINTS

    forAll(current_pointcloud_,i)
    {

			 
		

        //COMPUTING VELOCITY AT THE PROBE POINT

        vector probePoint=current_pointcloud_[i];               
        vector velocity(0.0,0.0,0.0);                                     //initialising probe point velocity
        
        const volVectorField& U =mesh_.lookupObject<volVectorField>("U");

		
        labelHashSet cell_neighbors=findNeighbourCells(probePoint);

	
			  	
        forAllConstIter(labelHashSet, cell_neighbors,iter) 
        {
            label celli = iter.key();
            vector cellicentre=mesh_.cellCentres()[celli];
            scalar W=weight(probePoint,cellicentre);

            //Info<<"U of neighbour cell"<<celli<<"is"<<U[celli]<<endl;
            //Info<<"Weight of neighbour cell"<<celli<<"is"<<W<<endl;

									   

            velocity=velocity+(W*U[celli]*hx*hy);

        }
				 
        Info<<"Velocity at the probe point is "<<"  "<<velocity<<" "<<"point"<<current_pointcloud_[i]<<endl;

        //MOVE THE POINT FOR THE FIRST TIME STEP//

        current_pointcloud_[i]=probePoint+ velocity*(delta);     
 

    }





    //
    vectorField force(pointcloud_.size(),Zero);           ///Initialising Force at each vertex

    forAll(area_,i)                                     ///For loop runs over each element
    {

        Info<<"element"<<i<<endl;

        //Coordinates of reference triangle:
        vector  s0=pointcloud_[connectivity_matrix_[i][0]];   
        vector  s1=pointcloud_[connectivity_matrix_[i][1]];
        vector  s2=pointcloud_[connectivity_matrix_[i][2]];

        //Coordinates of deformed triangle:
        vector  x0=current_pointcloud_[connectivity_matrix_[i][0]];   
        vector  x1=current_pointcloud_[connectivity_matrix_[i][1]];
        vector  x2=current_pointcloud_[connectivity_matrix_[i][2]];

        //difference between coordinates of reference nodes:
        vector  ds1=s1-s0;
        vector  ds2=s2-s0;
                 

        //difference between coordinates of reference nodes:
        vector  dx1=x1-x0;
        vector  dx2=x2-x0;

 
        //Deformation gradient:

        //First row:
        vector2D  a0( (dx1[0]*ds2[1]-dx2[0]*ds1[1]) , (-dx1[0]*ds2[0]+dx2[0]*ds1[0]) );
                    
        //Second row:
        vector2D  a1( (dx1[1]*ds2[1] - dx2[1]*ds1[1]) , (-dx1[1]*ds2[0] + dx2[1]*ds1[0]) );
                     
        //Tensor s:

        tensor2D s(ds1[0],ds2[0],ds1[1],ds2[1]);

        //Info<<"area"<<area_[i]<<endl;

        a0=(1.0/(det(s)))*a0;
        a1=(1.0/(det(s)))*a1;
                     

        Info<<"def_grad"<<a0<<endl;
        Info<<"def_grad"<<a1<<endl;

        //Piola Kirchoff Stress:

        //For simplicity we set deformationgradient=Piola Kirchoff stress

        vector2D P0=a0;
        vector2D P1=a1;
                    
        //Derivative of deformation gradient wrt position:
        vector2D dela_delx0( ds2[1]-ds1[1] , -ds2[0]+ds1[0]);
        vector2D dela_delx1( ds2[1] ,  -ds2[0] );
        vector2D dela_delx2(-ds1[1] ,  ds1[0] );

 
        dela_delx0 = -(1.0/(det(s)))*dela_delx0;
        dela_delx1 =  (1.0/(det(s)))*dela_delx1;
        dela_delx2 =  (1.0/(det(s)))*dela_delx2;
	         

                    
        vector f0( P0&dela_delx0, P1&dela_delx0, 0.0);
        vector f1( P0&dela_delx1, P1&dela_delx1, 0.0);
        vector f2( P0&dela_delx2, P1&dela_delx2, 0.0);
                   

        f0 = -area_[i]*f0;
        f1 = -area_[i]*f1;
        f2 = -area_[i]*f2;
                   
        //Info<<"f0"<<f0<<endl;
        //Info<<"f1"<<f1<<endl;
        //Info<<"f2"<<f2<<endl;
                  

        force[connectivity_matrix_[i][0]]=force[connectivity_matrix_[i][0]] + f0;
        force[connectivity_matrix_[i][1]]=force[connectivity_matrix_[i][1]] + f1;
        force[connectivity_matrix_[i][2]]=force[connectivity_matrix_[i][2]] + f2;

    }  

                    
    volVectorField& f = const_cast<volVectorField&>(mesh_.lookupObject<volVectorField>("f"));

    //FOR LOOP RUNNING OVER ALL CLOUD POINTS
 
    forAll(current_pointcloud_,i)
    {

			 
			 

        vector probePoint=current_pointcloud_[i];               
        vector probePoint_force=force[i];

			    
        labelHashSet cell_neighbors=findNeighbourCells(probePoint);

			  	
        forAllConstIter(labelHashSet, cell_neighbors,iter) 
        {
            label celli = iter.key();
            vector cellicentre=mesh_.cellCentres()[celli];
            scalar W=weight(cellicentre,probePoint);


            //Info<<"Weight of neighbour cell"<<celli<<"is"<<W<<endl;
								 	    
            f[celli]=f[celli]+W*probePoint_force;
            //Info<<"f of neighbour cell"<<celli<<"is"<<f[celli]<<endl;
							
        }
 

    }

  
    //Info<<"grid_force"<<f<<endl;

















 

 

      
    return true;
}



//TODO check whether the cloud point is in mesh or not

Foam::labelHashSet Foam::functionObjects::interpolateOnCloudOfPoints::findNeighbourCells(const vector probePoint) const
{
    
	
    // find cell containing this point
    label celli = mesh_.findCell(probePoint);


    // container for neighbours set by dumping the cell containing it
    labelHashSet neighbourCellSet(0);
    neighbourCellSet.set(celli);


    // number of layers
    int nLayers = 2;
    for (int n = 0; n < nLayers; n++)
    {
        // make a copy of marked cells
        labelHashSet markedNeighbours = neighbourCellSet;

        // loop over all marked cells
        forAllConstIter(labelHashSet, markedNeighbours,iter)
        {
            celli = iter.key();

            // get points of celli
            labelList celliPoints = mesh_.cellPoints()[celli];

            forAll(celliPoints,j)
            {
                // get neighbor cells of j th point
                labelList cellJNeighbours = mesh_.pointCells()[celliPoints[j]];

                // append these cells in neighbourCellSet
                forAll(cellJNeighbours, k)
                {
                    neighbourCellSet.set(cellJNeighbours[k]);
                }
        
            }
        }  
    }

    return neighbourCellSet;
}







bool Foam::functionObjects::interpolateOnCloudOfPoints::write()
{
    return true;
}


// ************************************************************************* //
