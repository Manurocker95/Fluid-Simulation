/*
	Manuel Rodríguez Matesanz
	URJC - Simulación Avanzada
	Práctica 3 - Fluidos
	Abril 2018

*/

#include <algorithm>
#include <math.h>
#include <stdio.h>

#include "scene.h"
#include "pcg_solver.h"

namespace
{
    //////////////////////////////////////////////
    // Add any custom classes or functions here! //
    //////////////////////////////////////////////
	
	// Dunno why this doesn't do what it should ://///
	// https://helloacm.com/cc-function-to-compute-the-bilinear-interpolation/
	// https://es.wikipedia.org/wiki/Interpolaci%C3%B3n_bilineal
	float BilinearInterpolation(float q11, float q12, float q21, float q22, float x1, float x2, float y1, float y2, float x, float y)
	{
		float x2x1, y2y1, x2x, y2y, yy1, xx1;
		x2x1 = x2 - x1;
		y2y1 = y2 - y1;
		x2x = x2 - x;
		y2y = y2 - y;
		yy1 = y - y1;
		xx1 = x - x1;

		float denom = (x2x1 * y2y1);

		
		// WHAT TO DO IF DENOM IS 0???????
		if (denom == 0)
			denom = 1.0f;

		float value = (1.0f / denom) * (
			q11 * x2x * y2y +
			q21 * xx1 * y2y +
			q12 * x2x * yy1 +
			q22 * xx1 * yy1
			);

		return value;
	}
}

/*
* 1) Copy grid velocity components  
* 2) For each one. getFacePos based on i,j to get the position
* 3) Interpolate for y between (u, (v1,v2,v3,v4)/4)
* 4) Pos - v*dt = pt-1 <- use bilinear func
*/


// advection
void Fluid2::fluidAdvection(const float dt)
{


	// ink advection
	{
		// We copy de ink grid
		Array2< float > ink_Copy;
		ink_Copy.copy(ink);

		unsigned int minimumX = 0;
		unsigned int maximumX = ink_Copy.getSize().x;
		unsigned int minimumY = 0;
		unsigned int maximumY = ink_Copy.getSize().y;


		// Current node id
		Index2 id(0, 0);
		// Index i+1
		Index2 id_nextX(0, 0);
		// Index j+1
		Index2 id_nextY(0, 0);

		Vec2 currentCellPos;
		Vec2 oldCellPos;

		Vec2 velocity;

		float vxInNode = 0;
		float vyInNode = 0;

		float vxInNextNode = 0;
		float vyInNextNode = 0;

		Vec2 oldCellPosIndex;

		float oldInkInterpolated = 0;

		// values to clamp interpolation
		Vec2 posMin;
		Vec2 posMax;

		// corner indices
		Index2 Q11;
		Index2 Q21;
		Index2 Q12;
		Index2 Q22;

		for (unsigned int i = 0; i < maximumX; i++)
		{
			for (unsigned int j = 0; j < maximumY; j++)
			{

				id.x = i;
				id.y = j;

				id_nextX.x = i + 1;
				id_nextX.y = j;

				id_nextY.x = i;
				id_nextY.y = j + 1;

				// velocity in i,j
				vxInNode = velocityX.getValue(id);
				vyInNode = velocityY.getValue(id);

				vxInNextNode = velocityX.getValue(id_nextX);
				vyInNextNode = velocityY.getValue(id_nextY);

				// Velocity used to integrate
				velocity.x = (vxInNode + vxInNextNode) / 2;
				velocity.y = (vyInNode + vyInNextNode) / 2;

				currentCellPos = grid.getCellPos(id);
				// We get old position (we can do x = x0 + vt for next pos and x0 = x - vt for previous one)
				oldCellPos = currentCellPos - dt * velocity;
				oldCellPosIndex = grid.getCellIndex(oldCellPos);

				
				// values to clamp 
				posMin.x = floorf(oldCellPosIndex.x);
				posMin.y = floorf(oldCellPosIndex.y);

				posMax.x = ceilf(oldCellPosIndex.x);
				posMax.y = ceilf(oldCellPosIndex.y);

				Q11.x = clamp((unsigned int)posMin.x, minimumX, maximumX - 1);
				Q11.y = clamp((unsigned int)posMin.y, minimumY, maximumY - 1);

				Q12.x = clamp((unsigned int)posMin.x, minimumX, maximumX - 1);
				Q12.y = clamp((unsigned int)posMax.y, minimumY, maximumY - 1);

				Q21.x = clamp((unsigned int)posMax.x, minimumX, maximumX - 1);
				Q21.y = clamp((unsigned int)posMin.y, minimumY, maximumY - 1);

				Q22.x = clamp((unsigned int)posMax.x, minimumX, maximumX - 1);
				Q22.y = clamp((unsigned int)posMax.y, minimumY, maximumY - 1);

				//(float q11, float q12, float q21, float q22, float x1, float x2, float y1, float y2, float x, float y)
				oldInkInterpolated = BilinearInterpolation(ink_Copy.getValue(Q11), ink_Copy.getValue(Q12), ink_Copy.getValue(Q21), ink_Copy.getValue(Q22), posMin.x, posMax.x, posMin.y, posMax.y, oldCellPosIndex.x, oldCellPosIndex.y);
				ink.setValue(id, oldInkInterpolated);
			}
		}
	}


	// velocity advection in X Axis
	{
		// We copy de velocty X grid
		Array2< float > velocityX_Copy;
		velocityX_Copy.copy(velocityX);

		unsigned int minimumX = 0;
		unsigned int maximumX = velocityX_Copy.getSize().x;
		unsigned int minimumY = 0;
		unsigned int maximumY = velocityX_Copy.getSize().y;

		Index2 id(0, 0);
		Index2 id1(0, 0);
		Index2 id2(0, 0);
		Index2 id3(0, 0);
		Index2 id4(0, 0);

		Vec2 currentFacePos(0, 0);
		Vec2 oldFacePos(0, 0);

		Vec2 velocity(0, 0);

		float vxInNode = 0;
		float vyInNode = 0;

		float interpolatedXVelocity = 0;
		Vec2 oldFacePosIndex;

		// values to clamp interpolation
		Vec2 posMin;
		Vec2 posMax;

		// corner indices
		Index2 Q11;
		Index2 Q21;
		Index2 Q12;
		Index2 Q22;

		/*

		--------0,2-------------1,2------
		|				 |				|
		0,1		0,1		1,1		1,1	   2,1
		|				 |				|
		------- 0,1-------------1,1------
		|				 |				|
		0,0		0,0		1,0		1,0	   2,0
		|				 |				|
		--------0,0-------------1,0------

			if ID = 1,0 (i,j)
			ID1 = 0,0 (i-1,j)
			ID2 = 1,0 (i,j)
			ID3 = 0,1 (i-1,j+1)
			ID4 = 1,1 (i, j+1)
		
		*/

		// We loop through each node (i,j)
		for (int i = 0; i < maximumX; i++)
		{
			for (int j = 0; j < maximumY; j++)
			{
				id.x = i;
				id.y = j;

				id1.x = i - 1;
				id1.y = j;
				id1.x = clamp(id1.x, minimumX, maximumX - 1);
				id1.y = clamp(id1.y, minimumY, maximumY - 1);

				id2.x = i;
				id2.y = j;
				id2.x = clamp(id2.x, minimumX, maximumX - 1);
				id2.y = clamp(id2.y, minimumY, maximumY - 1);

				id3.x = i - 1;
				id3.y = j + 1;
				id3.x = clamp(id3.x, minimumX, maximumX - 1);
				id3.y = clamp(id3.y, minimumY, maximumY - 1);

				id4.x = i;
				id4.y = j + 1;
				id4.x = clamp(id4.x, minimumX, maximumX - 1);
				id4.y = clamp(id4.y, minimumY, maximumY - 1);

				// (u, (v1,v2,v3,v4)/4)
				vxInNode = velocityX_Copy.getValue(id);
				vyInNode = (velocityY.getValue(id1) + velocityY.getValue(id2) + velocityY.getValue(id3) + velocityY.getValue(id4)) / 4;

				velocity.x = vxInNode;
				velocity.y = vyInNode;

				currentFacePos = grid.getFaceXPos(id);
				// We get old position (we can do x = x0 + vt for next pos and x0 = x - vt for previous one)
				oldFacePos = currentFacePos - dt * velocity;
				oldFacePosIndex = grid.getFaceIndex(oldFacePos, 0);

				// values to clamp 
				posMin.x = floorf(oldFacePosIndex.x);
				posMin.y = floorf(oldFacePosIndex.y);

				posMax.x = ceilf(oldFacePosIndex.x);
				posMax.y = ceilf(oldFacePosIndex.y);

				Q11.x = clamp((unsigned int)posMin.x, minimumX, maximumX - 1);
				Q11.y = clamp((unsigned int)posMin.y, minimumY, maximumY - 1);

				Q12.x = clamp((unsigned int)posMin.x, minimumX, maximumX - 1);
				Q12.y = clamp((unsigned int)posMax.y, minimumY, maximumY - 1);

				Q21.x = clamp((unsigned int)posMax.x, minimumX, maximumX - 1);
				Q21.y = clamp((unsigned int)posMin.y, minimumY, maximumY - 1);

				Q22.x = clamp((unsigned int)posMax.x, minimumX, maximumX - 1);
				Q22.y = clamp((unsigned int)posMax.y, minimumY, maximumY - 1);

				//(float q11, float q12, float q21, float q22, float x1, float x2, float y1, float y2, float x, float y)
				interpolatedXVelocity = BilinearInterpolation(velocityX_Copy.getValue(Q11), velocityX_Copy.getValue(Q12), velocityX_Copy.getValue(Q21), velocityX_Copy.getValue(Q22), posMin.x, posMax.x, posMin.y, posMax.y, oldFacePosIndex.x, oldFacePosIndex.y);
				velocityX.setValue(id, interpolatedXVelocity);
			}
		}
	}

	// velocity Y advection
	{
		// We copy de velocty Y grid
		Array2< float > velocityY_Copy;
		velocityY_Copy.copy(velocityY);

		unsigned int minimumX = 0;
		unsigned int maximumX = velocityY_Copy.getSize().x;
		unsigned int minimumY = 0;
		unsigned int maximumY = velocityY_Copy.getSize().y;

		Index2 id(0, 0);
		Index2 id1(0, 0);
		Index2 id2(0, 0);
		Index2 id3(0, 0);
		Index2 id4(0, 0);

		Vec2 currentFacePos(0, 0);
		Vec2 oldFacePos(0, 0);

		Vec2 velocity(0, 0);

		float vxInNode = 0;
		float vyInNode = 0;

		float interpolatedYVelocity = 0;
		Vec2 oldFacePosIndex;

		// values to clamp interpolation
		Vec2 posMin;
		Vec2 posMax;

		// corner indices
		Index2 Q11;
		Index2 Q21;
		Index2 Q12;
		Index2 Q22;

		/*

		--------0,2-------------1,2------
		|				 |				|
		0,1		0,1		1,1		1,1	   2,1
		|				 |				|
		------- 0,1-------------1,1------
		|				 |				|
		0,0		0,0		1,0		1,0	   2,0
		|				 |				|
		--------0,0-------------1,0------

		if ID = 0,1 (i,j)
		ID1 = 0,0 (i,j-1)
		ID2 = 1,0 (i+1,j-1)
		ID3 = 0,1 (i,j)
		ID4 = 1,1 (i+1, j)

		*/

		// We loop through each node (i,j)
		for (int i = 0; i < maximumX; i++)
		{
			for (int j = 0; j < maximumY; j++)
			{
				id.x = i;
				id.y = j;

				id1.x = i;
				id1.y = j - 1;
				id1.x = clamp(id1.x, minimumX, maximumX - 1);
				id1.y = clamp(id1.y, minimumY, maximumY - 1);

				id2.x = i + 1;
				id2.y = j - 1;
				id2.x = clamp(id2.x, minimumX, maximumX - 1);
				id2.y = clamp(id2.y, minimumY, maximumY - 1);

				id3.x = i;
				id3.y = j;
				id3.x = clamp(id3.x, minimumX, maximumX - 1);
				id3.y = clamp(id3.y, minimumY, maximumY - 1);

				id4.x = i + 1;
				id4.y = j;
				id4.x = clamp(id4.x, minimumX, maximumX - 1);
				id4.y = clamp(id4.y, minimumY, maximumY - 1);

				// ((u1,u2,u3,u4)/4, v)
				vxInNode = (velocityX.getValue(id1) + velocityX.getValue(id2) + velocityX.getValue(id3) + velocityX.getValue(id4)) / 4;
				vyInNode = velocityY_Copy.getValue(id);

				velocity.x = vxInNode;
				velocity.y = vyInNode;

				currentFacePos = grid.getFaceYPos(id);
				// We get old position (we can do x = x0 + vt for next pos and x0 = x - vt for previous one)
				oldFacePos = currentFacePos - dt * velocity;
				oldFacePosIndex = grid.getFaceIndex(oldFacePos, 1);

				// values to clamp 
				posMin.x = floorf(oldFacePosIndex.x);
				posMin.y = floorf(oldFacePosIndex.y);

				posMax.x = ceilf(oldFacePosIndex.x);
				posMax.y = ceilf(oldFacePosIndex.y);

				Q11.x = clamp((unsigned int)posMin.x, minimumX, maximumX - 1);
				Q11.y = clamp((unsigned int)posMin.y, minimumY, maximumY - 1);

				Q12.x = clamp((unsigned int)posMin.x, minimumX, maximumX - 1);
				Q12.y = clamp((unsigned int)posMax.y, minimumY, maximumY - 1);

				Q21.x = clamp((unsigned int)posMax.x, minimumX, maximumX - 1);
				Q21.y = clamp((unsigned int)posMin.y, minimumY, maximumY - 1);

				Q22.x = clamp((unsigned int)posMax.x, minimumX, maximumX - 1);
				Q22.y = clamp((unsigned int)posMax.y, minimumY, maximumY - 1);

				//(float q11, float q12, float q21, float q22, float x1, float x2, float y1, float y2, float x, float y)
				interpolatedYVelocity = BilinearInterpolation(velocityY_Copy.getValue(Q11), velocityY_Copy.getValue(Q12), velocityY_Copy.getValue(Q21), velocityY_Copy.getValue(Q22), posMin.x, posMax.x, posMin.y, posMax.y, oldFacePosIndex.x, oldFacePosIndex.y);
				velocityY.setValue(id, interpolatedYVelocity);
			}
		}	
	}
}

// emission
void Fluid2::fluidEmission()
{
    if( Scene::testcase >= Scene::SMOKE )
    {
        // emit source ink
        {
			
			for (unsigned int i = 2; i < ink.getSize().x / 4; ++i)
				for (unsigned int j = 2; j < ink.getSize().y / 4; ++j)
				{
					ink[Index2(i, j)] = 1.0f;
				}
		}

        // emit source velocity
        {
			// AS TEST_ADVECTION
			/*
			Array2< float >& u = velocityX;
			for (unsigned int i = 0; i < u.getSize().x; ++i)
			for (unsigned int j = 0; j < u.getSize().y; ++j)
			u[Index2(i, j)] = 2.0f;

			*/


			Array2< float >& v = velocityY;
			for (unsigned int i = 0; i < v.getSize().x; ++i)
				for (unsigned int j = 0; j < v.getSize().y; ++j)
					v[Index2(i, j)] = 2.0f;
	
        }
    }

}

// volume forces
void Fluid2::fluidVolumeForces( const float dt )
{
    if( Scene::testcase >= Scene::SMOKE )
    {
		Index2 id (0,0);
        // gravity
		for (int i = 0; i < velocityY.getSize().x; i++)
		{
			for (int j = 0; j < velocityY.getSize().y; j++)
			{
				id.x = 0;
				id.y = 0;
				velocityY.setValue(id, velocityY.getValue(id) + (Scene::kGravity*dt));
			}
		}

    }
}

// viscosity
void Fluid2::fluidViscosity( const float dt )
{
    if( Scene::testcase >= Scene::SMOKE )
    {
		// u*ij = uij + Dt / dens * coefVisc *((ui+1j - 2uij + ui-1,j)/(dx^2) + (uij+1 - 2uij + uij-1)/dy^2)
	
		// Velocity X
		{
			Array2< float > velocityX_Copy;
			velocityX_Copy.copy(velocityX);

			unsigned int minimumX = 0;
			unsigned int maximumX = velocityX_Copy.getSize().x;
			unsigned int minimumY = 0;
			unsigned int maximumY = velocityX_Copy.getSize().y;

			float dtViscosityDensity = (dt * Scene::kViscosity) / Scene::kDensity;
			float DX2 = pow(grid.getCellDx().x, 2);
			float DY2 = pow(grid.getCellDx().y, 2);

			Index2 id(0, 0);
			
			Index2 id1(0, 0);
			Index2 id2(0, 0);
			Index2 id3(0, 0);
			Index2 id4(0, 0);

			for (unsigned int i = 0; i < maximumX; i++) 
			{
				for (unsigned int j = 0; j < maximumY; j++)
				{
					id.x = i;
					id.y = j;

					id1.x = (clamp(i + 1, (unsigned int)0, maximumX - 1));
					id1.y = j;

					id2.x = (clamp(i - 1, (unsigned int)0, maximumX - 1));
					id2.y = j;

					id3.x = i;
					id3.y = clamp(j + 1, (unsigned int)0, maximumY - 1);

					id4.x = i;
					id4.y = clamp(j - 1, (unsigned int)0, maximumY - 1);

					float node1 = (velocityX_Copy.getValue(id1) - (2.0f * velocityX_Copy.getValue(id)) + velocityX_Copy.getValue(id2)) / DX2;
					float node2 = (velocityX_Copy.getValue(id3) - (2.0f * velocityX_Copy.getValue(id)) + velocityX_Copy.getValue(id4)) / DY2;

					velocityX.setValue(id, velocityX_Copy.getValue(id) + dtViscosityDensity * ((node1)+(node2)));
				}
			}

		}

		// Velocity Y
		{
			Array2< float > velocityY_Copy;
			velocityY_Copy.copy(velocityY);

			unsigned int minimumX = 0;
			unsigned int maximumX = velocityY_Copy.getSize().x;
			unsigned int minimumY = 0;
			unsigned int maximumY = velocityY_Copy.getSize().y;

			float dtViscosityDensity = (dt * Scene::kViscosity) / Scene::kDensity;
			float DX2 = pow(grid.getCellDx().x, 2);
			float DY2 = pow(grid.getCellDx().y, 2);

			Index2 id(0, 0);

			Index2 id1(0, 0);
			Index2 id2(0, 0);
			Index2 id3(0, 0);
			Index2 id4(0, 0);

			for (unsigned int i = 0; i < maximumX; i++)
			{
				for (unsigned int j = 0; j < maximumY; j++)
				{
					id.x = i;
					id.y = j;

					id1.x = (clamp(i + 1, (unsigned int)0, maximumX - 1));
					id1.y = j;

					id2.x = (clamp(i - 1, (unsigned int)0, maximumX - 1));
					id2.y = j;

					id3.x = i;
					id3.y = clamp(j + 1, (unsigned int)0, maximumY - 1);

					id4.x = i;
					id4.y = clamp(j - 1, (unsigned int)0, maximumY - 1);

					float node1 = (velocityY_Copy.getValue(id1) - (2.0f * velocityY_Copy.getValue(id)) + velocityY_Copy.getValue(id2)) / DX2;
					float node2 = (velocityY_Copy.getValue(id3) - (2.0f * velocityY_Copy.getValue(id)) + velocityY_Copy.getValue(id4)) / DY2;

					velocityX.setValue(id, velocityY_Copy.getValue(id) + dtViscosityDensity * ((node1)+(node2)));
				}
			}
		}   
    }
}

// pressure
void Fluid2::fluidPressureProjection( const float dt )
{
    if( Scene::testcase >= Scene::SMOKE )
    {
        // pressure
    }
}
