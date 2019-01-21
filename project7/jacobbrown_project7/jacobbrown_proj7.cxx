/*=========================================================================

				Project 7 for CIS 410 (W18)
				3D isosurfacing code
				implemented by Jacob Brown 2/22/2018

===========================================================================*/
/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include "TriangleList.h"

#define ISOVALUE 3.2f

// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    return dims[0]*dims[1]*dims[2];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    idx[0] = pointId%dims[0];
    idx[1] = (pointId/dims[0])%dims[1];
    idx[2] = pointId/(dims[0]*dims[1]);
}

// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    idx[2] = cellId/((dims[0]-1)*(dims[1]-1));
}

// Interpolation Function
float Interpolation(float value, float a, float b, float fa, float fb) {
	return (b - a) != 0 ? fa + ((value - a) / (b - a)) * (fb - fa) : 0;
}

int main()
{
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj6B.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
	float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

	TriangleList tl;

	// **Beginning of my code**
	int no_cells = GetNumberOfCells(dims);
	int no_triangles, case_id, v0_bin, v1_bin, v2_bin, v3_bin, v4_bin, v5_bin, v6_bin, v7_bin;
	int v0_log[3], v1_log[3], v2_log[3], v3_log[3], v4_log[3], v5_log[3], v6_log[3], v7_log[3];

	float v0_f, v1_f, v2_f, v3_f, v4_f, v5_f, v6_f, v7_f;

	int num_triangles[256];

	int lookup[256][16] = {
		{ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 0: aphilli9 */
		{ 0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 1: chatfiel */
		{ 0, 1,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 2: chesshir */
		{ 1, 3, 8, 1, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 3: cneugass */
		{ 2, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 4: cpalk */
		{ 0, 8, 10, 0, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 5: criegler */
		{ 3, 10, 2, 0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 6: cworkman */
		{ 1, 2, 10, 1, 9, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1 },  /* 7: ssane */
		{ 1, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 8: eewing */
		{ 1, 2, 11, 0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 9: gmorriso */
		{ 0, 2, 9, 2, 9, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 10: jbayes */
		{ 8, 2, 3, 8, 2, 11, 8, 9, 11, -1, -1, -1, -1, -1, -1, -1 },  /*11: Areej Alghamdi (3) 0000 1011*/
		{ 1, 3, 11, 3, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 12: jhorn3 */
		{ 0, 8, 10, 0, 1, 10, 1, 10, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 13: jlowen */
		{ 0, 3, 9, 3, 9, 10, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 14: pem */
		{ 8, 9, 11, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 15: kdawes */
		{ 7, 8, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 16: kpinto */
		{ 0, 3, 4, 3, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 17: mcmillan */
		{ 0, 1, 9, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 18: nivaldot */
		{ 1, 3, 7, 1, 4, 7, 1, 4, 9, -1, -1, -1, -1, -1, -1, -1 },  /* 19: sgrady2 */
		{ 7, 4, 8, 2, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 20: touermi */
		{ 0, 2, 4, 2, 4, 10, 4, 7, 10, -1, -1, -1, -1, -1, -1, -1 },  /* 21: tristanj */
		{ 0, 1, 9, 2, 3, 10, 4, 7, 9, -1, -1, -1, -1, -1, -1, -1 },  /* 22: aphilli9 */
		{ 1, 2, 9, 2, 9, 10, 4, 9, 10, 4, 7, 10, -1, -1, -1, -1 },  /* 23: chatfiel */
		{ 1, 2, 11, 4,  7,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 24: chesshir */
		{ 0, 3, 4, 0, 4, 7, 1, 2, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 25: cneugass */
		{ 0, 2, 9, 2, 9, 11,  4,  7,  8, -1, -1, -1, -1, -1, -1, -1 },  /* 26: cpalk */
		{ 2, 3, 7, 2, 9, 11, 4, 7, 9, 2, 7, 9, -1, -1, -1, -1 },  /* 27: criegler */
		{ 10, 3, 11, 4, 8, 10, 1, 11, 3, -1, -1, -1, -1, -1, -1, -1 },  /* 28: cworkman */
		{ 7,  4, 10,  0,  1,  4,  1,  4, 10,  1, 10, 11, -1, -1, -1, -1 },  /*29: David Stevens (4) */
		{ 4, 7, 8, 0, 3, 10, 0, 9, 10, 9, 10, 11, -1, -1, -1, -1 },  /* 30: eewing */
		{ 9, 10, 11, 4, 9, 10, 4, 7, 10, -1, -1, -1, -1, -1, -1, -1 },  /* 31: gmorriso */
		{ 4, 5, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 32: jbayes */
		{ 0, 3, 8, 4, 5, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 33: jbrawner */
		{ 0, 1, 5, 0, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 34: jhorn3 */
		{ 1, 3, 5, 3, 5, 8, 4, 5, 8, -1, -1, -1, -1, -1, -1, -1 },  /* 35: pem */
		{ 2, 3, 10, 4, 5, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 36: jnelson */
		{ 0, 2, 10, 0, 8, 10, 4, 5, 9, -1, -1, -1, -1, -1, -1, -1 },  /* 37: kdawes */
		{ 5, 4, 0, 2, 10, 3, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1 },  /* 38: kpinto */
		{ 1, 2, 5, 2, 5, 10, 5, 10, 4, 10, 4, 8, -1, -1, -1, -1 },  /* 39: mcmillan */
		{ 1, 2, 10, 4, 5, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 40: nivaldot */
		{ 0, 3, 8, 1, 2, 11, 4, 5, 9, -1, -1, -1, -1, -1, -1, -1 },  /* 41: sgrady2 */
		{ 0, 2, 4, 2, 5, 11, 2, 4, 5, -1, -1, -1, -1, -1, -1, -1 },  /* 42: touermi */
		{ 3, 4, 8, 3, 4, 5, 2, 3, 5, 2, 5, 11, -1, -1, -1, -1 },  /* 43: tristanj */
		{ 1, 2, 11, 2, 3, 10, 4, 5, 9, -1, -1, -1, -1, -1, -1, -1 },  /* 44: aphilli9 */
		{ 8, 10, 11, 4, 5, 9, 0, 1, 8, 1, 8, 11, -1, -1, -1, -1 },  /* 45: chatfiel */
		{ 0,  4,  5,  0,  3, 10,  0,  5, 10,  5, 10, 11, -1, -1, -1, -1 },  /* 46: chesshir */
		{ 5, 8, 11, 8, 10, 11, 4, 5, 8, -1, -1, -1, -1, -1, -1, -1 },  /* 47: cneugass */
		{ 5,  7,  9,  7,  8,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 48: cpalk */
		{ 0, 3, 9, 3, 5, 7, 3, 5, 9, -1, -1, -1, -1, -1, -1, -1 },  /* 49: criegler */
		{ 0, 8, 7, 0, 7, 1, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1 },  /* 50: cworkman */
		{ 1, 3, 5, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 51: hang */
		{ 2, 3, 10, 7, 8, 9, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1 },  /* 52: eewing */
		{ 5, 7, 9, 2, 7, 10, 0, 2, 9, 6, 7, 9, -1, -1, -1, -1 },  /* 53: gmorriso */
		{ 2, 3, 10, 7, 1, 5, 8, 1, 7, 0, 1, 8, -1, -1, -1, -1 },  /* 54: jbayes */
		{ 10, 2, 1, 7, 5, 1, 10, 1, 7, -1, -1, -1, -1, -1, -1, -1 },  /* 55: jbrawner */
		{ 7, 8, 9, 9, 5, 7, 1, 2, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 56: jhorn3 */
		{ 0, 1, 2, 0, 2, 3, 0, 1, 5, 2, 3, 7, -1, -1, -1, -1 },  /* 57: jlowen */
		{ 0, 2, 8, 2, 5, 8, 2, 5, 11, 5, 7, 8, -1, -1, -1, -1 }  ,  /* 58: jnelson */
		{ 2, 5, 11, 2, 3, 5, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1 },  /* 59: kdawes */
		{ 9, 5, 8, 8, 7, 5, 11, 10, 3, 1, 3, 11, -1, -1, -1, -1 },  /* 60: kpinto */
		{ 0, 1, 9, 5, 7, 11, 7, 10, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 61: mcmillan */
		{ 0, 3, 8, 5, 7, 11, 5, 10, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 62: nivaldot */
		{ 5, 7, 10, 5, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 63: sgrady2 */
		{ 6, 7, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 64: liuly */
		{ 0, 3, 8, 6, 7, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 65: tristanj */
		{ 0, 1, 9, 6, 7, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 66: aphilli9 */
		{ 6, 7, 10, 1, 3, 8, 1, 8, 9, -1, -1, -1, -1, -1, -1, -1 },  /* 67: chatfiel */
		{ 2,  3,  7,  2,  6,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 68: chesshir */
		{ 0, 2, 6, 0, 6, 7, 0, 7, 8, -1, -1, -1, -1, -1, -1, -1 },  /* 69: cneugass */
		{ 2,  3,  7,  2,  6,  7,  0,  1,  9, -1, -1, -1, -1, -1, -1, -1 },  /* 70: cpalk */
		{ 6, 7, 8, 1, 8, 9, 1, 2, 6, 1, 6, 8, -1, -1, -1, -1 },  /* 71: criegler */
		{ 1, 2, 11, 6, 7, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 72: cworkman */
		{ 6, 7, 10, 0, 3, 8, 1, 2, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 73: liuly */
		{ 0, 2, 9, 2, 9, 11, 6, 7, 10, -1, -1, -1, -1, -1, -1, -1 },  /* 74: eewing */
		{ 8, 9, 11, 6, 7, 10, 2, 3, 11, 3, 8, 11, -1, -1, -1, -1 },  /* 75: gmorriso */
		{ 7, 3, 1, 7, 1, 11, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 76: jbayes */
		{ 7, 1, 11, 7, 6, 11, 0, 8, 1, 1, 8, 7, -1, -1, -1, -1 },  /* 77: jbrawner */
		{ 0, 9, 11, 0, 3, 7, 0, 7, 11, 6, 7, 11, -1, -1, -1, -1 },  /* 78: jhorn3 */
		{ 8, 9, 11, 7, 8, 11, 6, 7, 11, -1, -1, -1, -1, -1, -1, -1 },  /*79: Stephanie Labasan (3) */
		{ 4, 8, 6, 10, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }  ,  /* 80: jnelson */
		{ 0, 3, 6, 0, 4, 6, 3, 6, 10, -1, -1, -1, -1, -1, -1, -1 },  /* 81: kdawes */
		{ 4, 8, 6, 0, 1, 9, 8, 10, 6, -1, -1, -1, -1, -1, -1, -1 },  /* 82: kpinto */
		{ 1, 3, 9, 3, 4, 9, 3, 4, 10, 4, 6, 10, -1, -1, -1, -1 },  /* 83: mcmillan */
		{ 3, 4, 8, 2, 3, 4, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1 },  /* 84: nivaldot */
		{ 0, 2, 4, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 85: sgrady2 */
		{ 2, 4, 6, 0, 1, 9, 3, 4, 8, 2, 3, 4, -1, -1, -1, -1 },  /* 86: touermi */
		{ 2, 4, 6, 1, 2, 4, 1, 4, 9, -1, -1, -1, -1, -1, -1, -1 },  /* 87: tristanj */
		{ 1, 2, 11, 4, 7, 8, 6, 7, 10, -1, -1, -1, -1, -1, -1, -1 },  /* 88: aphilli9 */
		{ 1, 2, 11, 0, 3, 10, 0, 6, 10, 0, 4, 6, -1, -1, -1, -1 },  /* 89: chatfiel */
		{ 8,  4, 10, 10,  4,  6,  2,  0,  9,  9,  2, 11, -1, -1, -1, -1 },  /* 90: chesshir */
		{ 3, 9, 11, 2, 3, 11, 3, 6, 10, 3, 4, 10, 3, 4, 9, -1 },  /* 91: cneugass */
		{ 1,  3,  8,  1,  6, 11,  1,  6,  8,  4,  6,  8, -1, -1, -1, -1 },  /* 92: cpalk */
		{ 0, 4, 6, 0, 1, 11, 0, 6, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 93: criegler */
		{ 3, 9, 11, 3, 4, 8, 3, 0, 9, 3, 4, 6, 3, 6, 11, -1 },  /* 94: cworkman */
		{ 4, 6, 11, 4, 9, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 95: pem */
		{ 4, 5, 9, 6, 7, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 96: roba */
		{ 0, 3, 8, 4, 5, 9, 6, 7, 10, -1, -1, -1, -1, -1, -1, -1 },  /* 97: gmorriso */
		{ 10, 6, 7, 0, 1, 5, 0, 4, 5, -1, -1, -1, -1, -1, -1, -1 },  /* 98: jbayes */
		{ 3, 1, 5, 4, 3, 5, 7, 6, 10, 8, 3, 4, -1, -1, -1, -1 },  /* 99: jbrawner */
		{ 4, 5, 9, 2, 3, 7, 2, 6, 7, -1, -1, -1, -1, -1, -1, -1 },  /* 100: jhorn3 */
		{ 0, 2, 9, 4, 7, 8, 4, 5, 6, -1, -1, -1, -1, -1, -1, -1 },  /* 101: jlowen */
		{ 0, 1, 5, 3, 6, 7, 2, 3, 6, 0, 4, 5, -1, -1, -1, -1 } ,  /* 102: jnelson */
		{ 1, 2, 8, 1, 8, 5, 2, 6, 8, 4, 5, 8, 7, 8, 6, -1 },  /* 103: kdawes */
		{ 4, 5, 9, 6, 7, 10, 1, 2, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 104: kpinto */
		{ 0, 1, 9, 2, 3, 10, 4, 7, 8, 5, 6, 11, -1, -1, -1, -1 },  /* 105: mcmillan */
		{ 0, 2, 4, 2, 4, 5, 2, 5, 10, 6, 7, 11, -1, -1, -1, -1 },  /* 106: nivaldot */
		{ 2, 3, 5, 2, 5, 11, 3, 4, 5, 3, 4, 8, 6, 7, 10, -1 },  /* 107: sgrady2 */
		{ 4, 5, 9, 1, 6, 11, 1, 6, 7, 1, 3, 7, -1, -1, -1, -1 },  /* 108: touermi */
		{ 0, 7, 8, 0, 1, 7, 1, 6, 7, 0, 6, 7, 4, 5, 9, -1 },  /* 109: tristanj */
		{ 0, 1, 9, 1, 2, 11, 2, 3, 10, 4, 5, 9, 6, 7, 10, -1 },  /* 110: aphilli9 */
		{ 4, 5, 11, 6, 7, 11, 4, 8, 11, 7, 8, 11, -1, -1, -1, -1 },  /* 111: chatfiel */
		{ 8,  9, 10,  9,  5,  6,  9,  6, 10, -1, -1, -1, -1, -1, -1, -1 },  /* 112: chesshir */
		{ 0, 5, 6, 0, 5, 9, 0, 3, 6, 3, 6, 10, -1, -1, -1, -1 },  /* 113: cneugass */
		{ 0,  1,  5,  5,  6, 10,  0,  5, 10,  0,  8, 10, -1, -1, -1, -1 },  /* 114: cpalk */
		{ 3, 5, 6, 3, 6, 10, 1, 3, 5, -1, -1, -1, -1, -1, -1, -1 },  /* 115: criegler */
		{ 2,  8,  9,  2,  3,  9,  2,  5,  9,  2,  5,  6, -1, -1, -1, -1 },  /* 116: bishara */
		{ 0, 2, 9, 2, 5, 9, 2, 5, 6, -1, -1, -1, -1, -1, -1, -1 },  /* 117: josh */
		{ 0,  3,  8,  1,  2,  5,  2,  5,  6, -1, -1, -1, -1, -1, -1, -1 },  /* 118: annag */
		{ 1, 2, 6, 1, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 119: gmorriso */
		{ 10, 9, 5, 1, 7, 11, 5, 6, 10, 8, 9, 10, -1, -1, -1, -1 },  /* 120: jbayes */
		{ 3, 0, 10, 2, 11, 1, 6, 5, 9, 0, 9, 6, 6, 10, 0, -1 },  /* 121: jbrawner */
		{ 0, 2, 5, 2, 5, 11, 5, 6, 10, 5, 8, 10, 0, 5, 8, -1 },  /* 122: jhorn3 */
		{ 2, 3, 10, 5, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 123: jlowen */
		{ 3, 8, 6, 9, 8, 6, 9, 5, 6, 1, 3, 6, -1, -1, -1, -1 } ,  /* 124: jnelson */
		{ 0, 1, 11, 0, 5, 6, 0, 5, 9, 0, 6, 11, -1, -1, -1, -1 },  /* 125: kdawes */
		{ 11, 5, 6, 0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 126: kpinto */
		{ 5, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 127: mcmillan */
		{ 5,  6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 128: bishara */
		{ 0, 3, 8, 5, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 129: sgrady2 */
		{ 0, 1, 9, 5, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 130: touermi */
		{ 5, 6, 11, 1, 3, 8, 1, 3, 9, -1, -1, -1, -1, -1, -1, -1 },  /* 131: tristanj */
		{ 2, 3, 10, 5, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 132: aphilli9 */
		{ 5, 6, 11, 0, 2, 10, 0, 8, 10, -1, -1, -1, -1, -1, -1, -1 },  /* 133: chatfiel */
		{ 0,  1,  9,  5,  6, 11,  2,  3, 10, -1, -1, -1, -1, -1, -1, -1 },  /* 134: chesshir */
		{ 5, 6, 11, 2, 9, 10, 8, 9, 10, 1, 2, 9, -1, -1, -1, -1 },  /* 135: cneugass */
		{ 1,  2,  6,  1,  5,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 136: cpalk */
		{ 0, 3, 8, 1, 2, 6, 1, 5, 6, -1, -1, -1, -1, -1, -1, -1 },  /* 137: criegler */
		{ 9, 5, 6, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1 },  /* 138: cworkman */
		{ 1, 2, 9, 2, 9, 10, 9, 10, 8, -1, -1, -1, -1, -1, -1, -1 },  /* 139: annag */
		{ 1, 3, 5, 3, 6, 10, 3, 5, 6, -1, -1, -1, -1, -1, -1, -1 },  /*140: James Kress (3) */
		{ 0, 1, 5, 0, 8, 10, 5, 6, 10, 0, 5, 10, -1, -1, -1, -1 },  /* 141: gmorriso */
		{ 0, 5, 9, 0, 5, 6, 0, 3, 6, 3, 6, 10, -1, -1, -1, -1 },  /* 142: liuly */
		{ 6, 9, 5, 8, 9, 10, 10, 9, 6, -1, -1, -1, -1, -1, -1, -1 },  /* 143: jbrawner */
		{ 4, 7, 8, 5, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 144: jhorn3 */
		{ 5, 6, 11, 0, 3, 4, 3, 4, 7, -1, -1, -1, -1, -1, -1, -1 },  /* 145: jlowen */
		{ 0, 1, 9, 4, 7, 8, 5, 6, 11, -1, -1, -1, -1, -1, -1, -1 } ,  /* 146: jnelson */
		{ 1, 3, 7, 5, 6, 11, 1, 7, 9, 4, 7, 9, -1, -1, -1, -1 },  /* 147: kdawes */
		{ 7, 8, 4, 11, 5, 6, 2, 3, 10, -1, -1, -1, -1, -1, -1, -1 },  /* 148: kpinto */
		{ 0, 2, 4, 2, 4, 5, 2, 5, 11, 7, 6, 10, -1, -1, -1, -1 },  /* 149: mcmillan */
		{ 0, 1, 9, 2, 3, 11, 4, 7, 8, 5, 6, 10, -1, -1, -1, -1 },  /* 150: nivaldot */
		{ 1, 2, 9, 2, 9, 10, 4, 7, 10, 4, 9, 10, 5, 6, 11, -1 },  /* 151: sgrady2 */
		{ 1, 2, 6, 1, 5, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1 },  /* 152: touermi */
		{ 0, 3, 4, 3, 4, 7, 1, 2, 5, 2, 5, 6, -1, -1, -1, -1 },  /* 153: tristanj */
		{ 0, 1, 9, 1, 2, 11, 4, 7, 8, 5, 6, 11, -1, -1, -1, -1 },  /* 154: aphilli9 */
		{ 5, 6, 9, 4, 7, 9, 2, 3, 9, 2, 6, 9, 3, 7, 9, -1 },  /* 155: chatfiel */
		{ 8,  4,  7,  3,  1,  5,  3, 10,  5,  5, 10,  6, -1, -1, -1, -1 },  /* 156: chesshir */
		{ 5, 6, 10, 4, 7, 10, 0, 4, 10, 1, 5, 10, 0, 1, 10, -1 },  /* 157: cneugass */
		{ 4,  7,  8,  0,  5,  6,  0,  5,  9,  0,  3,  6,  3,  6, 10, -1 },  /* 158: cpalk */
		{ 5, 6, 9, 6, 9, 10, 4, 7, 9, 7, 9, 10, -1, -1, -1, -1 },  /* 159: criegler */
		{ 4, 9, 11, 6, 4, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 160: cworkman */
		{ 0, 3, 8, 4, 6, 11, 4, 9, 11, -1, -1, -1, -1, -1, -1, -1 },  /*161: Monisha Balireddi (3) */
		{ 0, 1, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1 },  /* 162: hang */
		{ 4, 6, 8, 1, 3, 8, 1, 6, 11, 1, 6, 8, -1, -1, -1, -1 },  /* 163: gmorriso */
		{ 6, 11, 4, 9, 11, 4, 10, 2, 3, -1, -1, -1, -1, -1, -1, -1 },  /* 164: jbayes */
		{ 10, 2, 8, 8, 2, 0, 6, 11, 4, 4, 11, 9, -1, -1, -1, -1 },  /* 165: jbrawner */
		{ 2, 3, 10, 0, 1, 6, 1, 6, 11, 0, 4, 6, -1, -1, -1, -1 },  /* 166: jhorn3 */
		{ 1, 2, 11, 4, 6, 8, 6, 8, 10, -1, -1, -1, -1, -1, -1, -1 },  /* 167: jlowen */
		{ 1, 2, 4, 1, 4, 9, 2, 4, 6, -1,-1, -1, -1, -1, -1, -1 } ,  /* 168: jnelson */
		{ 0, 3, 8, 1, 2, 9, 2, 4, 6, 2, 4, 9, -1, -1, -1, -1 },  /* 169: kdawes */
		{ 0, 2, 4, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 170: kpinto */
		{ 3, 4, 8, 2, 3, 4, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1 },  /* 171: mcmillan */
		{ 1, 3, 9, 3, 4, 9, 3, 4, 11, 4, 6, 11, -1, -1, -1, -1 },  /* 172: nivaldot */
		{ 0, 1, 8, 1, 4, 6, 1, 4, 9, 1, 6, 10, 1, 8, 10, -1 },  /* 173: sgrady2 */
		{ 0, 4, 6, 0, 3, 6, 3, 6, 10, -1, -1, -1, -1, -1, -1, -1 },  /* 174: touermi */
		{ 4, 6, 8, 6, 8, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 175: tristanj */
		{ 8,  9, 11,  6,  8, 11,  6,  7,  8, -1, -1, -1, -1, -1, -1, -1 },  /* 176: bishara */
		{ 0, 3, 7, 0, 9, 11, 6, 7, 11, 0, 7, 11, -1, -1, -1, -1 },  /* 177: chatfiel */
		{ 0,  1,  8,  1,  7,  8,  1,  7, 11,  7, 11,  6, -1, -1, -1, -1 },  /* 178: chesshir */
		{ 6, 7, 11, 1, 3, 7, 1, 7, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 179: cneugass */
		{ 2,  3, 10,  6,  7,  8,  6,  8, 11,  8,  9, 11, -1, -1, -1, -1 },  /* 180: cpalk */
		{ 0, 2, 7, 2, 7, 10, 0, 7, 9, 6, 7, 11, 7, 9, 11, -1 },  /* 181: criegler */
		{ 2, 3, 10, 0, 1, 8, 6, 7, 11, 1, 8, 7, 1, 11, 7, -1 },  /* 182: cworkman */
		{ 1, 2, 10, 1, 6, 11, 1, 6, 7, 1, 7, 10, -1, -1, -1, -1 },  /* 183: hang */
		{ 6, 7, 8, 1, 2, 6, 1, 6, 8, 1, 8, 9, -1, -1, -1, -1 },  /* 184: liuly */
		{ 1, 2, 9, 0, 3, 9, 3, 7, 9, 6, 7, 9, 2, 6, 9, -1 },  /* 185: gmorriso */
		{ 0, 8, 7, 6, 7, 0, 2, 0, 6, -1, -1, -1, -1, -1, -1, -1 },  /* 186: jbayes */
		{ 7, 3, 2, 7, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 187: jbrawner */
		{ 1, 3, 6, 1, 6, 9, 6, 8, 9, 6, 7, 8, 3, 6, 10 , -1 },  /* 188: jhorn3 */
		{ 0, 1, 9, 6, 7, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 189: jlowen */
		{ 0, 3, 10, 0, 6, 7, 0, 6, 10, 0, 8, 7, -1, -1, -1, -1 } ,  /* 190: jnelson */
		{ 6, 7, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 191: kdawes */
		{ 11, 5, 10, 5, 7, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 192: kpinto */
		{ 0, 3, 8, 5, 7, 10, 5, 10, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 193: mcmillan */
		{ 0, 1, 9, 5, 7, 10, 7, 10, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 194: nivaldot */
		{ 1, 3, 8, 1, 8, 9, 5, 7, 11, 7, 10, 11, -1, -1, -1, -1 },  /* 195: sgrady2 */
		{ 3, 5, 7, 2, 3, 5, 2, 5, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 196: touermi */
		{ 0, 2, 8, 2, 5, 8, 5, 7, 8, 2, 5, 11, -1, -1, -1, -1 },  /* 197: tristanj */
		{ 0, 1, 9, 2, 3, 10, 7, 6, 10, 5, 6, 11, -1, -1, -1, -1 },  /* 198: aphilli9 */
		{ 2, 5, 11, 1, 2, 9, 2, 5, 7, 2, 7, 8, 2, 8, 9, -1 },  /* 199: chatfiel */
		{ 1,  5,  7,  1,  2, 10,  7, 10,  1, -1, -1, -1, -1, -1, -1, -1 },  /* 200: chesshir */
		{ 0, 3, 8, 1, 5, 7, 1, 2, 7, 2, 7, 10, -1, -1, -1, -1 },  /* 201: cneugass */
		{ 0,  2,  9,  2,  7,  9,  2,  7, 10,  5,  7,  9, -1, -1, -1, -1 },  /* 202: cpalk */
		{ 2, 3, 8, 2, 5, 9, 2, 5, 7, 2, 7, 10, 2, 8, 9, -1 },  /* 203: criegler */
		{ 5, 7, 3, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 204: cworkman */
		{ 1, 5, 7, 0, 7, 8, 0, 1, 7, -1, -1, -1, -1, -1, -1, -1 },  /* 205: ssane */
		{ 3,  5,  7,  0,  3,  5,  0,  5,  9, -1, -1, -1, -1, -1, -1, -1 },  /* 206: bishara */
		{ 5, 7, 9, 7, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 207: gmorriso */
		{ 8, 10, 11, 8, 5, 11, 8, 5, 4, -1, -1, -1, -1, -1, -1, -1 },  /* 208: jbayes */
		{ 3, 0, 10, 5, 11, 10, 4, 5, 0, 10, 0, 5, -1, -1, -1, -1 },  /* 209: jbrawner */
		{ 0, 1, 9, 4, 5, 11, 4, 8, 11, 8, 10, 11, -1, -1, -1, -1 },  /* 210: jhorn3 */
		{ 4, 5, 9, 1, 3, 11, 3, 10, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 211: jlowen */
		{ 2, 3, 5, 2, 5, 11, 3, 4, 5, 3, 4, 8,  -1, -1, -1, -1 } ,  /* 212: jnelson */
		{ 0, 2, 4, 2, 4, 5, 2, 5, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 213: kdawes */
		{ 3, 5, 8, 2, 3, 11, 0, 1, 9, 4, 5, 8, 3, 5, 11, -1 },  /* 214: kpinto */
		{ 1, 2, 11, 1, 5, 11, 1, 5, 9, 4, 5, 9, -1, -1, -1, -1 },  /* 215: mcmillan */
		{ 1, 2, 5, 2, 5, 11, 4, 5, 11, 4, 8, 11, -1, -1, -1, -1 },  /* 216: nivaldot */
		{ 0, 3, 10, 0, 4, 10, 1, 2, 10, 1, 5, 10, 4, 5, 10, -1 },  /* 217: sgrady2 */
		{ 2, 5, 10, 0, 5, 9, 0, 2, 5, 4, 5, 8, 5, 8, 10, -1 },  /* 218: touermi */
		{ 4, 5, 9, 2, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 219: tristanj */
		{ 1, 3, 5, 4, 5, 8, 3, 8, 5, -1, -1, -1, -1, -1, -1, -1 },  /*220: James Kress (3) */
		{ 0, 1, 5, 0, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 221: chatfiel */
		{ 0, 5, 9, 0, 3, 5, 4, 5, 8, 3, 5, 8, -1, -1, -1, -1 },  /* 222: hang */
		{ 4, 5, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 223: cneugass */
		{ 4,  7, 10,  9, 10, 11,  4,  9, 10, -1, -1, -1, -1, -1, -1, -1 },  /* 224: cpalk */
		{ 0, 3, 8, 9, 10, 11, 4, 7, 9, 7, 9, 10, -1, -1, -1, -1 },  /* 225: criegler */
		{ 11, 10, 1, 1, 10, 4, 1, 0, 4, 4, 7, 10, -1, -1, -1, -1 },  /* 226: cworkman */
		{ 4,  7,  8,  1,  3, 10,  1, 10, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 227: bishara */
		{ 3, 4, 7, 2, 3, 4, 2, 4, 9, 2, 9, 11, -1, -1, -1, -1 },  /* 228: josh */
		{ 7, 9, 11, 0, 7, 8, 0, 2, 7, 4, 7, 9, 2, 7, 11, -1 },  /* 229: gmorriso */
		{ 11, 2, 3, 11, 1, 0, 11, 7, 3, 11, 4, 7, 11, 4, 0, -1 },  /* 230: jbayes */
		{ 2, 1, 11, 7, 4, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 231: jbrawner */
		{ 1, 2, 9, 4, 7, 10, 4, 9, 10, 2, 9, 10, -1, -1, -1, -1 },  /* 232: jhorn3 */
		{ 0, 1, 9, 2, 3, 10, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1 },  /* 233: jlowen */
		{ 0, 2, 4, 2, 4, 10, 4, 7, 10, -1, -1, -1, -1, -1, -1, -1 } ,  /* 234: jnelson */
		{ 2, 4, 10, 2, 3, 4, 3, 4, 8, 4, 7, 10, -1, -1, -1, -1 },  /* 235: kdawes */
		{ 1, 4, 9, 1, 3, 7, 1, 4, 7, -1, -1, -1, -1, -1, -1, -1 },  /* 236: kpinto */
		{ 0, 1, 9, 0, 4, 9, 0, 4, 8, 4, 7, 8, -1, -1, -1, -1 },  /* 237: mcmillan */
		{ 0, 3, 4, 3, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 238: nivaldot */
		{ 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 239: sgrady2 */
		{ 8, 9, 11, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 240: touermi */
		{ 0, 3, 9, 3, 9, 10, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 241: tristanj */
		{ 0,  1, 11,  0, 11, 10,  0, 10,  8, -1, -1, -1, -1, -1, -1, -1 },  /*242: Brandon Hildreth (3)  1111 0010*/
		{ 1, 3, 11, 3, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 243: chatfiel */
		{ 8,  9, 11,  2,  3,  8,  2,  8, 11, -1, -1, -1, -1, -1, -1, -1 },  /* 244: chesshir */
		{ 2, 9, 11, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 245: cneugass */
		{ 0,  1,  8,  2,  3,  8,  1,  8, 11,  2,  8, 11, -1, -1, -1, -1 },  /* 246: cpalk */
		{ 1, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 247: criegler */
		{ 8, 9, 10, 9, 1, 10, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1 },  /* 248: cworkman */
		{ 0, 1, 9, 2, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 249: josh */
		{ 0, 2, 10, 0, 8, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /*250: Paul Elliott (2) */
		{ 2, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 251: gmorriso */
		{ 1, 3, 8, 1, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 252: jbayes */
		{ 0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 253: jbrawner */
		{ 0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 254: pem */
		{ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },  /* 255: jlowen */
	};

	// Determines how many triangles per case and stores 
	// the info in an array
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 16; j++) {
			if (lookup[i][j] == -1) {
				num_triangles[i] = j / 3;
				break;
			}
		}
	}

	// Runs through every cell
	for (int i = 0; i < no_cells; i++) {

		// Initializes case ID
		case_id = 0;

		// Determines the logical cell indices based on
		// current cell number
		GetLogicalCellIndex(v0_log, i, dims);

		v1_log[0] = v0_log[0] + 1;
		v1_log[1] = v0_log[1];
		v1_log[2] = v0_log[2];

		v2_log[0] = v0_log[0];
		v2_log[1] = v0_log[1] + 1;
		v2_log[2] = v0_log[2];

		v3_log[0] = v0_log[0] + 1;
		v3_log[1] = v0_log[1] + 1;
		v3_log[2] = v0_log[2];

		v4_log[0] = v0_log[0];
		v4_log[1] = v0_log[1];
		v4_log[2] = v0_log[2] + 1;

		v5_log[0] = v0_log[0] + 1;
		v5_log[1] = v0_log[1];
		v5_log[2] = v0_log[2] + 1;

		v6_log[0] = v0_log[0];
		v6_log[1] = v0_log[1] + 1;
		v6_log[2] = v0_log[2] + 1;

		v7_log[0] = v0_log[0] + 1;
		v7_log[1] = v0_log[1] + 1;
		v7_log[2] = v0_log[2] + 1;

		// Grabs the field values at every index
		v0_f = F[GetPointIndex(v0_log, dims)];
		v1_f = F[GetPointIndex(v1_log, dims)];
		v2_f = F[GetPointIndex(v2_log, dims)];
		v3_f = F[GetPointIndex(v3_log, dims)];
		v4_f = F[GetPointIndex(v4_log, dims)];
		v5_f = F[GetPointIndex(v5_log, dims)];
		v6_f = F[GetPointIndex(v6_log, dims)];
		v7_f = F[GetPointIndex(v7_log, dims)];

		// Compares field values to isovalue to determine
		// which case to use
		v0_f < ISOVALUE ? v0_bin = 0 : v0_bin = 1;
		v1_f < ISOVALUE ? v1_bin = 0 : v1_bin = 1;
		v2_f < ISOVALUE ? v2_bin = 0 : v2_bin = 1;
		v3_f < ISOVALUE ? v3_bin = 0 : v3_bin = 1;
		v4_f < ISOVALUE ? v4_bin = 0 : v4_bin = 1;
		v5_f < ISOVALUE ? v5_bin = 0 : v5_bin = 1;
		v6_f < ISOVALUE ? v6_bin = 0 : v6_bin = 1;
		v7_f < ISOVALUE ? v7_bin = 0 : v7_bin = 1;
		
		// "Binary" to find case number
		v0_bin == 1 ? case_id += 1 : case_id += 0;
		v1_bin == 1 ? case_id += 2 : case_id += 0;
		v2_bin == 1 ? case_id += 4 : case_id += 0;
		v3_bin == 1 ? case_id += 8 : case_id += 0;
		v4_bin == 1 ? case_id += 16 : case_id += 0;
		v5_bin == 1 ? case_id += 32 : case_id += 0;
		v6_bin == 1 ? case_id += 64 : case_id += 0;
		v7_bin == 1 ? case_id += 128 : case_id += 0;

		// Gets the number of triangles to draw
		no_triangles = num_triangles[case_id];

		// For each triangle in a specified case...
		for (int j = 0; j < no_triangles; j++) {

			int edge[3];
			float pt[3][3];

			// Grabs the edge number based on the case number
			edge[0] = lookup[case_id][3 * j];
			edge[1] = lookup[case_id][(3 * j) + 1];
			edge[2] = lookup[case_id][(3 * j) + 2];

			// Determines points based on edges and interpolation
			for (int k = 0; k < 3; k++) {
				if (edge[k] == 0) {
					pt[k][0] = Interpolation(ISOVALUE, v0_f, v1_f, X[v0_log[0]], X[v1_log[0]]);
					pt[k][1] = Y[v0_log[1]];
					pt[k][2] = Z[v0_log[2]];
				} else if (edge[k] == 1) {
					pt[k][0] = X[v1_log[0]];
					pt[k][1] = Interpolation(ISOVALUE, v1_f, v3_f, Y[v1_log[1]], Y[v3_log[1]]);
					pt[k][2] = Z[v1_log[2]];
				} else if (edge[k] == 2) {
					pt[k][0] = Interpolation(ISOVALUE, v2_f, v3_f, X[v2_log[0]], X[v3_log[0]]);
					pt[k][1] = Y[v2_log[1]];
					pt[k][2] = Z[v2_log[2]];
				} else if (edge[k] == 3) {
					pt[k][0] = X[v0_log[0]];
					pt[k][1] = Interpolation(ISOVALUE, v0_f, v2_f, Y[v0_log[1]], Y[v2_log[1]]);
					pt[k][2] = Z[v0_log[2]];
				} else if (edge[k] == 4) {
					pt[k][0] = Interpolation(ISOVALUE, v4_f, v5_f, X[v4_log[0]], X[v5_log[0]]);
					pt[k][1] = Y[v4_log[1]];
					pt[k][2] = Z[v4_log[2]];
				} else if (edge[k] == 5) {
					pt[k][0] = X[v5_log[0]];
					pt[k][1] = Interpolation(ISOVALUE, v5_f, v7_f, Y[v5_log[1]], Y[v7_log[1]]);
					pt[k][2] = Z[v5_log[2]];
				} else if (edge[k] == 6) {
					pt[k][0] = Interpolation(ISOVALUE, v6_f, v7_f, X[v6_log[0]], X[v7_log[0]]);
					pt[k][1] = Y[v6_log[1]];
					pt[k][2] = Z[v6_log[2]];
				} else if (edge[k] == 7) {
					pt[k][0] = X[v4_log[0]];
					pt[k][1] = Interpolation(ISOVALUE, v4_f, v6_f, Y[v4_log[1]], Y[v6_log[1]]);
					pt[k][2] = Z[v4_log[2]];
				} else if (edge[k] == 8) {
					pt[k][0] = X[v0_log[0]];
					pt[k][1] = Y[v0_log[1]];
					pt[k][2] = Interpolation(ISOVALUE, v0_f, v4_f, Z[v0_log[2]], Z[v4_log[2]]);
				} else if (edge[k] == 9) {
					pt[k][0] = X[v1_log[0]];
					pt[k][1] = Y[v1_log[1]];
					pt[k][2] = Interpolation(ISOVALUE, v1_f, v5_f, Z[v1_log[2]], Z[v5_log[2]]);
				} else if (edge[k] == 10) {
					pt[k][0] = X[v2_log[0]];
					pt[k][1] = Y[v2_log[1]];
					pt[k][2] = Interpolation(ISOVALUE, v2_f, v6_f, Z[v2_log[2]], Z[v6_log[2]]);
				} else if (edge[k] == 11) {
					pt[k][0] = X[v3_log[0]];
					pt[k][1] = Y[v3_log[1]];
					pt[k][2] = Interpolation(ISOVALUE, v3_f, v7_f, Z[v3_log[2]], Z[v7_log[2]]);
				}
			}

			// Adds triangle to be drawn
			tl.AddTriangle(pt[0][0], pt[0][1], pt[0][2], pt[1][0], pt[1][1], pt[1][2], pt[2][0], pt[2][1], pt[2][2]);
		}
	}

	vtkPolyData * pd = tl.MakePolyData();

	// ** End of my code **

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
