#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <string>
#include <iostream>
#include "gclgrid.h"

using namespace std;

vtkStructuredGrid*	convert_gcl3d_to_vtksg(GCLscalarfield3d &g, bool swap_xz = true)
{

	vtkStructuredGrid	*pGrid 			= vtkStructuredGrid::New();
	vtkPoints		*pPoints		= vtkPoints::New();
	vtkDoubleArray		*pData			= vtkDoubleArray::New();

/*
	if(swap_xz)	pGrid->SetDimensions(g.n3, g.n2, g.n1);
	else 		pGrid->SetDimensions(g.n1, g.n2, g.n3);
*/
	pGrid->SetDimensions(g.n1, g.n2, g.n3);

	pData->SetNumberOfTuples(g.n1 * g.n2 * g.n3);
	pData->SetNumberOfComponents(1);

/*
	for(int i = 0; i < g.n1; i++)
	{
		for(int j = 0; j < g.n2; j++)
		{
			for(int k = 0; k < g.n3; k++)
			{
				// insert the 'geometry'
				if(swap_xz)	pPoints->InsertNextPoint(g.x3[i][j][k], g.x2[i][j][k], g.x1[i][j][k]);
				else		pPoints->InsertNextPoint(g.x1[i][j][k], g.x2[i][j][k], g.x3[i][j][k]);

				// insert the 'data'
				// corrupted.  not the original
				pData->SetValue((k * g.n2 * g.n3) + (j * g.n3) + k, g.val[i][j][k]);
				
			}
		}
	}
*/
	int kk;
	for(int k=g.n3-1,kk=0;k>=0;k--,kk++)
	{
		for(int j = 0; j < g.n2; j++)
		{
			for(int i = 0; i < g.n1; i++)
			{
				pPoints->InsertNextPoint(g.x1[i][j][k], 
					g.x2[i][j][k], g.x3[i][j][k]);

				// insert the 'data'
				//pData->SetValue((i * g.n2 * g.n3) + (j * g.n3) + k, g.val[i][j][k]);
				pData->SetValue((kk*g.n1*g.n2)+(j*g.n1)+i,
					g.val[i][j][k]);
				
			}
		}
	}

	pGrid->SetPoints(pPoints);
	pGrid->GetPointData()->SetScalars(pData);

	pPoints->Delete();
	pData->Delete();

	return(pGrid);
}

template<typename W>
void write_vtksg(vtkStructuredGrid *pGrid, const char *filename)
{
	W*	pWriter = W::New();

	pWriter->SetFileName(filename);
	pWriter->SetInput(pGrid);
	pWriter->Write();
	pWriter->Delete();

}

void output_gcl3d_to_vtksg(GCLscalarfield3d &g, const string filename,
							bool use_xml = false, bool swap_xz = true)
{
	
	vtkStructuredGrid	*pGrid = 0; 	

	pGrid = convert_gcl3d_to_vtksg(g, swap_xz);

	if(use_xml)	write_vtksg<vtkXMLStructuredGridWriter>(pGrid, 
				filename.c_str());
	else		write_vtksg<vtkStructuredGridWriter>(pGrid, 
				filename.c_str());

	pGrid->Delete();
}



