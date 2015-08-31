#include "GCLgridError.h"
#include "gclgrid.h"
#include "dmatrix.h"
#include "GCLMasked.h"
/* First in this file - GCLMask code. */
GCLMask::GCLMask()
{
    n1_mask=0;
    n2_mask=0;
}
GCLMask::GCLMask(GCLgrid& parent)
{
    n1_mask=parent.n1;
    n2_mask=parent.n2;
    int ntotal=parent.n1*parent.n2;
    valid.reserve(ntotal);
    int k;
    for(k=0;k<ntotal;++k)valid.push_back(false);
}
GCLMask::GCLMask(const GCLMask& parent) : valid(parent.valid)
{
    n1_mask=parent.n1_mask;
    n2_mask=parent.n2_mask;
}
GCLMask& GCLMask::operator=(const GCLMask& parent)
{
    if(this != & parent)
    {
        this->valid=parent.valid;
        n1_mask=parent.n1_mask;
        n2_mask=parent.n2_mask;
    }
    return(*this);
}
bool GCLMask::mask_point(int i, int j)
{
    if(i<0 || (i>=n1_mask)) return false;
    if(j<0 || (j>=n2_mask)) return false;
    valid[this->voffset(i,j)] = false;
    return(true);
}
bool GCLMask::enable_point(int i, int j)
{
    if(i<0 || (i>=n1_mask)) return false;
    if(j<0 || (j>=n2_mask)) return false;
    valid[this->voffset(i,j)] = true;
    return(true);
}
bool GCLMask::point_is_valid(int i, int j)
{
    if(i<0 || (i>=n1_mask)) return false;
    if(j<0 || (j>=n2_mask)) return false;
    int k=voffset(i,j);
    return(valid[k]);
}



GCLMaskedGrid::GCLMaskedGrid(GCLgrid& parent) : GCLgrid(parent), GCLMask(parent)
{
}
GCLMaskedGrid::GCLMaskedGrid(GCLgrid& parent, dmatrix& mask, double maskmin)
                        : GCLgrid(parent)
{
    const string base_error("GCLMaskedGrid constructor with dmatrix mask:  ");
    if(mask.columns()!=n2) throw GCLgridError(base_error
            + "number of columns in mask array does not match grid");
    if(mask.rows()!=n1) throw GCLgridError(base_error
            + "number of columns in mask array does not match grid");
    int i,j;
    for(i=0;i<parent.n1;++i)
        for(j=0;j<parent.n2;++j)
        {
            if(mask(i,j)>maskmin)
                this->mask_point(i,j);
            else
                this->enable_point(i,j);
        }
}
GCLMaskedGrid::GCLMaskedGrid(GCLMaskedGrid& parent) 
    : GCLgrid(dynamic_cast<GCLgrid&>(parent)), GCLMask(dynamic_cast<GCLMask&>(parent))
{
}
GCLMaskedGrid& GCLMaskedGrid::operator=(const GCLMaskedGrid& parent)
{
    if(this != & parent)
    {
        this->GCLgrid::operator=(parent);
        this->GCLMask::operator=(parent);
    }
    return(*this);
}
GCLMaskedVectorField::GCLMaskedVectorField(GCLMaskedGrid& g, int nvsize)
    : GCLvectorfield(dynamic_cast<GCLgrid&>(g),nvsize), GCLMask(dynamic_cast<GCLMask&>(g))
{
}
GCLMaskedVectorField::GCLMaskedVectorField(const GCLMaskedVectorField& parent)
    : GCLvectorfield(dynamic_cast<const GCLvectorfield&>(parent)),
        GCLMask(dynamic_cast<const GCLMask&>(parent))
{
}
GCLMaskedVectorField& GCLMaskedVectorField::operator=(const GCLMaskedVectorField& parent)
{
    if(this != & parent)
    {
       this->GCLvectorfield::operator=(parent);
       this->GCLMask::operator=(parent);
    }
    return(*this);
}
GCLMaskedScalarField::GCLMaskedScalarField(GCLMaskedGrid& g) 
    : GCLscalarfield(dynamic_cast<GCLgrid&>(g)), GCLMask(dynamic_cast<GCLMask&>(g))
{
}
GCLMaskedScalarField::GCLMaskedScalarField(const GCLMaskedScalarField& parent)
    : GCLscalarfield(dynamic_cast<const GCLscalarfield&>(parent)), 
        GCLMask(dynamic_cast<const GCLMask&>(parent))
{
}
GCLMaskedScalarField& GCLMaskedScalarField::operator=(const GCLMaskedScalarField& parent)
{
    if(this != & parent)
    {
       this->GCLscalarfield::operator=(parent);
       this->GCLMask::operator=(parent);
    }
    return(*this);
}
