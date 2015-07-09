GCLMaskedGrid::GCLMaskedGrid() : GCLgrid()
{
    valid=NULL;
}
GCLMaskedGrid::GCLMaskedGrid(GCLgrid& parent) : GCLgrid(parent)
{
    int i;
    int ntotal=parent.n1*parent*n2;
    valid = new bool[ntotal];
    /* There is a way to do this with the new operator, but here
       I set all valid explicitly */
    for(i=0;i<ntotal;++i) valid[i]=true;
}
GCLMaskedGrid::GCLMaskedGrid(GCLgrid& parent, dmatrix& mask, double maskmin)
                        : GCLgrid(parent)
{
    const string base_error("GCLMaskedGrid constructor with dmatrix mask:  ");
    if(mask.columns()!=n2) throw GCLgridError(base_error
            + "number of columns in mask array does not match grid");
    if(mask.rows()!=n1) throw GCLgridError(base_error
            + "number of columns in mask array does not match grid");
    int ntotal=parent.n1*parent*n2;
    valid = new bool[ntotal];
    int i,j,k;
    for(i=0;i<parent.n1;++i)
        for(j=0;j<parent.n2;++j)
        {
            k=voffset(i,j);
            if(mask(i,j)>maskmin)
            {
                valid[k]=false;
            }
            else
                valid[k]=true;
        }
}
GCLMaskedGrid::GCLMaskedGrid(GCLMaskedGrid& parent) 
    : GCLgrid(dynamic_cast<GCLgrid&>(parent))
{
    int i;
    int ntotal=parent.n1*parent*n2;
    valid = new bool[ntotal];
    for(i=0;i<ntotal;++i) valid[i]=parent.valid[i];
}
GCLMaskedGrid::~GCLMaskedGrid()
{
    if(valid!=NULL) delete [] valid;
}
GCLMaskedGrid& GCLMaskedGrid::operator=(const GCLMaskedGrid& parent)
{
    if(this != & parent)
    {
        this->GCLgrid::operator=(parent);
        int ntotal=n1*n2;
        valid=new bool[ntotal];
        for(int k=0;k<ntotal;++i) valid=parent.valid[k];
    }
    return(*this);
}
bool GCLMaskedGrid::mask_point(int i, int j)
{
    if(i<0 || (i>=n1)) return false;
    if(j<0 || (j>=n2)) return false;
    valid[this->voffset(i,j)] = false;
    return(true);
}
bool GCLMaskedGrid::point_is_valid(int i, int j)
{
    if(i<0 || (i>=n1)) return false;
    if(j<0 || (j>=n2)) return false;
    int k=voffset(i,j);
    return(valid[k]);
}
GCLMaskedVectorField::GCLMaskedVectorField() : GCLgrid()
{
    nv=0;
    val=NULL;
}
GCLMaskedVectorField::GCLMaskedVectorField(GCLMaskedGrid& g, int nvsize)
{
    nv=nvsize;
    val=create_3dgrid_contiguous(g.n1,g.n2,nv);
    int i,j,k;
    for(i=0;i<n1;++i)
        for(j=0;j<n2;++j)
            for(k=0;k<nv;++k) val[i][j][k]=GCLFieldNullValue;
}
GCLMaskedVectorField::GCLMaskedVectorField(const GCLMaskedVectorField& parent)
    : GCLMaskedGrid(dynamic_cast<GCLMaskedGrid>(parent))
{
    nv=parent.nv;
    val=create_3dgrid_contiguous(n1,n2,nv);
    int i,j,k;
    for(i=0;i<n1;++i)
        for(j=0;j<n2;++j)
            for(k=0;k<nv;++k) val[i][j][k]=parent.val[i][j][k];
}
GCLMaskedVectorField::~GCLMaskedVectorField()
{
    if(val!=NULL) free_3dgrid_contiguous(val,n1,n2);
}
GCLMaskedVectorField& GCLMaskedVectorField::operator=(const GCLMaskedVectorField& parent)
{
    if(this != & parent)
    {
       this->GCLMaskedGrid::operator=(parent);
       val=create_3dgrid_contiguous(n1,n2,nv);
       int i,j,k;
       for(i=0;i<n1;++i)
           for(j=0;j<n2;++j)
               for(k=0;k<nv;++k) val[i][j][k]=parent.val[i][j][k];
    }
    return(*this);
}
GCLMaskedScalarField::GCLMaskedScalarField() : GCLgrid()
{
    nv=0;
    val=NULL;
}
GCLMaskedScalarField::GCLMaskedScalarField(GCLMaskedGrid& g)
{
    val=create_2dgrid_contiguous(g.n1,g.n2);
    int i,j;
    for(i=0;i<n1;++i)
        for(j=0;j<n2;++j) val[i][j]=GCLFieldNullValue;
}
GCLMaskedScalarField::GCLMaskedScalarField(const GCLMaskedScalarField& parent)
    : GCLMaskedGrid(dynamic_cast<GCLMaskedGrid>(parent))
{
    val=create_2dgrid_contiguous(n1,n2);
    int i,j;
    for(i=0;i<n1;++i)
        for(j=0;j<n2;++j)
            val[i][j]=parent.val[i][j];
}
GCLMaskedScalarField::~GCLMaskedScalarField()
{
    if(val!=NULL) free_2dgrid_contiguous(val,n1,n2);
}
GCLMaskedScalarField& GCLMaskedScalarField::operator=(const GCLMaskedScalarField& parent)
{
    if(this != & parent)
    {
       this->GCLMaskedGrid::operator=(parent);
       val=create_wdgrid_contiguous(n1,n2);
       int i,j;
       for(i=0;i<n1;++i)
           for(j=0;j<n2;++j) val[i][j]=parent.val[i][j];
    }
    return(*this);
}
