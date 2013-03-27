/*! Undefined values in a field are set to this magic value */
const double GCLFieldNullValue(-9.99999e-99);
class GCLMaskedGrid : public GCLgrid
{
    public:
        /* \brief Default constructor.  Minimal initializer */
        GCLMaskedGrid();
        /* \brief Construct with no original mask.

           Often we want to apply our own mask.  This constructor clones 
           g and defines all points to be visible.
           \param g parent grid to clone
           */
        GCLMaskedGrid(GCLgrid& g);
        /* \brief Construct with a preassembled mask defined by a matrix of doubles.

           \param parent is a grid that will be cloned.
           \param mask is a double array used as a mask 
           \param maskmin is a threshold used to avoid float ambiguity 
               of dmatrix values.  Any value in dmatrix greater than maskmin
               is used to say mask that point.   Normally set values to some
               sensible positive number like 1.0
           \exception GCLgridError is thrown if mask size is not consistent 
               with parent.
               */
        GCLMaskedGrid(GCLgrid& parent,dmatrix& mask,double maskmin=1.0e-13);
        /* \brief Standard copy constructor */
        GCLMaskedGrid(GCLMaskedGrid& parent);
        /* \brief standard destructor.

           Not trivial here as this has to deal with mask data which
           for this implementation requires freeing some dymamic memory.*/
        ~GCLMaskedGrid();
        /* \brief  Define point i,j of the grid as masked. 

           When building a mask point by point this is a core method.
           \param i is x1 index of point to mask
           \param j is x2 index of point to mask
           \return true on success.  false if failed (outside range implied)
           \exception GCLgridError thrown if index passed is outside grid range.
           */
        bool mask_point(int i, int j);
        /* \brief Query if a point is marked as visible.

           This is a core method.  Returns true if the point at grid
           position i,j is marked as visible(valid).  Note if the
           index is outside the grid silently returns false since
           by definition such a point is invisible. 
           \param i is x1 index of point to query
           \param j is x2 index of point to query

           */
        bool point_is_valid(int i, int j);
        GCLMaskedGrid& operator=(const GCLMaskedGrid& parent);
    private:
        /* Created as a single array, easily indexed 
           because GCLgrid is n1xn2 */
        bool *valid;
        /* convenience routine */
        int voffset(int i, int j)
        {
            int result;
            result=j*n1+i;
            return(result);
        }
};
/* Variant of GCLscalarfield but with a smaller interface */
class GCLMaskedScalarField : public GCLMaskedGrid
{
    public:
        int nv;
        GCLMaskedScalarField();
        GCLMaskedScalarField(GCLMaskedGrid& g)
        GCLMaskedScalarField(const GCLMaskedScalarField& parent);
        ~GCLMaskedScalarField();
        GCLMaskedScalarField& operator=(const GCLMaskedScalarField& parent);
};
/* This is a variant of the GCLvectorfield but with a smaller
   interface.  For now used only for slabmodel so unless I release 
   the newest code this will not be a major issue. */
class GCLMaskedVectorField : public GCLMaskedGrid
{
    public:
        int nv;
        GCLMaskedVectorField();
        GCLMaskedVectorField(GCLMaskedGrid& g, int nvsize);
        GCLMaskedVectorField(const GCLMaskedVectorField& parent);
        ~GCLMaskedVectorField();
        GCLMaskedVectorField& operator=(const GCLMaskedVectorField& parent);
};
class TiePoint
{
        public:
                    double radius;
                            int i_tie, j_tie;
};

