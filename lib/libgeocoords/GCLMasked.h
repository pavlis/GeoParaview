#ifndef _GCLMASK_H_
#define _GCLMASK_H_
#include "dmatrix.h"
#include "gclgrid.h"
/*! Undefined values in a field are set to this magic value */
const double GCLFieldNullValue(-9.99999e-99);
/*! \brief  Used to apply a mask to a GCLgrid object.

  GCLgrid objects have a matrix structure that limits their natural
  extent to a distorted square in 3 space.  A mask is one way to 
  build a more complex shape.   Grid cells are defined as either valid
  or invalid (true or false respectively).   
  */
class GCLMask
{
    public:
        /*! Default constructor.  zero initializer */
        GCLMask();
        /*! \brief Build an initial mask based on a pattern grid.

          A GCLMask object is used to turn components on or off of a GCLgrod
          or one of the field objects derived from it.  This constructor
          builds a pattern based on g and initializes all cells to
          invalid.   Use enable_point method to turn points back on. */
        GCLMask(GCLgrid& g);
        /*! Standard copy constructor. */
        GCLMask(const GCLMask& parent);
        GCLMask& operator=(const GCLMask& parent);
        /*! \brief  Define point i,j of the grid as masked. 

           When building a mask point by point this is a core method.
           \param i is x1 index of point to mask
           \param j is x2 index of point to mask
           \return true on success.  false if failed (outside range implied)
           */
        bool mask_point(int i, int j);
        /*! \brief  Define point i,j of the grid as valid.

           When building a mask point by point this is a core method.
           This is the inverse of mask_point as it makes a point valid.
           \param i is x1 index of point to mask
           \param j is x2 index of point to mask
           \return true on success.  false if failed (outside range implied)
           */
        bool enable_point(int i, int j);
        /*! \brief Query if a point is marked as visible.

           This is a core method.  Returns true if the point at grid
           position i,j is marked as visible(valid).  Note if the
           index is outside the grid silently returns false since
           by definition such a point is invisible. 
           \param i is x1 index of point to query
           \param j is x2 index of point to query

           */
        bool point_is_valid(int i, int j);
    protected:
        /* Created as a single array, easily indexed 
           because GCLgrid is n1xn2 */
        vector<bool> valid;
        /* These are the sizes used to create mask.  Constructors check these for 
           consistency when applied to a parent GCLgrid */
        int n1_mask, n2_mask;
    private:
        /* convenience routine */
        int voffset(int i, int j)
        {
            int result;
            result=j*n1_mask+i;
            return(result);
        }
};
class GCLMaskedGrid : public GCLgrid, public GCLMask
{
    public:
        /* \brief Default constructor.  */
        GCLMaskedGrid(){};
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
        GCLMaskedGrid& operator=(const GCLMaskedGrid& parent);
};
/* Variant of GCLscalarfield but with a smaller interface */
class GCLMaskedScalarField : public GCLscalarfield, public GCLMask
{
    public:
        GCLMaskedScalarField(GCLMaskedGrid& g);
        GCLMaskedScalarField(const GCLMaskedScalarField& parent);
        ~GCLMaskedScalarField();
        GCLMaskedScalarField& operator=(const GCLMaskedScalarField& parent);
};
/* This is a variant of the GCLvectorfield but with a smaller
   interface.  For now used only for slabmodel so unless I release 
   the newest code this will not be a major issue. */
class GCLMaskedVectorField : public GCLvectorfield, public GCLMask
{
    public:
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
#endif
