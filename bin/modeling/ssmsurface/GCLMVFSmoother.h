class GCLMVFSmoother
{
    public:
        GCLMVFSmoother();
        /* \brief  Define coefficients with a dmatrix.

           This defines a rectangular fir filter for the grid
           input with a dmatrix.  i will run on x1 axis while
           j will run over x2 axis.  
           \param firfilter contains the fir coefficients (normalization
             is not required)
           \param i0 offset in i (x1) direction for 0 lag.  C indexing.
           \param j0 offset in j (x2) direction for 0 lag.  C indexing.
           */
        GCLMVFSmoother(dmatrix& firfilter,int i0in, int j0in);
        GCLMaskedScalarField apply(GCLMaskedScalarField& parent);
    private:
        dmatrix firfilter;
        int nrow,ncol;
        int i0,j0;  //Offset of 0 value
};
