namespace mymex {

    double* set_field_double(mxArray* my_struct, const char *name, size_t ndim, size_t *dims){
        
        auto var = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL); 
        mxSetField(my_struct,0,name,var);
        return (double*) mxGetPr(var);    

    }

    size_t* set_field_int(mxArray* my_struct, const char *name, size_t ndim, size_t* dims){

        auto var = mxCreateNumericArray(ndim,dims,mxINT32_CLASS,mxREAL); 
        mxSetField(my_struct,0,name,var);
        return (size_t*) mxGetData(var);    

    }

    double** set_field_cell(mxArray * my_struct, const char *name, 
                            size_t ndim_cell, size_t *dims_cell, 
                            size_t *ndim, size_t **dims){
        
        // a. field
        auto cell = mxCreateCellArray(ndim_cell, dims_cell);
        mxSetField(my_struct,0,name,cell);

        // b. allocate memory
        auto out = new double*[dims_cell[0]];
        for(size_t t = 0; t < dims_cell[0]; t++){

            auto array  = mxCreateNumericArray(ndim[0],dims[0],mxDOUBLE_CLASS,mxREAL);
            mxSetCell(cell, t, array);
            out[t] = mxGetPr(array);

        }

        return out;

    }

}