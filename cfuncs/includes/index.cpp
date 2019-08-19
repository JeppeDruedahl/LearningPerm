namespace index {

size_t d2(size_t i_x1, size_t i_x2, size_t Nx2){

    size_t i_grid =   i_x1*Nx2
                 + i_x2;

    return i_grid;

}

size_t d3(size_t i_x1, size_t i_x2, size_t i_x3, size_t Nx2, size_t Nx3){

    size_t i_grid =   i_x1*Nx2*Nx3
                 + i_x2*Nx3
                 + i_x3;

    return i_grid;

}

} // index