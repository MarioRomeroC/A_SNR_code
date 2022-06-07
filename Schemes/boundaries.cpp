#ifndef BOUNDARIES_CPP
#define BOUNDARIES_CPP

/**
    Here you will find all functions related to ghost cells and boundary conditions
**/

///RIGHT BOUNDARY
//These are expected to work behind the curtains
void reflect_right(Cell& last,Cell& penult, Cell* grid, const int N=NCELLS){
    /***
        This function will create a reflection at the right boundary of grid (last values).

        All magnitudes are the same, except momentum which changes sign
        | grid[N-1] |   last    |
        |   <---    |   --->    |

        This is generally a good standard condition, but is not suitable for cosmological simulations
    ***/

    Cell* true_last;
    Cell* true_penult;

    Cell* finder;
    finder = grid;
    true_last = find_last(finder);
    true_penult = true_last->prev;

    ///Start computing the correct positions and widths
    my_float r_last  = true_last->location() + true_last->width();
    my_float dr_last = true_last->width();

    my_float r_penult  = true_penult->location() + 2.*dr_last + true_penult->width();
    my_float dr_penult = true_penult->width();

    ///End copying magnitudes to last and penult, except momentum
    my_float l_data[NUM_EQUATIONS];
    my_float p_data[NUM_EQUATIONS];

    for(int elem=0;elem<NUM_EQUATIONS;elem++){
        l_data[elem] = true_last->eq(elem);
        p_data[elem] = true_penult->eq(elem);
    }
    //Reverse momentums
    l_data[2] = - l_data[2];
    p_data[2] = - p_data[2];

    last = Cell(l_data,true_last->get_pressure(),true_last->get_temperature(),true_last->get_adiabatic_const(),
                true_last->get_metallicity(),r_last - 0.5*dr_last,r_last + 0.5*dr_last,true_last->get_system());

    penult = Cell(p_data,true_penult->get_pressure(),true_penult->get_temperature(),true_penult->get_adiabatic_const(),
                true_penult->get_metallicity(),r_penult - 0.5*dr_penult,r_penult + 0.5*dr_penult,true_penult->get_system());

    last.update_tcool(true_last->get_cooling_time());
    last.update_meanweight(true_last->get_meanweight());

    penult.update_tcool(true_penult->get_cooling_time());
    penult.update_meanweight(true_penult->get_meanweight());

    ///That's all
}

void copy_right(Cell& last,Cell& penult, Cell* grid, const int N=NCELLS){
    /***
        This function will create a copy at the right boundary of grid (last values).

        All magnitudes are the same, even momentum
        | grid[N-1] |   last    |
        |   <---    |   <---    |

        This is useful for Noh test problem, in which you have a radially inward velocity
    ***/
    //I bet that this function is more inefficient than 'reflect_right', but is far clearer

    Cell* true_last;
    Cell* true_penult;

    Cell* finder;
    finder = grid;

    true_last = find_last(finder);
    true_penult = true_last->prev;

    ///Start copying magnitudes to last and penult
    last   = *true_last;
    penult = *true_penult;

    ///Now compute the correct positions
    my_float r_last   = last.location() + last.width();
    my_float r_penult = penult.location() + 2.*last.width() + penult.width();

    ///Now rellocate last and penult
    last   = last.move_location(r_last);
    penult = penult.move_location(r_penult);

    last.update_tcool(true_last->get_cooling_time());
    last.update_meanweight(true_last->get_meanweight());

    penult.update_tcool(true_penult->get_cooling_time());
    penult.update_meanweight(true_penult->get_meanweight());

    ///That's all

}

///LEFT BOUNDARY
//These are expected to work behind the curtains
void reflect_left(Cell& first,Cell& second, Cell* grid, const int N=NCELLS){
    /***
        This function will create a reflection at the right boundary of grid (last values).

        All magnitudes are the same, except momentum which changes sign
        |   first   |  grid[0]  |
        |   <---    |   --->    |

        This is generally a good standard condition, but is not suitable for cosmological simulations
    ***/

    Cell* true_first;
    Cell* true_second;

    Cell* finder;
    finder = grid;
    true_first = find_first(finder);
    true_second = true_first->next;

    ///Start computing the correct positions and widths
    my_float r_first  = true_first->location() - true_first->width();
    my_float dr_first = true_first->width();

    my_float r_second  = true_second->location() - 2.*dr_first - true_second->width();
    my_float dr_second = true_second->width();

    ///End copying magnitudes to last and penult, except momentum
    my_float f_data[NUM_EQUATIONS];
    my_float s_data[NUM_EQUATIONS];

    for(int elem=0;elem<NUM_EQUATIONS;elem++){
        f_data[elem] = true_first->eq(elem);
        s_data[elem] = true_second->eq(elem);
    }
    //Reverse momentums
    f_data[2] = -f_data[2];
    s_data[2] = -s_data[2];

    first = Cell(f_data,true_first->get_pressure(),true_first->get_temperature(),true_first->get_adiabatic_const(),
                true_first->get_metallicity(),r_first - 0.5*dr_first,r_first + 0.5*dr_first,true_first->get_system());

    second = Cell(s_data,true_second->get_pressure(),true_second->get_temperature(),true_second->get_adiabatic_const(),
                true_second->get_metallicity(),r_second - 0.5*dr_second,r_second + 0.5*dr_second,true_second->get_system());

    first.update_tcool(true_first->get_cooling_time());
    first.update_meanweight(true_first->get_meanweight());

    second.update_tcool(true_second->get_cooling_time());
    second.update_meanweight(true_second->get_meanweight());

    ///That's all

}

void copy_left(Cell& first,Cell& second, Cell* grid, const int N=NCELLS){
    /***
        This function will create a copy at the right boundary of grid (last values).

        All magnitudes are the same, even momentum
        |   first   |  grid[0]  |
        |   <---    |   <---    |

        Can be an interesting condition for spherical coordinates if the center is a sink?
        Can't tell!
    ***/
    //I bet that this function is more inefficient than 'reflect_left', but is far clearer

    Cell* true_first;
    Cell* true_second;


    Cell* finder;
    finder = grid;
    true_first = find_first(finder);
    true_second = true_first->next;

    ///Start copying magnitudes to last and penult
    first  = *true_first;//grid[0];
    second = *true_second;//grid[1];

    ///Now compute the correct positions
    my_float r_first  = first.location() - first.width();
    my_float r_second = second.location() - 2.*first.width() - second.width();

    ///Now rellocate last and penult
    first  = first.move_location(r_first);
    second = second.move_location(r_second);

    first.update_tcool(true_first->get_cooling_time());
    first.update_meanweight(true_first->get_meanweight());

    second.update_tcool(true_second->get_cooling_time());
    second.update_meanweight(true_second->get_meanweight());

    ///That's all

}

///MAIN BOUNDARY CONDITIONS ROUTINE

void create_boundary(Cell& first, Cell& second, Cell& penult, Cell& last, Cell* grid, const int N){
    /***
        Self-explanatory. Creating boundary conditions according to parameters given in 'define.h'
    ***/

    ///Create right boundary condition

    switch(RBOUNDARY){

        case 0: //Make a reflection
            reflect_right(last,penult,grid,N);
            break;
        case 1: //Make a copy
            copy_right(last,penult,grid,N);
            break;
        default:
            throw std::runtime_error("Boundary condition not understood");
    }

    ///Create left boundary condition
    //Careful with cylindrical/spherical coordinates, it should be a reflection

    //Same as right boundary
    switch(LBOUNDARY){

        case 0: //Make a reflection
            reflect_left(first,second,grid,N);
            break;
        case 1: //Make a copy
            if(SYMMETRY > 0){
                std::cout<<"Warning: For cylindrical or spherical coordinates you should use a reflection for the left boundary..."<<endl;
            }
            copy_left(first,second,grid,N);
            break;
        default:
            throw std::runtime_error("Boundary condition not understood");
    }

}

#endif
