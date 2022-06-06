#ifndef AMR_MAIN_CPP
#define AMR_MAIN_CPP
#if AMR_ON == 1
/**
    Adaptive Mesh Refinement (AMR) main functions
**/

enum Interpols {Position, Density, Energy, Momentum, Entropy, Temperature, Pressure, Cooling, N_interpols }; //Cooling is cooling time.

Cell left_neighbour(my_float X, Cell* old_grid, const int curr_N){
    /*
        There are two cases (X1 and X2) to consider here
                        X1  X2
        |  j-2  |  j-1  |   j   |  j+1  |  j+2  |

        Condition for selecting neighbour is if(X < x_j) get x_{j-1}
        So in this sketch, X1 (in the middle of two cells) return j
        And X2 (Same position as x_j) return x_{j}
    */

    for(int j=0;j<curr_N;++j){
        if(X < old_grid[j].location()){ //Last cell will never fullfill this condition if X=old_grid[N-1].location !
            if(j != 0){ return old_grid[j-1];}
            else{ return old_grid[0];}
        }
    }

    if( X >= old_grid[curr_N-1].location()){
        return old_grid[curr_N-1];
    }else{
        std::cout<<old_grid[curr_N-1].location()<<" "<<old_grid[curr_N-2].location()<<" "<<X<<std::endl;
        throw std::runtime_error("Unable to find a left neighbour");
    }
}

Cell right_neighbour(my_float X, Cell* old_grid, const int curr_N){
    /*
        There are two cases (X1 and X2) to consider here
                        X1  X2
        |  j-2  |  j-1  |   j   |  j+1  |  j+2  |

        Condition for selecting neighbour is if(X >= x_j) get x_{j+1}
        So in this sketch, X1 (in the middle of two cells) return j
        And X2 (Same position as x_j) return x_{j+1} (You should expect x_j, but that equality condition returns next value)
        This is done in order to get an interval to interpolate. Interpolation does not care about that,
        however computers do (an interval of [x_j,x_j] returns a division by 0)
    */
    for(int j=curr_N-1;j>=0;--j){
        if(X > old_grid[j].location()){
            if(j != curr_N-1){ return old_grid[j+1];}
            else{ return old_grid[curr_N-1];}
        }
        else if( X == old_grid[j].location()){
            return old_grid[j+1];
        }
    }

    if( X <= old_grid[0].location()){
        return old_grid[0];
    }else{
        std::cout<<old_grid[0].location()<<" "<<old_grid[1].location()<<" "<<X<<std::endl;
        throw std::runtime_error("Unable to find a right neighbour");
    }

}

///Split-Merge routines
void Cell::split(Cell* grid, const int max_N){
    /*
        This routine splits the caller (this) into two Cells of equal widths.
        Caller (*this) becomes the LEFT cell
        this routine also find an unused cell in grid, becoming the RIGHT cell
    */

    //(1) Find an unused cell
    Cell* right = find_unused_cell(grid,max_N);
    //(2) Define auxilliary pointers
    Cell* next_cell;
    Cell* prev_cell;
    next_cell = this->next; //This link will be cut shortly, so better save a copy
    prev_cell = this->prev;
    Cell* ambient;
    ambient = find_last(grid);

    int old_flag_for_shock = this->flag_for_shock;

    Cell right_cell;
    if(next_cell != nullptr){
        right_cell = *next_cell;
    }else{
        right_cell = *this;
    }

    Cell left_cell;
    if(prev_cell != nullptr){
        left_cell = *prev_cell;
    }else{
        left_cell = *this;
    }

    //(3) Prepare the new cells
    my_float new_dr  = this->dr / 2.;
    my_float left_r  = this->r - this->dr /4. ;
    my_float right_r = this->r + this->dr /4. ;

    //(3.1) Right cell
    my_float R_data[NUM_EQUATIONS];
    my_float R_P, R_T, R_gamma, R_Z, R_tcool, R_meanweight;

    for(int elem=0;elem<NUM_EQUATIONS;++elem){
        R_data[elem] = select_interpolation(right_r, this->location(), this->y[elem], right_cell.location(), right_cell.y[elem]);
    }

    R_T = select_interpolation(right_r, this->location(), this->get_temperature(), right_cell.location(), right_cell.get_temperature());
    R_gamma = select_interpolation(right_r, this->location(), this->get_adiabatic_const(), right_cell.location(), right_cell.get_adiabatic_const());
    R_Z = select_interpolation(right_r, this->location(), this->get_metallicity(), right_cell.location(), right_cell.get_metallicity());

    //Pressure is trickier
    R_P = pressure(R_data[0],R_data[1],R_data[2]/R_data[0],R_gamma);
    if(this->system_E == false || R_P == FAILURE){
        //Dual energy again!
        R_P = pressure(R_data[0],R_data[3],R_gamma);
        *right = Cell(R_data, R_P, R_T, R_gamma, R_Z, right_r-0.5*new_dr, right_r+0.5*new_dr, false);
    }else{
        *right = Cell(R_data, R_P, R_T, R_gamma, R_Z, right_r-0.5*new_dr, right_r+0.5*new_dr, true);
    }
    //Cooling time is even more trickier, because *right will now have t_cool = -INFINITY
    R_tcool = - select_interpolation(right_r, this->location(), -this->get_cooling_time(), right_cell.location(), -right_cell.get_cooling_time());
    R_meanweight = select_interpolation(right_r, this->location(), this->get_meanweight(), right_cell.location(), right_cell.get_meanweight());

    right->update_tcool(R_tcool);
    right->update_meanweight(R_meanweight);

    //(3.2) Left cell
    //You can update directly!
    my_float old_r  = this->r;
    my_float old_dr = this->dr;
    this->r  = left_r;
    this->dr = new_dr;

    for(int elem=0;elem<NUM_EQUATIONS;++elem){
        this->y[elem] = select_interpolation(left_r, left_cell.location(), left_cell.y[elem], old_r, this->y[elem]);
    }

    this->T = select_interpolation(left_r, left_cell.location(), left_cell.get_temperature(), old_r, this->get_temperature());
    this->gamma = select_interpolation(left_r, left_cell.location(), left_cell.get_adiabatic_const(), old_r, this->get_adiabatic_const());
    this->Z = select_interpolation(left_r, left_cell.location(), left_cell.get_metallicity(), old_r, this->get_metallicity());

    this->P = pressure(this->y[0],this->y[1],this->y[2]/this->y[0],this->gamma);
    if(this->system_E == false || this->P == FAILURE){
        this->P = pressure(this->y[0],this->y[3],this->gamma);
        system_E = false;
    }
    this->t_cool = -select_interpolation(left_r, left_cell.location(), -left_cell.get_cooling_time(), old_r, -this->get_cooling_time());
    this->meanweight = select_interpolation(left_r, left_cell.location(), left_cell.get_meanweight(), old_r, this->get_meanweight());

    //(4) Linking new cells
    this->prev = prev_cell;
    if(this->prev != nullptr){
        prev_cell->next = this;
    }

    this->next = right;
    right->prev = this;

    right->next = next_cell;
    if(right->next != nullptr){
        next_cell->prev = right;
    }

    //Touch Level of Refinement (LoR)
    this->LoR += 1;
    right->LoR = this->LoR;

    //And reevaluate cells
    this->flag_for_merge = false; //What's the point of merging a cell that has been split? In fact, that leads to errors
    this->flag_for_split = false;
    this->flag_for_shock = old_flag_for_shock;
    
    right->flag_for_merge = false;
    right->flag_for_split = false;
    right->flag_for_shock = old_flag_for_shock;

    //That's all!
}

void Cell::split(Cell* grid, const int max_N, Spline* gridish){ //Size of gridish is N_interpols, which is part of the enum (a global variable)
    /*
        This routine splits the caller (this) into two Cells of equal widths.
        Caller (*this) becomes the LEFT cell
        this routine also find an unused cell in grid, becoming the RIGHT cell
    */

    //(1) Find an unused cell
    Cell* right = find_unused_cell(grid,max_N);
    //(2) Define auxilliary pointers
    Cell* next_cell;
    Cell* prev_cell;
    next_cell = this->next; //This link will be cut shortly, so better save a copy
    prev_cell = this->prev;
    Cell* ambient;
    ambient = find_last(grid);

    int old_flag_for_shock = this->flag_for_shock;

    Cell right_cell;
    if(next_cell != nullptr){
        right_cell = *next_cell;
    }else{
        right_cell = *this;
    }

    Cell left_cell;
    if(prev_cell != nullptr){
        left_cell = *prev_cell;
    }else{
        left_cell = *this;
    }

    #if INTERPOL_RULE == 0
        bool tricky_split = false;
    #elif INTERPOL_RULE == 1
        bool tricky_split = true;
    #else
        bool tricky_split = in_shock_or_ambient(&left_cell,this) | in_shock_or_ambient(this,&right_cell);
    #endif

    //(4) Prepare the new cells
    my_float new_dr  = this->dr / 2.;
    my_float left_r  = this->r - this->dr /4. ;
    my_float right_r = this->r + this->dr /4. ;

    //(4.1) Right cell
    my_float R_data[NUM_EQUATIONS];
    my_float R_T,R_tcool, R_P;

    if(tricky_split = false){
        //Data from splines.
        my_float R_data[NUM_EQUATIONS];
        R_data[0] = gridish[Density].get_value(right_r);
        R_data[1] = gridish[Energy].get_value(right_r);
        R_data[2] = gridish[Momentum].get_value(right_r);
        R_data[3] = gridish[Entropy].get_value(right_r);

        R_T = gridish[Temperature].get_value(right_r);
        R_tcool = gridish[Cooling].get_value(right_r);
        R_P = gridish[Pressure].get_value(right_r);
        }else{
        //Data from lineal imterpolation as well
        for(int elem=0;elem<NUM_EQUATIONS;++elem){
            R_data[elem] = select_interpolation(right_r, this->location(), this->y[elem], right_cell.location(), right_cell.y[elem]);
        }
        
        R_T = select_interpolation(right_r, this->location(), this->get_temperature(), right_cell.location(), right_cell.get_temperature());
        R_tcool = - select_interpolation(right_r, this->location(), -this->get_cooling_time(), right_cell.location(), -right_cell.get_cooling_time());
        R_P = select_interpolation(right_r, this->location(), this->get_pressure(), right_cell.location(), right_cell.get_pressure());
        }
        bool R_system = this->system_E;

        //Data from lineal interpolation
        my_float R_gamma = select_interpolation(right_r, this->location(), this->get_adiabatic_const(), right_cell.location(), right_cell.get_adiabatic_const());
        my_float R_Z = select_interpolation(right_r, this->location(), this->get_metallicity(), right_cell.location(), right_cell.get_metallicity());
        my_float R_meanweight = select_interpolation(right_r, this->location(), this->get_meanweight(), right_cell.location(), right_cell.get_meanweight());

        //Data from physics
        *right = Cell(R_data, R_P, R_T, R_gamma, R_Z, right_r-0.5*new_dr, right_r+0.5*new_dr, R_system);
        //Update cell tcool and meanweight
        right->update_tcool(R_tcool);
        right->update_meanweight(R_meanweight);

    //(4.2) Left cell
    //You can update directly!
    my_float old_r  = this->r;
    my_float old_dr = this->dr;
    this->r  = left_r;
    this->dr = new_dr;

    if(tricky_split = false){
        //Data from splines
        this->y[0] = gridish[Density].get_value(left_r);
        this->y[1] = gridish[Energy].get_value(left_r);
        this->y[2] = gridish[Momentum].get_value(left_r);
        this->y[3] = gridish[Entropy].get_value(left_r);

        this->T = gridish[Temperature].get_value(left_r);
        this->t_cool = gridish[Cooling].get_value(left_r);
        this->P = gridish[Pressure].get_value(left_r);
        
    }else{
        //Lineal interpolation
        for(int elem=0;elem<NUM_EQUATIONS;++elem){
            this->y[elem] = select_interpolation(left_r, left_cell.location(), left_cell.y[elem], old_r, this->y[elem]);
        }

            this->T = select_interpolation(left_r, left_cell.location(), left_cell.get_temperature(), old_r, this->get_temperature());
            this->t_cool = -select_interpolation(left_r, left_cell.location(), -left_cell.get_cooling_time(), old_r, -this->get_cooling_time());
            this->P = select_interpolation(left_r, left_cell.location(), left_cell.get_pressure(), old_r, this->get_pressure());
    }

    //Data from lineal interpolation
    this->gamma = select_interpolation(left_r, left_cell.location(), left_cell.get_adiabatic_const(), old_r, this->get_adiabatic_const());
    this->Z = select_interpolation(left_r, left_cell.location(), left_cell.get_metallicity(), old_r, this->get_metallicity());
    this->meanweight = select_interpolation(left_r, left_cell.location(), left_cell.get_meanweight(), old_r, this->get_meanweight());
        
    //(5) Linking new cells

    this->prev = prev_cell;
    if(this->prev != nullptr){
        prev_cell->next = this;
        //prev_cell->need_split(ambient);
        //prev_cell->need_merge(ambient);
    }

    this->next = right;
    right->prev = this;

    right->next = next_cell;
    if(right->next != nullptr){
        next_cell->prev = right;
    }

    //Touch LoR
    this->LoR += 1;
    right->LoR = this->LoR;

    //And reevaluate cells
    this->flag_for_merge = false; //What's the point of merging a cell that has been split? In fact, that leads to errors
    this->flag_for_split = false;
    this->flag_for_shock = old_flag_for_shock;
    
    right->flag_for_merge = false;
    right->flag_for_split = false;
    right->flag_for_shock = old_flag_for_shock;

    //That's all!

}


void Cell::merge(Cell* merged, Cell* old_grid, Cell* grid, const int curr_N){
    /*
        This function merges *this with *merged.
        The result will become *this
        And *merged will become unused

        Ideally, you want *merged to be this->next;
    */

    //(1) Prepare result
    my_float new_dr = this->width() + merged->width();
    my_float new_r  = 0.5*( this->location() + merged->location() );
    Cell* ambient;
    ambient = find_last(grid); //Warning here if the shock reaches the right simulation boundary!!!

    Cell right_cell = right_neighbour(new_r,old_grid,curr_N);
    Cell left_cell  = left_neighbour(new_r,old_grid,curr_N);

    //(2) Update *this directly
    this->r  = new_r;
    this->dr = new_dr;

    for(int elem=0;elem<NUM_EQUATIONS;++elem){
        this->y[elem] = select_interpolation(new_r, left_cell.location(), left_cell.y[elem], right_cell.location(), right_cell.y[elem]);
    }

    this->T = select_interpolation(new_r, left_cell.location(), left_cell.get_temperature(), right_cell.location(), right_cell.get_temperature());
    this->gamma = select_interpolation(new_r, left_cell.location(), left_cell.get_adiabatic_const(), right_cell.location(), right_cell.get_adiabatic_const());
    this->Z = select_interpolation(new_r, left_cell.location(), left_cell.get_metallicity(), right_cell.location(), right_cell.get_metallicity());

    this->P = pressure(this->y[0],this->y[1],this->y[2]/this->y[0],this->gamma);
    if(this->system_E == false || this->P == FAILURE){
        this->P = pressure(this->y[0],this->y[3],this->gamma);
        system_E = false;
    }
    this->t_cool = select_interpolation(new_r, left_cell.location(), left_cell.get_cooling_time(), right_cell.location(), right_cell.get_cooling_time());
    this->meanweight = select_interpolation(new_r, left_cell.location(), left_cell.get_meanweight(), right_cell.location(), right_cell.get_meanweight());

    //(3) Linking
    Cell* next_cell;
    next_cell = merged->next;
    Cell* prev_cell;
    prev_cell = this->prev;
    make_unused_cell(merged);

    this->next = next_cell;
    if(this->next != nullptr){
        next_cell->prev = this;
    }
    this->prev = prev_cell;
    if(this->prev != nullptr){
        prev_cell->next = this;
    }

    //(4) Touch LoR
    this->LoR -= 1;

    this->flag_for_split = false; //What's the point of splitting a cell that has been merged? this leads to errors
    this->flag_for_merge = false; //It can be reconsidered again, but I prefer to play safe

    //That's all!
}

void Cell::merge(Cell* merged, Spline* gridish){
    /*
        This function merges *this with *merged.
        The result will become *this
        And *merged will become unused

        Ideally, you want *merged to be this->next;
    */

    //(1) Prepare result
    my_float new_dr = this->width() + merged->width();
    my_float new_r  = 0.5*( this->location() + merged->location() );

    #if INTERPOL_RULE == 0
        bool tricky_cell = false;
    #elif INTERPOL_RULE == 1
        bool tricky_cell = true;
    #else
        bool tricky_cell = in_shock_or_ambient(this,merged);
    #endif

    Cell right_cell = *merged;
    Cell left_cell = *this;

    //(2) Update *this directly
    my_float old_r  = this->r;
    my_float old_dr = this->dr;
    this->r  = new_r;
    this->dr = new_dr;

    if(tricky_cell = false){
        //Data from splines
        this->y[0] = gridish[Density].get_value(new_r);
        this->y[1] = gridish[Energy].get_value(new_r);
        this->y[2] = gridish[Momentum].get_value(new_r);
        this->y[3] = gridish[Entropy].get_value(new_r);

        this->T = gridish[Temperature].get_value(new_r);
        this->t_cool = gridish[Cooling].get_value(new_r);
        this->P = gridish[Pressure].get_value(new_r);
    }else{
        //Lineal interpolation as usual
        for(int elem=0;elem<NUM_EQUATIONS;++elem){
            this->y[elem] = select_interpolation(new_r, left_cell.location(), left_cell.y[elem], right_cell.location(), right_cell.y[elem]);
        }

        this->T = select_interpolation(new_r, left_cell.location(), left_cell.get_temperature(), right_cell.location(), right_cell.get_temperature());
        this->t_cool = select_interpolation(new_r, left_cell.location(), left_cell.get_cooling_time(), right_cell.location(), right_cell.get_cooling_time());
        this->P = select_interpolation(new_r, left_cell.location(), left_cell.get_pressure(), right_cell.location(), right_cell.get_pressure());
    }

    //Lineal interpolations for the rest
    this->gamma = select_interpolation(new_r, left_cell.location(), left_cell.get_adiabatic_const(), right_cell.location(), right_cell.get_adiabatic_const());
    this->Z = select_interpolation(new_r, left_cell.location(), left_cell.get_metallicity(), right_cell.location(), right_cell.get_metallicity());
    this->meanweight = select_interpolation(new_r, left_cell.location(), left_cell.get_meanweight(), right_cell.location(), right_cell.get_meanweight());

    //(3) Linking
    Cell* next_cell;
    next_cell = merged->next;
    Cell* prev_cell;
    prev_cell = this->prev;
    make_unused_cell(merged);

    this->next = next_cell;
    if(this->next != nullptr){
        next_cell->prev = this;
    }
    
    this->prev = prev_cell;
    if(this->prev != nullptr){
        prev_cell->next = this;
    }

    //(4) Touch LoR
    this->LoR -= 1;

    this->flag_for_split = false; //What's the point of splitting a cell that has been merged? this leads to errors
    this->flag_for_merge = false; //It can be reconsidered again, but I prefer to play safe

    //That's all!
}


///Main AMR algorithm

void refine_grid(Cell* grid,const int max_N, const bool delay_time = false){
    /*
        This routine does:
        -Creates a copy of current grid (using recent merged/split cells for interpolations leads to oscillations)
        -Checks cells to be split
        -Checks cells to be merged
    */

    //(1) Copy old grid and auxiliaries
    Cell* curr;
    Cell* first;
    Cell* ambient;
    first = find_first(grid);
    curr = first;
    ambient = find_last(grid);

    ///NEXT LOOP TAKES THE FOLLOWING STEPS
    /**(1) Create a copy of previous grid in order to get an interpolating spline**/
    /**(2) Flag if a cell needs to be merged or split **/
    const int curr_N = count_cells(grid);

    my_float** table = new my_float* [N_interpols];
    for(int elem=0;elem<N_interpols;++elem){
        table[elem] = new my_float[curr_N];
    }

    for(int i=0;i<curr_N;++i){
        //Copy the grid into table
        table[Position][i]    = curr->location();
        table[Density][i]     = curr->eq(0);
        table[Energy][i]      = curr->eq(1);
        table[Momentum][i]    = curr->eq(2);
        table[Entropy][i]     = curr->eq(3);
        table[Temperature][i] = curr->get_temperature();
        table[Pressure][i]    = curr->get_pressure();
        table[Cooling][i]     = curr->get_cooling_time();

        // An important side effect I exploit in both split and merge routines above is that old_grid has their cells ordered (i= n is guaranteed to be 'n' position)
        // Another side effect is that old_grid[i].next is the location of next cell in grid, not old_grid. So better not to use!

        //Initialize cell flags
        curr->prepare_for_refinement();

        //Check if curr needs to be split/merged
        curr->check_refinement(ambient,delay_time);

        curr = curr->next;
    }

    //Initialize the Splines
    Spline* interpolator = new Spline [N_interpols];

    interpolator[Density] = Spline(table[Position],table[Density],curr_N);
    interpolator[Energy] = Spline(table[Position],table[Energy],curr_N);
    interpolator[Momentum] = Spline(table[Position],table[Momentum],curr_N);
    interpolator[Entropy] = Spline(table[Position],table[Entropy],curr_N);
    interpolator[Temperature] = Spline(table[Position],table[Temperature],curr_N);
    interpolator[Pressure] = Spline(table[Position],table[Pressure],curr_N);
    interpolator[Cooling] = Spline(table[Position],table[Cooling],curr_N);

    //Delete unneccesary memory (because it's now a duplicate)
    for(int elem=0;elem<N_interpols;++elem){
        delete[] table[elem];
    }
    delete[] table;

    ///SPLIT/MERGE CELLS BASED ON THEIR FLAGS
    curr = first;
    int i_test = 0;
    while(curr->next != nullptr){

        i_test++;
        //Does a cell need to be split?
        if( curr->should_be_split() == true ){ //In case of should_be_merged == true; split cell wins
            curr->split(grid,max_N,interpolator);
        }
        //Does a cell need to be merged? (This condition is far trickier!)
        else if( (curr->should_be_merged() & curr->next->should_be_merged() ) == true && curr->compare_LoRs(curr->next) == 0){
            curr->merge(curr->next,interpolator);
        }
        //else do nothing

        //Next lines are needed in case that the penultimate cell is merged with the last one. Penultimate cell becomes last, and last->next = nullptr
        if(curr->next != nullptr){ curr = curr->next; }
        else{ break; }
    }

    ///CHECK THE RESULT TO PREVENT FROM OSCILLATIONS
    //If an oscillation happens, split again the cell, due to we have 1 merge per cell, interpolation should revert to original values.
    curr = first->next;
    while(curr->next != nullptr){

        //LoR oscillations
        if(curr->prev->level_of_refinement() == curr->next->level_of_refinement() && curr->compare_LoRs(curr->next) > 0){
            curr->split(grid,max_N,interpolator);
        }
        else if(curr->prev->level_of_refinement() == curr->next->level_of_refinement() && curr->compare_LoRs(curr->next) < 0){
            curr->prev->split(grid,max_N,interpolator);
            curr->next->split(grid,max_N,interpolator);
        }
        //Jumps of 2 in LoRs
        else if( curr->compare_LoRs(curr->prev) < -1 ){
            curr->prev->split(grid,max_N,interpolator);
        }
        else if( curr->compare_LoRs(curr->next) < -1 ){
            curr->next->split(grid,max_N,interpolator);
        }

        curr = curr->next;
    }
    //Finished!
    delete[] interpolator;

}

#endif //Ends AMR_ON
#endif //Ends ifndef
