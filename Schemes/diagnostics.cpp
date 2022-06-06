#ifndef DIAGNOSTICS_CPP
#define DIAGNOSTICS_CPP

/**
    These functions are used for debugging
**/

bool warning_about_very_small_values_given = false;

//Grid diagnostics
void grid_diagnostics(Cell* grid, const int N = NCELLS * pow(2,MAX_REF - INI_REF)){
    std::cout<<"Counting grid usage..."<<std::endl;

    int used_cells = 0; //next and/or prev cell are linked
    int boundary_cells = 0; //next or prev cell are linked (only one of them)
    int null_cells = 0; //Unused cells, written as null cells to use a different word
    int bad_cells = 0; //Used cells with bad positions
    int i; //Grid position

    Cell* curr;
    curr = grid;

    for(i=0;i<N;++i){
        if(curr->next != nullptr && curr->prev != nullptr ){
            used_cells++;
            if( curr->location() >= curr->next->location() || curr->location() <= curr->prev->location() ){
                bad_cells++;
            }
        }
        else if(curr->next != nullptr || curr->prev != nullptr){ used_cells++; boundary_cells++; }
        else if(curr->next == nullptr && curr->prev == nullptr){ null_cells++; }
        curr++;
    }
    std::cout<<"Grid count   = "<<i<<std::endl;
    std::cout<<"Used cells   = "<<used_cells<<std::endl;
    std::cout<<"Boundaries   = "<<boundary_cells<<std::endl;
    std::cout<<"Unused cells = "<<null_cells<<std::endl;
    std::cout<<"Bad cells    = "<<bad_cells<<std::endl;
}

//Timestep diagnostics
bool exceeded_time(time_t starting_time){
    //Simulation took too long
    const my_float max_time = 7200.; //In seconds

    my_float current_time = difftime(time(NULL),starting_time);
    if(current_time > max_time){
        std::cout<<"Warning: Simulation took too much time. Stopping..."<<std::endl;
        return true;
    }
    else{
        return false;
    }

}

bool ambient_is_kept(Cell* grid, const int N, Cell ambient){
    //Returns true if last cell matches with the reference (usually, the ambient)
    Cell* last;
    last = find_last(grid);
    bool veredict = true;

    const my_float relative_error = 1.E-8;
    for(int elem=0;elem<NUM_EQUATIONS;++elem){
        if( std::abs( (last->prev->eq(elem) - ambient.eq(elem)) /ambient.eq(elem)) > relative_error && ambient.eq(elem) != 0.){ //Second condition is to cover against zeros.
            veredict = false;
            std::cout<<"Warning: Simulation reached a boundary. Stopping..."<<std::endl;
            break;
        }else if(warning_about_very_small_values_given == false && std::abs(last->prev->eq(elem)) > std::abs(ambient.eq(elem))){
            std::cout<<"Warning: Very small variation in ambient detected in ";
            switch(elem){
                case 0:
                    std::cout<<"Density ";
                    break;
                case 1:
                    std::cout<<"Energy ";
                    break;
                case 2:
                    std::cout<<"Momentum ";
                    break;
                case 3:
                    std::cout<<"Entropy ";
                    break;
                default:
                    std::cout<<"Unexpected ";
            }
            std::cout<<"variable. Ignoring..."<<std::endl;
            warning_about_very_small_values_given = true;
        }
    }

    return veredict;
}

#endif
