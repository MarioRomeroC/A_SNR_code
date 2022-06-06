#ifndef AMR_AUX_CPP
#define AMR_AUX_CPP

/**
    Adaptive Mesh Refinement (AMR) auxilliary functions
**/

///Working with adaptive chains
void make_unused_cell(Cell* candidate){
    //Breaks links of candidate, making it effectively unused
    if(candidate->next != nullptr){candidate->next->prev = nullptr;}
    if(candidate->prev != nullptr){candidate->prev->next = nullptr;}

    //Break candidate pointers
    candidate->next = nullptr;
    candidate->prev = nullptr;
}

Cell* find_first(Cell* finder){
    //As title says, return *first using a pointer to an arbitary cell

    while( finder->prev != nullptr ){
        finder = finder->prev;
    }

    return finder;
}

Cell* find_last(Cell* finder){

   while( finder->next != nullptr ){
        finder = finder->next;
    }

    return finder;

}

int count_cells(Cell* finder){
    //Count the number of active cells

    Cell* curr = find_first(finder);

    int i = 0;
    while( curr->next != nullptr ){
        curr = curr->next;
        i++;
    }

    return i;
}

void make_a_link(Cell* left, Cell* right){
    left->next  = right;
    right->prev = left;
}

Cell* find_unused_cell(Cell* grid, const int max_N){
    /*
        When grid is initialized in main.cpp, you're creating the worst case scenario: All cells with max refinement.
        Due to is very unlikely to have that during the simulation, some parts of grid will be used and some others not.
        NOT used parts are cells whose pointers are null for both 'next' and 'prev' and does not take part of the chain.
        This function finds the first one of these, and returns a pointer to it.

        Note that this function works as expected as long as you have initialized the chain first.
    */

    Cell* curr;
    curr = grid;

    int i = 0;

    while(curr->next !=nullptr || curr->prev != nullptr){
        //Breaks when both are null pointers
        curr++; //Move through the array
        i++;
        if(i >= max_N){
            //Something went wrong (you are looking outside the chain of cells
            std::cout<<"Error diagnostics!"<<std::endl;
            grid_diagnostics(grid,max_N);
            delete[] grid;
            throw runtime_error("Moving out of grid!");
        }
    }
    //I have found the first cell whose next and prev are null!
    return curr;
}

///Allowing auxiliary grids

void rearrange_links(Cell* grid, Cell* grid2, const int N = NCELLS * pow(2,MAX_REF - INI_REF)){
    //Given *grid, links *grid2. Does not work if lengths of both arrays are not equal

    Cell* curr;
    curr = grid; //This grid is good

    int j=0; //grid2 count
    while(curr->next != nullptr){
        grid2[j].Link_next(grid2[j+1]);
        grid2[j+1].Link_prev(grid2[j]);
        j++;
        curr = curr->next;
    }
    //That's all!
}

#endif
