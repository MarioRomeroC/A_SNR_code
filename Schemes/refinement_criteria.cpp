#ifndef REFINEMENT_CRITERIA
#define REFINEMENT_CRITERIA

#define MAX_GRADIENT 0.1
#define MIN_GRADIENT 0.005

bool check_gradient(my_float left_value, my_float right_value, my_float curr_value, my_float threesold){

    my_float gradient = std::abs( (right_value-left_value)/curr_value );

    if(gradient > threesold){ return true; }
    else{ return false; }
}


void Cell::needs_merge(Cell* ambient, bool extern_condition = false){
    //This function is not called if cell is flagged for splitting

    //Do NOT merge if a potential jump of two is created
    if(this->compare_LoRs(this->prev) > 0 || this->compare_LoRs(this->next) > 0){
        this->flag_for_merge = false;
        return;
    }
    else if(this->prev != nullptr && this->prev->flag_for_split == true){ //Note that we advance with curr = curr->next. The condition for next cannot be evaluated
        this->flag_for_merge = false;
        return;
    }

    //Do NOT merge if an external condition is true (and cell is not ambient
    if(extern_condition){
        if(this->y[2] != ambient->y[2]){
            this->flag_for_merge = false;
        }else{
            this->flag_for_merge = true;
        }
        return;
    }

    //To avoid segfaults...
    Cell* right;
    if(this->next != nullptr){ right = this->next; }
    else{ right = this; }
    Cell* left;
    if(this->prev != nullptr){ left = this->prev; }
    else{ left = this; }

    //MERGE if energy density gradient is low
    if(this->flag_for_gradient > 0){
        this->flag_for_merge = false;
        return;
    }
    if(check_gradient(left->y[1],right->y[1],this->y[1],MIN_GRADIENT) == false){
        this->flag_for_merge = true;
        return;
    }

    //That's all!

}

void Cell::needs_split(Cell* ambient, bool extern_condition = false){

    //Split if we are in a situation of a jump of 2
    if(this->compare_LoRs(this->prev) > 1 || this->compare_LoRs(this->next) > 1){
        this->flag_for_split = true;
        return;
    }
    else if(this->compare_LoRs(this->prev) == 1 && this->prev->flag_for_split == true){
        this->flag_for_split = true;
        return;
    }

    //Split if an external condition is true AND the cell is not an ambient one
    if(extern_condition == true){
        if(this->y[2] != ambient->y[2]){
            this->flag_for_split = true;
        }else{
            this->flag_for_split = false;
        }

        return;
    }


    //To avoid segfaults...
    Cell* right;
    if(this->next != nullptr){ right = this->next; }
    else{ right = this; }
    Cell* left;
    if(this->prev != nullptr){ left = this->prev; }
    else{ left = this; }

    if(this->flag_for_gradient > 0){
        this->flag_for_split = true;
        return;
    }
    if(check_gradient(left->y[1],right->y[1],this->y[1],MAX_GRADIENT) == true){
        this->flag_for_split = true;
        return;
    }

    //That's all!
}

void Cell::near_shock(){

    const int num_neighbors = MAX_REF;
    //Security first!
    Cell* right;
    if(this->next != nullptr){ right = this->next; }
    else{ right = this; }
    Cell* left;
    if(this->prev != nullptr){ left = this->prev; }
    else{ left = this; }

    if(in_shock(left,right,this) == true){
        //Flag this cell and neighbors
        this->flag_for_shock = 1;
        return;
    }
    if(left->prev != nullptr){
        //Look what you have behind
        if(left->flag_for_shock > 0 && left->flag_for_shock != num_neighbors){
            this->flag_for_shock = left->flag_for_shock + 1;
        }
    }
    if(right->next != nullptr){
        int i=1;
        //Look what you have next
        while(i<num_neighbors && right->next != nullptr){ //We advance in right, so repetition is necessary
            i++;
            if(in_shock(right->prev,right->next,right) == true){
                this->flag_for_shock = i;
                break;
            }else{
                right = right->next;
            }
        }
    }

    return;

}

void Cell::high_gradient(){

    const int num_neighbors = MAX_REF;
    const int min_ref = MAX_REF - 1; //We are checking HIGH gradients, so (1) You are already in MAX_REF and gradient is still high or
                                                                        //(2) you are in MAX_REF-1 and gonna get split to MAX_REF anyway
    //Security first!
    Cell* right;
    if(this->next != nullptr){ right = this->next; }
    else{ right = this; }
    Cell* left;
    if(this->prev != nullptr){ left = this->prev; }
    else{ left = this; }

    if( this->LoR >= min_ref && check_gradient(left->y[1],right->y[1],this->y[1],MAX_GRADIENT) == true){
        this->flag_for_gradient = 1;
        return;
    }
    if(left->prev != nullptr){
        //Look what you have behind
        if(left->flag_for_gradient > 0 && left->flag_for_gradient != num_neighbors){
            this->flag_for_gradient = left->flag_for_gradient + 1;
        }
    }
    if(right->next != nullptr){
        int i = 1;
        //Look what you have next
        while(i<num_neighbors && right->next != nullptr){
            i++;
            if( right->LoR >= min_ref && check_gradient(right->prev->y[1],right->next->y[1],right->y[1],MAX_GRADIENT) == true){
                this->flag_for_gradient = i;
                break;
            }else{
                right = right->next;
            }
        }
    }
    return;
}

void Cell::check_refinement(Cell* ambient, bool extern_condition = false){
    //OK, this is the function to check if a cell has to be split/merged or not
    //...and is a boolean nightmare of conditions. You've been warned

    this->near_shock();
    this->high_gradient();
    if(this->LoR < MAX_REF){ this->needs_split(ambient, extern_condition); }
    if(this->LoR > 0 && this->flag_for_split != true){ this->needs_merge(ambient,extern_condition); }

    return;
}

#endif
