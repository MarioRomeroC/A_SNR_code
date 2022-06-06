#ifndef SPLINE_CLASS_H
#define SPLINE_CLASS_H

//Enhanced natural cubic spline class, based on previous code.

class Spline{
    //Note: needs cmath and algorithm libraries to work

    private:
        //Your data
        my_float* x;
        my_float* y;

        int N;

        //Needed data
        my_float* d; //Second derivatives

        ///Private functions
        void find_interval(my_float X, int& i1, int& i2){
            //Find the subinterval [x[i1], x[i2]] in which X is located

            if(X < x[0]){//First case is X < x[0]
                i1 = 0;
                i2 = 1;
            }
            else if(X > x[N-1]){ //Second case is the opposite: X>x[N-1]
                i1 = N-2;
                i2 = N-1;
            }
            else{//Third, and most common case
                for(int i=1;i<N;++i){
                    if( x[i-1]<=X && X<=x[i] ){ //Got them!
                        i1 = i-1;
                        i2 = i;
                        break;
                    }
                }
            }
        }

        void find_nearest(my_float X, int& i1, int& i2){
            //Given x find nearest values (by their indexes i1, i2) in x*
            //i1 and i2 are your outputs, as long they are declared variables in the scope in which this function is called this works fine
            my_float dist =  1e6 * std::abs(x[N-1]-x[0]); //Arbitrary large distance
            my_float low1 = dist;
            my_float low2 = dist;
            i1 = 0;
            i2 = 1;
            //First nearest
            for(int i=0;i<N;i++){
                dist = std::abs(X-x[i]);
                if( dist < low1 ){
                    low1 = dist;
                    i1   = i;
                }
            }

            //Second nearest;
            for(int i=0;i<N;i++){
                if(i == i1){continue;} //Exclude first nearest for this loop
                dist = std::abs(X-x[i]);
                if( dist < low2 ){
                    low2 = dist;
                    i2   = i;
                }
            }

            //I need these indexes ordered as i1 = i, and i2 = i+1
            if(i1>i2){
                my_float temp = i2;
                i2 = i1;
                i1 = temp;
            }
        }

        void sort_values(my_float* X, my_float* Y){ //Both arrays' size is N, and their result is *y created
            //Sort *x is easy:
            // std::sort(x, x+N);
            // However we have to rearrange the positions of Y according to the sorted x array
            // so... I don't know if it's better to call sort x, and after rearrange y
            // or sorting with my code (slower), but at the same time assign the y.
            // I'm taking first option

            y = new my_float[N];

            bool *assigned = new bool[N];
            for(int j=0;j<N;++j){
                assigned[j] = false;
            }
            /* Imagine that I give you these arrays:
                x = [ 1, 3, 5, 6, 6, 6, 6, 6]
                y = [ 2, 4, 0, 1, 1, 3, 5, 7]
               You are interpolating something that is not a function.
               Having x[i] == X[j] does not guarantee saving latter points (If you indeed have a function, it doesn't matter, all y(x) is uniquely defined).
            */

            for(int i=0;i<N;++i){
                for(int j=0;j<N;++j){
                    if(x[i] == X[j] && assigned[j] == false){ //Shoudn't give an error of 1 != 1 if is called afterwards initialization and sort, no mathematical operation shoudn't have been done!
                        y[i] = Y[j];
                        assigned[j] = true;
                        break; //Right position found! No need to continue
                    }
                }
            }

            delete[] assigned;
        }

        public:
            Spline(){
                x = nullptr;
                y = nullptr;
                d = nullptr;
            }
            //C++ rule of three: If you have pointer members to dynamic memory, you need a destructor, a copy constructor, and assigment operator.
            ~Spline(){ //Destructor
                if( x != nullptr ){
                    delete[] x;
                    x = nullptr;
                }
                if( y != nullptr ){
                    delete[] y;
                    y = nullptr;
                }
                if( d != nullptr ){
                    delete[] d;
                    d = nullptr;
                }

            }
            Spline(const Spline& other){ //Copy Constructor

                //Create new dynamic memory for the copy
                x = new my_float[other.N];
                y = new my_float[other.N];
                d = new my_float[other.N];

                N = other.N;

                for(int i=0;i<N;++i){
                    x[i] = other.x[i];
                    y[i] = other.y[i];
                    d[i] = other.d[i];
                }

            }
            Spline& operator=(const Spline& other ){ //Assignment operator
                //This represent the operation *this = other. WE ARE OVERRIDING *this, with potential memory leaks if we are not careful!

                //Prepare new arrays
                my_float* new_x = new my_float[other.N];
                my_float* new_y = new my_float[other.N];
                my_float* new_d = new my_float[other.N];

                for(int i=0;i<other.N;++i){
                    new_x[i] = other.x[i];
                    new_y[i] = other.y[i];
                    new_d[i] = other.d[i];
                }

                //Now override *this
                if( x != nullptr ){ delete[] x; }
                if( y != nullptr ){ delete[] y; }
                if( d != nullptr ){ delete[] d; }
                //Unlike destructors, {x,y,f} won't be null pointers, but...
                x = new_x;
                y = new_y;
                d = new_d;

                N = other.N;

            }

            ///CONSTRUCTOR
            Spline(my_float* X, my_float* Y, const int n){
                //Arguments are two dynamic arrays, it does not need to be sorted

                N = n;
                x = new my_float[N];
                //y is 'new'ed in sort_values
                d = new my_float[N];


                //Init x
                for(int i=0;i<N;++i){
                    x[i] = X[i];
                }

                //Init y
                std::sort(x,x+N);
                sort_values(X,Y); //y is 'new'ed there.

                //Now the complicated part: obtain the second derivatives!
                //which leads to solving a tridiagonal matrix

                d[0]   = 0.;
                d[N-1] = 0.; //Natural spline condition: second derivatives are zero at the borders

                //Because we have already assign two derivatives, these arrays are N-2 in size
                my_float* alpha = new my_float[N-2];
                my_float* beta  = new my_float[N-2];
                my_float* gamma = new my_float[N-2];
                my_float* delta = new my_float[N-2];
                my_float temp1,temp2;

                int j;
                for(int i=1;i<N-1;i++){ //Fill values
                        j = i-1; //To fit correctly these new arrays
                        alpha[j] = (x[i] - x[i-1])/6.;
                        beta[j]  = (x[i+1] - x[i-1])/3.;

                        temp1 = (x[i+1] - x[i])/6.;
                        temp2 = (y[i+1] - y[i])/( x[i+1] - x[i] ) - (y[i] - y[i-1])/( x[i] - x[i-1] );

                        if(i==1){
                            gamma[j] = temp1 / beta[j];
                            delta[j] = temp2 / beta[j];
                        }else{
                            gamma[j] = temp1 / ( beta[j] - alpha[j]*gamma[j-1] );
                            delta[j] = ( temp2 - alpha[j]*delta[j-1] ) / ( beta[j] - alpha[j]*gamma[j-1] );
                        }
                }

                //Now solve the system
                int i;
                for(j=N-3; j>0;j--){
                    i=j+1;
                    if(j==N-3){
                        d[i] = delta[j];
                    }else{
                        d[i] = delta[j] - gamma[j]*d[i+1];
                    }
                }

                //Freeing auxilliary arrays
                delete[] alpha;
                delete[] beta;
                delete[] gamma;
                delete[] delta;

            }

            my_float get_value(my_float X){
                //Find the suitable points
                int i1;
                int i2;

                find_interval(X,i1,i2);

                my_float x1 = x[i1];
                my_float y1 = y[i1];
                my_float d1 = d[i1];

                my_float x2 = x[i2];
                my_float y2 = y[i2];
                my_float d2 = d[i2];

                //And now construct the cubic spline in that interval, and get the value
                my_float A = ( x2 - X )/( x2 - x1 );
                my_float B = 1.-A;
                my_float C = A*(A*A - 1.)*(x2-x1)*(x2-x1)/6.;
                my_float D = B*(B*B - 1.)*(x2-x1)*(x2-x1)/6.;

                return A*y1 + B*y2 + C*d1 + D*d2;
            }


};


#endif
