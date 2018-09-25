//
// Created by jmo510 on 16/08/18.
//


#include "time_dependence.h"





// Check for consistency with the flux balance equation - outputs the correct value of the albedo
double flux_convergence( double Sext, double K, double E, double a, double b, double original_albedo) {
    bool converged = false;
    double low = 0;
    double mid = Sext/2.0;
    double high = Sext;
    double PREC = 0.01;

    //Initial
   // int counter  = 0;
    while (true) {
   // counter +=1;
   // std::cout<<counter<<std::endl;
        if ((K * pow(E, a) * pow(low, b) + low - Sext) * (K * pow(E, a) * pow(mid, b) + mid - Sext) < 0) {
            //Then the root lies between

            high = mid;
            mid = (mid + low) / 2.0;



        } else if ((K * pow(E, a) * pow(mid, b) + mid - Sext) * (K * pow(E, a) * pow(high, b) + high - Sext) < 0) {
            //THe root lies between

            low = mid;
            mid = (mid + high) / 2.0;
        }


            //std::cout << "| " <<low<<" | "<<mid<<" | "<<high<<"\n";
        //std::cout <<std::scientific <<  " result= " << fabs(K * pow(E, a) * pow(mid, b) + mid - Sext) - PREC*Sext << std::endl;

        if (Sext < 0.0001){
            return original_albedo;
        }

        if (fabs(K * pow(E, a) * pow(mid, b) + mid - Sext) < PREC) {
            //std::cout <<"Done!!\n";
            //std::cout<<1.0 / (1.0 + (1.0 / K ) * pow(E, -a) * pow(mid, 1-b))<<std::endl;
            return 1.0 - mid/Sext; //albedo
        }


    }
}

