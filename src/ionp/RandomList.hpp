/*****************************************************************************
*Class to create lists with random numbers

*16/05/2022 Lauritz Klünder
******************************************************************************/

#ifndef RANDOMLIST_HPP_
#define RANDOMLIST_HPP_

#include "RandNumberBetween.hpp"

class RandomList
{
  public:
                    RandomList          ( double var0, double var1, bool var2 );

    void            generateList        (  );
    double          getNumber           (  );


  private:
    double low, high;
    bool gauss;
    unsigned int counter, listLen;
    std::vector<double> randomVector;
    RandNumberBetween RandNum_Normal;
    RandNumberGaussian RandNum_Gauss;


};

inline RandomList::RandomList( double var0, double var1, bool var2 ):
    low(var0),
    high(var1),
    gauss(var2),
    counter(0),
    listLen(10000),
    randomVector(),
    RandNum_Normal( var0, var1 ),
    RandNum_Gauss( var0, var1 ){

        for ( unsigned int i = 0; i < listLen; i++){
            randomVector.push_back(0.);
        } 
}

#endif