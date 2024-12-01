/*****************************************************************************
*Class to create lists with random numbers

*16/05/2022 Lauritz Klünder
******************************************************************************/

#include <vector>
#include <iostream>

#include "RandNumberBetween.hpp"
#include "RandomList.hpp"


void RandomList::generateList(  ){
    if ( gauss == false ){
        for (unsigned int i = 0; i < listLen; i++) {
		    randomVector[i] = RandNum_Normal();
	    }
    } 
    else{
        for (unsigned int i = 0; i < listLen; i++) {
		    randomVector[i] = RandNum_Gauss();
	    }
    } 
}

double RandomList::getNumber(  ){

    double randNumber;

    if ( counter == listLen ){
        counter = 0;
    } 

    if ( counter == 0 ){
        generateList();
        randNumber = randomVector[counter]; 
    } 
    else{
        randNumber = randomVector[counter]; 
    } 

    counter++;

    return randNumber;

}