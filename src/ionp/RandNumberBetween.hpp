/*****************************************************************************
*Class to define a mathematical a random number (header .hpp)
*
* how to: create an object with creator and give the boundaries eg. RandNumberBetween r(0,100)
* draw a random number using () operator e.g r(), output is different each time eg. 0<r()<100

*16/05/2022 Lauritz Klünder
******************************************************************************/
#ifndef RANDNUMBERBETWEEN_HPP_
#define RANDNUMBERBETWEEN_HPP_

#include <algorithm>
#include <iostream>
#include <random>
#include <vector>
#include <chrono> //random generator (from time)

#include "mptmacros.h"


class RandNumberBetween {
public:
	RandNumberBetween(double low, double high)
		: random_engine( SC_UI(std::chrono::system_clock::now().time_since_epoch().count()) ),
		distribution{ low, high } {}

	double operator() () {
		return distribution(random_engine);
	}

private:
	std::mt19937 random_engine;
	std::uniform_real_distribution<double> distribution;
};

class RandNumberGaussian {
public:

	RandNumberGaussian(double mean, double std)
		: random_engine( SC_UI(std::chrono::system_clock::now().time_since_epoch().count()) ),
		distribution{ mean, std } {}

	double operator() () {
		return distribution(random_engine);
	}

private:
	std::mt19937 random_engine;
	std::normal_distribution<double> distribution;
};

class RandNumberLognorm {
public:

	RandNumberLognorm(double mean, double std)
		: random_engine( SC_UI(std::chrono::system_clock::now().time_since_epoch().count()) ),
		distribution{ mean, std } {}

	double operator() () {
		return distribution(random_engine);
	}

private:
	std::mt19937 random_engine;
	std::lognormal_distribution<double> distribution;
};


class RandNumberHist {
public:

	RandNumberHist(const std::vector<double>& bins_var, const std::vector<int>& frequencies_var)
		: bins(bins_var),
		frequencies(frequencies_var) {}

	// Function to generate a random number based on the histogram
	double generateRandomNumber() {
	    // Calculate the total frequency
	    int totalFrequency = 0;
	    for (int freq : frequencies) {
	        totalFrequency += freq;
	    }

	    // Generate a random number between 1 and the total frequency
	    int randomNumber = rand() % totalFrequency + 1;

	    // Find the bin corresponding to the random number
	    int cumulativeFrequency = 0;
	    for (size_t i = 0; i < bins.size(); ++i) {
	        cumulativeFrequency += frequencies[i];
	        if (randomNumber <= cumulativeFrequency) {
	            return bins[i];
	        }
	    }

	    // This line should never be reached
	    return -1.0;
	}

	std::vector<double> randomVector(const int lenVec){
		srand(static_cast<unsigned int>(time(nullptr)));
		std::vector<double> Vec;
		for (int i = 0; i<lenVec; i++){
			double random = generateRandomNumber();
			Vec.push_back(random);
		}

		return Vec;
	}

private:
	std::vector<double> bins;
	std::vector<int> frequencies;
};


#endif
