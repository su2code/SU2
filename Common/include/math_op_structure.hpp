
#pragma once
#include <complex>
#include <iostream>
#include <fstream>
#include <valarray>
#include "../../Common/include/config_structure.hpp"


typedef std::complex<su2double> Complex;
typedef std::valarray<Complex> CArray;


class FFT {
  public:




	/*!
	 * \brief Constructor of the  class.
	 */
	FFT();

	/*!
	 * \brief Destructor of the class.
	 */
	~FFT(void);


        void fft_r2(CArray& x);
        void ifft_r2(CArray& x);


};

