#include <complex>
#include <iostream>
#include <fstream>
#include <valarray>
#include "../include/math_op_structure.hpp"






FFT::FFT() {


}


FFT::~FFT(void) {


}





 //Function, which computes the FFT using Cooley-Tukey method
void FFT::fft_r2(CArray& x)
{


    const size_t N = x.size();
    if (N <= 1) return;

    // Rearrange the input data into two parts (Divide):
    CArray even = x[std::slice(0, N/2, 2)];
    CArray  odd = x[std::slice(1, N/2, 2)];

    // Calculate the FFT recursively (Conquer):
    fft_r2(even);
    fft_r2(odd);

    // Finally combine the data:
    for (size_t k = 0; k < N/2; ++k)
    {
        Complex t = std::polar(su2double(1.0), su2double(-2 * PI_NUMBER * k / N)) * odd[k];
        x[k    ] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}

//Function, which computes the inverse FFT using Cooley-Tukey method
void FFT::ifft_r2(CArray& x)
{
    // First, conjugate the complex numbers
    x = x.apply(std::conj);

    // Second, apply the forward fft to the data
    fft_r2( x );

    // Third, conjugate the complex numbers once more
    x = x.apply(std::conj);

    // Finally, scale the numbers
    x /= su2double(x.size());
}

