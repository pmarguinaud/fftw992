/* Minimal MKL compatibility header for fft992-based mkl library */
#ifndef _MKL_H_
#define _MKL_H_

/* This is a stub header providing compatibility with ecbuild FindMKL.cmake.
   The actual implementation uses fft992 (Temperton FFT) wrapped in an
   FFTW3-compatible interface. */

#define MKL_FFTW3_COMPAT 1

#endif /* _MKL_H_ */
