blas = dependency('blas')
lapack = dependency('lapack')
omp = dependency('openmp')
fftw3 = dependency('fftw3')
fftw3_omp = declare_dependency(
  link_args: ['-lfftw3_omp' ]
)

executable('ml-tdse',
           files('main.f90',
           'adiabatic.f90',
           'nuclear_wv.f90',
           'propagation.f90',
           'pulse.f90',
           'setpot.f90'),
           dependencies: [ blas, lapack, omp, fftw3, fftw3_omp ],
           install : true
           )