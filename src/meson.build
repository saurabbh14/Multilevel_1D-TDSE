blas = dependency('blas')
lapack = dependency('lapack')
omp = dependency('openmp')
fftw3 = dependency('fftw3')
fftw3_omp = declare_dependency(
  link_args: ['-lfftw3_omp', 
    '-fexternal-blas']
)

subdir('input_modules')
prog_sources = files('main.f90',
           'nuclear_wv.f90',
           'propagation.f90',
           'initializer.f90')

incdir = include_directories('input_modules')
executable('ML-TDSE',
            mod_sources + prog_sources,
           dependencies: [ blas, lapack, omp, fftw3, fftw3_omp ],
           install : true,
           include_directories : incdir
           )
