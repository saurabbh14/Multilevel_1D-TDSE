blas = dependency('blas')
lapack = dependency('lapack')
omp = dependency('openmp')
fftw3 = dependency('fftw3')
fftw3_omp = declare_dependency(
  link_args: ['-lfftw3_omp', 
    '-fexternal-blas']
)

subdir('IO_modules')
subdir('pulse_gen')
subdir('fft')
subdir('performance_modules')
subdir('blas_module')
prog_sources = files('main.f90',
           'nuclear_wv.f90',
           'propagation.f90',
           'initializer.f90')

incdir = include_directories('IO_modules', 
                    'pulse_gen',
                    'fft')
executable('ML-TDSE',
            IO_mod_sources
            + fft_mod_sources
            + performance_mod_sources
            + blas_mod_sources
            + pulse_mod_sources 
            + prog_sources,
           dependencies: [ blas, lapack, omp, fftw3, fftw3_omp ],
           install : true,
           include_directories : incdir
           )
