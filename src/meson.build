blas = dependency('blas')
lapack = dependency('lapack')
omp = dependency('openmp')
fftw3 = dependency('fftw3')
fftw3_omp = declare_dependency(
  link_args: ['-lfftw3_omp', 
    '-fexternal-blas']
)

subdir('IO_modules')
subdir('libs')
subdir('processes')
prog_sources = files('main.f90',
           'nuclear_wv.f90',
           'propagation.f90')

incdir = include_directories('IO_modules', 
                    'pulse_gen',
                    'libs',
                    'processes')
executable('ML-TDSE',
            IO_mod_sources
            + libs_sources
            + processes_sources,
           dependencies: [ blas, lapack, omp, fftw3, fftw3_omp ],
           install : true,
           include_directories : incdir
           )
