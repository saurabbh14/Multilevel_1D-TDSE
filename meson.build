project('ML-TDSE', 'fortran',
  version : '0.1',
  default_options : ['warning_level=2',
    'optimization=0'
    ]
  )
add_project_arguments([ '-static',
                    '-fbounds-check', 
                    '-g',
                    '-fcheck=all',
                    '-fbacktrace'], language: 'fortran')
  
subdir('src')
