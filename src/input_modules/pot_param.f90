!> Parameters related to potentials and absorbers.
module pot_param
    use varprecision, only: dp
    real(dp):: R0     ! Grid-Parameter, start..
    real(dp)::Rend   !..and end
    real(dp),parameter:: cpmR=3.2d0*2*2 !*2 !absorber position from the end of R-grid
end module pot_param
