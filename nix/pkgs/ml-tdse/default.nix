{ stdenv
, lib
, gfortran
, meson
, ninja
, pkg-config
, blas
, lapack
, fftw
}:

stdenv.mkDerivation (finalAttrs: {
  pname = "ML-TDSE";
  version = "dev";

  src = lib.cleanSourceWith {
    filter = path: type: !(builtins.elem path [ "nix" ]);
    src = lib.cleanSource ../../../.;
  };

  nativeBuildInputs = [
    gfortran
    meson
    ninja
    pkg-config
  ];

  buildInputs = [
    blas
    lapack
    (lib.getLib fftw)
    (lib.getDev fftw)
  ];

  disableHardening = "all";

  meta.mainProgram = "ML-TDSE";
})
