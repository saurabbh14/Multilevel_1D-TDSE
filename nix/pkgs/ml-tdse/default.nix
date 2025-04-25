{ stdenv
, lib
, gfortran
, meson
, ninja
, pkg-config
, openblasCompat
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
    openblasCompat
    (lib.getLib fftw)
    (lib.getDev fftw)
  ];

  
  # Add a postInstall phase to ensure module files are copied
  postInstall = ''
    # Create include directory if it doesn't exist
    mkdir -p $out/include/ML-TDSE
    
    # Copy all mod files to the include directory
    find . -name "*.mod" -exec cp {} $out/include/ML-TDSE/ \;
  '';

  disableHardening = "all";

  meta.mainProgram = "ML-TDSE";
})
