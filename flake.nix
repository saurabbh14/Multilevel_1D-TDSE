{
  description = "Saurabh's Multilevel TDSE solver";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    nix-bundler.url = "github:nixos/bundlers";
  };

  outputs = { self, nixpkgs, flake-utils, nix-bundler }: flake-utils.lib.eachDefaultSystem (system:
    let
      pkgs = import nixpkgs {
        inherit system;
        overlays = [ (import ./nix/overlay.nix) ];
      };
    in
    {
      packages.default = pkgs.ml-tdse;

      devShell = with pkgs; mkShell {
        buildInputs = [ fortran-language-server ] ++
          ml-tdse.nativeBuildInputs ++
          ml-tdse.buildInputs;
      };

      formatter = pkgs.nixpkgs-fmt;

      bundlers = { inherit (nix-bundler.bundlers."${system}") toArx toDockerImage; };
    }
  );

}
