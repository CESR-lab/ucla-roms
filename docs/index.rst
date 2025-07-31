The Regional Ocean Modeling System (ROMS)
=========================================
ROMS is a free-surface, terrain following, primitive equation ocean model, developed by our collaborators at UCLA. UCLA-ROMS is highly optimized, supports state-of-the-art biogeochemical modelling, and is well-suited for nested regional modelling from ocean basins down to estuaries.

This fork aims to expand support for marine carbon dioxide removal (mCDR) research while adhering to community driven open-source practices.

Key technical information:
--------------------------

- Terrain-following S-coordinates (vertical) with parameterized stretching enhancing resolution near surface and bottom
- staggered Arakawa-C grid (horizontal) supporting spherical or Cartesian curvilinear coordinates
- Split-explicit time-stepping, separately advancing barotropic and baroclinic modes
- Third order upstream advection of momentum and tracers
- K-profile parameterization for vertical turbulence closure :cite:`large1994`, with optional alternative schemes
- Support for non-hydrostatic modelling via an included library :cite:`kanarska2007`
- Support for biogeochemical modelling via MARBL :cite:`long2021` or its predecessor, BEC :cite:`moore2004`


.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   Installing ROMS <installation/index>

.. toctree::
   :maxdepth: 1
   :caption: References

   References <references>
