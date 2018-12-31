# MPI-CRTM
Simple code for calculating WRF output files Brightness Temperature ( TBB ) observed by MTSAT-2.

By calling the Community Radiative Transfer Model ( CRTM https://www.jcsda.noaa.gov/projects_crtm.php )'s static libaray and Message Passing Interface ( Intel's MPI https://software.intel.com/en-us/mpi-library , Free of Charge because of Intel's Student discount ), excutables can be running in clusters, but the performance needs to be supercharged.

Community Gridpoint Statistical Interpolation ( GSI https://dtcenter.org/com-GSI/users/ ) system's cloud effective radius ( set_crtm_cloudmod.f90 ) is used in this module.

Example Figures:
