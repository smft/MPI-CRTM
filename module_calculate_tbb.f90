!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculating Geostationary Satellite Tbb
! Using CRTM Radiative Model
! Qi Zhang
! dg1628023@smail.nju.edu.cn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this is semi-full aware version !
! COPYED FROM GSI :} !
SUBROUTINE effective_radius_semi(nuclei_name,nz,theta,pressure,&
                                 qvapor,q_nuclei,eff_rad)
  implicit none
  character(len=20), intent(in) :: nuclei_name
  integer, intent(in) :: nz
  integer :: k,j,i
  real(kind=8), intent(in) :: theta(nz), pressure(nz),&
                              qvapor(nz), q_nuclei(nz)
  real(kind=8), intent(inout) :: eff_rad(nz)
  real(kind=8) :: tem4, indexw, indexi, indexr, temp, rho, rhos, &
                  pi, rhog, rhoh
  pi=3.141592654
  indexi=50.0
  indexr=350.0
  rhos=150.0
  rhog=400.0
  rhoh=900.0
  do k=1,nz,1
    temp=theta(k)*(pressure(k)/100000.)**(287.04/1003)
    rho=pressure(k)/(278.04*temp*(1+0.608*qvapor(k)))
    if ((273.15-temp)>0) then
      tem4=(273.15-temp)*0.05
    else
      tem4=0
    endif
    if(trim(nuclei_name)=="cloud") then
      if(tem4>1) then
        eff_rad(k)=10.0
      else
        eff_rad(k)=5+5*tem4
      endif
    else if(trim(nuclei_name)=="ice") then
      eff_rad(k)=1.5*indexi
    else if(trim(nuclei_name)=="rain") then
      eff_rad(k)=1.5*indexr
    else if(trim(nuclei_name)=="snow") then
      eff_rad(k)=1.5*(rho*q_nuclei(k)/(pi*rhos*20.0e3))**(1.0/3.0)*1.0e6
    else if(trim(nuclei_name)=="graupel") then
      eff_rad(k)=1.5*(rho*q_nuclei(k)/(pi*rhog*20.0e3))**(1.0/3.0)*1.0e6
    else if(trim(nuclei_name)=="hail") then
      eff_rad(k)=1.5*(rho*q_nuclei(k)/(pi*rhoh*20.0e3))**(1.0/3.0)*1.0e6
    endif
  end do
END SUBROUTINE effective_radius_semi

! this is the NASA GMAO Simple version :) !
! Static Effective Radius !
SUBROUTINE effective_radius_simple(nz,nuclei_name,q_nuclei,eff_rad)
  implicit none
  character(len=20), intent(in) :: nuclei_name
  integer, intent(in) :: nz
  integer :: k
  real(kind=8), intent(in) :: q_nuclei(nz)
  real(kind=8), intent(inout) :: eff_rad(nz)
  do k=1,nz,1
    if (q_nuclei(k)>0) then
      if ( trim(nuclei_name) == 'cloud' ) then
        eff_rad(k)=10.0
      else if ( trim(nuclei_name) == 'ice' ) then
        eff_rad(k)=30.0
      else if ( trim(nuclei_name) == 'hail' ) then
        eff_rad(k)=1200.0
      else if ( trim(nuclei_name) == 'graupel' ) then
        eff_rad(k)=800.0
      else if ( trim(nuclei_name) == 'rain' ) then
        eff_rad(k)=400.0
      else if ( trim(nuclei_name) == 'snow' ) then
        eff_rad(k)=600.0
      endif
    else
      eff_rad(k)=0.0
    endif
  end do
END SUBROUTINE effective_radius_simple 

SUBROUTINE calculate_tbb(nlevels,theta, qvapor, pressure, &
                          qcloud, qsnow, qgraupel, qice, qrain, &
                          qhail, lai, u10, v10, seaice, &
                          snowh, vegfrac, tsk, ivegtyp, smois, stemp, &
                          xland, lat, lon, ptop, coszen, channel, tbb)
  !NEED temperature, qvapor, pressure, qcloud, qsnow, qgraupel, qice, qrain, \
  !qhail, ercloud, erice, errain, ersnow, ergraupel, erhail for ATMOSPHERE
  !NEED lai, u10, v10, seaice, snowh, vegfrac, vegfrac, tsk, ivegtyp, xland, \
  !landuse, for SURFACE
  !ptop, mp_physics, lat, lon, nlevels for GEO
  !sensor, channel, coeff_path, request_var, crtm_out are mannually defined
  !Only calculating Single Profile, No need to call a 3D array
  use crtm_module
  implicit none
  !define input and output variables
  integer, intent(in) :: nlevels, channel
  real(kind = 8), intent(in) :: theta(nlevels), &
                                qvapor(nlevels), &
                                pressure(nlevels), &
                                qcloud(nlevels), &
                                qsnow(nlevels), &
                                qgraupel(nlevels), &
                                qice(nlevels), &
                                qrain(nlevels), &
                                qhail(nlevels)
  real(kind = 8) :: eff_rad_cloud(nlevels), &
                          eff_rad_ice(nlevels), &
                          eff_rad_rain(nlevels), &
                          eff_rad_snow(nlevels), &
                          eff_rad_graupel(nlevels), &
                          eff_rad_hail(nlevels)
  real(kind = 8), intent(in) :: lai, u10, v10, &
                                seaice, snowh, &
                                vegfrac, tsk, &
                                xland, lat, lon, &
                                ptop, coszen, &
                                smois, stemp, ivegtyp
  real(kind = 8), intent(inout) :: tbb
  character(len = 20) :: landuse="USGS", sensor="imgr_mt2"
  character(len = 250) :: coeff_path="/home/qzhang/CRTM/bin/"
  character(len=50) :: request_var="Brightness Temperature"
  character(len=50) :: nuclei_name
  ! define crtm lib use variables
  real :: sat_lat, sat_lon, rho
  integer :: n_channels,n_clouds,n_layers,nz
  integer, parameter :: n_sensors = 1, n_profiles = 1, &
                        n_absorbers = 2, n_aerosols = 0
  ! H20 and O3
  ! only single profiel
  character(len = 20) :: sensor_id(n_sensors)
  type(crtm_channelinfo_type) :: chinfo(n_sensors)
  type(crtm_geometry_type) :: geo(n_sensors)
  type(crtm_options_type) :: opt(n_sensors)
  type(crtm_atmosphere_type) :: atm(n_sensors)
  type(crtm_surface_type) :: sfc(n_sensors)
  type(crtm_rtsolution_type), allocatable :: rts(:,:)
  real :: o3, veg_coverage
  real :: lvl_pressure(nlevels+1), lay_pressure(nlevels),temp(nlevels)

  real, allocatable :: r_eff(:)
  integer :: alloc_stat
  integer :: err_stat
  integer :: k,j,i, vegtyp
  integer :: usgs_to_npoess(24), igbp_to_npoess(20)
  real :: sat_zenith, r_sat, vzenith
  real :: sun_zenith, sun_azimuth
  real :: land_coverage, water_coverage, snow_coverage, ice_coverage
  real, parameter :: rd = 287.05
  real, parameter :: cp = 1003.
  real, parameter :: r2d = 180./3.141592654
  real :: mass_lay
  nz = nlevels
  n_layers = nlevels
  n_clouds = 6
  nuclei_name="cloud"
  CALL effective_radius_semi(nuclei_name,nz,theta,pressure,&
                             qvapor,qcloud,eff_rad_cloud)
  nuclei_name="ice"
  CALL effective_radius_semi(nuclei_name,nz,theta,pressure,&
                             qvapor,qice,eff_rad_ice)
  nuclei_name="rain"
  CALL effective_radius_semi(nuclei_name,nz,theta,pressure,&
                             qvapor,qrain,eff_rad_rain)
  nuclei_name="snow"
  CALL effective_radius_semi(nuclei_name,nz,theta,pressure,&
                             qvapor,qsnow,eff_rad_snow)
  do k=1,nz,1
    eff_rad_hail(k)=0
    eff_rad_graupel(k)=0
  enddo
  ! start CRTM module
  ! initialize crtm
  sensor_id(1)=sensor
  err_stat = crtm_init(sensor_id, chinfo, &
                       Load_CloudCoeff = .TRUE., &
                       Load_AerosolCoeff = .TRUE., &
                       IRlandCoeff_File = &
                       'NPOESS.IRland.EmisCoeff.bin', &
                       IRwaterCoeff_File = &
                       'Nalli.IRwater.EmisCoeff.bin', &
                       MWwaterCoeff_File = &
                       'FASTEM5.MWwater.EmisCoeff.bin', &
                       VISlandCoeff_File = &
                       'NPOESS.VISland.EmisCoeff.bin', &
                       File_Path = trim(coeff_path), &
                       Quiet = .True.)
  err_stat = crtm_channelinfo_subset(chinfo(1), Channel_Subset = (/channel/))
  CALL crtm_atmosphere_create(atm,n_layers,n_absorbers,n_clouds,n_aerosols)
  atm(1)%n_layers = n_layers
  atm(1)%absorber_id(1) = H2O_ID
  atm(1)%absorber_id(2) = O3_ID
  atm(1)%absorber_units(1) = MASS_MIXING_RATIO_UNITS
  atm(1)%absorber_units(2) = VOLUME_MIXING_RATIO_UNITS
  atm(1)%cloud(1)%n_layers = n_layers
  atm(1)%cloud(1)%Type = WATER_CLOUD
  atm(1)%cloud(2)%n_layers = n_layers
  atm(1)%cloud(2)%Type = ICE_CLOUD
  atm(1)%cloud(3)%n_layers = n_layers
  atm(1)%cloud(3)%Type = RAIN_CLOUD
  atm(1)%cloud(4)%n_layers = n_layers
  atm(1)%cloud(4)%Type = SNOW_CLOUD
  atm(1)%cloud(5)%n_layers = n_layers
  atm(1)%cloud(5)%Type = GRAUPEL_CLOUD
  atm(1)%cloud(6)%n_layers = n_layers
  atm(1)%cloud(6)%Type = HAIL_CLOUD
  atm(1)%Level_Pressure(0) = TOA_PRESSURE
  do k=1,nz,1
    lay_pressure(k) = pressure(k)
    temp(k) = theta(k) * (lay_pressure(k)/100000.)**(rd/cp)
    ! Will define level 1 later
    if (k > 1) then
      lvl_pressure(k) = (lay_pressure(k) + lay_pressure(k-1)) * 0.5
    endif
    ! Make sure pressure is decreasing
    if (k > 1 .and. lay_pressure(k) > lay_pressure(k-1)) then
      lay_pressure(k-1) = lay_pressure(k) * 1.0001
    endif
    if (k > 2 .and. lvl_pressure(k) > lvl_pressure(k-1)) then
      lvl_pressure(k-1) = lvl_pressure(k) * 1.0001
    endif
  end do
    lvl_pressure(nlevels+1)=ptop
    lvl_pressure(1)=lay_pressure(1)+(lvl_pressure(2)-lay_pressure(2))
    n_channels = crtm_channelinfo_n_channels(chinfo(1))
    allocate(rts(n_channels,n_profiles))
    CALL crtm_surface_create(sfc,chinfo(1)%n_channels)
    CALL crtm_rtsolution_create(rts,n_layers)
    geo%sensor_zenith_angle = 0.
    geo%sensor_scan_angle = 0.
    geo%source_zenith_angle = acos(coszen) * r2d
    sfc(1)%sensordata%n_channels = chinfo(1)%n_channels
    sfc(1)%sensordata%sensor_id = chinfo(1)%sensor_id
    sfc(1)%sensordata%WMO_sensor_id = chinfo(1)%WMO_sensor_id
    sfc(1)%sensordata%WMO_Satellite_id = chinfo(1)%WMO_Satellite_id
    sfc(1)%sensordata%sensor_channel = chinfo(1)%sensor_channel
    if (xland > 1.5) then
      veg_coverage = 0.0
      if (snowh > 0.1) then
        snow_coverage = 1.0
        land_coverage = 0.0
        ice_coverage = 0.0
        water_coverage = 0.0
      else
        snow_coverage = 0.0
        land_coverage = 0.0
        ice_coverage = 0.0
        water_coverage = 1.0
      endif
    else
      veg_coverage = vegfrac/100.
      if (seaice > 0) then
        if (snowh > 0.1) then
          snow_coverage = 1.0
          land_coverage = 0.0
          ice_coverage = 0.0
          water_coverage = 0.0
        else
          snow_coverage = 0.0
          land_coverage = 0.0
          ice_coverage = 1.0
          water_coverage = 0.0
        endif
      else
        if (snowh > 0.1) then
          snow_coverage = 1.0
          land_coverage = 0.0
          ice_coverage = 0.0
          water_coverage = 0.0
        else
          snow_coverage = 0.0
          land_coverage = 1.0
          ice_coverage = 0.0
          water_coverage = 0.0
        endif
      endif
    endif
    sfc(1)%Land_Coverage = land_coverage
    sfc(1)%Water_Coverage = water_coverage
    sfc(1)%Snow_Coverage = snow_coverage
    sfc(1)%Ice_Coverage = ice_coverage
    if (trim(landuse) .eq. 'USGS') then
      usgs_to_npoess = (/15,  1,  5, 11,  6,  6,  6,  7, 13, &
                        6,  8,  9,  8,  9, 12,  1, 18, 18, &
                        5, 10, 10, 10, 10,  1 /)
      vegtyp = usgs_to_npoess(min(max(1.,ivegtyp),24.))
    else
      igbp_to_npoess = (/ 9,  8,  9,  8, 12,  7, 19, 17, 17, &
                         7, 17,  2, 15,  2,  1,  1,  1, 10, &
                         10, 10  /)
      vegtyp = igbp_to_npoess(min(max(1.,ivegtyp),20.))
    endif
    sfc(1)%Land_Type = vegtyp
    sfc(1)%Wind_Speed = sqrt(u10*u10+v10*v10)
    sfc(1)%Land_Temperature = tsk
    sfc(1)%Snow_temperature = min(tsk,273.15)
    sfc(1)%Water_temperature = max(tsk,273.15)
    sfc(1)%Ice_temperature = min(tsk,273.15)
    ! soil_moist(i,j,1) !Get from WRF in g cm-3
    sfc(1)%Soil_Moisture_Content = smois
    sfc(1)%Canopy_Water_Content = smois
    sfc(1)%Vegetation_Fraction = veg_coverage
    ! soil_temp(i,j,1) !Default degree K
    sfc(1)%Soil_Temperature = stemp
    ! Get from WRF in mm
    sfc(1)%Snow_Depth = snowh*1000.
    ! From WRF in m2/m2
    sfc(1)%LAI = lai
  do k = 1, nz, 1
    atm(1)%Level_Pressure(k) = lvl_pressure(n_layers+1-k)
    atm(1)%Pressure(k) = lay_pressure(n_layers+1-k)
    atm(1)%Temperature(k) = temp(n_layers+1-k)
    ! H2O g/kg
    atm(1)%Absorber(k,1) = max(0.,qvapor(n_layers+1-k)*1000.)
    atm(1)%Absorber(k,2) = o3
    ! Need water content in kg_water/m2 --> need kg_air/m2 factor
    rho=atm(1)%Level_Pressure(k)/(278.04*atm(1)%Temperature(k)*&
        (1+0.608*atm(1)%Absorber(k,1)/1000))
    mass_lay = (atm(1)%Level_Pressure(k)-atm(1)%Level_Pressure(k-1))/(rho*9.81)
    !write(*,*) rho,mass_lay
    ! calculate effective radius, pre_calculated with NASA formula (GSI)
    atm(1)%Cloud(1)%Water_Content(k) = max(0.,qcloud(n_layers+1-k)*mass_lay*rho)
    atm(1)%Cloud(2)%Water_Content(k) = max(0.,qice(n_layers+1-k)*mass_lay*rho)
    atm(1)%Cloud(3)%Water_Content(k) = max(0.,qrain(n_layers+1-k)*mass_lay*rho)
    atm(1)%Cloud(4)%Water_Content(k) = max(0.,qsnow(n_layers+1-k)*mass_lay*rho)
    atm(1)%Cloud(5)%Water_Content(k) = max(0.,qgraupel(n_layers+1-k)*mass_lay*rho)
    atm(1)%Cloud(6)%Water_Content(k) = max(0.,qhail(n_layers+1-k)*mass_lay*rho)
    atm(1)%Cloud(1)%Effective_Radius(k) = max(0.,eff_rad_cloud(n_layers+1-k))
    atm(1)%Cloud(2)%Effective_Radius(k) = max(0.,eff_rad_ice(n_layers+1-k))
    atm(1)%Cloud(3)%Effective_Radius(k) = max(0.,eff_rad_rain(n_layers+1-k))
    atm(1)%Cloud(4)%Effective_Radius(k) = max(0.,eff_rad_snow(n_layers+1-k))
    atm(1)%Cloud(5)%Effective_Radius(k) = max(0.,eff_rad_graupel(n_layers+1-k))
    atm(1)%Cloud(6)%Effective_Radius(k) = max(0.,eff_rad_hail(n_layers+1-k))
  end do
    ! start calculate tbb !
    !CALL crtm_surface_zero(sfc)
    err_stat = CRTM_Forward(atm, sfc, geo, chinfo(1:1), rts)
    if (trim(request_var) == 'Radiance') then
      tbb = rts(1,1)%Radiance
    else if (trim(request_var) == 'Brightness Temperature') then
      tbb = rts(1,1)%Brightness_Temperature
    else if (trim(request_var) == 'Up Radiance') then
      tbb = rts(1,1)%Up_Radiance
    else if (trim(request_var) == 'Down Radiance') then
      tbb = rts(1,1)%Down_Radiance
    else if (trim(request_var) == 'Down Solar Radiance') then
      tbb = rts(1,1)%Down_Solar_Radiance
    endif
    deallocate(rts,STAT = alloc_stat )
  !write(*,*) tbb
  CALL crtm_atmosphere_destroy(atm)
  CALL crtm_surface_destroy(sfc)
  err_stat = crtm_destroy(chinfo)
END SUBROUTINE calculate_tbb

