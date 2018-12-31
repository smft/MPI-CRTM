!!!!             Definition                !!!!
! target "!#" stands for chunk description    !
! target "!!" stands for block description    !
! target "! " stands for sentence description !
!           May God Bless Myself              !

!*******************************************************************************
subroutine calculate_temperature(theta,pressure,base_pressure,temperature,nz)
  use omp_lib
  implicit none
  real(kind=8), intent(in) :: theta(nz), pressure(nz), base_pressure(nz)
  real(kind=8), intent(inout) :: temperature(nz)
  integer, intent(in) :: nz
  integer :: i
  !$OMP DO
    do i=1,nz,1
      temperature(i)=(theta(i)+290)/((1000.0/((pressure(i)+base_pressure(i))/100))**(287.04/1003))
    enddo
  !$OMP END DO
end subroutine calculate_temperature
!*******************************************************************************

!*******************************************************************************
subroutine read_netcdf_3D(file_name,var_name,nx,ny,nz,data_)
  use omp_lib
  use netcdf
  implicit none
  character(len=100), intent(in) :: var_name,file_name
  real(kind=8), intent(inout) :: data_(nx*ny*nz)
  real(kind=8), allocatable, dimension(:,:,:,:) :: trans
  integer, intent(in) :: nx,ny,nz
  integer :: status_, i, j, k, count_, var_id, file_id
  status_=nf90_open(trim(file_name),nf90_nowrite,file_id)
  allocate(trans(nx,ny,nz,1))
  status_=nf90_inq_varid(file_id, trim(var_name), var_id)
  status_=nf90_get_var(file_id, var_id, trans, (/1,1,1,1/), (/nx,ny,nz,1/))
  status_=nf90_close(file_id)
  count_=1
  !$OMP DO
    do i=1,nx,1
      do j=1,ny,1
        !do k=1,nz,1
        data_(count_:count_+nz-1)=trans(i,j,:,1)
        count_=count_+nz
        !enddo
      enddo
    enddo
  !$OMP END DO
  deallocate(trans)
end subroutine read_netcdf_3D
!*******************************************************************************

!*******************************************************************************
subroutine read_netcdf_2D(file_name,var_name,nx,ny,data_)
  use omp_lib
  use netcdf
  implicit none
  character(len=100), intent(in) :: var_name,file_name
  real(kind=8), intent(inout) :: data_(nx*ny)
  real(kind=8), allocatable, dimension(:,:,:) :: trans
  integer, intent(in) :: nx,ny
  integer :: status_, i, j, count_, var_id, file_id
  status_=nf90_open(trim(file_name),nf90_nowrite,file_id)
  allocate(trans(nx,ny,1))
  status_=nf90_inq_varid(file_id, trim(var_name), var_id)
  status_=nf90_get_var(file_id, var_id, trans, (/1,1,1/), (/nx,ny,1/))
  status_=nf90_close(file_id)
  count_=1
  !$OMP DO
    do i=1,nx,1
      do j=1,ny,1
        data_(count_)=trans(i,j,1)
        count_=count_+1
      enddo
    enddo
  !$OMP END DO
  deallocate(trans)
end subroutine read_netcdf_2D
!*******************************************************************************

!*******************************************************************************
!*******************************************************************************
program test_mpi_scatter_bcast
  use mpi
  use netcdf
  use crtm_module

  implicit none
  !! define variables !!
  character (len=100) :: file_name, output_name
  integer :: local_mpi_err_flag, local_mpi_rank, local_mpi_size
  integer :: read_status, file_id, nz, ny, nx, nsoil, var_id
  integer :: i, j, k, count_, singular_size_to_add, normal_size, i_sendcounts, process_i
  real(kind=8) :: ptop
  real(kind=8), allocatable, dimension(:) :: lai, u10, v10, seaice, snowh, &
                                              vegfra, ivegtyp, xland, &
                                              xlat, xlong, coszen, &
                                              smois, stemp, tsk
  real(kind=8), allocatable, dimension(:) :: process_lai, process_u10, process_v10, &
                                              process_seaice, process_snowh, &
                                              process_vegfra, process_ivegtyp, & 
                                              process_xland, process_xlat, process_xlong, &
                                              process_coszen, process_smois, process_stemp, &
                                              process_tsk
  real(kind=8), allocatable, dimension(:) :: theta, pressure, base_pressure, &
                                              qvapor, qcloud, qrain, qice, &
                                              qsnow, temperature, tbb_1d
  real(kind=8), allocatable, dimension(:) :: process_theta,process_pressure,true_pressure,&
                                              process_base_pressure,process_temperature,&
                                              process_qvapor,process_qcloud,process_qrain,&
                                              process_qice,process_qsnow,process_qgraupel,&
                                              process_qhail,process_tbb,&
                                              eff_rad_cloud,eff_rad_rain,eff_rad_ice,&
                                              eff_rad_snow,eff_rad_graupel,eff_rad_hail
  real(kind=8), allocatable, dimension(:,:) :: tbb_2d
  !! start mpi !!
  call MPI_INIT(local_mpi_err_flag)
  call MPI_COMM_RANK(MPI_COMM_WORLD,local_mpi_rank,local_mpi_err_flag)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,local_mpi_size,local_mpi_err_flag)

  if (local_mpi_rank == 0) then
    !# acquire file name from namelist#!
    open(155,FILE="./namelist.crtm",FORM="FORMATTED",STATUS="OLD")
    read(155,*) file_name,output_name
    write(*,*) file_name,output_name
    !# reading netcdf variables in a serial mod #!
    ! open netcdf file in reading mod !
    read_status=nf90_open(trim(file_name),nf90_nowrite,file_id)
    ! get netcdf dimensions !
    read_status=nf90_inq_dimid(file_id,"south_north",var_id)
    read_status=nf90_inquire_dimension(file_id,var_id,len=ny)
    read_status=nf90_inq_dimid(file_id,"west_east",var_id)
    read_status=nf90_inquire_dimension(file_id,var_id,len=nx)
    read_status=nf90_inq_dimid(file_id,"bottom_top",var_id)
    read_status=nf90_inquire_dimension(file_id,var_id,len=nz)
    read_status=nf90_inq_dimid(file_id,"soil_layers_stag",var_id)
    read_status=nf90_inquire_dimension(file_id,var_id,len=nsoil)
    read_status=nf90_close(file_id)
    ! get the 3D data !
    allocate(theta(nx*ny*nz))
    call read_netcdf_3D(file_name,"T",nx,ny,nz,theta)
    allocate(pressure(nx*ny*nz))
    call read_netcdf_3D(file_name,"P",nx,ny,nz,pressure)
    allocate(base_pressure(nx*ny*nz))
    call read_netcdf_3D(file_name,"PB",nx,ny,nz,base_pressure)
    allocate(qvapor(nx*ny*nz))
    call read_netcdf_3D(file_name,"QVAPOR",nx,ny,nz,qvapor)
    allocate(qcloud(nx*ny*nz))
    call read_netcdf_3D(file_name,"QCLOUD",nx,ny,nz,qcloud)
    allocate(qrain(nx*ny*nz))
    call read_netcdf_3D(file_name,"QRAIN",nx,ny,nz,qrain)
    allocate(qice(nx*ny*nz))
    call read_netcdf_3D(file_name,"QICE",nx,ny,nz,qice)
    allocate(qsnow(nx*ny*nz))
    call read_netcdf_3D(file_name,"QSNOW",nx,ny,nz,qsnow)
    allocate(temperature(nx*ny*nz))
    ! get the 2D data !
    allocate(lai(nx*ny))
    call read_netcdf_2D(file_name,"LAI",nx,ny,lai)
    allocate(vegfra(nx*ny))
    call read_netcdf_2D(file_name,"VEGFRA",nx,ny,vegfra)
    allocate(ivegtyp(nx*ny))
    call read_netcdf_2D(file_name,"TSK",nx,ny,ivegtyp)
    allocate(u10(nx*ny))
    call read_netcdf_2D(file_name,"U10",nx,ny,u10)
    allocate(v10(nx*ny))
    call read_netcdf_2D(file_name,"V10",nx,ny,v10)
    allocate(seaice(nx*ny))
    call read_netcdf_2D(file_name,"SEAICE",nx,ny,seaice)
    allocate(snowh(nx*ny))
    call read_netcdf_2D(file_name,"SNOWH",nx,ny,snowh)
    allocate(tsk(nx*ny))
    call read_netcdf_2D(file_name,"TSK",nx,ny,tsk)
    allocate(smois(nx*ny*nsoil))
    call read_netcdf_3D(file_name,"SMOIS",nx,ny,nsoil,smois)
    allocate(stemp(nx*ny*nsoil))
    call read_netcdf_3D(file_name,"TSLB",nx,ny,nsoil,stemp)
    allocate(xland(nx*ny))
    call read_netcdf_2D(file_name,"XLAND",nx,ny,xland)
    allocate(xlat(nx*ny))
    call read_netcdf_2D(file_name,"XLAT",nx,ny,xlat)
    allocate(xlong(nx*ny))
    call read_netcdf_2D(file_name,"XLONG",nx,ny,xlong)
    allocate(coszen(nx*ny))
    call read_netcdf_2D(file_name,"COSZEN",nx,ny,coszen)
    allocate(tbb_1d(nx*ny))
    ! get ptop !
    ptop=5000.0
    write(*,*) "SUCCESS reading WRFOUT"
    ! get scatter size !
    i_sendcounts=((nx*ny)/local_mpi_size)*nz
    call MPI_Bcast(i_sendcounts,1,MPI_INTEGER,0,MPI_COMM_WORLD,local_mpi_err_flag)
    call MPI_Bcast(nz,1,MPI_INTEGER,0,MPI_COMM_WORLD,local_mpi_err_flag)
    call MPI_Bcast(nsoil,1,MPI_INTEGER,0,MPI_COMM_WORLD,local_mpi_err_flag)
    call MPI_Bcast(ptop,1,MPI_REAL8,0,MPI_COMM_WORLD,local_mpi_err_flag)
    write(*,*) "SUCCESS Broardcast Const on Processor ",local_mpi_rank
  else
    call MPI_Bcast(i_sendcounts,1,MPI_INTEGER,0,MPI_COMM_WORLD,local_mpi_err_flag)
    call MPI_Bcast(nz,1,MPI_INTEGER,0,MPI_COMM_WORLD,local_mpi_err_flag)
    call MPI_Bcast(nsoil,1,MPI_INTEGER,0,MPI_COMM_WORLD,local_mpi_err_flag)
    call MPI_Bcast(ptop,1,MPI_REAL8,0,MPI_COMM_WORLD,local_mpi_err_flag)
    write(*,*) "SUCCESS Broardcast Const on Processor ",local_mpi_rank
  endif
  ! scatter data gonna trasnmiting !
  allocate(process_theta(i_sendcounts))
  allocate(process_pressure(i_sendcounts))
  allocate(process_base_pressure(i_sendcounts))
  allocate(true_pressure(i_sendcounts))
  allocate(process_qvapor(i_sendcounts))
  allocate(process_qcloud(i_sendcounts))
  allocate(process_qrain(i_sendcounts))
  allocate(process_qice(i_sendcounts))
  allocate(process_qsnow(i_sendcounts))
  allocate(process_qgraupel(i_sendcounts))
  allocate(process_qhail(i_sendcounts))
  allocate(process_tbb(i_sendcounts))
  call MPI_Scatter(theta,i_sendcounts,MPI_REAL8,&
                    process_theta,i_sendcounts,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(pressure,i_sendcounts,MPI_REAL8,&
                    process_pressure,i_sendcounts,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(base_pressure,i_sendcounts,MPI_REAL8,&
                    process_base_pressure,i_sendcounts,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(qvapor,i_sendcounts,MPI_REAL8,&
                    process_qvapor,i_sendcounts,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(qcloud,i_sendcounts,MPI_REAL8,&
                    process_qcloud,i_sendcounts,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(qrain,i_sendcounts,MPI_REAL8,&
                    process_qrain,i_sendcounts,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(qice,i_sendcounts,MPI_REAL8,&
                    process_qice,i_sendcounts,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(qsnow,i_sendcounts,MPI_REAL8,&
                    process_qsnow,i_sendcounts,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  allocate(process_lai(i_sendcounts/nz))
  allocate(process_vegfra(i_sendcounts/nz))
  allocate(process_ivegtyp(i_sendcounts/nz))
  allocate(process_u10(i_sendcounts/nz))
  allocate(process_v10(i_sendcounts/nz))
  allocate(process_seaice(i_sendcounts/nz))
  allocate(process_snowh(i_sendcounts/nz))
  allocate(process_tsk(i_sendcounts/nz))
  allocate(process_smois((i_sendcounts/nz)*nsoil))
  allocate(process_stemp((i_sendcounts/nz)*nsoil))
  allocate(process_xland(i_sendcounts/nz))
  allocate(process_xlat(i_sendcounts/nz))
  allocate(process_xlong(i_sendcounts/nz))
  allocate(process_coszen(i_sendcounts/nz))
  call MPI_Scatter(lai,i_sendcounts/nz,MPI_REAL8,&
                    process_lai,i_sendcounts/nz,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(vegfra,i_sendcounts/nz,MPI_REAL8,&
                    process_vegfra,i_sendcounts/nz,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(ivegtyp,i_sendcounts/nz,MPI_REAL8,&
                    process_ivegtyp,i_sendcounts/nz,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(u10,i_sendcounts/nz,MPI_REAL8,&
                    process_u10,i_sendcounts/nz,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(v10,i_sendcounts/nz,MPI_REAL8,&
                    process_v10,i_sendcounts/nz,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(seaice,i_sendcounts/nz,MPI_REAL8,&
                    process_seaice,i_sendcounts/nz,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(snowh,i_sendcounts/nz,MPI_REAL8,&
                    process_snowh,i_sendcounts/nz,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(tsk,i_sendcounts/nz,MPI_REAL8,&
                    process_tsk,i_sendcounts/nz,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(smois,(i_sendcounts/nz)*nsoil,MPI_REAL8,&
                    process_smois,(i_sendcounts/nz)*nsoil,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(stemp,(i_sendcounts/nz)*nsoil,MPI_REAL8,&
                    process_stemp,(i_sendcounts/nz)*nsoil,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(xland,i_sendcounts/nz,MPI_REAL8,&
                    process_xland,i_sendcounts/nz,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(xlat,i_sendcounts/nz,MPI_REAL8,&
                    process_xlat,i_sendcounts/nz,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(xlong,i_sendcounts/nz,MPI_REAL8,&
                    process_xlong,i_sendcounts/nz,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Scatter(coszen,i_sendcounts/nz,MPI_REAL8,&
                    process_coszen,i_sendcounts/nz,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  write(*,*) "SUCCESS Scatter Array on Processor ",local_mpi_rank
  ! start calcluate pressure !
  do process_i=1,i_sendcounts,1
    true_pressure(process_i)=process_pressure(process_i)+process_base_pressure(process_i) 
    process_theta(process_i)=process_theta(process_i)+290.0
  enddo 
  write(*,*) "SUCCESS Calculate Pressure on Processor ",local_mpi_rank
  
  do process_i=1,i_sendcounts/nz,1
  !do process_i=1,10,1
     ! nlevels,theta, qvapor, pressure, &
     ! qcloud, qsnow, qgraupel, qice, qrain, &
     ! qhail, lai, u10, v10, seaice, &
     ! snowh, vegfrac, tsk, ivegtyp, smois, stemp, &
     ! xland, lat, lon, ptop, coszen, channel, tbb
     !write(*,*) process_theta((process_i-1)*nz+1:process_i*nz)
     !write(*,*) process_qvapor((process_i-1)*nz+1:process_i*nz)
     !write(*,*) true_pressure((process_i-1)*nz+1:process_i*nz)
     !write(*,*) process_qcloud((process_i-1)*nz+1:process_i*nz)
     !write(*,*) process_qsnow((process_i-1)*nz+1:process_i*nz)
     !write(*,*) process_qgraupel((process_i-1)*nz+1:process_i*nz)
     !write(*,*) process_qice((process_i-1)*nz+1:process_i*nz)
     !write(*,*) process_qrain((process_i-1)*nz+1:process_i*nz)
     !write(*,*) process_qhail((process_i-1)*nz+1:process_i*nz)
     call calculate_tbb(nz,&
                          process_theta((process_i-1)*nz+1:process_i*nz),&
                          process_qvapor((process_i-1)*nz+1:process_i*nz),&
                          true_pressure((process_i-1)*nz+1:process_i*nz), &
                          process_qcloud((process_i-1)*nz+1:process_i*nz),&
                          process_qsnow((process_i-1)*nz+1:process_i*nz),&
                          process_qgraupel((process_i-1)*nz+1:process_i*nz),&
                          process_qice((process_i-1)*nz+1:process_i*nz),&
                          process_qrain((process_i-1)*nz+1:process_i*nz),&
                          process_qhail((process_i-1)*nz+1:process_i*nz),&
                          process_lai(process_i),&
                          process_u10(process_i),&
                          process_v10(process_i),&
                          process_seaice(process_i),&
                          process_snowh(process_i),&
                          process_vegfra(process_i),&
                          process_tsk(process_i),&
                          process_ivegtyp(process_i),&
                          process_smois((process_i-1)*nsoil+1),&
                          process_stemp((process_i-1)*nsoil+1),&
                          process_xland(process_i),&
                          process_xlat(process_i),&
                          process_xlong(process_i),&
                          ptop,&
                          process_coszen(process_i),&
                          4,&
                          process_tbb(process_i))
  enddo 
  call MPI_Barrier(MPI_COMM_WORLD,local_mpi_err_flag)
  call MPI_Gather(process_tbb,i_sendcounts/nz,MPI_REAL8,&
                    tbb_1d,i_sendcounts/nz,MPI_REAL8,&
                    0,MPI_COMM_WORLD,local_mpi_err_flag)
  if (local_mpi_rank == 0) then
    allocate(tbb_2d(nx,ny))
    count_=1
    do i=1,nx,1
      do j=1,ny,1
        tbb_2d(i,j)=tbb_1d(count_)
        count_=count_+1
      enddo
    enddo
    write(*,*) count_,nx,ny
    open(156,file=trim(output_name),form="unformatted")
    write(156) tbb_2d
    close(156)
  endif
  !! end mpi !!
  call MPI_FINALIZE(local_mpi_err_flag)
end program
!*******************************************************************************
!*******************************************************************************
