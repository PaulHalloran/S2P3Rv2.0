!***************************************************************************************
!
!                                    S2P3 v7.0
!
!------------------Shelf Sea Physics and Primary Production v7.0------------------------
!
!   1-D MODEL OF THE EQUATION OF MOTION USING THE Canuto k-e TURBULENCE CLOSURE SCHEME
! with a simple single-species model of carbon fixation in response to light and nitrate
!
!***************************************************************************************
!
!                                Jonathan Sharples
!           University of Liverpool and NERC National Oceanography Centre
!
!                      In: J. H. Simpson and J. Sharples
!       Introduction to the Physical and Biological Oceanography of Shelf Seas
!                       Cambridge University Press 2012
!
!***************************************************************************************
!
!   Modified for regional application, compiled and executed in unix environments
!     (Marsh, Hickman and Sharples, Geoscientific Model Development, submitted)
!
!***************************************************************************************!
!


MODULE physics
  implicit none
  double precision :: lon,lat,kb,f0,f1,lambda,heat_shade,vismax,qsdt,stressx,stressy,windspeed,polarisn(6),rad_out, &
  dewT,stresss,stressb,Px,Py,slopex(5),slopey(5),tsteps,orient(5),semi_major(5),semi_minor(5), &
  first_temp,sx,sy,radsum,omega(5),major_phase(5),background_mol,u_bed,v_bed,speed_mag(200),mean_speed, &
  rad0,lw_rad_down0,speed,speed1,speed2,speed3,time,depth,newdepth,grav_accel,radiation_zerocloud(40000), &
  radiation_fullcloud(40000), &
  uac,uc,gac1,gac2,gac,gc1,gc2,gc,strat_jul,rad_sum,acount
  DOUBLE PRECISION :: wind_speed(40000),wind_dir(40000),std_dir(40000),std_speed(40000),&
  radiation(40000),humid(40000),airP(40000),airT(40000),rad_input(40000),lw_rad_down(40000), &
  cloud(40000),declination(40000)
  real :: tid,tid2,uamp(5),vamp(5),vphase(5),uphase(5)
  integer :: idmet,iday_front1,iday_front2,iday_front3,iday_front4
  data omega /1.405278d-4,1.454440d-4,1.378616d-4,6.759596d-5,7.293472e-5/
END MODULE
MODULE variables_all
 implicit none
 DOUBLE PRECISION :: c1,c2,c3,rad,time_step,dz,newdz
 double precision :: cnpar=0.5
 INTEGER :: N,newN,i,dataplot,ikey,icont,iquit,message_to_quit
 INTEGER :: imode
END MODULE
MODULE turbulence
 implicit none
 double precision :: NN(0:200),SS(0:200),lscale(0:200),h(200),NN_mean(0:200),SS_mean(0:200),DD_mean(0:200),RN_mean(0:200)
 double precision :: Kz(0:200),Nz(0:200),sm(0:200),KK_mean(0:200)
 double precision :: sh(0:200),Ri(0:200),P(0:200),B(0:200)
 double precision :: cmue1(0:200),cmue2(0:200),as(0:200),an(0:200)
 double precision, allocatable :: au(:),bu(:),cu(:),du(:),ru(:),qu(:),tke(:),tkeold(:),eps(:)
 double precision :: kappa=0.41
 double precision :: Kz_bg,Nz_bg				! Used in simple IW model
END MODULE
MODULE biology
  implicit none
  double precision :: x_old(200),x_new(200),ni_old(200),ni_new(200),quo(200), &
  net_grow(200),gross_grow(200),prod_net(200),prod_gross(200),prod_net_daily(200),prod_gross_daily(200),uptake(200),g0(40000)
  double precision :: sub_quota,alpha,alpha_persec,chl_carbon,max_quota, &
    Nrecycle,Pmax10,Pmax10_persec,Q10,gT0,uptake_max,uptake_max_persec, &
    half_saturation,graze_thresh,respiration,respiration_persec,rT0,rQ10, &
    swim_speed,swim_speed_persec,sink_speed, &
    sink_speed_persec,graze_min,graze_min_persec,graze_max,graze_max_persec,graze_daymax, &
    seed_biomass,chl_abscross,x_new_max,temp_x_new_max
  DOUBLE PRECISION :: s_old(200),s_new(200),din_source(200),radbar(200),par_percent, &
  din_rate,bed_din,par_atten,par_shade,graz(200),Pmax,Pm
  integer :: growth_model
  integer, parameter :: eppley=1
  integer, parameter :: Qten=2
END MODULE
MODULE file_names
  implicit none
  character(len=100) :: metfile_in,surface_out,phys_profile_day,phys_profile_hour, &
  bio_profile_day,bio_profile_hour,initialdata,monthly_out,user_metfile
END MODULE
! --------------------------------------------------------------------
! MODULE  MyTrigonometricFunctions:
!    This module provides the following functions and constants
!    (1) RadianToDegree()     - converts its argument in radian to
!                               degree
!    (2) DegreeToRadian()     - converts its argument in degree to
!                               radian
!    (3) MySIN()              - compute the sine of its argument in
!                               degree
!    (4) MyCOS()              - compute the cosine of its argument
!                               in degree
! --------------------------------------------------------------------

MODULE  MyTrigonometricFunctions
  implicit none
   REAL(8), PARAMETER :: PI        = 3.1415926
   REAL(8), PARAMETER :: Degree180 = 180.0
   REAL(8), PARAMETER :: R_to_D    = Degree180/PI
   REAL(8), PARAMETER :: D_to_R    = PI/Degree180

CONTAINS

! --------------------------------------------------------------------
! FUNCTION  RadianToDegree():
!    This function takes a REAL argument in radian and converts it to
! the equivalent degree.
! --------------------------------------------------------------------

   REAL FUNCTION  RadianToDegree(Radian)
      IMPLICIT  NONE
      REAL(8), INTENT(IN) :: Radian

      RadianToDegree = Radian * R_to_D
   END FUNCTION  RadianToDegree

! --------------------------------------------------------------------
! FUNCTION  DegreeToRadian():
!    This function takes a REAL argument in degree and converts it to
! the equivalent radian.
! --------------------------------------------------------------------

   REAL FUNCTION  DegreeToRadian(Degree)
      IMPLICIT  NONE
      REAL(8), INTENT(IN) :: Degree

      DegreeToRadian = Degree * D_to_R
   END FUNCTION  DegreeToRadian

! --------------------------------------------------------------------
! FUNCTION  MySIN():
!    This function takes a REAL argument in degree and computes its
! sine value.  It does the computation by converting its argument to
! radian and uses Fortran's sin().
! --------------------------------------------------------------------

   REAL FUNCTION  MySIN(x)
      IMPLICIT  NONE
      REAL(8), INTENT(IN) :: x

      MySIN = SIN(DegreeToRadian(x))
   END FUNCTION  MySIN

! --------------------------------------------------------------------
! FUNCTION  MySIN():
!    This function takes a REAL argument in degree and computes its
! cosine value.  It does the computation by converting its argument to
! radian and uses Fortran's cos().
! --------------------------------------------------------------------

   REAL FUNCTION  MyCOS(x)
      IMPLICIT  NONE
      REAL(8), INTENT(IN) :: x

      MyCOS = COS(DegreeToRadian(x))
   END FUNCTION  MyCOS

END MODULE  MyTrigonometricFunctions
!
MODULE VARIABLES
!
!   Shared variables for any routine with 'USE VARIABLES'
!
          IMPLICIT NONE
!
          INTEGER, PARAMETER            :: MAX_CHILD = 20  ! Max child windows
          INTEGER                       :: ICHILD    =  0  ! No. of children
          INTEGER, DIMENSION(MAX_CHILD) :: CHILDREN        ! Child handles
          LOGICAL                       :: PLOTT = .FALSE. ! Graphic plotted?
!
      DOUBLE PRECISION :: velx_old(200),velx_new(200),vely_old(200),  &
                          vely_new(200),ri_mean(200),Kz_mean(200), &
                          tke_mean(200), diss_mean(200),rad_mean(200), &
                          temp_old(200),temp_new(200),density(200), &
                          vel_mean(200),u_mean(200),v_mean(200),c1flux_mean(200), &
                          n1flux_mean(200),sflux_mean(200),hourly_net, &
                          hourly_gross,daily_net,daily_gross,prod_net2, &
                          prod_gross2,tpn1,tpn2,tpg1,tpg2,total_flux, &
                          surf_gross1,surf_net1,total_uptake,surf_uptake
real :: height(200),xtotal(200),wind_distribution(100)
INTEGER :: n_rad_choice,nhr_out,ndays,infile_error,iyear,met_type,metfile_error
integer :: output_start,output_end, output_hr_start, output_hr_end,ibmp(5),irun

      END MODULE VARIABLES
!
MODULE GRAPHICS_VARIABLES
  implicit none
  REAL, allocatable :: contour_data(:,:),xdata(:),ydata(:,:),profile(:,:)
  REAL :: r_height(200),r_dz,parmax,velmax,wmax
  real :: hu3a(5000),h_hu3(5000),t_hu3(5000),x1_hu3(5000),n1_hu3(5000),u1_hu3(5000),g1_hu3(5000), &
  x2_hu3(5000),n2_hu3(5000),u2_hu3(5000),g2_hu3(5000),s_hu3(5000),hkz_hu3(5000),Kz_hu3(5000), &
  hu3b(5000),sf_hu3(5000),c1f_hu3(5000),n1f_hu3(5000)
!
  real :: originx(15),originy(15),startx(15),starty(15),deltax(15),deltay(15),xincr(15),yincr(15), &
  endx(15),endy(15),xlen(15),ylen(15),ticklen,temp_min,temp_max,zmin,zmax1,zmax2,z_scale
  real :: charhx=0.006
  real :: charhy=0.0225
  INTEGER ncolour(15),ig(15),nticksx(15),nticksy(15),irunplot,screenx,screeny,irenew,black,white,JD_front,ifront
  character (len=15) :: labelx(15),labely(15),labelz(15),last_x(12),last_y(12)
  integer red(29),green(29),blue(29)
  data red /193,167,138,106,72,36,0,0,0,0,0,0,0,0,0,0,0,36,72,106,138,167,193,214,231,244,252, &
  255,252/
  data green /0,0,0,0,0,0,0,36,72,106,138,167,193,214,231,244,252,255,252,244,231,214,193,167,138, &
  106,72,36,0/
  data blue /138,167,193,214,231,244,252,255,252,244,231,214,193,167,138,106,72,36,0,0,0, &
  0,0,0,0,0,0,0,0/
END MODULE GRAPHICS_VARIABLES

!
!*****************************************************************************
!
      MODULE INTERFACES
          IMPLICIT NONE
          INTERFACE
              SUBROUTINE ProcessMenu2(IDENT,QUIT,unique_job_id)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN    ) :: IDENT
                  LOGICAL, INTENT(IN OUT) :: QUIT
                  character(len=36) :: unique_job_id
              END SUBROUTINE ProcessMenu2
!
              SUBROUTINE About()
                  IMPLICIT NONE
              END SUBROUTINE About
!
              SUBROUTINE Redraw()
                  IMPLICIT NONE
              END SUBROUTINE Redraw
!
              SUBROUTINE Plotit(IPLT)
                  IMPLICIT NONE
                  INTEGER, INTENT (IN    ) :: IPLT
              END SUBROUTINE Plotit
              SUBROUTINE WExit(QUIT,IWIN)
                  IMPLICIT NONE
                  LOGICAL,           INTENT (IN OUT) :: QUIT
                  INTEGER, OPTIONAL, INTENT (IN    ) :: IWIN
              END SUBROUTINE WExit
          END INTERFACE
      END MODULE INTERFACES
!
!*****************************************************************************
!
      PROGRAM phytowin
!
!   Program generated by WiDE Wizard at 09:17 on 06 Jun 2001.
!
      USE INTERFACES
      use turbulence
      use variables_all
      use variables
      use physics
      use biology
      use graphics_variables
      use file_names
!
      IMPLICIT NONE
!
!   Variable declarations
!
      LOGICAL           :: QUIT  = .FALSE.
      INTEGER           :: ITYPE
      INTEGER           :: IDENT
      integer :: run_year,start_year,iline
      character(len=12) :: lat_in_domain
      character(len=12) :: lon_in_domain
      character(len=36) :: unique_job_id
      character(len=300) :: met_data_location
      real :: woa_nutrient
      INTEGER           :: include_depth_output,include_temp_surface_output,include_temp_bottom_output,&
      include_chlorophyll_surface_output,include_phyto_biomass_surface_output,&
      include_phyto_biomass_bottom_output,include_PAR_surface_output,include_PAR_bottom_output,&
      include_windspeed_output,include_stressx_output,include_stressy_output,include_Etide_output,&
      include_Ewind_output,include_u_mean_surface_output,include_u_mean_bottom_output,&
      include_grow1_mean_surface_output,include_grow1_mean_bottom_output,&
      include_uptake1_mean_surface_output,include_uptake1_mean_bottom_output,include_tpn1_output,&
      include_tpg1_output,include_speed3_output

!
!   Initialise Winteracter
!
!      print*, 'getting physics defaults'
      call get_physics_defaults(lat_in_domain,lon_in_domain,run_year,start_year,unique_job_id,&
      met_data_location,iline,woa_nutrient,include_depth_output,include_temp_surface_output,&
      include_temp_bottom_output,include_chlorophyll_surface_output,include_phyto_biomass_surface_output,&
      include_phyto_biomass_bottom_output,include_PAR_surface_output,include_PAR_bottom_output,&
      include_windspeed_output,include_stressx_output,include_stressy_output,include_Etide_output,&
      include_Ewind_output,include_u_mean_surface_output,include_u_mean_bottom_output,&
      include_grow1_mean_surface_output,include_grow1_mean_bottom_output,include_uptake1_mean_surface_output,&
      include_uptake1_mean_bottom_output,include_tpn1_output,include_tpg1_output,include_speed3_output)
      depth=newdepth; N=newN; dz=newdz
!      print*, 'getting phyto defaults'
      call get_phyto_defaults()
!      print*, 'getting grazing defaults'
      call get_grazing_defaults()
      allocate (contour_data(1,1),xdata(1),ydata(1,1),profile(1,1))
      allocate (au(1),bu(1),cu(1),du(1),ru(1),qu(1),tke(1),tkeold(1),eps(1))
!      ifront=1; met_type=1; metfile_error=0
      ifront=1; met_type=0; metfile_error=0
!      print*, 'met_type = (1 for default)'
!      read(5,*) met_type
!      print*, 'getting met'
      call get_met(lat_in_domain,lon_in_domain,run_year,start_year,unique_job_id,met_data_location,iline)
      surface_out='surface'//unique_job_id//'.dat'; phys_profile_day='physday.dat'; bio_profile_day='biolday.dat';&
       monthly_out='monthly'//unique_job_id//'.dat'
      phys_profile_hour='physhour.dat'; bio_profile_hour='biolhour.dat'
      metfile_in='Model default meteorological data'
      user_metfile=''
      initialdata='Initialisation data not saved'
      output_start=0; output_end=0; output_hr_start=0; output_hr_end=0; JD_front=182
      temp_min=5.0; temp_max=20.0; zmin=0.0; zmax1=10.0; zmax2=10.0
      irunplot=-1; irun=0

!!! new calls to ProcessMenu (rma, 21/10/11)
!      print*, 'IDENT = '
      IDENT = 1
!      read(5,*) IDENT
      CALL ProcessMenu(IDENT,QUIT,unique_job_id)
!!! end new calls to ProcessMenu (rma, 21/10/11)

      STOP
      END PROGRAM phytowin
!
!*****************************************************************************
!
      SUBROUTINE ProcessMenu(IDENT,QUIT,unique_job_id)
!
!   This subroutine processes the menu selections
!
      USE INTERFACES
      use turbulence
      use biology
      use physics
      use graphics_variables
      use variables_all
      use variables
      use file_names
!
      IMPLICIT NONE
!
      integer nhu3a,nhu3b
      INTEGER, INTENT (IN)     :: IDENT
      LOGICAL, INTENT (IN OUT) :: QUIT
      ! character(len=4)  :: run_year,start_year,iline
      INTEGER  :: run_year,start_year,iline
      character(len=12) :: lat_in_domain
      character(len=12) :: lon_in_domain
      character(len=36) :: unique_job_id
      character(len=300) :: met_data_location
      real :: woa_nutrient
      INTEGER           :: include_depth_output,include_temp_surface_output,include_temp_bottom_output,&
      include_chlorophyll_surface_output,include_phyto_biomass_surface_output,include_phyto_biomass_bottom_output,&
      include_PAR_surface_output,include_PAR_bottom_output,include_windspeed_output,include_stressx_output,&
      include_stressy_output,include_Etide_output,include_Ewind_output,include_u_mean_surface_output,&
      include_u_mean_bottom_output,include_grow1_mean_surface_output,include_grow1_mean_bottom_output,&
      include_uptake1_mean_surface_output,include_uptake1_mean_bottom_output,include_tpn1_output,&
      include_tpg1_output,include_speed3_output


!
!   Branch depending on chosen menu item
!
              call save_work(lat_in_domain,lon_in_domain,run_year,start_year,unique_job_id,met_data_location,&
              iline,woa_nutrient,include_depth_output,include_temp_surface_output,include_temp_bottom_output,&
              include_chlorophyll_surface_output,include_phyto_biomass_surface_output,&
              include_phyto_biomass_bottom_output,include_PAR_surface_output,include_PAR_bottom_output,&
              include_windspeed_output,include_stressx_output,include_stressy_output,include_Etide_output,&
              include_Ewind_output,include_u_mean_surface_output,include_u_mean_bottom_output,&
              include_grow1_mean_surface_output,include_grow1_mean_bottom_output,include_uptake1_mean_surface_output,&
              include_uptake1_mean_bottom_output,include_tpn1_output,include_tpg1_output,include_speed3_output)
              deallocate (contour_data,xdata,ydata,profile,au,bu,cu,du,ru,qu,tke,tkeold,eps)
!
!       Open required files for results, and write headers
!
              open(21,file=surface_out,status='replace')
              if(output_end.gt.0.0)then
                open(22,file=bio_profile_day,status='replace')
                open(23,file=phys_profile_day,status='replace')
              end if
              if(output_hr_end.gt.0.0)then
                open(24,file=bio_profile_hour,status='replace')
                open(25,file=phys_profile_hour,status='replace')
              end if
              open(26,file=monthly_out,status='replace')

              if(ifront.eq.1)then
                open(27,file='front_data'//unique_job_id//'.dat',status='replace')
                write(27,fmt="('    hu3   height    temp    X1      N1      U1      G1      DIN      hKz      log10Kz')")
                ifront=0
              else
                i=1
                open(27,file='front_data'//unique_job_id//'.dat')
                read(27,*,end=322)
320			  read(27,*,end=322) hu3a(i),h_hu3(i),t_hu3(i),x1_hu3(i),n1_hu3(i),u1_hu3(i),g1_hu3(i), &
			    s_hu3(i),hkz_hu3(i),Kz_hu3(i)
			    i=i+1
			    goto 320
322		    close(27)
		      nhu3a=i-1
			    open(27,file='front_data'//unique_job_id//'.dat',status='replace')
			    write(27,fmt="('	hu3   height	temp	X1	N1	U1	G1	DIN	 hKz	  log10Kz')")
			    if(nhu3a.gt.1)then
			      do i=1,nhu3a
				write(27,fmt="(10f8.2)") hu3a(i),h_hu3(i),t_hu3(i),x1_hu3(i),n1_hu3(i),u1_hu3(i),g1_hu3(i),s_hu3(i),hkz_hu3(i),Kz_hu3(i)
			      end do
			    end if
			  end if
!
!       Write headers to output files
              WRITE(21,421)
              if(output_end.gt.0)then
                WRITE(22,422)
                write(23,423)
              end if
              if(output_hr_end.gt.0)then
                WRITE(24,424)
                WRITE(25,425)
              end if
              write(26,426)

421   FORMAT(1x,' JD     time       SH        Ts        Tb       Ts-Tb      PHI     &
                       & spd     stress     Qs       Qflux      CHLs      CHLt     &
                       & DINs    netp      grossp')
422   FORMAT(1x,'   time       hu3     height     PAR       chl1      n1        grow1    &
                       & uptake1   chl2      n2        grow2     uptake2   DIN')
423   FORMAT(1x,'    time      hu3      height     temp     sigmat       u         v     &
                       &  h_turb     Ri        logKz    Logdiss   Logtke')
424   FORMAT(1x,'   time       hu3     height     PAR       chl1      n1        grow1    &
                       &  uptake1   chl2      n2        grow2     uptake2   DIN')
425   FORMAT(1x,'    time      hu3      height     temp     sigmat     u         v       &
                       &  h_turb     N2         S2         Ri       logKz     Logdiss    Logtke')
426   FORMAT(1x,' month   SH        Ts        Tb      delta-T   Wstress    Chlt       Ct       Cgross    Cnet     accumC')
!

              CALL run_model(run_year,start_year,unique_job_id,iline,woa_nutrient,include_depth_output,&
			    include_temp_surface_output,include_temp_bottom_output,include_chlorophyll_surface_output,&
			    include_phyto_biomass_surface_output,include_phyto_biomass_bottom_output,include_PAR_surface_output,&
			    include_PAR_bottom_output,include_windspeed_output,include_stressx_output,include_stressy_output,&
			    include_Etide_output,include_Ewind_output,include_u_mean_surface_output,include_u_mean_bottom_output,&
          include_grow1_mean_surface_output,include_grow1_mean_bottom_output,include_uptake1_mean_surface_output,&
          include_uptake1_mean_bottom_output,include_tpn1_output,include_tpg1_output,include_speed3_output)

              surface_out='surface'//unique_job_id//'.dat'; phys_profile_day='physday.dat'; bio_profile_day='biolday.dat';&
               monthly_out='monthly'//unique_job_id//'.dat'

              close(21)
              if(output_end.gt.0.0)close(22)
              if(output_end.gt.0.0)close(23)
              if(output_hr_end.gt.0.0)close(24)
              if(output_hr_end.gt.0.0)close(25)
              close(26)
              close(27)
              close(28)
              close(50)
              close(55)
              close(56)
              close(57)
              close(54)

      RETURN
      END SUBROUTINE ProcessMenu
!******************************************************************************************
!!!!      SUBROUTINE SpawnChild(IHANDLE, in_flags,in_x,in_y,in_width,in_height,in_title)
!
!*****************************************************************************
!
      Subroutine save_work(lat_in_domain,lon_in_domain,run_year,start_year,unique_job_id,met_data_location,iline,woa_nutrient,&
			    include_depth_output,include_temp_surface_output,include_temp_bottom_output,include_chlorophyll_surface_output,&
          include_phyto_biomass_surface_output,include_phyto_biomass_bottom_output,include_PAR_surface_output,&
          include_PAR_bottom_output,include_windspeed_output,include_stressx_output,include_stressy_output,&
          include_Etide_output,include_Ewind_output,include_u_mean_surface_output,include_u_mean_bottom_output,&
          include_grow1_mean_surface_output,include_grow1_mean_bottom_output,include_uptake1_mean_surface_output,&
          include_uptake1_mean_bottom_output,include_tpn1_output,include_tpg1_output,include_speed3_output)
      USE VARIABLES
      use physics
      use turbulence
      use variables_all
      use biology
      use file_names

      implicit none

      integer :: iflags,testdate(8),idate(8),iclock(6),int_date,int_clock
      integer :: run_year,start_year,iline
      real :: real_clock
      real :: woa_nutrient_tmp, woa_nutrient
      real :: smaj1, smin1, smaj2, smin2, smaj3, smin3, smaj4, smin4, smaj5, smin5, alldepth
      character(len=6) :: met_time
      character(len=30) :: filter
      character (len=50) :: filein
      character(len=8) :: date
      character(len=10) :: clock,zone
      character(len=100) :: initialdata2
      character(len=1) :: ans
      character(len=3) :: type
      character(len=300) :: domain_file
      character(len=300) :: nutrient_file
      character(len=12) :: lat_in_domain
      character(len=12) :: lon_in_domain
      character(len=36) :: unique_job_id
      character(len=300) :: met_data_location
      integer :: include_depth_output,include_temp_surface_output,include_temp_bottom_output,include_chlorophyll_surface_output,&
			    include_phyto_biomass_surface_output,include_phyto_biomass_bottom_output,include_PAR_surface_output,&
          include_PAR_bottom_output,include_windspeed_output,include_stressx_output,include_stressy_output,&
          include_Etide_output,include_Ewind_output,include_u_mean_surface_output,include_u_mean_bottom_output,&
          include_grow1_mean_surface_output,include_grow1_mean_bottom_output,include_uptake1_mean_surface_output,&
          include_uptake1_mean_bottom_output,include_tpn1_output,include_tpg1_output,include_speed3_output



      if(initialdata.eq.'Initialisation data not saved')then
        ! initialdata2='initial.txt'
        initialdata2='initial'//unique_job_id//'.txt'
      else
        initialdata2=initialdata
      end if

      initialdata=initialdata2
      ! initialdata = 'initial'//unique_job_id//'.txt'
      open(20,file=initialdata,status='replace')
      call date_and_time(date,clock,zone,testdate)

      write(20,fmt="('Initialisation and driving parameters for S2P3 (Shelf Sea Physics and Primary Production) model')")
      write(20,fmt="('In: Simpson & Sharples, Introduction to the Physical and Biological Oceanography of Shelf Seas')")
      write(20,fmt="('    Cambridge University Press, 2012.')")
      write(20,*)
      write(20,fmt="('File generated:')")
      write(20,*)
      write(20,*)
      write(20,fmt="('Physics parameters:')")
      write(20,fmt="('Total depth (m)                                        =   ',f9.2)") depth
      write(20,fmt="('Number of depth cells                                  =   ',i9)") N
      write(20,fmt="('Depth cell thickness (m)                               =   ',f9.2)") dz
      write(20,fmt="('Time step (s)                                          =   ',f9.2)") time_step
      write(20,fmt="('Longitude (degrees, positive east)                     =   ',f9.2)") lon
      write(20,fmt="('Latitude (degrees, positive north)                     =   ',f9.2)") lat
      write(20,fmt="('Bottom quadratic drag coefficient                      =   ',f9.5)") kb
      write(20,fmt="('Maximum diffusivity and viscosity (m2 s-1)             =   ',f9.6)") vismax
      write(20,fmt="('Background viscosity (m2 s-1)                          =   ',f9.6)") Nz_bg
      write(20,fmt="('Background diffusivity (m2 s-1)                        =   ',f9.6)") Kz_bg
      write(20,fmt="('Initial water temperature (deg C)                      =   ',f9.2)") first_temp
      write(20,fmt="('Heat vertical attenuation (m-1)                        =   ',f9.3)") lambda
      write(20,fmt="('Chl effect on heat attenuation (m2 (mg Chl)-1)         =   ',f9.5)") heat_shade
      write(20,fmt="('PAR vertical attenuation (m-1)                         =   ',f9.3)") par_atten
      write(20,fmt="('Fraction of surface radiation that is PAR              =   ',f9.2)") par_percent
      write(20,fmt="('Maximum seabed dissolved inorganic N (mmol m-3)        =   ',f9.2)") bed_din
      write(20,fmt="('Maximum flux of inorganic N from seabed (mmol m-2 d-1) =   ',f9.2)") din_rate

      write(20,*)
      write(20,fmt="('Tidal parameters:')")
      write(20,fmt="('                                            M2       S2       N2       O1       K1')")
      write(20,fmt="('u amplitude (m s-1)                  : ',5f9.3)") (uamp(i),i=1,5)
      write(20,fmt="('u phase (radians)                    : ',5f9.3)") (uphase(i),i=1,5)
      write(20,fmt="('v amplitude (m s-1)                  : ',5f9.3)") (vamp(i),i=1,5)
      write(20,fmt="('v phase (radians)                    : ',5f9.3)") (vphase(i),i=1,5)

!	     Calculate ellipse parameters

      f0=(4.0*3.142/(24.0*3600.0))*dsin(lat*3.142/180.0); f1=-f0		! calculate Coriolis parameter.
      do i=1,5
         uac=0.5*sqrt(uamp(i)**2.0+vamp(i)**2.0+(2.0*uamp(i)*vamp(i)*sin(vphase(i)-uphase(i))))
         uc=0.5*sqrt(uamp(i)**2.0+vamp(i)**2.0-(2.0*uamp(i)*vamp(i)*sin(vphase(i)-uphase(i))))
         gac1=(vamp(i)*COS(vphase(i)))-(uamp(i)*SIN(uphase(i)))
         gac2=(uamp(i)*COS(uphase(i)))+(vamp(i)*SIN(vphase(i)))
         gc1=(vamp(i)*COS(vphase(i)))+(uamp(i)*SIN(uphase(i)))
         gc2=(uamp(i)*COS(uphase(i)))-(vamp(i)*SIN(vphase(i)))
         gac=dATAN2(gac1,gac2)
         gc=dATAN2(gc1,gc2)

         if (uamp(i).lt.0.0001.and.vamp(i).lt.0.0001) then
             orient(i)=0.0
             major_phase(i)=0.0
             semi_major(i)=0.0
             semi_minor(i)=0.0
             polarisn(i)=0.0
           else
             orient(i)=(gac+gc)/2.0
             major_phase(i)=1.57-(gac-gc)/2.0
             semi_major(i)=uac+uc
             semi_minor(i)=uac-uc
             polarisn(i)=semi_minor(i)/semi_major(i)
         end if

      END do

!      print*, 'over-ride M2, S2, N2 tides (y/n)?'
!      read(*,'(a1)') ans
!      if(ans.eq.'y') then
        ans = 'y'

!	open(1,file='s12_m2_s2_n2_h.dat',status='old')
!	open(1,file='s12_m2_s2_n2_h_sec.dat',status='old')
!	open(1,file='s12_m2_s2_n2_h_tim.dat',status='old')

! get tidal current amplitudes

! first read whether map, sec or tim
        read(5,'(i4)') start_year
        read(5,'(i4)') run_year
        read(5,'(a12)') lat_in_domain
        read(5,'(a12)') lon_in_domain
        read(5,'(a300)') domain_file
        read(5,'(a300)') nutrient_file
        read(5,'(a36)')  unique_job_id
        read(5,'(a300)') met_data_location
        read(5,'(a3)') type
        read(5,*) iline
        read(5,'(f6.1)') smaj1
        read(5,'(f6.1)') smin1
        read(5,'(f6.1)') smaj2
        read(5,'(f6.1)') smin2
        read(5,'(f6.1)') smaj3
        read(5,'(f6.1)') smin3
        read(5,'(f6.1)') smaj4
        read(5,'(f6.1)') smin4
        read(5,'(f6.1)') smaj5
        read(5,'(f6.1)') smin5
        read(5,'(f6.1)') woa_nutrient
        read(5,'(f6.1)') alldepth
        read(5,'(i1)') include_depth_output
        read(5,'(i1)') include_temp_surface_output
        read(5,'(i1)') include_temp_bottom_output
        read(5,'(i1)') include_chlorophyll_surface_output
        read(5,'(i1)') include_phyto_biomass_surface_output
        read(5,'(i1)') include_phyto_biomass_bottom_output
        read(5,'(i1)') include_PAR_surface_output
        read(5,'(i1)') include_PAR_bottom_output
        read(5,'(i1)') include_windspeed_output
        read(5,'(i1)') include_stressx_output
        read(5,'(i1)') include_stressy_output
        read(5,'(i1)') include_Etide_output
        read(5,'(i1)') include_Ewind_output
        read(5,'(i1)') include_u_mean_surface_output
        read(5,'(i1)') include_u_mean_bottom_output
        read(5,'(i1)') include_grow1_mean_surface_output
        read(5,'(i1)') include_grow1_mean_bottom_output
        read(5,'(i1)') include_uptake1_mean_surface_output
        read(5,'(i1)') include_uptake1_mean_bottom_output
        read(5,'(i1)') include_tpn1_output
        read(5,'(i1)') include_tpg1_output
        read(5,'(i1)') include_speed3_output
        imode=1

! 	do i = 1,iline
! 	read(1,'(16x,10f6.1)') smaj1, smin1, smaj2, smin2, smaj3, smin3, smaj4, smin4, smaj5, smin5
! !Note, reading in 1st 6 constituents of teh tides, but only using the first 5
! !	write(6,'(16x,6f6.1)') smaj1, smin1, smaj2, smin2, smaj3, smin3
! 	if(i.eq.iline) then
	semi_major(1)=smaj1*1e-2
	semi_minor(1)=smin1*1e-2
	semi_major(2)=smaj2*1e-2
	semi_minor(2)=smin2*1e-2
	semi_major(3)=smaj3*1e-2
	semi_minor(3)=smin3*1e-2
  semi_major(4)=smaj4*1e-2
	semi_minor(4)=smin4*1e-2
  semi_major(5)=smaj5*1e-2
	semi_minor(5)=smin5*1e-2
  !       endif
	! enddo

	close (1)

      do i=1,5
         polarisn(i)=semi_minor(i)/semi_major(i)
      enddo

      write(20,*)
      write(20,fmt="('tidal ellipse orientation (radians)  : ',5f9.3)") (orient(i),i=1,5)
      write(20,fmt="('tidal ellipse polarisation           : ',5f9.3)") (polarisn(i),i=1,5)
      write(20,fmt="('tidal ellipse semi-major axis (m s-1): ',5f9.3)") (semi_major(i),i=1,5)

      write(20,*)
      write(20,fmt="('Biology parameters:')")
      write(20,fmt="('Growth model (1=Eppley, 2=Q10)                                    =   ',i4)") growth_model

      if(growth_model.eq.2)then
      write(20,*)
      write(20,fmt="('    Parameters for Q10 growth model:')")
      write(20,fmt="('    Reference maximum growth rate (d-1) =  ',f9.2)") Pmax10
      write(20,fmt="('    Reference temperature (deg C)       =  ',f9.2)") gT0
      write(20,fmt="('    Q10 exponent for growth             =  ',f9.2)") Q10
      write(20,*)
      end if

      write(20,fmt="('Max light utilisation coefficient (mg C (mg Chl)-1 d-1 (W m-2)-1) =   ',f9.4)") alpha
      write(20,fmt="('Reference respiration rate (mg C (mg Chl)-1 d-1)                  =   ',f9.4)") respiration
      write(20,fmt="('Reference temperature for respiration rate (deg C)                =   ',f9.4)") rT0
      write(20,fmt="('Q10 exponent for respiration                                      =   ',f9.4)") rQ10
      write(20,fmt="('Chl:carbon (mg Chl (mg C)-1)                                      =   ',f9.4)") chl_carbon
      write(20,fmt="('Near-bed seed stock of phytoplankton (mg C m-3)                   =   ',f9.4)") seed_biomass
      write(20,fmt="('Pigment absorption cross-section (m2 (mg Chl)-1)                  =   ',f9.4)") chl_abscross

			write(20,fmt="('Maximum nitrate uptake rate (mmol (mg Chl)-1 d-1)		  =   ',f9.4)") uptake_max
			write(20,fmt="('Maximum cell nutrient quota (mmol N (mg Chl)-1)			  =   ',f9.4)") max_quota
			write(20,fmt="('Subsistence cell nutrient quota (mmol N (mg Chl)-1)		  =   ',f9.4)") sub_quota
			write(20,fmt="('Nitrate uptake half-saturation concentration (mmol m-3)		  =   ',f9.4)") half_saturation
			write(20,fmt="('Swimming speed (m d-1)						  =   ',f9.4)") swim_speed
			write(20,fmt="('Sinking speed (m d-1)						  =   ',f9.4)") sink_speed
			write(20,fmt="('Minimum grazing impact (d-1)					  =   ',f9.4)") graze_min
			write(20,fmt="('Amplitude of seasonal grazing impact (d-1)			  =   ',f9.4)") Graze_max
			write(20,fmt="('Year day on which maximum grazing impact is reached		  =   ',i9)") int(graze_daymax)
			write(20,fmt="('Biomass threshold for grazing (mg Chl m-3)			  =   ',f9.4)") graze_thresh
			write(20,fmt="('Proportion of grazed organic N recycled to ambient DIN (0.0-1.0)  =   ',f9.4)") Nrecycle

      write(20,*)
      write(20,fmt="('----End of initialisation data file----')")
      close(20)
      return
      end subroutine save_work

!*****************************************************************************
!
      Subroutine load_work()
      USE VARIABLES
      use physics
      use turbulence
      use variables_all
      use biology
      use file_names

      implicit none

      integer :: iflags,i_algae,inerror,npoints,igraze_daymax
      real :: dummy
      character(len=25) :: filter
      character(len=6) :: met_time
      character (len=8) :: ntxt
      character(len=100) :: initialdata2
      character (len=32) :: txt_dz,txt_dt
      character (len=6) :: txt_rdz,txt_rdt
      character (len=60) :: string1
      character (len=39) :: string2
      character (len=70) :: string3
      character (len=43) :: string4
      character(len=36) :: unique_job_id

      if(initialdata.eq.'Initialisation data not saved')then
        initialdata2='initial'//unique_job_id//'.txt'
      else
        initialdata2=initialdata
      end if

      initialdata=initialdata2
      ! initialdata = 'initial'//unique_job_id//'.txt'
      open(20,file=initialdata,status='old')

			read(20,*)
			read(20,*)
			read(20,*)
			read(20,*)
			read(20,*)
			read(20,*)
			read(20,*)
			read(20,*)
			read(20,*)
			read(20,*)
      read(20,fmt="(a60,f9.2)") string1,newdepth
      read(20,fmt="(a60,i9)") string1,newN
      read(20,fmt="(a60,f9.2)") string1,newdz
      read(20,fmt="(a60,f9.2)") string1,time_step
      read(20,fmt="(a60,f9.2)") string1,lon
      read(20,fmt="(a60,f9.2)") string1,lat
      read(20,fmt="(a60,f9.2)") string1,kb
      read(20,fmt="(a60,f9.6)") string1,vismax
      read(20,fmt="(a60,f9.6)") string1,Nz_bg
      read(20,fmt="(a60,f9.6)") string1,Kz_bg
      read(20,fmt="(a60,f9.2)") string1,first_temp
      read(20,fmt="(a60,f9.2)") string1,lambda
      read(20,fmt="(a60,f9.2)") string1,heat_shade
      read(20,fmt="(a60,f9.2)") string1,par_atten
      read(20,fmt="(a60,f9.2)") string1,par_percent
      read(20,fmt="(a60,f9.2)") string1,bed_din
      read(20,fmt="(a60,f9.2)") string1,din_rate
      read(20,*)
      read(20,*)
      read(20,*)
      read(20,fmt="(a39,5f9.3)") string2,(uamp(i),i=1,5)
      read(20,fmt="(a39,5f9.3)") string2,(uphase(i),i=1,5)
      read(20,fmt="(a39,5f9.3)") string2,(vamp(i),i=1,5)
      read(20,fmt="(a39,5f9.3)") string2,(vphase(i),i=1,5)
      read(20,*)
      read(20,*)
      read(20,*)
      read(20,*)
      read(20,*)
      read(20,*)
      read(20,fmt="(a70,i4)") string3,growth_model

      if(growth_model.eq.2)then
        read(20,*)
        read(20,*)
        read(20,fmt="(a43,f9.2)") string4,Pmax10
        read(20,fmt="(a43,f9.2)") string4,gT0
        read(20,fmt="(a43,f9.2)") string4,Q10
        read(20,*)
      end if

      read(20,fmt="(a70,f9.4)") string3,alpha
      read(20,fmt="(a70,f9.4)") string3,respiration
      read(20,fmt="(a70,f9.4)") string3,rT0
      read(20,fmt="(a70,f9.4)") string3,rQ10
      read(20,fmt="(a70,f9.4)") string3,chl_carbon
      read(20,fmt="(a70,f9.4)") string3,seed_biomass
      read(20,fmt="(a70,f9.4)") string3,chl_abscross
			read(20,fmt="(a70,f9.4)") string3,uptake_max
			read(20,fmt="(a70,f9.4)") string3,max_quota
			read(20,fmt="(a70,f9.4)") string3,sub_quota
			read(20,fmt="(a70,f9.4)") string3,half_saturation
			read(20,fmt="(a70,f9.4)") string3,swim_speed
			read(20,fmt="(a70,f9.4)") string3,sink_speed
			read(20,fmt="(a70,f9.4)") string3,graze_min
			read(20,fmt="(a70,f9.4)") string3,Graze_max
			read(20,fmt="(a70,i9)") string3,igraze_daymax
			graze_daymax=dble(igraze_daymax)
			read(20,fmt="(a70,f9.4)") string3,graze_thresh
			read(20,fmt="(a70,f9.4)") string3,Nrecycle

!	     Calculate ellipse parameters

      f0=(4.0*3.142/(24.0*3600.0))*dsin(lat*3.142/180.0); f1=-f0		! calculate Coriolis parameter.
      do i=1,5
         uac=0.5*sqrt(uamp(i)**2.0+vamp(i)**2.0+(2.0*uamp(i)*vamp(i)*sin(vphase(i)-uphase(i))))
         uc=0.5*sqrt(uamp(i)**2.0+vamp(i)**2.0-(2.0*uamp(i)*vamp(i)*sin(vphase(i)-uphase(i))))
         gac1=(vamp(i)*COS(vphase(i)))-(uamp(i)*SIN(uphase(i)))
         gac2=(uamp(i)*COS(uphase(i)))+(vamp(i)*SIN(vphase(i)))
         gc1=(vamp(i)*COS(vphase(i)))+(uamp(i)*SIN(uphase(i)))
         gc2=(uamp(i)*COS(uphase(i)))-(vamp(i)*SIN(vphase(i)))
         gac=dATAN2(gac1,gac2)
         gc=dATAN2(gc1,gc2)

         if (uamp(i).lt.0.0001.and.vamp(i).lt.0.0001) then
             orient(i)=0.0
             major_phase(i)=0.0
             semi_major(i)=0.0
             semi_minor(i)=0.0
             polarisn(i)=0.0
           else
             orient(i)=(gac+gc)/2.0
             major_phase(i)=1.57-(gac-gc)/2.0
             semi_major(i)=uac+uc
             semi_minor(i)=uac-uc
             polarisn(i)=semi_minor(i)/semi_major(i)
         end if

      END do
!
!             calculate maximum timestep allowed for model stability
!
      newdz=newdepth/dble(newN)
      time_step=dble(INT((newdz**2.0)/(2.0*vismax)))
      time_step=time_step-0.5d0                       ! this line prevents exact equality of the stability condition
      tsteps=3600.0/time_step
      do WHILE(dble(INT(tsteps)).ne.tsteps)
              time_step=time_step-1.0; tsteps=3600.0/time_step
      end do
      IF(time_step.lt.1.0)then
        time_step=0.01
        if(dble(INT((newdz**2.0)/(2.0*vismax))).gt.0.1)time_step=0.1
        if(dble(INT((newdz**2.0)/(2.0*vismax))).gt.0.2)time_step=0.2
        if(dble(INT((newdz**2.0)/(2.0*vismax))).gt.0.25)time_step=0.25
        if(dble(INT((newdz**2.0)/(2.0*vismax))).gt.0.5)time_step=0.5
        if(dble(INT((newdz**2.0)/(2.0*vismax))).gt.0.75)time_step=0.75
      end if
      if(time_step.lt.0.09)then
        print*, &
        'The time step is less than 0.1 second!!'//char(13)// &
        'You need to reduce the maximum diffusivity or'//char(13)// &
        'reduce the number of depth cells','Time step too small!'
	call physics_input()
      Else
        txt_dz='Depth resolution is '//txt_rdz//' metres'
        txt_dt='Time step is '//txt_rdt//' seconds'
        print*, &
        txt_dz//char(13)//txt_dt//char(13)// &
        ' ','Information'
      End if
      close(20)
3001  return
      end subroutine load_work

!********************************************************************************
!
Subroutine physics_input()	! Get physics driving parameters

!
      USE VARIABLES
      use physics
      use turbulence
      USE BIOLOGY
      use file_names
      use variables_all
!
      IMPLICIT NONE
      integer :: itype,idepth
      character (len=32) :: txt_dz,txt_dt
      character (len=6) :: txt_rdz,txt_rdt
      real :: rdepth,rlat,rk,rvismax,rfirst_temp,rlambda,rheat_shade, &
      rpar_atten,rpar_percent,rbed_din,rdin_rate

1001  rdepth=real(newdepth); rlat=real(lat); rk=real(kb)
      rvismax=real(vismax); rfirst_temp=real(first_temp); rlambda=real(lambda)
      rheat_shade=real(heat_shade); rpar_atten=real(par_atten); rpar_percent=real(par_percent); rbed_din=real(bed_din)
      rdin_rate=real(din_rate)
!
!   Load and show the modal dialog
!
              ! call get_physics_defaults()
              call tide_parameters(1)
!
!	     Calculate ellipse parameters

			f0=(4.0*3.142/(24.0*3600.0))*dsin(lat*3.142/180.0); f1=-f0								! calculate Coriolis parameter.
      do i=1,5
         uac=0.5*sqrt(uamp(i)**2.0+vamp(i)**2.0+(2.0*uamp(i)*vamp(i)*sin(vphase(i)-uphase(i))))
         uc=0.5*sqrt(uamp(i)**2.0+vamp(i)**2.0-(2.0*uamp(i)*vamp(i)*sin(vphase(i)-uphase(i))))
         gac1=(vamp(i)*COS(vphase(i)))-(uamp(i)*SIN(uphase(i)))
         gac2=(uamp(i)*COS(uphase(i)))+(vamp(i)*SIN(vphase(i)))
         gc1=(vamp(i)*COS(vphase(i)))+(uamp(i)*SIN(uphase(i)))
         gc2=(uamp(i)*COS(uphase(i)))-(vamp(i)*SIN(vphase(i)))
         gac=dATAN2(gac1,gac2)
         gc=dATAN2(gc1,gc2)

         if (uamp(i).lt.0.0001.and.vamp(i).lt.0.0001) then
             orient(i)=0.0
             major_phase(i)=0.0
             semi_major(i)=0.0
             semi_minor(i)=0.0
             polarisn(i)=0.0
           else
             orient(i)=(gac+gc)/2.0
             major_phase(i)=1.57-(gac-gc)/2.0
             semi_major(i)=uac+uc
             semi_minor(i)=uac-uc
             polarisn(i)=semi_minor(i)/semi_major(i)
         end if

      END do

!             calculate maximum timestep allowed for model stability
!
      newdz=newdepth/dble(newN)
      time_step=dble(INT((newdz**2.0)/(2.0*vismax)))
      time_step=time_step-0.5d0                       ! this line prevents exact equality of the stability condition
      tsteps=3600.0/time_step
      do WHILE(dble(INT(tsteps)).ne.tsteps)
              time_step=time_step-1.0; tsteps=3600.0/time_step
      end do
      IF(time_step.lt.1.0)then
        time_step=0.01
        if(dble(INT((newdz**2.0)/(2.0*vismax))).gt.0.1)time_step=0.1
        if(dble(INT((newdz**2.0)/(2.0*vismax))).gt.0.2)time_step=0.2
        if(dble(INT((newdz**2.0)/(2.0*vismax))).gt.0.25)time_step=0.25
        if(dble(INT((newdz**2.0)/(2.0*vismax))).gt.0.5)time_step=0.5
        if(dble(INT((newdz**2.0)/(2.0*vismax))).gt.0.75)time_step=0.75
      end if
      if(time_step.lt.0.09)then
        print*, &
        'The time step is less than 0.1 second!!'//char(13)// &
        'You need to reduce the maximum diffusivity or'//char(13)// &
        'reduce the number of depth cells','Time step too small!'
	goto 1001
      Else
!        call irealtostring(real(newdz),txt_rdz,'(f6.2)')
!        call irealtostring(real(time_step),txt_rdt,'(f6.2)')
        txt_dz='Depth resolution is '//txt_rdz//' metres'
        txt_dt='Time step is '//txt_rdt//' seconds'
        print*, &
        txt_dz//char(13)//txt_dt//char(13)// &
        ' ','Information'
      End if
      RETURN
      END SUBROUTINE physics_input
!
!*********************************************************************************************
!*********************************************************************************************
!
      subroutine tide_parameters(ntide)

      use variables
      use physics
      use variables_all
      IMPLICIT NONE
      integer :: itype,ntide
!
      RETURN
      END SUBROUTINE tide_parameters
!
!*********************************************************************************************
!**********************************************************************************************
!
subroutine get_physics_defaults(lat_in_domain,lon_in_domain,run_year,start_year,unique_job_id,met_data_location,iline,&
			    woa_nutrient,include_depth_output,include_temp_surface_output,include_temp_bottom_output,&
			    include_chlorophyll_surface_output,include_phyto_biomass_surface_output,include_phyto_biomass_bottom_output,&
			    include_PAR_surface_output,include_PAR_bottom_output,include_windspeed_output,include_stressx_output,&
			    include_stressy_output,include_Etide_output,include_Ewind_output,include_u_mean_surface_output,&
			    include_u_mean_bottom_output,include_grow1_mean_surface_output,include_grow1_mean_bottom_output,&
			    include_uptake1_mean_surface_output,include_uptake1_mean_bottom_output,include_tpn1_output,&
			    include_tpg1_output,include_speed3_output)

use physics
use turbulence
use biology
use variables_all
use variables

implicit none

real alldepth
integer :: run_year,start_year,iline
character(len=3) :: type
character(len=12) :: lat_in_domain
character(len=12) :: lon_in_domain
character(len=36) :: unique_job_id
character(len=300) :: domain_file
character(len=300) :: nutrient_file
character(len=300) :: met_data_location
real :: smaj1, smin1, smaj2, smin2, smaj3, smin3, smaj4, smin4, smaj5, smin5
real :: woa_nutrient_tmp, woa_nutrient
integer :: include_depth_output,include_temp_surface_output,include_temp_bottom_output,include_chlorophyll_surface_output,&
			    include_phyto_biomass_surface_output,include_phyto_biomass_bottom_output,include_PAR_surface_output,&
			    include_PAR_bottom_output,include_windspeed_output,include_stressx_output,include_stressy_output,&
          include_Etide_output,include_Ewind_output,include_u_mean_surface_output,include_u_mean_bottom_output,&
          include_grow1_mean_surface_output,include_grow1_mean_bottom_output,include_uptake1_mean_surface_output,&
          include_uptake1_mean_bottom_output,include_tpn1_output,include_tpg1_output,include_speed3_output

! get lon, lat, depth

! first read whether map, sec or tim
        read(5,'(i4)') start_year
        read(5,'(i4)') run_year
        read(5,'(a12)') lat_in_domain
        read(5,'(a12)') lon_in_domain
        read(5,'(a300)') domain_file
        read(5,'(a300)') nutrient_file
        read(5,'(a36)')  unique_job_id
        read(5,'(a300)') met_data_location
        read(5,'(a3)') type
        read(5,*) iline
        read(5,'(f6.1)') smaj1
        read(5,'(f6.1)') smin1
        read(5,'(f6.1)') smaj2
        read(5,'(f6.1)') smin2
        read(5,'(f6.1)') smaj3
        read(5,'(f6.1)') smin3
        read(5,'(f6.1)') smaj4
        read(5,'(f6.1)') smin4
        read(5,'(f6.1)') smaj5
        read(5,'(f6.1)') smin5
        read(5,'(f6.1)') woa_nutrient
        read(5,'(f6.1)') alldepth
        read(5,'(i1)') include_depth_output
        read(5,'(i1)') include_temp_surface_output
        read(5,'(i1)') include_temp_bottom_output
        read(5,'(i1)') include_chlorophyll_surface_output
        read(5,'(i1)') include_phyto_biomass_surface_output
        read(5,'(i1)') include_phyto_biomass_bottom_output
        read(5,'(i1)') include_PAR_surface_output
        read(5,'(i1)') include_PAR_bottom_output
        read(5,'(i1)') include_windspeed_output
        read(5,'(i1)') include_stressx_output
        read(5,'(i1)') include_stressy_output
        read(5,'(i1)') include_Etide_output
        read(5,'(i1)') include_Ewind_output
        read(5,'(i1)') include_u_mean_surface_output
        read(5,'(i1)') include_u_mean_bottom_output
        read(5,'(i1)') include_grow1_mean_surface_output
        read(5,'(i1)') include_grow1_mean_bottom_output
        read(5,'(i1)') include_uptake1_mean_surface_output
        read(5,'(i1)') include_uptake1_mean_bottom_output
        read(5,'(i1)') include_tpn1_output
        read(5,'(i1)') include_tpg1_output
        read(5,'(i1)') include_speed3_output
        imode=1

! do i = 2,iline+1
! read(1,'(f8.3,f8.3,62x,f6.1)') lon,lat,alldepth
! ! read(1,'(f8.3,f8.3,38x,f6.1)') lon,lat,alldepth
! ! write(6,'(f8.3,f8.3,38x,f6.1)') lon,lat,alldepth
! if(i.eq.iline+1) newdepth = alldepth
! enddo
! close (1)

read( lon_in_domain, * )  lon
read( lat_in_domain, * )  lat

newdepth = alldepth

! number of depth levels; bottom drag coefficient
!newN=20; kb=0.003
!original values above - just changed to 10 to speed up for demonstartion
!purposes

!newN=20; kb=0.0025
!if(alldepth.lt.30) newN = 10
!if(alldepth.lt.10) newN = 5
!if(alldepth.lt.5) newN = 3

!kb=0.0025
kb=0.0025
!All vertical levels set to 2m
if(alldepth.gt.0.0) then
newN=nint(alldepth/2.0)
newdepth=2.0*newN
else
newN=1
newdepth=0.0
endif

!set up equal depth intervals of 2m for sections:
if(imode.eq.2) then
if(alldepth.gt.0.0) then
newN=nint(alldepth/2.0)
newdepth=2.0*newN
else
newN=1
newdepth=0.0
endif
endif

!open(1,file='s12_m2_s2_n2_h.dat',status='old')
!open(1,file='s12_m2_s2_n2_h_sec.dat',status='old')
!open(1,file='s12_m2_s2_n2_h_tim.dat',status='old')
vismax=0.1; Nz_bg=1.0e-5; Kz_bg=1.0e-5; first_temp=25.0; lambda=0.1; heat_shade=0.012; bed_din=woa_nutrient; din_rate=10.0
par_atten=0.1; par_percent=0.45
uamp(1)=0.4; uamp(2:5)=0.0; uphase(1:5)=0.0; vphase(1:5)=0.0
vamp(1:5)=0.0
newdz=newdepth/dble(newN)
time_step=dble(INT((newdz**2.0)/(2.0*vismax)))
time_step=time_step-0.5d0                       ! this line prevents exact equality of the stability condition
tsteps=3600.0/time_step
do WHILE(dble(INT(tsteps)).ne.tsteps)
        time_step=time_step-1.0; tsteps=3600.0/time_step
end do
IF(time_step.lt.1.0)then
  time_step=0.01
  if(dble(INT((newdz**2.0)/(2.0*vismax))).gt.0.1)time_step=0.1
  if(dble(INT((newdz**2.0)/(2.0*vismax))).gt.0.2)time_step=0.2
  if(dble(INT((newdz**2.0)/(2.0*vismax))).gt.0.25)time_step=0.25
  if(dble(INT((newdz**2.0)/(2.0*vismax))).gt.0.5)time_step=0.5
  if(dble(INT((newdz**2.0)/(2.0*vismax))).gt.0.75)time_step=0.75
end if
!
!     Calculate ellipse parameters
!
f0=(4.0*3.142/(24.0*3600.0))*dsin(lat*3.142/180.0); f1=-f0								! calculate Coriolis parameter.
do i=1,5
   uac=0.5*sqrt(uamp(i)**2.0+vamp(i)**2.0+(2.0*uamp(i)*vamp(i)*sin(vphase(i)-uphase(i))))
   uc=0.5*sqrt(uamp(i)**2.0+vamp(i)**2.0-(2.0*uamp(i)*vamp(i)*sin(vphase(i)-uphase(i))))
   gac1=(vamp(i)*COS(vphase(i)))-(uamp(i)*SIN(uphase(i)))
   gac2=(uamp(i)*COS(uphase(i)))+(vamp(i)*SIN(vphase(i)))
   gc1=(vamp(i)*COS(vphase(i)))+(uamp(i)*SIN(uphase(i)))
   gc2=(uamp(i)*COS(uphase(i)))-(vamp(i)*SIN(vphase(i)))
   gac=dATAN2(gac1,gac2)
   gc=dATAN2(gc1,gc2)

   if (uamp(i).lt.0.0001.and.vamp(i).lt.0.0001) then
       orient(i)=0.0
       major_phase(i)=0.0
       semi_major(i)=0.0
       semi_minor(i)=0.0
       polarisn(i)=0.0
     else
       orient(i)=(gac+gc)/2.0
       major_phase(i)=1.57-(gac-gc)/2.0
       semi_major(i)=uac+uc
       semi_minor(i)=uac-uc
       polarisn(i)=semi_minor(i)/semi_major(i)
   end if

END do

return
end subroutine get_physics_defaults
!
!*****************************************************************************
!
Subroutine phyto_input()	! Get phytoplankton driving parameters


      use biology
      USE VARIABLES
      use file_names

      IMPLICIT NONE
      integer :: itype
      real :: r(17)
!
              initialdata='Initialisation data not saved'

              call get_phyto_defaults()
!
      RETURN
      END SUBROUTINE phyto_input
!
!**********************************************************************************************
!
subroutine get_phyto_defaults()
use biology
implicit none
sub_quota=0.2; alpha=4.0; chl_carbon=0.03; max_quota=1.0
uptake_max=2.0; half_saturation=0.3; respiration=3.5
rT0=15.0; rQ10=1.0
swim_speed=0.0; sink_speed=0.0
seed_biomass=0.0; chl_abscross=0.012
growth_model=1; Pmax10=1.2; gT0=15.0; Q10=1.0
return
end subroutine get_phyto_defaults
!
!***********************************************************************************************
!
Subroutine grazing()	! Get phytoplankton grazing parameters

      use biology
      USE VARIABLES
      use file_names
!
      IMPLICIT NONE
      integer :: itype
      real :: r(6)

              initialdata='Initialisation data not saved'

              call get_grazing_defaults()
!
      RETURN
      END SUBROUTINE grazing
!
!************************************************************************************************
!
      subroutine get_grazing_defaults()
      use biology
      implicit none

      graze_min=0.12; graze_max=0.0
      graze_daymax=180.0
      graze_thresh=0.1; Nrecycle=0.5

      return
      end subroutine get_grazing_defaults
!************************************************************************************************
!
      subroutine bio_rates
      use biology
      alpha_persec=alpha/(24.0*3600.0)
      uptake_max_persec=uptake_max/(24.0*3600.0)
      respiration_persec=respiration/(24.0*3600.0)
      swim_speed_persec=swim_speed/(24.0*3600.0)
      sink_speed_persec=sink_speed/(24.0*3600.0)
      Pmax10_persec=Pmax10/(24.0*3600.0)
      graze_min_persec=graze_min/(24.0*3600.0)
      graze_max_persec=graze_max/(24.0*3600.0)

      return
      end subroutine bio_rates
!************************************************************************************************
!
      subroutine get_met_type()

      use variables
      use file_names
      IMPLICIT NONE
      integer :: itype,ioption,iflags
      character :: filter

      RETURN
      END SUBROUTINE get_met_type
!
!*************************************************************************************************
!
      subroutine get_met(lat_in_domain,lon_in_domain,run_year,start_year,unique_job_id,met_data_location,iline)

      use physics
      use variables
      use file_names
      use MyTrigonometricFunctions
      IMPLICIT NONE

      integer :: i,idum
      integer :: run_year,start_year,iline
      double precision :: day_angle,direct
!      character(len=:), allocatable :: fileplace
      character(len=12) :: lat_in_domain
      character(len=12) :: lon_in_domain
      character(len=4) :: run_year_str
      character(len=36) :: unique_job_id
      character(len=300) :: met_data_location

      if(met_type.eq.1)then             ! DEfault Celtic Sea met data
        open(60,file='Celtic_met.dat',status='replace')
        ndays=365
        direct=30.0
        do i=1,ndays
          wind_speed(i)=dsqrt(62.9+26.8*dcos(dble(i)*0.017214191)+7.9*dsin(dble(i)*0.017214191))
          wind_dir(i)=direct
          direct=direct+72.0
          if(direct.gt.360.0)direct=direct-360.0
          humid(i)=81.4-2.3*cos(real(i)*0.017214191)+0.8*sin(real(i)*0.017214191)
          airT(i)=12.5-3.3*cos(real(i)*0.017214191)-2.5*sin(real(i)*0.017214191)
          airP(i)=1016.0 !*(1.014-0.0004*cos(real(i)*0.017214191)+0.0009*sin(real(i)*0.017214191))
          day_angle=360.0*real(i)/real(ndays)
          rad_input(i) = 200.0
                declination(i)=0.39637-22.9133*MyCOS(day_angle)+4.02543*MySIN(day_angle)-0.3872*MyCOS(2.0*day_angle)&
                +0.052*MySIN(2.0*day_angle)
                                cloud(i)=66.5+5.1*cos(real(i)*0.017214191)-4.1*sin(real(i)*0.017214191) !100.0*(0.15*cos(real(i)*0.017214191)+0.65)
                                write(60,fmt="(i4,7f10.2)") i,wind_speed(i),wind_dir(i),cloud(i),airT(i),airP(i),humid(i)&
                                ,rad_input(i)
!          print*, i,wind_speed(i),wind_dir(i),cloud(i),airT(i),airP(i),humid(i)
        end do
        close(60)
      else
!        open(60,file=user_metfile,status='old')
        ! fileplace = "../met/spatial_data/"
        write(run_year_str,3001) run_year
3001    format (I4)
        ! open(60,file=fileplace//"meterological_datalat"//TRIM(ADJUSTL(lat_in_domain))//"lon"//TRIM(ADJUSTL(lon_in_domain))//"_"//run_year_str//".dat",status='old')
        open(60,file=TRIM(ADJUSTL(met_data_location))//"meterological_datalat"//TRIM(ADJUSTL(lat_in_domain))&
        //"lon"//TRIM(ADJUSTL(lon_in_domain))//"_"//run_year_str//".dat",status='old')


        i=1
1301    metfile_error=1
        read(60,*,end=1302,err=1303) idum,wind_speed(i),wind_dir(i),cloud(i),airT(i),airP(i),humid(i),rad_input(i),lw_rad_down(i)
!        print*, idum,wind_speed(i),wind_dir(i),cloud(i),airT(i),airP(i),humid(i)
        metfile_error=0
        i=i+1
        goto 1301
1302    ndays=i-1
        close(60)
!        if(ndays.lt.365)then
        if(ndays.lt.366)then
          metfile_error=1
        else
          metfile_error=0
        end if
        do i=1,ndays
!          day_angle=(360.0*dble(i)/dble(ndays)) * (dble(ndays)/365.0)
!The change in the line above alows for simulations longer than 1 year, but assumes all years are 365 days long.
!                declination(i)=0.39637-22.9133*MyCOS(day_angle)+4.02543*MySIN(day_angle)-0.3872*&
!                MyCOS(2.0*day_angle)+0.052*MySIN(2.0*day_angle)
               end do
      end if

      wind_speed(ndays+1)=wind_speed(ndays)
      wind_dir(ndays+1)=wind_dir(ndays)
      humid(ndays+1)=humid(ndays)
      airP(ndays+1)=airP(ndays)
      airT(ndays+1)=airT(ndays)
      cloud(ndays+1)=cloud(ndays)
      rad_input(ndays+1)=rad_input(ndays)
      lw_rad_down(ndays+1)= lw_rad_down(ndays)
 !     declination(ndays+1)=declination(ndays)

1303  return
      end subroutine get_met
!
!************************************************************************************************
!
      subroutine front_set()

      use graphics_variables
      IMPLICIT NONE
      integer :: itype
!
      RETURN
      END SUBROUTINE front_set
!
!************************************************************************************************
!
      subroutine output()

      use variables
      use file_names
      use graphics_variables

      IMPLICIT NONE
      integer :: itype,iflags,out_error
      character :: filter

              if(output_end+output_start.gt.0.and.output_end-output_start.le.0)then
                print*,'For daily and hourly data end day must be greater than start day ','Information'
                out_error=1
              end if
              if(output_hr_end+output_hr_start.gt.0.and.output_hr_end-output_hr_start.le.0)then
                print*,'For daily and hourly data end day must be greater than start day ','Information'
                out_error=1
              end if
      RETURN
      END SUBROUTINE output
!
!*************************************************************************************************
!      SUBROUTINE About()
!
!*****************************************************************************
!      SUBROUTINE Save_Picture()
!
!******************************************************************************
!      SUBROUTINE Redraw()
!
!*****************************************************************************
!      SUBROUTINE Plotit(IPLT)
!
!********************************************************************
!********************************************************************
!********************************************************************
!         MAIN MODEL STARTS HERE..........
!********************************************************************
!********************************************************************
!********************************************************************

     subroutine run_model(run_year,start_year,unique_job_id,iline,woa_nutrient,include_depth_output,include_temp_surface_output,&
			    include_temp_bottom_output,include_chlorophyll_surface_output,include_phyto_biomass_surface_output,&
          include_phyto_biomass_bottom_output,include_PAR_surface_output,include_PAR_bottom_output,include_windspeed_output,&
          include_stressx_output,include_stressy_output,include_Etide_output,include_Ewind_output,include_u_mean_surface_output,&
          include_u_mean_bottom_output,include_grow1_mean_surface_output,include_grow1_mean_bottom_output,&
          include_uptake1_mean_surface_output,include_uptake1_mean_bottom_output,include_tpn1_output,&
          include_tpg1_output,include_speed3_output)

     use variables
     use variables_all
     use physics
     use biology
     use turbulence
     use file_names
     use graphics_variables
     use MyTrigonometricFunctions

     implicit none

     LOGICAL           :: QUIT  = .FALSE.
     integer :: itype,iday,j,julian_day,imonth,mon_day(12),montest

      INTEGER :: ihandle,nsteps_hour,in_flags,itarget,nsteps_count,itotal,nyear,nhr,nsteps_minu,nsteps_count2
      INTEGER :: id,itime,ndays_total,iday_length,n_cubed,iplot_count,iwind,n_wind,nturb,nturbd,ihu3,grad_i,itide
      integer :: run_year,start_year,iline
      real :: month_gross1,month_net1,month_dt,month_ts,month_tb,month_stress,month_chl
      real :: zsum,xm2,ym2,xs2,ys2,starttime,beninput,total_diss,grow1_mean(200),grow2_mean(200),u_cubed
      real :: uptake1_mean(200),uptake2_mean(200),cd,sx2,w_stress,w_dir,Zscm,Tchl,chl_max,rad_in,grad_max
      real :: wx,wy,Zn2,MaxN2,accumulated
      real :: wind_factor,day_tan,day_time,noon_time,rad_ideal,sin_solar_elev
      real :: woa_nutrient
      double precision :: rad_amp,day_angle,sunrise,sunset,day_seconds,S_plus,S_minus,NN_calc,Pr,angle1,angle2,phi,PEA, &
      r,s_phasex(5),s_phasey(5),bed_factor,alf,alf1,alf2,gusty,wind_mean,u3_mean,Etide,Ewind,Qflux,Kz5(200),loghu3

      character (len=20) :: in_title
      character (len=5) :: idstr
      double precision, DIMENSION(:), ALLOCATABLE :: recltest
      integer :: reclen
      character(len=36) :: unique_job_id
      INTEGER ::include_depth_output,include_temp_surface_output,include_temp_bottom_output,include_chlorophyll_surface_output,&
			    include_phyto_biomass_surface_output,include_phyto_biomass_bottom_output,include_PAR_surface_output,&
          include_PAR_bottom_output,include_windspeed_output,include_stressx_output,include_stressy_output,&
          include_Etide_output,include_Ewind_output,include_u_mean_surface_output,include_u_mean_bottom_output,&
          include_grow1_mean_surface_output,include_grow1_mean_bottom_output,include_uptake1_mean_surface_output,&
          include_uptake1_mean_bottom_output,include_tpn1_output,include_tpg1_output,include_speed3_output


!
!
u_cubed=0.0; n_cubed=0
N=newN; depth=newdepth; dz=newdz

! map output
!      if(imode.eq.1) then

! avoid shallow and deep ocean
!      if(depth.lt.2.0.or.depth.gt.200.0) then
!      iday = 1
!      strat_jul = -999.99
!      timeloop_init: do itime=0,ndays
!      write(6,fmt="(i4,2f9.3,4f8.2)")iday,lon,lat,depth,strat_jul,strat_jul,strat_jul
!      iday=iday+1
!      END do timeloop_init
!      go to 999
!the go to 999 jumps to 999 continue, so skip the main routine for these water depth, just writing out missing data
!      endif

!      endif
!
!       Allocate graphics and turbulence arrays
!
    allocate (contour_data(ndays+5,N),xdata(ndays+5),ydata(ndays+5,3),profile(6,N))
    allocate (au(0:N),bu(0:N),cu(0:N),du(0:N),ru(0:N),qu(0:N),tke(0:N),tkeold(0:N),eps(0:N))
!
!       Set physical arrays to initial values.
!
velx_old(1:200)=0.0; vely_old(1:200)=0.0; temp_old(1:200)=first_temp; velx_new(1:200)=0.1
vely_new(1:200)=0.1; temp_new(1:200)=first_temp; Nz(1:200)=vismax/2.0; Kz(1:200)=vismax/2.0
      tke(1:N) = 3.0e-6
      eps(1:N) = 5.0e-10
      rad_amp=radiation(1); radsum=0.0; wx=1.0; wy=1.0
!
!       ....and the biological arrays...
!
x_old(1:200)=0.1; ni_old(1:200)=0.1; ni_old(1:200)=0.1; s_old(1:200)=bed_din; x_new(1:200)=0.1
ni_new(1:200)=0.1; x_new(1:200)=0.1
x_old(1:200)=x_new(1:200); ni_new(1:200)=0.1; s_new(1:200)=bed_din
tpn1=0.0; tpn2=0.0; tpg1=0.0; tpg2=0.0; total_flux=0.0
hourly_net=0.0; hourly_gross=0.0; accumulated=0.0; daily_net=0.0; daily_gross=0.0; prod_net=0.0; &
prod_gross=0.0; prod_net_daily=0.0; prod_gross_daily=0.0
total_uptake=0.0; surf_uptake=0.0
call bio_rates()   ! This changes the daily bio rates into per sec

!Read variables from restart file if running a new year that is not the 1st year of the simulation
if(run_year.ne.start_year) then
allocate(recltest(6*200))
inquire(iolength=reclen) recltest
deallocate(recltest)
!open(1,file="restart.dat",form="unformatted",status="old",action="read")
open(1,file="restart"//unique_job_id//".dat",access='DIRECT',recl=reclen,form='UNFORMATTED')
!open(1,file="restart.dat",access='DIRECT',recl=1200*8,form='UNFORMATTED')
read(1,rec = iline) temp_old,velx_old,vely_old,ni_old,x_old,s_old
close (1)
endif



!	....and the array of days in the months....
!
if(mod(iyear,4).eq.0)then
  mon_day(1)=31; mon_day(2)=29; mon_day(3)=31; mon_day(4)=30; mon_day(5)=31; mon_day(6)=30; mon_day(7)=31; mon_day(8)=31
  mon_day(9)=30; mon_day(10)=31; mon_day(11)=30; mon_day(12)=31
else
  mon_day(1)=31; mon_day(2)=28; mon_day(3)=31; mon_day(4)=30; mon_day(5)=31; mon_day(6)=30; mon_day(7)=31; mon_day(8)=31
  mon_day(9)=30; mon_day(10)=31; mon_day(11)=30; mon_day(12)=31
end if
month_chl=0.0; month_gross1=0.0; month_net1=0.0; month_dt=0.0; month_ts=0.0
month_tb=0.0; month_stress=0.0
!
f0=(4.0*3.142/(24.0*3600.0))*dsin(lat*3.142/180.0); f1=-f0			! calculate Coriolis parameter.
alf=time_step*f0; alf1=1.0/(1.0+((alf**2.0)/4.0)); alf2=1.0-((alf**2.0)/4.0)	! constants used for semi-implicit Coriolis

parmax=real(int(1368.0*0.76*MySIN(abs(lat))*(1.0-0.06)*par_percent*DEXP(-(dz/2.0)*par_atten)))
velmax=REAL(INT(100.0*(sqrt(sum(semi_major(1:5)**2.0)))))			! max velocity for screen plot
idmet=1			! counter for meteorology arrays

!!!!!!!! open graphics window and draw axes !!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       calculate first density profile from temperature....
!
call EQN_of_STATE(temp_new,density)
!
!       SET UP LENGTHSCALE AND HEIGHT ARRAYS
!
zsum=0.0
do i=1,N
      r_height(i)=real(i)*dz-dz/2.0
      height(i)=real(i)*dz
      lscale(i)=0.41*r_height(i)*dsqrt(1.0-(r_height(i)/depth))
      h(i)=dz	! h(i) is array of depth cell thicknesses.
END do
!
c1=1.0/(dz**2.0); c3=c1; c2=kb/dz                       ! constants used in equation of motion
bed_factor=-0.0016*(dz**2.0)+0.0328*dz+0.9357		! this factor is used to correct i=1 currents to 1m above bed
wind_mean=0.0; n_wind=0					! used to calculated mean cubed wind speed
!
!       set up grazing impact array g0
!
do i=1,40000
!do i=1,365
  g0(i)=graze_min_persec+graze_max_persec*sin((3.1415927/366.0)*(real(i)+(182-graze_daymax)));
!  g0(i)=graze_min_persec+graze_max_persec*sin((3.1415927/365.0)*(real(i)+(182-graze_daymax)));
  if(g0(i).lt.graze_min_persec)g0(i)=graze_min_persec+(graze_min_persec-g0(i));
end do
!g0(367)=g0(1)
!g0(366)=g0(1)
!
!       set all counters and initial met variables
!	nsteps_count counts timesteps within individual hours,
!	nyear keeps track of years, nhr steps through whole hours in each day,
!       nhours keeps track of total number of hours stepped through,
!	itotal is the total no. of timesteps for the run,
!       tu3 continuously sums the depth-mean current speed cubed,
!
nsteps_hour=int(3600.0/time_step)					! nsteps_hour is the number of timesteps in 1 hour.
nsteps_minu=int(10.0*60.0/time_step)					! nsteps_minu is the number of timesteps in 10 minutes.
nsteps_count=0								! used to count time steps in each hour
nsteps_count2=0								! used to count time steps in each 10 minutes
nhr=0									! used to count hours through each day
! itotal=nsteps_hour*24*(ndays+3)						! total number of time steps in simulation, including 3 day spin-up
itotal=nsteps_hour*24*ndays					! total number of time steps in simulation
                                        ! note extra 3 days not included because starts from restart files
day_seconds=0.0								! used to track time in seconds through each day
iplot_count=0								! keeps track of the number of daily profiles plotted to the screen

rad0=radiation(1); windspeed=0.0; stressx=0.0; stressy=0.0	! set meteorlogical data for spin-up

! zero arrays used to hold mean turbulence, C and N flux profiles
SS_mean(0:N)=0.0; DD_mean(0:N)=0.0; KK_mean(0:N)=0.0; RN_mean(0:N)=0.0; nturb=0; nturbd=0
!
!	Calculate the tidal slope amplitudes....
!
do itide=1,5
  slopex(itide)=(omega(itide)+(polarisn(itide)*f0))*semi_major(itide)
  slopey(itide)=(f0+(polarisn(itide)*omega(itide)))*semi_major(itide)
!  print*, itide, slopex(itide), slopey(itide)
end do
!
!
!       ****************START MAIN TIME LOOP HERE****************
!
! time=0.0 is at 0000 hrs January 1st, loop starts on January 1st minus
!	72.0 hours to give the eqn of motion 3 days to spin up, including ramping
! up the tidal forcing.
!now removed because model starts from restart files
!
! starttime=-72.0*3600.0; idmet=1; iday=-3; julian_day=-3; imonth=1; montest=1
starttime=0.0; idmet=1; iday=1; julian_day=0; imonth=1; montest=1

!
        rad_sum=0.0
        acount=0.0
timeloop: do itime=1,itotal                     ! <A NAME="START OF TIME LOOP">
!
!       detect window resize/expose, and renew screen
!
	time=starttime+(dble(itime-1)*time_step)		! time in seconds
!
!       update counters....
!
	nsteps_count=nsteps_count+1
	nsteps_count2=nsteps_count2+1
	if(nsteps_count.eq.nsteps_hour+1)nsteps_count=1			! reset counter at each hour
	if(nsteps_count2.eq.nsteps_minu+1)nsteps_count2=1		! reset counter every 10 minutes
	day_seconds=day_seconds+time_step				! increase seconds through the day
!
!       CALCULATE THE CURRENT PROFILE
!
!       first the x and y sea surface slopes for tides in the equation of motion
!
  do itide=1,5
    Px=-slopex(itide)*dcos(omega(itide)*time-major_phase(itide))
    Py=-slopey(itide)*dsin(omega(itide)*time-major_phase(itide))
    angle1=datan2(Py,Px);																													! Next lines rotate slope forcing to account for tidal ellipse orientation
    angle2=angle1+orient(itide);																									!
    Pr=dsqrt(Px**2+Py**2);																												!
    Px=Pr*dcos(angle2);																														!
    Py=Pr*dsin(angle2);																														!
		call EQN_PRESSURE(Px,velx_old)																								! x-component pressure
		call EQN_PRESSURE(Py,vely_old)																								! y-component pressure
	end do

	call EQN_FRICTION(velx_old,velx_new,vely_old,vely_new,stressx,u_bed,v_bed)		! x-component surface and bed friction
	u_bed=velx_new(1)/bed_factor																									! bed x-velocity at h=1m
	call EQN_FRICTION(vely_old,vely_new,velx_old,velx_new,stressy,v_bed,u_bed)		! y-component surface and bed friction
	v_bed=vely_new(1)/bed_factor																									! bed y-velocity at h=1m
	call EQN_CORIOLIS(velx_new,vely_new,time_step,f0,N,alf,alf1,alf2)							! semi-implicit Coriolis
!
!	Calculate u-cubed
!
!	if(iday.gt.-2)then
	! if(iday.gt.2)then
	  speed_mag(1:N)=dsqrt(velx_new(1:N)**2.0+vely_new(1:N)**2.0)
	  mean_speed=sum(speed_mag(1:N))/real(N)
	  u_cubed=u_cubed+mean_speed**3.
	  n_cubed=n_cubed+1
!	end if
!
! Update radiation
!
!  if(iday.gt.1)idmet=iday
  idmet=iday
  noon_time=dble(idmet-1)+0.5
  day_time=(time/(24.0*3600.0))-noon_time
!  day_angle=(2.0*3.1415927)*day_time
!  sin_solar_elev=dcos(day_angle)*MyCOS(lat)*MyCOS(declination(idmet))+MySIN(lat)*MySIN(declination(idmet))
!now directly reading downwelling radiation in
  rad0 = rad_input(idmet)
  lw_rad_down0 = lw_rad_down(idmet)
!  rad_ideal=1368.0*0.76*sin_solar_elev*(1.0-0.06)										! clear sky irradiance W m-2
!  rad0=rad_ideal*(1.0-0.01*cloud(idmet)*0.4-0.000038*cloud(idmet)**2.0)					! irradiance with cloud
!	if(rad0.lt.0.0)rad0=0.0
!
! Heat the surface and distribute heat through the water column
!
	call SURFACE_HEAT()
!      if(mod(itime,24*nsteps_hour).eq.1) print*, itime, rad0, rad_out
!      print*, itime, rad0, temp_old(N), rad_out
  radsum=radsum+(rad0+rad_out)*time_step					! Cumulative total net heat supply
  temp_old(1:N)=temp_new(1:N)							! Update temperature array
!
!  call CONVECTION(temp_old,temp_new,N,nsteps_count,nhr,int(tid2))	    	! Convective mixing routine
!  temp_old(1:N)=temp_new(1:N)							! Update temperature array
!
!  Then mix the vertical thermal structure
!
!      print*, itime, temp_old(N)
  call scalar_diff(temp_new,temp_old,time_step,N,cnpar)
!      print*, itime, temp_new(N)
!
!       Calculate new density profile....
!
  call EQN_of_STATE(temp_new,density)
!
! TURBULENCE CLOSURE BIT.
!
!	surface and bottom stress boundary conditions...
!
	stresss=dsqrt(stressx**2.0+stressy**2.0)
	sx=kb*1025.0*(dsqrt(u_bed**2.0+v_bed**2.0)*u_bed)
	sy=kb*1025.0*(dsqrt(u_bed**2.0+v_bed**2.0)*v_bed)
	stressb=dsqrt(sx**2.0+sy**2.0)

!        if(mod(itime,30*24*nsteps_hour).eq.1) print*, itime/(30*24*nsteps_hour), temp_new(1), temp_new(N)
!
	call TURBULENCE_ke(iday,nsteps_count)
!
! Next lines add turbulence and flux profiles to produce daily mean profiles
!
	SS_mean(0:N)=SS_mean(0:N)+SS(0:N)
	DD_mean(0:N)=DD_mean(0:N)+eps(0:N)
	KK_mean(0:N)=KK_mean(0:N)+Kz(0:N)
	RN_mean(0:N)=RN_mean(0:N)+(NN(0:N)/SS(0:N))
	do i=1,N-1
    NN_calc=-(9.81/1025.0)*(density(i+1)-density(i))/dz
	  if(NN_calc.lt.1.0e-6)NN_calc=1.0e-6
	  NN_mean(i)=NN_mean(i)+NN_calc
	end do
	nturb=nturb+1						! This is the counter being updated for use with hourly mean profiles
!
! PRIMARY PRODUCTION BITS
!
! Profile of PAR
!
  call PAR_profile()
!
! algal production and update of DIN
!
	call PRIMARY_PROD(julian_day)
	call DIN_update()
!
! then turbulent mixing of the nutrient
!
  call scalar_diff(s_new,s_old,time_step,N,cnpar)
!
! input DIN from sediments to bottom depth cell...
!
  beninput=time_step*(din_rate/(24.0*3600.0))*(1.0-(s_old(1)/bed_din))/dz
	s_new(1)=s_new(1)+beninput
!
! update temperature and biochemical arrays
!
	velx_old(1:N)=velx_new(1:N); vely_old(1:N)=vely_new(1:N); temp_old(1:N)=temp_new(1:N); xtotal(1:N)=0.0
  x_old(1:N)=x_new(1:N); ni_old(1:N)=ni_new(1:N); s_old(1:N)=s_new(1:N)
!
!       DATA OUTPUT and Met data update
!
	if(nsteps_count.eq.1)then					! i.e. every hour
	  tid=real((time/3600.0)/24.0)		! tid is the total time in days
	  julian_day=max(1,int(tid+0.5))
    u3_mean=u_cubed/real(n_cubed)					! keep calculating mean u-cubed
!
! calculation of depth mean current speed
!
	  speed3=SUM(dsqrt(velx_new(1:N)**2.0+vely_new(1:N)**2.0))/DBLE(N)	! depth-mean current speed
    IF(speed1.lt.speed2.and.speed2.ge.speed3)speed=speed2							! maximum speed during this day
    speed1=speed2; speed2=speed3
	  NN_mean(0:N)=NN_mean(0:N)/dble(nturb)															! calculate hourly-averaged turbulence parameter profiles
	  SS_mean(0:N)=SS_mean(0:N)/dble(nturb)															!
	  DD_mean(0:N)=DD_mean(0:N)/dble(nturb)															!
	  KK_mean(0:N)=KK_mean(0:N)/dble(nturb)															!
	  RN_mean(0:N)=RN_mean(0:N)/dble(nturb)															!
!
		hourly_net=0.0; hourly_gross=0.0
	  do i=1,N
      hourly_net=hourly_net+dz*prod_net(i)												! Net total primary production / mg C m-2 h-1
      hourly_gross=hourly_gross+dz*prod_gross(i)									! Gross total primary production / mg C m-2 h-1
	  end do
	  daily_net=daily_net+hourly_net
	  daily_gross=daily_gross+hourly_gross

    u_mean(1:N)=u_mean(1:N)+velx_new(1:N)																							! Update profiles for daily means
    v_mean(1:N)=v_mean(1:N)+vely_new(1:N)																							! Update profiles for daily means
    vel_mean(1:N)=vel_mean(1:N)+dsqrt(vely_new(1:N)**2.0+velx_new(1:N)**2.0)					! Update profiles for daily means
	  diss_mean(1:N)=diss_mean(1:N)+eps(1:N)																						!
	  Ri_mean(1:N)=Ri_mean(1:N)+Ri(1:N)																									!
	  Kz_mean(1:N)=Kz_mean(1:N)+Kz(1:N)																									!
	  tke_mean(1:N)=tke_mean(1:N)+tke(1:N)																							!
	  rad_mean(1:N)=rad_mean(1:N)+radbar(1:N)																						!
	  grow1_mean(1:N)=grow1_mean(1:N)+net_grow(1:N)																			!
	  uptake1_mean(1:N)=uptake1_mean(1:N)+uptake(1:N)																		!
!
! output to screen and data files
!
! Hourly output of profiles between times set in graphics information file
!
	  ! if(tid.GE.real(output_hr_start).and.tid.LE.real(output_hr_end))then
	  ! grad_max=0.0
	  ! grad_i=1
    !    do i=1,N
	  !      write(24,fmt="(10f10.4)") tid,dlog10(depth/u3_mean),dble(i)*dz-dz/2.0,temp_new(i),radbar(i),x_new(i),ni_new(i), &
	  !      3600.0*net_grow(i),3600.0*uptake(i),s_new(i)
    !      IF(Ri(i).gt.10.0)Ri(i)=10.0
    !      IF(Ri(i).LT.-1.0)Ri(i)=-1.0
    !      IF(RN_mean(i).gt.10.0)RN_mean(i)=10.0
    !      IF(RN_mean(i).LT.-1.0)RN_mean(i)=-1.0
    !      write(25,fmt="(8f10.3,7f11.6)") tid,dlog10(depth/u3_mean),DBLE(i)*dz-dz/2.0,temp_new(i), &
    !      density(i)-1000.0,velx_new(i),vely_new(i),DBLE(i)*dz,NN_mean(i),SS_mean(i),log10(DD_mean(i)), &
    !      Ri(i),dlog10(Kz(i)),dlog10(eps(i)),log10(tke(i))
	  !    END do
    ! END if
    NN_mean(0:N)=0.0		! Reset arrays for hourly means
    SS_mean(0:N)=0.0		!
    DD_mean(0:N)=0.0		!
    KK_mean(0:N)=0.0		!
    RN_mean(0:N)=0.0		!
    nturb=0
!
!
	  call wind_components(wx,wy)												! Calculate hourly wind components
    windspeed=sqrt(wx**2+wy**2)												!	Calculate wind speed
    wind_mean=wind_mean+(windspeed**3.0)							! Keep track of cubed-wind speed (to allow calculation of wind power)
    n_wind=n_wind+1																		!

    cd=(0.63+0.066*windspeed)/density(N)              ! drag coefficient related to wind speed
	  stressx=1.3*cd*windspeed*wx												! x component of surface wind drag
    stressy=1.3*cd*windspeed*wy												! y component of surface wind drag
!
!	Ouput section for daily data....
!
	  if(nhr.eq.12)then				! At NOON
	    Etide=1000.0*0.003*kb*1025.0*(u3_mean)										! Mixing power in the tidal currents (assumes constant mixing efficiency 0.003)
	    Ewind=1000.0*0.023*cd*0.025*1.3*(wind_mean/real(n_wind))	! Mixing power in the wind (assumes constant mixing efficiency 0.023 and slippage factor=0.025)
      u_mean(1:N)=100.0*u_mean(1:N)/24.0												! daily mean u profile in cm/s !!
      v_mean(1:N)=100.0*v_mean(1:N)/24.0												! daily mean v profile in cm/s !!
      vel_mean(1:N)=100.0*vel_mean(1:N)/24.0										! daily mean v profile in cm/s !!
      diss_mean(1:N)=diss_mean(1:N)/24.0												! daily mean turbulent dissipation profile m2 s-3
      Ri_mean(1:N)=Ri_mean(1:N)/24.0												! daily mean Richardson number
      Kz_mean(1:N)=Kz_mean(1:N)/24.0												! daily mean diffusivity m2 s-1
      tke_mean(1:N)=tke_mean(1:N)/24.0												! daily mean turbulent kinetic energy m2 s-3
      rad_mean(1:N)=rad_mean(1:N)/24.0												! daily mean profile of PAR W m-2
      grow1_mean(1:N)=3600.0*grow1_mean(1:N)									! daily mean growth rate d-1
      uptake1_mean(1:N)=3600.0*uptake1_mean(1:N)								! daily mean DIN uptake rate mmol DIN (mg C)-1 d-1
!	Calculate monthly means
	    if(tid.ge.0.0)then
	      if(montest.lt.mon_day(imonth))then
		month_chl=month_chl+sum(x_new(1:N))*dz
		month_gross1=month_gross1+0.001*daily_gross
		month_net1=month_net1+0.001*daily_net
		month_ts=month_ts+temp_new(N)
		month_tb=month_tb+temp_new(1)
		month_dt=month_dt+(temp_new(N)-temp_new(1))
		month_stress=month_stress+sqrt(stressx**2.0+stressy**2.0)
		montest=montest+1
	      else
		month_chl=(month_chl+(sum(x_new(1:N))*dz))/real(montest)
		month_gross1=(month_gross1+0.001*daily_gross)
		month_net1=(month_net1+0.001*daily_net)
		month_ts=(month_ts+temp_new(N))/real(montest)
		month_tb=(month_tb+temp_new(1))/real(montest)
		month_dt=(month_dt+(temp_new(N)-temp_new(1)))/real(montest)
		month_stress=(month_stress+sqrt(stressx**2.0+stressy**2.0))/real(montest)
		write(26,fmt="(i6,f8.3,3f10.2,f10.4,5f10.2)") &
                      & imonth,dlog10(depth/u3_mean),month_ts,month_tb, &
                      & month_dt,month_stress,month_chl,0.001*month_chl/chl_carbon, &
                      & month_gross1,month_net1,month_net1+accumulated
!                if(imonth.eq.7) print*, 'Jul Ts-Tb = ',month_dt
                if(imonth.eq.7) strat_jul = month_dt
		accumulated=accumulated+month_net1
		month_gross1=0.0; month_net1=0.0; month_dt=0.0
		month_ts=0.0; month_tb=0.0; month_stress=0.0; month_chl=0.0
		montest=1
		imonth=imonth+1
	      end if
	    end if

      id=INT(tid)
!
!     Plotting the daily time series of surface conditions
!
      IF(tid.ge.0.0)then																						! i.e. output switches on out of spin-up time
!
!     Plotting blank profiles over last output date
!
!       store latest profiles so that it can be plotted this day and blanked next day
!
	profile(1,1:N)=REAL(temp_new(1:N))				! temperature
	profile(2,1:N)=REAL(x_new(1:N))					! chlorophyll
	profile(3,1:N)=REAL(s_new(1:N))					! dissolved inorganic nitrogen
        profile(4,1:N)=REAL(dlog10(Kz_mean(1:N)))     ! log Kz
        profile(5,1:N)=REAL(rad_mean(1:N))            ! PAR
        profile(6,1:N)=REAL(vel_mean(1:N))            ! speed
!
        do dataplot=4,9
		ig(dataplot)=1
		if(dataplot.ne.7)then
	    do i=1,N
!              call PLOT(profile(dataplot-3,i),r_height(i),last_x(dataplot),last_y(dataplot),dataplot,ncolour(dataplot))
            end do
          end if
	  if(dataplot.eq.7)then
	    do i=1,N-1
!              call PLOT(profile(dataplot-3,i),height(i),last_x(dataplot),last_y(dataplot),dataplot,ncolour(dataplot))
            end do
          end if
        END do
!
!       Put time series data into arrays used for screen renew events
!
        irenew=irenew+1
        xdata(irenew)=REAL(tid)
        ydata(irenew,1)=REAL(temp_new(N))
        ydata(irenew,2)=REAL(temp_new(1))
        ydata(irenew,3)=REAL(x_new(N))
        contour_data(irenew,1:N)=x_new(1:N)
!
	tpn1=tpn1+0.001*daily_net																		! Keep track of total water column net and gross production
	tpg1=tpg1+0.001*daily_gross																	!
!
	      w_stress=dsqrt(stressx**2.0+stressy**2.0)										! Calculate parameters for output to file
	      Tchl=sum(x_new(1:N)*dz)																			! column integrated biomass
	Qflux=radsum/(24.0*3600.0)																	! net heat flux (+ is into sea surface)
	      loghu3=dlog10(depth/u3_mean)																! log10 (h/u3)
	      phi=PEA(density,N,dz)

!
!				Output surface data file
!
	      write(21,fmt="(i4,f10.1,f10.2,3f10.3,5f10.2,3f10.2,2f10.3)")iday+1, &
                    & tid,loghu3,temp_new(n),temp_new(1),temp_new(n)-temp_new(1), &
                    & phi,speed,w_stress,rad0,Qflux,x_new(n),Tchl,s_new(n),0.001*daily_net, &
                    & 0.001*daily_gross

! section output
      if(imode.eq.2) then

!		      OUTPUT SECTION DATA
!                  *** CENTRED ON DAY 190 ***

              if(iday.eq.190)then
                do i=1,N
                  write(6,fmt="(12f8.3)") lon,lat,dlog10(depth/u3_mean),height(i)-dz/2.0, &
                  temp_new(i),x_new(i),ni_new(i),uptake1_mean(i),grow1_mean(i),s_new(i),height(i),dlog10(Kz_mean(i))
                end do
              end if

              end if

! time series output
if(imode.eq.3) then

!                       OUTPUT TEMPERATURE DATA AT A POINT

x_new_max = 0.0
temp_x_new_max = 0.0
do i=1,N
if(x_new(i).gt.x_new_max) then
x_new_max = x_new(i)
temp_x_new_max = temp_new(i)
endif
enddo

if(iday.ge.1.and.iday.le.365)then
write(6,fmt="(i4,7f8.3)") iday,temp_new(N),temp_new(1),x_new(N),s_new(N),x_new_max,temp_x_new_max,s_new(1)
end if

end if


! map output
!Note, I moved this to get output for each day
if(imode.eq.1) then
! include_depth_output=1
! include_temp_surface_output=1
!write(6,'(2f8.3,5f8.2)') lon,lat,depth,dlog10(depth/u3_mean),strat_jul,rad_sum/acount,month_net1+accumulated,bottom_phyto_biomass
! write(6,fmt="(i4,2f8.3,5f8.2)")iday,lon,lat,depth,temp_new(N),temp_new(1),x_new(N),x_new(1)/chl_carbon

write(6,fmt="(i4,2f8.3)",advance="no")iday,lon,lat
! to include additional variables add more if the template if statements copied bleow before write(6,fmt="()") and add to list where others defined etc
! then add to run_map_parallel by extending the lists containing existing output variables
! if(include_XXNEWXX_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")XXNEWXX
! end if

if(include_depth_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")depth
end if

if(include_temp_surface_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")temp_new(N)
end if

if(include_temp_bottom_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")temp_new(1)
end if

if(include_chlorophyll_surface_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")x_new(N)
end if

if(include_phyto_biomass_surface_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")x_new(N)/chl_carbon
end if

if(include_phyto_biomass_bottom_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")x_new(1)/chl_carbon
end if

if(include_PAR_surface_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")rad_mean(N)
end if

if(include_PAR_bottom_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")rad_mean(1)
end if

if(include_windspeed_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")windspeed
end if

if(include_stressx_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")stressx
end if

if(include_stressy_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")stressy
end if

if(include_Etide_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")Etide
end if

if(include_Ewind_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")Ewind
end if

if(include_u_mean_surface_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")u_mean(N)
end if

if(include_u_mean_bottom_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")u_mean(1)
end if

if(include_grow1_mean_surface_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")grow1_mean(N)
end if

if(include_grow1_mean_bottom_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")grow1_mean(1)
end if

if(include_uptake1_mean_surface_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")uptake1_mean(N)
end if

if(include_uptake1_mean_bottom_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")uptake1_mean(1)
end if

if(include_tpn1_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")tpn1
end if

if(include_tpg1_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")tpg1
end if

if(include_speed3_output.eq.1) then
  write(6,fmt="(1f8.2)",advance="no")speed3
end if

write(6,fmt="()")

! ,5f8.2)")iday,lon,lat,depth,temp_new(N),temp_new(1),x_new(N),x_new(1)/chl_carbon
! write (*,"(3f8.3)",advance="no") a,b,c

!note 2f8.3,8f8.2 means that there are two columns with 8 spaces for digits,
!three numbers below the decimal place, then 8 columns with 8 spaces for
!and two numbers behind the decimal place
endif

!
!         profiles between times set in graphics information file

        do i=1,N
		      if(tid.GE.real(output_start).and.tid.LE.real(output_end))then
				    write(22,fmt="(9f10.4)") tid,dlog10(depth/u3_mean), &
                                          & dble(i)*dz-dz/2.0,rad_mean(i),x_new(i),ni_new(i), &
                                          & grow1_mean(i),uptake1_mean(i),s_new(i)
            IF(Ri_mean(i).gt.10.0)Ri_mean(i)=10.0
            IF(Ri_mean(i).LT.-1.0)Ri_mean(i)=-1.0
            write(23,fmt="(12f10.3)") tid,dlog10(depth/u3_mean),DBLE(i)*dz-dz/2, &
                 & temp_new(i),density(i)-1000.0,0.01*u_mean(i),0.01*v_mean(i),DBLE(i)*dz, &
              Ri_mean(i),dlog10(Kz_mean(i)),dlog10(diss_mean(i)),dlog10(tke_mean(i))
          end if
	      END do

	    end if		! end of tid.ge.0.0 if-then

!         Reset daily mean turbulence and biology profiles

      diss_mean(1:N)=0.0; Kz_mean(1:N)=0.0; tke_mean(1:N)=0.0; Ri_mean(1:N)=0.0; rad_mean(1:N)=0.0
      grow1_mean(1:N)=0.0; uptake1_mean(1:N)=0.0; vel_mean(1:N)=0.0; u_mean(1:N)=0.0; v_mean(1:N)=0.0
      total_diss=0.0; radsum=0.0; daily_net=0.0; daily_gross=0.0

	  end if		  ! end of NOON selection

	  prod_gross=0.0; prod_net=0.0;
    wind_mean=0.0; n_wind=0
    nhr=nhr+1             ! update hours
	  if(nhr.eq.24)then
	    nhr=0
	    iday=iday+1
	  end if

	END if						!!!!!!! End of hourly selection

END do timeloop             ! END OF MAIN TIME LOOP

month_chl=(month_chl+(sum(x_new(1:N)))*dz)/real(montest)
month_gross1=(month_gross1+0.001*daily_gross)
month_net1=(month_net1+0.001*daily_net)
month_ts=(month_ts+temp_new(N))/real(montest)
month_tb=(month_tb+temp_new(1))/real(montest)
month_dt=(month_dt+(temp_new(N)-temp_new(1)))/real(montest)
month_stress=(month_stress+sqrt(stressx**2.0+stressy**2.0))/real(montest)

        write(26,fmt="(i6,f8.3,3f10.2,f10.4,5f10.2)") imonth, &
             & dlog10(depth/u3_mean),month_ts,month_tb,month_dt, &
             & month_stress,month_chl,0.001*month_chl/chl_carbon, &
             & month_gross1,month_net1,month_net1+accumulated
!print*, '   lon     lat   depth log(h/u3) strat  accumC '


! ! map output
! if(imode.eq.1) then
! !write(6,'(2f8.3,5f8.2)') lon,lat,depth,dlog10(depth/u3_mean),strat_jul,rad_sum/acount,month_net1+accumulated,bottom_phyto_biomass
! ! write(6,fmt="(i4,2f8.3,5f8.2)")iday,lon,lat,depth,temp_new(N),temp_new(1),x_new(N),x_new(1)/chl_carbon
! !note 2f8.3,8f8.2 means that there are two columns with 8 spaces for digits,
! !three numbers below the decimal place, then 8 columns with 8 spaces for
! !and two numbers behind the decimal place
!
! write(6,fmt="(i4,2f8.3)",advance="no")iday,lon,lat
!
!
! if(include_depth_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")depth
! end if
!
! if(include_temp_surface_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")temp_new(N)
! end if
!
! if(include_temp_bottom_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")temp_new(1)
! end if
!
! if(include_chlorophyll_surface_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")x_new(N)
! end if
!
! if(include_phyto_biomass_surface_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")x_new(N)/chl_carbon
! end if
!
! if(include_phyto_biomass_bottom_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")x_new(1)/chl_carbon
! end if
!
! if(include_PAR_surface_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")rad_mean(N)
! end if
!
! if(include_PAR_bottom_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")rad_mean(1)
! end if
!
! if(include_windspeed_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")windspeed
! end if
!
! if(include_stressx_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")stressx
! end if
!
! if(include_stressy_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")stressy
! end if
!
! if(include_Etide_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")Etide
! end if
!
! if(include_Ewind_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")Ewind
! end if
!
! if(include_u_mean_surface_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")u_mean(N)
! end if
!
! if(include_u_mean_bottom_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")u_mean(1)
! end if
!
! if(include_grow1_mean_surface_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")grow1_mean(N)
! end if
!
! if(include_grow1_mean_bottom_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")grow1_mean(1)
! end if
!
! if(include_uptake1_mean_surface_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")uptake1_mean(N)
! end if
!
! if(include_uptake1_mean_bottom_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")uptake1_mean(1)
! end if
!
! if(include_tpn1_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")tpn1
! end if
!
! if(include_tpg1_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")tpg1
! end if
!
! if(include_speed3_output.eq.1) then
!   write(6,fmt="(1f8.2)",advance="no")speed3
! end if
!
! write(6,fmt="()")
!
! endif

! time series output
if(imode.eq.3) then

x_new_max = 0.0
temp_x_new_max = 0.0
do i=1,N
if(x_new(i).gt.x_new_max) then
x_new_max = x_new(i)
temp_x_new_max = temp_new(i)
endif
enddo

write(6,fmt="(i4,7f8.3)") iday,temp_new(N),temp_new(1),x_new(N),s_new(N),x_new_max,temp_x_new_max,s_new(1)

endif

!
110 ndays_total=INT((time-starttime)/(3600.0*24.0))
!close(60)

tpn1=tpn1*real(ndays)/REAL(ndays_total); tpg1=tpg1*real(ndays)/REAL(ndays_total)
call results_report()
!call WBitmapPut(1,0,0)
!

!Write out a restart file to get read in when teh next year runs

!6 is the number of variables I think need to be written out and read back in for restart
!This approch to specifying recl seems overly complicated but should ensure compatability
!with ifortran and gfortran
allocate(recltest(6*200))
inquire(iolength=reclen) recltest
deallocate(recltest)
! open(1,file="restart.dat",access='DIRECT',recl=1200*8,form='UNFORMATTED')
open(1,file="restart"//unique_job_id//".dat",access='DIRECT',recl=reclen,form='UNFORMATTED')
write(1,rec = iline) temp_old,velx_old,vely_old,ni_old,x_old,s_old
close (1)




999 continue
return
end subroutine run_model                    !!!!!!!! END OF MAIN MODEL !!!!!!!!
!**********************************************************

!**********************************************************
!
!       FUNCTIONS AND SUBROUTINES.
!
!**********************************************************

!**********************************************************
!				Potential Energy Anomaly
!**********************************************************
!
FUNCTION PEA(density,N,dz)
implicit none

double precision :: density(200),pe1,pe2,dz,mean_rho,PEA
integer :: N,j

mean_rho=sum(density(1:N))/dble(N)
pe1=mean_rho*9.81*dz*dble(N)*dz*(dble(N)/2.0)
pe2=0.0
do j=1,N
  pe2=pe2+density(j)*9.81*(dble(j)-0.5)*dz*dz
end do
PEA=(pe1-pe2)/(dble(N)*dz)

return
end
!
!**********************************************************
!       Equation of state of seawater (quadratic in temperature,
!	      constant salinity of 35.00):
!**********************************************************
!
SUBROUTINE EQN_of_STATE(temp,density)
!
use variables_all
implicit none
DOUBLE PRECISION temp(200),density(200),t,rho
INTEGER it,irho
!
!
!	all scalars are rounded to 0.001 to avoid rounding errors
!	causing apparent instabilities in the turbulence routine
!
do i=1,N
	it=nint(temp(i)*1000.0)
	T=dble(it)/1000.0
	rho=(-5.29468d-3*T**2.0)-(6.24956d-2*T)+1028.11
	irho=nint(rho*1000.0)
	density(i)=dble(irho)/1000.0
END do
RETURN
end
!
!**********************************************************
!				Rotate current vector by angle
!**********************************************************
!
SUBROUTINE current_rotate(u,v,angle,N)

implicit none

double precision :: u(200),v(200),angle,angle1(200),angle2(200),u2(200),v2(200),r
integer :: N,i

do i=1,N
    angle1(i)=datan2(v(i),u(i));
    angle2(i)=angle1(i)+angle;
    r=dsqrt(u(i)**2+v(i)**2);
    u2(i)=r*dcos(angle2(i));
    v2(i)=r*dsin(angle2(i));
end do
v=v2
u=u2

return
end
!
!**********************************************************
!       Equation of motion subroutine - pressure term
!**********************************************************
!
SUBROUTINE EQN_PRESSURE(Pgrad,u_old)
!
use variables_all
implicit none

double precision :: Pgrad,u_old(200)

do i=1,N
  u_old(i)=u_old(i)+time_step*Pgrad
!  print*, i, time_step, Pgrad, u_old(i)
end do

return
end
!
!**********************************************************
!       Equation of motion subroutine - friction term
!**********************************************************
!
SUBROUTINE EQN_FRICTION(u_old,u_new,v_old,v_new,stress,u_bed,v_bed)
!
use variables_all
use turbulence
implicit none
!
double precision :: u_old(200),u_new(200),v_old(200),v_new(200),stress,dv,u_bed,v_bed
!
!       Intermediate points.
!
DO i=2,N-1
	dv=time_step*(c1*(Nz(i)*(u_old(i+1)-u_old(i))-Nz(i-1)* &
    (u_old(i)-u_old(i-1))))
	u_new(i)=u_old(i)+dv
END do
!
!       Bottom depth element with quadratic boundary condition.
!
dv=time_step*(c1*Nz(1)*(u_old(2)-u_old(1))- &
  c2*dsqrt(u_bed**2.0+v_bed**2.0)*u_bed)
u_new(1)=u_old(1)+dv
!
!       Surface depth element with wind stress boundary condition.
!
dv=time_step*(-c1*Nz(N-1)*(u_old(N)-u_old(N-1))+stress/(1025.0*dz))
u_new(N)=u_old(N)+dv
!
RETURN
end
!
!**********************************************************
!       Semi-implicit Coriolis.
!**********************************************************
!
SUBROUTINE EQN_CORIOLIS(u,v,dt,f0,N,alf,alf1,alf2)

implicit none
double precision :: u(200),v(200),dt,alf,alf1,alf2,f0
integer :: i,N

do i=1,N
  u(i)=alf1*(alf2*u(i)+alf*v(i))
  v(i)=alf1*(-alf*u(i)+alf2*v(i))
end do

return
end
!**********************************************************
!       Surface heating and vertical heat distribution.
!**********************************************************
!
SUBROUTINE SURFACE_HEAT()

use physics
use biology
use variables_all
use variables
implicit none
!
double precision :: qs,rad1,rad2,atten,dts1,dts2,hl,sk,ske
double precision :: sat_vap,vap,spec_hum1,spec_hum2,surf_temp,rad_in
integer :: ns

surf_temp=temp_old(N)+273.15
rad_in=rad0
sat_vap=10.0**((0.7859+0.03477*temp_old(N))/(1.0+0.00412*temp_old(N)))
vap=0.01*humid(idmet)*sat_vap
spec_hum1=0.62*sat_vap/(airP(idmet)-0.38*sat_vap)
spec_hum2=0.62*vap/(airP(idmet)-0.38*vap)
!hl=0.985*5.67d-8*(surf_temp**4.0)*(0.39-0.05*vap**0.5)*(1.0-0.6d-4*cloud(idmet)**2.0)	!original codes Longwave
!Net longwave radiation flux=IR_surface_emissivity_1_for_black_body(Stefan_Boltzmann_constant * T**4.8 - downwelling_longwave) Note T should really be skin temp
!hl=0.985*5.67d-8*(surf_temp**4.0)-(lw_rad_down0)	!Longwave
hl=0.98*5.67d-8*(surf_temp**4.0)-(lw_rad_down0*0.95)	!Longwave The 0.95 is
!playing around to account for lw being absorbed and lost in the skin rather
!than top meters
sk=1.45d-3*1.3*1004.0*wind_speed(idmet)*(temp_old(N)-airT(idmet))			!Sensible
ske=1.5d-3*1.3*wind_speed(idmet)*(spec_hum1-spec_hum2)*(2.5d6-2.3d3*temp_old(N))        !Evaporation
rad_out=-(hl+sk+ske)

rad_sum = rad_sum + rad_in+rad_out
acount = acount + 1

qs=rad_out*time_step
dts1=qs/(3900.0*1025.0*dz)	       	      ! Heat loss from sea surface
dts2=0.55*rad_in*time_step/(3900.0*1025.0*dz) ! 55% of heat (red end of spectrum) assumed absorbed in upper bgrid cell
rad1=0.45*rad_in                              ! 45% of heat available for distribution through water column
DO i=n,1,-1
	atten=lambda+heat_shade*x_new(i)	    ! local attenuation coefficient related to local total chl
	rad2=rad1*EXP(-atten*dz)
	temp_new(i)=temp_old(i)+time_step*(rad1-rad2)/(3900.0*1025.0*dz)
	rad1=rad2
!        print*, i,temp_new(i)
END do
temp_new(1)=temp_new(1)+time_step*rad1/(3900.0*1025.0*dz)
temp_new(N)=temp_new(N)+dts1+dts2

RETURN
end
!
!
!**********************************************************
!       Convective mixing.
!**********************************************************
!
SUBROUTINE CONVECTION(old,new,N,ns,nhr,iday)

implicit none
double precision :: old(200),new(200),t_save,meant
integer i,N,i_save,ns,nhr,iday

if(old(N).lt.old(1))then	! short cut for winter convection over whole water column
  t_save=old(N)
  do i=2,N
    new(i)=old(i-1)
  end do
  new(1)=t_save
  i_save=1
else
  i_save=N-1
  do while(old(N).lt.old(i_save))
    i_save=i_save-1
    if(i_save.eq.1)i_save=-1
  end do
  i_save=i_save+1	! This is the level that the surface temp needs to be moved to
  if(i_save.gt.0)then
    meant=sum(old(i_save:N))/dble(N-i_save+1)
    do i=i_save,N
      new(i)=meant
    end do
  end if
end if

return
end
!
!**********************************************************
!       Scalar diffusion - implicit version.
!**********************************************************
!
SUBROUTINE scalar_diff(new,old,time_step,N,cnpar)

use turbulence
implicit none
DOUBLE PRECISION :: old(200),new(200),time_step,cnpar,c,a
double precision :: Y(0:N)
integer :: N,i

do i=1,N
  Y(i)=old(i)
end do
Y(0)=Y(1)

do i=2,N-1
  c=2.*time_step*Kz(i)  /(h(i)+h(i+1))/h(i)
  a=2.*time_step*Kz(i-1)/(h(i)+h(i-1))/h(i)
  cu(i)=-cnpar*c						!i+1,n+1
  au(i)=-cnpar*a						!i-1,n+1
  bu(i)=1-au(i)-cu(i)						!i  ,n+1
  du(i)=Y(i)+(1-cnpar)*(a*Y(i-1)-(a+c)*Y(i)+c*Y(i+1))
end do

!  Surface
a=2.*time_step*Kz(N-1)/(h(N)+h(N-1))/h(N)
au(N)=-cnpar*a
bu(N)=1-au(N)
du(N)=Y(N)+(1-cnpar)*a*(Y(N-1)-Y(N))

!  Bottom
c=2.*time_step*Kz(1)/(h(1)+h(2))/h(1)
cu(1)=-cnpar*c
bu(1)=1-cu(1)
du(1)=Y(1)+(1-cnpar)*c*(Y(2)-Y(1))

call Tridiagonal(N,1,N,Y)

new(1:N)=Y(1:N)

return
end subroutine scalar_diff
!
!*********************************************************
!      Gaussian distribution for wind variability.
!*********************************************************
!
SUBROUTINE gaussian(dist)
!
implicit none
double precision :: average,std
real :: dist(100),dx,x1,twopi,amp,rprob(16),factor
integer :: i,i2,i3,prob(16),itest,xsum

twopi=sqrt(2.0*3.1415927)
amp=1.0/(twopi)
std=1.0
dx=0.25
x1=-2.0
do i=1,16
  rprob(i)=amp*exp(-0.5*(x1)**2.0)*dx
  prob(i)=nint(rprob(i)*100)
  x1=x1+dx
end do
i2=1
factor=-2.0
xsum=0
itest=0
do i=1,16
  do i3=i2,i2+prob(i)-1
    dist(i3)=factor*std
    xsum=xsum+1
    itest=itest+1
  end do
  i2=itest+1
  factor=factor+0.25
end do
dist(xsum:100)=0.0
!
RETURN
end
!
!**********************************************************
!       Wind components, including std variability.
!**********************************************************
!
SUBROUTINE wind_components(wx,wy)
!
use variables
use physics
implicit none
real :: dir2,wx,wy,ws,pi
integer :: iwind

pi=3.1415927

ws=wind_speed(idmet)
dir2=wind_dir(idmet)
dir2=dir2*3.142/180.0

if(dir2.ge.0.0.and.dir2.lt.pi/2.0)then
  wx=ws*sin(dir2)
  wy=ws*cos(dir2)
end if

if(dir2.ge.pi/2.0.and.dir2.lt.pi)then
  wx=ws*cos(dir2-pi/2.0)
  wy=-ws*sin(dir2-pi/2.0)
end if

if(dir2.ge.pi.and.dir2.lt.3.0*pi/2.0)then
  wx=-ws*sin(dir2-pi)
  wy=-ws*cos(dir2-pi)
end if

if(dir2.ge.3.0*pi/2.0.and.dir2.le.2.0*pi)then
  wx=-ws*cos(dir2-3.0*pi/2.0)
  wy=ws*sin(dir2-3.0*pi/2.0)
end if

return
end
!
!**********************************************************
!       Turbulence closure subroutine Canuto's k-e
!**********************************************************
!
SUBROUTINE TURBULENCE_ke(iday,ncount)

use physics
use variables_all
use variables
use turbulence
implicit none

double precision :: z0b,z0s,u_taus,u_taub,LLk,x,rich(200),pot,rich2,iwshear,rich_iw
double precision :: iw_sm,iw_sh,iw_cm,iw_ch,iw_an,iw_as

double precision :: mol_nu=1.3e-6
double precision :: h0b=0.05
double precision :: charnock=1400
double precision :: alpha_wave=0.0
double precision :: crit_rich=0.85

INTEGER iday,ncount

!   calculate friction velocities and roughness:
u_taub=sqrt(stressb/density(1))
u_taus=((stressx/density(N))**2+(stressy/density(N))**2.)**0.25
z0b=(0.1*mol_nu/u_taub)+0.03*h0b
z0s=charnock*u_taus**2./9.81
if(z0s.lt.0.02)z0s=0.02

!   calculate buoyancy frequency, shear, and the production and buoyancy terms
SS(0:N)=0.0
!!HACK!!
NN(0:N)=1.0d-6
SS(0:N)=1.0d-6
do i=1,N-1
  SS(i)=SS(i)+0.5*((cnpar*(velx_new(i+1)-velx_new(i))*(velx_new(i+1)-velx_new(i))+(1.-cnpar)* &
    (velx_new(i+1)-velx_new(i))*(velx_new(i+1)-velx_new(i)))/(0.5*(h(i+1)+h(i)))/h(i)+(cnpar* &
    (velx_new(i+1)-velx_new(i))*(velx_new(i+1)-velx_new(i))+(1.-cnpar)*(velx_new(i+1)-velx_new(i)) &
    *(velx_new(i+1)-velx_new(i)))/(0.5*(h(i+1)+h(i)))/h(i+1))
  SS(i)=SS(i)+0.5*((cnpar*(vely_new(i+1)-vely_new(i))*(vely_new(i+1)-vely_new(i))+(1.-cnpar)* &
    (vely_new(i+1)-vely_new(i))*(vely_new(i+1)-vely_new(i)))/(0.5*(h(i+1)+h(i)))/h(i)+(cnpar* &
    (vely_new(i+1)-vely_new(i))*(vely_new(i+1)-vely_new(i))+(1.-cnpar)*(vely_new(i+1)-vely_new(i)) &
    *(vely_new(i+1)-vely_new(i)))/(0.5*(h(i+1)+h(i)))/h(i+1))

  NN(i)=-0.00957*(density(i+1)-density(i))/dz
  P(i)=Nz(i)*SS(i)
  B(i)=-Kz(i)*NN(i)
!  print *,' i, NN, SS = ',i,NN(i),SS(i)
  Ri(i)=NN(i)/SS(i)
end do
SS(0)=SS(1); NN(0)=NN(1); P(0)=P(1); B(0)=P(1); NN(N)=0.0; SS(N)=SS(N-1); Ri(N)=0.0

call tke_calc(u_taus,u_taub,z0s)

call dissipation(z0b,z0s,u_taus,u_taub) !This gets dissipation eps(i) and lengthscale Lscale(i)

do i=0,N
  LLk=Lscale(i)*Lscale(i)/tke(i)
  as(i)=LLk*SS(i)
  an(i)=LLk*NN(i)
end do

call stability_funcs()

do i=1,N
  x=sqrt(tke(i))*Lscale(i)
  Nz(i)=cmue1(i)*x
  Kz(i)=cmue2(i)*x
end do
Nz(0)=kappa*u_taub*z0b
Nz(N)=kappa*u_taus*z0s
Kz(0)=kappa*u_taub*z0b
Kz(N)=kappa*u_taus*z0s

!   Next lines apply background mixing
where(Nz.lt.Nz_bg)Nz=Nz_bg
where(Kz.lt.Kz_bg)Kz=Kz_bg
! Add molecular viscosity
Nz(0:N)=Nz(0:N)+1.0e-6
Kz(0:N)=Kz(0:N)+1.0e-6

! Limit mixing at Ri<0 to avoid numerical instability
where(Nz.gt.vismax)Nz=vismax
where(Kz.gt.vismax)Kz=vismax

RETURN
end
!
!**********************************************************
!	Calculate dissipation and lengthscale
!**********************************************************
!
SUBROUTINE dissipation(z0b,z0s,u_taus,u_taub)

use physics
use variables_all
use turbulence

implicit none

double precision :: z0b,z0s,u_taus,u_taub,sig_e0,sig_e1,cde
double precision :: cm0=0.527
double precision :: ce1=1.44
double precision :: ce2=1.92
double precision :: galp=0.53
double precision :: avh(0:n),flux(0:n)
double precision :: pminus(0:n),pplus(0:n)
double precision :: prod,buoyan,diss,craig_m
double precision :: cee3,epslim
double precision :: peps,sig_e(0:N)
double precision :: eps_min=5e-10
double precision :: L_min=0.01
double precision :: cm_craig=0.73
logical :: length_lim=.true.

cde=cm0**3.0
craig_m=sqrt(1.5*cm_craig**2*1.0/kappa**2)
sig_e0=(4./3.*craig_m+1.)*(craig_m+1.)*kappa**2/(ce2*cm_craig**2)
sig_e1= kappa*kappa*cm0/(ce2-ce1)/cde
flux(1:N-1)=0.0

! Without surface wave contribution
sig_e(0:N)=sig_e1

do i=1,N
   avh(i)=0.5*(Nz(i-1)/sig_e(i-1)+Nz(i)/sig_e(i))
end do

flux(1)=avh(1)*cde*(tkeold(1)**1.5)/(kappa*(z0b+h(1))**2.)

! Without surface wave impact on flux
flux(N-1)=cmue1(N-1)*sqrt(tkeold(N-1))*kappa*(dz+z0s)/sig_e(N-1)*cde*tkeold(N-1)**1.5/(kappa*(z0s+dz)**2.)

avh(1)=0.0; avh(N)=0.0

do i=1,n-1
  if (B(i).gt.0) then
    cee3=1.
  else
    cee3=-0.629
  end if
  prod=ce1*eps(i)/tkeold(i)*P(i)
  buoyan=cee3*eps(i)/tkeold(i)*B(i)
  diss=ce2*eps(i)*eps(i)/tkeold(i)
  if (prod+buoyan.gt.0) then
    pplus(i)=prod+buoyan
    pminus(i)=diss
  else
    pplus(i)=prod
    pminus(i)=diss-buoyan
  end if
end do

do i=1,N-1
   au(i)=-2.*time_step*avh(i)/(h(i)+h(i+1))/h(i)
   cu(i)=-2.*time_step*avh(i+1)/(h(i)+h(i+1))/h(i+1)
   bu(i)=1.-au(i)-cu(i)+pminus(i)*time_step/eps(i)
   du(i)=(1+pplus(i)*time_step/eps(i))*eps(i)+flux(i)*time_step/(0.5*(h(i)+h(i+1)))
end do

call Tridiagonal(n,1,n-1,eps)
eps(0) = cde*sqrt(tke(0)*tke(0)*tke(0))/kappa/z0b
eps(n) = cde*sqrt(tke(n)*tke(n)*tke(n))/kappa/z0s

do i=0,n
  if ((NN(i).gt.0).and.(length_lim)) then
    epslim=cde/sqrt(2.)/galp*tke(i)*sqrt(NN(i))
  else
    epslim=eps_min
  end if
  if (eps(i).lt.epslim) eps(i)=epslim
  Lscale(i)=cde*sqrt(tke(i)*tke(i)*tke(i))/eps(i)
!  if ( NN(i).eq.0.0d0 ) print *,' duff NN for i = ',i
  if ( NN(i).eq.0.0d0 ) NN = 1.0d-6
  if ( Lscale(i).gt.(0.267*sqrt(2.0*tke(i)/NN(i))) ) Lscale(i)=0.267*sqrt(2.0*tke(i)/NN(i))
  if ( Lscale(i).lt.L_min) Lscale(i)=L_min
end do

RETURN
End
!**********************************************************
!       Calculate stability functions
!**********************************************************
!
SUBROUTINE stability_funcs()

use physics
use variables_all
use variables
use turbulence
implicit none

double precision :: d(6),s(6)
double precision, parameter:: cm0=0.5270
double precision, parameter:: tnmin=-12.27
double precision, parameter:: a2_cm03=2./cm0**3
double precision :: tsmax
double precision :: tn,ts,dd

data d/4.2483000e+002,2.7132000e+001,3.0498889e+000,2.3040000e-001,1.3866000e-001,-8.9538560e-004/
data s/2.2728405e+001,9.2456320e-001,-6.4200000e-003,2.3800000e+001,2.4000000e-001,4.6702411e-002/

do i=1,N-1
  tn=4./cm0**6*an(i)
  if(tn.lt.tnmin)tn=tnmin
  ts=4./cm0**6*as(i)
  tsmax=(d(1)+d(2)*tn+d(4)*tn*tn)/(d(3)+d(5)*tn)
  if(ts.gt.tsmax)ts=tsmax
  dd=d(1)+d(2)*tn+d(3)*ts+d(4)*tn*tn+d(5)*tn*ts+d(6)*ts*ts
  sm(i)=(s(1)+s(2)*tn+s(3)*ts)/dd
  sh(i)=(s(4)+s(5)*tn+s(6)*ts)/dd
  cmue1(i)=a2_cm03*sm(i)
  cmue2(i)=a2_cm03*sh(i)
end do
cmue1(0)=cmue1(1)
cmue1(n)=cmue1(n-1)
cmue2(0)=cmue2(1)
cmue2(n)=cmue2(n-1)

return
end
!**********************************************************
!       Calculate TKE
!**********************************************************
!
SUBROUTINE tke_calc(u_taus,u_taub,z0s)

use turbulence
use variables_all
use variables
implicit none

double precision :: u_taus,u_taub,z0s,numtke,craig_m
double precision :: avh(0:N)
double precision :: pminus(0:N),pplus(0:N)
double precision :: prod,buoyan,diss,cde
double precision :: cm0=0.527
double precision :: k_min=3e-6
double precision :: cm_craig=0.73

tkeold(0:N)=tke(0:N)
cde=cm0**3.
do i=2,N-1
  avh(i)=0.5*(Nz(i-1)+Nz(i))
end do

avh(1)=0
avh(N)=0

!avh(1)=0.5*(Nz(0)+Nz(1))
!avh(N)=0.5*(Nz(N)+Nz(N-1))
!avh(1)=u_taub**4*2./(eps(0)+eps(1))
!avh(N)=u_taus**4*2./(eps(N)+eps(N-1))

do i=N-1,1,-1
  prod=P(i)
  buoyan=B(i)
  diss=eps(i)
  if (prod+buoyan.gt.0) then
    pplus(i)=prod+buoyan
    pminus(i)=diss
  else
    pplus(i)=prod
    pminus(i)=diss-buoyan
  end if
end do

do i=1,N-1
  au(i)=-2.*time_step*avh(i)/(h(i)+h(i+1))/h(i)
  cu(i)=-2.*time_step*avh(i+1)/(h(i)+h(i+1))/h(i+1)
  bu(i)=1.-au(i)-cu(i)+pminus(i)*time_step/tke(i)
  du(i)=(1+pplus(i)*time_step/tke(i))*tke(i)
end do

call tridiagonal(N,1,N-1,tke)

tke(0)=u_taub*u_taub/sqrt(cm0*cde)
tke(N)=u_taus*u_taus/sqrt(cm0*cde)

where(tke.lt.k_min)tke=k_min

return
end
!**********************************************************
!       Tridiaginal matrix solution
!**********************************************************
!
SUBROUTINE tridiagonal(N,fi,lt,value)

use turbulence
implicit none

integer :: N,lt,fi,i
double precision :: value(0:N)

ru(lt)=au(lt)/bu(lt)
qu(lt)=du(lt)/bu(lt)

do i=lt-1,fi+1,-1
  ru(i)=au(i)/(bu(i)-cu(i)*ru(i+1))
  qu(i)=(du(i)-cu(i)*qu(i+1))/(bu(i)-cu(i)*ru(i+1))
end do

qu(fi)=(du(fi)-cu(fi)*qu(fi+1))/(bu(fi)-cu(fi)*ru(fi+1))

value(fi)=qu(fi)
do i=fi+1,lt
  value(i)=qu(i)-ru(i)*value(i-1)
end do

return
end
!**********************************************************
!       Profile of layer mean PAR irradiance
!**********************************************************
!
SUBROUTINE PAR_profile()

use biology
use variables_all
use physics
implicit none
double precision :: rad_surf,d
!

rad_surf=par_percent*rad0
do i=N,1,-1
	d=par_atten+chl_abscross*x_new(i)	       ! local attenuation coefficient related to local chl
	radbar(i)=(rad_surf/(d*dz))*(1.0-DEXP(-d*dz))
	rad_surf=rad_surf*DEXP(-d*dz)
END do
return
end
!
!**********************************************************
!       Primary production subroutine.
!**********************************************************
!
SUBROUTINE PRIMARY_PROD(jd)

use biology
use turbulence
use variables_all
use variables
use physics
implicit none
!
double precision :: mu1,mu2,downsource_ni,downsource_x,upsource_ni,upsource_x,loss_ni,loss_x,w(200)
integer :: jd

IF(x_old(1).lt.seed_biomass*chl_carbon)then   ! limit low values of biomass by seedstock value
  x_old(1)=seed_biomass*chl_carbon
  x_new(1)=seed_biomass*chl_carbon
  ni_old(1)=x_old(1)*sub_quota
  ni_new(1)=ni_old(1)
END if
!
!       MICROBIAL NUTRIENT QUOTA, vertical velocities, growth, uptake and grazing rates....
!
do i=1,n
	quo(i)=ni_new(i)/x_new(i)      ! calculate cell quota
!
  select case(growth_model)
    case(eppley)
      Pmax=0.59*dexp(0.0633*temp_old(i))/(24.0*3600.0)	!Eppley Curve
	  case(Qten)
	    Pmax=Pmax10_persec*Q10**((temp_old(i)-gT0)/10.0)
	end select
  Pm=Pmax*((quo(i)-sub_quota)/(max_quota-sub_quota))
  gross_grow(i)=Pm*(1.0-dexp(-alpha_persec*radbar(i)*chl_carbon/Pm))
  net_grow(i)=gross_grow(i)-(chl_carbon*(respiration_persec*rQ10**((temp_old(i)-rT0)/10.0)))

	uptake(i)=uptake_max_persec*(1.0-(quo(i)/max_quota))*(s_old(i))/(half_saturation+(s_old(i)))
	if(net_grow(i).lt.0.0)uptake(i)=uptake(i)+net_grow(i)*quo(i)
!
!       If it is light, then the upwards vertical velocity is applied
!
  IF(radbar(i).gt.0.0)then
    w(i)=swim_speed_persec
  else
    w(i)=0.0d0
  END if
!
!       then add the sinking speed
!
  w(i)=w(i)+sink_speed_persec
!
!       the next lines limit grazing if x is less than xthresh
!
	if(x_old(i).gt.graze_thresh)then
    graz(i)=g0(jd)
  else
		graz(i)=0.0
  end if
END do
!
!       vertical movement of cell chl and ni
!
do i=2,n-1
  IF(w(i+1).lt.0.0)then
    downsource_ni=-ni_old(i+1)*w(i+1)*time_step/dz
    downsource_x=-x_old(i+1)*w(i+1)*time_step/dz
  else
    downsource_ni=0.0d0
    downsource_x=0.0d0
  END if
  IF(w(i-1).ge.0.0)then
    upsource_ni=ni_old(i-1)*w(i-1)*time_step/dz
    upsource_x=x_old(i-1)*w(i-1)*time_step/dz
  else
    upsource_ni=0.0d0
    upsource_x=0.0d0
  END if
  loss_ni=ABS(w(i))*ni_old(i)*time_step/dz
  ni_new(i)=ni_new(i)+upsource_ni+downsource_ni-loss_ni
  loss_x=ABS(w(i))*x_old(i)*time_step/dz
  x_new(i)=x_new(i)+upsource_x+downsource_x-loss_x
END do
!
!       surface depth cell...
!
IF(w(N-1).ge.0.0)ni_new(N)=ni_new(N)+ni_old(n-1)*w(N-1)*time_step/dz
IF(w(N).lt.0.0)x_new(N)=x_new(N)+x_old(N)*w(N)*time_step/dz
IF(w(N-1).ge.0.0)x_new(N)=x_new(N)+x_old(n-1)*w(N-1)*time_step/dz
IF(w(N).lt.0.0)ni_new(N)=ni_new(N)+ni_old(N)*w(N)*time_step/dz
!
!       bottom depth cell
!
IF(w(2).lt.0.0)ni_new(1)=ni_new(1)-w(2)*ni_old(2)*time_step/dz
IF(w(1).ge.0.0)x_new(1)=x_new(1)-w(1)*x_old(1)*time_step/dz
IF(w(2).lt.0.0)x_new(1)=x_new(1)-w(2)*x_old(2)*time_step/dz
IF(w(1).ge.0.0)ni_new(1)=ni_new(1)-w(1)*ni_old(1)*time_step/dz
!
!       Mixing, growth, and grazing of phytoplankton:
!
!       vertical diffusion of algal cells
!
call scalar_diff(x_new,x_old,time_step,N,cnpar)
call scalar_diff(ni_new,ni_old,time_step,N,cnpar)
!
!       growth of biomass and uptake of nutrient
!
x_new(1:N)=x_new(1:N)+time_step*net_grow(1:N)*x_old(1:N)
ni_new(1:N)=ni_new(1:N)+time_step*uptake(1:N)*x_old(1:N)
!
!       grazing of algal biomass and internal nutrient
!
x_new(1:N)=x_new(1:N)-time_step*graz(1:N)*x_old(1:N)
ni_new(1:N)=ni_new(1:N)-time_step*graz(1:N)*quo(1:N)*x_old(1:N)
!
!       finally, keep track of the removal/addition of nitrate.
!
din_source(1:N)=time_step*(Nrecycle*graz(1:N)*quo(1:N)*x_old(1:N)-uptake(1:N)*x_old(1:N))

do i=1,N
	prod_net(i)=prod_net(i)+time_step*net_grow(i)*x_new(i)/chl_carbon
	prod_gross(i)=prod_gross(i)+time_step*gross_grow(i)*x_new(i)/chl_carbon
	prod_net_daily(i)=prod_net_daily(i)+time_step*net_grow(i)*x_new(i)/chl_carbon
	prod_gross_daily(i)=prod_gross_daily(i)+time_step*gross_grow(i)*x_new(i)/chl_carbon
end do

RETURN
end
!
!**********************************************************
!       Updating dissolved nutrient after phytoplankton utilisation
!**********************************************************
subroutine DIN_update()
!
use variables_all
use biology
use variables
implicit none

integer :: j
!
!       plus grazing and minus uptake
!
  s_new(1:N)=s_new(1:N)+din_source(1:N)
!
!       and then update....
!
s_old(1:N)=s_new(1:N)
din_source(1:N)=0.0
!
return
end
!
!**************************************************************************
subroutine results_report()

use variables

implicit none

     LOGICAL           :: QUIT  = .FALSE.
     integer :: itype

      Character (len=6) :: tpn1_txt,tpg1_txt
      Character (len=100) :: info1,info2
!
      info1='Mean annual gross production: '//tpg1_txt
      info2='Mean annual net production:   '//tpn1_txt
!
!
!   show the data info
!
      RETURN
      END SUBROUTINE results_report
!

!*****************************************************************************
!********************    END OF ENTIRE PROGRAMME ! ************************
!**************************************************************************
