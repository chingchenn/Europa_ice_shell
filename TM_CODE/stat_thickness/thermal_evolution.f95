!-----------------------------------------------------------------------------------------------|
!  Calculates ice layer properties and thermal evolution of icy bodies.                         |
!                                                                                               |
!  - First calculates ice layer properties as a function of ice shell thickness, with results   |
!    written in "_shell-properties.dat" files.                                                  |
!  - Then calculates thermal evolution with results written in "_thermal-evolution.dat" files.  |
!  - Finally, calculates temperature profile in ice layer at different timesteps, with results  |
!    written in _T-profiles.dat.                                                                |
!                                                                                               |
!  This version (15-04-2021) by Frederic Deschamps with heat diffusion solver by Henri Samuel   |
!                                                                                               |
!  Compile:   gfortran -O4 -o thermal_evolution thermal_evolution_mod.f95 thermal_evolution.f95 |
!  Execute:   thermal_evolution                                                                 |
!                                                                                               |
!  Input parameters in file thermal_evolution_param. See details therein.                       |
!-----------------------------------------------------------------------------------------------|
      program thermal_evolution

      use thermal_evolution_mod

      implicit real(8) (a-h,o-z)

      character(len=120) :: casename,filename
      integer, parameter :: nimax=10
      logical :: convecting,stat_thickness
      integer :: heat_param,heat_timevar
      integer :: convec
      real(8) :: km,kref,Lsolid
      real(8) :: rx(nimax),drxdt(nimax),rscal(nimax),rout(nimax)
      real(8) :: r_body(npzmax),T_body(npzmax)
      real(8) :: time_frame(npwmax),r_write(npwmax,npzmax),T_write(npwmax,npzmax)
      real(8), allocatable :: Htidal_r(:)

      internal_heating = .false.
      convecting       = .false.
      stat_thickness   = .false.

!---> Reads entry parameters from file thermal_evolution_param
      open(1,file='thermal_evolution_param',status='unknown')
      write(*,*) 'Satellite properties: radius (km), density (kg/m^3), surface temperature (K), gravity acceleration (m/s^2)'
      read(1,*) R_tot
      read(1,*) rho_m
      read(1,*) Tsurf
      read(1,*) g
      write(*,*) 'Silicate core properties: radius (km), density (kg/m^3), diffusivity (m^2/s), heat capacity (J/kg/K)'
      read(1,*) R_core
      read(1,*) rho_core
      read(1,*) kappa_core
      read(1,*) Cp_core
      write(*,*) 'Chondritic model for initial radiogenic concentration: '
      write(*,*) '     (1) Lodders (2003); (2) CI; (3) CM; (4) LL'
      read(1,*) chondritic_model
      write(*,*) 'Reference vistosity (Pas) and activation energy (kJ/mol)'
      read(1,*) etaref
      read(1,*) Eact
      write(*,*) 'Scaling laws for isoviscous (1) or stagnant-lid (2) fluid ?'
      read(1,*) scaling_law
      write(*,*) 'Temperature-dependent thermal conductivity (.false./.true.)' 
      read(1,*) variable_conductivity
      write(*,*) 'Reference conductivity (if constant conductivity: surface=1, interior=2, manual=3)'
      read(1,*) reference_conductivity
      read(1,*) k_input
      write(*,*) 'Scaling laws for 3D-Cartesian (1) or spherical (2) geometry ?'
      read(1,*) geometry
      write(*,*) 'Keep internal heating rate (1) or total power (2) constant with ice shell thickness'
      read(1,*) heat_param
      if (heat_param == 1) then 
         write(*,*) 'Internal heating rate (W/kg) dissipated throughout the body'
         read(1,*) H0_total
      else
         write(*,*) 'Total power (TW) dissipated throughout the body'
         read(1,*) P0_total
         P0_total  = P0_total*1.0d12
      endif
      write(*,*) 'Fraction of tidal heating within ice ;layer and within core (sum should be 1.0)' 
      read(1,*) Hfrac_ice,Hfrac_core
      if (heat_param == 1) then 
         H0_ice  = H0_total*Hfrac_ice
         H0_core = H0_total*Hfrac_core
      else
         P0_ice  = P0_total*Hfrac_ice
         P0_core = P0_total*Hfrac_core
      endif
      write(*,*) 'Tidal heating time-dependence :'
      write(*,*) '   - Constant . . . . . . . . . . . . . . . . . . . . 1'
      write(*,*) '   - Sinusoidal . . . . . . . . . . . . . . . . . . . 2'
      write(*,*) '   - Sinusoidal + exponential decrease  . . . . . . . 3'
      write(*,*) '   - Sinusoidal-square  . . . . . . . . . . . . . . . 4'
      write(*,*) '   - Sinusoidal-square + exponential decrease . . . . 5'
      read(1,*) heat_timevar
      write(*,*) 'Period for tidal heating time-variation in Gyr, regular decay constant'
      read(1,*) Theat_Gyr,decay_heat
      write(*,*) 'Eccentricity min and max (relative to current eccentricity)'
      read(1,*) emin,emax
      write(*,*) 'Initial fractions of ammonia (NH3) and methanol (CH3OH) in vol %'
      read(1,*) pNH3_init
      read(1,*) pCH3OH_init
      write(*,*) 'Number of point in iceI layer + ocean (max = ',npimax,')'
      read(1,*) npz
      write(*,*) 'Initial thickness of ice I layer (km)'
      read(1,*) z_bot_init
      write(*,*) 'Begining and finishing times (yr), number of timestep, number of output frames for T-profile'
      read(1,*) time_start,time_end,npt
      write(*,*) 'Timescale: (1) Linear; (2) Logarithmic'
      read(1,*) itimescale
      write(*,*) 'Calculate average, min and max in ice shell over given time-window (.faslse./.true.) ?'
      read(1,*) stat_thickness
      write(*,*) 'Starting and ending time of time window (in Gyr)'
      read(1,*) tstart_Gyr,tend_Gyr
      write(*,*) 'Number of output frames for T-profile (max = 51)'
      read(1,*) npw
      write(*,*) 'Core heat flux from Kirk and Stevenson (1987) (1) or by Solving heat diffusion equation (2)'
      read(1,*) core_heatflux_mode
      write(*,*) 'Integration method :'
      write(*,*) '   - 4th order Runge-Kutta . . . . . . 1'
      write(*,*) '   - Adaptative Runge-Kutta  . . . . . 2'
      write(*,*) '   - Bulirsh-Stoer . . . . . . . . . . 3'
      read(1,*) iint
      write(*,*) 'Case name:'
      read(1,'(a)') casename
      close(1)

!---> Various definitions and switches
      call nbrcarac(120,casename,nstc)

      npframe = nint(dble(npt - 1)/dble(npw - 1))

! Radii and volumes
      R_tot        = R_tot*1.0d3
      R_core       = R_core*1.0d3
      Volume_core  = 4.0d0*pi*R_core**3.0d0/3.0d0
      Surface_core = 4.0d0*pi*R_core**2.0d0
      z_bot_init   = z_bot_init*1.0d3
      r_bot_init   = R_tot - z_bot_init
      Volume_init  = 4.0d0*pi*(R_tot**3.0d0 - R_core**3.0d0)/3.0d0
      Volume_iceI  = 4.0d0*pi*(R_tot**3.0d0 - r_bot_init**3.0d0)/3.0d0

! Set activation energy to J/mol (input is kJ/mol)
      Eact = Eact*1.0d3

! Core and ice thermodynamic properties
      k_core      = kappa_core*Cp_core*rho_core
      k_iceI      = k_iceI_a/Tsurf + k_iceI_b
      Cp_iceI     = Cp_iceI_a*Tsurf + Cp_iceI_b
      kappa_iceI  = k_iceI/rho_iceI/Cp_iceI
      kappa_input = k_input/rho_iceI/Cp_iceI

      write(*,*)
      write(*,*) '|----------- >   Surface properties :'
      write(*,1000) Tsurf
      write(*,1005) g
      write(*,1010) k_iceI,Cp_iceI,kappa_iceI
      write(*,*)

! Precision & order of ODE for time integration
      eps = 1.0d-3
      n_int = 1

! Switch for ocean initial composition
      pNH3_init   = pNH3_init/1.0d2
      pCH3OH_init = pCH3OH_init/1.0d2
      if ((pNH3_init == 0.0d0).and.(pCH3OH_init == 0.0d0)) then
         composition = 1
         write(*,*) 'Pure H2O ocean'
         endif
      if ((pNH3_init /= 0.0d0).and.(pCH3OH_init == 0.0d0)) then 
         composition     = 2
         pCLA_init = pNH3_init
         pCLAmax   = pNH3max
         write(*,*) 'H2O-NH3 system'
      endif
      if ((pNH3_init == 0.0d0).and.(pCH3OH_init /= 0.0d0)) then
         composition     = 3
         pCLA_init = pCH3OH_init
         pCLAmax   = pCH3OHmax
      endif

!---> Calculates heat flux and other properties as a function of inner radius of the ice I layer
!     Can be removed if to calculate thermal evolution
      dr_bot           = (R_tot - R_core)/dble(npz-1)
      r_bot            = R_tot - 0.1d0*dr_bot
      Volume_IceI_init = 4.0d0*pi*(R_tot**3.0d0 - r_bot**3.0d0)/3.0d0

      filename=casename(1:nstc)//'_shell-properties.dat'
      open(11,file=filename,status='unknown')
      write(11,1100) R_tot/1.0d3,R_core/1.0d3,Tsurf,g
      write(11,1105)
      do 100 iz=1,npz

         r_bot = R_tot - dr_bot*dble(iz-1)
         z_bot = R_tot - r_bot
         Volume_iceI = 4.0d0*pi*(R_tot**3.0d0 - r_bot**3.0d0)/3.0d0

         if (heat_param == 1) then
            Hint_ice  = H0_ice
            Pint_ice  = H0_ice*rho_iceI*Volume_iceI
            Hint_core = H0_core
            Pint_core = H0_core*rho_core*Volume_core
         else
            Pint_ice  = P0_ice
            if (Volume_iceI <= 0.0d0) Volume_iceI = Volume_IceI_init
            Hint_ice  = P0_ice/rho_iceI/Volume_iceI
            Pint_core = P0_core
            Hint_core = P0_core/rho_core/Volume_core
         endif

         if (Hint_ice <= 0.0d0) then
            internal_heating = .false.
         else
            internal_heating = .true.
         endif

! Temperature at bottom of iceI layer, and temperature jump
         call Temperature_bottom_iceI_shell(z_bot,pres,pCLA,Tbot,dTbotdr,DeltaT)

! Heat flux, ice layer interior temperature, stagnant lid thickness
         call HeatFlux_bottom_iceI_shell(convecting,z_bot,Tbot,DeltaT,Tm,etam,km,alfam,Ra,Hadim,Hcrit,     &
                                                               Fcond,Fbot,Ftop,dlid,z_topTBL,z_botTBL)

!---> Test to check whether the ice layers melts at its bottom
         if (convecting) then
            d_botTBL = z_bot - z_botTBL
            if (d_botTBL >= 0.0d0) then
               zz    = z_bot - 0.5d0*d_botTBL
               dT_zz = (Tbot - Tm)/d_botTBL
               Tzz   = Tm + dT_zz*(zz - z_botTBL)
            else
               Tzz = Tm
            endif
         else
            zz    = z_bot - 2.0d0
            dT_zz = (Tbot - Tsurf)/z_bot
            Tzz   = Tbot + dT_zz*(zz - z_bot)
         endif
         call T_fusion_H2O(zz,pres,TH2O,dTbotdr)

         imelt = 0
         if (Tzz > TH2O) imelt = 1

         Pres_MPa    = pres*1.0d3
         zbot_km     = z_bot/1.0d3
         dlid_km     = dlid/1.0d3
         z_topTBL_km = z_topTBL/1.0d3
         z_botTBL_km = z_botTBL/1.0d3
         Fcond_mW    = Fcond*1.0d3
         Fbot_mW     = Fbot*1.0d3
         Ftop_mW     = Ftop*1.0d3

         f = r_bot/R_tot
         if (f /= 0.0) then
            Fcond_bot_mW = Fcond_mW/f
         else
            Fcond_bot_mW = 0.0d0
         endif
         Fcond_top_mW = Fcond_mW*f

         if (convecting) then
            convec = 1
         else
            convec = 0
         endif
         write(*,1020) zbot_km,Pres_MPa,Hint_ice,convec,imelt,pCLA,Tbot,Tm,etam,km,Fcond_bot_mW,Fcond_top_mW,Fbot_mW,     &
                                                                              Ftop_mW,dlid_km,z_topTBL_km,z_botTBL_km
         write(11,1110) zbot_km,Pres_MPa,Hint_ice,convec,imelt,pCLA,Tbot,Tm,etam,km,Fcond_bot_mW,Fcond_top_mW,Fbot_mW,    &
                                                                               Ftop_mW,dlid_km,z_topTBL_km,z_botTBL_km

100   continue
      close(11)
      write(*,*)

!---> Now calculates time evolution of radial structure
      if (itimescale == 1) then
         dtime     = (time_end - time_start)/dble(npt-1)
         dtime_sec = dtime*year_sec
      else
         time_start_log = dlog10(time_start)
         time_end_log   = dlog10(time_end)
         dtime_log      = (time_end_log - time_start_log)/dble(npt-1)
      endif

      r_bot     = r_bot_init
      dr_init   = 2.0d3
      r0_core   = 1.0d3
      npc       = nint((R_core - r0_core)/dr_init) + 1
      if (npc > npcmax) npc = npcmax
      dr_core      = (R_core - r0_core)/dble(npc-1)
      npr_TBL_min  = 5

      allocate(Htidal_r(npc))

      T_body(:)    = 0.0d0
      r_body(:)    = 0.0d0
      r_write(:,:) = R_tot
      T_write(:,:) = Tsurf

      icount       = 0
      iframe       = 0

      a0_ecc = 0.5d0*(emax + emin)
      a1_ecc = 0.5d0*(emax - emin)

      npt_stat  = 0
      time_stat = 0.0d0

      filename=casename(1:nstc)//'_thermal-evolution.dat'
      open(12,file=filename,status='unknown')
      write(12,1200) R_tot/1.0d3,R_core/1.0d3,Tsurf,g
      write(12,1205)
      do 200 i=1,npt

         if (itimescale == 1) then
            time = time_start + dtime*dble(i-1)
         else
            time_log = time_start_log + dtime_log*dble(i-1)
            time     = 10**time_log
            if (i == 1) then
               dtime = time
            else
               dtime = time - time_save
            endif
            dtime_sec = dtime*year_sec 
         endif
         time_save = time

         icount = icount + 1
         iwrite = 0
         if ((i == 1).or.(i == npt).or.(icount == npframe)) then
            iwrite             = 1            
            iframe             = iframe + 1
            time_frame(iframe) = time
            if (icount == npframe) icount = 0
         endif

         z_bot = R_tot - r_bot
         Volume_iceI = 4.0d0*pi*(R_tot**3.0d0 - r_bot**3.0d0)/3.0d0

         if (heat_param == 1) then
            P0_ice  = H0_ice*rho_iceI*Volume_iceI
            P0_core = H0_core*rho_core*Volume_core
         else
            H0_ice  = 0.0d0
            H0_core = 0.0d0
            if (Volume_iceI > 0.0d0) H0_ice = P0_ice/rho_iceI/Volume_iceI
            if (Volume_core > 0.0d0) H0_core = P0_core/rho_core/Volume_core
         endif

         call Time_dependent_tidal_heating(heat_timevar,time,Theat_Gyr,a0_ecc,a1_ecc,decay_heat,P0_ice,Pint_ice)
         call Time_dependent_tidal_heating(heat_timevar,time,Theat_Gyr,a0_ecc,a1_ecc,decay_heat,P0_core,Pint_core)
         Hint_ice  = 0.0d0
         Hint_core = 0.0d0
         if (Volume_iceI > 0.0d0) Hint_ice = Pint_ice/rho_iceI/Volume_iceI
         if (Volume_core > 0.0d0) Hint_core = Pint_core/rho_core/Volume_core

         if (Hint_ice <= 0.0d0) then
            internal_heating = .false.
         else
            internal_heating = .true.
         endif

! Temperature at bottom of iceI layer, and temperature jump
         call Temperature_bottom_iceI_shell(z_bot,pres,pCLA,Tbot,dTbotdr,DeltaT)

! Radiogenic heat in the core, and heat flux at the top of the core following Kirk and Stevenson (1987) (F_radio_KS)
         call Radiogenic_heat_core(time,H_core,F_radio_KS,DT_core)

! Radiogenic and tidal heat in the core in W/m^3, Hrad and Htidal
         Hrad   = H_core*rho_core
         Htidal = Hint_core*rho_core

! Radial profile of temperature in the core calculated by solving heat diffusion equation (T_core).
!    - Note: input heating rates Hrad and Htidal should be in W/m^3
!    - Initial condition is set to ice layer bottom temperature (Tbot) at t= time_start
!    - Top boundary condition is set to Tbot at each time step
         Tcs = Tbot
         if (i == 1) T_core(1:npc) = Tbot
         call compute_Tcore_profile(npc,dr_core,dtime_sec,Tcs,Hrad,Htidal,k_core,rho_core,Cp_core,r0_core,T_core(1:npc))

! Heat flux at the top of the core :
!    - Radiogenic heat flux calculated from Kirk and Stevenson (1987)
!    - Heat flux calculated from the temperature profile calculated by solving heat diffusion equation
         if (core_heatflux_mode == 1) then
            F_core = F_radio_KS + Htidal*R_core/3.0d0
         else
            F_core = k_core*(-T_core(npc-2) + 4.0d0*T_core(npc-1) - 3.0d0*T_core(npc))/2.0d0/dr_core
         endif
         P_core    = Surface_core*F_core
         Prad_TW   = Hrad*Volume_core/1.0d12
         Ptidal_TW = Htidal*Volume_core/1.0d12

! Adiabatic temperature increase in ocean
         Tadia    = Tbot*(1.0d0 - alfa_H2O*rhog*(R_core - r_bot)/rho_H2O/Cp_H2O)
         dTadiadr = -Tbot*alfa_H2O*rhog/rho_H2O/Cp_H2O

! Heat flux at the bottom of ice I layer
         call HeatFlux_bottom_iceI_shell(convecting,z_bot,Tbot,DeltaT,Tm,etam,km,alfam,Ra,Hadim,Hcrit,     &
                                                               Fcond,Fbot,Ftop,dlid,z_topTBL,z_botTBL)

!---> Test to check whether the ice layers melts at its bottom
         if (convecting) then
            d_botTBL = z_bot - z_botTBL
            if (d_botTBL >= 0.0d0) then
               zz    = z_bot - 0.5d0*d_botTBL
               dT_zz = (Tbot - Tm)/d_botTBL
               Tzz   = Tm + dT_zz*(zz - z_botTBL)
            else
               Tzz = Tm
            endif
         else
            zz    = z_bot - 2.0d0
            dT_zz = (Tbot - Tsurf)/z_bot
            Tzz   = Tbot + dT_zz*(zz - z_bot)
         endif
         call T_fusion_H2O(zz,pres,TH2O,dTbotdr)

         imelt = 0
!!!         if (Tzz > TH2O) then
!!!
!!!            imelt   = 1
!!!            dtime   = dtime_save
!!!            r_bot   = r_bot_save
!!!            drbotdt = drbotdt_save
!!!
!!!         else

            if (Tzz > TH2O) imelt = 1

            dtime_save   = dtime
            r_bot_save   = r_bot
            drbotdt_save = drbotdt

            f = r_bot/R_tot
            if (f /= 0.0) then
               Fcond_bot = Fcond/f
            else
               Fcond_bot = 0.0d0
            endif
            Fcond_top_mW = Fcond*f

            dPower     = r_bot**2.0d0*Fbot - R_core**2.0d0*F_core 
            cool_ocean = rho_H2O*Cp_H2O*(dTbotdr - dTadiadr)*(r_bot**3.0d0 - R_core**3.0d0)/3.0d0
            Lsolid     = rho_iceI*L_iceI*r_bot**2.0d0
            drbotdt    = dPower/(cool_ocean - Lsolid)*year_sec

            time_Gyr     = time/1.0d9
            dtime_Gyr    = dtime/1.0d9
            P_TW         = P_core/1.0d12
            Pint_ice_TW  = Pint_ice/1.0d12
            Pres_MPa     = pres*1.0d3
            zbot_km      = z_bot/1.0d3
            rbot_km      = r_bot/1.0d3
            Fcore_mW     = F_core*1.0d3
            Fcond_mW     = Fcond*1.0d3
            Fbot_mW      = Fbot*1.0d3
            Ftop_mW      = Ftop*1.0d3
            dlid_km      = dlid/1.0d3
            z_topTBL_km  = z_topTBL/1.0d3

!!!            write(*,1030) dTbotdr,dTadiadr,dPower,cool_ocean,Lsolid,drbotdt
!!!         endif

!---> Calculates anchor points for temperature profile

!     Reference conductivity
         if (variable_conductivity) then
            kref = k_iceI
         else
            if (reference_conductivity == 1) then
               kref = k_iceI
            else if (reference_conductivity == 2) then
               kref = km
            else if (reference_conductivity == 3) then
               kref = k_input
            endif
         endif

         rlid     = 0.0d0
         r_topTBL = 0.0d0
         r_botTBL = 0.0d0
         rbot_ice = rbot_km*1.0d3
         if (convecting) then
            rlid  = R_tot - dlid_km*1.0d3
            Fsurf = Ftop_mW/1.0d3
            Fbase = Fbot_mW/1.0d3
            if (.not.variable_conductivity) then
               tcond_TBL = kref
               if (geometry == 1) then
                  d_topTBL = tcond_TBL*(Tm - Tsurf)/Fsurf
                  r_topTBL = R_tot - d_topTBL
                  d_botTBL = tcond_TBL*(Tbot - Tm)/Fbase
                  r_botTBL = rbot_ice + d_botTBL
                  if (r_botTBL <= rbot_ice) r_botTBL = rbot_ice
               else
                  f        = r_bot/R_tot
                  Flin     = tcond_TBL*(Tm - Tsurf)/R_tot
                  r_topTBL = R_tot*Fsurf/(Fsurf + Flin)
                  Flin     = tcond_TBL*(Tbot - Tm)/rbot_ice
                  r_botTBL = rbot_ice*Fbase/(Fbase - Flin)
                  if (r_botTBL <= rbot_ice) r_botTBL = rbot_ice
               endif
            else
               tcond_topTBL = kref
               tcond_botTBL = k_iceI_a/Tm + k_iceI_b
               if (geometry == 1) then
                  d_topTBL = tcond_topTBL*Tsurf*dlog(Tm/Tsurf)/Fsurf
                  r_topTBL = R_tot - d_topTBL
                  d_botTBL = tcond_botTBL*Tm*dlog(Tbot/Tm)/Fbase
                  r_botTBL = rbot_ice + d_botTBL
               else
                  f        = r_bot/R_tot
                  Flin     = tcond_topTBL*(Tm - Tsurf)/R_tot
                  r_topTBL = R_tot*Fsurf/(Fsurf + Flin)
                  Flin     = tcond_botTBL*(Tbot - Tm)/rbot_ice
                  r_botTBL = rbot_ice*Fbase/(Fbase - Flin)
               endif
            endif
         endif
         rlid_km     = rlid/1.0d3
         r_topTBL_km = r_topTBL/1.0d3
         r_botTBL_km = r_botTBL/1.0d3

!---> Calculates temperature profile in planet/body

         rbot_ice = rbot_km*1.0d3
         Tice     = Tm
         Tice_bot = Tbot

!     Temperature profile within the silcate core
         r_body(1) = 0.0d0
         T_body(1) = T_core(1)
         if (iwrite == 1) then
            T_write(iframe,1) = T_body(1)
            r_write(iframe,1) = r_body(1)
         endif
         do 110 ic=1,npc
            r_body(ic+1) = r0_core + dr_core*dble(ic-1)
            T_body(ic+1) = T_core(ic)
            if (iwrite == 1) then
               T_write(iframe,ic+1) = T_body(ic+1)
               r_write(iframe,ic+1) = r_body(ic+1)
            endif
110      continue
         npr = npc + 1

!     Temperature within ice layer
         if (convecting) then

            npr_botTBL = nint((r_botTBL - rbot_ice)/dr_init) + 1
            if (npr_botTBL < npr_TBL_min) npr_botTBL = npr_TBL_min
            if (npr_botTBL > npimax) npr_botTBL = npimax
            npr_topTBL = nint((R_tot - r_topTBL)/dr_init) + 1
            if (npr_topTBL < npr_TBL_min) npr_topTBL = npr_TBL_min
            if (npr_topTBL > npimax) npr_topTBL = npimax

!     Temperature profile within bottom TBL
            h_TBL  = r_botTBL - rbot_ice
            DTemp  = Tice_bot - Tice
            dr_TBL = 0.0d0
            if (npr_botTBL > 1) dr_TBL = h_TBL/dble(npr_botTBL-1)
            do 120 ipt=1,npr_botTBL
               r_body(npr+ipt) = rbot_ice + dr_TBL*dble(ipt-1)
               xrad            = r_botTBL/r_body(npr+ipt)
               one_xrad        = 1.0d0 - xrad
               f_TBL           = 1.0d0/xrad
               aux_H_1         = f_TBL*(1.0d0 + f_TBL)*one_xrad
               aux_H_2         = 1.0d0 - f_TBL**2.0d0
               T_Hint          = rho_iceI*Hint_ice*R_tot**2.0d0*(aux_H_1 + aux_H_2)/6.0d0/kref
               T_body(npr+ipt) = Tice - DTemp*rbot_ice*(1.0d0 - xrad)/h_TBL + T_Hint
               if (iwrite == 1) then
                  T_write(iframe,npr+ipt) = T_body(npr+ipt)
                  r_write(iframe,npr+ipt) = r_body(npr+ipt)
               endif
120         continue
            npr = npr + npr_botTBL

!     Temperature profile within top TBL + stagnant lid
            h_TBL  = R_tot - r_topTBL
            f_TBL  = r_topTBL/R_tot
            DTemp  = Tice - Tsurf
            dr_TBL = 0.0d0
            if (npr_topTBL > 1) dr_TBL = h_TBL/dble(npr_topTBL - 1)
            do 130 ipt=1,npr_topTBL
               r_body(npr+ipt) = r_topTBL + dr_TBL*dble(ipt-1)
               xrad            = R_tot/r_body(npr+ipt)
               one_xrad        = 1.0d0 - xrad
               rxrad           = 1.0d0/xrad
               aux_H_1         = f_TBL*(1.0d0 + f_TBL)*one_xrad
               aux_H_2         = 1.0d0 - rxrad**2.0d0
               T_Hint          = rho_iceI*Hint_ice*R_tot**2.0d0*(aux_H_1 + aux_H_2)/6.0d0/kref
               if (.not.variable_conductivity) then
                  T_body(npr+ipt) = Tsurf - DTemp*r_topTBL*one_xrad/h_TBL + T_Hint
               else
                  aux_k_1         = -R_tot*f_TBL*one_xrad*dlog(Tice/Tsurf)/h_TBL
                  aux_k_2         = T_Hint/Tsurf
                  T_body(npr+ipt) = Tsurf*dexp(aux_k_1 + aux_k_2)
               endif
               if (iwrite == 1) then
                  T_write(iframe,npr+ipt) = T_body(npr+ipt)
                  r_write(iframe,npr+ipt) = r_body(npr+ipt)
               endif
130         continue
            npr_tot = npr + npr_topTBL

         else

!     If ice layer does not convect, temperature profile is conductive throughout the layer.
            npr_ice = nint((R_tot - rbot_ice)/dr_init) + 1
            if (npr_ice < npr_TBL_min) npr_ice = npr_TBL_min
            if (npr_ice > npimax) npr_ice = 2*npimax

            h_ice  = R_tot - rbot_ice
            f_ice  = rbot_ice/R_tot
            DTemp  = Tice_bot - Tsurf
            dr_ice = 0.0d0
            if (npr_ice > 1) dr_ice = h_ice/dble(npr_ice - 1)
            do 140 ipt=1,npr_ice
               r_body(npr+ipt) = rbot_ice + dr_ice*dble(ipt-1)
               xrad            = R_tot/r_body(npr+ipt)
               one_xrad        = 1.0d0 - xrad
               rxrad           = 1.0d0/xrad
               aux_H_1         = f_ice*(1.0d0 + f_ice)*one_xrad
               aux_H_2         = 1.0d0 - rxrad**2.0d0
               T_Hint          = rho_iceI*Hint_ice*R_tot**2.0d0*(aux_H_1 + aux_H_2)/6.0d0/kref
               if (.not.variable_conductivity) then
                  T_body(npr+ipt) = Tsurf - DTemp*rbot_ice*one_xrad/h_ice + T_Hint
               else
                  aux_k_1         = -R_tot*f_ice*one_xrad*dlog(Tice_bot/Tsurf)/h_ice
                  aux_k_2         = T_Hint/Tsurf
                  T_body(npr+ipt) = Tsurf*dexp(aux_k_1 + aux_k_2)
               endif
               if (iwrite == 1) then
                  T_write(iframe,npr+ipt) = T_body(npr+ipt)
                  r_write(iframe,npr+ipt) = r_body(npr+ipt)
               endif
140         continue
            npr_tot = npr + npr_ice

         endif

         if (i == 1) then
            npr_write = npr_tot
         else
            if (npr_tot > npr_write) npr_write = npr_tot
         endif

!---> Integrate bottom radius to next point according to chosen method
         rx(1)    = r_bot
         drxdt(1) = drbotdt

         call integrate(iint,n_int,eps,time,dtime,rx,drxdt)
         r_bot = rx(1)
         if (r_bot > R_core) then
            r_bot = rx(1)
         else
            r_bot = R_core
         endif

         if (convecting) then
            convec = 1
         else
            convec = 0
         endif

!---> Average, min and max in ice shell thickness within prescribed time-window
         if ((stat_thickness).and.(time_Gyr >= tstart_Gyr).and.(time_Gyr <= tend_Gyr)) then
            npt_stat  = npt_stat + 1
            time_stat = time_stat + dtime_Gyr
            if (npt_stat == 1) then
               zbot_km_min    = zbot_km
               zbot_km_max    = zbot_km
               zbot_km_avg    = 0.0d0
               dlid_km_min    = dlid_km
               dlid_km_max    = dlid_km
               dlid_km_avg    = 0.0d0
            else
               if (zbot_km < zbot_km_min) zbot_km_min = zbot_km
               if (zbot_km > zbot_km_max) zbot_km_max = zbot_km
               if (dlid_km < dlid_km_min) dlid_km_min = dlid_km
               if (dlid_km > dlid_km_max) dlid_km_max = dlid_km
            endif
            zbot_km_avg = zbot_km_avg + zbot_km*dtime_Gyr
            dlid_km_avg = dlid_km_avg + dlid_km*dtime_Gyr
         endif

         write(*,1040) time_Gyr,Prad_TW,Ptidal_TW,DT_core,Fcore_mW,Pint_ice_TW,Hint_ice,convec,imelt,Pres_MPa,zbot_km,pCLA,    &
                                                                                  Tbot,Tm,Fbot_mW,Ftop_mW,dlid_km,T_core(1)
         write(12,1210) time_Gyr,Prad_TW,Ptidal_TW,Fcore_mW,Pint_ice_TW,Hint_ice,convec,imelt,Pres_MPa,zbot_km,pCLA,Tbot,Tm,   &
                                                                                           Fbot_mW,Ftop_mW,dlid_km,T_core(1)
         
200   continue
      close(12)

      rnpt_stat   = 1.0d0/dble(npt_stat)
      zbot_km_avg = zbot_km_avg/time_stat
      dlid_km_avg = dlid_km_avg/time_stat
      write(*,*)
      write(*,1050) tstart_Gyr,tend_Gyr,zbot_km_avg,zbot_km_min,zbot_km_max
      write(*,1060) tstart_Gyr,tend_Gyr,dlid_km_avg,dlid_km_min,dlid_km_max

      if (iframe /= npw) npw = iframe
      filename=casename(1:nstc)//'_T-profiles.dat'
      open(15,file=filename,status='unknown')
      write(15,1500) R_tot/1.0d3,R_core/1.0d3,Tsurf,rho_core,k_core
      write(15,1505) (time_frame(iframe)/1.0d9,iframe=1,npw)
      do 250 ir = 1,npr_write
         write(15,1510) (r_write(iframe,ir)/1.0d3,T_write(iframe,ir),iframe=1,npw)
250   continue
      close(15)

1000  format(' Tsurf = ',f8.2,' K')
1005  format(' g = ',f6.3,' m/s^2')
1010  format(' k_I = ',f7.3,' W/m/K ;  Cp_I = ',f8.2,' J/kg/K ; kappa_I = ',e12.7,' m^2/s')
1020  format(2(1x,f7.2),1x,d11.5,2(1x,i2),1x,f8.5,2(1x,f7.2),1x,d11.5,1x,f6.3,4(1x,f8.3),3(1x,f7.2))
1030  format(6(1x,d11.5))
1040  format(1x,d11.5,2(1x,f8.4),1x,f8.2,1x,f8.3,1x,f8.4,1x,d11.5,2(1x,i2),1x,f7.2,1x,f8.3,1x,f8.5,2(1x,f7.2),            &
                                                                                   2(1x,f8.3),1x,f7.2,1x,f8.2)
1050  format('Average, min and max in ice shell thickness between ',f5.2,' and ',f5.2,' Gyr    :',3(1x,f7.2))
1060  format('Average, min and max in stagnant lid thickness between ',f5.2,' and ',f5.2,' Gyr :',3(1x,f7.2))
1100  format('R_tot =',f7.1,' R_core = ',f7.1,'  Tsurf = ',f8.2,'  g = ',f6.3)
1105  format(3x,'zbot',4x,'P(MPa)',1x,'Hint',5x,'conv',1x,'melt',1x,'pCLA',3x,'Tbot',5x,'Tm',5x,'etam',8x,'k',            &
                      3x,'Fcond_bot',1x,'Fcond_top',3x,'Fbot',5x,'Ftop',4x,'dlid',1x,'z_topTBL',1x,'z_botTBL')
1110  format(2(1x,f7.2),1x,d11.5,2(1x,i2),1x,f8.5,2(1x,f7.2),1x,d11.5,1x,f6.3,4(1x,f8.3),3(1x,f7.2))
1200  format('R_tot =',f7.1,' R_core = ',f7.1,'  Tsurf = ',f8.2,'  g = ',f6.3)
1205  format(2x,'time_Gyr',3x,'Prad',4x,'Ptidal',4x,'Fcore',5x,'Pint',5x,'Hint',5x,'conv melt',1x,'P(MPa)',3x,'zbot',     &
                                               5x,'%vol',4x,'Tbot',5x,'Tm',6x,'Fbot',5x,'Ftop',5x,'dlid',3x,'T_core')
1210  format(1x,d11.5,2(1x,f8.4),1x,f8.3,1x,f8.4,1x,d11.5,2(1x,i2),1x,f7.2,1x,f8.3,1x,f8.5,2(1x,f7.2),2(1x,f8.3),         &
                                                                                                 1x,f7.2,1x,f8.2)
1500  format('R_tot =',f7.1,' R_core = ',f7.1,'  Tsurf = ',f8.2,'  rho_core = ',f8.2,'  k_core = ',f7.2)
1505  format(1x,'Radius (km) and temperature (K) at ',*(1x,f7.3),' Gyr')
1510  format(*(1x,f8.2,1x,f8.2))

      stop

      end program thermal_evolution

!-----------------------------------------------------------------------------------------
!  Time derivative of ice shell inner radius                                             |
!-----------------------------------------------------------------------------------------
      subroutine derivative(n_int,time,rx,drxdt)

      use thermal_evolution_mod

      implicit real(8) (a-h,o-z)

      integer, parameter :: nimax=10
      logical :: convecting
      real(8) :: km,Lsolid
      real(8) :: rx(nimax),drxdt(nimax)
      real(8) :: T_core_aux(npcmax)

      r_bot = rx(1)
      z_bot = R_tot - r_bot
      dtime     = time - time_save
      dtime_sec = dtime*year_sec

!---> Temperature at bottom of iceI layer, and temperature jump
      call Temperature_bottom_iceI_shell(z_bot,pres,pCLA,Tbot,dTbotdr,DeltaT)

!---> Heat flux at the top of the core
      call Radiogenic_heat_core(time,H_core,F_radio_KS,DT_core)

!---> Radiogenic and tidal heat in the core in W/m^3, Hrad and Htidal
      Hrad   = H_core*rho_core
      Htidal = Hint_core*rho_core

!---> Radial profile of temperature in the core
      Tcs               = Tbot
      T_core_aux(1:npc) = T_core(1:npc)
      call compute_Tcore_profile(npc,dr_core,dtime_sec,Tcs,Hrad,Htidal,k_core,rho_core,Cp_core,r0_core,T_core_aux(1:npc))

!---> Heat flux at the top of the core :
!    - Radiogenic heat flux calculated from Kirk and Stevenson (1987)
!    - Heat flux calculated from the temperature profile calculated by solving heat diffusion equation
      if (core_heatflux_mode == 1) then
         F_core = F_radio_KS + Htidal*R_core/3.0d0
      else
         F_core = k_core*(-T_core(npc-2) + 4.0d0*T_core(npc-1) - 3.0d0*T_core(npc))/2.0d0/dr_core
      endif

!---> Derivative of adiabatic temperature increase in ocean
      dTadiadr = -Tbot*alfa_H2O*rhog/rho_H2O/Cp_H2O

!---> Heat flux at the bottom of ice I layer
      call HeatFlux_bottom_iceI_shell(convecting,z_bot,Tbot,DeltaT,Tm,etam,km,alfam,Ra,Hadim,Hcrit,     &
                                                            Fcond,Fbot,Ftop,dlid,z_topTBL,z_botTBL)

      dPower     = r_bot**2.0d0*Fbot - R_core**2.0d0*F_core 
      cool_ocean = rho_H2O*Cp_H2O*(dTbotdr - dTadiadr)*(r_bot**3.0d0 - R_core**3.0d0)/3.0d0
      Lsolid     = rho_iceI*L_iceI*r_bot**2.0d0

      drxdt(1)   = dPower/(cool_ocean - Lsolid)*year_sec

      return

      end subroutine derivative

!-----------------------------------------------------------------------------------------
!  Radiogenic heat and heat flux at the top of the core from Kirk & Stevenson (1987)     | 
!-----------------------------------------------------------------------------------------
      subroutine Radiogenic_heat_core(time,Hrad,flux,DTemp)

!---> Power is calculated from heat production rate, present-day concentration, and 
!     desintegration constant for each element. Four elements considered:
!     1 = U235, 2 = U238, 3 = Th232, 4 = K40

      use thermal_evolution_mod

      implicit real(8) (a-h,o-z)

      integer, parameter :: nele=4                   ! Number of elements
      real(8) :: H(nele),lambda(nele)                ! Heat release and decay cosntants
      real(8) :: C0(nele),C(nele)                    ! Initial and comptemporaneous concentrations 

      time_year = time
      time_sec  = time*year_sec

!---> Heat release in W/kg (Turcotte & Schubert (2001)
      H(1) = 5.687d-4
      H(2) = 9.465d-5
      H(3) = 2.638d-5
      H(4) = 2.917d-5

!---> Initial concentrations 
      select case (chondritic_model)

!  From chondritic meteorites (Lodders, 2003)
      case(1)
         C(1) =   5.4d-9
         C(2) =  19.9d-9
         C(3) =  38.7d-9
         C(4) = 738.0d-9

!  From CI meteorites (Lodders & Fegley, 1998)
      case(2)
         C(1) =   5.09d-9
         C(2) =  16.09d-9
         C(3) =  36.31d-9
         C(4) = 760.54d-9

!  From CM meteorites (Lodders & Fegley, 1998)
      case(3)
         C(1) =   7.63d-9
         C(2) =  24.31d-9
         C(3) =  51.34d-9
         C(4) = 511.64d-9

!  From LL meteorites (Lodders & Fegley, 1998)
      case(4)
         C(1) =    9.54d-9
         C(2) =   30.16d-9
         C(3) =   58.85d-9
         C(4) = 1216.90d-9

      case(:0, 5:)
         write(*,*) 'Selected initial concentration model not valid. Progam stops . . . '
         stop
      end select

!---> Decay constants (in yr^-1)
      lambda(1) = 9.8485d-10
      lambda(2) = 1.5514d-10
      lambda(3) = 4.9405d-11
      lambda(4) = 5.4279d-10

      timeref   = 0.0d0
      do 5 i=1,nele
         C0(i) = C(i)*dexp(lambda(i)*timeref)
5     continue

!---> Total heat per unit mass
      Hsum = 0.0d0
      Hrad = 0.0d0
      do 10 i=1,nele
         if (time_year == 0.0d0) time_year = 1.0d0
         a    = lambda(i)*time_year
         Hrad = Hrad + H(i)*C0(i)*dexp(-a)
         Hsum = Hsum + H(i)*C0(i)*(1.0d0 - dexp(-a))/a
10    continue

!---> Surface core flux and temperature acording to Kirk & Stevenson (1987)
      flux  = 2.0d0*dsqrt(kappa_core*time_sec/pi)*rho_core*Hsum
      DTemp = 0.5d0*dsqrt(kappa_core*time_sec*pi)*flux/k_core

      return

      end subroutine Radiogenic_heat_core

!-----------------------------------------------------------------------------------------
!  Time-dependent tidal heating                                                          |
!-----------------------------------------------------------------------------------------
      subroutine Time_dependent_tidal_heating(heat_timevar,time,period,a0,a1,a_expo,H0,H)

      use thermal_evolution_mod

      implicit real(8) (a-h,o-z)

      integer :: heat_timevar

      time_Gyr = time/1.0d9

      select case (heat_timevar)

!-    Heating constant with time
      case(1)
         H = H0

!-    Sinusoidal
      case(2)
         frq   = 2.0d0*pi/period
         vsin  = a0 + a1*dsin(time_Gyr*frq)
!!!         vsin  = (a0 + a1*dsin(time_Gyr*frq))/(a0 + a1)
         H     = H0*vsin

!-    Sinusoidal + regular exponential decrease
      case(3)
         frq   = 2.0d0*pi/period
         vsin  = a0 + a1*dsin(time_Gyr*frq)
!!!         vsin  = (a0 + a1*dsin(time_Gyr*frq))/(a0 + a1)
         expo  = dexp(-a_expo*time)
         H     = H0*expo*vsin   

!-    Sinusoidal^2
      case(4)
         frq   = 2.0d0*pi/period
         vsin  = a0 + a1*dsin(time_Gyr*frq)
!!!         vsin  = (a0 + a1*dsin(time_Gyr*frq))/(a0 + a1)
         H     = H0*vsin**2.0d0   

!-    Sinusoidal^2 + regular exponential decrease
      case(5)
         frq   = 2.0d0*pi/period
         vsin  = a0 + a1*dsin(time_Gyr*frq)
!!!         vsin  = (a0 + a1*dsin(time_Gyr*frq))/(a0 + a1)
         expo  = dexp(-a_expo*time)
         H     = H0*expo*vsin**2.0d0

!-    Undifined evolution
      case(:0, 6:)
         write(*,*) 'Undefined tidal heating time-dependence ; heating set to constant'
         H = H0
      end select

      end subroutine Time_dependent_tidal_heating

!-----------------------------------------------------------------------------------------
!  Temperature at bottom of iceI layer, and temperature jump                             |
!-----------------------------------------------------------------------------------------
      subroutine Temperature_bottom_iceI_shell(z_bot,pres,pCLA,Tbot,dTbotdr,DeltaT)

!---> Calculates temperature at bottom of iceI layer and derivative, and temperature jump 
!     from phase diagrams of water and volatiles

      use thermal_evolution_mod

      implicit real(8) (a-h,o-z)

      rhog = rho_iceI*g

      if (composition == 1) then
         call T_fusion_H2O(z_bot,pres,Tbot,dTbotdr)
         pCLA = 0.0d0
      else

         r_bot    = R_tot - z_bot
         Volume   = 4.0d0*pi*(r_bot**3.0d0 - R_core**3.0d0)/3.0d0
         pCLA_vol = pCLA_init*Volume_init/Volume
         if (pCLA_vol > 1.0d0) pCLA_vol = 1.0d0
         pCLA     = pCLA_vol*rho_NH3/(pCLA_vol*rho_NH3 + (1.0d0 - pCLA_vol)*rho_H2O)

         if (pCLA > pCLAmax) pCLA = pCLAmax
         call T_fusion_mixture(z_bot,pres,pCLA,TH2O,TNH3_15pc,TNH3_eutec,TCH3OH,Tbot)

         dr = 5.0d0
         if (z_bot == 0.0d0) then
            zm = z_bot
            zp = z_bot + dr
         else if (z_bot == (R_tot - R_core)) then
            zm = z_bot - dr
            zp = z_bot
         else
            zm = z_bot - dr
            zp = z_bot + dr
         endif         
         call T_fusion_mixture(zm,pres_zm,pCLA,TH2O,TNH3_15pc,TNH3_eutec,TCH3OH,T_zm)
         call T_fusion_mixture(zp,pres_zp,pCLA,TH2O,TNH3_15pc,TNH3_eutec,TCH3OH,T_zp)
         dTbotdr = -(T_zp - T_zm)/(zp - zm)

      endif
      DeltaT = Tbot - Tsurf

      return

      end subroutine Temperature_bottom_iceI_shell

!-----------------------------------------------------------------------------------------
!  Heat flux at bottom of iceI layer, and temperature jump                               |
!-----------------------------------------------------------------------------------------
      subroutine Heatflux_bottom_iceI_shell(convecting,z_bot,Tbot,DeltaT,Tm,etam,km,alfam,Ra,Hadim,Hcrit,     &
                                                                  Fcond,Fbot,Ftop,dlid,z_topTBL,z_botTBL)

!---> Calculates heat flux at bottom of iceI layer from appropriate scaling laws
!     of thermal convection or conductive profiles

      use thermal_evolution_mod

      implicit real(8) (a-h,o-z)

      logical :: convecting
      real(8) :: km,km_cond,kref,kappam

!---> Set ratio of inner-to-outer radii, f
      r_bot = R_tot - z_bot
      f     = r_bot/R_tot

!---> Set ratio of top-to-bottom thermal conducivity, rk, if applicable 
      if (variable_conductivity) then
         rk = Tsurf/Tbot
         rT = DeltaT/Tsurf
      else
         rk = 1.0d0
         rT = 0.0d0
      endif

!---> Average temperature for conductive layer, and viscosity and conductivity
!     at this temperature
      if (variable_conductivity) then
         Tm_cond = Tsurf*dexp(0.5d0*dlog(Tbot/Tsurf))
         km_cond = k_iceI*Tsurf/Tm_cond
      else

         if (geometry == 1) then
            Tm_cond = Tsurf + 0.5d0*(Tbot - Tsurf)
         else
            Tm_cond = Tsurf + (Tbot - Tsurf)/(1.0d0 + f)
         endif

         if (reference_conductivity == 1) then
             km_cond = k_iceI
         else if (reference_conductivity == 2) then
             if (Hint_ice > 0.0d0) then
                cte_k   = 1.0d0 - rho_iceI*Hint_ice*z_bot**2.0d0/8.0d0/k_iceI
                km_cond = cte_k*k_iceI*Tsurf/Tm_cond
             else
                km_cond = k_iceI*Tsurf/Tm_cond
             endif
         else if (reference_conductivity == 3) then
             km_cond = k_input
         endif

         if (Hint_ice > 0.0d0) then
            Tm_cond_H = rho_iceI*Hint_ice*z_bot**2.0d0/8.0d0/km_cond
            Tm_cond = Tm_cond + Tm_cond_H
         endif
         
      endif
      etam_cond = etaref*dexp(Eact*(TH2O/Tm_cond - 1.0d0)/R/TH2O)

!---> Temperature in the well-mixed interior, Tm, from stagnant-lid scaling
      if (scaling_law == 1) then
         Tm = 0.5d0*(Tsurf + Tbot)
      else
         a = c1*R/Eact/rk**c2
         if (geometry == 2) a = a/f**c3
         Tm = 0.5*(dsqrt(1.0d0 + 4.0d0*a*Tbot) - 1.0d0)/a
      endif

!---> Viscosity at temperature of well-mixed interior (if applicable)
      if (scaling_law == 1) then
         etam = etaref
      else
         call T_fusion_H2O(z_bot,pres_dum,TH2O,dTbotdr_dum)     
         etam = etaref*dexp(Eact*(TH2O/Tm - 1.0d0)/R/TH2O)
      endif

!---> Thermal conductivity at temperature of well-mixed interior.
!     Surface value of k is stored in k_iceI
      km     = k_iceI*Tsurf/Tm
      kappam = kappa_iceI*Tsurf/Tm

!---> Calculates effective Rayleigh number, at values of k, eta, T, etc ...
!     determined above. 
      alfam = alfa_iceI*Tm/250.0d0
      if (variable_conductivity) then
         Ra = alfam*rho_iceI*g*DeltaT*z_bot**3.0d0/etam/kappam
      else
         if (reference_conductivity == 1) then
            Ra = alfam*rho_iceI*g*DeltaT*z_bot**3.0d0/etam/kappa_iceI
         else if (reference_conductivity == 2) then
            Ra = alfam*rho_iceI*g*DeltaT*z_bot**3.0d0/etam/kappam
         else if (reference_conductivity == 3) then
            Ra = alfam*rho_iceI*g*DeltaT*z_bot**3.0d0/etam/kappa_input
!!!            Ra = alfam*rho_iceI*g*DeltaT*z_bot**3.0d0/etam/kappam
         endif
      endif
      Ra_etaref = Ra*etam/etaref

!---> Internal heat case: re-assess Tm, etam, km, and Ra taking with a scaling
!     law for temperature that takes into account internal heating.
!     scaling law for temperature
      if (internal_heating) then
         c4 = c4pos
         c5 = c5pos
         c6 = c6pos
         c7 = c7pos
         call newton_search(z_bot,DeltaT,rk,f,TH2O,Tm,etam,km,Ra,Hadim)
      else
         Hadim = 0.0d0
      endif

!---> Defines reference conductivity to rescale heat flux
      if (variable_conductivity) then
         kref = k_iceI
      else
         if (reference_conductivity == 1) then
             kref = k_iceI
         else if (reference_conductivity == 2) then
             kref = km
         else if (reference_conductivity == 3) then
             kref = k_input
         endif
      endif

!---> Characteristic and conductive heat fluxes. kref is reference thermal conductivity,
!     defined above.
      Fcarac = kref*DeltaT/z_bot
      
      Fcond        = 0.0d0
      Fcond_surf_H = 0.0d0
      Fcond_bot_H  = 0.0d0
      if (z_bot /= 0.0d0) then
         if (variable_conductivity) then
            Fcond_adim = dlog(1.0d0 + rT)/rT
            Fcond      = kref*Tsurf*dlog(Tbot/Tsurf)/z_bot
         else
!!!            Fcond_adim = 1.0d0
            Fcond_adim = f + (f + 2.0d0)*Hadim/6.0d0
            Fcond_surf_H = rho_iceI*z_bot*Hint_ice*(2.0d0 + f)/6.0d0
            Fcond_bot_H  = rho_iceI*z_bot*Hint_ice*(1.0d0 + 2.0d0*f)/f/6.0d0
            Fcond        = kref*DeltaT/z_bot
         endif
      endif
      Fcond_surf = Fcond*f + Fcond_surf_H
      Fcond_bot  = Fcond/f - Fcond_bot_H

!---> Convective heat flux at bottom and top of the shell (Fbot and Ftop) and Racrit
!     from Stengel et al. (1982) (see also Deschamps and Lin, 2014)
      if (scaling_law == 1) then
         Fbot_adim = ac*Ra**betac/f
         Ftop_adim = Fbot_adim*f**2.0d0
         Racrit    = 657.0d0
      else

         gama = Eact*DeltaT/R/Tm**2.0d0

         aparam    = aparam_pos
         Ftop_adim = aparam*Ra**beta*rk**ck/gama**cgam
         if (geometry == 2) then
            Ftop_adim = Ftop_adim*f**cf
         else
            Ftop_adim = Ftop_adim*f
         endif
         Fbot_adim = Ftop_adim/f**2.0d0

!---> Old formulation (GJI 2021)
!!!         Fbot_adim = aparam*Ra**beta*rk**ck/gama**cgam
!!!         if (geometry == 2) then
!!!            Fbot_adim = Fbot_adim/f**cf
!!!         else
!!!            Fbot_adim = Fbot_adim/f
!!!         endif
!!!         Ftop_adim = Fbot_adim*f**2.0d0

         if (internal_heating) then
            Fbot_adim = Fbot_adim - (1.0d0 + f + f**2.0d0)*Hadim/3.0d0/f**2.0d0

!---> Switches to scaling laws for Ur > 1 if bottom heat flux predicted by scaling for
!     Ur < 1 is lower than 0.
            if (Fbot_adim < 0.0d0) then
               c4     = c4neg
               c5     = c5neg
               c6     = c6neg
               c7     = c7neg
               aparam = aparam_neg

               call newton_search(z_bot,DeltaT,rk,f,TH2O,Tm,etam,km,Ra,Hadim)

               Ftop_adim = aparam*Ra**beta*rk**ck/gama**cgam
               if (geometry == 2) then
                  Ftop_adim = Ftop_adim*f**cf
               else
                  Ftop_adim = Ftop_adim*f
               endif
               Fbot_adim = Ftop_adim/f**2.0d0 - (1.0d0 + f + f**2.0d0)*Hadim/3.0d0/f**2.0d0

!     Sets bottom heat flux to zero if bottom heat flux predicted by scaling for
!     Ur > 1 is again larger than 0.
               if (Fbot_adim > 0.0d0) then

                  Fbot_adim = 0.0d0
                  Ftop_adim = (1.0d0 + f + f**2.0d0)*Hadim/3.0d0

                  aparam = 0.5d0*(aparam_pos + aparam_neg)
                  cte1 = R*Tbot/Eact/beta
                  cte2 = Ftop_adim/aparam/Ra_etaref**beta

                  aux = cte1*dlog(cte2*gama**cgam)
                  fx  = Tm*(1.0d0 - aux) - Tbot

                  nite = 0
                  do while ((dabs(fx) > eps).and.(nite <= nitemax))

                     nite = nite+1
                     dTm  = Tm*1.0d-3

                     Tmm   = Tm - 0.5d0*Tm
                     gamam = Eact*DeltaT/R/Tmm**2.0d0
                     auxm  = cte1*dlog(cte2*gamam**cgam)
                     fxm   = Tmm*(1.0d0 - auxm) - Tbot

                     Tmp   = Tm + 0.5d0*Tm
                     gamap = Eact*DeltaT/R/Tmp**2.0d0
                     auxp  = cte1*dlog(cte2*gamap**cgam)
                     fxp   = Tmp*(1.0d0 - auxp) - Tbot

                     a_fx = (fxp - fxm)/dTm
                     b_fx = fx - a_fx*Tm
                     Tm   = -b_fx/a_fx 

                     gama = Eact*DeltaT/R/Tm**2.0d0
                     aux  = cte1*dlog(cte2*gama**cgam)
                     fx   = Tm*(1.0d0 - aux) - Tbot

                  enddo

                  etam   = etaref*dexp(Eact*(TH2O/Tm - 1.0d0)/R/TH2O)
                  kappam = kappa_iceI*Tsurf/Tm
                  alfam  = alfa_iceI*Tm/250.0d0
                  Ra     = alfam*rho_iceI*g*DeltaT*z_bot**3.0d0/etam/kappam

               endif

            endif

            Hcrit = 0.184d0*dexp(-0.19d0*gama)*Ra_etaref**0.31d0

         endif

         Racrit = 6.1d0*gama**4.0d0
      endif

      convecting = .true.
      if (Ftop_adim >= Fcond_adim) then
!!!      if (Ra >= Racrit) then
         Fbot = Fbot_adim*Fcarac
         Ftop = Ftop_adim*Fcarac
      else
         convecting = .false.
         Fbot  = Fcond_bot
         Ftop  = Fcond_surf
!!!         Tm    = Tm_cond
!!!         etam  = etam_cond
!!!         km    = km_cond
      endif

!---> Conductive lid
      if (.not.convecting) then
         Tlid_adim = 1.0d0
         Tlid      = Tbot
         dlid_adim = 1.0d0
         dlid      = z_bot
      else   
         Tm_adim   = (Tm - Tsurf)/DeltaT
         if (Tm < Tbot) then
            Tlid_adim = 2.0d0*Tm_adim - 1.0d0
            Tlid      = 2.0d0*Tm - Tbot
         else
            Tlid_adim = Tm_adim
            Tlid      = Tm
         endif

         if (variable_conductivity) then
!-> The two lines below (used in Deschamps, GJI, 2021), is an approximation as it identifies
!   the curvature of the lid, flid =  (R - dlid)/R = 1 - (1 -f)*dlid/D, with the curvature
!   of the whole lid. It is now replaced with the next lines.
!!!            dlid_adim = f*dlog(1.0d0 + rT*Tlid_adim)/Ftop_adim/rT
!!!            dlid      = kref*Tsurf*f*dlog(Tlid/Tsurf)/Ftop
            aux1_dlid = dlog(1.0d0 + rT*Tlid_adim)/rT
            dlid_adim = aux1_dlid/(Ftop_adim + (1.0d0 - f)*aux1_dlid)
            aux2_dlid = kref*Tsurf*dlog(Tlid/Tsurf)
            dlid      = aux2_dlid/(Ftop + (1.0d0 - f)*aux2_dlid/z_bot)
         else
            dlid_adim = f*Tlid_adim/Ftop_adim
            if (internal_heating) then
!->  The two lines below assume Tlid = 2*Tm -1.0 (not correct for a mix-heated fluid and
!    is strictly speaking valid for f = 1 (no curvature) only.
!!!               aux_dlid  = 1.0d0 - 2.0d0*Hadim*Tlid_adim/Ftop_adim**2.0d0
!!!               dlid_adim = Ftop_adim*(1.0d0 - dsqrt(aux_dlid))/Hadim
!->  Calculates dlid from appropriate scaling
               aux_dlid = f*gama**1.21d0/Ra**0.27d0/rk**0.82d0
               if (Ftop <= 0.0d0) then
                  dlid_adim = 0.63d0*aux_dlid
               else
                  dlid_adim = 0.67d0*aux_dlid
               endif
               dlid = dlid_adim*z_bot
            else
!-> The line below (used in Deschamps, GJI, 2021), is an approximation as it identifies the
!   curvature of the lid, flid =  (R - dlid)/R = 1 - (1 -f)*dlid/D, with the curvature of
!   the whole lid. It is now replaced with the next line.
!!!               dlid_adim = f*Tlid_adim/Ftop_adim
               dlid_adim = Tlid_adim/(Ftop_adim + (1.0d0 - f)*Tlid_adim)
               dlid      = dlid_adim*z_bot
!!!            dlid      = kref*f*(Tlid - Tsurf)/Ftop
            endif
         endif
      endif
      if (dlid_adim > 1.0d0) dlid_adim = 1.0d0
      if (dlid_adim < 0.0d0) dlid_adim = 0.0d0
      if (dlid > z_bot) dlid = z_bot
      if (dlid < 0.0d0) dlid = 0.0d0

!---> Top and bottom TBL thickness
      dT_topTBL = Tm - Tsurf
      if (.not.convecting) then
         z_topTBL = z_bot
      else
         Hint_vol = Hint_ice*rho_iceI
         if (f == 1.0d0) then
            delta = 1.0d0 - 2.0d0*kref*Hint_vol*dT_topTBL/Ftop**2.0d0
            if (delta >= 0.0d0) then
               z_topTBL = Ftop*(1.0d0 - dsqrt(delta))/Hint_ice*rho_iceI
            else
               z_topTBL = Ftop/Hint_ice*rho_iceI
            endif
         else
            a0_pol = -dT_topTBL
            a1_pol = Ftop/kref + (1.0d0 - f)*dT_topTBL/z_bot
            a2_pol = -Hint_vol/2.0d0/kref
            a3_pol = (1.0d0 - f)*Hint_vol/6.0d0/kref/z_bot

            b0_pol = a1_pol
            b1_pol = 2.0d0*a2_pol
            b2_pol = 3.0d0*a3_pol
            
            z_topTBL = dlid
            gx       = a3_pol*z_topTBL**3.0d0 + a2_pol*z_topTBL**2.0d0 + a1_pol*z_topTBL + a0_pol

            nite = 0
            do while ((dabs(gx) > eps).and.(nite <= nitemax))

               nite = nite+1

               a_gx     = b2_pol*z_topTBL**2.0d0 + b1_pol*z_topTBL + b0_pol
               b_gx     = gx - a_gx*z_topTBL
               z_topTBL = -b_gx/a_gx

               gx       = a3_pol*z_topTBL**3.0d0 + a2_pol*z_topTBL**2.0d0 + a1_pol*z_topTBL + a0_pol

            end do

         endif
      endif

      dT_botTBL = Tbot - Tm
      Fbase_eps = 1.0d-5
      Fbase     = Fbot
      if ((dT_botTBL <= 0.0d0).and.(Fbase >= 0.0d0)) Fbase = 0.0d0
      if ((.not.convecting).or.(dabs(Fbase) <= Fbase_eps)) then
         z_botTBL = z_bot
      else
         if (.not.variable_conductivity) then
            tcond_botTBL = kref
            if (f == 1.d0) then
               d_botTBL = tcond_botTBL*dT_botTBL/Fbase
               z_botTBL = z_bot - d_botTBL
            else
               Flin     = tcond_botTBL*dT_botTBL/r_bot
               r_botTBL = r_bot*Fbase/(Fbase - Flin)
               z_botTBL = R_tot - r_botTBL
               if (z_botTBL >= z_bot) z_botTBL = z_bot
            endif
         else
            tcond_botTBL = k_iceI_a/Tm + k_iceI_b
            if (f == 1.0d0) then
               d_botTBL = tcond_botTBL*Tm*dlog(Tbot/Tm)/Fbase
               z_botTBL = z_bot - d_botTBL
            else
               Flin     = tcond_botTBL*dT_botTBL/r_bot
               r_botTBL = r_bot*Fbase/(Fbase - Flin)
               z_botTBL = R_tot - r_botTBL
               if (z_botTBL >= z_bot) z_botTBL = z_bot
            endif
         endif
      endif

      return

      end subroutine Heatflux_bottom_iceI_shell

!-----------------------------------------------------------------------------------------
!  Solve for internal temperature in case of system with internal heating following      |
!  appropriate parameterization                                                          |
!-----------------------------------------------------------------------------------------
      subroutine newton_search(z_bot,DeltaT,rk,f,TH2O,Tm,etam,km,Ra,Hadim)

      use thermal_evolution_mod

      implicit real(8) (a-h,o-z)

      real(8) :: km,kref,kappam

      nitemax = 100
      eps     = 1.0d-5

      Tbot = Tsurf + DeltaT
      if (geometry == 2) then
         cste1 = c1*R/Eact/rk**c2/f**c3
         ageom = (1.0d0 + f + f**2.0d0)/3.0d0
         cste2 = (c4 + c5*f)*ageom**c6
      else
         cste1 = c1*R/Eact/rk**c2
         ageom = 1.0d0
         cste2 = c4 + c5
      endif

!---> NB: Internal heating should be non-dimensionalized usid the characteristic heat 
!     flux, which is defined with the characteristic thermal conductivity. In the case
!     of temperature-dependent conductivity, kref = ksurf.
     if (variable_conductivity) then
         kref = k_iceI
      else
         if (reference_conductivity == 1) then
            kref = k_iceI
         else if (reference_conductivity == 2) then
            kref = km
         else if (reference_conductivity == 3) then
            kref = k_input
         endif
      endif
      Hadim = rho_iceI*Hint_ice*z_bot**2.0d0/kref/DeltaT
      fx    = Tm - Tbot + cste1*Tm**2.0d0 - cste2*DeltaT*Hadim**c6/Ra**c7

      nite = 0
      do while ((dabs(fx) > eps).and.(nite <= nitemax))

         nite = nite+1
         dTm  = Tm*1.0d-3

         Tmm = Tm - 0.5d0*dTm
         if ((variable_conductivity).and.(reference_conductivity == 2)) then
            km   = k_iceI*Tsurf/Tmm
            kref = km
         endif
         kappam = kappa_iceI*Tsurf/Tmm
         alfam  = alfa_iceI*Tmm/250.0d0
         etam   = etaref*dexp(Eact*(TH2O/Tmm - 1.0d0)/R/TH2O)
         Ra     = alfam*rho_iceI*g*DeltaT*z_bot**3.0d0/etam/kappam
         Hadim  = rho_iceI*Hint_ice*z_bot**2.0d0/kref/DeltaT
         fxm    = Tmm - Tbot + cste1*Tmm**2.0d0 - cste2*DeltaT*Hadim**c6/Ra**c7

         Tmp = Tm + 0.5d0*dTm
         if ((variable_conductivity).and.(reference_conductivity == 2)) then
            km   = k_iceI*Tsurf/Tmp
            kref = km
         endif
         kappam = kappa_iceI*Tsurf/Tmp
         alfam  = alfa_iceI*Tmp/250.0d0
         etam   = etaref*dexp(Eact*(TH2O/Tmp - 1.0d0)/R/TH2O)
         Ra     = alfam*rho_iceI*g*DeltaT*z_bot**3.0d0/etam/kappam
         Hadim  = rho_iceI*Hint_ice*z_bot**2.0d0/kref/DeltaT
         fxp    = Tmp - Tbot + cste1*Tmp**2.0d0 - cste2*DeltaT*Hadim**c6/Ra**c7

         a_fx = (fxp - fxm)/dTm
         b_fx = fx - a_fx*Tm
         Tm   = -b_fx/a_fx

         if ((variable_conductivity).and.(reference_conductivity == 2)) then
            km   = k_iceI*Tsurf/Tm
            kref = km
         endif
         etam   = etaref*dexp(Eact*(TH2O/Tm - 1.0d0)/R/TH2O)
         kappam = kappa_iceI*Tsurf/Tm
         alfam  = alfa_iceI*Tm/250.0d0
         Ra     = alfam*rho_iceI*g*DeltaT*z_bot**3.0d0/etam/kappam
         Hadim  = rho_iceI*Hint_ice*z_bot**2.0d0/kref/DeltaT

         fx = Tm - Tbot + cste1*Tm**2.0d0 - cste2*DeltaT*Hadim**c6/Ra**c7

      end do

      return

      end subroutine newton_search

!-----------------------------------------------------------------------------------------
!  Integrate to next point according to chosen method                                    |
!-----------------------------------------------------------------------------------------
      subroutine integrate(iint,n_int,eps,time,dtime,rx,drxdt)

      use thermal_evolution_mod

      implicit real(8) (a-h,o-z)

      integer, parameter :: nimax=10
      real(8) :: rx(nimax),drxdt(nimax),rscal(nimax),rout(nimax)

!---> 4th-order Runge-Kutta
      if (iint == 1) then
         call rk4(rx,drxdt,n_int,time,dtime,rout)
         do 10 j=1,n_int
            rx(j) = rout(j)
10       continue
      endif

!---> 4th-order Runge-Kutta with adaptative step
      if (iint == 2) then
         do 20 j=1,n_int
            rscal(j) = dabs(rx(j)) + dtime*dabs(drxdt(j))
20       continue
         call rkqc(rx,drxdt,n_int,time,dtime,eps,rscal,hdid,hnext)
      endif

!---> Bulirsh-Stoer
      if (iint == 3) then
         do 30 j=1,n_int
            rscal(j) = dabs(rx(j)) + dtime*dabs(drxdt(j))
30       continue
         call bsstep(rx,drxdt,n_int,time,dtime,eps,rscal,hdid,hnext)
      endif

      return

      end subroutine integrate

!-----------------------------------------------------------------------------------------
!  Quality controle Runge-Kutta, d'apres Press et al., 1986                              |
!  call subroutine rk4 (4th order Runge Kutta                                            |
!-----------------------------------------------------------------------------------------
      subroutine rkqc(y,dydx,n,x,htry,eps,yscal,hdid,hnext)

      implicit real(8) (a-h,o-z)

      integer, parameter :: nimax=10
      real(8) :: y(n),dydx(n),yscal(n),ytemp(nimax),ysav(nimax),dysav(nimax)

      pgrow  = -0.2d0
      pshrnk = -0.25d0
      fcor   = 1.0d0/1.5d1
      one    = 1.0d0
      safety = 0.9d0
      errcon = 6.0d-4

      xsav = x
      do 5 i=1,n
         ysav(i)  = y(i)
         dysav(i) = dydx(i)
5     continue

      h  = htry
10    hh = 0.5d0*h
      call rk4(ysav,dysav,n,xsav,hh,ytemp)
      x = xsav + hh
      call derivative(n,x,ytemp,dydx)
      call rk4(ytemp,dydx,n,x,hh,y)
      x = xsav + h
      if (x == xsav) then 
         write(*,*) 'Stepsize not significant in RKQC'
         stop
      endif 
      call rk4(ysav,dysav,n,xsav,h,ytemp)
      errmax = 0.0d0
      do 20 i=1,n
         ytemp(i) = y(i) - ytemp(i)
         drscal   = dabs(ytemp(i)/yscal(i))
         errmax   = max(errmax,drscal)
20    continue
      errmax = errmax/eps
      if (errmax > one) then
         h = safety*h*errmax**pshrnk
         goto 10
      else
         hdid = h
         if (errmax > errcon) then
            hnext = safety*h*errmax**pgrow
         else
            hnext = 4.0d0*h
         endif
      endif

      do 30 i=1,n
         y(i) = y(i) + ytemp(i)*fcor
30    continue

      return

      end subroutine rkqc

!-----------------------------------------------------------------------------------------
!   Runge-Kunta ordre 4, d'apres Press et al. (1986)                                     |
!-----------------------------------------------------------------------------------------
      subroutine rk4(y,dydx,n,x,h,yout)

      implicit real(8) (a-h,o-z)

      integer, parameter :: nimax=10
      real(8) y(n),dydx(n),yout(n),yt(nimax),dyt(nimax),dym(nimax)

      hh = h*0.5d0
      h6 = h/6.0d0
      xh = x + hh
      do 10 i=1,n
         yt(i) = y(i) + hh*dydx(i)
10    continue
      call derivative(n,xh,yt,dyt)
      do 20 i=1,n
         yt(i) = y(i) + hh*dyt(i)
20    continue
      call derivative(n,xh,yt,dym)
      do 30 i=1,n
         yt(i)  = y(i) + h*dym(i)
         dym(i) = dyt(i) + dym(i)
30    continue
      call derivative(n,x+h,yt,dyt)
      do 40 i=1,n
         yout(i) = y(i) + h6*(dydx(i) + 2.0d0*dym(i) + dyt(i))
40    continue

      return

      end subroutine rk4

!-----------------------------------------------------------------------------------------
!   Bulirsch-Stoer method, after Press et al. (1986) (pp. 563-568)                       |
!-----------------------------------------------------------------------------------------
      subroutine bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext)

      implicit real(8) (a-h,o-z)

      integer, parameter :: nimax=10,imax=11,nuse=7
      integer :: nseq(imax)
      real(8) :: y(nv),dydx(nv),yscal(nv)
      real(8) ::  yerr(nimax),ysav(nimax),dysav(nimax),yseq(nimax)

      data nseq/2,4,6,8,12,16,24,32,48,64,96/
      data one,shrink,grow/1.0d0,0.95d0,1.2d0/

      h    = htry
      xsav = x
      do 10 i=1,nv
         ysav(i)  = y(i)
         dysav(i) = dydx(i)
10    continue
20    do 40 i=1,imax
         call mmid(ysav,dysav,nv,xsav,h,nseq(i),yseq)
         xest = (h/dble(nseq(i)))**2.0d0
         call rzextr(i,xest,yseq,y,yerr,nv,nuse)
         if (i > 3) then
            errmax = 0.0d0
            do 30 j=1,nv
               aux    = dabs(yerr(j)/yscal(j))
               errmax = max(errmax,aux)
30          continue
            errmax = errmax/eps
            if (errmax < one) then
               x    = x + h
               hdid = h
               if (i == nuse) then
                  hnext = h*shrink
               else if (i == nuse-1) then
                  hnext = h*grow
               else
                  hnext = h*dble(nseq(nuse-1))/dble(nseq(i))
               endif
               return
            endif
         endif
40    continue

      expo = dble(imax-nuse)/2.0d0
      h    = 0.25d0*h/2.0d0**expo
      if (x+h == x) then
         write(*,*) 'Step size underflow'
         stop
      endif
      goto 20

      end subroutine bsstep

!-----------------------------------------------------------------------------------------
!   Modified Midpoint, after Press et al. (1986) (pp. 560-562)                           |
!-----------------------------------------------------------------------------------------
      subroutine mmid(y,dydx,nv,xs,htot,nstep,yout)

      implicit real(8) (a-h,o-z)

      integer, parameter :: nimax=10
      real*8 :: y(nv),dydx(nv),yout(nv),ym(nimax),yn(nimax)

      h = htot/dble(nstep)
      do 10 i=1,nv
         ym(i) = y(i)
         yn(i) = y(i) + h*dydx(i)
10    continue
      x = xs + h

      call derivative(nv,x,yn,yout)
      h2 = 2.0d0*h
      do 30 n=2,nstep
         do 20 i=1,nv
            aux   = ym(i) + h2*yout(i)
            ym(i) = yn(i)
            yn(i) = aux
20       continue
         x = x + h
         call derivative(nv,x,yn,yout)
30    continue

      do 40 i=1,nv
         yout(i) = 0.5d0*(ym(i) + yn(i) + h*yout(i))
40    continue

      return

      end subroutine mmid

!-----------------------------------------------------------------------------------------
!   Rational extrapolation, after Press et al. (1986) (pp. 566-567)                      |
!-----------------------------------------------------------------------------------------
      subroutine rzextr(iest,xest,yest,yz,dy,nv,nuse)

      implicit real(8) (a-h,o-z)

      integer, parameter :: nimax=10,imax=11,ncol=7
      real(8) :: yest(nv),yz(nv),dy(nv),x(imax),d(nimax,ncol),fx(ncol)

      x(iest) = xest
      if (iest == 1) then
         do 10 j=1,nv
            yz(j)  = yest(j)
            d(j,1) = yest(j)
            dy(j)  = yest(j)
10       continue
      else
         m1 = min(iest,nuse)
         do 20 k=1,m1-1
            fx(k+1) = x(iest-k)/xest
20       continue
         do 40 j=1,nv
            yy     = yest(j)
            v      = d(j,1)
            c      = yy
            d(j,1) = yy
            do 30 k=2,m1
               b1 = fx(k)*v
               b  = b1 - c
               if (b /= 0.0d0) then
                  b   = (c - v)/b
                  ddy = c*b
                  c   = b1*b
               else
                  ddy = v
               endif
               if (k /= m1) v = d(j,k)
               d(j,k) = ddy
               yy     = yy + ddy
30          continue
            dy(j) = ddy
            yz(j) = yy
40       continue
      endif

      return

      end subroutine rzextr

!-----------------------------------------------------------------------------------------
!   Polynomial extrapolation, after Press et al. (1986) (p. 568)                         |
!-----------------------------------------------------------------------------------------
      subroutine pzextr(iest,xest,yest,yz,dy,nv,nuse)

      implicit real(8) (a-h,o-z)

      integer, parameter :: nimax=10,imax=11,ncol=7
      real(8) :: yest(nv),yz(nv),dy(nv),x(imax),qcol(nimax,ncol),d(nimax)

      x(iest) = xest
      do 10 j=1,nv
         dy(j) = yest(j)
         yz(j) = yest(j)
10    continue

      if (iest == 1) then
         do 20 j=1,nv
            qcol(j,1) = yest(j)
20       continue
      else
         m1 = min(iest,nuse)
         do 30 j=1,nv
            d(j) = yest(j)
30       continue
         do 50 k1=1,m1-1
            delta = 1.0d0/(x(iest-k1) - xest)
            f1    = xest*delta
            f2    = x(iest-k1)*delta
            do 40 j=1,nv
               q          = qcol(j,k1)
               qcol(j,k1) = dy(j)
               delta      = d(j)-q
               dy(j)      = f1*delta
               d(j)       = f2*delta
               yz(j)      = yz(j)+dy(j)
40          continue
50       continue
         do 60 j=1,nv
            qcol(j,m1) = dy(j)
60       continue
      endif

      return

      end subroutine pzextr


!-----------------------------------------------------------------------------------------
      subroutine nbrcarac(ncar,chain,nst)

!---> Counts number of character in one chain of characters

      implicit real(8) (a-h,o-z)
      character(len=ncar) :: chain

      nst = 0
      do 5 i=1,ncar
         if (chain(i:i) /= ' ') nst = nst + 1
5     continue

      return

      end subroutine nbrcarac

!-----------------------------------------------------------------------------------------
! Solves heat diffusion equation. By Henri Samuel, Sept. 2020                            |
!-----------------------------------------------------------------------------------------
      SUBROUTINE compute_Tcore_profile(nc,dr,dt,Tcs,Hrad,Htidal,kc,rho,Cp,r0,T)

      IMPLICIT NONE

      INTEGER                  , INTENT(IN)  :: nc
      REAL(8)                  , INTENT(IN)  :: dr,dt,Tcs,Hrad,Htidal,kc,rho,Cp,r0
      REAL(8), DIMENSION(1:nc) , INTENT(OUT) :: T
      REAL(8), DIMENSION(1:nc)               :: cw,cc,ce,rhs,r,k,H
      REAL(8)                                :: ri,dri,dr2i,dti
      INTEGER                                :: ir

      T(nc) = Tcs        
      dri   = 1.0d0/dr
      dr2i  = dri*dri
      dti   = 0.0d0 ; IF (dt > 0.0d0) dti = 1.0d0/dt

!---> Conductivity and internal heating (radioactive & Tidal) profiles ------------------- !
!     NB: Here radioactive and tidal heating are assumed to be constants with r. 
!     However any distribution can be considered and need to be specified below. 
      DO ir=1,nc
         r(ir) = r0 + dr*(ir-1)
         H(ir) = Hrad + Htidal  
         k(ir) = kc
      ENDDO

      DO ir=2,nc-1
         ri      = 1.0d0/r(ir)
         cw( ir) = -k(ir)*dr2i + 0.5d0*dri*((k(ir+1) - k(ir-1))*0.5d0*dri + 2.0d0*k(ir)*ri)
         cc( ir) = rho*Cp*dti  + 2.0d0*k(ir)*dr2i
         ce( ir) = -k(ir)*dr2i - 0.5d0*dri*((k(ir+1) - k(ir-1))*0.5d0*dri + 2.0d0*k(ir)*ri)
         rhs(ir) = rho*Cp*T(ir)*dti + H(ir)
      ENDDO

!---> NB: first-order FD for the Neumann BC at r=r0  for now, Dirichlet at nc, second order FD everywhere else.
      cw( 1) =  0.0d0 ; cw( nc) = 0.0d0
      cc( 1) = +1.0d0 ; cc( nc) = 1.0d0
      ce( 1) = -1.0d0 ; ce( nc) = 0.0d0

      rhs(1) =  0.0d0 ; rhs(nc) = T(nc)

      CALL tridiag(cw,cc,ce,rhs,T(1:nc),nc)

      END SUBROUTINE compute_Tcore_profile

!-----------------------------------------------------------------------------------------
! Tridiagonal matrix resolution (from Numerical recipes                                  |
!-----------------------------------------------------------------------------------------
      SUBROUTINE tridiag(a,b,c,rhs,sol,n)

!---> From Num. recipes
      INTEGER,                 INTENT(IN)  :: n
      REAL(8), DIMENSION(1:n), INTENT(IN)  :: a,b,c,rhs
      REAL(8), DIMENSION(1:n), INTENT(OUT) :: sol
      INTEGER                              :: j
      REAL(8)                              :: bet
      REAL(8), DIMENSION(1:n)              :: gam

      IF (b(1) == 0.0d0) STOP 'Tridag: rewrite equations'
      bet    = b(1)
      sol(1) = rhs(1)/bet

      DO j=2,n
         gam(j) = c(j-1)/bet
         bet    = b(j)-a(j)*gam(j)
         IF (bet == 0.0d0) STOP 'Tridag failed'
         sol(j) = (rhs(j) - a(j)*sol(j-1))/bet
      ENDDO

      DO j=n-1,1,-1
         sol(j) = sol(j) - gam(j+1)*sol(j+1)
      ENDDO

      END SUBROUTINE tridiag

