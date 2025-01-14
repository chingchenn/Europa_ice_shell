!----------------------------------------------------------------------------------------------------|
!  Calculates thermal evolution of icy bodies                                                        |
!                                                                                                    |
!  This version (15-04-2021) by Frederic Deschamps with heat diffusion solver by Henri Samuel        |
!                                                                                                    |
!  Compile:   gfortran -O4 -o thermal_evolution_vs_H-and-NH3 thermal_evolution_vs_H-and-NH3_mod.f95  |
!                                                       thermal_evolution_vs_H-and-NH3.f95           |
!  Execute:   thermal_evolution_vs_H-and-NH3                                                         |
!  Input parameters in file thermal_evolution_vs_H-and-NH3_param                                     |
!----------------------------------------------------------------------------------------------------|
      program thermal_evolution

      use thermal_evolution_mod

      implicit real(8) (a-h,o-z)

      character(len=120) :: casename,filename
      integer, parameter :: nimax=10
      logical :: convecting,stat_thickness
      integer :: heat_param,heat_timevar,heat_scale
      integer :: convec
      real(8) :: km,Lsolid
      real(8) :: rx(nimax),drxdt(nimax),rscal(nimax),rout(nimax)
      real(8) :: time_frame(npwmax)

!---> Reads entry parameters from file thermal_evolution_param
      open(1,file='thermal_evolution_vs_H-and-NH3_param',status='unknown')
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
      write(*,*) 'Reference vistosity (Pas) min, max, and number of points'
      read(1,*) etaref_min,etaref_max,npetaref
      write(*,*) 'Activation energy (kJ/mol)'
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
         write(*,*) 'Internal heating rate (W/kg) in ice layer min, max, and number of points'
         read(1,*) Hint_min,Hint_max,npHint
      else
         write(*,*) 'Total power (TW) min, max, and number of points'
         read(1,*) Pint_min,Pint_max,npHint
         Pint_min = Pint_min*1.0d12
         Pint_max = Pint_max*1.0d12
      endif
      write(*,*) 'Fraction of tidal heating within ice ;layer and within core (sum shoulfd be 1.0)' 
      read(1,*) Hfrac_ice,Hfrac_core
      write(*,*) 'Tidal heating time-dependence (1 = constant; 2 = sinusoidal; ...)'
      read(1,*) heat_timevar
      write(*,*) 'Period for tidal heating time-variation in Gyr, regular decay constant'
      read(1,*) Theat_Gyr,decay_heat
      write(*,*) 'Eccentricity min and max (relative to current eccentricity)'
      read(1,*) emin,emax
      write(*,*) 'Heating rate or power scale linear (1) or logarithmic (2)'
      read(1,*) heat_scale
      write(*,*) 'Composition: (1) pure H2O; (2) H2O + NH3; (3) H2O + CH3OH'
      read(1,*) composition
      write(*,*) 'Initial fractions of volatile in vol % min, max, and number of points'
      read(1,*) pCLA_init_min,pCLA_init_max,npCLA
      write(*,*) 'Number of point in iceI layer + ocean (max = ',npimax,')'
      read(1,*) npz
      write(*,*) 'Initial thickness of ice I layer (km)'
      read(1,*) z_bot_init
      write(*,*) 'Begining and finishing times (yr), number of timestep'
      read(1,*) time_start,time_end,npt
      write(*,*) 'Timescale: (1) Linear; (2) Logarithmic'
      read(1,*) itimescale
      write(*,*) 'Calculate average, min and max in ice shell over given time-window (.faslse./.true.) ?'
      read(1,*) stat_thickness
      write(*,*) 'Starting and ending time of time window (in Gyr)'
      read(1,*) tstart_Gyr,tend_Gyr
      write(*,*) 'Core heat flux from Kirk and Stevenson (1987) (1) or by Solving heat diffusion equation (2)'
      read(1,*) core_heatflux_mode
      write(*,*) 'Integration method:'
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
      R_tot       = R_tot*1.0d3
      R_core      = R_core*1.0d3
      Volume_core = 4.0d0*pi*R_core**3.0d0/3.0d0
      z_bot_init  = z_bot_init*1.0d3
      r_bot_init  = R_tot - z_bot_init
      Volume_init = 4.0d0*pi*(R_tot**3.0d0 - R_core**3.0d0)/3.0d0
      Volume_iceI = 4.0d0*pi*(R_tot**3.0d0 - r_bot_init**3.0d0)/3.0d0

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

! Initiate parameterization for Hint and etaref

      if (heat_param == 1) then
         if (Hint_min == 0.0d0) then
            heat_scale = 1
            write(*,*) 'Heating rate scale switched to linear because Hint_min = 0'
            write(*,*)
         endif
      else
         if (Pint_min == 0.0d0) then
            heat_scale = 1
            write(*,*) 'Power scale switched to linear because Pint_min = 0'
            write(*,*)
         endif
      endif
      if (heat_param == 1) then
         if (heat_scale == 1) then
            dHint = 0.0d0
            if (npHint > 1) dHint = (Hint_max - Hint_min)/dble(npHint-1)
         else
            dHint_log = 0.0d0
            if (npHint > 1) dHint_log = (dlog10(Hint_max) - dlog10(Hint_min))/dble(npHint-1)
         endif
      else
         if (heat_scale == 1) then
            dPint = 0.0d0
            if (npHint > 1) dPint = (Pint_max - Pint_min)/dble(npHint - 1)
         else
            dPint_log = 0.0d0
            if (npHint > 1) dPint_log = (dlog10(Pint_max) - dlog10(Pint_min))/dble(npHint-1)
         endif
      endif

      detaref = 0.0d0
      if (npetaref > 1) detaref = (dlog10(etaref_max) - dlog10(etaref_min))/dble(npetaref-1)

! Switch for ocean initial composition
      if (composition == 1) then
         npCLA         = 0
         pCLA_init_min = 0.0d0
         pCLA_init_max = 0.0d0
      else
         pCLA_init_min = pCLA_init_min/1.0d2
         pCLA_init_max = pCLA_init_max/1.0d2
      endif
      if ((npCLA /= 1).and.((pCLA_init_min == 0.0d0).and.(pCLA_init_max == 0.0d0))) then
         composition = 1
         npCLA = 1
      endif
      dpCLA = 0.0d0
      if (npCLA > 1) dpCLA = (pCLA_init_max - pCLA_init_min)/dble(npCLA-1)

      if (itimescale == 1) then
         dtime     = (time_end - time_start)/dble(npt-1)
         dtime_sec = dtime*year_sec
      else
         time_start_log = dlog10(time_start)
         time_end_log   = dlog10(time_end)
         dtime_log      = (time_end_log - time_start_log)/dble(npt-1)
      endif

      r_bot   = r_bot_init
      dr_init = 2.0d3

      r0_core      = 1.0d3
      npc          = nint((R_core- r0_core)/dr_init) + 1
      if (npc > npcmax) npc = npcmax
      dr_core      = (R_core - r0_core)/dble(npc-1)

      a0_ecc = 0.5d0*(emax + emin)
      a1_ecc = 0.5d0*(emax - emin)

      if (heat_param == 1) then
         filename=casename(1:nstc)//'_Hcrit-vs-zbot.dat'
      else
         filename=casename(1:nstc)//'_Pcrit-vs-zbot.dat'
      endif
      open(11,file=filename,status='unknown')
      write(11,1100) R_tot/1.0d3,R_core/1.0d3,Tsurf,g
      if (heat_param == 1) then
         filename=casename(1:nstc)//'_Hmelt-vs-zbot.dat'
      else
         filename=casename(1:nstc)//'_Pmelt-vs-zbot.dat'     
      endif
      open(12,file=filename,status='unknown')
      write(12,1200) R_tot/1.0d3,R_core/1.0d3,Tsurf,g
      filename=casename(1:nstc)//'_thermal-evolution.dat'
      open(14,file=filename,status='unknown')
      write(14,1400) R_tot/1.0d3,R_core/1.0d3,Tsurf,g
      filename=casename(1:nstc)//'_stat-ice-thickness.dat'
      open(15,file=filename,status='unknown')
      write(15,1500) R_tot/1.0d3,R_core/1.0d3,Tsurf,g
      filename=casename(1:nstc)//'_stat-stagnant-lid.dat'
      open(16,file=filename,status='unknown')
      write(16,1600) R_tot/1.0d3,R_core/1.0d3,Tsurf,g
      if (heat_param == 1) then
         write(11,1105)
         write(12,1205)     
         write(14,1405)
         write(15,1505)
         write(16,1605)
      else
         write(11,1110)
         write(12,1210)
         write(14,1410)
         write(15,1510)
         write(16,1610)
      endif
      
      do 400 iHint=1,npHint

         if (heat_param == 1) then
            if (heat_scale == 1) then
               Hint_input = Hint_min + dble(iHint-1)*dHint
            else
               Hint_log   = dlog10(Hint_min) + dble(iHint-1)*dHint_log
               Hint_input = 10.0d0**Hint_log
            endif
            H0_ice  = Hint_input*Hfrac_ice
            H0_core = Hint_input*Hfrac_core
         else
            if (heat_scale == 1) then
               Pint_input = Pint_min + dble(iHint-1)*dPint
            else
               Pint_log   = dlog10(Pint_min) + dble(iHint-1)*dPint_log
               Pint_input = 10.0d0**Pint_log
            endif
            P0_ice  = Pint_input*Hfrac_ice
            P0_core = Pint_input*Hfrac_core
         endif

         do 350 ietaref=1,npetaref

            etaref_log = dlog10(etaref_min) + dble(ietaref-1)*detaref
            etaref     = 10**etaref_log

            do 300 iCLA=1,npCLA

               pCLA_init = pCLA_init_min + dble(iCLA-1)*dpCLA
               if (composition == 2) pCLAmax = pNH3max
               if (composition == 3) pCLAmax = pCH3OHmax

               dr_bot = (R_tot - R_core)/dble(npz-1)
               r_bot            = R_tot - 0.1d0*dr_bot
               Volume_IceI_init = 4.0d0*pi*(R_tot**3.0d0 - r_bot**3.0d0)/3.0d0

               imelt  = 0
               iFneg  = 0
               do 200 iz=1,npz

                  r_bot = R_tot - dr_bot*dble(iz-1)
                  z_bot = R_tot - r_bot
                  f     = r_bot/R_tot
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

                  TH_coef = rho_iceI*Hint_ice*R_tot**2.0d0/6.0d0/kref

                  npzz = nint((z_bot - z_topTBL)/dr_bot)
                  dr_zz = (z_bot - z_topTBL)/dble(npzz - 1)

                  ibase_melt = 0
                  zbase_melt = z_bot

                  if (convecting) then
                     if (z_botTBL /= z_bot) then
                        dT_zz = (Tbot - Tm)/(z_bot - z_botTBL)
                     else
                        dT_zz = 0.0d0
                     endif
                  else
                     dT_zz = (Tbot - Tsurf)/z_bot
                  endif

                  do 50 izz=1,npzz

                     zz = z_bot - dr_zz*dble(izz - 1)
                     rr = R_tot - zz

                     call T_fusion_H2O(zz,pres,TH2O,dTbotdr)

                     if (convecting) then
                        if (zz == z_bot) then
                           Tzz = TH2O
                        else if ((zz < z_bot).and.(zz > z_botTBL)) then
                           Tzz = Tm + dT_zz*(zz - z_botTBL)
                        else
                           Tzz = Tm   
                        endif
                     else
                        TH  = TH_coef*(f*(1 + f)*(1.0d0 - R_tot/rr) + (1.0d0 - (rr/R_tot)**2.0))
                        Tzz = Tsurf - (Tbot - Tsurf)*R_tot*f*(1.0d0 - R_tot/rr)/z_bot + TH
                     endif

                     if (Tzz > TH2O) then
                        if (ibase_melt == 0) then
                           zbase_melt = zz
                           ibase_melt = 1
                        endif
                     endif

50                continue

                  zbot_km  = z_bot/1.0d3
                  dlid_km  = dlid/1.0d3
                  Fcond_mW = Fcond*1.0d3
                  Fbot_mW  = Fbot*1.0d3
                  Ftop_mW  = Ftop*1.0d3

                  if (iz > 1) then

                     if ((Fbot <= 0.0d0).and.(iFneg == 0)) then
                        iFneg = 1

                        dxdz      = (Fbot_mW - Fbot_mW_s)/(zbot_km - zbot_km_s)
                        zbot_km_w = zbot_km_s - Fbot_mW_s/dxdz
                        dxdz      = (pCLA - pCLA_s)/(zbot_km - zbot_km_s)
                        pCLA_w    = pCLA_s + dxdz*(zbot_km_w - zbot_km_s)
                        dxdz      = (Tbot - Tbot_s)/(zbot_km - zbot_km_s)
                        Tbot_w    = Tbot_s + dxdz*(zbot_km_w - zbot_km_s)
                        dxdz      = (Tm - Tm_s)/(zbot_km - zbot_km_s)
                        Tm_w      = Tm_s + dxdz*(zbot_km_w - zbot_km_s)
                        dxdz      = (Fbot_mW - Fbot_mW_s)/(zbot_km - zbot_km_s)
                        Fbot_mW_w = Fbot_mW_s + dxdz*(zbot_km_w - zbot_km_s)
                        dxdz      = (Ftop_mW - Ftop_mW_s)/(zbot_km - zbot_km_s)
                        Ftop_mW_w = Ftop_mW_s + dxdz*(zbot_km_w - zbot_km_s)
                        dxdz      = (dlid_km - dlid_km_s)/(zbot_km - zbot_km_s)
                        dlid_km_w = dlid_km_s + dxdz*(zbot_km_w - zbot_km_s)

                        if (convecting) then
                           convec = 1
                        else
                           convec = 0
                        endif
                        if (heat_param == 1) then
                           write(11,1120) Hint_ice,Pint_ice,etaref,pCLA_init,convec,zbot_km_w,pCLA_w,Tbot_w,   &
                                                                          Tm_w,Fbot_mW_w,Ftop_mW_w,dlid_km_w
                        else
                           write(11,1120) Pint_ice,Hint_ice,etaref,pCLA_init,convec,zbot_km_w,pCLA_w,Tbot_w,   &
                                                                          Tm_w,Fbot_mW_w,Ftop_mW_w,dlid_km_w
                        endif

                     endif
                     zbot_km_s = zbot_km
                     pCLA_s    = pCLA
                     Tbot_s    = Tbot
                     Tm_s      = Tm
                     Fbot_mW_s = Fbot_mW
                     Ftop_mW_s = Ftop_mW
                     dlid_km_s = dlid_km

!!!                     if ((Tm >= TH2O).and.(imelt == 0)) then
                     if ((zbase_melt < z_bot).and.(imelt == 0)) then
                        imelt = 1
                        if (convecting) then
                           convec = 1
                        else
                           convec = 0
                        endif
                        if (heat_param == 1) then
                           write(12,1220) Hint_ice,Pint_ice,etaref,pCLA_init,convec,zbot_km,pCLA,Tbot,Tm,Fbot_mW,   &
                                                                                                  Ftop_mW,dlid_km
                        else
                           write(12,1220) Hint_ice,Pint_ice,etaref,pCLA_init,convec,zbot_km,pCLA,Tbot,Tm,Fbot_mW,   &
                                                                                                  Ftop_mW,dlid_km
                        endif

                     endif

                  endif

200            continue

               r_bot = r_bot_init

               npt_stat  = 0
               time_stat = 0.0d0

!---> Now calculates time evolution of radial structure
               do 250 i=1,npt

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

                  z_bot = R_tot - r_bot
                  Volume_iceI = 4.0d0*pi*(R_tot**3.0d0 - r_bot**3.0d0)/3.0d0

                  if (heat_param == 1) then
                     P0_ice  = H0_ice*rho_iceI*Volume_iceI
                     P0_core = H0_core*rho_core*Volume_core
                  else
                     H0_ice  = P0_ice/rho_iceI/Volume_iceI
                     H0_core = P0_core/rho_core/Volume_core
                  endif

                  call Time_dependent_tidal_heating(heat_timevar,time,Theat_Gyr,a0_ecc,a1_ecc,decay_heat,H0_ice,Hint_ice)
                  call Time_dependent_tidal_heating(heat_timevar,time,Theat_Gyr,a0_ecc,a1_ecc,decay_heat,H0_core,Hint_core)

                  Pint_ice  = Hint_ice*rho_iceI*Volume_iceI
                  Pint_core = Hint_core*rho_core*Volume_core


                  if (heat_param == 1) then
                     write(20,*) Hint_input,time,heat_timevar,Hint_ice,Pint_ice
                  else
                     write(20,*) Pint_input,time,heat_timevar,Hint_ice,Pint_ice
                  endif

                  Pint_ice  = Hint_ice*rho_iceI*Volume_iceI
                  Pint_core = Hint_core*rho_core*Volume_core

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
!!!                  if (Tzz > TH2O) then
!!!
!!!                     imelt   = 1
!!!                     r_bot   = r_bot_old
!!!                     drbotdt = drbotdt_old
!!!
!!!                  else

                     r_bot_old   = r_bot
                     drbotdt_old = drbotdt

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

                     time_Gyr  = time/1.0d9
                     dtime_Gyr = dtime/1.0d9
                     P_TW      = P_core/1.0d12
                     Pres_MPa  = pres*1.0d3
                     zbot_km   = z_bot/1.0d3
                     rbot_km   = r_bot/1.0d3
                     Fcore_mW  = F_core*1.0d3
                     Fcond_mW  = Fcond*1.0d3
                     Fbot_mW   = Fbot*1.0d3
                     Ftop_mW   = Ftop*1.0d3
                     dlid_km   = dlid/1.0d3

!!!                  write(*,1030) dTbotdr,dTadiadr,dPower,cool_ocean,Lsolid,drbotdt
!!!                  endif

! Integrate bottom radius to next point according to chosen method
                  rx(1)    = r_bot
                  drxdt(1) = drbotdt
                  call integrate(iint,n_int,eps,time,dtime,rx,drxdt)
                  r_bot = rx(1)
                  if (r_bot > R_core) then
                     r_bot = rx(1)
                  else
                     r_bot = R_core
                  endif

! Average, min and max in ice shell thickness within prescribed time-window
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

250            continue

               call T_fusion_H2O(z_bot,pres,TH2O,dTbotdr)

               rnpt_stat   = 1.0d0/dble(npt_stat)
               zbot_km_avg = zbot_km_avg/time_stat
               dlid_km_avg = dlid_km_avg/time_stat

               if (convecting) then
                  convec = 1
               else
                  convec = 0
               endif
               if (heat_param == 1) then
                  write(*,1040) Hint_input,etaref,pCLA_init,convec,Pres_MPa,zbot_km,pCLA,Tbot,Tm,TH2O,Fbot_mW,Ftop_mW,dlid_km,   &
                                                                                                                     T_core(1)
                  write(14,1420) Hint_input,etaref,pCLA_init,convec,Pres_MPa,zbot_km,pCLA,Tbot,Tm,TH2O,Fbot_mW,Ftop_mW,dlid_km,   &
                                                                                                                      T_core(1)
                  write(15,1520) Hint_input,etaref,pCLA_init,zbot_km_avg,zbot_km_min,zbot_km_max
                  write(16,1620) Hint_input,etaref,pCLA_init,dlid_km_avg,dlid_km_min,dlid_km_max
               else
                  write(*,1040) Pint_input,etaref,pCLA_init,convec,Pres_MPa,zbot_km,pCLA,Tbot,Tm,TH2O,Fbot_mW,Ftop_mW,dlid_km,   &
                                                                                                                     T_core(1)
                  write(14,1420) Pint_input,etaref,pCLA_init,convec,Pres_MPa,zbot_km,pCLA,Tbot,Tm,TH2O,Fbot_mW,Ftop_mW,dlid_km,   &
                                                                                                                      T_core(1)
                  write(15,1520) Pint_input,etaref,pCLA_init,zbot_km_avg,zbot_km_min,zbot_km_max
                  write(16,1620) Pint_input,etaref,pCLA_init,dlid_km_avg,dlid_km_min,dlid_km_max
               endif

300         continue
350      continue
400   continue
      close(11)
      close(12)
      close(14)
      close(15)
      close(16)

1000  format(' Tsurf = ',f8.2,' K')
1005  format(' g = ',f6.3,' m/s^2')
1010  format(' k_I = ',f7.3,' W/m/K ;  Cp_I = ',f8.2,' J/kg/K ; kappa_I = ',e12.7,' m^2/s')
1030  format(6(1x,d11.5))
1040  format(2(1x,d11.5),1x,f8.4,1x,i2,1x,f7.2,1x,f8.3,1x,f8.5,3(1x,f7.2),2(1x,f8.3),1x,f7.2,1x,f8.2)
1100  format('R_tot =',f7.1,' R_core = ',f7.1,'  Tsurf = ',f8.2,'  g = ',f6.3)
1105  format(5x,'Hint',7x,'Pint',7x,'etaref',6x,'pCLA',1x,'conv',3x,'zbot',5x,'%vol',4x,'Tbot',5x,'Tm' 6x,     &
                                                                               'Fbot',5x,'Ftop',5x,'dlid')
1110  format(5x,'Pint',7x,'Hint',7x,'etaref',6x,'pCLA',1x,'conv',3x,'zbot',5x,'%vol',4x,'Tbot',5x,'Tm' 6x,     &
                                                                               'Fbot',5x,'Ftop',5x,'dlid')
1120  format(3(1x,d11.5),1x,f8.4,1x,i2,1x,f8.3,1x,f8.5,2(1x,f7.2),2(1x,f8.3),1x,f7.2)
1200  format('R_tot =',f7.1,' R_core = ',f7.1,'  Tsurf = ',f8.2,'  g = ',f6.3)
1205  format(5x,'Hint',7x,'Pint',7x,'etaref',6x,'pCLA',1x,'conv',3x,'zbot',5x,'%vol',4x,'Tbot',5x,'Tm' 6x,     &
                                                                               'Fbot',5x,'Ftop',5x,'dlid')
1210  format(5x,'Pint',7x,'Hint',7x,'etaref',6x,'pCLA',1x,'conv',3x,'zbot',5x,'%vol',4x,'Tbot',5x,'Tm' 6x,     &
                                                                               'Fbot',5x,'Ftop',5x,'dlid')
1220  format(3(1x,d11.5),1x,f8.4,1x,i2,1x,f8.3,1x,f8.5,2(1x,f7.2),2(1x,f8.3),1x,f7.2)
1400  format('R_tot =',f7.1,' R_core = ',f7.1,'  Tsurf = ',f8.2,'  g = ',f6.3)
1405  format(5x,'Hint',7x,'etaref',6x,'pCLA',1x,'conv',1x,'P(MPa)',3x,'zbot',5x,'%vol',4x,'Tbot',5x,'Tm',      &
                                                     6x,'TH2O',6x,'Fbot',5x,'Ftop',5x,'dlid',3x,'T_core')
1410  format(5x,'Pint',7x,'etaref',6x,'pCLA',1x,'conv',1x,'P(MPa)',3x,'zbot',5x,'%vol',4x,'Tbot',5x,'Tm',      &
                                                     6x,'TH2O',6x,'Fbot',5x,'Ftop',5x,'dlid',3x,'T_core')
1420  format(2(1x,d11.5),1x,f8.4,1x,i2,1x,f7.2,1x,f8.3,1x,f8.5,3(1x,f7.2),2(1x,f8.3),1x,f7.2,1x,f8.2)
1500  format('R_tot =',f7.1,' R_core = ',f7.1,'  Tsurf = ',f8.2,'  g = ',f6.3)
1505  format(5x,'Hint',7x,'etaref',6x,'pCLA',3x,'zbot_avg',1x,'zbot_min',1x,'zbot_max')
1510  format(5x,'Pint',7x,'etaref',6x,'pCLA',3x,'zbot_avg',1x,'zbot_min',1x,'zbot_max')
1520  format(2(1x,d11.5),1x,f8.4,3(1x,f8.3))
1600  format('R_tot =',f7.1,' R_core = ',f7.1,'  Tsurf = ',f8.2,'  g = ',f6.3)
1605  format(5x,'Hint',7x,'etaref',6x,'pCLA',3x,'dlid_avg',1x,'dlid_min',1x,'dlid_max')
1610  format(5x,'Pint',7x,'etaref',6x,'pCLA',3x,'dlid_avg',1x,'dlid_min',1x,'dlid_max')
1620  format(2(1x,d11.5),1x,f8.4,3(1x,f8.3))

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

      r_bot     = rx(1)
      z_bot     = R_tot - r_bot
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

