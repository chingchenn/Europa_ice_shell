!************************************************************************************|
!  Modules needed for thermal evolution calculation                                  |
!************************************************************************************|

module thermal_evolution_mod

      implicit real(8) (a-h,o-z)

!---> Parameters
      
!   Size of arrays
      integer, parameter :: nmax = 2001, npcmax = 2001, npimax = 2001, npzmax = npcmax + 2*npimax, npwmax = 101

!   Mathematical and physical constants
      real(8), parameter :: pi = 4.0d0*datan(1.0d0), cVol = 4.0d0*pi/3.0d0                              !! Geometrical constants
      real(8), parameter :: R = 8.32d0                                                                  !! Physical constants
      real(8), parameter :: year_sec = 3.15576d7                                                        !! Number of seconds per year

!   Properties of ices and volatiles
      real(8), parameter :: rho_iceI = 9.20d2, alfa_iceI = 1.56d-4, L_iceI = 284.0d3                    !! Properties of water iceI
!!!      real(8), parameter :: k_iceI_a = 488.0d0, k_iceI_b = 0.4685d0                                     !! Conductivity from Kirk and Stevenson (1987)
      real(8), parameter :: k_iceI_a = 566.8d0, k_iceI_b = 0.0d0                                        !! Conductvity from Andersson and Suga (1994)
      real(8), parameter :: Cp_iceI_a = 7.037d0, Cp_iceI_b = 185.0d0                                    !! Heat capacity of ice I
      real(8), parameter :: alfa_H2O = 3.0d-4, rho_H2O = 1.0d3, Cp_H2O = 4.18d3                         !! Properties of liquid water
      real(8), parameter :: alfa_NH3 = 3.0d-4, rho_NH3 = 0.734d3, Cp_NH3 = 4.07d3                       !! Properties of liquid ammonia
      real(8), parameter :: alfa_CH3OH = 3.0d-4, rho_CH3OH = 0.791d3, Cp_CH3OH = 2.48d3                 !! Properties of liquid methanol
      real(8), parameter :: pNH3max = 0.321d0, pCH3OHmax = 0.88d0                                       !! Eutectic of H2O-NH3 and H2O-CH3OH

!   Parameters of heat-flux and temperature scaling laws (average and uncertainties)
      real(8), parameter :: ac = 0.258d0, betac = 0.321d0                                               !! Heat flux scaling for constant viscosity
      real(8), parameter :: aparam_pos = 1.46d0, aparam_neg = 1.57d0                                    !! Top heat flux scaling for T-dep. visco, 3D-Cartesian
!!!      real(8), parameter :: beta = 0.270d0, cgam = 1.21d0, ck = 0.82d0, cf = 0.22d0 
      real(8), parameter :: beta = 0.270d0, cgam = 1.21d0, ck = 0.82d0, cf = 0.0d0 
      real(8), parameter :: c1 = 1.23d0, c2 = 0.1d0, c3 = 1.5d0                                         !! Temperature scaling for T-dep. visco., 3D-Cartesian/spherical,
      real(8), parameter :: c4pos = 3.48d0, c5pos = -2.26d0, c6pos = 1.00d0, c7pos = 0.250d0            !! ... and internal heating
      real(8), parameter :: c4neg = 4.36d0, c5neg = -3.00d0, c6neg = 1.72d0, c7neg = 0.333d0
      real(8), parameter :: atau = 0.37d0, btau = 0.91d0, ctau = 0.16d0

!---> Common scalar variables
      logical :: internal_heating,variable_conductivity                                                 !! Presence of internal heating, variable conductivity
      integer :: core_heatflux_mode,chondritic_model                                                    !! Method of core heat flux calculation, chondritic composition model
      integer :: scaling_law,geometry                                                                   !! Mode of convection, geometry
      integer :: composition,reference_conductivity                                                     !! Ocean omposition, reference conductivity, 
      integer :: npc                                                                                    !! Number of radial nodes in the core
      integer :: ichond_model                                                                           !! Chondritic model
      real(8) :: r0_core                                                                                !! Core minimum radius
      real(8) :: k_iceI,Cp_iceI,kappa_iceI                                                              !! Properties of ice I at surface temperature
      real(8) :: kappa_input,k_input                                                                    !! Input thermal diffusivity and conductivity (if ikver=3)
      real(8) :: R_tot,Tsurf,Volume_init                                                                !! Total and core radii, surface temperature, initial volume of ocean
      real(8) :: R_core,kappa_core,rho_core,Cp_core,k_core,dr_core                                      !! Core properties
      real(8) :: g,rhog                                                                                 !! Surface gravity acceleration, density*g
      real(8) :: etaref,Eact                                                                            !! Reference viscosity, activation energy for viscosity
      real(8) :: Hint_ice,Hint_core                                                                     !! Internal heating rate in ice layer and core
      real(8) :: z_I_to_III                                                                             !! Depth of the ice I ---> ice III transition
      real(8) :: pCLA_init,pCLAmax
      real(8) :: alfa_mix,rho_mix,Cp_mix
      real(8) :: time_save
      real(8) :: aparam,c4,c5,c6,c7

!---> Common array variables
      real(8) :: T_core(npcmax)                                                                         !! Core temperature profile

contains

!-------------------------------------------------------------------------------------
! Calculate liquidus for a mix at a given depth                                      |
!-------------------------------------------------------------------------------------
      subroutine T_fusion_mixture(z,p,pCLA,TH2O,TNH3_15pc,TNH3_eutec,TCH3OH,T)

      implicit real(8) (a-h,o-z)
      real(8), dimension(nmax) :: pint,Tint

      call T_fusion_H2O(z,p,TH2O,dTH2Odr)
      call T_fusion_NH3_15pc(z,p,TNH3_15pc)
      call T_fusion_NH3_eutec(z,p,TNH3_eutec)
      call T_fusion_Methanol(z,p,TCH3OH)

!---> Water-ammonia system
      if (composition == 2) then

         if (pCLA == pNH3max) then
            T = TNH3_eutec
         else
            pint(1) = 0.0d0
            pint(2) = 1.50d-01
            pint(3) = 3.21d-01
            Tint(1) = TH2O
            Tint(2) = TNH3_15pc
            Tint(3) = TNH3_eutec
            x0 = pCLA
            call polint(3,pint,Tint,x0,y0)
            T = y0
         endif

!---> alpha, rho, & heat capacity of liquid ocean
         alfa_mix = (1.0d0 - pCLA)*alfa_H2O + pCLA*alfa_NH3
         rho_mix = (1.0d0 - pCLA)*rho_H2O + pCLA*rho_NH3
!!!         Cp_mix = (1.0d0 - pCLA)*Cp_H2O + pCLA*Cp_NH3
         Cp_mix = Cp_H2O*Cp_NH3/(pCLA*Cp_H2O + (1.0d0 - pCLA)*Cp_NH3)

      endif

!---> Water-methanol system (2D interpolation)
      if (composition == 3) then

         if (pCLA == pCH3OHmax) then
            call T_fusion_H2O(z,p,TH2O,dTH2Odr)
            call T_fusion_Methanol(z,p,TCH3OH)
            TCH3OH_eutec = pCH3OHmax*TCH3OH + (1.0d0 - pCH3OHmax)*TH2O
            T = TCH3OH_eutec
         else
            ax = (TH2O - TCH3OH)/(273.15d0 - 174.45d0)
            pCLA2 = pCLA*pCLA
            if (pCLA <= 0.46d0) T_p0 = 273.15d0 - 83.208d0*pCLA - 41.529d0*pCLA2
            if ((pCLA > 0.46d0).and.(pCLA <= 0.70d0)) T_p0 = 173.37d0 + 339.4d0*pCLA - 489.52d0*pCLA2
            if ((pCLA > 0.70d0).and.(pCLA <= 0.88d0)) T_p0 = -69.528d0 + 676.58d0*pCLA - 476.39d0*pCLA2
            if (pCLA > 0.88d0) T_p0 = -130.27d0 + 485.27d0*pCLA - 180.55d0*pCLA2
            T = TH2O + ax*(T_p0 - 273.15d0)
         endif

!---> alpha, rho, & heat capacity of liquid ocean
         alfa_mix = (1.0d0 - pCLA)*alfa_H2O + pCLA*alfa_CH3OH
         rho_mix = (1.0d0 - pCLA)*rho_H2O + pCLA*rho_CH3OH
!!!         Cp_mix = (1.0d0 - pCLA)*Cp_H2O + pCLA*Cp_CH3OH
         Cp_mix = Cp_H2O*Cp_CH3OH/(pCLA*Cp_H2O + (1.0d0 - pCLA)*Cp_CH3OH)

      endif

!---> Uncomment to take values of alpha, rho, & heat capacity for pure water
!!!      alfa_mix = alfa_H2O
!!!      rho_mix = rho_H2O
!!!      Cp_mix = Cp_H2O

      return

      end subroutine T_fusion_mixture

!-------------------------------------------------------------------------------------
!  Phase diagram of pure water                                                       |
!-------------------------------------------------------------------------------------
      subroutine T_fusion_H2O(z,p,T,dTdr)

      implicit real(8) (a-h,o-z)

      data p1,p2,p3,p4,p5/0.0d0,0.207d0,0.346d0,0.625d0,2.21d0/
      data T1,T2,T3,T4,T5/2.7316d2,2.5115d2,2.5615d2,2.7331d2,3.55d2/
      data a1,a2,a3,a4,a5/-0.3952d0,0.062d0,0.410d0,0.707d0,0.5342d0/
      data d1,d2,d3,d4,d5/9.0d0,6.0d1,8.1d0,4.46d0,5.52d0/

      p = rhog*z/1.0d9

      if (p < p2) then
         T = T1*(1.0d0 + (p - p1)/a1)**(1.0d0/d1)
         dTdp = T/a1/d1/(1.0d0 + (p - p1)/a1)
      endif    
      if ((p >= p2).and.(p <= p3)) then
         T = T2*(1.0d0 + (p - p2)/a2)**(1.0d0/d2)
         dTdp = T/a2/d2/(1.0d0 + (p - p2)/a2)
      endif   
      if ((p > p3).and.(p <= p4)) then
         T = T3*(1.0d0 + (p - p3)/a3)**(1.0d0/d3)
         dTdp = T/a3/d3/(1.0d0 + (p - p3)/a3)
      endif   
      if ((p > p4).and.(p <= p5)) then
         T = T4*(1.0d0 + (p - p4)/a4)**(1.0d0/d4)
         dTdp = T/a4/d4/(1.0d0 + (p - p4)/a4)
       endif   
      if (p > p5) Tbase = T5*(1.0d0 + (p - p5)/a5)**(1.0d0/c5)

      dTdr = -rhog*dTdp/1.0d9

      return

      end subroutine T_fusion_H2O

!-------------------------------------------------------------------------------------
!  Phase diagram of pure methanol                                                    |
!-------------------------------------------------------------------------------------
      subroutine T_fusion_Methanol(z,p,T)

      implicit real(8) (a-h,o-z)

      data p1,p2/0.00d0,0.25d0/
      data T1,T2/1.7574d2,1.8148d2/
      data a11,a12,a21,a22/1.4005d2,-7.8269d1,1.1099d2,-4.507d1/

      p=rhog*z/1.0d9

      if (p <= p2) T = T1 + a11*p + a12*p**2.d0
      if (p > p2) T = T2 + a21*p + a22*p**2.d0

      return

      end subroutine T_fusion_Methanol

!-------------------------------------------------------------------------------------
!  Phase diagram for water + 15% NH3                                                 |
!-------------------------------------------------------------------------------------
      subroutine T_fusion_NH3_15pc(z,p,T)

      implicit real(8) (a-h,o-z)

      data p1,p2/0.00d0,0.21d0/
      data T1,T2/2.53d2,1.935d2/
      data a11,a12,a21,a22/-8.51d1,-2.818d2,1.498d2,-5.05d1/

      p=rhog*z/1.d9

      if (p <= p2) T = T1 + a11*p + a12*p**2.d0
      if (p > p2) T = T2 + a21*p + a22*p**2.d0

      return

      end subroutine T_fusion_NH3_15pc

!-------------------------------------------------------------------------------------
!  Phase diagram of NH3 at eutectic (32.3% NH3)                                      |
!-------------------------------------------------------------------------------------
      subroutine T_fusion_NH3_eutec(z,p,T)

      implicit real(8) (a-h,o-z)

      data p1,p2/0.0d0,0.297d0/
      data T1,T2,T3/1.7629d02,1.769d02,1.457d02/
      data a11,a12,a21,a22/2.43d01,-7.7583d01,1.1658d2,-3.92d1/

      p=rhog*z/1.0d9
!!!      if (p.le.p2) T = T1 + a11*p + a12*p**2.d0
!!!      if (p.gt.p2) T = T3 + a21*p + a22*p**2.d0

      T = 1.74d2 + 1.5d1*p

      return

      end subroutine T_fusion_NH3_eutec

!-------------------------------------------------------------------------------------
! Polynomial interpolation                                                           |
!-------------------------------------------------------------------------------------
      subroutine polint(n,x,y,xint,yint)

      implicit real(8) (a-h,o-z)
      real(8), dimension(nmax) :: x,y
      real(8), dimension(nmax,nmax) :: p

      do 10 i=1,n
         p(1,i) = y(i)
10    continue
      do 30 i=2,n
         m = i - 1
         do 20 j=1,n-i+1
            a1 = xint - x(j+m)
            a2 = x(j) - xint
            a3 = x(j) - x(j+m)
            p(i,j) = (a1*p(i-1,j) + a2*p(i-1,j+1))/a3
20       continue
30    continue
      yint = p(n,1)

      return

      end subroutine polint

end module thermal_evolution_mod

!-------------------------------------------------------------------------------------
MODULE precision_mod

IMPLICIT NONE

INTEGER, PARAMETER :: my_kind = SELECTED_REAL_KIND(8) 

END MODULE precision_mod

