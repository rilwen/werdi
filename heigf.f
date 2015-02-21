      subroutine heigf(kx, ky, kz, C4_SI, exx, BG, thetaM, phiM, H6x6)

      implicit none

      double precision kx, ky, kz


      double precision C4_SI, exx, BG, thetaM, phiM

      double precision ialpha, hbar, c, e, m, gamma1, gamma2, gamma3,
     $     kappa, g0, Lambda, Eg, eta, hartree_eV, au_GPa, c_SI, C4,
     $     b, d, c11, c12, c44

      parameter ( ialpha = 137.0359991d0, hbar = 1, c = ialpha, e = 1,
     $     m = 1, hartree_eV = 27.2113845d0, au_GPa = 29421.012d0,
     $     gamma1 = 6.85d0, gamma2 = 2.1d0, gamma3 = 2.9d0,
     $     kappa = 1.2d0, g0 = 2.0023193043718d0,
     $     Lambda = 0.341d0/hartree_eV, Eg=1.519d0/hartree_eV,
     $     eta=Lambda/(Eg+Lambda), c_SI=299792458,
     $     b = -2.0d0/hartree_eV, d = -4.8d0/hartree_eV,
     $     c11 = 119d0/au_GPa, c12 = 53.8d0/au_GPa,
     $     c44 = 59.5d0/au_GPa )


      double precision cg(3, 3)
      double precision H6x6Re(6, 6), H6x6Im(6, 6)
      double complex H6x6(6, 6)


      integer i, j


      C4 = C4_SI/c_SI*c*hbar

      do i = 1, 6
         do j = 1, 6
            H6x6Re(i, j) = 0
            H6x6Im(i, j) = 0
         end do
      end do

      cg(1,1) = 4 * (5 * c11 - 49 * c12 + c44) * c44 * exx / (9 * c11 **
     # 2 + 101 * c44 * c11 + 9 * c11 * c12 - 34 * c44 * c12 + 4 * c44 **
     # 2 - 18 * c12 ** 2)
      cg(1,2) = -3 * (5 * c11 ** 2 + c44 * c11 + 5 * c11 * c12 + 2 * c44
     # * c12 - 10 * c12 ** 2) * exx / (9 * c11 ** 2 + 101 * c44 * c11 + 
     #9 * c11 * c12 - 34 * c44 * c12 + 4 * c44 ** 2 - 18 * c12 ** 2)
      cg(1,3) = -3 * (5 * c11 ** 2 + c44 * c11 + 5 * c11 * c12 + 2 * c44
     # * c12 - 10 * c12 ** 2) * exx / (9 * c11 ** 2 + 101 * c44 * c11 + 
     #9 * c11 * c12 - 34 * c44 * c12 + 4 * c44 ** 2 - 18 * c12 ** 2)
      cg(2,1) = -3 * (5 * c11 ** 2 + c44 * c11 + 5 * c11 * c12 + 2 * c44
     # * c12 - 10 * c12 ** 2) * exx / (9 * c11 ** 2 + 101 * c44 * c11 + 
     #9 * c11 * c12 - 34 * c44 * c12 + 4 * c44 ** 2 - 18 * c12 ** 2)
      cg(2,2) = 4 * (27 * c11 - 5 * c12 + c44) * c44 * exx / (9 * c11 **
     # 2 + 101 * c44 * c11 + 9 * c11 * c12 - 34 * c44 * c12 + 4 * c44 **
     # 2 - 18 * c12 ** 2)
      cg(2,3) = -(9 * c11 ** 2 - 7 * c44 * c11 + 9 * c11 * c12 - 14 * c4
     #4 * c12 - 18 * c12 ** 2) * exx / (9 * c11 ** 2 + 101 * c44 * c11 +
     # 9 * c11 * c12 - 34 * c44 * c12 + 4 * c44 ** 2 - 18 * c12 ** 2)
      cg(3,1) = -3 * (5 * c11 ** 2 + c44 * c11 + 5 * c11 * c12 + 2 * c44
     # * c12 - 10 * c12 ** 2) * exx / (9 * c11 ** 2 + 101 * c44 * c11 + 
     #9 * c11 * c12 - 34 * c44 * c12 + 4 * c44 ** 2 - 18 * c12 ** 2)
      cg(3,2) = -(9 * c11 ** 2 - 7 * c44 * c11 + 9 * c11 * c12 - 14 * c4
     #4 * c12 - 18 * c12 ** 2) * exx / (9 * c11 ** 2 + 101 * c44 * c11 +
     # 9 * c11 * c12 - 34 * c44 * c12 + 4 * c44 ** 2 - 18 * c12 ** 2)
      cg(3,3) = 4 * (27 * c11 - 5 * c12 + c44) * c44 * exx / (9 * c11 **
     # 2 + 101 * c44 * c11 + 9 * c11 * c12 - 34 * c44 * c12 + 4 * c44 **
     # 2 - 18 * c12 ** 2)

     
      H6x6Re(1,1) = -hbar ** 2 / m * gamma1 * kx ** 2 / 0.2D1 - hbar ** 
     #2 / m * gamma1 * ky ** 2 / 0.2D1 - hbar ** 2 / m * gamma1 * kz ** 
     #2 / 0.2D1 - hbar ** 2 / m * gamma2 * kx ** 2 / 0.2D1 - hbar ** 2 /
     # m * gamma2 * ky ** 2 / 0.2D1 + hbar ** 2 / m * gamma2 * kz ** 2 +
     # b * cg(1,1) / 0.2D1 + b * cg(2,2) / 0.2D1 - b * cg(3,3) + 0.3D1 *
     # BG * cos(thetaM) + 0.3D1 / 0.2D1 * C4 * cg(3,1) * kx - 0.3D1 / 0.
     #2D1 * C4 * cg(3,2) * ky
      H6x6Re(1,2) = hbar ** 2 / m * gamma3 * kx * kz * sqrt(0.3D1) - d *
     # cg(1,3) + BG * sin(thetaM) * cos(phiM) * sqrt(0.3D1) + C4 * sqrt(
     #0.3D1) * cg(1,2) * ky / 0.2D1 - C4 * sqrt(0.3D1) * cg(1,3) * kz / 
     #0.2D1
      H6x6Re(1,3) = hbar ** 2 / m * gamma2 * kx ** 2 * sqrt(0.3D1) / 0.2
     #D1 - hbar ** 2 / m * gamma2 * ky ** 2 * sqrt(0.3D1) / 0.2D1 - b * 
     #cg(1,1) * sqrt(0.3D1) / 0.2D1 + b * cg(2,2) * sqrt(0.3D1) / 0.2D1
      H6x6Re(1,5) = -hbar ** 2 / m * gamma3 * kx * kz * sqrt(0.2D1) * sq
     #rt(0.3D1) / 0.2D1 + d * cg(3,1) * sqrt(0.2D1) / 0.2D1 + BG * sin(t
     #hetaM) * cos(phiM) * sqrt(0.2D1) * sqrt(0.3D1) - C4 * sqrt(0.2D1) 
     #* sqrt(0.3D1) * cg(1,2) * ky / 0.4D1 + C4 * sqrt(0.2D1) * sqrt(0.3
     #D1) * cg(1,3) * kz / 0.4D1 + C4 * eta * sqrt(0.2D1) * sqrt(0.3D1) 
     #* cg(1,2) * ky / 0.8D1 - C4 * eta * sqrt(0.2D1) * sqrt(0.3D1) * cg
     #(1,3) * kz / 0.8D1
      H6x6Re(1,6) = -hbar ** 2 / m * gamma2 * kx ** 2 * sqrt(0.2D1) * sq
     #rt(0.3D1) / 0.2D1 + hbar ** 2 / m * gamma2 * ky ** 2 * sqrt(0.2D1)
     # * sqrt(0.3D1) / 0.2D1 + b * cg(1,1) * sqrt(0.2D1) * sqrt(0.3D1) /
     # 0.2D1 - b * cg(2,2) * sqrt(0.2D1) * sqrt(0.3D1) / 0.2D1
      H6x6Re(2,1) = hbar ** 2 / m * gamma3 * kx * kz * sqrt(0.3D1) - d *
     # cg(1,3) + BG * sin(thetaM) * cos(phiM) * sqrt(0.3D1) + C4 * sqrt(
     #0.3D1) * cg(1,2) * ky / 0.2D1 - C4 * sqrt(0.3D1) * cg(1,3) * kz / 
     #0.2D1
      H6x6Re(2,2) = -hbar ** 2 / m * gamma1 * kx ** 2 / 0.2D1 - hbar ** 
     #2 / m * gamma1 * ky ** 2 / 0.2D1 - hbar ** 2 / m * gamma1 * kz ** 
     #2 / 0.2D1 + hbar ** 2 / m * gamma2 * kx ** 2 / 0.2D1 + hbar ** 2 /
     # m * gamma2 * ky ** 2 / 0.2D1 - hbar ** 2 / m * gamma2 * kz ** 2 -
     # b * cg(1,1) / 0.2D1 - b * cg(2,2) / 0.2D1 + b * cg(3,3) + BG * co
     #s(thetaM) + C4 * cg(3,1) * kx / 0.2D1 - C4 * cg(3,2) * ky / 0.2D1
      H6x6Re(2,3) = 0.2D1 * BG * sin(thetaM) * cos(phiM) + C4 * cg(1,2) 
     #* ky - C4 * cg(1,3) * kz
      H6x6Re(2,4) = hbar ** 2 / m * gamma2 * kx ** 2 * sqrt(0.3D1) / 0.2
     #D1 - hbar ** 2 / m * gamma2 * ky ** 2 * sqrt(0.3D1) / 0.2D1 - b * 
     #cg(1,1) * sqrt(0.3D1) / 0.2D1 + b * cg(2,2) * sqrt(0.3D1) / 0.2D1
      H6x6Re(2,5) = -hbar ** 2 / m * gamma2 * kx ** 2 * sqrt(0.2D1) / 0.
     #2D1 - hbar ** 2 / m * gamma2 * ky ** 2 * sqrt(0.2D1) / 0.2D1 + hba
     #r ** 2 / m * gamma2 * kz ** 2 * sqrt(0.2D1) + b * cg(1,1) * sqrt(0
     #.2D1) / 0.2D1 + b * cg(2,2) * sqrt(0.2D1) / 0.2D1 - b * cg(3,3) * 
     #sqrt(0.2D1) - 0.2D1 * BG * cos(thetaM) * sqrt(0.2D1) + C4 * sqrt(0
     #.2D1) * cg(3,1) * kx / 0.2D1 - C4 * sqrt(0.2D1) * cg(3,2) * ky / 0
     #.2D1 - C4 * sqrt(0.2D1) * eta * cg(3,1) * kx / 0.4D1 + C4 * sqrt(0
     #.2D1) * eta * cg(3,2) * ky / 0.4D1
      H6x6Re(2,6) = 0.3D1 / 0.2D1 * hbar ** 2 / m * gamma3 * kx * kz * s
     #qrt(0.2D1) - d * sqrt(0.3D1) * cg(3,1) * sqrt(0.2D1) / 0.2D1 + BG 
     #* sin(thetaM) * cos(phiM) * sqrt(0.2D1) - C4 * sqrt(0.2D1) * cg(1,
     #2) * ky / 0.4D1 + C4 * sqrt(0.2D1) * cg(1,3) * kz / 0.4D1 + C4 * e
     #ta * sqrt(0.2D1) * cg(1,2) * ky / 0.8D1 - C4 * eta * sqrt(0.2D1) *
     # cg(1,3) * kz / 0.8D1
      H6x6Re(3,1) = hbar ** 2 / m * gamma2 * kx ** 2 * sqrt(0.3D1) / 0.2
     #D1 - hbar ** 2 / m * gamma2 * ky ** 2 * sqrt(0.3D1) / 0.2D1 - b * 
     #cg(1,1) * sqrt(0.3D1) / 0.2D1 + b * cg(2,2) * sqrt(0.3D1) / 0.2D1
      H6x6Re(3,2) = 0.2D1 * BG * sin(thetaM) * cos(phiM) + C4 * cg(1,2) 
     #* ky - C4 * cg(1,3) * kz
      H6x6Re(3,3) = -hbar ** 2 / m * gamma1 * kx ** 2 / 0.2D1 - hbar ** 
     #2 / m * gamma1 * ky ** 2 / 0.2D1 - hbar ** 2 / m * gamma1 * kz ** 
     #2 / 0.2D1 + hbar ** 2 / m * gamma2 * kx ** 2 / 0.2D1 + hbar ** 2 /
     # m * gamma2 * ky ** 2 / 0.2D1 - hbar ** 2 / m * gamma2 * kz ** 2 -
     # b * cg(1,1) / 0.2D1 - b * cg(2,2) / 0.2D1 + b * cg(3,3) - BG * co
     #s(thetaM) - C4 * cg(3,1) * kx / 0.2D1 + C4 * cg(3,2) * ky / 0.2D1
      H6x6Re(3,4) = -hbar ** 2 / m * gamma3 * kx * kz * sqrt(0.3D1) + d 
     #* cg(1,3) + BG * sin(thetaM) * cos(phiM) * sqrt(0.3D1) + C4 * sqrt
     #(0.3D1) * cg(1,2) * ky / 0.2D1 - C4 * sqrt(0.3D1) * cg(1,3) * kz /
     # 0.2D1
      H6x6Re(3,5) = 0.3D1 / 0.2D1 * hbar ** 2 / m * gamma3 * kx * kz * s
     #qrt(0.2D1) - d * sqrt(0.3D1) * cg(3,1) * sqrt(0.2D1) / 0.2D1 - BG 
     #* sin(thetaM) * cos(phiM) * sqrt(0.2D1) + C4 * sqrt(0.2D1) * cg(1,
     #2) * ky / 0.4D1 - C4 * sqrt(0.2D1) * cg(1,3) * kz / 0.4D1 - C4 * e
     #ta * sqrt(0.2D1) * cg(1,2) * ky / 0.8D1 + C4 * eta * sqrt(0.2D1) *
     # cg(1,3) * kz / 0.8D1
      H6x6Re(3,6) = hbar ** 2 / m * gamma2 * kx ** 2 * sqrt(0.2D1) / 0.2
     #D1 + hbar ** 2 / m * gamma2 * ky ** 2 * sqrt(0.2D1) / 0.2D1 - hbar
     # ** 2 / m * gamma2 * kz ** 2 * sqrt(0.2D1) - b * cg(1,1) * sqrt(0.
     #2D1) / 0.2D1 - b * cg(2,2) * sqrt(0.2D1) / 0.2D1 + b * cg(3,3) * s
     #qrt(0.2D1) - 0.2D1 * BG * cos(thetaM) * sqrt(0.2D1) + C4 * sqrt(0.
     #2D1) * cg(3,1) * kx / 0.2D1 - C4 * sqrt(0.2D1) * cg(3,2) * ky / 0.
     #2D1 - C4 * sqrt(0.2D1) * eta * cg(3,1) * kx / 0.4D1 + C4 * sqrt(0.
     #2D1) * eta * cg(3,2) * ky / 0.4D1
      H6x6Re(4,2) = hbar ** 2 / m * gamma2 * kx ** 2 * sqrt(0.3D1) / 0.2
     #D1 - hbar ** 2 / m * gamma2 * ky ** 2 * sqrt(0.3D1) / 0.2D1 - b * 
     #cg(1,1) * sqrt(0.3D1) / 0.2D1 + b * cg(2,2) * sqrt(0.3D1) / 0.2D1
      H6x6Re(4,3) = -hbar ** 2 / m * gamma3 * kx * kz * sqrt(0.3D1) + d 
     #* cg(1,3) + BG * sin(thetaM) * cos(phiM) * sqrt(0.3D1) + C4 * sqrt
     #(0.3D1) * cg(1,2) * ky / 0.2D1 - C4 * sqrt(0.3D1) * cg(1,3) * kz /
     # 0.2D1
      H6x6Re(4,4) = -hbar ** 2 / m * gamma1 * kx ** 2 / 0.2D1 - hbar ** 
     #2 / m * gamma1 * ky ** 2 / 0.2D1 - hbar ** 2 / m * gamma1 * kz ** 
     #2 / 0.2D1 - hbar ** 2 / m * gamma2 * kx ** 2 / 0.2D1 - hbar ** 2 /
     # m * gamma2 * ky ** 2 / 0.2D1 + hbar ** 2 / m * gamma2 * kz ** 2 +
     # b * cg(1,1) / 0.2D1 + b * cg(2,2) / 0.2D1 - b * cg(3,3) - 0.3D1 *
     # BG * cos(thetaM) - 0.3D1 / 0.2D1 * C4 * cg(3,1) * kx + 0.3D1 / 0.
     #2D1 * C4 * cg(3,2) * ky
      H6x6Re(4,5) = hbar ** 2 / m * gamma2 * kx ** 2 * sqrt(0.2D1) * sqr
     #t(0.3D1) / 0.2D1 - hbar ** 2 / m * gamma2 * ky ** 2 * sqrt(0.2D1) 
     #* sqrt(0.3D1) / 0.2D1 - b * cg(1,1) * sqrt(0.2D1) * sqrt(0.3D1) / 
     #0.2D1 + b * cg(2,2) * sqrt(0.2D1) * sqrt(0.3D1) / 0.2D1
      H6x6Re(4,6) = -hbar ** 2 / m * gamma3 * kx * kz * sqrt(0.2D1) * sq
     #rt(0.3D1) / 0.2D1 + d * cg(3,1) * sqrt(0.2D1) / 0.2D1 - BG * sin(t
     #hetaM) * cos(phiM) * sqrt(0.2D1) * sqrt(0.3D1) + C4 * sqrt(0.2D1) 
     #* sqrt(0.3D1) * cg(1,2) * ky / 0.4D1 - C4 * sqrt(0.2D1) * sqrt(0.3
     #D1) * cg(1,3) * kz / 0.4D1 - C4 * eta * sqrt(0.2D1) * sqrt(0.3D1) 
     #* cg(1,2) * ky / 0.8D1 + C4 * eta * sqrt(0.2D1) * sqrt(0.3D1) * cg
     #(1,3) * kz / 0.8D1
      H6x6Re(5,1) = -hbar ** 2 / m * gamma3 * kx * kz * sqrt(0.2D1) * sq
     #rt(0.3D1) / 0.2D1 + d * cg(3,1) * sqrt(0.2D1) / 0.2D1 + BG * sin(t
     #hetaM) * cos(phiM) * sqrt(0.2D1) * sqrt(0.3D1) - C4 * sqrt(0.2D1) 
     #* sqrt(0.3D1) * cg(1,2) * ky / 0.4D1 + C4 * sqrt(0.2D1) * sqrt(0.3
     #D1) * cg(1,3) * kz / 0.4D1 + C4 * eta * sqrt(0.2D1) * sqrt(0.3D1) 
     #* cg(1,2) * ky / 0.8D1 - C4 * eta * sqrt(0.2D1) * sqrt(0.3D1) * cg
     #(1,3) * kz / 0.8D1
      H6x6Re(5,2) = -hbar ** 2 / m * gamma2 * kx ** 2 * sqrt(0.2D1) / 0.
     #2D1 - hbar ** 2 / m * gamma2 * ky ** 2 * sqrt(0.2D1) / 0.2D1 + hba
     #r ** 2 / m * gamma2 * kz ** 2 * sqrt(0.2D1) + b * cg(1,1) * sqrt(0
     #.2D1) / 0.2D1 + b * cg(2,2) * sqrt(0.2D1) / 0.2D1 - b * cg(3,3) * 
     #sqrt(0.2D1) - 0.2D1 * BG * cos(thetaM) * sqrt(0.2D1) + C4 * sqrt(0
     #.2D1) * cg(3,1) * kx / 0.2D1 - C4 * sqrt(0.2D1) * cg(3,2) * ky / 0
     #.2D1 - C4 * sqrt(0.2D1) * eta * cg(3,1) * kx / 0.4D1 + C4 * sqrt(0
     #.2D1) * eta * cg(3,2) * ky / 0.4D1
      H6x6Re(5,3) = 0.3D1 / 0.2D1 * hbar ** 2 / m * gamma3 * kx * kz * s
     #qrt(0.2D1) - d * sqrt(0.3D1) * cg(3,1) * sqrt(0.2D1) / 0.2D1 - BG 
     #* sin(thetaM) * cos(phiM) * sqrt(0.2D1) + C4 * sqrt(0.2D1) * cg(1,
     #2) * ky / 0.4D1 - C4 * sqrt(0.2D1) * cg(1,3) * kz / 0.4D1 - C4 * e
     #ta * sqrt(0.2D1) * cg(1,2) * ky / 0.8D1 + C4 * eta * sqrt(0.2D1) *
     # cg(1,3) * kz / 0.8D1
      H6x6Re(5,4) = hbar ** 2 / m * gamma2 * kx ** 2 * sqrt(0.2D1) * sqr
     #t(0.3D1) / 0.2D1 - hbar ** 2 / m * gamma2 * ky ** 2 * sqrt(0.2D1) 
     #* sqrt(0.3D1) / 0.2D1 - b * cg(1,1) * sqrt(0.2D1) * sqrt(0.3D1) / 
     #0.2D1 + b * cg(2,2) * sqrt(0.2D1) * sqrt(0.3D1) / 0.2D1
      H6x6Re(5,5) = -Lambda - hbar ** 2 / m * gamma1 * kx ** 2 / 0.2D1 -
     # hbar ** 2 / m * gamma1 * ky ** 2 / 0.2D1 - hbar ** 2 / m * gamma1
     # * kz ** 2 / 0.2D1 - BG * cos(thetaM) + C4 * cg(3,1) * kx - C4 * c
     #g(3,2) * ky - C4 * eta * cg(3,1) * kx + C4 * eta * cg(3,2) * ky
      H6x6Re(5,6) = -BG * sin(thetaM) * cos(phiM) + C4 * cg(1,2) * ky - 
     #C4 * cg(1,3) * kz - C4 * eta * cg(1,2) * ky + C4 * eta * cg(1,3) *
     # kz
      H6x6Re(6,1) = -hbar ** 2 / m * gamma2 * kx ** 2 * sqrt(0.2D1) * sq
     #rt(0.3D1) / 0.2D1 + hbar ** 2 / m * gamma2 * ky ** 2 * sqrt(0.2D1)
     # * sqrt(0.3D1) / 0.2D1 + b * cg(1,1) * sqrt(0.2D1) * sqrt(0.3D1) /
     # 0.2D1 - b * cg(2,2) * sqrt(0.2D1) * sqrt(0.3D1) / 0.2D1
      H6x6Re(6,2) = 0.3D1 / 0.2D1 * hbar ** 2 / m * gamma3 * kx * kz * s
     #qrt(0.2D1) - d * sqrt(0.3D1) * cg(3,1) * sqrt(0.2D1) / 0.2D1 + BG 
     #* sin(thetaM) * cos(phiM) * sqrt(0.2D1) - C4 * sqrt(0.2D1) * cg(1,
     #2) * ky / 0.4D1 + C4 * sqrt(0.2D1) * cg(1,3) * kz / 0.4D1 + C4 * e
     #ta * sqrt(0.2D1) * cg(1,2) * ky / 0.8D1 - C4 * eta * sqrt(0.2D1) *
     # cg(1,3) * kz / 0.8D1
      H6x6Re(6,3) = hbar ** 2 / m * gamma2 * kx ** 2 * sqrt(0.2D1) / 0.2
     #D1 + hbar ** 2 / m * gamma2 * ky ** 2 * sqrt(0.2D1) / 0.2D1 - hbar
     # ** 2 / m * gamma2 * kz ** 2 * sqrt(0.2D1) - b * cg(1,1) * sqrt(0.
     #2D1) / 0.2D1 - b * cg(2,2) * sqrt(0.2D1) / 0.2D1 + b * cg(3,3) * s
     #qrt(0.2D1) - 0.2D1 * BG * cos(thetaM) * sqrt(0.2D1) + C4 * sqrt(0.
     #2D1) * cg(3,1) * kx / 0.2D1 - C4 * sqrt(0.2D1) * cg(3,2) * ky / 0.
     #2D1 - C4 * sqrt(0.2D1) * eta * cg(3,1) * kx / 0.4D1 + C4 * sqrt(0.
     #2D1) * eta * cg(3,2) * ky / 0.4D1
      H6x6Re(6,4) = -hbar ** 2 / m * gamma3 * kx * kz * sqrt(0.2D1) * sq
     #rt(0.3D1) / 0.2D1 + d * cg(3,1) * sqrt(0.2D1) / 0.2D1 - BG * sin(t
     #hetaM) * cos(phiM) * sqrt(0.2D1) * sqrt(0.3D1) + C4 * sqrt(0.2D1) 
     #* sqrt(0.3D1) * cg(1,2) * ky / 0.4D1 - C4 * sqrt(0.2D1) * sqrt(0.3
     #D1) * cg(1,3) * kz / 0.4D1 - C4 * eta * sqrt(0.2D1) * sqrt(0.3D1) 
     #* cg(1,2) * ky / 0.8D1 + C4 * eta * sqrt(0.2D1) * sqrt(0.3D1) * cg
     #(1,3) * kz / 0.8D1
      H6x6Re(6,5) = -BG * sin(thetaM) * cos(phiM) + C4 * cg(1,2) * ky - 
     #C4 * cg(1,3) * kz - C4 * eta * cg(1,2) * ky + C4 * eta * cg(1,3) *
     # kz
      H6x6Re(6,6) = -Lambda - hbar ** 2 / m * gamma1 * kx ** 2 / 0.2D1 -
     # hbar ** 2 / m * gamma1 * ky ** 2 / 0.2D1 - hbar ** 2 / m * gamma1
     # * kz ** 2 / 0.2D1 + BG * cos(thetaM) - C4 * cg(3,1) * kx + C4 * c
     #g(3,2) * ky + C4 * eta * cg(3,1) * kx - C4 * eta * cg(3,2) * ky

      H6x6Im(1,2) = -hbar ** 2 / m * gamma3 * ky * kz * sqrt(0.3D1) + d 
     #* cg(2,3) - BG * sin(thetaM) * sin(phiM) * sqrt(0.3D1) - C4 * sqrt
     #(0.3D1) * cg(2,3) * kz / 0.2D1 + C4 * sqrt(0.3D1) * cg(2,1) * kx /
     # 0.2D1
      H6x6Im(1,3) = -hbar ** 2 / m * gamma3 * kx * ky * sqrt(0.3D1) + d 
     #* cg(1,2)
      H6x6Im(1,5) = hbar ** 2 / m * gamma3 * ky * kz * sqrt(0.2D1) * sqr
     #t(0.3D1) / 0.2D1 - d * cg(2,3) * sqrt(0.2D1) / 0.2D1 - BG * sin(th
     #etaM) * sin(phiM) * sqrt(0.2D1) * sqrt(0.3D1) + C4 * sqrt(0.2D1) *
     # sqrt(0.3D1) * cg(2,3) * kz / 0.4D1 - C4 * sqrt(0.2D1) * sqrt(0.3D
     #1) * cg(2,1) * kx / 0.4D1 - C4 * eta * sqrt(0.2D1) * sqrt(0.3D1) *
     # cg(2,3) * kz / 0.8D1 + C4 * eta * sqrt(0.2D1) * sqrt(0.3D1) * cg(
     #2,1) * kx / 0.8D1
      H6x6Im(1,6) = hbar ** 2 / m * gamma3 * kx * ky * sqrt(0.2D1) * sqr
     #t(0.3D1) - d * cg(1,2) * sqrt(0.2D1)
      H6x6Im(2,1) = hbar ** 2 / m * gamma3 * ky * kz * sqrt(0.3D1) - d *
     # cg(2,3) + BG * sin(thetaM) * sin(phiM) * sqrt(0.3D1) + C4 * sqrt(
     #0.3D1) * cg(2,3) * kz / 0.2D1 - C4 * sqrt(0.3D1) * cg(2,1) * kx / 
     #0.2D1
      H6x6Im(2,3) = -0.2D1 * BG * sin(thetaM) * sin(phiM) - C4 * cg(2,3)
     # * kz + C4 * cg(2,1) * kx
      H6x6Im(2,4) = -hbar ** 2 / m * gamma3 * kx * ky * sqrt(0.3D1) + d 
     #* cg(1,2)
      H6x6Im(2,6) = -0.3D1 / 0.2D1 * hbar ** 2 / m * gamma3 * ky * kz * 
     #sqrt(0.2D1) + d * sqrt(0.3D1) * cg(2,3) * sqrt(0.2D1) / 0.2D1 - BG
     # * sin(thetaM) * sin(phiM) * sqrt(0.2D1) + C4 * sqrt(0.2D1) * cg(2
     #,3) * kz / 0.4D1 - C4 * sqrt(0.2D1) * cg(2,1) * kx / 0.4D1 - C4 * 
     #eta * sqrt(0.2D1) * cg(2,3) * kz / 0.8D1 + C4 * eta * sqrt(0.2D1) 
     #* cg(2,1) * kx / 0.8D1
      H6x6Im(3,1) = hbar ** 2 / m * gamma3 * kx * ky * sqrt(0.3D1) - d *
     # cg(1,2)
      H6x6Im(3,2) = 0.2D1 * BG * sin(thetaM) * sin(phiM) + C4 * cg(2,3) 
     #* kz - C4 * cg(2,1) * kx
      H6x6Im(3,4) = hbar ** 2 / m * gamma3 * ky * kz * sqrt(0.3D1) - d *
     # cg(2,3) - BG * sin(thetaM) * sin(phiM) * sqrt(0.3D1) - C4 * sqrt(
     #0.3D1) * cg(2,3) * kz / 0.2D1 + C4 * sqrt(0.3D1) * cg(2,1) * kx / 
     #0.2D1
      H6x6Im(3,5) = 0.3D1 / 0.2D1 * hbar ** 2 / m * gamma3 * ky * kz * s
     #qrt(0.2D1) - d * sqrt(0.3D1) * cg(2,3) * sqrt(0.2D1) / 0.2D1 - BG 
     #* sin(thetaM) * sin(phiM) * sqrt(0.2D1) + C4 * sqrt(0.2D1) * cg(2,
     #3) * kz / 0.4D1 - C4 * sqrt(0.2D1) * cg(2,1) * kx / 0.4D1 - C4 * e
     #ta * sqrt(0.2D1) * cg(2,3) * kz / 0.8D1 + C4 * eta * sqrt(0.2D1) *
     # cg(2,1) * kx / 0.8D1
      H6x6Im(4,2) = hbar ** 2 / m * gamma3 * kx * ky * sqrt(0.3D1) - d *
     # cg(1,2)
      H6x6Im(4,3) = -hbar ** 2 / m * gamma3 * ky * kz * sqrt(0.3D1) + d 
     #* cg(2,3) + BG * sin(thetaM) * sin(phiM) * sqrt(0.3D1) + C4 * sqrt
     #(0.3D1) * cg(2,3) * kz / 0.2D1 - C4 * sqrt(0.3D1) * cg(2,1) * kx /
     # 0.2D1
      H6x6Im(4,5) = hbar ** 2 / m * gamma3 * kx * ky * sqrt(0.2D1) * sqr
     #t(0.3D1) - d * cg(1,2) * sqrt(0.2D1)
      H6x6Im(4,6) = -hbar ** 2 / m * gamma3 * ky * kz * sqrt(0.2D1) * sq
     #rt(0.3D1) / 0.2D1 + d * cg(2,3) * sqrt(0.2D1) / 0.2D1 - BG * sin(t
     #hetaM) * sin(phiM) * sqrt(0.2D1) * sqrt(0.3D1) + C4 * sqrt(0.2D1) 
     #* sqrt(0.3D1) * cg(2,3) * kz / 0.4D1 - C4 * sqrt(0.2D1) * sqrt(0.3
     #D1) * cg(2,1) * kx / 0.4D1 - C4 * eta * sqrt(0.2D1) * sqrt(0.3D1) 
     #* cg(2,3) * kz / 0.8D1 + C4 * eta * sqrt(0.2D1) * sqrt(0.3D1) * cg
     #(2,1) * kx / 0.8D1
      H6x6Im(5,1) = -hbar ** 2 / m * gamma3 * ky * kz * sqrt(0.2D1) * sq
     #rt(0.3D1) / 0.2D1 + d * cg(2,3) * sqrt(0.2D1) / 0.2D1 + BG * sin(t
     #hetaM) * sin(phiM) * sqrt(0.2D1) * sqrt(0.3D1) - C4 * sqrt(0.2D1) 
     #* sqrt(0.3D1) * cg(2,3) * kz / 0.4D1 + C4 * sqrt(0.2D1) * sqrt(0.3
     #D1) * cg(2,1) * kx / 0.4D1 + C4 * eta * sqrt(0.2D1) * sqrt(0.3D1) 
     #* cg(2,3) * kz / 0.8D1 - C4 * eta * sqrt(0.2D1) * sqrt(0.3D1) * cg
     #(2,1) * kx / 0.8D1
      H6x6Im(5,3) = -0.3D1 / 0.2D1 * hbar ** 2 / m * gamma3 * ky * kz * 
     #sqrt(0.2D1) + d * sqrt(0.3D1) * cg(2,3) * sqrt(0.2D1) / 0.2D1 + BG
     # * sin(thetaM) * sin(phiM) * sqrt(0.2D1) - C4 * sqrt(0.2D1) * cg(2
     #,3) * kz / 0.4D1 + C4 * sqrt(0.2D1) * cg(2,1) * kx / 0.4D1 + C4 * 
     #eta * sqrt(0.2D1) * cg(2,3) * kz / 0.8D1 - C4 * eta * sqrt(0.2D1) 
     #* cg(2,1) * kx / 0.8D1
      H6x6Im(5,4) = -hbar ** 2 / m * gamma3 * kx * ky * sqrt(0.2D1) * sq
     #rt(0.3D1) + d * cg(1,2) * sqrt(0.2D1)
      H6x6Im(5,6) = BG * sin(thetaM) * sin(phiM) - C4 * cg(2,3) * kz + C
     #4 * cg(2,1) * kx + C4 * eta * cg(2,3) * kz - C4 * eta * cg(2,1) * 
     #kx
      H6x6Im(6,1) = -hbar ** 2 / m * gamma3 * kx * ky * sqrt(0.2D1) * sq
     #rt(0.3D1) + d * cg(1,2) * sqrt(0.2D1)
      H6x6Im(6,2) = 0.3D1 / 0.2D1 * hbar ** 2 / m * gamma3 * ky * kz * s
     #qrt(0.2D1) - d * sqrt(0.3D1) * cg(2,3) * sqrt(0.2D1) / 0.2D1 + BG 
     #* sin(thetaM) * sin(phiM) * sqrt(0.2D1) - C4 * sqrt(0.2D1) * cg(2,
     #3) * kz / 0.4D1 + C4 * sqrt(0.2D1) * cg(2,1) * kx / 0.4D1 + C4 * e
     #ta * sqrt(0.2D1) * cg(2,3) * kz / 0.8D1 - C4 * eta * sqrt(0.2D1) *
     # cg(2,1) * kx / 0.8D1
      H6x6Im(6,4) = hbar ** 2 / m * gamma3 * ky * kz * sqrt(0.2D1) * sqr
     #t(0.3D1) / 0.2D1 - d * cg(2,3) * sqrt(0.2D1) / 0.2D1 + BG * sin(th
     #etaM) * sin(phiM) * sqrt(0.2D1) * sqrt(0.3D1) - C4 * sqrt(0.2D1) *
     # sqrt(0.3D1) * cg(2,3) * kz / 0.4D1 + C4 * sqrt(0.2D1) * sqrt(0.3D
     #1) * cg(2,1) * kx / 0.4D1 + C4 * eta * sqrt(0.2D1) * sqrt(0.3D1) *
     # cg(2,3) * kz / 0.8D1 - C4 * eta * sqrt(0.2D1) * sqrt(0.3D1) * cg(
     #2,1) * kx / 0.8D1
      H6x6Im(6,5) = -BG * sin(thetaM) * sin(phiM) + C4 * cg(2,3) * kz - 
     #C4 * cg(2,1) * kx - C4 * eta * cg(2,3) * kz + C4 * eta * cg(2,1) *
     # kx

      do i = 1, 6
         do j = 1, 6
            H6x6(i, j) = dcmplx(H6x6Re(i, j), H6x6Im(i, j))
         end do
      end do

      end

