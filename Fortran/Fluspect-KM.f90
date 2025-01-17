!***************************************************************************************
! This code is updated from the Fluspect model (MATLAB version).
!   https://github.com/Christiaanvandertol/Fluspect
!
! An analytical algorithm is applied to accelerate Fluspect.
!
! Author: Hanyu Shi (shihanyu@mail.bnu.edu.cn).
!
! Reference:
! [1] Shi and Xiao (2025). Accelerating Fluspect With an Analytical Algorithm in
!       Simulating Mesophyll Fluorescence Matrices.
!       IEEE Geosci. Remote Sens. Lett. 22. doi:10.1109/LGRS.2025.3529029.
! [2] Vilfan et al. (2016). Fluspect-B: A model for leaf fluorescence, reflectance
!       and transmittance spectra.
!       Remote Sens. Environ. 186, 596–615. doi:10.1016/j.rse.2016.09.017
! [3] Vilfan et al. (2018). Extending Fluspect to simulate xanthophyll driven leaf
!       reflectance dynamics.
!       Remote Sens. Environ. 211, 345–356. doi:10.1016/j.rse.2018.04.012
!***************************************************************************************

module func_Fluspect_KM
    use fluspect2021_dataSpec, only : nw,KcaV,KcaZ,phi,k_Car,k_Car,k_Cab,k_Cm,k_Cw, &
        k_Brown,k_Anth,refractive,phiI,phiII,lambda,phi,k_CBC,k_Prot
    implicit none

    integer, private :: iloop
    integer, parameter :: nwlf=211, nwle=351
    integer, parameter :: Iwle_s=1, Iwle_e=351
    integer, parameter :: Iwlf_s=241, Iwlf_e=451
    real(8), parameter :: wlf(nwlf) = [(iloop*1.,iloop=640,850)]
    real(8), parameter :: wle(nwle) = [(iloop*1.,iloop=400,750)]
    real(8), parameter :: Ih(nwle) = 1.d0, Iv(nwlf) = 1.d0
contains

    ! This subroutine uses phiI and phiII.
    !   A combination of them (phi) can be also easily implemented.
    subroutine Fluspect_KM(Cab,Cca,V2Z,Cw,Cdm,Cs,Cant,Cbc,Cp,N,fqeI,fqeII,RT,FLUO_I,FLUO_II)
        implicit none

        real(8), intent(in) :: Cab,Cca,V2Z,Cw,Cdm,Cs,Cant,Cbc,Cp,N,fqeI,fqeII
        real(8), intent(out) :: RT(nw,10),FLUO_I(nwlf,nwle,2),FLUO_II(nwlf,nwle,2)

        real(8) :: kca(nw)
        real(8) :: tau(nw),rho(nw)
        real(8) :: r21(nw),t21(nw),ralf(nw),talf(nw),r12(nw),t12(nw)
        real(8) :: kChlrel(nw),kChl(nw),s(nw),k(nw)
        real(8) :: f(nwlf,nwle),g(nwlf,nwle),fn(nwlf,nwle),gn(nwlf,nwle)

        real(8) :: Kabs(nw)
        real(8) :: refl(nw),tran(nw)

        real(8), parameter :: ang = 59.d0


        if (V2Z < 0.d0) then
            Kca = k_Car
        else
            Kca = (1.d0-V2Z)*KcaV + V2Z*KcaZ
        endif

        Kabs = (Cab*k_Cab+Cca*Kca+Cdm*k_Cm+Cw*k_Cw+Cs*k_Brown+Cant*k_Anth+Cbc*k_CBC+Cp*k_Prot) / N


        kChlrel = 0.d0
        where (Kabs > 0)
            kChlrel = Cab*k_Cab / (Kabs*N)
        end where

        call M1(ang,N,refractive,Kabs,refl,tran,r21,t21,ralf,talf,r12,t12)
        call M2(refl,tran,r21,t21,ralf,talf,rho,tau,s,k)


        RT(:,2) = tran
        RT(:,1) = refl

        kChl = kChlrel * k


        call M3(fqeI,phiI,kChl,rho,tau,s,k,f,g)
        call M4(rho,tau,r21,t21,talf,g,f,gn,fn)
        FLUO_I(:,:,1) = fn  ! Forward, PS-I
        FLUO_I(:,:,2) = gn  ! Backward, PS-II


        call M3(fqeII,phiII,kChl,rho,tau,s,k,f,g)
        call M4(rho,tau,r21,t21,talf,g,f,gn,fn)
        FLUO_II(:,:,1) = fn  ! Forward, PS-I
        FLUO_II(:,:,2) = gn  ! Backward, PS-II

        return
    end subroutine


    subroutine M1(ang,N,rnnl,kabs,refl,tran,r21,t21,ralf,talf,r12,t12)
        implicit none
        real(8), intent(in) :: ang,N,rnnl(nw),kabs(nw)
        real(8), intent(out) :: r21(nw),t21(nw),ralf(nw),talf(nw),r12(nw),t12(nw)
        real(8), intent(out) :: refl(nw),tran(nw)

        real(8) :: tau(nw),r(nw),t(nw)
        real(8) :: denom(nw),Ra(nw),Ta(nw)
        real(8) :: d(nw),rq(nw),tq(nw),a(nw),b(nw)
        real(8) :: bNm1(nw),bN2(nw),a2(nw),Rsub(nw),Tsub(nw)


        call S13AAF(nw,kabs,tau)

        call tav_abs(90.d0,rnnl,t12)
        call tav_abs(ang,rnnl,talf)

        ralf    = 1.-talf
        r12     = 1.-t12
        t21     = t12/(rnnl**2)
        r21     = 1-t21
        ! top surface side
        denom   = 1-r21*r21*tau**2
        Ta      = talf*tau*t21/denom
        Ra      = ralf+r21*tau*Ta
        ! bottom surface side
        t       = t12*tau*t21/denom
        r       = r12+r21*tau*t

        ! Stokes
        D       = sqrt((1.+r+t)*(1.+r-t)*(1.-r+t)*(1.-r-t))
        rq      = r**2
        tq      = t**2
        a       = (1.+rq-tq+D)/(2*r)
        b       = (1.-rq+tq+D)/(2*t)

        bNm1    = b**(N-1)
        bN2     = bNm1**2
        a2      = a**2
        denom   = a2*bN2-1.
        Rsub    = a*(bN2-1.)/denom
        Tsub    = bNm1*(a2-1.)/denom

        ! Case of zero absorption
        where ((r+t) .ge. 1.0)
            Tsub    = t/(t+(1.-t)*(N-1))
            Rsub    = 1-Tsub
        end where

        ! Reflectance and transmittance of the leaf: combine top layer with next N-1 layers
        denom   = 1-Rsub*r
        tran    = Ta*Tsub/denom
        refl    = Ra+Ta*Rsub*t/denom
    end subroutine

    subroutine M2(refl,tran,r21,t21,ralf,talf,r,t,s,k)
        implicit none
        real(8), intent(in) :: refl(nw),tran(nw),r21(nw),t21(nw),ralf(nw),talf(nw)
        real(8), intent(out) :: r(nw),t(nw),s(nw),k(nw)

        real(8) :: Z(nw),Rb(nw),D(nw),a(nw),b(nw)


        Rb = (refl-ralf)/(talf*t21+(refl-ralf)*r21)
        Z = tran*(1-Rb*r21)/(talf*t21)
        r = (Rb-r21*Z**2) / (1-(r21*Z)**2)  ! rho
        t = (1-Rb*r21) / (1-(r21*Z)**2) * Z ! tau

        where (r < 0.d0)
            r = 0.d0
        endwhere

        where ((r+t) < 1.)
            D = sqrt((1.+r+t)*(1.+r-t)*(1.-r+t)*(1.-r-t))
            a = (1. + r**2 - t**2 + D) / (2.*r)
            b = (1. - r**2 + t**2 + D) / (2.*t)
        elsewhere
            a = 1.
            b = 1.
        endwhere

        s = r / t
        k = dlog(b)
        where ((a>1))
            s = 2.*a / (a**2-1.) * dlog(b)
            k = (a-1.) / (a+1.) * dlog(b)
        endwhere

    end subroutine

    subroutine M3(fqe,PS_phi,kChl,rho,tau,s,k,Mf,Mb)
        implicit none
        real(8), intent(in) :: fqe,PS_phi(nw),kChl(nw),rho(nw),tau(nw),s(nw),k(nw)
        real(8), intent(out) :: Mf(nwlf,nwle),Mb(nwlf,nwle)

        real(8) :: ae(nwle),se(nwle),mme(nwle),taue(nwle)
        real(8) :: af(nwlf),sf(nwlf),mmf(nwlf),rhof(nwlf),tauf(nwlf)
        real(8) :: A1(nwlf,nwle),A2(nwlf,nwle),A3(nwlf,nwle),A4(nwlf,nwle)
        real(8) :: H0(nwlf,nwle),H1(nwlf,nwle),H2(nwlf,nwle),H3(nwlf,nwle),H4(nwlf,nwle)

        real(8) :: tmp1(nwlf,nwle),tmp2(nwlf,nwle)

        real(8) :: mempe(nwlf,nwle),meppe(nwlf,nwle),mempf(nwlf,nwle),meppf(nwlf,nwle)
        real(8) :: WWf(nwlf),WWg(nwlf)
        real(8) :: sigmoid(nwlf,nwle),varphi(nwlf,nwle)

        integer :: i


        sigmoid = 1.d0 / (1.d0+vec2mat_matmul(nwlf,nwle,dexp(-wlf/10.d0),dexp(wle/10.d0)))
        varphi = fqe*vec2mat_matmul(nwlf,nwle,(.5*PS_phi(Iwlf_s:Iwlf_e)),kChl(Iwle_s:Iwle_e))*sigmoid


        se = s(Iwle_s:Iwle_e)
        sf = s(Iwlf_s:Iwlf_e)
        ae = k(Iwle_s:Iwle_e)+se
        af = k(Iwlf_s:Iwlf_e)+sf
        mme = dsqrt(ae**2 - se**2)
        mmf = dsqrt(af**2 - sf**2)
        tauf = tau(Iwlf_s:Iwlf_e)
        rhof = rho(Iwlf_s:Iwlf_e)
        taue = tau(Iwle_s:Iwle_e)

        mempe = vec2mat_matmul(nwlf,nwle,Iv,mme-(ae+se))
        meppe = vec2mat_matmul(nwlf,nwle,Iv,mme+(ae+se))
        mempf = vec2mat_matmul(nwlf,nwle,Iv,mme) - vec2mat_matmul(nwlf,nwle,af+sf,Ih)
        meppf = vec2mat_matmul(nwlf,nwle,Iv,mme) + vec2mat_matmul(nwlf,nwle,af+sf,Ih)

        A1 = -mempe*mempf
        A2 = mempe*meppf
        A3 = -meppe*meppf
        A4 = meppe*mempf

        tmp1 = vec2mat_matmul(nwlf,nwle,Iv,exp(mme))
        tmp2 = vec2mat_matmul(nwlf,nwle,Iv,exp(-mme))
        H1 = A3*tmp1 - A1*tmp2
        H2 = A4 - A2
        H3 = A3 - A1
        H4 = A4*tmp1 - A2*tmp2


        tmp1 = vec2mat_matmul(nwlf,nwle,Iv,mme**2) - vec2mat_matmul(nwlf,nwle,mmf**2,Ih)
        do i = 1, 111
            tmp1(i,i+240) = 1.d0
        enddo
        H0 = 0.5d0 * vec2mat_matmul(nwlf,nwle,Iv,taue/mme) / tmp1


        tmp1 = vec2mat_matmul(nwlf,nwle,tauf,Ih)
        tmp2 = vec2mat_matmul(nwlf,nwle,rhof,Ih)
        Mf = H0*(H3 - H1*tmp1 - H2*tmp2)
        Mb = H0*(H4 - H1*tmp2 - H2*tmp1)


        WWf = tauf * (af - rhof*(1.d0+sf)) / k(Iwlf_s:Iwlf_e)
        WWg = 0.5d0 * (1.d0 - rhof**2 - tauf**2*(1.d0+2.d0*sf)) / k(Iwlf_s:Iwlf_e)


        do i = 1, 111
            Mf(i,i+240) = WWf(i)
            Mb(i,i+240) = WWg(i)
        enddo

        Mf = Mf * varphi
        Mb = Mb * varphi

    end subroutine

    subroutine M4(rho,tau,r21,t21,talf,g,f,gn,fn)
        implicit none
        real(8), intent(in) :: rho(nw),tau(nw),r21(nw),t21(nw),talf(nw),f(nwlf,nwle),g(nwlf,nwle)
        real(8), intent(out) :: fn(nwlf,nwle),gn(nwlf,nwle)

        real(8) :: Rb(nw)
        real(8) :: XeXe(nwlf,nwle),XfXf(nwlf,nwle),YeYe(nwlf,nwle),YfYf(nwlf,nwle)
        real(8) :: tmp1(nwle),tmp2(nwlf),AA(nwlf,nwle),BB(nwlf,nwle)

        Rb = rho + tau**2*r21/(1.d0-rho*r21)

        tmp1 = talf(Iwle_s:Iwle_e) / (1.d0-r21(Iwle_s:Iwle_e)*Rb(Iwle_s:Iwle_e))
        XeXe = vec2mat_matmul(nwlf,nwle,Iv,tmp1)
        tmp2 = t21(Iwlf_s:Iwlf_e) / (1.d0-r21(Iwlf_s:Iwlf_e)*Rb(Iwlf_s:Iwlf_e))
        XfXf = vec2mat_matmul(nwlf,nwle,tmp2,Ih)

        tmp1 = tau(Iwle_s:Iwle_e)*r21(Iwle_s:Iwle_e) / (1.d0-rho(Iwle_s:Iwle_e)*r21(Iwle_s:Iwle_e))
        YeYe = vec2mat_matmul(nwlf,nwle,Iv,tmp1)
        tmp2 = tau(Iwlf_s:Iwlf_e)*r21(Iwlf_s:Iwlf_e) / (1.d0-rho(Iwlf_s:Iwlf_e)*r21(Iwlf_s:Iwlf_e))
        YfYf = vec2mat_matmul(nwlf,nwle,tmp2,Ih)

        AA = XeXe * (1.+YeYe*YfYf)*XfXf
        BB = XeXe * (YeYe+YfYf)*XfXf

        gn = AA * g + BB * f
        fn = AA * f + BB * g
    end subroutine


    subroutine tav_abs(theta,nr,tav)
        implicit none

        real(8), intent(in) :: theta, nr(nw)
        real(8), intent(out) :: tav(nw)

        real(8) :: pi,rd
        real(8) :: n2(nw),np(nw),nm(nw)
        real(8) :: a(nw),k(nw),sa(nw),b1(nw),b2(nw),b3(nw),b(nw),a3(nw)
        real(8) :: ts(nw),tp(nw),tp1(nw),tp2(nw),tp3(nw),tp4(nw),tp5(nw)


        pi  = atan(1.)*4.
        rd  = pi/180.
        n2  = nr**2.
        np  = n2+1.
        nm  = n2-1.
        a   = (nr+1)*(nr+1.)/2.
        k   = -(n2-1)*(n2-1.)/4.
        sa  = sin(theta*rd)

        if (theta.eq.90.) then
            b1=0.
        else
            b1  = dsqrt((sa**2-np/2)*(sa**2-np/2)+k)
        endif

        b2  = sa**2-np/2
        b   = b1-b2
        b3  = b**3
        a3  = a**3
        ts  = (k**2./(6*b3)+k/b-b/2)-(k**2./(6*a3)+k/a-a/2)

        tp1 = -2*n2*(b-a)/(np**2)
        tp2 = -2*n2*np*dlog(b/a)/(nm**2)
        tp3 = n2*(1./b-1./a)/2
        tp4 = 16*n2**2.*(n2**2+1)*dlog((2*np*b-nm**2)/(2*np*a-nm**2))/(np**3.*nm**2)
        tp5 = 16*n2**3.*(1./(2*np*b-nm**2)-1./(2*np*a-nm**2))/(np**3)
        tp  = tp1+tp2+tp3+tp4+tp5
        tav = (ts+tp)/(2*sa**2)

        return
    end subroutine

    subroutine S13AAF(nw,k,tau)
        ! exponential integral: S13AAF routine from the NAG library
        implicit none
        integer, intent(in) :: nw
        real(8), intent(in) :: k(nw)
        real(8), intent(out) :: tau(nw)

        real(8) :: xx(nw),yy(nw)

        where (k.le.0.0)
            tau = 1
        end where
        where (k.gt.0.0.and.k.le.4.0)
            xx  = 0.5*k-1.0
            yy  = (((((((((((((((-3.60311230482612224d-13 &
                *xx+3.46348526554087424d-12)*xx-2.99627399604128973d-11) &
                *xx+2.57747807106988589d-10)*xx-2.09330568435488303d-9) &
                *xx+1.59501329936987818d-8)*xx-1.13717900285428895d-7) &
                *xx+7.55292885309152956d-7)*xx-4.64980751480619431d-6) &
                *xx+2.63830365675408129d-5)*xx-1.37089870978830576d-4) &
                *xx+6.47686503728103400d-4)*xx-2.76060141343627983d-3) &
                *xx+1.05306034687449505d-2)*xx-3.57191348753631956d-2) &
                *xx+1.07774527938978692d-1)*xx-2.96997075145080963d-1
            yy  = (yy*xx+8.64664716763387311d-1)*xx+7.42047691268006429d-1
            yy  = yy-log(k)
            tau = (1.0-k)*dexp(-k)+k**2*yy
        end where
        where (k.gt.4.0.and.k.le.85.0)
            xx  = 14.5/(k+3.25)-1.0
            yy  = (((((((((((((((-1.62806570868460749d-12 &
                *xx-8.95400579318284288d-13)*xx-4.08352702838151578d-12) &
                *xx-1.45132988248537498d-11)*xx-8.35086918940757852d-11) &
                *xx-2.13638678953766289d-10)*xx-1.10302431467069770d-9) &
                *xx-3.67128915633455484d-9)*xx-1.66980544304104726d-8) &
                *xx-6.11774386401295125d-8)*xx-2.70306163610271497d-7) &
                *xx-1.05565006992891261d-6)*xx-4.72090467203711484d-6) &
                *xx-1.95076375089955937d-5)*xx-9.16450482931221453d-5) &
                *xx-4.05892130452128677d-4)*xx-2.14213055000334718d-3
            yy  = ((yy*xx-1.06374875116569657d-2)*xx-8.50699154984571871d-2)*xx+9.23755307807784058d-1
            yy  = exp(-k)*yy/k
            tau = (1.0-k)*dexp(-k)+k**2*yy
        end where
        where (k.gt.85.0)
            tau = 0
        end where
    end subroutine

    function vec2mat_matmul(N1,N2,V1,V2)
        implicit none
        integer, intent(in) :: N1,N2
        real(8), intent(in) :: V1(N1),V2(N2)
        real(8) :: vec2mat_matmul(N1,N2)
        real(8) :: shapeV1(N1,1),shapeV2(1,N2)

        shapeV1 = reshape(V1,(/N1,1/))
        shapeV2 = reshape(V2,(/1,N2/))
        vec2mat_matmul = matmul(shapeV1, shapeV2)

        return
    end function

end module
