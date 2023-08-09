module trajectory
    use kind_module
    implicit none

    integer, parameter :: mpnt = 100000


    integer nrefj(mpnt)
    !!common/refl/nrefj(mpnt)

    integer mbeg(mpnt),mend(mpnt),mbad(mpnt)
    real(wp) rbeg(mpnt) !sav2008
    real(wp) tetbeg(mpnt),xnrbeg(mpnt),xmbeg(mpnt),yn3beg(mpnt)
    !!common/viewdat/mbeg,mend,mbad,rbeg,tetbeg,xnrbeg,xmbeg,yn3beg   
    !data mbad /mpnt*0/
contains

subroutine init_trajectory
    use constants
    use driver_module
    implicit none
    nrefj = 0
    
    dland = zero
    dcoll = zero
    perpn = zero 
    dalf  = zero
    vel = zero
    jrad = zero
    iww = zero
    tetai = zero
    xnpar = zero
    izz = zero

    mbeg = zero
    mend = zero
    mbad = zero

    rbeg = zero

    tetbeg = zero
    xnrbeg = zero
    xmbeg = zero
    yn3beg = zero
end subroutine 

subroutine view(tview, ispectr,nnz,ntet) !sav2008
!!!writing trajectories into a file
    use constants
    use approximation
    use plasma
    use dispersion_module, only: zatukh
    use rt_parameters, only :  nr, itend0, kv, nmaxm    
    use spectrum1D, only: ynzm, pm 
    use dispersion_module
    use driver_module !, only: jrad, iww, izz, length
    implicit real(wp) (a-h,o-z)
    real(wp), intent(in) :: tview

    integer, intent(in) :: ispectr, nnz, ntet  !sav#
    !common /bcef/ ynz,ynpopq
    !common /vth/ vthc(length),poloidn(length)
    real(wp) vthcg,npoli
    !common /a0ghp/ vlf,vrt,dflf,dfrt
    
    integer i, n, itr, ntraj
    integer jrc,nturn,ib,ie,jr,ifast,idir,iv
    integer jznak,jdlt,mn,mm,jchek,itet,inz
    integer, parameter :: unit_bias = 10
    integer, parameter :: m=7
    real(wp), parameter :: pleft=1.d-10 !m may be chaged together with name(m)
    real(wp) ptt(m),pll(m),pcc(m),paa(m)
    real(wp) pt_c(m),pl_c(m),pc_c(m),pa_c(m)
    character(40) fname
    character*40 name(m)
    save name !sav#
    data name/'lhcd/out/1.dat','lhcd/out/2.dat','lhcd/out/3.dat','lhcd/out/4.dat','lhcd/out/5.dat','lhcd/out/rest.dat','lhcd/out/traj.dat'/
    !if(iview.eq.0) return
    print *, 'view_time=',tview
    print *, name(m)
    if (ispectr>0) then
        write(fname,'("lhcd/traj/pos/", f9.7,".dat")') tview
    else
        write(fname,'("lhcd/traj/neg/", f9.7,".dat")') tview
    endif
    print *, fname
    name(m) = fname
    print *, name(m)

    htet=zero
    h=1d0/dble(nr+1)
    if(ntet.ne.1) htet=(tet2-tet1)/(ntet-1)
    open(1,file='lhcd/out/lcms.dat')
    write(1,*)'     R(m)            Z(m)'
    write(1,*)
    xr=1.d0
    xdl=fdf(xr,cdl,ncoef,xdlp)
    xly=fdf(xr,cly,ncoef,xlyp)
    xgm=fdf(xr,cgm,ncoef,xgmp)
    do i=1,101
        th=dble(i-1)*pi2/dble(100)
        cotet=dcos(th)
        sitet=dsin(th)
        xx=-xdl+xr*cotet-xgm*sitet**2
        zz=xr*xly*sitet
        x=(r0+rm*xx)/1d2
        z=(z0+rm*zz)/1d2
        write(1,5) x,z
    end do
    close(1)
    open(1,file='lhcd/out/npar_crit.dat')
    write(1,*)'  Npar_crit=sqrt(50/Te(keV))'
    write(1,*)
    write(1,*)'   rho         Npar_strong absorption'
    write(1,*)
    do i=1,101
        xr=dble(i-1)/dble(100)
        tmp=ft(xr)/0.16d-8  !Te,  KeV
        parn_c=dsqrt(50d0/tmp)
        write(1,5) xr,parn_c
    end do
    close(1)

    do n=1,m
        open(n+unit_bias,file=name(n))
        write(n+unit_bias,*)
        close(n+unit_bias)
        open(n+unit_bias,file=name(n))
        write(n+unit_bias,3)
        write(n+unit_bias,*)
        ptt(n)=zero
        pll(n)=zero
        pcc(n)=zero
        paa(n)=zero
    end do
    ntraj=0 !sav2008
    do itr=1,nnz*ntet !sav2008
        pow=1.d0
        pl=zero
        pc=zero
        pa=zero
        pdec1=zero
        pdec1z=zero
        pdec3=zero
        pdec3z=zero
        pdecv=zero
        pintld=zero
        pintal=zero
        jrc=nr+1
        jznak=-1
        nturn=1
        if(mbad(itr).eq.0) then 
            ntraj=ntraj+1
            ib=mbeg(itr)
            ie=mend(itr)
10          continue
            do i=ib,ie
                v=vel(i)
                jr=jrad(i)
                refr=perpn(i)
                npoli=poloidn(i)
                ifast=iww(i)
                vthcg=vthc(i)
                idir=izz(i)
                dek3=zero
                th=tetai(i)
                parn=xnpar(i)
!         xn1=an1(i)
!         xn2=an2(i)
                if(itend0.gt.0) then
                    argum=clt/(refr*valfa)
                    dek3=zatukh(argum,abs(jr),vperp,kv)
                end if
                ! old variant
                ! call raspr(v,abs(jr),iv,df)
                !
                call distr(v,abs(jr),iv,df)
                if(jr.lt.0) then    !case of turn
                    jr=-jr
!variant          pintld=-dland(i)*df
!!          pintld=-dland(i)*(dflf+dfrt)/2d0
                    pintld=dabs(dland(i)*(dflf+dfrt)/2d0)
                    pdec2=dexp(-2d0*dcoll(i))
                    pintal=dabs(dalf(i)*dek3)
                else
                    pdec2=dcoll(i)
                    pdecv=dland(i)
!!          pdec1=-pdecv*df
                    pdec1=dabs(pdecv*df)
                    pdec3=dabs(dalf(i)*dek3)
                    pintld=(pdec1+pdec1z)/2d0*h
                    pintal=(pdec3+pdec3z)/2d0*h
                    pdec1z=pdec1
                    pdec3z=pdec3
                end if
                powpr=pow
                powd=pow*dexp(-2d0*pintld)
                powcol=powd*pdec2
                powal=powcol*dexp(-2d0*pintal)
                pow=powal
                pil=pintld
                pic=.5d0*dabs(dlog(pdec2))
                pia=pintal
                pt=1.d0-pow  !total absorbed power
                denom=pil+pic+pia
                powdamped=1.d0-dexp(-2.d0*denom)
                domin=powpr*powdamped
                if(denom.ne.zero) then
                    fff=domin/denom
                    pl=pl+dabs(pil*fff)  !el. Landau absorbed power
                    pc=pc+dabs(pic*fff)  !el. collisions absorbed power
                    pa=pa+dabs(pia*fff)  !alpha Landau absorbed power
                end if
                xr=h*dble(jr)
                cotet=dcos(th)
                sitet=dsin(th)
                xdl=fdf(xr,cdl,ncoef,xdlp)
                xly=fdf(xr,cly,ncoef,xlyp)
                xgm=fdf(xr,cgm,ncoef,xgmp)
                xx=-xdl+xr*cotet-xgm*sitet**2
                zz=xr*xly*sitet
                x=(r0+rm*xx)/1d2
                z=(z0+rm*zz)/1d2 
                jdlt=jr-jrc
                jrc=jr
                if(jdlt*jznak.lt.0.and.nturn.lt.m-1) then
                    nturn=nturn+1
                    jznak=-jznak
                end if
                mn = nturn + unit_bias
                write(mn, 7) x,z,xr,th,parn,npoli,pt,pl,pc,pa,ifast,idir,itr
                mm = m + unit_bias
                write(mm, 7) x,z,xr,th,parn,npoli,pt,pl,pc,vthcg,ifast,idir,itr
                do n=m,nturn,-1
                    pt_c(n)=pt
                    pl_c(n)=pl
                    pc_c(n)=pc
                    pa_c(n)=pa
                end do
                if(pt.ge.1d0-pleft) go to 11 !maximal absorbed power along a ray
            end do
            jchek=jrad(ie+1)
            if(jchek.ne.0) then  !continue this trajectory
                ib=idnint(dland(ie+1))
                ie=idnint(dcoll(ie+1))
                goto 10
            end if
11          continue
            do n=1,m
                ptt(n)=ptt(n)+pt_c(n)
                pll(n)=pll(n)+pl_c(n)
                pcc(n)=pcc(n)+pc_c(n)
                paa(n)=paa(n)+pa_c(n)
            end do
            if(itr.lt.nnz*ntet) then
                write(m+unit_bias,*)
                write(nturn+unit_bias,*)
            end if
        end if
    end do
    do n=1,m
        close(n+unit_bias)
    end do
    do n=1,m
        ptt(n)=ptt(n)/dble(ntraj)
        pll(n)=pll(n)/dble(ntraj)
        pcc(n)=pcc(n)/dble(ntraj)
        paa(n)=paa(n)/dble(ntraj)
    end do
    open(1,file='lhcd/out/info_traj.dat')
    write(1,20) tview
    write(1,*)
    do n=1,m-1
        if(n.lt.m-1) then
            write(1,8) n,ptt(n),pll(n),pcc(n),paa(n)
        else
            write(1,9) ptt(n),pll(n),pcc(n),paa(n)
        end if
        write(1,*)
    end do
    write(1,*)
    write(1,1)
    write(1,*)
    itr=0
    do itet=1,ntet
        tetin=tet1+htet*(itet-1)
        do inz=1,nnz
            itr=itr+1
            write(1,6) itr,mbad(itr),tetin,ynzm(inz),rbeg(itr)
        end do
    end do
    close(1)
    open(1,file='lhcd/out/absorp.dat')
    write(1,2)
    write(1,*)
    do n=1,m-1
        if(n.eq.1) then
            dpt=ptt(n)
            dpl=pll(n)
            dpc=pcc(n)
            dpa=paa(n)
        else
            dpt=ptt(n)-ptt(n-1)
            dpl=pll(n)-pll(n-1)
            dpc=pcc(n)-pcc(n-1)
            dpa=paa(n)-paa(n-1)
        end if
        write(1,4) n,ptt(n),pll(n),pcc(n),paa(n),dpt,dpl,dpc,dpa
    end do
    close(1)

1     format(2x,'N_traj',3x,'mbad',6x,'theta',9x,'Npar',9x,'rho_start')
2     format('R_pass',4x,'Ptot',6x,'Pland',6x,'Pcoll',8x,'Pa',7x,'dPtot',6x,'dPland',5x,'dPcoll',6x,'dPa')
3     format(5x,'R',10x,'Z',11x,'rho',8x,'theta',7x,'N_par',7x,'N_pol',6x,'P_tot',7x,'P_land',6x,'P_coll',6x,'vth',4x,'slow=1',4x,'out=1',2x,'N_traj',6x)
4     format(i3,5x,8(f6.3,5x))
5     format(6(e13.6,3x))
6     format(2(i6,2x),4(e13.6,1x))
7     format(10(e11.4,1x),i5,2x,i5,2x,i5)
8     format('after radial pass=',i3,2x,' P_tot=',f6.3,2x,' P_land=',f5.3,2x,' P_coll=',f6.3,2x,' P_a=',f6.3)
9     format('Total passes:           P_tot=',f6.3,2x,' P_land=',f5.3,2x,' P_coll=',f6.3,2x,' P_a=',f6.3)
20    format('written time slice (seconds) =',f9.3)
    end

    subroutine traj(xm0, tet0, xbeg, nmax, nb1, nb2, nomth, nomnz, pabs) !sav2009
        use constants
        use approximation
        use plasma
        use rt_parameters
        use dispersion_module
        use runge_kutta_module
        use driver_module
        real(wp), intent(in)    :: xm0
        real(wp), intent(in)    :: tet0
        real(wp), intent(inout) :: xbeg
        real(wp), intent(in)    :: pabs
        integer,  intent(inout) :: nmax        
        integer,  intent(inout) :: nb1, nb2        
        integer,  intent(in)    :: nomth, nomnz

        !use manager_mod  !, only: ivar, iroot
        !implicit real*8 (a-h,o-z)
        !external extd4


        !common /abc/ rzz,tetzz,xmzz,iznzz,iwzz,irszz
        !common /abcd/ irs
        !common /abcde/ izn!,iw
        !common /abcdg/ iabsorp
        !common /bcg/ hrad
        !common /bcef/ ynz,ynpopq
        !common /be1/ xnr1,xnr2,xnr3,xnr4
        !common /be2/ ider
        !common /bg/ im4

        integer :: nrefl
        integer :: irep
        integer :: irf, irf1
        integer :: ib2
        integer :: irs0        

        real(wp), parameter :: pgdop=0.02d0
        real(wp), parameter :: hmin=0.d-7 !sav2008, old hmin=1.d-7
        real(wp) :: eps0
        real(wp) :: rrange0, hdrob0, tet, xm, hr
        real(wp) :: xsav, xend,hsav, h1
        real(wp) :: ystart(2),yy(4)
        real(wp) :: rnew2
        real(wp) :: cotet, sitet
        real(wp) :: xdl, xly, xgm, xgmp, xx, zz, xxx, zzz
        real(wp) :: xdlp, xlyp, xnr, prt, prm
        real(wp) :: ynz0, x1, x2, rexi, tetnew
        real(wp) :: xmnew, rnew, xnrnew, xnrv
        real(wp) :: pg1, pg2, pg3, pg4, pg

        eps0=eps
        rrange0=rrange
        hdrob0=hdrob
        iroot=1
        nrefl=0
        ider=1
        im4=0
        nb1=0
        nb2=0
        irep=0
        tet=tet0
        xm=xm0
        hr=1.d0/dble(nr+1) !sav2008
        hrad=hr
        !---------------------------------------
        ! find saving point and define
        ! parameters, depending on direction
        !---------------------------------------
  
  10    irf1=idnint(xbeg/hr)
        if (dabs(irf1*hr-xbeg).lt.tin)  then
            xsav=hr*irf1
        else
            irf=int(xbeg/hr)
            if (irs.eq.1)  xsav=hr*irf
            if (irs.eq.-1) xsav=hr*(irf+1)
        end if
        xend= 0.5d0 - 0.5d0*irs + tin*irs
        if (ipri.gt.2) write (*,*) 'xbeg-xend',xbeg,xend
        hsav = -hr*irs
        h1 = hsav
        !---------------------------------------
        ! solve eqs. starting from xbeg
        !---------------------------------------
        ystart(1) = tet
        ystart(2) = xm
        call driver2(ystart,xbeg,xend,xsav,hmin,h1, pabs)
        tet = ystart(1)
        xm = ystart(2)
        ib2 = 0
            rnew2=ystart(3)
            cotet=dcos(tet)
            sitet=dsin(tet)
            xdl = fdf(rnew2,cdl,ncoef,xdlp)
            xly = fdf(rnew2,cly,ncoef,xlyp)
            xgm = fdf(rnew2,cgm,ncoef,xgmp)
            xx=-xdl+rnew2*cotet-xgm*sitet**2
            zz=rnew2*xly*sitet
            xxx=(r0+rm*xx)/1d2
            zzz=(z0+rm*zz)/1d2 
            ! open(33,file='lhcd/out/dots.dat',position="append")
            ! write(33,*)xxx, zzz, nomth, nomnz
            ! close(33)

        !---------------------------------------
        ! absorption
        !---------------------------------------
        if(iabsorp.ne.0) then
            if(ipri.gt.2) write (*,*)'in traj() iabsorp=',iabsorp
            nmax=nrefl
            return
        end if
        if (xend.eq.xbeg) nb1=nb1+1
        !sav2008 20    continue
  
        !--------------------------------------------------------
        !  pass turning point
        !--------------------------------------------------------
        irs0=irs
        ider=0
        call disp2(xend,xm,tet,xnr,prt,prm)
        ider=1
        ynz0 = ynz
  40    yy(1)=tet 
        yy(2)=xm
        yy(3)=xend
        yy(4)=xnr
        x1=0d0
        x2=1d+10
        rexi=xend
        call driver4(yy,x1,x2,rexi,hmin, extd4)
        if(iabsorp.eq.-1) return !failed to turn

        tetnew = yy(1)
        xmnew  = yy(2)
        rnew   = yy(3)
        xnrnew = yy(4)

            cotet = cos(tetnew)
            sitet = sin(tetnew)
            xdl = fdf(rnew,cdl,ncoef,xdlp)
            xly = fdf(rnew,cly,ncoef,xlyp)
            xgm = fdf(rnew,cgm,ncoef,xgmp)
            xx = -xdl+rnew*cotet-xgm*sitet**2
            zz = rnew*xly*sitet
            xxx = (r0+rm*xx)/1d2
            zzz = (z0+rm*zz)/1d2 
            ! open(33,file='lhcd/out/dots.dat',position="append")
            ! write(33,*)xxx, zzz, nomth, nomnz
            ! close(33)

        if(ipri.gt.2) write (*,*) 'from r=',rexi,'to r=',rnew
  
        !---------------------------------------
        ! find mode
        !---------------------------------------
        iroot=3
        ider=0
        xnrv=xnrnew
        call disp2(rnew,xmnew,tetnew,xnrv,prt,prm)
        ider=1
        iroot=1
        !ipric      if (ipri.gt.2) then
        !ipric       write (*,*)'nr check, r=',rnew,' tet=',tetnew
        !ipric       write (*,*)'iw=',iw,' izn=',izn
        !ipric       write (*,*) xnrnew,xnr1
        !ipric       write (*,*) xnr2,xnr3,xnr4
        !ipric       pause
        !ipric      end if
        pg1 = abs(xnrnew-xnr1)
        pg2 = abs(xnrnew-xnr2)
        pg3 = abs(xnrnew-xnr3)
        pg4 = abs(xnrnew-xnr4)

        pg = dmin1(pg1,pg2,pg3,pg4)
        if (dabs(pg/xnrnew).gt.pgdop) then
            !---------------------------------------------
            ! bad accuracy, continue with 4 equations
            !--------------------------------------------
            ib2=ib2+1
            nb2=nb2+1
            if (ib2.gt.4) then
                if (ipri.gt.1) write (*,*) 'error: cant leave 4 eqs'
                iabsorp=-1
                return
            end if
            eps=eps/5d0
            rrange=rrange*2d0
            hdrob=hdrob*2d0
            goto 40
        end if
        !-------------------------------------
        !          change wave type
        !-------------------------------------
        if (pg.ne.pg1) then
            if (pg.eq.pg2) izn=-izn
            if (pg.eq.pg3) iw=-iw
            if (pg.eq.pg4) iw=-iw
            if (pg.eq.pg4) izn=-izn
        end if
        if (irs0.ne.irs) nrefl=nrefl+1
        xbeg=rnew
        tet=tetnew
        xm=xmnew
        im4=1
        eps=eps0
        rrange=rrange0
        hdrob=hdrob0
        if(nrefl.lt.nmax) goto 10
        rzz=xbeg
        tetzz=tet
        xmzz=xm
        iznzz=izn
        iwzz=iw
        irszz=irs
    end  


    

end module trajectory