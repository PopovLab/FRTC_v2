module spectrum_mod
    use kind_module
    implicit none    

    type SpectrumPoint
        !real(wp) nz и ny не нужны
        ! 
        !real(wp) ny
        !
        real(wp) Ntor   
        !! Ntau=-Ntor   
        real(wp) Npol
        !! Ntet=Npol
        real(wp) power
        !! power

    contains
    end type SpectrumPoint

    type Spectrum
        integer size
        !! size of spectrum
        real(wp) input_power
        !! power of spectrum
        real(wp) power_ratio
        !! доля входной мощности
        real(wp) max_power
        !!
        real(wp) sum_power
        !! суммарная power        
        integer direction
        !! направление спектра   +1 или -1 или 0 - полный
        type(SpectrumPoint), allocatable ::  data(:)
        !! 
    contains
        procedure :: get_positive_part => get_positive_part_method
        procedure :: get_negative_part => get_negative_part_method
        procedure :: calc_max_power => calc_max_power_method
        procedure :: normalization => normalization_method
        procedure :: integral_trapez => integral_trapez_method
        procedure :: write => write_spectrum
    end type Spectrum

    interface Spectrum
        module procedure :: spectrum_constructor
        !module procedure :: read_spectrum
    end interface Spectrum    
contains
    function spectrum_constructor(size) result(this)
        !- конструктор для spectrum
        implicit none
        type(Spectrum) :: this
        integer, value :: size
        this%size = size
        this%input_power = 0
        this%sum_power = 0
        allocate(this%data(size))
    end function spectrum_constructor 

    subroutine calc_max_power_method(this)
        use constants, only: xsgs
        use rt_parameters, only : ntet
        implicit none
        class(Spectrum),  intent(inout) :: this
        type(SpectrumPoint) :: p           
        real(wp) max_power, pnorm
        integer i
        max_power = 0
        pnorm = this%input_power*xsgs/ntet
        print *, 'pnorm =', pnorm        
        do i = 1, this%size
            p = this%data(i)
            p%power = p%power*pnorm
            p%Ntor = this%direction * p%Ntor
            this%data(i) = p
            if (p%power>max_power)  max_power = p%power
        end do        
        
        this%max_power = max_power
        print *, 'this%max_power = ', this%max_power
    end subroutine

    function integral_trapez_method(this) result(sum)
        !! вычисление полной мощности спектра
        !! интегрирование методом трапеций
        implicit none
        class(Spectrum),  intent(in) :: this
        real(wp) :: sum
        integer :: i
        type(SpectrumPoint) :: p1, p2
        p1 = this%data(1)
        sum = 0
        do i = 2, this%size
            p2 = this%data(i)
            sum = sum + 0.5d0*(p2%power + p1%power)*(p2%Ntor-p1%Ntor)
            p1 = p2
        end do 
    end function
    
    subroutine normalization_method(this)
        use constants, only: xsgs
        use rt_parameters, only : ntet
        implicit none
        class(Spectrum),  intent(inout) :: this
        type(SpectrumPoint) :: p           
        real(wp) max_power, pnorm, p_sum
        integer i
        !p_sum = this%integral_trapez()
        p_sum = 0
        do i = 1, this%size
            p_sum = p_sum + this%data(i)%power
        end do
        max_power = 0
        pnorm = this%input_power*xsgs/ntet/p_sum
        print *, 'pnorm =', pnorm        
        do i = 1, this%size
            p = this%data(i)
            p%power = p%power*pnorm
            p%Ntor = this%direction * p%Ntor
            !p%Ntor = p%Ntor
            this%data(i) = p
            if (p%power>max_power)  max_power = p%power
        end do        
        
        this%max_power = max_power
        print *, 'this%max_power = ', this%max_power
    end subroutine

    function get_positive_part_method(this) result(spectr)
        !! 
        implicit none
        class(Spectrum), intent(in) :: this
        type(Spectrum) :: spectr, tmp_spectr
        type(SpectrumPoint) :: p        
        integer i, n
        print *, 'read positive'
        tmp_spectr = Spectrum(this%size)
        n = 0
        do i = 1, this%size
            p = this%data(i)
            if (p%Ntor>0) then
                
                n = n + 1                
                tmp_spectr%data(n) = p
                tmp_spectr%sum_power = tmp_spectr%sum_power + p%power
            end if
        end do
        tmp_spectr%size = n

        spectr = Spectrum(n)
        spectr%sum_power = tmp_spectr%sum_power
        do i = 1, n
            spectr%data(i) = tmp_spectr%data(i)
        end do
        spectr%size = n
        spectr%direction = +1
        spectr%power_ratio = spectr%sum_power/this%sum_power
        spectr%input_power = spectr%power_ratio * this%input_power
        print *, this%size, n        
        print *, 'sum_power ', this%sum_power, spectr%sum_power
        print *, 'power_ratio ', this%power_ratio, spectr%power_ratio
        print *, 'input_power ', this%input_power, spectr%input_power

    end function get_positive_part_method

    function get_negative_part_method(this) result(spectr)
        !! 
        implicit none
        class(Spectrum), intent(in) :: this
        type(Spectrum) :: spectr, tmp_spectr
        type(SpectrumPoint) :: p
        integer i, n
        print *, 'negative positive'
        tmp_spectr = Spectrum(this%size)
        n = 0
        do i = 1, this%size
            p = this%data(i)
            if (p%Ntor<0) then
                n = n + 1                
                p%Ntor = -p%Ntor
                tmp_spectr%data(n) = p
                tmp_spectr%sum_power = tmp_spectr%sum_power + p%power
            end if
        end do
        tmp_spectr%size = n

        spectr = Spectrum(n)
        spectr%sum_power = tmp_spectr%sum_power
        do i = 1, n
            spectr%data(i) = tmp_spectr%data(n + 1 - i)
        end do
        spectr%size = n
        spectr%direction = -1
        spectr%power_ratio = spectr%sum_power/this%sum_power
        spectr%input_power = spectr%power_ratio * this%input_power
        print *, this%size, n        
        print *, 'sum_power ', this%sum_power, spectr%sum_power
        print *, 'power_ratio ', this%power_ratio, spectr%power_ratio
        print *, 'input_power ', this%input_power, spectr%input_power

    end function get_negative_part_method

    function read_spectrum(file_name) result(spectr)
        !- чтение spectrum из файла
        implicit none
        type(Spectrum) :: spectr
        character (len = *), value :: file_name 
        logical                     :: res
        integer i,n,stat
        real(wp) sum_power
        !integer, value :: size
        print *, file_name      
        ! Check if the file exists
        inquire( file=trim(file_name), exist=res )        
        if (.not.res) then
            print *, 'spectrum file not exists'
            stop
        end if

        open(20,file=file_name)
        n=-1
        stat=0
        do while(stat == 0)
            n=n+1
            read (20,*,iostat=stat)
        enddo

        spectr%size = n
        spectr%input_power = 0
        spectr%sum_power = 0
        spectr%direction = 1
        spectr%power_ratio = 1
        sum_power = 0
        allocate(spectr%data(n))        
        print *,'Spectrum size = ',  n
        rewind(20)
        do i=1,n
            read (20,*) spectr%data(i)%Ntor, spectr%data(i)%Npol, spectr%data(i)%power
            sum_power = sum_power + spectr%data(i)%power
        enddo
        !sum_power
        !do i=1,n
        !    spectr%data(i)%power = spectr%data(i)%power/sum_power
        !enddo
        spectr%sum_power = sum_power
        close(20)

    end function read_spectrum         

    subroutine write_spectrum(this, spectrum_name)
        !! write spectrum to file
        implicit none
        class(Spectrum), intent(inout) :: this
        character(len=*), intent(in) :: spectrum_name
        character(len=256)  :: fname
        integer i, n
        integer, parameter :: iu = 21
        write(fname,'("lhcd/", A, ".dat")') spectrum_name
        print *,'write to:',  fname

        open(iu, file=fname)
        do i=1, this%size 
            write (iu,'(3(ES22.14))') this%data(i)%Ntor, this%data(i)%Npol, this%data(i)%power
        enddo
    end subroutine write_spectrum

    subroutine divide_spectrum(spectr, pos_spectr, neg_spectr)
        !! деление спектра на две части
        implicit none
        type(Spectrum), intent(in)  :: spectr
        type(Spectrum), intent(out) :: pos_spectr, neg_spectr
        type(Spectrum):: tmp_spectr
        type(SpectrumPoint) :: p
        integer i, pos_n, neg_n

        pos_spectr = Spectrum(spectr%size)
        tmp_spectr = Spectrum(spectr%size)
        pos_n = 0
        neg_n = 0
        do i = 1, spectr%size
            p = spectr%data(i)

            if (p%Ntor>0) then
                pos_n = pos_n + 1                
                pos_spectr%data(pos_n) = p
                pos_spectr%sum_power = pos_spectr%sum_power + p%power
            end if
            if (p%Ntor<0) then
                neg_n = neg_n + 1                
                p%Ntor = -p%Ntor
                tmp_spectr%data(neg_n) = p
                tmp_spectr%sum_power = tmp_spectr%sum_power + p%power
            endif
        end do
        pos_spectr%size = pos_n

        neg_spectr = Spectrum(neg_n)
        neg_spectr%sum_power = tmp_spectr%sum_power
        do i = 1, neg_n
            neg_spectr%data(i) = tmp_spectr%data(neg_n + 1 - i)
        end do
        neg_spectr%size = neg_n
        pos_spectr%direction = +1
        neg_spectr%direction = -1
        pos_spectr%power_ratio = pos_spectr%sum_power/spectr%sum_power
        neg_spectr%power_ratio = neg_spectr%sum_power/spectr%sum_power
        pos_spectr%input_power = pos_spectr%power_ratio * spectr%input_power
        neg_spectr%input_power = neg_spectr%power_ratio * spectr%input_power        
        print *, pos_n, neg_n        
        print *, 'sum_power ', spectr%sum_power, pos_spectr%sum_power, neg_spectr%sum_power
        print *, 'power_ratio ', pos_spectr%power_ratio, neg_spectr%power_ratio
        print *, 'input_power ', spectr%input_power, pos_spectr%input_power, neg_spectr%input_power
    end subroutine



    function make_spline_approximation(spectr) result(appx_spectr)
        !! approximation of input LH spectrum
            use constants, only: zero, xsgs
            use spline_module
            use rt_parameters, only: nnz, ntet, pabs0
            implicit none
            type(Spectrum), intent(in) :: spectr
            type(Spectrum) :: appx_spectr
            integer :: ispectr, ispl
            real(wp), allocatable :: ynzm0(:),pm0(:)
            real(wp), allocatable :: ynzm(:),pm(:)
            real(wp), allocatable :: yn2z(:),powinp(:)
            integer innz, i
            real(wp) dxx, xx0, xx1, xx2, yy1, yy2, pinp
            real(wp) dpw, dpower, pwcurr, ptot, dynn
            real(wp) pmax, pnorm, plaun
            ispectr = spectr%direction
            plaun = spectr%input_power
            ispl = spectr%size
            allocate(ynzm(nnz),pm(nnz))
            allocate(ynzm0(ispl),pm0(ispl))
            allocate(yn2z(ispl),powinp(ispl))
            do i = 1, spectr%size
                ynzm0(i) = spectr%data(i)%Ntor
                pm0(i) = spectr%data(i)%power
            end do
            call splne(ynzm0,pm0,ispl,yn2z)
            innz=100*ispl
            dxx=(ynzm0(ispl)-ynzm0(1))/innz
            xx2=ynzm0(1)
            yy2=pm0(1)
            pinp=0d0
            do i=1,innz
                  xx1=xx2
                  yy1=yy2
                  xx2=xx1+dxx
                  call splnt(ynzm0,pm0,yn2z,ispl,xx2,yy2,dynn)
                  dpw=.5d0*(yy2+yy1)*(xx2-xx1)
                  pinp=pinp+dpw
            end do
            dpower=pinp/dble(nnz)
            xx2=ynzm0(1)
            yy2=pm0(1)
            pwcurr= 0
            ptot=zero
            do i=1,nnz-1
                xx0=xx2
      11        continue
                xx1=xx2
                yy1=yy2
                xx2=xx1+dxx
                call splnt(ynzm0,pm0,yn2z,ispl,xx2,yy2,dynn)
                dpw=.5d0*(yy2+yy1)*(xx2-xx1)
                if(pwcurr+dpw.gt.dpower) then
                    xx2=xx1+dxx*(dpower-pwcurr)/dpw
                    call splnt(ynzm0,pm0,yn2z,ispl,xx2,yy2,dynn)
                    dpw=.5d0*(yy2+yy1)*(xx2-xx1)
                    pwcurr=pwcurr+dpw
                else
                    pwcurr=pwcurr+dpw
                    go to 11
                end if
                ynzm(i)=.5d0*(xx2+xx0)
                pm(i)=pwcurr
                ptot=ptot+pwcurr
                pwcurr=zero
            end do
            ynzm(nnz)=.5d0*(ynzm0(ispl)+xx2)
            pm(nnz)=pinp-ptot
            pnorm=plaun*xsgs/(pinp*ntet)
            pmax=-1d+10
            do i=1,nnz
                call splnt(ynzm0,pm0,yn2z,ispl,ynzm(i),powinp(i),dynn)
                pm(i)=pm(i)*pnorm
                if (pm(i).gt.pmax) pmax=pm(i)
                ynzm(i)=dble(ispectr)*ynzm(i) !sav2009
            end do
            !pabs=pabs0*pmax/1.d2
            appx_spectr = Spectrum(nnz)
            do i= 1, nnz
                appx_spectr%data(i) = SpectrumPoint(power = pm(i), Ntor = ynzm(i), Npol = 0)
            end do
            appx_spectr%input_power = plaun
            appx_spectr%max_power = pmax
            appx_spectr%direction = ispectr
            appx_spectr%power_ratio = spectr%power_ratio
        end function    
end module spectrum_mod

module spectrum1D
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    use spectrum_mod
    implicit none    

    type(Spectrum) :: full_spectrum
    type(Spectrum) :: pos_spectr, neg_spectr
    integer :: ispl
    !! size of spectrum
    real(wp) :: plaun
    !! power of spectrum
    real(wp) :: ynzm0(1001)
    !+
    real(wp) :: pm0(1001)
    !+ 
    real(wp) :: ynzm(1001), pm(1001)
    !!common /a0a1/ ynzm(1001),pm(1001)     
    real(wp) :: pabs
    !!common /a0gh/ pabs

    integer, parameter, private :: HEADER_LENGTH = 53

    contains
    subroutine read_positive_spectrum(file_name, p_in)
        implicit none
        character(*) file_name
        real(wp) :: p_in        
        integer, parameter :: iunit = 20        
        integer :: i, i1
        real(wp) :: anz,apz
        
        open(iunit, file= file_name)
        do i = 1, HEADER_LENGTH
            read(iunit,*)
        end do
        do i=1,10000
            read (iunit,*) anz,apz
            if(apz.eq.-88888.d0) then
                plaun=p_in*anz !input power in positive spectrum
                exit
            end if
            ynzm0(i)=anz
            pm0(i)=apz
            i1=i
        end do
        close(iunit)
        ispl=i1

        if(ispl.gt.4001) stop 'too many points in spectrum'

    end subroutine read_positive_spectrum        

    subroutine read_negative_spectrum(file_name, p_in)
        implicit none
        character(*) file_name
        real(wp) :: p_in
        integer, parameter :: iunit = 20        
        integer :: i, i1
        real(wp) :: anz, apz

        open(iunit, file= file_name)        
        do i = 1, HEADER_LENGTH
            read(iunit,*)
        end do
        apz=0.d0
        do while(apz.ne.-88888.d0)
            read (iunit,*) anz,apz
        end do
        read(iunit,*)
        plaun=p_in*(1.d0-anz) !input power in negative spectrum
        if (plaun > 0.d0) then
            do i=1,10000
                read (iunit,*,end=10) ynzm0(i),pm0(i)
            i1=i
            end do        
        end if
10      close(iunit)
        ispl=i1     

        if(ispl.gt.4001) stop 'too many points in spectrum'

    end subroutine read_negative_spectrum    

    subroutine spectrum_approximation(ispectr)
    !! approximation of input LH spectrum
        use constants, only: zero, xsgs
        use spline_module
        use rt_parameters, only: nnz, ntet, pabs0
        implicit none
        integer, intent(in) :: ispectr
        real(wp) yn2z(1001),powinp(1001)
        integer innz, i
        real(wp) dxx, xx0, xx1, xx2, yy1, yy2, pinp
        real(wp) dpw, dpower, pwcurr, ptot, dynn
        real(wp) pmax, pnorm
        call splne(ynzm0,pm0,ispl,yn2z)
        innz=100*ispl
        dxx=(ynzm0(ispl)-ynzm0(1))/innz
        xx2=ynzm0(1)
        yy2=pm0(1)
        pinp=0d0
        do i=1,innz
              xx1=xx2
              yy1=yy2
              xx2=xx1+dxx
              call splnt(ynzm0,pm0,yn2z,ispl,xx2,yy2,dynn)
              dpw=.5d0*(yy2+yy1)*(xx2-xx1)
              pinp=pinp+dpw
        end do
  
        dpower=pinp/dble(nnz)
        xx2=ynzm0(1)
        yy2=pm0(1)
        pwcurr=zero
        ptot=zero
        do i=1,nnz-1
            xx0=xx2
  11        continue
            xx1=xx2
            yy1=yy2
            xx2=xx1+dxx
            call splnt(ynzm0,pm0,yn2z,ispl,xx2,yy2,dynn)
            dpw=.5d0*(yy2+yy1)*(xx2-xx1)
            if(pwcurr+dpw.gt.dpower) then
                xx2=xx1+dxx*(dpower-pwcurr)/dpw
                call splnt(ynzm0,pm0,yn2z,ispl,xx2,yy2,dynn)
                dpw=.5d0*(yy2+yy1)*(xx2-xx1)
                pwcurr=pwcurr+dpw
            else
                pwcurr=pwcurr+dpw
                go to 11
            end if
            ynzm(i)=.5d0*(xx2+xx0)
            pm(i)=pwcurr
            ptot=ptot+pwcurr
            pwcurr=zero
        end do
        ynzm(nnz)=.5d0*(ynzm0(ispl)+xx2)
        pm(nnz)=pinp-ptot
        pnorm=plaun*xsgs/(pinp*ntet)
        pmax=-1d+10
        do i=1,nnz
            call splnt(ynzm0,pm0,yn2z,ispl,ynzm(i),powinp(i),dynn)
            pm(i)=pm(i)*pnorm
            if (pm(i).gt.pmax) pmax=pm(i)
            ynzm(i)=dble(ispectr)*ynzm(i) !sav2009
        end do
        pabs=pabs0*pmax/1.d2
    end subroutine

    subroutine copy_to_spectrum_1D(spectr)
        use spectrum_mod
        implicit none
        type(Spectrum) :: spectr
        type(SpectrumPoint) ::p
        integer i
        do i= 1, spectr%size
            p = spectr%data(i)
            ynzm0(i) = p%Ntor
            pm0(i) = p%power
        end do        
        plaun = spectr%input_power
        ispl = spectr%size
    end subroutine

    function create_spectrum() result(spectr)
        use spectrum_mod
        use rt_parameters, only: nnz
        implicit none
        type(Spectrum) :: spectr
        type(SpectrumPoint) :: p
        integer i 
        real(wp) :: pmax
        pmax = 0
        spectr = Spectrum(nnz)
        do i= 1, nnz
            p = SpectrumPoint(power = pm(i), Ntor = ynzm(i), Npol = 0)
            if (pm(i) > pmax) pmax=pm(i)
            spectr%data(i) = p
        end do
        spectr%max_power = pmax
    end function


end module spectrum1D