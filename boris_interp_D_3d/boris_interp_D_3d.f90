!  boris_interp_D_3d.f90 
!
!  FUNCTIONS:
!  boris_interp_D_3d - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: boris_interp_D_3d
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************
 SUBROUTINE ELAPSE_SECOND (ELAPSE) 
!
!  based on QueryPerformanceCounter 
!  
      use ifwin, only: T_LARGE_INTEGER,QueryPerformanceCounter, QueryPerformanceFrequency
! 
!     Returns the total elapsed time in seconds 
!     This is the fastest and most accurate timing routine 
! 
      real*8,   intent (out) :: elapse 
! 
      real*8    :: freq  = 1 
      logical*4 :: first = .true. 
      integer*8 :: start = 0 
      integer*8 :: num 
      logical*4 :: ll 
      type(T_LARGE_INTEGER) :: arg
! 
!   Calibrate this time using QueryPerformanceFrequency 
      if (first) then 
!         num   = 0 
         ll    = QueryPerformanceFrequency (arg) 
         num = transfer(arg,num)
         freq  = 1.0d0 / dble (num) 
!         start = 0 
         ll    = QueryPerformanceCounter (arg) 
         start = transfer(arg,start)
         first = .false. 
         WRITE (*,*) 'QueryPerformanceFrequency   :',num,' ticks per second' 
      end if 
! 
!      num    = 0 
      ll     = QueryPerformanceCounter (arg) 
      num = transfer(arg,num)
      elapse = dble (num-start) * freq 
      return 
    end 
    
    
    program boris_interp_D_3d

    implicit double precision(a-h,o-z)
    integer, parameter:: nx=701, ny=701, nz=701
    double precision B(3, 1), E(3, 1), speed(3, 1), matrix(3, 3), matrix_trans(3,3), x1(3, 1), x_0(3,1), x_d(nx), x_d1(nx+1), y_d(ny), y_d1(ny+1), z_d(nz), z_d1(nz+1), v_par(3, 1), v_ort(3, 1), E_par(3, 1), E_ort(3, 1)
    real, allocatable:: Bx(:,:,:), By(:,:,:), Bz(:,:,:), Ex(:,:,:), Ey(:,:,:), Ez(:,:,:)
    

100 format(f22.15, 1x, f22.15, 1x, f22.15, 1x, f22.15, 1x, f22.15, 1x, f22.15, 1x, f22.15, 1x, f22.15)
    
101 format(f14.11, 1x, f14.11, 1x, f14.11)
    q=1. 
    am=1. 
    c=1. 
    pi = 3.1415926535
    phi=pi/6.

    matrix(1, 1)=1.
    matrix(1, 2)=0.
    matrix(1, 3)=0.
    matrix(2, 1)=0.
    matrix(3, 1)=0.
    matrix(2, 2)=cos(phi)
    matrix(2, 3)=-sin(phi)
    matrix(3, 2)=sin(phi)
    matrix(3, 3)=cos(phi)
    matrix_trans=transpose(matrix)
    !x0=0.
    !y0=-1.
    !z0=-1.
    !xn=1.1
    !yn=0.1
    !zn=0.1
    
    x0=-1.1
    y0=-1.1
    z0=-1.1
    xn=1.1
    yn=1.1
    zn=1.1
    hx=(xn-x0)/(nx-1)
    hy=(yn-y0)/(ny-1)
    hz=(zn-z0)/(nz-1)
    !open(7, FILE="interp_d_3d_int_ht=0025_nt=40000.txt", access="append")

    t0=0.
    tn=1000.
    open(2, FILE="boris_interp_D_3d_phi=30_ht=0003125_nt=32000_var1_2.txt")! interp_3d_ht=pi160_ht=pi160_2 - плоский вариант в 3д
    ht=0.003125
    nt=32000
    k=64
    sum=0.
    sum_int=0.
    sum_fil=0.
    do  i=1, nx
        x_d(i)=x0+i*hx
        x_d1(i)=x0-hx/2+i*hx
    end do
    x_d1(nx+1)=xn+hx/2.
    do  i=1, ny
        y_d(i)=y0+i*hy
        y_d1(i)=y0-hy/2+i*hy
    end do
    y_d1(ny+1)=yn+hy/2.
    do  i=1, nz
        z_d(i)=x0+i*hz
        z_d1(i)=x0-hz/2+i*hz
    end do
    y_d1(nz+1)=zn+hz/2.

    !do jn=1, 100
        allocate(Bx(nx+1, ny, nz))
        allocate(By(nx, ny+1, nz))
        allocate(Bz(nx, ny, nz+1))
        allocate(Ex(nx, ny+1, nz+1))
        allocate(Ey(nx+1, ny, nz+1))
        allocate(Ez(nx+1, ny+1, nz))
        call ELAPSE_SECOND(start) 
        Bx0=0.
        By0=0.
        Ez0=0.
        do i=1, nx+1
            do j=1,ny
                do k=1, nz
                    x1=reshape((/x_d1(i), y_d(j), z_d(k)/), (/3, 1/))
                    x_0=matmul(matrix_trans, x1)
                    !Bz0=sqrt(x_d1(i)**2+y_d(j)**2)
                    Bz0=sqrt(x_0(1,1)**2+x_0(2,1)**2)
                    !Bz0=1.
                    Bx(i, j, k)=Bx0*matrix(1, 1)+By0*matrix(1, 2)+Bz0*matrix(1, 3)
                end do
            end do
        end do
    
        do  i=1, nx
            do  j=1, ny+1
                do k=1, nz
                    x1=reshape((/x_d(i), y_d1(j), z_d(k)/), (/3, 1/))
                    x_0=matmul(matrix_trans, x1)
                    Bz0=sqrt(x_0(1,1)**2+x_0(2,1)**2)
                    !Bz0=sqrt(x_d(i)**2+y_d1(j)**2)
                    !Bz0=1.
                    By(i, j, k)=Bx0*matrix(2, 1)+By0*matrix(2, 2)+Bz0*matrix(2, 3)
                end do
            end do
        end do
    
    
        do i=1, nx
            do j=1, ny
                do  k=1, nz+1
                    x1=reshape((/x_d(i), y_d(j), z_d1(k)/), (/3, 1/))
                    x_0=matmul(matrix_trans, x1)
                    Bz0=sqrt(x_0(1,1)**2+x_0(2,1)**2)
                    !Bz0=sqrt(x_d(i)**2+y_d(j)**2)
                    !Bz0=1.
                    Bz(i, j, k)=Bx0*matrix(3, 1)+By0*matrix(3, 2)+Bz0*matrix(3, 3)
                end do
            end do
        end do
            
        do i=1, nx
            do j=1, ny+1
                do k=1, nz+1
                    x1=reshape((/x_d(i), y_d1(j), z_d1(k)/), (/3, 1/))
                    x_0=matmul(matrix_trans, x1)
                    !Ex0=0.01*x_d(i)/(sqrt(x_d(i)**2+y_d1(j)**2))**3
                    !Ey0=0.01*y_d1(j)/(sqrt(x_d(i)**2+y_d1(j)**2))**3
                    !Ex0=0.
                    !Ey0=0.
                    Ex0=0.01*x_0(1,1)/(sqrt(x_0(1,1)**2+x_0(2,1)**2))**3
                    Ey0=0.01*x_0(2,1)/(sqrt(x_0(1,1)**2+x_0(2,1)**2))**3
                    Ex(i, j, k)=Ex0*matrix(1, 1)+Ey0*matrix(1, 2)+Ez0*matrix(1, 3)
                end do
            end do
        end do
            
        do i=1, nx+1
            do j=1, ny
                do k=1, nz+1
                    x1=reshape((/x_d1(i), y_d(j), z_d1(k)/), (/3, 1/))
                    x_0=matmul(matrix_trans, x1)
                    !Ex0=0.01*x_d1(i)/(sqrt(x_d1(i)**2+y_d(j)**2))**3
                    !Ey0=0.01*y_d(j)/(sqrt(x_d1(i)**2+y_d(j)**2))**3
                    !Ex0=0.
                    !Ey0=0.
                    Ex0=0.01*x_0(1,1)/(sqrt(x_0(1,1)**2+x_0(2,1)**2))**3
                    Ey0=0.01*x_0(2,1)/(sqrt(x_0(1,1)**2+x_0(2,1)**2))**3
                    Ey(i, j, k)=Ex0*matrix(2, 1)+Ey0*matrix(2, 2)+Ez0*matrix(2, 3)
                end do
            end do
        end do
            
        do i=1, nx+1
            do j=1, ny+1
                do k=1, nz
                    x1=reshape((/x_d1(i), y_d1(j), z_d(k)/), (/3, 1/))
                    x_0=matmul(matrix_trans, x1)
                    !Ex0=0.01*x_d1(i)/(sqrt(x_d1(i)**2+y_d1(j)**2))**3
                    !Ey0=0.01*y_d1(j)/(sqrt(x_d1(i)**2+y_d1(j)**2))**3
                    !Ex0=0.
                    !Ey0=0.
                    Ex0=0.01*x_0(1,1)/(sqrt(x_0(1,1)**2+x_0(2,1)**2))**3
                    Ey0=0.01*x_0(2,1)/(sqrt(x_0(1,1)**2+x_0(2,1)**2))**3
                    Ez(i, j, k)=Ex0*matrix(3, 1)+Ey0*matrix(3, 2)+Ez0*matrix(3, 3)
                end do
            end do
        end do
        call ELAPSE_SECOND(fin)
        sum_fil=sum_fil+fin-start
        x=0.9
        y=0.
        z=0.
        u_01=0.1
        v_01=0.
        w_01=0.
    
        u_0=u_01*matrix(1, 1)+v_01*matrix(1, 2)+w_01*matrix(1, 3)
        v_0=u_01*matrix(2, 1)+v_01*matrix(2, 2)+w_01*matrix(2, 3)
        w_0=u_01*matrix(3, 1)+v_01*matrix(3, 2)+w_01*matrix(3, 3)
    
        s2 = (x-x0) / hx
        i = idint(s2 + 1.d0)
        i1 = idint(s2 + 1.5d0)
        s1 = i - s2
        s2 = i1 - 0.5d0 - s2
        s4 = (y-y0) / hy
        l = idint(s4 + 1.d0)
        l1 = idint(s4 + 1.5d0)
        s3 = l - s4
        s4 = l1 - 0.5d0 - s4
        s6 = (z-z0) / hz
        k = idint(s6 + 1.d0)
        k1 = idint(s6 + 1.5d0)
        s5 = k - s6
        s6 = k1 - 0.5d0 - s6
        s11 = 1.d0 - s1
        s21 = 1.d0 - s2
        s31 = 1.d0 - s3
        s41 = 1.d0 - s4
        s51 = 1.d0 - s5
        s61 = 1.d0 - s6
    
        Ex_total = s1 * (s4 * (s6 * Ex(i, l1, k1) + s61 * Ex(i, l1, k1 + 1)) +&
	    +s41 * (s6 * Ex(i, l1 + 1, k1) + s61 * Ex(i, l1 + 1, k1 + 1))) +&
	    +s11 * (s4 * (s6 * Ex(i + 1, l1, k1) + s61 * Ex(i + 1, l1, k1 + 1)) +&
	    +s41 * (s6 * Ex(i + 1, l1 + 1, k1) + s61 * Ex(i + 1, l1 + 1, k1 + 1)))
        
        Ey_total = s2 * (s3 * (s6 * Ey(i1, l, k1) + s61 * Ey(i1, l, k1 + 1)) +&
	    +s31 * (s6 * Ey(i1, l + 1, k1) + s61 * Ey(i1, l + 1, k1 + 1))) +&
	    +s21 * (s3 * (s6 * Ey(i1 + 1, l, k1) + s61 * Ey(i1 + 1, l, k1 + 1)) +&
		    +s31 * (s6 * Ey(i1 + 1, l + 1, k1) + s61 * Ey(i1 + 1, l + 1, k1 + 1)))
    
	    Ez_total = s2 * (s4 * (s5 * Ez(i1, l1, k) + s51 * Ez(i1, l1, k + 1)) +&
	    +s41 * (s5 * Ez(i1, l1 + 1, k) + s51 * Ez(i1, l1 + 1, k + 1))) +&
	    +s21 * (s4 * (s5 * Ez(i1 + 1, l1, k) + s51 * Ez(i1 + 1, l1, k + 1)) +&
		    +s41 * (s5 * Ez(i1 + 1, l1 + 1, k) + s51 * Ez(i1 + 1, l1 + 1, k + 1)))
    
	    Bx_total = s2 * (s3 * (s5 * Bx(i1, l, k) + s51 * Bx(i1, l, k + 1)) +&
	    +s31 * (s5 * Bx(i1, l + 1, k) + s51 * Bx(i1, l + 1, k + 1))) +&
	    +s21 * (s3 * (s5 * Bx(i1 + 1, l, k) + s51 * Bx(i1 + 1, l, k + 1)) +&
		    +s31 * (s5 * Bx(i1 + 1, l + 1, k) + s51 * Bx(i1 + 1, l + 1, k + 1)))
    
	    By_total = s1 * (s4 * (s5 * By(i, l1, k) + s51 * By(i, l1, k + 1)) +&
	    +s41 * (s5 * By(i, l1 + 1, k) + s51 * By(i, l1 + 1, k + 1))) +&
	    +s11 * (s4 * (s5 * By(i + 1, l1, k) + s51 * By(i + 1, l1, k + 1)) +&
		    +s41 * (s5 * By(i + 1, l1 + 1, k) + s51 * By(i + 1, l1 + 1, k + 1)))
    
	    Bz_total = s1 * (s3 * (s61 * Bz(i, l, k1) + s6 * Bz(i, l, k1 + 1)) +&
	    +s31 * (s61 * Bz(i, l + 1, k1) + s6 * Bz(i, l + 1, k1 + 1))) +&
	    +s11 * (s3 * (s61 * Bz(i + 1, l, k1) + s6 * Bz(i + 1, l, k1 + 1)) +&
		    +s31 * (s61 * Bz(i + 1, l + 1, k1) + s6 * Bz(i + 1, l + 1, k1 + 1)))
    
        ux=u_0-ht*q/2.*(Ex_total+(v_0*Bz_total-w_0*By_total)/c)
        uy=v_0-ht*q/2.*(Ey_total+(w_0*Bx_total-u_0*Bz_total)/c)
        uz=w_0-ht*q/2.*(Ez_total+(u_0*By_total-v_0*Bx_total)/c)
    
        i_t=1
        t=t0+ht
        energy=am*(ux*ux+uy*uy+uz*uz)/2
        write(2, 100) t0, x, y, z, ux, uy, uz, energy
        s=0.d0
        p=q*ht/am/c/2.
        time=0.
        time_int=0.
        do n=1,nt
            call ELAPSE_SECOND(start)
            s2 = (x-x0) / hx
            i = idint(s2 + 1.d0)
            i1 = idint(s2 + 1.5d0)
            s1 = i - s2
            s2 = i1 - 0.5d0 - s2
            s4 = (y-y0) / hy
            l = idint(s4 + 1.d0)
            l1 = idint(s4 + 1.5d0)
            s3 = l - s4
            s4 = l1 - 0.5d0 - s4
            s6 = (z-z0) / hz
            k = idint(s6 + 1.d0)
            k1 = idint(s6 + 1.5d0)
            s5 = k - s6
            s6 = k1 - 0.5d0 - s6
            s11 = 1.d0 - s1 
            s21 = 1.d0 - s2
            s31 = 1.d0 - s3
            s41 = 1.d0 - s4
            s51 = 1.d0 - s5
            s61 = 1.d0 - s6
    
    
            Ex_total = s1 * (s4 * (s6 * Ex(i, l1, k1) + s61 * Ex(i, l1, k1 + 1)) +&
	        +s41 * (s6 * Ex(i, l1 + 1, k1) + s61 * Ex(i, l1 + 1, k1 + 1))) +&
	        +s11 * (s4 * (s6 * Ex(i + 1, l1, k1) + s61 * Ex(i + 1, l1, k1 + 1)) +&
		    +s41 * (s6 * Ex(i + 1, l1 + 1, k1) + s61 * Ex(i + 1, l1 + 1, k1 + 1)))
        
            Ey_total = s2 * (s3 * (s6 * Ey(i1, l, k1) + s61 * Ey(i1, l, k1 + 1)) +&
		    +s31 * (s6 * Ey(i1, l + 1, k1) + s61 * Ey(i1, l + 1, k1 + 1))) +&
		    +s21 * (s3 * (s6 * Ey(i1 + 1, l, k1) + s61 * Ey(i1 + 1, l, k1 + 1)) +&
			    +s31 * (s6 * Ey(i1 + 1, l + 1, k1) + s61 * Ey(i1 + 1, l + 1, k1 + 1)))
    
	        Ez_total = s2 * (s4 * (s5 * Ez(i1, l1, k) + s51 * Ez(i1, l1, k + 1)) +&
		    +s41 * (s5 * Ez(i1, l1 + 1, k) + s51 * Ez(i1, l1 + 1, k + 1))) +&
		    +s21 * (s4 * (s5 * Ez(i1 + 1, l1, k) + s51 * Ez(i1 + 1, l1, k + 1)) +&
			    +s41 * (s5 * Ez(i1 + 1, l1 + 1, k) + s51 * Ez(i1 + 1, l1 + 1, k + 1)))
    
	        Bx_total = s2 * (s3 * (s5 * Bx(i1, l, k) + s51 * Bx(i1, l, k + 1)) +&
		    +s31 * (s5 * Bx(i1, l + 1, k) + s51 * Bx(i1, l + 1, k + 1))) +&
		    +s21 * (s3 * (s5 * Bx(i1 + 1, l, k) + s51 * Bx(i1 + 1, l, k + 1)) +&
			    +s31 * (s5 * Bx(i1 + 1, l + 1, k) + s51 * Bx(i1 + 1, l + 1, k + 1)))
    
	        By_total = s1 * (s4 * (s5 * By(i, l1, k) + s51 * By(i, l1, k + 1)) +&
		    +s41 * (s5 * By(i, l1 + 1, k) + s51 * By(i, l1 + 1, k + 1))) +&
		    +s11 * (s4 * (s5 * By(i + 1, l1, k) + s51 * By(i + 1, l1, k + 1)) +&
			    +s41 * (s5 * By(i + 1, l1 + 1, k) + s51 * By(i + 1, l1 + 1, k + 1)))
    
	        Bz_total = s1 * (s3 * (s6 * Bz(i, l, k1) + s61 * Bz(i, l, k1 + 1)) +&
		    +s31 * (s6 * Bz(i, l + 1, k1) + s61 * Bz(i, l + 1, k1 + 1))) +&
		    +s11 * (s3 * (s6 * Bz(i + 1, l, k1) + s61 * Bz(i + 1, l, k1 + 1)) +&
			    +s31 * (s6 * Bz(i + 1, l + 1, k1) + s61 * Bz(i + 1, l + 1, k1 + 1)))
            call ELAPSE_SECOND(finish_int)
            b_0=sqrt(Bx_total**2+By_total**2+Bz_total**2)
            B=reshape((/Bx_total, By_total, Bz_total/), (/3, 1/))
            E=reshape((/Ex_total, Ey_total, Ez_total/), (/3, 1/))
            speed=reshape((/ux, uy, uz/), (/3, 1/))
        
            do j=1, 3
                B(j, 1)=B(j, 1)/b_0
            end do
            !Bx_total=Bx_total/b_0
            !By_total=By_total/b_0
            !Bz_total=Bz_total/b_0
            vb=0.
            eb=0.
            do j=1, 3
                vb=vb+speed(j, 1)*B(j, 1)
                eb=eb+E(j, 1)*B(j, 1)
            end do
            !vb=ux*Bx_total+uy*By_total+uz*Bz_total
            !eb=Ex_total*Bx_total+Ey_total*By_total+Ez_total*Bz_total
            !v_parx=vb*Bx_total
            !v_pary=vb*By_total
            !v_parz=vb*Bz_total
            !v_ortx=ux-v_parx
            !v_orty=uy-v_pary
            !v_ortz=uz-v_parz
            !E_parx=eb*Bx_total
            !E_pary=eb*By_total
            !E_parz=eb*Bz_total
            !E_ortx=Ex_total-E_parx
            !E_orty=Ey_total-E_pary
            !E_ortz=Ez_total-E_parz
            do j=1, 3
                v_par(j, 1)=vb*B(j, 1)
                v_ort(j, 1)=speed(j, 1)-v_par(j, 1)
                E_par(j, 1)=eb*B(j, 1)
                E_ort(j, 1)=(E(j, 1)-E_par(j, 1))*c/b_0
            end do
            theta=q*b_0/am/c
            cs=cos(theta*ht)
            cs1=1-cs
            sn=sin(theta*ht)
!            write(*,*) cs, cs1, sn
            ux=E_ort(2, 1)*B(3, 1)-E_ort(3, 1)*B(2, 1)+(v_ort(1, 1)-E_ort(2, 1)*B(3, 1)+E_ort(3, 1)*B(2, 1))*cs+(v_ort(2, 1)*B(3, 1)-v_ort(3, 1)*B(2, 1)+E_ort(1, 1))*sn+v_par(1, 1)
            uy=E_ort(3, 1)*B(1, 1)-E_ort(1, 1)*B(3, 1)+(v_ort(2, 1)-E_ort(3, 1)*B(1, 1)+E_ort(1, 1)*B(3, 1))*cs+(v_ort(3, 1)*B(1, 1)-v_ort(1, 1)*B(3, 1)+E_ort(2, 1))*sn+v_par(2, 1)
            uz=E_ort(1, 1)*B(2, 1)-E_ort(2, 1)*B(1, 1)+(v_ort(3, 1)-E_ort(1, 1)*B(2, 1)+E_ort(2, 1)*B(1, 1))*cs+(v_ort(1, 1)*B(2, 1)-v_ort(2, 1)*B(1, 1)+E_ort(3, 1))*sn+v_par(3, 1)
        
            !ux=E_orty*Bz_total*cs1-E_ortz*By_total*cs1+E_ortx*sn+v_ortx*cs+(v_orty*Bz_total-v_ortz*By_total)*sn+v_parx
            !uy=E_ortz*Bx_total*cs1-E_ortx*Bz_total*cs1+E_orty*sn+v_orty*cs+(v_ortz*Bx_total-v_ortx*Bz_total)*sn+v_pary
            !uz=E_ortx*By_total*cs1-E_orty*Bx_total*cs1+E_ortz*sn+v_ortz*cs+(v_ortx*By_total-v_orty*Bx_total)*sn+v_parz
    
            x=x+ht*ux
            y=y+ht*uy
            z=z+ht*uz
    
            call ELAPSE_SECOND(finish)
            !write(*,*) n, start, finish_int, finish
            energy=am*(ux**2+uy**2+uz**2)/2.
            time=time+finish-start
            time_int=time_int+finish_int-start
            write(2, 100) t, x, y, z, ux, uy, uz, energy
            t=i_t*ht+t0
            i_t=i_t+1
        end do
        sum=sum+time
        sum_int=sum_int+time_int
        !write(7,*) jn, sum, sum_int, sum_fil
        !write(*,*) jn, sum, sum_int, sum_fil
        deallocate(Bx)
        deallocate(By)
        deallocate(Bz)
        deallocate(Ex)
        deallocate(Ey)
        deallocate(Ez)
   ! end do
   ! sum=sum/100
   ! sum_int=sum_int/100
   ! sum_fil=sum_fil/100
   ! open(4, file="time_d_3d_interp.txt", access='append')
   ! write(4,*) nt, ht, sum, sum_int, sum_fil
    pause

    end program boris_interp_D_3d
!
!!  boris_interp_D_3d.f90 
!!
!!  FUNCTIONS:
!!  boris_interp_D_3d - Entry point of console application.
!!
!
!!****************************************************************************
!!
!!  PROGRAM: boris_interp_D_3d
!!
!!  PURPOSE:  Entry point for the console application.
!!
!!****************************************************************************
! SUBROUTINE ELAPSE_SECOND (ELAPSE) 
!!
!!  based on QueryPerformanceCounter 
!!  
!      use ifwin, only: T_LARGE_INTEGER,QueryPerformanceCounter, QueryPerformanceFrequency
!! 
!!     Returns the total elapsed time in seconds 
!!     This is the fastest and most accurate timing routine 
!! 
!      real*8,   intent (out) :: elapse 
!! 
!      real*8    :: freq  = 1 
!      logical*4 :: first = .true. 
!      integer*8 :: start = 0 
!      integer*8 :: num 
!      logical*4 :: ll 
!      type(T_LARGE_INTEGER) :: arg
!! 
!!   Calibrate this time using QueryPerformanceFrequency 
!      if (first) then 
!!         num   = 0 
!         ll    = QueryPerformanceFrequency (arg) 
!         num = transfer(arg,num)
!         freq  = 1.0d0 / dble (num) 
!!         start = 0 
!         ll    = QueryPerformanceCounter (arg) 
!         start = transfer(arg,start)
!         first = .false. 
!         WRITE (*,*) 'QueryPerformanceFrequency   :',num,' ticks per second' 
!      end if 
!! 
!!      num    = 0 
!      ll     = QueryPerformanceCounter (arg) 
!      num = transfer(arg,num)
!      elapse = dble (num-start) * freq 
!      return 
!    end 
!    
!    
!    program boris_interp_D_3d
!
!    implicit double precision(a-h,o-z)
!    integer, parameter:: nx=201, ny=201, nz=201
!    double precision B(3, 1), E(3, 1), speed(3, 1), matrix(3, 3), matrix_trans(3,3), x1(3, 1), x_0(3,1), x_d(nx), x_d1(nx+1), y_d(ny), y_d1(ny+1), z_d(nz), z_d1(nz+1), v_par(3, 1), v_ort(3, 1), E_par(3, 1), E_ort(3, 1)
!    real, allocatable:: Bx(:,:,:), By(:,:,:), Bz(:,:,:), Ex(:,:,:), Ey(:,:,:), Ez(:,:,:)
!    
!    allocate(Bx(nx+1, ny, nz))
!    allocate(By(nx, ny+1, nz))
!    allocate(Bz(nx, ny, nz+1))
!    allocate(Ex(nx, ny+1, nz+1))
!    allocate(Ey(nx+1, ny, nz+1))
!    allocate(Ez(nx+1, ny+1, nz))
!
!100 format(f22.15, 1x, f22.15, 1x, f22.15, 1x, f22.15, 1x, f22.15, 1x, f22.15, 1x, f22.15, 1x, f22.15)
!    
!101 format(f14.11, 1x, f14.11, 1x, f14.11)
!    q=1. 
!    am=1. 
!    c=1. 
!    pi = 3.1415926535
!    phi=pi/6.
!
!    matrix(1, 1)=1.
!    matrix(1, 2)=0.
!    matrix(1, 3)=0.
!    matrix(2, 1)=0.
!    matrix(3, 1)=0.
!    matrix(2, 2)=cos(phi)
!    matrix(2, 3)=-sin(phi)
!    matrix(3, 2)=sin(phi)
!    matrix(3, 3)=cos(phi)
!    matrix_trans=transpose(matrix)
!    !x0=0.
!    !y0=-1.
!    !z0=-1.
!    !xn=1.1
!    !yn=0.1
!    !zn=0.1
!    
!    x0=-1.1
!    y0=-1.1
!    z0=-1.1
!    xn=1.1
!    yn=1.1
!    zn=1.1
!    hx=(xn-x0)/(nx-1)
!    hy=(yn-y0)/(ny-1)
!    hz=(zn-z0)/(nz-1)
!    
!    do  i=1, nx
!        x_d(i)=x0+i*hx
!        x_d1(i)=x0-hx/2+i*hx
!    end do
!    x_d1(nx+1)=xn+hx/2.
!    do  i=1, ny
!        y_d(i)=y0+i*hy
!        y_d1(i)=y0-hy/2+i*hy
!    end do
!    y_d1(ny+1)=yn+hy/2.
!    do  i=1, nz
!        z_d(i)=x0+i*hz
!        z_d1(i)=x0-hz/2+i*hz
!    end do
!    y_d1(nz+1)=zn+hz/2.
!    
!    Bx0=0.
!    By0=0.
!    Ez0=0.
!    do i=1, nx+1
!        do j=1,ny
!            do k=1, nz
!                x1=reshape((/x_d1(i), y_d(j), z_d(k)/), (/3, 1/))
!                x_0=matmul(matrix_trans, x1)
!                !Bz0=sqrt(x_d1(i)**2+y_d(j)**2)
!                Bz0=sqrt(x_0(1,1)**2+x_0(2,1)**2)
!                !Bz0=1.
!                Bx(i, j, k)=Bx0*matrix(1, 1)+By0*matrix(1, 2)+Bz0*matrix(1, 3)
!            end do
!        end do
!    end do
!    
!    do  i=1, nx
!        do  j=1, ny+1
!            do k=1, nz
!                x1=reshape((/x_d(i), y_d1(j), z_d(k)/), (/3, 1/))
!                x_0=matmul(matrix_trans, x1)
!                Bz0=sqrt(x_0(1,1)**2+x_0(2,1)**2)
!                !Bz0=sqrt(x_d(i)**2+y_d1(j)**2)
!                !Bz0=1.
!                By(i, j, k)=Bx0*matrix(2, 1)+By0*matrix(2, 2)+Bz0*matrix(2, 3)
!            end do
!        end do
!    end do
!    
!    
!    do i=1, nx
!        do j=1, ny
!            do  k=1, nz+1
!                x1=reshape((/x_d(i), y_d(j), z_d1(k)/), (/3, 1/))
!                x_0=matmul(matrix_trans, x1)
!                Bz0=sqrt(x_0(1,1)**2+x_0(2,1)**2)
!                !Bz0=sqrt(x_d(i)**2+y_d(j)**2)
!                !Bz0=1.
!                Bz(i, j, k)=Bx0*matrix(3, 1)+By0*matrix(3, 2)+Bz0*matrix(3, 3)
!            end do
!        end do
!    end do
!            
!    do i=1, nx
!        do j=1, ny+1
!            do k=1, nz+1
!                x1=reshape((/x_d(i), y_d1(j), z_d1(k)/), (/3, 1/))
!                x_0=matmul(matrix_trans, x1)
!                !Ex0=0.01*x_d(i)/(sqrt(x_d(i)**2+y_d1(j)**2))**3
!                !Ey0=0.01*y_d1(j)/(sqrt(x_d(i)**2+y_d1(j)**2))**3
!                !Ex0=0.
!                !Ey0=0.
!                Ex0=0.01*x_0(1,1)/(sqrt(x_0(1,1)**2+x_0(2,1)**2))**3
!                Ey0=0.01*x_0(2,1)/(sqrt(x_0(1,1)**2+x_0(2,1)**2))**3
!                Ex(i, j, k)=Ex0*matrix(1, 1)+Ey0*matrix(1, 2)+Ez0*matrix(1, 3)
!            end do
!        end do
!    end do
!            
!    do i=1, nx+1
!        do j=1, ny
!            do k=1, nz+1
!                x1=reshape((/x_d1(i), y_d(j), z_d1(k)/), (/3, 1/))
!                x_0=matmul(matrix_trans, x1)
!                !Ex0=0.01*x_d1(i)/(sqrt(x_d1(i)**2+y_d(j)**2))**3
!                !Ey0=0.01*y_d(j)/(sqrt(x_d1(i)**2+y_d(j)**2))**3
!                !Ex0=0.
!                !Ey0=0.
!                Ex0=0.01*x_0(1,1)/(sqrt(x_0(1,1)**2+x_0(2,1)**2))**3
!                Ey0=0.01*x_0(2,1)/(sqrt(x_0(1,1)**2+x_0(2,1)**2))**3
!                Ey(i, j, k)=Ex0*matrix(2, 1)+Ey0*matrix(2, 2)+Ez0*matrix(2, 3)
!            end do
!        end do
!    end do
!            
!    do i=1, nx+1
!        do j=1, ny+1
!            do k=1, nz
!                x1=reshape((/x_d1(i), y_d1(j), z_d(k)/), (/3, 1/))
!                x_0=matmul(matrix_trans, x1)
!                !Ex0=0.01*x_d1(i)/(sqrt(x_d1(i)**2+y_d1(j)**2))**3
!                !Ey0=0.01*y_d1(j)/(sqrt(x_d1(i)**2+y_d1(j)**2))**3
!                !Ex0=0.
!                !Ey0=0.
!                Ex0=0.01*x_0(1,1)/(sqrt(x_0(1,1)**2+x_0(2,1)**2))**3
!                Ey0=0.01*x_0(2,1)/(sqrt(x_0(1,1)**2+x_0(2,1)**2))**3
!                Ez(i, j, k)=Ex0*matrix(3, 1)+Ey0*matrix(3, 2)+Ez0*matrix(3, 3)
!            end do
!        end do
!    end do
!    time=0.
!    do i =1, 100
!        start=0.
!        finish=0.
!        time=time+bor(Bx, By, Bz, Ex, Ey, Ez, hx, hy, hz, nx, ny, nz, start, finish)
!    end do
!    open(4, file="time_d_3d_interp.txt", access='append')
!    write(4,*) nt, ht, time
!    pause
!
!    end program boris_interp_D_3d
!
!    function bor(Bx, By, Bz, Ex, Ey, Ez, hx, hy, hz, nx, ny, nz, start, finish)
!    double precision bor
!    implicit double precision(a-h,o-z)
!   ! integer, parameter:: nx=201, ny=201, nz=201
!    double precision B(3, 1), E(3, 1), speed(3, 1), matrix(3, 3), matrix_trans(3,3), x1(3, 1), x_0(3,1), x_d(nx), x_d1(nx+1), y_d(ny), y_d1(ny+1), z_d(nz), z_d1(nz+1), v_par(3, 1), v_ort(3, 1), E_par(3, 1), E_ort(3, 1)
!    !real, allocatable:: Bx(:,:,:), By(:,:,:), Bz(:,:,:), Ex(:,:,:), Ey(:,:,:), Ez(:,:,:)
!    !
!    !allocate(Bx(nx+1, ny, nz))
!    !allocate(By(nx, ny+1, nz))
!    !allocate(Bz(nx, ny, nz+1))
!    !allocate(Ex(nx, ny+1, nz+1))
!    !allocate(Ey(nx+1, ny, nz+1))
!    !allocate(Ez(nx+1, ny+1, nz))
!
!    x0=-1.1
!    y0=-1.1
!    z0=-1.1
!    xn=1.1
!    yn=1.1
!    zn=1.1
!    t0=0.
!    tn=1000.
!    k=64
!    !open(2, FILE="boris_interp_D_3d_phi=30_ht=00125_nt=100_var1.txt")! interp_3d_ht=pi160_ht=pi160_2 - плоский вариант в 3д
!    ht=0.0125
!    nt=100
!
!    x=0.9
!    y=0.
!    z=0.
!    u_01=0.1
!    v_01=0.
!    w_01=0.
!
!    u_0=u_01*matrix(1, 1)+v_01*matrix(1, 2)+w_01*matrix(1, 3)
!    v_0=u_01*matrix(2, 1)+v_01*matrix(2, 2)+w_01*matrix(2, 3)
!    w_0=u_01*matrix(3, 1)+v_01*matrix(3, 2)+w_01*matrix(3, 3)
!    
!    s2 = (x-x0) / hx
!    i = idint(s2 + 1.d0)
!    i1 = idint(s2 + 1.5d0)
!    s1 = i - s2
!    s2 = i1 - 0.5d0 - s2
!    s4 = (y-y0) / hy
!    l = idint(s4 + 1.d0)
!    l1 = idint(s4 + 1.5d0)
!    s3 = l - s4
!    s4 = l1 - 0.5d0 - s4
!    s6 = (z-z0) / hz
!    k = idint(s6 + 1.d0)
!    k1 = idint(s6 + 1.5d0)
!    s5 = k - s6
!    s6 = k1 - 0.5d0 - s6
!    s11 = 1.d0 - s1
!    s21 = 1.d0 - s2
!    s31 = 1.d0 - s3
!    s41 = 1.d0 - s4
!    s51 = 1.d0 - s5
!    s61 = 1.d0 - s6
! 
!    Ex_total = s1 * (s4 * (s6 * Ex(i, l1, k1) + s61 * Ex(i, l1, k1 + 1)) +&
!	+s41 * (s6 * Ex(i, l1 + 1, k1) + s61 * Ex(i, l1 + 1, k1 + 1))) +&
!	+s11 * (s4 * (s6 * Ex(i + 1, l1, k1) + s61 * Ex(i + 1, l1, k1 + 1)) +&
!	+s41 * (s6 * Ex(i + 1, l1 + 1, k1) + s61 * Ex(i + 1, l1 + 1, k1 + 1)))
!        
!    Ey_total = s2 * (s3 * (s6 * Ey(i1, l, k1) + s61 * Ey(i1, l, k1 + 1)) +&
!	+s31 * (s6 * Ey(i1, l + 1, k1) + s61 * Ey(i1, l + 1, k1 + 1))) +&
!	+s21 * (s3 * (s6 * Ey(i1 + 1, l, k1) + s61 * Ey(i1 + 1, l, k1 + 1)) +&
!		+s31 * (s6 * Ey(i1 + 1, l + 1, k1) + s61 * Ey(i1 + 1, l + 1, k1 + 1)))
! 
!	Ez_total = s2 * (s4 * (s5 * Ez(i1, l1, k) + s51 * Ez(i1, l1, k + 1)) +&
!	+s41 * (s5 * Ez(i1, l1 + 1, k) + s51 * Ez(i1, l1 + 1, k + 1))) +&
!	+s21 * (s4 * (s5 * Ez(i1 + 1, l1, k) + s51 * Ez(i1 + 1, l1, k + 1)) +&
!		+s41 * (s5 * Ez(i1 + 1, l1 + 1, k) + s51 * Ez(i1 + 1, l1 + 1, k + 1)))
! 
!	Bx_total = s2 * (s3 * (s5 * Bx(i1, l, k) + s51 * Bx(i1, l, k + 1)) +&
!	+s31 * (s5 * Bx(i1, l + 1, k) + s51 * Bx(i1, l + 1, k + 1))) +&
!	+s21 * (s3 * (s5 * Bx(i1 + 1, l, k) + s51 * Bx(i1 + 1, l, k + 1)) +&
!		+s31 * (s5 * Bx(i1 + 1, l + 1, k) + s51 * Bx(i1 + 1, l + 1, k + 1)))
! 
!	By_total = s1 * (s4 * (s5 * By(i, l1, k) + s51 * By(i, l1, k + 1)) +&
!	+s41 * (s5 * By(i, l1 + 1, k) + s51 * By(i, l1 + 1, k + 1))) +&
!	+s11 * (s4 * (s5 * By(i + 1, l1, k) + s51 * By(i + 1, l1, k + 1)) +&
!		+s41 * (s5 * By(i + 1, l1 + 1, k) + s51 * By(i + 1, l1 + 1, k + 1)))
! 
!	Bz_total = s1 * (s3 * (s61 * Bz(i, l, k1) + s6 * Bz(i, l, k1 + 1)) +&
!	+s31 * (s61 * Bz(i, l + 1, k1) + s6 * Bz(i, l + 1, k1 + 1))) +&
!	+s11 * (s3 * (s61 * Bz(i + 1, l, k1) + s6 * Bz(i + 1, l, k1 + 1)) +&
!		+s31 * (s61 * Bz(i + 1, l + 1, k1) + s6 * Bz(i + 1, l + 1, k1 + 1)))
!   
!    ux=u_0-ht*q/2.*(Ex_total+(v_0*Bz_total-w_0*By_total)/c)
!    uy=v_0-ht*q/2.*(Ey_total+(w_0*Bx_total-u_0*Bz_total)/c)
!    uz=w_0-ht*q/2.*(Ez_total+(u_0*By_total-v_0*Bx_total)/c)
!
!    i_t=1
!    t=t0+ht
!    energy=am*(ux*ux+uy*uy+uz*uz)/2
!    !write(2, 100) t0, x, y, z, ux, uy, uz, energy
!    s=0.d0
!    p=q*ht/am/c/2.
!    time=0.
!    time_int=0.
!    do n=1,nt
!        call ELAPSE_SECOND(start)
!        s2 = (x-x0) / hx
!        i = idint(s2 + 1.d0)
!        i1 = idint(s2 + 1.5d0)
!        s1 = i - s2
!        s2 = i1 - 0.5d0 - s2
!        s4 = (y-y0) / hy
!        l = idint(s4 + 1.d0)
!        l1 = idint(s4 + 1.5d0)
!        s3 = l - s4
!        s4 = l1 - 0.5d0 - s4
!        s6 = (z-z0) / hz
!        k = idint(s6 + 1.d0)
!        k1 = idint(s6 + 1.5d0)
!        s5 = k - s6
!        s6 = k1 - 0.5d0 - s6
!        s11 = 1.d0 - s1 
!        s21 = 1.d0 - s2
!        s31 = 1.d0 - s3
!        s41 = 1.d0 - s4
!        s51 = 1.d0 - s5
!        s61 = 1.d0 - s6
!  
!  
!        Ex_total = s1 * (s4 * (s6 * Ex(i, l1, k1) + s61 * Ex(i, l1, k1 + 1)) +&
!	    +s41 * (s6 * Ex(i, l1 + 1, k1) + s61 * Ex(i, l1 + 1, k1 + 1))) +&
!	    +s11 * (s4 * (s6 * Ex(i + 1, l1, k1) + s61 * Ex(i + 1, l1, k1 + 1)) +&
!		+s41 * (s6 * Ex(i + 1, l1 + 1, k1) + s61 * Ex(i + 1, l1 + 1, k1 + 1)))
!        
!        Ey_total = s2 * (s3 * (s6 * Ey(i1, l, k1) + s61 * Ey(i1, l, k1 + 1)) +&
!		+s31 * (s6 * Ey(i1, l + 1, k1) + s61 * Ey(i1, l + 1, k1 + 1))) +&
!		+s21 * (s3 * (s6 * Ey(i1 + 1, l, k1) + s61 * Ey(i1 + 1, l, k1 + 1)) +&
!			+s31 * (s6 * Ey(i1 + 1, l + 1, k1) + s61 * Ey(i1 + 1, l + 1, k1 + 1)))
!  
!	    Ez_total = s2 * (s4 * (s5 * Ez(i1, l1, k) + s51 * Ez(i1, l1, k + 1)) +&
!		+s41 * (s5 * Ez(i1, l1 + 1, k) + s51 * Ez(i1, l1 + 1, k + 1))) +&
!		+s21 * (s4 * (s5 * Ez(i1 + 1, l1, k) + s51 * Ez(i1 + 1, l1, k + 1)) +&
!			+s41 * (s5 * Ez(i1 + 1, l1 + 1, k) + s51 * Ez(i1 + 1, l1 + 1, k + 1)))
!  
!	    Bx_total = s2 * (s3 * (s5 * Bx(i1, l, k) + s51 * Bx(i1, l, k + 1)) +&
!		+s31 * (s5 * Bx(i1, l + 1, k) + s51 * Bx(i1, l + 1, k + 1))) +&
!		+s21 * (s3 * (s5 * Bx(i1 + 1, l, k) + s51 * Bx(i1 + 1, l, k + 1)) +&
!			+s31 * (s5 * Bx(i1 + 1, l + 1, k) + s51 * Bx(i1 + 1, l + 1, k + 1)))
!  
!	    By_total = s1 * (s4 * (s5 * By(i, l1, k) + s51 * By(i, l1, k + 1)) +&
!		+s41 * (s5 * By(i, l1 + 1, k) + s51 * By(i, l1 + 1, k + 1))) +&
!		+s11 * (s4 * (s5 * By(i + 1, l1, k) + s51 * By(i + 1, l1, k + 1)) +&
!			+s41 * (s5 * By(i + 1, l1 + 1, k) + s51 * By(i + 1, l1 + 1, k + 1)))
!  
!	    Bz_total = s1 * (s3 * (s6 * Bz(i, l, k1) + s61 * Bz(i, l, k1 + 1)) +&
!		+s31 * (s6 * Bz(i, l + 1, k1) + s61 * Bz(i, l + 1, k1 + 1))) +&
!		+s11 * (s3 * (s6 * Bz(i + 1, l, k1) + s61 * Bz(i + 1, l, k1 + 1)) +&
!			+s31 * (s6 * Bz(i + 1, l + 1, k1) + s61 * Bz(i + 1, l + 1, k1 + 1)))
!        !call ELAPSE_SECOND(finish_int)
!        b_0=sqrt(Bx_total**2+By_total**2+Bz_total**2)
!        B=reshape((/Bx_total, By_total, Bz_total/), (/3, 1/))
!        E=reshape((/Ex_total, Ey_total, Ez_total/), (/3, 1/))
!        speed=reshape((/ux, uy, uz/), (/3, 1/))
!        
!        do j=1, 3
!            B(j, 1)=B(j, 1)/b_0
!        end do
!        vb=0.
!        eb=0.
!        do j=1, 3
!            vb=vb+speed(j, 1)*B(j, 1)
!            eb=eb+E(j, 1)*B(j, 1)
!        end do
!        do j=1, 3
!            v_par(j, 1)=vb*B(j, 1)
!            v_ort(j, 1)=speed(j, 1)-v_par(j, 1)
!            E_par(j, 1)=eb*B(j, 1)
!            E_ort(j, 1)=(E(j, 1)-E_par(j, 1))*c/b_0
!        end do
!        theta=q*b_0/am/c
!        cs=cos(theta*ht)
!        cs1=1-cs
!        sn=sin(theta*ht)
!        
!        
!        ux=E_ort(2, 1)*B(3, 1)*cs1-E_ort(3, 1)*B(2, 1)*cs1+E_ort(1, 1)*sn+v_ort(1, 1)*cs+(v_ort(2, 1)*B(3, 1)-v_ort(3, 1)*B(2, 1))*sn+v_par(1, 1)
!        uy=E_ort(3, 1)*B(1, 1)*cs1-E_ort(1, 1)*B(3, 1)*cs1+E_ort(2, 1)*sn+v_ort(2, 1)*cs+(v_ort(3, 1)*B(1, 1)-v_ort(1, 1)*B(3, 1))*sn+v_par(2, 1)
!        uz=E_ort(1, 1)*B(2, 1)*cs1-E_ort(2, 1)*B(1, 1)*cs1+E_ort(3, 1)*sn+v_ort(3, 1)*cs+(v_ort(1, 1)*B(2, 1)-v_ort(2, 1)*B(1, 1))*sn+v_par(3, 1)
!
!        x=x+ht*ux
!        y=y+ht*uy
!        z=z+ht*uz
!
!        call ELAPSE_SECOND(finish)
!      !  write(*,*) nt, finish, start, s
!        energy=am*(ux**2+uy**2+uz**2)/2.
!        time=time+finish-start
!       ! time_int=time_int+finish_int-start
!        !write(2, 100) t, x, y, z, ux, uy, uz, energy
!        t=i_t*ht+t0
!        i_t=i_t+1
!    end do
!    bor=time
!    end function bor
    


