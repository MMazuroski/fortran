!  CFD_Prova1_Questao1_MarinaMazuroski.f90 
!
!  FUNCTIONS:
!  CFD_Prova1_Questao1_MarinaMazuroski - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: CFD_Prova1_Questao1_MarinaMazuroski
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program CFD_Prova1_Questao1_MarinaMazuroski

    implicit none
    
    integer,parameter:: nv=5 !Volumes de controle
    integer:: i
    real:: L, dx, dt, f, tempo, tempof, pi
    real, dimension(nv):: Ap, Apo, Aw, Ae, Su, Sp
    real, dimension(nv):: T, x, Tpo, a, b, V
    real:: ro, rL, T0, h, Tinf, Cp, k, rho
    
    f=1
    dt=1
    pi=3.1415
    
    L=0.1
    ro=0.011
    rL=0.001
    T0=100
    h=2500
    Tinf=20
    Cp=466.2
    k=57.7
    rho=7840

    dx=L/nv
    write(*,*) 'Tempo final'
    read(*,*) tempof
    
    !Malha
    x(1)=0.5D0*dx
    do i=2,nv
        x(i)=x(i-1)+dx
    end do
    
    !Campo inicial de T
    do i=1,nv
        T(i)=T0
        Tpo(i)=T0
    end do
    
    tempo=0
    do while (tempo<tempof)   
    
    !Parametros
    do i=1,nv
        a(i)=x(i)-dx/2.0D0
        A(i)=((rL-ro)*a(i)/L+ro)
        b(i)=x(i)+dx/2.0D0
        B(i)=((rL-ro)*b(i)/L+ro)
        V(i)=pi*dx/3*(A(i)*A(i)+A(i)*B(i)+B(i)*B(i))
        if(i==1) then
            Aw(i)=0
            Ae(i)=k*pi*B(i)**2.0D0/dx
            Apo(i)=rho*Cp*V(i)/dt
            Su(i)=(2*k*pi*A(i)**2.0D0/dx)*T0+Apo(i)*Tpo(i)
            Sp(i)=-(2*k*pi*A(i)**2.0D0/dx)
            Ap(i)=Apo(i)+f*(Aw(i)+Ae(i))-Sp(i)
        else if((i>1).and.(i<nv)) then
            Aw(i)=k*pi*A(i)**2.0D0/dx
            Ae(i)=k*pi*B(i)**2.0D0/dx
            Apo(i)=rho*Cp*V(i)/dt
            Su(i)=Apo(i)*T(i)
            Sp(i)=0
            Ap(i)=Apo(i)+f*(Aw(i)+Ae(i))-Sp(i)
        else
            Aw(i)=k*pi*A(i)**2.0D0/dx
            Ae(i)=0
            Apo(i)=rho*Cp*V(i)/dt
            Su(i)=pi*B(i)**2.0D0*Tinf/(1.0D0/h+dx/(2.0D0*k))+Apo(i)*Tpo(i)
            Sp(i)=-pi*B(i)**2.0D0/(1.0D0/h+dx/(2.0D0*k))
            Ap(i)=Apo(i)+f*(Aw(i)+Ae(i))-Sp(i)
        end if
    end do
    
      
        call tdma1D(Aw, Ae, Su, Sp, Ap, nv, T)
        tempo=tempo+dt
        write(*,*) 'Tempo', tempo, 'segundos'
        write(*,*) 'Temperatura'
        write(*,*) T
        do i=1,nv
            Tpo(i)=T(i)
        end do
    end do
    
    

    pause
    !Tecplot
    150 format (2x,30(F20.12,1X))
    open(15,file="Transiente.plt")
    write(15,*) 'Title="Transiente"'
    write(15,*) "Variables='x','T'"
    write(15,*) 'ZONE T="Transiente" i=',nv
    
    do i=1,nv
        write(15,150) x(i), T(i)
    end do
    close(15)
    
    
    Contains
    SUBROUTINE tdma1D(Aw, Ae, Su, Sp, Ap, nv, T)
    
    integer:: i, k, it_final
    integer:: nv
    real, dimension(nv):: Ap, Aw, Ae, Su, Sp
    real, dimension(nv):: Am, Bm, Cm, Dm, P, Q
    real, dimension(nv):: T
    real, dimension(1:nv,1:2):: M

    do i=1,nv
        if (i==1) then
            Am(i)=Ap(i)
            Bm(i)=-Ae(i)
            Cm(i)=0
            Dm(i)=Su(i)
            P(i)=-Bm(i)/Am(i)
            Q(i)=Dm(i)/Am(i)
        else
            Am(i)=Ap(i)
            Bm(i)=-Ae(i)
            Cm(i)=-Aw(i)
            Dm(i)=Su(i)
            P(i)=-Bm(i)/(Am(i)+Cm(i)*P(i-1))
            Q(i)=(Dm(i)-Cm(i)*Q(i-1))/(Am(i)+Cm(i)*P(i-1))
        end if
    end do
    T(nv)=Q(nv)
    do i=nv,2,-1
        T(i-1)=P(i-1)*T(i)+Q(i-1)
    end do
        
    end subroutine tdma1D


    end program CFD_Prova1_Questao1_MarinaMazuroski

