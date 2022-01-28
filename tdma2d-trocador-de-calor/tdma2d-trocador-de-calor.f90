!  CFD_Trabalho1_Verificacao_MarinaMazuroski_2.f90 
!
!  FUNCTIONS:
!  CFD_Trabalho1_Verificacao_MarinaMazuroski_2 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: CFD_Trabalho1_Verificacao_MarinaMazuroski_2
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program CFD_Trabalho1_Verificacao_MarinaMazuroski_2

    implicit none

    integer, parameter:: nvi=70    
    integer, parameter:: nvj=70
    integer:: i, j, it
    real:: L, W, dx, dy, T_est
    real:: residuo, RMS, tol, tempo, tempof, dt
    real:: h, Tinf, rho, cp, k, alfa, Tini
    real, dimension(nvi, nvj):: Ap, Aw, Ae, An, As, Su, Sp, Apo
    real, dimension(nvi, nvj):: T, Tpo
    real, dimension(nvi):: x, R_final
    real, dimension(nvj):: y
    
    !h=10.0D0
    h=100.0D0
    !Tinf=0
    Tinf=0.0D0
    !k=202.4D0
    k=387.6D0
    !rho=2719.0D0
    rho=8978.0D0
    !cp=871.0D0
    cp=381.0D0
    alfa=k/(rho*cp)
    !Tini=50
    Tini=500.0D0
    
    !L=1.0D0
    L=0.05D0
    !W=1.0D0
    W=0.05D0
    dx=L/nvi
    dy=W/nvj
    
    dt=1
    it=0.0D0
    residuo=0.0D0
    RMS=0.0D0
    tol=9.0D-03
    
    !Malha
    x(1)=0.5D0*dx
    do i=2,nvi
        x(i)=x(i-1)+dx
    end do
    
    y(1)=0.5D0*dy
    do j=2, nvj
        y(j)=y(j-1)+dy
    end do
    
    !Campo inicial de T
    do j=1, nvj
        do i=1, nvi
            Tpo(i,j)=Tini
            T(i,j)=Tini
        end do
    end do
    !write(*,*) 'Campo inicial'
    !write(*,*) T
    !pause
    
    !Apo=0   
    write(*,*) 'Tempo final?'
    read(*,*) tempof
    tempo=0
    
    do while (tempo<tempof)
    
    !Parametros
    do i=1,nvi !Passagem de colunas    
        if(i==1) then !Primeira coluna
            do j=1,nvj !Elementos da primeira coluna
                
                if(j==1) then !Primeiro elemento da primeira coluna
                    Aw(i,j)= 0
                    Ae(i,j)= dy/dx
                    An(i,j)= dx/dy
                    As(i,j)= 0
                    Apo(i,j)= dx*dy/(alfa*dt)
                    Su(i,j)= Apo(i,j)*Tpo(i,j)
                    Sp(i,j)= 0
                    Ap(i,j)= Apo(i,j) + Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
                    !write(*,*) 'Primeiro elemento da primeira coluna'
                    !write(*,*) 'Aw  Ae  An  As'
                    !write(*,*) Aw(i,j), Ae(i,j), An(i,j), As(i,j)
                    !write(*,*) 'Su  Sp  Ap'
                    !write(*,*) Su(i,j), Sp(i,j), Ap(i,j)
                    !pause
        
                else if((j>1).and.(j<nvj)) then !Elementos intermediários da primeira coluna
                    Aw(i,j)= 0
                    Ae(i,j)= dy/dx
                    An(i,j)= dx/dy
                    As(i,j)= dx/dy
                    Apo(i,j)= dx*dy/(alfa*dt)
                    Su(i,j)= Apo(i,j)*Tpo(i,j)
                    Sp(i,j)= 0
                    Ap(i,j)= Apo(i,j) + Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
                    !write(*,*) 'Elementos intermediarios da primeira coluna'
                    !write(*,*) 'Aw  Ae  An  As'
                    !write(*,*) Aw(i,j), Ae(i,j), An(i,j), As(i,j)
                    !write(*,*) 'Su  Sp  Ap'
                    !write(*,*) Su(i,j), Sp(i,j), Ap(i,j)
                    !pause   
            
                else !Último elemento da primeira coluna
                    Aw(i,j)= 0
                    Ae(i,j)= dy/dx
                    An(i,j)= 0
                    As(i,j)= dx/dy
                    Apo(i,j)= dx*dy/(alfa*dt)
                    Su(i,j)= (dx/k)*Tinf/(1.0D0/h+dy/(2.0D0*k)) + Apo(i,j)*Tpo(i,j)
                    Sp(i,j)= -(dx/k)*1.0D0/(1.0D0/h+dy/(2.0D0*k))
                    Ap(i,j)= Apo(i,j) + Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
                    !write(*,*) 'Ultimo elemento da primeira coluna'
                    !write(*,*) 'Aw  Ae  An  As'
                    !write(*,*) Aw(i,j), Ae(i,j), An(i,j), As(i,j)
                    !write(*,*) 'Su  Sp  Ap'
                    !write(*,*) Su(i,j), Sp(i,j), Ap(i,j)
                    !pause
                end if
                
            end do !Encerra primeira coluna
                    
        else if ((i>1).and.(i<nvi)) then !Colunas intermediárias
            do j=1,nvj !Elementos das colunas intermediárias
                
                if (j==1) then !Primeira linha da coluna intermediária
                    Aw(i,j)= dy/dx
                    Ae(i,j)= dy/dx
                    An(i,j)= dx/dy
                    As(i,j)= 0
                    Apo(i,j)= dx*dy/(alfa*dt)
                    Su(i,j)= Apo(i,j)*Tpo(i,j)
                    Sp(i,j)= 0
                    Ap(i,j)= Apo(i,j) + Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
                    !write(*,*) 'Primeira linha da coluna intermediaria'
                    !write(*,*) 'Aw  Ae  An  As'
                    !write(*,*) Aw(i,j), Ae(i,j), An(i,j), As(i,j)
                    !write(*,*) 'Su  Sp  Ap'
                    !write(*,*) Su(i,j), Sp(i,j), Ap(i,j)
                    !pause

                else if ((j>1).and.(j<nvj)) then !Nós centrais
                    Aw(i,j)= dy/dx
                    Ae(i,j)= dy/dx
                    An(i,j)= dx/dy
                    As(i,j)= dx/dy
                    Apo(i,j)= dx*dy/(alfa*dt)
                    Su(i,j)= Apo(i,j)*Tpo(i,j)
                    Sp(i,j)= 0
                    Ap(i,j)= Apo(i,j) + Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
                    !write(*,*) 'No central'
                    !write(*,*) 'Aw  Ae  An  As'
                    !write(*,*) Aw(i,j), Ae(i,j), An(i,j), As(i,j)
                    !write(*,*) 'Su  Sp  Ap'
                    !write(*,*) Su(i,j), Sp(i,j), Ap(i,j)
                    !pause
    
                else !Última linha das colunas intermediárias
                    Aw(i,j)= dy/dx
                    Ae(i,j)= dy/dx
                    An(i,j)= 0
                    As(i,j)= dx/dy
                    Apo(i,j)= dx*dy/(alfa*dt)
                    Su(i,j)= (dx/k)*Tinf/(1.0D0/h+dy/(2.0D0*k)) + Apo(i,j)*Tpo(i,j)
                    Sp(i,j)= -(dx/k)*1.0D0/(1.0D0/h+dy/(2.0D0*k))
                    Ap(i,j)= Apo(i,j) + Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
                    !write(*,*) 'Ultima linha da coluna intermediaria'
                    !write(*,*) 'Aw  Ae  An  As'
                    !write(*,*) Aw(i,j), Ae(i,j), An(i,j), As(i,j)
                    !write(*,*) 'Su  Sp  Ap'
                    !write(*,*) Su(i,j), Sp(i,j), Ap(i,j)
                    !pause
                
                end if
                
            end do !Encerra colunas intermediárias
            
        else !Última coluna
            do j=1,nvj
                if (j==1) then !Primeiro elemento da última coluna
                    Aw(i,j)= dy/dx
                    Ae(i,j)= 0
                    An(i,j)= dx/dy
                    As(i,j)= 0
                    Apo(i,j)= dx*dy/(alfa*dt)
                    Su(i,j)= (dy/k)*Tinf/(1.0D0/h+dx/(2.0D0*k)) + Apo(i,j)*Tpo(i,j)
                    Sp(i,j)= -(dy/k)*1.0D0/(1.0D0/h+dx/(2.0D0*k))
                    Ap(i,j)= Apo(i,j) + Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
                    !write(*,*) 'Primeiro elemento da ultima coluna'
                    !write(*,*) 'Aw  Ae  An  As'
                    !write(*,*) Aw(i,j), Ae(i,j), An(i,j), As(i,j)
                    !write(*,*) 'Su  Sp  Ap'
                    !write(*,*) Su(i,j), Sp(i,j), Ap(i,j)
                    !pause
                
                else if ((j>1).and.(j<nvj)) then !Elementos intermediários da última coluna
                    Aw(i,j)= dy/dx
                    Ae(i,j)= 0
                    An(i,j)= dx/dy
                    As(i,j)= dx/dy
                    Apo(i,j)= dx*dy/(alfa*dt)
                    Su(i,j)= (dy/k)*Tinf/(1.0D0/h+dx/(2.0D0*k)) + Apo(i,j)*Tpo(i,j)
                    Sp(i,j)= -(dy/k)*1.0D0/(1.0D0/h+dx/(2.0D0*k))
                    Ap(i,j)= Apo(i,j) + Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
                    !write(*,*) 'Elemento intermediario da ultima coluna'
                    !write(*,*) 'Aw  Ae  An  As'
                    !write(*,*) Aw(i,j), Ae(i,j), An(i,j), As(i,j)
                    !write(*,*) 'Su  Sp  Ap'
                    !write(*,*) Su(i,j), Sp(i,j), Ap(i,j)
                    !pause
                    
                else !Último elemento da última coluna
                    Aw(i,j)= dy/dx
                    Ae(i,j)= 0
                    An(i,j)= 0
                    As(i,j)= dx/dy
                    Apo(i,j)= dx*dy/(alfa*dt)
                    Su(i,j)= (dx/k)*Tinf/(1.0D0/h+dy/(2.0D0*k)) + (dy/k)*Tinf/(1.0D0/h+dx/(2.0D0*k)) + Apo(i,j)*Tpo(i,j)
                    Sp(i,j)= -(dx/k)*1.0D0/(1.0D0/h+dy/(2.0D0*k)) + (dy/k)*1.0D0/(1.0D0/h+dx/(2.0D0*k))
                    Ap(i,j)= Apo(i,j) + Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
                    !write(*,*) 'Ultimo elemento da ultima coluna'
                    !write(*,*) 'Aw  Ae  An  As'
                    !write(*,*) Aw(i,j), Ae(i,j), An(i,j), As(i,j)
                    !write(*,*) 'Su  Sp  Ap'
                    !write(*,*) Su(i,j), Sp(i,j), Ap(i,j)
                    !pause
                    
                end if
            end do !Encerra última coluna
        end if !Encerra colunas
    end do !Encerra linhas 
    
    !Chamada TDMA 2D
    
    call tdma2d (Aw, Ae, An, As, Ap, Su, Sp, nvi, nvj, T, it, residuo, RMS, tol)
    tempo=tempo+dt
    write(*,*) 'Tempo:', tempo, 'segundos'
    write(*,*) 'Temperatura'
    write(*,*) T

    
    do j=1, nvj
        do i=1, nvi
            Tpo(i,j)=T(i,j)
        end do
    end do
    
    end do
    pause
    
    do i=1,nvi
        R_final(i)=T(i, nvj/2)
    end do
    
    !Resultados plotados no TecPlot:
    !Temperatura através de MVF:
    150 format(2x,30(F20.12,1X))
    open(15,file="Verificacao_100s_70.plt")
    write (15,*) 'Title ="Verificacao"'
    write (15,*) "Variables = 'x', 'T'"
    write (15,*) 'ZONE T = "CAMPO" i=',nvi
    do i=1,nvi
        write(15,150) x(i), R_final(i)
    end do
    close(15)
    
    !Campo 3D
    20 format(2x,30(F20.12,1X))
    open(2,file="CampoT_3D.plt")
    write(2,*) 'Title="CampoT"'
    write(2,*) "Variables='x','y','T'"
    write(2,*) 'ZONE T = "CAMPO" j=',nvj,' i=',nvi
    do i=1,nvi
    do j=1,nvj
        write(2,20) x(j), y(i), T(j,i)
    end do
    end do
    close(2)
    
    
    Contains
    SUBROUTINE tdma2d(Aw, Ae, An, As, Ap, Su, Sp, nvi, nvj, T, it, residuo, RMS, tol)
    
    integer:: nvi, nvj
    integer:: i, j, it
    real:: residuo, RMS, tol
    real, dimension (nvi, nvj):: Aw, Ae, An, As, Ap, Su, Sp
    real, dimension (nvi, nvj):: T
    real, dimension (nvi, nvj):: Am, Bm, Cm, Dm, P, Q, SOL
    
    160format(2x,30(F20.12,1X))
        open(16,file="R.plt")
        write(16,*) 'Title ="RMS"'
        write(16,*) "Variables= 'it', ‘RMS'"    
    
    !Varredura por linha, de baixo para cima
60     do j=1,nvj
        if (j==1) then !Primeira linha
            do i=1,nvi !Subindo nas linhas
                Am(i,j)= Ap(i,j)
                Bm(i,j)= -Ae(i,j)
                Cm(i,j)= -Aw(i,j)
                Dm(i,j)= An(i,j)*T(i,j+1) + Su(i,j)
                if (i==1) then
                    P(i,j)= -Bm(i,j) / Am(i,j)
                    Q(i,j)= Dm(i,j) / Am(i,j)
                else 
                    P(i,j)= -Bm(i,j) / (Am(i,j) + Cm(i,j)*P(i-1,j))
                    Q(i,j)= (Dm(i,j) - Cm(i,j)*Q(i-1,j)) / (Am(i,j) + Cm(i,j)*P(i-1,j))
                end if
            end do
           
        else if ((j>1).and.(j<nvj)) then !Linhas intermediárias 
            do i=1,nvi !Subindo nas linhas
                Am(i,j)= Ap(i,j)
                Bm(i,j)= -Ae(i,j)
                Cm(i,j)= -Aw(i,j)
                Dm(i,j)= An(i,j)*T(i,j+1) + As(i,j)*T(i,j-1) + Su(i,j)
                if (i==1) then
                    P(i,j)= -Bm(i,j) / Am(i,j)
                    Q(i,j)= Dm(i,j) / Am(i,j)
                else 
                    P(i,j)= -Bm(i,j) / (Am(i,j) + Cm(i,j)*P(i-1,j))
                    Q(i,j)= (Dm(i,j) - Cm(i,j)*Q(i-1,j)) / (Am(i,j) + Cm(i,j)*P(i-1,j))
                end if
            end do
            
        else !Última linha
            do i=1,nvi !Subindo nas linhas
                Am(i,j)= Ap(i,j)
                Bm(i,j)= -Ae(i,j)
                Cm(i,j)= -Aw(i,j)
                Dm(i,j)= As(i,j)*T(i,j-1) + Su(i,j)
                if (i==1) then
                    P(i,j)= -Bm(i,j) / Am(i,j)
                    Q(i,j)= Dm(i,j) / Am(i,j)
                else 
                    P(i,j)= -Bm(i,j) / (Am(i,j) + Cm(i,j)*P(i-1,j))
                    Q(i,j)= (Dm(i,j) - Cm(i,j)*Q(i-1,j)) / (Am(i,j) + Cm(i,j)*P(i-1,j))
                end if
            end do
        end if
            
        !Cálculo decrescente de T
        SOL(nvi,j)=Q(nvi,j)
        do i=nvi,2,-1
            SOL(i-1,j)= P(i-1,j)*SOL(i,j) + Q(i-1,j)
        end do
    end do
   
    T(:,:)=SOL(:,:)
    !Fim da varredura por linha de baixo para cima
   
    
    
    !Varredura por linha, de cima para baixo
     do j=nvj,1,-1
        if (j==1) then !Primeira linha
            do i=1,nvi !Subindo nas linhas
                Am(i,j)= Ap(i,j)
                Bm(i,j)= -Ae(i,j)
                Cm(i,j)= -Aw(i,j)
                Dm(i,j)= An(i,j)*T(i,j+1) + Su(i,j)
                if (i==1) then
                    P(i,j)= -Bm(i,j) / Am(i,j)
                    Q(i,j)= Dm(i,j) / Am(i,j)
                else 
                    P(i,j)= -Bm(i,j) / (Am(i,j) + Cm(i,j)*P(i-1,j))
                    Q(i,j)= (Dm(i,j) - Cm(i,j)*Q(i-1,j)) / (Am(i,j) + Cm(i,j)*P(i-1,j))
                end if
            end do
           
        else if ((j>1).and.(j<nvj)) then !Linhas intermediárias 
            do i=1,nvi !Subindo nas linhas
                Am(i,j)= Ap(i,j)
                Bm(i,j)= -Ae(i,j)
                Cm(i,j)= -Aw(i,j)
                Dm(i,j)= An(i,j)*T(i,j+1) + As(i,j)*T(i,j-1) + Su(i,j)
                if (i==1) then
                    P(i,j)= -Bm(i,j) / Am(i,j)
                    Q(i,j)= Dm(i,j) / Am(i,j)
                else 
                    P(i,j)= -Bm(i,j) / (Am(i,j) + Cm(i,j)*P(i-1,j))
                    Q(i,j)= (Dm(i,j) - Cm(i,j)*Q(i-1,j)) / (Am(i,j) + Cm(i,j)*P(i-1,j))
                end if
            end do
            
        else !Última linha
            do i=1,nvi !Subindo nas linhas
                Am(i,j)= Ap(i,j)
                Bm(i,j)= -Ae(i,j)
                Cm(i,j)= -Aw(i,j)
                Dm(i,j)= As(i,j)*T(i,j-1) + Su(i,j)
                if (i==1) then
                    P(i,j)= -Bm(i,j) / Am(i,j)
                    Q(i,j)= Dm(i,j) / Am(i,j)
                else 
                    P(i,j)= -Bm(i,j) / (Am(i,j) + Cm(i,j)*P(i-1,j))
                    Q(i,j)= (Dm(i,j) - Cm(i,j)*Q(i-1,j)) / (Am(i,j) + Cm(i,j)*P(i-1,j))
                end if
            end do
        end if
            
        !Cálculo decrescente de T
        SOL(nvi,j)=Q(nvi,j)
        do i=nvi,2,-1
            SOL(i-1,j)= P(i-1,j)*SOL(i,j) + Q(i-1,j)
        end do
    end do
   
    T(:,:)=SOL(:,:) 
    !Fim da varredura por linha de cima para baixo

    
    
    !Varredura por coluna, da esquerda para a direita
    do i=1,nvi
        if (i==1) then !Primeira coluna
            do j=1,nvj !Subindo na coluna
                Am(i,j)= Ap(i,j)
                Bm(i,j)= -An(i,j)
                Cm(i,j)= -As(i,j)
                Dm(i,j)= Ae(i,j)*T(i+1,j) + Su(i,j)
                if (j==1) then
                    P(i,j)= -Bm(i,j) / Am(i,j)
                    Q(i,j)= Dm(i,j) / Am(i,j)
                else 
                    P(i,j)= -Bm(i,j) / (Am(i,j) + Cm(i,j)*P(i,j-1))
                    Q(i,j)= (Dm(i,j) - Cm(i,j)*Q(i,j-1)) / (Am(i,j) + Cm(i,j)*P(i,j-1))
                end if
            end do
           
        else if ((i>1).and.(i<nvi)) then !Colunas intermediárias 
            do j=1,nvj
                Am(i,j)= Ap(i,j)
                Bm(i,j)= -An(i,j)
                Cm(i,j)= -As(i,j)
                Dm(i,j)= Ae(i,j)*T(i+1,j) + Aw(i,j)*T(i-1,j) + Su(i,j)
                if (j==1) then
                    P(i,j)= -Bm(i,j) / Am(i,j)
                    Q(i,j)= Dm(i,j) / Am(i,j)
                else 
                    P(i,j)= -Bm(i,j) / (Am(i,j) + Cm(i,j)*P(i,j-1))
                    Q(i,j)= (Dm(i,j) - Cm(i,j)*Q(i,j-1)) / (Am(i,j) + Cm(i,j)*P(i,j-1))
                end if
            end do
            
        else !Última coluna
            do j=1,nvj
                Am(i,j)= Ap(i,j)
                Bm(i,j)= -An(i,j)
                Cm(i,j)= -As(i,j)
                Dm(i,j)= Aw(i,j)*T(i-1,j) + Su(i,j)
                if (j==1) then
                    P(i,j)= -Bm(i,j) / Am(i,j)
                    Q(i,j)= Dm(i,j) / Am(i,j)
                else 
                    P(i,j)= -Bm(i,j) / (Am(i,j) + Cm(i,j)*P(i,j-1))
                    Q(i,j)= (Dm(i,j) - Cm(i,j)*Q(i,j-1)) / (Am(i,j) + Cm(i,j)*P(i,j-1))
                end if
            end do
        end if
            
        !Cálculo decrescente de T
        SOL(i,nvj)=Q(i,nvj)
        do j=nvj,2,-1
            SOL(i,j-1)= P(i,j-1)*SOL(i,j) + Q(i,j-1)
        end do
    end do
   
    T(:,:)=SOL(:,:)
    
    !Fim da varredura por coluna da esquerda para a direita

    
    
    
    !Varredura por coluna, da direita para a esquerda
    do i=nvi,1,-1
        if (i==1) then !Primeira coluna
            do j=1,nvj !Subindo na coluna
                Am(i,j)= Ap(i,j)
                Bm(i,j)= -An(i,j)
                Cm(i,j)= -As(i,j)
                Dm(i,j)= Ae(i,j)*T(i+1,j) + Su(i,j)
                if (j==1) then
                    P(i,j)= -Bm(i,j) / Am(i,j)
                    Q(i,j)= Dm(i,j) / Am(i,j)
                else 
                    P(i,j)= -Bm(i,j) / (Am(i,j) + Cm(i,j)*P(i,j-1))
                    Q(i,j)= (Dm(i,j) - Cm(i,j)*Q(i,j-1)) / (Am(i,j) + Cm(i,j)*P(i,j-1))
                end if
            end do
           
        else if ((i>1).and.(i<nvi)) then !Colunas intermediárias 
            do j=1,nvj
                Am(i,j)= Ap(i,j)
                Bm(i,j)= -An(i,j)
                Cm(i,j)= -As(i,j)
                Dm(i,j)= Ae(i,j)*T(i+1,j) + Aw(i,j)*T(i-1,j) + Su(i,j)
                if (j==1) then
                    P(i,j)= -Bm(i,j) / Am(i,j)
                    Q(i,j)= Dm(i,j) / Am(i,j)
                else 
                    P(i,j)= -Bm(i,j) / (Am(i,j) + Cm(i,j)*P(i,j-1))
                    Q(i,j)= (Dm(i,j) - Cm(i,j)*Q(i,j-1)) / (Am(i,j) + Cm(i,j)*P(i,j-1))
                end if
            end do
            
        else !Última coluna
            do j=1,nvj
                Am(i,j)= Ap(i,j)
                Bm(i,j)= -An(i,j)
                Cm(i,j)= -As(i,j)
                Dm(i,j)= Aw(i,j)*T(i-1,j) + Su(i,j)
                if (j==1) then
                    P(i,j)= -Bm(i,j) / Am(i,j)
                    Q(i,j)= Dm(i,j) / Am(i,j)
                else 
                    P(i,j)= -Bm(i,j) / (Am(i,j) + Cm(i,j)*P(i,j-1))
                    Q(i,j)= (Dm(i,j) - Cm(i,j)*Q(i,j-1)) / (Am(i,j) + Cm(i,j)*P(i,j-1))
                end if
            end do
        end if
            
        !Cálculo decrescente de T
        SOL(i,nvj)=Q(i,nvj)
        do j=nvj,2,-1
            SOL(i,j-1)= P(i,j-1)*SOL(i,j) + Q(i,j-1)
        end do
    end do
   
    T(:,:)=SOL(:,:)
    
    !Fim da varredura por coluna da direita para a esquerda

    
    
    !Resíduo
    RMS=0.0D0
    do j=1,nvj
    !write(*,*) 'Residuo'
        if (j==1) then !Primeira linha
            do i=1,nvi
                if (i==1) then
                    residuo= Ae(i,j)*SOL(i+1,j) + An(i,j)*SOL(i,j+1) + Su(i,j) - Ap(i,j)*SOL(i,j)
                    !write(*,*) residuo
                else if ((i>1).and.(i<nvi)) then
                    residuo= Aw(i,j)*SOL(i-1,j) + Ae(i,j)*SOL(i+1,j) + An(i,j)*SOL(i,j+1) + Su(i,j) - Ap(i,j)*SOL(i,j)
                    !write(*,*) residuo
                else
                    residuo= Aw(i,j)*SOL(i-1,j) + An(i,j)*SOL(i,j+1) + Su(i,j) - Ap(i,j)*SOL(i,j)
                    !write(*,*) residuo
                end if
            RMS=RMS+(residuo**2.0D0)
            end do
            
        else if ((j>1).and.(j<nvj)) then
            do i=1,nvi
                if (i==1) then
                    residuo= Ae(i,j)*SOL(i+1,j) + An(i,j)*SOL(i,j+1) + As(i,j)*SOL(i,j-1) + Su(i,j) - Ap(i,j)*SOL(i,j)
                    !write(*,*) residuo
                else if ((i>1).and.(i<nvi)) then
                    residuo= Aw(i,j)*SOL(i-1,j) + Ae(i,j)*SOL(i+1,j) + An(i,j)*SOL(i,j+1) + As(i,j)*SOL(i,j-1) + Su(i,j) - Ap(i,j)*SOL(i,j)
                    !write(*,*) residuo
                else
                    residuo= Aw(i,j)*SOL(i-1,j) + An(i,j)*SOL(i,j+1) + As(i,j)*SOL(i,j-1)+ Su(i,j) - Ap(i,j)*SOL(i,j)
                    !write(*,*) residuo
                end if
            RMS=RMS+(residuo**2.0D0)    
            end do
            
        else 
            do i=1,nvi
                if (i==1) then
                    residuo= Ae(i,j)*SOL(i+1,j) + As(i,j)*SOL(i,j-1) + Su(i,j) - Ap(i,j)*SOL(i,j)
                    !write(*,*) residuo
                else if ((i>1).and.(i<nvi)) then
                    residuo= Aw(i,j)*SOL(i-1,j) + Ae(i,j)*SOL(i+1,j) + As(i,j)*SOL(i,j-1) + Su(i,j) - Ap(i,j)*SOL(i,j)
                    !write(*,*) residuo
                else
                    residuo= Aw(i,j)*SOL(i-1,j) + As(i,j)*SOL(i,j-1)+ Su(i,j) - Ap(i,j)*SOL(i,j)
                    !write(*,*) residuo
                end if
            RMS=RMS+(residuo**2.0D0)
            end do
        end if
    end do
    RMS=sqrt(RMS)
    write(*,*) 'RMS e it'
    write(*,*) RMS, it
    
    write(16,*) it, RMS
    
    if (RMS<=tol) then
        write(*,*) 'Convergiu'
    else
        it=it+1
        go to 60
    end if
    
    end SUBROUTINE tdma2d            

    end program CFD_Trabalho1_Verificacao_MarinaMazuroski_2

