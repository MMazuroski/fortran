!  CFD_Prova1_Questao2_MarinaMazuroski.f90 
!
!  FUNCTIONS:
!  CFD_Prova1_Questao2_MarinaMazuroski - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: CFD_Prova1_Questao2_MarinaMazuroski
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program CFD_Prova1_Questao2_MarinaMazuroski

    implicit none

    integer, parameter:: nvi=4
    integer, parameter:: nvj=4
    integer:: i, j, kk
    real:: L, Ro, u, Tent, q, Cp, k, rho, dx, dr, pi, Few, Dns
    real:: residuo, RMS, tol
    real, dimension(nvi, nvj):: Ap, Aw, Ae, An, As, Su, Sp
    real, dimension(nvi, nvj):: T, c, d, Areaew, Arean, Areas
    real, dimension(nvi):: x
    real, dimension(nvj):: r
    
    L=1.6
    Ro=0.16
    u=0.0063
    Tent=20
    q=-1500
    Cp=4183
    k=0.586
    rho=998
        
    dx=L/nvi
    dr=Ro/nvj
    pi=3.1415
    Few=rho*u
    Dns=k/(Cp*dr)
    
    kk=0.0D0
    residuo=0.0D0
    RMS=0.0D0
    tol=9.0D-03
    
    !Malha
    x(1)=0.5D0*dx
    do i=2,nvi
        x(i)=x(i-1)+dx
    end do
    
    r(1)=0.5D0*dr
    do j=2, nvj
        r(j)=r(j-1)+dr
    end do
    
    !Campo inicial de T
    do j=1, nvj
        do i=1, nvi
            T(i,j)=Tent
        end do
    end do
    
    !Parametros
    
do i=1,nvi
    do j=1,nvj
        c(i,j)=r(j)-dr/2
        d(i,j)=r(j)+dr/2
        Areaew(i,j)=pi*(d(i,j)*d(i,j)-c(i,j)*c(i,j))
        Arean(i,j)=2*pi*d(i,j)*dx
        Areas(i,j)=2*pi*c(i,j)*dx
    end do
end do

    do i=1,nvi !Passagem de colunas    
        if(i==1) then !Primeira coluna
            do j=1,nvj !Elementos da primeira coluna
                if(j==1) then !Primeiro elemento da primeira coluna
                    Aw(i,j)= 0
                    Ae(i,j)= 0
                    An(i,j)= Dns*Arean(i,j)
                    As(i,j)= 0
                    Su(i,j)= Few*Areaew(i,j)*Tent
                    Sp(i,j)= -(Few*Areaew(i,j))
                    Ap(i,j)= Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
                    !write(*,*) 'Primeiro elemento da primeira coluna'
                    !write(*,*) 'Aw  Ae  An  As'
                    !write(*,*) Aw(i,j), Ae(i,j), An(i,j), As(i,j)
                    !write(*,*) 'Su  Sp  Ap'
                    !write(*,*) Su(i,j), Sp(i,j), Ap(i,j)
                    !pause
        
                else if((j>1).and.(j<nvj)) then !Elementos intermediários da primeira coluna
                    Aw(i,j)= 0
                    Ae(i,j)= 0
                    An(i,j)= Dns*Arean(i,j)
                    As(i,j)= Dns*Areas(i,j)
                    Su(i,j)= Few*Areaew(i,j)*Tent
                    Sp(i,j)= -(Few*Areaew(i,j))
                    Ap(i,j)= Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
                    !write(*,*) 'Elementos intermediarios da primeira coluna'
                    !write(*,*) 'Aw  Ae  An  As'
                    !write(*,*) Aw(i,j), Ae(i,j), An(i,j), As(i,j)
                    !write(*,*) 'Su  Sp  Ap'
                    !write(*,*) Su(i,j), Sp(i,j), Ap(i,j)
                    !pause   
            
                else !Último elemento da primeira coluna
                    Aw(i,j)= 0
                    Ae(i,j)= 0
                    An(i,j)= 0
                    As(i,j)= Dns*Areas(i,j)
                    Su(i,j)= Few*Areaew(i,j)*Tent - Arean(i,j)*q/Cp
                    Sp(i,j)= -(Few*Areaew(i,j))
                    Ap(i,j)= Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
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
                    Aw(i,j)= Few*Areaew(i,j)
                    Ae(i,j)= 0
                    An(i,j)= Dns*Arean(i,j)
                    As(i,j)= 0
                    Su(i,j)= 0
                    Sp(i,j)= 0
                    Ap(i,j)= Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
                    !write(*,*) 'Primeira linha da coluna intermediaria'
                    !write(*,*) 'Aw  Ae  An  As'
                    !write(*,*) Aw(i,j), Ae(i,j), An(i,j), As(i,j)
                    !write(*,*) 'Su  Sp  Ap'
                    !write(*,*) Su(i,j), Sp(i,j), Ap(i,j)
                    !pause

                else if ((j>1).and.(j<nvj)) then !Nós centrais
                    Aw(i,j)= Few*Areaew(i,j)
                    Ae(i,j)= 0
                    An(i,j)= Dns*Arean(i,j)
                    As(i,j)= Dns*Areas(i,j)
                    Su(i,j)= 0
                    Sp(i,j)= 0
                    Ap(i,j)= Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
                    !write(*,*) 'No central'
                    !write(*,*) 'Aw  Ae  An  As'
                    !write(*,*) Aw(i,j), Ae(i,j), An(i,j), As(i,j)
                    !write(*,*) 'Su  Sp  Ap'
                    !write(*,*) Su(i,j), Sp(i,j), Ap(i,j)
                    !pause
    
                else !Última linha das colunas intermediárias

                    Aw(i,j)= Few*Areaew(i,j)
                    Ae(i,j)= 0
                    An(i,j)= 0
                    As(i,j)= Dns*Areas(i,j)
                    Su(i,j)= -Arean(i,j)*q/Cp
                    Sp(i,j)= 0
                    Ap(i,j)= Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
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
                    Aw(i,j)= Few*Areaew(i,j)
                    Ae(i,j)= 0
                    An(i,j)= Dns*Arean(i,j)
                    As(i,j)= 0
                    Su(i,j)= 0
                    Sp(i,j)= 0
                    Ap(i,j)= Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
                    !write(*,*) 'Primeiro elemento da ultima coluna'
                    !write(*,*) 'Aw  Ae  An  As'
                    !write(*,*) Aw(i,j), Ae(i,j), An(i,j), As(i,j)
                    !write(*,*) 'Su  Sp  Ap'
                    !write(*,*) Su(i,j), Sp(i,j), Ap(i,j)
                    !pause
                
                else if ((j>1).and.(j<nvj)) then !Elementos intermediários da última coluna
                    Aw(i,j)= Few*Areaew(i,j)
                    Ae(i,j)= 0
                    An(i,j)= Dns*Arean(i,j)
                    As(i,j)= Dns*Areas(i,j)
                    Sp(i,j)= 0
                    Ap(i,j)= Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
                    !write(*,*) 'Elemento intermediario da ultima coluna'
                    !write(*,*) 'Aw  Ae  An  As'
                    !write(*,*) Aw(i,j), Ae(i,j), An(i,j), As(i,j)
                    !write(*,*) 'Su  Sp  Ap'
                    !write(*,*) Su(i,j), Sp(i,j), Ap(i,j)
                    !pause
                    
                else !Último elemento da última coluna
                    Aw(i,j)= Few*Areaew(i,j)
                    Ae(i,j)= 0
                    An(i,j)= 0
                    As(i,j)= Dns*Areas(i,j)
                    Su(i,j)= -Arean(i,j)*q/Cp
                    Sp(i,j)= 0
                    Ap(i,j)= Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) - Sp(i,j)
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
    write(*,*) 'Iniciar?'
    pause
    
    call tdma2d (Aw, Ae, An, As, Ap, Su, Sp, nvi, nvj, T, kk, residuo, RMS, tol)
    write(*,*) 'Temperatura'
    write(*,*) T
    pause
    
    Open(Unit=43,File="..\Temperaturas.dat",Status="Unknown") 
        
        Write(43, Fmt='(A15)', Advance = 'Yes') "Matrix_T"
        do j=nvj,1,-1
            do i=1,nvi
                Write(43, Fmt='(1PG15.7E2)', Advance = 'No') T(i,j)
            end do
            Write(43, Fmt='(A15)', Advance = 'Yes') " "
        end do
        Write(43, Fmt='(A15)', Advance = 'Yes') " "
        
    close(43)
    
   
    20format(2x,30(F20.12,1X))
    open(2,file="CampoT_3D.plt")
    write(2,*) 'Title="CampoT"'
    write(2,*) "Variables='x','y','T'"
    write(2,*) 'ZONE T = "CAMPO" j=',nvj,' i=',nvi
    do i=1,nvi
    do j=1,nvj
    write(2,20) x(i), r(j), T(i,j)
    end do
    end do
    close(2)
    
    
    Contains
    SUBROUTINE tdma2d(Aw, Ae, An, As, Ap, Su, Sp, nvi, nvj, T, k, residuo, RMS, tol)
    
    integer:: nvi, nvj
    integer:: i, j, k
    real:: residuo, RMS, tol
    real, dimension (nvi, nvj):: Aw, Ae, An, As, Ap, Su, Sp
    real, dimension (nvi, nvj):: T
    real, dimension (nvi, nvj):: Am, Bm, Cm, Dm, P, Q, SOL
    
    160format(2x,30(F20.12,1X))
        open(16,file="R.plt")
        write(16,*) 'Title ="RMS"'
        write(16,*) "Variables= 'k', ‘RMS'"    
    
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
    write(*,*) 'RMS e k'
    write(*,*) RMS, k
    
    write(16,*) k, RMS
    
    if (RMS<=tol) then
        write(*,*) 'Convergiu'
    else if (k<100) then
        k=k+1
        go to 60
    else
        write(*,*) 'Nao convergiu'
    end if
    
    end SUBROUTINE tdma2d            


    end program CFD_Prova1_Questao2_MarinaMazuroski

