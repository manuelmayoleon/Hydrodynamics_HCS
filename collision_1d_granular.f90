!final_version.f90
PROGRAM final_version

    !!!!!!!!!!!!!!!!!!Programa para calcular las colisiones!!!!!!!!!!!
    !!  Calculo de colisiones en 1d con condiciones periodicas      !!
    !!                                                              !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*****************************************************************************
  !  Modified:
  !
  !    19 Jule 2021
  !
  !  Author:
  !
  !    Manuel Mayo León 
  !
!! AVISOS URGENTES 
!** CONSEJOS 
!???? EXPLICACIONES 
!///  tachado


implicit none



    INTEGER:: i,j,k,l,m
    REAL(kind=8),ALLOCATABLE,DIMENSION(:):: r,v !vector con posiciones y velocidades
    REAL(kind=8),ALLOCATABLE,DIMENSION(:,:)::sumv ! suma de velocidades para cada colision
    REAL(kind=8),ALLOCATABLE,DIMENSION(:):: tmp !temperaturas en las dos direcciones del espacio
    REAL(kind=8),ALLOCATABLE,DIMENSION(:,:,:):: density !densidad para tiempo t en funcion de k
    REAL(kind=8)::temp,H,longy,sigma,rho !temperaturas en "y" y en "z", altura vertical y anchura, sigma==tamaño de la particula, rho=density
    REAL(kind=8)::alpha  !! coeficiente de restitucion y velocidad que se introduce a traves de la pared 
    LOGICAL :: granular,pos_para_t
    REAL(kind=8)::tcol,colt !tiempo de colision, tiempo para comparar y tiempo inicial
    INTEGER::nstep,tray,n,iseed !numero de repeticiones que se realizan (tiempo) y numero de iteraciones  (numero de copias)
    REAL(kind=8)::rab,vab !distancias y velocidades relativas
    INTEGER,DIMENSION(2)::ni !particulas que colisionan
    REAL,ALLOCATABLE,DIMENSION(:)::colisiones!numero de colisiones, tiempos de relajacion por repeticiones
    REAL(kind=8)::bij,qij,discr,t !bij=(ri-rj)*(vi-vj), discr es el discriminante de la solucion de segundo grado, t=tiempo de colision
    REAL(kind=8),ALLOCATABLE,DIMENSION(:)::tiempos !tiempos de colision
    REAL(kind=8), parameter :: pi = 4 * atan (1.0_8)
    ! REAL(kind=8) :: num_onda,ts,gamma,lin_dens, nu,eta,kapa
    REAL(kind=8) :: num_onda
    
    !!! para deteminar el tiempo de cálculo
    REAL(kind=4):: start, finish
    character(len=10)::alfa, colpp
    
    !!notar que consideramos KT=1
    !!inicializamos variables

       
    sigma=1.0e-6
    
  
    n=1000
   
    rho=0.03d00
   
    
    ! alpha=1.00
    alpha=0.991

    


   
    ! longy=REAL(n,8)/(rho)
    longy = 10

    num_onda=2*pi/longy

    ! eta = (3/4)*(1-alpha**2.00)
    ! nu=(1+alpha)*(3*(1+8)*(1-alpha))
    ! kapa = (1)/(nu+2*eta)
    



    !?? Temperaturas cercanas al equilibrio 
    ! temp=0.0+num_onda
    temp = 100
    !** 1 colision equivale  a 10^5 pasos temporales
    nstep=10 !para 1.5*sigma
    
    !! Numero de trayectorias en las que se promedia  
    tray=1

    



!! Allocateamos 

    
    
    
    ALLOCATE(r(n),v(n),sumv(tray,nstep),tmp(nstep),colisiones(tray),tiempos(nstep))
    ALLOCATE(density(tray,nstep,2))

    

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '            MD 1D SIMULATION                 '
    write ( *, '(a)' ) '            FORTRAN90 version                '
    
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Temperature  = ', temp
    
    write ( *, '(a,i10)' ) '  number of steps = ', nstep
    write ( *, '(a,i8)' ) &
      '  The number of iterations taken  = ', tray
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,g14.6)' ) '  diameter (sigma) = ', sigma
    write ( *, '(a,g14.6)' ) '  density (rho) = ', rho
   
    write ( *, '(a,g14.6)' ) '  alpha = ', alpha
  
     
    write ( *, '(a)' ) ' '

 !!! Convertir a caracteres alfa y epsilon !!!!!!!!
    WRITE(alfa,'(F10.3)')   alpha
    WRITE( colpp,'(F10.3)')   real(nstep)/n
    

   

  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! llamamos al generador de numeros aleatorios
    iseed=2312
    call dran_ini(iseed)
    !pongo el numero de colisiones a cero
    colisiones(:)=0.d00
    sumv = 0.d0

    !inicializo el tiempo para calcular el tiempo de computo
    call cpu_time(start)
    !Abro los archivos en los que voy a guardar los valores de las temperaturas, velocidades, posiciones, etc...
    
    !?? If granular eqv true, the particles lost energy on each collision (inelastic disks) 
    !??else, the collision is elastic between particles 
    granular=.TRUE.
    ! granular=.FALSE.
    
 
    !?? If boolean eqv true, the program saves positions for diferent times steps     
    ! pos_para_t = .TRUE.
    pos_para_t = .FALSE.
    

    !!!! para guardar los valores de las posiciones y velocidades iniciales!!!!!!

    call save_initial_distribution()
   
    density=0
    ! tiempo_relajacion=0

    DO i=1,tray

      

        !!inicializo los tiempos 
        t=0.0
        
        !!inicializo las posiciones y velocidades 
        
        CALL inicializacion1d(n,longy,temp,r,v)
        !! Guardo posiciones y velocidades iniciales 
        call save_initial_distribution()
        !! Miro si estan susperpuestas las posiciones 
        ! CALL superpuesto()
        !! Calculo la temp inicial
        print*, "temperatura inicial" 
        print*, sum(v(:)**2)/n     
             
        ! DO k = 1, n
        !         print*,"velocidad de la particula ",k,"", v(k) 
        !          print*,"pos de la particula ",k,"", r(k) 
        ! END DO    
        
      
        DO j=1,nstep
            ni=0
            colt = HUGE(1.0)
            
            DO k=1,n

                CALL calcular_tiempo_de_colision(k)
                ! CALL     calcular_tiempo_de_colision_point_particles(k) 
                 
            END DO

            ! print*, tcol

            !para saber si los tiempos son negativos de nuevo
            IF (colt<0) THEN 
                WRITE(*,*) 'TIEMPOS NEGATIVOS'
                PRINT*,  'iteracion', j
                STOP 
            END IF 
            !!!Hacer avanzar las posiciones hasta el momento de colisión
            
            t=t+colt
            tiempos(j)=t

            print*, colt

 
                 !obtenemos las densidades para todo t en el eje horizontal. r(:,1) :eje y. r(:,2) :eje z.  
            
            density(i,j,1)=sum(cos(num_onda*r(:)) )
            density(i,j,2)=sum(sin(num_onda*r(:)) )
    
         ! density(i,j,2)=sum(sin(num_onda*r(:,1)) )
           
           
            !! avanzo las posiciones hasta el momento de colisión
                  
                     
                            
            CALL evolve_positions()


            
            

            !colision entre particulas  
           
             
            CALL collide_point_particles(ni(1),ni(2))
            colisiones(i)=colisiones(i)+1
            print*,"Colisiones p.p"
            print*,colisiones/n

            ! print*,'particulas que colisionan', ni(1),ni(2)
            ! print*,'velociades de las particulas que colisionan', v(ni(1)),v(ni(2))

            ! print*,'posiciones de las particulas que colisionan', r(ni(1)),r(ni(2))
            ! CALL superpuesto() 

       


            ! print*, "temperatura tras", j , "colisiones" 
            !     print*, sum(v(:)**2)/n
            !! OBTENEMOS LOS VALORES DE LAS TEMPERATURA 
               DO l=1, n 
                
                sumv(i,j)=sumv(i,j)+  v(l)**2/n   


              
                        
               END DO

               Print*,' temp ' , sumv(i,j)

      

               

                      
                

            

            IF ((j>=nstep-10).AND. (pos_para_t .EQV. .TRUE.) ) THEN 
              CALL save_med_distribution(alfa,j)
            END IF 
             
               
         
                

            
  
        END DO

       
       
       !! Vemos si se superponen las particulas
        ! CALL superpuesto()
        !!Calculamos la temperatura final
        print*, "temperatura final" 
        print*, sum(v(:)**2)/n     
       

        ! Guardamos los valores de las velocidades para todas las trayectorias de forma consecutiva. Asi aumentamos la estadística

        ! OPEN(11,FILE='velocidad_' // trim(adjustl(alfa)) // '.txt',  FORM ='FORMATTED',STATUS='UNKNOWN',POSITION='APPEND'&
        ! ,ACTION='READWRITE')   
        !         DO l=1,n
        !             WRITE(11,*) v(l,1), v(l,2)
        !         END DO
        !     CLOSE(11)
       

        !PRINT*, "numero de colisiones en la iteración ",i ,":", colisiones(i)
       
    END DO
    
    tmp=0.0

    !!Calculamos la temperatura promedio 
    DO l=1,nstep
        DO m=1,tray
        tmp(l)=tmp(l) + sumv(l,m)/tray

        ! tmp(l,m)=2*tmp(l,m)/(temp+tempz)
        END DO
    END DO 








!! Guardo las temperaturas 
    OPEN(9,FILE='temperaturas_' // trim(adjustl(alfa)) // 'colpp_' // trim(adjustl(colpp)) // '.txt',STATUS='unknown')
    DO l=1,nstep
        WRITE(9,*) tmp(l)
    END DO 
    CLOSE(9)   




!! Guardo las posiciones finales 

    OPEN(10,FILE='pos_' // trim(adjustl(alfa)) // 'colpp_' // trim(adjustl(colpp)) // '.txt',STATUS='unknown')   
        DO l=1,n
         WRITE(10,*) r(l)
        END DO
    CLOSE(10) 
!! Guardo las velocidades finales 
        OPEN(10,FILE='vel_' // trim(adjustl(alfa)) //  'colpp_' // trim(adjustl(colpp)) // '.txt',STATUS='unknown')   
        DO l=1,n
         WRITE(10,*) v(l)
        END DO
    CLOSE(10) 



!! Guardo los tiempos de colision 

    OPEN(12,FILE='tiemposdecol_' // trim(adjustl(alfa)) //  'colpp_' // trim(adjustl(colpp)) // '.txt',STATUS='unknown') 
    DO l=1,nstep
        WRITE(12,*) tiempos(l)
    END DO
    CLOSE(12) 
    





 
   
    call save_data_file()
 
    !final del programa
    call cpu_time(finish)
    !calcula el tiempo de computo
    WRITE(*,*) '(Tiempo = ', finish-start , 'segundos.)'

    
    
     WRITE(*,*) '(Colisiones p.p = ', real(colisiones(:))/real(n) , ')'


    DEALLOCATE(r,v,sumv,tmp,colisiones,tiempos)



    CONTAINS
        !!!! para guardar los valores de las posiciones y velocidades iniciales!!!!!!
        subroutine save_initial_distribution()

            IMPLICIT NONE

            INTEGER:: kk


            OPEN(7,FILE='velocidad_init.txt',STATUS='unknown')                       
            OPEN(8,FILE='posiciones_init.txt',STATUS='unknown')                      
            ! CALL inicializacion1d(n,longy,temp,r,v)                      
            DO kk=1,n                                                                
                WRITE(7,*) v(kk)
                WRITE(8,*)  r(kk)
            END DO
            CLOSE(7)
            CLOSE(8)

        end subroutine save_initial_distribution

        subroutine save_med_distribution(aa,bb)
            
            character(len=10)::cc,aa
            INTEGER :: bb,ii
            ! REAL(kind=8) :: aa 
            WRITE(cc,'(I10)') bb 
            ! WRITE(dd,'(F8.3)') aa 

            OPEN(98,FILE='pos_' // trim(adjustl(aa)) // '_tiempo_'// trim(adjustl(cc)) // '.txt',STATUS='unknown')                       
            OPEN(99,FILE='vel_' // trim(adjustl(aa)) // '_tiempo_' // trim(adjustl(cc)) // '.txt',STATUS='unknown')                      
                                 
            DO ii=1,n                                                                
                WRITE(98,*) r(ii)
                WRITE(99,*)  v(ii)
            END DO
            CLOSE(98)
            CLOSE(99)

        end subroutine save_med_distribution

        subroutine save_data_file()
            implicit none 
            LOGICAL :: file_exists
            INQUIRE(FILE="data.txt", EXIST=file_exists)   ! file_exists will be TRUE if the file
            

                                               ! exists and FALSE otherwise
            OPEN(UNIT=35,FILE='data.txt', FORM ='FORMATTED',STATUS='UNKNOWN',POSITION='APPEND'&
            ,ACTION='READWRITE')
            IF (file_exists.EQV. .FALSE.) THEN
                    write(35,* ) ' N ',' T ', ' alpha    ' , ' rho ',  ' colisiones p.p. ',  ' tiempo '
                     write(35,* )  n , temp, alpha ,  rho ,  real(colisiones(:))/real(n) , nstep
                    write(35,* ) ' '
            ELSE
                write(35,* )  n , temp, alpha ,  rho ,  real(colisiones(:))/real(n) , nstep
                write(35,* ) ' '
            END IF 
        end subroutine save_data_file
   

        SUBROUTINE collide ( a, b)
            IMPLICIT NONE
            INTEGER, INTENT(in)  :: a, b   ! Colliding atom indices
            ! LOGICAL :: granular
            ! This routine implements collision dynamics, updating the velocities
            ! The colliding pair (i,j) is assumed to be in contact already

            REAL(kind=8) :: rij, vij
            REAL(kind=8)               :: factor, modulus

            rij = r(a) - r(b)
            ! print*, r(a), r(b)
            rij= rij - longy*ANINT( rij/(longy) ) ! Separation vector
            vij = v(a) - v(b)           ! Relative velocity

            factor =  rij *  vij !dot product of relative velocities and relative positions
            ! modulus = DOT_PRODUCT (rij,rij) !modulus of relative positions
            modulus = sigma !! El módulo es el diámetro de las partículas. 
            
            if(granular .EQV. .TRUE.) THEN 
                vij    = -((1.0d0+alpha)*factor * rij)/(2.0d0*modulus)
            ELSE
                vij    = -factor * rij / modulus
            END IF

            v(a) = v(a) + vij
            v(b) = v(b) - vij


                ! PRINT*, "velocidad particula ",a,"",v(a)
                ! PRINT*, "velocidad particula ",b,"",v(b)
        END SUBROUTINE collide

        subroutine calcular_tiempo_de_colision(c) 
            IMPLICIT NONE

            INTEGER, INTENT(in)  :: c  
            INTEGER :: ii
            DO  ii = 1, n

                IF(ii/=c) THEN 
                    
                            rab=r(c)-r(ii) ! calculamos posiciones relativas
                            rab=rab-longy*ANINT((rab-longy/2.0)/(longy)) ! condiciones periodicas
                            vab=v(c)-v(ii)   !calculamos velocidades relativas
                            bij    = rab * vab    ! obtenemos el producto escalar (ri-rj)*(vi-vj)
                            
                            !! FIRST WAY TO COMPUTE 

                                ! IF (bij<0.0d00 ) THEN !Same as sign(1.0,bij)<0
                                IF (sign(1.0d00,bij)<0) THEN 
                                    ! PRINT*, sign(1.0d00,bij)
                                        ! PRINT*, bij
                                    discr=bij**2-(rab**2-sigma**2)*vab**2
                                    IF( discr>0.0) THEN ! si colisiona con la sucesiva particula
                                        ! tcol = ( -bij - SQRT ( discr ) ) / ( SUM ( vab**2 ) )
                                        !! ALTERNATIVE WAY 
                                        
                                        ! qij=-(bij+sign(1.0d00,bij)*dsqrt(discr))    
                                        !  
                                        ! tcol=MIN(qij/abs((vab**2)),((rab**2)-sigma**2)/qij )
                                        
                                        tcol = ( - bij  - dSQRT ( discr ) ) /  ( vab**2 ) 
                                        !comprobar que los tiempos no son negativos
                                        IF (tcol<0) THEN 
                                            PRINT*, 'colisión:',c,ii,'tiempo',tcol
                                        END IF 
                                       
                                        IF (tcol<colt ) THEN
                                       
                                            ! PRINT*, 'colisión:',c,ii,'tiempo',tcol
                                         
                                            colt=tcol
                                            ni(1)=c
                                            ni(2)=ii
                                        END IF
                                    END IF
                                END IF
                     
                END IF
            END DO 
            


        end subroutine calcular_tiempo_de_colision

        
        SUBROUTINE collide_point_particles ( a, b)
            IMPLICIT NONE
            INTEGER, INTENT(in)  :: a, b   ! Colliding atom indices
            real(kind= 8):: v1,v2

                
                ! PRINT*, "velocidad particula  ",a,"antes",v(a)
                ! PRINT*, "velocidad particula ",b,"antes",v(b)
            v1 = v(a)
            v2 = v(b)
            v(a) = v1*(1-alpha)*0.5d0 + v2*(1+alpha)*0.5d0
            v(b) = v2*(1-alpha)*0.5d0 + v1*(1+alpha)*0.5d0


                ! PRINT*, "velocidad particula  ",a,"despues",v(a)
                ! PRINT*, "velocidad particula ",b,"despues",v(b)
        END SUBROUTINE collide_point_particles



        subroutine calcular_tiempo_de_colision_point_particles(c) 
            IMPLICIT NONE

            INTEGER, INTENT(in)  :: c  
            INTEGER :: ii
            DO  ii = c, n
                IF (c == 1) THEN 
                          rab=r(c)-r(n) ! calculamos posiciones relativas
                            rab=rab- longy*ANINT((rab-longy/2.0)/(longy))  ! condiciones periodicas
                            vab=v(c)-v(n)   !calculamos velocidades relativas
                            bij    = rab * vab    ! obtenemos el producto escalar (ri-rj)*(vi-vj)
                            
                            !! FIRST WAY TO COMPUTE 

                                !? IF (bij<0.0d00 ) THEN !Same as sign(1.0,bij)<0
                                IF (sign(1.0d00,vab)<0 ) THEN 
                                   
                                        
                                        tcol = (  abs(rab)  ) /  ( abs(vab) ) 
                                        !!comprobar que los tiempos no son negativos
                                        IF (tcol<0) THEN 
                                            ! PRINT*, 'colisión:',c,ii,'tiempo',tcol
                                        END IF 
                                        
                                        IF (tcol<colt ) THEN
                                       
                                            ! PRINT*, 'colisión:',c,ii,'tiempo',tcol
                                         
                                            colt=tcol
                                            ni(1)=c
                                            ni(2)=n
                                        END IF
                                    
                                END IF

                END IF 


                IF(ii/=c) THEN 
                    
                            rab=r(c)-r(ii) ! calculamos posiciones relativas
                            rab=rab-longy*ANINT(rab/(longy)) ! condiciones periodicas
                            vab=v(c)-v(ii)   !calculamos velocidades relativas
                            bij    = rab * vab    ! obtenemos el producto escalar (ri-rj)*(vi-vj)
                            
                            !! FIRST WAY TO COMPUTE 

                                !? IF (bij<0.0d00 ) THEN !Same as sign(1.0,bij)<0
                                IF (sign(1.0d00,vab)<0 ) THEN 
                                   
                                        
                                        tcol = (  abs(rab)  ) /  ( abs(vab) ) 
                                        !!comprobar que los tiempos no son negativos
                                        IF (tcol<0) THEN 
                                            ! PRINT*, 'colisión:',c,ii,'tiempo',tcol
                                        END IF 
                                        
                                        IF (tcol<colt ) THEN
                                       
                                            ! PRINT*, 'colisión:',c,ii,'tiempo',tcol
                                         
                                            colt=tcol
                                            ni(1)=c
                                            ni(2)=ii
                                        END IF
                                    
                                END IF
                     
                END IF
            END DO 
            


        end subroutine calcular_tiempo_de_colision_point_particles

        subroutine  evolve_positions()
            IMPLICIT NONE
          

                r(:)     = r(:) + colt * v(:) 
                r(:)     = r(:) - longy*ANINT((r(:)-longy/2.0)/(longy)) ! Apply periodic boundaries

           
          
        end subroutine evolve_positions
      




        SUBROUTINE superpuesto() 
        LOGICAL::super
        INTEGER::q

        super=.FALSE. !! 1 si es falso, 0 si es verdadero
        iloop:DO q=1,n
            IF(ABS(r(q+1)-r(q)) <( 0.999) ) THEN
                !PRINT*, 'particula ',q,'con',q+1,'superpuestas',ABS(r(q+1,1)-r(q,1))
            super=.TRUE.
            !PRINT*, 'particula',q,'distansia', ABS(r(q+1,1)-r(q,1))
                EXIT iloop !salir del loop 
            END IF
        END DO iloop

        IF (super.EQV. .TRUE.) THEN 
        PRINT*, "LAS PARTICULAS ESTAN SUPERPUESTAS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! :,("
        PRINT * , "particulas" , q, q+1
        PRINT * , "distancia" , abs(r(q+1)-r(q))
        ELSE
            PRINT*, "LAS PARTICULAS ESTAN BIEN :)"
        END IF

        END SUBROUTINE superpuesto




END PROGRAM final_version