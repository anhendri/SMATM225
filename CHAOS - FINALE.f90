!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ETUDE D'UN SYSTEME A N CORPS A L'AIDE DE L'INDICATEUR DE CHAOS DE MEGNO
!!-------------------------------------------------------------------------------------------------------------------------
!! Programme realise dans le cadre du cours "SMATM225 - Chaos et Determinisme" de l'UNAMUR
!! Derniere modification : Le 07/06/16, a Namur
!! Auteurs : P. CASTENETTO - A. HENDRICKX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MODULE COMMON_VALUE
!!-------------------------------------------------------------------------------------------------------------------------
!! R    : R_ij = 1/r_ij^alpha
!! a    : Matrice hessienne du systeme
!! ALPH : Alpha
!! NTLD : Ntilde
!! NVAR : Nombre de masse
module common_value
  implicit none
  real(kind=8), allocatable :: R(:,:), a(:,:)
  real(kind=8)              :: ALPH, NTLD
  integer                   :: NVAR
end module


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MODULE INTERFACES
!!-------------------------------------------------------------------------------------------------------------------------
module interfaces
  interface
    function calc_H(u0, NTOT) result(H)
      use common_value
      integer :: NTOT, POSX, POSY
      real(kind=8) :: u0(NTOT), H
    end function
    function ode45( NTOT, TMIN, TMAX, PAST, u0, f) result(VARS)
      integer      :: NTOT, NUMT, POSX
      external     :: f
      real(kind=8) :: PAST, TMAX, f0(NTOT), f1(NTOT), f2(NTOT), f3(NTOT), TMIN, u0(NTOT)
      real(kind=8), allocatable :: VARS(:,:)
    end function
    function eulerS( NTOT, TMIN, TMAX, PAST, u0, f) result(VARS)
      integer      :: NTOT, NUMT, POSX
      external     :: f
      real(kind=8) :: PAST, TMAX, f0(NTOT), f1(NTOT), f2(NTOT), f3(NTOT), TMIN, u0(NTOT)
      real(kind=8), allocatable :: VARS(:,:)
    end function
    function verlet( NTOT, TMIN, TMAX, PAST, u0, f) result(VARS)
      integer      :: NTOT, NUMT, POSX
      external     :: f
      real(kind=8) :: PAST, TMAX, f0(NTOT), f1(NTOT), f2(NTOT), f3(NTOT), TMIN, u0(NTOT)
      real(kind=8), allocatable :: VARS(:,:)
    end function
  end interface
end module


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROGRAMME PRINCIPAL
!!-------------------------------------------------------------------------------------------------------------------------
!! NTOT : Nombre d'equation totale (4*NVAR + 1)
!! TMIN : Temps initial d'integration
!! TMAX : Temps final d'integration
!! PAST : Pas d'integration pour le temps
!! E0   : Energie moyenne par masse
program main
  use common_value
  use interfaces
  implicit none

  external                       :: calc
  integer                        :: INDX, POSX, POSY, cas, NTOT, NUMT, NTWK
  real(kind=8)                   :: TMIN, TMAX, PAST, TINI, H0, H1, PASF, PASP, E0
  real(kind=8), allocatable,save :: u0(:), u1(:), VARS(:,:)

    call INIT_RANDOM_SEED()                                                                                                ! INITIALISATION DU GENERATEUR ALEATOIRE

    open(10,file='INPUT.in')                                                                                               ! OUVERTURE DU FICHIER D'ENTREE
    read(10,*) NVAR, ALPH, cas, TMIN, TMAX, PAST, PASP, PASF, E0, NTWK                                                     ! LECTURE DES PARAMETRES DANS LE FICHIER D'ENTREE
    close(10)                                                                                                              ! FERMETURE DU FICHIER D'ENTREE
    NTOT = 4*NVAR+2                                                                                                        ! CALCUL DU NOMBRE TOTAL D'EQUATION
    NUMT = floor((TMAX-TMIN)/PAST)+1                                                                                       ! CALCUL DU NOMBRE D'INTERVALLES D'INTEGRATION

    allocate(R(NVAR,NVAR), a(NVAR,NVAR),u0(NTOT),u1(NTOT),VARS(NUMT,NTOT))                                                 ! ALLOCATION DE LA MEMOIRE

    if(cas==1) then                                                                                                        ! CALCUL DE Ntilde ET DE R EN FONCTION DU CAS (V. ARTICLE)
        NTLD = sum(1./(/((/(POSY-POSX, POSY=POSX+1,NVAR+1)/),POSX=0,NVAR)/)**ALPH)/NVAR
    elseif(cas==2) then
        NTLD = sum(2./(/(POSX,POSX=1,NVAR/2-1)/)**ALPH)
    endif

    if(NTWK==1) then
        if(cas==1) then                                                                                                    ! CALCUL DE Ntilde ET DE R EN FONCTION DU CAS (V. ARTICLE)
            forall(POSX=1:NVAR,POSY=1:NVAR) R(POSX,POSY) = abs(POSX-POSY)
        elseif(cas==2) then
            forall(POSX=1:NVAR,POSY=1:NVAR) R(POSX,POSY) = min(abs(POSX-POSY),NVAR-abs(POSX-POSY))
        endif
        WHERE(R==0) R=1.                                                                                                   ! CETTE LIGNE SERT A EVITER LES INFINIS DANS LE CALCUL DE 1/R
        R = R**(-ALPH)                                                                                                     ! CALCUL DE 1/R** ALPHA (POUR EVITER DE REFAIRE LA CALCUL A CHAQUE FOIS)
    else
        open(10,file='Matrix.dat')
        read(*,*) R
        close(10)
    endif

    TINI = TMIN                                                                                                            ! SAUVEGARDE DU TEMPS INITIAL

    u0 = 0.                                                                                                                ! CREATION DES CONDITIONS INITIALES
    call random_number(u0(NVAR+1:2*NVAR))
    u0(NVAR+1:2*NVAR-1) = 2.*sqrt(2.*E0)*(u0(NVAR+1:2*NVAR-1)-.5)
    u0(2*NVAR) = sqrt(2.*E0*NVAR - sum(u0(NVAR+1:2*NVAR-1)**2.))                                                           ! AJUSTEMENT DE L'ENERGIE TOTALE INITIALE
    call random_number(u0(3*NVAR+1:4*NVAR))
    u1 = u0                                                                                                                ! SAUVEGARDE DES CONDITIONS INITIALES

    open(100,file = 'INITIALE.dat')                                                                                        ! OUVERTURE DU FICHIER DE SAUVEGARDE DES CONDITIONS INITIALES
    write(100,'(1e15.7 )') u0                                                                                              ! ECRITURE DES CONDITIONS INITIALES DANS LE FICHIER
    close(100)                                                                                                             ! FERMETURE DU FICHIER

    open(101,file = 'ENERGIE.dat')                                                                                         ! OUVERTURE DU FICHIER DE SAUVEGARDE DES ENERGIES
    write(101,'(1e15.7)') calc_H(u0, NTOT)                                                                                 ! CALCUL DE L'ENERGIE INITIALE

    VARS = ode45(NTOT,TMIN,TMAX,PAST,u0,calc)                                                                              ! INTEGRATION PAR ODE45
    open(10,file='OUTPUT1.out')                                                                                            ! OUVERTURE DU FICHIER DE SAUVEGARDE DES Nu/P EN FONCTION DU TEMPS
    do INDX = 1,NUMT                                                                                                       ! POUR CHAQUE VALEUR DU TEMPS ...
        write(10,*) TINI+INDX*PAST, VARS(INDX,1:4*NVAR)                                                                    ! ECRITURE DANS LE FICHIER DES ANGLES ET VITESSES CORRESPONDANTES
    enddo
    write(101,'(1e15.7)') calc_H(VARS(NUMT,:), NTOT)                                                                       ! SAUVEGARDE DE L'ENERGIE APRES INTEGRATION

    TMIN = TINI                                                                                                            ! REINITIALISATION DES CONDITIONS INITIALES
    u0 = u1
    VARS = eulerS(NTOT,TMIN,TMAX,PAST,u0,calc)                                                                             ! INTEGRATION PAR LA METHODE D'EULER
    open(10,file='OUTPUT2.out')                                                                                            ! OUVERTURE DU FICHIER DE SAUVEGARDE DES Nu/P EN FONCTION DU TEMPS
    do INDX = 1,NUMT                                                                                                       ! POUR CHAQUE VALEUR DU TEMPS ...
        write(10,*) TINI+INDX*PAST, VARS(INDX,1:4*NVAR)                                                                    ! ECRITURE DANS LE FICHIER DES ANGLES ET VITESSES CORRESPONDANTES
    enddo
    write(101,'(1e15.7)') calc_H(VARS(NUMT,:), NTOT)                                                                       ! SAUVEGARDE DE L'ENERGIE APRES INTEGRATION

    TMIN = TINI                                                                                                            ! REINITIALISATION DES CONDITIONS INITIALES
    u0 = u1
    VARS = verlet(NTOT,TMIN,TMAX,PAST,u0,calc)                                                                             ! INTEGRATION PAR LA METHODE DE VERLET
    open(10,file='OUTPUT3.out')                                                                                            ! OUVERTURE DU FICHIER DE SAUVEGARDE DES Nu/P EN FONCTION DU TEMPS
    do INDX = 1,NUMT                                                                                                       ! POUR CHAQUE VALEUR DU TEMPS ...
        write(10,*) TINI+INDX*PAST, VARS(INDX,1:4*NVAR)                                                                    ! ECRITURE DANS LE FICHIER DES ANGLES ET VITESSES CORRESPONDANTES
    enddo
    write(101,'(1e15.7)') calc_H(VARS(NUMT,:), NTOT)                                                                       ! SAUVEGARDE DE L'ENERGIE APRES INTEGRATION
    close(101)                                                                                                             ! FERMETURE

    open(100,file = 'MEGNO.dat')                                                                                           ! FICHIER CONTENANT L'EVOLUTION DE MEGNO EN FONCTION DU TEMPS
    write(100,'(4e15.7)') (TINI+INDX*PAST, VARS(INDX,NTOT-1:)/(TINI+INDX*PAST), 2*VARS(INDX,NTOT)/(TINI+INDX*PAST)/TMAX,  &
                           INDX = 1,NUMT)
    close(100)

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FONCTION CALC_H : Fonction permettant de calculer l'hamiltonien a base des coordonnees
!!-------------------------------------------------------------------------------------------------------------------------
!! u0   : Coordonnees pour lesquelles on veut calculer H
!! H    : Valeur de l'energie pour les coordonnees u0
function calc_H(u0,NTOT) result(H)
use common_value
implicit none
integer      :: NTOT, POSX, POSY
real(kind=8) :: u0(NTOT), H

    H = 0.                                                                                                                 ! INITIALISATION DE L'HAMILTONIEN
    do POSX = 1,NVAR
        H = H + 1/2.*u0(NVAR+POSX)**2.                                                                                     ! CALCUL DE L'ENERGIE CINETIQUE
        do POSY = 1,NVAR
            if(POSX.NE.POSY) H = H + 1./(2.*NTLD)*(1.-cos(u0(POSX)-u0(POSY)))*R(POSX,POSY)                                 ! CALCUL DE L'ENERGIE POTENTIELLE
        enddo
    enddo
end function



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FUNCTION CALC : Fonction permettant de calculer les derivees a partir des coordonnees
!!-------------------------------------------------------------------------------------------------------------------------
!! UDOT : Vecteurs contenant les derivees evaluee au temps TMPS aux coordonnees U
!! U    : Coordonnees auxquelles on veut evaluer la derivee
!! TMPS : Temps auquel on veut evaluer la derivee (utile pour MEGNO)
!! NTOT : Nombre d'equations
function calc( TMPS, NTOT, u) result(UDOT)
use common_value
  implicit none
  integer      :: NTOT, POSX, POSY
  real(kind=8) :: TMPS, u(NTOT), UDOT(NTOT)

    a = 0.                                                                                                                 ! CALCUL DE LA MATRICE A (MATRICE HESSIENNE)
    do POSX=1,NVAR
        do POSY = 1,NVAR
            if(POSX.NE.POSY) then
                a(POSX,POSX) = a(POSX,POSX) - cos(u(POSX)-u(POSY))*R(POSX,POSY)                                            ! FORMULE POUR Aii
                a(POSX,POSY) = cos(u(POSX)-u(POSY))*R(POSX,POSY)                                                           ! FORMULE POUR Aij
            endif
        enddo
    enddo
    a = a/NTLD;

    UDOT = 0.                                                                                                              ! INITIALISATION DES DERIVEES
    UDOT(:NVAR) = u(NVAR+1:2*NVAR);                                                                                        ! CALCUL DES Nu POINTS
    do POSX=1,NVAR                                                                                                         ! CALCUL DES P POINTS
        UDOT(POSX+NVAR) = -sum(sin(u(POSX)-u(:POSX-1))*R(POSX,:POSX-1))/NTLD -sum(sin(u(POSX)-u(POSX+1:NVAR)) &
                        * R(POSX,POSX+1:NVAR))/NTLD
    enddo
    UDOT(2*NVAR+1:3*NVAR) = u(3*NVAR+1:4*NVAR);                                                                            ! CALCUL DES Delta Nu POINT
    UDOT(3*NVAR+1:4*NVAR) = matmul(a,UDOT(2*NVAR+1:3*NVAR))                                                                ! CALCUL DES Delta P POINT
    UDOT(4*NVAR+1) = 2*sum(UDOT(2*NVAR+1:4*NVAR)*u(2*NVAR+1:4*NVAR))/sum(u(2*NVAR+1:4*NVAR)**2)*TMPS                       ! CALCUL POUR MEGNO
    if(TMPS>0) UDOT(4*NVAR+2) = U(4*NVAR+1)/TMPS                                                                           ! CALCUL POUR MEGNO MOYEN
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FUNCTION ODE45 : Fonction permettant l'integration par Runge Kutta d'ordre 4-5
!! ------------------------------------------------------------------------------------------------------------------------
!! NTOT : Nombre d'equation totale (4*NVAR + 1)
!! TMIN : Temps initial d'integration
!! TMAX : Temps final d'integration
!! PAST : Pas d'integration pour le temps
!! u0   : Conditions initiales du système
!! f    : Fonction permettant d'évaluer la dérivée
!! VARS : Matrice contenant, pour chaque temps, l'état du système
function ode45(NTOT, TMIN, TMAX, PAST, u0, f) result(VARS)
  implicit none
  integer                   :: NTOT, NUMT, INDX
  real(kind=8)              :: PAST, TMAX, f0(NTOT), f1(NTOT), f2(NTOT), f3(NTOT), TMIN, u0(NTOT)
  real(kind=8), allocatable :: VARS(:,:)

  interface
    function f(TMPS, NTOT, u) result(UDOT)
      integer      :: NTOT
      real(kind=8) :: TMPS, u(NTOT), UDOT(NTOT)
    end function
  end interface

    NUMT = ceiling((TMAX-TMIN)/PAST)                                                                                       ! CALCUL DU NOMBRE D'INTERVALLES DE TEMPS
    allocate(VARS(NUMT,NTOT))                                                                                              ! ALLOCATION DE LA MEMOIRE

    do INDX = 1,NUMT                                                                                                       ! BOUCLE SUR LE TEMPS
        f0           = f(TMIN       , NTOT, u0          )
        f1           = f(TMIN+PAST/2, NTOT, u0+PAST*f0/2)
        f2           = f(TMIN+PAST/2, NTOT, u0+PAST*f1/2)
        f3           = f(TMIN+PAST  , NTOT, u0+PAST*f2  )
        u0           = u0 + PAST*(f0 + 2.*f1 + 2.*f2 + f3)/6.                                                              ! CALCUL DES NOUVELLES COORDONNEES
        TMIN         = TMIN + PAST                                                                                         ! CALCUL DU TEMPS SUIVANT
        VARS(INDX,:) = u0                                                                                                  ! SAUVEGARDER DES COORDONNEES
    end do
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FUNCTION EULER : Fonction permettant l'integration par la méthode d'Euler
!! ------------------------------------------------------------------------------------------------------------------------
!! NTOT : Nombre d'equation totale (4*NVAR + 1)
!! TMIN : Temps initial d'integration
!! TMAX : Temps final d'integration
!! PAST : Pas d'integration pour le temps
!! u0   : Conditions initiales du systeme
!! f    : Fonction permettant d'evaluer la derivee
!! VARS : Matrice contenant, pour chaque temps, l'etat du systeme
function eulerS(NTOT, TMIN, TMAX, PAST, u0, f) result(VARS)
  implicit none
  integer                   :: NTOT, NUMT, INDX, NVAR
  real(kind=8)              :: PAST, TMAX, f0(NTOT), f1(NTOT), f2(NTOT), f3(NTOT), TMIN, u0(NTOT)
  real(kind=8), allocatable :: VARS(:,:)

  interface
    function f(TMPS, NTOT, u) result(UDOT)
      integer      :: NTOT
      real(kind=8) :: TMPS, u(NTOT), UDOT(NTOT)
    end function
  end interface

    NUMT = ceiling((TMAX-TMIN)/PAST)                                                                                       ! CALCUL DU NOMBRE D'INTERVALLES D'INTEGRATION
    NVAR = (NTOT-2)/4                                                                                                      ! CALCUL DU NOMBRE DE MASSE
    allocate(VARS(NUMT,NTOT))                                                                                              ! ALLOCATION DE LA MEMOIRE

    do INDX = 1,NUMT                                                                                                       ! BOUCLE SUR LE TEMPS
        f0                  = f(TMIN, NTOT, u0)
        u0(:NVAR)           = u0(1:NVAR) + PAST*f0(:NVAR)
        u0(2*NVAR+1:3*NVAR) = u0(2*NVAR+1:3*NVAR) + PAST*f0(2*NVAR+1:3*NVAR)
        f0                  = f(TMIN, NTOT, u0)
        u0(NVAR+1:2*NVAR)   = u0(NVAR+1:2*NVAR) + PAST*f0(NVAR+1:2*NVAR)
        u0(3*NVAR+1:)       = u0(3*NVAR+1:) + PAST*f0(3*NVAR+1:)                                                           ! CALCUL DES NOUVELLES COORDONNEES
        TMIN         = TMIN + PAST                                                                                         ! CALCUL DU TEMPS SUIVANT
        VARS(INDX,:) = u0                                                                                                  ! SAUVEGARDER DES COORDONNEES
    end do
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FUNCTION VERLET : Fonction permettant l'integration par la méthode de Verlet
!! ------------------------------------------------------------------------------------------------------------------------
!! NTOT : Nombre d'equation totale (4*NVAR + 1)
!! TMIN : Temps initial d'integration
!! TMAX : Temps final d'integration
!! PAST : Pas d'integration pour le temps
!! u0   : Conditions initiales du systeme
!! f    : Fonction permettant d'evaluer la derivee
!! VARS : Matrice contenant, pour chaque temps, l'etat du systeme
function verlet(NTOT, TMIN, TMAX, PAST, u0, f) result(VARS)
  implicit none
  integer                   :: NTOT, NUMT, INDX, NVAR
  real(kind=8)              :: PAST, TMAX, f0(NTOT), f1(NTOT), f2(NTOT), f3(NTOT), TMIN, u0(NTOT)
  real(kind=8), allocatable :: VARS(:,:)

  interface
    function f(TMPS, NTOT, u) result(UDOT)
      integer      :: NTOT
      real(kind=8) :: TMPS, u(NTOT), UDOT(NTOT)
    end function
  end interface

    NUMT = ceiling((TMAX-TMIN)/PAST)                                                                                       ! CALCUL DU NOMBRE D'INTERVALLES D'INTEGRATION
    NVAR = (NTOT-2)/4                                                                                                      ! CALCUL DU NOMBRE DE MASSE
    allocate(VARS(NUMT,NTOT))                                                                                              ! ALLOCATION DE LA MEMOIRE

    do INDX = 1,NUMT                                                                                                       ! BOUCLE SUR LE TEMPS
        f0                  = f(TMIN, NTOT, u0)
        u0(NVAR+1:2*NVAR)   = u0(NVAR+1:2*NVAR) + PAST/2.*f0(NVAR+1:2*NVAR)
        u0(3*NVAR+1:)       = u0(3*NVAR+1:) + PAST/2.*f0(3*NVAR+1:)
        f0                  = f(TMIN, NTOT, u0)
        u0(:NVAR)           = u0(:NVAR) + PAST*f0(:NVAR)
        u0(2*NVAR+1:3*NVAR) = u0(2*NVAR+1:3*NVAR) + PAST*f0(2*NVAR+1:3*NVAR)
        f0                  = f(TMIN, NTOT, u0)
        u0(NVAR+1:2*NVAR)   = u0(NVAR+1:2*NVAR) + PAST/2.*f0(NVAR+1:2*NVAR)
        u0(3*NVAR+1:)       = u0(3*NVAR+1:) + PAST/2.*f0(3*NVAR+1:)                                                        ! CALCUL DES NOUVELLES COORDONNEES
        TMIN         = TMIN + PAST                                                                                         ! CALCUL DU TEMPS SUIVANT
        VARS(INDX,:) = u0                                                                                                  ! SAUVEGARDER DES COORDONNEES
    end do
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SOUS-ROUTINE INIT_RANDOM_SEED
subroutine Init_Random_Seed()
  integer              :: INDX, NMBR, clock
  integer, allocatable :: seed(:)

    call random_seed(size = NMBR)
    allocate(seed(NMBR))
    call system_clock(COUNT=clock)
    seed = clock + 37*(/(INDX, INDX = 0, NMBR-1)/)
    call random_seed(PUT = seed)
    deallocate(seed)
end subroutine Init_Random_Seed