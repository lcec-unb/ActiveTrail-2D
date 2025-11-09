program active_trail2d
  implicit none

  !===============================
  ! Variáveis
  !===============================
  integer :: i,j,k,l
  integer :: n,m, Nf, Nl
  real    :: xmin,xmax,ymin,ymax
  real    :: dx, dy, dt, mass, friction
  real, allocatable :: S(:,:), F(:,:)
  integer, allocatable :: SAUX(:)        ! índice linear -> índice da partícula
  real, allocatable :: V1(:),V2(:),X1(:),X2(:)
  real, allocatable :: RAND1(:), RAND2(:)
  real, allocatable :: FORCA1(:), FORCA2(:)
  real, allocatable :: FICA(:)
  real              :: EVAP(1)
  integer           :: npast, out_stride
  real              :: p,q,time
  logical           :: magpart

  ! ===== Novos para EM + Repulsão =====
  integer, allocatable :: head(:,:), next(:), cell_of(:)
  real,    allocatable :: FX(:), FY(:)     ! força determinística (repulsão) por partícula
  real :: d0, k_rep                         ! raio efetivo e constante de repulsão
  real :: sigma                             ! intensidade do ruído (EM): aceleração estocástica

  character(42) texto
  !===============================
  ! Parâmetros da simulação 
  ! Lidos do config.dat
  !===============================

505   FORMAT(1X,A40,1X,E11.4E2)
507   FORMAT(1X,A40,1X,I6)
508   FORMAT(1X,A40,L10)

open(3,file='config.dat')
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,508) texto,magpart
      READ(3,507) texto,Nf
      READ(3,507) texto,n
      READ(3,507) texto,m
      READ(3,505) texto,xmin
      READ(3,505) texto,xmax
      READ(3,505) texto,ymin
      READ(3,505) texto,ymax      
      READ(3,505) texto,dt
      READ(3,505) texto,time
      READ(3,507) texto,out_stride
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ(3,505) texto,mass
      READ(3,505) texto,friction
      READ(3,505) texto,d0
      READ(3,505) texto,k_rep
      READ(3,505) texto,sigma
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto  
      READ(3,505) texto,p
      READ(3,505) texto,q     
  
  npast = int(time/dt)
  dx  = (xmax - xmin)/real(n-1)
  dy  = (ymax - ymin)/real(m-1)

  !===============================
  ! Alocações
  !===============================
  allocate(S(n,m));          S = 0.0
  allocate(F(n,m));          F = 0.0
  allocate(SAUX(Nl));        SAUX = 0

  allocate(V1(Nf));          V1 = 0.0
  allocate(V2(Nf));          V2 = 0.0
  allocate(X1(Nf));          X1 = 0.0
  allocate(X2(Nf));          X2 = 0.0

  allocate(RAND1(Nf))
  allocate(RAND2(Nf))
  allocate(FORCA1(Nf));      FORCA1 = 0.0
  allocate(FORCA2(Nf));      FORCA2 = 0.0
  allocate(FICA(Nf))

  allocate(head(n,m));       head = 0
  allocate(next(Nf));        next = 0
  allocate(cell_of(Nf));     cell_of = 0
  allocate(FX(Nf));          FX = 0.0
  allocate(FY(Nf));          FY = 0.0

  call banner()

  !===============================
  ! Distribuição inicial
  !===============================
  call randomica(xmin,xmax,RAND1,Nf,1)
  call randomica(ymin,ymax,RAND2,Nf,2)

  do i=1,Nf
    X1(i) = RAND1(i)
    X2(i) = RAND2(i)
  end do

  ! Mapeia partículas -> células (O(Nf))
  call rebuild_cell_occupancy(n,m,xmin,ymin,dx,dy,Nf,X1,X2,S,SAUX)

  write(*,*) 'Initial distribution of particles - ok'

  !===============================
  ! Rastro inicial aleatório (seed leve)
  !===============================
  call seed_traces(n,m, int(0.4*real(n)), F)  ! espalha ~40%*n marcas iniciais

  !===============================
  ! Arquivos de saída (CSV)
  !===============================
  open (2,file='trajectories.csv', status='replace')
  write(2,*) 't, id, x, y, u, v'

  open (3,file='traces.csv', status='replace')
  write(3,*) 't, x, y'

  !===============================
  ! Loop temporal principal
  !===============================
  do k=1,npast
    write(*,*) 'Current time-step', k

    ! --------- RUÍDO GAUSSIANO (EM): aceleração estocástica ---------
    call randn_gauss(Nf, sigma*sqrt(dt), FORCA1)  ! ax ~ N(0, sigma^2) * sqrt(dt)
    call randn_gauss(Nf, sigma*sqrt(dt), FORCA2)  ! ay ~ N(0, sigma^2) * sqrt(dt)

    ! --------- FORÇA DETERMINÍSTICA DE REPULSÃO (curto alcance) -----
    head = 0; next = 0; cell_of = 0
    call build_cell_lists(n,m,xmin,ymin,dx,dy,Nf,X1,X2, head,next, cell_of)
    FX = 0.0; FY = 0.0
    call compute_repulsion(n,m,dx,dy, Nf, X1,X2, head,next, d0, k_rep, FX,FY)

    ! --------- INTEGRADOR EULER–MARUYAMA ----------------------------
    do i=1,Nf
      ! dV = [ (F_rep - friction*V)/mass ]*dt + sigma*sqrt(dt)*N(0,1)
      V1(i) = V1(i) + ( (FX(i) - friction*V1(i))/mass )*dt + FORCA1(i)
      V2(i) = V2(i) + ( (FY(i) - friction*V2(i))/mass )*dt + FORCA2(i)

      ! X_{n+1} = X_n + V_{n+1}*dt
      X1(i) = X1(i) + V1(i)*dt
      X2(i) = X2(i) + V2(i)*dt

      ! Condição de contorno reflexiva simples
      if (X1(i) <= xmin) X1(i) = xmin + 0.1*dx
      if (X1(i) >= xmax) X1(i) = xmax - 0.1*dx
      if (X2(i) <= ymin) X2(i) = ymin + 0.1*dy
      if (X2(i) >= ymax) X2(i) = ymax - 0.1*dy
    end do

    ! --------- Reconstroi ocupação de células (zera e marca) ---------
    S = 0.0
    SAUX = 0
    call rebuild_cell_occupancy_unique(n,m,xmin,ymin,dx,dy, &
                                       Nf,X1,X2,FORCA1,FORCA2,S,SAUX)

    ! --------- Regras ORIGINAIS de decisão/evaporação (inalteradas) --
    call randomica(0.0, 1.0, FICA, Nf, 222222)
    call apply_stay_or_kick(n,m,dx,dy,q,S,SAUX,F,FICA,FORCA1,FORCA2,X1,X2)

    call update_traces(n,m,p,S,F,EVAP, 333333 + k)

    ! Saída (amostrada por stride)
    if (mod(k,out_stride)==0) then
      do i=1,Nf
        write(2,'(F12.6,1x,I8,1x,F12.6,1x,F12.6,1x,F12.6,1x,F12.6)') &
             k*dt, i, X1(i), X2(i), V1(i), V2(i)
      end do
      do l=1,n
        do j=1,m
          if (F(l,j) == 1.0) then
            write(3,'(F12.6,1x,F12.6,1x,F12.6)') k*dt, (l-0.5)*dx + xmin, (j-0.5)*dy + ymin
          end if
        end do
      end do
    end if
  end do

  close(2)
  close(3)

contains

  subroutine banner()
    implicit none
    write(*,*) '==============================================================='
    write(*,*) '* ActiveTrail-2D  |  v0.1 + (EM noise + short-range repulsion) *'
    write(*,*) '==============================================================='
  end subroutine banner

  subroutine rebuild_cell_occupancy(n,m,xmin,ymin,dx,dy,Nf,X1,X2,S,SAUX)
    implicit none
    integer, intent(in) :: n,m,Nf
    real,    intent(in) :: xmin,ymin,dx,dy
    real,    intent(in) :: X1(:), X2(:)
    real,    intent(inout) :: S(:,:)
    integer, intent(inout) :: SAUX(:)
    integer :: i, ix, iy, idx

    do i=1,Nf
      ix = 1 + int( floor( (X1(i)-xmin)/dx ) )
      iy = 1 + int( floor( (X2(i)-ymin)/dy ) )
      if (ix < 1) ix = 1; if (ix > n) ix = n
      if (iy < 1) iy = 1; if (iy > m) iy = m
      S(ix,iy) = S(ix,iy) + 1.0
      idx = (ix-1)*m + iy
      SAUX(idx) = i
    end do
  end subroutine rebuild_cell_occupancy

  subroutine rebuild_cell_occupancy_unique(n,m,xmin,ymin,dx,dy, &
                                           Nf,X1,X2,FORCA1,FORCA2,S,SAUX)
    ! Mantida do v0.1: “empurrão” leve para evitar mais de 1 por célula
    implicit none
    integer, intent(in) :: n,m,Nf
    real,    intent(in) :: xmin,ymin,dx,dy
    real,    intent(inout) :: X1(:), X2(:)
    real,    intent(in) :: FORCA1(:), FORCA2(:)
    real,    intent(inout) :: S(:,:)
    integer, intent(inout) :: SAUX(:)
    integer :: i, ix, iy, idx, tries

    do i=1,Nf
      tries = 0
      do
        ix = 1 + int( floor( (X1(i)-xmin)/dx ) )
        iy = 1 + int( floor( (X2(i)-ymin)/dy ) )
        if (ix < 1) ix = 1; if (ix > n) ix = n
        if (iy < 1) iy = 1; if (iy > m) iy = m
        idx = (ix-1)*m + iy

        if (S(ix,iy) < 1.0) then
          S(ix,iy) = S(ix,iy) + 1.0
          SAUX(idx) = i
          exit
        else
          X1(i) = X1(i) + 0.05*FORCA1(i)
          X2(i) = X2(i) + 0.05*FORCA2(i)
          if (X1(i) <= xmin) X1(i) = xmin + 0.1*dx
          if (X1(i) >= xmax) X1(i) = xmax - 0.1*dx
          if (X2(i) <= ymin) X2(i) = ymin + 0.1*dy
          if (X2(i) >= ymax) X2(i) = ymax - 0.1*dy
          tries = tries + 1
          if (tries >= 5) then
            S(ix,iy) = S(ix,iy) + 1.0
            SAUX(idx) = i
            exit
          end if
        end if
      end do
    end do
  end subroutine rebuild_cell_occupancy_unique

  subroutine apply_stay_or_kick(n,m,dx,dy,q,S,SAUX,F,FICA,FORCA1,FORCA2,X1,X2)
    implicit none
    integer, intent(in) :: n,m
    real,    intent(in) :: dx,dy,q
    real,    intent(in) :: S(:,:), F(:,:)
    integer, intent(in) :: SAUX(:)
    real,    intent(in) :: FICA(:), FORCA1(:), FORCA2(:)
    real,    intent(inout) :: X1(:), X2(:)
    integer :: l,j, idx, pid
    do l=1,n
      do j=1,m
        if (S(l,j) == 1.0) then
          if (F(l,j) == 0.0) then
            idx = (l-1)*m + j
            pid = SAUX(idx)
            if (pid > 0) then
              if (FICA(pid) > q) then
                X1(pid) = X1(pid) + 0.1*FORCA1(pid)
                X2(pid) = X2(pid) + 0.1*FORCA2(pid)
              end if
            end if
          end if
        end if
      end do
    end do
  end subroutine apply_stay_or_kick

  subroutine update_traces(n,m,p,S,F,EVAP, seed_hint)
    implicit none
    integer, intent(in) :: n,m, seed_hint
    real,    intent(in) :: p
    real,    intent(in) :: S(:,:)
    real,    intent(inout) :: F(:,:)
    real,    intent(inout) :: EVAP(1)
    integer :: l,j

    do l=1,n
      do j=1,m
        if (S(l,j) == 1.0) then
          F(l,j) = 1.0
        else
          if (F(l,j) == 1.0) then
            call randomica(0.0,1.0,EVAP,1, seed_hint + 13*(l+j))
            if (EVAP(1) < p) F(l,j) = 0.0
          end if
        end if
      end do
    end do
  end subroutine update_traces

  subroutine seed_traces(n,m, nmarks, F)
    implicit none
    integer, intent(in) :: n,m, nmarks
    real,    intent(inout) :: F(:,:)
    real, allocatable :: rx(:), ry(:)
    integer :: i, ix, iy

    allocate(rx(nmarks), ry(nmarks))
    call randomica(1.0, real(n), rx, nmarks, 7777)
    call randomica(1.0, real(m), ry, nmarks, 8888)

    do i=1,nmarks
      ix = max(1, min(n, int(rx(i))))
      iy = max(1, min(m, int(ry(i))))
      F(ix,iy) = 1.0
    end do

    deallocate(rx,ry)
  end subroutine seed_traces

  !===============================
  ! Integradores (RK4, mantidos — não usados no EM)
  !===============================
  subroutine resvel(a,b,c,d,e)
    implicit none
    real, intent(inout) :: a
    real, intent(in)    :: b
    real, intent(in)    :: c
    real, intent(in)    :: d
    real, intent(in)    :: e
    real :: k1,k2,k3,k4
    k1 = b*(c-d*a)/e
    k2 = b*(c-d*(a+(0.5*k1)))/e
    k3 = b*(c-d*(a+(0.5*k2)))/e
    k4 = b*(c-d*(a+k3))/e
    a = a + (1.0/6.0)*(k1+2.0*k2+2.0*k3+k4)
  end subroutine resvel

  subroutine respos(a,b,c)
    implicit none
    real, intent(inout) :: a
    real, intent(in)    :: b
    real, intent(in)    :: c
    real :: k1,k2,k3,k4
    k1 = b*(c)
    k2 = b*((0.5*k1)+c)
    k3 = b*((0.5*k2)+c)
    k4 = b*((k3)+c)
    a = a + (1.0/6.0)*(k1+2.0*k2+2.0*k3+k4)
  end subroutine respos

  !===============================
  ! RNG: semear uma única vez
  !===============================
  subroutine randomica(a,b,c,n,d)
    implicit none
    real,    intent(in)  :: a,b        ! intervalo
    integer, intent(in)  :: n
    real,    intent(out) :: c(n)
    integer, intent(in)  :: d          ! "hint" para perturbar, sem reseedar sempre

    integer :: i, m, e
    integer :: f(8)
    integer, allocatable :: seed(:)
    logical, save :: seeded = .false.

    call random_seed(size = m)
    if (.not. seeded) then
      allocate (seed(m))
      call date_and_time(values=f)
      call system_clock(count=e)
      do i=1,m
        seed(i) =  47 + f(8)*i*12000 + e*(3+i)
      end do
      call random_seed(put = seed)
      deallocate(seed)
      seeded = .true.
    end if

    do i=1, max(0, mod(d,17))
      call random_number(c(1:1))
    end do

    call random_number(c)
    c = a + (b-a)*c
  end subroutine randomica

  !===============================
  ! GAUSSIANO (Box–Muller) p/ EM
  !===============================
  subroutine randn_gauss(N, scale, G)
    implicit none
    integer, intent(in) :: N
    real,    intent(in) :: scale
    real,    intent(out):: G(N)
    real, allocatable :: U1(:), U2(:)
    integer :: i, M
    M = (N+1)/2
    allocate(U1(M), U2(M))
    call randomica(1.0e-6, 1.0-1.0e-6, U1, M, 901)
    call randomica(1.0e-6, 1.0-1.0e-6, U2, M, 902)
    do i=1,M
      G(2*i-1) = scale * sqrt(-2.0*log(U1(i))) * cos(2.0*3.141592653589793*U2(i))
      if (2*i <= N) then
        G(2*i)   = scale * sqrt(-2.0*log(U1(i))) * sin(2.0*3.141592653589793*U2(i))
      end if
    end do
    deallocate(U1,U2)
  end subroutine randn_gauss

  !===============================
  ! Cell-lists e Repulsão curta
  !===============================
  subroutine build_cell_lists(n,m,xmin,ymin,dx,dy,Nf,X1,X2, head,next, cell_of)
    implicit none
    integer, intent(in) :: n,m,Nf
    real,    intent(in) :: xmin,ymin,dx,dy
    real,    intent(in) :: X1(:), X2(:)
    integer, intent(inout) :: head(n,m), next(:), cell_of(:)
    integer :: i, ix, iy
    do i=1,Nf
      ix = 1 + int( floor( (X1(i)-xmin)/dx ) )
      iy = 1 + int( floor( (X2(i)-ymin)/dy ) )
      if (ix < 1) ix = 1; if (ix > n) ix = n
      if (iy < 1) iy = 1; if (iy > m) iy = m
      next(i)     = head(ix,iy)
      head(ix,iy) = i
      cell_of(i)  = (ix-1)*m + iy
    end do
  end subroutine build_cell_lists

  subroutine compute_repulsion(n,m,dx,dy, Nf, X1,X2, head,next, d0,k_rep, FX,FY)
    implicit none
    integer, intent(in) :: n,m,Nf
    real,    intent(in) :: dx,dy, d0, k_rep
    real,    intent(in) :: X1(:), X2(:)
    integer, intent(in) :: head(n,m), next(:)
    real,    intent(inout) :: FX(:), FY(:)

    integer :: ix,iy, i, j, cix, ciy
    real :: rx, ry, r, overlap, invr

    do ix=1,n
      do iy=1,m
        i = head(ix,iy)
        do while (i /= 0)
          do cix=max(1,ix-1), min(n,ix+1)
            do ciy=max(1,iy-1), min(m,iy+1)
              j = head(cix,ciy)
              do while (j /= 0)
                if (j > i) then
                  rx = X1(i) - X1(j)
                  ry = X2(i) - X2(j)
                  r  = sqrt(rx*rx + ry*ry)
                  if (r > 1.0e-12) then
                    if (r < d0) then
                      overlap = d0 - r
                      invr = 1.0/r
                      FX(i) = FX(i) + k_rep*overlap*(rx*invr)
                      FY(i) = FY(i) + k_rep*overlap*(ry*invr)
                      FX(j) = FX(j) - k_rep*overlap*(rx*invr)
                      FY(j) = FY(j) - k_rep*overlap*(ry*invr)
                    end if
                  end if
                end if
                j = next(j)
              end do
            end do
          end do
          i = next(i)
        end do
      end do
    end do
  end subroutine compute_repulsion

end program active_trail2d
