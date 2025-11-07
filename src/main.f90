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

  !===============================
  ! Parâmetros da simulação
  !===============================
  n   = 50
  m   = 50
  Nl  = n*m
  Nf  = 10000
  xmin= 0.0
  ymin= 0.0
  xmax= 10.0
  ymax= 10.0
  dt  = 0.01
  mass= 0.1
  friction = 0.1
  p   = 0.6     ! prob. de evaporação do rastro
  q   = 0.8     ! prob. de permanecer quando NÃO há rastro (se houver rastro, prob. é maior que q)
  time= 10.0
  npast = int(time/dt)
  out_stride = 5        ! grava a cada 10 passos (ajuste livre)
  magpart = .false.

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
  allocate(FORCA1(Nf))
  allocate(FORCA2(Nf))
  allocate(FICA(Nf))

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

    ! Forças aleatórias (ruído) por partícula
    call randomica(-dx, dx, FORCA1, Nf, 1000 + 2*k)
    call randomica(-dy, dy, FORCA2, Nf, 1000 + 2*k + 1)

    ! Integra velocidades e posições
    do i=1,Nf
      call resvel(V1(i), dt, FORCA1(i), friction, mass)  ! x
      call resvel(V2(i), dt, FORCA2(i), friction, mass)  ! y
      call respos(X1(i), dt, V1(i))                      ! x
      call respos(X2(i), dt, V2(i))                      ! y

      ! Condição de contorno reflexiva simples
      if (X1(i) <= xmin) X1(i) = xmin + 0.1*dx
      if (X1(i) >= xmax) X1(i) = xmax - 0.1*dx
      if (X2(i) <= ymin) X2(i) = ymin + 0.1*dy
      if (X2(i) >= ymax) X2(i) = ymax - 0.1*dy
    end do

    ! Reconstroi ocupação de células (zera e marca)
    S = 0.0
    SAUX = 0
    call rebuild_cell_occupancy_unique(n,m,xmin,ymin,dx,dy, &
                                       Nf,X1,X2,FORCA1,FORCA2,S,SAUX)

    ! Decisão de ficar/sair quando NÃO há rastro (usa prob. q)
    call randomica(0.0, 1.0, FICA, Nf, 222222)
    call apply_stay_or_kick(n,m,dx,dy,q,S,SAUX,F,FICA,FORCA1,FORCA2,X1,X2)

    ! Deixa/evapora rastro
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
    write(*,*) '* ActiveTrail-2D  |  v0.1 (higienização básica)               *'
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
    ! Sem varredura da grade: mapeia cada partícula e aplica "empurrão"
    ! leve se a célula já estiver ocupada — mantém 1 partícula por célula.
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
          ! célula ocupada → empurra levemente e tenta de novo (máx 5 vezes)
          X1(i) = X1(i) + 0.05*FORCA1(i)
          X2(i) = X2(i) + 0.05*FORCA2(i)
          if (X1(i) <= xmin) X1(i) = xmin + 0.1*dx
          if (X1(i) >= xmax) X1(i) = xmax - 0.1*dx
          if (X2(i) <= ymin) X2(i) = ymin + 0.1*dy
          if (X2(i) >= ymax) X2(i) = ymax - 0.1*dy
          tries = tries + 1
          if (tries >= 5) then
            ! desiste do empurrão, registra mesmo assim para não travar
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
    ! Se célula ocupada e NÃO há rastro, partícula fica com prob. q; senão “chute” leve
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
  ! Integradores (RK4, mantidos)
  !===============================
  subroutine resvel(a,b,c,d,e)
    implicit none
    real, intent(inout) :: a          ! velocidade
    real, intent(in)    :: b          ! dt
    real, intent(in)    :: c          ! força (ruído)
    real, intent(in)    :: d          ! atrito
    real, intent(in)    :: e          ! massa
    real :: k1,k2,k3,k4
    k1 = b*(c-d*a)/e
    k2 = b*(c-d*(a+(0.5*k1)))/e
    k3 = b*(c-d*(a+(0.5*k2)))/e
    k4 = b*(c-d*(a+k3))/e
    a = a + (1.0/6.0)*(k1+2.0*k2+2.0*k3+k4)
  end subroutine resvel

  subroutine respos(a,b,c)
    implicit none
    real, intent(inout) :: a      ! posição
    real, intent(in)    :: b      ! dt
    real, intent(in)    :: c      ! velocidade
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

    ! Não reseeda: apenas consome alguns números extra se quiser variar sequência
    do i=1, max(0, mod(d,17))
      call random_number(c(1:1))
    end do

    call random_number(c)
    c = a + (b-a)*c
  end subroutine randomica

end program active_trail2d
