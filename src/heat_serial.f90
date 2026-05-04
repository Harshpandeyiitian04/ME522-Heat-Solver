!=============================================================================
! heat_serial.f90
!
! Serial 2D Steady-State Heat Conduction Solver
! Governing equation : d2T/dx2 + d2T/dy2 = 0   (2-D Laplace equation)
! Domain             : [0,1] x [0,1]
! Boundary conditions: T = 100 on top edge, T = 0 on other three edges
! Discretisation     : Uniform N x N interior grid, Finite Difference Method
! Iterative solver   : Jacobi method with max-norm convergence test
!
! Usage  : ./heat_serial  N  [cache]
!            N     = number of interior grid points per side
!            cache = 1  => also run loop-ordering comparison (optional)
!
! Outputs:
!   data/serial_timing.csv  – one CSV line:  N, threads, wall_time, iterations
!   data/temp_NxN.dat       – full temperature field for plotting
!=============================================================================
program heat_serial
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  !-- Problem parameters
  integer  :: N, max_iter, iter, ios, do_cache_flag, i, j
  real(dp) :: tol, err, t_start, t_end, elapsed
  real(dp), allocatable :: T(:,:), Tnew(:,:)
  character(len=20) :: arg
  character(len=60) :: fname

  !-- Defaults
  N              = 256
  max_iter       = 200000
  tol            = 1.0d-6
  do_cache_flag  = 0

  !-- Command-line arguments
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg);  read(arg,*) N
  end if
  if (command_argument_count() >= 2) then
    call get_command_argument(2, arg);  read(arg,*) do_cache_flag
  end if

  write(*,'(A)')           repeat('=',56)
  write(*,'(A,I4,A,I4)')   '  Serial 2D Heat Solver   grid = ', N,' x ',N
  write(*,'(A)')           repeat('=',56)

  !-- Allocate: indices 0..N+1 include ghost boundary points
  allocate(T(0:N+1, 0:N+1), Tnew(0:N+1, 0:N+1), stat=ios)
  if (ios /= 0) stop '*** Allocation failed ***'

  !-- Initialise to zero
  T    = 0.0d0
  Tnew = 0.0d0

  !-- Dirichlet BCs:  top edge T=100,  bottom/left/right T=0
  T(0:N+1, N+1) = 100.0d0
  Tnew(0:N+1, N+1) = 100.0d0

  !-- Jacobi iteration
  call cpu_time(t_start)

  iter = 0
  err  = huge(1.0d0)

  do while (err > tol .and. iter < max_iter)
    iter = iter + 1
    err  = 0.0d0

    ! Cache-FRIENDLY loop order: i is the INNER loop
    ! Fortran stores arrays in column-major order, so T(i,j) and T(i+1,j)
    ! are adjacent in memory.  Iterating i fastest maximises cache reuse.
    do j = 1, N
      do i = 1, N
        Tnew(i,j) = 0.25d0 * ( T(i-1,j) + T(i+1,j) &
                              + T(i,j-1) + T(i,j+1) )
        err = max(err, abs(Tnew(i,j) - T(i,j)))
      end do
    end do

    T(1:N, 1:N) = Tnew(1:N, 1:N)     ! swap interior only (BCs preserved)
  end do

  call cpu_time(t_end)
  elapsed = t_end - t_start

  write(*,'(A,I8)')        '  Iterations : ', iter
  write(*,'(A,ES12.4)')    '  Residual   : ', err
  write(*,'(A,F12.6,A)')   '  CPU time   : ', elapsed, ' s'

  !-- Append to timing CSV
  call ensure_data_dir()
  open(unit=20, file='data/serial_timing.csv', position='append', &
       status='unknown', iostat=ios)
  write(20,'(I0,",",I0,",",F18.8,",",I0)') N, 1, elapsed, iter
  close(20)

  !-- Write temperature field
  write(fname,'(A,I0,A,I0,A)') 'data/temp_', N, 'x', N, '.dat'
  call write_field(T, N, trim(fname))

  !-- Optional cache-ordering comparison
  if (do_cache_flag == 1) call run_cache_test(N, min(500, max_iter), tol)

  deallocate(T, Tnew)

contains

  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  subroutine ensure_data_dir()
    integer :: stat
    call execute_command_line('mkdir -p data', exitstat=stat)
  end subroutine

  !--------------------------------------------------------------------------
  subroutine write_field(T, N, fname)
    integer,          intent(in) :: N
    real(dp),         intent(in) :: T(0:N+1, 0:N+1)
    character(len=*), intent(in) :: fname
    integer :: i, j, fu

    open(newunit=fu, file=fname, status='replace', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') '  WARNING: could not open ', fname
      return
    end if

    do j = 0, N+1
      do i = 0, N+1
        if (i < N+1) then
          write(fu,'(F10.5,1X)', advance='no') T(i,j)
        else
          write(fu,'(F10.5)') T(i,j)
        end if
      end do
    end do
    close(fu)
    write(*,'(A,A)') '  Field written to ', fname
  end subroutine write_field

  !--------------------------------------------------------------------------
  ! Demonstrates the effect of loop ordering on cache performance.
  ! In Fortran (column-major) iterating i as the inner loop is cache-friendly;
  ! swapping the loops causes strided access and reduces throughput.
  !--------------------------------------------------------------------------
  subroutine run_cache_test(N, nstep, tol_in)
    integer,  intent(in) :: N, nstep
    real(dp), intent(in) :: tol_in
    real(dp), allocatable :: A(:,:), B(:,:)
    real(dp) :: e, c1, c2, c3, c4, slowdown
    integer  :: ii, jj, kk, fu

    allocate(A(0:N+1,0:N+1), B(0:N+1,0:N+1))

    write(*,'(A)') ''
    write(*,'(A)') repeat('-',56)
    write(*,'(A,I4,A,I4,A,I0,A)') &
      '  Cache-Ordering Test  grid=',N,'x',N,'  steps=',nstep,' (fixed)'
    write(*,'(A)') repeat('-',56)

    ! --- Cache-FRIENDLY: i-inner (column-major order) ---
    A = 0.0d0; B = 0.0d0
    A(0:N+1, N+1) = 100.0d0;  B(0:N+1, N+1) = 100.0d0
    call cpu_time(c1)
    do kk = 1, nstep
      e = 0.0d0
      do jj = 1, N
        do ii = 1, N
          B(ii,jj) = 0.25d0*(A(ii-1,jj)+A(ii+1,jj)+A(ii,jj-1)+A(ii,jj+1))
          e = max(e, abs(B(ii,jj)-A(ii,jj)))
        end do
      end do
      A(1:N,1:N) = B(1:N,1:N)
      if (e < tol_in) exit
    end do
    call cpu_time(c2)

    ! --- Cache-UNFRIENDLY: j-inner (strided access in memory) ---
    A = 0.0d0; B = 0.0d0
    A(0:N+1, N+1) = 100.0d0;  B(0:N+1, N+1) = 100.0d0
    call cpu_time(c3)
    do kk = 1, nstep
      e = 0.0d0
      do ii = 1, N          !  i is OUTER — causes stride-N access on A(ii,jj)
        do jj = 1, N
          B(ii,jj) = 0.25d0*(A(ii-1,jj)+A(ii+1,jj)+A(ii,jj-1)+A(ii,jj+1))
          e = max(e, abs(B(ii,jj)-A(ii,jj)))
        end do
      end do
      A(1:N,1:N) = B(1:N,1:N)
      if (e < tol_in) exit
    end do
    call cpu_time(c4)

    if (c2-c1 > 0.0d0) then
      slowdown = (c4-c3)/(c2-c1)
    else
      slowdown = 0.0d0
    end if
    write(*,'(A,F10.6,A)') '  Friendly  (i-inner): ', c2-c1, ' s'
    write(*,'(A,F10.6,A)') '  Unfriendly(j-inner): ', c4-c3, ' s'
    write(*,'(A,F6.2,A)')  '  Slowdown factor    : ', slowdown, 'x'

    ! Append to cache CSV
    open(newunit=fu, file='data/cache_timing.csv', position='append', &
         status='unknown')
    write(fu,'(I0,",",F14.8,",",F14.8)') N, c2-c1, c4-c3
    close(fu)

    deallocate(A, B)
  end subroutine run_cache_test

end program heat_serial
