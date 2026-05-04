!=============================================================================
! heat_parallel.f90
!
! OpenMP-Parallel 2D Steady-State Heat Conduction Solver
! Same problem as heat_serial.f90 but with Jacobi loop parallelised.
!
! Key OpenMP directives used:
!   !$OMP PARALLEL DO  – parallelises the outer j-loop across threads
!   PRIVATE(i)         – each thread owns its own copy of loop variable i
!   REDUCTION(max:err) – combines per-thread max-residuals into a global max
!   SCHEDULE(static)   – static block partitioning for load balance
!
! Usage: ./heat_parallel  N  nthreads
!=============================================================================
program heat_parallel
  use omp_lib
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  integer  :: N, nthreads, max_iter, iter, ios
  real(dp) :: tol, err, t_start, t_end, elapsed
  real(dp), allocatable :: T(:,:), Tnew(:,:)
  character(len=20) :: arg
  integer :: i, j

  !-- Defaults
  N        = 256
  nthreads = 1
  max_iter = 200000
  tol      = 1.0d-6

  !-- Command-line arguments
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg);  read(arg,*) N
  end if
  if (command_argument_count() >= 2) then
    call get_command_argument(2, arg);  read(arg,*) nthreads
  end if

  call omp_set_num_threads(nthreads)

  write(*,'(A)')              repeat('=',56)
  write(*,'(A,I4,A,I4,A,I2)') &
    '  Parallel 2D Heat Solver  grid=',N,'x',N,'  threads=',nthreads
  write(*,'(A)')              repeat('=',56)

  !-- Allocate
  allocate(T(0:N+1, 0:N+1), Tnew(0:N+1, 0:N+1), stat=ios)
  if (ios /= 0) stop '*** Allocation failed ***'

  !-- Initialise
  T    = 0.0d0
  Tnew = 0.0d0

  !-- Dirichlet BCs
  T(0:N+1, N+1)    = 100.0d0
  Tnew(0:N+1, N+1) = 100.0d0

  !-- Jacobi iteration – timed with omp_get_wtime for wall-clock time
  t_start = omp_get_wtime()

  iter = 0
  err  = huge(1.0d0)

  do while (err > tol .and. iter < max_iter)
    iter = iter + 1
    err  = 0.0d0

    !--------------------------------------------------------------------------
    ! Outer j-loop distributed across threads.
    ! Each thread writes to Tnew(i, j_chunk) while reading T — no race hazard
    ! because Tnew and T are different arrays (pure Jacobi, no Gauss-Seidel).
    ! REDUCTION(max:err) ensures the per-thread max-residuals are combined.
    !--------------------------------------------------------------------------
    !$OMP PARALLEL DO PRIVATE(i) REDUCTION(max:err) SCHEDULE(static)
    do j = 1, N
      do i = 1, N
        Tnew(i,j) = 0.25d0 * ( T(i-1,j) + T(i+1,j) &
                              + T(i,j-1) + T(i,j+1) )
        err = max(err, abs(Tnew(i,j) - T(i,j)))
      end do
    end do
    !$OMP END PARALLEL DO

    ! Copy updated values back to T (parallelise the memory-copy too)
    !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static)
    do j = 1, N
      do i = 1, N
        T(i,j) = Tnew(i,j)
      end do
    end do
    !$OMP END PARALLEL DO
  end do

  t_end   = omp_get_wtime()
  elapsed = t_end - t_start

  write(*,'(A,I8)')       '  Iterations : ', iter
  write(*,'(A,ES12.4)')   '  Residual   : ', err
  write(*,'(A,F12.6,A)')  '  Wall time  : ', elapsed, ' s'
  write(*,'(A,F8.2)')     '  Threads    : ', real(nthreads, dp)

  !-- Append to timing CSV: N, threads, wall_time, iterations
  call execute_command_line('mkdir -p data')
  open(unit=20, file='data/parallel_timing.csv', position='append', &
       status='unknown', iostat=ios)
  write(20,'(I0,",",I0,",",F18.8,",",I0)') N, nthreads, elapsed, iter
  close(20)

  deallocate(T, Tnew)
end program heat_parallel
