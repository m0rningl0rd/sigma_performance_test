! A very simple program to test performance of different singular value
! decomposition methods.
! It is used for sigma subgrid model in LES, see doi:10.1063/1.3623274
! Each method is executed n_trials times, and total time is printed for each.
! Method 1: Using LAPACK sgesvd to directly compute singular values
! Method 2: Using LAPACK ssyev to compute eigenvalues, sigma = sqrt(lambda)
! Method 3: Using LAPACK sgesdd to directly compute singular values
! Method 4: Using self-contained algebraic method given in doi:10.1006/jmre.2001.2400

program sigma_performance_test
  implicit none
  integer, parameter :: seed = 31, n_trials = 1000000
  real :: A(3,3) = 0., A1(3,3) = 0., A2(3,3) = 0., sigma(3) = 0., start, finish
  integer :: k, l_work

  ! Either use a pre-defined matrix
  A1 = reshape( (/1., 2., 3.,     &
                  4., 5., 6.,     &
                  7., 8., 9./),   & 
                  shape(A1),      & 
                  order=(/2,1/) )

  ! Or create a random one based on pre-defined seed
  call random_number(A2)
  
  ! Computations will use whatever is contained in "A"
  A = A2

  write(*,*) "*** Computing singular values of matrix A_ij: ***"
  write(*,*) "---------------------------------------------------"
  call print_matrix(A)
  write(*,*) "---------------------------------------------------"
  write(*,'(i,1x,a)') n_trials, "times execution for each method"
  write(*,*) "---------------------------------------------------"
  write(*,'(a)') "*** METHOD 1: Directly with LAPACK sgesvd"

  l_work = sgesvd_l_work(3)
  call cpu_time(start)
  do k = 1, n_trials
    call method_1(A, sigma, l_work)
  enddo
  call print_vector(sigma)
  call cpu_time(finish)
  write(*,'(a,1x,f7.1, 1x, a)') "Time elapsed:", 1.e3*(finish - start), "ms."

  write(*,*) 
  write(*,'(a)') "*** METHOD 2: From sqrt of eigenvalues via LAPACK ssyev"
  l_work = ssyev_l_work(3)
  call cpu_time(start)
  do k = 1, n_trials
    call method_2(A, sigma, l_work)
  enddo
  call print_vector(sigma)
  call cpu_time(finish)
  write(*,'(a,1x,f7.1, 1x, a)') "Time elapsed:", 1.e3*(finish - start), "ms."

  write(*,*) 
  write(*,'(a)') "*** METHOD 3: Directly with LAPACK sgesdd"
  l_work = sgesdd_l_work(3)
  call cpu_time(start)
  do k = 1, n_trials
    call method_3(A, sigma, l_work)
  enddo
  call print_vector(sigma)
  call cpu_time(finish)
  write(*,'(a,1x,f7.1, 1x, a)') "Time elapsed:", 1.e3*(finish - start), "ms."

  write(*,*) 
  write(*,'(a)') "*** METHOD 4: Self contained algebraic method"
  call cpu_time(start)
  do k = 1, n_trials
    call method_4(A, sigma)
  enddo
  call print_vector(sigma)
  call cpu_time(finish)
  write(*,'(a,1x,f7.3, 1x, a)') "Time elapsed:", 1.e3*(finish - start), "ms."

contains

  subroutine print_matrix(a)
    implicit none
    real, dimension(:,:), intent(in) :: a
    integer :: i, j, m, n
    m = size(a,1)
    n = size(a,2)
    do i = 1, m
      write(*,'("|",<m>(e15.8,1x,"|"))') (a(i,j), j = 1,n)
    enddo

  end subroutine print_matrix

  subroutine print_vector(v)
    implicit none
    real, dimension(:), intent(in) :: v
    integer                        :: i, m

    m  = size(v)
    write(*,'("[",<m-1>(e15.8,","),e15.8,"]")') ( v(i), i = 1,m )

  end subroutine print_vector

  function sgesvd_l_work(n) result (l_work)
    implicit none
    integer             :: l_work
    integer, intent(in) :: n
    real                :: a(n,n), s(n), work
    integer             :: info

    call sgesvd('n', 'n', n, n, a, n, s, a, n, a, n, work, -1, info)
    l_work = nint(work)

  end function sgesvd_l_work

  function sgesdd_l_work(n) result (l_work)
    implicit none
    integer             :: l_work
    integer, intent(in) :: n
    real                :: a(n,n), s(n), work
    integer             :: info

    call sgesdd('n', n, n, a, n, s, a, n, a, n, work, -1, l_work, info)
    l_work = nint(work)

  end function sgesdd_l_work

  function ssyev_l_work(n) result (l_work)
    implicit none
    integer             :: l_work
    integer, intent(in) :: n
    real                :: a(n,n), s(n), work
    integer             :: info

    call ssyev('n', 'upper', 3, a, 3, s, work, -1, info)
    l_work = nint(work)

  end function ssyev_l_work

  subroutine method_1(g, sigma, l_work)
    implicit none
    real, intent(in)    :: g(3,3)
    real, intent(out)   :: sigma(3)
    integer, intent(in) :: l_work
    real                :: a(3,3), work(l_work)
    integer             :: info

    sigma(:) = 0.
    a(:,:) = g(:,:)
    call sgesvd('n', 'n', 3, 3, a, 3, sigma, a, 3, a, 3, work, l_work, info)
    if (info > 0) then
      write(*,*) 'sgesvd failed!'
    endif

  end subroutine method_1

  subroutine method_2(g, sigma, l_work)
    implicit none
    real, intent(in)    :: g(3,3)
    real, intent(out)   :: sigma(3)
    integer, intent(in) :: l_work
    real                :: work(l_work)
    integer             :: info

    call ssyev('n', 'upper', 3, matmul(transpose(g),g), 3, sigma, work, l_work, info)
    sigma(1:3) = sigma(3:1:-1)
    sigma = sqrt(sigma)

  end subroutine method_2

  subroutine method_3(g, sigma, l_work)
    implicit none
    real, intent(in)    :: g(3,3)
    real, intent(out)   :: sigma(3)
    integer, intent(in) :: l_work
    real                :: work(l_work)
    real                :: a(3,3)
    integer             :: info, iwork(24)

    sigma(:) = 0.
    a(:,:) = g(:,:)
    call sgesdd('n', 3, 3, a, 3, sigma, a, 3, a, 3, work, l_work, iwork, info)
    if (info > 0) then
      write(*,*) 'sgesdd failed!'
    endif

  end subroutine method_3

  subroutine method_4(g, sigma)
    implicit none
    real, intent(in)  :: g(3,3)
    real, intent(out) :: sigma(3)
    real, parameter   :: pi = 3.141593
    real              :: gg(3,3), gg2(3,3), I1, I2, I3, alpha1, alpha2, alpha3

    gg  = matmul(transpose(g),g)
    gg2 = matmul(gg,gg)

    I1  = gg(1,1) + gg(2,2) + gg(3,3)
    I2  = 0.5*(I1*I1 - (gg2(1,1) + gg2(2,2) + gg2(3,3)))
    I3  = m33det(gg)

    alpha1 = I1*I1/9. - I2/3.
    alpha2 = I1*I1*I1/27. - I1*I2/6. + I3/2.
    alpha3 = (acos(alpha2/alpha1**(3./2.)))/3.

    ! Below causes NaN with matrix A1, you can test
    ! sigma(1) = sqrt(I1/3. + 2.*sqrt(alpha1)*cos(alpha3))
    ! sigma(2) = sqrt(I1/3. - 2.*sqrt(alpha1)*cos(pi/3. + alpha3))
    ! sigma(3) = sqrt(max(0., I1/3. - 2.*sqrt(alpha1)*cos(pi/3. - alpha3)))

    sigma(1) = sqrt(max(0., I1/3. + 2.*sqrt(alpha1)*cos(alpha3)))
    sigma(2) = sqrt(max(0., I1/3. - 2.*sqrt(alpha1)*cos(pi/3. + alpha3)))
    sigma(3) = sqrt(max(0., I1/3. - 2.*sqrt(alpha1)*cos(pi/3. - alpha3)))

  end subroutine method_4

  function m33det(a) result (det)
    implicit none
    real, intent(in) :: a(3,3)
    real             :: det

    det = a(1,1)*a(2,2)*a(3,3) &
        - a(1,1)*a(2,3)*a(3,2) &
        - a(1,2)*a(2,1)*a(3,3) &
        + a(1,2)*a(2,3)*a(3,1) &
        + a(1,3)*a(2,1)*a(3,2) &
        - a(1,3)*a(2,2)*a(3,1)

  end function m33det

end program sigma_performance_test
