
module bous_module

    ! Module parameters:

    logical      :: bouss  ! Turn on the dispersive terms
    real(kind=8) :: B_param 
    logical      :: use_bous_sw_thresh 
    real(kind=8) :: bous_sw_thresh  ! Threshold for the transition to SWE
    real(kind=8) :: sw_depth
    real(kind = 8 ), allocatable, dimension (:,:) :: A_inverse

    save
    
contains

      !----------------------
      subroutine matvec_triad ( n, x, y, nelt, ia, ja, a, isym )                       !
      !
      !! MATVEC_TRIAD computes A*X for a matrix A stored in SLAP Triad form.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    21 July 2004
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input, integer N, the number of elements in the vectors.
      !
      !    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
      !
      !    Output, real ( kind = 8 ) Y(N), the product A * X.
      !
      !    Input, integer NELT, the number of nonzero entries in A.
      !
      !    Input, integer IA(NELT), JA(NELT), real ( kind = 8 ) A(NELT), the data
      !    structure storing the sparse matrix.
      !
      !    Input, integer ISYM, is 0 if all nonzero entries of the matrix
      !    are stored, and 1 if only the diagonal and upper or lower triangle
      !    are stored.
      !
        implicit none

        integer ( kind = 4 ) n
        integer ( kind = 4 ) nelt

        real ( kind = 8 ) a(nelt)
        integer ( kind = 4 ) ia(nelt)
        integer ( kind = 4 ) isym
        integer ( kind = 4 ) ja(nelt)
        integer ( kind = 4 ) k
        real ( kind = 8 ) x(n)
        real ( kind = 8 ) y(n)

        y(1:n) = 0.0D+00

        do k = 1, nelt
          y(ia(k)) = y(ia(k)) + a(k) * x(ja(k))
        end do

        return
      end
      
!====================     
      subroutine msolve_identity ( n, r, z, nelt, ia, ja, a, isym, rwork, iwork) !*****************************************************************************80
      !
      !! MSOLVE_IDENTITY applies the identity matrix preconditioner.
      !
      !  Discussion:
      !
      !    Most SLAP solver routines require a preconditioner routine
      !    that can solve M * Z = R.  If no preconditioning is required,
      !    then you can simply behave as though the preconditioning matrix
      !    M was the identity matrix.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    21 July 2004
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) N, the number of elements in the vectors.
      !
      !    Input, real ( kind = 8 ) R(N), the right hand side.
      !
      !    Output, real ( kind = 8 ) Z(N), the solution of M * Z = R.
      !
      !    Input, integer ( kind = 4 ) NELT, the number of nonzero entries in A.
      !
      !    Input, integer ( kind = 4 ) IA(NELT), JA(NELT), real ( kind = 8 ) A(NELT), 
      !    the data structure storing the sparse matrix.
      !
      !    Input, integer ( kind = 4 ) ISYM, is 0 if all nonzero entries of the matrix
      !    are stored, and 1 if only the diagonal and upper or lower triangle
      !    are stored.
      !
      !    Input, real ( kind = 8 ) RWORK(*), a real array that
      !    can be used to pass information to the preconditioner.
      !
      !    Input, integer ( kind = 4 ) IWORK(*), an integer array that
      !    can be used to pass information to the preconditioner.
      !
        implicit none

        integer ( kind = 4 ) n
        integer ( kind = 4 ) nelt

        real ( kind = 8 ) a(nelt)
        integer ( kind = 4 ) ia(nelt)
        integer ( kind = 4 ) isym
        integer ( kind = 4 ) iwork(*)
        integer ( kind = 4 ) ja(nelt)
        real ( kind = 8 ) r(n)
        real ( kind = 8 ) rwork(*)
        real ( kind = 8 ) z(n)

        z(1:n) = r(1:n)

        return
      end

end module bous_module
