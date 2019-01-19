module math_utilities
!*Description:* This module contains solvers required for blood vessel problems
!
  use arrays, only: dp
  implicit none
  private
  public ax_cr,diagonal_pointer_cr,ilu_cr,lus_cr,mult_givens,rearrange_cr
  public unit_vector,vector_length
  public mean,variance

contains
!
!###########################################################################
!
  subroutine ax_cr ( n, ia, ja, a, x, w )
  !*Description:* Computes A*x for a matrix stored in sparse compressed row form
    implicit none

    integer ( kind = 4 ) n !the order of the system
    integer ( kind = 4 ) ia(*) !ia(n+1) row indices
    integer ( kind = 4 ) ja(*) !ja(nz_num) column indices
    real ( kind = 8 ) a(*) !a(nz_num) Matrix values
    real ( kind = 8 ) x(*) !x(n) Vector to be multiplied by A
    real ( kind = 8 ) w(*) !w(n) Value of A*x

    integer ( kind = 4 ) i
    integer ( kind = 4 ) k1
    integer ( kind = 4 ) k2

    w(1:n) = 0.0D+00

    do i = 1, n
       k1 = ia(i)
       k2 = ia(i+1) - 1
       w(i) = w(i) + dot_product ( a(k1:k2), x(ja(k1:k2)) )
    end do

    return
  end subroutine ax_cr
!
!##############################################################################
!
  subroutine ilu_cr ( n, nz_num, ia, ja, a, ua, l )
  ! *Description:* computes the incomplete LU factorization of a matrix. For a matrix
  ! stored in compressed row format.
  ! Input, integer ( kind = 4 ) UA(N), the index of the diagonal element of each row.
  ! Output, real ( kind = 8 ) L(NZ_NUM), the ILU factorization of A.  
    integer ( kind = 4 ) n
    integer ( kind = 4 ) nz_num
    integer ( kind = 4 ) ia(*) !ia(n+1)
    integer ( kind = 4 ) ja(*) !ja(nz_num)
    real ( kind = 8 ) a(*) !a(nz_num)
    integer ( kind = 4 ) ua(*) !ua(n)
    real ( kind = 8 ) l(*) !l(nz_num)

    integer ( kind = 4 ) i
    integer ( kind = 4 ) iw(n)
    integer ( kind = 4 ) j
    integer ( kind = 4 ) jj
    integer ( kind = 4 ) jrow
    integer ( kind = 4 ) jw
    integer ( kind = 4 ) k
    real ( kind = 8 ) tl


    !  Copy A.
    l(1:nz_num) = a(1:nz_num)

    do i = 1, n ! for each row, up to max number of rows
       !  IW points to the nonzero entries in row I.
       iw(1:n) = -1
       do k = ia(i), ia(i+1) - 1 !for each
          iw(ja(k)) = k
       end do
       do j = ia(i), ia(i+1) - 1
          jrow = ja(j)
          if ( i <= jrow ) then
             exit
          end if
          tl = l(j) * l(ua(jrow))
          l(j) = tl
          do jj = ua(jrow) + 1, ia(jrow+1) - 1
             jw = iw(ja(jj))
             if ( jw /= -1 ) then
                l(jw) = l(jw) - tl * l(jj)
             end if
          end do
       end do
       ua(i) = j
       if ( jrow /= i ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'ILU_CR - Fatal error!'
          write ( *, '(a)' ) '  JROW ~= I'
          write ( *, '(a,i8)' ) '  JROW = ', jrow
          write ( *, '(a,i8)' ) '  I    = ', i
          stop
       end if
       if ( l(j) == 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'ILU_CR - Fatal error!'
          write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
          write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
          stop
       end if
       l(j) = 1.0D+00 / l(j)
    end do

    l(ua(1:n)) = 1.0D+00 / l(ua(1:n))

    return
  end subroutine ilu_cr
!
!##############################################################################
!
subroutine diagonal_pointer_cr ( n, ia, ja, ua )
!*Description:* finds diagonal entries in a sparse compressed row matrix.
! The array UA can be used to locate the diagonal elements of the matrix.
! It is assumed that every row of the matrix includes a diagonal element,
! and that the elements of each row have been ascending sorted.
    integer ( kind = 4 ) n
    integer ( kind = 4 ) ia(*) !ia(n+1)
    integer ( kind = 4 ) ja(*) !ja(nz_num)
    integer ( kind = 4 ) ua(*) !ua(n)

    integer ( kind = 4 ) i
    integer ( kind = 4 ) k

    ua(1:n) = -1

    do i = 1, n
       do k = ia(i), ia(i+1) - 1
          if ( ja(k) == i ) then
             ua(i) = k
          end if
       end do
    end do
    return
  end subroutine diagonal_pointer_cr

  !*****************************************************************************80

  subroutine lus_cr ( n, ia, ja, l, ua, r, z )
!*Description:* this subroutine applies the incomplete LU preconditioner.
! The linear system M * Z = R is solved for Z.  M is the incomplete
! LU preconditioner matrix, and R is a vector supplied by the user.
! So essentially, we're solving L * U * Z = R.
! Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
! of each row.
! Input, real ( kind = 8 ) R(N), the right hand side.
! Output, real ( kind = 8 ) Z(N), the solution of the system M * Z = R.
    implicit none

    integer ( kind = 4 ) n
    integer ( kind = 4 ) ia(*) !ia(n+1)
    integer ( kind = 4 ) ja(*) !ja(nz_num)
    real ( kind = 8 ) l(*) !l(nz_num)
    integer ( kind = 4 ) ua(*) !ua(n)
    real ( kind = 8 ) r(*) !r(n)

    integer ( kind = 4 ) i
    integer ( kind = 4 ) j
    real ( kind = 8 ) w(n)
    real ( kind = 8 ) z(n)

    !  Copy R in.
    w(1:n) = r(1:n)

    !  Solve L * w = w where L is unit lower triangular.
    do i = 2, n
       do j = ia(i), ua(i) - 1
          w(i) = w(i) - l(j) * w(ja(j))
       end do
    end do

    !  Solve U * w = w, where U is upper triangular.
    do i = n, 1, -1
       do j = ua(i) + 1, ia(i+1) - 1
          w(i) = w(i) - l(j) * w(ja(j))
       end do
       w(i) = w(i) / l(ua(i))
    end do

    !  Copy Z out.
    z(1:n) = w(1:n)

    return
  end subroutine lus_cr

  !*****************************************************************************80
  subroutine mult_givens ( c, s, k, g )
!*Description:* This subroutine applies a Givens rotation to two successive entries of a vector.
! In order to make it easier to compare this code with the Original C, the vector indexing is 0-based.
! Input, real ( kind = 8 ) C, S, the cosine and sine of a Givens rotation.
!
! Input, integer ( kind = 4 ) K, indicates the location of the first vector entry.
!
! Input/output, real ( kind = 8 ) G(1:K+1), the vector to be modified.
! On output, the Givens rotation has been applied to entries G(K) and G(K+1).

    implicit none

    real ( kind = 8 ) c
    real ( kind = 8 ) s
    integer ( kind = 4 ) k
    real ( kind = 8 ) g(*) !g(1:k+1)

    real ( kind = 8 ) g1
    real ( kind = 8 ) g2

    g1 = c * g(k) - s * g(k+1)
    g2 = s * g(k) + c * g(k+1)

    g(k)   = g1
    g(k+1) = g2

    return
  end subroutine mult_givens


  !*****************************************************************************80
  subroutine rearrange_cr ( n, ia, ja, a )
!*Description:* This subroutine sorts a sparse compressed row matrix.
! It guarantees that the entries in the CR matrix are properly sorted.
!
! After the sorting, the entries of the matrix are rearranged in such
! a way that the entries of each column are listed in ascending order
! of their column values.
! Input, integer ( kind = 4 ) N, the order of the system.
!
! Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
! Input, integer ( kind = 4 ) IA(N+1), the compressed row indices.
!
! Input/output, integer ( kind = 4 ) JA(NZ_NUM), the column indices.
! On output, these may have been rearranged by the sorting.
!
! Input/output, real ( kind = 8 ) A(NZ_NUM), the matrix values.  On output,
!    the matrix values may have been moved somewhat because of the sorting.
!
    implicit none

    integer ( kind = 4 ) n
    integer ( kind = 4 ) ia(*) !ia(n+1)
    integer ( kind = 4 ) ja(*) !ja(nz_num)
    real ( kind = 8 ) a(*) !a(nz_num)

    integer ( kind = 4 ) i
    integer ( kind = 4 ) i4temp
    integer ( kind = 4 ) k
    integer ( kind = 4 ) l
    real ( kind = 8 ) r8temp

    do i = 1, n

       do k = ia(i), ia(i+1) - 2
          do l = k + 1, ia(i+1) - 1

             if ( ja(l) < ja(k) ) then
                i4temp = ja(l)
                ja(l)  = ja(k)
                ja(k)  = i4temp

                r8temp = a(l)
                a(l)   = a(k)
                a(k)   = r8temp
             end if

          end do
       end do

    end do

    return
  end subroutine rearrange_cr
  !!!###############################################################

  function unit_vector(A)

    !###  Calculates the unit vector for an arbitrary 3x1 vector

    real(dp),intent(in) :: A(*)
    real(dp) :: length_a,unit_vector(3)

    length_a = vector_length(A)
    if(length_a.gt.1.0e-6_dp)then
       unit_vector(1:3) = A(1:3)/length_a
    else
       WRITE(*,*) ' >>WARNING: Cannot normalise a zero length vector'
       WRITE(*,*) ' We recommend debugging, but hit enter to continue'
       read(*,*)
    endif

  end function unit_vector

  !!!##################################################

  function vector_length(A)

    !###  Calculates the length of a 3x1 vector

    real(dp),intent(in) :: A(*)
    real(dp) :: vector_length
    integer :: i

    vector_length = 0.0_dp
    do i=1,3
       vector_length = vector_length + A(i)*A(i)
    enddo
    vector_length = dsqrt(vector_length)

  end function vector_length

    !!!##################################################

  function mean(vec_len,vector)

    !###  Calculates the mean value of entries in a vector

    integer, intent(in) :: vec_len
    real(dp),intent(in) :: vector(vec_len)
    real(dp) :: mean

    mean = sum(vector)/vec_len

  end function mean

      !!!##################################################

  function variance(vec_len,vector,mean_value)

    !###  Calculates the mean value of entries in a vector

    integer, intent(in) :: vec_len
    real(dp),intent(in) :: vector(vec_len)
    real(dp), intent(in) :: mean_value

    real(dp) :: variance
    integer :: i

    variance = 0.0_dp
    do i = 1,vec_len
       variance = variance + (vector(i) - mean_value)**2.0_dp
    enddo
    variance = variance / (vec_len - 1)

  end function variance


end module math_utilities
