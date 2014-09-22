! 
!     Copyright 2013 Chris Pardy <cpardy@unsw.edu.au>
! 
!     This file is part of the mpmi R package.
! 
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, version 3.
! 
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.
!


! Epanechnikov kernel
! Pairwise only
! NO Jackknife
subroutine cmipwnjk(v1, v2, lv, h1, h2, ans)
    use iface
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Input variables:
    
    ! Length of vectors 
    integer, intent(in) :: lv
    ! Data vectors
    real(kind=rdble), dimension(lv), intent(in) :: v1, v2
    ! Smoothing bandwidths in each dimension
    ! (corresponding to v1 and v2 respectively)
    real(kind=rdble), intent(in) :: h1, h2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Output variables:
    ! ans = raw MI
    real(kind=rdble), intent(out) :: ans 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local variables:
    ! Loop indices
    integer :: i, j, k

    ! Temporary variables for calculating kernel matrix
    real(kind=rdble) :: t1, t2

    ! Sums of kernel distances for each point (of lv points)
    !
    ! s1 & s2 hold sums of kernel distances from each point
    ! to all other points
    !
    ! s12 holds the sums of product kernels for each point
    real(kind=rdble), dimension(lv) :: s1, s2, s12

    ! Kernel matrices for vectors 1 & 2
    real(kind=rdble), dimension(lv, lv) :: kmat1, kmat2 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ans = 0.0

    ! Pre-calculate kernel distances
    ! Inefficient matrix of kernel distances (should probably pack into vector)
    kmat1 = 0.0 
    kmat2 = 0.0 
    t1 = 0.0
    t2 = 0.0

    ! Separate loops hopefully help cache locality
    ! Vector 1:
    do i = 1, lv
        do j = i + 1, lv
            ! Epanechnikov kernel
            t1 = (v1(j) - v1(i)) / h1
            if (abs(t1) .ge. 1.0) then
                t1 = 0.0
            else
                t1 = 1.0 - (t1 * t1) 
            end if
            kmat1(i, j) = t1

            ! Symmetrise
            kmat1(j, i) = kmat1(i, j)
        end do
        kmat1(i, i) = kmat1(i, i) + 1.0
    end do
    ! Vector 2:
    do i = 1, lv
        do j = i + 1, lv
            ! Epanechnikov kernel
            t2 = (v2(j) - v2(i)) / h2
            if (abs(t2) .ge. 1.0) then
                t2 = 0.0
            else
                t2 = 1.0 - (t2 * t2) 
            end if
            kmat2(i, j) = t2

            ! Symmetrise
            kmat2(j, i) = kmat2(i, j)
        end do
        kmat2(i, i) = kmat2(i, i) + 1.0
    end do

    s1 = 0.0
    s2 = 0.0
    s12 = 0.0

    ! N.B., this uses the simple 'product kernel'
    ! approach for 2D kernel density estimation
    do i = 1, lv
        do j = i + 1, lv

            s1(i) = s1(i) + kmat1(i,j)
            s2(i) = s2(i) + kmat2(i,j)

            ! Use product kernel for joint distribution
            s12(i) = s12(i) + kmat1(i,j) * kmat2(i,j)

            ! Using kernel symmetry
            s1(j) = s1(j) + kmat1(i,j) 
            s2(j) = s2(j) + kmat2(i,j) 
            s12(j) = s12(j) + kmat1(i,j) * kmat2(i,j)
        end do

        ! For when i == j
        s1(i) = s1(i) + 1.0
        s2(i) = s2(i) + 1.0
        s12(i) = s12(i) + 1.0
        
        ! Accumulate raw MI value
        ans = ans + log(s12(i) / (s1(i) * s2(i)))
    end do
    ans = ans / lv + log(dble(lv))
    
end subroutine

subroutine cmimnjk(cdat, nrc, ncc, mis, h)
    use iface
    implicit none

    ! Input variables
    integer, intent(in) :: nrc, ncc
    real(kind=rdble), dimension(nrc, ncc), intent(in) :: cdat
    real(kind=rdble), dimension(ncc), intent(in) :: h

    ! Output matrices
    real(kind=rdble), dimension(ncc, ncc), intent(out) :: mis

    ! Arrays to hold non-missing observations only
    ! Reuse 'static' arrays for speed
    real(kind=rdble), dimension(nrc) :: cvec, svec

    ! Local variables
    integer :: i, j, nok, k
    logical, dimension(nrc) :: ok

    ! R function to check real missing values
    integer :: rfinite

    !$omp parallel do default(none) shared(ncc, nrc, cdat, &
    !$omp h, mis)  &
    !$omp private(ok, nok, cvec, svec, i, j) &
    !$omp schedule(dynamic)
    do i = 1, ncc
        do j = i, ncc
            ! Remove missing observations pairwise
            do k = 1, nrc
                if (rfinite(cdat(k,i)) == 1 .and. rfinite(cdat(k,j)) == 1) then
                    ok(k) = .true.
                else
                    ok(k) = .false.
                end if
            end do

            nok = count(ok)

            ! Only perform calculation if there are non-missing values
            ! in both input vectors (set to 3 for no real reason)
            if (nok > 2) then
                ! Pack non-missing values
                cvec = pack(cdat(:,i), mask = ok)
                svec = pack(cdat(:,j), mask = ok)

                ! Call pairwise continuous MI subroutine.
                call cmipwnjk(cvec(1:nok), svec(1:nok), nok, h(i), h(j), mis(i,j))
            else
                ! Set all results to zero
                mis(i, j) = 0.0
            end if

            ! Symmetrise result matrix
            if (i .ne. j) then
                mis(j, i) = mis(i, j)
            end if
        end do
    end do
    !$omp end parallel do
end subroutine


