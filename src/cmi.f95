! Epanechnikov kernel
! Pairwise only
! Jackknife
subroutine cmipw(v1, v2, lv, h1, h2, ans, mps, zvalue)
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
    ! mps = jackknife bias corrected MI
    ! zvalue = z value for hypothesis that mps == 0
    real(kind=rdble), intent(out) :: ans, mps, zvalue
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
    
    ! Jackknife replication MI values
    real(kind=rdble), dimension(lv) :: ansjk

    ! Jackknife pseudo-values
    real(kind=rdble), dimension(lv) :: ps

    ! Temporary variables for jackknife
    real(kind=rdble) :: ts1, ts2, ts12

    ! SD of pseudo-values
    real(kind=rdble) :: sdps

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

    ! Get jackknife estimates
    ansjk = 0.0
    do k = 1, lv
        do i = 1, lv
            ! Exclude kth observation
            if (i .ne. k) then
                ! Subtract kernel distances corresponding 
                ! to kth (excluded) observation
                ts1 = s1(i) - kmat1(k, i)
                ts2 = s2(i) - kmat2(k, i)
                ts12 = s12(i) - kmat1(k, i) * kmat2(k, i)

                ! Accumulate jackknife MI values
                ansjk(k) = ansjk(k) + log(ts12 / (ts1 * ts2))
            end if
        end do
    end do
    ansjk = ansjk / (dble(lv) - 1.0) + log(dble(lv) - 1.0)

    ! Tukey's jackknife pseudo values
    ps = dble(lv) * ans - (dble(lv) - 1.0) * ansjk

    ! Bias corrected MI
    mps = sum(ps) / dble(lv)
    ! Get z-value for hypothesis that mps == 0
    sdps = sqrt(sum((ps - mps) * (ps - mps)) / (dble(lv) - 1.0))
    zvalue = sqrt(dble(lv)) * mps / sdps
end subroutine

subroutine cmim(cdat, nrc, ncc, mis, bcmis, zmat, h)
    use iface
    implicit none

    ! Input variables
    integer, intent(in) :: nrc, ncc
    real(kind=rdble), dimension(nrc, ncc), intent(in) :: cdat
    real(kind=rdble), dimension(ncc), intent(in) :: h
    ! integer, dimension(nrs, ncs), intent(in) :: sdat

    ! Output matrices
    real(kind=rdble), dimension(ncc, ncc), intent(out) :: mis
    real(kind=rdble), dimension(ncc, ncc), intent(out) :: bcmis
    ! Matrix of z-values
    real(kind=rdble), dimension(ncc, ncc), intent(out) :: zmat

    ! Arrays to hold non-missing observations only
    ! Reuse 'static' arrays for speed
    real(kind=rdble), dimension(nrc) :: cvec, svec

    ! Local variables
    integer :: i, j, nok, k
    logical, dimension(nrc) :: ok

    ! R function to check real missing values
    integer :: rfinite

    !$omp parallel do default(none) shared(ncc, nrc, cdat, &
    !$omp h, mis, bcmis, zmat)  &
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
                call cmipw(cvec(1:nok), svec(1:nok), nok, h(i), h(j), mis(i,j), bcmis(i,j), zmat(i,j))
            else
                ! Set all results to zero
                mis(i, j) = 0.0
                bcmis(i, j) = 0.0
                zmat(i, j) = 0.0
            end if

            ! Symmetrise result matrix
            if (i .ne. j) then
                mis(j, i) = mis(i, j)
                bcmis(j, i) = bcmis(i, j)
                zmat(j, i) = zmat(i, j)
            end if
        end do
    end do
    !$omp end parallel do
end subroutine


