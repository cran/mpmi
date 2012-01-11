! Subroutine to calculate mixed-pair MI value with
! jackknife bias correction and z-value.
subroutine mmipwnjk(cts, lc, disc, h, ans) !, mps, zvalue)
    use iface
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Input variables:
    ! Length of input vectors
    integer, intent(in) :: lc

    ! Vector for discrete variable (group membership)
    integer, dimension(lc), intent(in) :: disc

    ! Bandwidth for kernel density estimation
    real(kind=rdble), intent(in) :: h

    ! Continuous input vector
    real(kind=rdble), dimension(lc), intent(in) :: cts
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Output variables:
    ! ans = raw MI value
    ! mps = bias corrected mi value
    ! zvalue = approximate z value for hypothesis that 
    ! mps == 0
    real(kind=rdble), intent(out) :: ans !, mps, zvalue
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local variables:
    ! Sums of kernel distances (see below)
    real(kind=rdble), dimension(lc) :: t2, t3

    ! Temporary variable for calculating kernel distances
    real(kind=rdble) t1 

    ! Loop indices
    integer :: i, j, k 
    ! Number of groups in discrete variable
    integer :: levd

    ! Kernel distance matrix
    real(kind=rdble), dimension(lc, lc) :: kmat

    ! Jackknifed MI scores
    ! real(kind=rdble), dimension(lc) :: ansjk

    ! Jackknife pseudo values
    ! real(kind=rdble), dimension(lc) :: ps

    ! Standard deviation of pseudo values
    ! real(kind=rdble) :: sdps

    ! Sums of kernel distances within jackknife
    ! (I.e., with kth observation removed)
    ! real(kind=rdble) :: t22, t32

    ! Dynamic arrays:
    ! Table of discrete variable
    integer, dimension(:), allocatable :: tab

    ! Above table scaled to probabilities / relative frequencies
    real(kind=rdble), dimension(:), allocatable :: ptab
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ans = 0.0

    ! Number of levels of factor (coded with integers)
    levd = maxval(disc)

    ! Initialise table
    allocate(tab(levd))
    tab = 0

    ! Tabulate discrete variable
    do i = 1, lc
        tab(disc(i)) = tab(disc(i)) + 1
    end do

    ! Probability table
    allocate(ptab(levd))
    do i = 1, levd
        ptab(i) = dble(tab(i)) / dble(lc)
    end do

    ! Matrix of kernel distances (inefficient - should probably pack into vector)
    kmat = 0.0 
    t1 = 0.0
    do i = 1, lc
        do j = i + 1, lc
            ! Epanechnikov kernel
            t1 = (cts(j) - cts(i)) / h
            if (abs(t1) .ge. 1.0) then
                t1 = 0.0
            else
                t1 = 1.0 - (t1 * t1) 
            end if
            kmat(i, j) = t1

            ! Symmetrise
            kmat(j, i) = kmat(i, j)
        end do
        kmat(i, i) = kmat(i, i) + 1.0
    end do

    ! Sum of kernel distances from cts(i) to all points
    t2 = 0.0
    ! Sum of kernel distances from cts(i) to all points in same subgroup
    t3 = 0.0

    ! Evaluate non-jackknifed MI and fill t2 and t3
    do i = 1, lc
        do j = i + 1, lc
            t2(i) = t2(i) + kmat(i, j)

            ! Using kernel symmetry
            t2(j) = t2(j) + kmat(i, j)

            if (disc(j) == disc(i)) then
                t3(i) = t3(i) + kmat(i, j)
                ! Using kernel symmetry
                t3(j) = t3(j) + kmat(i, j)
            end if
        end do

        ! For when i == j
        t2(i) = t2(i) + 1.0
        t3(i) = t3(i) + 1.0

        ! Accumulate MI
        ans = ans + ptab(disc(i)) * log(lc * t3(i) / (tab(disc(i)) * t2(i))) / tab(disc(i))
    end do

    ! ansjk = 0.0
    ! do k = 1, lc
    !     ! Remove kth observation from table of counts
    !     ! (Because tab(disc(i)) may equal tab(disc(k)) below)
    !     tab(disc(k)) = tab(disc(k)) - 1

    !     do i = 1, lc
    !         ! Exclude kth observation
    !         if (i .ne. k) then
    !             ! Subtracting kernel distances to kth observation (as per
    !             ! jackknife)
    !             t22 = t2(i) - kmat(k, i)

    !             if (disc(i) == disc(k)) then
    !                 t32 = t3(i) - kmat(k, i)
    !             else
    !                 t32 = t3(i)
    !             end if
    !    
    !             ! Accumulate kth jackknife MI value
    !             ansjk(k) = ansjk(k) + ((tab(disc(i)))/(lc-1.0)) * log((lc-1.0) * t32 / ((tab(disc(i))) * t22)) / (tab(disc(i)))
    !         end if
    !     end do

    !     ! Put kth observation back in table of counts
    !     tab(disc(k)) = tab(disc(k)) + 1
    ! end do

    ! Get bias corrected MI value and z-value from the jackknife
    ! using Tukey's pseudo-value approach.
    ! (There are probably more efficient ways to do this.)
    ! ps = dble(lc) * ans - (dble(lc) - 1.0) * ansjk

    ! mps = sum(ps) / dble(lc)
    ! sdps = sqrt(sum((ps - mps) * (ps - mps)) / (dble(lc) - 1.0))
    ! zvalue = sqrt(dble(lc)) * mps / sdps

    deallocate(tab)
    deallocate(ptab)
end subroutine

subroutine mmimnjk(cdat, nrc, ncc, sdat, nrs, ncs, mis, h)
    use iface
    implicit none

    ! Input variables
    integer, intent(in) :: nrc, ncc, nrs, ncs !, switch
    ! integer, intent(in), optional :: m
    real(kind=rdble), dimension(nrc, ncc), intent(in) :: cdat
    real(kind=rdble), dimension(ncc), intent(in), optional :: h
    integer, dimension(nrs, ncs), intent(in) :: sdat

    ! Output matrix (now holds bias corrected estimates)
    real(kind=rdble), dimension(ncc, ncs), intent(out) :: mis
    ! real(kind=rdble), dimension(ncc, ncs), intent(out) :: bcmis
    ! Matrix of z-values
    ! real(kind=rdble), dimension(ncc, ncs), intent(out) :: zmat

    ! Arrays to hold non-missing observations only
    ! Reuse 'static' arrays for speed
    real(kind=rdble), dimension(nrc) :: cvec
    integer, dimension(nrs) :: svec

    ! Local variables
    integer :: i, j, k, nok
    logical, dimension(nrc) :: ok

    ! Function to get R's code for missing integers
    integer :: rnaint
    ! Holds R's coding for missing integers (INT_MIN)
    integer :: naint

    ! R function to check real missing values
    integer :: rfinite

    naint = rnaint() ! Asks R for its missing integer coding

    !$omp parallel do default(none) shared(ncc, ncs, cdat, sdat, &
    !$omp nrc, naint, h, mis)  &
    !$omp private(ok, nok, cvec, svec, i, j) &
    !$omp schedule(dynamic)
    do i = 1, ncc
        do j = 1, ncs
            if (i <= ncc .and. j <= ncs) then
                ! Remove missing values pairwise
                do k = 1, nrc
                    if (rfinite(cdat(k,i)) == 1 .and. sdat(k,j) /= naint) then
                        ok(k) = .true.
                    else
                        ok(k) = .false.
                    end if
                end do

                nok = count(ok)

                ! Pack non-missing values
                cvec = pack(cdat(:,i), mask = ok)
                svec = pack(sdat(:,j), mask = ok)

                ! Call pairwise mixed MI subroutine
                call mmipwnjk(cvec(1:nok), nok, svec(1:nok), h(i), mis(i,j))
            end if
        end do
    end do
    !$omp end parallel do
end subroutine

