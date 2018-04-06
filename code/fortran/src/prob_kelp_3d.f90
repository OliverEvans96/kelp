! -*- f90 -*-
subroutine probkelp(arr, n)
    
    double precision, dimension(n,n) :: arr
    integer i, j, n
    
    !f2py intent(in) n
    !f2py intent(inout) arr
    
    do i=1,n
        do j=1,n
            if (mod(i+j, 2) .eq. 0) then
                arr(i,j) = 3
            else
                arr(i,j) = 1
            end if
        end do
    end do
end subroutine
