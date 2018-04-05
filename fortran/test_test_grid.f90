program test_test_grid
    use test_grid
    integer ntheta, nphi
    double precision, dimension(:), allocatable :: theta_p, phi_p, area_p
    double precision, dimension(:,:), allocatable :: vsf_integral

    ntheta = 10
    nphi = 10
    nomega = ntheta * (nphi - 2) + 2

    allocate(theta_p(nomega))
    allocate(phi_p(nomega))
    allocate(area_p(nomega))
    allocate(vsf_integral(nomega,nomega))

    call make_vsf(ntheta, nphi, theta_p, phi_p, area_p, vsf_integral)

    deallocate(vsf_integral)
    deallocate(area_p)
    deallocate(phi_p)
    deallocate(theta_p)
end program

