module fortepianoWrapper
    use precision
    use variables
    use fpCosmology
    implicit None

    contains

    !functions for const
    subroutine storeNuDensMat(vec)
        real(kind=kind(0.d0)), dimension(:), intent(in) :: vec
        call vec_2_densMat(vec)
    end subroutine storeNuDensMat

    !functions for cosmology
    function w_photonNumberDensity(z) result(r)
        real(kind=kind(0.d0)), intent(in) :: z
        real(kind=kind(0.d0)) :: r
        r=photonNumberDensity(z)
    end function w_photonNumberDensity

    function w_photonDensity(z) result(r)
        real(kind=kind(0.d0)), intent(in) :: z
        real(kind=kind(0.d0)) :: r
        r=photonDensity(z)
    end function w_photonDensity
end module fortepianoWrapper
