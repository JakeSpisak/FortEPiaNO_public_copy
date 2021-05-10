module fortepianoWrapper
    use precision
    use variables
    use fpCosmology
    implicit None

    contains

    subroutine storeNuDensMat(vec)
        real(kind=1.d0), dimension(:), intent(in) :: vec
        vec_2_densMat(vec)
    end subroutine storeNuDensMat
end module fortepianoWrapper
