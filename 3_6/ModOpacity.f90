module ModOpacity

    use ModLookUpTable, only: LookUpTables,LookUpTable,nLookUpTables

    implicit none

    type(LookUpTable),pointer   ::  Opacity_table
    integer                     ::  nlogT_opacity,nlogR_opacity
    real(8),allocatable         ::  logT_opacity(:)
    real(8),allocatable         ::  logR_opacity(:)
    contains

    subroutine ModOpacity_init
        implicit none
        integer ::  iLookUpTable

        do iLookUpTable=1,nLookUpTables
            if (LookUpTables(iLookUpTable)%name == 'OPACITY') then
                Opacity_table=>LookUpTables(iLookUpTable)

                ! Get sizes
                nlogT_opacity = Opacity_table%sizes(1)
                nlogR_opacity = Opacity_table%sizes(2)

                ! Get logT and logR grids of the opacity table.
                allocate(logT_opacity(nlogT_opacity))
                allocate(logR_opacity(nlogR_opacity))
                logT_opacity = Opacity_table%data_3D(:,1,1)
                logR_opacity = Opacity_table%data_3D(1,:,2)
                exit
            end if
        end do
    end subroutine ModOpacity_init
end module ModOpacity