program test_lookuptable

    use ModLookUpTable
    implicit none

    call ModLookUpTable_Read('/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02.dat',logical_unit=1)
    call ModLookUpTable_Read('/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02_entropy_contour_9.252752.dat',logical_unit=1)
    call ModLookUpTable_Read('/Users/xianyu/Documents/research/fortran/Falcon/3_6/RadCoolPhoto_8.0.dat',logical_unit=1)
    call ModLookUpTable_Read('/Users/xianyu/Documents/research/fortran/Falcon/3_6/opacity_x70_y28_z2.dat',logical_unit=1)

    write(*,*) 'LookUpTable_ptr%sizes = ', LookUpTables(4)%sizes
    write(*,*) 'LookUpTable_ptr%dimension = ', LookUpTables(4)%dimension
    write(*,*) size(LookUpTables(4)%data_3D)
    print *,LookUpTables(4)%data_3D(:,1,1)
    print *,LookUpTables(4)%data_3D(1,:,2)
    print *,LookUpTables(4)%name

    !print *,1

    !print *,lookuptables(1)%data_3D(:,100,1)
end program test_lookuptable