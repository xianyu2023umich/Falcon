program test2

    use ModLinReconstruct


    implicit none

    integer             ::  nr=100,nt=100,np=100,ng=2,nvar=5
    integer             ::  ir,it,ip,test
    real,allocatable    ::  primitive(:,:,:,:),d_primitive(:,:,:,:)

    allocate(primitive(-ng+1:nr+ng,-ng+1:nt+ng,-ng+1:np+ng,nvar))
    allocate(d_primitive(-ng+2:nr+ng-1,-ng+2:nt+ng-1,-ng+2:np+ng-1,5))

    call random_number(primitive)
    print *,1
    call ModLinReconstruct_minmod(nvar,nr,nt,np,ng,1,primitive,d_primitive)
    print *,2
    print *,d_primitive(10,10,10,:)
    print *,3
    call ModLinReconstruct_minmod_new(nvar,nr,nt,np,ng,1,primitive,d_primitive)
    print *,4
    print *,d_primitive(10,10,10,:)
    print *,5
    


    

end program test2
