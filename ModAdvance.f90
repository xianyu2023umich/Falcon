module ModAdvance

    use ieee_arithmetic
    use ModEquation
    use ModGC
    use ModCommunication
    use ModBoundary
    use ModTimeStep

    contains

    subroutine ModAdvance_rk4(Tree,CFL_ad,CFL_df,if_all_same_sizes,MpiSize,MpiRank)

        implicit none

        type(OcTree),intent(inout),target   ::  Tree                        ! Tree
        real,intent(in)                     ::  CFL_ad,CFL_df               ! CFL
        logical,intent(in)                  ::  if_all_same_sizes           ! if all local blocks same size
        integer,intent(in)                  ::  MpiSize,MpiRank

        type(Block),pointer                 ::  Block1                      ! block pointer
        integer                             ::  iLocalBlock                 ! i of block
        integer                             ::  rk_index                    ! runge-kutta index
        real,pointer                        ::  primitive_i(:,:,:,:),&      ! primitive pointers
                                                primitive_i_1(:,:,:,:)
        real,pointer                        ::  EQN_update_Ri(:,:,:,:)      ! EQN-update_R pointer
        real,allocatable,target             ::  EQN_update_R1(:,:,:,:),&    ! the following four used for rk
                                                EQN_update_R2(:,:,:,:),&
                                                EQN_update_R3(:,:,:,:),&
                                                EQN_update_R4(:,:,:,:)
        
        logical                             ::  if_top,if_bottom            ! if block at top or bottom
        real                                ::  dt_local,dt                 ! dt
        integer :: i,j,k

        call ModTimeStep_Dynamo_HD(Tree,CFL_ad,CFL_df,dt_local)
        call ModCommunication_GlobalTimeStep(dt_local,dt)
          
        ! see if all the local blocks in the tree
        ! have the same block size. If yes then it
        ! there is no need to allocate and deallocate
        ! every time.

        if (if_all_same_sizes) then


            ! allocate the Runge-kutta R1 to R4 according to the Block1 size.
            !
            Block1 => Tree%LocalBlocks(1)
            allocate(EQN_update_R1(1:Block1%nvar,1:Block1%ni,1:Block1%nj,1:Block1%nk))
            !allocate(EQN_update_R2(1:Block1%nvar,1:Block1%ni,1:Block1%nj,1:Block1%nk))
            !allocate(EQN_update_R3(1:Block1%nvar,1:Block1%ni,1:Block1%nj,1:Block1%nk))
            !allocate(EQN_update_R4(1:Block1%nvar,1:Block1%ni,1:Block1%nj,1:Block1%nk))

            ! loop from runge-kutta index of 1 to 4
            !
            do rk_index=1,4
                ! loop all the blocks
                !
                do iLocalBlock=1,size(Tree%LocalBlocks)

                    ! assign block pointer and see if it's at top or bottom.
                    !
                    Block1=>Tree%LocalBlocks(iLocalBlock)
                    call ModBundary_if_top_bottom(Block1,Tree%xijk_range,if_top,if_bottom)

                    ! first assign the pointers to simplify the code
                    !
                    select case(rk_index)
                    case(1)
                        primitive_i=>Block1%primitive   ; primitive_i_1=>Block1%primitive_k2
                        !EQN_update_Ri=>EQN_update_R1
                    case(2)
                        primitive_i=>Block1%primitive_k2; primitive_i_1=>Block1%primitive_k3
                        !EQN_update_Ri=>EQN_update_R2
                    case(3)
                        primitive_i=>Block1%primitive_k3; primitive_i_1=>Block1%primitive_k4
                        !EQN_update_Ri=>EQN_update_R3
                    case(4)
                        primitive_i=>Block1%primitive_k4; primitive_i_1=>Block1%primitive
                        !EQN_update_Ri=>EQN_update_R4
                    end select
                    
                    ! Get the EQN_update_Ri
                    !
                    call ModEquation_Dynamo_HD(primitive_i,Block1%ni,Block1%nj,Block1%nk,Block1%ng,&
                        Block1%dxi,Block1%dxj,Block1%dxk,Block1%xk,'cartesian',EQN_update_R1,if_top,if_bottom)
                    
                    
                    ! Get the next RK
                    !
                    primitive_i_1(:,1:Block1%ni,1:Block1%nj,1:Block1%nk)=&
                        Block1%primitive(:,1:Block1%ni,1:Block1%nj,1:Block1%nk)+dt*EQN_update_R1/(5.0-rk_index)

                    do i=1,Block1%ni
                        do j=1,Block1%nj
                            do k=1,Block1%nk
                                if (ieee_is_nan(EQN_update_R1(2,i,j,k)) .or. &
                                    ieee_is_nan(EQN_update_R1(3,i,j,k)) .or. &
                                    ieee_is_nan(EQN_update_R1(4,i,j,k))) then
                                    
                                    print *,'start debugging'
                                    print *,EQN_update_R1(:,i,j,k)
                                    stop 1
                                end if
                            end do
                        end do
                    end do
                    
                end do

                if (rk_index<4) call ModAdvance_CommunicateAll(Tree,rk_index+1,MpiSize,MpiRank)

            end do
        else
            

            ! loop from runge-kutta index of 1 to 4
            !
            do rk_index=1,4

                ! loop all the blocks
                !
                do iLocalBlock=1,size(Tree%LocalBlocks)

                    ! allocate the Runge-kutta R1 to R4
                    !
                    allocate(EQN_update_R1(1:Block1%nvar,1:Block1%ni,1:Block1%nj,1:Block1%nk))
                    allocate(EQN_update_R2(1:Block1%nvar,1:Block1%ni,1:Block1%nj,1:Block1%nk))
                    allocate(EQN_update_R3(1:Block1%nvar,1:Block1%ni,1:Block1%nj,1:Block1%nk))
                    allocate(EQN_update_R4(1:Block1%nvar,1:Block1%ni,1:Block1%nj,1:Block1%nk))

                    ! assign the pointers to simplify the code
                    !
                    select case(rk_index)
                    case(1)
                        primitive_i=>Block1%primitive   ; primitive_i_1=>Block1%primitive_k2
                        EQN_update_Ri=>EQN_update_R1
                    case(2)
                        primitive_i=>Block1%primitive_k2; primitive_i_1=>Block1%primitive_k3
                        EQN_update_Ri=>EQN_update_R2
                    case(3)
                        primitive_i=>Block1%primitive_k3; primitive_i_1=>Block1%primitive_k4
                        EQN_update_Ri=>EQN_update_R3
                    case(4)
                        primitive_i=>Block1%primitive_k4; primitive_i_1=>Block1%primitive
                        EQN_update_Ri=>EQN_update_R4
                    end select

                    ! assign block pointer and see if it's at top or bottom.
                    !
                    Block1=>Tree%LocalBlocks(iLocalBlock)
                    call ModBundary_if_top_bottom(Block1,Tree%xijk_range,if_top,if_bottom)
                    call ModEquation_Dynamo_HD(primitive_i,Block1%ni,Block1%nj,Block1%nk,Block1%ng,&
                        Block1%dxi,Block1%dxj,Block1%dxk,Block1%xk,'cartesian',EQN_update_R1,if_top,if_bottom)
                    
                    primitive_i_1(:,1:Block1%ni,1:Block1%nj,1:Block1%nk)=&
                        primitive_i(:,1:Block1%ni,1:Block1%nj,1:Block1%nk)+dt*EQN_update_Ri/(5.0-rk_index)
                end do

                if (rk_index<4) call ModAdvance_CommunicateAll(Tree,rk_index+1,MpiSize,MpiRank)
            end do
        end if

        call ModAdvance_CommunicateAll(Tree,1,MpiSize,MpiRank)

    end subroutine ModAdvance_rk4

    subroutine ModAdvance_CommunicateAll(Tree,rk_index,MpiSize,MpiRank)
        implicit none

        type(OcTree),intent(inout),target   ::  Tree            ! Tree
        integer,intent(in)                  ::  rk_index        ! rk index
        integer,intent(in)                  ::  MpiSize,MpiRank ! for mpi
        
        call ModCommunication_SendRecvGCAll(Tree,rk_index,MpiSize,MpiRank)
        !call ModCommunication_RecvGCAll(Tree,rk_index,MpiSize,MpiRank)
        
        
        call ModGC_CommunicateGCLocal(Tree,MpiRank,rk_index)
    end subroutine ModAdvance_CommunicateAll

end module ModAdvance