module ModAdvance

    use ieee_arithmetic
    use ModEquation
    use ModGC
    use ModBoundary
    use ModTimeStep
    use ModCommunication
    use ModParameters,  only:   ni,nj,nk,nvar

    contains

    subroutine ModAdvance_rk4(Tree,CFL_ad,if_all_same_sizes,MpiSize,MpiRank,dt)
        implicit none
        type(YYTree),intent(inout),target   ::  Tree                        ! Tree
        real,intent(in)                     ::  CFL_ad                      ! CFL
        logical,intent(in)                  ::  if_all_same_sizes           ! if all local blocks same size
        integer,intent(in)                  ::  MpiSize,MpiRank

        type(Block),pointer                 ::  Block1                      ! block pointer
        integer                             ::  iLocalBlock                 ! i of block
        integer                             ::  rk_index                    ! runge-kutta index
        logical                             ::  if_rk_input,if_rk_output    ! If use rk for R i/o
        real,allocatable                    ::  EQN_update_R(:,:,:,:)       ! the following four used for rk

        real                                ::  dt_local,dt                 ! dt
        !integer :: i,j,k

        call ModTimeStep_Dynamo_HD(Tree,CFL_ad,dt_local)
        call ModCommunication_GlobalTimeStep(dt_local,dt)
          
        ! see if all the local blocks in the tree
        ! have the same block size. If yes then it
        ! there is no need to allocate and deallocate
        ! every time.

        if (if_all_same_sizes) then

            ! allocate the Runge-kutta R1 to R4 according to the Block1 size.
            Block1 => Tree%LocalBlocks(1)
            allocate(EQN_update_R(1:ni,1:nj,1:nk,1:nvar))

            ! loop from runge-kutta index of 1 to 4
            ! loop all the blocks
            do rk_index=1,4
                do iLocalBlock=1,size(Tree%LocalBlocks)
                    ! assign block pointer and see if it's at top or bottom.
                    Block1=>Tree%LocalBlocks(iLocalBlock)
                    !if (Block1%iBlock==6) print *,Block1%primitive_rk(32,-1:3,1,3)
                    
                    ! Decide if rk or ot for the input and output primitives
                    if_rk_input=(rk_index>1)
                    if_rk_output=(rk_index<4)

                    ! Get the EQN_update_R
                    call ModEquation_Dynamo_HD(Block1,if_rk_input,EQN_update_R)

                    ! Get the next RK
                    if (if_rk_output) then
                        Block1%primitive_rk(1:ni,1:nj,1:nk,:)=&
                            Block1%primitive(1:ni,1:nj,1:nk,:)+dt*EQN_update_R/(5.0-rk_index)
                    else
                        Block1%primitive(1:ni,1:nj,1:nk,:)=&
                            Block1%primitive(1:ni,1:nj,1:nk,:)+dt*EQN_update_R/(5.0-rk_index)
                    end if
                end do
                !if(MpiRank==2) print *,dt*EQN_update_R(32,62:64,1,3)/(5.0-rk_index)
                !if(MpiRank==2) print *,Tree%LocalBlocks(2)%primitive(32,62:66,1,3)
                !if(MpiRank==2) print *,Tree%LocalBlocks(2)%primitive_rk(32,62:66,1,3)
                ! Communicate the GC for the current primitives
                call ModAdvance_CommunicateAll(Tree,if_rk_output,MpiSize,MpiRank)
                !call ModAdvance_CommunicateAll(Tree,if_rk_output,MpiSize,MpiRank)
                !if(MpiRank==2) print *,Tree%LocalBlocks(2)%primitive(32,62:66,1,3)
                !if(MpiRank==2) print *,Tree%LocalBlocks(2)%primitive_rk(32,62:66,1,3)
                !if(MpiRank==2) print *,1
                !if(MpiRank==2) print *,1
            end do
        else
        end if
    end subroutine ModAdvance_rk4

    subroutine ModAdvance_CommunicateAll(Tree,if_rk,MpiSize,MpiRank)
        implicit none

        type(YYTree),intent(inout),target   ::  Tree            ! Tree
        logical,intent(in)                  ::  if_rk           ! if_rk
        integer,intent(in)                  ::  MpiSize,MpiRank ! for mpi
        
        call ModCommunication_SendRecvGCAll(Tree,if_rk,MpiSize,MpiRank)
        call ModGC_CommunicateGCLocal(Tree,MpiRank,if_rk)
    end subroutine ModAdvance_CommunicateAll
end module ModAdvance