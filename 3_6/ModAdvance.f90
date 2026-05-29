module ModAdvance

    use ModBlock,           only:   BlockType
    use ModYinYangTree,     only:   YYTree
    use ModEquation,        only:   ModEquation_Dynamo
    use ModDiffusion,       only:   ModDiffusion_Aritificial_1
    use ModTimeStep,        only:   ModTimeStep_Dynamo
    use ModWaveSpeed,       only:   ModWaveSpeed_Dynamo
    use ModCommunication,   only:   ModCommunication_SendRecvGC,ModCommunication_SendRecvGC_new
    use ModParameters,      only:   ni,nj,nk,nvar,MpiRank,iEquation,CFL,DivB_option
    use ModCheck,           only:   ModCheck_primitive
    use ModBoundary,        only:   ModBoundary_Dynamo_primitives
    use ModDivB,            only:   ModDivB_GLM
    use MPI

    contains

    subroutine ModAdvance_rk4(Tree,if_all_same_sizes,dt_global)
        implicit none
        type(YYTree),intent(inout),target   ::  Tree                        ! Tree
        logical,intent(in)                  ::  if_all_same_sizes           ! if all local blocks same size

        type(BlockType),pointer             ::  Block1                      ! block pointer
        integer                             ::  iLocalBlock                 ! i of block
        integer                             ::  rk_index                    ! runge-kutta index
        logical                             ::  if_rk_input,if_rk_output    ! If use rk for R i/o
        real(8)                             ::  dt_local,dt_global          ! dt
        integer                             ::  ierr

        ! Compute wave speed for all local blocks before timestep and equation
        do iLocalBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iLocalBlock)
            select case(iEquation)
            case(0)
                call ModWaveSpeed_Dynamo(Block1)
            case(1)
                call ModWaveSpeed_Dynamo(Block1)
            end select
        end do

        select case(iEquation)
        case(0)
            call ModTimeStep_Dynamo(Tree,CFL,dt_local)
        case(1)
            call ModTimeStep_Dynamo(Tree,CFL,dt_local)
        end select

        call MPI_AllReduce(dt_local,dt_global,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD,ierr)
          
        ! see if all the local blocks in the tree
        ! have the same block size. If yes then it
        ! there is no need to allocate and deallocate
        ! every time.

        if (if_all_same_sizes) then

            ! loop from runge-kutta index of 1 to 4
            ! loop all the blocks
            do rk_index=1,4
                do iLocalBlock=1,size(Tree%LocalBlocks)
                    Block1=>Tree%LocalBlocks(iLocalBlock)
                    
                    ! Determine which is used as input or output.
                    if_rk_input=(rk_index>1)
                    if_rk_output=(rk_index<4)

                    ! Set primitive pointer for this RK stage
                    if (if_rk_input) then
                        Block1%primitive=>Block1%primitive_rk_IV
                    else
                        Block1%primitive=>Block1%primitive_IV
                    end if

                    select case(iEquation)
                    case(0)
                        call ModBoundary_Dynamo_primitives(Block1)
                        call ModEquation_Dynamo(Block1)
                    case(1)
                        call ModBoundary_Dynamo_primitives(Block1)
                        call ModEquation_Dynamo(Block1)

                        if (Block1%if_involve_B .and. DivB_option==1) then
                            call ModDivB_GLM(Block1)
                        end if
                    end select

                    if (if_rk_output) then
                        Block1%primitive_rk_IV(1:ni,1:nj,1:nk,:)=&
                            Block1%primitive_IV(1:ni,1:nj,1:nk,:)+dt_global*Block1%EQN_update_R_IV/(5.0-rk_index)
                    else
                        Block1%primitive_IV(1:ni,1:nj,1:nk,:)=&
                            Block1%primitive_IV(1:ni,1:nj,1:nk,:)+dt_global*Block1%EQN_update_R_IV/(5.0-rk_index)
                    end if
                end do

                call ModCommunication_SendRecvGC_new(Tree,if_rk_output)                
            end do

            ! Do the artificial diffusion after the rk loop

            do iLocalBlock=1,size(Tree%LocalBlocks)
               Block1=>Tree%LocalBlocks(iLocalBlock)

               Block1%EQN_update_R_IV=0.0
               call ModDiffusion_Aritificial_1(Block1,Block1%EQN_update_R_IV,2,.false.)
               Block1%primitive_IV(1:ni,1:nj,1:nk,:)=&
                   Block1%primitive_IV(1:ni,1:nj,1:nk,:)+0.5*dt_global*Block1%EQN_update_R_IV
            end do
            call ModCommunication_SendRecvGC_new(Tree,.false.)
        else
        end if
    end subroutine ModAdvance_rk4
end module ModAdvance