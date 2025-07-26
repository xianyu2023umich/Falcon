module ModPFSS
    ! This module setups the PFSS linear problem and calls the linear solver.
    ! The production of magnetogram is done in ModMagnetogram.f90. Here we only
    ! do the PFSS.

    use ModMath,                   only: ModMath_rebin_2D_factor2, ModMath_prolongate_3D_factor2
    use ModSpherical,              only: ModSpherical_Grad_f_O2
    use ModYinYang,                only: ModYinYang_CoordConv_0D
    use ModYinYangTree,            only: YYTree
    use ModBlock,                  only: BlockType,Multigrid_level
    use ModParameters,             only: ni,nj,nk,ng,nvar,MpiRank,MpiSize,Multigrid_nLevels
    use ModLinearSolver,           only: ModLinearSolver_Multigrid_CG, ModLinearSolver_JCB_CG
    !use ModCommunication,          only: ModCommunication_SendRecvHC

    implicit none

    contains

    ! The setup sets the b and A (Aii, Ai, Aj, Ak) for the PFSS problem.

    subroutine ModPFSS_setup(Tree)
        implicit none
        type(YYTree),target             ::  Tree

        type(Multigrid_level),pointer   ::  MGL1
        type(BlockType),pointer         ::  Block1
        integer                         ::  iBlock,iLevel
        integer                         ::  j,k

        ! First setup b
        ! b should be Br (magnetogram) times the area of the face.
        ! Here we substitute the magnetogram for all the
        ! MGL levels using the original map and rebin.

        do iBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iBlock)

            do iLevel=1,Multigrid_nLevels
                MGL1=>Block1%Multigrid_levels(iLevel)
                MGL1%b_LLL=0.0

                if (Block1%if_bottom .and. allocated(Block1%magnetogram_LL)) then                
                    if (iLevel==1) then
                        MGL1%b_LLL(1,:,:)=Block1%magnetogram_LL*MGL1%Si_FLL(1,:,:)
                    else
                        MGL1%b_LLL(1,:,:)=MGL1%Si_FLL(1,:,:)*&
                            ModMath_rebin_2D_factor2(MGL1%MGL_finer%b_LLL(1,:,:)/MGL1%MGL_finer%Si_FLL(1,:,:),&
                            MGL1%MGL_finer%mj,MGL1%MGL_finer%mk)
                    end if
                end if
            end do
            
        end do

        ! For all the other points the b is always zero.

        ! Matrix setup
        do iLevel=1,Multigrid_nLevels
            do iBlock=1,Tree%nLocalBlocks
                ! Get the matrix
                Block1=>Tree%LocalBlocks(iBlock)
                MGL1=>Block1%Multigrid_levels(iLevel)

                ! Deal with the boundary conditions
                ! For the bottom boundary, the bottom face is considered with b_LLL(1,:,:).
                ! So we need to set the non-diagonal A_i_FLL(1,:,:) to be zero. Other things
                ! still the same.

                MGL1%Ai_FLL=-MGL1%Si_FLL/MGL1%Di_FLL
                MGL1%Aj_LFL=-MGL1%Sj_LFL/MGL1%Dj_LFL
                MGL1%Ak_LLF=-MGL1%Sk_LLF/MGL1%Dk_LLF

                if (Block1%if_bottom) MGL1%Ai_FLL(1,:,:)=0.0

                MGL1%Aii_LLL=-( MGL1%Ai_FLL(1:MGL1%mi,1:MGL1%mj,1:MGL1%mk)+&
                        MGL1%Ai_FLL(2:MGL1%mi+1,1:MGL1%mj,1:MGL1%mk)+&
                        MGL1%Aj_LFL(1:MGL1%mi,1:MGL1%mj,1:MGL1%mk)+&
                        MGL1%Aj_LFL(1:MGL1%mi,2:MGL1%mj+1,1:MGL1%mk)+&
                        MGL1%Ak_LLF(1:MGL1%mi,1:MGL1%mj,1:MGL1%mk)+&
                        MGL1%Ak_LLF(1:MGL1%mi,1:MGL1%mj,2:MGL1%mk+1))
                
                if (Block1%if_top) then
                    MGL1%Aii_LLL(MGL1%mi,:,:)=MGL1%Aii_LLL(MGL1%mi,:,:)-MGL1%Ai_FLL(MGL1%mi+1,:,:)
                    MGL1%Ai_FLL(MGL1%mi+1,:,:)=0.0
                end if
            end do
        end do

        ! Initial guess
        ! test only
        !do iBlock=1,Tree%nLocalBlocks
        !    Block1=>Tree%LocalBlocks(iBlock)
        !    MGL1=>Block1%Multigrid_levels(1)
        !    do j=1,MGL1%mj
        !        do k=1,MGL1%mk
        !            MGL1%x_III(:,j,k)=[0.99,0.2,0.6]!0.333798230   ,   0.154550985     ,  3.98327261E-02]
        !        end do
        !    end do
        !end do

    end subroutine ModPFSS_setup

    subroutine ModPFSS_solve(Tree,max_iter,tol,i_option)
        implicit none
        type(YYTree),target             ::  Tree
        integer,intent(in)              ::  max_iter
        real,intent(in)                 ::  tol
        integer,intent(in)              ::  i_option

        type(BlockType),pointer         ::  Block1
        integer                         ::  iBlock,iLevel

        ! Solve using Multigrid preconditioned conjugate gradient.

        select case (i_option)
        case (1)
            call ModLinearSolver_Multigrid_CG(Tree,1,max_iter,tol)
        case (2)
            call ModLinearSolver_JCB_CG(Tree,1,max_iter,tol)
        case (3)
            call ModLinearSolver_JCB_CG(Tree,2,max_iter,tol)
            do iBlock=1,Tree%nLocalBlocks
                Block1=>Tree%LocalBlocks(iBlock)
                
                Block1%Multigrid_levels(1)%x_III=ModMath_prolongate_3D_factor2(&
                    Block1%Multigrid_levels(2)%x_III(1:Block1%Multigrid_levels(2)%mi,&
                    1:Block1%Multigrid_levels(2)%mj,1:Block1%Multigrid_levels(2)%mk),&
                    Block1%Multigrid_levels(1)%mi,Block1%Multigrid_levels(1)%mj,Block1%Multigrid_levels(1)%mk)
            end do
        end select
        
        ! Now we get Phi as the solution (x). We need to convert it to B.

        do iBlock=1,Tree%nLocalBlocks

            ! First get Phi
            Block1=>Tree%LocalBlocks(iBlock)
            Block1%Phi_LLL=Block1%Multigrid_levels(1)%x_III

            ! Then use derivatives to get B
            ! B = -grad(Phi)

            Block1%B0_IV=ModSpherical_Grad_f_O2(ni,nj,nk,1,&
                Block1%xi_I,Block1%xj_I,Block1%dxi,Block1%dxj,Block1%dxk,&
                Block1%Phi_LLL)
        end do

        call ModPFSS_finalize(Tree)
    end subroutine ModPFSS_solve

    ! Set the arrays to zero.

    subroutine ModPFSS_finalize(Tree)
        implicit none
        type(YYTree),target             ::  Tree

        type(Multigrid_level),pointer   ::  MGL1
        type(BlockType),pointer         ::  Block1
        integer                         ::  iBlock,iLevel

        do iBlock=1,Tree%nLocalBlocks
            do iLevel=1,Multigrid_nLevels
                MGL1=>Block1%Multigrid_levels(iLevel)

                ! Set everything to zero
                MGL1%Aii_LLL=0.0
                MGL1%Ai_FLL=0.0
                MGL1%Aj_LFL=0.0
                MGL1%Ak_LLF=0.0
                MGL1%x_III=0.0
                MGL1%b_LLL=0.0

                ! Set krylov arrays to zero
                MGL1%p_III=0.0
                MGL1%r_LLL=0.0
                MGL1%z_III=0.0
                MGL1%res_LLL=0.0
                MGL1%Ap_LLL=0.0
                MGL1%Az_LLL=0.0
            end do
        end do
    end subroutine ModPFSS_finalize

end module ModPFSS