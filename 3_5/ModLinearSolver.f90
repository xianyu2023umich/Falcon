module ModLinearSolver

    use ModMath,                   only: ModMath_rebin_3D_factor2,ModMath_prolongate_3D_factor2
    use ModYinYangTree,            only: YYTree
    use ModBlock,                  only: BlockType,Multigrid_level,GC_target
    use ModParameters,             only: ni,nj,nk,ng,nvar,MpiRank,MpiSize,Multigrid_nLevels
    use ModMultigrid,              only: ModMultigrid_Jacobi_Ax_b,ModMultigrid_Jacobi_Mz_r
    use ModCommunication,          only: ModCommunication_SendRecvHC
    use MPI

    implicit none

    contains

    subroutine ModLinearSolver_A_mul_z(MGL1)
        implicit none
        type(Multigrid_level)           ::  MGL1
        integer                         ::  mi,mj,mk

        mi=MGL1%mi
        mj=MGL1%mj
        mk=MGL1%mk

        MGL1%Az_LLL=&
                MGL1%Aii_LLL            *MGL1%z_III(1:mi,1:mj,1:mk)        &
            +   MGL1%Ai_FLL(1:mi,:,:)   *MGL1%z_III(0:mi-1,1:mj,1:mk)      &
            +   MGL1%Ai_FLL(2:mi+1,:,:) *MGL1%z_III(2:mi+1,1:mj,1:mk)      &
            +   MGL1%Aj_LFL(:,1:mj,:)   *MGL1%z_III(1:mi,0:mj-1,1:mk)      &
            +   MGL1%Aj_LFL(:,2:mj+1,:) *MGL1%z_III(1:mi,2:mj+1,1:mk)      &
            +   MGL1%Ak_LLF(:,:,1:mk)   *MGL1%z_III(1:mi,1:mj,0:mk-1)      &
            +   MGL1%Ak_LLF(:,:,2:mk+1) *MGL1%z_III(1:mi,1:mj,2:mk+1)
    end subroutine ModLinearSolver_A_mul_z

    subroutine ModLinearSolver_A_mul_p(MGL1)
        implicit none
        type(Multigrid_level)           ::  MGL1
        integer                         ::  mi,mj,mk

        mi=MGL1%mi
        mj=MGL1%mj
        mk=MGL1%mk

        MGL1%Ap_LLL=&
                MGL1%Aii_LLL            *MGL1%p_III(1:mi,1:mj,1:mk)        &
            +   MGL1%Ai_FLL(1:mi,:,:)   *MGL1%p_III(0:mi-1,1:mj,1:mk)      &
            +   MGL1%Ai_FLL(2:mi+1,:,:) *MGL1%p_III(2:mi+1,1:mj,1:mk)      &
            +   MGL1%Aj_LFL(:,1:mj,:)   *MGL1%p_III(1:mi,0:mj-1,1:mk)      &
            +   MGL1%Aj_LFL(:,2:mj+1,:) *MGL1%p_III(1:mi,2:mj+1,1:mk)      &
            +   MGL1%Ak_LLF(:,:,1:mk)   *MGL1%p_III(1:mi,1:mj,0:mk-1)      &
            +   MGL1%Ak_LLF(:,:,2:mk+1) *MGL1%p_III(1:mi,1:mj,2:mk+1)
    end subroutine ModLinearSolver_A_mul_p

    subroutine ModLinearSolver_Get_r0_All(Tree,iLevel)
        implicit none
        type(YYTree),target             ::  Tree
        integer,intent(in)              ::  iLevel

        type(Multigrid_level),pointer   ::  MGL1
        integer                         ::  iBlock

        do iBlock=1,Tree%nLocalBlocks
            MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
            call ModLinearSolver_Get_r0(MGL1)
        end do
    end subroutine ModLinearSolver_Get_r0_All

    subroutine ModLinearSolver_Get_r0(MGL1)
        implicit none
        type(Multigrid_level)           ::  MGL1
        integer                         ::  mi,mj,mk

        mi=MGL1%mi
        mj=MGL1%mj
        mk=MGL1%mk

        MGL1%r_LLL=MGL1%b_LLL&
            -   MGL1%Aii_LLL           *MGL1%x_III(1:mi,1:mj,1:mk)      & 
            -   MGL1%Ai_FLL(1:mi,:,:)  *MGL1%x_III(0:mi-1,1:mj,1:mk)    &
            -   MGL1%Ai_FLL(2:mi+1,:,:)*MGL1%x_III(2:mi+1,1:mj,1:mk)    &
            -   MGL1%Aj_LFL(:,1:mj,:)  *MGL1%x_III(1:mi,0:mj-1,1:mk)    &
            -   MGL1%Aj_LFL(:,2:mj+1,:)*MGL1%x_III(1:mi,2:mj+1,1:mk)    &
            -   MGL1%Ak_LLF(:,:,1:mk)  *MGL1%x_III(1:mi,1:mj,0:mk-1)    &
            -   MGL1%Ak_LLF(:,:,2:mk+1)*MGL1%x_III(1:mi,1:mj,2:mk+1)
    end subroutine ModLinearSolver_Get_r0

    subroutine ModLinearSolver_JCB_CG(Tree,iLevel,max_iter,tol)
        implicit none
        type(YYTree),target             ::  Tree
        integer,intent(in)              ::  iLevel
        integer,intent(in)              ::  max_iter
        real,intent(in)                 ::  tol

        type(Multigrid_level),pointer   ::  MGL1
        integer                         ::  iBlock
        real                            ::  pAp_local,rz_sum_local,rz_next_sum_local,alpha,beta
        real                            ::  pAp_global,rz_sum_global,rz_next_sum_global
        integer                         ::  iter,ierr
        integer                         ::  mi,mj,mk

        mi=Tree%LocalBlocks(1)%Multigrid_levels(iLevel)%mi
        mj=Tree%LocalBlocks(1)%Multigrid_levels(iLevel)%mj
        mk=Tree%LocalBlocks(1)%Multigrid_levels(iLevel)%mk

        ! Communicate x
        call ModCommunication_SendRecvHC(Tree,iLevel,i_option=2)

        ! Get the r0
        call ModLinearSolver_Get_r0_All(Tree,iLevel)
        ! Get the z0
        rz_sum_local=0.0
        do iBlock=1,Tree%nLocalBlocks
            MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
            MGL1%z_III(1:mi,1:mj,1:mk)=MGL1%r_LLL*MGL1%V_LLL
            rz_sum_local=rz_sum_local+sum(MGL1%r_LLL*MGL1%z_III(1:mi,1:mj,1:mk))
        end do

        ! Get the global sum
        call MPI_Allreduce(rz_sum_local,rz_sum_global,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

        ! Get the p0
        do iBlock=1,Tree%nLocalBlocks
            MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
            MGL1%p_III=MGL1%z_III
        end do
        ! Communicate the p0

        call ModCommunication_SendRecvHC(Tree,iLevel,i_option=1)
        do iter=1,max_iter
            ! Get the Ap
            pAp_local=0.0
            do iBlock=1,Tree%nLocalBlocks
                MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
                call ModLinearSolver_A_mul_p(MGL1)
                pAp_local=pAp_local+sum(MGL1%p_III(1:mi,1:mj,1:mk)*MGL1%Ap_LLL)
            end do

            ! Get the global sum of pAp
            call MPI_Allreduce(pAp_local,pAp_global,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ! Get the alpha
            alpha=rz_sum_global/pAp_global
            
            ! Update x
            ! Update the r and z
            ! No need to communicate the r and z, because the r and z do not contain the ghost cells.
            do iBlock=1,Tree%nLocalBlocks
                MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
                MGL1%x_III(1:mi,1:mj,1:mk)=&
                    MGL1%x_III(1:mi,1:mj,1:mk)+alpha*MGL1%p_III(1:mi,1:mj,1:mk)
                MGL1%r_LLL=MGL1%r_LLL-alpha*MGL1%Ap_LLL
                MGL1%z_III(1:mi,1:mj,1:mk)=MGL1%r_LLL*MGL1%V_LLL
            end do

            ! Get the global sum of rz
            rz_next_sum_local=0.0
            do iBlock=1,Tree%nLocalBlocks
                MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
                rz_next_sum_local=rz_next_sum_local+sum(MGL1%r_LLL*MGL1%z_III(1:mi,1:mj,1:mk))
            end do
            call MPI_Allreduce(rz_next_sum_local,rz_next_sum_global,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

            ! Exit if accurate enough
            if(rz_next_sum_global<tol) exit
            
            ! Get the beta
            beta=rz_next_sum_global/rz_sum_global
            do iBlock=1,Tree%nLocalBlocks   
                MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
                MGL1%p_III(1:mi,1:mj,1:mk)=MGL1%z_III(1:mi,1:mj,1:mk)+beta*MGL1%p_III(1:mi,1:mj,1:mk)
            end do

            ! Communicate the p
            call ModCommunication_SendRecvHC(Tree,iLevel,i_option=1)

            rz_sum_global=rz_next_sum_global

            print *,iter,rz_sum_global
        end do
        !if (MpiRank==0) print *,iter
    end subroutine ModLinearSolver_JCB_CG

    ! Multigrid V-cycle preconditioner
    ! What a multigrid V cycle does is to solve Mz=r, where M is a approximate of A.
    ! Consider two levels h (finer )and 2h (coarser), then the multigrid does following steps:
    !
    ! 1. Initial guess (e.g. z_h=0)
    ! 2. Use GS or jacobi smoothing to get z_h
    ! 3. Compute the residual r^res_h=r_h-A_h*z_h
    ! 4. Restrict r^res_h to r^res_2h: r^res_2h=R * r^res_h
    ! 5. Solve M_2h (or A_2h if possible) z_2h=r^res_2h
    ! 6. Prolongate z_2h to P z_h:
    ! 7. Update z_h+=P z_2h
    !
    ! In total we use z_h = z_h + P  A_2h^{-1} R (r_h - A_h z_h) to get an estimation of
    ! the solution of M_h z_h = r_h. Obivously, in this process we need A_2h^{-1} -- which
    ! can be also estimated by Multigrid V cycle. Therefore this is a recursive process.
    !
    ! If the level is the coarsest level, then we use CG to solve the linear system.
    ! In this particular version we use Jacobi smoothing twice before and once after.

    recursive subroutine ModLinearSolver_Multigrid_V_cycle(Tree,iLevel,max_iter,tol)
        implicit none
        type(YYTree),target             ::  Tree
        integer,intent(in)              ::  iLevel
        integer,intent(in)              ::  max_iter
        real,intent(in)                 ::  tol

        type(Multigrid_level),pointer   ::  MGL1,MGL_next
        integer                         ::  iBlock
        integer                         ::  iSmoothing

        print *,111,iLevel

        ! Initialize the z_III to 0.0
        do iBlock=1,Tree%nLocalBlocks
            MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
            MGL1%z_III=0.0
        end do

        print *,222,iLevel

        ! First use Jacobi smoothing twice.
        ! Since Jacobi takes the neighbor cells, we need mpi communication.

        do iSmoothing=1,3
            call ModCommunication_SendRecvHC(Tree,iLevel,i_option=3)
            do iBlock=1,Tree%nLocalBlocks ! First smoothing
                MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
                call ModMultigrid_Jacobi_Mz_r(MGL1)
            end do
        end do

        print *,333,iLevel

        ! Compute the residual
        do iBlock=1,Tree%nLocalBlocks
            MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
            call ModLinearSolver_A_mul_z(MGL1)
            MGL1%res_LLL=MGL1%r_LLL-MGL1%Az_LLL
        end do

        print *,444,iLevel

        ! Two options again: if next level is coarsest then we will pass
        ! residual to b_2h in A_2h x_2h=b_2h and solve for x_2h. Otherwise we will pass it
        ! to r_2h in A_2h z_2h=r_2h and solve for z_2h.

        if (Multigrid_nLevels/=1) then

            print *,555,iLevel
            if (iLevel+1==Multigrid_nLevels) then
                ! Pass it to b
                do iBlock=1,Tree%nLocalBlocks
                    MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
                    MGL_next=>MGL1%MGL_coarser
                    MGL_next%b_LLL=ModMath_rebin_3D_factor2(MGL1%res_LLL,MGL1%mi,MGL1%mj,MGL1%mk)
                end do
    
                ! call CG solver
    
                call ModLinearSolver_JCB_CG(Tree,iLevel+1,max_iter,tol)

    
                ! Prolongate the solution to x_h
                do iBlock=1,Tree%nLocalBlocks
                    MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
                    MGL_next=>MGL1%MGL_coarser
                    MGL1%x_III=ModMath_prolongate_3D_factor2(MGL_next%x_III,MGL_next%mi,MGL_next%mj,MGL_next%mk)
                end do

            else

                print *,666,iLevel
                ! Restrct the res_h to r_2h
                do iBlock=1,Tree%nLocalBlocks
                    MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
                    MGL_next=>MGL1%MGL_coarser
                    MGL_next%r_LLL=ModMath_rebin_3D_factor2(MGL1%res_LLL,MGL1%mi,MGL1%mj,MGL1%mk)
                end do
    
                ! call V cycle
                call ModLinearSolver_Multigrid_V_cycle(Tree,iLevel+1,max_iter,tol)

                print *,12345,iLevel
    
                ! Prolongate the z_2h solution to P z_2h, and add it to z_h
                do iBlock=1,Tree%nLocalBlocks
                    MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
                    MGL_next=>MGL1%MGL_coarser
                    MGL1%z_III(1:ni,1:nj,1:nk)=MGL1%z_III(1:ni,1:nj,1:nk)+&
                        ModMath_prolongate_3D_factor2(MGL_next%z_III(1:ni,1:nj,1:nk),MGL_next%mi,MGL_next%mj,MGL_next%mk)
                end do

                print *,123456,iLevel
            end if
        end if

        print *,777,iLevel

        ! No matther which one we need one more smoothing
        !call ModCommunication_SendRecvHC(Tree,iLevel,i_option=3)
        do iSmoothing=1,3
            call ModCommunication_SendRecvHC(Tree,iLevel,i_option=3)
            do iBlock=1,Tree%nLocalBlocks ! First smoothing
                MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
                call ModMultigrid_Jacobi_Mz_r(MGL1)
            end do
        end do

        print *,888,iLevel
    end subroutine ModLinearSolver_Multigrid_V_cycle

    subroutine ModLinearSolver_Multigrid_CG(Tree,iLevel,max_iter,tol)
        implicit none
        type(YYTree),target             ::  Tree
        integer,intent(in)              ::  iLevel
        integer,intent(in)              ::  max_iter
        real,intent(in)                 ::  tol

        type(Multigrid_level),pointer   ::  MGL1
        integer                         ::  iBlock
        real                            ::  pAp_local,rz_sum_local,rz_next_sum_local,alpha,beta
        real                            ::  pAp_global,rz_sum_global,rz_next_sum_global
        integer                         ::  iter,ierr


        ! Communicate x
        call ModCommunication_SendRecvHC(Tree,iLevel,i_option=2)

        ! Get the r0
        call ModLinearSolver_Get_r0_All(Tree,iLevel)

        ! Get the z0
        rz_sum_local=0.0

        call ModLinearSolver_Multigrid_V_cycle(Tree,iLevel,max_iter,tol*10.0)

        ! The following is turning off preconditioner test-- M=I so that z=r.
        !
        !do iBlock=1,Tree%nLocalBlocks
        !    MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
        !    MGL1%z_III(1:ni,1:nj,1:nk)=MGL1%r_LLL
        !end do

        ! Get the r0
        call ModLinearSolver_Get_r0_All(Tree,iLevel)

        do iBlock=1,Tree%nLocalBlocks
            MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
            rz_sum_local=rz_sum_local+sum(MGL1%r_LLL*MGL1%z_III(1:ni,1:nj,1:nk))
        end do

        ! Get the global sum
        call MPI_Allreduce(rz_sum_local,rz_sum_global,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

        print *,rz_sum_global,1111

        ! Get the p0
        do iBlock=1,Tree%nLocalBlocks
            MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
            MGL1%p_III=MGL1%z_III
        end do
        ! Communicate the p0
        call ModCommunication_SendRecvHC(Tree,iLevel,i_option=1)

        print *,2222
        do iter=1,max_iter
            
            ! Get the Ap
            pAp_local=0.0
            do iBlock=1,Tree%nLocalBlocks
                MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
                call ModLinearSolver_A_mul_p(MGL1)
                pAp_local=pAp_local+sum(MGL1%p_III(1:ni,1:nj,1:nk)*MGL1%Ap_LLL)
            end do

            ! Get the global sum of pAp
            call MPI_Allreduce(pAp_local,pAp_global,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ! Get the alpha
            alpha=rz_sum_global/pAp_global
            
            ! Update x
            ! Update the r and z
            ! No need to communicate the r and z, because the r and z do not contain the ghost cells.
            do iBlock=1,Tree%nLocalBlocks
                MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
                MGL1%x_III(1:ni,1:nj,1:nk)=&
                    MGL1%x_III(1:ni,1:nj,1:nk)+alpha*MGL1%p_III(1:ni,1:nj,1:nk)
                MGL1%r_LLL=MGL1%r_LLL-alpha*MGL1%Ap_LLL
            end do

            call ModLinearSolver_Multigrid_V_cycle(Tree,iLevel,max_iter,tol)

            ! The following is turning off preconditioner test-- M=I so that z=r.
            !
            !do iBlock=1,Tree%nLocalBlocks
            !    MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
            !    MGL1%z_III(1:ni,1:nj,1:nk)=MGL1%r_LLL
            !end do

            ! Get the global sum of rz
            rz_next_sum_local=0.0
            do iBlock=1,Tree%nLocalBlocks
                MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
                rz_next_sum_local=rz_next_sum_local+sum(MGL1%r_LLL*MGL1%z_III(1:ni,1:nj,1:nk))
            end do
            call MPI_Allreduce(rz_next_sum_local,rz_next_sum_global,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

            ! Exit if accurate enough
            if(rz_next_sum_global<tol) exit
            
            ! Get the beta
            beta=rz_next_sum_global/rz_sum_global
            do iBlock=1,Tree%nLocalBlocks   
                MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
                MGL1%p_III(1:ni,1:nj,1:nk)=MGL1%z_III(1:ni,1:nj,1:nk)+beta*MGL1%p_III(1:ni,1:nj,1:nk)
            end do

            ! Communicate the p
            call ModCommunication_SendRecvHC(Tree,iLevel,i_option=1)
            rz_sum_global=rz_next_sum_global
            print *,rz_sum_global
        end do
        if (MpiRank==0) print *,iter
        print *,1111
    end subroutine ModLinearSolver_Multigrid_CG

end module ModLinearSolver