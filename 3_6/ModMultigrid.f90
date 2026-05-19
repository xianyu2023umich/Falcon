module ModMultigrid

    ! This is a relatively basic module so it shouldn't use too many other modules.
    ! It doesn't use ModGC. ModGC uses it. 
    ! It doesn't include multigrid methods. It will be done in some solver modules.

    use ModParameters,      only:   ni,nj,nk,Multigrid_nLevels,iGeometry
    use ModBlock,           only:   GC_target,Multigrid_level,BlockType
    use ModYinYangTree,     only:   YYTree

    implicit none

    contains

    ! Do everything below

    subroutine ModMultigrid_InitAll(Tree,nLevelmax,if_krylov)
        implicit none
        type(YYTree),target             ::  Tree
        integer,intent(in)              ::  nLevelmax
        logical,intent(in)              ::  if_krylov

        integer                         ::  iBlock
        type(BlockType),pointer         ::  Block1

        call ModMultigrid_get_nLevels(nLevelmax)

        do iBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iBlock)
            call ModMultigrid_Init(Block1,if_krylov)
        end do
    end subroutine ModMultigrid_InitAll

    subroutine ModMultigrid_Init(Block1,if_krylov)
        implicit none
        type(BlockType),target          ::  Block1
        logical,intent(in)              ::  if_krylov
        integer                         ::  iLevel
        type(Multigrid_level),pointer   ::  MGL1

        allocate(Block1%Multigrid_levels(Multigrid_nLevels))

        do iLevel=1,Multigrid_nLevels
            MGL1=>Block1%Multigrid_levels(iLevel)
            MGL1%iBlock=Block1%iBlock
            MGL1%iLevel=iLevel
            
            ! The finer level
            if(iLevel==1) then
                MGL1%MGL_finer=>null()
            else
                MGL1%MGL_finer=>Block1%Multigrid_levels(iLevel-1)
            end if

            ! The coarser level 
            if(iLevel==Multigrid_nLevels) then
                MGL1%MGL_coarser=>null()
            else
                MGL1%MGL_coarser=>Block1%Multigrid_levels(iLevel+1)
            end if

            ! Initialize the grid
            call ModMultigrid_InitGrid(MGL1,Block1%xijk_range,Block1%if_yin)

            ! Initialize the linear solver
            MGL1%if_krylov=if_krylov
            call ModMultigrid_InitLinearSolver(MGL1,if_krylov)
        end do
    end subroutine ModMultigrid_Init

    ! First I want a subroutine to get how many levels to be set.
    ! It should be determined by the factors of ni, nj and nk.
    ! For example, if nijk=8, then I can have levels in which mi,mj,mk are 1,2,4,8.
    ! Right now I only consider factor of 2. I don't want things like 3 or 5.
    
    subroutine ModMultigrid_get_nLevels(nLevelmax)
        implicit none
        integer,intent(in)              ::  nLevelmax
        integer                         ::  iLevel

        do iLevel=1,nLevelmax
            if (mod(ni,2**(iLevel-1))==0 .and. mod(nj,2**(iLevel-1))==0 .and. mod(nk,2**(iLevel-1))==0) then
                Multigrid_nLevels=iLevel
            else
                exit
            end if 
        end do
    end subroutine ModMultigrid_get_nLevels

    ! Now I want to set the levels.
    ! Allocate the multigrid_levels and allocate the arrays depending on the needs.
    ! mg is always 1, since multigrid is only used for parabolic problem and those
    ! problems only need neighbouring cells (there's even no need to set mg cuz it's
    ! never used.)
    ! HC sources and targets are done in ModGC.f90.

    subroutine ModMultigrid_InitGrid(MGL1,xijk_range,if_yin)
        implicit none
        type(Multigrid_level)           ::  MGL1
        real,intent(in)                 ::  xijk_range(3,2)
        logical,intent(in)              ::  if_yin
        integer                         ::  i,j,k
        integer                         ::  mi,mj,mk

        MGL1%if_yin=if_yin
        MGL1%mi=ni/(2**(MGL1%iLevel-1))
        MGL1%mj=nj/(2**(MGL1%iLevel-1))
        MGL1%mk=nk/(2**(MGL1%iLevel-1))
        mi=MGL1%mi
        mj=MGL1%mj
        mk=MGL1%mk

        ! Allocate HC_iBlock_III
        allocate(MGL1%HC_iBlock_III(0:mi+1,0:mj+1,0:mk+1))  
        MGL1%HC_iBlock_III=-777                             

        ! allocate the grid
        allocate(MGL1%xi_I(0:mi+1))
        allocate(MGL1%xj_I(0:mj+1))
        allocate(MGL1%xk_I(0:mk+1))
        allocate(MGL1%xi_F(1:mi+1))
        allocate(MGL1%xj_F(1:mj+1))
        allocate(MGL1%xk_F(1:mk+1))                                        

        ! Get the positions of the cell centers
        do i=0,mi+1; MGL1%xi_I(i)=(xijk_range(1,1)*(mi-i+0.5)+xijk_range(1,2)*(i-0.5))/mi; end do
        do j=0,mj+1; MGL1%xj_I(j)=(xijk_range(2,1)*(mj-j+0.5)+xijk_range(2,2)*(j-0.5))/mj; end do
        do k=0,mk+1; MGL1%xk_I(k)=(xijk_range(3,1)*(mk-k+0.5)+xijk_range(3,2)*(k-0.5))/mk; end do
        
        ! Get dxi, dxj, dxk
        MGL1%dxi=(xijk_range(1,2)-xijk_range(1,1))/mi
        MGL1%dxj=(xijk_range(2,2)-xijk_range(2,1))/mj
        MGL1%dxk=(xijk_range(3,2)-xijk_range(3,1))/mk

        ! Get the positions of the cell faces
        do i=1,mi+1; MGL1%xi_F(i)=(xijk_range(1,1)*(mi+1-i)+xijk_range(1,2)*(i-1))/mi; end do
        do j=1,mj+1; MGL1%xj_F(j)=(xijk_range(2,1)*(mj+1-j)+xijk_range(2,2)*(j-1))/mj; end do
        do k=1,mk+1; MGL1%xk_F(k)=(xijk_range(3,1)*(mk+1-k)+xijk_range(3,2)*(k-1))/mk; end do
        
        ! Set the face areas and cell volumes based on the geometry
        allocate(   MGL1%Si_FLL(1:mi+1,1:mj,1:mk),&
                    MGL1%Sj_LFL(1:mi,1:mj+1,1:mk),&
                    MGL1%Sk_LLF(1:mi,1:mj,1:mk+1))
        allocate(   MGL1%Di_FLL(1:mi+1,1:mj,1:mk),&
                    MGL1%Dj_LFL(1:mi,1:mj+1,1:mk),&
                    MGL1%Dk_LLF(1:mi,1:mj,1:mk+1))
        allocate(   MGL1%V_LLL(1:mi,1:mj,1:mk))

        select case(iGeometry)
        case(0)
            do k=1,mk
                do j=1,mj
                    do i=1,mi+1
                        MGL1%Si_FLL(i,j,k)=(MGL1%xj_F(j+1)-MGL1%xj_F(j))*(MGL1%xk_F(k+1)-MGL1%xk_F(k))
                    end do
                end do 
            end do
            do k=1,mk
                do j=1,mj+1
                    do i=1,mi
                        MGL1%Sj_LFL(i,j,k)=(MGL1%xi_F(i+1)-MGL1%xi_F(i))*(MGL1%xk_F(k+1)-MGL1%xk_F(k))
                    end do
                end do
            end do
            do k=1,mk+1
                do j=1,mj
                    do i=1,mi
                        MGL1%Sk_LLF(i,j,k)=(MGL1%xj_F(j+1)-MGL1%xj_F(j))*(MGL1%xi_F(i+1)-MGL1%xi_F(i))
                    end do
                end do
            end do
            do k=1,mk
                do j=1,mj
                    do i=1,mi
                        MGL1%V_LLL(i,j,k)=(MGL1%xi_F(i+1)-MGL1%xi_F(i))*&
                                    (MGL1%xj_F(j+1)-MGL1%xj_F(j))*(MGL1%xk_F(k+1)-MGL1%xk_F(k))
                    end do
                end do
            end do
            MGL1%Di_FLL=MGL1%dxi
            MGL1%Dj_LFL=MGL1%dxj
            MGL1%Dk_LLF=MGL1%dxk
        case(1)
            do k=1,mk
                do j=1,mj
                    do i=1,mi+1
                        MGL1%Si_FLL(i,j,k)=MGL1%xi_F(i)**2*&
                                    (cos(MGL1%xj_F(j))-cos(MGL1%xj_F(j+1)))*&
                                    (MGL1%xk_F(k+1)-MGL1%xk_F(k))
                    end do
                end do
            end do
            do k=1,mk
                do j=1,mj+1
                    do i=1,mi
                        MGL1%Sj_LFL(i,j,k)=0.5*(MGL1%xi_F(i+1)**2-MGL1%xi_F(i)**2)*(MGL1%xk_F(k+1)-MGL1%xk_F(k))*&
                                        sin(MGL1%xj_F(j))
                    end do
                end do
            end do
            do k=1,mk+1
                do j=1,mj
                    do i=1,mi
                        MGL1%Sk_LLF(i,j,k)=0.5*(MGL1%xi_F(i+1)**2-MGL1%xi_F(i)**2)*(MGL1%xj_F(j+1)-MGL1%xj_F(j))
                    end do
                end do
            end do
            do k=1,mk
                do j=1,mj
                    do i=1,mi
                        MGL1%V_LLL(i,j,k)=1.0/3.0*(MGL1%xi_F(i+1)**3-MGL1%xi_F(i)**3)*&
                                    (cos(MGL1%xj_F(j))-cos(MGL1%xj_F(j+1)))*&
                                    (MGL1%xk_F(k+1)-MGL1%xk_F(k))
                    end do
                end do
            end do
            MGL1%Di_FLL=MGL1%dxi
            do k=1,mk
                do j=1,mj+1
                    do i=1,mi
                        MGL1%Dj_LFL(i,j,k)=MGL1%xi_I(i)*MGL1%dxj
                    end do
                end do
            end do
            do k=1,mk+1
                do j=1,mj
                    do i=1,mi
                        MGL1%Dk_LLF(i,j,k)=MGL1%xi_I(i)*MGL1%dxk*sin(MGL1%xj_I(j))
                    end do
                end do
            end do
        end select
    end subroutine ModMultigrid_InitGrid

    subroutine ModMultigrid_InitLinearSolver(MGL1,if_krylov)
        implicit none
        type(Multigrid_level)           ::  MGL1
        logical,intent(in)              ::  if_krylov

        ! Allocate the arrays

        allocate(MGL1%Aii_LLL(1:MGL1%mi,1:MGL1%mj,1:MGL1%mk))
        allocate(MGL1%Ai_FLL(1:MGL1%mi+1,1:MGL1%mj,1:MGL1%mk))
        allocate(MGL1%Aj_LFL(1:MGL1%mi,1:MGL1%mj+1,1:MGL1%mk))
        allocate(MGL1%Ak_LLF(1:MGL1%mi,1:MGL1%mj,1:MGL1%mk+1))
        allocate(MGL1%x_III(0:MGL1%mi+1,0:MGL1%mj+1,0:MGL1%mk+1))
        allocate(MGL1%b_LLL(1:MGL1%mi,1:MGL1%mj,1:MGL1%mk))

        ! Initialize the arrays
        
        MGL1%Aii_LLL=0.0
        MGL1%Ai_FLL=0.0
        MGL1%Aj_LFL=0.0
        MGL1%Ak_LLF=0.0
        MGL1%x_III=0.0
        MGL1%b_LLL=0.0

        ! Initialize the linear solver

        if(if_krylov) then
            allocate(MGL1%p_III(0:MGL1%mi+1,0:MGL1%mj+1,0:MGL1%mk+1))
            allocate(MGL1%r_LLL(1:MGL1%mi,1:MGL1%mj,1:MGL1%mk))
            allocate(MGL1%z_III(0:MGL1%mi+1,0:MGL1%mj+1,0:MGL1%mk+1))
            allocate(MGL1%res_LLL(1:MGL1%mi,1:MGL1%mj,1:MGL1%mk))
            allocate(MGL1%Ap_LLL(1:MGL1%mi,1:MGL1%mj,1:MGL1%mk))
            allocate(MGL1%Az_LLL(1:MGL1%mi,1:MGL1%mj,1:MGL1%mk))

            MGL1%p_III=0.0
            MGL1%r_LLL=0.0
            MGL1%z_III=0.0
            MGL1%res_LLL=0.0
            MGL1%Ap_LLL=0.0
            MGL1%Az_LLL=0.0
        end if
    end subroutine ModMultigrid_InitLinearSolver

    ! Jacobi iteration
    ! For an Ax=b problem we expand it to A_ii * x_i = b_i - sum_j A_ij * x_j,
    ! where A_ij (i/=j) are the off-diagonal elements given by Ai, Aj, and Ak at
    ! the cell faces.
    subroutine ModMultigrid_Jacobi_Ax_b(MGL1)
        implicit none
        type(Multigrid_level)           ::  MGL1
        integer                         ::  mi,mj,mk

        mi=MGL1%mi
        mj=MGL1%mj
        mk=MGL1%mk

        ! Jacobi iteration
        ! Vectorized version

        MGL1%x_III(1:mi,1:mj,1:mk)=MGL1%b_LLL(1:mi,1:mj,1:mk)           - &
            MGL1%Ai_FLL(1:mi,1:mj,1:mk)*MGL1%x_III(0:mi-1,1:mj,1:mk)    - &
            MGL1%Ai_FLL(2:mi+1,1:mj,1:mk)*MGL1%x_III(2:mi+1,1:mj,1:mk)  - &
            MGL1%Aj_LFL(1:mi,1:mj,1:mk)*MGL1%x_III(1:mi,0:mj-1,1:mk)    - &
            MGL1%Aj_LFL(1:mi,2:mj+1,1:mk)*MGL1%x_III(1:mi,2:mj+1,1:mk)  - &
            MGL1%Ak_LLF(1:mi,1:mj,1:mk)*MGL1%x_III(1:mi,1:mj,0:mk-1)    - &
            MGL1%Ak_LLF(1:mi,1:mj,2:mk+1)*MGL1%x_III(1:mi,1:mj,2:mk+1)
        
        MGL1%x_III(1:mi,1:mj,1:mk)=MGL1%x_III(1:mi,1:mj,1:mk)/MGL1%Aii_LLL(1:mi,1:mj,1:mk)
    end subroutine ModMultigrid_Jacobi_Ax_b

    subroutine ModMultigrid_Jacobi_Mz_r(MGL1)
        implicit none
        type(Multigrid_level)           ::  MGL1
        integer                         ::  mi,mj,mk

        mi=MGL1%mi
        mj=MGL1%mj
        mk=MGL1%mk

        ! Jacobi iteration
        ! Vectorized version

        MGL1%z_III(1:mi,1:mj,1:mk)=MGL1%r_LLL(1:mi,1:mj,1:mk)           - &
            MGL1%Ai_FLL(1:mi,1:mj,1:mk)*MGL1%z_III(0:mi-1,1:mj,1:mk)    - &
            MGL1%Ai_FLL(2:mi+1,1:mj,1:mk)*MGL1%z_III(2:mi+1,1:mj,1:mk)  - &
            MGL1%Aj_LFL(1:mi,1:mj,1:mk)*MGL1%z_III(1:mi,0:mj-1,1:mk)    - &
            MGL1%Aj_LFL(1:mi,2:mj+1,1:mk)*MGL1%z_III(1:mi,2:mj+1,1:mk)  - &
            MGL1%Ak_LLF(1:mi,1:mj,1:mk)*MGL1%z_III(1:mi,1:mj,0:mk-1)    - &
            MGL1%Ak_LLF(1:mi,1:mj,2:mk+1)*MGL1%z_III(1:mi,1:mj,2:mk+1)
        
        MGL1%z_III(1:mi,1:mj,1:mk)=MGL1%z_III(1:mi,1:mj,1:mk)/MGL1%Aii_LLL(1:mi,1:mj,1:mk)
    end subroutine ModMultigrid_Jacobi_Mz_r
end module ModMultigrid