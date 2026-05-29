module ModInitiation

    use ModBlock,       only:   BlockType
    use ModYinYang,     only:   ModYinyang_CoordConv_0D,&
                                ModYinYang_VecConv_0D
    use ModYinYangTree, only:   YYTree
    use ModParameters,  only:   r_range,ni,nj,nk,ng,&
                                Initiation_type_index,randVelocity_rms,&
                                InitiationB_type_index,Bphi_uniform,&
                                if_involve_B
    use ModConst,       only:   dpi
    contains

    subroutine ModInitiation_DoAll(Tree)
        implicit none
        type(YYTree),target     ::  Tree

        select case (Initiation_type_index)
        case (1)
            call ModInitiation_harmonic(Tree)
        case (2)
            call ModInitiation_all_ones(Tree)
        case (3)
            call ModInitiation_rand_velocity(Tree, randVelocity_rms)
        case default
        end select

        ! Only initialize B if if_involve_B is true.

        if (if_involve_B) then
            select case (InitiationB_type_index)
        case (1)
            call ModInitiation_B_uniform_Bph(Tree, Bphi_uniform)
        end select
        end if
    end subroutine

    subroutine ModInitiation_harmonic(Tree)
        
        implicit none
        type(YYTree),target     ::  Tree
        real(8)                 ::  vec(3)
        real(8)                 ::  coord(3)
        type(BlockType),pointer ::  Block1
        integer                 ::  iLocalBlock,ir,it,ip

        
        do iLocalBlock=1,Tree%nLocalBlocks
            
            Block1=>Tree%LocalBlocks(iLocalBlock)

            do ip=-ng+1,ng+nk
                do it=-ng+1,ng+nj
                    do ir=-ng+1,ng+ni
                        Block1%primitive_IV(ir,it,ip,:)=0.000000
                        coord=[Block1%xi_I(ir),Block1%xj_I(it),Block1%xk_I(ip)]
                        if (.not. Block1%if_yin) coord=ModYinyang_CoordConv_0D(coord)
                        vec=sin(dpi*(coord(1)-r_range(1))/(r_range(2)-r_range(1)))*[1.,1.,1.]*1.e-2*&
                            sin(12.*coord(3))*sin(coord(2))**4
                        Block1%primitive_IV(ir,it,ip,Block1%vr_:Block1%vp_)=vec
                        if (.not. Block1%if_yin) then
                            Block1%primitive_IV(ir,it,ip,Block1%vr_:Block1%vp_)=&
                            ModYinYang_VecConv_0D(coord,vec)
                        end if
                        
                    end do
                    
                end do
                
            end do
            Block1%primitive_rk_IV=Block1%primitive_IV
        end do
    end subroutine ModInitiation_harmonic

    subroutine ModInitiation_all_ones(Tree)
        implicit none
        type(YYTree),target     ::  Tree
        type(BlockType),pointer ::  Block1
        integer                 ::  iLocalBlock,ir,it,ip

        do iLocalBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iLocalBlock)
            do ip=-ng+1,ng+nk
                do it=-ng+1,ng+nj
                    do ir=-ng+1,ng+ni
                        Block1%primitive_IV(ir,it,ip,:)=1.0d0
                    end do
                end do
            end do
            Block1%primitive_rk_IV=Block1%primitive_IV
        end do
    end subroutine ModInitiation_all_ones

    subroutine ModInitiation_rand_velocity(Tree, v_rms)
        implicit none
        type(YYTree), target   :: Tree
        real(8), intent(in)    :: v_rms

        type(BlockType), pointer :: Block1
        integer :: iLocalBlock, ir, it, ip
        real(8) :: coord(3), vec(3)

        do iLocalBlock = 1, Tree%nLocalBlocks
            Block1 => Tree%LocalBlocks(iLocalBlock)
            do ip = 1, nk
                do it = 1, nj
                    do ir = 1, ni
                        call random_number(vec)
                        vec = (vec - 0.5d0) * 2.0d0 * v_rms

                        coord = [Block1%xi_I(ir), Block1%xj_I(it), Block1%xk_I(ip)]
                        if (.not. Block1%if_yin) coord = ModYinyang_CoordConv_0D(coord)

                        Block1%primitive_IV(ir, it, ip, Block1%vr_:Block1%vp_) = vec
                        if (.not. Block1%if_yin) &
                            Block1%primitive_IV(ir, it, ip, Block1%vr_:Block1%vp_) = &
                                ModYinYang_VecConv_0D(coord, vec)
                    end do
                end do
            end do
            Block1%primitive_rk_IV = Block1%primitive_IV
        end do
    end subroutine ModInitiation_rand_velocity

    ! -------------------------------------------------------
    subroutine ModInitiation_B(Tree, B_option, B0)
    ! Set the magnetic field initial condition on top of an existing HD state.
    ! B_option 1: uniform Bph = B0 [G] everywhere (Br = Bt = 0) in the global Yin frame.
    ! -------------------------------------------------------
        implicit none
        type(YYTree), target :: Tree
        integer, intent(in)  :: B_option
        real(8), intent(in)  :: B0

        select case(B_option)
        case(1)
            call ModInitiation_B_uniform_Bph(Tree, B0)
        end select
    end subroutine ModInitiation_B


    ! -------------------------------------------------------
    subroutine ModInitiation_B_uniform_Bph(Tree, B0)
    ! Br = Bt = 0, Bp = B0 (in Yin frame) for all cells including ghost cells.
    ! Yang blocks receive the appropriately rotated components.
    ! psi (GLM scalar) is initialised to zero.
    ! -------------------------------------------------------
        implicit none
        type(YYTree), target     :: Tree
        real(8), intent(in)      :: B0
        type(BlockType), pointer :: Block1
        integer :: iLocalBlock, ir, it, ip
        real(8) :: coord(3), B_vec(3)

        B_vec = [0.0d0, 0.0d0, B0]   ! uniform Bph in Yin frame

        do iLocalBlock = 1, Tree%nLocalBlocks
            Block1 => Tree%LocalBlocks(iLocalBlock)

            do ip = -ng+1, ng+nk
                do it = -ng+1, ng+nj
                    do ir = -ng+1, ng+ni

                        ! Default: set field in Yin frame
                        Block1%primitive_IV(ir,it,ip,Block1%br_:Block1%bp_) = B_vec

                        ! Yang block: rotate B_vec from Yin frame to Yang frame
                        if (.not. Block1%if_yin) then
                            coord = [Block1%xi_I(ir), Block1%xj_I(it), Block1%xk_I(ip)]
                            coord = ModYinyang_CoordConv_0D(coord)   ! Yang → Yin
                            Block1%primitive_IV(ir,it,ip,Block1%br_:Block1%bp_) = &
                                ModYinYang_VecConv_0D(coord, B_vec)
                        end if

                        ! GLM scalar starts at zero
                        if (Block1%psi_ > 0) &
                            Block1%primitive_IV(ir,it,ip,Block1%psi_) = 0.0d0
                    end do
                end do
            end do

            ! Sync RK copy for B variables
            Block1%primitive_rk_IV(:,:,:,Block1%br_:Block1%bp_) = &
                Block1%primitive_IV(:,:,:,Block1%br_:Block1%bp_)
            if (Block1%psi_ > 0) &
                Block1%primitive_rk_IV(:,:,:,Block1%psi_) = 0.0d0
        end do
    end subroutine ModInitiation_B_uniform_Bph

end module ModInitiation
