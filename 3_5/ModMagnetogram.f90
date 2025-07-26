module ModMagnetogram
    ! This module contains the subroutines for the magnetogram.

    use ModYinYang,                only: ModYinYang_CoordConv_0D
    use ModYinYangTree,            only: YYTree
    use ModBlock,                  only: BlockType
    use ModParameters,             only: ni,nj,nk,ng,nvar,MpiRank,MpiSize,Multigrid_nLevels

    implicit none

    contains

    subroutine ModMagnetogram_dipole_magnetogram_ALL(Tree,n_dipoles,xyz_dipoles,q_dipoles,l_dipoles)
        implicit none
        type(YYTree),target         ::  Tree
        integer,intent(in)          ::  n_dipoles
        real,intent(in)             ::  xyz_dipoles(3,n_dipoles)
        real,intent(in)             ::  q_dipoles(n_dipoles)
        real,intent(in)             ::  l_dipoles(3,n_dipoles)

        integer                     ::  iBlock,iDipole
        integer                     ::  j,k
        real                        ::  xyz_positive(3),xyz_negative(3)
        real                        ::  rtp_face_center(3),xyz_face_center(3)
        real                        ::  xyz_positive_to_face(3),xyz_negative_to_face(3)
        real                        ::  r_positive,r_negative
        type(BlockType),pointer     ::  Block1

        do iBlock = 1, Tree%nLocalBlocks
            Block1 => Tree%LocalBlocks(iBlock)
            

            ! Only do if it's PFSS and bottom
            if (Block1%if_PFSS .and. Block1%if_bottom) then

                ! Allocate the magnetogram
                allocate(Block1%magnetogram_LL(1:nj,1:nk))
                Block1%magnetogram_LL=0.0

                ! Get the magnetogram by looping over the dipoles
                do iDipole = 1, n_dipoles
                    ! Get the position of the positive and negative poles
                    xyz_positive=xyz_dipoles(:,iDipole)+l_dipoles(:,iDipole)/2.0
                    xyz_negative=xyz_dipoles(:,iDipole)-l_dipoles(:,iDipole)/2.0

                    ! Loop over the bottom boundary
                    do j=1,nj
                        do k=1,nk
                            ! Get the position of this face center
                            ! Since it's at the bottom face, r is xi_F(1), j is xj_F(), k is xk_I(1)
                            rtp_face_center=[Block1%xi_F(1),Block1%xj_F(j),Block1%xk_I(k)]
                            if(Block1%if_yin) rtp_face_center=ModYinYang_CoordConv_0D(rtp_face_center)
                            xyz_face_center=rtp_face_center(1)*&
                                [cos(rtp_face_center(3))*sin(rtp_face_center(2)),&
                                sin(rtp_face_center(3))*sin(rtp_face_center(2)),&
                                cos(rtp_face_center(2))]

                            ! The vector from positive/negative pole to face center
                            xyz_positive_to_face=xyz_face_center-xyz_positive
                            xyz_negative_to_face=xyz_face_center-xyz_negative

                            ! Get the distance between the face center and the positive and negative poles
                            r_positive=sqrt(sum(xyz_positive_to_face**2))
                            r_negative=sqrt(sum(xyz_negative_to_face**2))

                            ! Get the Br
                            ! Projection factor is sum(xyz_positive_to_face*xyz_face_center)/r_positive
                            ! The factor 1/Block1%xi_F(1) is for normalization of xyz_face_center.
                            ! The additional 1/r_positive(negative) is for normalization of the vector.
                            Block1%magnetogram_LL(j,k)=Block1%magnetogram_LL(j,k)+&
                                q_dipoles(iDipole)*&
                                sum(xyz_positive_to_face*xyz_face_center/Block1%xi_F(1))/r_positive**3-&
                                q_dipoles(iDipole)*&
                                sum(xyz_negative_to_face*xyz_face_center/Block1%xi_F(1))/r_negative**3
                        end do ! j
                    end do ! k
                end do ! iDipole
            end if
        end do ! iBlock
    end subroutine ModMagnetogram_dipole_magnetogram_ALL
end module ModMagnetogram