module ModHeatConduction

    ! This module contains the linear setup for heat conduction
    ! and calls the linear solver.
    ! Right now we only consider the Spitzer heat flow.

    use ModYinYangTree,            only: YYTree
    use ModBlock,                  only: BlockType,Multigrid_level
    use ModVariables,              only: rho_,te_
    use ModParameters,             only: ni,nj,nk,ng,nvar,MpiRank,MpiSize,Multigrid_nLevels
    use ModCommunication,          only: ModCommunication_SendRecvHC

    implicit none

    contains
end module ModHeatConduction