Module Cutcell
    !    
    USE PrecisionVar
    USE Mesh
    USE StateVariables, only : ight,jght,R1,xc,yc,pi
    !
    Implicit none
    !
    Private
    !
    Real(kind=dp),parameter,public :: epsi = 1.d-3
    !
    Type,public:: Cell
      Integer(kind=it4b),dimension(:,:,:),allocatable   :: Cell_Type
                          ! CellType: 0 fluid cell and not target cell
                          !           1 boudary cell
                          !           2 solid cell
      Integer(kind=it4b),dimension(:,:,:),allocatable   :: MoExCell
      ! Cell number for Poisson solving equation only for Pressure cell
      Integer(kind=it4b),dimension(:,:,:),allocatable   :: PosNu
      Integer(kind=it4b):: ExtCell
      Integer(kind=it4b),dimension(:,:,:,:),allocatable :: MsCe
      ! Cell number of Implicit solving for diffusive term for u-velocity and v-velocity
      ! Fluid volume fraction and levelset function
      Real(kind=dp),dimension(:,:,:),allocatable   :: Vof,Phi,VofL,PhiL
      ! Normal vector of boundary
      Real(kind=dp),dimension(:,:,:),allocatable   :: nx,ny,nz,nxL,nyL,nzL
      ! The center of cell
      Real(kind=dp),dimension(:,:,:,:),allocatable :: Cell_Cent
      ! The area of edge in x,y,z direction
      Real(kind=dp),dimension(:,:,:),allocatable   :: EEArea
      Real(kind=dp),dimension(:,:,:),allocatable   :: NEArea,TEArea
      ! The cellface center of velocity cell
      Real(kind=dp),dimension(:,:,:),allocatable   :: WlLh,delh
      Real(kind=dp),dimension(:,:,:),allocatable   :: WePr
      Real(kind=dp),dimension(:,:,:,:),allocatable :: FCE,FCN,FCT
      Real(kind=dp),dimension(:,:,:),allocatable   :: EtaE,EtaN,EtaT,DAlE,DAlN,DAlT
      Real(kind=dp),dimension(:,:,:),allocatable   :: AlE,AlN,AlT,SxE,SyN,SzT
    End type
    !
    Public:: GridPreprocess,DefineMomentumExchangeCell
    Public:: NumberExternalCell,NewCellFace
    !
    Interface GridPreprocess
      Module procedure GridPreprocess
    End interface GridPreprocess
    !
    Interface DefineMomentumExchangeCell
      Module procedure DefineMomentumExchangeCell
    End interface
    !
    Interface NumberExternalCell
      Module procedure NumberExternalCell
    End interface
    !
    Interface NewCellFace
      Module procedure NewCellFace
    End interface
    !
    Contains
    !
      Subroutine GridPreprocess(itt,PGrid,UGrid,VGrid,WGrid, PCell,UCell,VCell,WCell)
        !                                                     
        Implicit none
        !
        Integer(kind=it8b),intent(in):: itt
        Type(Grid),intent(in):: PGrid,UGrid,VGrid,WGrid
        Type(Cell),intent(inout):: PCell,UCell,VCell,WCell

        Integer(kind=it4b):: i,j,k
        Real(kind=dp):: Face(3),minVof,CellFaceArea
        !
        ! defines the type of cell (liquid, gas (solid ?), or cut cell)
        !
        Call DefineCell(PCell)
        Call DefineCell(UCell)
        Call DefineCell(VCell)
        Call DefineCell(WCell)
        !
        ! compute the face centers, face area covered by fluid
        !
        Call CellGeoCal(PGrid,PCell)
        Call CellGeoCal(UGrid,UCell)
        Call CellGeoCal(VGrid,VCell)
        Call CellGeoCal(WGrid,WCell)
        !
        ! Correct Cell face area of U,V,W Cell
        !
        minVof=1.d0
        !
        Do j = 1,Jmax
          Do k = 1,Kmax
            Do i = 1,Imax-1
              !
              If(UCell%Cell_Type(i,j,k)/=1.and.UCell%Cell_Type(i+1,j,k)/=1)then
                !      
                If(dabs(UGrid%dx(i,j,k)-UGrid%dx(i+1,j,k))>1.d-10) then
                  !      
                  UCell%FCE(i,j,k,1) = 0.5d0*PGrid%dx(i+1,j,k)
                  UCell%FCN(i,j,k,1) = 0.d0
                  UCell%FCT(i,j,k,1) = 0.d0
                  !
                End if
                !
              End if
              !
            End do
            !
            UCell%FCN(Imax,j,k,1) = 0.d0
            UCell%FCT(Imax,j,k,1) = 0.d0
            !
          End do
        End do
        !
        minVof=1.d0
        !
        Do i = 1,Imax
          Do k = 1,Kmax
            Do j = 1,Jmax-1
              !
              If(VCell%Cell_Type(i,j,k)/=1.and.VCell%Cell_Type(i,j+1,k)/=1)then
                !      
                If(dabs(VGrid%dy(i,j,k)-VGrid%dy(i,j+1,k))>1.d-10) then
                  !      
                  VCell%FCN(i,j,k,2) = 0.5d0*PGrid%dy(i,j+1,k)
                  VCell%FCE(i,j,k,2) = 0.d0
                  VCell%FCT(i,j,k,2) = 0.d0
                  !
                End if
                !
              End if
              !
            End do
            !
            VCell%FCE(i,Jmax,k,2) = 0.d0
            VCell%FCT(i,Jmax,k,2) = 0.d0
            !
          End do
        End do
        !
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 1,Kmax-1
              !
              If(WCell%Cell_Type(i,j,k)/=1.and.WCell%Cell_Type(i,j,k+1)/=1)then
                !      
                If(dabs(WGrid%dz(i,j,k)-WGrid%dz(i,j,k+1))>1.d-10) then
                  !      
                  WCell%FCT(i,j,k,3) = 0.5d0*PGrid%dz(i,j,k+1)
                  WCell%FCE(i,j,k,3) = 0.d0
                  WCell%FCN(i,j,k,3) = 0.d0
                  !
                End if
                !
              End if
              !
            End do
            !
            WCell%FCE(i,j,Kmax,3) = 0.d0
            WCell%FCN(i,j,Kmax,3) = 0.d0
            !
          End do
        End do
        !
        ! Correct Cell Center of U,V,W Cell
        !
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 1,Kmax
              !
              PCell%Cell_Cent(i,j,k,1) = 0.d0
              PCell%Cell_Cent(i,j,k,2) = 0.d0
              PCell%Cell_Cent(i,j,k,3) = 0.d0
              !
              If(PCell%vof(i,j,k)<1.d0.and.PCell%vof(i,j,k)>0.d0) then
                !        
                ! East Cell Face
                Call FaceGeoCal(PGrid%dy(i,j,k),PGrid%dz(i,j,k),PCell%ny(i,j,k),&
                                PCell%nz(i,j,k),PCell%nx(i,j,k),               &
                                0.d0,PCell%phi(i,j,k),                         &
                                UCell%EEArea(i-1,j,k),FaCe)
                !        
                UCell%FCE(i-1,j,k,1) = FaCe(3)+PGrid%dx(i,j,k)/2.d0
                UCell%FCE(i-1,j,k,2) = FaCe(1);UCell%FCE(i-1,j,k,3) = FaCe(2)
                !
                ! North Cell Face
                Call FaceGeoCal(PGrid%dz(i,j,k),PGrid%dx(i,j,k),PCell%nz(i,j,k),&
                                PCell%nx(i,j,k),PCell%ny(i,j,k),               &
                                0.d0,PCell%phi(i,j,k),                         &
                                VCell%NEArea(i,j-1,k),FaCe)
                !        
                VCell%FCN(i,j-1,k,2) = FaCe(3)+PGrid%dy(i,j,k)/2.d0
                VCell%FCN(i,j-1,k,1) = FaCe(2);VCell%FCN(i,j-1,k,3) = FaCe(1)
                !
                ! Top Cell Face
                Call FaceGeoCal(PGrid%dx(i,j,k),PGrid%dy(i,j,k),PCell%nx(i,j,k),&
                                PCell%ny(i,j,k),PCell%nz(i,j,k),               &
                                0.d0,PCell%phi(i,j,k),                         &
                                WCell%TEArea(i,j,k-1),FaCe)
                !        
                WCell%FCT(i,j,k-1,3) = FaCe(3)+PGrid%dz(i,j,k)/2.d0
                WCell%FCT(i,j,k-1,1) = FaCe(1);WCell%FCT(i,j,k-1,2) = FaCe(2)
                !
              Elseif(PCell%vof(i,j,k)>0.5d0) then
                !       
                UCell%FCE(i-1,j,k,1) = PGrid%dx(i,j,k)/2.d0
                UCell%FCE(i-1,j,k,2) = 0.d0;UCell%FCE(i-1,j,k,3) = 0.d0
                UCell%EEArea(i-1,j,k) = 1.d0
                !
                VCell%FCN(i,j-1,k,2) = PGrid%dy(i,j,k)/2.d0
                VCell%FCN(i,j-1,k,1) = 0.d0;VCell%FCN(i,j-1,k,3) = 0.d0
                VCell%NEArea(i,j-1,k) = 1.d0
                !
                WCell%FCT(i,j,k-1,3) = PGrid%dz(i,j,k)/2.d0
                WCell%FCT(i,j,k-1,1) = 0.d0;WCell%FCT(i,j,k-1,2) = 0.d0
                WCell%TEArea(i,j,k-1) = 1.d0
                !
              Else
                !
                UCell%EEArea(i-1,j,k) = 0.d0;UCell%FCE(i-1,j,k,:) = 0.d0
                VCell%NEArea(i,j-1,k) = 0.d0;VCell%FCN(i,j-1,k,:) = 0.d0
                WCell%TEArea(i,j,k-1) = 0.d0;WCell%FCT(i,j,k-1,:) = 0.d0
                !
              End if
              !
            End do
          End do
        End do
        !
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 1,Kmax
              !
              ! The New Cell Center for U
              If(UCell%vof(i,j,k)<1.d0.and.UCell%vof(i,j,k)>0.d0) then
                !      
                Call FaceGeoCal(UGrid%dy(i,j,k),UGrid%dz(i,j,k),UCell%ny(i,j,k),&
                                UCell%nz(i,j,k),UCell%nx(i,j,k),               &
                                0.d0,UCell%phi(i,j,k),                         &
                                CellFaceArea,FaCe)
                !        
                PCell%FCE(i,j,k,1) = FaCe(3)+PGrid%dx(i,j,k)/2.d0
                PCell%FCE(i,j,k,2) = FaCe(1);PCell%FCE(i,j,k,3) = FaCe(2)
                UCell%Cell_Cent(i,j,k,1) = FaCe(3)
                UCell%Cell_Cent(i,j,k,2) = FaCe(1)
                UCell%Cell_Cent(i,j,k,3) = FaCe(2)
                !
              Elseif(UCell%vof(i,j,k)>0.5d0) then
                !      
                PCell%EEArea(i,j,k) = 1.d0
                PCell%FCE(i,j,k,1) = PGrid%dx(i,j,k)/2.d0
                PCell%FCE(i,j,k,2) = 0.d0
                PCell%FCE(i,j,k,3) = 0.d0
                UCell%Cell_Cent(i,j,k,:) = 0.d0
                !
              Else
                !      
                PCell%EEArea(i,j,k) = 0.d0
                PCell%FCE(i,j,k,:) = 0.d0
                UCell%Cell_Cent(i,j,k,:) = 0.d0
                !
              End if
              !
              ! New Cell Center for V
              If(VCell%vof(i,j,k)<1.d0.and.VCell%vof(i,j,k)>0.d0) then
                !      
                Call FaceGeoCal(VGrid%dz(i,j,k),VGrid%dx(i,j,k),VCell%nz(i,j,k),&
                                VCell%nx(i,j,k),VCell%ny(i,j,k),               &
                                0.d0,VCell%phi(i,j,k),                         &
                                CellFaceArea,FaCe)
                !        
                PCell%FCN(i,j,k,1) = FaCe(2);PCell%FCN(i,j,k,3) = FaCe(1)
                PCell%FCN(i,j,k,2) = FaCe(3)+PGrid%dy(i,j,k)/2.d0
                VCell%Cell_Cent(i,j,k,1) = FaCe(2)
                VCell%Cell_Cent(i,j,k,2) = FaCe(3)
                VCell%Cell_Cent(i,j,k,3) = FaCe(1)
                !
              Elseif(VCell%vof(i,j,k)>0.5d0) then
                !      
                PCell%NEArea(i,j,k) = 1.d0
                PCell%FCN(i,j,k,1) = 0.d0;PCell%FCN(i,j,k,3) = 0.d0
                PCell%FCN(i,j,k,2) = PGrid%dy(i,j,k)/2.d0
                VCell%Cell_Cent(i,j,k,:) = 0.d0
                !
              Else
                !      
                PCell%NEArea(i,j,k) = 0.d0
                PCell%FCN(i,j,k,:) = 0.d0
                VCell%Cell_Cent(i,j,k,:) = 0.d0
                !
              End if
              !
              ! New Cell Center for WCell
              If(WCell%vof(i,j,k)<1.d0.and.WCell%vof(i,j,k)>0.d0) then
                !      
                Call FaceGeoCal(WGrid%dx(i,j,k),WGrid%dy(i,j,k),WCell%nx(i,j,k),&
                                WCell%ny(i,j,k),WCell%nz(i,j,k),               &
                                0.d0,WCell%phi(i,j,k),                         &
                                CellFaceArea,FaCe)
                !           
                PCell%FCT(i,j,k,1) = FaCe(1);PCell%FCT(i,j,k,2) = FaCe(2)
                PCell%FCT(i,j,k,3) = FaCe(3)+PGrid%dz(i,j,k)/2.d0
                WCell%Cell_Cent(i,j,k,1) = FaCe(1)
                WCell%Cell_Cent(i,j,k,2) = FaCe(2)
                WCell%Cell_Cent(i,j,k,3) = FaCe(3)
                !
              Elseif(WCell%vof(i,j,k)>0.5d0) then
                !      
                PCell%TEArea(i,j,k) = 1.d0
                PCell%FCT(i,j,k,1) = 0.d0;PCell%FCT(i,j,k,2) = 0.d0
                PCell%FCT(i,j,k,3) = PGrid%dz(i,j,k)/2.d0
                WCell%Cell_Cent(i,j,k,:) = 0.d0
                !
              Else
                !      
                PCell%TEArea(i,j,k) = 0.d0
                PCell%FCT(i,j,k,:) = 0.d0
                WCell%Cell_Cent(i,j,k,:) = 0.d0
                !
              End if
              !
            End do
          End do
        End do
        !
      End subroutine GridPreprocess
      !
      ! defines the type of cell (liquid, gas, or cut cell)
      !
      Subroutine DefineCell(TCell)
        !      
        Implicit none
        !
        Type(Cell),intent(inout):: TCell

        Integer(kind=it4b):: i,j,k
        !
        ! Define cell basing the volume fraction of fluid
        !
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 1,Kmax
              !
              if(TCell%vof(i,j,k)>=1.d0-epsi) then ! liquid phase 
                !      
                TCell%Cell_Type(i,j,k) = 0  ! External cell
                TCell%vofL(i,j,k)=TCell%VofL(i,j,k)/TCell%vof(i,j,k)
                TCell%vof(i,j,k)=1.d0
                !
              end if  
              !
              If(TCell%vof(i,j,k)<=epsi) then ! gas phase
                TCell%Cell_Type(i,j,k) = 2  ! Internal Cell
                TCell%vof(i,j,k)=0.d0
                TCell%vofL(i,j,k)=0.d0
              End if
              !
              If(TCell%vof(i,j,k)>=epsi.and.TCell%vof(i,j,k)<=1.d0-epsi)then ! Cut cell 
                TCell%Cell_Type(i,j,k) = 1  
              End if
              !                    
            End do
          End do
        End do
        !
      End subroutine DefineCell
      !
      ! Calculate the Edge Area of Each Cell
      !
      Subroutine CellGeoCal(TGrid, TCell)
        !      
        Implicit none
        !
        Type(Grid),intent(in):: TGrid
        Type(Cell),intent(inout):: TCell
        !
        Integer(kind=it4b):: i,j,k
        Real(kind=dp):: MaxFace,AverageArea
        Real(kind=dp):: tol,FaCe(3)
        Real(kind=dp),dimension(:,:,:,:),allocatable::CeE,CeW,CeN,CeS,CeT,CeB
        Real(kind=dp),dimension(:,:,:),allocatable:: WEArea,SEArea,BEArea
        !
        Allocate(WEArea(Imax,Jmax,Kmax))
        Allocate(SEArea(Imax,Jmax,Kmax))
        Allocate(BEArea(Imax,Jmax,Kmax))
        Allocate(CeE(Imax,Jmax,Kmax,3))
        Allocate(CeW(Imax,Jmax,Kmax,3))
        Allocate(CeN(Imax,Jmax,Kmax,3))
        Allocate(CeS(Imax,Jmax,Kmax,3))
        Allocate(CeT(Imax,Jmax,Kmax,3))
        Allocate(CeB(Imax,Jmax,Kmax,3))
        !
        tol = 1.d-20
        !
        WEArea = 1.d0
        SEArea = 1.d0
        BEArea = 1.d0
        !
        TCell%WlLh = 0.d0
        TCell%EEArea(:,:,:) = 1.d0
        TCEll%NEArea(:,:,:) = 1.d0
        TCell%TEArea(:,:,:) = 1.d0
        !
        ! compute the cell center and the covered area
        !
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 1,Kmax
              !
              TCell%Cell_Cent(i,j,k,1) = 0.d0
              TCell%Cell_Cent(i,j,k,2) = 0.d0
              TCell%Cell_Cent(i,j,k,3) = 0.d0
              !
              ! Top Cell Face
              !
              Call FaceGeoCal(TGrid%dx(i,j,k),TGrid%dy(i,j,k),TCell%nx(i,j,k), &
                              TCell%ny(i,j,k),TCell%nz(i,j,k),                 &
                              TGrid%dz(i,j,k)/2.d0,TCell%phi(i,j,k),           &
                              TCell%TEArea(i,j,k),FaCe)
              CeT(i,j,k,1) = FaCe(1)
              CeT(i,j,k,2) = FaCe(2)
              CeT(i,j,k,3) = FaCe(3)
              !                
              ! Bottom Cell Face
              !
              Call FaceGeoCal(TGrid%dx(i,j,k),TGrid%dy(i,j,k),TCell%nx(i,j,k), &
                              TCell%ny(i,j,k),TCell%nz(i,j,k),                 &
                              -TGrid%dz(i,j,k)/2.d0,TCell%phi(i,j,k),          &
                              BEArea(i,j,k),FaCe)
              CeB(i,j,k,1) = FaCe(1)
              CeB(i,j,k,2) = FaCe(2)
              CeB(i,j,k,3) = FaCe(3)
              !                
              ! East Cell Face
              !
              Call FaceGeoCal(TGrid%dy(i,j,k),TGrid%dz(i,j,k),TCell%ny(i,j,k), &
                              TCell%nz(i,j,k),TCell%nx(i,j,k),                 &
                              TGrid%dx(i,j,k)/2.d0,TCell%phi(i,j,k),           &
                              TCell%EEArea(i,j,k),FaCe)
              CeE(i,j,k,1) = FaCe(3)
              CeE(i,j,k,2) = FaCe(1)
              CeE(i,j,k,3) = FaCe(2)
              !                
              ! West Cell Face
              !
              Call FaceGeoCal(TGrid%dy(i,j,k),TGrid%dz(i,j,k),TCell%ny(i,j,k), &
                              TCell%nz(i,j,k),TCell%nx(i,j,k),                 &
                              -TGrid%dx(i,j,k)/2.d0,TCell%phi(i,j,k),          &
                              WEArea(i,j,k),FaCe)
              CeW(i,j,k,1) = FaCe(3)
              CeW(i,j,k,2) = FaCe(1)
              CeW(i,j,k,3) = FaCe(2)
              !                
              ! North Cell Face
              !
              Call FaceGeoCal(TGrid%dz(i,j,k),TGrid%dx(i,j,k),TCell%nz(i,j,k), &
                              TCell%nx(i,j,k),TCell%ny(i,j,k),                 &
                              TGrid%dy(i,j,k)/2.d0,TCell%phi(i,j,k),           &
                              TCell%NEArea(i,j,k),FaCe)
              CeN(i,j,k,1) = FaCe(2)
              CeN(i,j,k,2) = FaCe(3)
              CeN(i,j,k,3) = FaCe(1)
              !                
              ! South Cell Face
              !
              Call FaceGeoCal(TGrid%dz(i,j,k),TGrid%dx(i,j,k),TCell%nz(i,j,k), &
                              TCell%nx(i,j,k),TCell%ny(i,j,k),                 &
                              -TGrid%dy(i,j,k)/2.d0,TCell%phi(i,j,k),          &
                              SEArea(i,j,k),FaCe)
              CeS(i,j,k,1) = FaCe(2)
              CeS(i,j,k,2) = FaCe(3)
              CeS(i,j,k,3) = FaCe(1)
              !                
              ! Calculate the face area of internal cell
              ! An edge sharing with internal cell will be set up at 0
              ! The other edges will be set up at actual value.
              !
              TCell%WlLh(i,j,k) = dsqrt(((TCell%EEArea(i,j,k)-WEArea(i,j,k))   &
                                *TGrid%dy(i,j,k)*TGrid%dz(i,j,k))**2.d0+       &
                                ((TCell%NEArea(i,j,k)-SEArea(i,j,k))           &
                                *TGrid%dx(i,j,k)*TGrid%dz(i,j,k))**2.d0+       &
                                ((TCell%TEArea(i,j,k)-BEArea(i,j,k))           &
                                *TGrid%dx(i,j,k)*TGrid%dy(i,j,k))**2.d0)
            End do
          End do
        End do
        !
        ! Adjusting the area of sharing face between 2 cells
        !
        ! East 
        Do j = 1,Jmax
          Do k = 1,Kmax
            Do i = 2,Imax
              !
              AverageArea = (WEArea(i,j,k)**2.d0+TCell%EEArea(i-1,j,k)**2.d0)  &
                           /(WEArea(i,j,k)+TCell%EEArea(i-1,j,k)+tol)
              TCell%FCE(i-1,j,k,2) = (CeW(i,j,k,2)*WEArea(i,j,k)+CeE(i-1,j,k,2)&
                                     *TCell%EEArea(i-1,j,k))/(WEArea(i,j,k)+   &
                                      TCell%EEArea(i-1,j,k)+tol)
              TCell%FCE(i-1,j,k,3) = (CeW(i,j,k,3)*WEArea(i,j,k)+CeE(i-1,j,k,3)&
                                     *TCell%EEArea(i-1,j,k))/(WEArea(i,j,k)+   &
                                      TCell%EEArea(i-1,j,k)+tol)
              TCell%FCE(i-1,j,k,1) = CeE(i-1,j,k,1)
              TCell%EEArea(i-1,j,k) = AverageArea
              !
            End do
          End do
        End do
        !
        ! North
        Do i = 1,Imax
          Do k = 1,Kmax
            Do j = 2,Jmax
              !
              AverageArea = (SEArea(i,j,k)**2.d0+TCell%NEArea(i,j-1,k)**2.d0)  &
                           /(SEArea(i,j,k)+TCell%NEArea(i,j-1,k)+tol)
              TCell%FCN(i,j-1,k,1) = (CeS(i,j,k,1)*SEArea(i,j,k)+CeN(i,j-1,k,1)&
                                     *TCell%NEArea(i,j-1,k))/(SEArea(i,j,k)+   &
                                      TCell%NEArea(i,j-1,k)+tol)
              TCell%FCN(i,j-1,k,3) = (CeS(i,j,k,3)*SEArea(i,j,k)+CeN(i,j-1,k,3)&
                                     *TCell%NEArea(i,j-1,k))/(SEArea(i,j,k)+   &
                                      TCell%NEArea(i,j-1,k)+tol)
              TCell%FCN(i,j-1,k,2) = CeN(i,j-1,k,2)
              TCell%NEArea(i,j-1,k) = AverageArea
              !
            End do
          End do
        End do
        !
        ! Top
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 2,Kmax
              !
              AverageArea = (BEArea(i,j,k)**2.d0+TCell%TEArea(i,j,k-1)**2.d0)  &
                           /(BEArea(i,j,k)+TCell%TEArea(i,j,k-1)+tol)
              TCell%FCT(i,j,k-1,1) = (CeB(i,j,k,1)*BEArea(i,j,k)+CeT(i,j,k-1,1)&
                                     *TCell%TEArea(i,j,k-1))/(BEArea(i,j,k)+   &
                                      TCell%TEArea(i,j,k-1)+tol)
              TCell%FCT(i,j,k-1,2) = (CeB(i,j,k,1)*BEArea(i,j,k)+CeT(i,j,k-1,2)&
                                     *TCell%TEArea(i,j,k-1))/(BEArea(i,j,k)+   &
                                      TCell%TEArea(i,j,k-1)+tol)
              TCell%FCT(i,j,k-1,3) = CeT(i,j,k-1,3)
              TCell%TEArea(i,j,k-1) = AverageArea
              !
            End do
          End do
        End do
        !
        Do i = 2,Imax
          Do j = 2,Jmax
            Do k = 2,Kmax
              !
              MaxFace = dmax1(TCell%EEArea(i-1,j,k),TCell%EEArea(i,j,k),       &
                              TCell%NEArea(i,j-1,k),TCell%NEArea(i,j,k),       &
                              TCell%TEArea(i,j,k-1),TCell%TEArea(i,j,k))
              !        
              If(.not.allocated(TCell%MoExCell)) then
                If(MaxFace<1.d-2) TCell%Cell_Type(i,j,k)=2
              End if
              !
            End do
          End do
        End do
        !
        do i = 1,Imax
          do j = 1,Jmax
            do k = 1,Kmax
              !
              if(TCell%Cell_Type(i,j,k)==2) then ! Set every cell face of solid cell to 0
                TCell%EEArea(i-1,j,k)=0.d0; TCell%EEArea(i,j,k)=0.d0
                TCell%NEArea(i,j-1,k)=0.d0; TCell%NEArea(i,j,k)=0.d0
                TCell%TEArea(i,j,k-1)=0.d0; TCell%TEArea(i,j,k)=0.d0
              end if
              !
              if(TCell%Cell_Type(i,j,k)==1) then
                if(i+1<=Imax) then
                  if(TCell%Cell_Type(i+1,j,k)==2) then
                    TCell%EEArea(i,j,k) = 0.d0
                  end if
                end if
                if(j+1<=Jmax) then
                  if(TCell%Cell_Type(i,j+1,k)==2) then
                    TCell%NEArea(i,j,k) = 0.d0
                  end if
                end if
                if(k+1<=Kmax) then
                  if(TCell%Cell_Type(i,j,k+1)==2) then
                    TCell%TEArea(i,j,k) = 0.d0
                  end if
                end if
                !
              end if
              !
            end do
          end do
        end do
        !
        ! Boundary conditions
        !
        TCell%FCE(0,:,:,1) = TCell%FCE(1,:,:,1)
        TCell%FCE(0,:,:,2) = TCell%FCE(1,:,:,2)
        TCell%FCE(0,:,:,3) = TCell%FCE(1,:,:,3)
        TCell%FCN(0,:,:,1) = TCell%FCN(1,:,:,1)
        TCell%FCN(0,:,:,2) = TCell%FCN(1,:,:,2)
        TCell%FCN(0,:,:,3) = TCell%FCN(1,:,:,3)
        TCell%FCT(0,:,:,1) = TCell%FCT(1,:,:,1)
        TCell%FCT(0,:,:,2) = TCell%FCT(1,:,:,2)
        TCell%FCT(0,:,:,3) = TCell%FCT(1,:,:,3)

        TCell%FCE(Imax,:,:,1) = TCell%FCE(Imax-1,:,:,1)
        TCell%FCE(Imax,:,:,2) = TCell%FCE(Imax-1,:,:,2)
        TCell%FCE(Imax,:,:,3) = TCell%FCE(Imax-1,:,:,3)

        TCell%FCE(:,0,:,1) = TCell%FCE(:,1,:,1)
        TCell%FCE(:,0,:,2) = TCell%FCE(:,1,:,2)
        TCell%FCE(:,0,:,3) = TCell%FCE(:,1,:,3)
        TCell%FCN(:,0,:,1) = TCell%FCN(:,1,:,1)
        TCell%FCN(:,0,:,2) = TCell%FCN(:,1,:,2)
        TCell%FCN(:,0,:,3) = TCell%FCN(:,1,:,3)
        TCell%FCT(:,0,:,1) = TCell%FCT(:,1,:,1)
        TCell%FCT(:,0,:,2) = TCell%FCT(:,1,:,2)
        TCell%FCT(:,0,:,3) = TCell%FCT(:,1,:,3)

        TCell%FCN(:,Jmax,:,1) = TCell%FCN(:,Jmax-1,:,1)
        TCell%FCN(:,Jmax,:,2) = TCell%FCN(:,Jmax-1,:,2)
        TCell%FCN(:,Jmax,:,3) = TCell%FCN(:,Jmax-1,:,3)

        TCell%FCE(:,:,0,1) = TCell%FCE(:,:,1,1)
        TCell%FCE(:,:,0,2) = TCell%FCE(:,:,1,2)
        TCell%FCE(:,:,0,3) = TCell%FCE(:,:,1,3)
        TCell%FCN(:,:,0,1) = TCell%FCN(:,:,1,1)
        TCell%FCN(:,:,0,2) = TCell%FCN(:,:,1,2)
        TCell%FCN(:,:,0,3) = TCell%FCN(:,:,1,3)
        TCell%FCT(:,:,0,1) = TCell%FCT(:,:,1,1)
        TCell%FCT(:,:,0,2) = TCell%FCT(:,:,1,2)
        TCell%FCT(:,:,0,3) = TCell%FCT(:,:,1,3)

        TCell%FCT(:,:,Kmax,1) = TCell%FCT(:,:,Kmax-1,1)
        TCell%FCT(:,:,Kmax,2) = TCell%FCT(:,:,Kmax-1,2)
        TCell%FCT(:,:,Kmax,3) = TCell%FCT(:,:,Kmax-1,3)
        !
        Do i = 2,Imax
          Do j = 2,Jmax
            Do k = 2,Kmax
              !
              If(dabs(TCell%vof(i,j,k)-1.d0)<epsi) then 
              ! When the cell is normal but one of its cell face is 0.  
                If(TCell%EEArea(i-1,j,k)==0.d0.or.TCell%EEArea(i,j,k)==0.d0)   &
                                                                           then
                  TCell%WlLh(i,j,k) = TGrid%dy(i,j,k)*TGrid%dz(i,j,k)
                  TCell%delh(i,j,k) = TGrid%dx(i,j,k)/2.d0
                End if
                If(TCell%NEArea(i,j-1,k)==0.d0.or.TCell%NEARea(i,j,k)==0.d0)   &
                                                                           then
                  TCell%WlLh(i,j,k) = TGrid%dx(i,j,k)*TGrid%dz(i,j,k)
                  TCell%delh(i,j,k) = TGrid%dy(i,j,k)/2.d0
                End if
                If(TCell%TEArea(i,j,k-1)==0.d0.or.TCell%TEArea(i,j,k)==0.d0)   &
                                                                           then
                  TCell%WlLh(i,j,k) = TGrid%dx(i,j,k)*TGrid%dy(i,j,k)
                  TCell%delh(i,j,k) = TGrid%dz(i,j,k)/2.d0
                End if
              End if
              !
            End do
          End do
        End do
        !
        Deallocate(CeE,CeW,CeN,CeS,CeT,CeB,WEArea,BEArea,SEArea)
        !
      End subroutine CellGeoCal
      !
      Subroutine DefineMomentumExchangeCell(PCell, UCell,VCell,WCell)
        !      
        Implicit none
        !
        Type(Cell),intent(in):: PCell
        Type(Cell),intent(inout):: UCell,VCell,WCell

        Integer(kind=it4b):: i,j,k
        !
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 1,Kmax
              !
              ! Define the U-cell applying the momentum exchange term
              !
              If(UCell%Cell_Type(i,j,k)==1) then ! cut cell
                !      
                If(i<Imax) then
                  !      
                  If(Pcell%Cell_Type(i,j,k)==2.or.PCell%Cell_Type(i+1,j,k)==2) &
                                                                           then
                    UCell%MoExCell(i,j,k) = 1
                    !
                  End if
                  !
                Else
                  !      
                  If(PCell%Cell_Type(i,j,k)==2) UCell%MoExCell(i,j,k) = 1
                  !
                End if
                !
              End if
              !
              ! Define the V-cell applying the momentum exchange term
              !
              If(VCell%Cell_Type(i,j,k)==1) then ! cut cell
                !      
                If(j<Jmax) then
                  !      
                  If(Pcell%Cell_Type(i,j,k)==2.or.Pcell%Cell_Type(i,j+1,k)==2) &
                                                                           then
                    VCell%MoExCell(i,j,k) = 1
                    !
                  End if
                  !
                Else
                  !
                  If(PCell%Cell_Type(i,j,k)==2) VCell%MoExCell(i,j,k) = 1
                  !
                End if
                !
              End if
              !
              ! Define the W-cell applying the momentum exchange term
              !
              If(WCell%Cell_Type(i,j,k)==1) then ! cut cell
                !       
                If(k<Kmax) then
                  !      
                  If(PCell%Cell_Type(i,j,k)==2.or.PCell%Cell_Type(i,j,k+1)==2) &
                                                                           then
                    WCell%MoExCell(i,j,k) = 1
                  End if
                  !
                Else
                  !
                  If(PCell%Cell_Type(i,j,k)==2) WCell%MoExCell(i,j,k) = 1
                  !
                End if
                !
              End if
              !
            End do
          End do
        End do
        !
      End subroutine DefineMomentumExchangeCell
      !
      ! numbering the pressure cell for poisson equation
      !
      Subroutine NumberExternalCell(iu,iv,iw, TCell)
        !      
        Implicit none
        !
        Integer(kind=it4b),intent(in):: iu,iv,iw
        Type(Cell),intent(inout):: TCell

        Integer(kind=it4b):: i,j,k
        Integer(kind=it4b):: ctr
        !
        ctr = 0
        !
        Do i = 1,Imax-iu
          Do j = 1,Jmax-iv
            Do k = 1,Kmax-iw
              !
              If(TCell%Cell_Type(i,j,k)/=2) then
                !      
                TCell%Posnu(i,j,k) = ctr
                ctr = ctr+1
                !
              Else
                !      
                TCell%Posnu(i,j,k) = -1
                !
              End if
              !
            End do
          End do
        End do
        !
        TCell%ExtCell = ctr-1
        !
      End subroutine NumberExternalCell
      !
      Subroutine NewCellFace(PGrid,UGrid,VGrid,WGrid, PCell,UCell,VCell,WCell)
        !      
        Implicit none
        !
        Type(Grid),intent(in):: PGrid,UGrid,VGrid,WGrid
        Type(Cell),intent(inout):: PCell,UCell,VCell,WCell

        Integer(kind=it4b):: i,j,k
        !
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 1,Kmax
              !
              If(Pcell%EEArea(i,j,k)>0.d0) then
                UCell%Cell_Cent(i,j,k,1) = 0.d0
                UCell%Cell_Cent(i,j,k,2) = PCell%FCE(i,j,k,2)
                UCell%Cell_Cent(i,j,k,3) = PCell%FCE(i,j,k,3)
              End if
              !
              If(PCell%NEArea(i,j,k)>0.d0) then
                VCell%Cell_Cent(i,j,k,1) = PCell%FCN(i,j,k,1)
                VCell%Cell_Cent(i,j,k,2) = 0.0d0
                VCell%Cell_Cent(i,j,k,3) = PCell%FCN(i,j,k,3)
              End if
              !
              If(PCell%TEArea(i,j,k)>0.d0) then
                WCell%Cell_Cent(i,j,k,1) = PCell%FCT(i,j,k,1)
                WCell%Cell_Cent(i,j,k,2) = PCell%FCT(i,j,k,2)
                WCell%Cell_Cent(i,j,k,3) = 0.d0
              End if
              !
              ! find the cell center for pressure cell
              !
              If(PCell%Cell_Type(i,j,k)==1) then

              End if
              !
            End do
          End do
        End do
        !
        Do i = 1,Imax-1
          Do j = 1,Jmax-1
            Do k = 1,Kmax-1
              !
              If(UCell%MoExCell(i,j,k)==1) then ! for very small UCell which connects to only one PCell
                Call CellLinking(i,j,k,UGrid, UCell)
              End if
              !
              If(VCell%MoExCell(i,j,k)==1) then ! for very small VCell which connects to only one PCell
                Call CellLinking(i,j,k,VGrid, VCell)
              End if
              !
              If(WCell%MoExCell(i,j,k)==1) then ! for very small WCell Which connects to only one PCell
                Call CellLinking(i,j,k,WGrid, WCell)
              End if
              !
            End do
          End do
        End do
        !
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 1,Kmax
              !
              UCell%delh(i,j,k)=dabs(UCell%Cell_Cent(i,j,k,1)*UCell%nx(i,j,k)+ &
                           UCell%Cell_Cent(i,j,k,2)*UCell%ny(i,j,k)+           &
                      UCell%Cell_Cent(i,j,k,3)*UCell%nz(i,j,k)+UCell%phi(i,j,k))
              VCell%delh(i,j,k)=dabs(VCell%Cell_Cent(i,j,k,1)*VCell%nx(i,j,k)+ &
                           VCell%Cell_Cent(i,j,k,2)*VCell%ny(i,j,k)+           &
                      VCell%Cell_Cent(i,j,k,3)*VCell%nz(i,j,k)+VCell%phi(i,j,k))
              WCell%delh(i,j,k)=dabs(WCell%Cell_Cent(i,j,k,1)*WCell%nx(i,j,k)+ &
                           WCell%Cell_Cent(i,j,k,2)*WCell%ny(i,j,k)+           &
                      WCell%Cell_Cent(i,j,k,3)*WCell%nz(i,j,k)+WCell%phi(i,j,k))
              !
            End do
          End do
        End do
        !
        ! Calculate all information needed for interpolation procedure
        !
        Call EastFaceInterpolationInf(1,UGrid,PGrid, Ucell)
        Call NorthFaceInterpolationInf(0,UGrid,VGrid, Ucell)
        Call TopFaceInterpolationInf(0,UGrid,WGrid, UCell)

        Call EastFaceInterpolationInf(0,VGrid,UGrid, VCell)
        Call NorthFaceInterpolationInf(1,VGrid,PGrid, Vcell)
        Call TopFaceInterpolationInf(0,VGrid,WGrid, Vcell)

        Call EastFaceInterpolationInf(0,WGrid,UGrid, WCell)
        Call NorthFaceInterpolationInf(0,WGrid,VGrid, WCell)
        Call TopFaceInterpolationInf(1,WGrid,PGrid, WCell)
        !
      End subroutine NewCellFace
      !
      Subroutine FaceGeoCal(d1,d2,n1,n2,n3,d3,phi, Area,Cent)
        !      
        Implicit none
        !
        Real(kind=dp),intent(in):: d1,d2,n1,n2,n3,d3,phi
        Real(kind=dp),intent(inout):: Area,Cent(3)

        Integer(kind=it4b):: temp,k
        Real(kind=dp):: d11,d22
        Real(kind=dp):: dpt(4),vol
        Real(kind=dp),dimension(:,:):: node(6,2)
        !
        d11 = d1/2.d0
        d22 = d2/2.d0
        !
        temp = 1
        !
        node = 0.d0
        !
        dpt(1) =  d11*n1-d22*n2+d3*n3+phi
        dpt(2) = -d11*n1-d22*n2+d3*n3+phi
        dpt(3) = -d11*n1+d22*n2+d3*n3+phi
        dpt(4) =  d11*n1+d22*n2+d3*n3+phi
        !
        If(dpt(1)>0.d0) then
          node(temp,1) = d11
          node(temp,2) = -d22
          temp = temp+1
        End if
        !
        If(dpt(1)*dpt(2)<0.d0) then
          node(temp,1) = (d22*n2-d3*n3-phi)/n1
          node(temp,2) = -d22
          temp = temp+1
        End if
        !
        If(dpt(2)>0.d0) then
          node(temp,1) = -d11
          node(temp,2) = -d22
          temp = temp+1
        End if
        !
        If(dpt(2)*dpt(3)<0.d0) then
          node(temp,1) = -d11
          node(temp,2) = (d11*n1-d3*n3-phi)/n2
          temp = temp+1
        End if
        !
        If(dpt(3)>0.d0) then
          node(temp,1) = -d11
          node(temp,2) = d22
          temp = temp+1
        End if
        !
        If(dpt(3)*dpt(4)<0.d0) then
          node(temp,1) = (-d22*n2-d3*n3-phi)/n1
          node(temp,2) = d22
          temp = temp+1
        End if
        !
        If(dpt(4)>0.d0) then
          node(temp,1) = d11
          node(temp,2) = d22
          temp = temp+1
        End if
        !
        If(dpt(4)*dpt(1)<0.d0) then
          node(temp,1) = d11
          node(temp,2) = (-d11*n1-d3*n3-phi)/n2
          temp = temp+1
        End if
        !
        node(temp,1) = node(1,1)
        node(temp,2) = node(1,2)
        !
        vol = 0.d0
        !
        Do k = 1,temp-1
          vol = vol+0.5d0*(node(k,1)*node(k+1,2)-node(k+1,1)*node(k,2))
        End do
        !
        Area = dabs(vol/(d1*d2))
        !
        Cent(:) = 0.d0
        !
        If(Area>0.d0) then
          !      
          Do k = 1,temp-1
            Cent(1) = Cent(1)+1.d0/6.d0/vol*(node(k,1)+node(k+1,1))*           &
                               (node(k,1)*node(k+1,2)-node(k+1,1)*node(k,2))
            Cent(2) = Cent(2)+1.d0/6.d0/vol*(node(k,2)+node(k+1,2))*           &
                               (node(k,1)*node(k+1,2)-node(k+1,1)*node(k,2))
          End do
          !
        End if
        !
        Cent(3) = d3
        !
      End subroutine FaceGeoCal
      !
      Subroutine EastFaceInterpolationInf(iu,TGrid,BGrid, TCell)
        !      
        Implicit none
        !
        Integer(kind=it4b),intent(in):: iu
        Type(Grid),intent(in):: TGrid,BGrid
        Type(Cell),intent(inout):: TCell

        Integer(kind=it4b):: i,j,k
        Real(kind=dp):: Sy,Sz,xf,yf,zf,nxf,nyf,nzf
        Real(kind=dp):: delh,delh1,delh2,delhec,delhec1,delhec2
        !
        Do j = 1,Jmax
          Do k = 1,Kmax
            Do i = 1,Imax-1
              !
              if(TCell%Cell_Type(i,j,k)/=2.and.TCell%Cell_Type(i+1,j,k)/=2) then
                !      
                TCell%SxE(i,j,k) = TCell%Cell_Cent(i+1,j,k,1)+BGrid%dx(i+iu,j,k)-&
                                   TCell%Cell_Cent(i,j,k,1)
                Sy = TCell%Cell_Cent(i+1,j,k,2)-TCell%Cell_Cent(i,j,k,2)
                Sz = TCell%Cell_Cent(i+1,j,k,3)-TCell%Cell_Cent(i,j,k,3)
                TCell%EtaE(i,j,k) = dabs((TCell%FCE(i,j,k,1)-                  &
                                TCell%Cell_Cent(i,j,k,1))/TCell%SxE(i,j,k))
                !        
                If(TCell%EtaE(i,j,k)>1.d0) TCell%EtaE(i,j,k) = 0.5d0
                !
                xf = (1.d0-TCell%EtaE(i,j,k))*TCell%Cell_Cent(i,j,k,1)+        &
                         TCell%EtaE(i,j,k)*(TCell%Cell_Cent(i+1,j,k,1)+        &
                                                          BGrid%dx(i+iu,j,k))
                yf = (1.d0-TCell%EtaE(i,j,k))*TCell%Cell_Cent(i,j,k,2)+        &
                         TCell%EtaE(i,j,k)*TCell%Cell_Cent(i+1,j,k,2)
                zf = (1.d0-TCell%EtaE(i,j,k))*TCell%Cell_Cent(i,j,k,3)+        &
                         TCell%EtaE(i,j,k)*TCell%Cell_Cent(i+1,j,k,3)
                ! 
                If(dabs(Sy)>1.d-4*TGrid%dy(i,j,k).or.                          &
                                       dabs(Sz)>1.d-4*TGrid%dz(i,j,k)) then
                  !                     
                  nxf = 0.5d0*(TCell%nx(i,j,k)+TCell%nx(i+1,j,k))
                  nyf = 0.5d0*(TCell%ny(i,j,k)+TCell%ny(i+1,j,k))
                  nzf = 0.5d0*(TCell%nz(i,j,k)+TCell%nz(i+1,j,k))
                  delh1 = xf*TCell%nx(i,j,k)+yf*TCell%ny(i,j,k)+               &
                                            zf*TCell%nz(i,j,k)+TCell%phi(i,j,k)
                  delh2 = (xf-BGrid%dx(i+iu,j,k))*TCell%nx(i+1,j,k)+           &
                         yf*TCell%ny(i+1,j,k)+zf*TCell%nz(i+1,j,k)+            &
                                                             TCell%phi(i+1,j,k)
                  delh = 0.5d0*(delh1+delh2)

                  delhec1 = TCell%FCE(i,j,k,1)*TCell%nx(i,j,k)+TCell%FCE(i,j,k,2)&
                         *TCell%ny(i,j,k)+TCell%FCE(i,j,k,3)*TCell%nz(i,j,k)+  &
                                                             TCell%phi(i,j,k)
                  delhec2=(TCell%FCE(i,j,k,1)-BGrid%dx(i+iu,j,k))*TCell%nx(i+1,j,k)&
                       +TCell%FCE(i,j,k,2)*TCell%ny(i+1,j,k)+TCell%FCE(i,j,k,3)&
                                       *TCell%nz(i+1,j,k)+TCell%phi(i+1,j,k)
                  delhec = 0.5d0*(delhec1+delhec2)
                  !
                  if(delh<epsi.and.TCell%EEArea(i,j,k)>0.d0) then
                    delh=delhec
                  end if
                  !
                  TCell%AlE(i,j,k) = dabs(delhec)/dabs(delh)
                  TCell%DAlE(i,j,k) = (Sy*nyf+Sz*nzf)/(TCell%SxE(i,j,k)*dabs(delh))
                  !
                Else
                  !      
                  TCell%DAlE(i,j,k) = 0.d0
                  TCell%AlE(i,j,k) = 1.d0
                  !
                End if
                !
              else
                !      
                TCell%AlE(i,j,k)=1.d0
                TCell%DalE(i,j,k)=0.d0
                TCell%SxE(i,j,k)=BGrid%dx(i+iu,j,k)
                TCell%EtaE(i,j,k)=0.5d0
                !
              end if  
              !
            End do
            !
            TCell%SxE(i,j,k) = TCell%SxE(i-1,j,k)
            TCell%EtaE(i,j,k) = TCell%EtaE(i-1,j,k)
            TCell%AlE(i,j,k) = TCell%AlE(i-1,j,k)
            TCell%DAlE(i,j,k) = TCell%DAlE(i-1,j,k)
            !
          End do
        End do
        !
        TCell%SxE(0,:,:) = TCell%SxE(1,:,:)
        TCell%EtaE(0,:,:) = TCell%EtaE(1,:,:)
        TCell%AlE(0,:,:) = TCell%AlE(1,:,:)
        TCell%DAlE(0,:,:) = Tcell%DAlE(1,:,:)
        !
      End Subroutine EastFaceInterpolationInf
      !
      Subroutine NorthFaceInterpolationInf(iv,TGrid,BGrid, TCell)
        !      
        Implicit none
        !
        Integer(kind=it4b),intent(in):: iv
        Type(Grid),intent(in):: TGrid,BGrid
        Type(Cell),intent(inout):: TCell

        Integer:: i,j,k
        Real(kind=dp):: Sx,Sz,xf,yf,zf,nxf,nyf,nzf
        Real(kind=dp):: delh,delh1,delh2,delhec,delhec1,delhec2
        !
        do i = 1,Imax
          do k = 1,Kmax
            do j = 1,Jmax-1
              !
              if(TCell%Cell_Type(i,j,k)/=2.and.TCell%Cell_Type(i,j+1,k)/=2) then
                !      
                TCell%SyN(i,j,k)=TCell%Cell_Cent(i,j+1,k,2)+BGrid%dy(i,j+iv,k)-&
                                 TCell%Cell_Cent(i,j,k,2)
                Sx = TCell%Cell_Cent(i,j+1,k,1)-TCell%Cell_Cent(i,j,k,1)
                Sz = TCell%Cell_Cent(i,j+1,k,3)-TCell%Cell_Cent(i,j,k,3)
                TCell%EtaN(i,j,k) = dabs((TCell%FCN(i,j,k,2)-                  &
                                    TCell%Cell_Cent(i,j,k,2))/TCell%SyN(i,j,k))
                !            
                if(TCell%EtaN(i,j,k)>1.d0) TCell%EtaN(i,j,k) = 0.5d0
                !
                xf=(1.d0-TCell%EtaN(i,j,k))*TCell%Cell_Cent(i,j,k,1)+          &
                         TCell%EtaN(i,j,k)*TCell%Cell_Cent(i,j+1,k,1)
                yf=(1.d0-TCell%EtaN(i,j,k))*TCell%Cell_Cent(i,j,k,2)+          &
                         TCell%EtaN(i,j,k)*(TCell%Cell_Cent(i,j+1,k,2)+        &
                                                           BGrid%dy(i,j+iv,k))
                zf=(1.d0-TCell%EtaN(i,j,k))*TCell%Cell_Cent(i,j,k,3)+          &
                         TCell%EtaN(i,j,k)*TCell%Cell_Cent(i,j+1,k,3)
                ! 
                if(dabs(Sx)>=1.d-4*TGrid%dx(i,j,k).or.dabs(Sz)>=1.d-4*         &
                                                         TGrid%dz(i,j,k)) then
                  !                                       
                  nxf = 0.5d0*(TCell%nx(i,j,k)+TCell%nx(i,j+1,k))
                  nyf = 0.5d0*(TCell%ny(i,j,k)+TCell%ny(i,j+1,k))
                  nzf = 0.5d0*(TCell%nz(i,j,k)+TCell%nz(i,j+1,k))
                  delh1 = xf*TCell%nx(i,j,k)+yf*TCell%ny(i,j,k)+               &
                        zf*TCell%nz(i,j,k)+TCell%phi(i,j,k)
                  delh2 = xf*TCell%nx(i,j+1,k)+(yf-BGrid%dy(i,j+iv,k))*        &
                           TCell%ny(i,j+1,k)+zf*TCell%nz(i,j+1,k)+             &
                                                            TCell%phi(i,j+1,k)
                  delh = 0.5d0*(delh1+delh2)
                  delhec1=TCell%FCN(i,j,k,1)*TCell%nx(i,j,k)+TCell%FCN(i,j,k,2)&
                         *TCell%ny(i,j,k)+TCell%nz(i,j,k)*TCell%FCN(i,j,k,3)+  &
                          TCell%phi(i,j,k)
                  delhec2=TCell%FCN(i,j,k,1)*TCell%nx(i,j+1,k)+(TCell%FCN(i,j,k,2)&
                       -BGrid%dy(i,j+1,k))*TCell%ny(i,j+1,k)+TCell%nz(i,j+1,k)*&
                                         TCell%FCN(i,j,k,3)+TCell%phi(i,j+1,k)
                  delhec = 0.5d0*(delhec1+delhec2)
                  !
                  if(delh<epsi.and.TCell%NEArea(i,j,k)>0.d0) then
                    delh = delhec
                  end if
                  !
                  TCell%AlN(i,j,k) = dabs(delhec)/dabs(delh)
                  TCell%DAlN(i,j,k) = (Sx*nxf+Sz*nzf)/(TCell%SyN(i,j,k)*dabs(delh))
                  !
                else
                  !
                  TCell%DAlN(i,j,k) = 0.d0
                  TCell%AlN(i,j,k) = 1.d0
                  !
                end if
                !
              else
                !
                TCell%EtaN(i,j,k)=0.5d0
                TCell%AlN(i,j,k)=1.d0
                TCell%DAlN(i,j,k)=0.d0
                TCell%SyN(i,j,k)=BGrid%dy(i,j+iv,k)
                !
              end if
              !
            end do
            !
            TCell%SyN(i,j,k) = TCell%SyN(i,j-1,k)
            TCell%EtaN(i,j,k) = TCell%EtaN(i,j-1,k)
            TCell%AlN(i,j,k) = TCell%AlN(i,j-1,k)
            TCell%DAlN(i,j,k) = TCell%DAlN(i,j-1,k)
            !
          end do
        end do
        !
        TCell%SyN(:,0,:) = TCell%SyN(:,1,:)
        TCell%EtaN(:,0,:) = TCell%EtaN(:,1,:)
        TCell%AlN(:,0,:) = TCell%AlN(:,1,:)
        TCell%DAlN(:,0,:) = TCell%DAlN(:,1,:)
        !
      End subroutine NorthFaceInterpolationInf
      !
      Subroutine TopFaceInterpolationInf(iw,TGrid,BGrid, TCell)
        !      
        Implicit none
        !
        Integer(kind=it4b),intent(in):: iw
        Type(Grid),intent(in):: TGrid,BGrid
        Type(Cell),intent(inout):: TCell

        Integer:: i,j,k
        Real(kind=dp):: Sx,Sy,xf,yf,zf,nxf,nyf,nzf
        Real(kind=dp):: delh,delh1,delh2,delhec,delhec1,delhec2
        !
        Do i = 1,Imax
          Do j = 1,Jmax
            Do k = 1,Kmax-1
              !
              if(TCell%Cell_Type(i,j,k)/=2.and.TCell%Cell_Type(i,j,k+1)/=2) then
                !      
                Sx = TCell%Cell_Cent(i,j,k+1,1)-TCell%Cell_Cent(i,j,k,1)
                Sy = TCell%Cell_Cent(i,j,k+1,2)-TCell%Cell_Cent(i,j,k,2)
                TCell%SzT(i,j,k)=TCell%Cell_Cent(i,j,k+1,3)+BGrid%dz(i,j,k+iw)-&
                                 TCell%Cell_Cent(i,j,k,3)
                TCell%EtaT(i,j,k)=dabs((TCell%FCT(i,j,k,3)-                    &
                                 TCell%Cell_Cent(i,j,k,3))/TCell%SzT(i,j,k))
                !         
                if(TCell%EtaT(i,j,k)>1.d0) TCell%EtaT(i,j,k)=0.5d0
                !
                xf=(1.d0-TCell%EtaT(i,j,k))*TCell%Cell_Cent(i,j,k,1)+          &
                         TCell%EtaT(i,j,k)*TCell%Cell_Cent(i,j,k+1,1)
                yf=(1.d0-TCell%EtaT(i,j,k))*TCell%Cell_Cent(i,j,k,2)+          &
                         TCell%EtaT(i,j,k)*TCell%Cell_Cent(i,j,k+1,2)
                zf=(1.d0-TCell%EtaT(i,j,k))*TCell%Cell_Cent(i,j,k,3)+          &
                         TCell%EtaT(i,j,k)*(TCell%Cell_Cent(i,j,k+1,3)+        &
                                                  BGrid%dz(i,j,k+iw))
                !                          
                if(dabs(Sx)>=1.d-4*TGrid%dx(i,j,k).or.dabs(Sy)>=1.d-4*         &
                                                         TGrid%dy(i,j,k)) then
                  !                                       
                  nxf=0.5d0*(TCell%nx(i,j,k)+TCell%nx(i,j,k+1))
                  nyf=0.5d0*(TCell%ny(i,j,k)+TCell%ny(i,j,k+1))
                  nzf=0.5d0*(TCell%nz(i,j,k)+TCell%nz(i,j,k+1))
                  delh1=xf*TCell%nx(i,j,k)+yf*TCell%ny(i,j,k)+                 &
                        zf*TCell%nz(i,j,k)+TCell%phi(i,j,k)
                  delh2=xf*TCell%nx(i,j,k+1)+yf*TCell%ny(i,j,k+1)+             &
                       (zf-BGrid%dz(i,j,k+iw))*TCell%nz(i,j,k+1)+              &
                                                            TCell%phi(i,j,k+1)
                  delh=0.5d0*(delh1+delh2)
                  delhec1=TCell%FCT(i,j,k,1)*TCell%nx(i,j,k)+TCell%FCT(i,j,k,2)&
                         *TCell%ny(i,j,k)+TCell%nz(i,j,k)*TCell%FCT(i,j,k,3)+  &
                          TCell%phi(i,j,k)
                  delhec2=TCell%FCT(i,j,k,1)*TCell%nx(i,j,k+1)+TCell%FCT(i,j,k,2)&
                       *TCell%ny(i,j,k+1)+TCell%nz(i,j,k+1)*(TCell%FCT(i,j,k,3)&
                       -BGrid%dz(i,j,k+iw))+TCell%phi(i,j,k+1)
                  delhec = 0.5d0*(delhec1+delhec2)
                  !  
                  if(delh<epsi.and.TCell%TEArea(i,j,k)>0.d0) then
                    delh = delhec
                  end if
                  !
                  TCell%AlT(i,j,k)=dabs(delhec)/dabs(delh)
                  TCell%DAlT(i,j,k)=(Sx*nxf+Sy*nyf)/(TCell%SzT(i,j,k)*dabs(delh))
                  !
                else
                  !
                  TCell%DAlT(i,j,k) = 0.d0
                  TCell%AlT(i,j,k) = 1.d0
                  !
                end if
                !
              else
                !
                TCell%EtaT(i,j,k)=0.5d0
                TCell%AlT(i,j,k)=1.d0
                TCell%DAlT(i,j,k)=0.d0
                TCell%SzT(i,j,k)=BGrid%dz(i,j,k+iw)
                !
              end if  
              !
            end do
            !
            TCell%SzT(i,j,k) = TCell%SzT(i,j,k-1)
            TCell%EtaT(i,j,k) = TCell%EtaT(i,j,k-1)
            TCell%DAlT(i,j,k) = TCell%DAlT(i,j,k-1)
            TCell%AlT(i,j,k) = TCell%AlT(i,j,k-1)
            !
          end do
        end do
        !
        TCell%SzT(:,:,0) = TCell%SzT(:,:,1)
        TCell%EtaT(:,:,0) = TCell%EtaT(:,:,1)
        TCell%AlT(:,:,0) = TCell%AlT(:,:,1)
        TCell%DAlT(:,:,0) = TCell%DAlT(:,:,1)
        !
      End Subroutine TopFaceInterpolationInf
      !
      Subroutine CellLinking(i,j,k,TGrid, TCell)
        !      
        Implicit none
        !
        Integer(kind=it4b),intent(in):: i,j,k
        Type(Grid),intent(in):: TGrid
        Type(Cell),intent(inout):: TCell

        Real(kind=dp):: nmax,tol
        !
        tol = 1.d-8
        !
        nmax = dmax1(dabs(TCell%nx(i,j,k)),dabs(TCell%ny(i,j,k)),              &
                     dabs(TCell%nz(i,j,k)))
        !
        If(dabs(dabs(TCell%nx(i,j,k))-nmax)<1.d-10) then
          !      
          If(TCell%nx(i,j,k)<0.d0) then
            !      
            TCell%MsCe(i,j,k,1) = i-1
            TCell%Cell_Cent(i,j,k,1) = TCell%Cell_Cent(i-1,j,k,1)-             &
                                       TGrid%dx(i-1,j,k)+tol*TGrid%dx(i-1,j,k)
            TCell%Cell_Cent(i,j,k,2) = TCell%Cell_Cent(i-1,j,k,2)
            TCell%Cell_Cent(i,j,k,3) = TCell%Cell_Cent(i-1,j,k,3)
            !
          Else
            !
            TCell%MsCe(i,j,k,1) = i+1
            TCell%Cell_Cent(i,j,k,1) = TCell%Cell_Cent(i+1,j,k,1)+             &
                                       TGrid%dx(i+1,j,k)-tol*TGrid%dx(i+1,j,k)
            TCell%Cell_Cent(i,j,k,2) = TCell%Cell_Cent(i+1,j,k,2)
            TCell%Cell_Cent(i,j,k,3) = TCell%Cell_Cent(i+1,j,k,3)
            !
          End if
          !
          TCell%MsCe(i,j,k,2) = j
          TCell%MsCe(i,j,k,3) = k
          !
        ElseIf(dabs(dabs(TCell%ny(i,j,k))-nmax)<1.d-10) then
          !      
          If(TCell%ny(i,j,k)<0.d0) then
            !      
            TCell%MsCe(i,j,k,2) = j-1
            TCell%Cell_Cent(i,j,k,2) = TCell%Cell_Cent(i,j-1,k,2)-             &
                                       TGrid%dy(i,j-1,k)+tol*TGrid%dy(i,j-1,k)
            TCell%Cell_Cent(i,j,k,1) = TCell%Cell_Cent(i,j-1,k,1)
            TCell%Cell_Cent(i,j,k,3) = TCell%Cell_Cent(i,j-1,k,3)
            !
          Else
            !      
            TCell%MsCe(i,j,k,2) = j+1
            TCell%Cell_Cent(i,j,k,2) = TCell%Cell_Cent(i,j+1,k,2)+             &
                                       TGrid%dy(i,j+1,k)-tol*TGrid%dy(i,j+1,k)
            TCell%Cell_Cent(i,j,k,1) = TCell%Cell_Cent(i,j+1,k,1)
            TCell%Cell_Cent(i,j,k,3) = TCell%Cell_Cent(i,j+1,k,3)
            !
          End if
          !
          TCell%MsCe(i,j,k,1) = i
          TCell%MsCe(i,j,k,3) = k
          !
        ElseIf(dabs(dabs(TCell%nz(i,j,k))-nmax)<1.d-10) then
          !      
          If(TCell%nz(i,j,k)<0.d0) then
            !     
            TCell%MsCe(i,j,k,3) = k-1
            TCell%Cell_Cent(i,j,k,3) = TCell%Cell_Cent(i,j,k-1,3)-             &
                                       TGrid%dz(i,j,k-1)+tol*TGrid%dz(i,j,k-1)
            TCell%Cell_Cent(i,j,k,1) = TCell%Cell_Cent(i,j,k-1,1)
            TCell%Cell_Cent(i,j,k,2) = TCell%Cell_Cent(i,j,k-1,2)
            !
          Else
            !
            TCell%MsCe(i,j,k,3) = k+1
            TCell%Cell_Cent(i,j,k,3) = TCell%Cell_Cent(i,j,k+1,3)+             &
                                       TGrid%dz(i,j,k+1)-tol*TGrid%dz(i,j,k+1)
            TCell%Cell_Cent(i,j,k,1) = TCell%Cell_Cent(i,j,k+1,1)
            TCell%Cell_Cent(i,j,k,2) = TCell%Cell_Cent(i,j,k+1,2)
            !
          End if
          !
          TCell%MsCe(i,j,k,1) = i
          TCell%MsCe(i,j,k,2) = j
          !
        End if
        !
      End subroutine CellLinking
      !
End Module CutCell
