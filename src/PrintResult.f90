Module PrintResult
    USE PrecisionVar
    USE Mesh
    USE StateVariables,ONLY : ight,jght,kght,variables,pi,Ktw,rey,xc,yc
    USE Clsvof
    USE Cutcell
    USE VTK
    USE VTR

    Implicit none
    Private

    character(len=1), parameter :: newline=achar(10)

    Public:: PrintResultTecplotPCent,PrintResultTecplotUCent,                  &
             PrintResultTecplotVCent,PrintResultTecplotWCent,                  &
             PrintResultTecplotPCentXY,PrintResultVTK,PrintResultVTR3D,	       &
             PrintResultTecplotPCentXZ	

    Interface PrintResultTecplotPCent
       Module procedure PrintResultTecplotPCent
    End Interface PrintResultTecplotPCent

    Interface PrintResultTecplotUCent
       Module procedure PrintResultTecplotUCent
    End Interface PrintResultTecplotUCent

    Interface PrintResultTecplotVCent
       Module procedure PrintResultTecplotVCent
    End Interface PrintResultTecplotVCent

    Interface PrintResultTecplotWCent
       Module procedure PrintResultTecplotWCent
    End Interface PrintResultTecplotWCent

    Interface PrintResultTecplotPCentXY
       Module procedure PrintResultTecplotPCentXY
    End Interface PrintResultTecplotPCentXY

    Interface PrintResultVTK
      module procedure PrintResultVTK
    End interface PrintResultVTK

    Interface PrintResultVTR3D
      module procedure PrintResultVTR3D
    End interface PrintResultVTR3D

    Contains

    Subroutine PrintResultTecplotPCent(TGrid,TVar,TCell,iter)
      Implicit none
      type(Grid),intent(in):: TGrid
      type(Variables),intent(inout):: TVar
      type(Cell),intent(in):: TCell
      Integer(kind=it8b),intent(in):: iter
      Integer(kind=it4b) i,j,k
      Real(kind=dp),dimension(:,:,:),allocatable:: p
      Character(15) curd

      Allocate(p(Imax,Jmax,Kmax))
      Do i = 1,Imax
        Do j = 1,Jmax
          Do k = 1,Kmax
            p(i,j,k) = TVar%p(i,j,k)
          End do
        End do
      End do
      Open(unit=11,file='text.txt',action='write')
      Write(11,*) iter
      Close(11)
      Open(unit=12,file='text.txt',action='read')
      Read(12,*) curd
      Close(12,status='delete')
      Open(unit=5,file=trim(curd)//'Pressure.dat',action='write')
      Write(5,*) 'title = "flow"'
      Write(5,*) 'variables ="x","y","z","p","u","v","w","phi","vof"'!//',"Pres","Mres"'
      Write(5,1110) Imax,Jmax,Kmax
      Write(5,"(f18.10)")(((TGrid%x(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TGrid%y(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TGrid%z(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((p(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      write(5,"(f18.10)")(((0.5d0*(TVar%u(i-1,j,k)+TVar%u(i,j,k)),i=1,Imax),   &
                                                            j=1,Jmax),k=1,Kmax)
      write(5,"(f18.10)")(((0.5d0*(TVar%v(i,j-1,k)+TVar%v(i,j,k)),i=1,Imax),   &
                                                            j=1,Jmax),k=1,Kmax)
      write(5,"(f18.10)")(((0.5d0*(TVar%w(i,j,k-1)+TVar%w(i,j,k)),i=1,Imax),   &
                                                            j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%phiL(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%vofL(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
  !    Write(5,"(f18.10)")(((TVar%pres(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
  !    Write(5,"(f18.10)")(((TVar%mres(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Close(5)
 1110 format(1x,'zone i=',i5,',j=',i5,',k=',i5,'f=block')
      Deallocate(p)
    End Subroutine PrintResultTecplotPCent

    Subroutine PrintResultVTR3D(TGrid,TVar,TCell,itt)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN)         :: TGrid
      TYPE(Variables),INTENT(IN)    :: TVar
      TYPE(Cell),INTENT(IN)         :: TCell
      INTEGER(kind=it8b),INTENT(IN) :: itt
      TYPE(VTR_file_handle)         :: fd
      CHARACTER(len=70)              :: dir
      dir = trim("/home/sontd/code/CutCell3DGFMCLSVOF/Result/")

      call VTR_open_file(Prefix="FlowField",dir=dir, itera=itt,FD=fd)
    ! use keyword argument due to huge number of optional dummy argument
    ! so we need keyword argument to specify the location of actual argument
      call VTR_write_mesh(fd,TGrid%x(:,1,1),TGrid%y(1,:,1),TGrid%z(1,1,:))
      call VTR_write_var(fd,"Velocity",TVar%u(1:IMax,1:JMax,1:KMax),		&
                                       TVar%v(1:IMax,1:JMax,1:Kmax),		&
                                       TVar%w(1:Imax,1:Jmax,1:Kmax))
      call VTR_write_var(fd,"Pressure",TVar%p(1:IMax,1:JMax,1:KMax))
      call VTR_write_var(fd,"FluidLvs",TCell%phi(1:Imax,1:Jmax,1:Kmax))
      call VTR_write_var(fd,"LiquidLvs",TCell%phiL(1:Imax,1:Jmax,1:Kmax))
      call VTR_Write_var(fd,"FluidVof",TCell%vof(1:Imax,1:Jmax,1:Kmax))
      call VTR_Write_var(fd,"LiquidVof",TCell%vofL(1:Imax,1:Jmax,1:Kmax))
      call VTR_Write_var(fd,"Mass error",Tvar%mres(1:Imax,1:Jmax,1:Kmax))
      call VTR_close_file(fd)
    end subroutine PrintResultVTR3D

    Subroutine PrintResultTecplotPCentXY(TGrid,TVar,TCell,iter)
      Implicit none
      type(Grid),intent(in):: TGrid
      type(Variables),intent(inout):: TVar
      type(Cell),intent(in):: TCell
      Integer(kind=it8b),intent(in):: iter
      Integer(kind=it4b) i,j,k,kp
      Real(kind=dp),dimension(:,:,:),allocatable:: p
      Character(15) curd
      Allocate(p(Imax,Jmax,Kmax))
      kp = int(Kmax/2)
      kp = kp
      Do i = 1,Imax
        Do j = 1,Jmax
          Do k = 1,Kmax
            p(i,j,k) = TVar%p(i,j,k)
          End do
        End do
      End do
      Open(unit=11,file='text.txt',action='write')
      Write(11,*) iter
      Close(11)
      Open(unit=12,file='text.txt',action='read')
      Read(12,*) curd
      Close(12,status='delete')
      Open(unit=5,file=trim(curd)//'PresXY.dat',action='write')
      Write(5,*) 'title = "flow"'
      Write(5,*) 'variables ="x","y","p","u","v","phi","vof","Pres","Mres"'
      Write(5,1110) Imax,Jmax
      Write(5,"(f18.10)")((TGrid%x(i,j,kp),i=1,Imax),j=1,Jmax)
      Write(5,"(f18.10)")((TGrid%y(i,j,kp),i=1,Imax),j=1,Jmax)
      Write(5,"(f18.10)")((p(i,j,kp),i=1,Imax),j=1,Jmax)
      write(5,"(f18.10)")((0.5d0*(TVar%u(i-1,j,kp)+TVar%u(i,j,kp)),i=1,Imax),  &
                                                            j=1,Jmax)
      write(5,"(f18.10)")((0.5d0*(TVar%v(i,j-1,kp)+TVar%v(i,j,kp)),i=1,Imax),  &
                                                            j=1,Jmax)
      Write(5,"(f18.10)")((TCell%phi(i,j,kp),i=1,Imax),j=1,Jmax)
      Write(5,"(f18.10)")((TCell%vof(i,j,kp),i=1,Imax),j=1,Jmax)
      Write(5,"(f18.10)")((TVar%pres(i,j,kp),i=1,Imax),j=1,Jmax)
      Write(5,"(f18.10)")((TVar%mres(i,j,kp),i=1,Imax),j=1,Jmax)
      Close(5)
 1110 format(1x,'zone i=',i5,',j=',i5,'f=block')
      Deallocate(p)
    End Subroutine PrintResultTecplotPCentXY
    
    Subroutine PrintResultTecplotPCentXZ(TGrid,TVar,TCell,iter)
      Implicit none
      type(Grid),intent(in):: TGrid
      type(Variables),intent(inout):: TVar
      type(Cell),intent(in):: TCell
      Integer(kind=it8b),intent(in):: iter
      Integer(kind=it4b) i,j,k,jp
      Real(kind=dp),dimension(:,:,:),allocatable:: p
      Character(15) curd
      Allocate(p(Imax,Jmax,Kmax))
      jp = int(Jmax/2)
      Do i = 1,Imax
        Do j = 1,Jmax
          Do k = 1,Kmax
            p(i,j,k) = TVar%p(i,j,k)
          End do
        End do
      End do
      Open(unit=11,file='text.txt',action='write')
      Write(11,*) iter
      Close(11)
      Open(unit=12,file='text.txt',action='read')
      Read(12,*) curd
      Close(12,status='delete')
      Open(unit=5,file=trim(curd)//'PresXZ.dat',action='write')
      Write(5,*) 'title = "flow"'
      Write(5,*) 'variables ="x","z","p","u","v","phi","vof"'!//',"Pres","Mres"'
      Write(5,1110) Imax,kmax
      Write(5,"(f18.10)")((TGrid%x(i,jp,k),i=1,Imax),k=1,kmax)
      Write(5,"(f18.10)")((TGrid%z(i,jp,k),i=1,Imax),k=1,kmax)
      Write(5,"(f18.10)")((p(i,jp,k),i=1,Imax),k=1,kmax)
      write(5,"(f18.10)")((0.5d0*(TVar%u(i-1,jp,k)+TVar%u(i,jp,k)),i=1,Imax),  &
                                                            k=1,kmax)
      write(5,"(f18.10)")((0.5d0*(TVar%w(i,jp,k-1)+TVar%v(i,jp,k)),i=1,Imax),  &
                                                            k=1,kmax)
      Write(5,"(f18.10)")((TCell%phil(i,jp,k),i=1,Imax),k=1,kmax)
      Write(5,"(f18.10)")((TCell%vofl(i,jp,k),i=1,Imax),k=1,kmax)
   !   Write(5,"(f18.10)")((TVar%pres(i,jp,k),i=1,Imax),k=1,kmax)
   !   Write(5,"(f18.10)")((TVar%mres(i,jp,k),i=1,Imax),k=1,kmax)
      Close(5)
 1110 format(1x,'zone i=',i5,',j=',i5,'f=block')
      Deallocate(p)
    End Subroutine PrintResultTecplotPCentXZ
 
    Subroutine PrintResultTecplotUCent(TGrid,TVar,TCell,itt)
      Implicit none
      type(Grid),intent(in):: TGrid
      type(Variables),intent(in):: TVar
      type(Cell),intent(in):: TCell
      Integer(kind=it8b),intent(in):: itt
      Integer(kind=it4b) i,j,k
      Character(15) curd
      Open(unit=11,file='text.txt',action='write')
      Write(11,*) itt
      Close(11)
      Open(unit=12,file='text.txt',action='read')
      Read(12,*) curd
      Close(12,status='delete')
      Open(unit=5,file=trim(curd)//'Uvelocity.dat',action='write')
      Write(5,*) 'title = "flow"'
      Write(5,*) 'variables ="x","y","z","u","Gp","phi","vof","ures"'
      Write(5,1110) Imax,Jmax,Kmax
      Write(5,"(f18.10)")(((TGrid%x(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TGrid%y(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TGrid%z(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%u(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%Gpu(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%phi(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%vof(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%ures(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Close(5)
 1110 format(1x,'zone i=',i5,',j=',i5,',k=',i5,'f=block')
    End Subroutine PrintResultTecplotUCent

    Subroutine PrintResultTecplotVCent(TGrid,TVar,TCell,itt)
      Implicit none
      type(Grid),intent(in):: TGrid
      type(Variables),intent(in):: TVar
      type(Cell),intent(in):: TCell
      Integer(kind=it8b),intent(in):: itt
      Integer(kind=it4b) i,j,k
      Character(15) curd
      Open(unit=11,file='text.txt',action='write')
      Write(11,*) itt
      Close(11)
      Open(unit=12,file='text.txt',action='read')
      Read(12,*) curd
      Close(12,status='delete')
      Open(unit=5,file=trim(curd)//'Vvelocity.dat',action='write')
      Write(5,*) 'title = "flow"'
      Write(5,*) 'variables ="x","y","z","v","Gp","phi","vof","vres"'
      Write(5,1110) Imax,Jmax,Kmax
      Write(5,"(f18.10)")(((TGrid%x(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TGrid%y(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TGrid%z(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%v(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%Gpv(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%phi(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%vof(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%vres(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Close(5)
 1110 format(1x,'zone i=',i5,',j=',i5,',k=',i5,'f=block')
    End Subroutine PrintResultTecplotVCent

    Subroutine PrintResultTecplotWCent(TGrid,TVar,TCell,itt)
      Implicit none
      type(Grid),intent(in):: TGrid
      type(Variables),intent(in):: TVar
      type(Cell),intent(in):: TCell
      Integer(kind=it8b),intent(in):: itt
      Integer(kind=it4b) i,j,k
      Character(15) curd
      Open(unit=11,file='text.txt',action='write')
      Write(11,*) itt
      Close(11)
      Open(unit=12,file='text.txt',action='read')
      Read(12,*) curd
      Close(12,status='delete')
      Open(unit=5,file=trim(curd)//'Wvelocity.dat',action='write')
      Write(5,*) 'title = "flow"'
      Write(5,*) 'variables ="x","y","z","w","Gp","phi","vof","wres"'
      Write(5,1110) Imax,Jmax,Kmax
      Write(5,"(f18.10)")(((TGrid%x(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TGrid%y(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TGrid%z(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%w(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%Gpw(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%phi(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%vof(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%wres(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Close(5)
 1110 format(1x,'zone i=',i5,',j=',i5,',k=',i5,'f=block')
    End Subroutine PrintResultTecplotWCent

    subroutine PrintResultVTK(TGrid,TVar,TCell,itt)
      use penf
      use vtk_fortran, only : vtk_file
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      implicit none
      type(vtk_file)                	 :: a_vtk_file                  !< A VTK file.
      type(Grid),INTENT(IN)         	 :: TGrid
      type(Variables),INTENT(IN)   	 :: TVar
      type(Cell),INTENT(IN)         	 :: TCell
      integer(kind=it8b),INTENT(IN) 	 :: itt 			!< Number of elements.
      real(R4P),allocatable,dimension(:) :: x,y,z                       !< X,Y,Z coordinates.
      integer(I4P)                  	 :: error                       !< Status error.
      integer(I4P)                  	 :: i,j,k
      character(15)                      :: curd
      !------------------------------------------------------------------------
      ! Print result for 3D data
      !------------------------------------------------------------------------
      allocate(x(IMax+1))
      allocate(y(JMax+1))
      allocate(z(Kmax+1))

      do i=1,IMax
        x(i)=TGrid%x(i,1,1)-TGrid%dx(i,1,1)/2.d0
      enddo
      x(Imax+1)=TGrid%x(Imax,1,1)+TGrid%dx(Imax,1,1)/2.d0
      do j=1,JMax
        y(j)=TGrid%y(1,j,1)-TGrid%dy(1,j,1)/2.d0
      enddo
      y(Jmax+1)=TGrid%y(1,Jmax,1)+TGrid%dy(1,Jmax,1)/2.d0
      do k=1,KMax
        z(k)=TGrid%z(1,1,k)-TGrid%dz(1,1,k)/2.d0
      end do
      z(Kmax+1)=TGrid%z(1,1,Kmax)+TGrid%dz(1,1,Kmax)/2.d0
      write(curd,'(i8.8)') itt
      error=a_vtk_file%initialize(format='BINARY', 			       &
      filename='PCell_Binary_3D'//trim(curd)//'.vtr', 		       	       &
                mesh_topology='RectilinearGrid',nx1=1,nx2=IMax+1,ny1=1,        &
                ny2=JMax+1,nz1=1,nz2=KMax+1)

   !   error = a_vtk_file%xml_writer%write_fielddata(action='open')
   !   error = a_vtk_file%xml_writer%write_fielddata(x=0._R8P, data_name='TIME')
   !   error = a_vtk_file%xml_writer%write_fielddata(x=1_I8P, data_name='CYCLE')
   !   error = a_vtk_file%xml_writer%write_fielddata(action='close')

      error=a_vtk_file%xml_writer%write_piece(nx1=1,nx2=IMax+1,ny1=1,          &
                                              ny2=JMax+1,nz1=1,nz2=KMax+1)
      error=a_vtk_file%xml_writer%write_geo(x=x,y=y,z=z)
      error=a_vtk_file%xml_writer%write_dataarray(location='cell',action='open')
      call WriteDataArrayVTK(a_vtk_file,TVar%p,"p")
      call WriteDataArrayVTK(a_vtk_file,TVar%u,"u")
      call WriteDataArrayVTK(a_vtk_file,TVar%v,"v")
      call WriteDataArrayVTK(a_vtk_file,TVar%w,"w")

      call WriteDataArrayVTK(a_vtk_file,TCell%phiL,"LiquidLvs")
      call WriteDataArrayVTK(a_vtk_file,TCell%phi,"FluidLvs")
      call WriteDataArrayVTK(a_vtk_file,TCell%vofL,"LiquidVof")
      call WriteDataArrayVTK(a_vtk_file,TCell%vof,"FluidVof")
      error=a_vtk_file%xml_writer%write_dataarray(location='cell',action='close')
      error=a_vtk_file%xml_writer%write_piece()
      error=a_vtk_file%finalize()
      deallocate(x,y,z)
    end subroutine PrintResultVTK

    subroutine WriteDataArrayVTK(a_vtk_file,InVar,VarName)
      use penf
      use vtk_fortran, only : vtk_file

      implicit none
      type(vtk_file),intent(inout)                      :: a_vtk_file
      real(R8P),dimension(:,:,:),allocatable,intent(in) :: InVar
      character(*),intent(in) 			        :: VarName
      integer(I4P)                  	   	        :: error                 !< Status error.
      integer(I4P)                  	   	        :: i,j,k
      integer(I4P)                  	   	        :: n
      real(R4P),allocatable,dimension(:)   	        :: v

      n=1
      allocate(v(1:IMax*JMax*KMax))
      v=0.d0
      do k=1,Kmax
        do j=1,Jmax
          do i=1,Imax
            v(n)=InVar(i,j,k)
            n=n+1
          end do
        end do
      end do
      error=a_vtk_file%xml_writer%write_dataarray(data_name=VarName,x=v)
      deallocate(v)
    end subroutine WriteDataArrayVTK
End module PrintResult
