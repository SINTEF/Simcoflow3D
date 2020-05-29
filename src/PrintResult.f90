Module PrintResult
    USE PrecisionVar
    USE Mesh
    USE StateVariables,ONLY : ight,jght,kght,variables,pi,Ktw,rey,xc,yc
    USE Clsvof
    USE Cutcell
    Implicit none
    Private
    character(len=1), parameter :: newline=achar(10)
    Public:: PrintResultTecplotPCent,PrintResultTecplotUCent,                  &
             PrintResultTecplotVCent,PrintResultTecplotWCent,                  &
             PrintResultTecplotPCentXY,PrintResultVTK
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
      Write(5,*) 'variables ="x","y","z","p","u","v","w","phi","vof","Pres","Mres"'
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
      Write(5,"(f18.10)")(((TCell%phi(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TCell%vof(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%pres(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Write(5,"(f18.10)")(((TVar%mres(i,j,k),i=1,Imax),j=1,Jmax),k=1,Kmax)
      Close(5)
 1110 format(1x,'zone i=',i5,',j=',i5,',k=',i5,'f=block')
      Deallocate(p)
    End Subroutine PrintResultTecplotPCent

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
      real(R4P),allocatable,dimension(:) :: x,y,z                       !< X coordinates.
      integer(I4P)                  	 :: error                       !< Status error.
      integer(I4P)                  	 :: i,j,k                       
      character(15)                      :: curd
!------------------------------------------------------------------------------
      ! Print result for 3D data  
!------------------------------------------------------------------------------
      allocate(x(0:IMax-1))
      allocate(y(0:JMax-1))
      allocate(z(0:Kmax-1))
      
      do i=0,IMax-1
        x(i)=TGrid%x(i+1,1,1)
      enddo
      do j=0,JMax-1
        y(j)=TGrid%y(1,j+1,1)
      enddo
      do k=0,KMax-1
        z(k)=TGrid%z(1,1,k+1)
      end do
      write(curd,'(i8.8)') itt
      error=a_vtk_file%initialize(format='BINARY', 			       &
      filename='PCell_Binary_3D'//trim(curd)//'.vtr', 		       	       &
                mesh_topology='RectilinearGrid',nx1=0,nx2=IMax-1,ny1=0,        &
                ny2=JMax-1,nz1=0,nz2=KMax-1) 

   !   error = a_vtk_file%xml_writer%write_fielddata(action='open')
   !   error = a_vtk_file%xml_writer%write_fielddata(x=0._R8P, data_name='TIME')
   !   error = a_vtk_file%xml_writer%write_fielddata(x=1_I8P, data_name='CYCLE')
   !   error = a_vtk_file%xml_writer%write_fielddata(action='close')
      
      error=a_vtk_file%xml_writer%write_piece(nx1=0,nx2=IMax-1,ny1=0,          &
                                              ny2=JMax-1,nz1=0,nz2=KMax-1)
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
      !-------------------------------------------------------------------------
      ! Print result for 2D
      !-------------------------------------------------------------------------
      deallocate(z)
      allocate(z(1))
      z(1)=0.d0
      error=a_vtk_file%initialize(format='BINARY', 			       &
      filename='PCell_Binary_2D'//trim(curd)//'.vtr', 		       	       &
                mesh_topology='RectilinearGrid',nx1=0,nx2=IMax-1,ny1=0,        &
                ny2=JMax-1,nz1=1,nz2=1) 
   !   test_passed = .true. ! nothing to test yet

   !   print "(A,L1)", new_line('a')//'Are all tests passed? ', all(test_passed)
   !   stop
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
      do i=1,Imax
        do j=1,Jmax
          do k=1,KMax
            v(n)=InVar(i,j,k)
            n=n+1
          end do
        end do
      end do  
      error=a_vtk_file%xml_writer%write_dataarray(data_name=VarName,x=v)
      deallocate(v)	       		
    end subroutine WriteDataArrayVTK	
End module PrintResult
