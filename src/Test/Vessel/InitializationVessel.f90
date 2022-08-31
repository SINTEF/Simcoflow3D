Module InitializationVessel
      use PrecisionVar
      Use Mesh
      Use StateVariables

      implicit none
      private
      public InitialGridVessel
      contains


      subroutine InitialGridVessel(Lref,Spoint,Epoint, TGrid)
        implicit none
        type(Point),intent(in)   :: SPoint,EPoint
        real(kind=dp),intent(in) :: Lref
        type(Grid),intent(inout) :: TGrid
        Integer                  :: i,j,k

        TGrid%Lref = Lref
        TGrid%dx(:,:,:) = (EPoint%x/Lref-SPoint%x/Lref)/dble(Imax)
        TGrid%dy(:,:,:) = (EPoint%y/Lref-SPoint%y/Lref)/dble(Jmax)
        TGrid%dz(:,:,:) = (EPoint%z/Lref-SPoint%z/Lref)/dble(Kmax)

        do i = 1,Imax
          do j = 1,Jmax
            do k = 1,Kmax
              TGrid%x(i,j,k) = SPoint%x/Lref+TGrid%dx(i,j,k)*(dble(i)-0.5d0)
              TGrid%y(i,j,k) = SPoint%y/Lref+TGrid%dy(i,j,k)*(dble(j)-0.5d0)
              TGrid%z(i,j,k) = SPoint%z/Lref+TGrid%dz(i,j,k)*(dble(k)-0.5d0)
            end do
          end do
        end do
      end subroutine InitialGridVessel
!
end module InitializationVessel
