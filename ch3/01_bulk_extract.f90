program extract
  implicit none
  integer    :: i,j,k
  real(8)    :: kx,ky,kz
  complex(8) :: ham(80,80,80,3,3)
  real(8)    :: temp(12)
  complex(8) :: ci = (0.d0, 1.d0)

  ! read in
  open(unit=100, file='SVOt2gWANN80.hk')
  read(100,*)
  do i=1,512000
    read(100,*) kx,ky,kz
    do j=1,3
      read(100,*) (temp(k), k=1,6)
      do k=1,3
        ham( nint(kx*80)+1, nint(ky*80)+1, nint(kz*80)+1,j,k) =  temp(2*k-1) + ci * temp(2*k)
      enddo
    enddo
  enddo
  close(unit=100)

  write(*,*) 'Success.'

  ! write out

  open(unit=101, file='bulk_extracted_path.hk')
  do i=1,41 ! 0 0 0 to 0 pi 0 (X)
    write(101,*) (real(ham(1,i,1,j,j)), j=1,3)
  enddo
  do i=2,41 ! 0 pi 0 to pi pi 0 (M)
    write(101,*) (real(ham(i,41,1,j,j)), j=1,3)
  enddo
  do i=40,1,-1 ! pi pi 0 to 0 0 0
    write(101,*) (real(ham(i,i,1,j,j)), j=1,3)
  enddo
  do i=2,41 ! 0 0 0 to pi pi pi (R)
    write(101,*) (real(ham(i,i,i,j,j)), j=1,3)
  enddo
  close(unit=101)



end program extract
