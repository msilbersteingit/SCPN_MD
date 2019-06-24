	program Finter
C parameter
C a
C boxlength
C s's

	parameter(num_atoms=54000)
	real a(num_atoms,6)
	real force(2000,3)
	real net_force(3)
	real tmp_force(3)
	real boxlength(3,2)
	real xyz(3)
	real ref(3)
	real vec(3)
	real dist
	character(len=10) s1(2), s2, s3(4), s4, s5(5), s9(8)

	open(1,file='dump.0')
	read(1,*) s1
	read(1,*) s2
	read(1,*) s3
	read(1,*) s4
	read(1,*) s5
	read(1,*) boxlength(1,1), boxlength(1,2)
	read(1,*) boxlength(2,1), boxlength(2,2)
	read(1,*) boxlength(3,1), boxlength(3,2)
	read(1,*) s9
	do j = 1, num_atoms, 1
		read(1,*) x1, x2, x3, x4, x5, x6
		a(j,1)=x1
		a(j,2)=x2
		a(j,3)=x3
		a(j,4)=x4
		a(j,5)=x5
		a(j,6)=x6
	enddo
	close(1)

	call sortrow(a,num_atoms) ! sort rows with the first column ascending

	do i = 1, 2000, 1
		net_force=(/ 0, 0, 0 /)

		do j = 1, 2000, 1
			tmp_force=(/ 0, 0, 0 /)

			if (abs(i-j) .gt. 3) then
				xyz=a(j,4:)
				ref=a(i,4:)
				call pbc(xyz,ref,boxlength)

				vec=xyz-ref

				dist=distcal(vec)
				if ((dist .le. 10.5) .and. (dist .ge. 3.0)) then
					tmp_force=-vec/dist*forcecal(dist,0.112,4.01)
				endif

				net_force=net_force+tmp_force
			endif
		enddo

		force(i,:)=net_force
	enddo

	open(2,file='force_intra.txt')
	do i = 1, 2000, 1
		write (2,*) force(i,:)
	enddo
	close(2)

	stop
	end
C -----------------------------------------------------------

	subroutine sortrow(A,num_rows)
	real A(num_rows,6), buf(6)
	integer irow, krow

	do irow = 1, num_rows, 1
		krow=minloc(A(irow:num_rows,1),dim=1)+irow-1
		buf(:)=A(irow,:)
		A(irow,:)=A(krow,:)
		A(krow,:)=buf(:)
	enddo

	return
	end

C -----------------------------------------------------------
	subroutine pbc(xyz,ref,boxlength)
	real xyz(3)
	real ref(3)
	real boxlength(3,2)
	real reflen(3)
	real vec(3)
	real bl(3)

	do j = 1, 3, 1
		bl(j)=boxlength(j,2)-boxlength(j,1)
		reflen(j)=bl(j)-10.0
	enddo

	do j = 1, 3, 1
		vec(j)=xyz(j)-ref(j)

		if (vec(j) .lt. reflen(j)*(-1)) then
			vec(j)=vec(j)+bl(j)
		endif
		if (vec(j) .gt. reflen(j)) then
			vec(j)=vec(j)-bl(j)
		endif
	enddo

	xyz=ref+vec

	return
	end

C -----------------------------------------------------------

	real function distcal(a)
	real a(3)

	distcal=sqrt(a(1)**2+a(2)**2+a(3)**2)

	return
	end

C -----------------------------------------------------------

	real function forcecal(r,epsilon,sigma)
	real r
	real epsilon
	real sigma
	real sr

	sr=sigma/r
	forcecal=-4*epsilon/r*((-12)*sr**12+6*sr**6)

	return
	end
