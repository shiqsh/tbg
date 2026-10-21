! Save completed scan points independently of the final, integrated output.
! Files are never overwritten; resume creates a new output/checkpoint pair.
module CHI_CHECKPOINT
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: LOAD_CHI_CHECKPOINT, CREATE_CHI_CHECKPOINT, APPEND_CHI_POINT
contains
    subroutine CHECK_IO(ios)
        integer, intent(in) :: ios
        if (ios /= 0) error stop 'CHI checkpoint: I/O error (file missing, exists, truncated or disk full)'
    end subroutine

    subroutine LOAD_CHI_CHECKPOINT(path, scan, ints, reals, eps_grid, phi_grid, values, completed)
        character(*), intent(in) :: path, scan
        integer, intent(in) :: ints(:)
        real(8), intent(in) :: reals(:), eps_grid(:), phi_grid(:)
        real(8), intent(out) :: values(:,:)
        integer, intent(out) :: completed
        integer :: unit, ios, index, j, ntokens
        integer :: saved_ints(size(ints))
        real(8) :: saved_reals(size(reals)), saved_eps(size(eps_grid)), saved_phi(size(phi_grid)), row(8)
        character(1024) :: line
        logical :: in_token

        if (size(values,1) /= 8) error stop 'CHI checkpoint: expected eight result columns'
        values = 0.d0
        completed = 0
        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        call CHECK_IO(ios)
        read(unit,'(A)',iostat=ios) line
        call CHECK_IO(ios)
        if (trim(line) /= 'TBG_CHI_CHECKPOINT_V1') error stop 'CHI checkpoint: wrong format'
        read(unit,'(A)',iostat=ios) line
        call CHECK_IO(ios)
        if (trim(line) /= trim(scan)) error stop 'CHI checkpoint: different scan parameter'
        read(unit,*,iostat=ios) saved_ints
        call CHECK_IO(ios)
        if (any(saved_ints /= ints)) error stop 'CHI checkpoint: different dimensions or index ranges'
        read(unit,*,iostat=ios) saved_reals
        call CHECK_IO(ios)
        if (any(saved_reals /= reals)) error stop 'CHI checkpoint: different physical/numerical parameters'
        read(unit,*,iostat=ios) saved_eps
        call CHECK_IO(ios)
        read(unit,*,iostat=ios) saved_phi
        call CHECK_IO(ios)
        if (any(saved_eps /= eps_grid) .or. any(saved_phi /= phi_grid)) &
            error stop 'CHI checkpoint: different parameter grids'
        do
            read(unit,'(A)',iostat=ios) line
            if (ios == iostat_end) exit
            call CHECK_IO(ios)
            ! Require exactly one integer and eight reals on each record.
            ntokens = 0
            in_token = .false.
            do j = 1, len_trim(line)
                if (line(j:j) == ' ' .or. line(j:j) == achar(9)) then
                    in_token = .false.
                else if (.not. in_token) then
                    ntokens = ntokens + 1
                    in_token = .true.
                endif
            enddo
            if (ntokens /= 9) error stop 'CHI checkpoint: incomplete or malformed point record'
            read(line,*,iostat=ios) index, row
            call CHECK_IO(ios)
            if (index /= completed+1 .or. index > size(values,2)) &
                error stop 'CHI checkpoint: nonconsecutive/out-of-range point'
            if (.not. all(ieee_is_finite(row))) error stop 'CHI checkpoint: nonfinite result'
            if (any(abs(row(7:8)-(row(3:4)+row(5:6))) > &
                    1.d-12*max(1.d0, maxval(abs(row))))) error stop 'CHI checkpoint: inconsistent chi total'
            values(:,index) = row
            completed = index
        enddo
        close(unit, iostat=ios)
        call CHECK_IO(ios)
    end subroutine

    subroutine CREATE_CHI_CHECKPOINT(path, scan, ints, reals, eps_grid, phi_grid, unit)
        character(*), intent(in) :: path, scan
        integer, intent(in) :: ints(:)
        real(8), intent(in) :: reals(:), eps_grid(:), phi_grid(:)
        integer, intent(out) :: unit
        integer :: ios
        open(newunit=unit, file=path, status='new', action='write', iostat=ios)
        call CHECK_IO(ios)
        write(unit,'(A)',iostat=ios) 'TBG_CHI_CHECKPOINT_V1'
        call CHECK_IO(ios)
        write(unit,'(A)',iostat=ios) trim(scan)
        call CHECK_IO(ios)
        write(unit,'(*(I0,1X))',iostat=ios) ints
        call CHECK_IO(ios)
        write(unit,'(*(ES26.17E3,1X))',iostat=ios) reals
        call CHECK_IO(ios)
        write(unit,'(*(ES26.17E3,1X))',iostat=ios) eps_grid
        call CHECK_IO(ios)
        write(unit,'(*(ES26.17E3,1X))',iostat=ios) phi_grid
        call CHECK_IO(ios)
        flush(unit, iostat=ios)
        call CHECK_IO(ios)
    end subroutine

    subroutine APPEND_CHI_POINT(unit, index, row)
        integer, intent(in) :: unit, index
        real(8), intent(in) :: row(8)
        integer :: ios
        if (.not. all(ieee_is_finite(row))) error stop 'CHI checkpoint: refusing nonfinite result'
        write(unit,'(I0,1X,8(ES26.17E3,1X))',iostat=ios) index, row
        call CHECK_IO(ios)
        flush(unit, iostat=ios)
        call CHECK_IO(ios)
    end subroutine
end module CHI_CHECKPOINT
