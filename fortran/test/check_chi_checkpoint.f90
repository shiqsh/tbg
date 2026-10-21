program check_chi_checkpoint
    use CHI_CHECKPOINT
    implicit none
    integer :: unit, count, clock_value, saved
    real(8) :: values(8,2), expected(8,2)
    character(80) :: path
    call system_clock(count=clock_value)
    write(path,'(A,I0,A)') 'check_chi_checkpoint_',clock_value,'.tmp'
    expected(:,1) = [1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,8.d0,10.d0]
    expected(:,2) = 2.d0*expected(:,1)
    call CREATE_CHI_CHECKPOINT(trim(path),'strain_epsilon',[2,3,0], &
        [1.d-200,1.d200], [0.d0,1.d-4], [0.5d0], unit)
    call APPEND_CHI_POINT(unit,1,expected(:,1))
    call APPEND_CHI_POINT(unit,2,expected(:,2))
    close(unit)
    call LOAD_CHI_CHECKPOINT(trim(path),'strain_epsilon',[2,3,0], &
        [1.d-200,1.d200], [0.d0,1.d-4], [0.5d0], values, saved)
    if (saved /= 2 .or. any(values /= expected)) error stop 'Checkpoint roundtrip failed'
    open(newunit=unit,file=trim(path),status='old')
    close(unit,status='delete')
    print *, 'Checkpoint roundtrip passed'
end program
