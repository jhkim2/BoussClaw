subroutine setprob()

    use bous_module

    implicit none

    integer :: iunit
    character(len=25) fname

    iunit = 7
    fname = 'setprob.data'
!   # open the unit with new routine from Clawpack 4.4 to skip over
!   # comment lines starting with #:
    call opendatafile(iunit, fname)

    read(7,*) bouss 
    read(7,*) B_param 
    read(7,*) use_bous_sw_thresh 
    read(7,*) bous_sw_thresh 
    read(7,*) sw_depth

end subroutine setprob
