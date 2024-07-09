module lung_mechanics_c
!*Brief Description:* This module wraps part of the lung_mechanics module that require a c interface
!
!*LICENSE:*
!
!
!*Contributor(s):* Merryn Tawhai
!
!*Full Description:*
!
!This module wraps part of the lung_mechanics module that require a c interface
  use arrays,only: dp

  implicit none
  !Interfaces
  private
  
contains
  !
  !######################################################################
  !
  !>
  subroutine deform_tissue_in_cavity_c(nsteps, posture, posture_len, filename, filename_len) &
       bind(C, name="deform_tissue_in_cavity_c")
    
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_STRING_LEN
    use lung_mechanics, only: deform_tissue_in_cavity
    implicit none
    
    integer,intent(in) :: nsteps, posture_len, filename_len
    type(c_ptr), value, intent(in) :: posture, filename
    character(len=MAX_STRING_LEN) :: posture_f, filename_f
    
    call strncpy(posture_f, posture, posture_len)
    call strncpy(filename_f, filename, filename_len)
    call deform_tissue_in_cavity(nsteps, posture_f, filename_f)
    
  end subroutine deform_tissue_in_cavity_c
  
  
end module lung_mechanics_c
