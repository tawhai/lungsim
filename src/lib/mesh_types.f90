module mesh_types
  !*Brief Description:* This module defines types for the fundamental components of the lung model.
  !
  !*LICENSE:*
  !
  !
  !*Contributor(s):* Merryn Tawhai
  !
  !*Full Description:*
  !
  !This module defines types for the fundamental components of the lung model: airways, arteries, veins, tissue.
  ! Other types for lung-specific structures should be added as required. Parameters that are spatially-varying
  ! (e.g. different for each airway) should be included here, whereas parameters that apply to the whole model
  ! (e.g. cardiac output) should be defined in parameter_types. 
  
  use precision
  
  implicit none

  type airway_tree
     integer :: n_elems
     integer :: n_nodes
     integer, allocatable :: elems(:)
     integer, allocatable :: nodes(:)
  end type airway_tree

  type(airway_tree) :: airways

end module mesh_types
