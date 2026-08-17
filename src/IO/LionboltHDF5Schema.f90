
!  ================================================================================
!    DO NOT MODIFY THIS FILE DIRECTLY!                                             
!    Any changes you make will NOT be saved, because this file is                  
!    generated upon compilation by [NittanyPhysics]/core/tools/SchemaGenerator.py  
!                                                                                  
!    In order to add keys to the HDF5 file, you must modify the schema.yaml file:  
!        src/IO/schema.yaml                                                        
!  ================================================================================

Module LionboltHDF5Schema
    Implicit None

    Character (3), Parameter :: SCHEMA_VERSION = '1.0'

    ! MAIN_ATTRIBUTES
    Character (12), Parameter :: ATTR_PROBTYPE = 'problem_type'

    ! XSLIBRARY
    Character (10), Parameter :: GROUP_XSLIBRARY = 'XS_library'

    ! GENERAL
    Character (11), Parameter :: ATTR_OBJTYPE = 'object_type'

    ! MESH
    Character (4), Parameter :: GROUP_MESH = 'mesh'
    Character (12), Parameter :: DATASET_R = 'global_nodes'
    Character (12), Parameter :: DATASET_CONNECTIVITY = 'connectivity'
    Character (6), Parameter :: DATASET_OFFSET = 'offset'
    Character (15), Parameter :: DATASET_VOLS = 'element_volumes'
    Character (19), Parameter :: ATTR_NUMMATS_M = 'number_of_materials'
    Character (18), Parameter :: ATTR_NUMELS = 'number_of_elements'
    Character (22), Parameter :: ATTR_NUMSD = 'number_of_spatial_dofs'
    Character (22), Parameter :: ATTR_NUMKG = 'number_of_global_nodes'
    Character (23), Parameter :: GROUP_MAT2SD = 'material_to_spatial_dof'

    ! ANGULAR
    Character (22), Parameter :: ATTR_NUMANGLES = 'number_of_angular_dofs'
    Character (15), Parameter :: DATASET_ANGWEIGHTS = 'angular_weights'
    Character (9), Parameter :: DATASET_ABSCISSAE = 'abscissae'

    ! ENERGY
    Character (21), Parameter :: ATTR_NUMENERGIES = 'number_of_energy_dofs'
    Character (11), Parameter :: DATASET_EGRID = 'energy_grid'

    ! PARTICLES
    Character (8), Parameter :: IS_PARTICLE = 'particle'
    Character (7), Parameter :: GROUP_FL = 'fluence'
    Character (18), Parameter :: GROUP_FL_UNC = 'uncollided_fluence'
    Character (15), Parameter :: GROUP_ANG_FL = 'angular_fluence'
    Character (26), Parameter :: GROUP_ANG_FL_UNC = 'uncollided_angular_fluence'

    ! PATTERNS
    Character (7), Parameter :: PATTERN_ENERGY = 'energy_'
    Character (6), Parameter :: PATTERN_ANGLE = 'angle_'
    Character (9), Parameter :: PATTERN_MATERIALS = 'material_'

    ! POSTPROC
    Character (11), Parameter :: PATTERN_DEP = '_deposition'
    Character (6), Parameter :: DATASET_EDEP = 'energy'
    Character (4), Parameter :: DATASET_DDEP = 'dose'
    Character (6), Parameter :: DATASET_CDEP = 'charge'
    Character (12), Parameter :: ATTR_EDEP = 'total_energy'
    Character (10), Parameter :: ATTR_DDEP = 'total_dose'
    Character (12), Parameter :: ATTR_CDEP = 'total_charge'

End Module