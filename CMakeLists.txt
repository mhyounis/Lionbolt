#  ==================
#    CMake settings  
#  ==================

CMAKE_MINIMUM_REQUIRED (VERSION 3.20)
PROJECT (Lionbolt LANGUAGES Fortran)

OPTION (BUILD_LB_EXEC      "Build Lionbolt program"           ON)
OPTION (ENABLE_DEV_COMPILE "Enable developer compile options" OFF)
# # # These are not currently in use, subject to change
OPTION (USE_MPI            "Use MPI"                          OFF)
OPTION (USE_OPENMP         "Use OpenMP"                       ON )
# # #

#  ====================
#    Fortran settings  
#  ====================

SET (CMAKE_Fortran_STANDARD 2008)

#  ========================================
#    Add NittanyPhysics as a subdirectory  
#  ========================================

# Make sure the executable isn't built
SET (BUILD_NP_EXEC OFF CACHE BOOL "" FORCE)

# If you want to override the GitHub download with a local version, then set
# the environment variables:
# export LB_NP_DEV_OVERRIDE=1
# export NITTANY=/path/to/your/NittanyPhysics
IF (DEFINED ENV{LB_NP_DEV_OVERRIDE} AND "$ENV{LB_NP_DEV_OVERRIDE}" STREQUAL "1")
    MESSAGE(STATUS "Using local NittanyPhysics: $ENV{NITTANY}")
    
    ADD_SUBDIRECTORY (
        $ENV{NITTANY}
        ${CMAKE_BINARY_DIR}/NittanyPhysics
    )
ELSE ()
    INCLUDE (FetchContent)
    
    MESSAGE (STATUS "Downloading NittanyPhysics from GitHub.")
    
    FetchContent_Declare (
        NittanyPhysics
        GIT_REPOSITORY https://github.com/mhyounis/NittanyPhysics.git
        GIT_TAG v1.0.0
    )
    
    FetchContent_MakeAvailable (NittanyPhysics)
ENDIF ()

#  ===================
#    Schema handling  
#  ===================

FIND_PACKAGE (Python REQUIRED COMPONENTS Interpreter)

SET (LIONBOLT_SCHEMA_YAML ${CMAKE_CURRENT_SOURCE_DIR}/src/IO/schema.yaml)
SET (LIONBOLT_SCHEMA      ${CMAKE_CURRENT_SOURCE_DIR}/src/IO/LionboltHDF5Schema.f90)

ADD_CUSTOM_COMMAND (
    OUTPUT  ${LIONBOLT_SCHEMA}
    COMMAND ${Python_EXECUTABLE} -u ${CORE_LIB_DIR}/tools/SchemaGenerator.py ${LIONBOLT_SCHEMA_YAML} ${LIONBOLT_SCHEMA} Lionbolt
    DEPENDS ${LIONBOLT_SCHEMA_YAML} ${CORE_LIB_DIR}/tools/SchemaGenerator.py
    VERBATIM
)

ADD_CUSTOM_TARGET (LionboltSchema DEPENDS ${LIONBOLT_SCHEMA})

#  =============================
#    Build the private library  
#  =============================

# First do a recursive glob to get all files and order them properly
FILE (GLOB_RECURSE PRIVATE_SOURCES
    src/*.f90
    third_party/SPARSKIT/*.f
)

# Append the schema file
LIST (APPEND PRIVATE_SOURCES ${LIONBOLT_SCHEMA})

# Define the executable sources (just for removal, this won't be used in the executable, for
# the sake of organization)
FILE (GLOB_RECURSE REMOVED_SOURCES
    src/main/*
    src/input/*
)

# Remove the main/ and input/ folders, as they are not part of the private library.
# As well as some temporary/under-development/other files.
LIST (REMOVE_ITEM PRIVATE_SOURCES
    ${REMOVED_SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/spaceangle/.depr/TransportDiagonals.f90 # Deprecated subroutine. Likely to make a comeback though
    ${CMAKE_CURRENT_SOURCE_DIR}/src/spaceangle/.depr/SweepArrays.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sources/.depr/PrepSource.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sources/.depr/BoundaryFluence.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sources/.depr/BoundarySource.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sources/.depr/TraceAngularFluence.f90
)

ADD_LIBRARY (LionboltPrivate OBJECT ${PRIVATE_SOURCES})

# Set the Lionbolt install directory to a variable for IO  
TARGET_COMPILE_DEFINITIONS (LionboltPrivate PRIVATE LIONBOLT_BASE_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}\")

# Make the schema build function a dependency
ADD_DEPENDENCIES (LionboltPrivate LionboltSchema)

# Link LionboltPrivate to core and NittanyAPI
TARGET_LINK_LIBRARIES (LionboltPrivate PUBLIC  core)
TARGET_LINK_LIBRARIES (LionboltPrivate PRIVATE NittanyAPI)

# Include the core and NittanyAPI module directory
TARGET_INCLUDE_DIRECTORIES (LionboltPrivate PRIVATE
    ${CMAKE_BINARY_DIR}/mod/core
    ${CMAKE_BINARY_DIR}/mod/NittanyAPI
)

# Set the module directory
SET_TARGET_PROPERTIES (LionboltPrivate PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/LionboltPrivate
)

# Set compile options
IF (ENABLE_DEV_COMPILE)
    TARGET_COMPILE_OPTIONS (LionboltPrivate PRIVATE
        -O3
        -march=native
        -finit-real=snan
        -cpp
        -fcheck=all
        -Wall
        -Wextra
    )
ELSE ()
    TARGET_COMPILE_OPTIONS (LionboltPrivate PRIVATE
        -O3
        -march=native
        -finit-real=snan
        -cpp
        -w
    )
ENDIF ()

# Link LionboltPrivate to BLAS and LAPACK
FIND_PACKAGE (BLAS REQUIRED)
FIND_PACKAGE (LAPACK REQUIRED)
TARGET_LINK_LIBRARIES  (LionboltPrivate PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

#  ========================
#    Build the public API  
#  ========================

ADD_LIBRARY (LionboltAPI STATIC src/LionboltAPI.f90)

# Link LionboltAPI to NittanyAPI
TARGET_LINK_LIBRARIES (LionboltAPI PRIVATE NittanyAPI)

# Link LionboltAPI to LionboltPrivate
TARGET_LINK_LIBRARIES (LionboltAPI PRIVATE 
    $<TARGET_OBJECTS:LionboltPrivate>
    ${BLAS_LIBRARIES} 
    ${LAPACK_LIBRARIES}
    ${OpenMP_Fortran_LIBRARIES}
)

# Set the module directory
SET_TARGET_PROPERTIES (LionboltAPI PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/LionboltAPI
)

# Include the LionboltPrivate module directory
TARGET_INCLUDE_DIRECTORIES (LionboltAPI PRIVATE ${CMAKE_BINARY_DIR}/mod/LionboltPrivate)

# Set compile options
IF (ENABLE_DEV_COMPILE)
    TARGET_COMPILE_OPTIONS (LionboltAPI PRIVATE
        -O3
        -march=native
        -finit-real=snan
        -cpp
        -fcheck=all
        -Wall
        -Wextra
    )
ELSE ()
    TARGET_COMPILE_OPTIONS (LionboltAPI PRIVATE
        -O3
        -march=native
        -finit-real=snan
        -cpp
        -w
    )
ENDIF ()

IF (BUILD_LB_EXEC)
    
    #  =================================
    #    Build the Lionbolt executable  
    #  =================================
    
    # If this is not a driver, use the main and input sources.
    # Otherwise, use the driver file
    IF (NOT DEFINED LBDRIVER)
        # Get the source files. For Lionbolt, the only non-API/internal source files are the entireties of main/ and input/
        FILE (GLOB_RECURSE LIONBOLT_SOURCES
            src/input/*
            src/main/*
        )
    ELSE ()
        SET (LIONBOLT_SOURCES ${LBDRIVER})
    ENDIF ()
    
    ADD_EXECUTABLE (Lionbolt ${LIONBOLT_SOURCES})
    
    # Set compile options
    IF (ENABLE_DEV_COMPILE)
        TARGET_COMPILE_OPTIONS (Lionbolt PRIVATE
            -O3
            -march=native
            -finit-real=snan
            -cpp
            -fcheck=all
            -Wall
            -Wextra
        )
    ELSE ()
        TARGET_COMPILE_OPTIONS (Lionbolt PRIVATE
            -O3
            -march=native
            -finit-real=snan
            -cpp
            -w
        )
    ENDIF ()
    
    # Link to core, NittanyAPI, and LionboltAPI
    TARGET_LINK_LIBRARIES (Lionbolt PRIVATE
        core
        NittanyAPI
        LionboltAPI
    )
    
    TARGET_INCLUDE_DIRECTORIES (Lionbolt PRIVATE
        ${CMAKE_BINARY_DIR}/mod/core
        ${CMAKE_BINARY_DIR}/mod/NittanyAPI
        ${CMAKE_BINARY_DIR}/mod/LionboltAPI
    )
ENDIF ()