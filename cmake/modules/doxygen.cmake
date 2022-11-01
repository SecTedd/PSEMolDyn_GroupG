find_package(Doxygen)

# make doc_doxygen optional if someone does not have / like doxygen –– Default is ON
option(BUILD_DOC "Build documentation (requires Doxygen)" OFF)

if(BUILD_DOC)
    if(NOT DOXYGEN_FOUND)
        message("Missing doxygen, required to generate the documentation")
    endif()

    # TODO: create CMake build option for the target.
    # TODO: Add a custom target for building the documentation.
    add_custom_target( 
        doc_doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in 
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

endif(BUILD_DOC)
