#--////////////////////////////////////////////////////////////////////////////
#--
#--  PLEASE PROVIDE COPYRIGHT INFORMATION HERE.
#--
#--////////////////////////////////////////////////////////////////////////////

set(PLUGIN_NAME "FiberToolbox")
set(${PLUGIN_NAME}Test_SOURCE_DIR ${${PLUGIN_NAME}_SOURCE_DIR}/Test)
set(${PLUGIN_NAME}Test_BINARY_DIR ${${PLUGIN_NAME}_BINARY_DIR}/Test)

include(${SIMPLProj_SOURCE_DIR}/Source/SIMPLib/Testing/SIMPLTestFunctions.cmake)

include_directories(${${PLUGIN_NAME}_SOURCE_DIR}/Test)
include_directories(${${PLUGIN_NAME}_BINARY_DIR})
include_directories(${SIMPLProj_SOURCE_DIR}/Source)
include_directories(${SIMPLProj_BINARY_DIR})

#------------------------------------------------------------------------------
# Set any additional include directories here

set(TESTFILES_SRC_DIR ${${PLUGIN_NAME}Test_SOURCE_DIR}/TestFiles)
set(TESTFILES_BINARY_DIR ${${PLUGIN_NAME}Test_BINARY_DIR}/TestFiles)
set(TEST_TEMP_DIR ${${PLUGIN_NAME}Test_BINARY_DIR}/Temp)

# Make sure the directory is created during CMake time
file(MAKE_DIRECTORY ${TEST_TEMP_DIR})


configure_file(${${PLUGIN_NAME}_SOURCE_DIR}/Test/TestFileLocations.h.in
               ${${PLUGIN_NAME}_BINARY_DIR}/Test/${PLUGIN_NAME}TestFileLocations.h @ONLY IMMEDIATE)

set(${PLUGIN_NAME}_Link_Libs Qt5::Core H5Support SIMPLib)

#------------------------------------------------------------------------------
# List all the source files here. They will NOT be compiled but instead
# be directly included in the main test source file. We list them here so that
# they will show up in IDEs
set(TEST_NAMES
  DetectEllipsoids

)

set( ${PLUGIN_NAME}_TEST_SRCS )
set(FilterTestIncludes "")
set(TestMainFunctors "")

foreach(name  ${TEST_NAMES})
  set( ${PLUGIN_NAME}_TEST_SRCS
    ${${PLUGIN_NAME}_TEST_SRCS}
    "${${PLUGIN_NAME}Test_SOURCE_DIR}/${name}Test.cpp"
    )
  string(CONCAT
    FilterTestIncludes
    ${FilterTestIncludes}
    "#include \"${name}Test.cpp\"\n"
    )

    string(CONCAT
      TestMainFunctors
     ${TestMainFunctors}
     "  ${name}Test()()|\n")
endforeach()

STRING(REPLACE "|" ";" TestMainFunctors ${TestMainFunctors}   )

configure_file(${SIMPLProj_SOURCE_DIR}/Source/SIMPLib/Testing/TestMain.cpp.in
               ${${PLUGIN_NAME}_BINARY_DIR}/${PLUGIN_NAME}UnitTest.cpp @ONLY)


# Set the source files properties on each source file.
foreach(f ${${PLUGIN_NAME}_TEST_SRCS})
  set_source_files_properties( ${f} PROPERTIES HEADER_FILE_ONLY TRUE)
endforeach()

AddSIMPLUnitTest(TESTNAME ${PLUGIN_NAME}UnitTest
                    SOURCES ${${PLUGIN_NAME}_BINARY_DIR}/${PLUGIN_NAME}UnitTest.cpp ${${PLUGIN_NAME}_TEST_SRCS}
                    FOLDER "${PLUGIN_NAME}Plugin/Test"
                    LINK_LIBRARIES ${${PLUGIN_NAME}_Link_Libs})

target_include_directories(${PLUGIN_NAME}UnitTest
                          PUBLIC ${PLUGINS_SOURCE_DIR}
                          PUBLIC ${PLUGINS_BINARY_DIR}
                          PUBLIC ${PLUGINS_BINARY_DIR}/${PLUGIN_NAME}/Test
)

if(MSVC)
  set_source_files_properties(${${PLUGIN_NAME}Test_BINARY_DIR}/${PLUGIN_NAME}UnitTest.cpp PROPERTIES COMPILE_FLAGS /bigobj)
endif()

#-------------------------------------------------------------------------------
#- This copies all the Prebuilt Pipeline files into the Build directory so the help
#- works from the Build Tree
add_custom_target(TestFileCopy ALL
            COMMAND ${CMAKE_COMMAND} -E copy_directory ${TESTFILES_SRC_DIR} ${TESTFILES_BINARY_DIR}
            COMMENT "Copying ${PLUGIN_NAME} Test Files to Build directory...")
set_target_properties(TestFileCopy PROPERTIES FOLDER ZZ_COPY_FILES)
