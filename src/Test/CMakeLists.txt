project(simcoflowtest Fortran)
set (CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(
    ${Simcoflow_INCLUDE_DIR}
    )
set(HYPRE_INCLUDE_DIR /usr/include/hypre/)
set(${PROJECT_NAME}_INCLUDE_DIR)
set(HYPRE_LIBRARIES /usr/lib/x86_64-linux-gnu/libHYPRE_core.so /usr/lib/x86_64-linux-gnu/libHYPRE.so /usr/lib/x86_64-linux-gnu/libHYPRE_FEI.so)
set (LINK_LIBS ${LINK_LIBS})
macro(my_add_test Test)
    add_test(NAME ${Test} WORKING_DIRECTORY ${PROJECT_BINARY_DIR} COMMAND
        ${Test})
endmacro()
#aux_source_directory(../ SRC_LIST)
list(APPEND SRC_LIST2
    ../VTKFortran/src/third_party/BeFoR64/src/lib/befor64.F90
    ../VTKFortran/src/third_party/BeFoR64/src/lib/befor64_pack_data_m.F90
    ../VTKFortran/src/third_party/FoXy/src/lib/foxy.f90
    ../VTKFortran/src/third_party/FoXy/src/lib/foxy_xml_file.f90
    ../VTKFortran/src/third_party/FoXy/src/lib/foxy_xml_tag.F90
    ../VTKFortran/src/third_party/PENF/src/lib/penf_b_size.F90
    ../VTKFortran/src/third_party/PENF/src/lib/penf.F90
    ../VTKFortran/src/third_party/PENF/src/lib/penf_global_parameters_variables.F90
    ../VTKFortran/src/third_party/PENF/src/lib/penf_stringify.F90
    ../VTKFortran/src/third_party/StringiFor/src/lib/stringifor.F90
    ../VTKFortran/src/third_party/StringiFor/src/lib/stringifor_string_t.F90
    ../VTKFortran/src/lib/vtk_fortran_dataarray_encoder.f90
    ../VTKFortran/src/lib/vtk_fortran.f90
    ../VTKFortran/src/lib/vtk_fortran_parameters.f90
    ../VTKFortran/src/lib/vtk_fortran_pvtk_file.f90
    ../VTKFortran/src/lib/vtk_fortran_vtk_file.f90
    ../VTKFortran/src/lib/vtk_fortran_vtk_file_xml_writer_abstract.f90
    ../VTKFortran/src/lib/vtk_fortran_vtk_file_xml_writer_appended.f90
    ../VTKFortran/src/lib/vtk_fortran_vtk_file_xml_writer_ascii_local.f90
    ../VTKFortran/src/lib/vtk_fortran_vtk_file_xml_writer_binary_local.f90
    ../VTKFortran/src/lib/vtk_fortran_vtm_file.F90 
    ../Solver.f90
    ../Precision.f90
    ../MPI.f90
    ../Mesh.f90
    ../Matrix.f90
    ../CutCell.f90 
    ../Clsvof.f90
    ../PrintResult.f90 
    ../PredictUVW.f90
    ../ProjectionP.f90 
    ../ComputePUVW.f90
    ../StateVariables.f90
    ../Constants.f90
    ../BoundaryInterface.f90
    ../BoundaryFunction.f90
    ../VTK_mod.f90
    ../VTR_mod.f90  	
)
list(APPEND DAMBREAKLIST
    ${SRC_LIST2}
    Dambreak/InitialVof.f90
    Dambreak/Dambreak.f90)
add_executable(Dambreak ${DAMBREAKLIST})
target_link_libraries(Dambreak ${LINK_LIBS})
my_add_test(testdambreak)
list(APPEND FLOWSPHERELIST
    ${SRC_LIST2}
    FlowSphere/FlowSphere.f90)
add_executable(FlowSphere ${FLOWSPHERELIST})
target_link_libraries(FlowSphere ${LINK_LIBS})
my_add_test(testflowsphere)
