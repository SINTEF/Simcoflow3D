project(Simcoflow Fortran)
find_package(MPI REQUIRED)
# add_definitions(${MPI_Fortran_COMPILE_FLAGS})
# Hypre,
# This should be ectracted out to a FindHYPRE.cmake
# The paths here are platform dependent. Please edit
set(Fortran_COMPILER_FLAGS -w)
set(CMAKE_REQUIRED_FLAGS ${Fortran_COMPILER_FLAGS})
set(HYPRE_INCLUDE_DIR /usr/include/hypre)
set(${PROJECT_NAME}_INCLUDE_DIR)
set(HYPRE_LIBRARIES /usr/lib/libHYPRE.so /usr/lib/libHYPRE_FEI.so /usr/lib/libHYPRE_struct_mv.so /usr/lib/libHYPRE_IJ_mv.so /usr/lib/libHYPRE_parcsr_ls.so)
set(HYPRE_LIBRARIES /usr/lib/x86_64-linux-gnu/libHYPRE_core.so /usr/lib/x86_64-linux-gnu/libHYPRE.so /usr/lib/x86_64-linux-gnu/libHYPRE_FEI.so)
include_directories(${MPI_Fortran_INCLUDE_PATH} ${HYPRE_INCLUDE_DIR} ${Simcoflow_INCLUDE_DIR})
link_directories(${MPI_Fortran_LIBRARIES})
#aux_source_directory(. SRC_LIST)
list(APPEND SRC_LIST
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/third_party/BeFoR64/src/lib/befor64.F90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/third_party/BeFoR64/src/lib/befor64_pack_data_m.F90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/third_party/FoXy/src/lib/foxy.f90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/third_party/FoXy/src/lib/foxy_xml_file.f90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/third_party/FoXy/src/lib/foxy_xml_tag.F90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/third_party/PENF/src/lib/penf_b_size.F90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/third_party/PENF/src/lib/penf.F90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/third_party/PENF/src/lib/penf_global_parameters_variables.F90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/third_party/PENF/src/lib/penf_stringify.F90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/third_party/StringiFor/src/lib/stringifor.F90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/third_party/StringiFor/src/lib/stringifor_string_t.F90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/lib/vtk_fortran_dataarray_encoder.f90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/lib/vtk_fortran.f90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/lib/vtk_fortran_parameters.f90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/lib/vtk_fortran_pvtk_file.f90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/lib/vtk_fortran_vtk_file.f90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/lib/vtk_fortran_vtk_file_xml_writer_abstract.f90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/lib/vtk_fortran_vtk_file_xml_writer_appended.f90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/lib/vtk_fortran_vtk_file_xml_writer_ascii_local.f90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/lib/vtk_fortran_vtk_file_xml_writer_binary_local.f90
   ${CMAKE_CURRENT_LIST_DIR}/VTKFortran/src/lib/vtk_fortran_vtm_file.F90 
   ${CMAKE_CURRENT_LIST_DIR}/Solver.f90
   ${CMAKE_CURRENT_LIST_DIR}/Precision.f90
   ${CMAKE_CURRENT_LIST_DIR}/MPI.f90
   ${CMAKE_CURRENT_LIST_DIR}/Mesh.f90
   ${CMAKE_CURRENT_LIST_DIR}/Matrix.f90
   ${CMAKE_CURRENT_LIST_DIR}/CutCell.f90 
   ${CMAKE_CURRENT_LIST_DIR}/Clsvof.f90
   ${CMAKE_CURRENT_LIST_DIR}/PrintResult.f90 
   ${CMAKE_CURRENT_LIST_DIR}/PredictUVW.f90
   ${CMAKE_CURRENT_LIST_DIR}/ProjectionP.f90 
   ${CMAKE_CURRENT_LIST_DIR}/ComputePUVW.f90
   ${CMAKE_CURRENT_LIST_DIR}/StateVariables.f90
   ${CMAKE_CURRENT_LIST_DIR}/Constants.f90
   ${CMAKE_CURRENT_LIST_DIR}/BoundaryInterface.f90
   ${CMAKE_CURRENT_LIST_DIR}/BoundaryFunction.f90
   ${CMAKE_CURRENT_LIST_DIR}/VTK_mod.f90
   ${CMAKE_CURRENT_LIST_DIR}/VTR_mod.f90)

set (LINK_LIBS ${LINK_LIBS} 
    ${MPI_Fortran_LIBRARIES}
    ${HYPRE_LIBRARIES}
    )
list(APPEND SRC_LIST_MAIN
    ${SRC_LIST}
    Main.f90)
add_executable(${PROJECT_NAME} ${SRC_LIST_MAIN})
target_link_libraries(${PROJECT_NAME} ${LINK_LIBS})
    
add_subdirectory(Test)
