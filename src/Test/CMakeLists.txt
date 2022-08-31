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

list(APPEND ParallelFlowLIST
    ${SRC_LIST}
    ParallelFlowTest/InitialVof.f90
    ParallelFlowTest/ParallelFlowTest.f90)
add_executable(ParallelFlowTest ${ParallelFlowLIST})
target_link_libraries(ParallelFlowTest ${LINK_LIBS})
my_add_test(testparallelflow)

list(APPEND ADVECTIONLIST
    ${SRC_LIST}
    AdvectionTest/InitialVof.f90
    AdvectionTest/AdvectionTest.f90)
add_executable(AdvectionTest ${ADVECTIONLIST})
target_link_libraries(AdvectionTest ${LINK_LIBS})
my_add_test(testadvection)

list(APPEND PlanarLIST
    ${SRC_LIST}
    PlanarTest/InitialVof.f90
    PlanarTest/PlanarTest.f90)
add_executable(PlanarTest ${PlanarLIST})
target_link_libraries(PlanarTest ${LINK_LIBS})
my_add_test(testplanar)

list(APPEND DAMBREAKLIST
    ${SRC_LIST}
    Dambreak/InitialVof.f90
    Dambreak/Dambreak.f90)
add_executable(Dambreak ${DAMBREAKLIST})
target_link_libraries(Dambreak ${LINK_LIBS})
my_add_test(testdambreak)

list(APPEND FLOWSPHERELIST
    ${SRC_LIST}
    FlowSphere/FlowSphere.f90)
add_executable(FlowSphere ${FLOWSPHERELIST})
target_link_libraries(FlowSphere ${LINK_LIBS})
my_add_test(testflowsphere)

list(APPEND VESSELLIST
    ${SRC_LIST}
    Vessel/Vessel.f90
    Vessel/InitializationVessel.f90)
add_executable(Vessel ${VESSELLIST})
target_link_libraries(Vessel ${LINK_LIBS})
my_add_test(testvessel)

list(APPEND FLOWSPHERE2PHASESLIST
    ${SRC_LIST}
    FlowSphere2Phases/FlowSphere2Phases.f90
    FlowSphere2Phases/InitialVof.f90)
add_executable(FlowSphere2Phases ${FLOWSPHERE2PHASESLIST})
target_link_libraries(FlowSphere2Phases ${LINK_LIBS})
my_add_test(testflowsphere2phases)

list(APPEND WAVELIST
    ${SRC_LIST}
    WavesSimulation/WavesMain.f90
    WavesSimulation/InitialVof.f90)
    add_executable(WavesSimulation ${WAVELIST})
    target_link_libraries(WavesSimulation ${LINK_LIBS})
my_add_test(testwave)

