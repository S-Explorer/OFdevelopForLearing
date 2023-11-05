add_executable(kernelPhaseTransitionFOAM 
    kernelPhaseChangeFOAM.C
    fluid/compressibleCourantNo.C
    solid/solidRegionDiffNo.C
)

target_include_directories(kernelPhaseTransitionFOAM PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/fluid
    ${CMAKE_CURRENT_SOURCE_DIR}/solid
    ${OpenFOAM_SRC}/OSspecific/POSIX/lnInclude
    ${OpenFOAM_SRC}/OpenFOAM/lnInclude
    ${OpenFOAM_SRC}/finiteVolume/lnInclude 
    ${OpenFOAM_SRC}/meshTools/lnInclude 
    ${OpenFOAM_SRC}/sampling/lnInclude 
    ${OpenFOAM_SRC}/transportModels/compressible/lnInclude 
    ${OpenFOAM_SRC}/thermophysicalModels/basic/lnInclude 
    ${OpenFOAM_SRC}/thermophysicalModels/specie/lnInclude 
    ${OpenFOAM_SRC}/thermophysicalModels/reactionThermo/lnInclude 
    ${OpenFOAM_SRC}/thermophysicalModels/solidThermo/lnInclude 
    ${OpenFOAM_SRC}/thermophysicalModels/chemistryModel/lnInclude 
    ${OpenFOAM_SRC}/ODE/lnInclude 
    ${OpenFOAM_SRC}/combustionModels/lnInclude 
    ${OpenFOAM_SRC}/TurbulenceModels/turbulenceModels/lnInclude 
    ${OpenFOAM_SRC}/TurbulenceModels/compressible/lnInclude 
    ${OpenFOAM_SRC}/thermophysicalModels/radiation/lnInclude 
    ${OpenFOAM_SRC}/regionModels/regionModel/lnInclude
)

target_link_directories(kernelPhaseTransitionFOAM PRIVATE 
    ${OpenFOAM_LIB}
    ${OpenFOAM_LIB}/sys-openmpi
    ${OpenFOAM_LIB}/dummy
)

target_link_libraries(kernelPhaseTransitionFOAM PRIVATE
    # OSspecific
    Pstream
    OpenFOAM
    finiteVolume 
    fvOptions 
    meshTools 
    sampling 
    compressibleTransportModels 
    fluidThermophysicalModels 
    specie 
    reactionThermophysicalModels 
    solidThermo 
    chemistryModel 
    ODE 
    combustionModels 
    turbulenceModels 
    compressibleTurbulenceModels 
    radiationModels 
    regionModels
)

install(TARGETS kernelPhaseTransitionFOAM DESTINATION ${OpenFOAM_USR_BIN})
