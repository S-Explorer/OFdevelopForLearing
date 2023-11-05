add_library(MassTransferCoupled SHARED ./massTransferCoupledFvPatchScalarField.C)

target_include_directories(MassTransferCoupled PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/lnInclude
    ${OpenFOAM_SRC}/OSspecific/POSIX/lnInclude
    ${OpenFOAM_SRC}/OpenFOAM/lnInclude
    ${OpenFOAM_SRC}/TurbulenceModels/compressible/lnInclude 
    ${OpenFOAM_SRC}/transportModels/compressible/lnInclude 
    ${OpenFOAM_SRC}/thermophysicalModels/basic/lnInclude 
    ${OpenFOAM_SRC}/thermophysicalModels/specie/lnInclude 
    ${OpenFOAM_SRC}/thermophysicalModels/reactionThermo/lnInclude 
    ${OpenFOAM_SRC}/thermophysicalModels/solidThermo/lnInclude 
    ${OpenFOAM_SRC}/thermophysicalModels/solidSpecie/lnInclude 
    ${OpenFOAM_SRC}/finiteVolume/lnInclude 
    ${OpenFOAM_SRC}/meshTools/lnInclude 
    ${OpenFOAM_SRC}/thermoTools/lnInclude 
    ${OpenFOAM_SRC}/thermophysicalModels/thermophysicalProperties/lnInclude 
    ${OpenFOAM_SRC}/fvOptions/lnInclude
)

target_link_directories(MassTransferCoupled PRIVATE 
    ${OpenFOAM_LIB}
)

target_link_libraries(MassTransferCoupled PRIVATE
    compressibleTransportModels 
    fluidThermophysicalModels 
    solidThermo 
    solidSpecie 
    thermoTools 
    turbulenceModels 
    specie 
    finiteVolume 
    fvOptions 
    meshTools
)

install(TARGETS MassTransferCoupled DESTINATION ${OpenFOAM_USR_LIB})