__precompile__(true)

module UNSflow

import Dierckx: Spline1D, derivative, evaluate

import ForwardDiff

import NLsolve: nlsolve, not_in_place

import DelimitedFiles

import Serialization

import PyCall
PyCall.pygui(:tk)

import PyPlot: plot, scatter, figure, xlabel, ylabel, xlim, ylim,
xticks, yticks, subplot, legend, axis, savefig, close

import SpecialFunctions: erf

import LaTeXStrings: @L_str

export
    # kinematics types and funtions
    MotionDef,
    KinemPar,
    KinemDef,
    EldUpDef,
    EldUptstartDef,
    ConstDef,
    EldRampReturnDef,
    EldUpIntDef,
    EldUpInttstartDef,
    SinDef,
    CosDef,
    StepGustDef,

    # 2D low-order solver types
    TwoDSurf,
    TwoDOFPar,
    KinemPar2DOF,
    TwoDSurf2DOF,
    TwoDVort,
    TwoDFlowField,
    KelvinCondition,
    KelvinCondition2DOF,
    KelvinKutta,
    KelvinKutta2DOF,

    # vortex count control utility
    delVortDef,
    delNone,
    delSpalart,

    # 3D low-order solver types
    ThreeDSurfSimple,
    KinemDef3D,
    ThreeDFieldSimple,
    KelvinConditionLLT,

    # utility functions
    simpleTrapz,
    camber_calc,
    find_tstep,
    simpleInterp,
    cleanWrite,

    # 2D low-order solver function
    update_boundpos,
    update_kinem,
    update_indbound,
    update_downwash,
    update_a0anda1,
    place_tev,
    place_lev,
    update_a2toan,
    mutual_ind,
    update_a2a3adot,
    update_bv,
    ind_vel,
    wakeroll,
    update_adot,
    update_externalvel,
    controlVortCount,
    update_kinem,

    #2D low-order solver methods
    lautat,
    lautatRoll,
    ldvm,
    ldvmLin,

    # Postprocessing functions
    calc_forces,
    writeStamp,

    # 2D plotting functions
    viewVort2D,
    viewVortConnect2D,

    # 2D plot output functions
    makeForcePlots,
    makeVortPlots2D,

    # 3D low-order solver methods
    QSLLTlautat,
    QSLLTlautatRoll,
    QSLLTldvm,
    QSLLTlautatLin,
    QSLLTldvmLin,
    StripldvmLin,

    # 2D low-order solver function
    calc_a0a13d,
    calc_a2toan3d,

    # 3D plot output functions
    makeForcePlots3Dstrip,
    makeVortPlots3Dstrip,
    makeTevstrPlots3Dstrip,

    # 3D Geometry
    DiscreteGeometry3D,
    Line2,
    Point3D,
    PolyLine2,
    Vector3D,

    # 3D Vorticity
    Vorticity3D,
    Vorticity3DCollector,
    Vorticity3DSimpleCollector,
    Vorticity3DAdaptive,
    VortexParticleFilamentAdaptive,
    VortexParticleVolumeAdaptive,
    VortexParticle3D,
    VortexRing,
    StraightVortexFilament

### source files

# kinematic types
include("kinem.jl")

# utility functions
include("utils.jl")

# vortex count control utility
include("delVort.jl")

# low-order 2D solvers
include("lowOrder2D/typedefs.jl")            # type definitions
include("lowOrder2D/calcs.jl")               # calculation functions
include("lowOrder2D/solvers.jl")             # solver methods
include("lowOrder2D/postprocess.jl")         # postprocessing functions

# low-order 3D solvers
# include("lowOrder3D/typedefs.jl")            # type definitions
# include("lowOrder3D/calcs.jl")               # calculation functions
# include("lowOrder3D/solvers.jl")             # solver methods
# include("lowOrder3D/postprocess.jl")         # postprocessing functions
include("lowOrder3D/ThreadWork.jl")
include("lowOrder3D/Vector3D.jl")
include("lowOrder3D/DiscreteGeometry3D.jl")
include("lowOrder3D/Point3D.jl")
include("lowOrder3D/Line2.jl")
include("lowOrder3D/PolyLine2.jl")
include("lowOrder3D/DiscreteGeometry3DToVTK.jl")
include("lowOrder3D/Vortex3DRegularisationFunctions.jl")
include("lowOrder3D/RedistributionScheme.jl")
include("lowOrder3D/Vorticity3D.jl")
include("lowOrder3D/Vorticity3DCollector.jl")
include("lowOrder3D/Vorticity3DSimpleCollector.jl")
include("lowOrder3D/Vorticity3DAdaptive.jl")
include("lowOrder3D/StraightVortexFilament.jl")
include("lowOrder3D/VortexParticle3D.jl")
include("lowOrder3D/VortexRing.jl")
include("lowOrder3D/SplineVortexParticleFilament.jl")
include("lowOrder3D/VortexParticleFilamentAdaptive.jl")
include("lowOrder3D/VortexParticleVolumeAdaptive.jl")

# 2D plotting functions
include("plots/plots2D.jl")

end
