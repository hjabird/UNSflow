#===============================================================================
    ThreeDStraightVortexFilament

    Initial code: HJAB 2018
------------------------------------------------------------------------------=#
type ThreeDStraightVortexFilament
    start_coord :: ThreeDVector
    end_coord :: ThreeDVector
    vorticity :: Float64

    function ThreeDStraightVortexFilament(
        start_coord :: ThreeDVector,
        end_coord :: ThreeDVector,
        vorticity :: T
        ) where T <: Real
        new(start_coord, end_coord, Float64(vorticity))
    end
end

function ThreeDStraightVortexFilament(
    start_coord :: ThreeDVector,
    end_coord :: ThreeDVector
    )
    return ThreeDStraightVortexFilament(start_coord, end_coord, 0.0)
end

ThreeDStraightVortexFilament() =
    ThreeDStraightVortexFilament(
        ThreeDVector([0., 0., 0.]),
        ThreeDVector([1., 1., 1.]),
        0.0)
#= END ThreeDStraightVortexFilament ------------------------------------------=#
