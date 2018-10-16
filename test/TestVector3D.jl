push!(LOAD_PATH, "../src/")
import UNSflow
import Test

let
    vec1 = UNSflow.Vector3D(1,1,1)
    vecm1 = UNSflow.Vector3D(-1,-1,-1)
    vec1b = UNSflow.Vector3D(1,1,1)
    vec0 = UNSflow.Vector3D(0,0,0)
    vec0b = UNSflow.Vector3D(0,0,0)
    vec2 = UNSflow.Vector3D(2,2,2)
    vecx = UNSflow.Vector3D(1,0,0)
    vecxb = UNSflow.Vector3D(1,0,0)
    vecy = UNSflow.Vector3D(0,1,0)
    vecz = UNSflow.Vector3D(0,0,1)

    # Comparisons
    Test.@test vec1 == vec1
    Test.@test vec1 == vec1b
    Test.@test vec0 == vec0b
    Test.@test vec1 != vec0
    Test.@test vecx == vecxb
    Test.@test vecy != vecx
    Test.@test vecx != vec0
    Test.@test vecx != vec1

    # Conversions
    Test.@test convert(Vector{Float64}, vec1) == [1.,1.,1.]
    Test.@test convert(Vector{Float64}, vecy) == [0.,1.,0.]
    Test.@test convert(UNSflow.Vector3D, [1,1,1]) == vec1
    Test.@test convert(UNSflow.Vector3D, [1,0,0]) == vecx
    Test.@test convert(Vector{Float64}, [vec1, vec0]) == [1.,1.,1.,0.,0.,0.]
    Test.@test convert(Vector{UNSflow.Vector3D}, [1,0,0,0,1,0]) == [vecx, vecy]
    Test.@test convert(Matrix{Float64}, [vec1, vec0]) == [1 0; 1 0; 1 0]
    Test.@test convert(Vector{UNSflow.Vector3D}, [1 0;1 0;1 0]) == [vec1, vec0]

    # Arithmatic

    Test.@test vec0 + vec0 == vec0
    Test.@test vec1 + vec0 == vec1
    Test.@test vecx + vecy + vecz == vec1

    Test.@test vecx - vecx == vec0
    Test.@test vec1 - vecx - vecy - vecz == vec0

    Test.@test -vec1 == vecm1

    Test.@test vec0 * 10 == vec0
    Test.@test vec1 * 2 == vec2
    Test.@test 2 * vec1 == vec2
    Test.@test [1 0 0; 0 1 0; 0 0 1] * vec1 == vec1
    Test.@test [1 1 1; 1 1 1; 1 1 1] * vecx == vec1

    Test.@test vec2 / 2 == vec1

    Test.@test UNSflow.cross(vecx, vecy) == vecz

    Test.@test UNSflow.dot(vecx, vec1) == 1
    Test.@test UNSflow.dot(vecx, vec0) == 0
    Test.@test UNSflow.dot(vec1, vec0) == 0
    Test.@test UNSflow.dot(vecx, vec2) == 2
    Test.@test UNSflow.dot(vecx, vecy) == 0
    Test.@test UNSflow.dot(vecy, [1 0; 0 1; 1 0]) == [0, 1]
    Test.@test UNSflow.dot(vec1 - vecy, [1 0; 0 1; 1 0]) == [2, 0]
    Test.@test UNSflow.dot(vecx, [1 0; 0 0; 1 0]) == [1, 0]
    Test.@test UNSflow.dot([1 0 1; 0 1 0], vecy) == [0, 1]
    Test.@test UNSflow.dot([1 0 1; 0 1 0], vec1 - vecy) == [2, 0]
    Test.@test UNSflow.dot([1 0 1; 0 0 0], vecx) == [1, 0]

    Test.@test abs(vecx) == 1
    Test.@test abs(vec1) == sqrt(3)

    tmp = UNSflow.Vector3D(2,3,4)
    UNSflow.zero!(tmp)
    Test.@test tmp == vec0

    Test.@test UNSflow.unit(vecx) == vecx
    Test.@test UNSflow.unit(vec1) == vec1 / sqrt(3)
    
    tmp = UNSflow.Vector3D(1,1,1)
    UNSflow.unit!(tmp)
    Test.@test tmp == vec1 / sqrt(3)

    Test.@test UNSflow.iszero(vec1) == false
    Test.@test UNSflow.iszero(vec0) == true

    Test.@test vec1[1] == 1
    Test.@test vecy[1] == 0
    Test.@test vecy[2] == 1
    Test.@test vecy[3] == 0

    Test.@test size(vecy) == (3,)
    Test.@test length(vecy) == 3

end
