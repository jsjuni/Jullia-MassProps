using TestItems

@testsnippet Setup begin
    using Artifacts
    using DataFrames
    using CSV
    using MassProps
    using MetaGraphsNext
    using Graphs

    data_dir = artifact"mass_props_data"

    read_data(filename) = CSV.read(joinpath(data_dir, filename), DataFrame; missingstring = "NA",delim = "\t")

    tree_from_edgelist(el, id, pid) = begin
        tree = MetaGraphsNext.MetaGraph(Graphs.SimpleDiGraph(), label_type = String)
        for row in eachrow(el)
            Graphs.add_vertex!(tree, row[id])
            if !ismissing(row[pid])
                Graphs.add_vertex!(tree, row[pid])
                Graphs.add_edge!(tree, row[id], row[pid])
            end
        end
        tree
    end

    mp_table = read_data("mp_table.tsv")
    mp_edges = read_data("mp_edges.tsv")
    mp_tree = tree_from_edgelist(mp_edges, "child", "parent")
    
    sawe_table = read_data("sawe_table.tsv")
    sawe_edges = read_data("sawe_edges.tsv")
    sawe_tree = tree_from_edgelist(sawe_edges, "child", "parent")
    
    test_table = read_data("test_table.tsv")
    test_edges = read_data("test_edges.tsv")
    test_tree = tree_from_edgelist(test_edges, "id", "parent")

    id_pos = "C.1.2.2.3.1.2.3"
    row_pos = mp_table[findfirst(mp_table.id .== id_pos), :]
    mass_pos = row_pos.mass
    center_mass_pos = [row_pos.Cx; row_pos.Cy; row_pos.Cz]
    inertia_pos = [  row_pos.Ixx -row_pos.Ixy -row_pos.Ixz;
                    -row_pos.Ixy  row_pos.Iyy -row_pos.Iyz;
                    -row_pos.Ixz -row_pos.Iyz  row_pos.Izz
                    ]
    poi_conv_pos = row_pos.POIconv
    point_pos = row_pos.Ipoint

    sigma_mass_pos = row_pos.sigma_mass
    sigma_center_mass_pos = [row_pos.sigma_Cx; row_pos.sigma_Cy; row_pos.sigma_Cz]
    sigma_inertia_pos = [ row_pos.sigma_Ixx row_pos.sigma_Ixy row_pos.sigma_Ixz;
                          row_pos.sigma_Ixy row_pos.sigma_Iyy row_pos.sigma_Iyz;
                          row_pos.sigma_Ixz row_pos.sigma_Iyz row_pos.sigma_Izz
                        ]
    mp_pos = (
        mass = mass_pos,
        center_mass = center_mass_pos,
        inertia = inertia_pos,
        poi_conv = poi_conv_pos,
        point = point_pos
    )

    mpu_pos = (
        sigma_mass = sigma_mass_pos,
        sigma_center_mass = sigma_center_mass_pos,
        sigma_inertia = sigma_inertia_pos
    )

    id_neg = "C.1.2.2.3.2.1.1"
    row_neg = mp_table[findfirst(mp_table.id .== id_neg), :]
    mass_neg = row_neg.mass
    center_mass_neg = [row_neg.Cx; row_neg.Cy; row_neg.Cz]
    inertia_neg = [  row_neg.Ixx row_neg.Ixy row_neg.Ixz;
                     row_neg.Ixy row_neg.Iyy row_neg.Iyz;
                     row_neg.Ixz row_neg.Iyz row_neg.Izz
                    ]
    poi_conv_neg = row_neg.POIconv
    point_neg = row_neg.Ipoint

    sigma_mass_neg = row_neg.sigma_mass
    sigma_center_mass_neg = [row_neg.sigma_Cx; row_neg.sigma_Cy; row_neg.sigma_Cz]
    sigma_inertia_neg = [row_neg.sigma_Ixx row_neg.sigma_Ixy row_neg.sigma_Ixz;
                          row_neg.sigma_Ixy row_neg.sigma_Iyy row_neg.sigma_Iyz;
                          row_neg.sigma_Ixz row_neg.sigma_Iyz row_neg.sigma_Izz
                        ]

     mp_neg = (
        mass = mass_neg,
        center_mass = center_mass_neg,
        inertia = inertia_neg,
        poi_conv = poi_conv_neg,
        point = point_neg
    )

   mpu_neg = (
        sigma_mass = sigma_mass_neg,
        sigma_center_mass = sigma_center_mass_neg,
        sigma_inertia = sigma_inertia_neg
    )

    target_id = "C.1"
end

@testitem "get_mass_props() for positive convention" setup = [Setup] begin

    result = MassProps.get_mass_props(mp_table, id_pos)

    @test result isa NamedTuple

    @test result.mass == mass_pos

    @test result.center_mass isa Vector{Float64}
    @test result.center_mass == center_mass_pos

    @test result.poi_conv == poi_conv_pos
    @test result.poi_conv == "+"

    @test result.point == point_pos

    @test result.inertia isa Matrix{Float64}
    @test result.inertia == inertia_pos

end
    
@testitem "get_mass_props() for negative convention" setup = [Setup] begin

    result = MassProps.get_mass_props(mp_table, id_neg)
 
    @test result isa NamedTuple

    @test result.mass == mass_neg

    @test result.center_mass isa Vector{Float64}
    @test result.center_mass == center_mass_neg

    @test result.poi_conv == poi_conv_neg
    @test result.poi_conv == "-"

    @test result.point == point_neg

    @test result.inertia isa Matrix{Float64}
    @test result.inertia == inertia_neg

end

@testitem "get_mass_props_unc()" setup = [Setup] begin

    result = MassProps.get_mass_props_unc(mp_table, id_pos)
 
    @test result isa NamedTuple

    @test result.sigma_mass == sigma_mass_pos
    
    @test result.sigma_center_mass isa Vector{Float64}
    @test result.sigma_center_mass == sigma_center_mass_pos   

    @test result.sigma_inertia isa Matrix{Float64}
    @test result.sigma_inertia == sigma_inertia_pos
end

@testitem "get_mass_props_and_unc() for positive convention" setup = [Setup] begin

    result = MassProps.get_mass_props_and_unc(mp_table, id_pos)
    mp = MassProps.get_mass_props(mp_table, id_pos)
    uc = MassProps.get_mass_props_unc(mp_table, id_pos)
    
    @test result isa NamedTuple

    @test result == merge(mp, uc)

end

@testitem "get_mass_props_and_unc() for negative convention" setup = [Setup] begin

    result = MassProps.get_mass_props_and_unc(mp_table, id_neg)
    mp = MassProps.get_mass_props(mp_table, id_neg)
    uc = MassProps.get_mass_props_unc(mp_table, id_neg)
    
    @test result isa NamedTuple

    @test result == merge(mp, uc)

end

@testitem "set_mass_props() for positive convention" setup = [Setup] begin
    ot = MassProps.set_mass_props(mp_table, target_id, mp_pos)
    row = ot[findfirst(ot.id .== target_id), :]

    @test row.mass == mp_pos.mass

    @test row.Cx == mp_pos.center_mass[1]
    @test row.Cy == mp_pos.center_mass[2]
    @test row.Cz == mp_pos.center_mass[3]

    @test row.Ixx == mp_pos.inertia[1, 1]
    @test row.Iyy == mp_pos.inertia[2, 2]
    @test row.Izz == mp_pos.inertia[3, 3]

    @test row.Ixy == -mp_pos.inertia[1, 2]
    @test row.Ixz == -mp_pos.inertia[1, 3]
    @test row.Iyz == -mp_pos.inertia[2, 3]

    @test row.Ipoint == mp_pos.point
    @test row.POIconv == mp_pos.poi_conv
end

@testitem "set_mass_props() for negative convention" setup = [Setup] begin
    ot = MassProps.set_mass_props(mp_table, target_id, mp_neg)
    row = ot[findfirst(ot.id .== target_id), :]

    @test row.mass == mp_neg.mass

    @test row.Cx == mp_neg.center_mass[1]
    @test row.Cy == mp_neg.center_mass[2]
    @test row.Cz == mp_neg.center_mass[3]

    @test row.Ixx == mp_neg.inertia[1, 1]
    @test row.Iyy == mp_neg.inertia[2, 2]
    @test row.Izz == mp_neg.inertia[3, 3]

    @test row.Ixy == mp_neg.inertia[1, 2]
    @test row.Ixz == mp_neg.inertia[1, 3]
    @test row.Iyz == mp_neg.inertia[2, 3]

    @test row.Ipoint == mp_neg.point
    @test row.POIconv == mp_neg.poi_conv
end

@testitem "set_mass_props() with invalid POI convention" setup = [Setup] begin
    mp_invalid = merge(mp_pos, (poi_conv = "invalid",))
    @test_throws ErrorException MassProps.set_mass_props(mp_table, target_id, mp_invalid)
end

@testitem "set_mass_props_unc()" setup = [Setup] begin

    ot = MassProps.set_mass_props_unc(mp_table, target_id, mpu_pos)
    row = ot[findfirst(ot.id .== target_id), :]

    @test row.sigma_mass == mpu_pos.sigma_mass

    @test row.sigma_Cx == mpu_pos.sigma_center_mass[1]
    @test row.sigma_Cy == mpu_pos.sigma_center_mass[2]
    @test row.sigma_Cz == mpu_pos.sigma_center_mass[3]

    @test row.sigma_Ixx == mpu_pos.sigma_inertia[1, 1]
    @test row.sigma_Iyy == mpu_pos.sigma_inertia[2, 2]
    @test row.sigma_Izz == mpu_pos.sigma_inertia[3, 3]

    @test row.sigma_Ixy == mpu_pos.sigma_inertia[1, 2]
    @test row.sigma_Ixz == mpu_pos.sigma_inertia[1, 3]
    @test row.sigma_Iyz == mpu_pos.sigma_inertia[2, 3]

end

@testitem "set_mass_props_and_unc() for positive convention" setup = [Setup] begin

    ot = MassProps.set_mass_props_and_unc(mp_table, target_id, merge(mp_pos, mpu_pos))
    row = ot[findfirst(ot.id .== target_id), :]

    @test row.mass == mp_pos.mass

    @test row.Cx == mp_pos.center_mass[1]
    @test row.Cy == mp_pos.center_mass[2]
    @test row.Cz == mp_pos.center_mass[3]

    @test row.Ixx == mp_pos.inertia[1, 1]
    @test row.Iyy == mp_pos.inertia[2, 2]
    @test row.Izz == mp_pos.inertia[3, 3]

    @test row.Ixy == -mp_pos.inertia[1, 2]
    @test row.Ixz == -mp_pos.inertia[1, 3]
    @test row.Iyz == -mp_pos.inertia[2, 3]

    @test row.Ipoint == mp_pos.point
    @test row.POIconv == mp_pos.poi_conv

    @test row.sigma_mass == mpu_pos.sigma_mass

    @test row.sigma_Cx == mpu_pos.sigma_center_mass[1]
    @test row.sigma_Cy == mpu_pos.sigma_center_mass[2]
    @test row.sigma_Cz == mpu_pos.sigma_center_mass[3]

    @test row.sigma_Ixx == mpu_pos.sigma_inertia[1, 1]
    @test row.sigma_Iyy == mpu_pos.sigma_inertia[2, 2]
    @test row.sigma_Izz == mpu_pos.sigma_inertia[3, 3]

    @test row.sigma_Ixy == mpu_pos.sigma_inertia[1, 2]
    @test row.sigma_Ixz == mpu_pos.sigma_inertia[1, 3]
    @test row.sigma_Iyz == mpu_pos.sigma_inertia[2, 3]

end

@testitem "set_mass_props_and_unc() for negative convention" setup = [Setup] begin

    ot = MassProps.set_mass_props_and_unc(mp_table, target_id, merge(mp_neg, mpu_neg))
    row = ot[findfirst(ot.id .== target_id), :]

    @test row.mass == mp_neg.mass

    @test row.Cx == mp_neg.center_mass[1]
    @test row.Cy == mp_neg.center_mass[2]
    @test row.Cz == mp_neg.center_mass[3]

    @test row.Ixx == mp_neg.inertia[1, 1]
    @test row.Iyy == mp_neg.inertia[2, 2]
    @test row.Izz == mp_neg.inertia[3, 3]

    @test row.Ixy == mp_neg.inertia[1, 2]
    @test row.Ixz == mp_neg.inertia[1, 3]
    @test row.Iyz == mp_neg.inertia[2, 3]

    @test row.Ipoint == mp_neg.point
    @test row.POIconv == mp_neg.poi_conv

    @test row.sigma_mass == mpu_neg.sigma_mass

    @test row.sigma_Cx == mpu_neg.sigma_center_mass[1]
    @test row.sigma_Cy == mpu_neg.sigma_center_mass[2]
    @test row.sigma_Cz == mpu_neg.sigma_center_mass[3]

    @test row.sigma_Ixx == mpu_neg.sigma_inertia[1, 1]
    @test row.sigma_Iyy == mpu_neg.sigma_inertia[2, 2]
    @test row.sigma_Izz == mpu_neg.sigma_inertia[3, 3]

    @test row.sigma_Ixy == mpu_neg.sigma_inertia[1, 2]
    @test row.sigma_Ixz == mpu_neg.sigma_inertia[1, 3]
    @test row.sigma_Iyz == mpu_neg.sigma_inertia[2, 3]

end

@testitem "combine_mass_props() for non-point masses" setup = [Setup] begin

    leaves = collect(test_table[map(!ismissing, test_table[:, :mass]), :id])
 
    mpl = map(id -> MassProps.get_mass_props(test_table, id), leaves)

    mpc = MassProps.combine_mass_props(mpl)

    @test mpc.mass == 21.0
    @test mpc.center_mass == [0.0, 0.0, 0.0]
    @test mpc.inertia ≈ [144.0 -4.8 -24.8; -4.8 144.0 -23.2; -24.8 -23.2 139.0]

    @test mpc.point == false
    
end

@testitem "combine_mass_props_unc() for non-point masses" setup = [Setup] begin

    amp = MassProps.get_mass_props_and_unc(sawe_table, "Combined")

    leaves = MetaGraphsNext.inneighbor_labels(sawe_tree, "Combined")
    mpul = map(id -> MassProps.get_mass_props_and_unc(sawe_table, id), leaves)
    
    mpuc = MassProps.combine_mass_props_unc(mpul, amp)

    @test isapprox(mpuc.sigma_mass, amp.sigma_mass, rtol = 1e-5)
    @test isapprox(mpuc.sigma_center_mass, amp.sigma_center_mass, rtol = 5e-3)
    @test isapprox(mpuc.sigma_inertia, amp.sigma_inertia, rtol = 2.1e-3) # published values are not accurate

end

@testitem "combine_mass_props_and_unc()" setup = [Setup] begin

    leaves = MetaGraphsNext.inneighbor_labels(sawe_tree, "Combined")
    mpul = map(id -> MassProps.get_mass_props_and_unc(sawe_table, id), leaves)

    @test isequal(MassProps.combine_mass_props_and_unc(mpul), merge(MassProps.combine_mass_props_unc(mpul, MassProps.combine_mass_props(mpul))))

end

@testitem "set_poi_conv_plus()" setup = [Setup] begin

    mp_plus = MassProps.set_poi_conv_plus(nothing, nothing, mp_neg)
    @test mp_plus.poi_conv == "+"

end

@testitem "set_poi_conv_minus()" setup = [Setup] begin

    mp_minus = MassProps.set_poi_conv_minus(nothing, nothing, mp_pos)
    @test mp_minus.poi_conv == "-"

end

@testitem "set_poi_conv_from_target()" setup = [Setup] begin

    mp_plus = MassProps.set_poi_conv_from_target(mp_table, id_pos, mp_neg)
    @test mp_plus.poi_conv == mp_pos.poi_conv
 
    mp_minus = MassProps.set_poi_conv_from_target(mp_table, id_neg, mp_pos)
    @test mp_minus.poi_conv == mp_neg.poi_conv

end

@testitem "validate_mass_props()" setup = [Setup] begin

    mp_valid = mp_pos
    @test MassProps.validate_mass_props(mp_pos) == true

    mp_invalid_mass = merge(mp_valid, (mass = -1.0,))
    # @test_throws ErrorException MassProps.validate_mass_props(mp_invalid_mass)

    mp_invalid_inertia = merge(mp_valid, (inertia = [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0],))
    # @test_throws ErrorException MassProps.validate_mass_props(mp_invalid_inertia)

    mp_invalid_poi_conv = merge(mp_valid, (poi_conv = "invalid",))
    # @test_throws ErrorException MassProps.validate_mass_props(mp_invalid_poi_conv)

end

@testitem "rollup_mass_props()" setup = [Setup] begin

    mp = MassProps.rollup_mass_props(test_tree, test_table)
    @test mp isa DataFrame

    top_row = mp[findfirst(mp.id .== "A.1"), :]

    @test top_row.mass == 21.0

    @test top_row.Cx == 0.0
    @test top_row.Cy == 0.0
    @test top_row.Cz == 0.0

    @test top_row.Ixx == 144.0
    @test top_row.Iyy == 144.0
    @test top_row.Izz == 139.0
    @test top_row.Ixy == -4.8
    @test top_row.Ixz == -24.8
    @test top_row.Iyz == -23.2

    @test top_row.POIconv == "-"

    @test top_row.Ipoint == false

end

@testitem "rollup_mass_props_and_unc()" setup = [Setup] begin

    expected = sawe_table[findfirst(sawe_table.id .== "Combined"), :]

    mp = MassProps.rollup_mass_props_and_unc(sawe_tree, sawe_table)
    actual = mp[findfirst(mp.id .== "Combined"), :]

    @test actual.mass == expected.mass

    cm_rtol = 2e-2
    @test isapprox(actual.Cx, expected.Cx, rtol = cm_rtol)
    @test isapprox(actual.Cy, expected.Cy, rtol = cm_rtol)
    @test isapprox(actual.Cz, expected.Cz, rtol = cm_rtol)

    it_rtol = 2e-2
    @test isapprox(actual.Ixx, expected.Ixx, rtol = it_rtol)
    @test isapprox(actual.Iyy, expected.Iyy, rtol = it_rtol)
    @test isapprox(actual.Izz, expected.Izz, rtol = it_rtol)
    @test isapprox(actual.Ixy, expected.Ixy, rtol = it_rtol)
    @test isapprox(actual.Ixz, expected.Ixz, rtol = it_rtol)
    @test isapprox(actual.Iyz, expected.Iyz, rtol = it_rtol)

    @test actual.POIconv == expected.POIconv
    @test actual.Ipoint == expected.Ipoint

    sigma_m_rtol = 1e-5
    @test isapprox(actual.sigma_mass, expected.sigma_mass, rtol = sigma_m_rtol)

    sigma_cm_rtol = 2e-2
    @test isapprox(actual.sigma_Cx, expected.sigma_Cx, rtol = sigma_cm_rtol)
    @test isapprox(actual.sigma_Cy, expected.sigma_Cy, rtol = sigma_cm_rtol)
    @test isapprox(actual.sigma_Cz, expected.sigma_Cz, rtol = sigma_cm_rtol)

    sigma_it_rtol = 2e-2
    @test isapprox(actual.sigma_Ixx, expected.sigma_Ixx, rtol = sigma_it_rtol)
    @test isapprox(actual.sigma_Iyy, expected.sigma_Iyy, rtol = sigma_it_rtol)
    @test isapprox(actual.sigma_Izz, expected.sigma_Izz, rtol = sigma_it_rtol)
    @test isapprox(actual.sigma_Ixy, expected.sigma_Ixy, rtol = sigma_it_rtol)
    @test isapprox(actual.sigma_Ixz, expected.sigma_Ixz, rtol = sigma_it_rtol)
    @test isapprox(actual.sigma_Iyz, expected.sigma_Iyz, rtol = sigma_it_rtol)

end