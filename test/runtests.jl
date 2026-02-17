using TestItems

@testsnippet Setup begin
    using Artifacts
    using DataFrames
    using CSV
    using MassProps

    data_dir = artifact"mass_props_data"

    read_data(filename) = CSV.read(joinpath(data_dir, filename), DataFrame; missingstring = "NA",delim = "\t")

    mp_table = read_data("mp_table.tsv")
    mp_edges = CSV.read(joinpath(data_dir, "mp_edges.tsv"), DataFrame; delim = "\t")
    
    sawe_table = CSV.read(joinpath(data_dir, "sawe_table.tsv"), DataFrame; delim = "\t")
    sawe_edges = CSV.read(joinpath(data_dir, "sawe_edges.tsv"), DataFrame; delim = "\t")
    
    test_table = CSV.read(joinpath(data_dir, "test_table.tsv"), DataFrame; delim = "\t")
    test_edges = CSV.read(joinpath(data_dir, "test_edges.tsv"), DataFrame; delim = "\t")

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
