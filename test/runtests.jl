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
end

@testitem "get_mass_props() for positive convention" setup = [Setup] begin

    result = MassProps.get_mass_props(mp_table, "C.1.2.2.3.1.2.3")
    row = mp_table[findfirst(mp_table.id .== "C.1.2.2.3.1.2.3"), :]

    @test result isa NamedTuple

    @test result.mass == row.mass

    @test result.center_mass isa Vector{Float64}
    @test result.center_mass == [row.Cx, row.Cy, row.Cz]

    @test result.poi_conv == row.POIconv
    @test result.poi_conv == "+"

    @test result.inertia isa Vector{Float64}
    @test result.inertia == [ row.Ixx, -row.Ixy, -row.Ixz,
                             -row.Ixy,  row.Iyy, -row.Iyz,
                             -row.Ixz, -row.Iyz,  row.Izz
                            ]

    
                                
end
    
@testitem "get_mass_props() for negative convention" setup = [Setup] begin

    result = MassProps.get_mass_props(mp_table, "C.1.2.2.3.2.1.1")
    row = mp_table[findfirst(mp_table.id .== "C.1.2.2.3.2.1.1"), :]

    @test result isa NamedTuple

    @test result.mass == row.mass

    @test result.center_mass isa Vector{Float64}
    @test result.center_mass == [row.Cx, row.Cy, row.Cz]

    @test result.poi_conv == row.POIconv
    @test result.poi_conv == "-"

    @test result.inertia isa Vector{Float64}
    @test result.inertia == [row.Ixx, row.Ixy, row.Ixz,
                             row.Ixy, row.Iyy, row.Iyz,
                             row.Ixz, row.Iyz, row.Izz
                            ]

    
                                
end
