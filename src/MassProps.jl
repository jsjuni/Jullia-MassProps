module MassProps

    using RollupTree

    get_mass_props(table, id) = begin
        row = RollupTree.df_get_row_by_id(table, id)
        poi_factor = row.POIconv == "-" ? 1.0 : -1.0
 
        (
            mass = row.mass,

            center_mass = [row.Cx, row.Cy, row.Cz],

            poi_conv = row.POIconv,

            inertia = [             row.Ixx, poi_factor * row.Ixy, poi_factor * row.Ixz,
                       poi_factor * row.Ixy,        row.Iyy, poi_factor * row.Iyz,
                       poi_factor * row.Ixz, poi_factor * row.Iyz,              row.Izz
            ],

            point = row.Ipoint
        )
    end

    get_mass_props_unc(table, id) = begin
        row = RollupTree.df_get_row_by_id(table, id)

        (
            sigma_mass = row.sigma_mass,

            sigma_center_mass = [row.sigma_Cx, row.sigma_Cy, row.sigma_Cz],

            sigma_inertia = [row.sigma_Ixx, row.sigma_Ixy, row.sigma_Ixz,
                             row.sigma_Ixy, row.sigma_Iyy, row.sigma_Iyz,
                             row.sigma_Ixz, row.sigma_Iyz, row.sigma_Izz
            ]
        )
    end

    get_mass_props_and_unc(table, id) = merge(get_mass_props(table, id), get_mass_props_unc(table, id))
    
end
