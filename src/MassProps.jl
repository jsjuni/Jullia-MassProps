module MassProps

    using RollupTree
    using LinearAlgebra

    get_mass_props(table, id) = begin
        row = RollupTree.df_get_row_by_id(table, id)
        poi_factor = row.POIconv == "-" ? 1.0 : -1.0
 
        (
            mass = row.mass,

            center_mass = [row.Cx; row.Cy; row.Cz],

            inertia = [             row.Ixx poi_factor * row.Ixy poi_factor * row.Ixz;
                       poi_factor * row.Ixy              row.Iyy poi_factor * row.Iyz;
                       poi_factor * row.Ixz poi_factor * row.Iyz              row.Izz
            ],

            poi_conv = row.POIconv,
            point = row.Ipoint
        )

    end

    get_mass_props_unc(table, id) = begin
        row = RollupTree.df_get_row_by_id(table, id)

        (
            sigma_mass = row.sigma_mass,

            sigma_center_mass = [row.sigma_Cx; row.sigma_Cy; row.sigma_Cz],

            sigma_inertia = [row.sigma_Ixx row.sigma_Ixy row.sigma_Ixz;
                             row.sigma_Ixy row.sigma_Iyy row.sigma_Iyz;
                             row.sigma_Ixz row.sigma_Iyz row.sigma_Izz
            ]
        )
    end

    get_mass_props_and_unc(table, id) = merge(get_mass_props(table, id), get_mass_props_unc(table, id))

    set_mass_props(table, id, mp) = begin
        
        cm = mp.center_mass
        it = (mp.inertia + transpose(mp.inertia)) / 2

        poi_factor = mp.poi_conv == "+" ? -1.0 : (mp.poi_conv == "-" ? 1.0 : error("invalid POI convention"))

        values = (
            mass = mp.mass,

            Cx = cm[1],
            Cy = cm[2],
            Cz = cm[3],

            Ixx = it[1, 1],
            Iyy = it[2, 2],
            Izz = it[3, 3],

            Ixy = poi_factor * it[1, 2],
            Ixz = poi_factor * it[1, 3],
            Iyz = poi_factor * it[2, 3],

            Ipoint = mp.point,
            POIconv = mp.poi_conv
        )

        RollupTree.df_set_row_by_id(table, id, values)
    end

    set_mass_props_unc(table, id, mp_unc) = begin
        values = (
            sigma_mass = mp_unc.sigma_mass,

            sigma_Cx = mp_unc.sigma_center_mass[1],
            sigma_Cy = mp_unc.sigma_center_mass[2],
            sigma_Cz = mp_unc.sigma_center_mass[3],

            sigma_Ixx = mp_unc.sigma_inertia[1, 1],
            sigma_Iyy = mp_unc.sigma_inertia[2, 2],
            sigma_Izz = mp_unc.sigma_inertia[3, 3],

            sigma_Ixy = mp_unc.sigma_inertia[1, 2],
            sigma_Ixz = mp_unc.sigma_inertia[1, 3],
            sigma_Iyz = mp_unc.sigma_inertia[2, 3]
        )

        RollupTree.df_set_row_by_id(table, id, values)
    end

    set_mass_props_and_unc(table, id, mpu) = set_mass_props_unc(set_mass_props(table, id, mpu), id, mpu)

    combine_mass_props(mpl) = begin
        
        mass = sum(mp.mass for mp in mpl)

        center_mass = sum(mp.mass .* mp.center_mass for mp in mpl) ./ mass

        inertia = sum(mp.inertia .+ mp.mass .* (dot(mp.center_mass, mp.center_mass) * I - mp.center_mass * transpose(mp.center_mass)) for mp in mpl)
        inertia = sum(
            map(mp -> begin
                d = mp.center_mass - center_mass
                mp.inertia .+ mp.mass .* (dot(d, d) * I - d * transpose(d))
            end, mpl)
        )

        (
            mass = mass,
            center_mass = center_mass,
            inertia = inertia,
            point = false
        )
    end

    combine_mass_props_unc(mpul, amp) = begin

        sigma_mass = sqrt(sum(mpu.sigma_mass^2 for mpu in mpul))

        sigma_center_mass = sqrt.(
            sum(
                (
                    (mpu.mass .* mpu.sigma_center_mass).^2 .+
                    (mpu.sigma_mass .* (mpu.center_mass - amp.center_mass)).^2
                ) for mpu in mpul
            )
        ) ./ amp.mass

        sigma_inertia = sqrt.(
            sum(
                map(mpu -> begin
                    d = mpu.center_mass - amp.center_mass
                    P = d .* mpu.sigma_center_mass'
                    p = diag(P)
                    Q = d .* d'

                    M1 = P  - diagm(p - 2 .* [p[2]; p[1]; p[1]])
                    M2 = P' - diagm(p - 2 .* [p[3]; p[3]; p[2]])
                    M3 = Q  - tr(Q) * I
                    M4 = mpu.mass^2 .* (M1.^2 .+ M2.^2) .+ (mpu.sigma_mass .* M3).^2
                    mpu.point ? M4 : mpu.sigma_inertia.^2 .+ M4
                end, mpul)
            )
        )

        merge(amp,
            (
                sigma_mass = sigma_mass,
                sigma_center_mass = sigma_center_mass,
                sigma_inertia = sigma_inertia
            )
        )

    end

    combine_mass_props_and_unc(mpul) = combine_mass_props_unc(mpul, combine_mass_props(mpul))

    set_poi_conv_plus(df, target, mp) = merge(mp, (poi_conv = "+",))
    set_poi_conv_minus(df, target, mp) = merge(mp, (poi_conv = "-",))
    set_poi_conv_from_target(df, target, mp) = merge(mp, (poi_conv = RollupTree.df_get_by_id(df, target, :POIconv),))

    update_mass_props(df, target, sources, override = set_poi_conv_from_target) = begin
        RollupTree.update_prop(
            df,
            target,
            sources,
            set_mass_props,
            get_mass_props,
            combine_mass_props,
            override
        )
    end

    update_mass_props_unc(df, target, sources) = begin
        RollupTree.update_prop(
            df,
            target,
            sources,
            set_mass_props_unc,
            get_mass_props_and_unc,
            l -> combine_mass_props_unc(l, get_mass_props(df, target))
        )
    end

    update_mass_props_and_unc(df, target, sources, override = set_poi_conv_from_target) = begin
           RollupTree.update_prop(
            df,
            target,
            sources,
            set_mass_props_and_unc,
            get_mass_props_and_unc,
            combine_mass_props_and_unc,
            override
        )
    end

    validate_mass_props(mp) = begin
        
        ismissing(mp.mass) && error("mass is missing")
        isnothing(mp.mass) && error("mass is nothing")
        mp.mass isa Real || error("mass must be a real number")
        mp.mass > 0.0 || error("mass must be positive")

        any(ismissing.(mp.center_mass)) && error("center mass component is missing")
        any(isnothing.(mp.center_mass)) && error("center mass component is nothing")
        mp.center_mass isa AbstractVector{<:Real} || error("center mass must be a vector of real numbers")
        length(mp.center_mass) == 3 || error("center mass must have three components")

        any(ismissing.(mp.inertia)) && error("inertia component is missing")
        any(isnothing.(mp.inertia)) && error("inertia component is nothing")
        mp.inertia isa AbstractMatrix{<:Real} || error("inertia must be a matrix of real numbers")
        size(mp.inertia) == (3, 3) || error("inertia must be a 3x3 matrix")

        ev = eigen(mp.inertia)
        all(isreal.(ev.values)) || error("inertia matrix must have real eigenvalues")
        any(ev.values .<= 0.0) && error("inertia matrix must be positive definite")
        all([
            ev.values[1] <= ev.values[2] + ev.values[3],
            ev.values[2] <= ev.values[1] + ev.values[3],
            ev.values[3] <= ev.values[1] + ev.values[2]
        ]) || error("inertia matrix must satisfy triangle inequalities")

        mp.poi_conv in ("+", "-") || error("POI convention must be '+' or '-'")

        mp.point isa Bool || error("point must be a boolean")

        return true

    end

    validate_mass_props_unc(mpu) = true

    validate_mass_props_and_unc(mpu) = validate_mass_props(mpu) && validate_mass_props_unc(mpu)

    validate_mass_props_table(tree, df) = RollupTree.validate_ds(tree, df, RollupTree.df_get_ids, get_mass_props, validate_mass_props)

    validate_mass_props_and_unc_table(tree, df) = RollupTree.validate_ds(tree, df, RollupTree.df_get_ids, get_mass_props_and_unc, validate_mass_props_and_unc)

    rollup_mass_props(tree, df, validate_df = validate_mass_props_table) = RollupTree.rollup(tree, df, update_mass_props, validate_df)

    rollup_mass_props_unc(tree, df, validate_df = validate_mass_props_and_unc_table) = RollupTree.rollup(tree, df, update_mass_props_unc, validate_df)

    rollup_mass_props_and_unc(tree, df, validate_df = validate_mass_props_and_unc_table) = RollupTree.rollup(tree, df, update_mass_props_and_unc, validate_df)

end