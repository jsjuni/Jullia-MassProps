module MassProps

    using RollupTree
    using LinearAlgebra

    export get_mass_props, get_mass_props_unc, get_mass_props_and_unc,
            set_mass_props, set_mass_props_unc, set_mass_props_and_unc,
            combine_mass_props, combine_mass_props_unc, combine_mass_props_and_unc,
            set_poi_conv_plus, set_poi_conv_minus, set_poi_conv_from_target,
            update_mass_props, update_mass_props_unc, update_mass_props_and_unc,
            validate_mass_props, validate_mass_props_unc, validate_mass_props_and_unc,
            validate_mass_props_table, validate_mass_props_and_unc_table,
            rollup_mass_props, rollup_mass_props_unc, rollup_mass_props_and_unc,
            X, Y, Z

    const X = 1
    const Y = 2
    const Z = 3

    get_mass_props(table, id) = begin
        row = df_get_row_by_id(table, id)
        poi_factor = row.POIconv == "-" ? 1.0 : -1.0

        (
            mass = row.mass,

            center_mass = [row.Cx; row.Cy; row.Cz],

            inertia = [             row.Ixx poi_factor * row.Ixy poi_factor * row.Ixz;
                       poi_factor * row.Ixy              row.Iyy poi_factor * row.Iyz;
                       poi_factor * row.Ixz poi_factor * row.Iyz              row.Izz
            ],

           point = row.Ipoint
        )

    end

    get_mass_props_unc(table, id) = begin
        row = df_get_row_by_id(table, id)

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

            Cx = cm[X],
            Cy = cm[Y],
            Cz = cm[Z],

            Ixx = it[X, X],
            Iyy = it[Y, Y],
            Izz = it[Z, Z],

            Ixy = poi_factor * it[X, Y],
            Ixz = poi_factor * it[X, Z],
            Iyz = poi_factor * it[Y, Z],

            Ipoint = mp.point,
            POIconv = mp.poi_conv
        )

        df_set_row_by_id(table, id, values)
    end

    set_mass_props_unc(table, id, mp_unc) = begin
        values = (
            sigma_mass = mp_unc.sigma_mass,

            sigma_Cx = mp_unc.sigma_center_mass[X],
            sigma_Cy = mp_unc.sigma_center_mass[Y],
            sigma_Cz = mp_unc.sigma_center_mass[Z],

            sigma_Ixx = mp_unc.sigma_inertia[X, X],
            sigma_Iyy = mp_unc.sigma_inertia[Y, Y],
            sigma_Izz = mp_unc.sigma_inertia[Z, Z],

            sigma_Ixy = mp_unc.sigma_inertia[X, Y],
            sigma_Ixz = mp_unc.sigma_inertia[X, Z],
            sigma_Iyz = mp_unc.sigma_inertia[Y, Z]
        )

        df_set_row_by_id(table, id, values)
    end

    set_mass_props_and_unc(table, id, mpu) = set_mass_props_unc(set_mass_props(table, id, mpu), id, mpu)

    combine_mass_props(mpl) = begin
        
        mass = sum(mp.mass for mp in mpl)

        center_mass = sum(mp.mass .* mp.center_mass for mp in mpl) ./ mass

        inertia = sum(
            map(mp -> begin
                d = mp.center_mass - center_mass
                Q = d .* d'
                M = mp.mass .* (tr(Q) * I - Q)
                mp.point ? M : mp.inertia .+ M
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

                    M1 = P  - diagm(p - 2 .* view(p, [Y, X, X]))
                    M2 = P' - diagm(p - 2 .* view(p, [Z, Z, Y]))
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
    set_poi_conv_from_target(df, target, mp) = merge(mp, (poi_conv = df_get_by_id(df, target, :POIconv),))

    update_mass_props(df, target, sources, override = set_poi_conv_from_target) = begin
        update_prop(
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
        update_prop(
            df,
            target,
            sources,
            set_mass_props_unc,
            get_mass_props_and_unc,
            l -> combine_mass_props_unc(l, get_mass_props(df, target))
        )
    end

    update_mass_props_and_unc(df, target, sources, override = set_poi_conv_from_target) = begin
           update_prop(
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

        mp.point isa Bool || error("point must be a boolean")

        return true

    end

    validate_mass_props_unc(mpu) = begin

        mpu.sigma_mass isa Real || error("mass uncertainty must be a real number")
        mpu.sigma_mass >= 0.0 || error("mass uncertainty must be non-negative")

        any(ismissing.(mpu.sigma_center_mass)) && error("center mass uncertainty component is missing")
        any(isnothing.(mpu.sigma_center_mass)) && error("center mass uncertainty component is nothing")
        mpu.sigma_center_mass isa AbstractVector{<:Real} || error("center mass uncertainty must be a vector of real numbers")
        length(mpu.sigma_center_mass) == 3 || error("center mass uncertainty must have three components")
        any(mpu.sigma_center_mass .< 0.0) && error("center mass uncertainty must be non-negative")

        any(ismissing.(mpu.sigma_inertia)) && error("inertia uncertainty component is missing")
        any(isnothing.(mpu.sigma_inertia)) && error("inertia uncertainty component is nothing")
        mpu.sigma_inertia isa AbstractMatrix{<:Real} || error("inertia uncertainty must be a matrix of real numbers")
        size(mpu.sigma_inertia) == (3, 3) || error("inertia uncertainty must be a 3x3 matrix")
        any(mpu.sigma_inertia .< 0.0) && error("inertia uncertainty must be non-negative")

        return true

    end

    validate_mass_props_and_unc(mpu) = validate_mass_props(mpu) && validate_mass_props_unc(mpu)

    validate_mass_props_table(tree, df) = validate_ds(tree, df, df_get_ids, get_mass_props, validate_mass_props)

    validate_mass_props_and_unc_table(tree, df) = validate_ds(tree, df, df_get_ids, get_mass_props_and_unc, validate_mass_props_and_unc)

    rollup_mass_props(tree, df, validate_df = validate_mass_props_table) = rollup(tree, df, update_mass_props, validate_df)

    rollup_mass_props_unc(tree, df, validate_df = validate_mass_props_and_unc_table) = rollup(tree, df, update_mass_props_unc, validate_df)

    rollup_mass_props_and_unc(tree, df, validate_df = validate_mass_props_and_unc_table) = rollup(tree, df, update_mass_props_and_unc, validate_df)

end