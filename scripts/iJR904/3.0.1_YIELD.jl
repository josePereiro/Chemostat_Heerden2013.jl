let
    objider = iJR.BIOMASS_IDER
    glcider = "EX_glc_LPAREN_e_RPAREN__REV"
    iterator = Hd.val(:D) |> enumerate |> collect 
    for (exp, D) in iterator
        try 
            # prepare model
            model = base_model(exp)
            M, N = size(model)
            model = ChU.fix_dims(model)
            ChU.invert_bkwds!(model)
            exglc_idx = ChU.rxnindex(model, glcider)
            objidx = ChU.rxnindex(model, objider)
            exp_growth = Hd.val(:D, exp)
            dgrowth = 0.01
            
            ChU.ub!(model, objider, exp_growth * (1 + dgrowth))

            # fba
            fbaout = ChLP.fba(model, objider)
            fba_growth = ChU.av(model, fbaout, objider)

            # yield max
            model.c .= 0.0; model.c[objidx] = 1.0
            d = zeros(N); d[exglc_idx] = 1.0
            status, yflxs, yield = ChLP.yLP(model, d)
            status != JuMP.MOI.OPTIMAL && @warn(status)
            ymax_growth = yflxs[objidx]
            LPDAT[YIELD, :status, exp] = status
            LPDAT[YIELD, :yield, exp] = yield
            LPDAT[YIELD, :yout, exp] = ChU.FBAout(yflxs, yflxs[objidx], objider, nothing)
            LPDAT[YIELD, :d, exp] = d
            LPDAT[YIELD, :model, exp] = model

            @info("Yield Maximization", 
                exp, status, yield,
                fba_growth, ymax_growth, exp_growth
            ); println()

        catch err
            @warn("Error", err, exp); println()
        end
    end
end
