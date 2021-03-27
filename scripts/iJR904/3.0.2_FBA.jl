let
    objider = iJR.BIOMASS_IDER
    costider = iJR.COST_IDER
    exglcider = iJR.EX_GLC_IDER
    max_sense = 1.0
    min_sense = 1.0

    iterator = Hd.val(:D) |> enumerate |> collect 
    for (exp, D) in iterator

        @info("Doing ", exp); println()

        # FBA_OPEN
        let
            model = base_model(exp)
            fbaout = ChLP.fba(model, objider, costider)
            
            LPDAT[FBA_OPEN, :model, exp] = model
            LPDAT[FBA_OPEN, :fbaout, exp] = fbaout
        end

        # FBA_MIN_COST
        let
            model = base_model(exp)
            exp_growth = Hd.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            fbaout = ChLP.fba(model, objider, costider)

            LPDAT[FBA_MIN_COST, :model, exp] = model
            LPDAT[FBA_MIN_COST, :fbaout, exp] = fbaout
        end

        # FBA_MAX_VG_YIELD
        let
            model = base_model(exp)
            exp_growth = Hd.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            fbaout1 = ChLP.fba(model, exglcider; sense = max_sense)
            exglc = ChU.av(model, fbaout1, exglcider)
            ChU.bounds!(model, exglcider, exglc, exglc)
            fbaout = ChLP.fba(model, costider; sense = min_sense)

            LPDAT[FBA_MAX_VG_YIELD, :model, exp] = model
            LPDAT[FBA_MAX_VG_YIELD, :fbaout, exp] = fbaout
        end
        
    end
end
