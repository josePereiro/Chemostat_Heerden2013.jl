# Collect
DAT = ChU.DictTree()
let 
    
    # CACHE
    DATfile = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, "2.1_DAT.jls")
    if isfile(DATfile) 
        global DAT = deserialize(DATfile) 
        @info("DAT CACHE LOADED")
        return
    end

    objider = iJR.BIOMASS_IDER
    DAT[:CONC_IDERS] = CONC_IDERS
    DAT[:FLX_IDERS] = FLX_IDERS
    DAT[:EXPS] = []

    # Find exps
    for exp in Hd.EXPS
        ok = false
        for method in ALL_METHODS
            ok = haskey(INDEX, method, :DFILE, exp) &&
                INDEX[method, :DFILE, exp] != :unfeasible
            !ok && break
        end
        !ok && continue
        push!(DAT[:EXPS], exp)
    end
    max_model = iJR.load_model("max_model"; uncompress = false)

    for exp in DAT[:EXPS], method in ALL_METHODS
            
        datfile = INDEX[method, :DFILE, exp]
        dat = deserialize(datfile)
        
        model = dat[:model]
        objidx = ChU.rxnindex(model, objider)
        exp_beta = dat[:exp_beta] 
        epout = dat[:epouts][exp_beta]
        exp_xi = Hd.val(:xi, exp)
        fva_model = iJR.load_model("fva_models", exp; uncompress = false)

        println()
        @info("Doing", 
            exp, method, exp_beta, 
            length(dat[:epouts]), 
            epout.iter, datfile
        ); 

        # Biomass
        ep_biom = ChU.av(model, epout, objidx)
        ep_std = sqrt(ChU.va(model, epout, objidx))
        Hd_biom = Hd.val("D", exp)
        max_lb, max_ub = ChU.bounds(max_model, objidx)
        fva_lb, fva_ub = ChU.bounds(fva_model, objidx)
        lb = max(max_lb, fva_lb)
        ub = min(max_ub, fva_ub)
        
        # store
        DAT[method, :ep   , :flx, objider, exp] = ep_biom
        DAT[method, :eperr, :flx, objider, exp] = ep_std
        DAT[method, :Hd   , :flx, objider, exp] = Hd_biom
        DAT[:Hd   , :flx, objider, exp] = Hd_biom
        DAT[method, :bounds, :flx, objider, exp] = (lb, ub)
        
        # mets
        for Hd_met in FLX_IDERS

                model_met = Hd_mets_map[Hd_met]
                model_exch = Hd_rxns_map[Hd_met]
                model_exchi = ChU.rxnindex(model, model_exch)

                # fuxes
                proj = ChLP.projection2D(model, objider, model_exchi; l = 50)
                ep_av = ChU.av(model, epout, model_exchi)
                ep_std = sqrt(ChU.va(model, epout, model_exchi))
                Hd_flx = Hd.val("u$Hd_met", exp)
                max_lb, max_ub = ChU.bounds(max_model, model_exchi)
                fva_lb, fva_ub = ChU.bounds(fva_model, model_exchi)
                lb = max(max_lb, fva_lb)
                ub = min(max_ub, fva_ub)
                
                DAT[method, :Hd, :flx, Hd_met, exp] = Hd_flx
                DAT[:Hd, :flx, Hd_met, exp] = Hd_flx
                DAT[method, :ep, :proj, Hd_met, exp] = proj
                DAT[method, :ep, :flx, Hd_met, exp] = ep_av
                DAT[method, :eperr, :flx, Hd_met, exp] = ep_std
                DAT[method, :bounds, :flx, Hd_met, exp] = (lb, ub)

        end

        for Hd_met in CONC_IDERS

            ep_av = DAT[method, :ep, :flx, Hd_met, exp]
            ep_std = DAT[method, :eperr, :flx, Hd_met, exp] 

            # conc (s = c + u*xi)
            c = Hd.val("c$Hd_met", exp, 0.0)
            ep_conc = max(c + ep_av * exp_xi, 0.0)
            Hd_conc = Hd.val("s$Hd_met", exp)
            
            DAT[method, :Hd, :conc, Hd_met, exp] = Hd_conc
            DAT[:Hd, :conc, Hd_met, exp] = Hd_conc
            DAT[method, :ep, :conc, Hd_met, exp] = ep_conc
            DAT[method, :eperr, :conc, Hd_met, exp] = ep_std * exp_xi
        end

    end # for exp in EXPS

    DAT[:EXPS] |> unique! |> sort!
    
    # saving
    serialize(DATfile, DAT)
end

# Inter project comunication
let
    CORR_DAT = isfile(iJR.CORR_DAT_FILE) ? ChU.load_data(iJR.CORR_DAT_FILE) : Dict()
    CORR_DAT[:MAXENT_EP] = DAT
    ChU.save_data(iJR.CORR_DAT_FILE, CORR_DAT)
end