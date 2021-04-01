# Collect
DAT = ChU.DictTree()
let 
    
    # CACHE
    DATfile = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, "2.1_DAT.jls")
    # if isfile(DATfile) 
    #     global DAT = deserialize(DATfile) 
    #     @info("DAT CACHE LOADED")
    #     return
    # end

    objider = iJR.BIOMASS_IDER
    DAT[:FLX_IDERS] = FLX_IDERS
    DAT[:EXPS] = []

    # Find exps
    for exp in Hd.EXPS
        ok = false
        for method in ALL_METHODS
            ok = haskey(ME_INDEX, method, :DFILE, exp) &&
                ME_INDEX[method, :DFILE, exp] != :unfeasible
            !ok && break
        end
        !ok && continue
        push!(DAT[:EXPS], exp)
    end
    max_model = iJR.load_model("max_model"; uncompress = false)

    for exp in DAT[:EXPS], method in ALL_METHODS
            
        datfile = ME_INDEX[method, :DFILE, exp]
        dat = deserialize(datfile)
        
        # ep data
        ep_model = dat[:model]
        exp_beta = dat[:exp_beta] 
        epout = dat[:epouts][exp_beta]
        exp_xi = Hd.val(:xi, exp)
        fva_model = iJR.load_model("fva_models", exp; uncompress = false)

        # fba
        fba_model = LPDAT[method, :model, exp]
        fbaout = LPDAT[method, :fbaout, exp]
        
        
        println()
        @info("Doing", 
            exp, method, exp_beta, 
            length(dat[:epouts]), 
            epout.iter, datfile
        ); 

        # ep
        ep_biom = ChU.av(ep_model, epout, objider)
        ep_std = sqrt(ChU.va(ep_model, epout, objider))
        Hd_biom = Hd.val("D", exp)

        # global
        max_lb, max_ub = ChU.bounds(max_model, objider)
        fva_lb, fva_ub = ChU.bounds(fva_model, objider)
        lb = max(max_lb, fva_lb)
        ub = min(max_ub, fva_ub)

        # # TODO collect fba aswell
        # # fba
        # fba_biom = ChU.av(fba_model, fbaout, objider)
        # DAT[method, :Hd, :flx, objider, exp] = Hd_flx
        # DAT[method, :lp, :flx, objider, exp] = fba_flx

        
        # store
        DAT[method, :ep     , :flx, objider, exp] = ep_biom
        DAT[method, :eperr  , :flx, objider, exp] = ep_std
        DAT[method, :Hd     , :flx, objider, exp] = Hd_biom
        DAT[method, :bounds , :flx, objider, exp] = (lb, ub)
        DAT[        :Hd     , :flx, objider, exp] = Hd_biom
        
        # mets
        for Hd_met in FLX_IDERS

                model_met = Hd_mets_map[Hd_met]
                model_exch = Hd_rxns_map[Hd_met]
                model_exchi = ChU.rxnindex(ep_model, model_exch)

                # fuxes
                proj = ChLP.projection2D(model, objider, model_exchi; l = 50)
                ep_av = ChU.av(ep_model, epout, model_exchi)
                ep_std = sqrt(ChU.va(ep_model, epout, model_exchi))
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