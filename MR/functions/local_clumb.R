local_clumb <- function (dat, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, 
    clump_p2 = 1, pop = "EUR") 
{
    pval_column <- "pval.exposure"
    if (!is.data.frame(dat)) {
        stop("Expecting data frame returned from format_data")
    }
    if ("pval.exposure" %in% names(dat) & "pval.outcome" %in% 
        names(dat)) {
        message("pval.exposure and pval.outcome columns present. Using pval.exposure for clumping.")
    }
    else if (!"pval.exposure" %in% names(dat) & "pval.outcome" %in% 
        names(dat)) {
        message("pval.exposure column not present, using pval.outcome column for clumping.")
        pval_column <- "pval.outcome"
    }
    else if (!"pval.exposure" %in% names(dat)) {
        message("pval.exposure not present, setting clumping p-value to 0.99 for all variants")
        dat$pval.exposure <- 0.99
    }
    else {
        pval_column <- "pval.exposure"
    }
    if (!"id.exposure" %in% names(dat)) {
        dat$id.exposure <- random_string(1)
    }
    d <- data.frame(rsid = dat$SNP, pval = dat[[pval_column]], 
        id = dat$id.exposure)
    out <- ieugwasr::ld_clump(d, clump_kb = clump_kb, clump_r2 = clump_r2, 
        clump_p = clump_p1,
    plink_bin =  "/mnt/data/lijincheng/miniconda3/envs/mamba/envs/jupyter/bin/plink",  #plink 
    bfile = paste0("/mnt/data/lijincheng/database/twosampleMR_1000gLD_reference/",pop)
)
        
    keep <- paste(dat$SNP, dat$id.exposure) %in% paste(out$rsid, 
        out$id)
    return(dat[keep, ])
}
