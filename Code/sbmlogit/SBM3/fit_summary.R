fit_summary <- function(sbm_fit)
{
    labs = get_labels(fit = sbm_fit)
    kS = seq_len(sbm_fit$ngroups)
    gams = col_summarize(as.matrix(sbm_fit$gamma))
    gams = c("Ref.",gams)
}