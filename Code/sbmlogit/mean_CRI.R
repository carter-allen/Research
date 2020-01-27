mean_CRI <- function(y, dig = 2)
{
    m = round(mean(y, na.rm = TRUE),dig)
    cri = round(quantile(y, probs = c(0.025,0.95)), dig)
    ret = paste0(m," (",cri[1],", ",cri[2],")")
    return(ret)
}