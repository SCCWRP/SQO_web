# format species inputs ---------------------------------------------------

#' format species inputs
#'
#' @param inps reactive inputs
#' @param biota table
formsppinp <- function(inps, biota){
  
  # format names and lips as tibble
  frminps <- reactiveValuesToList(inps) %>% 
    enframe('Biota', 'value') %>% 
    filter(!grepl('selectized', Biota) & grepl('indic[0-9]lip', Biota)) %>% 
    mutate(
      var = case_when(
        grepl('lip', Biota) ~ 'lipid', 
        T ~ 'spp'
      ), 
      Biota = gsub('lip', '', Biota)
    ) %>% 
    spread(var, value) %>% 
    unnest
  
  # add frminps user input to biota
  out <- biota %>% 
    left_join(frminps, by = 'Biota') %>% 
    mutate(lipid.x = ifelse(is.na(lipid.x), lipid.y, lipid.x)) %>% 
    rename(lipid = lipid.x) %>% 
    select(-lipid.y)
  
  return(out)

}

# format site/env constants inputs ---------------------------------------------------

#' format site/env constants inputs
#'
#' @param inps reactive inputs
#' @param constants table
formcnsinp <- function(inps, constants){

  # format names and lips as tibble
  frminps <- reactiveValuesToList(inps) %>% 
    enframe('Constant', 'Value') %>% 
    filter(Constant %in% c('SA', 'SL', 'Cox', 'T', 'salinity', 'ocsed', 'vss', 'xpoc', 'xdoc')) %>% 
    unnest
  
  # add frminps user input to constants
  out <- constants %>% 
    left_join(frminps, by = 'Constant') %>% 
    mutate(Value.x = ifelse(is.na(Value.x), Value.y, Value.x)) %>% 
    rename(Value = Value.x) %>% 
    select(-Value.y)
  
  return(out)
  
}

# format contaminant inputs ---------------------------------------------------

#' format contaminant inputs
#'
#' @param inps reactive inputs
#' @param contam table
formcntinp <- function(inps, contam){

  # format names and contams as tibble
  frminps <- reactiveValuesToList(inps) %>% 
    enframe('Chem', 'Value') %>% 
    filter(grepl('^alpha|^gamma|^Oxy|^Dieldrin|^op\\-|^pp\\-|^PCB', Chem)) %>% 
    unnest %>% 
    separate(Chem, c('Chem', 'loc'), sep = '_') %>% 
    spread(loc, Value) %>% 
    rename(
      `cs_ng.g` = sed,  # sediment concentration
      `cd_ng.g` = dis,   # dissolved surface water concentration
      `cp_ng.g` = por    # porewater concentration
    )
  
  # add frminps user input to contaminants
  out <- contam %>% 
    left_join(frminps, by = 'Chem') %>% 
    mutate(
      cs_ng.g.x = ifelse(is.na(cs_ng.g.x), cs_ng.g.y, cs_ng.g.x),
      cd_ng.g.x = ifelse(is.na(cd_ng.g.x), cd_ng.g.y, cd_ng.g.x),
      cp_ng.g.x = ifelse(is.na(cp_ng.g.x), cp_ng.g.y, cp_ng.g.x)
      ) %>% 
    rename(
      cs_ng.g = cs_ng.g.x,
      cd_ng.g = cd_ng.g.x,
      cp_ng.g = cp_ng.g.x
      ) %>% 
    select(-cs_ng.g.y, -cd_ng.g.y, -cp_ng.g.y)
  
  return(out)
  
}

# format mcs inputs -------------------------------------------------------

#' Format mcs inputs
#'
#' @param inps reactive inputs
#' @param contam table
formmcsinp <- function(inps, mcsparms){

  # format names and input mean/sd as tibble
  frminps <- reactiveValuesToList(inps) %>% 
    enframe('MCSvar', 'Value') %>% 
    filter(grepl('X$|SD$', MCSvar)) %>% 
    unnest
  
  # default parameters from table
  mcsparms <- mcsparms %>% 
    dplyr::select(name, value) %>% 
    rename(
      MCSvar = name, 
      Value = value
    )
  
  # combine
  out <- mcsparms %>% 
    bind_rows(frminps)
  
  return(out)
  
}
    
# format seafood proportions in diet from user inputs ---------------------

#' Format seafood proportions in diet from user inputs
#'
#' @param inps reactive input list
#'
formpropseaf <- function(inps){
  
  out <- reactiveValuesToList(inps) %>% 
    enframe('Biota', 'value') %>% 
    filter(!grepl('selectized', Biota) & grepl('indic[0-9]seaf$', Biota)) %>% 
    unnest %>% 
    mutate(
      Biota = gsub('seaf$', '', Biota),
      value = ifelse(is.na(value), 0, value)
      ) %>% 
    arrange(Biota) %>% 
    pull(value)
  
  return(out)

}

# weighted average observed tissue concentration (ng/g), from empirical data -------------------

#' Calculate weighted average observed tissue concentration (ng/g), from empirical data
#'
#' @param mcsparms input mcsparms data frame, observed average concentrations extracted
#' @param inps shiny reactives, extracts proportion seafood
#'
wgt_avg_fun <- function(mcsparms, inps){
  
  # propseaf for guild species
  propseaf <- formpropseaf(inps) 
  
  # observed contaminants from user input, mean only
  contobs <- mcsparms %>% 
    filter(grepl('^indic.*X$', MCSvar)) %>% 
    mutate(
      contam = gsub('^indic[0-9](.*)X$', '\\1', MCSvar), 
      MCSvar = gsub('(^indic[0-9]).*$', '\\1', MCSvar), 
      Value = case_when(
        is.na(Value) ~ 0, 
        T ~ Value
      )
    ) %>% 
    arrange(contam, MCSvar)
  
  # weighted average observed tissue conc
  wgt_avg <- contobs %>% 
    group_by(contam) %>%
    summarise(
      wgt_obs = Value %*% propseaf
    )
  
  return(wgt_avg)
  
}

# generate log normal variables from mean and sd inputs -------------------

#' Generate log normal variables from mean and sd inputs
#'
#' @param nsim number of simulations
#' @param X mean value from user input for single contaminant
#' @param SD SD value from user input for single contaminant
#' 
#' @details http://yasai.rutgers.edu/yasai-guide-27.html
genlognorm_fun <- function(nsim, X, SD){
  
  # # genlognormal, see link for doc
  # sims <- suppressWarnings(rlnorm(nsim, meanlog = X, sdlog = SD)) %>% 
  #   log(.) %>% 
  #   pmax(0, .)
  
  # genlognormal, see link for doc
  mu <- log(X) - 0.5 * log(1 + SD ^ 2 / X ^ 2)
  sigma <- (log(SD ^ 2 / X ^ 2 + 1)) ^ 0.5
  sims <- exp(rnorm(nsim, mu, sigma))
  
  simi <- seq(1:nsim)
  out <- data.frame(i = simi, sims = sims)
  return(out)
  
}

# mcs function for modelled tissue concentration --------------------------

#' mcs function for modelled tissue concentration
#'
#' @param nsim number of simulations
#' @param meanse mean and se values from user input for each guild species and contaminant class
#' @param propseaf proportion of seafood diet, output from formpropseaf
#'
modtiscon_mcs_fun <- function(nsim, meanse, propseaf){
  
  # simulated tissue concentrations across guild species, all sims
  sims <- meanse %>%  
    group_by(MCSvar, contam) %>% 
    mutate(
      ests = purrr::pmap(list(nsim, X, SD), genlognorm_fun)
    ) %>% 
    dplyr::select(-SD, -X) %>% 
    unnest
  
  # weighted tissue concentrations across guilds for each contam, all sims
  out <- sims %>% 
    dplyr::group_by(i) %>% 
    nest() %>% 
    mutate(
      wgtave = purrr::map(data, function(x){
        
        sumprod <- x %>% 
          arrange(contam, MCSvar) %>% 
          mutate(
            sims = case_when(
              is.na(sims) ~ 0, 
              T ~ sims
            )
          ) %>% 
          group_by(contam) %>% 
          summarise(
            wgtave = sims %*% propseaf
          )
        
        return(sumprod)
        
      })
    ) %>% 
    dplyr::select(-data)
  
  return(out)

}

# mcs function for sediment concentration ---------------------------------

#' mcs function for sediment concentration
#'
#' @param nsim number of simulations
#' @param sedmeanse sediment mean and se values from user input for each contaminant class
#' @param propseaf proportion of seafood diet, output from formpropseaf
#' @param SUF site use factor
#' @param CVBAF bioaccumulation factor sd/mean from mcsparms
#' @param indic_sum indicator guild concentrations sums across contaminants
#'
modsedcon_mcs_fun <- function(nsim, sedmeanse, propseaf, SUF, CVBAF, indic_sum){

  # simulated sediment concentrations, all contams
  sedsims <- sedmeanse %>%  
    group_by(contam) %>% 
    mutate(
      ests = purrr::pmap(list(nsim, X, SD), genlognorm_fun) 
    ) %>% 
    dplyr::select(-SD, -X) %>% 
    unnest %>% 
    dplyr::ungroup()
   
  # bioaccumulation sims
  biosims <- indic_sum %>% 
    tidyr::gather('contam', 'val', -species) %>% 
    filter(grepl('\\_calc$', contam)) %>% 
    mutate(
      contam = gsub('\\_calc$', '', contam)
    ) %>% 
    dplyr::group_by(species, contam) %>% 
    mutate(
      ests = purrr::pmap(list(nsim, val, CVBAF * val), genlognorm_fun)
    ) %>% 
    dplyr::select(-val) %>% 
    dplyr::ungroup()
    
  # combine sediment sims with biosims and SUF
  estcncsims <- biosims %>% 
    unnest %>% 
    full_join(sedsims, ., by = c('contam', 'i')) %>% 
    full_join(SUF, by = c('i', 'species')) %>% 
    mutate(
      estcnc = sims.x * suf * sims.y
    ) %>% 
    dplyr::select(-MCSvar, -sims.x, -sims.y, -suf)
  
  # weighted sediment concentrations across guilds for each contam, all sims
  out <- estcncsims %>% 
    dplyr::group_by(i) %>% 
    nest() %>% 
    mutate(
      wgtave = purrr::map(data, function(x){
        
        sumprod <- x %>% 
          arrange(contam, species) %>% 
          mutate(
            estcnc = case_when(
              is.na(estcnc) ~ 0, 
              T ~ estcnc
            )
          ) %>% 
          group_by(contam) %>% 
          summarise(
            wgtave = estcnc %*% propseaf
          )
        
        return(sumprod)
        
      })
    ) %>% 
    dplyr::select(-data)

  return(out)
  
}


# site use factor sims ----------------------------------------------------

suf_mcs_fun <- function(nsim, constants, mcsparms){
  
  # site area and length
  SA <- constants %>% 
    filter(Constant %in% 'SA') %>% 
    pull(Value)
  SL <- constants %>% 
    filter(Constant %in% 'SL') %>%
    pull(Value)
  
  # home range mean and sd for guild species
  hrvals <- mcsparms %>% 
    filter(grepl('^HR[0-9]', MCSvar)) %>% 
    rename(species = MCSvar) %>% 
    mutate(
      var = case_when(
        grepl('X$', species) ~ 'X', 
        grepl('SD$', species) ~ 'SD'
      ),
      species = gsub('^HR', 'indic', species),
      species = gsub('X$|SD$', '', species)
    ) %>% 
    spread(var, Value)
  
  # home range sims
  sufsims <- hrvals %>% 
    group_by(species) %>% 
    mutate(
      suf = purrr::map(list(species), function(...){

        # indic1, indic8, indic9
        if(grepl('1$|8$|9$', species))
          out <- genlognorm_fun(nsim, X, SD) %>% 
            mutate(
              sims = SL / sims,
              sims = ifelse(is.infinite(sims), 0, sims)
              )
        
        # indic2, indic3, indic4, indic5, indic7
        if(grepl('2$|3$|4$|5$|7$', species))
          out <- genlognorm_fun(nsim, X, SD) %>% 
            mutate(
              sims = SA / sims,
              sims = ifelse(is.infinite(sims), 0, sims)
            )
        
        # indic6
        if(grepl('6$', species)){
          out <- (SL * 1000) / pgamma(runif(nsim, 0, 1), shape = X, scale = SD) 
          simi <- seq(1:nsim)
          out <- tibble(i = simi, sims = out)
        }
        
        return(out)
        
      })
    ) %>% 
    dplyr::select(-SD, -X) %>% 
    unnest %>% 
    mutate(
      sims = pmin(1, sims)
    ) %>% 
    rename(suf = sims)
  
  return(sufsims)
  
}

# MCS function --------------------------------------------------------------
 
#' MCS function
#'
#' @param inps reactive inputs
#' @param nsim number of MC sims
#' @param indic_sum output from indic_sum_fun
#' @param mcsparms MCS parameter inputs
#' @param constants constants inputs
#'
mcs_fun <- function(inps, nsim, indic_sum, mcsparms, constants){
  
  ##
  # inputs 
  
  # CVBAF
  CVBAF <- mcsparms %>% 
    filter(MCSvar == 'CVBAF') %>% 
    pull
  
  # seafood diet proportion
  propseaf <- formpropseaf(inps) 

  # mean and se values from observed contaminants, from user inputs
  meanse <- mcsparms %>% 
    filter(grepl('^indic', MCSvar)) %>% 
    mutate(
      var = case_when(
        grepl('X$', MCSvar) ~ 'X', 
        grepl('SD$', MCSvar) ~ 'SD'
      ),
      contam = gsub('^indic[0-9]|X$|SD$', '', MCSvar), 
      MCSvar = gsub('(^indic[0-9]).*$', '\\1', MCSvar)
    ) %>% 
    spread(var, Value)
  
  # mean and se values for sediment contaminants, from user inputs
  sedmeanse <- mcsparms %>% 
    filter(grepl('^sed', MCSvar)) %>% 
    mutate(
      var = case_when(
        grepl('X$', MCSvar) ~ 'X', 
        grepl('SD$', MCSvar) ~ 'SD'
      ),
      contam = gsub('^sed|X$|SD$', '', MCSvar), 
      MCSvar = gsub('(^sed).*$', '\\1', MCSvar)
    ) %>% 
    spread(var, Value)
  
  ##
  # modeled tissue concentration for consumption risk, mcs
  # returns weighted concentrations across all sims
  modtiscon <- modtiscon_mcs_fun(nsim, meanse, propseaf)
  
  ##
  # site use function sims
  SUF <- suf_mcs_fun(nsim, constants, mcsparms)
  
  ##
  # modeled sediment contribution to tissue concentration, mcs
  # returns weighted concentrations across all sims
  modsedcon <- modsedcon_mcs_fun(nsim, sedmeanse, propseaf, SUF, CVBAF, indic_sum)
  
  ## 
  # combine modeled tissue and sediment concentrations to get site linkages
  out <- modtiscon %>% 
    full_join(modsedcon, by = 'i') %>% 
    tidyr::unnest() %>% 
    mutate(sitsedlnk = wgtave1 / wgtave) %>% 
    dplyr::select(-wgtave, -contam1, -wgtave1)
    
  return(out)

}

# MCS summary function ----------------------------------------------------

#' Summarize MCS results, compare with observed
#'
#' @param mcsres MC results, output from \code{mcs_fun}
#'
#' @return
#' @export
#'
#' @examples
mcs_sum_fun <- function(mcsres){
  
  # get percentiles
  persitsed <- mcsres %>% 
    group_by(contam) %>% 
    nest %>% 
    mutate(
      percnt = purrr::map(data, function(x){
        
        prc <- quantile(x$sitsedlnk, c(0, .01, .05, .1, 0.25, .5, .75, 0.9, .95, .99, 1)) %>% 
          enframe 
        
        return(prc)
        
      })
    ) %>% 
    dplyr::select(-data) %>% 
    unnest %>% 
    mutate(name = factor(name, levels = c('0%', '1%', '5%', '10%', '25%', '50%', '75%', '90%', '95%', '99%', '100%'))) %>% 
    rename(percentile = name)
  
  return(persitsed)
  
}

# SQO assessment summary --------------------------------------------------

#' SQO assessment summary
#'
#' @param wgtavg weighted average observed tissue concentrations from input, by contaminant category, output from \code{wgt_avg_fun}
#' @param MCSsum summary results by percentiles of MC analyses, output from \code{mcs_sum_fun}
#' @param tischmthr lookup table for tissue chemistry thresholds
#' @param constants constants from user inputs and lookup table, only SCT is used (sediment linkage threshold)
#' @param finalsiteassess final site assessment lookup table
#'
#' @return
#' @export
#'
#' @examples
sqo_sum_fun <- function(wgtavg, MCSsum, tischmthr, constants, finalsiteassess){

  # category scores and labels, final labels
  levs <- c('1', '2', '3', '4', '5')
  labs <- c('Very Low', 'Low', 'Moderate', 'High', 'Very High')
  flabs <- c('Unimpacted', 'Likely Unimpacted', 'Possibly Impacted', 'Likely Impacted', 'Clearly Impacted')
  
  # sediment linkage threshold
  SCT <- constants %>% 
    filter(Constant %in% 'SCT') %>% 
    pull(Value)
  
  # thresholds in nested format for join
  tischmthr <- tischmthr %>% 
    group_by(contam) %>% 
    nest(.key = 'thr')
  
  # quartiles from MCSsum
  mcsres <- MCSsum %>% 
    filter(grepl('25|50|75', percentile))

  # combined data to get category outcomes
  cmb <- wgtavg %>% 
    full_join(mcsres, by = 'contam') %>% 
    spread(percentile, value) %>% 
    full_join(tischmthr, by = 'contam') 
  
  # get category outcomes
  # chmscr/chmlab - chemical exposure category score
  # lnkscr/lnklab - site linkage category score
  # sitscr/sitlab - final site assessment category score
  sums <- cmb %>% 
    mutate(
      wgt_est = wgt_obs * `50%`,
      chmscr = purrr::pmap(list(wgt_obs, thr), function(wgt_obs, thr){
   
        val <- thr %>% pull(val)
        scr <- 1 + findInterval(wgt_obs, val)
        
        return(scr)
        
      }),
      chmlab = factor(as.character(chmscr), levels = levs, labels = labs), 
      chmlab = as.character(chmlab)
    ) %>% 
    rowwise() %>% 
    mutate(
      lnkscr = 4 - findInterval(SCT, c(`25%`, `50%`, `75%`)), 
      lnklab = factor(as.character(lnkscr), levels = levs, labels = labs), 
      lnklab = as.character(lnklab)
    ) %>% 
    unite('cmbscr', chmscr, lnkscr, sep = ', ', remove = F) %>% 
    mutate(
      sitscr = factor(cmbscr, levels = finalsiteassess[[1]], labels = finalsiteassess[[2]]), 
      sitscr = as.numeric(as.character(sitscr)), 
      sitlab = factor(sitscr, levels = levs, labels = flabs), 
      sitlab = as.character(sitlab)
    )
  
  # final formatting (no calcs)
  out <- sums %>% 
    select(-thr) %>% 
    unnest %>% 
    select(contam, wgt_obs, `25%`, `50%`, `75%`, wgt_est, chmscr, chmlab, lnkscr, lnklab, sitscr, sitlab) %>% 
    rename(
      Compound = contam,
      `Weighted observed tissue conc. (ng/g)` = wgt_obs,
      `Weighted estimated tissue conc. (ng/g)` = wgt_est,
      `Chemical exposure score` = chmscr, 
      `Chemical exposure category` = chmlab, 
      `Site linkage score` = lnkscr, 
      `Site linkage category` = lnklab, 
      `Site assessment score` = sitscr, 
      `Site assessment category` = sitlab
    )
  
  return(out)
  
}