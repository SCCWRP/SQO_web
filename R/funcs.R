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
    filter(grepl('X$|SD$|indic[0-9]seaf$', MCSvar)) %>% 
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
