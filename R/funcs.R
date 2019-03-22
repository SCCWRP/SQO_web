# bioaccumulation function ------------------------------------------------

#FoodWeb_SQO
# Purpose: To implement the Gobas/Arnot Food Web Bioaccumulation Model
#          of San Francisco Bay for various species and contaminants.
# 6.17.2010
# Michelle Lent
# edited 6.23.10 & 6.30.10 by BG and ML
FoodWeb_SQO <- function(NumSim, csed, cwater, KowTS, Kow, EdA, EdB, xdoc, 
                              ddoc, xpoc, dpoc, alphapoc, alphadoc, ocsed, ds, taxa, A, B, T, lipid, nloc, 
                              nlom, wc, beta, betap, mo, mp, kM, Wb, Cox, vss, scav, preyprop, cbiota, 
                              vld, vcd, vnd, vwd, GR, assimEff_1, assimEff_2, assimEff_3){ #updated 6.30.2010
 
  ## initialize the output variables
  temp <- data.frame(cbiota = NA, cprey = NA, k1=NA, k2=NA, GR=NA, Gv=NA, Gd=NA, 
                     Gf=NA, vlg=NA, vcg=NA, vng=NA, vwg=NA, kgb=NA, ke=NA, kd=NA, Ew=NA, Ed=NA, phi=NA, 
                     cpw=NA, assimEff_1=NA, assimEff_2=NA, assimEff_3=NA);
  Output <- temp[rep(1:nrow(temp), times=NumSim), ]
  
  lden <- 0.9; # lipid density - added 6.30.2010

  ##
  ## Calculate contaminant dependant parameters
  ##
  #Ew <- 1/(1.85+1.55/KowTS); # gill chemical uptake efficiency
  Ew <- 1/(1.85+155/KowTS); # gill chemical uptake efficiency #fixed 6.30.2010
  Ed <- 1/(EdA*Kow + EdB); # dietary chemical transfer efficiency (gut uptake efficiency)
  
  # freely dissolved contaminant fraction in overlying water column
  #phi <- 1/(1 + xpoc*dpoc*alphapoc*KowTS + xdoc*ddoc*alphapoc*alphadoc*KowTS);
  phi <- 1/(1 + xpoc*dpoc*alphapoc*KowTS + xdoc*ddoc*alphadoc*KowTS); #fixed 6.30.2010
  
  # [contaminant] in porewater
  #cpw = csed*ds/(ocsed*0.35*KowTS);
  #corrected to remove density of organic carbon in suspended sediment as per email from Jon Arnot 2/1/06
  #cpw <- csed/(ocsed*0.35*KowTS);
  cpw <- csed*ds/(ocsed*alphapoc*KowTS); #fixed 6.30.2010
  
  ## assign contaminant specific outputs
  Output$Ew <- Ew;
  Output$Ed <- Ed;
  Output$phi <- phi;
  Output$cpw <- cpw;
  
  ##
  ## Now proceed with the biota-specific calculations.
  ##
  ## Calculations are taxa dependant. The main IF loop here represents that dependancy.

  if (taxa==0) {Output$cbiota <- csed}; # This is a dummy taxa used by sediment
  
  if (taxa>=1 & taxa<2)
  {
    # plants (i.e., phytoplankton) only have selected parameters
    
    # aqueous uptake rate constant different for phytoplankton
    k1 <- ((A + (B/KowTS))^-1);
    
    # biota-water partition coefficient for phytoplankton; i.e., the bioconcentration factor
    #kbw <- (lipid*KowTS + nlom*betap*KowTS + wc); 		
    kbw <- (lipid*KowTS/lden + nloc*betap*KowTS + nlom*beta*KowTS + wc); #updated 6.30.2010
    
    # elimination rate constant for phytoplankton ******AM changed from kbwp********
    k2 <- k1/kbw; 
    kG <- GR
    # calculate concentration in organism
    # based on uptake/loss pseudoequilibrium
    # main model equation (only includes plant parameters)
    Output$cbiota <- k1*(mo*phi*cwater + mp*cpw)/(k2 + kG + kM);
    # assign output variables
    Output$k1 <- k1;
    Output$k2 <- k2;
    #		Output$kG <- kG;
  }
  
  if (taxa>=2)  # ML- CHANGED THIS FROM taxa>1 b/c code was trying calculated fecal egestion rate for taxa 1.1 which is sea lettuce (ulva)!!! However, should check that this doesn't cause other problems. ###
  {
    if (taxa==2|taxa==3|taxa==4|taxa==5) 
    {
      Gv <- (1400*Wb^0.65)/Cox;  # Calculate Gv - gill ventilation rate
      k1 <- Ew*Gv/Wb;  # k1 = aqueous uptake rate constant (L/kg*d)            
      #kbw <- (lipid*KowTS + nlom*beta*KowTS + wc);  # kbw = biota-water partition coefficient
      kbw <- (lipid*KowTS/lden + nloc*betap*KowTS + nlom*beta*KowTS + wc);  #updated 6/30/2010 & fixed 7/1/2010
      
      k2 <- k1/kbw;	 # k2 = gill elimination rate constant
      
      # calculate feeding rate
      if (taxa==3|taxa==5) 
      {
        Gd <- 0.022 * (Wb^0.85) * exp(0.06*T); # non filter feeders
      }
      if (taxa==2|taxa==4) 
      {
        Gd <- Gv*vss*scav; # filter feeders (kg/d)
      }
      
      kd <- Ed*Gd/Wb;	# calculate kd dietary uptake rate constant (kg/kg*d)
    }
    if (taxa>5) 
    { # seals and birds -- !! NOT READY YET !!
      #fprintf(' DIETCALC : Model not implemented for taxa > 5 ...\n');
      Output$cbiota <- NA; 
    }	
    
    kG <- GR * Wb^-0.2;
    #-----------------
    
    # calculate fecal egestion rate (kg/d)
    Gf <- Gd*((1-assimEff_1)*vld + (1-assimEff_2)*(vnd+vcd) + (1-assimEff_3)*vwd); #updated 6/30/2010
    
    # lipid (vlg), nloc (vcg), nlom (vng), and water (vwg) fraction of gut contents (kg/kg), respectively
    ## Updated 6/30/2010 ##
    vlg <- (1-assimEff_1)*vld/((1-assimEff_1)*vld+(1-assimEff_2)*(vnd+vcd)+(1-assimEff_3)*vwd);
    vng <- (1-assimEff_2)*vnd/((1-assimEff_1)*vld+(1-assimEff_2)*(vnd+vcd)+(1-assimEff_3)*vwd); 
    vcg <- (1-assimEff_2)*vcd/((1-assimEff_1)*vld+(1-assimEff_2)*(vnd+vcd)+(1-assimEff_3)*vwd);
    vwg <- (1-assimEff_3)*vwd/((1-assimEff_1)*vld+(1-assimEff_2)*(vnd+vcd)+(1-assimEff_3)*vwd);
    
    kgb <- (vlg*KowTS/lden + vcg*betap*KowTS + vng*beta*KowTS + vwg)/
      (lipid*KowTS/lden + nlom*beta*KowTS + wc); # gut-biota partition coefficient (unitless) - updated 6/30/2010
    
    ke <- Gf*Ed*kgb/Wb; # fecal egestion rate constant (1/d)
    
    # calculate total concentration available from prey for dietary uptake
    Output$cprey <- as.vector(as.matrix(preyprop) %*% cbiota);
    
    # finally, calculate concentration in organism based on uptake/loss pseudo-equilibrium.
    Output$cbiota <- (k1*(mo*phi*cwater + mp*cpw)+kd*Output$cprey[1])/(k2 + ke + kG + kM);	
    
    # assign results to out variable
    Output$k1  <- k1;
    Output$k2  <- k2;
    #Output$kG  <- kG;
    Output$GR  <- GR;
    Output$Gv  <- Gv;
    Output$Gd  <- Gd;
    Output$Gf  <- Gf;
    Output$assimEff_1 <- assimEff_1;
    Output$assimEff_2 <- assimEff_2;	
    Output$assimEff_3 <- assimEff_3;	
    Output$vlg <- vlg;
    Output$vcg <- vcg;
    Output$vng <- vng;
    Output$vwg <- vwg;
    Output$kgb <- kgb;
    Output$ke  <- ke;
    Output$kd  <- kd;
    
  }
  return(Output)
  
}


# bioaccumulation batch ---------------------------------------------------

#' Run bioaccumulation model function FoodWeb_SQO in batch for all species, contaminants
#'
#' @param biota input biological data
#' @param contam input contaminant data
#' @param biota_preyprop input prey proportions
#' @param constants input constants
#'
bioaccum_batch <- function(biota, contam, biota_preyprop, constants){
  
  ### HAVE HARD-CODED THE PARAMETER NAMES
  theParamNames <- c("k1","k2", "GR", "Gv", "Gd", "Gf", "vlg", "vcg", "vng", "vwg", "kgb", 
                     "ke", "kd", "Ew", "Ed", "phi", "cpw", "assimEff_1", "assimEff_2", "assimEff_3");

  nspecies = nrow(biota_preyprop); 
  nprey = ncol(biota_preyprop);
  ncontam = length(contam$Chem); 
  nparams = length(theParamNames);
  
  #dietary composition (lipid, organic carbon, organic matter, water)
  biota_vld  = as.matrix(biota_preyprop) %*% biota$lipid;
  biota_vcd  = as.matrix(biota_preyprop) %*% biota$nloc;
  biota_vnd  = as.matrix(biota_preyprop) %*% biota$nlom;
  biota_wc = (1 - biota$lipid - biota$nloc - biota$nlom); #wc = (1 -lipid -nloc -nlom)
  biota_vwd  = as.matrix(biota_preyprop) %*% biota_wc; 
  #respiratory uptake of overlying water
  biota_mo = 1 - biota$mp;
  
  ## ASSIGN CONSTANT PARAMETERS from input data##
  EdA = constants$Value[which(constants$Constant=='EdA')];
  EdB = constants$Value[which(constants$Constant=='EdB')];
  xdoc = constants$Value[which(constants$Constant=='xdoc')];
  ddoc = constants$Value[which(constants$Constant=='ddoc')];
  xpoc = constants$Value[which(constants$Constant=='xpoc')];
  dpoc = constants$Value[which(constants$Constant=='dpoc')];
  alphapoc = constants$Value[which(constants$Constant=='alphapoc')];
  alphadoc = constants$Value[which(constants$Constant=='alphadoc')];
  ocsed = constants$Value[which(constants$Constant=='ocsed')]; 
  ds = constants$Value[which(constants$Constant=='ds')];
  A = constants$Value[which(constants$Constant=='A')];
  B = constants$Value[which(constants$Constant=='B')];
  T = constants$Value[which(constants$Constant=='T')];
  Cox = constants$Value[which(constants$Constant=='Cox')];
  vss = constants$Value[which(constants$Constant=='vss')];
  scav = constants$Value[which(constants$Constant=='scav')];
  
  #some descriptions of parameters
  ## that are used in model
  #cbiota = 'contaminant concentration in biota';
  #cprey  = 'contaminant concentration in prey diet';
  #k1 = 'aqueous uptake rate constant';
  #k2 = 'elimination rate constant';
  #GR = 'growth rate coefficient for fish or invertebrates; growth rate for phytoplankton';
  #Gv = 'gill ventilation rate';
  #Gd = 'feeding rate';
  #Gf = 'fecal egestion rate';
  #vlg = 'lipid fraction of gut';
  #vcg = 'nloc fraction of gut';
  #vng = 'nlom fraction of gut';
  #vwg = 'water fraction of gut';
  #kgb = 'gut-biota partition coefficient';
  #ke =  'fecal egestion rate constant (1/d)';
  #kd =  'dietary uptake rate constant';
  #Ew =  'gill chemical uptake efficiency';
  #Ed =  'dietary chemical transfer efficiency (also called gut uptake efficiency)';
  #phi = 'freely dissolved contaminant fraction in overlying water column';
  #cpw = 'contaminant concentration in porewater';
  #assimEff_1  = 'assimilation efficiency for lipid';
  #assimEff_2  = 'assimilation efficiency for nlom or nloc';
  #assimEff_3  = 'assimilation efficiency for water';
  
  ##Initialize output variables data objects (matrices and an array)
  #matrix that stores organism contaminant concentrations
  CBIOTA <- matrix(numeric(nspecies*ncontam), nrow=nspecies,dimnames=list(biota$Biota,contam$Chem));
  #contaminant concentrations of organisms' prey
  CPREY  <- matrix(numeric(nspecies*ncontam), nrow=nspecies,dimnames=list(biota$Biota,contam$Chem));
  #calculated biota sediment accumulation factor (BSAF)
  BSAF <-  matrix(numeric(nspecies*ncontam), nrow=nspecies,dimnames=list(biota$Biota,contam$Chem));
  #array stores each parameter specific to species and contaminant type
  PARAMS <- array(NA, c(nspecies,ncontam,nparams),dimnames=list(biota$Biota,contam$Chem,theParamNames));
  
  #PART 1
  ####Deterministic version###########################################################3
  
  #N.B. the use of for loops is rather inefficient computationally but was easy to code
  #and moderately easy to debug.  Perhaps there is a more efficient way.
  
  # loops over ALL contaminants
  # wrap the loop execution in withProgress
  withProgress(
    message='Please wait',
    detail='Doing important stuff...',
    value = 0, {
      for (icontam in 1:ncontam){
        
        # loop over species
        for (ispecies in 1:nspecies) {
          
          # log
          txt <- paste0('Contaminant', icontam, 'of', ncontam, '\n\tspecies', ispecies, 'of', nspecies)
          
          ## ASSIGN VARIABLES ##
          csed <- contam$cs_ng.g[icontam];  #sediment contaminant concentration 
          Kow <- contam$Kow[icontam]; #octanol-water partitioning coefficient (affects contaminant partitioning)
          KowTS <- contam$KowTS[icontam]; #temperature and salinity corrected octanol water partitioning
          Koc <- alphapoc*KowTS; # from Gobas and Arnot, 2010 Supp. pg. 5 (cit. Seth et at., 1999)
          
          ## Multiple water conc options ##
          #cwater <- contam$cd_ng.g[icontam]; # water concentration from optional input
          #cwater <- csed/KowTS; # simulated water concentration using Kow only
          cwater <- (csed/ocsed)/Koc; # simulated water concentration using %TOC in sediment 
          
          #other parameters, most reviewed above; described in Gobas and Arnot 2010
          beta <- contam$beta[icontam]; 
          betap <- contam$betap[icontam]; 
          kM <- contam$kM[icontam];
          taxa <- biota$taxa[ispecies];
          lipid <- biota$lipid[ispecies];
          nloc <- biota$nloc[ispecies];
          nlom <- biota$nlom[ispecies];
          wc <- biota_wc[ispecies]; 
          mp <- biota$mp[ispecies];
          mo <- biota_mo[ispecies];
          Wb <- biota$Wb[ispecies];
          GR <- biota$GR[ispecies];
          preyprop <- biota_preyprop[ispecies,]; # a vector of prey proportions
          cbiota <- CBIOTA[,icontam];
          vld <- biota_vld[ispecies]; # lipid proportions in diet
          vcd <- biota_vcd[ispecies]; # nonlipid O.C. proportions in diet
          vnd <- biota_vnd[ispecies]; # nonlipid O.M. proportions in diet
          vwd <- biota_vwd[ispecies]; # water proportions in diet
          assimEff_1 <- biota$assimEff_1[ispecies]; # assimilation efficiency for lipid
          assimEff_2 <- biota$assimEff_2[ispecies]; # assimilation efficiency for non-lipid crabon (nloc) or org matter (nlom)
          assimEff_3 <- biota$assimEff_3[ispecies]; # assimilation efficiency for water
          
          ###########################
          ### CALL FOOD WEB MODEL ###
 
          #call the function, using all the various model input parameters
          Results <- FoodWeb_SQO(NumSim=1, csed, cwater, KowTS, Kow, EdA, EdB, xdoc, ddoc, xpoc, dpoc, alphapoc, 
                                       alphadoc, ocsed, ds, taxa, A, B, T, lipid, nloc, nlom, wc, beta, betap, mo, mp, kM, Wb, Cox, 
                                       vss, scav, preyprop, cbiota, vld, vcd, vnd, vwd, GR, assimEff_1, assimEff_2, assimEff_3) 
          
          #extract the biota and prey contaminant information
          CBIOTA[ispecies,icontam] <- Results$cbiota;
          CPREY[ispecies,icontam] <- Results$cprey;
          # check column headers on Results match dimnames on PARAMS 
          for (o in 1:(length(Results)-2)) 
          { if (dimnames(PARAMS)[[3]][o]==colnames(Results)[o+2]) {PARAMS[ispecies,icontam,o] <- Results[[o+2]]} else { print('Error in ordering of parameters')} 
          };  
        } #end of loop over species
      } #end of loop over contaminants
      incProgress(amount=1)
    })
  
  ### CALCULATE BSAF (simple measure of biota versus sediment contamination) ###
  for (i in 1:nspecies) {
    BSAF[i,] <- CBIOTA[i,]/contam$cs_ng.g
  }
  
  out <- list(BSAF = BSAF, CBIOTA = CBIOTA, contam = contam)
  
  return(out)

}



# format species inputs ---------------------------------------------------

#' format species inputs
#'
#' @param inps reactive inputs
#' @param biota table
formsppinp <- function(inps, biota){
  
  # format names and lips as tibble
  frminps <- reactiveValuesToList(inps) %>% 
    enframe('Biota', 'value') %>% 
    filter(!grepl('selectized', Biota) & grepl('indic', Biota)) %>% 
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
    select(-lipid.y, -spp)
  
  return(out)

}

# format constants inputs ---------------------------------------------------

#' format constants inputs
#'
#' @param inps reactive inputs
#' @param constants table
formcnsinp <- function(inps, constants){

  # format names and lips as tibble
  frminps <- reactiveValuesToList(inps) %>% 
    enframe('Constant', 'Value') %>% 
    filter(Constant %in% c('Cox', 'T', 'salinity', 'MCS', 'ocsed', 'vss', 'xpoc', 'xdoc')) %>% 
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
      `cw_ng.g` = sur,   # surface water concentration
      `cp_ng.g` = por    # porewater concentration
    )
  
  # add frminps user input to contaminants
  out <- contam %>% 
    left_join(frminps, by = 'Chem') %>% 
    mutate(
      cs_ng.g.x = ifelse(is.na(cs_ng.g.x), cs_ng.g.y, cs_ng.g.x),
      cd_ng.g.x = ifelse(is.na(cd_ng.g.x), cd_ng.g.y, cd_ng.g.x),
      cw_ng.g.x = ifelse(is.na(cw_ng.g.x), cw_ng.g.y, cw_ng.g.x),
      cp_ng.g.x = ifelse(is.na(cp_ng.g.x), cp_ng.g.y, cp_ng.g.x)
      ) %>% 
    rename(
      cs_ng.g = cs_ng.g.x,
      cd_ng.g = cd_ng.g.x,
      cw_ng.g = cw_ng.g.x,
      cp_ng.g = cp_ng.g.x
      ) %>% 
    select(-cs_ng.g.y, -cd_ng.g.y, -cw_ng.g.y, -cp_ng.g.y)
  
  return(out)
  
}

# format contaminant concentration data as calculated ---------------------

#' Format calculated contaminant concentrations
#'
#' @param contam input contaminants from formatted user inputs
#' @param constants input constants from formatted user inputs
cntcalc <- function(contam, constants){

  # total organic carbon input (ocsed)
  ocsed <- constants %>% 
    filter(Constant %in% 'ocsed') %>% 
    pull(Value)
  
  # salinity input (Sal)
  Sal <- constants %>% 
    filter(Constant %in% 'salinity') %>% 
    pull(Value)
  
  # temperature
  Temp <- constants %>% 
    filter(Constant %in% 'T') %>% 
    pull(Value)
  
  # poc concentration in water
  xpoc <- constants %>% 
    filter(Constant %in% 'xpoc') %>% 
    pull(Value)
  
  # doc concentration in water
  xdoc <- constants %>% 
    filter(Constant %in% 'xdoc') %>% 
    pull(Value)
  
  # disequalibrium factor for poc partitioning
  dpoc <- constants %>% 
    filter(Constant %in% 'dpoc') %>% 
    pull(Value)
  
  # disequalibrium factor for doc partitioning
  ddoc <- constants %>% 
    filter(Constant %in% 'ddoc') %>% 
    pull(Value)
  
  # proportionality constant describing phase partitioning of poc
  alphapoc <- constants %>% 
    filter(Constant %in% 'alphapoc') %>% 
    pull(Value)
  
  # proportionality constant describing phase partitioning of doc
  alphadoc <- constants %>% 
    filter(Constant %in% 'alphadoc') %>% 
    pull(Value)
  
  # log kow temp corrected
  # log kow sal corrected (whic includes temp correction)
  # phi
  # calculated dissolved surface water concentration
  # free dissolved surface water concentration
  # calculated poreweqter concentration
  # log koc
  contam <- contam %>% 
    mutate(
      logkow_tempcor = log10(Kow) - ((delt_uow / (log(10) * 0.0083145)) * ((1 / (273 + Temp)) - (1 / 298))),
      logkow_salcor = log10(10 ^ logkow_tempcor* (10 ^ (0.0018 * LeBas_Molar_Volume * (0.5 * Sal / 35)))),
      phi = 1 / (1 + (xpoc * dpoc * alphapoc * 10 ^ logkow_salcor) + (xdoc * ddoc * alphadoc * 10 ^ logkow_salcor)),
      calc_cd_pg.l = ifelse(!is.na(cs_ng.g), 1e6 * (cs_ng.g / (ocsed * (0.35 * 10 ^ logkow_salcor)) / 8), 0),
      free_cd_ng.ml = ifelse(is.na(cd_ng.g), ifelse(is.na(cw_ng.g), calc_cd_pg.l, ifelse(cw_ng.g * phi <= calc_cd_pg.l, cw_ng.g * phi, calc_cd_pg.l)), ifelse(cd_ng.g <= calc_cd_pg.l, cd_ng.g, calc_cd_pg.l)) / 1e6,
      calc_cp_ng.ml = ifelse(cp_ng.g > 0 & !is.na(cp_ng.g), cp_ng.g / 1e6, ((cs_ng.g / ocsed) / (0.35 * 10 ^ logkow_salcor))), 
      log_koc = log10(0.35 * 10 ^ logkow_salcor), 
      KowTS = logkow_salcor
    )

  return(contam)
  
}
