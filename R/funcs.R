#FoodWeb_SQO_v2.R
# Purpose: To implement the Gobas/Arnot Food Web Bioaccumulation Model
#          of San Francisco Bay for various species and contaminants.
# 6.17.2010
# Michelle Lent
# edited 6.23.10 & 6.30.10 by BG and ML


FoodWeb_SQO_v2MCS <- function(NumSim, csed, cwater, KowTS, Kow, EdA, EdB, xdoc, 
                              ddoc, xpoc, dpoc, alphapoc, alphadoc, ocsed, ds, taxa, A, B, T, lipid, nloc, 
                              nlom, wc, beta, betap, mo, mp, kM, Wb, Cox, vss, scav, preyprop, cbiota, 
                              vld, vcd, vnd, vwd, GR, assimEff_1, assimEff_2, assimEff_3) #updated 6.30.2010
{
  
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