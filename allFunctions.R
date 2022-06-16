library(CHNOSZ)
library(stringr)
library(knitr)
library(ggplot2)
internalEnergy <- function(dq,dw,dv, p)
{
      # Energy change in a system
      ΔE <- dq - dw       # change in energy by work or heat
      ΔE <- dq - p*dv     # change in work by change in volume
}

systemEnthalpy <- function(dq,dw,dp,dv,p,v)
{
      # find ΔE
      ΔE <- internalEnergy(dq,dw,dv,p)
      
      # Enthalpy change in a system
      ΔH <- ΔE + (p*dv) + (v*dp)
}

systemEntropy <- function(dq,dw,dp,dv,p,v,t)
{
      if(length(v)==0){
            #using ΔE
            ΔE <- internalEnergy(dq,dw,dv,p)
            ΔS <- (ΔE+(p*dv))/t
      }
      else{
            #using ΔH
            ΔH <- systemEnthalpy(dq,dw,dp,dv,p,v)
            ΔS <- (ΔH-(v*dp))/t
      }
}

systemConditions <- function(dq,dw,dv,dp,q,w,v,p,t) 
{
      ΔE <- internalEnergy(dq, dw, dv, p)
      print(paste("The change in energy is", ΔE)) 
      
      ΔH <- systemEnthalpy(dq, dw, dp, dv, p, v)
      print(paste("The change in enthalpy is", ΔH))
      
      ΔS <- systemEntropy(dq,dw,dp,dv,p,v,t)
      print(paste("The change in entropy is", ΔS)) 
}

GibbsFreeEnergy <- function(dq,dw,dv,dp,q,w,v,p,t)
{
      ΔE <- internalEnergy(dq, dw, dv, p)
      ΔH <- systemEnthalpy(dq, dw, dp, dv, p, v)
      ΔS <- systemEntropy(dq,dw,dp,dv,p,v,t)
      
      ΔG <- ΔH - tΔS
      print(paste("Therefore, the change in Gibbs free energy is", ΔG, "."))
      
      if(ΔG < 0){print("This process is spontaneous.")}
      if(ΔG > 0){print("This process is not spontaneous.")}
      if(ΔG == 0){print("This process is at equilibrium.")}
}

activityCoefficient <- function(a,m)
{
      γ=a/m  
}

chemicalPotential <- function(reactant,molesR,productA,molesA,productB,molesB,productC,molesC,Temp,Pres)
{
      E.units("J")
      
      reactant <- data.frame(subcrt(reactant, T = Temp, P = Pres))
      print(reactant)
      nameReactant <- paste("reactant$out.", reactant$species.name, ".G",sep ="")
      ΔGr <- eval(parse(text = nameReactant))/1000
      
      a <- data.frame(subcrt(productA, T = Temp, P = Pres))
      nameA <- paste("a$out.", str_replace(productA, "[\\+\\-\\*]", "."), ".G",sep ="")
      ΔGa <- eval(parse(text = nameA))/1000
      
      b <- data.frame(subcrt(productB, T = Temp, P = Pres))
      nameB <- paste("b$out.", str_replace(productB, "[\\+\\-\\*]", "."), ".G",sep ="")
      ΔGb <- eval(parse(text = nameB))/1000
      
      c <- data.frame(subcrt(productC, T = Temp, P = Pres))
      nameC <- paste("c$out.", str_replace(productC, "[\\+\\-\\*]", "."), ".G",sep ="")
      ΔGc <- eval(parse(text = nameC))/1000
      
      print(paste("Reactants","⇌", "Products"))
      print(paste(reactant$species.name, "⇌", molesA, productA,"+", molesB, productB,"+", molesC, productC,sep=" "))
      print("Products - Reactants")
      print("(ΔGa + ΔGb + ΔGc) - (ΔGr)")
      print(paste("(", ΔGa, "+", ΔGb, "+", ΔGc, ") - (", ΔGr, ")", sep=""))
      ΔGtotal <- ((molesA*ΔGa)+(molesB*ΔGb)+(molesC*ΔGc))-(molesR*ΔGr)
      print(paste("ΔG =", ΔGtotal, "kJ/mol"))
      return(ΔGtotal)
}

equilibriumConstant <- function(reactant,molesR,productA,molesA,productB,molesB,productC,molesC,Temp,Pres)
{
      ΔG <- chemicalPotential(reactant,molesR,productA,molesA,productB,molesB,productC,molesC,Temp,Pres)
      
      logKsp <- (-1*ΔG)/5.708
      print(paste("Ksp = 10^", logKsp))
}

ionicStrength <- function(species = data.frame())
{
      IS <- is.integer(0)
      for(i in 1:nrow(species))
      {
            m <- species[i, "concentration"]
            z <- species[i, "charge"]
            output <- m * z^2
            IS <- IS + output
      }
      solutionIonicStrength <- 0.5*IS
      cat("The solution ionic strength is", solutionIonicStrength)
}

activity <- function(y,m)
{
      a <- y*m  
}

printEverything <- function(...)
{
      ## Mineral dissolving in water
      cat("MINERALS DISSOLVING IN WATER\n")
      cat("Reactants ⇌ Products\n")
      cat("  ",reactant, "⇌", molesA, productA,"+", molesB, productB,"+", molesC, productC, "\n",sep=" ")
      cat("\nGibbs Free Energy Calculations\n")
      cat("  ΔG = Products - Reactants")
      cat("\n  ΔG = (ΔGa + ΔGb + ΔGc) - (ΔGr)\n")
      cat("  ΔG = (", ΔGa, "+", ΔGb, "+", ΔGc, ") - (", ΔGr, ")", sep="")
      cat("  ΔG =", ΔGtotalSolid, "kJ/mol\n")
      if(ΔGtotalSolid < 0){cat("   This process is spontaneous.\n")}
      if(ΔGtotalSolid > 0){cat("   This process is not spontaneous.\n")}
      if(ΔGtotalSolid == 0){cat("   This process is at equilibrium.\n")}
      cat("\nK Equilibrium Calculation\n")
      cat("  Ksp = 10^", logKsp,"\n")
      
      ## Aqueous complex
      cat("\nAQUEOUS COMPLEX STABILITY\n")
      cat("Reactants ⇌ Products\n")
      cat("  ",cation," + ", anion, "⇌", aqComplex, "\n",sep=" ")
      cat("\nGibbs Free Energy Calculations")
      cat("\n  ΔG = Products - Reactants")
      cat("\n  ΔG = (ΔGaq)-(ΔGcation + ΔGanion)\n")
      cat("  ΔG = (", ΔGaq, ") - (", ΔGcation, "+", ΔGanion, ")", sep="")
      cat("  ΔG =", ΔGtotalAq, "kJ/mol\n")
      if(ΔGtotalAq < 0){cat("   This process is spontaneous.\n")}
      if(ΔGtotalAq > 0){cat("   This process is not spontaneous.\n")}
      if(ΔGtotalAq == 0){cat("   This process is at equilibrium.\n")}
      cat("\nAqueous Complex Stability")
      cat("\n Kstab = 10^", logKstab)
      cat("\n New solubility = ", solubility)
}
