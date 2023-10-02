# ==============================================================================
# === openMSE Exercise 2: Constructing Operating Models  =======================
# ==============================================================================
#
# In this practical exercise, we will explore four ways to construct openMSE
# operating models (OMs):
#
# 1. Using built-in Examples
#
# 2. Specify OM with expert judgement
#
# 3. Import a stock assessment
#
# 4. Condition OM on available fishery data
#


# === Setup (only required if new R session) ===================================
#
# Load the openMSE package
library(openMSE)

# Set the working directory to the Exercises folder.
#
# Replace with the path to the Exercises folder on your machine:
setwd('G:\\R Exercises')
getw

# === Section 1: Using built-in Examples =======================================

# We've already seen the OMs in the openMSE library:
library(MSEextra)
avail('OM')

# In Exercise 1 we also built an OM from pre-built components using the `new`
# function:

AlbOM <- new('OM', Albacore, FlatE_NDom, Generic_Obs, Perfect_Imp)

# This can be a useful approach for rapidly building operating models, especially
# for data-limited fisheries where there is insufficient data to condition the
# model.

# This OM is built from four sub-components:

# 1. An object of class Stock
class(Albacore)

# 2. An object of class Fleet
class(FlatE_NDom)

# 3. An object of class Obs
class(Generic_Obs)

# 4. An object of class Imp
class(Perfect_Imp)

# The Stock object contains all the information related to the biological
# characteristics of the stock:
Albacore
slotNames(Albacore)

# Most of the slots in a Stock object require two values: lower and upper bounds
# of a uniform distribution

# For example, the natural mortality rate for the Albacore stock is a uniform
# distribution between 0.35 and 0.45:
Albacore@M

# A value will be sampled from this range for each simulation in the model.
# So for example, our OM has 48 simulations:
AlbOM@nsim

# that means 48 values of M will be drawn from this distribution.

# The same thing applies to the other biological parameters
Albacore@Linf
Albacore@K

# Later, we'll learn how to generate correlated life-history parameters, or
# over-ride the default assumption of a uniform distribution.

# It's easy to change the values:

AlbOM@M <- c(0.2, 0.2) # make M values constant across simulations
AlbOM@M <- c(0.1, 0.4) # increase the range (more uncertainty)

# The same principle applies to the other objects:
FlatE_NDom@L5 # first length at 5% selection
FlatE_NDom@LFS # first length at full selection
FlatE_NDom@Vmaxlen # vulnerability of largest size class

Generic_Obs@Cobs  # SD for log-normal observation error for catches

Perfect_Imp@TACFrac # Fraction of TAC overages (none)

# There are a number of objects of each class available in openMSE:
avail('Stock')
avail('Fleet')
avail('Obs')
avail('Imp')

# The individual objects can be plotted:
plot(Albacore)
plot(DecE_HDom, Albacore)
plot(IncE_HDom, Albacore)

# As we've seen before you can also plot the OM, which will produce a report
# with plots of the individual components and the simulated fishery dynamics:
plot(AlbOM)


# You may be now confused by a large number of plots that include
# names you do not understand. In most cases these plots refer to
# operating model parameters that are slots in the objects. You
# can find out more about these objects and what the slots mean by using
# the package documentation that comes built-in to openMSE
#
# For example to learn more about objects of class Stock you can use
# the question mark operator:

class?Stock

# You can also find this information by browsing the MSEtool documentation page:
# e.g.: https://msetool.openmse.com/reference/Stock-class.html

# Help documentation for the other objects is available in the same way:
class?Fleet
class?Obs
class?Imp


# These objects can be stitched together to quickly build operating models.

# This can be a useful way to get started with building operating models,
# especially for data-limited fisheries.

# You can select objects that most closely match your fishery, and modify the
# individual parameters based on your knowledge of the fishery.

# Blue shark with decreasing effort with heavily dome-shaped selectivity:
Blue_shark_1 <- new('OM', Blue_shark, DecE_HDom, Imprecise_Unbiased, Overages)

# Now with increasing effort and same selectivity pattern:
Blue_shark_2 <- new('OM', Blue_shark, IncE_HDom, Imprecise_Unbiased, Overages)


# === Section 2: Specify OM with expert judgement ======================
#
# Another way to create an Operating Model is to import the parameters from an
# Excel workbook (CSV files also work).
#
# Again, this is most useful for data-limited fisheries where you are not
# able to condition an OM with data.
#
# For example, the Exercises/Data directory contains an OM for the Lane Snapper
# (Lutjanus synagris) from the Gulf of Mexico.
#
# Open the Excel workbook and inspect the values in the Stock, Fleet, Obs, Imp,
# and OM tabs. Notice how the names in each tab match up to the slot names in
# the corresponding objects.
#
# The XL2OM function is used to import the OM from Excel:

Lane_Snapper <- XL2OM('Data/Lane_Snapper/Lane_Snapper_GOM_NOAA')

avail('OM')
# Confirm the parameters in the Lane_Snapper OM object match those in the
# spreadsheet; e.g.:
# Stock parameters
Lane_Snapper@M

# Fleet parameters
Lane_Snapper@nyears

# If time allows, you could make a copy of the Lane_Snapper_GOM_NOAA.xlsx file
# and modify some of the parameters in it and import the new OM into the R
# session.

# You will notice that the Data/Lane_Snapper directory also contains a HTML
# document. This document is a companion OM report that describes the rationale
# behind the parameter choices for the operating model.
#
# It is important that operating models are described in detail. This
# provides a full record of the decisions and choices that went into creating
# the OM, and provides a transparent and reproducible approach that others can
# follow.

# You can generate a blank Excel OM file with the `OMinit` function:
dir.create('MyOM')
OMinit('MyOM', dir='MyOM')

# This function creates `MyOM.xlsx` and `MyOM.rmd` in `MyOM` directory.

# `MyOM.xlsx` is a spreadsheet with sheets corresponding to the components of an
# OM: `Stock`, `Fleet`, `Obs`, and `Imp`, `and` OM worksheets.

# To populate the OM, open `MyOM.xlsx` in Excel or any other suitable software,
# and fill in the values for each parameter.
# Once complete, you can import the OM using `XL2OM` as above.

# The `MyOM.rmd` file can be opened in any text editor or RStudio, and
# contains a skeleton for the OM documentation.
#
# The `OMinit` function also creates several folders in the working directory:
# `data`, `docs`, `images`, and `robustness.` These sub-directories can be used to
# store data, documents, images, and other information that is reference in
# the OM Report.

# You can learn more about generating data-limited OMs with expert judgement
# on the openMSE website:
# https://openmse.com/om-data-limited/


# === Section 3: Import a stock assessment ===========================

# When an age-structured stock assessment model exists for a fishery, it is
# possible to import this assessment model directly into openMSE as an operating
# model for closed-loop simulation testing.

# openMSE includes several functions for importing assessment models:

# ASAP (Age-structured Assessment Model; NOAA):
?ASAP2OM

# Awatea (version of Coleraine; University of Washington):
?Awatea2OM

# iSCAM (integrated statistical catch-at-age model); UBC, DFO)
?iSCAM2OM

# Stock Synthesis 3 (NOAA)
?SS2OM

# WHAM (Wood's Hole Assessment Model)
?WHAM2OM

# There is also a generic `Assess2OM` function.
# This function can be used to import the results of any stock assessment,
# provided you have estimates of numbers-at-age, F-at-age, and the life-history
# parameters:

?Assess2OM



# The Data folder has a SS3 model for the Strait of Gibraltar Blackspot Seabream
# (Pagellus bogaraveo)

# This model has two stocks (females and males) and two fishing fleets (Morocco
# and Spain)

# We can import this assessment into openMSE with the `SS2MOM` function
BSB_MOM <- SS2MOM("Data/Blackspot_Seabream/SBR_ss3ref")

# `BSB_MOM` is a object of class `MOM` (multi-operating model)
class(BSB_MOM)

# It has two stocks:
length(BSB_MOM@Stocks)
# and two fleets:
length(BSB_MOM@Fleets[[1]])

# Let's compare the simulated dynamics from the imported operating model,
# to the predictions in the SS3 model:

plot_SS2MOM(BSB_MOM,
            SSdir="Data/Blackspot_Seabream/SBR_ss3ref",
            filename='BSB_MOM',
            dir="Data/Blackspot_Seabream")

# The plots in the report show that the simulated fishery exactly matches the
# values from the SS3 model. That's reassuring!


# `MOM` objects are used for multi-stock and/or multi-fleet MSE.
# They are essentially made up of nested lists of Stock and Fleet objects,
# with some extra stuff to deal with interactions between stocks and/or fleets.

# We're not going to cover multi-stock/fleet MSE in this course. You can learn
# more about the multi-MSE features of openMSE here:
# https://openmse.com/features-multimse/


# The `SS2OM` function will import a SS3 model and collapse it down to a single
# stock and fleet:
BSB_OM <- SS2OM("Data/Blackspot_Seabream/SBR_ss3ref")

class(BSB_OM)

# This approach calculates the overall selectivity and exploitation pattern,
# and describes this as a single fleet.

# The biological parameters from the two stocks (growth, natural mortality, etc)
# are averaged into a single stock.

# Let's compare the simulated dynamics from the single-stock OM with the SS3
# model:
plot_SS2OM(BSB_OM,
            SSdir="Data/Blackspot_Seabream/SBR_ss3ref",
            filename='BSB_OM',
            dir="Data/Blackspot_Seabream")

# The report shows that the fishery dynamics from the single stock/fleet OM
# don't quite match the predictions from SS3 (main differences is absolute
# scale of the catches due to different empirical weight-at-age between sexes).
#
# But it's close enough for our purposes!


# Let's examine the `BSB_OM` object in a bit more detail.

# The number of simulations is set in the call to `SS2OM` (default 48)
BSB_OM@nsim

# Maximum age
BSB_OM@maxage

# Ages start at age-0, so 20 age classes:
0:BSB_OM@maxage

# 39 historical years
BSB_OM@nyears

# 50 projection years
BSB_OM@proyears



# Instead of using the individual slots in the `Stock` and `Fleet`
# sub-components of the `OM` object, the fishery dynamics from imported OMs
# use the Custom Parameters (cpars) feature.

# Custom Parameters is a powerful feature that allows you to directly specify
# all the internal workings of the operating model.This allows lots of
# flexibility for more complex dynamics, such as time-varying M or
# changes in the historical selectivity or retention pattern.

# The `validcpars` function produces a searchable table with all the parameters
# that can be passed in through `cpars`
validcpars()

# We'll look at a few of these now, but won't have time to go through them
# all.

# You can learn more about Custom Parameters here:
# https://openmse.com/features-custom-parameters/

# We can examine the names of the variables in the `cpars` slot:
names(BSB_OM@cpars)

# M-at-Age array:
dim(BSB_OM@cpars$M_ageArray)
# 48 simulations, 20 age classes, 89 years

# Recruitment devations:
dim(BSB_OM@cpars$Perr_y)
# 48 simulations, nyears+proyears+maxage years

# As we've seen before, the `OM` object contains all of the information for the
# historical fishery and the projection period:
Years <- c(rev(seq(BSB_OM@CurrentYr, by=-1, length.out=BSB_OM@nyears+BSB_OM@maxage)),
           seq(BSB_OM@CurrentYr+1, by=1, length.out=BSB_OM@proyears))

rec_devs <- data.frame(Year=rep(Years,each=BSB_OM@nsim),
                       Sim=1:BSB_OM@nsim,
                       Dev=as.vector(BSB_OM@cpars$Perr_y))

rec_devs$Sim <- factor(rec_devs$Sim)

library(ggplot2)
ggplot(rec_devs, aes(x=Year, y=Dev, color=Sim)) +
  geom_line(alpha=0.5) +
  theme_bw() +
  guides(color='none') +
  labs(y='Recruitment Deviation')

# Notice that the recruitment deviations are identical for the historical period
# and are generated for the projection period.

# The recruitment deviations were generated from a log-normal distribution using
# the standard deviation and auto-correlation from the historical
# recruitment deviations

# plot a single simulation
library(dplyr)
sim <- 10
ggplot(rec_devs %>% filter(Sim==sim), aes(x=Year, y=Dev)) +
  geom_line(alpha=0.5) +
  geom_vline(xintercept=BSB_OM@CurrentYr, linetype=3) +
  theme_bw() +
  guides(color='none') +
  labs(y='Recruitment Deviation')


# SD - sigmaR (specified in the SS3 assessment):
BSB_OM@Perr

# lag-1 auto-correlation (calculated from the historical deviations)
BSB_OM@AC


# Note that the `SS2OM` function has arguments to specify the observation and
# implementation parameters:
args(SS2OM)

?SS2OM

# By default, these are:
# Obs = Generic_Obs
# Imp = Perfect_Imp

# You can manually change these values, e.g., assume no observation
# error on catches:
BSB_OM@Cobs
BSB_OM@Cobs <- c(0.01, 0.1)
BSB_OM@Cbiascv
BSB_OM@Cbiascv <- 0.001

# Or replace the Observation object with something else:
BSB_OM <- Replace(BSB_OM, Perfect_Info)
BSB_OM@Cobs

# Later we'll see how to condition the observation parameters using the
# available fishery data.

# === Section 4: Condition OM on available fishery data=========================

# The final way to generate operating models is to use the available fishery
# data to condition an OM using the Rapid Conditioning Model (RCM).

# Documentation for the RCM is available here:
# https://openmse.com/tutorial-rcm/

# We're going to spend most of the last day of the course looking at RCM in
# detail.

# For now, we'll have a quick look just to see how it works.

# To use `RCM` we need two objects:
# 1. A fishery data object
# 2. A OM object

# The fishery data object can either be an object of class `Data`, which we've
# looked at before, or a `RCMdata` object.
# `RCMdata` objects have a bit more flexibility for dealing with data from
# multiple fleets.

# We'll import an RCMdata object containing the data for the Blackspot Seabream:
BSB_data <- readRDS('Data/Blackspot_Seabream/BSB.rcmdata.gz')

# This object contains catch, indices, catch-at-length data, etc for each fleet:
slotNames(BSB_data)

matplot(BSB_data@Chist, type='l', xlab='Year', ylab='Catch')
matplot(BSB_data@Index, type='l', xlab='Year', ylab='Index')

class?RCMdata

# The OM object containing the values for the life-history parameters (Linf, K,
# M, etc) is needed.

# We'll use the OM we generated earlier when we imported the SS3 model.

# Create a new blank OM:
BSB_OM_RCM <- new('OM')

# Replace all the Stock parameters from those in the OM generated from SS3:
BSB_OM_RCM <- Replace(BSB_OM_RCM, BSB_OM, 'Stock')

# We also need to provide the Observation and Implementation objects to the OM
# These aren't used in the RCM fitting and we can replace them later if we like:
BSB_OM_RCM <- Replace(BSB_OM_RCM, Generic_Obs)
BSB_OM_RCM <- Replace(BSB_OM_RCM, Perfect_Imp)

# Give our OM an informative name (used for plots and reports):
BSB_OM_RCM@Name <- 'Blackspot seabream - RCM'

# This OM object now only contains our assumed life-history parameters:
# Growth:
BSB_OM_RCM@Linf
BSB_OM_RCM@K
BSB_OM_RCM@t0

# M is a single value:
BSB_OM_RCM@M

# We'll add the age-dependent M back in with cpars:
BSB_OM_RCM@cpars$M_ageArray <- BSB_OM@cpars$M_ageArray

# The maturity parameters are not populated:
BSB_OM_RCM@L50
BSB_OM_RCM@L50 + BSB_OM_RCM@L50_95

# Add the Maturity-at-Age back in to `cpars`:
BSB_OM_RCM@cpars$Mat_age <- BSB_OM@cpars$Mat_age

# Assumed sigma R:
BSB_OM_RCM@Perr

# Finally, we need to provide the number of historical years:
BSB_OM_RCM@nyears <- BSB_OM@nyears
BSB_OM_RCM@nyears

# And some starting values for the fishery selectivity:

# Plot overall length composition:
barplot(apply(BSB_data@CAL, 2, sum), names.arg=BSB_data@length_bin)

# Set initial values for selectivity parameters:
BSB_OM_RCM@L5 <- c(22,22)
BSB_OM_RCM@LFS <- c(34,34)
BSB_OM_RCM@Vmaxlen <- c(1,1)

# Now that we have our fishery data and our OM containing life-history parameters
# we can run RCM:
# (here we are mirroring the selectivity for the indices to the fishing fleets,
# and saying the selectivity for Fleet 1 is logisitic and Fleet 2 may be
# dome-shaped)

BSB_RCMfit <- RCM(BSB_OM_RCM, BSB_data, s_selectivity=c(1,2),
              resample=F, selectivity=c("logistic","dome"))

# RCM returns an object of class `RCModel`:
class(BSB_RCMfit)

# Generate a report for the fitted RCM model:
plot(BSB_RCMfit, filename='BSB_RCM', dir = 'Data/Blackspot_Seabream')


# The conditioned OM is available in the @OM slot of the `RCModel` object:
BSB_OM_RCM <- BSB_RCMfit@OM

# Let's spool-up the historical fishery dynamics for the two models and compare
# the results:

BSB_SS_Hist <- Simulate(BSB_OM) # the SS3 assessment
BSB_RCM_Hist <- Simulate(BSB_OM_RCM) # our RCM-conditioned OM

nyears <- BSB_OM@nyears
Years <- seq(BSB_OM@CurrentYr, by=-1, length.out=nyears) %>% rev()

DF <- data.frame(Year=Years,
                   F=c(BSB_SS_Hist@SampPars$Fleet$Find[1,],
                       BSB_RCM_Hist@SampPars$Fleet$Find[1,]),
                   B=c(rowSums(BSB_SS_Hist@TSdata$Biomass[1,,]),
                       rowSums(BSB_RCM_Hist@TSdata$Biomass[1,,])),
                   SB=c(rowSums(BSB_SS_Hist@TSdata$SBiomass[1,,]),
                       rowSums(BSB_RCM_Hist@TSdata$SBiomass[1,,])),
                   SB0=c(rep(BSB_SS_Hist@Ref$ReferencePoints$SSB0[1], nyears),
                         rep(BSB_RCM_Hist@Ref$ReferencePoints$SSB0[1], nyears)),
                   SBMSY=c(rep(BSB_SS_Hist@Ref$ReferencePoints$SSBMSY[1], nyears),
                          rep(BSB_RCM_Hist@Ref$ReferencePoints$SSBMSY[1], nyears)),
                   Model=c(rep('SS3', nyears), rep('RCM', nyears)))

DF <- DF %>% mutate(SB_SB0=SB/SB0,
                    SB_SBMSY=SB/SBMSY)
DF$SB0 <- DF$SBMSY <- NULL
DF <- DF %>% tidyr::pivot_longer(., cols=c(2:4, 6:7))

ggplot(DF, aes(x=Year, y=value, color=Model)) +
  facet_wrap(~name, scales='free') +
  geom_line() +
  theme_bw()

# Note that the RCM model assumes a single-sex population

# The `compare_RCM` function can be used to compare runs of RCM with different
# assumptions/data:
?compare_RCM()

# We'll explore this and more features of RCM on the last day of the course.



