## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE---------------------------------------------------
library(D3mirt)
knitr::knit_hooks$set(webgl = hook_webgl)

## ---- message = FALSE, results = 'hide', eval = FALSE-------------------------
#  # Load data
#  data("anes0809offwaves")
#  x <- anes0809offwaves
#  x <- x[,3:22] # Remove columns for age and gender

## -----------------------------------------------------------------------------
# Optional: Load the EFA data for this example directly from the package file
load(system.file("extdata/efa.Rdata", package = "D3mirt"))

## -----------------------------------------------------------------------------
# Call to modid()
a <- modid(x, efa = FALSE)
summary(a)

## -----------------------------------------------------------------------------
# Call to modid with increased lower and higher bound
b <- modid(x, efa = FALSE, lower = 1, upper = 1 )
summary(b)

## -----------------------------------------------------------------------------
# Override factor order by reversing columns in the original data frame
c <- modid(x, efa = FALSE, fac.order = c(3, 2, 1))
summary(c)

## ---- message = FALSE, results = 'hide', eval = FALSE-------------------------
#  # Load data
#  data("anes0809offwaves")
#  x <- anes0809offwaves
#  x <- x[,3:22] # Remove columns for age and gender
#  
#  # Fit a three-dimensional graded response model with orthogonal factors
#  # Example below uses Likert items from the built-in data set "anes0809offwaves"
#  # Item W7Q3 and item W7Q20 was selected with `modid()`
#  # The model specification set all items in the data set (1-20)
#  # to load on all three factors (F1-F3)
#  # The START and FIXED commands are used on the two items to identify the DMIRT model
#  spec <-  ' F1 = 1-20
#             F2 = 1-20
#             F3 = 1-20
#  
#             START=(W7Q3,a2,0)
#             START=(W7Q3,a3,0)
#  
#             START=(W7Q20,a3,0)
#  
#             FIXED=(W7Q3,a2)
#             FIXED=(W7Q3,a3)
#  
#             FIXED=(W7Q20,a3) '
#  
#  mod1 <- mirt::mirt(x,
#                     spec,
#                     itemtype = 'graded',
#                     SE = TRUE,
#                     method = 'QMCEM')

## -----------------------------------------------------------------------------
# Optional: Load the mod1 as a data frame directly from the package file
load(system.file("extdata/mod1.Rdata", package = "D3mirt"))

## -----------------------------------------------------------------------------
# Call D3mirt()
g <- D3mirt(mod1)
summary(g) # Show summary of results

## -----------------------------------------------------------------------------
# Call to D3mirt(), including optional nested lists for three constructs
# Item W7Q16 is not included in any construct because of model violations
# The model violations for the W7Q16 item can be seen when plotting the model
c <- list(list(1,2,3,4,5,6,7,8,9,10),
          list(11,12,13,14),
          list(15,17,18,19,20))
g <- D3mirt(mod1, c)
summary(g)

## ---- testg1, webgl = TRUE, fig.width = 7, fig.height = 7---------------------
# Plot RGL device
plot(g, view = c(15, 20, 0.6))

## ---- testg2, webgl=TRUE, fig.width = 7, fig.height = 7-----------------------
# Plot RGL device with constructs visible and named
plot(g, constructs = TRUE, 
        construct.lab = c("Compassion", "Fairness", "Conformity"), 
        view = c(15, 20, 0.6))

## ---- testg3, webgl=TRUE, fig.width = 7, fig.height = 7-----------------------
# Item W7Q16 has location 16 in the data set (gender and age excluded)
# The item is plotted together with constructs to aid the visual interpretation
plot(g, constructs = TRUE, 
        items = 16, 
        construct.lab = c("Compassion", "Fairness", "Conformity"), 
        view = c(15, 20, 0.6))

## ---- testg4, webgl = TRUE, fig.width = 7, fig.height = 7---------------------
# Plot RGL device on item difficulty level 5
plot(g, diff.level = 5, 
        view = c(15, 20, 0.6))

## ---- testg5, webgl=TRUE, fig.width = 7, fig.height = 7-----------------------
# A selection of Conformity items from the model plotted with constructs
plot(g, constructs = TRUE, 
        items = c(15,17,18,19,20), 
        construct.lab = c("Compassion", "Fairness", "Conformity"), 
        view = c(15, 20, 0.6))

## ---- testg6, webgl=TRUE, fig.width = 7, fig.height = 7-----------------------
# Plot RGL device with items in uniform length and constructs visible and named
plot(g, scale = TRUE, 
        constructs = TRUE, 
        construct.lab = c("Compassion", "Fairness", "Conformity"), 
        view = c(15, 20, 0.6))

## ---- eval = FALSE------------------------------------------------------------
#  # Extract respondent factor scores from mod1 with fscores() function
#  f <- mirt::fscores(mod1,
#                     method="EAP",
#                     full.scores = TRUE,
#                     full.scores.SE = FALSE, QMC = TRUE)

## -----------------------------------------------------------------------------
# Optional: Load the respondent factor scores for this example directly from the package file
load(system.file("extdata/fscores.Rdata", package = "D3mirt"))

## -----------------------------------------------------------------------------
# Attach f to the gender variable (column 2 from anes0809offwaves data set; "W3XGENDER")
# Use cbind with fscores() output attached first
data("anes0809offwaves")
x <- anes0809offwaves
z <- data.frame(cbind(f, x[,2]))

## ---- testg7, webgl=TRUE, fig.width = 7, fig.height = 7-----------------------
# Plot profiles with item vector arrows hidden
# Score levels: 1 = Blue ("male") and 2 = Red ("female")
plot(g, hide = TRUE, 
        profiles = z, 
        levels = z[,4], 
        sphere.col = c("blue", "red"), 
        x.lab = "Compassion", 
        y.lab="Conformity", 
        z.lab="Fairness", 
        view = c(16, 20, 0.6))

## -----------------------------------------------------------------------------
# Column bind fscores() with age variable ("W3Xage")
y <- data.frame(cbind(f, x[,1]))

# Subset data frame y conditioned on age <= 30
z1 <- subset(y, y[,4] <= 30)

# Subset data frame y conditioned on age >= 70
z2 <- subset(y, y[,4] >= 70)

# Row bind z1 and z2
z <- rbind(z1,z2)

# Check number of factor levels with nlevels() and as.factor()
nlevels(as.factor(z1[,4]))
nlevels(as.factor(z2[,4]))

# Use rep() to create a color vector to color groups based on the nlevels() output
# z1 has 14 factor levels and z2 has 16 factor levels
# z1 respondents are colored red and z2 are colored blue
colvec <- c(rep("red", 14), 
            rep("blue", 16))

## ---- testg8, webgl=TRUE, fig.width = 7, fig.height = 7-----------------------
# Call plot() with profile data on age with item vector arrows hidden
plot(g, hide = TRUE, 
     profiles = z, 
     levels = z[,4], 
     sphere.col = colvec, 
     x.lab = "Compassion", 
     y.lab="Conformity", 
     z.lab="Fairness", 
     view = c(15, 20, 0.6))

## -----------------------------------------------------------------------------
# Column bind fscores() with age variable ("W3Xage")
y <- data.frame(cbind(f, x[,1]))

# Subset data frame y conditioned on age <= 30
z1 <- subset(y, y[,4] <= 30)

# Use rep() to create a color vector to color groups based on the nlevels() output
# z1 has 14 factor levels
colvec <- c(rep("red", 14))

## ---- testg9, webgl=TRUE, fig.width = 7, fig.height = 7-----------------------
# Call plot() with profile data on age with item vector arrows hidden
plot(g, hide = TRUE, 
     profiles = z1, 
     levels = z1[,4], 
     sphere.col = colvec, 
     x.lab = "Compassion", 
     y.lab="Conformity", 
     z.lab="Fairness", 
     ci = TRUE, 
     ci.level = 0.95, 
     ellipse.col = "orange",
     view = c(15, 20, 0.6))

## ---- eval = FALSE------------------------------------------------------------
#  # Export an open RGL device to the console that can be saved as an html or image file
#  plot(g, constructs = TRUE)
#  s <- scene3d()
#  rgl::rglwidget(s,
#                 width = 1040,
#                 height = 1040)
#  
#  # Export a snap shoot of an open RGL device directly to file
#  plot(g, constructs = TRUE)
#  rgl::rgl.snapshot('RGLdevice.png',
#                      fmt = 'png')

