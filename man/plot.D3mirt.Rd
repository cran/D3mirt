% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotD3mirt.R
\name{plot.D3mirt}
\alias{plot.D3mirt}
\title{Plot Method for Objects of Class \code{D3mirt}}
\usage{
\method{plot}{D3mirt}(
  x,
  scale = FALSE,
  hide = FALSE,
  ind.scores = FALSE,
  diff.level = NULL,
  items = NULL,
  item.names = TRUE,
  item.lab = NULL,
  constructs = FALSE,
  construct.lab = NULL,
  adjust.lab = c(0.5, -0.8),
  x.lab = "X",
  y.lab = "Y",
  z.lab = "Z",
  title = "",
  line = -5,
  font = 1,
  cex = 1,
  font.col = "black",
  axis.scalar = 1.1,
  axis.length = NULL,
  axis.points = "black",
  axis.ticks = TRUE,
  nticks = 4,
  width.rgl.x = 1040,
  width.rgl.y = 1040,
  view = c(15, 20, 0.6),
  show.plane = TRUE,
  plane.col = "grey80",
  background = "white",
  type = "rotation",
  col = c("black", "grey20", "grey40", "grey60", "grey80"),
  arrow.width = 0.6,
  n = 20,
  theta = 0.2,
  barblen = 0.03,
  c.scalars = c(1, 1),
  c.type = "rotation",
  c.col = "black",
  c.arrow.width = 0.6,
  c.n = 20,
  c.theta = 0.2,
  c.barblen = 0.03,
  profiles = NULL,
  levels = NULL,
  sphere.col = c("black", "grey20", "grey40", "grey60", "grey80"),
  spheres.r = 0.05,
  ci = FALSE,
  ci.level = 0.95,
  ellipse.col = "grey80",
  ellipse.alpha = 0.2,
  ...
)
}
\arguments{
\item{x}{A S3 object of class \code{D3mirt}.}

\item{scale}{Logical, if item vector arrow length should visualize the MDISC. If set to FALSE, the vector arrow length will be of one unit length. The default is \code{scale = FALSE}.}

\item{hide}{Logical, if items should be plotted. The default is \code{hide = FALSE}.}

\item{ind.scores}{Logical, should output plot respondents trait scores. The default is \code{ind.scores = FALSE}.}

\item{diff.level}{Optional. Plotting of a single level of difficulty indicated by an integer.}

\item{items}{Optional. The user can input a list of integers indicating what item vector arrows will be visible while the remaining item vector arrows are hidden.}

\item{item.names}{Logical, if item labels should be plotted. The default is \code{item.names = TRUE}.}

\item{item.lab}{Optional. String vector of item names that will override row names extracted from the data frame. Note that row names are not overwritten.
Instead, the string vector in \code{item.lab} prints item labels on the item vector arrows currently displayed following the order of item vector arrows in the graphical output.
For example, when plotting in the default mode (plotting all item vectors) the labels will follow the order of the items in the data frame.
If a selection of items is plotted with \code{items}, e.g., \code{items = c(24,34,25)}, then the item labels will be displayed following the order of the vector in \code{items} left to right.
In this case, item label 1 will be printed on item 24, item label 2 on item 34, and item label 3 on item 25, and so on.}

\item{constructs}{Logical, if construct vector arrows should be plotted. The default is \code{constructs = FALSE}.}

\item{construct.lab}{Optional. String vector of names for constructs, similar to \code{item.lab}.}

\item{adjust.lab}{Vector of parameters for the position of the item and construct labels for the \code{text3d} function.
The first value is for horizontal adjustment, and the second is for vertical adjustment. The default is \code{adjust.lab = c(0.5, -0.8)}.}

\item{x.lab}{Labels for x-axis, the default is \code{x.lab = "X"}.}

\item{y.lab}{Labels for y-axis, the default is \code{y.lab = "Y"}.}

\item{z.lab}{Labels for y-axis, the default is \code{z.lab = "Z"}.}

\item{title}{The main title for the graphical device, plotted with the \code{title3d()} function. The default is no title.}

\item{line}{Title placement for \code{title3d()}. The default is \code{line = -5}.}

\item{font}{A numeric font number from 1 to 5, the default is \code{font = 1}. See \link[rgl:texts]{rgl::text3d} for more on font options.}

\item{cex}{A numeric character expansion value to adjust font size, the default is \code{cex = 1}.}

\item{font.col}{Color of axes, numbers, and fonts. The default is \code{font.col = "black"}.}

\item{axis.scalar}{Scalar multiple for adjusting the length of all axes (x, y, z) in the 3D model proportionally. The default is \code{axis.scalar = 1.1}.}

\item{axis.length}{Optional. For adjusting the length of the axis manually by entering a numeric or a numeric vector.
For instance, c(3,2,4,3,3,2) indicate axis coordinates x = 3, -x = 3, y = 4, -y = 3, z = 3, -z = 2.
Note that a symmetric model can be created easily by adding a single numeric in the \code{axis.length} argument (e.g., \code{axis.length = 4}) because the function repeats the last value in the vector to cover all axes points.
The default is \code{axis.length = NULL}.}

\item{axis.points}{Color of axis points for the \code{points3d()} function. The default is \code{axis.points = "black"}.}

\item{axis.ticks}{Logical, if axis ticks from the \code{axis3d()} function should be plotted. The default is \code{axis.ticks = TRUE}.}

\item{nticks}{Number of ticks for \code{axis3d()}.
The function repeats the last numeric value in the vector to cover all axis.
The user can, therefore, adjust the number of ticks with one numeric value (e.g., \code{nticks = 6}) or up to three (e.g., \code{nticks = c(6,4,8)} corresponding to the for the x, y, and z axes respectively.
The default is \code{nticks = 4}.}

\item{width.rgl.x}{Width in the x direction for \code{par3d()}. The default is \code{width.rgl.x = 1040}.}

\item{width.rgl.y}{Width in the y direction for \code{par3d()}. The default is \code{width.rgl.y = 1040}.}

\item{view}{Vector with polar coordinates and zoom factor for the \code{view3d} function. The default is \code{view = c(15,20, 1)}.}

\item{show.plane}{Logical, if xz-plane should be visible in the graphical device. The default is \code{show.plane = TRUE}.}

\item{plane.col}{Color of the plane, the default is \code{plane.col = "grey80"}.}

\item{background}{Set background color for the graphical device, the default is \code{background = "white"}.}

\item{type}{Type of vector arrow for items, the default is \code{type = "rotation"}. See \link[rgl:arrow3d]{rgl::arrow3d} for more options regarding arrow types.}

\item{col}{Vector of colors representing difficulty levels for item response functions used in \code{arrow3d()}. The default is \code{col = c("black", "grey20", "grey40", "grey60", "grey80")}.}

\item{arrow.width}{Width of vector arrows for \code{arrow3d()}. The default is \code{arrow.width = 0.6}.}

\item{n}{Number of barbs for the vector arrows from \code{arrow3d()}. The default is \code{n = 20}.}

\item{theta}{Opening angle of barbs for vector arrows from \code{arrow3d()}. The default is \code{theta = 0.2}.}

\item{barblen}{The length of the barbs for vector arrows from \code{arrow3d()}. The default is \code{barblen = 0.03}.}

\item{c.scalars}{Set of scalars for adjusting construct arrow length proportionally.
The first numeric adjusts the length in the positive direction, and the second numeric the length in the negative direction. The default is \code{c.scalars = c(1,1)}.}

\item{c.type}{Type of vector arrow for constructs. See \link[rgl:arrow3d]{rgl::arrow3d} for more options regarding arrow types. The default is \code{c.type = "rotation"}.}

\item{c.col}{Color of construct vector arrows from \code{arrow3d()}, the default is \code{c.col = "black"}.}

\item{c.arrow.width}{Width of construct vector arrows for \code{arrow3d()}. The default is \code{c.arrow.width = 0.6}.}

\item{c.n}{Number of barbs for the construct vector arrows from the \code{arrow3d()} function. The default is \code{c.n = 20}.}

\item{c.theta}{Opening angle of barbs for construct vector arrows from \code{arrow3d()}. The default is \code{c.theta = 0.2}.}

\item{c.barblen}{The length of the barbs for construct vector arrows from \code{arrow3d()}. The default is \code{c.barblen = 0.03}.}

\item{profiles}{Data frame with coordinates for spheres representing a subset of respondent scores. The default is \code{profiles = NULL}.}

\item{levels}{Optional. A column with values indicating levels for sphere colors from the \code{sphere.col} vector. The default is \code{levels = NULL}.}

\item{sphere.col}{Color vector for \code{spheres3d()}. The default is \code{sphere.col = c("black", "grey20", "grey40", "grey60", "grey80")}.}

\item{spheres.r}{Radius of the spheres for \code{spheres3d()}. The default is \code{spheres.r = 0.05}.}

\item{ci}{Logical, if spheres should include an ellipsoid outlining a confidence region returned from the \code{ellipse3d()} function. The default is \code{ci = FALSE}.}

\item{ci.level}{Level of confidence for \code{ellipse3d()}, the default is \code{ci.level = 0.95}.}

\item{ellipse.col}{Color of the ellipse from \code{ellipse3d()}. The default is \code{ellipse.col = "grey80"}.}

\item{ellipse.alpha}{Opacity for the confidence region from \code{ellipse3d()}. The default is \code{ellipse.alpha = 0.20}.}

\item{...}{Additional arguments passed to RGL or methods.}
}
\value{
A RGL graphical device.
}
\description{
For graphing of objects of class \code{D3mirt} from the \code{\link[=D3mirt]{D3mirt()}} function using the rgl 3D visualization device system (Adler & Murdoch, 2022).
}
\details{
The plotting function allows plotting of all items, a selection of items as well as plotting a single item.
Length of the vector arrows can be set to one unit length across all item vector arrows by setting \code{scale = TRUE}.
This removes the visualization of the MDISC parameter.
Note that when scaling items with \code{scale = TRUE}, the \code{plot()} function does not change the length of the model axes.
This often means that the axes of the model may need to be adjusted, which can be achieved proportionally with \code{axis.scalar} or manually with \code{axis.length}.

The user has the option of adding constructs to the graphical output with \code{constructs = TRUE} (see the documentation for \link{D3mirt} or the package vignette regarding constructs).
Other options include plotting one level of difficulty at a time with the \code{diff.level} argument if polytomous items are used in the model.
Item names are displayed by default, but the user has the option of adding new item labels for the items with \code{item.lab}, as well as labeling constructs with \code{construct.lab}.

Regarding the interpretation of results, the angle of the vector arrows indicates what traits, located along the orthogonal axes, an item can be said to describe (Reckase, 2009, 1985, Reckase & McKinley, 1991).
For instance, an item located at 0 degrees seen from the x-axis, and 90 degrees as seen from the y and z-axis, only describes trait x.
Such an item is unidimensional since its direction vector lies parallel and on the x-axis.
In contrast, an item located at 45 degrees between all three axes in a three-dimensional model describes all three traits in the model equally well.
Such an item is within-multidimensional with respect to all three latent traits used in the analysis because its direction vector points in a neutral direction in the model.

When plotting the \code{D3mirt} model with \code{plot()}, it is possible to visually observe statistical violations in the graphical output returned.
For instance, shorter vector arrows indicate weaker discrimination and, therefore, higher amounts of statistical violations.
Moreover, if a polytomous item struggles or even fails to describe any of the latent variables in the model, it can often lead to an extreme stretch of the MDIFF range.
This is comparable to trace lines turning horizontal in a unidimensional item response theory model.

The plot function can also display respondent scores in the three-dimensional model space, represented as spheres whose coordinates are derived from the respondent's factor scores.
This allows for a profile analysis in which respondent rows are separated or selected conditioned on some external criteria.
To do this, the user must first extract respondent factor scores with \link[mirt:fscores]{mirt::fscores} (Chalmers, 2012) and then use some selection process to separate or subset respondent rows.
The resulting data frame is used in the \code{profiles} argument.
If desired, a confidence interval can be added to the spheres by setting \code{ci = TRUE}.
A general advice is to hide vector arrows with \code{hide = TRUE} when analyzing respondent profiles to avoid visual cluttering.
For more on profile analysis (e.g., preparation and examples), see package vignette.

The returned RGL device can, for example, be exported to the R console and saved as an interactive HTML file or as a still shoot (see examples below).
In the the latter case, the model perspective in the still shoot can be manually adjusted by changing the \code{view} argument for the function.
}
\examples{
\donttest{
# Load data
data("anes0809offwaves")
x <- anes0809offwaves
x <- x[, 3:22] # Remove columns for age and gender

# Call D3mirt() with constructs assigned to con
con <- list(c(1,2,3,4,5,6,7,8,9,10),
            c(11,12,13,14),
            c(15,17,18,19,20))
mod <- D3mirt(x, modid = c("W7Q3", "W7Q20"), con.items = con)

# Plot RGL device with constructs invisible
plot(mod)

# Plot RGL device with constructs visible and named
plot(mod, constructs = TRUE,
     construct.lab = c("Compassion", "Fairness", "Conformity"))

plot(mod, constructs = TRUE,
     items = c(15, 17, 18, 19, 20),
     construct.lab = c("Compassion", "Fairness", "Conformity"))

# Item W7Q16 has location 16 in the data set (gender and age excluded)
# Below, the item is plotted together with construct to aid the visual interpretation
plot(mod, constructs = TRUE,
     items = 16,
     construct.lab = c("Compassion", "Fairness", "Conformity"))

# Plot RGL device on item difficulty level 5
plot(mod, diff.level = 5)

# Plot RGL device with scaled items and constructs visible and named
plot(mod, scale = TRUE,
     constructs = TRUE,
     construct.lab = c("Compassion", "Fairness", "Conformity"))

# Profile Analysis
# Plot respondents trait scores with item vectors hidden and no constructs
plot(mod, hide = TRUE, ind.scores = TRUE)

# Plot respondents separated on gender
# Start by assigning the gender variable to a data frame
# In this example, this is done by sub-setting the gender column
x <- anes0809offwaves

# Call plot() and use the gender variable column in the levels argument
# Respondent data on gender is in column two, x[, 2]
# In the function call below, both items and constructs are hidden
# Score levels: 1 = Blue ("male") and 2 = Red ("female")
plot(mod, hide = TRUE, ind.scores = TRUE,
    levels = x[, 2],
    sphere.col = c("blue", "red"),
    x.lab = "Compassion",
    y.lab="Conformity",
    z.lab="Fairness")

# Add a 95\% CI to respondent factor scores on <= 30 y.o.
# Column bind trait scores with the age variable "W3Xage"
z <- data.frame(cbind(mod$fscores, x[, 1]))

# Subset data frame z conditioned on age <= 30
z1 <- subset(z, z[, 4] <= 30)

# Use rep() to create a color vector to color groups based on the nlevels() output
# z1 has 14 factor levels
colvec <- c(rep("red", 14))

# Call plot() with profile data on age with item vector arrows hidden
# Use the profiles argument for the data frame containing the subset to be plotted
plot(mod, hide = TRUE,
    profiles = z1,
    levels = z1[, 4],
    sphere.col = colvec,
    x.lab = "Compassion",
    y.lab="Conformity",
    z.lab="Fairness",
    ci = TRUE,
    ci.level = 0.95,
    ellipse.col = "orange")
}
\dontrun{
# Export an open RGL device to the console to be saved as HTML or image file
plot(mod, constructs = TRUE)
s <- rgl::scene3d()
rgl::rglwidget(s,
               width = 1040,
               height = 1040)

# Export a snapshoot of an open RGL device directly to file
plot(mod, constructs = TRUE)
rgl::rgl.snapshot('RGLdevice.png',
                    fmt = 'png')
}
}
\references{
Adler, D., & Murdoch, D. (2022). \emph{Rgl: 3d Visualization Using OpenGL} Computer software.

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.\cr
https://doi.org/10.18637/jss.v048.i06

Reckase, M. D. (2009). \emph{Multidimensional Item Response Theory}. Springer.

Reckase, M. D. (1985). The Difficulty of Test Items That Measure More Than One Ability. \emph{Applied Psychological Measurement, 9}(4), 401-412.\cr
https://doi.org/10.1177/014662168500900409

Reckase, M. D., & McKinley, R. L. (1991). The Discriminating Power of Items That Measure More Than One Dimension. \emph{Applied Psychological Measurement, 15}(4), 361-373.\cr
https://doi.org/10.1177/014662169101500407
}
\author{
Erik Forsberg
}
