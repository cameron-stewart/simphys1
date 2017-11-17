from numpy import *
from matplotlib.pyplot import *

# Create 100 equidistant points between 0 and 2*pi
xs = linspace(0, 2.*pi, 100)

# Create a title for the plot
title("How to make plots")
# Create labels for the axes
# Note that in the label, you can use TeX-strings!
xlabel("$x$")
ylabel("$f(x)$")
# Now create two plots of different functions in a single figure
plot(xs, sin(xs), 'o-', label="$f(x)=\sin(x)$")
plot(xs, sin(xs)**2, '*:', label="$f(x)=\sin^2(x)$")
# Now create the legend box
legend()

# Create a second plot in loglog scale
figure()
title("A second plot, this time in loglog scale")
loglog(xs, xs**2, 'o-', label="$f(x)=x^2$")
loglog(xs, xs**0.5, '*:', label="$f(x)=\sqrt{x}$")
legend()

# Create a third plot with four subplots
figure("A title for the window")
# Select a subplot
# '221' describes the number of subplots in the figure and which subplot to select
# '22' tells to make a grid of 2 by 2 subplots
# '1' tells to select the first plot
subplot(221, title="$\sin(x)$")
plot(xs, sin(xs))
subplot(222, title="$\cos(x)$")
plot(xs, cos(xs))
subplot(223, title="$\sin^2(x)$")
plot(xs, sin(xs)**2)
subplot(224, title="$\cos^2(x)$")
plot(xs, cos(xs)**2)

show()
