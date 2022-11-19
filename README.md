# spline-interpolation-NumericalAnalysis
Solution of Problems 7 and 8, Section 6.4 from the book Numerical Analysis: Mathematics of Scientific Computation

Section 6.4, Problem 7:
Draw a script letter. Then reproduce it with the aid of cubic splines and a plotter. Proceed as follows: Select a modest number of points on the curve, say n = 11. Label these t = 1, 2 , . . . , n. For each point, obtain the corresponding x- and y-coordinates. Then fit x = S x (t) and y = S y (t), using cubic spline interpolating functions S x and S y . This will produce a parametric representation of the original curve. Compute a large number of values of S x (t) and S y (t) to give to the plotter

<Section 6.4, Problem 8:>
Interpret the results of the following numerical experiment and draw some conclusions.
-a. Define p to be the polynomial of degree 20 that interpolates the function f(x) = (1 + 6x 2 )~ l at 21 equally spaced nodes in the interval [—1, 1]. Include the endpoints as nodes. Print a table of f(x), p(x), and f(x) — p(x) at 41 equally spaced points on the interval.

-b. Repeat the experiment using the Chebyshev nodes given by xi = cos[(i - 1)JT/20] (1 < i < 21) c. With 21 equally spaced knots, repeat the experiment using a cubic interpolating spline.

-c. With 21 equally spaced knots, repeat the experiment using a cubic interpolating spline.
