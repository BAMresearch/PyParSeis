# PyParSeis
Python code for Parallel Seismic test valuation (determination of foundation depth)

# Purpose
The Parallel seimic (PS) method ist used to determine the length of foundation piles since decades. A hammer ist used to provide an impact on top of the foundation, generating a stress wave propagating downwards in side the pile. A part of the energy leaks into the soil. Hydrophones or geophones are used to detect this seismic wave in a parallel borehole close to the pile.
From the recording, the first arrival time is determined and plotted against deetctor depth. In the past, just the depth of the change of slope was used to provide an estimate of pile length (below the pile waves are propagating at a mostly much slower speed in the soil). Unfortunately, this approach leads to an overestimation of pile length, potentially leading to unsafe load capacity extimates. Improved methods have been developed in the early 2000s, e.g. to incorporate the distance between pile and borehole.
In my PhD thesis , I have deloped a method using curve fitting, using a model, which in addition to wave speed and pile length incorporates pile diameter, distance pile-borehole and pile inclination (actually relative inclination between pile and borehole). See references below. The coding was done first in C, then in Matlab. When my group at BAM slowly moved from Matlab to Python to avoid license cost, I thought it to be worthy to do an exercise on my own and started to recode my results. In 2021, I have provided a version to an engineering firm, who used it in practice and found it useful. Meanwhile, I have contrbuted to an ASTM standard on PS, which mentions the method.  

Feel free to use, comment and improve, as long it goes along with the GPL 3.0 license, but on your own risk. No warranties on whatever.

# Content (so far)
- test data (synthetic)
- python scripts for producing syntehic datasets (first arrival times)(
- python scripts for evaluating fiste arrival data (without inclination, yet)

# References
My PhD-Thesis: Niederleithinger, Ernst, 2010: "https://publishup.uni-potsdam.de/opus4-ubp/frontdoor/deliver/index/docId/4708/file/niederleithinger_diss.pdf". Dissertation, Universität Potsdam, https://publishup.uni-potsdam.de/opus4-ubp/frontdoor/deliver/index/docId/4708/file/niederleithinger_diss.pdf
Main paper: Niederleithinger, Ernst, 2012. „Improvement and extension of the parallel seismic method for foundation depth measurement“. SOILS AND FOUNDATIONS 52 (6): 1093–101. https://doi.org/10.1016/j.sandf.2012.11.023.
ASTM D8381/D8381M-21: Standard Test Methods for Measuring the Depth of Deep Foundations by Parallel Seismic Logging

Ernst Niederleithinger, 2025, ernst.niederleithinger@bam.de
 
