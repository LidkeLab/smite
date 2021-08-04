% SOFTWARE included here:
%
% FRCresolution_software
%    This software is distributed as accompanying software for the article
%    "Measuring Image Resolution in Optical Nanoscopy" by R.J.P.
%    Nieuewenhuizen, K.A. Lidke, M. Bates, D. Leyton Puig, D. Gru\"nwald,
%    S. Stallinga, B. Rieger, Nature Methods, 2013 doi:10.1038/nmeth.2448.
%    Download via: http://www.diplib.org/add-ons
%    NOTE: cfsmooth in qcorrection_{ims,los}.m is the MATLAB code for smooth
%          from the Curve Fitting Toolbox.  This is needed because DIPimage
%          also has a smooth function which will typically shadow MATLAB's
%          smooth.
%
%    Quantitative Imaging Group Faculty of Applied Sciences,
%    Delft University of Technology
%    Lorentzweg 1, 2628 CJ Delft, The Netherlands
%    contact: Bernd Rieger, b.rieger@tudelft.nl
%
% PlotSpread
%    plotSpread allows creating "beeswarm plots", i.e. point distributions
%    where jitter has been added to the data points to avoid overlap.  It
%    further allows specifying groups within the data to show the distribution
%    of the groups inside a distribution.  plotSpread is most suited to
%    visualizing distributions with small amounts of data points.  If the
%    points become too dense, it becomes difficult to appreciate the relative
%    importance of modes of a distribution, in which case "distributionPlot"
%    should be used.
%
%    plotSpread uses the excellent "distinguishable_colors" to choose default
%    colors for different categories.
%
%    In addition, the .zip file contains "myErrorbar" (a modification of the
%    built-in errorbar), "repeatEntries" for easy repetition of entries in a
%    list, and "isEven" to test whether a number is even.  For both
%    "repeatEntries" and "isEven" there are better alternatives on the File
%    Exchange.
%
%    Cite As
%    Jonas (2020). plot spread points (beeswarm plot)
%    (https://www.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot),
%    MATLAB Central File Exchange. Retrieved September 17, 2020.
%
% uipickfiles
%    This is a GUI application that allows multiple files or directories to be
%    selected and allows you to manage the list (remove files, reorder, etc.)
%    before returning.  It has basic filtering as well as regular expression
%    filtering and navigation of the file system is easy.  The output is
%    configurable (cell, struct or char arrays).  It is written entirely in M
%    and so is platform independent.
%
%    Cite As
%       Douglas Schwarz (2021). uipickfiles: uigetfile on steroids
%       (https://www.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles-uigetfile-on-steroids),
%       MATLAB Central File Exchange. Retrieved August 4, 2021.
% -----------------------------------------------------------------------------
% setupExternalSoftware
%    function for startup.m to setup required paths for smite external software
