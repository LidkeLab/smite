% Examples of how to call Clustering routines.

% If true, save plots produced into the ResultsDir defined below.
Saving = false;

% --- 2D ---

SMF = smi_core.SingleMoleculeFitting();
SMF.Data.ResultsDir = 'Results';

if Saving
   if ~exist(SMF.Data.ResultsDir, 'dir')
      mkdir(SMF.Data.ResultsDir);
   end
end

% Simulate a sparse collection of hextets.
SIM = smi_sim.SimSMLM();
kmer = 6;                   % Hextets
PixelSize = 100;            % nm in a pixel
radius_kTet = 25/PixelSize; % Radius of the k-tets (pixel)
SZ = 100;                   % linear size of image (pixel)
receptorDensity = 1;        % Receptors/um^2
SIM.Rho = receptorDensity / (1000 / PixelSize)^2;
SIM.StartState = 'Equib';% If Equib, particle starts randomly on or off
SIM.SZ = SZ;
SIM.simkTets(kmer, radius_kTet);
SMD_Data = SIM.genNoisySMD(SIM.SMD_Model);

SZnm = SZ * PixelSize;

Npts = 1000;
% Generate a random cover of points.
xy  = SZnm * rand(2*Npts, 2);
xyz = SZnm * rand(3*Npts, 3);

% Create arrays of (x, y) or (x, y, z) coordinates (nm).
XY1 = zeros(Npts, 2);
XY1(:, 1) = xy(1:2:end, 1);
XY1(:, 2) = xy(2:2:end, 2);

XY3 = zeros(Npts, 3);
XY3(:, 1) = xyz(1:3:end, 1);
XY3(:, 2) = xyz(2:3:end, 1);
XY3(:, 3) = xyz(3:3:end, 1);

% Make SMD structures (pixel).
SMD1.X = XY1(:, 1) ./ PixelSize;
SMD1.Y = XY1(:, 2) ./ PixelSize;

SMD3.X = XY3(:, 1) ./ PixelSize;
SMD3.Y = XY3(:, 2) ./ PixelSize;
SMD3.Z = XY3(:, 3) ./ PixelSize;

SMD2 = SMD_Data;
indx = SMD2.X < 0 | SMD2.X > SZ | SMD2.Y < 0 | SMD2.Y > SZ;
SMD2.X(indx) = [];
SMD2.Y(indx) = [];
SMD2.X_SE(indx) = [];
SMD2.Y_SE(indx) = [];

% Create an array of (x, y) coordinates (nm).
XY2 = zeros(numel(SMD2.X), 2);
XY2(:, 1) = SMD2.X .* PixelSize;
XY2(:, 2) = SMD2.Y .* PixelSize;

fprintf('SMD1 [xmin, xmax, ymin, ymax] = [%5.1f, %5.1f, %5.1f, %5.1f] px\n',...
       min(SMD1.X), max(SMD1.X), min(SMD1.Y), max(SMD1.Y)); 
fprintf('SMD2 [xmin, xmax, ymin, ymax] = [%5.1f, %5.1f, %5.1f, %5.1f] px\n',...
       min(SMD2.X), max(SMD2.X), min(SMD2.Y), max(SMD2.Y)); 

% ========== Set up the smi_cluster.Clustering class (c) =========

% Note that some of the examples below use (x, y) coordinate arrays and some
% use SMD structures to indicate that both types of input may be accepted.
c = smi_cluster.Clustering(SMF);
c.PixelSize = PixelSize;
c.Timing = false;
ROI = [0, SZnm, 0, SZnm];   % nm
A_ROI = (ROI(2) - ROI(1)) * (ROI(4) - ROI(3));   % ROI area (nm^2)

% Compare each ROI nearest neighbor (nn) distances to a random distribution of
% points with the same density.
h = c.nn_ROIrandom(SMD1, A_ROI, 'Random Cover');
if Saving
   saveas(h, fullfile(SMF.Data.ResultsDir, 'RND2_nn'), 'png');
   close
end
h = c.nn_ROIrandom(SMD2, A_ROI, 'Hextets');
if Saving
   saveas(h, fullfile(SMF.Data.ResultsDir, 'HEX2_nn'), 'png');
   close
end

% ========== 2D ==========

fprintf('\n2D examples:\n\n');

% ---------- Random ----------

fprintf('Random:\n\n');
E = 200;      % minimum separation of points between clusters AND
              % maximum separation of points within  clusters (nm)
minPts = 3;   % minimum number of points required to form a cluster
% Ratio of local density / overall density for a point's Voronoi
% region to be considered sufficiently dense for clustering purposes.
c.Alpha = 1.2;

% ---[ XY ]---

% Call with (x, y) coordinate array (nm).  XY1 is a random cover.
algorithm = 'DBSCAN';
[nC, C, centers, ptsI] = c.cluster(algorithm, XY1, E, minPts);
fprintf('%s (E = %g, minPts = %d) number of clusters = %d\n', ...
        algorithm, E, minPts, nC);

% Compute statistics from the results of clustering.
results = c.clusterStats(XY1, C, centers);
fprintf('percentage clustered = %.1f\n', ...
        100 * results.n_clustered/results.n_points);

% Plot the clusters computed.
h = c.plotClusters(XY1, C, centers, ptsI, algorithm);
if Saving
   saveas(h, fullfile(SMF.Data.ResultsDir, ...
                sprintf('RND2_%s_E=%d,N=%d', algorithm, E, minPts)), 'png');
else
   figure(h)
end

% ---[ SMD ]---

% Call with SMD structure (pixel).  SMD1 is a random cover.
for algorithm_range = {'DBSCAN', 'Hierarchal', 'Voronoi'}
   algorithm = algorithm_range{1};
   if strcmp(algorithm, 'Voronoi')
      % Note that Voronoi uses Alpha (and optionally E) for its criteria.
      [nC, C, centers, ptsI] = c.cluster(algorithm, SMD1, [], minPts);
   else
      [nC, C, centers, ptsI] = c.cluster(algorithm, SMD1, E, minPts);
   end
   fprintf('%s (E = %g, minPts = %d) number of clusters = %d\n', ...
           algorithm, E, minPts, nC);

   % Compute statistics from the results of clustering.
   results = c.clusterStats(SMD1, C, centers);
   fprintf('percentage clustered = %.1f\n', ...
           100 * results.n_clustered/results.n_points);

   % Plot the clusters computed.
   h = c.plotClusters(SMD1, C, centers, ptsI, algorithm);
   if Saving
      saveas(h, fullfile(SMF.Data.ResultsDir, ...
                   sprintf('RND2_%s_E=%d,N=%d', algorithm, E, minPts)), 'png');
   else
      figure(h)
   end
end % for algorithm_range
fprintf('\n');

% ---------- Hextets ----------

fprintf('Hextets:\n\n');
E = 50;       % minimum separation of points between clusters AND
              % maximum separation of points within  clusters (nm)
minPts = 3;   % minimum number of points required to form a cluster
% Ratio of local density / overall density for a point's Voronoi
% region to be considered sufficiently dense for clustering purposes.
c.Alpha = 1.2;

% ---[ XY ]---

% Call with (x, y) coordinate array (nm).
algorithm = 'DBSCAN';
[nC, C, centers, ptsI] = c.cluster(algorithm, XY2, E, minPts);
fprintf('%s (E = %g, minPts = %d) number of clusters = %d\n', ...
        algorithm, E, minPts, nC);

% Compute statistics from the results of clustering.
results = c.clusterStats(XY2, C, centers);
fprintf('percentage clustered = %.1f\n', ...
        100 * results.n_clustered/results.n_points);

% Plot the clusters computed.
h = c.plotClusters(XY2, C, centers, ptsI, algorithm);
if Saving
   saveas(h, fullfile(SMF.Data.ResultsDir, ...
                sprintf('HEX2_%s_E=%d,N=%d', algorithm, E, minPts)), 'png');
else
   figure(h)
end

% Plot the clusters computed with the SE of the localizations
h = c.plotClustersSE(SMD2, C, centers, ptsI, algorithm);
if Saving
   saveas(h, fullfile(SMF.Data.ResultsDir, ...
                sprintf('RND2_%s_E=%d,N=%d_SE', algorithm, E, minPts)), 'png');
else
   figure(h)
end

% ---[ SMD ]---

% Call with SMD structure (pixel).  SMD2 is a scattering of hextets.
for algorithm_range = {'DBSCAN', 'Hierarchal', 'Voronoi', 'H-SET'}
   algorithm = algorithm_range{1};
   if strcmp(algorithm, 'Voronoi')
      % Note that Voronoi uses Alpha (and optionally E) for its criteria.
      [nC, C, centers, ptsI] = c.cluster(algorithm, SMD2, [], minPts);
   else
      % H-SET needs to use an SMD structure with X_SE, Y_SE defined.
      [nC, C, centers, ptsI] = c.cluster(algorithm, SMD2, E, minPts);
   end
   fprintf('%s (E = %g, minPts = %d) number of clusters = %d\n', ...
           algorithm, E, minPts, nC);

   % Compute statistics from the results of clustering.
   results = c.clusterStats(SMD2, C, centers);
   fprintf('percentage clustered = %.1f\n', ...
           100 * results.n_clustered/results.n_points);

   % Plot the clusters computed.
   h = c.plotClusters(SMD2, C, centers, ptsI, algorithm);
   if Saving
      saveas(h, fullfile(SMF.Data.ResultsDir, ...
                   sprintf('HEX2_%s_E=%d,N=%d', algorithm, E, minPts)), 'png');
   else
      figure(h)
   end
end % for algorithm_range

% ========== 3D ==========

fprintf('\n3D examples:\n\n');

% ---------- Random ----------

fprintf('Random:\n\n');
E = 500;      % minimum separation of points between clusters AND
              % maximum separation of points within  clusters (nm)
minPts = 3;   % minimum number of points required to form a cluster
% Ratio of local density / overall density for a point's Voronoi
% region to be considered sufficiently dense for clustering purposes.
c.Alpha = 1.2;

% ---[ SMD ]---

% Call with SMD structure (pixel).  SMD3 is a random cover.
for algorithm_range = {'DBSCAN', 'Hierarchal', 'Voronoi'}
   algorithm = algorithm_range{1};
   if strcmp(algorithm, 'Voronoi')
      % Note that Voronoi uses Alpha (and optionally E) for its criteria.
      [nC, C, centers, ptsI] = c.cluster(algorithm, SMD3, [], minPts);
   else
      [nC, C, centers, ptsI] = c.cluster(algorithm, SMD3, E, minPts);
   end
   fprintf('%s (E = %g, minPts = %d) number of clusters = %d\n', ...
           algorithm, E, minPts, nC);

   % Compute statistics from the results of clustering.
   results = c.clusterStats(SMD3, C, centers);
   fprintf('percentage clustered = %.1f\n', ...
           100 * results.n_clustered/results.n_points);

   % Plot the clusters computed.
   h = c.plotClusters3(SMD3, C, centers, ptsI, algorithm);
   view(-37.5, 30);
   if Saving
      saveas(h, fullfile(SMF.Data.ResultsDir, ...
                   sprintf('RND3_%s_E=%d,N=%d', algorithm, E, minPts)), 'png');
   else
      figure(h)
   end
end % for algorithm_range
fprintf('\n');

fprintf('Done.\n');
