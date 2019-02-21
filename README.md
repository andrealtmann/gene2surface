# gene2surface
Visualizing gene expression on the FreeSurfer surface

Parts of the pipeline require a working installation of
 - R
 - MATLAB
 - FreeSurfer

The visualization is based on kernel density smoothing (KDS) using either eucledian or geodesic distance. For geodesic distance FreeSurfer's implementation was used to compute shortest paths from source vertices (training data) to all vertices on the surface. There are larger (preprocessed) files needed to execute the pipeline:
 - RData for cortical samples from the Allen Data
 - pre-computed geodesic distances

In it's default version the pipleine uses updated MNI coordinates obtained from:
https://github.com/chrisfilo/alleninf

In order to run the pipeline simply run the following commands:
> R --no-save -q < src/run_pipeline.R
(alterantively open R and source('src/run_pipeline.R')

Run matlab to produce the new curv file:
> matlab -nodisplay -nojvm -nosplash -nodesktop -r "run('./src/array2curv.m'); exit;" 

Run freeview to display
> freeview -f ./Data/lh.smoothwm:overlay=./results/lh.KCTD4_geod --viewport 3d

