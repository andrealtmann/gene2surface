%this converts the output from R into a freesurfer file
addpath('/Applications/freesurfer/matlab/');

%gene expression file
subj='mean';
dist='geod';
%probe='A_23_P48325';
%genename='KCTD4';
probe='A_23_P300600';
genename='NEFH';
exprfile=['results/' subj '_' dist '_' probe '.txt'];

%%
%freesurfer file_path
freepath='/Users/andre/work/MRI_templates/freesurfer_mni/MNI152_2mm/surf/';
freefile='lh.thickness';

%%
[dummy, faces] = read_curv([freepath freefile]);
expr = csvread(exprfile);

%%
write_curv([freepath 'lh.' genename '_' dist], expr, faces);

