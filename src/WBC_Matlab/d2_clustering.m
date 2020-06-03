
% %%%%%%%%%%%%%%%%%%%%%%%%% %
% Timing Data               %
% %%%%%%%%%%%%%%%%%%%%%%%%% %

% To run this experiment, uncomment all of these:

numOfSamples = 4;
s_modalities = 1;
d_modalities = [2];
num_clusters = 2;

% and one of these:

filename='../../data/timing/data25.d2';
%filename='../../data/timing/data50.d2';
%filename='../../data/timing/data75.d2';
%filename='../../data/timing/data100.d2';

% Note that the number is the number of points per dataset
% There are four datasets per experiment

% %%%%%%%%%%%%%%%%%%%%%%%%% %
% Cubic Structures Data     %
% %%%%%%%%%%%%%%%%%%%%%%%%% %

% To run this experiment, uncomment all of these:

%numOfSamples = 6;
%s_modalities = 1;
%d_modalities = [3];
%num_clusters = 2;

% and one of these:

% No transformation
%filename='../../data/cubic_structures/transformed_data/data0.d2';
% Reflection
%filename='../../data/cubic_structures/transformed_data/data1.d2';
% Rotation
%filename='../../data/cubic_structures/transformed_data/data2.d2';
% Translation
%filename='../../data/cubic_structures/transformed_data/data3.d2';

% %%%%%%%%%%%%%%%%%%%%%%%%% %
% Carbon Allotropes Data    %
% %%%%%%%%%%%%%%%%%%%%%%%%% %

% To run this experiment, uncomment all of these:

%numOfSamples = 6;
%s_modalities = 1;
%d_modalities = [3];
%num_clusters = 2;

% and one of these:

% No transformation
%filename='../../data/carbon_allotropes/transformed_data/data0.d2';
% Reflection
%filename='../../data/carbon_allotropes/transformed_data/data1.d2';
% Rotation
%filename='../../data/carbon_allotropes/transformed_data/data2.d2';
% Translation
%filename='../../data/carbon_allotropes/transformed_data/data3.d2';

% %%%%%%%%%%%%%%%%%%%%%%%%% %
% The rest runs the algo... %
% %%%%%%%%%%%%%%%%%%%%%%%%% %

db = loaddata(numOfSamples, s_modalities, d_modalities, filename);
global statusIterRec;

tic
max_stride = max(cellfun(@(x) max(x.stride), db));
kantorovich_prepare(max_stride, max_stride);
[clusters, labels] = d2clusters(db, num_clusters);
toc

disp("labels")
disp(labels)