function wavecnv_segmenter( big_h5_file, level, threshmethod, hdf_datafield, output_file_string, slice_toggle, start_string, end_string,scaling_string, p_value_string, chr_string,bpoint_only_string )

%function wavecnv_hdf5_complete( big_h5_file, level, threshmethod, hdf_datafield, output_file_string, slice_toggle, start_string, end_string, scaling_string, p_value_string, chr_string, bpoint_only_string )

%%%%%%%%%%%%%%%thr_method comment%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% thr_method = 'sqtwolog';
%'rigrsure' 'heursure' 'sqtwolog' 'minimaxi' are possibilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bpoint_only = str2num( bpoint_only_string );

if bpoint_only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Grab, precondition, and compute SWT on data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data = hdfinput( big_h5_file, hdf_datafield, slice_toggle, start_string, end_string, scaling_string ); %grab 2-column data from hdf5 file
data_fid = fopen(big_h5_file);
data = cell2mat(textscan(data_fid, '%f %f', 'delimiter', ','));
fclose(data_fid);
%data = mat2cell(data);
%celldisp(data);

scaling_toggle2=str2num(scaling_string);
if scaling_toggle2
raw_coverage_factor = sum( data(:,2) )./length( data(:,2) );
diploid_coverage_factor = 2./raw_coverage_factor;
data(:,2) = data(:,2).*diploid_coverage_factor;
end 

level = str2num(level);
pow = 2^level;
ext_length = ceil(length(data(:,2))/pow)*pow - length(data(:,2)) 
testRDdata=wextend( '1', 'per', data(:,2), ext_length -1  ,'r');

if ~(rem(length(testRDdata),2^level)) 
testRDdata=wextend( '1', 'per', data(:,2), ext_length-1  ,'r');
else
    testRDdata=wextend( '1', 'per', data(:,2), ext_length  ,'r');
end

%implement level by level thresholding on swt

SWA = zeros( level, length(testRDdata) );
SWD = zeros(level, length(testRDdata) );
tic;[SWA,SWD] = swt( testRDdata, level, 'haar' );toc; %get whole SWT structure; subsequently reconstruct SWC = [SWD; SWA(end,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BREAKPOINTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zerocrossings_ind = find( diff( [0;sign( (2^(-1))*(SWA(3, 1:end-ext_length)').^2 - (2^-(level))*(SWA(level, 1:end-ext_length)').^2 ) ] ) ~= 0);  
														%conjecture: zero crossings of "sufficiently distant" 
														%approximation coefficient bands represent physically meaningful breakpoints

%the issue now is how to pair these putative breakpoints, but first we filter them out using two concepts
%1) remove all points separated by fewer than 10 indices and save them to a separate file
%2) for what clusters remain that are adjacent, perform a two-sided kstest to determine if they are statistically equivalent. If they are, merge them and repeat until all are reported as distinct.

new_set = unique( [ 1 ; zerocrossings_ind ; length(SWA(1,1:end-ext_length)) ] ); %insert beginning and end points

zerocrossings_cluster = [new_set(1:end-1) new_set(2:end)];           %classic-pair
zerocrossings_cluster(2:end, 1) = zerocrossings_cluster(2:end, 1)+1; %double-pair

%note that the kstest p-value returned by matlab is only really asymptotically accurate when n1*n2 / (n1+n2) ~> 4, i.e. n1 ~ n2 ~ 10, say.

%grab those elements of the matrix that have a high p-value (some cutoff) 

%%%%%%%%%start merger%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

p_val = str2num( p_value_string );

datacluster3 = shortmaker( zerocrossings_cluster, [data(:,1) data(:,2)] );

while sum( datacluster3(:,7) > p_val ) > 0 %given any p-values above p_value_string, repeat this procedure of cluster agglomerization

merge_regions = consecutiveindex( find( datacluster3(:,7) > p_val ) ); %find and print out consecutive indices of high-p clusters
blist=find( cellfun(@length, merge_regions)==1 ); %tag singleton high-p regions and add rightward neighbour

for i = 1:length( blist )
merge_regions{1, blist(i)} = [ merge_regions{1,blist(i)} merge_regions{1,blist(i)}+1];
end

%merge high-p regions via vectorized cellfun
merged_datacluster = [ datacluster3( cellfun(@min, merge_regions)' ,1 ) datacluster3( cellfun(@max, merge_regions)', 2 ) ]; %what are the extrema of the high-p regions
                                                                                                                            %which now have a minimum width of 2
%what are the other regions (without the padding)                                                                                                                            
rest=[ datacluster3( setdiff( [1:length(datacluster3)],[merge_regions{:,:}]' )', 1) datacluster3( setdiff( [1:length(datacluster3)],[merge_regions{:,:}]' )', 2) ];
                                                                                                                                                                    %complement
%form new short table by merging them
total_new = [merged_datacluster; rest];
[ss,ssidx] = sort( total_new(:,1), 'ascend' );
%sort and assemble
total_new_sorted = total_new(ssidx, :);
%recompute p-values of new clusters
datacluster3 = shortmaker(  total_new_sorted, data );
end
%%%%%end_merger%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%END BREAKPOINTS

%write breakpoints%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chr = str2num(chr_string);

datacluster3 = [repmat(chr, length( total_new_sorted), 1) datacluster3];

%breakpoint table
fid = fopen( [output_file_string '_short' ], 'w' );
fprintf( fid, '%lu\t %lu\t %lu\t %lu\t %lu\t %.4f\t %.4f\t %.4f\n', datacluster3' );
fclose(fid);

% $$$ total_hdf_string = big_h5_file;
% $$$ 
% $$$ master_short_table = mastershortmaker_final( datacluster3, total_hdf_string );
% $$$ 
% $$$ fid = fopen( [output_file_string '_short_annotated' ], 'w' );
% $$$ fprintf( fid, [repmat( '%u\t', 1, 5) '%.4f\t' '%.4f\t' '%e\t' '%.4f\t' '%e\t' '%.4f\t' repmat('%u\t', 1, 3) '%.4f\t' repmat( '%.4f\t', 1, 6) ...
% $$$                                                                        repmat('%u\t',1,3) '%.4f\t' '%.4f\n'], master_short_table' );
% $$$ fclose(fid);



else   %%%%%%%%%%%%end brekpoint only condition


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BREAKPOINT PLUS DENOISING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = hdfinput( big_h5_file, hdf_datafield, slice_toggle, start_string, end_string, scaling_string ); %grab data from hdf5 file

level = str2num(level);

%periodic right-padding--later removed-- to ensure divisibility by 2^level for filter banks

%extend rightward boundary to avoid artificial discontinuities and also simulataneously satisfy the 2^level divisibility requirement

pow = 2^level;
ext_length = ceil(length(data(:,2))/pow)*pow - length(data(:,2))
testRDdata=wextend( '1', 'per', data(:,2), ext_length -1  ,'r');

if ~(rem(length(testRDdata),2^level))
testRDdata=wextend( '1', 'per', data(:,2), ext_length-1  ,'r');
else
    testRDdata=wextend( '1', 'per', data(:,2), ext_length  ,'r');
end

%implement level by level thresholding on swt

SWA = zeros( level, length(testRDdata) );
SWD = zeros(level, length(testRDdata) );
tic;[SWA,SWD] = swt( testRDdata, level, 'haar' );toc; %get whole SWT structure; subsequently reconstruct SWC = [SWD; SWA(end,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BREAKPOINTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zerocrossings_ind = find( diff( [0;sign( (2^(-1))*(SWA(3, 1:end-ext_length)').^2 - (2^-(level))*(SWA(level, 1:end-ext_length)').^2 ) ] ) ~= 0);
                                                                                                                %conjecture: zero crossings of "sufficiently distant"
                                                                                                                %approximation coefficient bands represent physically meaningful breakpoints

%the issue now is how to pair these putative breakpoints, but first we filter them out using two concepts
%1) remove all points separated by fewer than 10 indices and save them to a separate file
%2) for what clusters remain that are adjacent, perform a two-sided kstest to determine if they are statistically equivalent. If they are, merge them and repeat until all are reported as distinct.

new_set = unique( [ 1 ; zerocrossings_ind ; length(SWA(1,1:end-ext_length)) ] ); %insert beginning and end points

zerocrossings_cluster = [new_set(1:end-1) new_set(2:end)];
zerocrossings_cluster(2:end, 1) = zerocrossings_cluster(2:end, 1)+1;

%note that the kstest p-value returned by matlab is only really asymptotically accurate when n1*n2 / (n1+n2) ~> 4, i.e. n1 ~ n2 ~ 10, say.

%grab those elements of the matrix that have a high p-value (some cutoff)

%%%%%%%%%start merger%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

p_val = str2num( p_value_string );

datacluster3 = shortmaker( zerocrossings_cluster, [data(:,1) data(:,2)] );

while sum( datacluster3(:,7) > p_val ) > 0 %given any p-values above p_value_string, repeat this procedure of cluster agglomerization

merge_regions = consecutiveindex( find( datacluster3(:,7) > p_val ) ); %find and print out consecutive indices of high-p clusters
blist=find( cellfun(@length, merge_regions)==1 ); %tag singleton high-p regions and add rightward neighbour

for i = 1:length( blist )
merge_regions{1, blist(i)} = [ merge_regions{1,blist(i)} merge_regions{1,blist(i)}+1];
end

%merge high-p regions via vectorized cellfun
merged_datacluster = [ datacluster3( cellfun(@min, merge_regions)' ,1 ) datacluster3( cellfun(@max, merge_regions)', 2 ) ]; %what are the extrema of the high-p regions
                                                                                                                            %which now have a minimum width of 2
%what are the other regions (without the padding)
rest=[ datacluster3( setdiff( [1:length(datacluster3)],[merge_regions{:,:}]' )', 1) datacluster3( setdiff( [1:length(datacluster3)],[merge_regions{:,:}]' )', 2) ];
                                                                                                                                                                    %complement
%form new short table by merging them
total_new = [merged_datacluster; rest];
[ss,ssidx] = sort( total_new(:,1), 'ascend' );
%sort and assemble
total_new_sorted = total_new(ssidx, :);
%recompute p-values of new clusters
datacluster3 = shortmaker(  total_new_sorted, data );
end

%%%%%%%%%%%%%%%%%%%%%END BREAKPOINTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HAAR DENOISING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SWC = zeros( level + 1, length(testRDdata) );

SWC = [SWD;SWA(end,:)]; %construct SWC

cfactor=repmat(.6745, level, 1);

% thr = wthrmngr('sw1ddenoLVL','penalhi', SWC,2);

for i = 1:size( SWC, 1) - 1
        thr(i) = mad( abs(SWC(i,:)),0 )./cfactor(i);   %threshold detail coefficients according to mad filter with Gaussian golden number .6745; kept as factor to recall that it $
end

 [thr ; max( abs(SWC(1:end-1,:)), [], 2)']    %thresholded value and maximum value of detail coefficients, as function of level

dswd = zeros( level, length( testRDdata) ); %preallocate thresholded coefficients

for i = 1:length(thr)
dswd(i,:) = wthresh(SWC(i,:), 'h', thr(i)*sqrt(2*log(length(SWC))));  %hard thresholding
end

den_swt_coeffs=[dswd; SWC(end,:)];  %form thresholded detail coefficients and last level of approx. coefficients to perform inverse transform

denoised_swt = zeros( length(SWC),1 ); %preallocate denoised signal in lieu of inverse transform

denoised_swt = iswt(den_swt_coeffs,'haar')';   %perform inverse transformation

[ 'detail coefficients thresholded with ' threshmethod]

denoised_swt_true = [data denoised_swt(1:end-ext_length,:)]; %remove padding

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%denoise with no detail coefficients (i.e. threshold everything )%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

denoisedA_swt = iswt( [zeros( level, length( testRDdata) ); SWC(end,:)  ],'haar')';
denoisedA_swt_true = [denoisedA_swt(1:end-ext_length,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SWA_bmat = [ SWA(3,1:end-ext_length)' SWC(end,1:end-ext_length)' ]; %record the L3 treetops' coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%grab associated qualities%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qdata = hdfinput( big_h5_file, 'GposMedianQ', slice_toggle, start_string, end_string, '0' ); %same dimension as read depth column, so unambiguous

%kmerdata = hdfinput( big_h5_file, GposMedianQ, slice_toggle, start_string, end_string, 0 ); idx is gpos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chr = str2num(chr_string);

master_output = [repmat( chr, length( denoised_swt_true), 1) denoised_swt_true denoisedA_swt_true Qdata(:,2) SWA_bmat]; %form chr | gpos | raw RD | SWD denoised | SWA-only denoised | med Q | SWA_1 | SWA_end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%write general output

fid = fopen(  [output_file_string '_denoised' ], 'w' );
fprintf( fid, '%lu\t %lu\t %.1f\t %.4f\t %.4f\t %lu\t %.4f\t %.4f\n',  master_output' );  %note the transpose necessary for fprintf
fclose(fid);


datacluster3 = [repmat(chr, length( total_new_sorted), 1) datacluster3];

%breakpoint table
fid = fopen( [output_file_string '_short' ], 'w' );
fprintf( fid, '%lu\t %lu\t %lu\t %lu\t %lu\t %.4f\t %.4f\t %.4f\n', datacluster3' );
fclose(fid);

total_hdf_string = big_h5_file;

master_short_table = mastershortmaker_final( datacluster3, total_hdf_string );

fid = fopen( [output_file_string '_short_annotated' ], 'w' );
fprintf( fid, [repmat( '%u\t', 1, 5) '%.4f\t' '%.4f\t' '%e\t' '%.4f\t' '%e\t' '%.4f\t' repmat('%u\t', 1, 3) '%.4f\t' repmat( '%.4f\t', 1, 6) ...
                                                                       repmat('%u\t',1,3) '%.4f\t' '%.4f\n'], master_short_table' );
fclose(fid);



end   %end total condition (breakpoints plus denoised signal )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDEND%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%AUXILLIARY FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function appended_short_table = shortmaker( paired_breakpoints, signal )

for i = 1:length(paired_breakpoints)
paired_breakpoints(i,3) = signal(paired_breakpoints(i,1),1);
paired_breakpoints(i,4) = signal(paired_breakpoints(i,2),1);
paired_breakpoints(i,5) = mad( signal(paired_breakpoints(i,1):paired_breakpoints(i,2),2),1 );
paired_breakpoints(i,6) = median( signal(paired_breakpoints(i,1):paired_breakpoints(i,2),2) );
datacluster{i,1} = signal(paired_breakpoints(i,1):paired_breakpoints(i,2),2);
end

%create general cell object in which to store the data in the putative clusters and perform the kstest on the cell object; append

datacluster2 = { {datacluster{1:end-1,1}}' {datacluster{2:end, 1}}' };
a = cellfun(@(x) x(:,1), datacluster2{1,1}, 'UniformOutput', false );  %pull objective column from 1st cell
b = cellfun(@(x) x(:,1), datacluster2{1,2}, 'UniformOutput', false );  %second cell
a( cellfun(@length, a) > 50  ) = cellfun( @(x) x(end-50:end,1), a(  cellfun(@length, a) > 50 ), 'UniformOutput', false );
b( cellfun(@length, b) > 50  ) = cellfun( @(x) x(1:50,1), b(  cellfun(@length, b) > 50 ), 'UniformOutput', false );
d = mat2cell([a b],ones(size(a,1),1),2);
[h p k] = cellfun(@(x)kstest2(x{1,1},x{1,2},.00001),d);

appended_short_table = [ paired_breakpoints [p;0] [k;0] ];


%%%%
function [conseuctivelist] = consecutiveindex( input )

conseuctivelist= mat2cell(input',1,diff([0,find(diff(input') ~= 1),length(input')]));

%%%%%%%%%%%%%%%%%%%%%%%%%read hdf5 funciton%%%%%%%%%%%%%%%%%%%%%%%%

function data = hdfinput( big_h5_file, hdf_datafield, slice_toggle, start_string, end_string, scaling_toggle )

%%%%HDF5 manipulations
fileID = H5F.open( big_h5_file, 'H5F_ACC_RDONLY','H5P_DEFAULT');

datasetID = H5D.open( fileID, hdf_datafield );

% Get dataspace
dataspaceID = H5D.get_space( datasetID );
[rank_GposRd dims_GposRd] = H5S.get_simple_extent_dims(dataspaceID)

stride=[];
block =[];

slice_toggle2=logical( str2num(slice_toggle) );
scaling_toggle2=logical( str2num(scaling_toggle));
if slice_toggle2   %slice_toggle controls if we want a slice--in idx space -- of the data instead of the whole thing

start_index = str2num( start_string);
end_index = str2num( end_string);
dimsRequired = [  end_index dims_GposRd(2)];
start = [start_index 0];

else
dimsRequired = [  dims_GposRd(1) dims_GposRd(2)]
start = [1 0]

end
count=dimsRequired - start;
H5S.select_hyperslab(dataspaceID,  'H5S_SELECT_SET', start, stride, count, block);

% Define the memory dataspace.
memspaceID = H5S.create_simple(rank_GposRd, count, []);
%output region
data_extract = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT')';

%%END HDF5 manipulations

if scaling_toggle2
raw_coverage_factor = sum( data_extract(:,2) )./length( data_extract(:,2) );
diploid_coverage_factor = 2./raw_coverage_factor;
data(:,2) = data_extract(:,2).*diploid_coverage_factor;
disp(['base scaling factor is ' num2str( raw_coverage_factor ) ] );
data(:,1) = data_extract(:,1);
else
data = data_extract;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FINAL SHORT TABLE ANNOTATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [master_short_table] = mastershortmaker_final( short_table,total_hdf_string )

master_short_table = [short_table zeros(length(short_table), 19) ];

disp('short table read and master_short preallocation complete' )

%slice data to cells in one shot

rawsigregions = hdfreader2_vec_input_2( total_hdf_string,short_table(:,2), short_table(:, 3),...
                                       short_table(:,4), short_table(:, 5),...
                                       short_table(:,4), short_table(:,5),...
                                       floor(short_table(:,4)/200), ceil(short_table(:, 5)/200)  );


disp( 'total HDF read and cell creation complete' )

rawsigregions2 = { {rawsigregions{1:end-1,1}}' {rawsigregions{2:end, 1}}' };  %stagger data for kstest2 cluster comparison, minding boundaries

disp( 'staggered raw region cells formed' )

aaa=cell2mat(rawsigregions2{1,1}(1));
bbb=cell2mat(rawsigregions2{1,2}(2));

aaa(1,:)
bbb(1,:)

a = cellfun(@(x) x(:,2), rawsigregions2{1,1}, 'UniformOutput', false );  %pull objective column from 1st cell
b = cellfun(@(x) x(:,2), rawsigregions2{1,2}, 'UniformOutput', false );  %second cell

disp( 'form the obective staggered columns for KSTEST2:1' )

%slice boundary points in logical way

a( cellfun(@length, a) > 50  ) = cellfun( @(x) x(end-50:end,1), a(  cellfun(@length, a) > 50 ), 'UniformOutput', false );
b( cellfun(@length, b) > 50  ) = cellfun( @(x) x(1:50,1), b(  cellfun(@length, b) > 50 ), 'UniformOutput', false );


% Reorganize the data such that all the inputs for each call to KSTEST2 are in a separate cell
d = mat2cell([a b],ones(size(a,1),1),2);
disp( 'form the obective staggered columns for KSTEST2:2' )
% Use CELLFUN to apply a function on every element of the cell array
[h p k] = cellfun(@(x)kstest2(x{1,1},x{1,2}),d);   %perform vectorized kstest2

disp( 'KStest2 performed' )


Qmed = cellfun(@(x) median(x(:,2)), {rawsigregions{:,2} } )';
Qmin = cellfun(@(x) min(x(:,2)),  {rawsigregions{:,2} } )';
Qmax = cellfun(@(x) max(x(:,2)),  {rawsigregions{:,2} } )';
Qmad = cellfun(@(x) mad(x(:,2)),  {rawsigregions{:,2} } )';
alpha = 1 - (Qmed - Qmax)./(Qmin - Qmax);

disp( 'Qmed,Qmin,Qmax,alpha formed' )

kstest_mat = [p k;[0 0]];

max_kmer = cellfun(@(x) max(x(:,2)),  {rawsigregions{:,3} } )';
min_kmer = cellfun(@(x) min(x(:,2)), {rawsigregions{:,3} } )';
med_kmer = cellfun(@(x) median(x(:,2)),  {rawsigregions{:,3} } )';
mad_kmer = cellfun(@(x) mad(x(:,2)),  {rawsigregions{:,3} } )';

max_GC_ratio = cellfun(@(x) max(x(:,7)),  {rawsigregions{:,4} } )';
min_GC_ratio = cellfun(@(x) min(x(:,7)),  {rawsigregions{:,4} } )' ;
med_GC_ratio = cellfun(@(x) median(x(:,7)),  {rawsigregions{:,4} } )';
mad_GC_ratio = cellfun(@(x) mad(x(:,7)),  {rawsigregions{:,4} } )';

max_GC_ref_ratio = cellfun(@(x) max(x(:,9)), {rawsigregions{:,4} } )';
med_GC_ref_ratio = cellfun(@(x) median(x(:,9)),  {rawsigregions{:,4} } )';
master_short_table = [ short_table kstest_mat max_kmer min_kmer med_kmer mad_kmer max_GC_ratio min_GC_ratio med_GC_ratio...
                                                 mad_GC_ratio max_GC_ref_ratio med_GC_ref_ratio Qmax Qmin Qmed Qmad alpha ];


disp('master_short_table formed, writing started' )

%fid = fopen(  [ short_table_string 'annotated_KSSEG' ], 'w' );
% header_array = {'Idx_start' ' Idx_end' ' Gpos_start' ' Gpos_end' ' RawRD'' ScaledRD' ' KS_h' ' KS_p',' KS_stat',...
%                  ' kmer-21_max' ' kmer-21_med' ' kmer-21_min' ' kmer-21_mad' ' GC_max' ' GC_med' ' GC_min' ' GC_mad' ' Q_med',...
%                 ' Q_max' 'Q_min' ' Q_mad' ' alpha_q'};
% dlmwrite( [short_table_string '_full_exp'], header_array ,  'delimiter', '' );
%fprintf( fid, [repmat( '%u\t', 1, 5) '%.4f\t' '%.4f\t' '%e\t' '%.4f\t' '%e\t' '%.4f\t' repmat('%u\t', 1, 3) '%.4f\t' repmat( '%.4f\t', 1, 6) ...
%                                                                       repmat('%u\t',1,3) '%.4f\t' '%.4f\n'], master_short_table' );

%fclose(fid);

disp('writing complete' )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data_extract] = hdfreader2_vec_input_2( big_h5_file, start_index1, end_index1, start_index2, end_index2,...
                                                                                    start_index3, end_index3, start_index4, end_index4 )

% Open the file for read-only access

fileID = H5F.open( big_h5_file, 'H5F_ACC_RDONLY','H5P_DEFAULT');

% Get the dataset of interest
datasetID_GposRd = H5D.open( fileID, 'GposRd' );
datasetID_GposKmer21freq = H5D.open( fileID, 'GposKmer' );
datasetID_GposMedianQ = H5D.open( fileID, 'GposMedianQ' );
datasetID_GcTable = H5D.open( fileID, 'GCtable' );

% Get dataspace

dataspaceID_GposRd = H5D.get_space( datasetID_GposRd );
[rank_GposRd dims_GposRd] = H5S.get_simple_extent_dims(dataspaceID_GposRd);

dataspaceID_GposKmer21freq = H5D.get_space( datasetID_GposKmer21freq );
[rank_GposKmer21freq dims_GposKmer21freq] = H5S.get_simple_extent_dims(dataspaceID_GposKmer21freq);

dataspaceID_GposMedianQ = H5D.get_space( datasetID_GposMedianQ );
[rank_GposMedianQ dims_GposMedianQ] = H5S.get_simple_extent_dims(dataspaceID_GposMedianQ);

dataspaceID_GcTable = H5D.get_space( datasetID_GcTable );
[rank_GcTable dims_GcTable] = H5S.get_simple_extent_dims(dataspaceID_GcTable);


%select region
data_extract = cell( length(start_index1), 4);

stride=[];
block = [];

%size of loop index must be the same for all classes of data

for i = 1:length(start_index1)

%GPOS_RD
dimsRequired = [ min([end_index1(i,1) dims_GposRd(1)]) dims_GposRd(2)]; %do not exceed the maximum file dimension
start = [start_index1(i,1)-1 0];
count=dimsRequired - start;
H5S.select_hyperslab(dataspaceID_GposRd,  'H5S_SELECT_SET', start, stride, count, block);

% Define the memory dataspace.
memspaceID = H5S.create_simple(rank_GposRd, count, []);
%output region
data_extract{i,1} = H5D.read(datasetID_GposRd, 'H5ML_DEFAULT', memspaceID, dataspaceID_GposRd, 'H5P_DEFAULT')';

%Q
% dimsRequired = [ min([end_index2(i,1) dims_GposMedianQ(1)]) dims_GposMedianQ(2)] %do not exceed the maximum file dimension
dimsRequired = [ end_index2(i,1) dims_GposMedianQ(2)];
start = [start_index2(i,1)-1 0];
count=dimsRequired - start;
if dims_GposMedianQ(1) > dimsRequired(1)
H5S.select_hyperslab(dataspaceID_GposMedianQ,  'H5S_SELECT_SET', start, stride, count, block);
% Define the memory dataspace.
memspaceID = H5S.create_simple(rank_GposMedianQ, count, []);
%output region
data_extract{i,2} = H5D.read(datasetID_GposMedianQ, 'H5ML_DEFAULT', memspaceID, dataspaceID_GposMedianQ, 'H5P_DEFAULT')';
else
data_extract{i,2} = [ 0 0 ];
end

%GC

if( start_index4(i,1) == end_index4(i,1) )
end_index4(i,1) = min([ start_index4(i,1)+1 dims_GcTable(1)]);   %assure that there is a minimum granularity in the GC table after compensation by ceil/200
end

dimsRequired = [ min([end_index4(i,1) dims_GcTable(1)]) dims_GcTable(2)]; %do not exceed the maximum file dimension
start = [ max( [1 start_index4(i,1)-1]) 0];
count=dimsRequired - start;
if count(1) == 0
count(1) = 1;
end
H5S.select_hyperslab(dataspaceID_GcTable,  'H5S_SELECT_SET', start, stride, count, block);
% Define the memory dataspace.
memspaceID = H5S.create_simple(rank_GcTable, count, []);
%output region
data_extract{i,4} = H5D.read(datasetID_GcTable, 'H5ML_DEFAULT', memspaceID, dataspaceID_GcTable, 'H5P_DEFAULT')';

%kmer21
dimsRequired = [ min([end_index3(i,1) dims_GposKmer21freq(1)]) dims_GposKmer21freq(2)]; %do not exceed the maximum file dimension
start = [start_index3(i,1)-1 0];
count=dimsRequired - start;
H5S.select_hyperslab(dataspaceID_GposKmer21freq,  'H5S_SELECT_SET', start, stride, count, block);
% Define the memory dataspace.
memspaceID = H5S.create_simple(rank_GposKmer21freq, count, []);
%output region
data_extract{i,3} = H5D.read(datasetID_GposKmer21freq, 'H5ML_DEFAULT', memspaceID, dataspaceID_GposKmer21freq, 'H5P_DEFAULT')';


end

