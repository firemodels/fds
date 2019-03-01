% Emanuele Gissi
% 6-7-2017
% geom_positive_errors.m

% This script catches errors produced during FDS setup of bad geometries cases.
% The success condition is the existance of the error message. The error message has
% its 'ERROR' label replaced by 'SUCCESS' label, thanks to the MISC parameter
% POSITIVE_ERROR_TEST=.TRUE. in the bad geometries cases.
% The absence of the 'ERROR' label makes current Firebot to consider the run a success.

% Clean up

close all
clear all

% Cases and error strings

infile{1} = '../../Verification/Complex_Geometry/geom_bad_inconsistent_normals.err';
errstring{1} = 'Face normals are probably pointing in the wrong direction';

infile{2} = '../../Verification/Complex_Geometry/geom_bad_open_surface.err';
errstring{2} = 'Non manifold geometry, open surface at edge';

infile{3} = '../../Verification/Complex_Geometry/geom_bad_self_intersection.err';
errstring{3} = 'self-intersection';

infile{4} = '../../Verification/Complex_Geometry/geom_bad_non_manifold_edge.err';
errstring{4} = 'Non manifold geometry, more than two triangles share edge';

infile{5} = '../../Verification/Complex_Geometry/geom_bad_intersection.err';
errstring{5} = 'intersection with other GEOM line geometry';

infile{6} = '../../Verification/Complex_Geometry/geom_bad_inverted_normals.err';
errstring{6} = 'Face normals are probably pointing in the wrong direction';

% FIXME This error condition is currently not caught
% It will be caught when boolean operations are performed
%infile{7} = '../../Verification/Complex_Geometry/geom_bad_non_manifold_vert.err';
%errstring{7} = 'Non manifold geometry at vertex'

% Check existance of errstring in each infile

for n = 1:6
    % infile exists?
    if ~exist(infile{n},'file')
        display(['Error: File ',infile{n},' does not exist. Skipping case.'])
        continue
    end

    % errstring in infile?
    errfile = fileread(infile{n});
    if isempty(strfind(errfile, errstring{n}));
        display(['Error: File ',infile{n},':'])
        display(['  Does not contain the following positive error message:'])
        display(['  <',errstring{n},'>'])
    end
end
