% Zeitableitung der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh4m1OL
% Use Code from Maple symbolic Code Generation
% 
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% 
% Input:
% qJ [8x1]
%   Generalized joint coordinates (joint angles)
% qJD [8x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,CB,CE,EP,OT,TA,TD]';
% 
% Output:
% Jg [6x8]
%   Zeitableitung der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 23:04
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD = palh4m1OL_jacobigD_sym_varpar(qJ, qJD, link_index, r_i_i_C, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(8,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [8 1]), ...
  'palh4m1OL_jacobigD_sym_varpar: qJ has to be [8x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [8 1]), ...
  'palh4m1OL_jacobigD_sym_varpar: qJD has to be [8x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh4m1OL_jacobigD_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1OL_jacobigD_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'palh4m1OL_jacobigD_sym_varpar: pkin has to be [7x1] (double)');

% Function calls
JaD_transl = palh4m1OL_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin);
JgD_rot = palh4m1OL_jacobigD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin);

JgD = [JaD_transl; JgD_rot];
end