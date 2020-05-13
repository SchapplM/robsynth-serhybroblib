% Zeitableitung der geometrischen Jacobi-Matrix für beliebiges Segment von
% picker2Dm2DE2
% Use Code from Maple symbolic Code Generation
% 
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% 
% Output:
% Jg [6x2]
%   Zeitableitung der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:02
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD = picker2Dm2DE2_jacobigD_sym_varpar(qJ, qJD, link_index, r_i_i_C, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2DE2_jacobigD_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'picker2Dm2DE2_jacobigD_sym_varpar: qJD has to be [2x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'picker2Dm2DE2_jacobigD_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'picker2Dm2DE2_jacobigD_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2DE2_jacobigD_sym_varpar: pkin has to be [9x1] (double)');

% Function calls
JaD_transl = picker2Dm2DE2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin);
JgD_rot = picker2Dm2DE2_jacobigD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin);

JgD = [JaD_transl; JgD_rot];
end