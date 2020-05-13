% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% fourbar1DE1
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% JgD_rot [3x1]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = fourbar1DE1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_jacobigD_rot_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1DE1_jacobigD_rot_sym_varpar: qJD has to be [1x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1DE1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_jacobigD_rot_sym_varpar: pkin has to be [4x1] (double)');
JgD_rot=NaN(3,1);
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:57:34
	% EndTime: 2020-04-24 19:57:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0; 0; 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:57:34
	% EndTime: 2020-04-24 19:57:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0; 0; 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:57:34
	% EndTime: 2020-04-24 19:57:34
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (145->21), mult. (198->34), div. (3->4), fcn. (51->4), ass. (0->21)
	t127 = pkin(1) ^ 2;
	t125 = cos(qJ(1));
	t133 = pkin(2) * t125;
	t131 = -0.2e1 * pkin(1) * t133 + t127;
	t137 = -pkin(3) - pkin(4);
	t117 = (pkin(2) - t137) * (pkin(2) + t137) + t131;
	t136 = -pkin(3) + pkin(4);
	t118 = (pkin(2) - t136) * (pkin(2) + t136) + t131;
	t124 = sin(qJ(1));
	t134 = pkin(2) * t124;
	t130 = qJD(1) * t134;
	t126 = pkin(2) ^ 2;
	t121 = t126 + t131;
	t135 = pkin(1) / t121;
	t138 = (-t117 - t118) * pkin(1) * t130 * t135;
	t132 = t117 * t118;
	t128 = sqrt(-t132);
	t116 = 0.1e1 / t128;
	t122 = pkin(1) - t133;
	t119 = -pkin(3) ^ 2 + pkin(4) ^ 2 + t121;
	t1 = [0; 0; (-(0.2e1 * t126 * t124 ^ 2 * pkin(1) + (t119 * t125 + t124 * t128) * pkin(2)) * qJD(1) * t135 + (-0.1e1 / t132 * t138 + 0.2e1 * t127 / t121 ^ 2 * t130) * (t119 * t134 + t122 * t128) - t122 * t116 * t138) * t116;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:57:34
	% EndTime: 2020-04-24 19:57:34
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (145->21), mult. (198->33), div. (3->4), fcn. (51->4), ass. (0->22)
	t126 = pkin(1) ^ 2;
	t124 = cos(qJ(1));
	t133 = pkin(1) * t124;
	t137 = -2 * pkin(2);
	t130 = t133 * t137 + t126;
	t136 = (-pkin(3) - pkin(4));
	t116 = ((pkin(2) - t136) * (pkin(2) + t136)) + t130;
	t135 = (-pkin(3) + pkin(4));
	t117 = ((pkin(2) - t135) * (pkin(2) + t135)) + t130;
	t123 = sin(qJ(1));
	t134 = pkin(1) * t123;
	t129 = qJD(1) * t134;
	t125 = pkin(2) ^ 2;
	t120 = t125 + t130;
	t132 = pkin(2) / t120;
	t138 = (-t116 - t117) * pkin(2) * t129 * t132;
	t131 = t116 * t117;
	t127 = sqrt(-t131);
	t115 = 0.1e1 / t127;
	t121 = -pkin(2) + t133;
	t118 = pkin(3) ^ 2 - pkin(4) ^ 2 + t120;
	t1 = [0; 0; (-(t126 * t123 ^ 2 * t137 + (-t118 * t124 - t123 * t127) * pkin(1)) * qJD(1) * t132 + (-0.1e1 / t131 * t138 + 0.2e1 * t125 / t120 ^ 2 * t129) * (-t118 * t134 + t121 * t127) - t121 * t115 * t138) * t115;];
	JgD_rot = t1;
end