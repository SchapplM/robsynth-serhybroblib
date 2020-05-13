% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% fourbar1turnTE
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% JR_rot [9x2]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = fourbar1turnTE_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_jacobiR_rot_sym_varpar: qJ has to be [2x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1turnTE_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_jacobiR_rot_sym_varpar: pkin has to be [5x1] (double)');
JR_rot=NaN(9,2);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:25
	% EndTime: 2020-04-12 19:20:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:25
	% EndTime: 2020-04-12 19:20:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0; t9, 0; 0, 0; -t9, 0; -t8, 0; 0, 0; 0, 0; 0, 0; 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:25
	% EndTime: 2020-04-12 19:20:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t7 = sin(qJ(2));
	t8 = sin(qJ(1));
	t14 = t8 * t7;
	t9 = cos(qJ(2));
	t13 = t8 * t9;
	t10 = cos(qJ(1));
	t12 = t10 * t7;
	t11 = t10 * t9;
	t1 = [-t13, -t12; t11, -t14; 0, t9; t14, -t11; -t12, -t13; 0, -t7; t10, 0; t8, 0; 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:26
	% EndTime: 2020-04-12 19:20:26
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (1294->42), mult. (1856->85), div. (88->4), fcn. (562->6), ass. (0->48)
	t120 = pkin(1) ^ 2;
	t116 = cos(qJ(2));
	t141 = pkin(1) * t116;
	t149 = -2 * pkin(2);
	t131 = t141 * t149 + t120;
	t110 = (pkin(2) ^ 2) + t131;
	t108 = 0.1e1 / t110;
	t114 = sin(qJ(2));
	t113 = t114 ^ 2;
	t130 = pkin(1) * pkin(2) / t110 ^ 2;
	t107 = pkin(3) ^ 2 - pkin(4) ^ 2 + t110;
	t142 = pkin(1) * t114;
	t104 = t107 * t142;
	t111 = -pkin(2) + t141;
	t146 = -pkin(3) - pkin(4);
	t105 = (pkin(2) - t146) * (pkin(2) + t146) + t131;
	t145 = -pkin(3) + pkin(4);
	t106 = (pkin(2) - t145) * (pkin(2) + t145) + t131;
	t121 = sqrt(-t105 * t106);
	t100 = -t111 * t121 + t104;
	t133 = t116 * t100;
	t135 = t114 * t121;
	t99 = -pkin(1) * t135 - t111 * t107;
	t123 = (t113 * t99 + t114 * t133) * t130;
	t137 = 0.1e1 / t121 * (-t105 - t106) * pkin(2) * t142;
	t93 = t104 + (-t121 * t116 + (t111 * t149 - t137) * t114) * pkin(1);
	t148 = -t93 / 0.2e1;
	t127 = t100 / 0.2e1 + t148;
	t147 = -t99 / 0.2e1;
	t94 = -t111 * t137 + 0.2e1 * t120 * t113 * pkin(2) + (t116 * t107 + t135) * pkin(1);
	t129 = t147 - t94 / 0.2e1;
	t151 = (t127 * t114 + t129 * t116) * t108 + t123;
	t144 = t114 / 0.2e1;
	t143 = t116 / 0.2e1;
	t115 = sin(qJ(1));
	t118 = 0.1e1 / pkin(3);
	t136 = t108 * t118;
	t126 = t136 * t143;
	t128 = t114 * t136;
	t140 = (t99 * t126 - t100 * t128 / 0.2e1) * t115;
	t117 = cos(qJ(1));
	t139 = (t100 * t126 + t99 * t128 / 0.2e1) * t117;
	t138 = t116 * t99;
	t134 = t115 * t118;
	t132 = t117 * t118;
	t124 = (-t100 * t113 + t114 * t138) * t130;
	t122 = (t124 + (-t129 * t114 + t127 * t116) * t108) * t118;
	t1 = [t140, ((t116 * t148 + t94 * t144) * t108 + t124) * t132 + t139; (-t138 / 0.2e1 + t100 * t144) * t108 * t132, t115 * t122; 0, t151 * t118; (-t133 / 0.2e1 + t114 * t147) * t108 * t134, -t151 * t132; t139, ((t94 * t143 + t93 * t144) * t108 - t123) * t134 + t140; 0, t122; t117, 0; t115, 0; 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:26
	% EndTime: 2020-04-12 19:20:26
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (510->27), mult. (728->49), div. (32->4), fcn. (211->6), ass. (0->32)
	t90 = sin(qJ(2));
	t110 = pkin(2) * t90;
	t104 = pkin(1) * t110;
	t92 = cos(qJ(2));
	t109 = pkin(2) * t92;
	t105 = (-0.2e1 * t109 + pkin(1)) * pkin(1);
	t95 = pkin(2) ^ 2;
	t87 = t95 + t105;
	t103 = 0.1e1 / t87 ^ 2 * t104;
	t112 = -pkin(3) - pkin(4);
	t82 = (pkin(2) - t112) * (pkin(2) + t112) + t105;
	t111 = -pkin(3) + pkin(4);
	t83 = (pkin(2) - t111) * (pkin(2) + t111) + t105;
	t97 = sqrt(-t82 * t83);
	t106 = t90 * t97;
	t108 = 0.1e1 / t97 * (-t82 - t83) * t104;
	t113 = 0.2e1 * pkin(1);
	t84 = -pkin(3) ^ 2 + pkin(4) ^ 2 + t87;
	t81 = t84 * t110;
	t88 = pkin(1) - t109;
	t77 = t88 * t97 + t81;
	t85 = 0.1e1 / t87;
	t94 = 0.1e1 / pkin(4);
	t98 = ((t88 * t108 + t95 * t90 ^ 2 * t113 + (t92 * t84 + t106) * pkin(2)) * t85 / 0.2e1 - t77 * t103) * t94;
	t107 = t85 * t94;
	t91 = sin(qJ(1));
	t102 = t91 * t107 / 0.2e1;
	t93 = cos(qJ(1));
	t101 = -t93 * t107 / 0.2e1;
	t76 = -pkin(2) * t106 + t88 * t84;
	t99 = (-(t81 + (-t92 * t97 + (t88 * t113 - t108) * t90) * pkin(2)) * t85 / 0.2e1 + t76 * t103) * t94;
	t1 = [t76 * t102, t93 * t99; t76 * t101, t91 * t99; 0, t98; t77 * t102, -t93 * t98; t77 * t101, -t91 * t98; 0, t99; t93, 0; t91, 0; 0, 0;];
	JR_rot = t1;
end