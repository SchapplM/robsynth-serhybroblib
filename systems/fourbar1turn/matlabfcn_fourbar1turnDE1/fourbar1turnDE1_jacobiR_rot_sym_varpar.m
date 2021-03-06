% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% fourbar1turnDE1
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
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = fourbar1turnDE1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_jacobiR_rot_sym_varpar: qJ has to be [2x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1turnDE1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_jacobiR_rot_sym_varpar: pkin has to be [5x1] (double)');
JR_rot=NaN(9,2);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:28:03
	% EndTime: 2020-04-12 19:28:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:28:03
	% EndTime: 2020-04-12 19:28:03
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
	% StartTime: 2020-04-12 19:28:03
	% EndTime: 2020-04-12 19:28:03
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
	% StartTime: 2020-04-12 19:28:04
	% EndTime: 2020-04-12 19:28:05
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (5454->40), mult. (7584->73), div. (360->8), fcn. (2162->11), ass. (0->47)
	t141 = pkin(1) ^ 2;
	t136 = cos(qJ(2));
	t165 = pkin(1) * t136;
	t168 = -2 * pkin(2);
	t154 = t165 * t168 + t141;
	t131 = (pkin(2) ^ 2) + t154;
	t143 = pkin(3) ^ 2;
	t128 = -pkin(4) ^ 2 + t131 + t143;
	t132 = -pkin(2) + t165;
	t134 = sin(qJ(2));
	t167 = -pkin(3) - pkin(4);
	t126 = (pkin(2) - t167) * (pkin(2) + t167) + t154;
	t166 = pkin(4) - pkin(3);
	t127 = (pkin(2) - t166) * (pkin(2) + t166) + t154;
	t142 = sqrt(-t126 * t127);
	t158 = t134 * t142;
	t120 = -pkin(1) * t158 - t132 * t128;
	t164 = t134 * pkin(1);
	t125 = t128 * t164;
	t121 = -t132 * t142 + t125;
	t129 = 0.1e1 / t131;
	t138 = 0.1e1 / pkin(3);
	t161 = t129 * t138;
	t116 = atan2(t121 * t161, t120 * t161);
	t113 = sin(t116);
	t114 = cos(t116);
	t170 = t136 * t113 + t134 * t114;
	t144 = t120 ^ 2;
	t118 = 0.1e1 / t144;
	t119 = t121 ^ 2;
	t130 = 0.1e1 / t131 ^ 2;
	t153 = t130 * t168;
	t162 = 0.1e1 / t142 * (-t126 - t127) * pkin(2) * t164;
	t163 = t121 * t134;
	t151 = (((0.2e1 * t141 * t134 ^ 2 * pkin(2) - t132 * t162) * t129 + ((t136 * t128 + t158) * t129 + t153 * t163) * pkin(1)) / t120 - (t125 * t129 + (-t142 * t136 * t129 + ((t132 * t168 - t162) * t129 + t120 * t153) * t134) * pkin(1)) * t121 * t118) * pkin(3) / (t119 * t118 + 0.1e1) * t131 * t138 * ((t119 + t144) / t143 * t130) ^ (-0.1e1 / 0.2e1) * t161;
	t169 = (t120 * t136 - t163) * t151;
	t160 = t134 * t113;
	t156 = t136 * t114;
	t137 = cos(qJ(1));
	t155 = t170 * t137;
	t149 = -t156 + t160;
	t147 = (t120 * t134 + t121 * t136) * t151;
	t146 = -t160 + t169;
	t145 = t147 + t170;
	t135 = sin(qJ(1));
	t110 = t135 * t156;
	t1 = [-t135 * t160 + t110, t137 * t147 + t155; t149 * t137, t145 * t135; 0, t149 - t169; -t170 * t135, (t146 + t156) * t137; t155, t135 * t146 + t110; 0, t145; t137, 0; t135, 0; 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:28:04
	% EndTime: 2020-04-12 19:28:04
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (2338->34), mult. (3232->68), div. (144->8), fcn. (906->11), ass. (0->37)
	t135 = -pkin(3) - pkin(4);
	t134 = -pkin(3) + pkin(4);
	t112 = sin(qJ(2));
	t133 = pkin(2) * t112;
	t114 = cos(qJ(2));
	t132 = pkin(2) * t114;
	t128 = (-0.2e1 * t132 + pkin(1)) * pkin(1);
	t104 = (pkin(2) - t135) * (pkin(2) + t135) + t128;
	t105 = (pkin(2) - t134) * (pkin(2) + t134) + t128;
	t120 = sqrt(-t104 * t105);
	t127 = pkin(1) * t133;
	t131 = (-t104 - t105) * t127 / t120;
	t118 = pkin(2) ^ 2;
	t109 = t118 + t128;
	t107 = 0.1e1 / t109;
	t116 = 0.1e1 / pkin(4);
	t130 = t107 * t116;
	t129 = t112 * t120;
	t121 = pkin(4) ^ 2;
	t106 = -pkin(3) ^ 2 + t109 + t121;
	t103 = t106 * t133;
	t108 = 0.1e1 / t109 ^ 2;
	t110 = pkin(1) - t132;
	t98 = -pkin(2) * t129 + t110 * t106;
	t122 = t98 ^ 2;
	t96 = 0.1e1 / t122;
	t99 = t110 * t120 + t103;
	t97 = t99 ^ 2;
	t125 = 0.2e1 * (-((t110 * t131 + (t114 * t106 + t129) * pkin(2)) * t107 / 0.2e1 + (-pkin(2) * t108 * t99 + t107 * t118 * t112) * t112 * pkin(1)) / t98 - (-(t103 + (-t112 * t131 - t114 * t120) * pkin(2)) * t107 / 0.2e1 + (-t107 * t110 + t108 * t98) * t127) * t99 * t96) * pkin(4) * t109 * t116 / (t97 * t96 + 0.1e1) * ((t122 + t97) / t121 * t108) ^ (-0.1e1 / 0.2e1) * t130;
	t113 = sin(qJ(1));
	t124 = t113 * t125;
	t115 = cos(qJ(1));
	t123 = t115 * t125;
	t93 = atan2(t99 * t130 / 0.2e1, -t98 * t130 / 0.2e1);
	t92 = cos(t93);
	t91 = sin(t93);
	t1 = [-t113 * t92, -t99 * t123; t115 * t92, -t99 * t124; 0, -t98 * t125; t113 * t91, t98 * t123; -t115 * t91, t98 * t124; 0, -t99 * t125; t115, 0; t113, 0; 0, 0;];
	JR_rot = t1;
end