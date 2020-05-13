% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% mg10hlDE1
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [17x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,AE,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 12:53
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = mg10hlDE1_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'mg10hlDE1_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'mg10hlDE1_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [17 1]), ...
  'mg10hlDE1_jacobig_rot_sym_varpar: pkin has to be [17x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:05
	% EndTime: 2020-04-11 12:36:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:05
	% EndTime: 2020-04-11 12:36:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:05
	% EndTime: 2020-04-11 12:36:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->20)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = t1;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = -t2;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:08
	% EndTime: 2020-04-11 12:36:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (726->38), mult. (1084->68), div. (28->7), fcn. (696->8), ass. (0->58)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = cos(qJ(2));
	t3 = cos(pkin(15));
	t5 = sin(qJ(2));
	t6 = sin(pkin(15));
	t8 = t3 * t2 - t6 * t5;
	t9 = t8 * pkin(2);
	t12 = t6 * t2 + t3 * t5;
	t15 = 0.2e1 * pkin(2) * pkin(1) * t12;
	t16 = pkin(1) ^ 2;
	t20 = -t15 + t16 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t24 = -t15 + t16 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t26 = sqrt(-t24 * t20);
	t27 = t26 * t9;
	t29 = -t12 * pkin(2) + pkin(1);
	t30 = 0.1e1 / t26;
	t32 = pkin(1) * t8;
	t36 = pkin(1) * pkin(2);
	t38 = t24 * pkin(2) * t32 + t36 * t8 * t20;
	t41 = -t12 * pkin(2);
	t42 = pkin(2) ^ 2;
	t43 = pkin(5) ^ 2;
	t44 = pkin(4) ^ 2;
	t45 = -t15 + t16 + t42 - t43 + t44;
	t47 = t8 ^ 2;
	t52 = 0.1e1 / pkin(4);
	t54 = -t15 + t16 + t42;
	t55 = 0.1e1 / t54;
	t58 = t45 * t9;
	t59 = t26 * t29 + t58;
	t61 = t54 ^ 2;
	t62 = 0.1e1 / t61;
	t64 = pkin(2) * t32;
	t69 = t45 * t29 - t27;
	t73 = t59 ^ 2;
	t74 = t69 ^ 2;
	t75 = 0.1e1 / t74;
	t78 = 0.1e1 / (t75 * t73 + 0.1e1);
	t101 = t78 * t54 * pkin(4) / t69 * (t55 * t52 * (-0.2e1 * pkin(1) * t47 * t42 + t38 * t30 * t29 + t45 * t41 - t27) + 0.2e1 * t64 * t62 * t52 * t59) - t78 * t75 * t54 * pkin(4) * t59 * (t55 * t52 * (-0.2e1 * t36 * t8 * t29 - t38 * t30 * t9 - t26 * t41 - t58) + 0.2e1 * t64 * t62 * t52 * t69);
	t104 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t1 * t101 + t1);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = (-t104 * t101 - t104);
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:10
	% EndTime: 2020-04-11 12:36:11
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (995->64), mult. (1226->104), div. (56->14), fcn. (710->10), ass. (0->79)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = cos(qJ(2));
	t3 = cos(pkin(15));
	t5 = sin(qJ(2));
	t6 = sin(pkin(15));
	t8 = t3 * t2 - t6 * t5;
	t9 = t8 * pkin(2);
	t12 = t6 * t2 + t3 * t5;
	t15 = 0.2e1 * pkin(2) * pkin(1) * t12;
	t16 = pkin(1) ^ 2;
	t20 = -t15 + t16 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t24 = -t15 + t16 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t26 = sqrt(-t24 * t20);
	t27 = t26 * t9;
	t29 = -t12 * pkin(2) + pkin(1);
	t30 = 0.1e1 / t26;
	t32 = pkin(1) * t8;
	t36 = pkin(1) * pkin(2);
	t38 = t24 * pkin(2) * t32 + t36 * t8 * t20;
	t41 = -t12 * pkin(2);
	t42 = pkin(2) ^ 2;
	t43 = pkin(5) ^ 2;
	t44 = pkin(4) ^ 2;
	t45 = -t15 + t16 + t42 - t43 + t44;
	t47 = t8 ^ 2;
	t52 = 0.1e1 / pkin(4);
	t54 = -t15 + t16 + t42;
	t55 = 0.1e1 / t54;
	t58 = t45 * t9;
	t59 = t26 * t29 + t58;
	t61 = t54 ^ 2;
	t62 = 0.1e1 / t61;
	t64 = pkin(2) * t32;
	t69 = t45 * t29 - t27;
	t73 = t59 ^ 2;
	t74 = t69 ^ 2;
	t75 = 0.1e1 / t74;
	t78 = 0.1e1 / (t73 * t75 + 0.1e1);
	t101 = t78 * t54 * pkin(4) / t69 * (t55 * t52 * (-0.2e1 * pkin(1) * t47 * t42 + t38 * t30 * t29 + t45 * t41 - t27) + 0.2e1 * t64 * t62 * t52 * t59) - t78 * t75 * t54 * pkin(4) * t59 * (t55 * t52 * (-0.2e1 * t36 * t8 * t29 - t38 * t30 * t9 - t26 * t41 - t58) + 0.2e1 * t64 * t62 * t52 * t69);
	t104 = 0.1e1 / pkin(7);
	t105 = qJ(6) + pkin(8);
	t106 = t105 ^ 2;
	t108 = 0.1e1 / t106 * t104;
	t109 = pkin(6) - pkin(7) - pkin(8) - qJ(6);
	t110 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t111 = t110 * t109;
	t112 = pkin(6) + pkin(7) - pkin(8) - qJ(6);
	t113 = pkin(6) + pkin(7) + pkin(8) + qJ(6);
	t114 = t113 * t112;
	t116 = sqrt(-t114 * t111);
	t119 = 0.1e1 / t105 * t104;
	t133 = (pkin(6) ^ 2);
	t134 = (pkin(7) ^ 2);
	t135 = (pkin(8) ^ 2);
	t138 = (qJ(6) ^ 2);
	t139 = -2 * pkin(8) * qJ(6) + t133 - t134 - t135 - t138;
	t142 = t139 ^ 2;
	t143 = 1 / t142;
	t147 = 0.1e1 / (-t143 * t114 * t111 + 0.1e1);
	t158 = -t147 / t139 * t105 * pkin(7) * (-t116 * t108 + (-t113 * t112 * t109 + t113 * t112 * t110 - t112 * t111 + t113 * t111) / t116 * t119 / 0.2e1) + t147 * t143 * t116 * t105 * pkin(7) * (-0.2e1 * t105 * t119 - t139 * t108);
	t160 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t1 * t101 + t1);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = (t1 * t158);
	unknown(2,1) = 0;
	unknown(2,2) = (-t160 * t101 - t160);
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = -(t160 * t158);
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:16
	% EndTime: 2020-04-11 12:36:16
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (2924->89), mult. (2520->186), div. (404->20), fcn. (852->12), ass. (0->110)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = cos(qJ(2));
	t3 = cos(pkin(15));
	t5 = sin(qJ(2));
	t6 = sin(pkin(15));
	t8 = t3 * t2 - t6 * t5;
	t9 = t8 * pkin(2);
	t12 = t6 * t2 + t3 * t5;
	t15 = 0.2e1 * pkin(2) * pkin(1) * t12;
	t16 = pkin(1) ^ 2;
	t20 = -t15 + t16 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t24 = -t15 + t16 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t26 = sqrt(-t24 * t20);
	t27 = t26 * t9;
	t29 = -t12 * pkin(2) + pkin(1);
	t30 = 0.1e1 / t26;
	t32 = pkin(1) * t8;
	t36 = pkin(1) * pkin(2);
	t38 = t24 * pkin(2) * t32 + t36 * t8 * t20;
	t41 = -t12 * pkin(2);
	t42 = pkin(2) ^ 2;
	t43 = pkin(5) ^ 2;
	t44 = pkin(4) ^ 2;
	t45 = -t15 + t16 + t42 - t43 + t44;
	t47 = t8 ^ 2;
	t52 = 0.1e1 / pkin(4);
	t54 = -t15 + t16 + t42;
	t55 = 0.1e1 / t54;
	t58 = t45 * t9;
	t59 = t26 * t29 + t58;
	t61 = t54 ^ 2;
	t62 = 0.1e1 / t61;
	t64 = pkin(2) * t32;
	t69 = t45 * t29 - t27;
	t73 = t59 ^ 2;
	t74 = t69 ^ 2;
	t75 = 0.1e1 / t74;
	t78 = 0.1e1 / (t73 * t75 + 0.1e1);
	t101 = t78 * t54 * pkin(4) / t69 * (t55 * t52 * (-0.2e1 * pkin(1) * t47 * t42 + t38 * t30 * t29 + t45 * t41 - t27) + 0.2e1 * t64 * t62 * t52 * t59) - t78 * t75 * t54 * pkin(4) * t59 * (t55 * t52 * (-0.2e1 * t36 * t8 * t29 - t38 * t30 * t9 - t26 * t41 - t58) + 0.2e1 * t64 * t62 * t52 * t69);
	t104 = 0.1e1 / pkin(7);
	t105 = qJ(6) + pkin(8);
	t106 = t105 ^ 2;
	t107 = 0.1e1 / t106;
	t108 = t107 * t104;
	t109 = pkin(6) - pkin(7) - pkin(8) - qJ(6);
	t110 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t111 = t110 * t109;
	t112 = pkin(6) + pkin(7) - pkin(8) - qJ(6);
	t113 = pkin(6) + pkin(7) + pkin(8) + qJ(6);
	t114 = t113 * t112;
	t116 = sqrt(-t114 * t111);
	t119 = 0.1e1 / t105 * t104;
	t120 = 0.1e1 / t116;
	t121 = t112 * t110;
	t127 = -t113 * t112 * t109 - t112 * t111 + t113 * t111 + t113 * t121;
	t133 = (pkin(6) ^ 2);
	t134 = (pkin(7) ^ 2);
	t135 = (pkin(8) ^ 2);
	t137 = 2 * pkin(8) * qJ(6);
	t138 = qJ(6) ^ 2;
	t139 = t133 - t134 - t135 - t137 - t138;
	t142 = t139 ^ 2;
	t143 = 1 / t142;
	t147 = 0.1e1 / (-t143 * t114 * t111 + 0.1e1);
	t158 = -t147 / t139 * t105 * pkin(7) * (-t116 * t108 + t127 * t120 * t119 / 0.2e1) + t147 * t143 * t116 * t105 * pkin(7) * (-0.2e1 * t105 * t119 - t139 * t108);
	t160 = 0.1e1 / pkin(6);
	t163 = 0.1e1 / t106 / t105 * t160;
	t164 = t104 * t116;
	t165 = t139 * t164;
	t168 = t107 * t160;
	t177 = t133 - t134 + t135 + t137 + t138;
	t178 = t104 * t177;
	t179 = t116 * t178;
	t182 = 0.2e1 * t104 * t105;
	t191 = t165 * t163 / 0.2e1 - t127 * t139 * t104 * t120 * t168 / 0.8e1 + t105 * t164 * t168 / 0.2e1 - t179 * t163 / 0.2e1 + t116 * t182 * t168 / 0.4e1 + t127 * t120 * t104 * t177 * t168 / 0.8e1;
	t192 = cos(pkin(16));
	t194 = t139 * t178;
	t205 = t104 * t113 * t121;
	t209 = t104 * t114;
	t212 = t109 * t168;
	t222 = -t194 * t163 / 0.2e1 + t139 * t182 * t168 / 0.4e1 - t105 * t178 * t168 / 0.2e1 + t205 * t109 * t163 / 0.2e1 + t209 * t110 * t168 / 0.4e1 - t209 * t212 / 0.4e1 + t104 * t113 * t110 * t212 / 0.4e1 - t104 * t121 * t212 / 0.4e1;
	t223 = sin(pkin(16));
	t228 = t194 * t168 - t205 * t212;
	t232 = -t165 * t168 + t179 * t168;
	t234 = -t192 * t228 / 0.4e1 + t223 * t232 / 0.4e1;
	t239 = t192 * t232 / 0.4e1 + t223 * t228 / 0.4e1;
	t240 = t239 ^ 2;
	t241 = t234 ^ 2;
	t242 = 0.1e1 / t241;
	t245 = 0.1e1 / (t242 * t240 + 0.1e1);
	t253 = t245 / t234 * (t192 * t191 + t223 * t222) - t245 * t242 * t239 * (t223 * t191 - t192 * t222);
	t256 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t1 * t101 + t1);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = (t1 * t158 + t1 * t253);
	unknown(2,1) = 0;
	unknown(2,2) = (-t256 * t101 - t256);
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = (-t256 * t158 - t256 * t253);
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:49
	% EndTime: 2020-04-11 12:36:49
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (8501->111), mult. (9922->247), div. (1004->23), fcn. (5020->21), ass. (0->135)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = cos(qJ(2));
	t3 = cos(pkin(15));
	t5 = sin(qJ(2));
	t6 = sin(pkin(15));
	t8 = t3 * t2 - t6 * t5;
	t9 = t8 * pkin(2);
	t12 = t6 * t2 + t3 * t5;
	t15 = 0.2e1 * pkin(2) * pkin(1) * t12;
	t16 = pkin(1) ^ 2;
	t20 = -t15 + t16 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t24 = -t15 + t16 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t26 = sqrt(-t24 * t20);
	t27 = t26 * t9;
	t29 = -t12 * pkin(2) + pkin(1);
	t30 = 0.1e1 / t26;
	t32 = pkin(1) * t8;
	t36 = pkin(1) * pkin(2);
	t38 = t24 * pkin(2) * t32 + t36 * t8 * t20;
	t41 = -t12 * pkin(2);
	t42 = pkin(2) ^ 2;
	t43 = pkin(5) ^ 2;
	t44 = pkin(4) ^ 2;
	t45 = -t15 + t16 + t42 - t43 + t44;
	t47 = t8 ^ 2;
	t52 = 0.1e1 / pkin(4);
	t54 = -t15 + t16 + t42;
	t55 = 0.1e1 / t54;
	t58 = t45 * t9;
	t59 = t26 * t29 + t58;
	t61 = t54 ^ 2;
	t62 = 0.1e1 / t61;
	t64 = pkin(2) * t32;
	t69 = t45 * t29 - t27;
	t73 = t59 ^ 2;
	t74 = t69 ^ 2;
	t75 = 0.1e1 / t74;
	t78 = 0.1e1 / (t75 * t73 + 0.1e1);
	t101 = t78 * t54 * pkin(4) / t69 * (t55 * t52 * (-0.2e1 * pkin(1) * t47 * t42 + t38 * t30 * t29 + t45 * t41 - t27) + 0.2e1 * t64 * t62 * t52 * t59) - t78 * t75 * t54 * pkin(4) * t59 * (t55 * t52 * (-0.2e1 * t36 * t8 * t29 - t38 * t30 * t9 - t26 * t41 - t58) + 0.2e1 * t64 * t62 * t52 * t69);
	t104 = cos(qJ(1));
	t105 = t5 * t104;
	t108 = 0.1e1 / t44;
	t114 = sqrt(t62 * t108 * t73 + t62 * t108 * t74);
	t116 = 0.1e1 / t114 * t55 * t52;
	t120 = -t116 * t69 * t3 + t116 * t59 * t6;
	t122 = t2 * t104;
	t127 = -t116 * t59 * t3 - t116 * t69 * t6;
	t129 = -t120 * t105 - t127 * t122;
	t130 = cos(pkin(17));
	t131 = 1 / pkin(7);
	t132 = qJ(6) + pkin(8);
	t134 = 0.1e1 / t132 * t131;
	t135 = (pkin(6) - pkin(7) - pkin(8) - qJ(6));
	t136 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t137 = t136 * t135;
	t138 = (pkin(6) + pkin(7) - pkin(8) - qJ(6));
	t139 = (pkin(6) + pkin(7) + pkin(8) + qJ(6));
	t140 = t139 * t138;
	t142 = sqrt(-t140 * t137);
	t144 = (pkin(6) ^ 2);
	t145 = (pkin(7) ^ 2);
	t146 = (pkin(8) ^ 2);
	t148 = 2 * pkin(8) * qJ(6);
	t149 = qJ(6) ^ 2;
	t150 = t144 - t145 - t146 - t148 - t149;
	t152 = atan2(t142 * t134, t150 * t134);
	t153 = t152 + pkin(16);
	t154 = cos(t153);
	t156 = sin(pkin(17));
	t157 = sin(t153);
	t159 = -t154 * t130 - t157 * t156;
	t163 = t127 * t105 - t120 * t122;
	t166 = t157 * t130 - t154 * t156;
	t169 = 0.1e1 / pkin(6);
	t170 = t132 ^ 2;
	t171 = 1 / t170;
	t172 = (t171 * t169);
	t173 = t144 - t145 + t146 + t148 + t149;
	t174 = t131 * t173;
	t175 = t150 * t174;
	t177 = (t135 * t172);
	t178 = (t138 * t136);
	t180 = t131 * t139 * t178;
	t182 = t175 * t172 - t180 * t177;
	t183 = cos(pkin(16));
	t185 = t131 * t142;
	t186 = t150 * t185;
	t188 = t142 * t174;
	t190 = -t186 * t172 + t188 * t172;
	t191 = sin(pkin(16));
	t193 = -t183 * t182 / 0.4e1 + t191 * t190 / 0.4e1;
	t197 = t183 * t190 / 0.4e1 + t191 * t182 / 0.4e1;
	t198 = t197 ^ 2;
	t199 = t193 ^ 2;
	t201 = sqrt(t198 + t199);
	t202 = 0.1e1 / t201;
	t210 = t171 * t131;
	t212 = 0.1e1 / t142;
	t218 = -(t139 * t138 * t135) - t138 * t137 + t139 * t137 + (t139 * t178);
	t226 = t150 ^ 2;
	t227 = 1 / t226;
	t231 = 0.1e1 / (-t227 * t140 * t137 + 0.1e1);
	t242 = -t231 / t150 * t132 * pkin(7) * (-t142 * t210 + t218 * t212 * t134 / 0.2e1) + t231 * t227 * t142 * t132 * pkin(7) * (-0.2e1 * t132 * t134 - (t150 * t210));
	t246 = 0.1e1 / t170 / t132 * t169;
	t259 = 0.2e1 * t131 * t132;
	t268 = t186 * t246 / 0.2e1 - t218 * t150 * t131 * t212 * t172 / 0.8e1 + t132 * t185 * t172 / 0.2e1 - t188 * t246 / 0.2e1 + t142 * t259 * t172 / 0.4e1 + t218 * t212 * t131 * t173 * t172 / 0.8e1;
	t282 = t131 * t140;
	t294 = -t175 * t246 / 0.2e1 + t150 * t259 * t172 / 0.4e1 - t132 * t174 * t172 / 0.2e1 + t180 * t135 * t246 / 0.2e1 + t282 * t136 * t172 / 0.4e1 - t282 * t177 / 0.4e1 + t131 * t139 * t136 * t177 / 0.4e1 - (t131 * t178 * t177) / 0.4e1;
	t299 = 0.1e1 / t199;
	t302 = 0.1e1 / (t299 * t198 + 0.1e1);
	t310 = t302 / t193 * (t183 * t268 + t191 * t294) - t302 * t299 * t197 * (-t183 * t294 + t191 * t268);
	t315 = t5 * t1;
	t317 = t2 * t1;
	t319 = -t120 * t315 - t127 * t317;
	t323 = -t120 * t317 + t127 * t315;
	t339 = t120 * t2 - t127 * t5;
	t343 = -t120 * t5 - t127 * t2;
	unknown(1,1) = 0;
	unknown(1,2) = (t1 * t101 + t1);
	unknown(1,3) = (-t202 * t193 * (t159 * t129 + t166 * t163) - t202 * t197 * (-t166 * t129 + t159 * t163));
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = (t1 * t242 + t1 * t310);
	unknown(2,1) = 0;
	unknown(2,2) = (-t104 * t101 - t104);
	unknown(2,3) = (-t202 * t193 * (t159 * t319 + t166 * t323) - t202 * t197 * (t159 * t323 - t166 * t319));
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = (-t104 * t242 - t104 * t310);
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = (-t202 * t193 * (t159 * t339 + t166 * t343) - t202 * t197 * (t159 * t343 - t166 * t339));
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobig_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:37:05
	% EndTime: 2020-04-11 12:37:06
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (14081->117), mult. (17329->264), div. (1604->23), fcn. (9195->23), ass. (0->143)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = cos(qJ(2));
	t3 = cos(pkin(15));
	t5 = sin(qJ(2));
	t6 = sin(pkin(15));
	t8 = t3 * t2 - t6 * t5;
	t9 = t8 * pkin(2);
	t12 = t6 * t2 + t3 * t5;
	t15 = 0.2e1 * pkin(2) * pkin(1) * t12;
	t16 = pkin(1) ^ 2;
	t20 = -t15 + t16 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t24 = -t15 + t16 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t26 = sqrt(-t24 * t20);
	t27 = t26 * t9;
	t29 = -t12 * pkin(2) + pkin(1);
	t30 = 0.1e1 / t26;
	t32 = pkin(1) * t8;
	t36 = pkin(1) * pkin(2);
	t38 = t24 * pkin(2) * t32 + t36 * t8 * t20;
	t41 = -t12 * pkin(2);
	t42 = pkin(2) ^ 2;
	t43 = pkin(5) ^ 2;
	t44 = pkin(4) ^ 2;
	t45 = -t15 + t16 + t42 - t43 + t44;
	t47 = t8 ^ 2;
	t52 = 0.1e1 / pkin(4);
	t54 = -t15 + t16 + t42;
	t55 = 0.1e1 / t54;
	t58 = t45 * t9;
	t59 = t26 * t29 + t58;
	t61 = t54 ^ 2;
	t62 = 0.1e1 / t61;
	t64 = pkin(2) * t32;
	t69 = t45 * t29 - t27;
	t73 = t59 ^ 2;
	t74 = t69 ^ 2;
	t75 = 0.1e1 / t74;
	t78 = 0.1e1 / (t75 * t73 + 0.1e1);
	t101 = t78 * t54 * pkin(4) / t69 * (t55 * t52 * (-0.2e1 * pkin(1) * t47 * t42 + t38 * t30 * t29 + t45 * t41 - t27) + 0.2e1 * t64 * t62 * t52 * t59) - t78 * t75 * t54 * pkin(4) * t59 * (t55 * t52 * (-0.2e1 * t36 * t8 * t29 - t38 * t30 * t9 - t26 * t41 - t58) + 0.2e1 * t64 * t62 * t52 * t69);
	t104 = cos(qJ(1));
	t105 = t5 * t104;
	t108 = 0.1e1 / t44;
	t114 = sqrt(t62 * t108 * t73 + t62 * t108 * t74);
	t116 = 0.1e1 / t114 * t55 * t52;
	t120 = -t116 * t69 * t3 + t116 * t59 * t6;
	t122 = t2 * t104;
	t127 = -t116 * t59 * t3 - t116 * t69 * t6;
	t129 = -t120 * t105 - t127 * t122;
	t130 = cos(pkin(17));
	t131 = 1 / pkin(7);
	t132 = qJ(6) + pkin(8);
	t134 = 0.1e1 / t132 * t131;
	t135 = (pkin(6) - pkin(7) - pkin(8) - qJ(6));
	t136 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t137 = t136 * t135;
	t138 = (pkin(6) + pkin(7) - pkin(8) - qJ(6));
	t139 = (pkin(6) + pkin(7) + pkin(8) + qJ(6));
	t140 = t139 * t138;
	t142 = sqrt(-t140 * t137);
	t144 = (pkin(6) ^ 2);
	t145 = (pkin(7) ^ 2);
	t146 = (pkin(8) ^ 2);
	t148 = 2 * pkin(8) * qJ(6);
	t149 = qJ(6) ^ 2;
	t150 = t144 - t145 - t146 - t148 - t149;
	t152 = atan2(t142 * t134, t150 * t134);
	t153 = t152 + pkin(16);
	t154 = cos(t153);
	t156 = sin(pkin(17));
	t157 = sin(t153);
	t159 = -t154 * t130 - t157 * t156;
	t163 = t127 * t105 - t120 * t122;
	t166 = t157 * t130 - t154 * t156;
	t168 = t159 * t129 + t166 * t163;
	t169 = 0.1e1 / pkin(6);
	t170 = t132 ^ 2;
	t171 = 1 / t170;
	t172 = (t171 * t169);
	t173 = t144 - t145 + t146 + t148 + t149;
	t174 = t131 * t173;
	t175 = t150 * t174;
	t177 = (t135 * t172);
	t178 = (t138 * t136);
	t180 = t131 * t139 * t178;
	t182 = t175 * t172 - t180 * t177;
	t183 = cos(pkin(16));
	t185 = t131 * t142;
	t186 = t150 * t185;
	t188 = t142 * t174;
	t190 = -t186 * t172 + t188 * t172;
	t191 = sin(pkin(16));
	t193 = -t183 * t182 / 0.4e1 + t191 * t190 / 0.4e1;
	t197 = t183 * t190 / 0.4e1 + t191 * t182 / 0.4e1;
	t198 = t197 ^ 2;
	t199 = t193 ^ 2;
	t201 = sqrt(t198 + t199);
	t202 = 0.1e1 / t201;
	t206 = -t166 * t129 + t159 * t163;
	t215 = sin(qJ(3));
	t217 = cos(qJ(3));
	t220 = t171 * t131;
	t222 = 0.1e1 / t142;
	t228 = -(t139 * t138 * t135) - t138 * t137 + t139 * t137 + (t139 * t178);
	t236 = t150 ^ 2;
	t237 = 1 / t236;
	t241 = 0.1e1 / (-t237 * t140 * t137 + 0.1e1);
	t252 = -t241 / t150 * t132 * pkin(7) * (-t142 * t220 + t228 * t222 * t134 / 0.2e1) + t241 * t237 * t142 * t132 * pkin(7) * (-0.2e1 * t132 * t134 - (t150 * t220));
	t256 = 0.1e1 / t170 / t132 * t169;
	t269 = 0.2e1 * t132 * t131;
	t278 = t186 * t256 / 0.2e1 - t228 * t150 * t131 * t222 * t172 / 0.8e1 + t132 * t185 * t172 / 0.2e1 - t188 * t256 / 0.2e1 + t142 * t269 * t172 / 0.4e1 + t228 * t222 * t131 * t173 * t172 / 0.8e1;
	t292 = t131 * t140;
	t304 = -t175 * t256 / 0.2e1 + t150 * t269 * t172 / 0.4e1 - t132 * t174 * t172 / 0.2e1 + t180 * t135 * t256 / 0.2e1 + t292 * t136 * t172 / 0.4e1 - t292 * t177 / 0.4e1 + t131 * t139 * t136 * t177 / 0.4e1 - (t131 * t178 * t177) / 0.4e1;
	t309 = 0.1e1 / t199;
	t312 = 0.1e1 / (t309 * t198 + 0.1e1);
	t320 = t312 / t193 * (t183 * t278 + t191 * t304) - t312 * t309 * t197 * (-t183 * t304 + t191 * t278);
	t325 = t5 * t1;
	t327 = t2 * t1;
	t329 = -t120 * t325 - t127 * t327;
	t333 = -t120 * t327 + t127 * t325;
	t335 = t159 * t329 + t166 * t333;
	t340 = t159 * t333 - t166 * t329;
	t357 = t120 * t2 - t127 * t5;
	t361 = -t120 * t5 - t127 * t2;
	t363 = t159 * t357 + t166 * t361;
	t368 = t159 * t361 - t166 * t357;
	unknown(1,1) = 0;
	unknown(1,2) = (t1 * t101 + t1);
	unknown(1,3) = (-t202 * t193 * t168 - t202 * t197 * t206);
	unknown(1,4) = (-t215 * (t202 * t197 * t168 - t202 * t193 * t206) + t217 * t1);
	unknown(1,5) = 0;
	unknown(1,6) = (t1 * t252 + t1 * t320);
	unknown(2,1) = 0;
	unknown(2,2) = (-t104 * t101 - t104);
	unknown(2,3) = (-t202 * t193 * t335 - t202 * t197 * t340);
	unknown(2,4) = (-t215 * (-t202 * t193 * t340 + t202 * t197 * t335) - t217 * t104);
	unknown(2,5) = 0;
	unknown(2,6) = (-t104 * t252 - t104 * t320);
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = (-t202 * t193 * t363 - t202 * t197 * t368);
	unknown(3,4) = -(t215 * (-t202 * t193 * t368 + t202 * t197 * t363));
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobig_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:37:58
	% EndTime: 2020-04-11 12:37:59
	% DurationCPUTime: 0.53s
	% Computational Cost: add. (25240->122), mult. (32144->275), div. (2804->23), fcn. (17544->25), ass. (0->151)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = cos(qJ(2));
	t3 = cos(pkin(15));
	t5 = sin(qJ(2));
	t6 = sin(pkin(15));
	t8 = t3 * t2 - t6 * t5;
	t9 = t8 * pkin(2);
	t12 = t6 * t2 + t3 * t5;
	t15 = 0.2e1 * pkin(2) * pkin(1) * t12;
	t16 = pkin(1) ^ 2;
	t20 = -t15 + t16 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t24 = -t15 + t16 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t26 = sqrt(-t24 * t20);
	t27 = t26 * t9;
	t29 = -t12 * pkin(2) + pkin(1);
	t30 = 0.1e1 / t26;
	t32 = pkin(1) * t8;
	t36 = pkin(1) * pkin(2);
	t38 = t24 * pkin(2) * t32 + t36 * t8 * t20;
	t41 = -t12 * pkin(2);
	t42 = pkin(2) ^ 2;
	t43 = pkin(5) ^ 2;
	t44 = pkin(4) ^ 2;
	t45 = -t15 + t16 + t42 - t43 + t44;
	t47 = t8 ^ 2;
	t52 = 0.1e1 / pkin(4);
	t54 = -t15 + t16 + t42;
	t55 = 0.1e1 / t54;
	t58 = t45 * t9;
	t59 = t26 * t29 + t58;
	t61 = t54 ^ 2;
	t62 = 0.1e1 / t61;
	t64 = pkin(2) * t32;
	t69 = t45 * t29 - t27;
	t73 = t59 ^ 2;
	t74 = t69 ^ 2;
	t75 = 0.1e1 / t74;
	t78 = 0.1e1 / (t75 * t73 + 0.1e1);
	t101 = t78 * t54 * pkin(4) / t69 * (t55 * t52 * (-0.2e1 * pkin(1) * t47 * t42 + t38 * t30 * t29 + t45 * t41 - t27) + 0.2e1 * t64 * t62 * t52 * t59) - t78 * t75 * t54 * pkin(4) * t59 * (t55 * t52 * (-0.2e1 * t36 * t8 * t29 - t38 * t30 * t9 - t26 * t41 - t58) + 0.2e1 * t64 * t62 * t52 * t69);
	t104 = cos(qJ(1));
	t105 = t5 * t104;
	t108 = 0.1e1 / t44;
	t114 = sqrt(t62 * t108 * t73 + t62 * t108 * t74);
	t116 = 0.1e1 / t114 * t55 * t52;
	t120 = -t116 * t69 * t3 + t116 * t59 * t6;
	t122 = t2 * t104;
	t127 = -t116 * t59 * t3 - t116 * t69 * t6;
	t129 = -t120 * t105 - t127 * t122;
	t130 = cos(pkin(17));
	t131 = 1 / pkin(7);
	t132 = qJ(6) + pkin(8);
	t134 = 0.1e1 / t132 * t131;
	t135 = (pkin(6) - pkin(7) - pkin(8) - qJ(6));
	t136 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t137 = t136 * t135;
	t138 = (pkin(6) + pkin(7) - pkin(8) - qJ(6));
	t139 = (pkin(6) + pkin(7) + pkin(8) + qJ(6));
	t140 = t139 * t138;
	t142 = sqrt(-t140 * t137);
	t144 = (pkin(6) ^ 2);
	t145 = (pkin(7) ^ 2);
	t146 = (pkin(8) ^ 2);
	t148 = 2 * pkin(8) * qJ(6);
	t149 = qJ(6) ^ 2;
	t150 = t144 - t145 - t146 - t148 - t149;
	t152 = atan2(t142 * t134, t150 * t134);
	t153 = t152 + pkin(16);
	t154 = cos(t153);
	t156 = sin(pkin(17));
	t157 = sin(t153);
	t159 = -t154 * t130 - t157 * t156;
	t163 = t127 * t105 - t120 * t122;
	t166 = t157 * t130 - t154 * t156;
	t168 = t159 * t129 + t166 * t163;
	t169 = 0.1e1 / pkin(6);
	t170 = t132 ^ 2;
	t171 = 1 / t170;
	t172 = (t171 * t169);
	t173 = t144 - t145 + t146 + t148 + t149;
	t174 = t131 * t173;
	t175 = t150 * t174;
	t177 = (t135 * t172);
	t178 = (t138 * t136);
	t180 = t131 * t139 * t178;
	t182 = t175 * t172 - t180 * t177;
	t183 = cos(pkin(16));
	t185 = t131 * t142;
	t186 = t150 * t185;
	t188 = t142 * t174;
	t190 = -t186 * t172 + t188 * t172;
	t191 = sin(pkin(16));
	t193 = -t183 * t182 / 0.4e1 + t191 * t190 / 0.4e1;
	t197 = t183 * t190 / 0.4e1 + t191 * t182 / 0.4e1;
	t198 = t197 ^ 2;
	t199 = t193 ^ 2;
	t201 = sqrt(t198 + t199);
	t202 = 0.1e1 / t201;
	t206 = -t166 * t129 + t159 * t163;
	t209 = -t202 * t193 * t168 - t202 * t197 * t206;
	t214 = t202 * t197 * t168 - t202 * t193 * t206;
	t215 = sin(qJ(3));
	t217 = cos(qJ(3));
	t223 = sin(qJ(4));
	t225 = cos(qJ(4));
	t228 = t171 * t131;
	t230 = 0.1e1 / t142;
	t236 = -(t139 * t138 * t135) - t138 * t137 + t139 * t137 + (t139 * t178);
	t244 = t150 ^ 2;
	t245 = 1 / t244;
	t249 = 0.1e1 / (-t245 * t140 * t137 + 0.1e1);
	t260 = -t249 / t150 * t132 * pkin(7) * (-t142 * t228 + t236 * t230 * t134 / 0.2e1) + t249 * t245 * t142 * t132 * pkin(7) * (-0.2e1 * t132 * t134 - (t150 * t228));
	t264 = 0.1e1 / t170 / t132 * t169;
	t277 = 0.2e1 * t131 * t132;
	t286 = t186 * t264 / 0.2e1 - t236 * t150 * t131 * t230 * t172 / 0.8e1 + t132 * t185 * t172 / 0.2e1 - t188 * t264 / 0.2e1 + t142 * t277 * t172 / 0.4e1 + t236 * t230 * t131 * t173 * t172 / 0.8e1;
	t300 = t131 * t140;
	t312 = -t175 * t264 / 0.2e1 + t150 * t277 * t172 / 0.4e1 - t132 * t174 * t172 / 0.2e1 + t180 * t135 * t264 / 0.2e1 + t300 * t136 * t172 / 0.4e1 - t300 * t177 / 0.4e1 + t131 * t139 * t136 * t177 / 0.4e1 - (t131 * t178 * t177) / 0.4e1;
	t317 = 0.1e1 / t199;
	t320 = 0.1e1 / (t317 * t198 + 0.1e1);
	t328 = t320 / t193 * (t183 * t286 + t191 * t312) - t320 * t317 * t197 * (-t183 * t312 + t191 * t286);
	t333 = t5 * t1;
	t335 = t2 * t1;
	t337 = -t120 * t333 - t127 * t335;
	t341 = -t120 * t335 + t127 * t333;
	t343 = t159 * t337 + t166 * t341;
	t348 = t159 * t341 - t166 * t337;
	t351 = -t202 * t193 * t343 - t202 * t197 * t348;
	t356 = -t202 * t193 * t348 + t202 * t197 * t343;
	t371 = t120 * t2 - t127 * t5;
	t375 = -t120 * t5 - t127 * t2;
	t377 = t159 * t371 + t166 * t375;
	t382 = t159 * t375 - t166 * t371;
	t385 = -t202 * t193 * t377 - t202 * t197 * t382;
	t390 = -t202 * t193 * t382 + t202 * t197 * t377;
	unknown(1,1) = 0;
	unknown(1,2) = (t1 * t101 + t1);
	unknown(1,3) = t209;
	unknown(1,4) = (t217 * t1 - t215 * t214);
	unknown(1,5) = (t223 * (t215 * t1 + t217 * t214) + t225 * t209);
	unknown(1,6) = (t1 * t260 + t1 * t328);
	unknown(2,1) = 0;
	unknown(2,2) = (-t104 * t101 - t104);
	unknown(2,3) = t351;
	unknown(2,4) = (-t217 * t104 - t215 * t356);
	unknown(2,5) = (t223 * (-t215 * t104 + t217 * t356) + t225 * t351);
	unknown(2,6) = (-t104 * t260 - t104 * t328);
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = t385;
	unknown(3,4) = -(t215 * t390);
	unknown(3,5) = (t223 * t217 * t390 + t225 * t385);
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobig_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:15
	% EndTime: 2020-04-11 12:36:15
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (1264->63), mult. (1368->102), div. (84->14), fcn. (724->10), ass. (0->79)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = cos(qJ(2));
	t3 = cos(pkin(15));
	t5 = sin(qJ(2));
	t6 = sin(pkin(15));
	t8 = t3 * t2 - t6 * t5;
	t9 = t8 * pkin(2);
	t12 = t6 * t2 + t3 * t5;
	t15 = 0.2e1 * pkin(2) * pkin(1) * t12;
	t16 = pkin(1) ^ 2;
	t20 = -t15 + t16 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t24 = -t15 + t16 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t26 = sqrt(-t24 * t20);
	t27 = t26 * t9;
	t29 = -t12 * pkin(2) + pkin(1);
	t30 = 0.1e1 / t26;
	t32 = pkin(1) * t8;
	t36 = pkin(1) * pkin(2);
	t38 = t24 * pkin(2) * t32 + t36 * t8 * t20;
	t41 = -t12 * pkin(2);
	t42 = pkin(2) ^ 2;
	t43 = pkin(5) ^ 2;
	t44 = pkin(4) ^ 2;
	t45 = -t15 + t16 + t42 - t43 + t44;
	t47 = t8 ^ 2;
	t52 = 0.1e1 / pkin(4);
	t54 = -t15 + t16 + t42;
	t55 = 0.1e1 / t54;
	t58 = t45 * t9;
	t59 = t26 * t29 + t58;
	t61 = t54 ^ 2;
	t62 = 0.1e1 / t61;
	t64 = pkin(2) * t32;
	t69 = t45 * t29 - t27;
	t73 = t59 ^ 2;
	t74 = t69 ^ 2;
	t75 = 0.1e1 / t74;
	t78 = 0.1e1 / (t75 * t73 + 0.1e1);
	t101 = t78 * t54 * pkin(4) / t69 * (t55 * t52 * (-0.2e1 * pkin(1) * t47 * t42 + t38 * t30 * t29 + t45 * t41 - t27) + 0.2e1 * t64 * t62 * t52 * t59) - t78 * t75 * t54 * pkin(4) * t59 * (t55 * t52 * (-0.2e1 * t36 * t8 * t29 - t38 * t30 * t9 - t26 * t41 - t58) + 0.2e1 * t64 * t62 * t52 * t69);
	t104 = 0.1e1 / pkin(7);
	t105 = qJ(6) + pkin(8);
	t106 = t105 ^ 2;
	t108 = 0.1e1 / t106 * t104;
	t109 = pkin(6) - pkin(7) - pkin(8) - qJ(6);
	t110 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t111 = t110 * t109;
	t112 = pkin(6) + pkin(7) - pkin(8) - qJ(6);
	t113 = pkin(6) + pkin(7) + pkin(8) + qJ(6);
	t114 = t113 * t112;
	t116 = sqrt(-t114 * t111);
	t119 = 0.1e1 / t105 * t104;
	t133 = (pkin(6) ^ 2);
	t134 = (pkin(7) ^ 2);
	t135 = (pkin(8) ^ 2);
	t138 = (qJ(6) ^ 2);
	t139 = -2 * pkin(8) * qJ(6) + t133 - t134 - t135 - t138;
	t142 = t139 ^ 2;
	t143 = 1 / t142;
	t147 = 0.1e1 / (-t143 * t114 * t111 + 0.1e1);
	t158 = -t147 / t139 * t105 * pkin(7) * (-t116 * t108 + (-t113 * t112 * t109 + t113 * t112 * t110 - t112 * t111 + t113 * t111) / t116 * t119 / 0.2e1) + t147 * t143 * t116 * t105 * pkin(7) * (-0.2e1 * t105 * t119 - t139 * t108);
	t162 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t1 * t101 + t1);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = (-t162 * t101 - t162);
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobig_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:07
	% EndTime: 2020-04-11 12:36:07
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (725->36), mult. (1084->69), div. (28->7), fcn. (694->8), ass. (0->57)
	unknown=NaN(3,6);
	t1 = cos(qJ(2));
	t2 = cos(pkin(15));
	t4 = sin(qJ(2));
	t5 = sin(pkin(15));
	t7 = t2 * t1 - t5 * t4;
	t10 = t5 * t1 + t2 * t4;
	t11 = pkin(1) * t10;
	t13 = 0.2e1 * pkin(2) * t11;
	t14 = pkin(1) ^ 2;
	t18 = -t13 + t14 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t22 = -t13 + t14 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t24 = sqrt(-t22 * t18);
	t26 = pkin(1) * t24 * t7;
	t27 = -t11 + pkin(2);
	t28 = 0.1e1 / t24;
	t30 = pkin(1) * t7;
	t34 = pkin(1) * pkin(2);
	t36 = pkin(2) * t22 * t30 + t18 * t34 * t7;
	t40 = pkin(2) ^ 2;
	t41 = pkin(5) ^ 2;
	t42 = pkin(4) ^ 2;
	t43 = -t13 + t14 + t40 + t41 - t42;
	t45 = t7 ^ 2;
	t50 = 0.1e1 / pkin(5);
	t52 = -t13 + t14 + t40;
	t53 = 0.1e1 / t52;
	t56 = t43 * t30;
	t57 = t24 * t27 + t56;
	t59 = t52 ^ 2;
	t60 = 0.1e1 / t59;
	t62 = pkin(2) * t30;
	t67 = t43 * t27 - t26;
	t71 = t57 ^ 2;
	t72 = t67 ^ 2;
	t73 = 0.1e1 / t72;
	t76 = 0.1e1 / (t73 * t71 + 0.1e1);
	t101 = t76 * t52 * pkin(5) / t67 * (t53 * t50 * (-pkin(1) * t10 * t43 - 0.2e1 * pkin(2) * t14 * t45 + t27 * t28 * t36 - t26) + 0.2e1 * t62 * t60 * t50 * t57) - t76 * t73 * t52 * pkin(5) * t57 * (t53 * t50 * (-t36 * pkin(1) * t28 * t7 + pkin(1) * t24 * t10 - 0.2e1 * t34 * t7 * t27 - t56) + 0.2e1 * t62 * t60 * t50 * t67);
	t102 = sin(qJ(1));
	t104 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t102 * t101);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = -(t104 * t101);
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 11
	%% Symbolic Calculation
	% From jacobig_rot_11_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:28
	% EndTime: 2020-04-11 12:36:28
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (1464->71), mult. (1640->114), div. (96->17), fcn. (870->10), ass. (0->86)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = cos(qJ(2));
	t3 = cos(pkin(15));
	t5 = sin(qJ(2));
	t6 = sin(pkin(15));
	t8 = t3 * t2 - t6 * t5;
	t9 = t8 * pkin(2);
	t12 = t6 * t2 + t3 * t5;
	t15 = 0.2e1 * pkin(2) * pkin(1) * t12;
	t16 = pkin(1) ^ 2;
	t20 = -t15 + t16 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t24 = -t15 + t16 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t25 = t24 * t20;
	t26 = sqrt(-t25);
	t27 = t26 * t9;
	t29 = -t12 * pkin(2) + pkin(1);
	t30 = 0.1e1 / t26;
	t32 = pkin(1) * t8;
	t36 = pkin(1) * pkin(2);
	t38 = t24 * pkin(2) * t32 + t36 * t8 * t20;
	t41 = -t12 * pkin(2);
	t42 = pkin(2) ^ 2;
	t43 = pkin(5) ^ 2;
	t44 = pkin(4) ^ 2;
	t45 = -t15 + t16 + t42 - t43 + t44;
	t47 = t8 ^ 2;
	t52 = 0.1e1 / pkin(4);
	t54 = -t15 + t16 + t42;
	t55 = 0.1e1 / t54;
	t58 = t45 * t9;
	t59 = t26 * t29 + t58;
	t61 = t54 ^ 2;
	t62 = 0.1e1 / t61;
	t64 = pkin(2) * t32;
	t69 = t45 * t29 - t27;
	t73 = t59 ^ 2;
	t74 = t69 ^ 2;
	t75 = 0.1e1 / t74;
	t78 = 0.1e1 / (t75 * t73 + 0.1e1);
	t82 = 0.2e1 * t38 * t30;
	t101 = t78 * t54 * pkin(4) / t69 * (t55 * t52 * (-0.2e1 * pkin(1) * t47 * t42 + t38 * t30 * t29 + t45 * t41 - t27) + 0.2e1 * t64 * t62 * t52 * t59) - t78 * t75 * t54 * pkin(4) * t59 * (t55 * t52 * (-t26 * t41 - t82 * t9 / 0.2e1 - t58 - 0.2e1 * t36 * t8 * t29) + 0.2e1 * t64 * t62 * t52 * t69);
	t103 = t15 - t16 - t42 + t43 + t44;
	t105 = t103 ^ 2;
	t106 = 0.1e1 / t105;
	t109 = 0.1e1 / (-t106 * t25 + 0.1e1);
	t117 = t109 / t103 * t82 / 0.2e1 - 0.2e1 * t109 * t106 * t26 * t64;
	t120 = 0.1e1 / pkin(7);
	t121 = qJ(6) + pkin(8);
	t122 = t121 ^ 2;
	t124 = 0.1e1 / t122 * t120;
	t125 = pkin(6) - pkin(7) - pkin(8) - qJ(6);
	t126 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t127 = t126 * t125;
	t128 = pkin(6) + pkin(7) - pkin(8) - qJ(6);
	t129 = pkin(6) + pkin(7) + pkin(8) + qJ(6);
	t130 = t129 * t128;
	t132 = sqrt(-t130 * t127);
	t135 = 0.1e1 / t121 * t120;
	t149 = (pkin(6) ^ 2);
	t150 = (pkin(7) ^ 2);
	t151 = (pkin(8) ^ 2);
	t154 = (qJ(6) ^ 2);
	t155 = -2 * pkin(8) * qJ(6) + t149 - t150 - t151 - t154;
	t158 = t155 ^ 2;
	t159 = 1 / t158;
	t163 = 0.1e1 / (-t159 * t130 * t127 + 0.1e1);
	t174 = -t163 / t155 * t121 * pkin(7) * (-t132 * t124 + (-t129 * t128 * t125 + t129 * t128 * t126 - t128 * t127 + t129 * t127) / t132 * t135 / 0.2e1) + t163 * t159 * t132 * t121 * pkin(7) * (-0.2e1 * t121 * t135 - t155 * t124);
	t178 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t1 * t101 + t1 * t117 + t1);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = (-t178 * t101 - t178 * t117 - t178);
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 12
	%% Symbolic Calculation
	% From jacobig_rot_12_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:09
	% EndTime: 2020-04-11 12:36:09
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (1670->78), mult. (1714->156), div. (208->16), fcn. (744->10), ass. (0->94)
	unknown=NaN(3,6);
	t1 = cos(qJ(2));
	t2 = cos(pkin(15));
	t4 = sin(qJ(2));
	t5 = sin(pkin(15));
	t7 = t2 * t1 - t5 * t4;
	t10 = t5 * t1 + t2 * t4;
	t11 = pkin(1) * t10;
	t13 = 0.2e1 * pkin(2) * t11;
	t14 = pkin(1) ^ 2;
	t18 = -t13 + t14 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t22 = -t13 + t14 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t24 = sqrt(-t22 * t18);
	t26 = pkin(1) * t24 * t7;
	t27 = -t11 + pkin(2);
	t28 = 0.1e1 / t24;
	t30 = pkin(1) * t7;
	t34 = pkin(1) * pkin(2);
	t36 = t22 * pkin(2) * t30 + t34 * t7 * t18;
	t40 = pkin(2) ^ 2;
	t41 = pkin(5) ^ 2;
	t42 = pkin(4) ^ 2;
	t43 = -t13 + t14 + t40 + t41 - t42;
	t45 = t7 ^ 2;
	t50 = 0.1e1 / pkin(5);
	t52 = -t13 + t14 + t40;
	t53 = 0.1e1 / t52;
	t56 = t43 * t30;
	t57 = t24 * t27 + t56;
	t59 = t52 ^ 2;
	t60 = 0.1e1 / t59;
	t62 = pkin(2) * t30;
	t67 = t43 * t27 - t26;
	t71 = t57 ^ 2;
	t72 = t67 ^ 2;
	t73 = 0.1e1 / t72;
	t76 = 0.1e1 / (t73 * t71 + 0.1e1);
	t101 = t76 * t52 * pkin(5) / t67 * (t53 * t50 * (-t43 * pkin(1) * t10 - 0.2e1 * pkin(2) * t45 * t14 + t36 * t28 * t27 - t26) + 0.2e1 * t62 * t60 * t50 * t57) - t76 * t73 * t52 * pkin(5) * t57 * (t53 * t50 * (-t36 * pkin(1) * t28 * t7 + pkin(1) * t24 * t10 - 0.2e1 * t34 * t7 * t27 - t56) + 0.2e1 * t62 * t60 * t50 * t67);
	t102 = sin(qJ(1));
	t104 = 0.1e1 / pkin(6);
	t105 = qJ(6) + pkin(8);
	t106 = t105 ^ 2;
	t109 = 0.1e1 / t106 / t105 * t104;
	t110 = pkin(6) - pkin(7) - pkin(8) - qJ(6);
	t111 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t112 = t111 * t110;
	t113 = pkin(6) + pkin(7) - pkin(8) - qJ(6);
	t114 = pkin(6) + pkin(7) + pkin(8) + qJ(6);
	t115 = t114 * t113;
	t117 = sqrt(-t115 * t112);
	t118 = 0.1e1 / pkin(7);
	t119 = t118 * t117;
	t120 = (pkin(6) ^ 2);
	t121 = (pkin(7) ^ 2);
	t122 = (pkin(8) ^ 2);
	t124 = 2 * pkin(8) * qJ(6);
	t125 = qJ(6) ^ 2;
	t126 = t120 - t121 - t122 - t124 - t125;
	t127 = t126 * t119;
	t131 = 0.1e1 / t106 * t104;
	t132 = 0.1e1 / t117;
	t135 = t113 * t111;
	t141 = -t114 * t113 * t110 - t113 * t112 + t114 * t112 + t114 * t135;
	t148 = t120 - t121 + t122 + t124 + t125;
	t149 = t118 * t148;
	t150 = t117 * t149;
	t153 = 0.2e1 * t118 * t105;
	t163 = t126 * t149;
	t165 = t110 * t131;
	t167 = t118 * t114 * t135;
	t169 = t163 * t131 - t167 * t165;
	t174 = -t127 * t131 + t150 * t131;
	t177 = 0.16e2 / t169 ^ 2;
	t180 = 0.1e1 / (0.1e1 + t177 * t174 ^ 2 / 0.16e2);
	t194 = t118 * t115;
	t210 = 0.4e1 * t180 / t169 * (t127 * t109 / 0.2e1 - t141 * t126 * t118 * t132 * t131 / 0.8e1 + t105 * t119 * t131 / 0.2e1 - t150 * t109 / 0.2e1 + t117 * t153 * t131 / 0.4e1 + t141 * t132 * t118 * t148 * t131 / 0.8e1) - t180 * t177 * t174 * (-t163 * t109 / 0.2e1 + t126 * t153 * t131 / 0.4e1 - t105 * t149 * t131 / 0.2e1 + t167 * t110 * t109 / 0.2e1 + t194 * t111 * t131 / 0.4e1 - t194 * t165 / 0.4e1 + t118 * t114 * t111 * t165 / 0.4e1 - t118 * t135 * t165 / 0.4e1) / 0.4e1;
	t212 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t102 * t101);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = (t102 * t210);
	unknown(2,1) = 0;
	unknown(2,2) = -(t212 * t101);
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = -(t212 * t210);
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
else
	Jg_rot=NaN(3,6);
end