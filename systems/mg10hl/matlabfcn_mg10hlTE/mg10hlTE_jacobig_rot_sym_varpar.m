% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% mg10hlTE
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
% Datum: 2020-04-11 12:29
% Revision: f5729c120b58d1b0137ade1aa9321d1ea2b3cc0a (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = mg10hlTE_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'mg10hlTE_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'mg10hlTE_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [17 1]), ...
  'mg10hlTE_jacobig_rot_sym_varpar: pkin has to be [17x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:25:21
	% EndTime: 2020-04-11 12:25:21
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
	% StartTime: 2020-04-11 12:25:21
	% EndTime: 2020-04-11 12:25:21
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
	% StartTime: 2020-04-11 12:25:21
	% EndTime: 2020-04-11 12:25:21
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
	% StartTime: 2020-04-11 12:25:22
	% EndTime: 2020-04-11 12:25:22
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (726->38), mult. (1084->68), div. (28->7), fcn. (696->8), ass. (0->58)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = cos(qJ(2));
	t3 = cos(pkin(15));
	t5 = sin(qJ(2));
	t6 = sin(pkin(15));
	t8 = t2 * t3 - t5 * t6;
	t9 = pkin(2) * t8;
	t12 = t2 * t6 + t5 * t3;
	t15 = 0.2e1 * t12 * pkin(1) * pkin(2);
	t16 = pkin(1) ^ 2;
	t20 = -t15 + t16 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t24 = -t15 + t16 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t26 = sqrt(-t20 * t24);
	t27 = t9 * t26;
	t29 = -pkin(2) * t12 + pkin(1);
	t30 = 0.1e1 / t26;
	t32 = t8 * pkin(1);
	t36 = pkin(1) * pkin(2);
	t38 = t32 * pkin(2) * t24 + t20 * t8 * t36;
	t41 = -pkin(2) * t12;
	t42 = pkin(2) ^ 2;
	t43 = pkin(5) ^ 2;
	t44 = pkin(4) ^ 2;
	t45 = -t15 + t16 + t42 - t43 + t44;
	t47 = t8 ^ 2;
	t52 = 0.1e1 / pkin(4);
	t54 = -t15 + t16 + t42;
	t55 = 0.1e1 / t54;
	t58 = t9 * t45;
	t59 = t29 * t26 + t58;
	t61 = t54 ^ 2;
	t62 = 0.1e1 / t61;
	t64 = t32 * pkin(2);
	t69 = t29 * t45 - t27;
	t73 = t59 ^ 2;
	t74 = t69 ^ 2;
	t75 = 0.1e1 / t74;
	t78 = 0.1e1 / (t73 * t75 + 0.1e1);
	t101 = ((-0.2e1 * t42 * t47 * pkin(1) + t29 * t30 * t38 + t41 * t45 - t27) * t52 * t55 + 0.2e1 * t59 * t52 * t62 * t64) / t69 * pkin(4) * t54 * t78 - ((-0.2e1 * t29 * t8 * t36 - t9 * t30 * t38 - t41 * t26 - t58) * t52 * t55 + 0.2e1 * t69 * t52 * t62 * t64) * t59 * pkin(4) * t54 * t75 * t78;
	t104 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t101 * t1 + t1);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = (-t101 * t104 - t104);
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
	% StartTime: 2020-04-11 12:25:22
	% EndTime: 2020-04-11 12:25:23
	% DurationCPUTime: 0.09s
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
	t78 = 0.1e1 / (t75 * t73 + 0.1e1);
	t101 = t78 * t54 * pkin(4) / t69 * (t55 * t52 * (-0.2e1 * pkin(1) * t47 * t42 + t38 * t30 * t29 + t45 * t41 - t27) + 0.2e1 * t64 * t52 * t59 * t62) - t78 * t75 * t54 * pkin(4) * t59 * (t55 * t52 * (-0.2e1 * t36 * t8 * t29 - t38 * t30 * t9 - t26 * t41 - t58) + 0.2e1 * t64 * t62 * t52 * t69);
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
	t158 = -t147 / t139 * t105 * pkin(7) * (-t108 * t116 + (-t113 * t112 * t109 + t113 * t112 * t110 - t112 * t111 + t113 * t111) / t116 * t119 / 0.2e1) + t147 * t143 * t116 * t105 * pkin(7) * (-0.2e1 * t105 * t119 - t139 * t108);
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
	% StartTime: 2020-04-11 12:25:24
	% EndTime: 2020-04-11 12:25:24
	% DurationCPUTime: 0.10s
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
	t78 = 0.1e1 / (t75 * t73 + 0.1e1);
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
	% StartTime: 2020-04-11 12:25:35
	% EndTime: 2020-04-11 12:25:36
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (5123->109), mult. (5146->247), div. (668->20), fcn. (2344->17), ass. (0->131)
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
	t107 = t55 * t52;
	t111 = -t107 * t69 * t3 + t107 * t59 * t6;
	t113 = t2 * t104;
	t118 = -t107 * t59 * t3 - t107 * t69 * t6;
	t120 = -t111 * t105 / 0.2e1 - t118 * t113 / 0.2e1;
	t121 = cos(pkin(17));
	t122 = 1 / pkin(7);
	t123 = qJ(6) + pkin(8);
	t125 = 0.1e1 / t123 * t122;
	t126 = (pkin(6) - pkin(7) - pkin(8) - qJ(6));
	t127 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t128 = t127 * t126;
	t129 = (pkin(6) + pkin(7) - pkin(8) - qJ(6));
	t130 = (pkin(6) + pkin(7) + pkin(8) + qJ(6));
	t131 = t130 * t129;
	t133 = sqrt(-t131 * t128);
	t135 = (pkin(6) ^ 2);
	t136 = (pkin(7) ^ 2);
	t137 = (pkin(8) ^ 2);
	t139 = 2 * pkin(8) * qJ(6);
	t140 = qJ(6) ^ 2;
	t141 = t135 - t136 - t137 - t139 - t140;
	t143 = atan2(t133 * t125, t141 * t125);
	t144 = t143 + pkin(16);
	t145 = cos(t144);
	t147 = sin(pkin(17));
	t148 = sin(t144);
	t150 = -t145 * t121 - t148 * t147;
	t154 = t118 * t105 / 0.2e1 - t111 * t113 / 0.2e1;
	t157 = t148 * t121 - t145 * t147;
	t160 = 0.1e1 / pkin(6);
	t161 = t123 ^ 2;
	t162 = 1 / t161;
	t163 = (t162 * t160);
	t164 = t135 - t136 + t137 + t139 + t140;
	t165 = t122 * t164;
	t166 = t141 * t165;
	t168 = (t126 * t163);
	t169 = (t129 * t127);
	t171 = t122 * t130 * t169;
	t173 = t166 * t163 - t171 * t168;
	t174 = cos(pkin(16));
	t176 = t122 * t133;
	t177 = t141 * t176;
	t179 = t133 * t165;
	t181 = -t177 * t163 + t179 * t163;
	t182 = sin(pkin(16));
	t184 = t174 * t173 / 0.4e1 - t182 * t181 / 0.4e1;
	t191 = -t174 * t181 / 0.4e1 - t182 * t173 / 0.4e1;
	t194 = t162 * t122;
	t196 = 0.1e1 / t133;
	t202 = -(t130 * t129 * t126) - t129 * t128 + t130 * t128 + (t130 * t169);
	t210 = t141 ^ 2;
	t211 = 1 / t210;
	t215 = 0.1e1 / (-t211 * t131 * t128 + 0.1e1);
	t226 = -t215 / t141 * t123 * pkin(7) * (-t133 * t194 + t202 * t196 * t125 / 0.2e1) + t215 * t211 * t133 * t123 * pkin(7) * (-0.2e1 * t123 * t125 - (t141 * t194));
	t230 = 0.1e1 / t161 / t123 * t160;
	t243 = 0.2e1 * t122 * t123;
	t252 = t177 * t230 / 0.2e1 - t202 * t141 * t122 * t196 * t163 / 0.8e1 + t123 * t176 * t163 / 0.2e1 - t179 * t230 / 0.2e1 + t133 * t243 * t163 / 0.4e1 + t202 * t196 * t122 * t164 * t163 / 0.8e1;
	t266 = t122 * t131;
	t278 = -t166 * t230 / 0.2e1 + t141 * t243 * t163 / 0.4e1 - t123 * t165 * t163 / 0.2e1 + t171 * t126 * t230 / 0.2e1 + t266 * t127 * t163 / 0.4e1 - t266 * t168 / 0.4e1 + t122 * t130 * t127 * t168 / 0.4e1 - (t122 * t169 * t168) / 0.4e1;
	t283 = t191 ^ 2;
	t284 = t184 ^ 2;
	t285 = 0.1e1 / t284;
	t288 = 0.1e1 / (t285 * t283 + 0.1e1);
	t296 = -t288 / t184 * (t174 * t252 + t182 * t278) + t288 * t285 * t191 * (-t174 * t278 + t182 * t252);
	t301 = t5 * t1;
	t303 = t2 * t1;
	t305 = -t111 * t301 / 0.2e1 - t118 * t303 / 0.2e1;
	t309 = t118 * t301 / 0.2e1 - t111 * t303 / 0.2e1;
	t323 = t111 * t2 / 0.2e1 - t118 * t5 / 0.2e1;
	t327 = -t118 * t2 / 0.2e1 - t111 * t5 / 0.2e1;
	unknown(1,1) = 0;
	unknown(1,2) = (t1 * t101 + t1);
	unknown(1,3) = (t184 * (t150 * t120 + t157 * t154) + t191 * (-t157 * t120 + t150 * t154));
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = (t1 * t226 + t1 * t296);
	unknown(2,1) = 0;
	unknown(2,2) = (-t104 * t101 - t104);
	unknown(2,3) = (t184 * (t150 * t305 + t157 * t309) + t191 * (t150 * t309 - t157 * t305));
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = (-t104 * t226 - t104 * t296);
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = (t184 * (t150 * t323 + t157 * t327) + t191 * (t150 * t327 - t157 * t323));
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobig_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:25:43
	% EndTime: 2020-04-11 12:25:44
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (7325->115), mult. (7777->258), div. (932->20), fcn. (3843->19), ass. (0->139)
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
	t107 = t55 * t52;
	t111 = -t107 * t69 * t3 + t107 * t59 * t6;
	t113 = t2 * t104;
	t118 = -t107 * t59 * t3 - t107 * t69 * t6;
	t120 = -t111 * t105 / 0.2e1 - t118 * t113 / 0.2e1;
	t121 = cos(pkin(17));
	t122 = 1 / pkin(7);
	t123 = qJ(6) + pkin(8);
	t125 = 0.1e1 / t123 * t122;
	t126 = (pkin(6) - pkin(7) - pkin(8) - qJ(6));
	t127 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t128 = t127 * t126;
	t129 = (pkin(6) + pkin(7) - pkin(8) - qJ(6));
	t130 = (pkin(6) + pkin(7) + pkin(8) + qJ(6));
	t131 = t130 * t129;
	t133 = sqrt(-t131 * t128);
	t135 = (pkin(6) ^ 2);
	t136 = (pkin(7) ^ 2);
	t137 = (pkin(8) ^ 2);
	t139 = 2 * pkin(8) * qJ(6);
	t140 = qJ(6) ^ 2;
	t141 = t135 - t136 - t137 - t139 - t140;
	t143 = atan2(t133 * t125, t141 * t125);
	t144 = t143 + pkin(16);
	t145 = cos(t144);
	t147 = sin(pkin(17));
	t148 = sin(t144);
	t150 = -t145 * t121 - t148 * t147;
	t154 = t118 * t105 / 0.2e1 - t111 * t113 / 0.2e1;
	t157 = t148 * t121 - t145 * t147;
	t159 = t150 * t120 + t157 * t154;
	t160 = 0.1e1 / pkin(6);
	t161 = t123 ^ 2;
	t162 = 1 / t161;
	t163 = (t162 * t160);
	t164 = t135 - t136 + t137 + t139 + t140;
	t165 = t122 * t164;
	t166 = t141 * t165;
	t168 = (t126 * t163);
	t169 = (t129 * t127);
	t171 = t122 * t130 * t169;
	t173 = t166 * t163 - t171 * t168;
	t174 = cos(pkin(16));
	t176 = t122 * t133;
	t177 = t141 * t176;
	t179 = t133 * t165;
	t181 = -t177 * t163 + t179 * t163;
	t182 = sin(pkin(16));
	t184 = t174 * t173 / 0.4e1 - t182 * t181 / 0.4e1;
	t188 = -t157 * t120 + t150 * t154;
	t191 = -t174 * t181 / 0.4e1 - t182 * t173 / 0.4e1;
	t197 = sin(qJ(3));
	t199 = cos(qJ(3));
	t202 = t162 * t122;
	t204 = 0.1e1 / t133;
	t210 = -(t130 * t129 * t126) - t129 * t128 + t130 * t128 + (t130 * t169);
	t218 = t141 ^ 2;
	t219 = 1 / t218;
	t223 = 0.1e1 / (-t219 * t131 * t128 + 0.1e1);
	t234 = -t223 / t141 * t123 * pkin(7) * (-t133 * t202 + t210 * t204 * t125 / 0.2e1) + t223 * t219 * t133 * t123 * pkin(7) * (-0.2e1 * t123 * t125 - (t141 * t202));
	t238 = 0.1e1 / t161 / t123 * t160;
	t251 = 0.2e1 * t122 * t123;
	t260 = t177 * t238 / 0.2e1 - t210 * t141 * t122 * t204 * t163 / 0.8e1 + t123 * t176 * t163 / 0.2e1 - t179 * t238 / 0.2e1 + t133 * t251 * t163 / 0.4e1 + t210 * t204 * t122 * t164 * t163 / 0.8e1;
	t274 = t122 * t131;
	t286 = -t166 * t238 / 0.2e1 + t141 * t251 * t163 / 0.4e1 - t123 * t165 * t163 / 0.2e1 + t171 * t126 * t238 / 0.2e1 + t274 * t127 * t163 / 0.4e1 - t274 * t168 / 0.4e1 + t122 * t130 * t127 * t168 / 0.4e1 - (t122 * t169 * t168) / 0.4e1;
	t291 = t191 ^ 2;
	t292 = t184 ^ 2;
	t293 = 0.1e1 / t292;
	t296 = 0.1e1 / (t293 * t291 + 0.1e1);
	t304 = -t296 / t184 * (t174 * t260 + t182 * t286) + t296 * t293 * t191 * (-t174 * t286 + t182 * t260);
	t309 = t5 * t1;
	t311 = t2 * t1;
	t313 = -t111 * t309 / 0.2e1 - t118 * t311 / 0.2e1;
	t317 = t118 * t309 / 0.2e1 - t111 * t311 / 0.2e1;
	t319 = t150 * t313 + t157 * t317;
	t323 = t150 * t317 - t157 * t313;
	t337 = t111 * t2 / 0.2e1 - t118 * t5 / 0.2e1;
	t341 = -t118 * t2 / 0.2e1 - t111 * t5 / 0.2e1;
	t343 = t150 * t337 + t157 * t341;
	t347 = t150 * t341 - t157 * t337;
	unknown(1,1) = 0;
	unknown(1,2) = (t1 * t101 + t1);
	unknown(1,3) = (t184 * t159 + t191 * t188);
	unknown(1,4) = (-t197 * (-t191 * t159 + t184 * t188) + t199 * t1);
	unknown(1,5) = 0;
	unknown(1,6) = (t1 * t234 + t1 * t304);
	unknown(2,1) = 0;
	unknown(2,2) = (-t104 * t101 - t104);
	unknown(2,3) = (t184 * t319 + t191 * t323);
	unknown(2,4) = (-t197 * (t184 * t323 - t191 * t319) - t199 * t104);
	unknown(2,5) = 0;
	unknown(2,6) = (-t104 * t234 - t104 * t304);
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = (t184 * t343 + t191 * t347);
	unknown(3,4) = -(t197 * (t184 * t347 - t191 * t343));
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobig_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:25:58
	% EndTime: 2020-04-11 12:25:59
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (11728->120), mult. (13040->269), div. (1460->20), fcn. (6840->21), ass. (0->147)
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
	t107 = t55 * t52;
	t111 = -t107 * t69 * t3 + t107 * t59 * t6;
	t113 = t2 * t104;
	t118 = -t107 * t59 * t3 - t107 * t69 * t6;
	t120 = -t111 * t105 / 0.2e1 - t118 * t113 / 0.2e1;
	t121 = cos(pkin(17));
	t122 = 1 / pkin(7);
	t123 = qJ(6) + pkin(8);
	t125 = 0.1e1 / t123 * t122;
	t126 = (pkin(6) - pkin(7) - pkin(8) - qJ(6));
	t127 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t128 = t127 * t126;
	t129 = (pkin(6) + pkin(7) - pkin(8) - qJ(6));
	t130 = (pkin(6) + pkin(7) + pkin(8) + qJ(6));
	t131 = t130 * t129;
	t133 = sqrt(-t131 * t128);
	t135 = (pkin(6) ^ 2);
	t136 = (pkin(7) ^ 2);
	t137 = (pkin(8) ^ 2);
	t139 = 2 * pkin(8) * qJ(6);
	t140 = qJ(6) ^ 2;
	t141 = t135 - t136 - t137 - t139 - t140;
	t143 = atan2(t133 * t125, t141 * t125);
	t144 = t143 + pkin(16);
	t145 = cos(t144);
	t147 = sin(pkin(17));
	t148 = sin(t144);
	t150 = -t145 * t121 - t148 * t147;
	t154 = t118 * t105 / 0.2e1 - t111 * t113 / 0.2e1;
	t157 = t148 * t121 - t145 * t147;
	t159 = t150 * t120 + t157 * t154;
	t160 = 0.1e1 / pkin(6);
	t161 = t123 ^ 2;
	t162 = 1 / t161;
	t163 = (t162 * t160);
	t164 = t135 - t136 + t137 + t139 + t140;
	t165 = t122 * t164;
	t166 = t141 * t165;
	t168 = (t126 * t163);
	t169 = (t129 * t127);
	t171 = t122 * t130 * t169;
	t173 = t166 * t163 - t171 * t168;
	t174 = cos(pkin(16));
	t176 = t122 * t133;
	t177 = t141 * t176;
	t179 = t133 * t165;
	t181 = -t177 * t163 + t179 * t163;
	t182 = sin(pkin(16));
	t184 = t174 * t173 / 0.4e1 - t182 * t181 / 0.4e1;
	t188 = -t157 * t120 + t150 * t154;
	t191 = -t174 * t181 / 0.4e1 - t182 * t173 / 0.4e1;
	t193 = t184 * t159 + t191 * t188;
	t196 = -t191 * t159 + t184 * t188;
	t197 = sin(qJ(3));
	t199 = cos(qJ(3));
	t205 = sin(qJ(4));
	t207 = cos(qJ(4));
	t210 = t162 * t122;
	t212 = 0.1e1 / t133;
	t218 = -(t130 * t129 * t126) - t129 * t128 + t130 * t128 + (t130 * t169);
	t226 = t141 ^ 2;
	t227 = 1 / t226;
	t231 = 0.1e1 / (-t227 * t131 * t128 + 0.1e1);
	t242 = -t231 / t141 * t123 * pkin(7) * (-t133 * t210 + t218 * t212 * t125 / 0.2e1) + t231 * t227 * t133 * t123 * pkin(7) * (-0.2e1 * t123 * t125 - (t141 * t210));
	t246 = 0.1e1 / t161 / t123 * t160;
	t259 = 0.2e1 * t122 * t123;
	t268 = t177 * t246 / 0.2e1 - t218 * t141 * t122 * t212 * t163 / 0.8e1 + t123 * t176 * t163 / 0.2e1 - t179 * t246 / 0.2e1 + t133 * t259 * t163 / 0.4e1 + t218 * t212 * t122 * t164 * t163 / 0.8e1;
	t282 = t122 * t131;
	t294 = -t166 * t246 / 0.2e1 + t141 * t259 * t163 / 0.4e1 - t123 * t165 * t163 / 0.2e1 + t171 * t126 * t246 / 0.2e1 + t282 * t127 * t163 / 0.4e1 - t282 * t168 / 0.4e1 + t122 * t130 * t127 * t168 / 0.4e1 - (t122 * t169 * t168) / 0.4e1;
	t299 = t191 ^ 2;
	t300 = t184 ^ 2;
	t301 = 0.1e1 / t300;
	t304 = 0.1e1 / (t301 * t299 + 0.1e1);
	t312 = -t304 / t184 * (t174 * t268 + t182 * t294) + t304 * t301 * t191 * (-t174 * t294 + t182 * t268);
	t317 = t5 * t1;
	t319 = t2 * t1;
	t321 = -t111 * t317 / 0.2e1 - t118 * t319 / 0.2e1;
	t325 = t118 * t317 / 0.2e1 - t111 * t319 / 0.2e1;
	t327 = t150 * t321 + t157 * t325;
	t331 = t150 * t325 - t157 * t321;
	t333 = t184 * t327 + t191 * t331;
	t336 = t184 * t331 - t191 * t327;
	t351 = t111 * t2 / 0.2e1 - t118 * t5 / 0.2e1;
	t355 = -t118 * t2 / 0.2e1 - t111 * t5 / 0.2e1;
	t357 = t150 * t351 + t157 * t355;
	t361 = t150 * t355 - t157 * t351;
	t363 = t184 * t357 + t191 * t361;
	t366 = t184 * t361 - t191 * t357;
	unknown(1,1) = 0;
	unknown(1,2) = (t1 * t101 + t1);
	unknown(1,3) = t193;
	unknown(1,4) = (t199 * t1 - t197 * t196);
	unknown(1,5) = (t205 * (t197 * t1 + t199 * t196) + t207 * t193);
	unknown(1,6) = (t1 * t242 + t1 * t312);
	unknown(2,1) = 0;
	unknown(2,2) = (-t104 * t101 - t104);
	unknown(2,3) = t333;
	unknown(2,4) = (-t199 * t104 - t197 * t336);
	unknown(2,5) = (t205 * (-t197 * t104 + t199 * t336) + t207 * t333);
	unknown(2,6) = (-t104 * t242 - t104 * t312);
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = t363;
	unknown(3,4) = -(t197 * t366);
	unknown(3,5) = (t205 * t199 * t366 + t207 * t363);
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobig_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:25:24
	% EndTime: 2020-04-11 12:25:25
	% DurationCPUTime: 0.06s
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
	t101 = t78 * t54 * pkin(4) / t69 * (t55 * t52 * (-0.2e1 * pkin(1) * t47 * t42 + t38 * t30 * t29 + t45 * t41 - t27) + 0.2e1 * t64 * t52 * t59 * t62) - t78 * t75 * t54 * pkin(4) * t59 * (t55 * t52 * (-0.2e1 * t36 * t8 * t29 - t38 * t30 * t9 - t26 * t41 - t58) + 0.2e1 * t64 * t62 * t52 * t69);
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
	% StartTime: 2020-04-11 12:25:21
	% EndTime: 2020-04-11 12:25:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (725->36), mult. (1084->69), div. (28->7), fcn. (694->8), ass. (0->57)
	unknown=NaN(3,6);
	t1 = cos(qJ(2));
	t2 = cos(pkin(15));
	t4 = sin(qJ(2));
	t5 = sin(pkin(15));
	t7 = t1 * t2 - t4 * t5;
	t10 = t1 * t5 + t4 * t2;
	t11 = t10 * pkin(1);
	t13 = 0.2e1 * t11 * pkin(2);
	t14 = pkin(1) ^ 2;
	t18 = -t13 + t14 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t22 = -t13 + t14 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t24 = sqrt(-t18 * t22);
	t26 = t7 * t24 * pkin(1);
	t27 = -t11 + pkin(2);
	t28 = 0.1e1 / t24;
	t30 = t7 * pkin(1);
	t34 = pkin(1) * pkin(2);
	t36 = t30 * pkin(2) * t22 + t18 * t7 * t34;
	t40 = pkin(2) ^ 2;
	t41 = pkin(5) ^ 2;
	t42 = pkin(4) ^ 2;
	t43 = -t13 + t14 + t40 + t41 - t42;
	t45 = t7 ^ 2;
	t50 = 0.1e1 / pkin(5);
	t52 = -t13 + t14 + t40;
	t53 = 0.1e1 / t52;
	t56 = t30 * t43;
	t57 = t27 * t24 + t56;
	t59 = t52 ^ 2;
	t60 = 0.1e1 / t59;
	t62 = t30 * pkin(2);
	t67 = t27 * t43 - t26;
	t71 = t57 ^ 2;
	t72 = t67 ^ 2;
	t73 = 0.1e1 / t72;
	t76 = 0.1e1 / (t71 * t73 + 0.1e1);
	t101 = ((-t10 * pkin(1) * t43 - 0.2e1 * t14 * t45 * pkin(2) + t27 * t28 * t36 - t26) * t50 * t53 + 0.2e1 * t57 * t50 * t60 * t62) / t67 * pkin(5) * t52 * t76 - ((-t7 * t28 * pkin(1) * t36 + t10 * t24 * pkin(1) - 0.2e1 * t27 * t7 * t34 - t56) * t50 * t53 + 0.2e1 * t67 * t50 * t60 * t62) * t57 * pkin(5) * t52 * t73 * t76;
	t102 = sin(qJ(1));
	t104 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t101 * t102);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = -(t101 * t104);
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
	% StartTime: 2020-04-11 12:25:27
	% EndTime: 2020-04-11 12:25:27
	% DurationCPUTime: 0.07s
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
	% StartTime: 2020-04-11 12:25:22
	% EndTime: 2020-04-11 12:25:22
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (1670->78), mult. (1714->156), div. (208->16), fcn. (744->10), ass. (0->94)
	unknown=NaN(3,6);
	t1 = cos(qJ(2));
	t2 = cos(pkin(15));
	t4 = sin(qJ(2));
	t5 = sin(pkin(15));
	t7 = t1 * t2 - t4 * t5;
	t10 = t1 * t5 + t4 * t2;
	t11 = t10 * pkin(1);
	t13 = 0.2e1 * t11 * pkin(2);
	t14 = pkin(1) ^ 2;
	t18 = -t13 + t14 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t22 = -t13 + t14 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t24 = sqrt(-t18 * t22);
	t26 = t7 * t24 * pkin(1);
	t27 = -t11 + pkin(2);
	t28 = 0.1e1 / t24;
	t30 = t7 * pkin(1);
	t34 = pkin(1) * pkin(2);
	t36 = t30 * pkin(2) * t22 + t18 * t7 * t34;
	t40 = pkin(2) ^ 2;
	t41 = pkin(5) ^ 2;
	t42 = pkin(4) ^ 2;
	t43 = -t13 + t14 + t40 + t41 - t42;
	t45 = t7 ^ 2;
	t50 = 0.1e1 / pkin(5);
	t52 = -t13 + t14 + t40;
	t53 = 0.1e1 / t52;
	t56 = t30 * t43;
	t57 = t27 * t24 + t56;
	t59 = t52 ^ 2;
	t60 = 0.1e1 / t59;
	t62 = t30 * pkin(2);
	t67 = t27 * t43 - t26;
	t71 = t57 ^ 2;
	t72 = t67 ^ 2;
	t73 = 0.1e1 / t72;
	t76 = 0.1e1 / (t71 * t73 + 0.1e1);
	t101 = ((-t10 * pkin(1) * t43 - 0.2e1 * t14 * t45 * pkin(2) + t27 * t28 * t36 - t26) * t50 * t53 + 0.2e1 * t57 * t50 * t60 * t62) / t67 * pkin(5) * t52 * t76 - ((-t7 * t28 * pkin(1) * t36 + t10 * t24 * pkin(1) - 0.2e1 * t27 * t7 * t34 - t56) * t50 * t53 + 0.2e1 * t67 * t50 * t60 * t62) * t57 * pkin(5) * t52 * t73 * t76;
	t102 = sin(qJ(1));
	t104 = 0.1e1 / pkin(6);
	t105 = qJ(6) + pkin(8);
	t106 = t105 ^ 2;
	t109 = t104 / t106 / t105;
	t110 = pkin(6) - pkin(7) - pkin(8) - qJ(6);
	t111 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t112 = t110 * t111;
	t113 = pkin(6) + pkin(7) - pkin(8) - qJ(6);
	t114 = pkin(6) + pkin(7) + pkin(8) + qJ(6);
	t115 = t113 * t114;
	t117 = sqrt(-t112 * t115);
	t118 = 0.1e1 / pkin(7);
	t119 = t117 * t118;
	t120 = (pkin(6) ^ 2);
	t121 = (pkin(7) ^ 2);
	t122 = (pkin(8) ^ 2);
	t124 = 2 * pkin(8) * qJ(6);
	t125 = qJ(6) ^ 2;
	t126 = t120 - t121 - t122 - t124 - t125;
	t127 = t119 * t126;
	t131 = t104 / t106;
	t132 = 0.1e1 / t117;
	t135 = t111 * t113;
	t141 = -t110 * t113 * t114 - t112 * t113 + t112 * t114 + t135 * t114;
	t148 = t120 - t121 + t122 + t124 + t125;
	t149 = t148 * t118;
	t150 = t149 * t117;
	t153 = 0.2e1 * t105 * t118;
	t163 = t149 * t126;
	t165 = t131 * t110;
	t167 = t135 * t114 * t118;
	t169 = t131 * t163 - t165 * t167;
	t174 = -t131 * t127 + t131 * t150;
	t177 = 0.16e2 / t169 ^ 2;
	t180 = 0.1e1 / (0.1e1 + t174 ^ 2 * t177 / 0.16e2);
	t194 = t115 * t118;
	t210 = 0.4e1 * (t109 * t127 / 0.2e1 - t131 * t132 * t118 * t126 * t141 / 0.8e1 + t131 * t119 * t105 / 0.2e1 - t109 * t150 / 0.2e1 + t131 * t153 * t117 / 0.4e1 + t131 * t148 * t118 * t132 * t141 / 0.8e1) / t169 * t180 - (-t109 * t163 / 0.2e1 + t131 * t153 * t126 / 0.4e1 - t131 * t149 * t105 / 0.2e1 + t109 * t110 * t167 / 0.2e1 + t131 * t111 * t194 / 0.4e1 - t165 * t194 / 0.4e1 + t165 * t111 * t114 * t118 / 0.4e1 - t165 * t135 * t118 / 0.4e1) * t174 * t177 * t180 / 0.4e1;
	t212 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t101 * t102);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = (t210 * t102);
	unknown(2,1) = 0;
	unknown(2,2) = -(t101 * t212);
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = -(t210 * t212);
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