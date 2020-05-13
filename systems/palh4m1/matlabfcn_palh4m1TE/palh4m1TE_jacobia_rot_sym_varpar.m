% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh4m1TE
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in palh4m1TE_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AD,CB,CE,EP,HC,OT,TA,TD]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 21:48
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = palh4m1TE_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1TE_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1TE_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1TE_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->15)
	unknown=NaN(3,5);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	Ja_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (6->3), div. (5->2), fcn. (6->2), ass. (0->21)
	unknown=NaN(3,5);
	t1 = sin(qJ(1));
	t2 = t1 ^ 2;
	t3 = cos(qJ(1));
	t4 = t3 ^ 2;
	t6 = t2 / t4;
	t8 = 0.1e1 / (0.1e1 + t6);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(3,1) = (t6 * t8 + t8);
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	Ja_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->15)
	unknown=NaN(3,5);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	Ja_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (6375->73), mult. (6989->226), div. (399->13), fcn. (2035->9), ass. (0->91)
	unknown=NaN(3,5);
	t1 = cos(qJ(1));
	t2 = cos(qJ(5));
	t3 = sin(qJ(5));
	t4 = t3 * pkin(1);
	t6 = 0.2e1 * t4 * pkin(2);
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t11 * t15);
	t19 = t2 * t17 * pkin(1);
	t20 = t4 - pkin(2);
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t27 = -t20 * t25 - t19;
	t28 = t1 * t27;
	t30 = t2 * pkin(1);
	t31 = t30 * t25;
	t32 = -t20 * t17 + t31;
	t33 = 0.1e1 / t32;
	t34 = sin(qJ(1));
	t35 = t34 ^ 2;
	t36 = t27 ^ 2;
	t38 = t32 ^ 2;
	t39 = 0.1e1 / t38;
	t42 = 0.1e1 / (t35 * t36 * t39 + 0.1e1);
	t43 = t33 * t42;
	t45 = 0.1e1 / t17;
	t46 = t2 * t45;
	t51 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t15 - 0.2e1 * t11 * (-qJ(2) - pkin(6) + pkin(3));
	t56 = -t46 * pkin(1) * t51 / 0.2e1 - 0.2e1 * t20 * t22;
	t58 = 0.1e1 / t22;
	t59 = -t6 + t7 + t21;
	t60 = 0.1e1 / t59;
	t61 = t58 * t60;
	t63 = t34 * t27;
	t64 = 0.1e1 / t23;
	t65 = t64 * t60;
	t69 = t22 * t59;
	t70 = t69 * t42;
	t72 = -t20 * t45;
	t76 = t72 * t51 / 0.2e1 + 0.2e1 * t30 * t22;
	t85 = t69 * t39 * t42;
	t87 = (-t34 * t56 * t61 + t63 * t65) * t33 * t70 / 0.2e1 + (-t32 * t64 * t60 + t76 * t58 * t60) * t34 * t27 * t85 / 0.2e1;
	t93 = pkin(1) * pkin(2);
	t95 = t30 * pkin(2) * t15 + t11 * t2 * t93;
	t102 = t3 * t17 * pkin(1) - t46 * pkin(1) * t95 + 0.2e1 * t20 * t2 * t93 - t31;
	t106 = t63 * t58;
	t107 = t59 ^ 2;
	t108 = 0.1e1 / t107;
	t110 = t108 * t2 * t93;
	t118 = t2 ^ 2;
	t122 = -0.2e1 * t7 * t118 * pkin(2) - t4 * t25 + t72 * t95 - t19;
	t126 = t32 * t58;
	t134 = (-t34 * t102 * t61 / 0.2e1 - t106 * t110) * t33 * t70 + (t122 * t58 * t60 / 0.2e1 + t126 * t108 * t30 * pkin(2)) * t34 * t27 * t85;
	t137 = t126 * t60;
	t139 = atan2(-t63 * t61 / 0.2e1, t137 / 0.2e1);
	t140 = cos(t139);
	t141 = t140 * t32;
	t143 = sin(t139);
	t144 = t143 * t34;
	t145 = t27 * t58;
	t146 = t145 * t60;
	t148 = t141 * t61 - t144 * t146;
	t149 = 0.2e1 / t148;
	t151 = t1 ^ 2;
	t155 = 0.4e1 / t148 ^ 2;
	t160 = 0.1e1 / (0.1e1 + t151 * t36 * t64 * t108 * t155 / 0.4e1);
	t178 = t61 * t155 * t160;
	t240 = 0.1e1 / t151;
	t242 = t39 * t23;
	t247 = 0.1e1 / (0.4e1 * t35 * t240 * t242 * t107 + 0.1e1);
	t258 = t1 * t32;
	t264 = t242 * t107 * t247;
	unknown(1,1) = -t28 * t43;
	unknown(1,2) = 0.2e1 * t87;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.2e1 * t134;
	unknown(2,1) = -t106 * t60 * t149 * t160 / 0.2e1 - (t1 * t36 * t43 * t140 * t34 * t61 + t28 * t42 * t143 * t58 * t60 - t143 * t1 * t146) * t1 * t27 * t178 / 0.4e1;
	unknown(2,2) = (t1 * t56 * t61 - t28 * t65) * t149 * t160 / 0.2e1 - (-0.2e1 * t87 * t140 * t34 * t146 + t144 * t27 * t64 * t60 - t144 * t56 * t58 * t60 - 0.2e1 * t87 * t143 * t137 + t140 * t76 * t61 - t141 * t65) * t1 * t27 * t178 / 0.4e1;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = (t1 * t102 * t61 / 0.2e1 + t28 * t58 * t110) * t149 * t160 - (-t134 * t143 * t137 + t140 * t122 * t61 / 0.2e1 + t141 * t58 * t110 - t134 * t140 * t34 * t146 - t144 * t102 * t58 * t60 / 0.2e1 - t144 * t145 * t110) * t1 * t27 * t178 / 0.2e1;
	unknown(3,1) = 0.2e1 * t35 * t33 * t22 * t59 * t240 * t247 + 0.2e1 * t33 * t22 * t59 * t247;
	unknown(3,2) = -0.2e1 * (t1 * t76 * t61 - t258 * t65) * t34 * t240 * t264;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = -0.4e1 * (t1 * t122 * t61 / 0.2e1 + t258 * t58 * t110) * t34 * t240 * t264;
	Ja_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->15)
	unknown=NaN(3,5);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	Ja_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->15)
	unknown=NaN(3,5);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	Ja_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:19
	% EndTime: 2020-04-11 21:48:21
	% DurationCPUTime: 0.93s
	% Computational Cost: add. (55060->158), mult. (59112->429), div. (4135->16), fcn. (16323->13), ass. (0->165)
	unknown=NaN(3,5);
	t1 = cos(qJ(1));
	t2 = sin(qJ(5));
	t3 = t2 * pkin(1);
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * t3 * pkin(2);
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t11 * t15);
	t19 = cos(qJ(5));
	t20 = pkin(1) * t19;
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t20 * t25;
	t27 = t4 * t17 + t26;
	t28 = t1 * t27;
	t29 = 0.1e1 / t23;
	t30 = t28 * t29;
	t31 = -t6 + t7 + t21;
	t32 = 0.1e1 / t31;
	t33 = 0.1e1 / pkin(3);
	t34 = t32 * t33;
	t35 = t34 * t17;
	t37 = t19 * t17;
	t38 = t37 * pkin(1);
	t40 = t4 * t25 - t38;
	t41 = t1 * t40;
	t42 = t41 * t29;
	t43 = t6 - t7 - t21 + t23 + t24;
	t45 = t32 * t43 * t33;
	t47 = -t30 * t35 + t42 * t45;
	t48 = cos(qJ(3));
	t52 = -t30 * t45 - t42 * t35;
	t53 = sin(qJ(3));
	t55 = -t47 * t48 / 0.4e1 - t52 * t53 / 0.4e1;
	t56 = t40 * t29;
	t58 = t27 * t29;
	t60 = t56 * t35 + t58 * t45;
	t64 = -t58 * t35 + t56 * t45;
	t66 = t60 * t48 / 0.4e1 + t64 * t53 / 0.4e1;
	t67 = 0.1e1 / t66;
	t68 = t55 * t67;
	t69 = sin(qJ(1));
	t70 = t69 * t27;
	t71 = t70 * t29;
	t73 = t69 * t40;
	t74 = t73 * t29;
	t76 = -t71 * t35 + t74 * t45;
	t80 = -t74 * t35 - t71 * t45;
	t82 = -t76 * t48 / 0.4e1 - t80 * t53 / 0.4e1;
	t83 = t82 ^ 2;
	t84 = t66 ^ 2;
	t85 = 0.1e1 / t84;
	t88 = 0.1e1 / (t83 * t85 + 0.1e1);
	t90 = 0.1e1 / t17;
	t91 = t4 * t90;
	t96 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t15 - 0.2e1 * t11 * (-qJ(2) - pkin(6) + pkin(3));
	t100 = t91 * t96 / 0.2e1 + 0.2e1 * t20 * t22;
	t102 = t69 * t100 * t29;
	t106 = 0.1e1 / t23 / t22;
	t107 = t70 * t106;
	t111 = t34 * t90 * t96;
	t114 = t19 * t90;
	t119 = -t114 * pkin(1) * t96 / 0.2e1 + 0.2e1 * t4 * t22;
	t121 = t69 * t119 * t29;
	t124 = t73 * t106;
	t128 = 0.2e1 * t32 * t22 * t33;
	t147 = -(-t102 * t35 / 0.4e1 + t107 * t35 / 0.2e1 - t71 * t111 / 0.8e1 + t121 * t45 / 0.4e1 - t124 * t45 / 0.2e1 + t74 * t128 / 0.4e1) * t48 - (-t102 * t45 / 0.4e1 + t107 * t45 / 0.2e1 - t71 * t128 / 0.4e1 - t121 * t35 / 0.4e1 + t124 * t35 / 0.2e1 - t74 * t111 / 0.8e1) * t53;
	t150 = t119 * t29;
	t153 = t40 * t106;
	t156 = t56 * t32;
	t157 = t33 * t90;
	t158 = t157 * t96;
	t161 = t100 * t29;
	t164 = t27 * t106;
	t181 = t58 * t32;
	t186 = (t150 * t35 / 0.4e1 - t153 * t35 / 0.2e1 + t156 * t158 / 0.8e1 + t161 * t45 / 0.4e1 - t164 * t45 / 0.2e1 + t58 * t128 / 0.4e1) * t48 + (t150 * t45 / 0.4e1 - t153 * t45 / 0.2e1 + t56 * t128 / 0.4e1 - t161 * t35 / 0.4e1 + t164 * t35 / 0.2e1 - t181 * t158 / 0.8e1) * t53;
	t188 = t85 * t88;
	t190 = t147 * t67 * t88 - t186 * t82 * t188;
	t193 = t76 * t53 / 0.4e1 - t80 * t48 / 0.4e1;
	t198 = -t60 * t53 / 0.4e1 + t64 * t48 / 0.4e1;
	t201 = -t198 * t82 * t188 + t193 * t67 * t88;
	t205 = pkin(1) * pkin(2);
	t207 = t20 * pkin(2) * t15 + t11 * t19 * t205;
	t211 = t19 ^ 2;
	t215 = -0.2e1 * t7 * t211 * pkin(2) + t91 * t207 - t3 * t25 - t38;
	t217 = t69 * t215 * t29;
	t220 = t31 ^ 2;
	t221 = 0.1e1 / t220;
	t222 = t29 * t221;
	t223 = t70 * t222;
	t225 = t20 * pkin(2);
	t226 = t33 * t17 * t225;
	t230 = 0.2e1 * t34 * t90 * t207;
	t241 = -t114 * pkin(1) * t207 + t2 * t17 * pkin(1) - 0.2e1 * t4 * t19 * t205 - t26;
	t243 = t69 * t241 * t29;
	t246 = t73 * t222;
	t248 = t43 * t33 * t225;
	t251 = t29 * t32;
	t254 = t20 * pkin(2) * t33;
	t274 = -(-t217 * t35 / 0.4e1 - t223 * t226 / 0.2e1 - t71 * t230 / 0.8e1 + t243 * t45 / 0.4e1 + t246 * t248 / 0.2e1 + t73 * t251 * t254 / 0.2e1) * t48 - (-t217 * t45 / 0.4e1 - t223 * t248 / 0.2e1 - t70 * t251 * t254 / 0.2e1 - t243 * t35 / 0.4e1 - t246 * t226 / 0.2e1 - t74 * t230 / 0.8e1) * t53;
	t277 = t241 * t29;
	t280 = t221 * t33;
	t282 = t37 * t205;
	t285 = 0.2e1 * t157 * t207;
	t288 = t215 * t29;
	t291 = t221 * t43;
	t315 = (t277 * t35 / 0.4e1 + t56 * t280 * t282 / 0.2e1 + t156 * t285 / 0.8e1 + t288 * t45 / 0.4e1 + t58 * t291 * t254 / 0.2e1 + t181 * t254 / 0.2e1) * t48 + (t277 * t45 / 0.4e1 + t56 * t291 * t254 / 0.2e1 + t156 * t254 / 0.2e1 - t288 * t35 / 0.4e1 - t58 * t280 * t282 / 0.2e1 - t181 * t285 / 0.8e1) * t53;
	t318 = -t315 * t82 * t188 + t274 * t67 * t88;
	t322 = atan2(t82, t66);
	t323 = cos(t322);
	t325 = sin(t322);
	t327 = t323 * t66 + t325 * t82;
	t328 = 0.1e1 / t327;
	t330 = t55 ^ 2;
	t331 = t327 ^ 2;
	t332 = 0.1e1 / t331;
	t335 = 0.1e1 / (t330 * t332 + 0.1e1);
	t345 = t332 * t335;
	t349 = t1 * t100 * t29;
	t352 = t28 * t106;
	t358 = t1 * t119 * t29;
	t361 = t41 * t106;
	t366 = -t349 * t35 / 0.4e1 + t352 * t35 / 0.2e1 - t30 * t111 / 0.8e1 + t358 * t45 / 0.4e1 - t361 * t45 / 0.2e1 + t42 * t128 / 0.4e1;
	t380 = -t349 * t45 / 0.4e1 + t352 * t45 / 0.2e1 - t30 * t128 / 0.4e1 - t358 * t35 / 0.4e1 + t361 * t35 / 0.2e1 - t42 * t111 / 0.8e1;
	t397 = -t47 * t53 / 0.4e1 + t52 * t48 / 0.4e1;
	t411 = t1 * t215 * t29;
	t414 = t28 * t222;
	t420 = t1 * t241 * t29;
	t423 = t41 * t222;
	t429 = -t411 * t35 / 0.4e1 - t414 * t226 / 0.2e1 - t30 * t230 / 0.8e1 + t420 * t45 / 0.4e1 + t423 * t248 / 0.2e1 + t41 * t251 * t254 / 0.2e1;
	t444 = -t411 * t45 / 0.4e1 - t414 * t248 / 0.2e1 - t28 * t251 * t254 / 0.2e1 - t420 * t35 / 0.4e1 - t423 * t226 / 0.2e1 - t42 * t230 / 0.8e1;
	t461 = t76 * t53 / 0.4e1 - t80 * t48 / 0.4e1;
	t462 = sin(qJ(4));
	t464 = cos(qJ(4));
	t469 = -t397 * t464 - t69 * t462;
	t470 = 0.1e1 / t469;
	t474 = -t397 * t462 + t69 * t464;
	t475 = t474 ^ 2;
	t476 = t469 ^ 2;
	t477 = 0.1e1 / t476;
	t480 = 0.1e1 / (t475 * t477 + 0.1e1);
	t486 = t477 * t480;
	t491 = -t366 * t53 + t380 * t48;
	t493 = t470 * t480;
	t497 = t474 * t477 * t480;
	t510 = -t429 * t53 + t444 * t48;
	unknown(1,1) = t68 * t88;
	unknown(1,2) = t190;
	unknown(1,3) = t201;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = t318;
	unknown(2,1) = (-t76 * t48 / 0.4e1 - t80 * t53 / 0.4e1) * t328 * t335 + (t68 * t88 * t323 * t82 - t55 * t88 * t325 + t325 * t55) * t55 * t345;
	unknown(2,2) = (t366 * t48 + t380 * t53) * t328 * t335 + (t190 * t323 * t82 - t190 * t325 * t66 + t325 * t147 + t323 * t186) * t55 * t345;
	unknown(2,3) = t397 * t328 * t335 + (t201 * t323 * t82 - t201 * t325 * t66 + t325 * t193 + t323 * t198) * t55 * t345;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = (t429 * t48 + t444 * t53) * t328 * t335 + (t318 * t323 * t82 - t318 * t325 * t66 + t325 * t274 + t323 * t315) * t55 * t345;
	unknown(3,1) = (t1 * t464 - t461 * t462) * t470 * t480 - (-t1 * t462 - t461 * t464) * t474 * t486;
	unknown(3,2) = -t491 * t462 * t493 + t491 * t464 * t497;
	unknown(3,3) = -t55 * t462 * t493 + t55 * t464 * t497;
	unknown(3,4) = t474 ^ 2 * t486 + t480;
	unknown(3,5) = -t510 * t462 * t493 + t510 * t464 * t497;
	Ja_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->15)
	unknown=NaN(3,5);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	Ja_rot = unknown;
else
	Ja_rot=NaN(3,5);
end