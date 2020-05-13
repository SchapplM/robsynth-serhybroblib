% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh4m1DE1
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
%   Wie in palh4m1DE1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AD,CB,CE,EP,HC,OT,TA,TD]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 22:26
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = palh4m1DE1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1DE1_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1DE1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1DE1_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:49
	% EndTime: 2020-04-11 22:24:49
	% DurationCPUTime: 0.02s
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
	% StartTime: 2020-04-11 22:24:49
	% EndTime: 2020-04-11 22:24:49
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
	% StartTime: 2020-04-11 22:24:50
	% EndTime: 2020-04-11 22:24:50
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
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:50
	% EndTime: 2020-04-11 22:24:50
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (22821->95), mult. (25564->276), div. (1449->21), fcn. (7155->11), ass. (0->116)
	unknown=NaN(3,5);
	t1 = cos(qJ(1));
	t2 = cos(qJ(5));
	t3 = sin(qJ(5));
	t4 = pkin(1) * t3;
	t6 = 0.2e1 * pkin(2) * t4;
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t15 * t11);
	t19 = pkin(1) * t17 * t2;
	t20 = t4 - pkin(2);
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t27 = -t25 * t20 - t19;
	t28 = t27 * t1;
	t30 = t2 * pkin(1);
	t31 = t25 * t30;
	t32 = -t17 * t20 + t31;
	t33 = 0.1e1 / t32;
	t34 = sin(qJ(1));
	t35 = t34 ^ 2;
	t36 = t27 ^ 2;
	t38 = t32 ^ 2;
	t39 = 0.1e1 / t38;
	t42 = 0.1e1 / (t39 * t36 * t35 + 0.1e1);
	t43 = t42 * t33;
	t45 = 0.1e1 / t17;
	t46 = t45 * t2;
	t51 = -0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3)) - 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t56 = -t51 * pkin(1) * t46 / 0.2e1 - 0.2e1 * t22 * t20;
	t58 = 0.1e1 / t22;
	t59 = -t6 + t7 + t21;
	t60 = 0.1e1 / t59;
	t61 = t60 * t58;
	t62 = 0.1e1 / t23;
	t63 = t62 * t38;
	t64 = t59 ^ 2;
	t65 = 0.1e1 / t64;
	t67 = t62 * t36;
	t69 = t65 * t63 + t65 * t67;
	t70 = sqrt(t69);
	t71 = 0.1e1 / t70;
	t72 = t71 * t61;
	t74 = t27 * t34;
	t76 = t71 * t60 * t62;
	t78 = t58 * t74;
	t80 = 0.1e1 / t70 / t69;
	t81 = t80 * t60;
	t82 = t62 * t32;
	t83 = -t45 * t20;
	t87 = t51 * t83 / 0.2e1 + 0.2e1 * t22 * t30;
	t91 = 0.1e1 / t23 / t22;
	t94 = t62 * t27;
	t99 = -t65 * t91 * t36 - t65 * t91 * t38 + t56 * t65 * t94 + t87 * t65 * t82;
	t100 = 0.2e1 * t99 * t81;
	t106 = t70 * t59;
	t107 = t42 * t106;
	t110 = t71 * t60;
	t113 = t58 * t32;
	t118 = t22 * t27;
	t121 = t42 * t39 * t106;
	t123 = t107 * t22 * t33 * (-t72 * t56 * t34 + t76 * t74 + t100 * t78 / 0.2e1) + t121 * t118 * t34 * (t110 * t58 * t87 - t110 * t82 - t100 * t113 / 0.2e1);
	t129 = pkin(1) * pkin(2);
	t131 = t15 * pkin(2) * t30 + t129 * t2 * t11;
	t138 = -t131 * pkin(1) * t46 + pkin(1) * t17 * t3 + 0.2e1 * t129 * t2 * t20 - t31;
	t141 = t65 * t58;
	t144 = t129 * t2 * t71;
	t150 = t2 ^ 2;
	t154 = -0.2e1 * pkin(2) * t150 * t7 + t131 * t83 - t25 * t4 - t19;
	t159 = 0.1e1 / t64 / t59;
	t161 = pkin(2) * t30;
	t170 = 0.2e1 * t138 * t65 * t94 + 0.2e1 * t154 * t65 * t82 + 0.4e1 * t161 * t159 * t63 + 0.4e1 * t161 * t159 * t67;
	t171 = t170 * t81;
	t189 = t107 * t22 * t33 * (-t72 * t138 * t34 - 0.2e1 * t144 * t141 * t74 + t171 * t78 / 0.2e1) + t121 * t118 * t34 * (t110 * t58 * t154 + 0.2e1 * t144 * t65 * t113 - t171 * t113 / 0.2e1);
	t192 = atan2(-t72 * t74, t110 * t113);
	t193 = cos(t192);
	t194 = t32 * t193;
	t196 = sin(t192);
	t197 = t34 * t196;
	t198 = t27 * t197;
	t200 = t72 * t194 - t72 * t198;
	t201 = 0.1e1 / t200;
	t202 = t1 ^ 2;
	t207 = t200 ^ 2;
	t208 = 0.1e1 / t207;
	t212 = 0.1e1 / (0.1e1 + t208 / t69 * t65 * t62 * t36 * t202);
	t230 = t58 * t27;
	t233 = t212 * t208 * t110;
	t239 = t58 * t28;
	t251 = t58 * t194;
	t256 = t110 * t230;
	t310 = 0.1e1 / t202;
	t313 = t64 * t23;
	t317 = 0.1e1 / (t69 * t313 * t39 * t310 * t35 + 0.1e1);
	t328 = t32 * t1;
	t330 = t58 * t328;
	t335 = t39 * t310;
	t338 = t317 * t69 * t313;
	unknown(1,1) = -t43 * t28;
	unknown(1,2) = t123;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = t189;
	unknown(2,1) = -t212 * t201 * t110 * t78 - t233 * t230 * t1 * (t72 * t34 * t193 * t43 * t36 * t1 + t110 * t58 * t196 * t42 * t28 - t72 * t27 * t1 * t196);
	unknown(2,2) = t212 * t201 * (t72 * t56 * t1 - t76 * t28 - t100 * t239 / 0.2e1) - t233 * t230 * t1 * (-t72 * t32 * t196 * t123 + t72 * t87 * t193 - t76 * t194 - t100 * t251 / 0.2e1 - t256 * t34 * t193 * t123 - t72 * t56 * t197 + t76 * t198 + t99 * t80 * t61 * t198);
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = t212 * t201 * (t72 * t138 * t1 + 0.2e1 * t144 * t141 * t28 - t171 * t239 / 0.2e1) - t233 * t230 * t1 * (-t72 * t32 * t196 * t189 + t72 * t154 * t193 + 0.2e1 * t144 * t141 * t194 - t171 * t251 / 0.2e1 - t256 * t34 * t193 * t189 - t72 * t138 * t197 - 0.2e1 * t161 * t71 * t65 * t230 * t197 + t170 * t80 * t61 * t198 / 0.2e1);
	unknown(3,1) = t317 * t310 * t106 * t22 * t33 * t35 + t317 * t106 * t22 * t33;
	unknown(3,2) = -t338 * t335 * t34 * (t72 * t87 * t1 - t76 * t328 - t100 * t330 / 0.2e1);
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = -t338 * t335 * t34 * (t72 * t154 * t1 + 0.2e1 * t144 * t141 * t328 - t171 * t330 / 0.2e1);
	Ja_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:50
	% EndTime: 2020-04-11 22:24:50
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
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:51
	% EndTime: 2020-04-11 22:24:51
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
	% StartTime: 2020-04-11 22:24:52
	% EndTime: 2020-04-11 22:24:55
	% DurationCPUTime: 2.66s
	% Computational Cost: add. (203952->220), mult. (225038->590), div. (18199->25), fcn. (58827->17), ass. (0->227)
	unknown=NaN(3,5);
	t1 = cos(qJ(1));
	t2 = sin(qJ(5));
	t3 = pkin(1) * t2;
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * pkin(2) * t3;
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t16 = t15 * t11;
	t17 = sqrt(-t16);
	t19 = cos(qJ(5));
	t20 = t19 * pkin(1);
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t25 * t20;
	t27 = t17 * t4 + t26;
	t28 = t27 * t1;
	t29 = 0.1e1 / t23;
	t30 = -t6 + t7 + t21;
	t31 = 0.1e1 / t30;
	t32 = t31 * t29;
	t33 = t32 * t28;
	t34 = t27 ^ 2;
	t35 = t29 * t34;
	t36 = t30 ^ 2;
	t37 = 0.1e1 / t36;
	t40 = pkin(1) * t17 * t19;
	t42 = t25 * t4 - t40;
	t43 = t42 ^ 2;
	t44 = t29 * t43;
	t46 = t37 * t35 + t37 * t44;
	t47 = sqrt(t46);
	t48 = 0.1e1 / t47;
	t49 = 0.1e1 / pkin(3);
	t50 = t49 * t48;
	t51 = 0.1e1 / t24;
	t52 = t51 * t29;
	t54 = t6 - t7 - t21 + t23 + t24;
	t55 = t54 ^ 2;
	t58 = t51 * t29 * t55 - t16 * t52;
	t59 = sqrt(t58);
	t60 = 0.1e1 / t59;
	t61 = t60 * t17;
	t62 = t61 * t50;
	t64 = t42 * t1;
	t65 = t32 * t64;
	t66 = t54 * t48;
	t67 = t60 * t49;
	t68 = t67 * t66;
	t70 = -t62 * t33 + t68 * t65;
	t71 = cos(qJ(3));
	t75 = -t68 * t33 - t62 * t65;
	t76 = sin(qJ(3));
	t78 = -t71 * t70 - t76 * t75;
	t79 = t29 * t42;
	t80 = t31 * t79;
	t82 = t29 * t27;
	t83 = t31 * t82;
	t85 = t62 * t80 + t68 * t83;
	t89 = -t62 * t83 + t68 * t80;
	t91 = t71 * t85 + t76 * t89;
	t92 = 0.1e1 / t91;
	t93 = t92 * t78;
	t94 = sin(qJ(1));
	t95 = t27 * t94;
	t96 = t32 * t95;
	t98 = t42 * t94;
	t99 = t32 * t98;
	t101 = -t62 * t96 + t68 * t99;
	t105 = -t62 * t99 - t68 * t96;
	t107 = -t71 * t101 - t76 * t105;
	t108 = t107 ^ 2;
	t109 = t91 ^ 2;
	t110 = 0.1e1 / t109;
	t113 = 0.1e1 / (t110 * t108 + 0.1e1);
	t115 = 0.1e1 / t17;
	t116 = t115 * t4;
	t118 = 0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3));
	t120 = 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t121 = -t118 - t120;
	t125 = t121 * t116 / 0.2e1 + 0.2e1 * t22 * t20;
	t127 = t32 * t125 * t94;
	t130 = 0.1e1 / t23 / t22;
	t131 = t31 * t130;
	t132 = t131 * t95;
	t136 = 0.1e1 / t47 / t46;
	t137 = t49 * t136;
	t142 = t115 * t19;
	t147 = -t121 * pkin(1) * t142 / 0.2e1 + 0.2e1 * t22 * t4;
	t152 = t125 * t37 * t82 - t37 * t130 * t34 - t37 * t130 * t43 + t147 * t37 * t79;
	t154 = 0.2e1 * t152 * t61 * t137;
	t157 = t60 * t115;
	t159 = t121 * t157 * t50;
	t163 = 0.1e1 / t59 / t58;
	t164 = t163 * t17;
	t170 = t29 * t54;
	t177 = 0.2e1 * t16 * t51 * t130 - 0.2e1 * t51 * t130 * t55 + 0.4e1 * t22 * t51 * t170 - t118 * t52 - t120 * t52;
	t179 = t177 * t164 * t50;
	t183 = t32 * t147 * t94;
	t185 = t131 * t98;
	t188 = t54 * t136;
	t190 = 0.2e1 * t152 * t67 * t188;
	t194 = 0.2e1 * t67 * t22 * t48;
	t196 = t163 * t49;
	t198 = t177 * t196 * t66;
	t222 = -t71 * (-t62 * t127 + 0.2e1 * t62 * t132 + t154 * t96 / 0.2e1 - t159 * t96 / 0.2e1 + t179 * t96 / 0.2e1 + t68 * t183 - 0.2e1 * t68 * t185 - t190 * t99 / 0.2e1 + t194 * t99 - t198 * t99 / 0.2e1) - t76 * (-t68 * t127 + 0.2e1 * t68 * t132 + t190 * t96 / 0.2e1 - t194 * t96 + t198 * t96 / 0.2e1 - t62 * t183 + 0.2e1 * t62 * t185 + t154 * t99 / 0.2e1 - t159 * t99 / 0.2e1 + t179 * t99 / 0.2e1);
	t226 = t31 * t29 * t147;
	t229 = t31 * t130 * t42;
	t232 = t136 * t31;
	t233 = t232 * t79;
	t234 = t17 * t49;
	t235 = 0.2e1 * t152 * t60;
	t236 = t235 * t234;
	t239 = t48 * t31;
	t240 = t239 * t79;
	t241 = t115 * t49;
	t243 = t121 * t60 * t241;
	t246 = t177 * t163;
	t247 = t246 * t234;
	t251 = t31 * t29 * t125;
	t254 = t31 * t130 * t27;
	t257 = t232 * t82;
	t258 = t49 * t54;
	t259 = t235 * t258;
	t263 = t239 * t82;
	t264 = t246 * t258;
	t288 = t71 * (t62 * t226 - 0.2e1 * t62 * t229 - t236 * t233 / 0.2e1 + t243 * t240 / 0.2e1 - t247 * t240 / 0.2e1 + t68 * t251 - 0.2e1 * t68 * t254 - t259 * t257 / 0.2e1 + t194 * t83 - t264 * t263 / 0.2e1) + t76 * (t68 * t226 - 0.2e1 * t68 * t229 - t259 * t233 / 0.2e1 + t194 * t80 - t264 * t240 / 0.2e1 - t62 * t251 + 0.2e1 * t62 * t254 + t236 * t257 / 0.2e1 - t243 * t263 / 0.2e1 + t247 * t263 / 0.2e1);
	t290 = t113 * t110;
	t292 = -t290 * t107 * t288 + t113 * t92 * t222;
	t295 = t76 * t101 - t71 * t105;
	t300 = t71 * t89 - t76 * t85;
	t303 = -t290 * t107 * t300 + t113 * t92 * t295;
	t307 = pkin(1) * pkin(2);
	t309 = t15 * pkin(2) * t20 + t307 * t19 * t11;
	t313 = t19 ^ 2;
	t317 = -0.2e1 * pkin(2) * t313 * t7 + t309 * t116 - t25 * t3 - t40;
	t319 = t32 * t317 * t94;
	t322 = t48 * t37 * t29;
	t323 = t322 * t95;
	t325 = pkin(2) * t20;
	t326 = t325 * t60 * t234;
	t333 = 0.1e1 / t36 / t30;
	t345 = -t309 * pkin(1) * t142 + pkin(1) * t17 * t2 - 0.2e1 * t307 * t19 * t4 - t26;
	t352 = 0.2e1 * t317 * t37 * t82 + 0.4e1 * t325 * t333 * t35 + 0.4e1 * t325 * t333 * t44 + 0.2e1 * t345 * t37 * t79;
	t354 = t352 * t61 * t137;
	t358 = 0.2e1 * t309 * t157 * t50;
	t371 = 0.2e1 * t15 * t307 * t19 * t52 + 0.2e1 * t325 * t11 * t52 + 0.4e1 * t325 * t51 * t170;
	t373 = t371 * t164 * t50;
	t377 = t32 * t345 * t94;
	t379 = t322 * t98;
	t381 = t325 * t60 * t258;
	t385 = t352 * t67 * t188;
	t388 = t48 * t32;
	t392 = t60 * t49 * pkin(2) * t20;
	t396 = t371 * t196 * t66;
	t422 = -t71 * (-t62 * t319 - 0.2e1 * t326 * t323 + t354 * t96 / 0.2e1 - t358 * t96 / 0.2e1 + t373 * t96 / 0.2e1 + t68 * t377 + 0.2e1 * t381 * t379 - t385 * t99 / 0.2e1 + 0.2e1 * t392 * t388 * t98 - t396 * t99 / 0.2e1) - t76 * (-t68 * t319 - 0.2e1 * t381 * t323 + t385 * t96 / 0.2e1 - 0.2e1 * t392 * t388 * t95 + t396 * t96 / 0.2e1 - t62 * t377 - 0.2e1 * t326 * t379 + t354 * t99 / 0.2e1 - t358 * t99 / 0.2e1 + t373 * t99 / 0.2e1);
	t426 = t31 * t29 * t345;
	t428 = t48 * t37;
	t429 = t49 * t428;
	t431 = t325 * t61;
	t434 = t352 * t60;
	t435 = t434 * t234;
	t439 = 0.2e1 * t309 * t60 * t241;
	t442 = t371 * t163;
	t443 = t442 * t234;
	t447 = t31 * t29 * t317;
	t449 = t54 * t428;
	t453 = t434 * t258;
	t458 = t442 * t258;
	t485 = t71 * (t62 * t426 + 0.2e1 * t431 * t429 * t79 - t435 * t233 / 0.2e1 + t439 * t240 / 0.2e1 - t443 * t240 / 0.2e1 + t68 * t447 + 0.2e1 * t392 * t449 * t82 - t453 * t257 / 0.2e1 + 0.2e1 * t392 * t263 - t458 * t263 / 0.2e1) + t76 * (t68 * t426 + 0.2e1 * t392 * t449 * t79 - t453 * t233 / 0.2e1 + 0.2e1 * t392 * t240 - t458 * t240 / 0.2e1 - t62 * t447 - 0.2e1 * t431 * t429 * t82 + t435 * t257 / 0.2e1 - t439 * t263 / 0.2e1 + t443 * t263 / 0.2e1);
	t488 = -t290 * t107 * t485 + t113 * t92 * t422;
	t492 = atan2(t107, t91);
	t493 = cos(t492);
	t495 = sin(t492);
	t497 = t107 * t495 + t91 * t493;
	t498 = 0.1e1 / t497;
	t500 = t78 ^ 2;
	t501 = t497 ^ 2;
	t502 = 0.1e1 / t501;
	t505 = 0.1e1 / (t502 * t500 + 0.1e1);
	t515 = t505 * t502;
	t519 = t32 * t125 * t1;
	t521 = t131 * t28;
	t531 = t32 * t147 * t1;
	t533 = t131 * t64;
	t541 = -t62 * t519 + 0.2e1 * t62 * t521 + t154 * t33 / 0.2e1 - t159 * t33 / 0.2e1 + t179 * t33 / 0.2e1 + t68 * t531 - 0.2e1 * t68 * t533 - t190 * t65 / 0.2e1 + t194 * t65 - t198 * t65 / 0.2e1;
	t560 = -t68 * t519 + 0.2e1 * t68 * t521 + t190 * t33 / 0.2e1 - t194 * t33 + t198 * t33 / 0.2e1 - t62 * t531 + 0.2e1 * t62 * t533 + t154 * t65 / 0.2e1 - t159 * t65 / 0.2e1 + t179 * t65 / 0.2e1;
	t577 = -t76 * t70 + t71 * t75;
	t591 = t32 * t317 * t1;
	t593 = t322 * t28;
	t603 = t32 * t345 * t1;
	t605 = t322 * t64;
	t615 = -t62 * t591 - 0.2e1 * t326 * t593 + t354 * t33 / 0.2e1 - t358 * t33 / 0.2e1 + t373 * t33 / 0.2e1 + t68 * t603 + 0.2e1 * t381 * t605 - t385 * t65 / 0.2e1 + 0.2e1 * t392 * t388 * t64 - t396 * t65 / 0.2e1;
	t636 = -t68 * t591 - 0.2e1 * t381 * t593 + t385 * t33 / 0.2e1 - 0.2e1 * t392 * t388 * t28 + t396 * t33 / 0.2e1 - t62 * t603 - 0.2e1 * t326 * t605 + t354 * t65 / 0.2e1 - t358 * t65 / 0.2e1 + t373 * t65 / 0.2e1;
	t653 = t76 * t101 - t71 * t105;
	t654 = sin(qJ(4));
	t656 = cos(qJ(4));
	t661 = -t656 * t577 - t654 * t94;
	t662 = 0.1e1 / t661;
	t666 = -t654 * t577 + t656 * t94;
	t667 = t666 ^ 2;
	t668 = t661 ^ 2;
	t669 = 0.1e1 / t668;
	t672 = 0.1e1 / (t669 * t667 + 0.1e1);
	t678 = t672 * t669;
	t683 = -t76 * t541 + t71 * t560;
	t685 = t672 * t662;
	t689 = t672 * t669 * t666;
	t702 = -t76 * t615 + t71 * t636;
	unknown(1,1) = t113 * t93;
	unknown(1,2) = t292;
	unknown(1,3) = t303;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = t488;
	unknown(2,1) = t505 * t498 * (-t71 * t101 - t76 * t105) + t515 * t78 * (t107 * t493 * t113 * t93 - t495 * t113 * t78 + t78 * t495);
	unknown(2,2) = t505 * t498 * (t71 * t541 + t76 * t560) + t515 * t78 * (t107 * t493 * t292 - t91 * t495 * t292 + t222 * t495 + t288 * t493);
	unknown(2,3) = t505 * t498 * t577 + t515 * t78 * (t107 * t493 * t303 - t91 * t495 * t303 + t295 * t495 + t300 * t493);
	unknown(2,4) = 0.0e0;
	unknown(2,5) = t505 * t498 * (t71 * t615 + t76 * t636) + t515 * t78 * (t107 * t493 * t488 - t91 * t495 * t488 + t422 * t495 + t485 * t493);
	unknown(3,1) = t672 * t662 * (t656 * t1 - t654 * t653) - t678 * t666 * (-t654 * t1 - t656 * t653);
	unknown(3,2) = -t685 * t654 * t683 + t689 * t656 * t683;
	unknown(3,3) = -t685 * t654 * t78 + t689 * t656 * t78;
	unknown(3,4) = t678 * t666 ^ 2 + t672;
	unknown(3,5) = -t685 * t654 * t702 + t689 * t656 * t702;
	Ja_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:49
	% EndTime: 2020-04-11 22:24:49
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
else
	Ja_rot=NaN(3,5);
end