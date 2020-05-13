% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh4m1DE2
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
%   Wie in palh4m1DE2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AD,CB,CE,EP,HC,OT,TA,TD]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 22:54
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = palh4m1DE2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1DE2_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1DE2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1DE2_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:32
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
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:32
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
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:32
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
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:33
	% DurationCPUTime: 0.42s
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
	% StartTime: 2020-04-11 22:53:33
	% EndTime: 2020-04-11 22:53:33
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
	% StartTime: 2020-04-11 22:53:33
	% EndTime: 2020-04-11 22:53:33
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
	% StartTime: 2020-04-11 22:53:33
	% EndTime: 2020-04-11 22:53:35
	% DurationCPUTime: 1.25s
	% Computational Cost: add. (92870->151), mult. (97320->389), div. (8047->25), fcn. (28018->16), ass. (0->167)
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
	t29 = 0.1e1 / t22;
	t30 = t29 * t28;
	t31 = -t6 + t7 + t21;
	t32 = 0.1e1 / t31;
	t33 = t27 ^ 2;
	t34 = 0.1e1 / t23;
	t35 = t34 * t33;
	t36 = t31 ^ 2;
	t37 = 0.1e1 / t36;
	t40 = pkin(1) * t17 * t19;
	t42 = t25 * t4 - t40;
	t43 = t42 ^ 2;
	t44 = t34 * t43;
	t46 = t37 * t35 + t37 * t44;
	t47 = sqrt(t46);
	t48 = 0.1e1 / t47;
	t49 = t48 * t32;
	t50 = 0.1e1 / pkin(3);
	t51 = t50 * t29;
	t53 = t6 - t7 - t21 + t23 + t24;
	t56 = atan2(t17 * t51, t50 * t29 * t53);
	t57 = t56 + qJ(3);
	t58 = sin(t57);
	t59 = t58 * t49;
	t61 = t42 * t1;
	t62 = t29 * t61;
	t63 = cos(t57);
	t64 = t63 * t49;
	t66 = t59 * t30 - t64 * t62;
	t67 = t29 * t42;
	t69 = t29 * t27;
	t71 = t59 * t67 + t64 * t69;
	t72 = 0.1e1 / t71;
	t73 = t72 * t66;
	t74 = sin(qJ(1));
	t75 = t27 * t74;
	t76 = t29 * t75;
	t78 = t42 * t74;
	t79 = t29 * t78;
	t81 = t59 * t76 - t64 * t79;
	t82 = t81 ^ 2;
	t83 = t71 ^ 2;
	t84 = 0.1e1 / t83;
	t87 = 0.1e1 / (t84 * t82 + 0.1e1);
	t89 = 0.1e1 / t17;
	t90 = t89 * t4;
	t95 = -0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3)) - 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t99 = t95 * t90 / 0.2e1 + 0.2e1 * t22 * t20;
	t106 = 0.1e1 / t47 / t46;
	t107 = t106 * t32;
	t108 = t34 * t27;
	t112 = 0.1e1 / t23 / t22;
	t115 = t34 * t42;
	t116 = t89 * t19;
	t121 = -t95 * pkin(1) * t116 / 0.2e1 + 0.2e1 * t22 * t4;
	t126 = t99 * t37 * t108 - t37 * t112 * t33 - t37 * t112 * t43 + t121 * t37 * t115;
	t128 = 0.2e1 * t126 * t58 * t107;
	t137 = 0.1e1 / t53;
	t140 = t53 ^ 2;
	t141 = 0.1e1 / t140;
	t144 = 0.1e1 / (-t141 * t16 + 0.1e1);
	t155 = t144 * t141 * t17;
	t157 = t144 * pkin(3) * t22 * t137 * (-t17 * t50 * t34 + t95 * t89 * t51 / 0.2e1) - t155 * pkin(3) * t22 * (0.2e1 * t50 * t29 * t22 - t50 * t34 * t53);
	t159 = t63 * t157 * t49;
	t167 = 0.2e1 * t126 * t63 * t107;
	t171 = t58 * t157 * t49;
	t173 = t59 * t29 * t99 * t74 - t59 * t34 * t75 - t128 * t76 / 0.2e1 + t159 * t76 - t64 * t29 * t121 * t74 + t64 * t34 * t78 + t167 * t79 / 0.2e1 + t171 * t79;
	t179 = t32 * t67;
	t180 = t58 * t106;
	t184 = t157 * t48;
	t190 = t32 * t69;
	t191 = t63 * t106;
	t197 = t59 * t29 * t121 - t126 * t180 * t179 - t126 * t191 * t190 + t63 * t184 * t179 - t58 * t184 * t190 + t64 * t29 * t99 - t64 * t108 - t59 * t115;
	t199 = t87 * t84;
	t201 = t87 * t72 * t173 - t199 * t81 * t197;
	t204 = t59 * t79 + t64 * t76;
	t209 = -t59 * t69 + t64 * t67;
	t212 = -t199 * t81 * t209 + t87 * t72 * t204;
	t216 = pkin(1) * pkin(2);
	t218 = t15 * pkin(2) * t20 + t216 * t19 * t11;
	t222 = t19 ^ 2;
	t226 = -0.2e1 * pkin(2) * t222 * t7 + t218 * t90 - t25 * t3 - t40;
	t230 = t37 * t29;
	t233 = pkin(2) * t20;
	t234 = t233 * t58 * t48;
	t241 = 0.1e1 / t36 / t31;
	t253 = -t218 * pkin(1) * t116 + pkin(1) * t17 * t2 - 0.2e1 * t216 * t19 * t4 - t26;
	t260 = 0.2e1 * t226 * t37 * t108 + 0.2e1 * t253 * t37 * t115 + 0.4e1 * t233 * t241 * t35 + 0.4e1 * t233 * t241 * t44;
	t262 = t260 * t58 * t107;
	t271 = t144 * t137 * t218 * t89 - 0.2e1 * t155 * t233;
	t273 = t63 * t271 * t49;
	t280 = t233 * t63 * t48;
	t284 = t260 * t63 * t107;
	t288 = t58 * t271 * t49;
	t290 = t59 * t29 * t226 * t74 + 0.2e1 * t234 * t230 * t75 - t262 * t76 / 0.2e1 + t273 * t76 - t64 * t29 * t253 * t74 - 0.2e1 * t280 * t230 * t78 + t284 * t79 / 0.2e1 + t288 * t79;
	t295 = t48 * t37;
	t304 = t271 * t48;
	t319 = t59 * t29 * t253 + 0.2e1 * t216 * t19 * t58 * t295 * t67 - t260 * t180 * t179 / 0.2e1 + t63 * t304 * t179 + t64 * t29 * t226 + 0.2e1 * t216 * t19 * t63 * t295 * t69 - t260 * t191 * t190 / 0.2e1 - t58 * t304 * t190;
	t322 = -t199 * t81 * t319 + t87 * t72 * t290;
	t323 = atan2(t81, t71);
	t324 = cos(t323);
	t326 = sin(t323);
	t328 = t71 * t324 + t81 * t326;
	t329 = 0.1e1 / t328;
	t331 = t66 ^ 2;
	t332 = t328 ^ 2;
	t333 = 0.1e1 / t332;
	t336 = 0.1e1 / (t333 * t331 + 0.1e1);
	t346 = t336 * t333;
	t350 = t29 * t99 * t1;
	t352 = t34 * t28;
	t358 = t29 * t121 * t1;
	t360 = t34 * t61;
	t380 = -t64 * t30 - t59 * t62;
	t394 = t29 * t226 * t1;
	t396 = t230 * t28;
	t403 = t29 * t253 * t1;
	t405 = t230 * t61;
	t424 = sin(qJ(4));
	t426 = cos(qJ(4));
	t431 = -t426 * t380 - t424 * t74;
	t432 = 0.1e1 / t431;
	t436 = -t424 * t380 + t426 * t74;
	t437 = t436 ^ 2;
	t438 = t431 ^ 2;
	t439 = 0.1e1 / t438;
	t442 = 0.1e1 / (t439 * t437 + 0.1e1);
	t448 = t442 * t439;
	t461 = -t64 * t350 + t64 * t352 + t167 * t30 / 0.2e1 + t171 * t30 - t59 * t358 + t59 * t360 + t128 * t62 / 0.2e1 - t159 * t62;
	t463 = t442 * t432;
	t467 = t442 * t439 * t436;
	t490 = -t64 * t394 - 0.2e1 * t280 * t396 + t284 * t30 / 0.2e1 + t288 * t30 - t59 * t403 - 0.2e1 * t234 * t405 + t262 * t62 / 0.2e1 - t273 * t62;
	unknown(1,1) = t87 * t73;
	unknown(1,2) = t201;
	unknown(1,3) = t212;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = t322;
	unknown(2,1) = t336 * t329 * t81 + t346 * t66 * (t81 * t324 * t87 * t73 - t326 * t87 * t66 + t66 * t326);
	unknown(2,2) = t336 * t329 * (-t59 * t350 + t59 * t352 + t128 * t30 / 0.2e1 - t159 * t30 + t64 * t358 - t64 * t360 - t167 * t62 / 0.2e1 - t171 * t62) + t346 * t66 * (t81 * t324 * t201 - t71 * t326 * t201 + t173 * t326 + t197 * t324);
	unknown(2,3) = t336 * t329 * t380 + t346 * t66 * (t81 * t324 * t212 - t71 * t326 * t212 + t204 * t326 + t209 * t324);
	unknown(2,4) = 0.0e0;
	unknown(2,5) = t336 * t329 * (-t59 * t394 - 0.2e1 * t234 * t396 + t262 * t30 / 0.2e1 - t273 * t30 + t64 * t403 + 0.2e1 * t280 * t405 - t284 * t62 / 0.2e1 - t288 * t62) + t346 * t66 * (t81 * t324 * t322 - t71 * t326 * t322 + t290 * t326 + t319 * t324);
	unknown(3,1) = t442 * t432 * (t426 * t1 - t424 * t204) - t448 * t436 * (-t424 * t1 - t426 * t204);
	unknown(3,2) = -t463 * t424 * t461 + t467 * t426 * t461;
	unknown(3,3) = -t463 * t424 * t66 + t467 * t426 * t66;
	unknown(3,4) = t448 * t436 ^ 2 + t442;
	unknown(3,5) = -t463 * t424 * t490 + t467 * t426 * t490;
	Ja_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:32
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