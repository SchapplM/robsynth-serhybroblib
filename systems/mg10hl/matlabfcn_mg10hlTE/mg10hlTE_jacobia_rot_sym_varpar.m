% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% mg10hlTE
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in mg10hlTE_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [17x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,AE,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 12:29
% Revision: f5729c120b58d1b0137ade1aa9321d1ea2b3cc0a (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = mg10hlTE_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'mg10hlTE_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'mg10hlTE_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [17 1]), ...
  'mg10hlTE_jacobia_rot_sym_varpar: pkin has to be [17x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
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
	Ja_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:25:21
	% EndTime: 2020-04-11 12:25:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (6->3), div. (5->2), fcn. (6->2), ass. (0->24)
	unknown=NaN(3,6);
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
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = (t6 * t8 + t8);
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Ja_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:25:21
	% EndTime: 2020-04-11 12:25:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:25:21
	% EndTime: 2020-04-11 12:25:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:25:22
	% EndTime: 2020-04-11 12:25:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:25:24
	% EndTime: 2020-04-11 12:25:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:25:26
	% EndTime: 2020-04-11 12:25:35
	% DurationCPUTime: 3.33s
	% Computational Cost: add. (198541->183), mult. (232096->397), div. (22569->23), fcn. (129836->22), ass. (0->198)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = sin(qJ(2));
	t3 = t2 * t1;
	t4 = cos(pkin(15));
	t5 = cos(qJ(2));
	t7 = sin(pkin(15));
	t9 = -t7 * t2 + t4 * t5;
	t10 = t9 * pkin(2);
	t13 = t4 * t2 + t7 * t5;
	t16 = 0.2e1 * pkin(2) * pkin(1) * t13;
	t17 = pkin(1) ^ 2;
	t21 = -t16 + t17 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t25 = -t16 + t17 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t27 = sqrt(-t25 * t21);
	t28 = t27 * t10;
	t30 = -t13 * pkin(2) + pkin(1);
	t31 = pkin(2) ^ 2;
	t32 = pkin(5) ^ 2;
	t33 = pkin(4) ^ 2;
	t34 = -t16 + t17 + t31 - t32 + t33;
	t36 = t34 * t30 - t28;
	t37 = t36 * t4;
	t38 = 0.1e1 / pkin(4);
	t39 = -t16 + t17 + t31;
	t41 = 0.1e1 / t39 * t38;
	t44 = t34 * t10;
	t45 = t27 * t30 + t44;
	t46 = t45 * t7;
	t48 = -t41 * t37 + t41 * t46;
	t49 = t48 * t3 / 0.2e1;
	t50 = t5 * t1;
	t51 = t36 * t7;
	t53 = t45 * t4;
	t55 = -t41 * t51 - t41 * t53;
	t57 = -t49 - t55 * t50 / 0.2e1;
	t58 = cos(pkin(17));
	t59 = 1 / pkin(7);
	t60 = qJ(6) + pkin(8);
	t62 = 0.1e1 / t60 * t59;
	t63 = (pkin(6) - pkin(7) - pkin(8) - qJ(6));
	t64 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t65 = t64 * t63;
	t66 = (pkin(6) + pkin(7) - pkin(8) - qJ(6));
	t67 = (pkin(6) + pkin(7) + pkin(8) + qJ(6));
	t68 = t67 * t66;
	t70 = sqrt(-t68 * t65);
	t72 = (pkin(6) ^ 2);
	t73 = (pkin(7) ^ 2);
	t74 = (pkin(8) ^ 2);
	t76 = 2 * pkin(8) * qJ(6);
	t77 = qJ(6) ^ 2;
	t78 = t72 - t73 - t74 - t76 - t77;
	t80 = atan2(t70 * t62, t78 * t62);
	t81 = t80 + pkin(16);
	t82 = cos(t81);
	t84 = sin(pkin(17));
	t85 = sin(t81);
	t87 = -t82 * t58 - t85 * t84;
	t90 = t48 * t50 / 0.2e1;
	t91 = t55 * t3 / 0.2e1 - t90;
	t94 = t85 * t58 - t82 * t84;
	t96 = t87 * t57 + t94 * t91;
	t97 = 0.1e1 / pkin(6);
	t98 = t60 ^ 2;
	t99 = 1 / t98;
	t100 = (t99 * t97);
	t101 = t72 - t73 + t74 + t76 + t77;
	t102 = t59 * t101;
	t103 = t78 * t102;
	t105 = (t63 * t100);
	t106 = (t66 * t64);
	t108 = t59 * t67 * t106;
	t110 = t103 * t100 - t108 * t105;
	t111 = cos(pkin(16));
	t113 = t59 * t70;
	t114 = t78 * t113;
	t116 = t70 * t102;
	t118 = -t114 * t100 + t116 * t100;
	t119 = sin(pkin(16));
	t121 = t111 * t110 / 0.4e1 - t119 * t118 / 0.4e1;
	t125 = -t94 * t57 + t87 * t91;
	t128 = -t111 * t118 / 0.4e1 - t119 * t110 / 0.4e1;
	t130 = -t121 * t96 - t128 * t125;
	t131 = t48 * t5 / 0.2e1;
	t133 = t131 - t55 * t2 / 0.2e1;
	t136 = t48 * t2 / 0.2e1;
	t137 = -t55 * t5 / 0.2e1 - t136;
	t139 = t87 * t133 + t94 * t137;
	t143 = -t94 * t133 + t87 * t137;
	t145 = t121 * t139 + t128 * t143;
	t146 = 0.1e1 / t145;
	t147 = t146 * t130;
	t148 = sin(qJ(1));
	t149 = t2 * t148;
	t150 = t48 * t149 / 0.2e1;
	t151 = t5 * t148;
	t153 = -t150 - t55 * t151 / 0.2e1;
	t156 = t48 * t151 / 0.2e1;
	t157 = t55 * t149 / 0.2e1 - t156;
	t159 = t87 * t153 + t94 * t157;
	t163 = -t94 * t153 + t87 * t157;
	t165 = -t121 * t159 - t128 * t163;
	t166 = t165 ^ 2;
	t167 = t145 ^ 2;
	t168 = 0.1e1 / t167;
	t171 = 0.1e1 / (t168 * t166 + 0.1e1);
	t173 = -t13 * pkin(2);
	t175 = 0.1e1 / t27;
	t180 = pkin(1) * pkin(2);
	t182 = t25 * pkin(2) * pkin(1) * t9 + t180 * t9 * t21;
	t189 = -t182 * t175 * t10 - 0.2e1 * t180 * t9 * t30 - t27 * t173 - t44;
	t194 = t39 ^ 2;
	t197 = t180 * t9 / t194;
	t203 = t9 ^ 2;
	t207 = -0.2e1 * pkin(1) * t203 * t31 + t182 * t175 * t30 + t34 * t173 - t28;
	t213 = -t41 * t189 * t4 / 0.2e1 - t197 * t38 * t37 + t41 * t207 * t7 / 0.2e1 + t197 * t38 * t46;
	t226 = -t41 * t189 * t7 / 0.2e1 - t197 * t38 * t51 - t41 * t207 * t4 / 0.2e1 - t197 * t38 * t53;
	t228 = -t156 - t213 * t149 + t55 * t149 / 0.2e1 - t226 * t151;
	t233 = t55 * t151 / 0.2e1 + t226 * t149 + t150 - t213 * t151;
	t241 = -t121 * (t87 * t228 + t94 * t233) - t128 * (-t94 * t228 + t87 * t233);
	t247 = -t136 + t213 * t5 - t55 * t5 / 0.2e1 - t226 * t2;
	t252 = t55 * t2 / 0.2e1 - t226 * t5 - t131 - t213 * t2;
	t260 = t121 * (t87 * t247 + t94 * t252) + t128 * (-t94 * t247 + t87 * t252);
	t262 = t171 * t168;
	t264 = t171 * t146 * t241 - t262 * t165 * t260;
	t265 = t99 * t59;
	t267 = 0.1e1 / t70;
	t273 = -(t67 * t66 * t63) + (t67 * t106) - t66 * t65 + t67 * t65;
	t281 = t78 ^ 2;
	t282 = 1 / t281;
	t286 = 0.1e1 / (-t282 * t68 * t65 + 0.1e1);
	t297 = t286 / t78 * t60 * pkin(7) * (-t70 * t265 + t273 * t267 * t62 / 0.2e1) - t286 * t282 * t70 * t60 * pkin(7) * (-(t78 * t265) - 0.2e1 * t60 * t62);
	t298 = t297 * t58;
	t300 = t297 * t84;
	t302 = t85 * t298 - t82 * t300;
	t306 = t82 * t298 + t85 * t300;
	t312 = 0.1e1 / t98 / t60 * t97;
	t315 = 0.2e1 * t59 * t60;
	t326 = t59 * t68;
	t338 = -t103 * t312 / 0.2e1 + t78 * t315 * t100 / 0.4e1 - t60 * t102 * t100 / 0.2e1 + t108 * t63 * t312 / 0.2e1 + t326 * t64 * t100 / 0.4e1 - t326 * t105 / 0.4e1 + t59 * t67 * t64 * t105 / 0.4e1 - (t59 * t106 * t105) / 0.4e1;
	t360 = t114 * t312 / 0.2e1 - t273 * t78 * t59 * t267 * t100 / 0.8e1 + t60 * t113 * t100 / 0.2e1 - t116 * t312 / 0.2e1 + t70 * t315 * t100 / 0.4e1 + t273 * t267 * t59 * t101 * t100 / 0.8e1;
	t362 = t111 * t338 - t119 * t360;
	t370 = -t111 * t360 - t119 * t338;
	t372 = -t121 * (t302 * t153 + t306 * t157) - t362 * t159 - t128 * (-t306 * t153 + t302 * t157) - t370 * t163;
	t385 = t121 * (t302 * t133 + t306 * t137) + t362 * t139 + t128 * (-t306 * t133 + t302 * t137) + t370 * t143;
	t388 = t171 * t146 * t372 - t262 * t165 * t385;
	t391 = -t87 * t153 - t94 * t157;
	t395 = t94 * t153 - t87 * t157;
	t398 = atan2(t165, t145);
	t399 = cos(t398);
	t401 = sin(t398);
	t403 = t145 * t399 + t165 * t401;
	t404 = 0.1e1 / t403;
	t406 = t130 ^ 2;
	t407 = t403 ^ 2;
	t408 = 0.1e1 / t407;
	t411 = 0.1e1 / (t408 * t406 + 0.1e1);
	t421 = t411 * t408;
	t427 = -t90 - t213 * t3 + t55 * t3 / 0.2e1 - t226 * t50;
	t432 = t55 * t50 / 0.2e1 + t226 * t3 + t49 - t213 * t50;
	t434 = t87 * t427 + t94 * t432;
	t438 = -t94 * t427 + t87 * t432;
	t455 = t302 * t57 + t306 * t91;
	t460 = t302 * t91 - t306 * t57;
	t478 = t121 * t395 - t128 * t391;
	t479 = sin(qJ(3));
	t481 = cos(qJ(3));
	t486 = t121 * t125 - t128 * t96;
	t489 = t479 * t148 + t481 * t486;
	t490 = 0.1e1 / t489;
	t494 = -t481 * t148 + t479 * t486;
	t495 = t494 ^ 2;
	t496 = t489 ^ 2;
	t497 = 0.1e1 / t496;
	t500 = 0.1e1 / (t497 * t495 + 0.1e1);
	t506 = t500 * t497;
	t511 = t121 * t438 - t128 * t434;
	t513 = t500 * t490;
	t517 = t500 * t497 * t494;
	t527 = t121 * t460 + t362 * t125 - t128 * t455 - t370 * t96;
	unknown(1,1) = t171 * t147;
	unknown(1,2) = t264;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = t388;
	unknown(2,1) = t411 * t404 * (t121 * t391 + t128 * t395) + t421 * t130 * (t165 * t399 * t171 * t147 - t401 * t171 * t130 + t130 * t401);
	unknown(2,2) = t411 * t404 * (t121 * t434 + t128 * t438) + t421 * t130 * (-t145 * t401 * t264 + t165 * t399 * t264 + t241 * t401 + t260 * t399);
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = t411 * t404 * (t121 * t455 + t370 * t125 + t128 * t460 + t362 * t96) + t421 * t130 * (-t145 * t401 * t388 + t165 * t399 * t388 + t372 * t401 + t385 * t399);
	unknown(3,1) = t500 * t490 * (-t481 * t1 + t479 * t478) - t506 * t494 * (t479 * t1 + t481 * t478);
	unknown(3,2) = t513 * t479 * t511 - t517 * t481 * t511;
	unknown(3,3) = t506 * t494 ^ 2 + t500;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = t513 * t479 * t527 - t517 * t481 * t527;
	Ja_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:25:30
	% EndTime: 2020-04-11 12:25:43
	% DurationCPUTime: 4.89s
	% Computational Cost: add. (289770->201), mult. (340406->449), div. (33314->25), fcn. (191367->24), ass. (0->218)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = sin(qJ(2));
	t3 = t2 * t1;
	t4 = cos(pkin(15));
	t5 = cos(qJ(2));
	t7 = sin(pkin(15));
	t9 = -t2 * t7 + t4 * t5;
	t10 = t9 * pkin(2);
	t13 = t2 * t4 + t5 * t7;
	t16 = 0.2e1 * pkin(2) * pkin(1) * t13;
	t17 = pkin(1) ^ 2;
	t21 = -t16 + t17 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t25 = -t16 + t17 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t27 = sqrt(-t25 * t21);
	t28 = t27 * t10;
	t30 = -t13 * pkin(2) + pkin(1);
	t31 = pkin(2) ^ 2;
	t32 = pkin(5) ^ 2;
	t33 = pkin(4) ^ 2;
	t34 = -t16 + t17 + t31 - t32 + t33;
	t36 = t30 * t34 - t28;
	t37 = t36 * t4;
	t38 = 0.1e1 / pkin(4);
	t39 = -t16 + t17 + t31;
	t41 = 0.1e1 / t39 * t38;
	t44 = t34 * t10;
	t45 = t27 * t30 + t44;
	t46 = t45 * t7;
	t48 = -t37 * t41 + t41 * t46;
	t49 = t48 * t3 / 0.2e1;
	t50 = t5 * t1;
	t51 = t36 * t7;
	t53 = t45 * t4;
	t55 = -t41 * t51 - t41 * t53;
	t57 = -t49 - t55 * t50 / 0.2e1;
	t58 = cos(pkin(17));
	t59 = 0.1e1 / pkin(7);
	t60 = qJ(6) + pkin(8);
	t62 = 0.1e1 / t60 * t59;
	t63 = pkin(6) - pkin(7) - pkin(8) - qJ(6);
	t64 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t65 = t64 * t63;
	t66 = pkin(6) + pkin(7) - pkin(8) - qJ(6);
	t67 = pkin(6) + pkin(7) + pkin(8) + qJ(6);
	t68 = t67 * t66;
	t70 = sqrt(-t68 * t65);
	t72 = (pkin(6) ^ 2);
	t73 = (pkin(7) ^ 2);
	t74 = (pkin(8) ^ 2);
	t76 = 2 * pkin(8) * qJ(6);
	t77 = qJ(6) ^ 2;
	t78 = t72 - t73 - t74 - t76 - t77;
	t80 = atan2(t70 * t62, t78 * t62);
	t81 = t80 + pkin(16);
	t82 = cos(t81);
	t84 = sin(pkin(17));
	t85 = sin(t81);
	t87 = -t58 * t82 - t84 * t85;
	t90 = t48 * t50 / 0.2e1;
	t91 = t55 * t3 / 0.2e1 - t90;
	t94 = t58 * t85 - t82 * t84;
	t96 = t57 * t87 + t91 * t94;
	t97 = 0.1e1 / pkin(6);
	t98 = t60 ^ 2;
	t99 = 0.1e1 / t98;
	t100 = t99 * t97;
	t101 = t59 * t70;
	t102 = t78 * t101;
	t104 = t72 - t73 + t74 + t76 + t77;
	t105 = t59 * t104;
	t106 = t70 * t105;
	t108 = -t100 * t102 + t100 * t106;
	t109 = cos(pkin(16));
	t111 = t78 * t105;
	t113 = t63 * t100;
	t114 = t66 * t64;
	t116 = t59 * t67 * t114;
	t118 = t100 * t111 - t113 * t116;
	t119 = sin(pkin(16));
	t121 = t109 * t108 / 0.4e1 + t119 * t118 / 0.4e1;
	t125 = -t57 * t94 + t87 * t91;
	t128 = t109 * t118 / 0.4e1 - t119 * t108 / 0.4e1;
	t130 = t121 * t96 + t125 * t128;
	t131 = sin(qJ(3));
	t133 = sin(qJ(1));
	t134 = cos(qJ(3));
	t136 = t130 * t131 - t133 * t134;
	t137 = t48 * t5 / 0.2e1;
	t139 = t137 - t55 * t2 / 0.2e1;
	t142 = t48 * t2 / 0.2e1;
	t143 = -t55 * t5 / 0.2e1 - t142;
	t145 = t139 * t87 + t143 * t94;
	t149 = -t139 * t94 + t143 * t87;
	t151 = t121 * t145 + t128 * t149;
	t152 = 0.1e1 / t151;
	t153 = t152 * t136;
	t154 = 0.1e1 / t131;
	t155 = t2 * t133;
	t156 = t48 * t155 / 0.2e1;
	t157 = t5 * t133;
	t159 = -t156 - t55 * t157 / 0.2e1;
	t162 = t48 * t157 / 0.2e1;
	t163 = t55 * t155 / 0.2e1 - t162;
	t165 = t159 * t87 + t163 * t94;
	t169 = -t159 * t94 + t163 * t87;
	t171 = t121 * t165 + t128 * t169;
	t173 = t134 * t1;
	t174 = t131 * t171 + t173;
	t175 = t174 ^ 2;
	t176 = t151 ^ 2;
	t177 = 0.1e1 / t176;
	t179 = t131 ^ 2;
	t180 = 0.1e1 / t179;
	t183 = 0.1e1 / (t175 * t177 * t180 + 0.1e1);
	t184 = t183 * t154;
	t186 = -t13 * pkin(2);
	t188 = 0.1e1 / t27;
	t193 = pkin(1) * pkin(2);
	t195 = pkin(1) * pkin(2) * t25 * t9 + t193 * t21 * t9;
	t202 = -t10 * t188 * t195 - 0.2e1 * t193 * t30 * t9 - t186 * t27 - t44;
	t207 = t39 ^ 2;
	t210 = t193 * t9 / t207;
	t216 = t9 ^ 2;
	t220 = -0.2e1 * pkin(1) * t216 * t31 + t188 * t195 * t30 + t186 * t34 - t28;
	t226 = -t41 * t202 * t4 / 0.2e1 - t210 * t38 * t37 + t41 * t220 * t7 / 0.2e1 + t210 * t38 * t46;
	t239 = -t41 * t202 * t7 / 0.2e1 - t210 * t38 * t51 - t41 * t220 * t4 / 0.2e1 - t210 * t38 * t53;
	t241 = -t162 - t226 * t155 + t55 * t155 / 0.2e1 - t239 * t157;
	t246 = t55 * t157 / 0.2e1 + t239 * t155 + t156 - t226 * t157;
	t254 = t121 * (t241 * t87 + t246 * t94) + t128 * (-t241 * t94 + t246 * t87);
	t260 = -t142 + t226 * t5 - t55 * t5 / 0.2e1 - t239 * t2;
	t265 = t55 * t2 / 0.2e1 - t239 * t5 - t137 - t226 * t2;
	t273 = t121 * (t260 * t87 + t265 * t94) + t128 * (-t260 * t94 + t265 * t87);
	t276 = t183 * t177 * t174;
	t278 = -t152 * t183 * t254 + t154 * t273 * t276;
	t280 = t131 * t1;
	t281 = t134 * t171 - t280;
	t288 = t134 * t152 * t174 * t180 * t183 - t152 * t184 * t281;
	t289 = t99 * t59;
	t291 = 0.1e1 / t70;
	t297 = -t63 * t66 * t67 + t114 * t67 - t65 * t66 + t65 * t67;
	t305 = t78 ^ 2;
	t306 = 1 / t305;
	t310 = 0.1e1 / (-t306 * t65 * t68 + 0.1e1);
	t321 = t310 / t78 * t60 * pkin(7) * (-t70 * t289 + t297 * t291 * t62 / 0.2e1) - t310 * t306 * t70 * t60 * pkin(7) * (-t289 * t78 - 0.2e1 * t60 * t62);
	t322 = t321 * t58;
	t324 = t321 * t84;
	t326 = t322 * t85 - t324 * t82;
	t330 = t322 * t82 + t324 * t85;
	t336 = 0.1e1 / t98 / t60 * t97;
	t349 = 0.2e1 * t59 * t60;
	t358 = t102 * t336 / 0.2e1 - t297 * t78 * t59 * t291 * t100 / 0.8e1 + t60 * t101 * t100 / 0.2e1 - t106 * t336 / 0.2e1 + t70 * t349 * t100 / 0.4e1 + t297 * t291 * t59 * t104 * t100 / 0.8e1;
	t372 = t59 * t68;
	t384 = -t111 * t336 / 0.2e1 + t78 * t349 * t100 / 0.4e1 - t60 * t105 * t100 / 0.2e1 + t116 * t63 * t336 / 0.2e1 + t372 * t64 * t100 / 0.4e1 - t372 * t113 / 0.4e1 + t59 * t67 * t64 * t113 / 0.4e1 - t59 * t114 * t113 / 0.4e1;
	t386 = t109 * t358 + t119 * t384;
	t394 = t109 * t384 - t119 * t358;
	t396 = t121 * (t159 * t326 + t163 * t330) + t386 * t165 + t128 * (-t159 * t330 + t163 * t326) + t394 * t169;
	t409 = t121 * (t139 * t326 + t143 * t330) + t386 * t145 + t128 * (-t139 * t330 + t143 * t326) + t394 * t149;
	t412 = -t152 * t183 * t396 + t154 * t276 * t409;
	t415 = -t159 * t87 - t163 * t94;
	t419 = t159 * t94 - t163 * t87;
	t421 = t121 * t415 + t128 * t419;
	t424 = t131 * t151;
	t425 = atan2(t174, -t424);
	t426 = cos(t425);
	t427 = t151 * t426;
	t429 = sin(t425);
	t431 = -t131 * t427 + t174 * t429;
	t432 = 0.1e1 / t431;
	t434 = t136 ^ 2;
	t435 = t431 ^ 2;
	t436 = 0.1e1 / t435;
	t439 = 0.1e1 / (t434 * t436 + 0.1e1);
	t450 = t439 * t436;
	t456 = -t90 - t226 * t3 + t55 * t3 / 0.2e1 - t239 * t50;
	t461 = t55 * t50 / 0.2e1 + t239 * t3 + t49 - t226 * t50;
	t463 = t456 * t87 + t461 * t94;
	t467 = -t456 * t94 + t461 * t87;
	t469 = t121 * t463 + t128 * t467;
	t471 = t439 * t432;
	t487 = -t130 * t134 - t131 * t133;
	t502 = t326 * t57 + t330 * t91;
	t507 = t326 * t91 - t330 * t57;
	t510 = t121 * t502 + t125 * t394 + t128 * t507 + t386 * t96;
	t526 = t134 * t421 + t280;
	t527 = cos(qJ(4));
	t531 = -t121 * t419 + t128 * t415;
	t532 = sin(qJ(4));
	t538 = -t121 * t125 + t128 * t96;
	t540 = -t487 * t532 + t527 * t538;
	t541 = 0.1e1 / t540;
	t545 = t487 * t527 + t532 * t538;
	t546 = t545 ^ 2;
	t547 = t540 ^ 2;
	t548 = 0.1e1 / t547;
	t551 = 0.1e1 / (t546 * t548 + 0.1e1);
	t557 = t551 * t548;
	t560 = t134 * t469;
	t564 = -t121 * t467 + t128 * t463;
	t586 = t134 * t510;
	t592 = -t121 * t507 - t125 * t386 + t128 * t502 + t394 * t96;
	unknown(1,1) = -t184 * t153;
	unknown(1,2) = t278;
	unknown(1,3) = t288;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = t412;
	unknown(2,1) = t439 * t432 * (-t131 * t421 + t173) + t450 * t136 * (-t153 * t154 * t174 * t183 * t426 - t136 * t183 * t429 + t136 * t429);
	unknown(2,2) = -t471 * t131 * t469 + t450 * t136 * (t131 * t254 * t429 - t131 * t273 * t426 + t174 * t278 * t426 + t278 * t424 * t429);
	unknown(2,3) = t439 * t432 * t487 + t450 * t136 * (t174 * t288 * t426 + t288 * t424 * t429 - t134 * t427 + t281 * t429);
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = -t471 * t131 * t510 + t450 * t136 * (t131 * t396 * t429 - t131 * t409 * t426 + t174 * t412 * t426 + t412 * t424 * t429);
	unknown(3,1) = t551 * t541 * (-t526 * t527 + t531 * t532) - t557 * t545 * (t526 * t532 + t527 * t531);
	unknown(3,2) = t551 * t541 * (-t527 * t560 + t532 * t564) - t557 * t545 * (t527 * t564 + t532 * t560);
	unknown(3,3) = t136 * t532 * t545 * t548 * t551 + t136 * t527 * t541 * t551;
	unknown(3,4) = t545 ^ 2 * t557 + t551;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = t551 * t541 * (-t527 * t586 + t532 * t592) - t557 * t545 * (t527 * t592 + t532 * t586);
	Ja_rot = unknown;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:25:35
	% EndTime: 2020-04-11 12:25:58
	% DurationCPUTime: 10.89s
	% Computational Cost: add. (656504->233), mult. (773627->521), div. (75945->23), fcn. (435842->26), ass. (0->238)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = sin(qJ(2));
	t3 = t2 * t1;
	t4 = cos(pkin(15));
	t5 = cos(qJ(2));
	t7 = sin(pkin(15));
	t9 = -t2 * t7 + t4 * t5;
	t10 = t9 * pkin(2);
	t13 = t2 * t4 + t5 * t7;
	t16 = 0.2e1 * pkin(2) * pkin(1) * t13;
	t17 = pkin(1) ^ 2;
	t21 = -t16 + t17 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t25 = -t16 + t17 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t27 = sqrt(-t25 * t21);
	t28 = t27 * t10;
	t30 = -t13 * pkin(2) + pkin(1);
	t31 = pkin(2) ^ 2;
	t32 = pkin(5) ^ 2;
	t33 = pkin(4) ^ 2;
	t34 = -t16 + t17 + t31 - t32 + t33;
	t36 = t30 * t34 - t28;
	t37 = t36 * t4;
	t38 = 0.1e1 / pkin(4);
	t39 = -t16 + t17 + t31;
	t41 = 0.1e1 / t39 * t38;
	t44 = t34 * t10;
	t45 = t27 * t30 + t44;
	t46 = t45 * t7;
	t48 = -t37 * t41 + t41 * t46;
	t49 = t48 * t3 / 0.2e1;
	t50 = t5 * t1;
	t51 = t36 * t7;
	t53 = t45 * t4;
	t55 = -t41 * t51 - t41 * t53;
	t57 = -t49 - t55 * t50 / 0.2e1;
	t58 = cos(pkin(17));
	t59 = 0.1e1 / pkin(7);
	t60 = qJ(6) + pkin(8);
	t62 = 0.1e1 / t60 * t59;
	t63 = pkin(6) - pkin(7) - pkin(8) - qJ(6);
	t64 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t65 = t64 * t63;
	t66 = pkin(6) + pkin(7) - pkin(8) - qJ(6);
	t67 = pkin(6) + pkin(7) + pkin(8) + qJ(6);
	t68 = t67 * t66;
	t70 = sqrt(-t68 * t65);
	t72 = (pkin(6) ^ 2);
	t73 = (pkin(7) ^ 2);
	t74 = (pkin(8) ^ 2);
	t76 = 2 * pkin(8) * qJ(6);
	t77 = qJ(6) ^ 2;
	t78 = t72 - t73 - t74 - t76 - t77;
	t80 = atan2(t70 * t62, t78 * t62);
	t81 = t80 + pkin(16);
	t82 = cos(t81);
	t84 = sin(pkin(17));
	t85 = sin(t81);
	t87 = -t58 * t82 - t84 * t85;
	t90 = t48 * t50 / 0.2e1;
	t91 = t55 * t3 / 0.2e1 - t90;
	t94 = t58 * t85 - t82 * t84;
	t96 = t57 * t87 + t91 * t94;
	t97 = 0.1e1 / pkin(6);
	t98 = t60 ^ 2;
	t99 = 0.1e1 / t98;
	t100 = t99 * t97;
	t101 = t59 * t70;
	t102 = t78 * t101;
	t104 = t72 - t73 + t74 + t76 + t77;
	t105 = t59 * t104;
	t106 = t70 * t105;
	t108 = -t100 * t102 + t100 * t106;
	t109 = cos(pkin(16));
	t111 = t78 * t105;
	t113 = t63 * t100;
	t114 = t66 * t64;
	t116 = t59 * t67 * t114;
	t118 = t100 * t111 - t113 * t116;
	t119 = sin(pkin(16));
	t121 = t109 * t108 / 0.4e1 + t119 * t118 / 0.4e1;
	t125 = -t57 * t94 + t87 * t91;
	t128 = t109 * t118 / 0.4e1 - t119 * t108 / 0.4e1;
	t130 = t121 * t96 + t125 * t128;
	t131 = cos(qJ(3));
	t133 = sin(qJ(1));
	t134 = sin(qJ(3));
	t136 = t130 * t131 + t133 * t134;
	t137 = sin(qJ(4));
	t141 = -t121 * t125 + t128 * t96;
	t142 = cos(qJ(4));
	t144 = -t136 * t137 - t141 * t142;
	t145 = t48 * t5 / 0.2e1;
	t147 = t145 - t55 * t2 / 0.2e1;
	t150 = t48 * t2 / 0.2e1;
	t151 = -t55 * t5 / 0.2e1 - t150;
	t153 = t147 * t87 + t151 * t94;
	t157 = -t147 * t94 + t151 * t87;
	t159 = t121 * t153 + t128 * t157;
	t160 = t131 * t159;
	t164 = -t121 * t157 + t128 * t153;
	t166 = t137 * t160 + t142 * t164;
	t167 = 0.1e1 / t166;
	t168 = t167 * t144;
	t169 = t2 * t133;
	t170 = t48 * t169 / 0.2e1;
	t171 = t5 * t133;
	t173 = -t170 - t55 * t171 / 0.2e1;
	t176 = t48 * t171 / 0.2e1;
	t177 = t55 * t169 / 0.2e1 - t176;
	t179 = t173 * t87 + t177 * t94;
	t183 = -t173 * t94 + t177 * t87;
	t185 = t121 * t179 + t128 * t183;
	t187 = t134 * t1;
	t188 = t131 * t185 - t187;
	t192 = -t121 * t183 + t128 * t179;
	t194 = -t137 * t188 - t142 * t192;
	t195 = t194 ^ 2;
	t196 = t166 ^ 2;
	t197 = 0.1e1 / t196;
	t200 = 0.1e1 / (t195 * t197 + 0.1e1);
	t202 = -t13 * pkin(2);
	t204 = 0.1e1 / t27;
	t209 = pkin(1) * pkin(2);
	t211 = pkin(1) * pkin(2) * t25 * t9 + t209 * t21 * t9;
	t218 = -t10 * t204 * t211 - 0.2e1 * t209 * t30 * t9 - t202 * t27 - t44;
	t223 = t39 ^ 2;
	t226 = t209 * t9 / t223;
	t232 = t9 ^ 2;
	t236 = -0.2e1 * pkin(1) * t232 * t31 + t204 * t211 * t30 + t202 * t34 - t28;
	t242 = -t41 * t218 * t4 / 0.2e1 - t226 * t38 * t37 + t41 * t236 * t7 / 0.2e1 + t226 * t38 * t46;
	t255 = -t41 * t218 * t7 / 0.2e1 - t226 * t38 * t51 - t41 * t236 * t4 / 0.2e1 - t226 * t38 * t53;
	t257 = -t176 - t242 * t169 + t55 * t169 / 0.2e1 - t255 * t171;
	t262 = t55 * t171 / 0.2e1 + t255 * t169 + t170 - t242 * t171;
	t264 = t257 * t87 + t262 * t94;
	t268 = -t257 * t94 + t262 * t87;
	t277 = -t137 * t131 * (t121 * t264 + t128 * t268) - t142 * (-t121 * t268 + t128 * t264);
	t283 = -t150 + t242 * t5 - t55 * t5 / 0.2e1 - t255 * t2;
	t288 = t55 * t2 / 0.2e1 - t255 * t5 - t145 - t242 * t2;
	t290 = t283 * t87 + t288 * t94;
	t294 = -t283 * t94 + t288 * t87;
	t303 = t137 * t131 * (t121 * t290 + t128 * t294) + t142 * (-t121 * t294 + t128 * t290);
	t305 = t200 * t197;
	t307 = t167 * t200 * t277 - t194 * t303 * t305;
	t309 = t131 * t1;
	t310 = -t134 * t185 - t309;
	t319 = t134 * t137 * t159 * t194 * t197 * t200 - t137 * t167 * t200 * t310;
	t322 = t137 * t192 - t142 * t188;
	t327 = -t137 * t164 + t142 * t160;
	t330 = t167 * t200 * t322 - t194 * t305 * t327;
	t331 = t99 * t59;
	t333 = 0.1e1 / t70;
	t339 = -t63 * t66 * t67 + t114 * t67 - t65 * t66 + t65 * t67;
	t347 = t78 ^ 2;
	t348 = 1 / t347;
	t352 = 0.1e1 / (-t348 * t65 * t68 + 0.1e1);
	t363 = t352 / t78 * t60 * pkin(7) * (-t70 * t331 + t339 * t333 * t62 / 0.2e1) - t352 * t348 * t70 * t60 * pkin(7) * (-t331 * t78 - 0.2e1 * t60 * t62);
	t364 = t363 * t58;
	t366 = t363 * t84;
	t368 = t364 * t85 - t366 * t82;
	t372 = t364 * t82 + t366 * t85;
	t374 = t173 * t368 + t177 * t372;
	t378 = 0.1e1 / t98 / t60 * t97;
	t391 = 0.2e1 * t59 * t60;
	t400 = t102 * t378 / 0.2e1 - t339 * t78 * t59 * t333 * t100 / 0.8e1 + t60 * t101 * t100 / 0.2e1 - t106 * t378 / 0.2e1 + t70 * t391 * t100 / 0.4e1 + t339 * t333 * t59 * t104 * t100 / 0.8e1;
	t414 = t59 * t68;
	t426 = -t111 * t378 / 0.2e1 + t78 * t391 * t100 / 0.4e1 - t60 * t105 * t100 / 0.2e1 + t116 * t63 * t378 / 0.2e1 + t414 * t64 * t100 / 0.4e1 - t414 * t113 / 0.4e1 + t59 * t67 * t64 * t113 / 0.4e1 - t59 * t114 * t113 / 0.4e1;
	t428 = t109 * t400 + t119 * t426;
	t432 = -t173 * t372 + t177 * t368;
	t436 = t109 * t426 - t119 * t400;
	t447 = -t137 * t131 * (t121 * t374 + t128 * t432 + t179 * t428 + t183 * t436) - t142 * (-t121 * t432 + t128 * t374 + t179 * t436 - t183 * t428);
	t452 = t147 * t368 + t151 * t372;
	t457 = -t147 * t372 + t151 * t368;
	t469 = t137 * t131 * (t121 * t452 + t128 * t457 + t153 * t428 + t157 * t436) + t142 * (-t121 * t457 + t128 * t452 + t153 * t436 - t157 * t428);
	t472 = t167 * t200 * t447 - t194 * t305 * t469;
	t475 = -t173 * t87 - t177 * t94;
	t479 = t173 * t94 - t177 * t87;
	t481 = t121 * t475 + t128 * t479;
	t483 = t131 * t481 + t187;
	t487 = -t121 * t479 + t128 * t475;
	t490 = atan2(t194, t166);
	t491 = cos(t490);
	t493 = sin(t490);
	t495 = t166 * t491 + t194 * t493;
	t496 = 0.1e1 / t495;
	t498 = t144 ^ 2;
	t499 = t495 ^ 2;
	t500 = 0.1e1 / t499;
	t503 = 0.1e1 / (t498 * t500 + 0.1e1);
	t513 = t503 * t500;
	t519 = -t90 - t242 * t3 + t55 * t3 / 0.2e1 - t255 * t50;
	t524 = t55 * t50 / 0.2e1 + t255 * t3 + t49 - t242 * t50;
	t526 = t519 * t87 + t524 * t94;
	t530 = -t519 * t94 + t524 * t87;
	t532 = t121 * t526 + t128 * t530;
	t533 = t131 * t532;
	t537 = -t121 * t530 + t128 * t526;
	t554 = -t130 * t134 + t131 * t133;
	t573 = t136 * t142 - t137 * t141;
	t588 = t368 * t57 + t372 * t91;
	t593 = t368 * t91 - t372 * t57;
	t596 = t121 * t588 + t125 * t436 + t128 * t593 + t428 * t96;
	t597 = t131 * t596;
	t603 = -t121 * t593 - t125 * t428 + t128 * t588 + t436 * t96;
	t620 = -t137 * t487 + t142 * t483;
	t621 = sin(qJ(5));
	t624 = -t134 * t481 + t309;
	t625 = cos(qJ(5));
	t630 = t554 * t621 + t573 * t625;
	t631 = 0.1e1 / t630;
	t635 = -t554 * t625 + t573 * t621;
	t636 = t635 ^ 2;
	t637 = t630 ^ 2;
	t638 = 0.1e1 / t637;
	t641 = 0.1e1 / (t636 * t638 + 0.1e1);
	t647 = t641 * t638;
	t652 = -t137 * t537 + t142 * t533;
	t654 = t134 * t532;
	t665 = t142 * t554;
	t690 = -t137 * t603 + t142 * t597;
	t692 = t134 * t596;
	unknown(1,1) = t200 * t168;
	unknown(1,2) = t307;
	unknown(1,3) = t319;
	unknown(1,4) = t330;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = t472;
	unknown(2,1) = t503 * t496 * (t137 * t483 + t142 * t487) + t513 * t144 * (t168 * t194 * t200 * t491 - t144 * t200 * t493 + t144 * t493);
	unknown(2,2) = t503 * t496 * (t137 * t533 + t142 * t537) + t513 * t144 * (-t166 * t307 * t493 + t194 * t307 * t491 + t277 * t493 + t303 * t491);
	unknown(2,3) = t503 * t496 * t137 * t554 + t513 * t144 * (-t134 * t137 * t159 * t491 - t137 * t310 * t493 - t166 * t319 * t493 + t194 * t319 * t491);
	unknown(2,4) = t503 * t496 * t573 + t513 * t144 * (-t166 * t330 * t493 + t194 * t330 * t491 + t322 * t493 + t327 * t491);
	unknown(2,5) = 0.0e0;
	unknown(2,6) = t503 * t496 * (t137 * t597 + t142 * t603) + t513 * t144 * (-t166 * t472 * t493 + t194 * t472 * t491 + t447 * t493 + t469 * t491);
	unknown(3,1) = t641 * t631 * (t620 * t621 - t624 * t625) - t647 * t635 * (t620 * t625 + t621 * t624);
	unknown(3,2) = t641 * t631 * (t621 * t652 + t625 * t654) - t647 * t635 * (-t621 * t654 + t625 * t652);
	unknown(3,3) = t641 * t631 * (t136 * t625 + t621 * t665) - t647 * t635 * (-t136 * t621 + t625 * t665);
	unknown(3,4) = -t144 * t625 * t635 * t638 * t641 + t144 * t621 * t631 * t641;
	unknown(3,5) = t635 ^ 2 * t647 + t641;
	unknown(3,6) = t641 * t631 * (t621 * t690 + t625 * t692) - t647 * t635 * (-t621 * t692 + t625 * t690);
	Ja_rot = unknown;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobia_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:25:24
	% EndTime: 2020-04-11 12:25:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobia_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:25:21
	% EndTime: 2020-04-11 12:25:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
elseif link_index == 11
	%% Symbolic Calculation
	% From jacobia_rot_11_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:25:27
	% EndTime: 2020-04-11 12:25:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
elseif link_index == 12
	%% Symbolic Calculation
	% From jacobia_rot_12_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:25:22
	% EndTime: 2020-04-11 12:25:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
else
	Ja_rot=NaN(3,6);
end