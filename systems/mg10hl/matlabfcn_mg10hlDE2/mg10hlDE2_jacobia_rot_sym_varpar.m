% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% mg10hlDE2
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
%   Wie in mg10hlDE2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [17x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,AE,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 13:01
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = mg10hlDE2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'mg10hlDE2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'mg10hlDE2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [17 1]), ...
  'mg10hlDE2_jacobia_rot_sym_varpar: pkin has to be [17x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:00:12
	% EndTime: 2020-04-11 13:00:12
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
	% StartTime: 2020-04-11 13:00:12
	% EndTime: 2020-04-11 13:00:12
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
	% StartTime: 2020-04-11 13:00:12
	% EndTime: 2020-04-11 13:00:12
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
	% StartTime: 2020-04-11 13:00:13
	% EndTime: 2020-04-11 13:00:13
	% DurationCPUTime: 0.07s
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
	% StartTime: 2020-04-11 13:00:14
	% EndTime: 2020-04-11 13:00:14
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
	% StartTime: 2020-04-11 13:00:15
	% EndTime: 2020-04-11 13:00:15
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
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:00:15
	% EndTime: 2020-04-11 13:00:17
	% DurationCPUTime: 0.81s
	% Computational Cost: add. (55640->120), mult. (50637->281), div. (8573->29), fcn. (18655->22), ass. (0->158)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = sin(qJ(2));
	t3 = cos(pkin(15));
	t5 = cos(qJ(2));
	t6 = sin(pkin(15));
	t8 = t3 * t2 + t6 * t5;
	t10 = -t8 * pkin(2) + pkin(1);
	t13 = 0.2e1 * pkin(2) * pkin(1) * t8;
	t14 = pkin(1) ^ 2;
	t18 = -t13 + t14 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t22 = -t13 + t14 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t24 = sqrt(-t22 * t18);
	t28 = -t6 * t2 + t3 * t5;
	t29 = t28 * pkin(2);
	t30 = pkin(2) ^ 2;
	t31 = pkin(5) ^ 2;
	t32 = pkin(4) ^ 2;
	t33 = -t13 + t14 + t30 - t31 + t32;
	t34 = t33 * t29;
	t35 = t24 * t10 + t34;
	t36 = 0.1e1 / pkin(4);
	t37 = t36 * t35;
	t38 = -t13 + t14 + t30;
	t39 = 0.1e1 / t38;
	t41 = t24 * t29;
	t43 = t33 * t10 - t41;
	t44 = t36 * t43;
	t46 = atan2(t39 * t37, t39 * t44);
	t47 = 0.1e1 / pkin(7);
	t48 = qJ(6) + pkin(8);
	t50 = 0.1e1 / t48 * t47;
	t51 = pkin(6) - pkin(7) - pkin(8) - qJ(6);
	t52 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t53 = t52 * t51;
	t54 = pkin(6) + pkin(7) - pkin(8) - qJ(6);
	t55 = pkin(6) + pkin(7) + pkin(8) + qJ(6);
	t56 = t55 * t54;
	t58 = sqrt(-t56 * t53);
	t60 = (pkin(6) ^ 2);
	t61 = (pkin(7) ^ 2);
	t62 = (pkin(8) ^ 2);
	t64 = 2 * pkin(8) * qJ(6);
	t65 = qJ(6) ^ 2;
	t66 = t60 - t61 - t62 - t64 - t65;
	t68 = atan2(t58 * t50, t66 * t50);
	t69 = 0.1e1 / pkin(6);
	t70 = t48 ^ 2;
	t71 = 0.1e1 / t70;
	t72 = t71 * t69;
	t73 = t47 * t58;
	t74 = t66 * t73;
	t76 = t60 - t61 + t62 + t64 + t65;
	t77 = t47 * t76;
	t78 = t58 * t77;
	t80 = -t74 * t72 + t78 * t72;
	t81 = cos(pkin(16));
	t83 = t66 * t77;
	t85 = t51 * t72;
	t86 = t54 * t52;
	t88 = t47 * t55 * t86;
	t90 = t83 * t72 - t88 * t85;
	t91 = sin(pkin(16));
	t93 = t81 * t80 / 0.4e1 + t91 * t90 / 0.4e1;
	t96 = -t81 * t90 / 0.4e1 + t91 * t80 / 0.4e1;
	t97 = atan2(t93, t96);
	t98 = qJ(2) + pkin(15) + t46 + pkin(17) - t68 - pkin(16) + t97;
	t99 = sin(t98);
	t100 = t99 * t1;
	t101 = cos(t98);
	t102 = 0.1e1 / t101;
	t103 = sin(qJ(1));
	t104 = t103 ^ 2;
	t105 = t99 ^ 2;
	t107 = t101 ^ 2;
	t108 = 0.1e1 / t107;
	t111 = 0.1e1 / (t104 * t105 * t108 + 0.1e1);
	t114 = 0.1e1 / t24;
	t116 = pkin(1) * t28;
	t120 = pkin(1) * pkin(2);
	t122 = pkin(2) * t116 * t22 + t120 * t18 * t28;
	t125 = -t8 * pkin(2);
	t127 = t28 ^ 2;
	t134 = t38 ^ 2;
	t135 = 0.1e1 / t134;
	t137 = pkin(2) * t116;
	t144 = t35 ^ 2;
	t145 = t43 ^ 2;
	t146 = 0.1e1 / t145;
	t149 = 0.1e1 / (t146 * t144 + 0.1e1);
	t171 = 0.1e1 + t149 * t38 * pkin(4) / t43 * (t39 * t36 * (-0.2e1 * pkin(1) * t127 * t30 + t10 * t114 * t122 + t33 * t125 - t41) + 0.2e1 * t137 * t135 * t37) - t149 * t146 * t38 * pkin(4) * t35 * (t39 * t36 * (-0.2e1 * t10 * t120 * t28 - t114 * t122 * t29 - t24 * t125 - t34) + 0.2e1 * t137 * t135 * t44);
	t176 = t111 * t108 * t103;
	t178 = t103 * t111 * t171 + t105 * t171 * t176;
	t179 = t71 * t47;
	t181 = 0.1e1 / t58;
	t187 = -t51 * t54 * t55 - t54 * t53 + t55 * t53 + t55 * t86;
	t195 = t66 ^ 2;
	t196 = 1 / t195;
	t200 = 0.1e1 / (-t196 * t53 * t56 + 0.1e1);
	t213 = 0.1e1 / t70 / t48 * t69;
	t226 = 0.2e1 * t47 * t48;
	t235 = t74 * t213 / 0.2e1 - t187 * t66 * t47 * t181 * t72 / 0.8e1 + t48 * t73 * t72 / 0.2e1 - t78 * t213 / 0.2e1 + t58 * t226 * t72 / 0.4e1 + t187 * t181 * t47 * t76 * t72 / 0.8e1;
	t249 = t47 * t56;
	t261 = -t83 * t213 / 0.2e1 + t66 * t226 * t72 / 0.4e1 - t48 * t77 * t72 / 0.2e1 + t88 * t51 * t213 / 0.2e1 + t249 * t52 * t72 / 0.4e1 - t249 * t85 / 0.4e1 + t47 * t55 * t52 * t85 / 0.4e1 - t47 * t86 * t85 / 0.4e1;
	t266 = t93 ^ 2;
	t267 = t96 ^ 2;
	t268 = 0.1e1 / t267;
	t271 = 0.1e1 / (t268 * t266 + 0.1e1);
	t279 = -t200 / t66 * t48 * pkin(7) * (-t58 * t179 + t187 * t181 * t50 / 0.2e1) + t200 * t196 * t58 * t48 * pkin(7) * (-t66 * t179 - 0.2e1 * t48 * t50) + t271 / t96 * (t81 * t235 + t91 * t261) - t271 * t268 * t93 * (t91 * t235 - t81 * t261);
	t284 = t103 * t111 * t279 + t105 * t176 * t279;
	t285 = t99 * t103;
	t286 = atan2(-t285, -t101);
	t287 = cos(t286);
	t289 = sin(t286);
	t290 = t103 * t289;
	t292 = -t101 * t287 - t99 * t290;
	t293 = 0.1e1 / t292;
	t294 = t1 ^ 2;
	t296 = t292 ^ 2;
	t297 = 0.1e1 / t296;
	t300 = 0.1e1 / (t105 * t294 * t297 + 0.1e1);
	t315 = t300 * t297 * t99;
	t318 = t171 * t1;
	t320 = t300 * t293 * t101;
	t334 = t279 * t1;
	t348 = t101 * t103;
	t349 = sin(qJ(3));
	t351 = cos(qJ(3));
	t354 = t101 * t1;
	t357 = t349 * t103 + t351 * t354;
	t358 = 0.1e1 / t357;
	t362 = -t351 * t103 + t349 * t354;
	t363 = t362 ^ 2;
	t364 = t357 ^ 2;
	t365 = 0.1e1 / t364;
	t368 = 0.1e1 / (t365 * t363 + 0.1e1);
	t374 = t368 * t365;
	t377 = t99 * t318;
	t379 = t368 * t358 * t349;
	t382 = t374 * t362 * t351;
	t388 = t99 * t334;
	unknown(1,1) = t111 * t102 * t100;
	unknown(1,2) = t178;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = t284;
	unknown(2,1) = -t300 * t293 * t285 - t315 * t1 * (-t1 * t102 * t103 * t105 * t111 * t287 - t1 * t289 * t99 + t100 * t111 * t289);
	unknown(2,2) = t320 * t318 - t315 * t1 * (-t101 * t171 * t290 + t101 * t178 * t289 + t171 * t287 * t99 - t178 * t285 * t287);
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = t320 * t334 - t315 * t1 * (-t101 * t279 * t290 + t101 * t284 * t289 + t279 * t287 * t99 - t284 * t285 * t287);
	unknown(3,1) = t368 * t358 * (-t351 * t1 - t349 * t348) - t374 * t362 * (t349 * t1 - t351 * t348);
	unknown(3,2) = -t379 * t377 + t382 * t377;
	unknown(3,3) = t362 ^ 2 * t374 + t368;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = -t379 * t388 + t382 * t388;
	Ja_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:00:16
	% EndTime: 2020-04-11 13:00:19
	% DurationCPUTime: 1.26s
	% Computational Cost: add. (83054->134), mult. (76202->329), div. (12914->31), fcn. (28656->24), ass. (0->173)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = sin(qJ(2));
	t3 = cos(pkin(15));
	t5 = cos(qJ(2));
	t6 = sin(pkin(15));
	t8 = t3 * t2 + t6 * t5;
	t10 = -t8 * pkin(2) + pkin(1);
	t13 = 0.2e1 * pkin(2) * pkin(1) * t8;
	t14 = pkin(1) ^ 2;
	t18 = -t13 + t14 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t22 = -t13 + t14 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t24 = sqrt(-t22 * t18);
	t28 = -t6 * t2 + t3 * t5;
	t29 = t28 * pkin(2);
	t30 = pkin(2) ^ 2;
	t31 = pkin(5) ^ 2;
	t32 = pkin(4) ^ 2;
	t33 = -t13 + t14 + t30 - t31 + t32;
	t34 = t33 * t29;
	t35 = t24 * t10 + t34;
	t36 = 0.1e1 / pkin(4);
	t37 = t36 * t35;
	t38 = -t13 + t14 + t30;
	t39 = 0.1e1 / t38;
	t41 = t24 * t29;
	t43 = t33 * t10 - t41;
	t44 = t36 * t43;
	t46 = atan2(t39 * t37, t39 * t44);
	t47 = 0.1e1 / pkin(7);
	t48 = qJ(6) + pkin(8);
	t50 = 0.1e1 / t48 * t47;
	t51 = pkin(6) - pkin(7) - pkin(8) - qJ(6);
	t52 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t53 = t52 * t51;
	t54 = pkin(6) + pkin(7) - pkin(8) - qJ(6);
	t55 = pkin(6) + pkin(7) + pkin(8) + qJ(6);
	t56 = t55 * t54;
	t58 = sqrt(-t56 * t53);
	t60 = (pkin(6) ^ 2);
	t61 = (pkin(7) ^ 2);
	t62 = (pkin(8) ^ 2);
	t64 = 2 * pkin(8) * qJ(6);
	t65 = qJ(6) ^ 2;
	t66 = t60 - t61 - t62 - t64 - t65;
	t68 = atan2(t58 * t50, t66 * t50);
	t69 = 0.1e1 / pkin(6);
	t70 = t48 ^ 2;
	t71 = 0.1e1 / t70;
	t72 = t71 * t69;
	t73 = t47 * t58;
	t74 = t66 * t73;
	t76 = t60 - t61 + t62 + t64 + t65;
	t77 = t47 * t76;
	t78 = t58 * t77;
	t80 = -t74 * t72 + t78 * t72;
	t81 = cos(pkin(16));
	t83 = t66 * t77;
	t85 = t51 * t72;
	t86 = t54 * t52;
	t88 = t47 * t55 * t86;
	t90 = t83 * t72 - t88 * t85;
	t91 = sin(pkin(16));
	t93 = t81 * t80 / 0.4e1 + t91 * t90 / 0.4e1;
	t96 = -t81 * t90 / 0.4e1 + t91 * t80 / 0.4e1;
	t97 = atan2(t93, t96);
	t98 = qJ(2) + pkin(15) + t46 + pkin(17) - t68 - pkin(16) + t97;
	t99 = cos(t98);
	t100 = t99 * t1;
	t101 = sin(qJ(3));
	t103 = sin(qJ(1));
	t104 = cos(qJ(3));
	t106 = t101 * t100 - t104 * t103;
	t107 = sin(t98);
	t108 = 0.1e1 / t107;
	t109 = t108 * t106;
	t110 = 0.1e1 / t101;
	t111 = t99 * t103;
	t114 = t104 * t1 + t101 * t111;
	t115 = t114 ^ 2;
	t116 = t107 ^ 2;
	t117 = 0.1e1 / t116;
	t119 = t101 ^ 2;
	t120 = 0.1e1 / t119;
	t123 = 0.1e1 / (t120 * t117 * t115 + 0.1e1);
	t124 = t123 * t110;
	t126 = 0.1e1 / t24;
	t128 = pkin(1) * t28;
	t132 = pkin(1) * pkin(2);
	t134 = t22 * pkin(2) * t128 + t132 * t28 * t18;
	t137 = -t8 * pkin(2);
	t139 = t28 ^ 2;
	t146 = t38 ^ 2;
	t147 = 0.1e1 / t146;
	t149 = pkin(2) * t128;
	t156 = t35 ^ 2;
	t157 = t43 ^ 2;
	t158 = 0.1e1 / t157;
	t161 = 0.1e1 / (t158 * t156 + 0.1e1);
	t183 = 0.1e1 + t161 * t38 * pkin(4) / t43 * (t39 * t36 * (-0.2e1 * pkin(1) * t139 * t30 + t134 * t126 * t10 + t33 * t137 - t41) + 0.2e1 * t149 * t147 * t37) - t161 * t158 * t38 * pkin(4) * t35 * (t39 * t36 * (-0.2e1 * t132 * t28 * t10 - t134 * t126 * t29 - t24 * t137 - t34) + 0.2e1 * t149 * t147 * t44);
	t189 = t123 * t117 * t114;
	t191 = t189 * t110 * t99 * t183 + t123 * t183 * t103;
	t194 = -t101 * t1 + t104 * t111;
	t201 = t123 * t120 * t114 * t104 * t108 - t124 * t108 * t194;
	t202 = t71 * t47;
	t204 = 0.1e1 / t58;
	t210 = -t55 * t54 * t51 - t54 * t53 + t55 * t53 + t55 * t86;
	t218 = t66 ^ 2;
	t219 = 1 / t218;
	t223 = 0.1e1 / (-t219 * t56 * t53 + 0.1e1);
	t236 = 0.1e1 / t70 / t48 * t69;
	t249 = 0.2e1 * t47 * t48;
	t258 = t74 * t236 / 0.2e1 - t210 * t66 * t47 * t204 * t72 / 0.8e1 + t48 * t73 * t72 / 0.2e1 - t78 * t236 / 0.2e1 + t58 * t249 * t72 / 0.4e1 + t210 * t204 * t47 * t76 * t72 / 0.8e1;
	t272 = t47 * t56;
	t284 = -t83 * t236 / 0.2e1 + t66 * t249 * t72 / 0.4e1 - t48 * t77 * t72 / 0.2e1 + t88 * t51 * t236 / 0.2e1 + t272 * t52 * t72 / 0.4e1 - t272 * t85 / 0.4e1 + t47 * t55 * t52 * t85 / 0.4e1 - t47 * t86 * t85 / 0.4e1;
	t289 = t93 ^ 2;
	t290 = t96 ^ 2;
	t291 = 0.1e1 / t290;
	t294 = 0.1e1 / (t291 * t289 + 0.1e1);
	t302 = -t223 / t66 * t48 * pkin(7) * (-t58 * t202 + t210 * t204 * t50 / 0.2e1) + t223 * t219 * t58 * t48 * pkin(7) * (-t66 * t202 - 0.2e1 * t48 * t50) + t294 / t96 * (t81 * t258 + t91 * t284) - t294 * t291 * t93 * (t91 * t258 - t81 * t284);
	t308 = t189 * t110 * t99 * t302 + t123 * t302 * t103;
	t309 = t101 * t107;
	t310 = atan2(t114, -t309);
	t311 = cos(t310);
	t312 = t107 * t311;
	t314 = sin(t310);
	t316 = -t101 * t312 + t114 * t314;
	t317 = 0.1e1 / t316;
	t319 = t106 ^ 2;
	t320 = t316 ^ 2;
	t321 = 0.1e1 / t320;
	t324 = 0.1e1 / (t321 * t319 + 0.1e1);
	t335 = t324 * t321;
	t338 = t183 * t1;
	t341 = t324 * t317 * t101;
	t346 = t101 * t99;
	t350 = t103 * t314;
	t360 = -t104 * t100 - t101 * t103;
	t373 = t302 * t1;
	t389 = cos(qJ(4));
	t391 = t107 * t103;
	t392 = sin(qJ(4));
	t396 = t107 * t1;
	t398 = -t392 * t360 + t389 * t396;
	t399 = 0.1e1 / t398;
	t403 = t389 * t360 + t392 * t396;
	t404 = t403 ^ 2;
	t405 = t398 ^ 2;
	t406 = 0.1e1 / t405;
	t409 = 0.1e1 / (t406 * t404 + 0.1e1);
	t415 = t409 * t406;
	t418 = t104 * t107;
	t419 = t389 * t418;
	t421 = t392 * t99;
	t426 = t392 * t418;
	t428 = t389 * t99;
	unknown(1,1) = -t124 * t109;
	unknown(1,2) = t191;
	unknown(1,3) = t201;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = t308;
	unknown(2,1) = t324 * t317 * t114 + t335 * t106 * (-t114 * t311 * t123 * t110 * t109 - t314 * t123 * t106 + t106 * t314);
	unknown(2,2) = t341 * t107 * t338 + t335 * t106 * (-t101 * t107 * t183 * t350 + t114 * t311 * t191 - t346 * t183 * t311 + t309 * t314 * t191);
	unknown(2,3) = t324 * t317 * t360 + t335 * t106 * (t114 * t311 * t201 + t309 * t314 * t201 - t104 * t312 + t194 * t314);
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = t341 * t107 * t373 + t335 * t106 * (-t101 * t107 * t302 * t350 + t114 * t311 * t308 - t346 * t302 * t311 + t309 * t314 * t308);
	unknown(3,1) = t409 * t399 * (t389 * t194 - t392 * t391) - t415 * t403 * (-t392 * t194 - t389 * t391);
	unknown(3,2) = t409 * t399 * (t419 * t338 + t421 * t338) - t415 * t403 * (-t426 * t338 + t428 * t338);
	unknown(3,3) = t409 * t406 * t403 * t392 * t106 + t409 * t399 * t389 * t106;
	unknown(3,4) = t415 * t403 ^ 2 + t409;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = t409 * t399 * (t419 * t373 + t421 * t373) - t415 * t403 * (-t426 * t373 + t428 * t373);
	Ja_rot = unknown;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:00:17
	% EndTime: 2020-04-11 13:00:22
	% DurationCPUTime: 2.71s
	% Computational Cost: add. (191220->156), mult. (176168->374), div. (29729->29), fcn. (66565->26), ass. (0->189)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = sin(qJ(2));
	t3 = cos(pkin(15));
	t5 = cos(qJ(2));
	t6 = sin(pkin(15));
	t8 = t3 * t2 + t6 * t5;
	t10 = -t8 * pkin(2) + pkin(1);
	t13 = 0.2e1 * pkin(2) * pkin(1) * t8;
	t14 = pkin(1) ^ 2;
	t18 = -t13 + t14 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t22 = -t13 + t14 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t24 = sqrt(-t22 * t18);
	t28 = -t6 * t2 + t3 * t5;
	t29 = t28 * pkin(2);
	t30 = pkin(2) ^ 2;
	t31 = pkin(5) ^ 2;
	t32 = pkin(4) ^ 2;
	t33 = -t13 + t14 + t30 - t31 + t32;
	t34 = t33 * t29;
	t35 = t24 * t10 + t34;
	t36 = 0.1e1 / pkin(4);
	t37 = t36 * t35;
	t38 = -t13 + t14 + t30;
	t39 = 0.1e1 / t38;
	t41 = t24 * t29;
	t43 = t33 * t10 - t41;
	t44 = t36 * t43;
	t46 = atan2(t39 * t37, t39 * t44);
	t47 = 0.1e1 / pkin(7);
	t48 = qJ(6) + pkin(8);
	t50 = 0.1e1 / t48 * t47;
	t51 = pkin(6) - pkin(7) - pkin(8) - qJ(6);
	t52 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t53 = t52 * t51;
	t54 = pkin(6) + pkin(7) - pkin(8) - qJ(6);
	t55 = pkin(6) + pkin(7) + pkin(8) + qJ(6);
	t56 = t55 * t54;
	t58 = sqrt(-t56 * t53);
	t60 = (pkin(6) ^ 2);
	t61 = (pkin(7) ^ 2);
	t62 = (pkin(8) ^ 2);
	t64 = 2 * pkin(8) * qJ(6);
	t65 = qJ(6) ^ 2;
	t66 = t60 - t61 - t62 - t64 - t65;
	t68 = atan2(t58 * t50, t66 * t50);
	t69 = 0.1e1 / pkin(6);
	t70 = t48 ^ 2;
	t71 = 0.1e1 / t70;
	t72 = t71 * t69;
	t73 = t47 * t58;
	t74 = t66 * t73;
	t76 = t60 - t61 + t62 + t64 + t65;
	t77 = t47 * t76;
	t78 = t58 * t77;
	t80 = -t74 * t72 + t78 * t72;
	t81 = cos(pkin(16));
	t83 = t66 * t77;
	t85 = t51 * t72;
	t86 = t54 * t52;
	t88 = t47 * t55 * t86;
	t90 = t83 * t72 - t88 * t85;
	t91 = sin(pkin(16));
	t93 = t81 * t80 / 0.4e1 + t91 * t90 / 0.4e1;
	t96 = -t81 * t90 / 0.4e1 + t91 * t80 / 0.4e1;
	t97 = atan2(t93, t96);
	t98 = qJ(2) + pkin(15) + t46 + pkin(17) - t68 - pkin(16) + t97;
	t99 = cos(t98);
	t100 = t99 * t1;
	t101 = cos(qJ(3));
	t103 = sin(qJ(1));
	t104 = sin(qJ(3));
	t106 = t101 * t100 + t104 * t103;
	t107 = sin(qJ(4));
	t109 = sin(t98);
	t110 = t109 * t1;
	t111 = cos(qJ(4));
	t113 = -t107 * t106 - t111 * t110;
	t114 = t101 * t109;
	t115 = t107 * t114;
	t116 = t111 * t99;
	t117 = t115 - t116;
	t118 = 0.1e1 / t117;
	t119 = t118 * t113;
	t120 = t99 * t103;
	t123 = -t104 * t1 + t101 * t120;
	t125 = t109 * t103;
	t126 = t111 * t125;
	t127 = -t107 * t123 - t126;
	t128 = t127 ^ 2;
	t129 = t117 ^ 2;
	t130 = 0.1e1 / t129;
	t133 = 0.1e1 / (t130 * t128 + 0.1e1);
	t135 = 0.1e1 / t24;
	t137 = pkin(1) * t28;
	t141 = pkin(1) * pkin(2);
	t143 = t22 * pkin(2) * t137 + t141 * t28 * t18;
	t146 = -t8 * pkin(2);
	t148 = t28 ^ 2;
	t155 = t38 ^ 2;
	t156 = 0.1e1 / t155;
	t158 = pkin(2) * t137;
	t165 = t35 ^ 2;
	t166 = t43 ^ 2;
	t167 = 0.1e1 / t166;
	t170 = 0.1e1 / (t167 * t165 + 0.1e1);
	t192 = 0.1e1 + t170 * t38 * pkin(4) / t43 * (t39 * t36 * (-0.2e1 * pkin(1) * t148 * t30 + t143 * t135 * t10 + t33 * t146 - t41) + 0.2e1 * t158 * t156 * t37) - t170 * t167 * t38 * pkin(4) * t35 * (t39 * t36 * (-0.2e1 * t141 * t28 * t10 - t143 * t135 * t29 - t24 * t146 - t34) + 0.2e1 * t158 * t156 * t44);
	t193 = t192 * t103;
	t196 = t115 * t193 - t116 * t193;
	t200 = t107 * t101;
	t204 = t111 * t109 * t192 + t200 * t99 * t192;
	t206 = t133 * t130;
	t208 = t133 * t118 * t196 - t206 * t127 * t204;
	t211 = -t101 * t1 - t104 * t120;
	t215 = t104 * t109;
	t220 = t133 * t130 * t127 * t107 * t215 - t133 * t118 * t107 * t211;
	t222 = t107 * t125;
	t223 = -t111 * t123 + t222;
	t226 = t111 * t114;
	t227 = t107 * t99;
	t228 = t226 + t227;
	t231 = t133 * t118 * t223 - t206 * t127 * t228;
	t232 = t71 * t47;
	t234 = 0.1e1 / t58;
	t240 = -t55 * t54 * t51 - t54 * t53 + t55 * t53 + t55 * t86;
	t248 = t66 ^ 2;
	t249 = 1 / t248;
	t253 = 0.1e1 / (-t249 * t56 * t53 + 0.1e1);
	t266 = 0.1e1 / t70 / t48 * t69;
	t279 = 0.2e1 * t47 * t48;
	t288 = t74 * t266 / 0.2e1 - t240 * t66 * t47 * t234 * t72 / 0.8e1 + t48 * t73 * t72 / 0.2e1 - t78 * t266 / 0.2e1 + t58 * t279 * t72 / 0.4e1 + t240 * t234 * t47 * t76 * t72 / 0.8e1;
	t302 = t47 * t56;
	t314 = -t83 * t266 / 0.2e1 + t66 * t279 * t72 / 0.4e1 - t48 * t77 * t72 / 0.2e1 + t88 * t51 * t266 / 0.2e1 + t302 * t52 * t72 / 0.4e1 - t302 * t85 / 0.4e1 + t47 * t55 * t52 * t85 / 0.4e1 - t47 * t86 * t85 / 0.4e1;
	t319 = t93 ^ 2;
	t320 = t96 ^ 2;
	t321 = 0.1e1 / t320;
	t324 = 0.1e1 / (t321 * t319 + 0.1e1);
	t332 = -t253 / t66 * t48 * pkin(7) * (-t58 * t232 + t240 * t234 * t50 / 0.2e1) + t253 * t249 * t58 * t48 * pkin(7) * (-t66 * t232 - 0.2e1 * t48 * t50) + t324 / t96 * (t81 * t288 + t91 * t314) - t324 * t321 * t93 * (t91 * t288 - t81 * t314);
	t333 = t332 * t103;
	t336 = t115 * t333 - t116 * t333;
	t343 = t111 * t109 * t332 + t200 * t99 * t332;
	t346 = t133 * t118 * t336 - t206 * t127 * t343;
	t349 = atan2(t127, t117);
	t350 = cos(t349);
	t352 = sin(t349);
	t354 = t117 * t350 + t127 * t352;
	t355 = 0.1e1 / t354;
	t357 = t113 ^ 2;
	t358 = t354 ^ 2;
	t359 = 0.1e1 / t358;
	t362 = 0.1e1 / (t359 * t357 + 0.1e1);
	t372 = t362 * t359;
	t375 = t192 * t1;
	t393 = -t104 * t100 + t101 * t103;
	t412 = t111 * t106 - t107 * t110;
	t425 = t332 * t1;
	t442 = -t111 * t123 + t222;
	t443 = sin(qJ(5));
	t445 = cos(qJ(5));
	t450 = t443 * t393 + t445 * t412;
	t451 = 0.1e1 / t450;
	t455 = -t445 * t393 + t443 * t412;
	t456 = t455 ^ 2;
	t457 = t450 ^ 2;
	t458 = 0.1e1 / t457;
	t461 = 0.1e1 / (t458 * t456 + 0.1e1);
	t467 = t461 * t458;
	t472 = -t226 * t375 - t227 * t375;
	t474 = t445 * t215;
	t480 = t443 * t215;
	t486 = t111 * t393;
	t511 = -t226 * t425 - t227 * t425;
	unknown(1,1) = t133 * t119;
	unknown(1,2) = t208;
	unknown(1,3) = t220;
	unknown(1,4) = t231;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = t346;
	unknown(2,1) = t362 * t355 * (-t107 * t123 - t126) + t372 * t113 * (t127 * t350 * t133 * t119 - t352 * t133 * t113 + t113 * t352);
	unknown(2,2) = t362 * t355 * (-t115 * t375 + t116 * t375) + t372 * t113 * (-t117 * t352 * t208 + t127 * t350 * t208 + t196 * t352 + t204 * t350);
	unknown(2,3) = t362 * t355 * t107 * t393 + t372 * t113 * (-t107 * t104 * t109 * t350 - t107 * t211 * t352 - t117 * t352 * t220 + t127 * t350 * t220);
	unknown(2,4) = t362 * t355 * t412 + t372 * t113 * (-t117 * t352 * t231 + t127 * t350 * t231 + t223 * t352 + t228 * t350);
	unknown(2,5) = 0.0e0;
	unknown(2,6) = t362 * t355 * (-t115 * t425 + t116 * t425) + t372 * t113 * (-t117 * t352 * t346 + t127 * t350 * t346 + t336 * t352 + t343 * t350);
	unknown(3,1) = t461 * t451 * (t445 * t211 + t443 * t442) - t467 * t455 * (-t443 * t211 + t445 * t442);
	unknown(3,2) = t461 * t451 * (-t474 * t375 + t443 * t472) - t467 * t455 * (t480 * t375 + t445 * t472);
	unknown(3,3) = t461 * t451 * (t445 * t106 + t443 * t486) - t467 * t455 * (-t443 * t106 + t445 * t486);
	unknown(3,4) = -t461 * t458 * t455 * t445 * t113 + t461 * t451 * t443 * t113;
	unknown(3,5) = t467 * t455 ^ 2 + t461;
	unknown(3,6) = t461 * t451 * (-t474 * t425 + t443 * t511) - t467 * t455 * (t480 * t425 + t445 * t511);
	Ja_rot = unknown;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobia_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:00:13
	% EndTime: 2020-04-11 13:00:13
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
	% StartTime: 2020-04-11 13:00:14
	% EndTime: 2020-04-11 13:00:14
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
	% StartTime: 2020-04-11 13:00:14
	% EndTime: 2020-04-11 13:00:14
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
	% StartTime: 2020-04-11 13:00:14
	% EndTime: 2020-04-11 13:00:14
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
else
	Ja_rot=NaN(3,6);
end