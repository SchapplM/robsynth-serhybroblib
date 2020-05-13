% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% mg10hlOL
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in mg10hlOL_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
% 
% Output:
% Ja_rot [3x13]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 13:06
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = mg10hlOL_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),uint8(0),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'mg10hlOL_jacobia_rot_sym_varpar: qJ has to be [13x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'mg10hlOL_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'mg10hlOL_jacobia_rot_sym_varpar: pkin has to be [16x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->39)
	unknown=NaN(3,13);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(1,9) = 0;
	unknown(1,10) = 0;
	unknown(1,11) = 0;
	unknown(1,12) = 0;
	unknown(1,13) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(2,7) = 0;
	unknown(2,8) = 0;
	unknown(2,9) = 0;
	unknown(2,10) = 0;
	unknown(2,11) = 0;
	unknown(2,12) = 0;
	unknown(2,13) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	unknown(3,9) = 0;
	unknown(3,10) = 0;
	unknown(3,11) = 0;
	unknown(3,12) = 0;
	unknown(3,13) = 0;
	Ja_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (6->3), div. (5->2), fcn. (6->2), ass. (0->45)
	unknown=NaN(3,13);
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
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(1,9) = 0;
	unknown(1,10) = 0;
	unknown(1,11) = 0;
	unknown(1,12) = 0;
	unknown(1,13) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(2,7) = 0;
	unknown(2,8) = 0;
	unknown(2,9) = 0;
	unknown(2,10) = 0;
	unknown(2,11) = 0;
	unknown(2,12) = 0;
	unknown(2,13) = 0;
	unknown(3,1) = (t6 * t8 + t8);
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	unknown(3,9) = 0;
	unknown(3,10) = 0;
	unknown(3,11) = 0;
	unknown(3,12) = 0;
	unknown(3,13) = 0;
	Ja_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->39)
	unknown=NaN(3,13);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(1,7) = NaN;
	unknown(1,8) = NaN;
	unknown(1,9) = NaN;
	unknown(1,10) = NaN;
	unknown(1,11) = NaN;
	unknown(1,12) = NaN;
	unknown(1,13) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(2,7) = NaN;
	unknown(2,8) = NaN;
	unknown(2,9) = NaN;
	unknown(2,10) = NaN;
	unknown(2,11) = NaN;
	unknown(2,12) = NaN;
	unknown(2,13) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	unknown(3,7) = NaN;
	unknown(3,8) = NaN;
	unknown(3,9) = NaN;
	unknown(3,10) = NaN;
	unknown(3,11) = NaN;
	unknown(3,12) = NaN;
	unknown(3,13) = NaN;
	Ja_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->39)
	unknown=NaN(3,13);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(1,7) = NaN;
	unknown(1,8) = NaN;
	unknown(1,9) = NaN;
	unknown(1,10) = NaN;
	unknown(1,11) = NaN;
	unknown(1,12) = NaN;
	unknown(1,13) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(2,7) = NaN;
	unknown(2,8) = NaN;
	unknown(2,9) = NaN;
	unknown(2,10) = NaN;
	unknown(2,11) = NaN;
	unknown(2,12) = NaN;
	unknown(2,13) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	unknown(3,7) = NaN;
	unknown(3,8) = NaN;
	unknown(3,9) = NaN;
	unknown(3,10) = NaN;
	unknown(3,11) = NaN;
	unknown(3,12) = NaN;
	unknown(3,13) = NaN;
	Ja_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:47
	% EndTime: 2020-04-11 13:05:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->39)
	unknown=NaN(3,13);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(1,7) = NaN;
	unknown(1,8) = NaN;
	unknown(1,9) = NaN;
	unknown(1,10) = NaN;
	unknown(1,11) = NaN;
	unknown(1,12) = NaN;
	unknown(1,13) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(2,7) = NaN;
	unknown(2,8) = NaN;
	unknown(2,9) = NaN;
	unknown(2,10) = NaN;
	unknown(2,11) = NaN;
	unknown(2,12) = NaN;
	unknown(2,13) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	unknown(3,7) = NaN;
	unknown(3,8) = NaN;
	unknown(3,9) = NaN;
	unknown(3,10) = NaN;
	unknown(3,11) = NaN;
	unknown(3,12) = NaN;
	unknown(3,13) = NaN;
	Ja_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->39)
	unknown=NaN(3,13);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(1,7) = NaN;
	unknown(1,8) = NaN;
	unknown(1,9) = NaN;
	unknown(1,10) = NaN;
	unknown(1,11) = NaN;
	unknown(1,12) = NaN;
	unknown(1,13) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(2,7) = NaN;
	unknown(2,8) = NaN;
	unknown(2,9) = NaN;
	unknown(2,10) = NaN;
	unknown(2,11) = NaN;
	unknown(2,12) = NaN;
	unknown(2,13) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	unknown(3,7) = NaN;
	unknown(3,8) = NaN;
	unknown(3,9) = NaN;
	unknown(3,10) = NaN;
	unknown(3,11) = NaN;
	unknown(3,12) = NaN;
	unknown(3,13) = NaN;
	Ja_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:49
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (10059->105), mult. (15974->221), div. (145->9), fcn. (23264->19), ass. (0->159)
	unknown=NaN(3,13);
	t1 = cos(qJ(1));
	t2 = sin(qJ(2));
	t3 = t1 * t2;
	t4 = cos(pkin(14));
	t5 = cos(qJ(3));
	t7 = sin(pkin(14));
	t8 = sin(qJ(3));
	t10 = -t4 * t5 + t7 * t8;
	t11 = t3 * t10;
	t12 = cos(qJ(2));
	t13 = t1 * t12;
	t16 = -t4 * t8 - t7 * t5;
	t18 = -t13 * t16 - t11;
	t19 = cos(pkin(16));
	t20 = -qJ(4) + pkin(15);
	t21 = cos(t20);
	t23 = sin(pkin(16));
	t24 = sin(t20);
	t26 = -t19 * t21 - t23 * t24;
	t29 = t13 * t10;
	t30 = t3 * t16 - t29;
	t33 = t19 * t24 - t23 * t21;
	t35 = t18 * t26 + t30 * t33;
	t36 = cos(qJ(5));
	t39 = t30 * t26;
	t40 = -t18 * t33 + t39;
	t41 = sin(qJ(5));
	t42 = t40 * t41;
	t43 = t35 * t36 + t42;
	t44 = t12 * t10;
	t46 = -t2 * t16 + t44;
	t49 = t2 * t10;
	t50 = -t12 * t16 - t49;
	t52 = t46 * t26 + t50 * t33;
	t55 = t50 * t26;
	t56 = -t46 * t33 + t55;
	t58 = -t52 * t36 - t56 * t41;
	t59 = 0.1e1 / t58;
	t60 = t43 * t59;
	t61 = sin(qJ(1));
	t62 = t61 * t2;
	t63 = t62 * t10;
	t64 = t61 * t12;
	t66 = -t64 * t16 - t63;
	t69 = t64 * t10;
	t70 = t62 * t16 - t69;
	t72 = t66 * t26 + t70 * t33;
	t75 = t70 * t26;
	t76 = -t66 * t33 + t75;
	t78 = t72 * t36 + t76 * t41;
	t79 = t78 ^ 2;
	t80 = t58 ^ 2;
	t81 = 0.1e1 / t80;
	t84 = 0.1e1 / (t79 * t81 + 0.1e1);
	t87 = t62 * t16 - t69;
	t89 = -t64 * t16;
	t90 = -t89 + t63;
	t98 = (t87 * t26 + t90 * t33) * t36 + (t90 * t26 - t87 * t33) * t41;
	t102 = -t12 * t16 - t49;
	t104 = -t2 * t16;
	t105 = -t104 - t44;
	t113 = -(t102 * t26 + t105 * t33) * t36 - (-t102 * t33 + t105 * t26) * t41;
	t115 = t81 * t84;
	t117 = -t113 * t78 * t115 + t98 * t59 * t84;
	t119 = t62 * t10 - t89;
	t123 = -t70 * t33;
	t127 = (t119 * t33 + t75) * t36 + (t119 * t26 + t123) * t41;
	t131 = -t12 * t10 - t104;
	t135 = -t50 * t33;
	t139 = -(t131 * t33 + t55) * t36 - (t131 * t26 + t135) * t41;
	t142 = -t139 * t78 * t115 + t127 * t59 * t84;
	t143 = t76 * t36;
	t145 = -t66 * t26 + t123;
	t146 = t145 * t41;
	t147 = t143 + t146;
	t150 = t56 * t36;
	t154 = -t150 - (-t46 * t26 + t135) * t41;
	t157 = -t154 * t78 * t115 + t147 * t59 * t84;
	t159 = -t72 * t41 + t143;
	t163 = t52 * t41 - t150;
	t166 = -t163 * t78 * t115 + t159 * t59 * t84;
	t170 = -t70 * t26 + t66 * t33;
	t173 = atan2(t78, t58);
	t174 = cos(t173);
	t176 = sin(t173);
	t178 = t174 * t58 + t176 * t78;
	t179 = 0.1e1 / t178;
	t181 = t43 ^ 2;
	t182 = t178 ^ 2;
	t183 = 0.1e1 / t182;
	t186 = 0.1e1 / (t181 * t183 + 0.1e1);
	t196 = t183 * t186;
	t200 = t3 * t16 - t29;
	t202 = -t13 * t16;
	t203 = -t202 + t11;
	t205 = t200 * t26 + t203 * t33;
	t209 = -t200 * t33 + t203 * t26;
	t225 = t3 * t10 - t202;
	t227 = t225 * t33 + t39;
	t229 = -t30 * t33;
	t231 = t225 * t26 + t229;
	t246 = t40 * t36;
	t248 = -t18 * t26 + t229;
	t264 = t35 * t41 - t246;
	t278 = -t170 * t36 + t146;
	t279 = sin(qJ(6));
	t281 = cos(qJ(6));
	t286 = t264 * t281 + t61 * t279;
	t287 = 0.1e1 / t286;
	t291 = t264 * t279 - t61 * t281;
	t292 = t291 ^ 2;
	t293 = t286 ^ 2;
	t294 = 0.1e1 / t293;
	t297 = 0.1e1 / (t292 * t294 + 0.1e1);
	t303 = t294 * t297;
	t308 = t205 * t41 - t209 * t36;
	t310 = t287 * t297;
	t314 = t291 * t294 * t297;
	t319 = t227 * t41 - t231 * t36;
	t326 = -t248 * t36 + t42;
	unknown(1,1) = t60 * t84;
	unknown(1,2) = t117;
	unknown(1,3) = t142;
	unknown(1,4) = t157;
	unknown(1,5) = t166;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(1,9) = 0.0e0;
	unknown(1,10) = 0.0e0;
	unknown(1,11) = 0.0e0;
	unknown(1,12) = 0.0e0;
	unknown(1,13) = 0.0e0;
	unknown(2,1) = (-t145 * t36 - t170 * t41) * t179 * t186 + (t60 * t84 * t174 * t78 - t43 * t84 * t176 + t176 * t43) * t43 * t196;
	unknown(2,2) = (-t205 * t36 - t209 * t41) * t179 * t186 + (t117 * t174 * t78 - t117 * t176 * t58 + t174 * t113 + t176 * t98) * t43 * t196;
	unknown(2,3) = (-t227 * t36 - t231 * t41) * t179 * t186 + (t142 * t174 * t78 - t142 * t176 * t58 + t176 * t127 + t174 * t139) * t43 * t196;
	unknown(2,4) = (-t248 * t41 - t246) * t179 * t186 + (t157 * t174 * t78 - t157 * t176 * t58 + t176 * t147 + t174 * t154) * t43 * t196;
	unknown(2,5) = t264 * t179 * t186 + (t166 * t174 * t78 - t166 * t176 * t58 + t176 * t159 + t174 * t163) * t43 * t196;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(2,9) = 0.0e0;
	unknown(2,10) = 0.0e0;
	unknown(2,11) = 0.0e0;
	unknown(2,12) = 0.0e0;
	unknown(2,13) = 0.0e0;
	unknown(3,1) = (-t1 * t281 + t278 * t279) * t287 * t297 - (t1 * t279 + t278 * t281) * t291 * t303;
	unknown(3,2) = t308 * t279 * t310 - t308 * t281 * t314;
	unknown(3,3) = t319 * t279 * t310 - t319 * t281 * t314;
	unknown(3,4) = t326 * t279 * t310 - t326 * t281 * t314;
	unknown(3,5) = t43 * t279 * t310 - t43 * t281 * t314;
	unknown(3,6) = t291 ^ 2 * t303 + t297;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	unknown(3,9) = 0.0e0;
	unknown(3,10) = 0.0e0;
	unknown(3,11) = 0.0e0;
	unknown(3,12) = 0.0e0;
	unknown(3,13) = 0.0e0;
	Ja_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:49
	% EndTime: 2020-04-11 13:05:50
	% DurationCPUTime: 0.57s
	% Computational Cost: add. (13962->127), mult. (22488->287), div. (228->11), fcn. (32607->21), ass. (0->182)
	unknown=NaN(3,13);
	t1 = cos(qJ(1));
	t2 = sin(qJ(2));
	t3 = t1 * t2;
	t4 = cos(pkin(14));
	t5 = cos(qJ(3));
	t7 = sin(pkin(14));
	t8 = sin(qJ(3));
	t10 = -t4 * t5 + t7 * t8;
	t11 = t3 * t10;
	t12 = cos(qJ(2));
	t13 = t1 * t12;
	t16 = -t4 * t8 - t7 * t5;
	t18 = -t13 * t16 - t11;
	t19 = cos(pkin(16));
	t20 = -qJ(4) + pkin(15);
	t21 = cos(t20);
	t23 = sin(pkin(16));
	t24 = sin(t20);
	t26 = -t19 * t21 - t23 * t24;
	t29 = t13 * t10;
	t30 = t3 * t16 - t29;
	t33 = t19 * t24 - t23 * t21;
	t35 = t18 * t26 + t30 * t33;
	t36 = sin(qJ(5));
	t39 = t30 * t26;
	t40 = -t18 * t33 + t39;
	t41 = cos(qJ(5));
	t42 = t40 * t41;
	t43 = t35 * t36 - t42;
	t44 = sin(qJ(6));
	t46 = sin(qJ(1));
	t47 = cos(qJ(6));
	t49 = t43 * t44 - t46 * t47;
	t50 = t12 * t10;
	t52 = -t2 * t16 + t50;
	t55 = t2 * t10;
	t56 = -t12 * t16 - t55;
	t58 = t52 * t26 + t56 * t33;
	t61 = t56 * t26;
	t62 = -t52 * t33 + t61;
	t64 = t58 * t36 - t62 * t41;
	t65 = 0.1e1 / t64;
	t66 = t49 * t65;
	t67 = 0.1e1 / t44;
	t68 = t46 * t2;
	t69 = t68 * t10;
	t70 = t46 * t12;
	t72 = -t70 * t16 - t69;
	t75 = t70 * t10;
	t76 = t68 * t16 - t75;
	t78 = t72 * t26 + t76 * t33;
	t81 = t76 * t26;
	t82 = -t72 * t33 + t81;
	t84 = t78 * t36 - t82 * t41;
	t86 = t1 * t47;
	t87 = t84 * t44 + t86;
	t88 = t87 ^ 2;
	t89 = t64 ^ 2;
	t90 = 0.1e1 / t89;
	t92 = t44 ^ 2;
	t93 = 0.1e1 / t92;
	t96 = 0.1e1 / (t88 * t90 * t93 + 0.1e1);
	t97 = t67 * t96;
	t100 = t68 * t16 - t75;
	t102 = -t70 * t16;
	t103 = -t102 + t69;
	t111 = (t100 * t26 + t103 * t33) * t36 - (-t100 * t33 + t103 * t26) * t41;
	t115 = -t12 * t16 - t55;
	t117 = -t2 * t16;
	t118 = -t117 - t50;
	t126 = (t115 * t26 + t118 * t33) * t36 - (-t115 * t33 + t118 * t26) * t41;
	t129 = t87 * t90 * t96;
	t131 = -t111 * t65 * t96 + t126 * t67 * t129;
	t133 = t68 * t10 - t102;
	t137 = -t76 * t33;
	t141 = (t133 * t33 + t81) * t36 - (t133 * t26 + t137) * t41;
	t145 = -t12 * t10 - t117;
	t149 = -t56 * t33;
	t153 = (t145 * t33 + t61) * t36 - (t145 * t26 + t149) * t41;
	t156 = t153 * t67 * t129 - t141 * t65 * t96;
	t157 = t82 * t36;
	t159 = -t72 * t26 + t137;
	t160 = t159 * t41;
	t161 = t157 - t160;
	t164 = t62 * t36;
	t168 = t164 - (-t52 * t26 + t149) * t41;
	t171 = t168 * t67 * t129 - t161 * t65 * t96;
	t173 = t78 * t41 + t157;
	t177 = t58 * t41 + t164;
	t180 = t177 * t67 * t129 - t173 * t65 * t96;
	t182 = t1 * t44;
	t183 = t84 * t47 - t182;
	t190 = t65 * t47 * t87 * t93 * t96 - t183 * t65 * t97;
	t194 = -t76 * t26 + t72 * t33;
	t196 = t159 * t36 - t194 * t41;
	t199 = t64 * t44;
	t200 = atan2(t87, -t199);
	t201 = cos(t200);
	t202 = t201 * t64;
	t204 = sin(t200);
	t206 = -t202 * t44 + t204 * t87;
	t207 = 0.1e1 / t206;
	t209 = t49 ^ 2;
	t210 = t206 ^ 2;
	t211 = 0.1e1 / t210;
	t214 = 0.1e1 / (t209 * t211 + 0.1e1);
	t225 = t211 * t214;
	t229 = t3 * t16 - t29;
	t231 = -t13 * t16;
	t232 = -t231 + t11;
	t234 = t229 * t26 + t232 * t33;
	t238 = -t229 * t33 + t232 * t26;
	t240 = t234 * t36 - t238 * t41;
	t242 = t207 * t214;
	t257 = t3 * t10 - t231;
	t259 = t257 * t33 + t39;
	t261 = -t30 * t33;
	t263 = t257 * t26 + t261;
	t265 = t259 * t36 - t263 * t41;
	t280 = t40 * t36;
	t282 = -t18 * t26 + t261;
	t284 = -t282 * t41 + t280;
	t300 = t35 * t41 + t280;
	t317 = -t43 * t47 - t46 * t44;
	t331 = t196 * t47 + t182;
	t332 = cos(qJ(7));
	t335 = -t194 * t36 - t160;
	t336 = sin(qJ(7));
	t341 = -t300 * t332 - t317 * t336;
	t342 = 0.1e1 / t341;
	t346 = -t300 * t336 + t317 * t332;
	t347 = t346 ^ 2;
	t348 = t341 ^ 2;
	t349 = 0.1e1 / t348;
	t352 = 0.1e1 / (t347 * t349 + 0.1e1);
	t358 = t349 * t352;
	t361 = t240 * t47;
	t365 = -t234 * t41 - t238 * t36;
	t376 = t265 * t47;
	t380 = -t259 * t41 - t263 * t36;
	t391 = t284 * t47;
	t394 = -t282 * t36 - t42;
	t405 = t300 * t47;
	unknown(1,1) = -t66 * t97;
	unknown(1,2) = t131;
	unknown(1,3) = t156;
	unknown(1,4) = t171;
	unknown(1,5) = t180;
	unknown(1,6) = t190;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(1,9) = 0.0e0;
	unknown(1,10) = 0.0e0;
	unknown(1,11) = 0.0e0;
	unknown(1,12) = 0.0e0;
	unknown(1,13) = 0.0e0;
	unknown(2,1) = (-t196 * t44 + t86) * t207 * t214 + (-t66 * t67 * t96 * t201 * t87 - t49 * t96 * t204 + t204 * t49) * t49 * t225;
	unknown(2,2) = -t240 * t44 * t242 + (t204 * t111 * t44 - t201 * t126 * t44 + t131 * t204 * t199 + t131 * t201 * t87) * t49 * t225;
	unknown(2,3) = -t265 * t44 * t242 + (t204 * t141 * t44 - t201 * t153 * t44 + t156 * t204 * t199 + t156 * t201 * t87) * t49 * t225;
	unknown(2,4) = -t284 * t44 * t242 + (t204 * t161 * t44 - t201 * t168 * t44 + t171 * t204 * t199 + t171 * t201 * t87) * t49 * t225;
	unknown(2,5) = -t300 * t44 * t242 + (t204 * t173 * t44 - t201 * t177 * t44 + t180 * t204 * t199 + t180 * t201 * t87) * t49 * t225;
	unknown(2,6) = t317 * t207 * t214 + (t190 * t204 * t199 + t190 * t201 * t87 + t204 * t183 - t202 * t47) * t49 * t225;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(2,9) = 0.0e0;
	unknown(2,10) = 0.0e0;
	unknown(2,11) = 0.0e0;
	unknown(2,12) = 0.0e0;
	unknown(2,13) = 0.0e0;
	unknown(3,1) = (-t331 * t332 + t335 * t336) * t342 * t352 - (t331 * t336 + t335 * t332) * t346 * t358;
	unknown(3,2) = (-t361 * t332 + t365 * t336) * t342 * t352 - (t365 * t332 + t361 * t336) * t346 * t358;
	unknown(3,3) = (-t376 * t332 + t380 * t336) * t342 * t352 - (t380 * t332 + t376 * t336) * t346 * t358;
	unknown(3,4) = (-t391 * t332 + t394 * t336) * t342 * t352 - (t394 * t332 + t391 * t336) * t346 * t358;
	unknown(3,5) = (-t405 * t332 + t43 * t336) * t342 * t352 - (t43 * t332 + t405 * t336) * t346 * t358;
	unknown(3,6) = t49 * t336 * t346 * t349 * t352 + t49 * t332 * t342 * t352;
	unknown(3,7) = t346 ^ 2 * t358 + t352;
	unknown(3,8) = 0.0e0;
	unknown(3,9) = 0.0e0;
	unknown(3,10) = 0.0e0;
	unknown(3,11) = 0.0e0;
	unknown(3,12) = 0.0e0;
	unknown(3,13) = 0.0e0;
	Ja_rot = unknown;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:49
	% EndTime: 2020-04-11 13:05:50
	% DurationCPUTime: 1.14s
	% Computational Cost: add. (30624->165), mult. (49299->370), div. (205->9), fcn. (71398->23), ass. (0->210)
	unknown=NaN(3,13);
	t1 = cos(qJ(1));
	t2 = sin(qJ(2));
	t3 = t1 * t2;
	t4 = cos(pkin(14));
	t5 = cos(qJ(3));
	t7 = sin(pkin(14));
	t8 = sin(qJ(3));
	t10 = -t4 * t5 + t7 * t8;
	t11 = t3 * t10;
	t12 = cos(qJ(2));
	t13 = t1 * t12;
	t16 = -t4 * t8 - t7 * t5;
	t18 = -t13 * t16 - t11;
	t19 = cos(pkin(16));
	t20 = -qJ(4) + pkin(15);
	t21 = cos(t20);
	t23 = sin(pkin(16));
	t24 = sin(t20);
	t26 = -t19 * t21 - t23 * t24;
	t29 = t13 * t10;
	t30 = t3 * t16 - t29;
	t33 = t19 * t24 - t23 * t21;
	t35 = t18 * t26 + t30 * t33;
	t36 = sin(qJ(5));
	t39 = t30 * t26;
	t40 = -t18 * t33 + t39;
	t41 = cos(qJ(5));
	t42 = t40 * t41;
	t43 = t35 * t36 - t42;
	t44 = cos(qJ(6));
	t46 = sin(qJ(1));
	t47 = sin(qJ(6));
	t49 = t43 * t44 + t46 * t47;
	t50 = sin(qJ(7));
	t53 = t40 * t36;
	t54 = -t35 * t41 - t53;
	t55 = cos(qJ(7));
	t57 = -t49 * t50 - t54 * t55;
	t58 = t12 * t10;
	t60 = -t2 * t16 + t58;
	t63 = t2 * t10;
	t64 = -t12 * t16 - t63;
	t66 = t60 * t26 + t64 * t33;
	t69 = t64 * t26;
	t70 = -t60 * t33 + t69;
	t71 = t70 * t41;
	t72 = t66 * t36 - t71;
	t73 = t72 * t44;
	t76 = t70 * t36;
	t77 = -t66 * t41 - t76;
	t79 = t73 * t50 + t77 * t55;
	t80 = 0.1e1 / t79;
	t81 = t57 * t80;
	t82 = t46 * t2;
	t83 = t82 * t10;
	t84 = t46 * t12;
	t86 = -t84 * t16 - t83;
	t89 = t84 * t10;
	t90 = t82 * t16 - t89;
	t92 = t86 * t26 + t90 * t33;
	t95 = t90 * t26;
	t96 = -t86 * t33 + t95;
	t97 = t96 * t41;
	t98 = t92 * t36 - t97;
	t100 = t1 * t47;
	t101 = t98 * t44 - t100;
	t104 = t96 * t36;
	t105 = -t92 * t41 - t104;
	t107 = -t101 * t50 - t105 * t55;
	t108 = t107 ^ 2;
	t109 = t79 ^ 2;
	t110 = 0.1e1 / t109;
	t113 = 0.1e1 / (t108 * t110 + 0.1e1);
	t116 = t82 * t16 - t89;
	t118 = -t84 * t16;
	t119 = -t118 + t83;
	t121 = t116 * t26 + t119 * t33;
	t125 = -t116 * t33 + t119 * t26;
	t134 = -(t121 * t36 - t125 * t41) * t44 * t50 - (-t121 * t41 - t125 * t36) * t55;
	t138 = -t12 * t16 - t63;
	t140 = -t2 * t16;
	t141 = -t140 - t58;
	t143 = t138 * t26 + t141 * t33;
	t147 = -t138 * t33 + t141 * t26;
	t156 = (t143 * t36 - t147 * t41) * t44 * t50 + (-t143 * t41 - t147 * t36) * t55;
	t158 = t110 * t113;
	t160 = -t156 * t107 * t158 + t134 * t80 * t113;
	t162 = t82 * t10 - t118;
	t164 = t162 * t33 + t95;
	t166 = -t90 * t33;
	t168 = t162 * t26 + t166;
	t177 = -(t164 * t36 - t168 * t41) * t44 * t50 - (-t164 * t41 - t168 * t36) * t55;
	t181 = -t12 * t10 - t140;
	t183 = t181 * t33 + t69;
	t185 = -t64 * t33;
	t187 = t181 * t26 + t185;
	t196 = (t183 * t36 - t187 * t41) * t44 * t50 + (-t183 * t41 - t187 * t36) * t55;
	t199 = -t196 * t107 * t158 + t177 * t80 * t113;
	t201 = -t86 * t26 + t166;
	t202 = t201 * t41;
	t206 = t201 * t36;
	t209 = -(t104 - t202) * t44 * t50 - (-t97 - t206) * t55;
	t213 = -t60 * t26 + t185;
	t221 = (-t213 * t41 + t76) * t44 * t50 + (-t213 * t36 - t71) * t55;
	t224 = -t221 * t107 * t158 + t209 * t80 * t113;
	t228 = t105 * t44 * t50 - t98 * t55;
	t234 = -t77 * t44 * t50 + t72 * t55;
	t237 = -t234 * t107 * t158 + t228 * t80 * t113;
	t239 = t1 * t44;
	t240 = -t98 * t47 - t239;
	t249 = t72 * t47 * t50 * t107 * t110 * t113 - t240 * t50 * t80 * t113;
	t252 = -t101 * t55 + t105 * t50;
	t257 = -t77 * t50 + t73 * t55;
	t260 = -t257 * t107 * t158 + t252 * t80 * t113;
	t263 = -t90 * t26 + t86 * t33;
	t265 = -t263 * t41 + t206;
	t267 = t265 * t44 + t100;
	t270 = -t263 * t36 - t202;
	t273 = atan2(t107, t79);
	t274 = cos(t273);
	t276 = sin(t273);
	t278 = t276 * t107 + t274 * t79;
	t279 = 0.1e1 / t278;
	t281 = t57 ^ 2;
	t282 = t278 ^ 2;
	t283 = 0.1e1 / t282;
	t286 = 0.1e1 / (t281 * t283 + 0.1e1);
	t296 = t283 * t286;
	t300 = t3 * t16 - t29;
	t302 = -t13 * t16;
	t303 = -t302 + t11;
	t305 = t300 * t26 + t303 * t33;
	t309 = t303 * t26 - t300 * t33;
	t311 = t305 * t36 - t309 * t41;
	t312 = t311 * t44;
	t316 = -t305 * t41 - t309 * t36;
	t332 = t3 * t10 - t302;
	t334 = t332 * t33 + t39;
	t336 = -t30 * t33;
	t338 = t332 * t26 + t336;
	t340 = t334 * t36 - t338 * t41;
	t341 = t340 * t44;
	t345 = -t334 * t41 - t338 * t36;
	t361 = -t18 * t26 + t336;
	t363 = -t361 * t41 + t53;
	t364 = t363 * t44;
	t367 = -t361 * t36 - t42;
	t382 = -t54 * t44;
	t400 = -t43 * t47 + t46 * t44;
	t419 = t49 * t55 - t54 * t50;
	t434 = t267 * t55 - t270 * t50;
	t435 = sin(qJ(8));
	t438 = -t265 * t47 + t239;
	t439 = cos(qJ(8));
	t444 = t400 * t435 + t419 * t439;
	t445 = 0.1e1 / t444;
	t449 = -t400 * t439 + t419 * t435;
	t450 = t449 ^ 2;
	t451 = t444 ^ 2;
	t452 = 0.1e1 / t451;
	t455 = 0.1e1 / (t450 * t452 + 0.1e1);
	t461 = t452 * t455;
	t466 = t312 * t55 - t316 * t50;
	t468 = t311 * t47;
	t481 = t341 * t55 - t345 * t50;
	t483 = t340 * t47;
	t496 = t364 * t55 - t367 * t50;
	t498 = t363 * t47;
	t511 = t382 * t55 - t43 * t50;
	t513 = -t54 * t47;
	t524 = t400 * t55;
	unknown(1,1) = t81 * t113;
	unknown(1,2) = t160;
	unknown(1,3) = t199;
	unknown(1,4) = t224;
	unknown(1,5) = t237;
	unknown(1,6) = t249;
	unknown(1,7) = t260;
	unknown(1,8) = 0.0e0;
	unknown(1,9) = 0.0e0;
	unknown(1,10) = 0.0e0;
	unknown(1,11) = 0.0e0;
	unknown(1,12) = 0.0e0;
	unknown(1,13) = 0.0e0;
	unknown(2,1) = (t267 * t50 + t270 * t55) * t279 * t286 + (t81 * t113 * t274 * t107 - t57 * t113 * t276 + t276 * t57) * t57 * t296;
	unknown(2,2) = (t312 * t50 + t316 * t55) * t279 * t286 + (t160 * t274 * t107 - t160 * t276 * t79 + t276 * t134 + t274 * t156) * t57 * t296;
	unknown(2,3) = (t341 * t50 + t345 * t55) * t279 * t286 + (t199 * t274 * t107 - t199 * t276 * t79 + t276 * t177 + t274 * t196) * t57 * t296;
	unknown(2,4) = (t364 * t50 + t367 * t55) * t279 * t286 + (t224 * t274 * t107 - t224 * t276 * t79 + t276 * t209 + t274 * t221) * t57 * t296;
	unknown(2,5) = (t382 * t50 + t43 * t55) * t279 * t286 + (t237 * t274 * t107 - t237 * t276 * t79 + t276 * t228 + t274 * t234) * t57 * t296;
	unknown(2,6) = t400 * t50 * t279 * t286 + (-t274 * t72 * t47 * t50 + t249 * t274 * t107 - t276 * t240 * t50 - t249 * t276 * t79) * t57 * t296;
	unknown(2,7) = t419 * t279 * t286 + (t260 * t274 * t107 - t260 * t276 * t79 + t276 * t252 + t274 * t257) * t57 * t296;
	unknown(2,8) = 0.0e0;
	unknown(2,9) = 0.0e0;
	unknown(2,10) = 0.0e0;
	unknown(2,11) = 0.0e0;
	unknown(2,12) = 0.0e0;
	unknown(2,13) = 0.0e0;
	unknown(3,1) = (t434 * t435 - t438 * t439) * t445 * t455 - (t434 * t439 + t438 * t435) * t449 * t461;
	unknown(3,2) = (t466 * t435 + t468 * t439) * t445 * t455 - (-t468 * t435 + t466 * t439) * t449 * t461;
	unknown(3,3) = (t481 * t435 + t483 * t439) * t445 * t455 - (-t483 * t435 + t481 * t439) * t449 * t461;
	unknown(3,4) = (t496 * t435 + t498 * t439) * t445 * t455 - (-t498 * t435 + t496 * t439) * t449 * t461;
	unknown(3,5) = (t511 * t435 + t513 * t439) * t445 * t455 - (-t513 * t435 + t511 * t439) * t449 * t461;
	unknown(3,6) = (t524 * t435 + t49 * t439) * t445 * t455 - (-t49 * t435 + t524 * t439) * t449 * t461;
	unknown(3,7) = -t57 * t439 * t449 * t452 * t455 + t57 * t435 * t445 * t455;
	unknown(3,8) = t449 ^ 2 * t461 + t455;
	unknown(3,9) = 0.0e0;
	unknown(3,10) = 0.0e0;
	unknown(3,11) = 0.0e0;
	unknown(3,12) = 0.0e0;
	unknown(3,13) = 0.0e0;
	Ja_rot = unknown;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobia_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->39)
	unknown=NaN(3,13);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(1,7) = NaN;
	unknown(1,8) = NaN;
	unknown(1,9) = NaN;
	unknown(1,10) = NaN;
	unknown(1,11) = NaN;
	unknown(1,12) = NaN;
	unknown(1,13) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(2,7) = NaN;
	unknown(2,8) = NaN;
	unknown(2,9) = NaN;
	unknown(2,10) = NaN;
	unknown(2,11) = NaN;
	unknown(2,12) = NaN;
	unknown(2,13) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	unknown(3,7) = NaN;
	unknown(3,8) = NaN;
	unknown(3,9) = NaN;
	unknown(3,10) = NaN;
	unknown(3,11) = NaN;
	unknown(3,12) = NaN;
	unknown(3,13) = NaN;
	Ja_rot = unknown;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobia_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->39)
	unknown=NaN(3,13);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(1,7) = NaN;
	unknown(1,8) = NaN;
	unknown(1,9) = NaN;
	unknown(1,10) = NaN;
	unknown(1,11) = NaN;
	unknown(1,12) = NaN;
	unknown(1,13) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(2,7) = NaN;
	unknown(2,8) = NaN;
	unknown(2,9) = NaN;
	unknown(2,10) = NaN;
	unknown(2,11) = NaN;
	unknown(2,12) = NaN;
	unknown(2,13) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	unknown(3,7) = NaN;
	unknown(3,8) = NaN;
	unknown(3,9) = NaN;
	unknown(3,10) = NaN;
	unknown(3,11) = NaN;
	unknown(3,12) = NaN;
	unknown(3,13) = NaN;
	Ja_rot = unknown;
elseif link_index == 11
	%% Symbolic Calculation
	% From jacobia_rot_11_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->39)
	unknown=NaN(3,13);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(1,7) = NaN;
	unknown(1,8) = NaN;
	unknown(1,9) = NaN;
	unknown(1,10) = NaN;
	unknown(1,11) = NaN;
	unknown(1,12) = NaN;
	unknown(1,13) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(2,7) = NaN;
	unknown(2,8) = NaN;
	unknown(2,9) = NaN;
	unknown(2,10) = NaN;
	unknown(2,11) = NaN;
	unknown(2,12) = NaN;
	unknown(2,13) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	unknown(3,7) = NaN;
	unknown(3,8) = NaN;
	unknown(3,9) = NaN;
	unknown(3,10) = NaN;
	unknown(3,11) = NaN;
	unknown(3,12) = NaN;
	unknown(3,13) = NaN;
	Ja_rot = unknown;
elseif link_index == 12
	%% Symbolic Calculation
	% From jacobia_rot_12_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->39)
	unknown=NaN(3,13);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(1,7) = NaN;
	unknown(1,8) = NaN;
	unknown(1,9) = NaN;
	unknown(1,10) = NaN;
	unknown(1,11) = NaN;
	unknown(1,12) = NaN;
	unknown(1,13) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(2,7) = NaN;
	unknown(2,8) = NaN;
	unknown(2,9) = NaN;
	unknown(2,10) = NaN;
	unknown(2,11) = NaN;
	unknown(2,12) = NaN;
	unknown(2,13) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	unknown(3,7) = NaN;
	unknown(3,8) = NaN;
	unknown(3,9) = NaN;
	unknown(3,10) = NaN;
	unknown(3,11) = NaN;
	unknown(3,12) = NaN;
	unknown(3,13) = NaN;
	Ja_rot = unknown;
else
	Ja_rot=NaN(3,13);
end