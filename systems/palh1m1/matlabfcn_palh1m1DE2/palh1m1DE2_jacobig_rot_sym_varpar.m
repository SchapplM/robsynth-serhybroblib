% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh1m1DE2
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% Jg_rot [3x4]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = palh1m1DE2_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_jacobig_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1DE2_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_jacobig_rot_sym_varpar: pkin has to be [23x1] (double)');
Jg_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:34
	% EndTime: 2020-04-15 18:49:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:35
	% EndTime: 2020-04-15 18:49:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:34
	% EndTime: 2020-04-15 18:49:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)), 0, 0; 0, -cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:35
	% EndTime: 2020-04-15 18:49:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t18 = cos(qJ(1));
	t17 = sin(qJ(1));
	t1 = [0, t17, t17, 0; 0, -t18, -t18, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:50:09
	% EndTime: 2020-04-15 18:50:18
	% DurationCPUTime: 5.28s
	% Computational Cost: add. (92684->87), mult. (140076->161), div. (6088->11), fcn. (88848->16), ass. (0->104)
	t259 = -2 * pkin(1);
	t258 = -2 * pkin(5);
	t257 = pkin(4) * pkin(5);
	t256 = -pkin(8) - pkin(3);
	t255 = -pkin(8) + pkin(3);
	t254 = (-pkin(9) - pkin(11));
	t253 = (-pkin(9) + pkin(11));
	t196 = sin(pkin(21));
	t252 = t196 / 0.2e1;
	t201 = cos(qJ(3));
	t251 = -t201 / 0.2e1;
	t205 = pkin(7) ^ 2;
	t199 = sin(qJ(2));
	t203 = cos(pkin(19));
	t247 = sin(pkin(19));
	t248 = cos(qJ(2));
	t190 = t199 * t203 - t247 * t248;
	t244 = pkin(7) * t190;
	t228 = (pkin(1) ^ 2) + t244 * t259;
	t187 = t205 + t228;
	t184 = pkin(3) ^ 2 - pkin(8) ^ 2 + t187;
	t188 = pkin(1) - t244;
	t191 = t199 * t247 + t203 * t248;
	t182 = (pkin(7) - t256) * (pkin(7) + t256) + t228;
	t183 = (pkin(7) - t255) * (pkin(7) + t255) + t228;
	t211 = sqrt(-t183 * t182);
	t243 = pkin(7) * t191;
	t227 = pkin(1) * t243;
	t236 = 0.1e1 / t211 * (t182 + t183) * t227;
	t170 = (t190 * t211 + (t188 * t259 - t184 - t236) * t191) * pkin(7);
	t234 = t191 * t211;
	t171 = t188 * t236 + t205 * t191 ^ 2 * t259 + (-t190 * t184 - t234) * pkin(7);
	t185 = 0.1e1 / t187;
	t198 = sin(qJ(3));
	t208 = 0.1e1 / pkin(3);
	t220 = 0.1e1 / t187 ^ 2 * t227;
	t178 = t184 * t243 + t188 * t211;
	t230 = t201 * t178;
	t177 = -pkin(7) * t234 + t188 * t184;
	t233 = t198 * t177;
	t156 = ((-t198 * t170 / 0.2e1 + t171 * t251) * t185 + (-t230 - t233) * t220) * t208;
	t231 = t201 * t177;
	t232 = t198 * t178;
	t157 = ((t170 * t251 + t198 * t171 / 0.2e1) * t185 + (-t231 + t232) * t220) * t208;
	t195 = pkin(23) + pkin(22);
	t193 = sin(t195);
	t194 = cos(t195);
	t152 = t194 * t156 + t193 * t157;
	t235 = t185 * t208;
	t173 = (t233 / 0.2e1 + t230 / 0.2e1) * t235;
	t174 = (-t231 / 0.2e1 + t232 / 0.2e1) * t235;
	t166 = t194 * t173 - t193 * t174;
	t246 = pkin(4) * t166;
	t229 = (pkin(5) ^ 2) - t246 * t258;
	t158 = ((pkin(4) - t254) * (pkin(4) + t254)) + t229;
	t159 = ((pkin(4) - t253) * (pkin(4) + t253)) + t229;
	t216 = (t158 + t159) * t257;
	t148 = 0.2e1 * t152 * t216;
	t153 = -t193 * t156 + t194 * t157;
	t210 = sqrt(-t159 * t158);
	t207 = pkin(4) ^ 2;
	t163 = t207 + t229;
	t160 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t163;
	t164 = pkin(5) + t246;
	t221 = t164 * t258 - t160;
	t155 = 0.1e1 / t210;
	t217 = t193 * t173 + t194 * t174;
	t239 = t155 * t217;
	t139 = (-t153 * t210 - t148 * t239 / 0.2e1 + t221 * t152) * pkin(4);
	t222 = t207 * t217 * t258;
	t240 = t155 * t164;
	t140 = t148 * t240 / 0.2e1 + t152 * t222 + (-t152 * t210 + t153 * t160) * pkin(4);
	t245 = pkin(4) * t217;
	t149 = t164 * t160 - t210 * t245;
	t150 = t160 * t245 + t164 * t210;
	t197 = cos(pkin(21));
	t161 = 0.1e1 / t163;
	t204 = 0.1e1 / pkin(11);
	t237 = t161 * t204;
	t147 = (-t149 * t197 / 0.2e1 + t150 * t252) * t237;
	t144 = 0.1e1 / t147;
	t226 = 0.1e1 / t163 ^ 2 * t257;
	t218 = t197 * t226;
	t213 = t152 * t218;
	t219 = t196 * t226;
	t215 = t152 * t219;
	t238 = t161 * t197;
	t223 = t238 / 0.2e1;
	t224 = -t238 / 0.2e1;
	t225 = t161 * t252;
	t145 = 0.1e1 / t147 ^ 2;
	t146 = (t149 * t252 + t150 * t197 / 0.2e1) * t237;
	t241 = t145 * t146;
	t242 = 0.1e1 / (t145 * t146 ^ 2 + 0.1e1) * t204;
	t250 = ((t139 * t225 + t140 * t223 + t149 * t215 + t150 * t213) * t144 - (t139 * t224 + t140 * t225 - t149 * t213 + t150 * t215) * t241) * t242 + 0.1e1;
	t151 = t217 * t216;
	t142 = (-t151 * t239 - t166 * t210 + t217 * t221) * pkin(4);
	t143 = t151 * t240 + t217 * t222 + (t166 * t160 - t210 * t217) * pkin(4);
	t212 = t217 * t218;
	t214 = t217 * t219;
	t249 = ((t142 * t225 + t143 * t223 + t149 * t214 + t150 * t212) * t144 - (t142 * t224 + t143 * t225 - t149 * t212 + t150 * t214) * t241) * t242 + 0.1e1;
	t202 = cos(qJ(1));
	t200 = sin(qJ(1));
	t1 = [0, t250 * t200, t249 * t200, 0; 0, -t250 * t202, -t249 * t202, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:51:36
	% EndTime: 2020-04-15 18:51:45
	% DurationCPUTime: 5.75s
	% Computational Cost: add. (101085->90), mult. (152678->163), div. (6688->11), fcn. (96812->19), ass. (0->106)
	t438 = -2 * pkin(1);
	t437 = -2 * pkin(5);
	t436 = pkin(4) * pkin(5);
	t435 = -pkin(8) - pkin(3);
	t434 = -pkin(8) + pkin(3);
	t433 = (-pkin(9) - pkin(11));
	t432 = (-pkin(9) + pkin(11));
	t375 = sin(pkin(21));
	t431 = t375 / 0.2e1;
	t380 = cos(qJ(3));
	t430 = -t380 / 0.2e1;
	t384 = pkin(7) ^ 2;
	t378 = sin(qJ(2));
	t382 = cos(pkin(19));
	t426 = sin(pkin(19));
	t427 = cos(qJ(2));
	t369 = t378 * t382 - t426 * t427;
	t423 = pkin(7) * t369;
	t407 = (pkin(1) ^ 2) + t423 * t438;
	t366 = t384 + t407;
	t363 = pkin(3) ^ 2 - pkin(8) ^ 2 + t366;
	t367 = pkin(1) - t423;
	t370 = t378 * t426 + t382 * t427;
	t361 = (pkin(7) - t435) * (pkin(7) + t435) + t407;
	t362 = (pkin(7) - t434) * (pkin(7) + t434) + t407;
	t390 = sqrt(-t362 * t361);
	t422 = pkin(7) * t370;
	t406 = pkin(1) * t422;
	t415 = 0.1e1 / t390 * (t361 + t362) * t406;
	t349 = (t369 * t390 + (t367 * t438 - t363 - t415) * t370) * pkin(7);
	t413 = t370 * t390;
	t350 = t367 * t415 + t384 * t370 ^ 2 * t438 + (-t369 * t363 - t413) * pkin(7);
	t364 = 0.1e1 / t366;
	t377 = sin(qJ(3));
	t387 = 0.1e1 / pkin(3);
	t399 = 0.1e1 / t366 ^ 2 * t406;
	t357 = t363 * t422 + t367 * t390;
	t409 = t380 * t357;
	t356 = -pkin(7) * t413 + t367 * t363;
	t412 = t377 * t356;
	t335 = ((-t377 * t349 / 0.2e1 + t350 * t430) * t364 + (-t409 - t412) * t399) * t387;
	t410 = t380 * t356;
	t411 = t377 * t357;
	t336 = ((t349 * t430 + t377 * t350 / 0.2e1) * t364 + (-t410 + t411) * t399) * t387;
	t374 = pkin(23) + pkin(22);
	t372 = sin(t374);
	t373 = cos(t374);
	t331 = t373 * t335 + t372 * t336;
	t414 = t364 * t387;
	t352 = (t412 / 0.2e1 + t409 / 0.2e1) * t414;
	t353 = (-t410 / 0.2e1 + t411 / 0.2e1) * t414;
	t345 = t373 * t352 - t372 * t353;
	t425 = pkin(4) * t345;
	t408 = (pkin(5) ^ 2) - t425 * t437;
	t337 = ((pkin(4) - t433) * (pkin(4) + t433)) + t408;
	t338 = ((pkin(4) - t432) * (pkin(4) + t432)) + t408;
	t395 = (t337 + t338) * t436;
	t327 = 0.2e1 * t331 * t395;
	t332 = -t372 * t335 + t373 * t336;
	t389 = sqrt(-t338 * t337);
	t386 = pkin(4) ^ 2;
	t342 = t386 + t408;
	t339 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t342;
	t343 = pkin(5) + t425;
	t400 = t343 * t437 - t339;
	t334 = 0.1e1 / t389;
	t396 = t372 * t352 + t373 * t353;
	t418 = t334 * t396;
	t316 = (-t332 * t389 - t327 * t418 / 0.2e1 + t400 * t331) * pkin(4);
	t401 = t386 * t396 * t437;
	t419 = t334 * t343;
	t317 = t327 * t419 / 0.2e1 + t331 * t401 + (-t331 * t389 + t332 * t339) * pkin(4);
	t424 = pkin(4) * t396;
	t328 = t343 * t339 - t389 * t424;
	t329 = t339 * t424 + t343 * t389;
	t376 = cos(pkin(21));
	t340 = 0.1e1 / t342;
	t383 = 0.1e1 / pkin(11);
	t416 = t340 * t383;
	t326 = (-t328 * t376 / 0.2e1 + t329 * t431) * t416;
	t323 = 0.1e1 / t326;
	t405 = 0.1e1 / t342 ^ 2 * t436;
	t397 = t376 * t405;
	t392 = t331 * t397;
	t398 = t375 * t405;
	t394 = t331 * t398;
	t417 = t340 * t376;
	t402 = t417 / 0.2e1;
	t403 = -t417 / 0.2e1;
	t404 = t340 * t431;
	t324 = 0.1e1 / t326 ^ 2;
	t325 = (t328 * t431 + t329 * t376 / 0.2e1) * t416;
	t420 = t324 * t325;
	t421 = 0.1e1 / (t325 ^ 2 * t324 + 0.1e1) * t383;
	t429 = ((t316 * t404 + t317 * t402 + t328 * t394 + t329 * t392) * t323 - (t316 * t403 + t317 * t404 - t328 * t392 + t329 * t394) * t420) * t421 + 0.1e1;
	t330 = t396 * t395;
	t321 = (-t330 * t418 - t345 * t389 + t396 * t400) * pkin(4);
	t322 = t330 * t419 + t396 * t401 + (t345 * t339 - t389 * t396) * pkin(4);
	t391 = t396 * t397;
	t393 = t396 * t398;
	t428 = ((t321 * t404 + t322 * t402 + t328 * t393 + t329 * t391) * t323 - (t321 * t403 + t322 * t404 - t328 * t391 + t329 * t393) * t420) * t421 + 0.1e1;
	t381 = cos(qJ(1));
	t379 = sin(qJ(1));
	t319 = qJ(2) + qJ(3) + atan2(t325, t326);
	t318 = sin(t319);
	t1 = [0, t429 * t379, t428 * t379, t381 * t318; 0, -t429 * t381, -t428 * t381, t379 * t318; 1, 0, 0, -cos(t319);];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:36
	% EndTime: 2020-04-15 18:49:37
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (1521->36), mult. (2350->67), div. (100->7), fcn. (1510->10), ass. (0->42)
	t106 = -2 * pkin(7);
	t84 = cos(pkin(18));
	t105 = t84 / 0.2e1;
	t104 = (-pkin(8) - pkin(3));
	t103 = (-pkin(8) + pkin(3));
	t100 = cos(qJ(2));
	t81 = sin(qJ(2));
	t83 = cos(pkin(19));
	t99 = sin(pkin(19));
	t78 = -t100 * t99 + t81 * t83;
	t102 = pkin(1) * t78;
	t79 = t100 * t83 + t81 * t99;
	t101 = pkin(1) * t79;
	t87 = pkin(1) ^ 2;
	t94 = t102 * t106 + t87;
	t70 = ((pkin(7) - t104) * (pkin(7) + t104)) + t94;
	t71 = ((pkin(7) - t103) * (pkin(7) + t103)) + t94;
	t88 = sqrt(-t71 * t70);
	t93 = pkin(7) * t101;
	t98 = 0.1e1 / t88 * (t70 + t71) * t93;
	t75 = (pkin(7) ^ 2) + t94;
	t73 = 0.1e1 / t75;
	t82 = sin(pkin(18));
	t97 = t73 * t82;
	t85 = 0.1e1 / pkin(8);
	t96 = t73 * t85;
	t95 = t79 * t88;
	t92 = t73 * t105;
	t91 = 0.1e1 / t75 ^ 2 * t93;
	t90 = t82 * t91;
	t89 = t84 * t91;
	t76 = -pkin(7) + t102;
	t72 = -pkin(3) ^ 2 + pkin(8) ^ 2 + t75;
	t66 = t72 * t101 - t76 * t88;
	t65 = -pkin(1) * t95 - t76 * t72;
	t64 = (t66 * t105 + t65 * t82 / 0.2e1) * t96;
	t63 = (t65 * t105 - t82 * t66 / 0.2e1) * t96;
	t62 = 0.1e1 / t63 ^ 2;
	t61 = -t76 * t98 + t87 * t79 ^ 2 * t106 + (-t78 * t72 - t95) * pkin(1);
	t60 = (t78 * t88 + (0.2e1 * t76 * pkin(7) - t72 - t98) * t79) * pkin(1);
	t58 = ((t61 * t92 + t66 * t89 + t60 * t97 / 0.2e1 + t65 * t90) / t63 - (t60 * t92 + t65 * t89 - t61 * t97 / 0.2e1 - t66 * t90) * t64 * t62) / (t64 ^ 2 * t62 + 0.1e1) * t85;
	t1 = [0, t58 * sin(qJ(1)), 0, 0; 0, -t58 * cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobig_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:36
	% EndTime: 2020-04-15 18:49:37
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (1522->37), mult. (2350->66), div. (100->7), fcn. (1512->10), ass. (0->42)
	t111 = -2 * pkin(1);
	t83 = sin(pkin(23));
	t110 = t83 / 0.2e1;
	t109 = -pkin(8) - pkin(3);
	t108 = -pkin(8) + pkin(3);
	t103 = sin(pkin(19));
	t104 = cos(qJ(2));
	t85 = sin(qJ(2));
	t88 = cos(pkin(19));
	t80 = -t104 * t103 + t85 * t88;
	t106 = pkin(7) * t80;
	t98 = (pkin(1) ^ 2) + t106 * t111;
	t72 = (pkin(7) - t109) * (pkin(7) + t109) + t98;
	t73 = (pkin(7) - t108) * (pkin(7) + t108) + t98;
	t92 = sqrt(-t73 * t72);
	t81 = t85 * t103 + t104 * t88;
	t105 = pkin(7) * t81;
	t97 = pkin(1) * t105;
	t102 = 0.1e1 / t92 * (t72 + t73) * t97;
	t89 = pkin(7) ^ 2;
	t77 = t89 + t98;
	t74 = pkin(3) ^ 2 - pkin(8) ^ 2 + t77;
	t78 = pkin(1) - t106;
	t62 = (t80 * t92 + (t78 * t111 - t102 - t74) * t81) * pkin(7);
	t100 = t81 * t92;
	t63 = t78 * t102 + t89 * t81 ^ 2 * t111 + (-t80 * t74 - t100) * pkin(7);
	t75 = 0.1e1 / t77;
	t90 = 0.1e1 / pkin(3);
	t101 = t75 * t90;
	t67 = -pkin(7) * t100 + t78 * t74;
	t68 = t74 * t105 + t78 * t92;
	t84 = cos(pkin(23));
	t65 = (-t84 * t67 / 0.2e1 + t68 * t110) * t101;
	t64 = 0.1e1 / t65 ^ 2;
	t66 = (t84 * t68 / 0.2e1 + t67 * t110) * t101;
	t95 = 0.1e1 / t77 ^ 2 * t97;
	t93 = t84 * t95;
	t94 = t83 * t95;
	t96 = t75 * t110;
	t99 = t84 * t75;
	t107 = ((t63 * t99 / 0.2e1 + t68 * t93 + t62 * t96 + t67 * t94) / t65 - (-t62 * t99 / 0.2e1 - t67 * t93 + t63 * t96 + t68 * t94) * t66 * t64) / (t66 ^ 2 * t64 + 0.1e1) * t90 + 0.1e1;
	t1 = [0, t107 * sin(qJ(1)), 0, 0; 0, -t107 * cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobig_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:36
	% EndTime: 2020-04-15 18:49:36
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (726->31), mult. (1088->50), div. (28->6), fcn. (696->8), ass. (0->32)
	t73 = pkin(6) ^ 2;
	t67 = sin(pkin(20));
	t68 = cos(pkin(20));
	t69 = sin(qJ(3));
	t71 = cos(qJ(3));
	t65 = t71 * t67 + t69 * t68;
	t82 = pkin(6) * t65;
	t78 = 0.2e1 * pkin(1) * t82 + t73;
	t62 = pkin(1) ^ 2 + t78;
	t60 = 0.1e1 / t62;
	t87 = t60 / 0.2e1;
	t86 = -pkin(2) - pkin(13);
	t85 = -pkin(2) + pkin(13);
	t66 = t69 * t67 - t71 * t68;
	t84 = pkin(1) * t66;
	t59 = pkin(2) ^ 2 - pkin(13) ^ 2 + t62;
	t63 = -pkin(1) - t82;
	t57 = (pkin(1) - t86) * (pkin(1) + t86) + t78;
	t58 = (pkin(1) - t85) * (pkin(1) + t85) + t78;
	t76 = sqrt(-t58 * t57);
	t81 = t66 * pkin(6);
	t53 = t59 * t81 - t63 * t76;
	t83 = pkin(6) * t53;
	t80 = 0.1e1 / t76 * (t57 + t58) * pkin(1) * t81;
	t79 = t66 * t76;
	t72 = cos(qJ(1));
	t70 = sin(qJ(1));
	t61 = 0.1e1 / t62 ^ 2;
	t52 = -pkin(6) * t79 - t63 * t59;
	t51 = 0.1e1 / t52 ^ 2;
	t49 = 0.2e1 * (((-t63 * t80 + (t65 * t59 - t79) * pkin(6)) * t87 + (-t60 * t73 * t66 + t61 * t83) * t84) / t52 - ((-t65 * t76 + (-t59 - t80) * t66) * t87 + (t52 * t61 + t60 * t63) * t84) * t51 * t83) / (t53 ^ 2 * t51 + 0.1e1) * t62;
	t1 = [0, t70, t49 * t70, 0; 0, -t72, -t49 * t72, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobig_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:36
	% EndTime: 2020-04-15 18:49:37
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (925->37), mult. (1360->56), div. (40->9), fcn. (842->8), ass. (0->41)
	t102 = pkin(6) ^ 2;
	t111 = pkin(1) ^ 2 + t102;
	t95 = sin(pkin(20));
	t96 = cos(pkin(20));
	t97 = sin(qJ(3));
	t99 = cos(qJ(3));
	t93 = t95 * t99 + t96 * t97;
	t118 = t93 * pkin(6);
	t110 = pkin(1) * t118;
	t92 = 0.2e1 * t110;
	t90 = t92 + t111;
	t88 = 0.1e1 / t90;
	t123 = t88 / 0.2e1;
	t122 = -pkin(2) - pkin(13);
	t121 = -pkin(2) + pkin(13);
	t94 = t95 * t97 - t96 * t99;
	t120 = pkin(1) * t94;
	t112 = t102 + t92;
	t83 = (pkin(1) - t122) * (pkin(1) + t122) + t112;
	t84 = (pkin(1) - t121) * (pkin(1) + t121) + t112;
	t115 = t84 * t83;
	t106 = sqrt(-t115);
	t117 = t94 * pkin(6);
	t103 = pkin(2) ^ 2;
	t108 = -pkin(13) ^ 2 + t111;
	t87 = t103 + t92 + t108;
	t91 = -pkin(1) - t118;
	t78 = -t106 * t91 + t87 * t117;
	t119 = pkin(6) * t78;
	t116 = 0.1e1 / t106 * (t83 + t84) * pkin(1) * t117;
	t113 = t106 * t94;
	t109 = pkin(6) * t113;
	t77 = -t87 * t91 - t109;
	t76 = 0.1e1 / t77 ^ 2;
	t86 = t103 - t108 - 0.2e1 * t110;
	t85 = 0.1e1 / t86 ^ 2;
	t89 = 0.1e1 / t90 ^ 2;
	t114 = 0.2e1 * (((-t91 * t116 + (t87 * t93 - t113) * pkin(6)) * t123 + (-t102 * t88 * t94 + t89 * t119) * t120) / t77 - ((-t93 * t106 + (-t87 - t116) * t94) * t123 + (t77 * t89 + t88 * t91) * t120) * t76 * t119) / (t76 * t78 ^ 2 + 0.1e1) * t90 + (0.1e1 / t86 * t116 - 0.2e1 * pkin(1) * t85 * t109) / (-t85 * t115 + 0.1e1);
	t100 = cos(qJ(1));
	t98 = sin(qJ(1));
	t1 = [0, t98, t114 * t98, 0; 0, -t100, -t114 * t100, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobig_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:56
	% EndTime: 2020-04-15 18:50:05
	% DurationCPUTime: 2.84s
	% Computational Cost: add. (45553->86), mult. (68854->160), div. (2940->13), fcn. (43692->16), ass. (0->92)
	t206 = pkin(5) ^ 2;
	t205 = pkin(7) ^ 2;
	t199 = sin(qJ(2));
	t203 = cos(pkin(19));
	t242 = sin(pkin(19));
	t243 = cos(qJ(2));
	t190 = t199 * t203 - t242 * t243;
	t238 = pkin(7) * t190;
	t253 = -2 * pkin(1);
	t226 = (pkin(1) ^ 2) + t238 * t253;
	t187 = t205 + t226;
	t184 = pkin(3) ^ 2 - pkin(8) ^ 2 + t187;
	t188 = pkin(1) - t238;
	t252 = -pkin(8) - pkin(3);
	t182 = (pkin(7) - t252) * (pkin(7) + t252) + t226;
	t251 = -pkin(8) + pkin(3);
	t183 = (pkin(7) - t251) * (pkin(7) + t251) + t226;
	t211 = sqrt(-t183 * t182);
	t191 = t199 * t242 + t203 * t243;
	t237 = pkin(7) * t191;
	t178 = t184 * t237 + t188 * t211;
	t201 = cos(qJ(3));
	t228 = t201 * t178;
	t232 = t191 * t211;
	t177 = -pkin(7) * t232 + t184 * t188;
	t198 = sin(qJ(3));
	t231 = t198 * t177;
	t185 = 0.1e1 / t187;
	t208 = 0.1e1 / pkin(3);
	t233 = t185 * t208;
	t173 = (t231 / 0.2e1 + t228 / 0.2e1) * t233;
	t229 = t201 * t177;
	t230 = t198 * t178;
	t174 = (-t229 / 0.2e1 + t230 / 0.2e1) * t233;
	t195 = pkin(23) + pkin(22);
	t193 = sin(t195);
	t194 = cos(t195);
	t162 = t173 * t194 - t174 * t193;
	t240 = t162 * pkin(5);
	t227 = 0.2e1 * pkin(4) * t240 + t206;
	t250 = -pkin(9) - pkin(11);
	t154 = (pkin(4) - t250) * (pkin(4) + t250) + t227;
	t249 = -pkin(9) + pkin(11);
	t155 = (pkin(4) - t249) * (pkin(4) + t249) + t227;
	t210 = sqrt(-t155 * t154);
	t254 = 0.1e1 / t210 * pkin(4) * pkin(5) * (t154 + t155);
	t159 = pkin(4) ^ 2 + t227;
	t157 = 0.1e1 / t159;
	t248 = t157 / 0.2e1;
	t225 = pkin(1) * t237;
	t235 = 0.1e1 / t211 * (t182 + t183) * t225;
	t167 = (t190 * t211 + (t188 * t253 - t184 - t235) * t191) * pkin(7);
	t247 = -t167 / 0.2e1;
	t168 = t188 * t235 + t205 * t191 ^ 2 * t253 + (-t184 * t190 - t232) * pkin(7);
	t246 = t168 / 0.2e1;
	t196 = sin(pkin(23));
	t245 = t196 / 0.2e1;
	t244 = -t201 / 0.2e1;
	t156 = pkin(9) ^ 2 - pkin(11) ^ 2 + t159;
	t160 = -pkin(4) - t240;
	t220 = t173 * t193 + t194 * t174;
	t239 = pkin(5) * t220;
	t146 = t156 * t239 - t160 * t210;
	t241 = pkin(5) * t146;
	t236 = t220 * t254;
	t197 = cos(pkin(23));
	t234 = t185 * t197;
	t221 = 0.1e1 / t187 ^ 2 * t225;
	t152 = ((t168 * t244 + t198 * t247) * t185 + (-t228 - t231) * t221) * t208;
	t153 = ((t167 * t244 + t198 * t246) * t185 + (-t229 + t230) * t221) * t208;
	t148 = t152 * t194 + t153 * t193;
	t149 = -t152 * t193 + t153 * t194;
	t171 = (-t197 * t177 / 0.2e1 + t178 * t245) * t233;
	t170 = 0.1e1 / t171 ^ 2;
	t172 = (t197 * t178 / 0.2e1 + t177 * t245) * t233;
	t158 = 0.1e1 / t159 ^ 2;
	t212 = pkin(4) * (-t157 * t206 * t220 + t158 * t241);
	t145 = -t156 * t160 - t210 * t239;
	t144 = 0.1e1 / t145 ^ 2;
	t219 = 0.1e1 / (t144 * t146 ^ 2 + 0.1e1) * t159;
	t213 = t144 * t219 * t241;
	t214 = pkin(4) * (t145 * t158 + t157 * t160);
	t215 = 0.1e1 / t145 * t219;
	t216 = t178 * t221;
	t217 = t177 * t221;
	t222 = t185 * t245;
	t223 = t148 * t254;
	t224 = 0.2e1 * ((-t160 * t223 + (-t148 * t210 + t149 * t156) * pkin(5)) * t248 + t148 * t212) * t215 - 0.2e1 * ((-t148 * t156 - t149 * t210 - t220 * t223) * t248 + t148 * t214) * t213 + ((t167 * t222 + t196 * t217 + t197 * t216 + t234 * t246) / t171 - (t168 * t222 + t196 * t216 - t197 * t217 + t234 * t247) * t172 * t170) / (t170 * t172 ^ 2 + 0.1e1) * t208 + 0.1e1;
	t202 = cos(qJ(1));
	t200 = sin(qJ(1));
	t139 = 0.2e1 * ((-t160 * t236 + (t162 * t156 - t210 * t220) * pkin(5)) * t248 + t220 * t212) * t215 - 0.2e1 * (-t162 * t210 * t248 + ((-t156 - t236) * t248 + t214) * t220) * t213;
	t1 = [0, t224 * t200, t139 * t200, 0; 0, -t224 * t202, -t139 * t202, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
end