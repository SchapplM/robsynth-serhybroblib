% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh3m1DE2
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
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% Jg_rot [3x4]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = palh3m1DE2_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_jacobig_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1DE2_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_jacobig_rot_sym_varpar: pkin has to be [19x1] (double)');
Jg_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:41
	% EndTime: 2020-04-20 16:20:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:42
	% EndTime: 2020-04-20 16:20:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:42
	% EndTime: 2020-04-20 16:20:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)), 0, 0; 0, -cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:41
	% EndTime: 2020-04-20 16:20:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t25 = cos(qJ(1));
	t24 = sin(qJ(1));
	t1 = [0, t24, t24, 0; 0, -t25, -t25, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:21:19
	% EndTime: 2020-04-20 16:21:30
	% DurationCPUTime: 5.40s
	% Computational Cost: add. (92684->87), mult. (140076->159), div. (6088->11), fcn. (88848->16), ass. (0->104)
	t259 = -2 * pkin(1);
	t258 = -2 * pkin(4);
	t257 = pkin(3) * pkin(4);
	t256 = -pkin(6) - pkin(2);
	t255 = -pkin(6) + pkin(2);
	t254 = (-pkin(8) - pkin(10));
	t253 = (-pkin(8) + pkin(10));
	t196 = sin(pkin(17));
	t252 = t196 / 0.2e1;
	t198 = sin(qJ(3));
	t251 = t198 / 0.2e1;
	t205 = pkin(5) ^ 2;
	t199 = sin(qJ(2));
	t201 = sin(pkin(16));
	t247 = cos(qJ(2));
	t248 = cos(pkin(16));
	t190 = t199 * t201 - t247 * t248;
	t244 = pkin(5) * t190;
	t228 = (pkin(1) ^ 2) + t244 * t259;
	t187 = t205 + t228;
	t184 = pkin(2) ^ 2 - pkin(6) ^ 2 + t187;
	t188 = pkin(1) - t244;
	t191 = t199 * t248 + t247 * t201;
	t182 = (pkin(5) - t256) * (pkin(5) + t256) + t228;
	t183 = (pkin(5) - t255) * (pkin(5) + t255) + t228;
	t211 = sqrt(-t183 * t182);
	t243 = pkin(5) * t191;
	t227 = pkin(1) * t243;
	t232 = 0.1e1 / t211 * (t182 + t183) * t227;
	t171 = (t190 * t211 + (t188 * t259 - t184 - t232) * t191) * pkin(5);
	t230 = t191 * t211;
	t172 = t188 * t232 + t205 * t191 ^ 2 * t259 + (-t184 * t190 - t230) * pkin(5);
	t185 = 0.1e1 / t187;
	t202 = cos(qJ(3));
	t208 = 0.1e1 / pkin(2);
	t220 = 0.1e1 / t187 ^ 2 * t227;
	t178 = t184 * t243 + t188 * t211;
	t233 = t178 * t202;
	t177 = -pkin(5) * t230 + t184 * t188;
	t236 = t177 * t198;
	t157 = ((t172 * t202 / 0.2e1 + t171 * t251) * t185 + (t233 + t236) * t220) * t208;
	t234 = t178 * t198;
	t235 = t177 * t202;
	t158 = ((-t171 * t202 / 0.2e1 + t172 * t251) * t185 + (t234 - t235) * t220) * t208;
	t195 = pkin(18) + pkin(19);
	t193 = sin(t195);
	t194 = cos(t195);
	t154 = -t157 * t193 - t158 * t194;
	t231 = t185 * t208;
	t173 = (-t235 / 0.2e1 + t234 / 0.2e1) * t231;
	t174 = (t233 / 0.2e1 + t236 / 0.2e1) * t231;
	t168 = t173 * t194 + t174 * t193;
	t245 = pkin(3) * t168;
	t229 = (pkin(4) ^ 2) - t245 * t258;
	t159 = ((pkin(3) - t254) * (pkin(3) + t254)) + t229;
	t160 = ((pkin(3) - t253) * (pkin(3) + t253)) + t229;
	t217 = (t159 + t160) * t257;
	t149 = t154 * t217;
	t153 = -t157 * t194 + t158 * t193;
	t210 = sqrt(-t160 * t159);
	t207 = pkin(3) ^ 2;
	t164 = t207 + t229;
	t161 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t164;
	t165 = pkin(4) + t245;
	t221 = t165 * t258 - t161;
	t156 = 0.1e1 / t210;
	t212 = t173 * t193 - t174 * t194;
	t239 = t156 * t212;
	t140 = (-t149 * t239 - t153 * t210 + t221 * t154) * pkin(3);
	t222 = t212 * t207 * t258;
	t240 = t156 * t165;
	t141 = t149 * t240 + t154 * t222 + (t153 * t161 - t154 * t210) * pkin(3);
	t246 = pkin(3) * t212;
	t150 = t161 * t165 - t210 * t246;
	t151 = t161 * t246 + t165 * t210;
	t197 = cos(pkin(17));
	t162 = 0.1e1 / t164;
	t204 = 0.1e1 / pkin(10);
	t237 = t162 * t204;
	t147 = (-t150 * t197 / 0.2e1 + t151 * t252) * t237;
	t145 = 0.1e1 / t147;
	t226 = 0.1e1 / t164 ^ 2 * t257;
	t219 = t154 * t226;
	t214 = t151 * t219;
	t216 = t150 * t219;
	t238 = t162 * t197;
	t223 = t238 / 0.2e1;
	t224 = -t238 / 0.2e1;
	t225 = t162 * t252;
	t146 = 0.1e1 / t147 ^ 2;
	t148 = (t151 * t197 / 0.2e1 + t150 * t252) * t237;
	t241 = t146 * t148;
	t242 = 0.1e1 / (t146 * t148 ^ 2 + 0.1e1) * t204;
	t250 = ((t140 * t225 + t141 * t223 + t196 * t216 + t197 * t214) * t145 - (t140 * t224 + t141 * t225 + t196 * t214 - t197 * t216) * t241) * t242 + 0.1e1;
	t152 = t212 * t217;
	t143 = (-t152 * t239 - t168 * t210 + t212 * t221) * pkin(3);
	t144 = t152 * t240 + t212 * t222 + (t161 * t168 - t210 * t212) * pkin(3);
	t218 = t212 * t226;
	t213 = t151 * t218;
	t215 = t150 * t218;
	t249 = ((t143 * t225 + t144 * t223 + t196 * t215 + t197 * t213) * t145 - (t143 * t224 + t144 * t225 + t196 * t213 - t197 * t215) * t241) * t242 + 0.1e1;
	t203 = cos(qJ(1));
	t200 = sin(qJ(1));
	t1 = [0, t250 * t200, t249 * t200, 0; 0, -t250 * t203, -t249 * t203, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:23:20
	% EndTime: 2020-04-20 16:23:33
	% DurationCPUTime: 5.95s
	% Computational Cost: add. (101086->91), mult. (152678->161), div. (6688->11), fcn. (96812->19), ass. (0->106)
	t438 = -2 * pkin(1);
	t437 = -2 * pkin(4);
	t436 = pkin(3) * pkin(4);
	t435 = -pkin(6) - pkin(2);
	t434 = -pkin(6) + pkin(2);
	t433 = (-pkin(8) - pkin(10));
	t432 = (-pkin(8) + pkin(10));
	t375 = sin(pkin(17));
	t431 = t375 / 0.2e1;
	t377 = sin(qJ(3));
	t430 = t377 / 0.2e1;
	t384 = pkin(5) ^ 2;
	t378 = sin(qJ(2));
	t380 = sin(pkin(16));
	t426 = cos(qJ(2));
	t427 = cos(pkin(16));
	t369 = t378 * t380 - t426 * t427;
	t423 = pkin(5) * t369;
	t407 = (pkin(1) ^ 2) + t423 * t438;
	t366 = t384 + t407;
	t363 = pkin(2) ^ 2 - pkin(6) ^ 2 + t366;
	t367 = pkin(1) - t423;
	t370 = t378 * t427 + t380 * t426;
	t361 = (pkin(5) - t435) * (pkin(5) + t435) + t407;
	t362 = (pkin(5) - t434) * (pkin(5) + t434) + t407;
	t390 = sqrt(-t362 * t361);
	t422 = pkin(5) * t370;
	t406 = pkin(1) * t422;
	t411 = 0.1e1 / t390 * (t361 + t362) * t406;
	t350 = (t369 * t390 + (t367 * t438 - t363 - t411) * t370) * pkin(5);
	t409 = t370 * t390;
	t351 = t367 * t411 + t384 * t370 ^ 2 * t438 + (-t369 * t363 - t409) * pkin(5);
	t364 = 0.1e1 / t366;
	t381 = cos(qJ(3));
	t387 = 0.1e1 / pkin(2);
	t399 = 0.1e1 / t366 ^ 2 * t406;
	t357 = t363 * t422 + t367 * t390;
	t412 = t357 * t381;
	t356 = -pkin(5) * t409 + t367 * t363;
	t415 = t356 * t377;
	t336 = ((t351 * t381 / 0.2e1 + t350 * t430) * t364 + (t412 + t415) * t399) * t387;
	t413 = t357 * t377;
	t414 = t356 * t381;
	t337 = ((-t350 * t381 / 0.2e1 + t351 * t430) * t364 + (t413 - t414) * t399) * t387;
	t374 = pkin(18) + pkin(19);
	t372 = sin(t374);
	t373 = cos(t374);
	t333 = -t372 * t336 - t373 * t337;
	t410 = t364 * t387;
	t352 = (-t414 / 0.2e1 + t413 / 0.2e1) * t410;
	t353 = (t412 / 0.2e1 + t415 / 0.2e1) * t410;
	t347 = t373 * t352 + t372 * t353;
	t424 = pkin(3) * t347;
	t408 = (pkin(4) ^ 2) - t424 * t437;
	t338 = ((pkin(3) - t433) * (pkin(3) + t433)) + t408;
	t339 = ((pkin(3) - t432) * (pkin(3) + t432)) + t408;
	t396 = (t338 + t339) * t436;
	t328 = t333 * t396;
	t332 = -t373 * t336 + t372 * t337;
	t389 = sqrt(-t339 * t338);
	t386 = pkin(3) ^ 2;
	t343 = t386 + t408;
	t340 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t343;
	t344 = pkin(4) + t424;
	t400 = t344 * t437 - t340;
	t335 = 0.1e1 / t389;
	t391 = t372 * t352 - t373 * t353;
	t418 = t335 * t391;
	t317 = (-t328 * t418 - t332 * t389 + t333 * t400) * pkin(3);
	t401 = t386 * t391 * t437;
	t419 = t335 * t344;
	t318 = t328 * t419 + t333 * t401 + (t332 * t340 - t333 * t389) * pkin(3);
	t425 = pkin(3) * t391;
	t329 = t344 * t340 - t389 * t425;
	t330 = t340 * t425 + t344 * t389;
	t376 = cos(pkin(17));
	t341 = 0.1e1 / t343;
	t383 = 0.1e1 / pkin(10);
	t416 = t341 * t383;
	t326 = (-t329 * t376 / 0.2e1 + t330 * t431) * t416;
	t324 = 0.1e1 / t326;
	t405 = 0.1e1 / t343 ^ 2 * t436;
	t397 = t376 * t405;
	t393 = t333 * t397;
	t398 = t375 * t405;
	t395 = t333 * t398;
	t417 = t341 * t376;
	t402 = t417 / 0.2e1;
	t403 = -t417 / 0.2e1;
	t404 = t341 * t431;
	t325 = 0.1e1 / t326 ^ 2;
	t327 = (t330 * t376 / 0.2e1 + t329 * t431) * t416;
	t420 = t325 * t327;
	t421 = 0.1e1 / (t327 ^ 2 * t325 + 0.1e1) * t383;
	t429 = ((t317 * t404 + t318 * t402 + t329 * t395 + t330 * t393) * t324 - (t317 * t403 + t318 * t404 - t329 * t393 + t330 * t395) * t420) * t421 + 0.1e1;
	t331 = t391 * t396;
	t322 = (-t331 * t418 - t347 * t389 + t391 * t400) * pkin(3);
	t323 = t331 * t419 + t391 * t401 + (t347 * t340 - t389 * t391) * pkin(3);
	t392 = t391 * t397;
	t394 = t391 * t398;
	t428 = ((t322 * t404 + t323 * t402 + t329 * t394 + t330 * t392) * t324 - (t322 * t403 + t323 * t404 - t329 * t392 + t330 * t394) * t420) * t421 + 0.1e1;
	t382 = cos(qJ(1));
	t379 = sin(qJ(1));
	t320 = qJ(2) + qJ(3) + atan2(t327, t326);
	t319 = sin(t320);
	t1 = [0, t429 * t379, t428 * t379, -t382 * t319; 0, -t429 * t382, -t428 * t382, -t379 * t319; 1, 0, 0, cos(t320);];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:43
	% EndTime: 2020-04-20 16:20:44
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (1521->36), mult. (2350->67), div. (100->7), fcn. (1510->10), ass. (0->42)
	t106 = -2 * pkin(5);
	t84 = cos(pkin(15));
	t105 = t84 / 0.2e1;
	t104 = (-pkin(6) - pkin(2));
	t103 = (-pkin(6) + pkin(2));
	t100 = cos(pkin(16));
	t81 = sin(qJ(2));
	t82 = sin(pkin(16));
	t99 = cos(qJ(2));
	t78 = -t99 * t100 + t81 * t82;
	t102 = pkin(1) * t78;
	t79 = t81 * t100 + t99 * t82;
	t101 = pkin(1) * t79;
	t87 = pkin(1) ^ 2;
	t94 = t102 * t106 + t87;
	t70 = ((pkin(5) - t104) * (pkin(5) + t104)) + t94;
	t71 = ((pkin(5) - t103) * (pkin(5) + t103)) + t94;
	t88 = sqrt(-t71 * t70);
	t93 = pkin(5) * t101;
	t98 = 0.1e1 / t88 * (t70 + t71) * t93;
	t75 = (pkin(5) ^ 2) + t94;
	t73 = 0.1e1 / t75;
	t83 = sin(pkin(15));
	t97 = t73 * t83;
	t85 = 0.1e1 / pkin(6);
	t96 = t73 * t85;
	t95 = t79 * t88;
	t92 = t73 * t105;
	t91 = 0.1e1 / t75 ^ 2 * t93;
	t90 = t83 * t91;
	t89 = t84 * t91;
	t76 = -pkin(5) + t102;
	t72 = -pkin(2) ^ 2 + pkin(6) ^ 2 + t75;
	t66 = t72 * t101 - t76 * t88;
	t65 = -pkin(1) * t95 - t76 * t72;
	t64 = (t66 * t105 - t65 * t83 / 0.2e1) * t96;
	t63 = (t65 * t105 + t66 * t83 / 0.2e1) * t96;
	t62 = 0.1e1 / t63 ^ 2;
	t61 = -t76 * t98 + t87 * t79 ^ 2 * t106 + (-t78 * t72 - t95) * pkin(1);
	t60 = (t78 * t88 + (0.2e1 * t76 * pkin(5) - t72 - t98) * t79) * pkin(1);
	t58 = ((t61 * t92 + t66 * t89 - t60 * t97 / 0.2e1 - t65 * t90) / t63 - (t60 * t92 + t65 * t89 + t61 * t97 / 0.2e1 + t66 * t90) * t64 * t62) / (t64 ^ 2 * t62 + 0.1e1) * t85;
	t1 = [0, t58 * sin(qJ(1)), 0, 0; 0, -t58 * cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobig_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:43
	% EndTime: 2020-04-20 16:20:44
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (1522->37), mult. (2350->66), div. (100->7), fcn. (1512->10), ass. (0->42)
	t111 = -2 * pkin(1);
	t83 = sin(pkin(19));
	t110 = t83 / 0.2e1;
	t109 = -pkin(6) - pkin(2);
	t108 = -pkin(6) + pkin(2);
	t89 = pkin(5) ^ 2;
	t103 = cos(qJ(2));
	t104 = cos(pkin(16));
	t85 = sin(qJ(2));
	t87 = sin(pkin(16));
	t80 = -t103 * t104 + t85 * t87;
	t106 = pkin(5) * t80;
	t98 = (pkin(1) ^ 2) + t106 * t111;
	t77 = t89 + t98;
	t75 = 0.1e1 / t77;
	t84 = cos(pkin(19));
	t101 = t75 * t84;
	t72 = (pkin(5) - t109) * (pkin(5) + t109) + t98;
	t73 = (pkin(5) - t108) * (pkin(5) + t108) + t98;
	t92 = sqrt(-t73 * t72);
	t81 = t103 * t87 + t104 * t85;
	t105 = pkin(5) * t81;
	t97 = pkin(1) * t105;
	t102 = 0.1e1 / t92 * (t72 + t73) * t97;
	t74 = pkin(2) ^ 2 - pkin(6) ^ 2 + t77;
	t78 = pkin(1) - t106;
	t62 = (t80 * t92 + (t78 * t111 - t102 - t74) * t81) * pkin(5);
	t99 = t81 * t92;
	t63 = t78 * t102 + t89 * t81 ^ 2 * t111 + (-t80 * t74 - t99) * pkin(5);
	t90 = 0.1e1 / pkin(2);
	t100 = t75 * t90;
	t67 = -pkin(5) * t99 + t78 * t74;
	t68 = t105 * t74 + t78 * t92;
	t65 = (-t67 * t84 / 0.2e1 + t68 * t110) * t100;
	t64 = 0.1e1 / t65 ^ 2;
	t66 = (t68 * t84 / 0.2e1 + t67 * t110) * t100;
	t95 = 0.1e1 / t77 ^ 2 * t97;
	t93 = t84 * t95;
	t94 = t83 * t95;
	t96 = t75 * t110;
	t107 = ((t63 * t101 / 0.2e1 + t68 * t93 + t62 * t96 + t67 * t94) / t65 - (-t62 * t101 / 0.2e1 - t67 * t93 + t63 * t96 + t68 * t94) * t66 * t64) / (t66 ^ 2 * t64 + 0.1e1) * t90 + 0.1e1;
	t1 = [0, t107 * sin(qJ(1)), 0, 0; 0, -t107 * cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobig_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:21:03
	% EndTime: 2020-04-20 16:21:14
	% DurationCPUTime: 2.83s
	% Computational Cost: add. (45553->86), mult. (68854->163), div. (2940->13), fcn. (43692->16), ass. (0->95)
	t251 = -2 * pkin(1);
	t250 = -pkin(6) - pkin(2);
	t249 = -pkin(6) + pkin(2);
	t248 = -pkin(8) - pkin(10);
	t247 = -pkin(8) + pkin(10);
	t204 = pkin(4) ^ 2;
	t203 = pkin(5) ^ 2;
	t197 = sin(qJ(2));
	t199 = sin(pkin(16));
	t240 = cos(qJ(2));
	t241 = cos(pkin(16));
	t188 = t197 * t199 - t240 * t241;
	t236 = pkin(5) * t188;
	t223 = (pkin(1) ^ 2) + t236 * t251;
	t185 = t203 + t223;
	t183 = 0.1e1 / t185;
	t206 = 0.1e1 / pkin(2);
	t226 = t183 * t206;
	t182 = pkin(2) ^ 2 - pkin(6) ^ 2 + t185;
	t186 = pkin(1) - t236;
	t180 = (pkin(5) - t250) * (pkin(5) + t250) + t223;
	t181 = (pkin(5) - t249) * (pkin(5) + t249) + t223;
	t209 = sqrt(-t181 * t180);
	t189 = t197 * t241 + t240 * t199;
	t235 = pkin(5) * t189;
	t176 = t182 * t235 + t186 * t209;
	t196 = sin(qJ(3));
	t230 = t176 * t196;
	t225 = t189 * t209;
	t175 = -pkin(5) * t225 + t182 * t186;
	t200 = cos(qJ(3));
	t231 = t175 * t200;
	t171 = (-t231 / 0.2e1 + t230 / 0.2e1) * t226;
	t229 = t176 * t200;
	t232 = t175 * t196;
	t172 = (t229 / 0.2e1 + t232 / 0.2e1) * t226;
	t193 = pkin(18) + pkin(19);
	t191 = sin(t193);
	t192 = cos(t193);
	t162 = t171 * t192 + t172 * t191;
	t237 = t162 * pkin(4);
	t224 = 0.2e1 * pkin(3) * t237 + t204;
	t158 = pkin(3) ^ 2 + t224;
	t156 = 0.1e1 / t158;
	t246 = t156 / 0.2e1;
	t222 = pkin(1) * t235;
	t228 = 0.1e1 / t209 * (t180 + t181) * t222;
	t166 = (t188 * t209 + (t186 * t251 - t182 - t228) * t189) * pkin(5);
	t245 = -t166 / 0.2e1;
	t167 = t186 * t228 + t203 * t189 ^ 2 * t251 + (-t188 * t182 - t225) * pkin(5);
	t244 = t167 / 0.2e1;
	t194 = sin(pkin(19));
	t243 = t194 / 0.2e1;
	t242 = t196 / 0.2e1;
	t155 = pkin(8) ^ 2 - pkin(10) ^ 2 + t158;
	t159 = -pkin(3) - t237;
	t153 = (pkin(3) - t248) * (pkin(3) + t248) + t224;
	t154 = (pkin(3) - t247) * (pkin(3) + t247) + t224;
	t208 = sqrt(-t154 * t153);
	t213 = t171 * t191 - t172 * t192;
	t238 = pkin(4) * t213;
	t145 = t155 * t238 - t159 * t208;
	t239 = pkin(4) * t145;
	t150 = 0.1e1 / t208;
	t234 = t150 * t159;
	t233 = t150 * t213;
	t195 = cos(pkin(19));
	t227 = t183 * t195;
	t219 = 0.1e1 / t185 ^ 2 * t222;
	t151 = ((t166 * t242 + t200 * t244) * t183 + (t229 + t232) * t219) * t206;
	t152 = ((t167 * t242 + t200 * t245) * t183 + (t230 - t231) * t219) * t206;
	t148 = -t151 * t191 - t152 * t192;
	t217 = pkin(3) * pkin(4) * (t153 + t154);
	t140 = t148 * t217;
	t147 = -t151 * t192 + t152 * t191;
	t169 = (-t175 * t195 / 0.2e1 + t176 * t243) * t226;
	t168 = 0.1e1 / t169 ^ 2;
	t170 = (t176 * t195 / 0.2e1 + t175 * t243) * t226;
	t157 = 0.1e1 / t158 ^ 2;
	t210 = pkin(3) * (-t156 * t204 * t213 + t157 * t239);
	t144 = -t155 * t159 - t208 * t238;
	t143 = 0.1e1 / t144 ^ 2;
	t218 = 0.1e1 / (t143 * t145 ^ 2 + 0.1e1) * t158;
	t211 = t143 * t218 * t239;
	t212 = pkin(3) * (t144 * t157 + t156 * t159);
	t214 = 0.1e1 / t144 * t218;
	t215 = t176 * t219;
	t216 = t175 * t219;
	t220 = t183 * t243;
	t221 = 0.2e1 * ((-t140 * t234 + (t147 * t155 - t148 * t208) * pkin(4)) * t246 + t148 * t210) * t214 - 0.2e1 * ((-t140 * t233 - t147 * t208 - t148 * t155) * t246 + t148 * t212) * t211 + ((t166 * t220 + t194 * t216 + t195 * t215 + t227 * t244) / t169 - (t167 * t220 + t194 * t215 - t195 * t216 + t227 * t245) * t170 * t168) / (t168 * t170 ^ 2 + 0.1e1) * t206 + 0.1e1;
	t201 = cos(qJ(1));
	t198 = sin(qJ(1));
	t146 = t213 * t217;
	t138 = 0.2e1 * ((-t146 * t234 + (t162 * t155 - t208 * t213) * pkin(4)) * t246 + t213 * t210) * t214 - 0.2e1 * ((-t146 * t233 - t155 * t213 - t162 * t208) * t246 + t213 * t212) * t211;
	t1 = [0, t221 * t198, t138 * t198, 0; 0, -t221 * t201, -t138 * t201, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
end