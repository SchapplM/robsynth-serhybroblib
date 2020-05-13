% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh1m1TE
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
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = palh1m1TE_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_jacobig_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1TE_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_jacobig_rot_sym_varpar: pkin has to be [23x1] (double)');
Jg_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)), 0, 0; 0, -cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t26 = cos(qJ(1));
	t25 = sin(qJ(1));
	t1 = [0, t25, t25, 0; 0, -t26, -t26, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:18:11
	% EndTime: 2020-04-13 14:18:17
	% DurationCPUTime: 5.22s
	% Computational Cost: add. (92684->87), mult. (140076->161), div. (6088->11), fcn. (88848->16), ass. (0->104)
	t294 = -2 * pkin(1);
	t293 = -2 * pkin(5);
	t292 = pkin(4) * pkin(5);
	t291 = -pkin(8) - pkin(3);
	t290 = -pkin(8) + pkin(3);
	t289 = (-pkin(9) - pkin(11));
	t288 = (-pkin(9) + pkin(11));
	t231 = sin(pkin(21));
	t287 = t231 / 0.2e1;
	t236 = cos(qJ(3));
	t286 = -t236 / 0.2e1;
	t240 = pkin(7) ^ 2;
	t234 = sin(qJ(2));
	t238 = cos(pkin(19));
	t282 = sin(pkin(19));
	t283 = cos(qJ(2));
	t225 = t234 * t238 - t283 * t282;
	t279 = pkin(7) * t225;
	t263 = (pkin(1) ^ 2) + t279 * t294;
	t222 = t240 + t263;
	t219 = pkin(3) ^ 2 - pkin(8) ^ 2 + t222;
	t223 = pkin(1) - t279;
	t226 = t234 * t282 + t283 * t238;
	t217 = (pkin(7) - t291) * (pkin(7) + t291) + t263;
	t218 = (pkin(7) - t290) * (pkin(7) + t290) + t263;
	t246 = sqrt(-t218 * t217);
	t278 = pkin(7) * t226;
	t262 = pkin(1) * t278;
	t271 = 0.1e1 / t246 * (t217 + t218) * t262;
	t205 = (t225 * t246 + (t223 * t294 - t219 - t271) * t226) * pkin(7);
	t269 = t226 * t246;
	t206 = t223 * t271 + t240 * t226 ^ 2 * t294 + (-t225 * t219 - t269) * pkin(7);
	t220 = 0.1e1 / t222;
	t233 = sin(qJ(3));
	t243 = 0.1e1 / pkin(3);
	t255 = 0.1e1 / t222 ^ 2 * t262;
	t213 = t219 * t278 + t223 * t246;
	t265 = t236 * t213;
	t212 = -pkin(7) * t269 + t223 * t219;
	t268 = t233 * t212;
	t191 = ((-t233 * t205 / 0.2e1 + t206 * t286) * t220 + (-t265 - t268) * t255) * t243;
	t266 = t236 * t212;
	t267 = t233 * t213;
	t192 = ((t205 * t286 + t233 * t206 / 0.2e1) * t220 + (-t266 + t267) * t255) * t243;
	t230 = pkin(23) + pkin(22);
	t228 = sin(t230);
	t229 = cos(t230);
	t187 = t229 * t191 + t228 * t192;
	t270 = t220 * t243;
	t208 = (t268 / 0.2e1 + t265 / 0.2e1) * t270;
	t209 = (-t266 / 0.2e1 + t267 / 0.2e1) * t270;
	t201 = t229 * t208 - t228 * t209;
	t281 = pkin(4) * t201;
	t264 = (pkin(5) ^ 2) - t281 * t293;
	t193 = ((pkin(4) - t289) * (pkin(4) + t289)) + t264;
	t194 = ((pkin(4) - t288) * (pkin(4) + t288)) + t264;
	t251 = (t193 + t194) * t292;
	t183 = 0.2e1 * t187 * t251;
	t188 = -t228 * t191 + t229 * t192;
	t245 = sqrt(-t194 * t193);
	t242 = pkin(4) ^ 2;
	t198 = t242 + t264;
	t195 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t198;
	t199 = pkin(5) + t281;
	t256 = t199 * t293 - t195;
	t190 = 0.1e1 / t245;
	t252 = t228 * t208 + t229 * t209;
	t274 = t190 * t252;
	t174 = (-t188 * t245 - t183 * t274 / 0.2e1 + t256 * t187) * pkin(4);
	t257 = t242 * t252 * t293;
	t275 = t190 * t199;
	t175 = t183 * t275 / 0.2e1 + t187 * t257 + (-t187 * t245 + t188 * t195) * pkin(4);
	t280 = pkin(4) * t252;
	t184 = t199 * t195 - t245 * t280;
	t185 = t195 * t280 + t199 * t245;
	t232 = cos(pkin(21));
	t196 = 0.1e1 / t198;
	t239 = 0.1e1 / pkin(11);
	t272 = t196 * t239;
	t182 = (-t184 * t232 / 0.2e1 + t185 * t287) * t272;
	t179 = 0.1e1 / t182;
	t261 = 0.1e1 / t198 ^ 2 * t292;
	t253 = t232 * t261;
	t248 = t187 * t253;
	t254 = t231 * t261;
	t250 = t187 * t254;
	t273 = t196 * t232;
	t258 = t273 / 0.2e1;
	t259 = -t273 / 0.2e1;
	t260 = t196 * t287;
	t180 = 0.1e1 / t182 ^ 2;
	t181 = (t184 * t287 + t185 * t232 / 0.2e1) * t272;
	t276 = t180 * t181;
	t277 = 0.1e1 / (t181 ^ 2 * t180 + 0.1e1) * t239;
	t285 = ((t174 * t260 + t175 * t258 + t184 * t250 + t185 * t248) * t179 - (t174 * t259 + t175 * t260 - t184 * t248 + t185 * t250) * t276) * t277 + 0.1e1;
	t186 = t252 * t251;
	t177 = (-t186 * t274 - t201 * t245 + t252 * t256) * pkin(4);
	t178 = t186 * t275 + t252 * t257 + (t201 * t195 - t245 * t252) * pkin(4);
	t247 = t252 * t253;
	t249 = t252 * t254;
	t284 = ((t177 * t260 + t178 * t258 + t184 * t249 + t185 * t247) * t179 - (t177 * t259 + t178 * t260 - t184 * t247 + t185 * t249) * t276) * t277 + 0.1e1;
	t237 = cos(qJ(1));
	t235 = sin(qJ(1));
	t1 = [0, t285 * t235, t284 * t235, 0; 0, -t285 * t237, -t284 * t237, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:19:11
	% EndTime: 2020-04-13 14:19:16
	% DurationCPUTime: 5.67s
	% Computational Cost: add. (101087->91), mult. (152702->171), div. (6688->11), fcn. (96836->16), ass. (0->107)
	t493 = -2 * pkin(1);
	t492 = -2 * pkin(5);
	t491 = pkin(4) * pkin(5);
	t490 = -pkin(8) - pkin(3);
	t489 = -pkin(8) + pkin(3);
	t488 = (-pkin(9) - pkin(11));
	t487 = (-pkin(9) + pkin(11));
	t427 = sin(pkin(21));
	t486 = t427 / 0.2e1;
	t432 = cos(qJ(3));
	t485 = -t432 / 0.2e1;
	t437 = pkin(7) ^ 2;
	t430 = sin(qJ(2));
	t433 = cos(qJ(2));
	t435 = cos(pkin(19));
	t482 = sin(pkin(19));
	t421 = t430 * t435 - t433 * t482;
	t479 = pkin(7) * t421;
	t463 = (pkin(1) ^ 2) + t479 * t493;
	t418 = t437 + t463;
	t415 = pkin(3) ^ 2 - pkin(8) ^ 2 + t418;
	t419 = pkin(1) - t479;
	t422 = t430 * t482 + t433 * t435;
	t413 = (pkin(7) - t490) * (pkin(7) + t490) + t463;
	t414 = (pkin(7) - t489) * (pkin(7) + t489) + t463;
	t443 = sqrt(-t414 * t413);
	t478 = pkin(7) * t422;
	t462 = pkin(1) * t478;
	t471 = 0.1e1 / t443 * (t413 + t414) * t462;
	t401 = (t421 * t443 + (t419 * t493 - t415 - t471) * t422) * pkin(7);
	t469 = t422 * t443;
	t402 = t419 * t471 + t437 * t422 ^ 2 * t493 + (-t421 * t415 - t469) * pkin(7);
	t416 = 0.1e1 / t418;
	t429 = sin(qJ(3));
	t440 = 0.1e1 / pkin(3);
	t455 = 0.1e1 / t418 ^ 2 * t462;
	t409 = t415 * t478 + t419 * t443;
	t465 = t432 * t409;
	t408 = -pkin(7) * t469 + t419 * t415;
	t468 = t429 * t408;
	t387 = ((-t429 * t401 / 0.2e1 + t402 * t485) * t416 + (-t465 - t468) * t455) * t440;
	t466 = t432 * t408;
	t467 = t429 * t409;
	t388 = ((t401 * t485 + t429 * t402 / 0.2e1) * t416 + (-t466 + t467) * t455) * t440;
	t426 = pkin(23) + pkin(22);
	t424 = sin(t426);
	t425 = cos(t426);
	t383 = t425 * t387 + t424 * t388;
	t470 = t416 * t440;
	t404 = (t468 / 0.2e1 + t465 / 0.2e1) * t470;
	t405 = (-t466 / 0.2e1 + t467 / 0.2e1) * t470;
	t397 = t425 * t404 - t424 * t405;
	t481 = pkin(4) * t397;
	t464 = (pkin(5) ^ 2) - t481 * t492;
	t389 = ((pkin(4) - t488) * (pkin(4) + t488)) + t464;
	t390 = ((pkin(4) - t487) * (pkin(4) + t487)) + t464;
	t451 = (t389 + t390) * t491;
	t377 = 0.2e1 * t383 * t451;
	t384 = -t424 * t387 + t425 * t388;
	t442 = sqrt(-t390 * t389);
	t439 = pkin(4) ^ 2;
	t394 = t439 + t464;
	t391 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t394;
	t395 = pkin(5) + t481;
	t456 = t395 * t492 - t391;
	t386 = 0.1e1 / t442;
	t452 = t424 * t404 + t425 * t405;
	t474 = t386 * t452;
	t368 = (-t384 * t442 - t377 * t474 / 0.2e1 + t456 * t383) * pkin(4);
	t457 = t452 * t439 * t492;
	t475 = t386 * t395;
	t369 = t377 * t475 / 0.2e1 + t383 * t457 + (-t383 * t442 + t384 * t391) * pkin(4);
	t480 = pkin(4) * t452;
	t380 = t395 * t391 - t442 * t480;
	t381 = t391 * t480 + t395 * t442;
	t428 = cos(pkin(21));
	t392 = 0.1e1 / t394;
	t436 = 0.1e1 / pkin(11);
	t472 = t392 * t436;
	t376 = (-t380 * t428 / 0.2e1 + t381 * t486) * t472;
	t373 = 0.1e1 / t376;
	t461 = 0.1e1 / t394 ^ 2 * t491;
	t453 = t428 * t461;
	t448 = t383 * t453;
	t454 = t427 * t461;
	t450 = t383 * t454;
	t473 = t392 * t428;
	t458 = t473 / 0.2e1;
	t459 = -t473 / 0.2e1;
	t460 = t392 * t486;
	t374 = 0.1e1 / t376 ^ 2;
	t375 = (t380 * t486 + t381 * t428 / 0.2e1) * t472;
	t476 = t374 * t375;
	t477 = 0.1e1 / (t375 ^ 2 * t374 + 0.1e1) * t436;
	t484 = ((t368 * t460 + t369 * t458 + t380 * t450 + t381 * t448) * t373 - (t368 * t459 + t369 * t460 - t380 * t448 + t381 * t450) * t476) * t477 + 0.1e1;
	t382 = t452 * t451;
	t371 = (-t382 * t474 - t397 * t442 + t452 * t456) * pkin(4);
	t372 = t382 * t475 + t452 * t457 + (t397 * t391 - t442 * t452) * pkin(4);
	t447 = t452 * t453;
	t449 = t452 * t454;
	t483 = ((t371 * t460 + t372 * t458 + t380 * t449 + t381 * t447) * t373 - (t371 * t459 + t372 * t460 - t380 * t447 + t381 * t449) * t476) * t477 + 0.1e1;
	t446 = t433 * t429 + t430 * t432;
	t445 = -t430 * t429 + t433 * t432;
	t444 = t445 * t375 + t446 * t376;
	t434 = cos(qJ(1));
	t431 = sin(qJ(1));
	t1 = [0, t484 * t431, t483 * t431, t444 * t434; 0, -t484 * t434, -t483 * t434, t444 * t431; 1, 0, 0, t446 * t375 - t445 * t376;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:59
	% EndTime: 2020-04-13 14:17:59
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (1521->36), mult. (2350->67), div. (100->7), fcn. (1510->10), ass. (0->42)
	t102 = -2 * pkin(7);
	t80 = cos(pkin(18));
	t101 = t80 / 0.2e1;
	t100 = (-pkin(8) - pkin(3));
	t99 = (-pkin(8) + pkin(3));
	t77 = sin(qJ(2));
	t79 = cos(pkin(19));
	t95 = sin(pkin(19));
	t96 = cos(qJ(2));
	t74 = t77 * t79 - t96 * t95;
	t98 = pkin(1) * t74;
	t75 = t77 * t95 + t96 * t79;
	t97 = pkin(1) * t75;
	t83 = pkin(1) ^ 2;
	t90 = t98 * t102 + t83;
	t66 = ((pkin(7) - t100) * (pkin(7) + t100)) + t90;
	t67 = ((pkin(7) - t99) * (pkin(7) + t99)) + t90;
	t84 = sqrt(-t67 * t66);
	t89 = pkin(7) * t97;
	t94 = 0.1e1 / t84 * (t66 + t67) * t89;
	t71 = (pkin(7) ^ 2) + t90;
	t69 = 0.1e1 / t71;
	t78 = sin(pkin(18));
	t93 = t69 * t78;
	t81 = 0.1e1 / pkin(8);
	t92 = t69 * t81;
	t91 = t75 * t84;
	t88 = t69 * t101;
	t87 = 0.1e1 / t71 ^ 2 * t89;
	t86 = t78 * t87;
	t85 = t80 * t87;
	t72 = -pkin(7) + t98;
	t68 = -pkin(3) ^ 2 + pkin(8) ^ 2 + t71;
	t62 = t68 * t97 - t72 * t84;
	t61 = -pkin(1) * t91 - t72 * t68;
	t60 = (t62 * t101 + t61 * t78 / 0.2e1) * t92;
	t59 = (t61 * t101 - t78 * t62 / 0.2e1) * t92;
	t58 = 0.1e1 / t59 ^ 2;
	t57 = -t72 * t94 + t83 * t75 ^ 2 * t102 + (-t74 * t68 - t91) * pkin(1);
	t56 = (t74 * t84 + (0.2e1 * t72 * pkin(7) - t68 - t94) * t75) * pkin(1);
	t54 = ((t57 * t88 + t62 * t85 + t56 * t93 / 0.2e1 + t61 * t86) / t59 - (t56 * t88 + t61 * t85 - t57 * t93 / 0.2e1 - t62 * t86) * t60 * t58) / (t60 ^ 2 * t58 + 0.1e1) * t81;
	t1 = [0, t54 * sin(qJ(1)), 0, 0; 0, -t54 * cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobig_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:59
	% EndTime: 2020-04-13 14:18:00
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (1522->37), mult. (2350->66), div. (100->7), fcn. (1512->10), ass. (0->42)
	t117 = -2 * pkin(1);
	t89 = sin(pkin(23));
	t116 = t89 / 0.2e1;
	t115 = -pkin(8) - pkin(3);
	t114 = -pkin(8) + pkin(3);
	t109 = sin(pkin(19));
	t110 = cos(qJ(2));
	t91 = sin(qJ(2));
	t94 = cos(pkin(19));
	t87 = t91 * t109 + t110 * t94;
	t111 = pkin(7) * t87;
	t103 = pkin(1) * t111;
	t86 = -t110 * t109 + t91 * t94;
	t112 = pkin(7) * t86;
	t104 = (pkin(1) ^ 2) + t112 * t117;
	t95 = pkin(7) ^ 2;
	t83 = t95 + t104;
	t101 = 0.1e1 / t83 ^ 2 * t103;
	t100 = t89 * t101;
	t81 = 0.1e1 / t83;
	t102 = t81 * t116;
	t90 = cos(pkin(23));
	t105 = t90 * t81;
	t78 = (pkin(7) - t115) * (pkin(7) + t115) + t104;
	t79 = (pkin(7) - t114) * (pkin(7) + t114) + t104;
	t98 = sqrt(-t79 * t78);
	t108 = 0.1e1 / t98 * (t78 + t79) * t103;
	t80 = pkin(3) ^ 2 - pkin(8) ^ 2 + t83;
	t84 = pkin(1) - t112;
	t68 = (t86 * t98 + (t84 * t117 - t108 - t80) * t87) * pkin(7);
	t106 = t87 * t98;
	t69 = t84 * t108 + t95 * t87 ^ 2 * t117 + (-t86 * t80 - t106) * pkin(7);
	t96 = 0.1e1 / pkin(3);
	t107 = t81 * t96;
	t73 = -pkin(7) * t106 + t84 * t80;
	t74 = t80 * t111 + t84 * t98;
	t71 = (-t90 * t73 / 0.2e1 + t74 * t116) * t107;
	t70 = 0.1e1 / t71 ^ 2;
	t72 = (t90 * t74 / 0.2e1 + t73 * t116) * t107;
	t99 = t90 * t101;
	t113 = ((t69 * t105 / 0.2e1 + t74 * t99 + t68 * t102 + t73 * t100) / t71 - (-t68 * t105 / 0.2e1 - t73 * t99 + t69 * t102 + t74 * t100) * t72 * t70) / (t72 ^ 2 * t70 + 0.1e1) * t96 + 0.1e1;
	t1 = [0, t113 * sin(qJ(1)), 0, 0; 0, -t113 * cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobig_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:59
	% EndTime: 2020-04-13 14:17:59
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (726->31), mult. (1088->50), div. (28->6), fcn. (696->8), ass. (0->32)
	t83 = pkin(6) ^ 2;
	t77 = sin(pkin(20));
	t78 = cos(pkin(20));
	t79 = sin(qJ(3));
	t81 = cos(qJ(3));
	t75 = t81 * t77 + t79 * t78;
	t92 = pkin(6) * t75;
	t88 = 0.2e1 * pkin(1) * t92 + t83;
	t72 = pkin(1) ^ 2 + t88;
	t70 = 0.1e1 / t72;
	t97 = t70 / 0.2e1;
	t96 = -pkin(2) - pkin(13);
	t95 = -pkin(2) + pkin(13);
	t76 = t79 * t77 - t81 * t78;
	t94 = pkin(1) * t76;
	t69 = pkin(2) ^ 2 - pkin(13) ^ 2 + t72;
	t73 = -pkin(1) - t92;
	t67 = (pkin(1) - t96) * (pkin(1) + t96) + t88;
	t68 = (pkin(1) - t95) * (pkin(1) + t95) + t88;
	t86 = sqrt(-t68 * t67);
	t91 = t76 * pkin(6);
	t63 = t69 * t91 - t73 * t86;
	t93 = pkin(6) * t63;
	t90 = 0.1e1 / t86 * (t67 + t68) * pkin(1) * t91;
	t89 = t76 * t86;
	t82 = cos(qJ(1));
	t80 = sin(qJ(1));
	t71 = 0.1e1 / t72 ^ 2;
	t62 = -pkin(6) * t89 - t73 * t69;
	t61 = 0.1e1 / t62 ^ 2;
	t59 = 0.2e1 * (((-t73 * t90 + (t75 * t69 - t89) * pkin(6)) * t97 + (-t70 * t83 * t76 + t71 * t93) * t94) / t62 - ((-t75 * t86 + (-t69 - t90) * t76) * t97 + (t62 * t71 + t70 * t73) * t94) * t61 * t93) / (t63 ^ 2 * t61 + 0.1e1) * t72;
	t1 = [0, t80, t59 * t80, 0; 0, -t82, -t59 * t82, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobig_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:59
	% EndTime: 2020-04-13 14:18:00
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (925->37), mult. (1360->56), div. (40->9), fcn. (842->8), ass. (0->41)
	t144 = -pkin(2) - pkin(13);
	t143 = -pkin(2) + pkin(13);
	t116 = sin(pkin(20));
	t117 = cos(pkin(20));
	t118 = sin(qJ(3));
	t120 = cos(qJ(3));
	t114 = t120 * t116 + t118 * t117;
	t139 = pkin(6) * t114;
	t131 = pkin(1) * t139;
	t113 = 0.2e1 * t131;
	t124 = pkin(2) ^ 2;
	t123 = pkin(6) ^ 2;
	t132 = pkin(1) ^ 2 + t123;
	t129 = -pkin(13) ^ 2 + t132;
	t108 = t113 + t124 + t129;
	t112 = -pkin(1) - t139;
	t133 = t113 + t123;
	t104 = (pkin(1) - t144) * (pkin(1) + t144) + t133;
	t105 = (pkin(1) - t143) * (pkin(1) + t143) + t133;
	t135 = t105 * t104;
	t127 = sqrt(-t135);
	t115 = t118 * t116 - t120 * t117;
	t138 = t115 * pkin(6);
	t99 = t108 * t138 - t112 * t127;
	t142 = pkin(6) * t99;
	t111 = t113 + t132;
	t109 = 0.1e1 / t111;
	t141 = t109 / 0.2e1;
	t140 = pkin(1) * t115;
	t107 = t124 - t129 - 0.2e1 * t131;
	t106 = 0.1e1 / t107 ^ 2;
	t110 = 0.1e1 / t111 ^ 2;
	t134 = t115 * t127;
	t130 = pkin(6) * t134;
	t136 = 0.1e1 / t127 * (t104 + t105) * pkin(1) * t138;
	t98 = -t112 * t108 - t130;
	t97 = 0.1e1 / t98 ^ 2;
	t137 = 0.2e1 * (((-t112 * t136 + (t114 * t108 - t134) * pkin(6)) * t141 + (-t109 * t123 * t115 + t110 * t142) * t140) / t98 - ((-t114 * t127 + (-t108 - t136) * t115) * t141 + (t109 * t112 + t110 * t98) * t140) * t97 * t142) * t111 / (t99 ^ 2 * t97 + 0.1e1) + (0.1e1 / t107 * t136 - 0.2e1 * pkin(1) * t106 * t130) / (-t106 * t135 + 0.1e1);
	t121 = cos(qJ(1));
	t119 = sin(qJ(1));
	t1 = [0, t119, t137 * t119, 0; 0, -t121, -t137 * t121, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobig_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:18:16
	% EndTime: 2020-04-13 14:18:19
	% DurationCPUTime: 2.90s
	% Computational Cost: add. (45553->86), mult. (68854->160), div. (2940->13), fcn. (43692->16), ass. (0->92)
	t231 = pkin(5) ^ 2;
	t230 = pkin(7) ^ 2;
	t224 = sin(qJ(2));
	t228 = cos(pkin(19));
	t267 = sin(pkin(19));
	t268 = cos(qJ(2));
	t215 = t224 * t228 - t268 * t267;
	t263 = pkin(7) * t215;
	t278 = -2 * pkin(1);
	t251 = (pkin(1) ^ 2) + t263 * t278;
	t212 = t230 + t251;
	t209 = pkin(3) ^ 2 - pkin(8) ^ 2 + t212;
	t213 = pkin(1) - t263;
	t277 = -pkin(8) - pkin(3);
	t207 = (pkin(7) - t277) * (pkin(7) + t277) + t251;
	t276 = -pkin(8) + pkin(3);
	t208 = (pkin(7) - t276) * (pkin(7) + t276) + t251;
	t236 = sqrt(-t208 * t207);
	t216 = t224 * t267 + t268 * t228;
	t262 = pkin(7) * t216;
	t203 = t209 * t262 + t213 * t236;
	t226 = cos(qJ(3));
	t253 = t226 * t203;
	t258 = t216 * t236;
	t202 = -pkin(7) * t258 + t213 * t209;
	t223 = sin(qJ(3));
	t256 = t223 * t202;
	t210 = 0.1e1 / t212;
	t233 = 0.1e1 / pkin(3);
	t259 = t210 * t233;
	t198 = (t256 / 0.2e1 + t253 / 0.2e1) * t259;
	t254 = t226 * t202;
	t255 = t223 * t203;
	t199 = (-t254 / 0.2e1 + t255 / 0.2e1) * t259;
	t220 = pkin(23) + pkin(22);
	t218 = sin(t220);
	t219 = cos(t220);
	t187 = t219 * t198 - t218 * t199;
	t265 = t187 * pkin(5);
	t252 = 0.2e1 * pkin(4) * t265 + t231;
	t275 = -pkin(9) - pkin(11);
	t179 = (pkin(4) - t275) * (pkin(4) + t275) + t252;
	t274 = -pkin(9) + pkin(11);
	t180 = (pkin(4) - t274) * (pkin(4) + t274) + t252;
	t235 = sqrt(-t180 * t179);
	t279 = 0.1e1 / t235 * pkin(4) * pkin(5) * (t179 + t180);
	t184 = pkin(4) ^ 2 + t252;
	t182 = 0.1e1 / t184;
	t273 = t182 / 0.2e1;
	t250 = pkin(1) * t262;
	t260 = 0.1e1 / t236 * (t207 + t208) * t250;
	t192 = (t215 * t236 + (t213 * t278 - t209 - t260) * t216) * pkin(7);
	t272 = -t192 / 0.2e1;
	t193 = t213 * t260 + t230 * t216 ^ 2 * t278 + (-t215 * t209 - t258) * pkin(7);
	t271 = t193 / 0.2e1;
	t221 = sin(pkin(23));
	t270 = t221 / 0.2e1;
	t269 = -t226 / 0.2e1;
	t181 = pkin(9) ^ 2 - pkin(11) ^ 2 + t184;
	t185 = -pkin(4) - t265;
	t245 = t218 * t198 + t219 * t199;
	t264 = pkin(5) * t245;
	t171 = t181 * t264 - t185 * t235;
	t266 = pkin(5) * t171;
	t261 = t245 * t279;
	t222 = cos(pkin(23));
	t257 = t222 * t210;
	t246 = 0.1e1 / t212 ^ 2 * t250;
	t177 = ((t193 * t269 + t223 * t272) * t210 + (-t253 - t256) * t246) * t233;
	t178 = ((t192 * t269 + t223 * t271) * t210 + (-t254 + t255) * t246) * t233;
	t173 = t219 * t177 + t218 * t178;
	t174 = -t218 * t177 + t219 * t178;
	t196 = (-t222 * t202 / 0.2e1 + t203 * t270) * t259;
	t195 = 0.1e1 / t196 ^ 2;
	t197 = (t222 * t203 / 0.2e1 + t202 * t270) * t259;
	t183 = 0.1e1 / t184 ^ 2;
	t237 = pkin(4) * (-t182 * t231 * t245 + t183 * t266);
	t170 = -t185 * t181 - t235 * t264;
	t169 = 0.1e1 / t170 ^ 2;
	t244 = 0.1e1 / (t171 ^ 2 * t169 + 0.1e1) * t184;
	t238 = t169 * t244 * t266;
	t239 = pkin(4) * (t170 * t183 + t182 * t185);
	t240 = 0.1e1 / t170 * t244;
	t241 = t222 * t246;
	t242 = t221 * t246;
	t247 = t210 * t270;
	t248 = t173 * t279;
	t249 = 0.2e1 * ((-t185 * t248 + (-t173 * t235 + t174 * t181) * pkin(5)) * t273 + t173 * t237) * t240 - 0.2e1 * ((-t173 * t181 - t174 * t235 - t245 * t248) * t273 + t173 * t239) * t238 + ((t192 * t247 + t202 * t242 + t203 * t241 + t257 * t271) / t196 - (t193 * t247 - t202 * t241 + t203 * t242 + t257 * t272) * t197 * t195) / (t197 ^ 2 * t195 + 0.1e1) * t233 + 0.1e1;
	t227 = cos(qJ(1));
	t225 = sin(qJ(1));
	t164 = 0.2e1 * ((-t185 * t261 + (t187 * t181 - t235 * t245) * pkin(5)) * t273 + t245 * t237) * t240 - 0.2e1 * (-t187 * t235 * t273 + ((-t181 - t261) * t273 + t239) * t245) * t238;
	t1 = [0, t249 * t225, t164 * t225, 0; 0, -t249 * t227, -t164 * t227, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
end