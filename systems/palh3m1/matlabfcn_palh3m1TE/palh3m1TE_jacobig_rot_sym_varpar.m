% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh3m1TE
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
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = palh3m1TE_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_jacobig_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1TE_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_jacobig_rot_sym_varpar: pkin has to be [19x1] (double)');
Jg_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)), 0, 0; 0, -cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
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
	% StartTime: 2020-04-18 09:52:20
	% EndTime: 2020-04-18 09:52:29
	% DurationCPUTime: 5.38s
	% Computational Cost: add. (92684->87), mult. (140076->159), div. (6088->11), fcn. (88848->16), ass. (0->104)
	t289 = -2 * pkin(1);
	t288 = -2 * pkin(4);
	t287 = pkin(3) * pkin(4);
	t286 = -pkin(6) - pkin(2);
	t285 = -pkin(6) + pkin(2);
	t284 = (-pkin(8) - pkin(10));
	t283 = (-pkin(8) + pkin(10));
	t226 = sin(pkin(17));
	t282 = t226 / 0.2e1;
	t228 = sin(qJ(3));
	t281 = t228 / 0.2e1;
	t235 = pkin(5) ^ 2;
	t229 = sin(qJ(2));
	t231 = sin(pkin(16));
	t277 = cos(qJ(2));
	t278 = cos(pkin(16));
	t220 = t229 * t231 - t277 * t278;
	t274 = pkin(5) * t220;
	t258 = (pkin(1) ^ 2) + t274 * t289;
	t217 = t235 + t258;
	t214 = pkin(2) ^ 2 - pkin(6) ^ 2 + t217;
	t218 = pkin(1) - t274;
	t221 = t229 * t278 + t277 * t231;
	t212 = (pkin(5) - t286) * (pkin(5) + t286) + t258;
	t213 = (pkin(5) - t285) * (pkin(5) + t285) + t258;
	t241 = sqrt(-t213 * t212);
	t273 = pkin(5) * t221;
	t257 = pkin(1) * t273;
	t262 = 0.1e1 / t241 * (t212 + t213) * t257;
	t201 = (t220 * t241 + (t218 * t289 - t214 - t262) * t221) * pkin(5);
	t260 = t221 * t241;
	t202 = t218 * t262 + t235 * t221 ^ 2 * t289 + (-t220 * t214 - t260) * pkin(5);
	t215 = 0.1e1 / t217;
	t232 = cos(qJ(3));
	t238 = 0.1e1 / pkin(2);
	t250 = 0.1e1 / t217 ^ 2 * t257;
	t208 = t214 * t273 + t218 * t241;
	t263 = t208 * t232;
	t207 = -pkin(5) * t260 + t214 * t218;
	t266 = t207 * t228;
	t187 = ((t202 * t232 / 0.2e1 + t201 * t281) * t215 + (t263 + t266) * t250) * t238;
	t264 = t208 * t228;
	t265 = t207 * t232;
	t188 = ((-t201 * t232 / 0.2e1 + t202 * t281) * t215 + (t264 - t265) * t250) * t238;
	t225 = pkin(18) + pkin(19);
	t223 = sin(t225);
	t224 = cos(t225);
	t184 = -t187 * t223 - t188 * t224;
	t261 = t215 * t238;
	t203 = (-t265 / 0.2e1 + t264 / 0.2e1) * t261;
	t204 = (t263 / 0.2e1 + t266 / 0.2e1) * t261;
	t198 = t203 * t224 + t204 * t223;
	t275 = pkin(3) * t198;
	t259 = (pkin(4) ^ 2) - t275 * t288;
	t189 = ((pkin(3) - t284) * (pkin(3) + t284)) + t259;
	t190 = ((pkin(3) - t283) * (pkin(3) + t283)) + t259;
	t247 = (t189 + t190) * t287;
	t179 = t184 * t247;
	t183 = -t187 * t224 + t188 * t223;
	t240 = sqrt(-t190 * t189);
	t237 = pkin(3) ^ 2;
	t194 = t237 + t259;
	t191 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t194;
	t195 = pkin(4) + t275;
	t251 = t195 * t288 - t191;
	t186 = 0.1e1 / t240;
	t242 = t203 * t223 - t204 * t224;
	t269 = t186 * t242;
	t170 = (-t179 * t269 - t183 * t240 + t251 * t184) * pkin(3);
	t252 = t237 * t242 * t288;
	t270 = t186 * t195;
	t171 = t179 * t270 + t184 * t252 + (t183 * t191 - t184 * t240) * pkin(3);
	t276 = pkin(3) * t242;
	t180 = t191 * t195 - t240 * t276;
	t181 = t191 * t276 + t195 * t240;
	t227 = cos(pkin(17));
	t192 = 0.1e1 / t194;
	t234 = 0.1e1 / pkin(10);
	t267 = t192 * t234;
	t177 = (-t180 * t227 / 0.2e1 + t181 * t282) * t267;
	t175 = 0.1e1 / t177;
	t256 = 0.1e1 / t194 ^ 2 * t287;
	t249 = t184 * t256;
	t245 = t181 * t249;
	t246 = t180 * t249;
	t268 = t192 * t227;
	t253 = t268 / 0.2e1;
	t254 = -t268 / 0.2e1;
	t255 = t192 * t282;
	t176 = 0.1e1 / t177 ^ 2;
	t178 = (t181 * t227 / 0.2e1 + t180 * t282) * t267;
	t271 = t176 * t178;
	t272 = 0.1e1 / (t176 * t178 ^ 2 + 0.1e1) * t234;
	t280 = ((t170 * t255 + t171 * t253 + t226 * t246 + t227 * t245) * t175 - (t170 * t254 + t171 * t255 + t226 * t245 - t227 * t246) * t271) * t272 + 0.1e1;
	t182 = t242 * t247;
	t173 = (-t182 * t269 - t198 * t240 + t242 * t251) * pkin(3);
	t174 = t182 * t270 + t242 * t252 + (t198 * t191 - t240 * t242) * pkin(3);
	t248 = t242 * t256;
	t243 = t227 * t248;
	t244 = t226 * t248;
	t279 = ((t173 * t255 + t174 * t253 + t180 * t244 + t181 * t243) * t175 - (t173 * t254 + t174 * t255 - t180 * t243 + t181 * t244) * t271) * t272 + 0.1e1;
	t233 = cos(qJ(1));
	t230 = sin(qJ(1));
	t1 = [0, t280 * t230, t279 * t230, 0; 0, -t280 * t233, -t279 * t233, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:53:57
	% EndTime: 2020-04-18 09:54:06
	% DurationCPUTime: 5.80s
	% Computational Cost: add. (101087->88), mult. (152702->165), div. (6688->11), fcn. (96836->16), ass. (0->104)
	t441 = pkin(3) ^ 2;
	t439 = pkin(5) ^ 2;
	t432 = sin(qJ(2));
	t434 = sin(pkin(16));
	t436 = cos(qJ(2));
	t484 = cos(pkin(16));
	t423 = t432 * t434 - t436 * t484;
	t481 = pkin(5) * t423;
	t495 = -2 * pkin(1);
	t465 = (pkin(1) ^ 2) + t481 * t495;
	t420 = t439 + t465;
	t418 = 0.1e1 / t420;
	t442 = 0.1e1 / pkin(2);
	t468 = t418 * t442;
	t417 = pkin(2) ^ 2 - pkin(6) ^ 2 + t420;
	t421 = pkin(1) - t481;
	t492 = -pkin(6) - pkin(2);
	t415 = (pkin(5) - t492) * (pkin(5) + t492) + t465;
	t491 = -pkin(6) + pkin(2);
	t416 = (pkin(5) - t491) * (pkin(5) + t491) + t465;
	t445 = sqrt(-t416 * t415);
	t424 = t432 * t484 + t436 * t434;
	t480 = pkin(5) * t424;
	t411 = t417 * t480 + t421 * t445;
	t431 = sin(qJ(3));
	t471 = t411 * t431;
	t467 = t424 * t445;
	t410 = -pkin(5) * t467 + t421 * t417;
	t435 = cos(qJ(3));
	t472 = t410 * t435;
	t406 = (-t472 / 0.2e1 + t471 / 0.2e1) * t468;
	t470 = t411 * t435;
	t473 = t410 * t431;
	t407 = (t470 / 0.2e1 + t473 / 0.2e1) * t468;
	t428 = pkin(18) + pkin(19);
	t426 = sin(t428);
	t427 = cos(t428);
	t401 = t406 * t427 + t407 * t426;
	t482 = pkin(3) * t401;
	t494 = -2 * pkin(4);
	t466 = (pkin(4) ^ 2) - t482 * t494;
	t397 = t441 + t466;
	t394 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t397;
	t398 = pkin(4) + t482;
	t490 = -pkin(8) - pkin(10);
	t392 = (pkin(3) - t490) * (pkin(3) + t490) + t466;
	t489 = -pkin(8) + pkin(10);
	t393 = (pkin(3) - t489) * (pkin(3) + t489) + t466;
	t444 = sqrt(-t393 * t392);
	t449 = t406 * t426 - t407 * t427;
	t483 = pkin(3) * t449;
	t383 = t394 * t398 - t444 * t483;
	t384 = t394 * t483 + t398 * t444;
	t430 = cos(pkin(17));
	t395 = 0.1e1 / t397;
	t438 = 0.1e1 / pkin(10);
	t474 = t395 * t438;
	t429 = sin(pkin(17));
	t488 = t429 / 0.2e1;
	t378 = (-t383 * t430 / 0.2e1 + t384 * t488) * t474;
	t376 = 0.1e1 / t378;
	t493 = pkin(3) * pkin(4);
	t463 = 0.1e1 / t397 ^ 2 * t493;
	t455 = t430 * t463;
	t456 = t429 * t463;
	t377 = 0.1e1 / t378 ^ 2;
	t379 = (t384 * t430 / 0.2e1 + t383 * t488) * t474;
	t478 = t377 * t379;
	t496 = (t383 * t456 + t384 * t455) * t376 - (-t383 * t455 + t384 * t456) * t478;
	t487 = t431 / 0.2e1;
	t464 = pkin(1) * t480;
	t469 = 0.1e1 / t445 * (t415 + t416) * t464;
	t404 = (t423 * t445 + (t421 * t495 - t417 - t469) * t424) * pkin(5);
	t405 = t421 * t469 + t439 * t424 ^ 2 * t495 + (-t417 * t423 - t467) * pkin(5);
	t457 = 0.1e1 / t420 ^ 2 * t464;
	t390 = ((t405 * t435 / 0.2e1 + t404 * t487) * t418 + (t470 + t473) * t457) * t442;
	t391 = ((-t404 * t435 / 0.2e1 + t405 * t487) * t418 + (t471 - t472) * t457) * t442;
	t387 = -t390 * t426 - t391 * t427;
	t454 = (t392 + t393) * t493;
	t380 = t387 * t454;
	t386 = -t390 * t427 + t391 * t426;
	t458 = t398 * t494 - t394;
	t389 = 0.1e1 / t444;
	t476 = t389 * t449;
	t371 = (-t380 * t476 - t386 * t444 + t458 * t387) * pkin(3);
	t459 = t449 * t441 * t494;
	t477 = t389 * t398;
	t372 = t380 * t477 + t387 * t459 + (t386 * t394 - t387 * t444) * pkin(3);
	t475 = t395 * t430;
	t460 = t475 / 0.2e1;
	t461 = -t475 / 0.2e1;
	t462 = t395 * t488;
	t479 = 0.1e1 / (t377 * t379 ^ 2 + 0.1e1) * t438;
	t486 = ((t371 * t462 + t372 * t460) * t376 - (t371 * t461 + t372 * t462) * t478 + t496 * t387) * t479 + 0.1e1;
	t385 = t449 * t454;
	t374 = (-t385 * t476 - t401 * t444 + t449 * t458) * pkin(3);
	t375 = t385 * t477 + t449 * t459 + (t394 * t401 - t444 * t449) * pkin(3);
	t485 = ((t374 * t462 + t375 * t460) * t376 - (t374 * t461 + t375 * t462) * t478 + t496 * t449) * t479 + 0.1e1;
	t448 = t431 * t436 + t432 * t435;
	t447 = t431 * t432 - t435 * t436;
	t446 = -t448 * t378 + t447 * t379;
	t437 = cos(qJ(1));
	t433 = sin(qJ(1));
	t1 = [0, t486 * t433, t485 * t433, t446 * t437; 0, -t486 * t437, -t485 * t437, t446 * t433; 1, 0, 0, -t447 * t378 - t448 * t379;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:03
	% EndTime: 2020-04-18 09:52:04
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (1521->36), mult. (2350->67), div. (100->7), fcn. (1510->10), ass. (0->42)
	t103 = -2 * pkin(5);
	t81 = cos(pkin(15));
	t102 = t81 / 0.2e1;
	t101 = (-pkin(6) - pkin(2));
	t100 = (-pkin(6) + pkin(2));
	t78 = sin(qJ(2));
	t79 = sin(pkin(16));
	t96 = cos(qJ(2));
	t97 = cos(pkin(16));
	t75 = t78 * t79 - t96 * t97;
	t99 = pkin(1) * t75;
	t76 = t78 * t97 + t96 * t79;
	t98 = pkin(1) * t76;
	t84 = pkin(1) ^ 2;
	t91 = t99 * t103 + t84;
	t67 = ((pkin(5) - t101) * (pkin(5) + t101)) + t91;
	t68 = ((pkin(5) - t100) * (pkin(5) + t100)) + t91;
	t85 = sqrt(-t68 * t67);
	t90 = pkin(5) * t98;
	t95 = 0.1e1 / t85 * (t67 + t68) * t90;
	t72 = (pkin(5) ^ 2) + t91;
	t70 = 0.1e1 / t72;
	t80 = sin(pkin(15));
	t94 = t70 * t80;
	t82 = 0.1e1 / pkin(6);
	t93 = t70 * t82;
	t92 = t76 * t85;
	t89 = t70 * t102;
	t88 = 0.1e1 / t72 ^ 2 * t90;
	t87 = t80 * t88;
	t86 = t81 * t88;
	t73 = -pkin(5) + t99;
	t69 = -pkin(2) ^ 2 + pkin(6) ^ 2 + t72;
	t63 = t69 * t98 - t73 * t85;
	t62 = -pkin(1) * t92 - t73 * t69;
	t61 = (t63 * t102 - t62 * t80 / 0.2e1) * t93;
	t60 = (t62 * t102 + t63 * t80 / 0.2e1) * t93;
	t59 = 0.1e1 / t60 ^ 2;
	t58 = -t73 * t95 + t84 * t76 ^ 2 * t103 + (-t75 * t69 - t92) * pkin(1);
	t57 = (t75 * t85 + (0.2e1 * t73 * pkin(5) - t69 - t95) * t76) * pkin(1);
	t55 = ((t58 * t89 + t63 * t86 - t57 * t94 / 0.2e1 - t62 * t87) / t60 - (t57 * t89 + t62 * t86 + t58 * t94 / 0.2e1 + t63 * t87) * t61 * t59) / (t61 ^ 2 * t59 + 0.1e1) * t82;
	t1 = [0, t55 * sin(qJ(1)), 0, 0; 0, -t55 * cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobig_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:04
	% EndTime: 2020-04-18 09:52:04
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (1522->37), mult. (2350->66), div. (100->7), fcn. (1512->10), ass. (0->42)
	t113 = -2 * pkin(1);
	t85 = sin(pkin(19));
	t112 = t85 / 0.2e1;
	t111 = -pkin(6) - pkin(2);
	t110 = -pkin(6) + pkin(2);
	t105 = cos(qJ(2));
	t106 = cos(pkin(16));
	t87 = sin(qJ(2));
	t89 = sin(pkin(16));
	t82 = -t105 * t106 + t87 * t89;
	t108 = pkin(5) * t82;
	t100 = (pkin(1) ^ 2) + t108 * t113;
	t91 = pkin(5) ^ 2;
	t79 = t91 + t100;
	t77 = 0.1e1 / t79;
	t86 = cos(pkin(19));
	t103 = t77 * t86;
	t74 = (pkin(5) - t111) * (pkin(5) + t111) + t100;
	t75 = (pkin(5) - t110) * (pkin(5) + t110) + t100;
	t94 = sqrt(-t75 * t74);
	t83 = t105 * t89 + t87 * t106;
	t107 = pkin(5) * t83;
	t99 = pkin(1) * t107;
	t104 = 0.1e1 / t94 * (t74 + t75) * t99;
	t76 = pkin(2) ^ 2 - pkin(6) ^ 2 + t79;
	t80 = pkin(1) - t108;
	t64 = (t82 * t94 + (t80 * t113 - t104 - t76) * t83) * pkin(5);
	t101 = t83 * t94;
	t65 = t80 * t104 + t91 * t83 ^ 2 * t113 + (-t82 * t76 - t101) * pkin(5);
	t92 = 0.1e1 / pkin(2);
	t102 = t77 * t92;
	t69 = -pkin(5) * t101 + t80 * t76;
	t70 = t76 * t107 + t80 * t94;
	t67 = (-t69 * t86 / 0.2e1 + t70 * t112) * t102;
	t66 = 0.1e1 / t67 ^ 2;
	t68 = (t70 * t86 / 0.2e1 + t69 * t112) * t102;
	t97 = 0.1e1 / t79 ^ 2 * t99;
	t95 = t86 * t97;
	t96 = t85 * t97;
	t98 = t77 * t112;
	t109 = ((t65 * t103 / 0.2e1 + t70 * t95 + t64 * t98 + t69 * t96) / t67 - (-t64 * t103 / 0.2e1 - t69 * t95 + t65 * t98 + t70 * t96) * t68 * t66) / (t68 ^ 2 * t66 + 0.1e1) * t92 + 0.1e1;
	t1 = [0, t109 * sin(qJ(1)), 0, 0; 0, -t109 * cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobig_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:26
	% EndTime: 2020-04-18 09:52:31
	% DurationCPUTime: 2.84s
	% Computational Cost: add. (45553->86), mult. (68854->163), div. (2940->13), fcn. (43692->16), ass. (0->95)
	t282 = -2 * pkin(1);
	t281 = -pkin(6) - pkin(2);
	t280 = -pkin(6) + pkin(2);
	t279 = -pkin(8) - pkin(10);
	t278 = -pkin(8) + pkin(10);
	t235 = pkin(4) ^ 2;
	t234 = pkin(5) ^ 2;
	t228 = sin(qJ(2));
	t230 = sin(pkin(16));
	t271 = cos(qJ(2));
	t272 = cos(pkin(16));
	t219 = t228 * t230 - t271 * t272;
	t267 = pkin(5) * t219;
	t254 = (pkin(1) ^ 2) + t267 * t282;
	t216 = t234 + t254;
	t214 = 0.1e1 / t216;
	t237 = 0.1e1 / pkin(2);
	t257 = t214 * t237;
	t213 = pkin(2) ^ 2 - pkin(6) ^ 2 + t216;
	t217 = pkin(1) - t267;
	t211 = (pkin(5) - t281) * (pkin(5) + t281) + t254;
	t212 = (pkin(5) - t280) * (pkin(5) + t280) + t254;
	t240 = sqrt(-t212 * t211);
	t220 = t228 * t272 + t271 * t230;
	t266 = pkin(5) * t220;
	t207 = t213 * t266 + t217 * t240;
	t227 = sin(qJ(3));
	t261 = t207 * t227;
	t256 = t220 * t240;
	t206 = -pkin(5) * t256 + t217 * t213;
	t231 = cos(qJ(3));
	t262 = t206 * t231;
	t202 = (-t262 / 0.2e1 + t261 / 0.2e1) * t257;
	t260 = t207 * t231;
	t263 = t206 * t227;
	t203 = (t260 / 0.2e1 + t263 / 0.2e1) * t257;
	t224 = pkin(18) + pkin(19);
	t222 = sin(t224);
	t223 = cos(t224);
	t193 = t223 * t202 + t222 * t203;
	t268 = t193 * pkin(4);
	t255 = 0.2e1 * pkin(3) * t268 + t235;
	t189 = pkin(3) ^ 2 + t255;
	t187 = 0.1e1 / t189;
	t277 = t187 / 0.2e1;
	t253 = pkin(1) * t266;
	t259 = 0.1e1 / t240 * (t211 + t212) * t253;
	t197 = (t219 * t240 + (t217 * t282 - t213 - t259) * t220) * pkin(5);
	t276 = -t197 / 0.2e1;
	t198 = t217 * t259 + t234 * t220 ^ 2 * t282 + (-t219 * t213 - t256) * pkin(5);
	t275 = t198 / 0.2e1;
	t225 = sin(pkin(19));
	t274 = t225 / 0.2e1;
	t273 = t227 / 0.2e1;
	t186 = pkin(8) ^ 2 - pkin(10) ^ 2 + t189;
	t190 = -pkin(3) - t268;
	t184 = (pkin(3) - t279) * (pkin(3) + t279) + t255;
	t185 = (pkin(3) - t278) * (pkin(3) + t278) + t255;
	t239 = sqrt(-t185 * t184);
	t244 = t222 * t202 - t223 * t203;
	t269 = pkin(4) * t244;
	t176 = t186 * t269 - t190 * t239;
	t270 = pkin(4) * t176;
	t181 = 0.1e1 / t239;
	t265 = t181 * t190;
	t264 = t181 * t244;
	t226 = cos(pkin(19));
	t258 = t214 * t226;
	t250 = 0.1e1 / t216 ^ 2 * t253;
	t182 = ((t197 * t273 + t231 * t275) * t214 + (t260 + t263) * t250) * t237;
	t183 = ((t198 * t273 + t231 * t276) * t214 + (t261 - t262) * t250) * t237;
	t179 = -t222 * t182 - t223 * t183;
	t248 = pkin(3) * pkin(4) * (t184 + t185);
	t171 = t179 * t248;
	t178 = -t223 * t182 + t222 * t183;
	t200 = (-t206 * t226 / 0.2e1 + t207 * t274) * t257;
	t199 = 0.1e1 / t200 ^ 2;
	t201 = (t207 * t226 / 0.2e1 + t206 * t274) * t257;
	t188 = 0.1e1 / t189 ^ 2;
	t241 = pkin(3) * (-t187 * t235 * t244 + t188 * t270);
	t175 = -t190 * t186 - t239 * t269;
	t174 = 0.1e1 / t175 ^ 2;
	t249 = 0.1e1 / (t176 ^ 2 * t174 + 0.1e1) * t189;
	t242 = t174 * t249 * t270;
	t243 = pkin(3) * (t175 * t188 + t187 * t190);
	t245 = 0.1e1 / t175 * t249;
	t246 = t226 * t250;
	t247 = t225 * t250;
	t251 = t214 * t274;
	t252 = 0.2e1 * ((-t171 * t265 + (t178 * t186 - t179 * t239) * pkin(4)) * t277 + t179 * t241) * t245 - 0.2e1 * ((-t171 * t264 - t178 * t239 - t179 * t186) * t277 + t179 * t243) * t242 + ((t197 * t251 + t206 * t247 + t207 * t246 + t258 * t275) / t200 - (t198 * t251 - t206 * t246 + t207 * t247 + t258 * t276) * t201 * t199) / (t201 ^ 2 * t199 + 0.1e1) * t237 + 0.1e1;
	t232 = cos(qJ(1));
	t229 = sin(qJ(1));
	t177 = t244 * t248;
	t169 = 0.2e1 * ((-t177 * t265 + (t193 * t186 - t239 * t244) * pkin(4)) * t277 + t244 * t241) * t245 - 0.2e1 * ((-t177 * t264 - t186 * t244 - t193 * t239) * t277 + t244 * t243) * t242;
	t1 = [0, t252 * t229, t169 * t229, 0; 0, -t252 * t232, -t169 * t232, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
end