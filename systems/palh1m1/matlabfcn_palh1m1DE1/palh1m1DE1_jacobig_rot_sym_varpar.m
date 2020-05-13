% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh1m1DE1
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
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = palh1m1DE1_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_jacobig_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1DE1_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_jacobig_rot_sym_varpar: pkin has to be [23x1] (double)');
Jg_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:57
	% EndTime: 2020-04-14 18:42:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:58
	% EndTime: 2020-04-14 18:42:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:58
	% EndTime: 2020-04-14 18:42:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)), 0, 0; 0, -cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:58
	% EndTime: 2020-04-14 18:42:58
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
	% StartTime: 2020-04-14 18:44:06
	% EndTime: 2020-04-14 18:44:14
	% DurationCPUTime: 5.27s
	% Computational Cost: add. (92684->87), mult. (140076->161), div. (6088->11), fcn. (88848->16), ass. (0->104)
	t292 = -2 * pkin(1);
	t291 = -2 * pkin(5);
	t290 = pkin(4) * pkin(5);
	t289 = -pkin(8) - pkin(3);
	t288 = -pkin(8) + pkin(3);
	t287 = (-pkin(9) - pkin(11));
	t286 = (-pkin(9) + pkin(11));
	t229 = sin(pkin(21));
	t285 = t229 / 0.2e1;
	t234 = cos(qJ(3));
	t284 = -t234 / 0.2e1;
	t238 = pkin(7) ^ 2;
	t232 = sin(qJ(2));
	t236 = cos(pkin(19));
	t280 = sin(pkin(19));
	t281 = cos(qJ(2));
	t223 = t232 * t236 - t281 * t280;
	t277 = pkin(7) * t223;
	t261 = (pkin(1) ^ 2) + t277 * t292;
	t220 = t238 + t261;
	t217 = pkin(3) ^ 2 - pkin(8) ^ 2 + t220;
	t221 = pkin(1) - t277;
	t224 = t232 * t280 + t281 * t236;
	t215 = (pkin(7) - t289) * (pkin(7) + t289) + t261;
	t216 = (pkin(7) - t288) * (pkin(7) + t288) + t261;
	t244 = sqrt(-t216 * t215);
	t276 = pkin(7) * t224;
	t260 = pkin(1) * t276;
	t269 = 0.1e1 / t244 * (t215 + t216) * t260;
	t203 = (t223 * t244 + (t221 * t292 - t217 - t269) * t224) * pkin(7);
	t267 = t224 * t244;
	t204 = t221 * t269 + t238 * t224 ^ 2 * t292 + (-t223 * t217 - t267) * pkin(7);
	t218 = 0.1e1 / t220;
	t231 = sin(qJ(3));
	t241 = 0.1e1 / pkin(3);
	t253 = 0.1e1 / t220 ^ 2 * t260;
	t211 = t217 * t276 + t221 * t244;
	t263 = t234 * t211;
	t210 = -pkin(7) * t267 + t221 * t217;
	t266 = t231 * t210;
	t189 = ((-t231 * t203 / 0.2e1 + t204 * t284) * t218 + (-t263 - t266) * t253) * t241;
	t264 = t234 * t210;
	t265 = t231 * t211;
	t190 = ((t203 * t284 + t231 * t204 / 0.2e1) * t218 + (-t264 + t265) * t253) * t241;
	t228 = pkin(23) + pkin(22);
	t226 = sin(t228);
	t227 = cos(t228);
	t185 = t227 * t189 + t226 * t190;
	t268 = t218 * t241;
	t206 = (t266 / 0.2e1 + t263 / 0.2e1) * t268;
	t207 = (-t264 / 0.2e1 + t265 / 0.2e1) * t268;
	t199 = t227 * t206 - t226 * t207;
	t279 = pkin(4) * t199;
	t262 = (pkin(5) ^ 2) - t279 * t291;
	t191 = ((pkin(4) - t287) * (pkin(4) + t287)) + t262;
	t192 = ((pkin(4) - t286) * (pkin(4) + t286)) + t262;
	t249 = (t191 + t192) * t290;
	t181 = 0.2e1 * t185 * t249;
	t186 = -t226 * t189 + t227 * t190;
	t243 = sqrt(-t192 * t191);
	t240 = pkin(4) ^ 2;
	t196 = t240 + t262;
	t193 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t196;
	t197 = pkin(5) + t279;
	t254 = t197 * t291 - t193;
	t188 = 0.1e1 / t243;
	t250 = t226 * t206 + t227 * t207;
	t272 = t188 * t250;
	t172 = (-t186 * t243 - t181 * t272 / 0.2e1 + t254 * t185) * pkin(4);
	t255 = t240 * t250 * t291;
	t273 = t188 * t197;
	t173 = t181 * t273 / 0.2e1 + t185 * t255 + (-t185 * t243 + t186 * t193) * pkin(4);
	t278 = pkin(4) * t250;
	t182 = t193 * t197 - t243 * t278;
	t183 = t193 * t278 + t197 * t243;
	t230 = cos(pkin(21));
	t194 = 0.1e1 / t196;
	t237 = 0.1e1 / pkin(11);
	t270 = t194 * t237;
	t180 = (-t182 * t230 / 0.2e1 + t183 * t285) * t270;
	t177 = 0.1e1 / t180;
	t259 = 0.1e1 / t196 ^ 2 * t290;
	t251 = t230 * t259;
	t246 = t185 * t251;
	t252 = t229 * t259;
	t248 = t185 * t252;
	t271 = t194 * t230;
	t256 = t271 / 0.2e1;
	t257 = -t271 / 0.2e1;
	t258 = t194 * t285;
	t178 = 0.1e1 / t180 ^ 2;
	t179 = (t182 * t285 + t183 * t230 / 0.2e1) * t270;
	t274 = t178 * t179;
	t275 = 0.1e1 / (t178 * t179 ^ 2 + 0.1e1) * t237;
	t283 = ((t172 * t258 + t173 * t256 + t182 * t248 + t183 * t246) * t177 - (t172 * t257 + t173 * t258 - t182 * t246 + t183 * t248) * t274) * t275 + 0.1e1;
	t184 = t250 * t249;
	t175 = (-t184 * t272 - t199 * t243 + t250 * t254) * pkin(4);
	t176 = t184 * t273 + t250 * t255 + (t199 * t193 - t243 * t250) * pkin(4);
	t245 = t250 * t251;
	t247 = t250 * t252;
	t282 = ((t175 * t258 + t176 * t256 + t182 * t247 + t183 * t245) * t177 - (t175 * t257 + t176 * t258 - t182 * t245 + t183 * t247) * t274) * t275 + 0.1e1;
	t235 = cos(qJ(1));
	t233 = sin(qJ(1));
	t1 = [0, t283 * t233, t282 * t233, 0; 0, -t283 * t235, -t282 * t235, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:47:43
	% EndTime: 2020-04-14 18:47:52
	% DurationCPUTime: 6.18s
	% Computational Cost: add. (109481->91), mult. (165302->171), div. (7288->11), fcn. (104804->19), ass. (0->110)
	t532 = -2 * pkin(1);
	t531 = -2 * pkin(5);
	t530 = pkin(4) * pkin(5);
	t529 = -pkin(8) - pkin(3);
	t528 = -pkin(8) + pkin(3);
	t527 = (-pkin(9) - pkin(11));
	t526 = (-pkin(9) + pkin(11));
	t466 = sin(pkin(21));
	t525 = t466 / 0.2e1;
	t471 = cos(qJ(3));
	t524 = -t471 / 0.2e1;
	t476 = pkin(7) ^ 2;
	t469 = sin(qJ(2));
	t472 = cos(qJ(2));
	t474 = cos(pkin(19));
	t521 = sin(pkin(19));
	t460 = t469 * t474 - t472 * t521;
	t518 = pkin(7) * t460;
	t502 = (pkin(1) ^ 2) + t518 * t532;
	t457 = t476 + t502;
	t454 = pkin(3) ^ 2 - pkin(8) ^ 2 + t457;
	t458 = pkin(1) - t518;
	t461 = t469 * t521 + t472 * t474;
	t452 = (pkin(7) - t529) * (pkin(7) + t529) + t502;
	t453 = (pkin(7) - t528) * (pkin(7) + t528) + t502;
	t482 = sqrt(-t453 * t452);
	t517 = pkin(7) * t461;
	t501 = pkin(1) * t517;
	t510 = 0.1e1 / t482 * (t452 + t453) * t501;
	t440 = (t460 * t482 + (t458 * t532 - t454 - t510) * t461) * pkin(7);
	t508 = t461 * t482;
	t441 = t458 * t510 + t476 * t461 ^ 2 * t532 + (-t460 * t454 - t508) * pkin(7);
	t455 = 0.1e1 / t457;
	t468 = sin(qJ(3));
	t479 = 0.1e1 / pkin(3);
	t494 = 0.1e1 / t457 ^ 2 * t501;
	t448 = t454 * t517 + t458 * t482;
	t504 = t471 * t448;
	t447 = -pkin(7) * t508 + t458 * t454;
	t507 = t468 * t447;
	t426 = ((-t468 * t440 / 0.2e1 + t441 * t524) * t455 + (-t504 - t507) * t494) * t479;
	t505 = t471 * t447;
	t506 = t468 * t448;
	t427 = ((t440 * t524 + t468 * t441 / 0.2e1) * t455 + (-t505 + t506) * t494) * t479;
	t465 = pkin(23) + pkin(22);
	t463 = sin(t465);
	t464 = cos(t465);
	t422 = t464 * t426 + t463 * t427;
	t509 = t455 * t479;
	t443 = (t507 / 0.2e1 + t504 / 0.2e1) * t509;
	t444 = (-t505 / 0.2e1 + t506 / 0.2e1) * t509;
	t436 = t464 * t443 - t463 * t444;
	t520 = pkin(4) * t436;
	t503 = (pkin(5) ^ 2) - t520 * t531;
	t428 = ((pkin(4) - t527) * (pkin(4) + t527)) + t503;
	t429 = ((pkin(4) - t526) * (pkin(4) + t526)) + t503;
	t490 = (t428 + t429) * t530;
	t418 = 0.2e1 * t422 * t490;
	t423 = -t463 * t426 + t464 * t427;
	t481 = sqrt(-t429 * t428);
	t478 = pkin(4) ^ 2;
	t433 = t478 + t503;
	t430 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t433;
	t434 = pkin(5) + t520;
	t495 = t434 * t531 - t430;
	t425 = 0.1e1 / t481;
	t491 = t463 * t443 + t464 * t444;
	t513 = t425 * t491;
	t406 = (-t423 * t481 - t418 * t513 / 0.2e1 + t495 * t422) * pkin(4);
	t496 = t478 * t491 * t531;
	t514 = t425 * t434;
	t407 = t418 * t514 / 0.2e1 + t422 * t496 + (-t422 * t481 + t423 * t430) * pkin(4);
	t519 = pkin(4) * t491;
	t419 = t434 * t430 - t481 * t519;
	t420 = t430 * t519 + t434 * t481;
	t467 = cos(pkin(21));
	t431 = 0.1e1 / t433;
	t475 = 0.1e1 / pkin(11);
	t511 = t431 * t475;
	t417 = (-t419 * t467 / 0.2e1 + t420 * t525) * t511;
	t414 = 0.1e1 / t417;
	t500 = 0.1e1 / t433 ^ 2 * t530;
	t492 = t467 * t500;
	t487 = t422 * t492;
	t493 = t466 * t500;
	t489 = t422 * t493;
	t512 = t431 * t467;
	t497 = t512 / 0.2e1;
	t498 = -t512 / 0.2e1;
	t499 = t431 * t525;
	t415 = 0.1e1 / t417 ^ 2;
	t416 = (t419 * t525 + t420 * t467 / 0.2e1) * t511;
	t515 = t415 * t416;
	t516 = 0.1e1 / (t415 * t416 ^ 2 + 0.1e1) * t475;
	t523 = ((t406 * t499 + t407 * t497 + t419 * t489 + t420 * t487) * t414 - (t406 * t498 + t407 * t499 - t419 * t487 + t420 * t489) * t515) * t516 + 0.1e1;
	t421 = t491 * t490;
	t412 = (-t421 * t513 - t436 * t481 + t491 * t495) * pkin(4);
	t413 = t421 * t514 + t491 * t496 + (t436 * t430 - t481 * t491) * pkin(4);
	t486 = t491 * t492;
	t488 = t491 * t493;
	t522 = ((t412 * t499 + t413 * t497 + t419 * t488 + t420 * t486) * t414 - (t412 * t498 + t413 * t499 - t419 * t486 + t420 * t488) * t515) * t516 + 0.1e1;
	t485 = t472 * t468 + t469 * t471;
	t484 = -t469 * t468 + t472 * t471;
	t411 = atan2(t416, t417);
	t408 = sin(t411);
	t409 = cos(t411);
	t483 = t484 * t408 + t485 * t409;
	t473 = cos(qJ(1));
	t470 = sin(qJ(1));
	t1 = [0, t523 * t470, t522 * t470, t483 * t473; 0, -t523 * t473, -t522 * t473, t483 * t470; 1, 0, 0, t485 * t408 - t484 * t409;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:59
	% EndTime: 2020-04-14 18:43:00
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
	% StartTime: 2020-04-14 18:43:01
	% EndTime: 2020-04-14 18:43:01
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (1522->37), mult. (2350->66), div. (100->7), fcn. (1512->10), ass. (0->42)
	t123 = -2 * pkin(1);
	t95 = sin(pkin(23));
	t122 = t95 / 0.2e1;
	t121 = -pkin(8) - pkin(3);
	t120 = -pkin(8) + pkin(3);
	t102 = 0.1e1 / pkin(3);
	t100 = cos(pkin(19));
	t115 = sin(pkin(19));
	t116 = cos(qJ(2));
	t97 = sin(qJ(2));
	t93 = t116 * t100 + t97 * t115;
	t117 = pkin(7) * t93;
	t109 = pkin(1) * t117;
	t101 = pkin(7) ^ 2;
	t92 = t97 * t100 - t116 * t115;
	t118 = pkin(7) * t92;
	t110 = (pkin(1) ^ 2) + t118 * t123;
	t89 = t101 + t110;
	t107 = 0.1e1 / t89 ^ 2 * t109;
	t96 = cos(pkin(23));
	t105 = t96 * t107;
	t106 = t95 * t107;
	t87 = 0.1e1 / t89;
	t108 = t87 * t122;
	t113 = t96 * t87;
	t84 = (pkin(7) - t121) * (pkin(7) + t121) + t110;
	t85 = (pkin(7) - t120) * (pkin(7) + t120) + t110;
	t104 = sqrt(-t85 * t84);
	t114 = 0.1e1 / t104 * (t84 + t85) * t109;
	t86 = pkin(3) ^ 2 - pkin(8) ^ 2 + t89;
	t90 = pkin(1) - t118;
	t74 = (t92 * t104 + (t90 * t123 - t114 - t86) * t93) * pkin(7);
	t111 = t104 * t93;
	t75 = t90 * t114 + t101 * t93 ^ 2 * t123 + (-t92 * t86 - t111) * pkin(7);
	t112 = t102 * t87;
	t79 = -pkin(7) * t111 + t90 * t86;
	t80 = t104 * t90 + t86 * t117;
	t77 = (-t96 * t79 / 0.2e1 + t80 * t122) * t112;
	t76 = 0.1e1 / t77 ^ 2;
	t78 = (t96 * t80 / 0.2e1 + t79 * t122) * t112;
	t119 = ((t75 * t113 / 0.2e1 + t80 * t105 + t74 * t108 + t79 * t106) / t77 - (-t74 * t113 / 0.2e1 - t79 * t105 + t75 * t108 + t80 * t106) * t78 * t76) / (t78 ^ 2 * t76 + 0.1e1) * t102 + 0.1e1;
	t1 = [0, t119 * sin(qJ(1)), 0, 0; 0, -t119 * cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobig_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:59
	% EndTime: 2020-04-14 18:43:00
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (726->31), mult. (1088->50), div. (28->6), fcn. (696->8), ass. (0->32)
	t78 = pkin(6) ^ 2;
	t72 = sin(pkin(20));
	t73 = cos(pkin(20));
	t74 = sin(qJ(3));
	t76 = cos(qJ(3));
	t70 = t76 * t72 + t74 * t73;
	t87 = pkin(6) * t70;
	t83 = 0.2e1 * pkin(1) * t87 + t78;
	t67 = pkin(1) ^ 2 + t83;
	t65 = 0.1e1 / t67;
	t92 = t65 / 0.2e1;
	t91 = -pkin(2) - pkin(13);
	t90 = -pkin(2) + pkin(13);
	t71 = t74 * t72 - t76 * t73;
	t89 = pkin(1) * t71;
	t64 = pkin(2) ^ 2 - pkin(13) ^ 2 + t67;
	t68 = -pkin(1) - t87;
	t62 = (pkin(1) - t91) * (pkin(1) + t91) + t83;
	t63 = (pkin(1) - t90) * (pkin(1) + t90) + t83;
	t81 = sqrt(-t63 * t62);
	t86 = t71 * pkin(6);
	t58 = t64 * t86 - t68 * t81;
	t88 = pkin(6) * t58;
	t85 = 0.1e1 / t81 * (t62 + t63) * pkin(1) * t86;
	t84 = t71 * t81;
	t77 = cos(qJ(1));
	t75 = sin(qJ(1));
	t66 = 0.1e1 / t67 ^ 2;
	t57 = -t84 * pkin(6) - t68 * t64;
	t56 = 0.1e1 / t57 ^ 2;
	t54 = 0.2e1 * (((-t68 * t85 + (t70 * t64 - t84) * pkin(6)) * t92 + (-t65 * t78 * t71 + t66 * t88) * t89) / t57 - ((-t70 * t81 + (-t64 - t85) * t71) * t92 + (t57 * t66 + t65 * t68) * t89) * t56 * t88) / (t58 ^ 2 * t56 + 0.1e1) * t67;
	t1 = [0, t75, t54 * t75, 0; 0, -t77, -t54 * t77, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobig_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:43:02
	% EndTime: 2020-04-14 18:43:03
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (925->37), mult. (1360->56), div. (40->9), fcn. (842->8), ass. (0->41)
	t141 = -pkin(2) - pkin(13);
	t140 = -pkin(2) + pkin(13);
	t113 = sin(pkin(20));
	t114 = cos(pkin(20));
	t115 = sin(qJ(3));
	t117 = cos(qJ(3));
	t111 = t113 * t117 + t114 * t115;
	t136 = pkin(6) * t111;
	t128 = pkin(1) * t136;
	t110 = 0.2e1 * t128;
	t121 = pkin(2) ^ 2;
	t120 = pkin(6) ^ 2;
	t129 = pkin(1) ^ 2 + t120;
	t126 = -pkin(13) ^ 2 + t129;
	t105 = t110 + t121 + t126;
	t109 = -pkin(1) - t136;
	t130 = t110 + t120;
	t101 = (pkin(1) - t141) * (pkin(1) + t141) + t130;
	t102 = (pkin(1) - t140) * (pkin(1) + t140) + t130;
	t132 = t102 * t101;
	t124 = sqrt(-t132);
	t112 = t113 * t115 - t114 * t117;
	t135 = pkin(6) * t112;
	t96 = t105 * t135 - t109 * t124;
	t139 = pkin(6) * t96;
	t108 = t110 + t129;
	t106 = 0.1e1 / t108;
	t138 = t106 / 0.2e1;
	t137 = pkin(1) * t112;
	t104 = t121 - t126 - 0.2e1 * t128;
	t103 = 0.1e1 / t104 ^ 2;
	t107 = 0.1e1 / t108 ^ 2;
	t131 = t112 * t124;
	t127 = pkin(6) * t131;
	t133 = 0.1e1 / t124 * (t101 + t102) * pkin(1) * t135;
	t95 = -t105 * t109 - t127;
	t94 = 0.1e1 / t95 ^ 2;
	t134 = 0.2e1 * (((-t109 * t133 + (t111 * t105 - t131) * pkin(6)) * t138 + (-t106 * t112 * t120 + t107 * t139) * t137) / t95 - ((-t111 * t124 + (-t105 - t133) * t112) * t138 + (t106 * t109 + t107 * t95) * t137) * t94 * t139) * t108 / (t94 * t96 ^ 2 + 0.1e1) + (0.1e1 / t104 * t133 - 0.2e1 * pkin(1) * t103 * t127) / (-t103 * t132 + 0.1e1);
	t118 = cos(qJ(1));
	t116 = sin(qJ(1));
	t1 = [0, t116, t134 * t116, 0; 0, -t118, -t134 * t118, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobig_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:44:14
	% EndTime: 2020-04-14 18:44:18
	% DurationCPUTime: 2.84s
	% Computational Cost: add. (45553->86), mult. (68854->160), div. (2940->13), fcn. (43692->16), ass. (0->92)
	t237 = pkin(5) ^ 2;
	t236 = pkin(7) ^ 2;
	t230 = sin(qJ(2));
	t234 = cos(pkin(19));
	t273 = sin(pkin(19));
	t274 = cos(qJ(2));
	t221 = t230 * t234 - t274 * t273;
	t269 = pkin(7) * t221;
	t284 = -2 * pkin(1);
	t257 = (pkin(1) ^ 2) + t269 * t284;
	t218 = t236 + t257;
	t215 = pkin(3) ^ 2 - pkin(8) ^ 2 + t218;
	t219 = pkin(1) - t269;
	t283 = -pkin(8) - pkin(3);
	t213 = (pkin(7) - t283) * (pkin(7) + t283) + t257;
	t282 = -pkin(8) + pkin(3);
	t214 = (pkin(7) - t282) * (pkin(7) + t282) + t257;
	t242 = sqrt(-t214 * t213);
	t222 = t230 * t273 + t274 * t234;
	t268 = pkin(7) * t222;
	t209 = t215 * t268 + t219 * t242;
	t232 = cos(qJ(3));
	t259 = t232 * t209;
	t263 = t222 * t242;
	t208 = -pkin(7) * t263 + t215 * t219;
	t229 = sin(qJ(3));
	t262 = t229 * t208;
	t216 = 0.1e1 / t218;
	t239 = 0.1e1 / pkin(3);
	t264 = t216 * t239;
	t204 = (t262 / 0.2e1 + t259 / 0.2e1) * t264;
	t260 = t232 * t208;
	t261 = t229 * t209;
	t205 = (-t260 / 0.2e1 + t261 / 0.2e1) * t264;
	t226 = pkin(23) + pkin(22);
	t224 = sin(t226);
	t225 = cos(t226);
	t193 = t204 * t225 - t205 * t224;
	t271 = t193 * pkin(5);
	t258 = 0.2e1 * pkin(4) * t271 + t237;
	t281 = -pkin(9) - pkin(11);
	t185 = (pkin(4) - t281) * (pkin(4) + t281) + t258;
	t280 = -pkin(9) + pkin(11);
	t186 = (pkin(4) - t280) * (pkin(4) + t280) + t258;
	t241 = sqrt(-t186 * t185);
	t285 = 0.1e1 / t241 * pkin(4) * pkin(5) * (t185 + t186);
	t190 = pkin(4) ^ 2 + t258;
	t188 = 0.1e1 / t190;
	t279 = t188 / 0.2e1;
	t256 = pkin(1) * t268;
	t266 = 0.1e1 / t242 * (t213 + t214) * t256;
	t198 = (t221 * t242 + (t219 * t284 - t215 - t266) * t222) * pkin(7);
	t278 = -t198 / 0.2e1;
	t199 = t219 * t266 + t236 * t222 ^ 2 * t284 + (-t221 * t215 - t263) * pkin(7);
	t277 = t199 / 0.2e1;
	t227 = sin(pkin(23));
	t276 = t227 / 0.2e1;
	t275 = -t232 / 0.2e1;
	t187 = pkin(9) ^ 2 - pkin(11) ^ 2 + t190;
	t191 = -pkin(4) - t271;
	t251 = t204 * t224 + t225 * t205;
	t270 = pkin(5) * t251;
	t177 = t187 * t270 - t191 * t241;
	t272 = pkin(5) * t177;
	t267 = t251 * t285;
	t228 = cos(pkin(23));
	t265 = t216 * t228;
	t252 = 0.1e1 / t218 ^ 2 * t256;
	t183 = ((t199 * t275 + t229 * t278) * t216 + (-t259 - t262) * t252) * t239;
	t184 = ((t198 * t275 + t229 * t277) * t216 + (-t260 + t261) * t252) * t239;
	t179 = t183 * t225 + t184 * t224;
	t180 = -t183 * t224 + t184 * t225;
	t202 = (-t228 * t208 / 0.2e1 + t209 * t276) * t264;
	t201 = 0.1e1 / t202 ^ 2;
	t203 = (t228 * t209 / 0.2e1 + t208 * t276) * t264;
	t189 = 0.1e1 / t190 ^ 2;
	t243 = pkin(4) * (-t188 * t237 * t251 + t189 * t272);
	t176 = -t187 * t191 - t241 * t270;
	t175 = 0.1e1 / t176 ^ 2;
	t250 = 0.1e1 / (t175 * t177 ^ 2 + 0.1e1) * t190;
	t244 = t175 * t250 * t272;
	t245 = pkin(4) * (t176 * t189 + t188 * t191);
	t246 = 0.1e1 / t176 * t250;
	t247 = t209 * t252;
	t248 = t208 * t252;
	t253 = t216 * t276;
	t254 = t179 * t285;
	t255 = 0.2e1 * ((-t191 * t254 + (-t179 * t241 + t180 * t187) * pkin(5)) * t279 + t179 * t243) * t246 - 0.2e1 * ((-t179 * t187 - t180 * t241 - t251 * t254) * t279 + t179 * t245) * t244 + ((t198 * t253 + t227 * t248 + t228 * t247 + t265 * t277) / t202 - (t199 * t253 + t227 * t247 - t228 * t248 + t265 * t278) * t203 * t201) / (t201 * t203 ^ 2 + 0.1e1) * t239 + 0.1e1;
	t233 = cos(qJ(1));
	t231 = sin(qJ(1));
	t170 = 0.2e1 * ((-t191 * t267 + (t193 * t187 - t241 * t251) * pkin(5)) * t279 + t251 * t243) * t246 - 0.2e1 * (-t193 * t241 * t279 + ((-t187 - t267) * t279 + t245) * t251) * t244;
	t1 = [0, t255 * t231, t170 * t231, 0; 0, -t255 * t233, -t170 * t233, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
end