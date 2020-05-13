% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh3m1DE1
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
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = palh3m1DE1_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_jacobig_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1DE1_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_jacobig_rot_sym_varpar: pkin has to be [19x1] (double)');
Jg_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)), 0, 0; 0, -cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
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
	% StartTime: 2020-04-19 18:19:02
	% EndTime: 2020-04-19 18:19:11
	% DurationCPUTime: 5.29s
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
	t221 = t229 * t278 + t231 * t277;
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
	t207 = -pkin(5) * t260 + t218 * t214;
	t266 = t207 * t228;
	t187 = ((t202 * t232 / 0.2e1 + t201 * t281) * t215 + (t263 + t266) * t250) * t238;
	t264 = t208 * t228;
	t265 = t207 * t232;
	t188 = ((-t201 * t232 / 0.2e1 + t202 * t281) * t215 + (t264 - t265) * t250) * t238;
	t225 = pkin(18) + pkin(19);
	t223 = sin(t225);
	t224 = cos(t225);
	t184 = -t223 * t187 - t224 * t188;
	t261 = t215 * t238;
	t203 = (-t265 / 0.2e1 + t264 / 0.2e1) * t261;
	t204 = (t263 / 0.2e1 + t266 / 0.2e1) * t261;
	t198 = t224 * t203 + t223 * t204;
	t275 = pkin(3) * t198;
	t259 = (pkin(4) ^ 2) - t275 * t288;
	t189 = ((pkin(3) - t284) * (pkin(3) + t284)) + t259;
	t190 = ((pkin(3) - t283) * (pkin(3) + t283)) + t259;
	t247 = (t189 + t190) * t287;
	t179 = t184 * t247;
	t183 = -t224 * t187 + t223 * t188;
	t240 = sqrt(-t190 * t189);
	t237 = pkin(3) ^ 2;
	t194 = t237 + t259;
	t191 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t194;
	t195 = pkin(4) + t275;
	t251 = t195 * t288 - t191;
	t186 = 0.1e1 / t240;
	t242 = t223 * t203 - t224 * t204;
	t269 = t186 * t242;
	t170 = (-t179 * t269 - t183 * t240 + t184 * t251) * pkin(3);
	t252 = t237 * t242 * t288;
	t270 = t186 * t195;
	t171 = t179 * t270 + t184 * t252 + (t183 * t191 - t184 * t240) * pkin(3);
	t276 = pkin(3) * t242;
	t180 = t195 * t191 - t240 * t276;
	t181 = t191 * t276 + t195 * t240;
	t227 = cos(pkin(17));
	t192 = 0.1e1 / t194;
	t234 = 0.1e1 / pkin(10);
	t267 = t192 * t234;
	t177 = (-t180 * t227 / 0.2e1 + t181 * t282) * t267;
	t175 = 0.1e1 / t177;
	t256 = 0.1e1 / t194 ^ 2 * t287;
	t248 = t227 * t256;
	t244 = t184 * t248;
	t249 = t226 * t256;
	t246 = t184 * t249;
	t268 = t192 * t227;
	t253 = t268 / 0.2e1;
	t254 = -t268 / 0.2e1;
	t255 = t192 * t282;
	t176 = 0.1e1 / t177 ^ 2;
	t178 = (t181 * t227 / 0.2e1 + t180 * t282) * t267;
	t271 = t176 * t178;
	t272 = 0.1e1 / (t176 * t178 ^ 2 + 0.1e1) * t234;
	t280 = ((t170 * t255 + t171 * t253 + t180 * t246 + t181 * t244) * t175 - (t170 * t254 + t171 * t255 - t180 * t244 + t181 * t246) * t271) * t272 + 0.1e1;
	t182 = t242 * t247;
	t173 = (-t182 * t269 - t198 * t240 + t242 * t251) * pkin(3);
	t174 = t182 * t270 + t242 * t252 + (t198 * t191 - t240 * t242) * pkin(3);
	t243 = t242 * t248;
	t245 = t242 * t249;
	t279 = ((t173 * t255 + t174 * t253 + t180 * t245 + t181 * t243) * t175 - (t173 * t254 + t174 * t255 - t180 * t243 + t181 * t245) * t271) * t272 + 0.1e1;
	t233 = cos(qJ(1));
	t230 = sin(qJ(1));
	t1 = [0, t280 * t230, t279 * t230, 0; 0, -t280 * t233, -t279 * t233, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:22:35
	% EndTime: 2020-04-19 18:22:46
	% DurationCPUTime: 6.32s
	% Computational Cost: add. (109481->91), mult. (165302->169), div. (7288->11), fcn. (104804->19), ass. (0->110)
	t534 = -2 * pkin(1);
	t533 = -2 * pkin(4);
	t532 = pkin(3) * pkin(4);
	t531 = -pkin(6) - pkin(2);
	t530 = -pkin(6) + pkin(2);
	t529 = (-pkin(8) - pkin(10));
	t528 = (-pkin(8) + pkin(10));
	t468 = sin(pkin(17));
	t527 = t468 / 0.2e1;
	t470 = sin(qJ(3));
	t526 = t470 / 0.2e1;
	t478 = pkin(5) ^ 2;
	t471 = sin(qJ(2));
	t473 = sin(pkin(16));
	t475 = cos(qJ(2));
	t523 = cos(pkin(16));
	t462 = t471 * t473 - t475 * t523;
	t520 = pkin(5) * t462;
	t504 = (pkin(1) ^ 2) + t520 * t534;
	t459 = t478 + t504;
	t456 = pkin(2) ^ 2 - pkin(6) ^ 2 + t459;
	t460 = pkin(1) - t520;
	t463 = t471 * t523 + t475 * t473;
	t454 = (pkin(5) - t531) * (pkin(5) + t531) + t504;
	t455 = (pkin(5) - t530) * (pkin(5) + t530) + t504;
	t484 = sqrt(-t455 * t454);
	t519 = pkin(5) * t463;
	t503 = pkin(1) * t519;
	t508 = 0.1e1 / t484 * (t454 + t455) * t503;
	t443 = (t462 * t484 + (t460 * t534 - t456 - t508) * t463) * pkin(5);
	t506 = t463 * t484;
	t444 = t460 * t508 + t478 * t463 ^ 2 * t534 + (-t462 * t456 - t506) * pkin(5);
	t457 = 0.1e1 / t459;
	t474 = cos(qJ(3));
	t481 = 0.1e1 / pkin(2);
	t496 = 0.1e1 / t459 ^ 2 * t503;
	t450 = t456 * t519 + t460 * t484;
	t509 = t450 * t474;
	t449 = -pkin(5) * t506 + t460 * t456;
	t512 = t449 * t470;
	t429 = ((t444 * t474 / 0.2e1 + t443 * t526) * t457 + (t509 + t512) * t496) * t481;
	t510 = t450 * t470;
	t511 = t449 * t474;
	t430 = ((-t443 * t474 / 0.2e1 + t444 * t526) * t457 + (t510 - t511) * t496) * t481;
	t467 = pkin(18) + pkin(19);
	t465 = sin(t467);
	t466 = cos(t467);
	t426 = -t465 * t429 - t466 * t430;
	t507 = t457 * t481;
	t445 = (-t511 / 0.2e1 + t510 / 0.2e1) * t507;
	t446 = (t509 / 0.2e1 + t512 / 0.2e1) * t507;
	t440 = t466 * t445 + t465 * t446;
	t521 = pkin(3) * t440;
	t505 = (pkin(4) ^ 2) - t521 * t533;
	t431 = ((pkin(3) - t529) * (pkin(3) + t529)) + t505;
	t432 = ((pkin(3) - t528) * (pkin(3) + t528)) + t505;
	t493 = (t431 + t432) * t532;
	t421 = t426 * t493;
	t425 = -t466 * t429 + t465 * t430;
	t483 = sqrt(-t432 * t431);
	t480 = pkin(3) ^ 2;
	t436 = t480 + t505;
	t433 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t436;
	t437 = pkin(4) + t521;
	t497 = t437 * t533 - t433;
	t428 = 0.1e1 / t483;
	t488 = t465 * t445 - t466 * t446;
	t515 = t428 * t488;
	t409 = (-t421 * t515 - t425 * t483 + t497 * t426) * pkin(3);
	t498 = t488 * t480 * t533;
	t516 = t428 * t437;
	t410 = t421 * t516 + t426 * t498 + (t425 * t433 - t426 * t483) * pkin(3);
	t522 = pkin(3) * t488;
	t422 = t437 * t433 - t483 * t522;
	t423 = t433 * t522 + t437 * t483;
	t469 = cos(pkin(17));
	t434 = 0.1e1 / t436;
	t477 = 0.1e1 / pkin(10);
	t513 = t434 * t477;
	t419 = (-t422 * t469 / 0.2e1 + t423 * t527) * t513;
	t417 = 0.1e1 / t419;
	t502 = 0.1e1 / t436 ^ 2 * t532;
	t495 = t426 * t502;
	t490 = t423 * t495;
	t492 = t422 * t495;
	t514 = t434 * t469;
	t499 = t514 / 0.2e1;
	t500 = -t514 / 0.2e1;
	t501 = t434 * t527;
	t418 = 0.1e1 / t419 ^ 2;
	t420 = (t423 * t469 / 0.2e1 + t422 * t527) * t513;
	t517 = t418 * t420;
	t518 = 0.1e1 / (t420 ^ 2 * t418 + 0.1e1) * t477;
	t525 = ((t409 * t501 + t410 * t499 + t468 * t492 + t469 * t490) * t417 - (t409 * t500 + t410 * t501 + t468 * t490 - t469 * t492) * t517) * t518 + 0.1e1;
	t424 = t488 * t493;
	t415 = (-t424 * t515 - t440 * t483 + t488 * t497) * pkin(3);
	t416 = t424 * t516 + t488 * t498 + (t440 * t433 - t483 * t488) * pkin(3);
	t494 = t488 * t502;
	t489 = t423 * t494;
	t491 = t422 * t494;
	t524 = ((t415 * t501 + t416 * t499 + t468 * t491 + t469 * t489) * t417 - (t415 * t500 + t416 * t501 + t468 * t489 - t469 * t491) * t517) * t518 + 0.1e1;
	t487 = t475 * t470 + t471 * t474;
	t486 = t471 * t470 - t475 * t474;
	t414 = atan2(t420, t419);
	t411 = sin(t414);
	t412 = cos(t414);
	t485 = t486 * t411 - t487 * t412;
	t476 = cos(qJ(1));
	t472 = sin(qJ(1));
	t1 = [0, t525 * t472, t524 * t472, t485 * t476; 0, -t525 * t476, -t524 * t476, t485 * t472; 1, 0, 0, -t487 * t411 - t486 * t412;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:15
	% EndTime: 2020-04-19 18:18:16
	% DurationCPUTime: 0.23s
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
	t72 = -pkin(2) ^ 2 + pkin(6) ^ 2 + t75;
	t76 = -pkin(5) + t102;
	t65 = -pkin(1) * t95 - t72 * t76;
	t90 = t65 * t91;
	t66 = t72 * t101 - t76 * t88;
	t89 = t66 * t91;
	t64 = (t66 * t105 - t65 * t83 / 0.2e1) * t96;
	t63 = (t65 * t105 + t66 * t83 / 0.2e1) * t96;
	t62 = 0.1e1 / t63 ^ 2;
	t61 = -t76 * t98 + t87 * t79 ^ 2 * t106 + (-t72 * t78 - t95) * pkin(1);
	t60 = (t78 * t88 + (0.2e1 * t76 * pkin(5) - t72 - t98) * t79) * pkin(1);
	t58 = ((t61 * t92 + t84 * t89 - t60 * t97 / 0.2e1 - t83 * t90) / t63 - (t60 * t92 + t84 * t90 + t61 * t97 / 0.2e1 + t83 * t89) * t64 * t62) / (t62 * t64 ^ 2 + 0.1e1) * t85;
	t1 = [0, t58 * sin(qJ(1)), 0, 0; 0, -t58 * cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobig_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:16
	% EndTime: 2020-04-19 18:18:16
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (1522->37), mult. (2350->66), div. (100->7), fcn. (1512->10), ass. (0->42)
	t118 = -2 * pkin(1);
	t90 = sin(pkin(19));
	t117 = t90 / 0.2e1;
	t116 = -pkin(6) - pkin(2);
	t115 = -pkin(6) + pkin(2);
	t110 = cos(qJ(2));
	t111 = cos(pkin(16));
	t92 = sin(qJ(2));
	t94 = sin(pkin(16));
	t88 = t110 * t94 + t92 * t111;
	t112 = pkin(5) * t88;
	t104 = pkin(1) * t112;
	t87 = -t110 * t111 + t92 * t94;
	t113 = pkin(5) * t87;
	t105 = (pkin(1) ^ 2) + t113 * t118;
	t96 = pkin(5) ^ 2;
	t84 = t96 + t105;
	t102 = 0.1e1 / t84 ^ 2 * t104;
	t91 = cos(pkin(19));
	t100 = t91 * t102;
	t101 = t90 * t102;
	t82 = 0.1e1 / t84;
	t103 = t82 * t117;
	t108 = t82 * t91;
	t79 = (pkin(5) - t116) * (pkin(5) + t116) + t105;
	t80 = (pkin(5) - t115) * (pkin(5) + t115) + t105;
	t99 = sqrt(-t80 * t79);
	t109 = 0.1e1 / t99 * (t79 + t80) * t104;
	t81 = pkin(2) ^ 2 - pkin(6) ^ 2 + t84;
	t85 = pkin(1) - t113;
	t69 = (t87 * t99 + (t85 * t118 - t109 - t81) * t88) * pkin(5);
	t106 = t88 * t99;
	t70 = t85 * t109 + t96 * t88 ^ 2 * t118 + (-t87 * t81 - t106) * pkin(5);
	t97 = 0.1e1 / pkin(2);
	t107 = t82 * t97;
	t74 = -pkin(5) * t106 + t85 * t81;
	t75 = t81 * t112 + t85 * t99;
	t72 = (-t74 * t91 / 0.2e1 + t75 * t117) * t107;
	t71 = 0.1e1 / t72 ^ 2;
	t73 = (t75 * t91 / 0.2e1 + t74 * t117) * t107;
	t114 = ((t70 * t108 / 0.2e1 + t75 * t100 + t69 * t103 + t74 * t101) / t72 - (-t69 * t108 / 0.2e1 - t74 * t100 + t70 * t103 + t75 * t101) * t73 * t71) / (t73 ^ 2 * t71 + 0.1e1) * t97 + 0.1e1;
	t1 = [0, t114 * sin(qJ(1)), 0, 0; 0, -t114 * cos(qJ(1)), 0, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobig_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:19:11
	% EndTime: 2020-04-19 18:19:16
	% DurationCPUTime: 2.82s
	% Computational Cost: add. (45553->86), mult. (68854->163), div. (2940->13), fcn. (43692->16), ass. (0->95)
	t279 = -2 * pkin(1);
	t278 = -pkin(6) - pkin(2);
	t277 = -pkin(6) + pkin(2);
	t276 = -pkin(8) - pkin(10);
	t275 = -pkin(8) + pkin(10);
	t232 = pkin(4) ^ 2;
	t231 = pkin(5) ^ 2;
	t225 = sin(qJ(2));
	t227 = sin(pkin(16));
	t268 = cos(qJ(2));
	t269 = cos(pkin(16));
	t216 = t225 * t227 - t268 * t269;
	t264 = pkin(5) * t216;
	t251 = (pkin(1) ^ 2) + t264 * t279;
	t213 = t231 + t251;
	t211 = 0.1e1 / t213;
	t234 = 0.1e1 / pkin(2);
	t254 = t211 * t234;
	t210 = pkin(2) ^ 2 - pkin(6) ^ 2 + t213;
	t214 = pkin(1) - t264;
	t208 = (pkin(5) - t278) * (pkin(5) + t278) + t251;
	t209 = (pkin(5) - t277) * (pkin(5) + t277) + t251;
	t237 = sqrt(-t209 * t208);
	t217 = t225 * t269 + t268 * t227;
	t263 = pkin(5) * t217;
	t204 = t210 * t263 + t214 * t237;
	t224 = sin(qJ(3));
	t258 = t204 * t224;
	t253 = t217 * t237;
	t203 = -pkin(5) * t253 + t210 * t214;
	t228 = cos(qJ(3));
	t259 = t203 * t228;
	t199 = (-t259 / 0.2e1 + t258 / 0.2e1) * t254;
	t257 = t204 * t228;
	t260 = t203 * t224;
	t200 = (t257 / 0.2e1 + t260 / 0.2e1) * t254;
	t221 = pkin(18) + pkin(19);
	t219 = sin(t221);
	t220 = cos(t221);
	t190 = t199 * t220 + t200 * t219;
	t265 = t190 * pkin(4);
	t252 = 0.2e1 * pkin(3) * t265 + t232;
	t186 = pkin(3) ^ 2 + t252;
	t184 = 0.1e1 / t186;
	t274 = t184 / 0.2e1;
	t250 = pkin(1) * t263;
	t256 = 0.1e1 / t237 * (t208 + t209) * t250;
	t194 = (t216 * t237 + (t214 * t279 - t210 - t256) * t217) * pkin(5);
	t273 = -t194 / 0.2e1;
	t195 = t214 * t256 + t231 * t217 ^ 2 * t279 + (-t216 * t210 - t253) * pkin(5);
	t272 = t195 / 0.2e1;
	t222 = sin(pkin(19));
	t271 = t222 / 0.2e1;
	t270 = t224 / 0.2e1;
	t183 = pkin(8) ^ 2 - pkin(10) ^ 2 + t186;
	t187 = -pkin(3) - t265;
	t181 = (pkin(3) - t276) * (pkin(3) + t276) + t252;
	t182 = (pkin(3) - t275) * (pkin(3) + t275) + t252;
	t236 = sqrt(-t182 * t181);
	t241 = t199 * t219 - t200 * t220;
	t266 = pkin(4) * t241;
	t173 = t183 * t266 - t187 * t236;
	t267 = pkin(4) * t173;
	t178 = 0.1e1 / t236;
	t262 = t178 * t187;
	t261 = t178 * t241;
	t223 = cos(pkin(19));
	t255 = t211 * t223;
	t247 = 0.1e1 / t213 ^ 2 * t250;
	t179 = ((t194 * t270 + t228 * t272) * t211 + (t257 + t260) * t247) * t234;
	t180 = ((t195 * t270 + t228 * t273) * t211 + (t258 - t259) * t247) * t234;
	t176 = -t179 * t219 - t180 * t220;
	t245 = pkin(3) * pkin(4) * (t181 + t182);
	t168 = t176 * t245;
	t175 = -t179 * t220 + t180 * t219;
	t197 = (-t203 * t223 / 0.2e1 + t204 * t271) * t254;
	t196 = 0.1e1 / t197 ^ 2;
	t198 = (t204 * t223 / 0.2e1 + t203 * t271) * t254;
	t185 = 0.1e1 / t186 ^ 2;
	t238 = pkin(3) * (-t184 * t232 * t241 + t185 * t267);
	t172 = -t183 * t187 - t236 * t266;
	t171 = 0.1e1 / t172 ^ 2;
	t246 = 0.1e1 / (t171 * t173 ^ 2 + 0.1e1) * t186;
	t239 = t171 * t246 * t267;
	t240 = pkin(3) * (t172 * t185 + t184 * t187);
	t242 = 0.1e1 / t172 * t246;
	t243 = t223 * t247;
	t244 = t222 * t247;
	t248 = t211 * t271;
	t249 = 0.2e1 * ((-t168 * t262 + (t175 * t183 - t176 * t236) * pkin(4)) * t274 + t176 * t238) * t242 - 0.2e1 * ((-t168 * t261 - t175 * t236 - t176 * t183) * t274 + t176 * t240) * t239 + ((t194 * t248 + t203 * t244 + t204 * t243 + t255 * t272) / t197 - (t195 * t248 - t203 * t243 + t204 * t244 + t255 * t273) * t198 * t196) / (t196 * t198 ^ 2 + 0.1e1) * t234 + 0.1e1;
	t229 = cos(qJ(1));
	t226 = sin(qJ(1));
	t174 = t241 * t245;
	t166 = 0.2e1 * ((-t174 * t262 + (t190 * t183 - t236 * t241) * pkin(4)) * t274 + t241 * t238) * t242 - 0.2e1 * ((-t174 * t261 - t183 * t241 - t190 * t236) * t274 + t241 * t240) * t239;
	t1 = [0, t249 * t226, t166 * t226, 0; 0, -t249 * t229, -t166 * t229, 0; 1, 0, 0, 0;];
	Jg_rot = t1;
end