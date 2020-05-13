% Calculate vector of centrifugal and Coriolis load on the joints for
% palh1m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh1m2TE_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_coriolisvecJ_fixb_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2TE_coriolisvecJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2TE_coriolisvecJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2TE_coriolisvecJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:43:18
% EndTime: 2020-05-01 20:43:23
% DurationCPUTime: 4.42s
% Computational Cost: add. (1248->312), mult. (1857->409), div. (0->0), fcn. (406->88), ass. (0->223)
t190 = cos(pkin(18));
t169 = t190 ^ 2;
t192 = mrSges(6,1) * pkin(9);
t183 = sin(qJ(4));
t199 = qJD(1) ^ 2;
t288 = t183 * t199;
t107 = t288 * t192;
t174 = Ifges(6,4) * t199;
t187 = cos(qJ(4));
t267 = t187 ^ 2 * t174;
t222 = -t174 + 0.2e1 * t267;
t285 = t187 * t199;
t181 = Ifges(6,1) - Ifges(6,2);
t289 = t181 * t183;
t328 = mrSges(6,2) * pkin(9);
t334 = t222 - t107 - (-t289 + t328) * t285;
t349 = t169 * t334;
t178 = mrSges(4,2) + mrSges(11,2);
t184 = sin(qJ(3));
t287 = t184 * t178;
t188 = cos(qJ(3));
t193 = -m(6) - m(5);
t224 = t193 * pkin(5) - mrSges(11,1) - mrSges(4,1);
t295 = t224 * t188;
t348 = t295 + t287;
t347 = 2 * qJD(1);
t165 = qJ(3) + pkin(19);
t201 = 2 * qJ(2);
t136 = t201 + t165;
t311 = mrSges(10,2) * sin(t136);
t313 = mrSges(10,1) * cos(t136);
t346 = (-t313 + t311) * pkin(2);
t159 = pkin(22) + pkin(21);
t248 = pkin(18) - t159;
t212 = -pkin(20) + t248;
t220 = -qJ(2) + t248;
t92 = -qJ(3) + t220;
t219 = qJ(2) + t248;
t93 = qJ(3) + t219;
t345 = (cos(t92) + cos(t93)) * mrSges(11,2);
t170 = qJ(3) + qJ(2);
t124 = 0.2e1 * t170;
t344 = (t193 * pkin(5) ^ 2 + Ifges(11,1) + Ifges(4,1) - Ifges(11,2) - Ifges(4,2)) * sin(t124);
t300 = mrSges(6,1) * pkin(11) - Ifges(6,5);
t293 = t300 * t183;
t148 = -mrSges(6,2) * pkin(11) + Ifges(6,6);
t294 = t148 * t187;
t342 = (-t293 + t294) * t199;
t152 = pkin(17) - pkin(18) + qJ(2);
t106 = 0.2e1 * t152;
t341 = (Ifges(3,4) + Ifges(10,4)) * cos(t201) + Ifges(7,4) * cos(t106);
t171 = t201 + qJ(3);
t291 = t178 * sin(t171);
t296 = t224 * cos(t171);
t340 = t296 + t291;
t339 = t183 * mrSges(6,1) + t187 * mrSges(6,2);
t137 = qJ(2) + t165;
t95 = 0.2e1 * t137;
t338 = (Ifges(11,4) + Ifges(4,4)) * cos(t124) + Ifges(9,4) * cos(t95);
t162 = qJD(3) + qJD(2);
t310 = mrSges(10,2) * sin(t165);
t312 = mrSges(10,1) * cos(t165);
t335 = pkin(2) * (t310 + t312);
t331 = pkin(1) / 0.2e1;
t330 = pkin(2) * mrSges(10,3);
t329 = m(11) * pkin(4);
t326 = pkin(4) * t199;
t197 = (qJD(4) ^ 2);
t132 = (2 * t197) + t199;
t325 = pkin(5) * t132;
t135 = -qJD(4) + t162;
t324 = pkin(5) * t135;
t323 = pkin(5) * t184;
t322 = pkin(5) * t188;
t321 = pkin(5) * t199;
t318 = qJD(2) / 0.2e1;
t317 = qJD(3) / 0.2e1;
t316 = -qJD(4) / 0.2e1;
t315 = mrSges(3,1) + mrSges(10,1);
t314 = pkin(5) * qJD(1);
t309 = pkin(15) * t199;
t308 = (m(10) * pkin(2) ^ 2 - Ifges(9,1) + Ifges(9,2)) * sin(t95);
t306 = t169 * t342;
t305 = (Ifges(7,1) - Ifges(7,2)) * sin(t106);
t129 = m(11) + m(4) + m(8) - t193;
t302 = (t129 * pkin(1) ^ 2 - Ifges(3,1) - Ifges(10,1) + Ifges(3,2) + Ifges(10,2)) * sin(t201);
t158 = t162 ^ 2;
t189 = cos(qJ(2));
t198 = qJD(2) ^ 2;
t301 = (pkin(1) * t198 + t158 * t323) * t189;
t299 = mrSges(6,2) * qJD(1);
t298 = pkin(14) * qJD(1);
t297 = pkin(15) * qJD(1);
t175 = sin(pkin(20));
t176 = cos(pkin(20));
t292 = t175 * t176;
t186 = sin(pkin(18));
t286 = t186 * t190;
t284 = -t197 + t198;
t155 = qJD(2) + t317;
t281 = qJD(3) * t155;
t280 = qJD(1) * qJD(4);
t277 = -pkin(22) + pkin(18);
t185 = sin(qJ(2));
t276 = t185 * t322;
t6 = t267 + (t289 / 0.2e1 - t328 / 0.2e1) * t285 - t107 / 0.2e1 - t174 / 0.2e1;
t275 = t6 * t286;
t254 = t300 * t288;
t255 = t148 * t285;
t274 = ((-t294 / 0.2e1 + t293 / 0.2e1) * t199 * t169 + t275 + t255 / 0.4e1 - t254 / 0.4e1) * t292;
t268 = (-t294 / 0.4e1 + t293 / 0.4e1) * t199 * t286;
t265 = -t326 / 0.2e1;
t264 = t325 / 0.4e1;
t263 = -t321 / 0.2e1;
t262 = -t321 / 0.4e1;
t261 = t321 / 0.2e1;
t256 = pkin(5) * t280;
t253 = qJD(1) * t331;
t252 = -t314 / 0.2e1;
t249 = t280 / 0.2e1;
t247 = t158 * t276;
t243 = 0.2e1 * t334;
t134 = qJD(4) + t162;
t234 = -t134 * t324 / 0.2e1;
t163 = qJD(2) + qJD(4);
t230 = t163 * t253;
t164 = qJD(2) - qJD(4);
t229 = t164 * t253;
t228 = t134 * t314 / 0.2e1;
t227 = t135 * t252;
t226 = 0.2e1 * t212;
t210 = t243 * t169 - 0.8e1 * t268 - t334;
t89 = -qJ(2) + t212;
t86 = qJ(2) + t212;
t85 = -qJ(4) + t212;
t84 = qJ(4) + t212;
t207 = t162 * t276 + t189 * (pkin(1) * qJD(2) + t162 * t323);
t73 = -qJ(2) + t85;
t72 = qJ(2) + t85;
t71 = -qJ(2) + t84;
t70 = qJ(2) + t84;
t111 = sin(t137);
t115 = cos(t137);
t150 = qJ(4) + t170;
t118 = sin(t150);
t151 = -qJ(4) + t170;
t119 = sin(t151);
t121 = cos(t150);
t122 = cos(t151);
t127 = pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3);
t140 = sin(t170);
t144 = cos(t170);
t153 = pkin(9) * m(6) + mrSges(5,1);
t154 = pkin(2) * m(10) + mrSges(9,1);
t61 = qJ(3) + t70;
t34 = sin(t61);
t62 = -qJ(3) + t71;
t35 = sin(t62);
t63 = qJ(3) + t72;
t36 = sin(t63);
t64 = -qJ(3) + t73;
t37 = sin(t64);
t38 = cos(t61);
t39 = cos(t62);
t40 = cos(t63);
t41 = cos(t64);
t74 = qJ(3) + t86;
t51 = sin(t74);
t75 = -qJ(3) + t89;
t52 = sin(t75);
t57 = cos(t74);
t58 = cos(t75);
t76 = sin(t92);
t77 = sin(t93);
t204 = t265 * t345 + (-t344 / 0.2e1 + t308 / 0.2e1) * t199 + (t52 * t261 + t51 * t263) * t153 + (t57 * t261 + t58 * t263) * t127 + (t76 * t326 / 0.2e1 + t77 * t265) * mrSges(11,1) + (mrSges(9,2) * t115 + t154 * t111 - t224 * t140 + t178 * t144) * t309 + (t38 * t262 + t40 * t321 / 0.4e1 + t39 * t264 - t41 * t325 / 0.4e1 + (t122 - t121) * t256) * mrSges(6,2) + ((-t118 - t119) * t256 + (t35 + t37) * t264 + (t36 + t34) * t262) * mrSges(6,1);
t200 = 0.2e1 * qJ(4);
t182 = mrSges(3,2) + mrSges(10,2);
t173 = qJ(2) - qJ(4);
t172 = qJ(2) + qJ(4);
t161 = t176 ^ 2;
t157 = mrSges(6,1) * t297;
t156 = mrSges(6,2) * t297;
t147 = cos(t173);
t146 = cos(t172);
t143 = sin(t173);
t142 = sin(t172);
t139 = -qJ(2) + t277;
t138 = qJ(2) + t277;
t126 = cos(t159);
t125 = sin(t159);
t123 = cos(t152);
t120 = sin(t152);
t117 = cos(t139);
t116 = cos(t138);
t113 = sin(t139);
t112 = sin(t138);
t105 = t158 + t199;
t101 = t300 - t328;
t100 = t300 + t328;
t99 = t148 + t192;
t98 = -t148 + t192;
t91 = cos(t220);
t90 = cos(t219);
t88 = -qJ(4) + t226;
t87 = qJ(4) + t226;
t69 = cos(t89);
t68 = cos(t86);
t67 = sin(t89);
t66 = sin(t86);
t60 = 0.2e1 * t85;
t59 = 0.2e1 * t84;
t56 = cos(t73);
t55 = cos(t72);
t54 = cos(t71);
t53 = cos(t70);
t50 = sin(t73);
t49 = sin(t72);
t48 = sin(t71);
t47 = sin(t70);
t3 = t339 * (-t105 * t189 * t322 + t185 * (t105 * t323 + pkin(1) * (t198 + t199)) - t309);
t2 = t186 * t3 + t339 * t190 * (t247 + t301);
t1 = t3 * t190 - 0.2e1 * t186 * t339 * (t247 / 0.2e1 + t301 / 0.2e1);
t4 = [-(t181 * sin(t200) + t98 * sin(t87) + t100 * cos(t88)) * t280 / 0.2e1 - 0.2e1 * (((Ifges(3,6) / 0.2e1 + Ifges(10,6) / 0.2e1) * qJD(2) + (t129 * pkin(1) + t315) * t297) * t189 + ((mrSges(7,2) * t298) + Ifges(7,5) * t318) * t120) * qJD(2) - 0.2e1 * (((Ifges(4,6) / 0.2e1 + Ifges(11,6) / 0.2e1) * qJD(2) + (Ifges(11,6) + Ifges(4,6)) * t317 - t224 * t297) * t140 + ((t330 / 0.2e1 - Ifges(9,5) / 0.2e1) * qJD(2) + mrSges(9,2) * t297 + (-Ifges(9,5) + t330) * t317) * t115 + (t154 * t297 + (t318 + t317) * Ifges(9,6)) * t111) * t162 + (0.2e1 * t182 * t297 + ((mrSges(11,3) + mrSges(4,3) + mrSges(5,3) + mrSges(8,3)) * pkin(1) - Ifges(3,5) - Ifges(10,5)) * qJD(2)) * qJD(2) * t185 + 0.2e1 * ((mrSges(7,1) * t298) - qJD(2) * Ifges(7,6) / 0.2e1) * qJD(2) * t123 + (t100 * t316 + t156) * qJD(4) * cos(t85) - qJD(4) * (t99 * t316 + t157) * sin(t85) + qJD(4) * (t98 * t316 + t157) * sin(t84) - (t163 * t47 + t164 * t48) * pkin(1) * t299 / 0.2e1 + (-sin(t60) / 0.4e1 + sin(t59) / 0.4e1) * t181 * t280 + (t99 * sin(t88) + t101 * cos(t87)) * t249 - qJD(4) * t339 * ((pkin(9) * qJD(1)) - qJD(4) * pkin(15)) + (t156 + qJD(4) * t101 / 0.2e1) * qJD(4) * cos(t84) + (-0.2e1 * t178 * t297 + (-mrSges(5,3) * pkin(5) + Ifges(11,5) + Ifges(4,5)) * t162) * t162 * t144 + (t340 * pkin(1) + t346) * t155 * t347 + (-t348 * pkin(1) + t335) * qJD(3) * qJD(1) + ((cos(t60) + cos(t59)) * t249 - cos(t200) * t280) * Ifges(6,4) + ((t122 + t121) * t234 + (t40 + t39) * t227 + (t41 + t38) * t228 + (t142 + t143) * t284 * t331 + t49 * t229 + t50 * t230) * mrSges(6,2) + ((t55 + t54) * t229 + (t56 + t53) * t230 + t134 * t37 * t252 + t118 * t234 + t35 * t227 + t34 * t228 + (t147 / 0.2e1 - t146 / 0.2e1) * pkin(1) * t284 + (qJD(1) * t36 + t134 * t119) * t324 / 0.2e1) * mrSges(6,1) + (t302 - t305 - 0.2e1 * t341 + ((t90 + t91) * t329 + (t68 + t69) * t153 + (t66 + t67) * t127 + (-t112 - t113) * mrSges(8,2) + (t116 + t117) * mrSges(8,1)) * pkin(1)) * qJD(2) * qJD(1) + (-t308 + t344 + 0.2e1 * t338 + (t345 + (-t76 + t77) * mrSges(11,1)) * pkin(4) + ((t51 - t52) * t153 + (-t57 + t58) * t127) * pkin(5)) * qJD(1) * t162; t204 + (-t302 / 0.2e1 + t305 / 0.2e1 + (-t182 * t185 + t315 * t189) * pkin(15) + (-mrSges(7,1) * t123 + mrSges(7,2) * t120) * pkin(14) - t346 + (t129 * pkin(15) * t189 + (-t69 / 0.2e1 - t68 / 0.2e1) * t153 + (-t66 / 0.2e1 - t67 / 0.2e1) * t127 + (t112 / 0.2e1 + t113 / 0.2e1) * mrSges(8,2) + (t47 / 0.4e1 - t49 / 0.4e1) * mrSges(6,2) + (-t116 / 0.2e1 - t117 / 0.2e1) * mrSges(8,1) + (-t53 / 0.4e1 - t55 / 0.4e1) * mrSges(6,1) + (-t90 / 0.2e1 - t91 / 0.2e1) * t329 - t340) * pkin(1) - t338 + t341) * t199 + 0.2e1 * t281 * t335 + (-0.2e1 * t348 * t281 + ((t48 / 0.4e1 - t50 / 0.4e1) * t132 + (t142 - t143) * t280) * mrSges(6,2) + ((-t54 / 0.4e1 - t56 / 0.4e1) * t132 + (-t146 - t147) * t280) * mrSges(6,1)) * pkin(1); t204 + ((-t310 / 0.2e1 - t312 / 0.2e1) * pkin(2) + (t287 / 0.2e1 + t295 / 0.2e1) * pkin(1)) * (0.2e1 * t198 + t199) + ((-t311 / 0.2e1 + t313 / 0.2e1) * pkin(2) + (-t291 / 0.2e1 - t296 / 0.2e1) * pkin(1) - t338) * t199; (-t1 * t175 + t2 * t176) * t125 + t210 * t161 + 0.4e1 * t274 - t349 + 0.4e1 * t268 + (mrSges(6,1) * t207 * t347 + t181 * t288) * t187 - 0.2e1 * t207 * t183 * t299 + t222 + (((-0.8e1 * t275 + 0.4e1 * t306 - 0.2e1 * t342) * t161 + 0.8e1 * (t6 * t169 - 0.2e1 * t268 - t267 / 0.2e1 + (-t289 / 0.4e1 + t328 / 0.4e1) * t285 + t107 / 0.4e1 + t174 / 0.4e1) * t292 - 0.2e1 * t306 + 0.4e1 * t275 + t255 - t254) * t125 + t1 * t176 + t2 * t175 + ((0.16e2 * t268 + t243 - 0.4e1 * t349) * t161 - 0.8e1 * t274 + t210) * t126) * t126;];
tauc = t4(:);
