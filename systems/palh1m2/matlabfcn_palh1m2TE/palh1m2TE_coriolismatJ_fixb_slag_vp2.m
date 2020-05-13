% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh1m2TE_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_coriolismatJ_fixb_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2TE_coriolismatJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2TE_coriolismatJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2TE_coriolismatJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:26:24
% EndTime: 2020-05-01 20:26:37
% DurationCPUTime: 3.87s
% Computational Cost: add. (1630->357), mult. (1766->474), div. (0->0), fcn. (470->105), ass. (0->226)
t166 = pkin(22) + pkin(21);
t241 = pkin(18) - t166;
t121 = -pkin(20) + t241;
t197 = sin(qJ(3));
t302 = pkin(5) * t197;
t131 = pkin(1) + t302;
t136 = -pkin(11) * m(6) + mrSges(5,2) - mrSges(6,3);
t163 = pkin(9) * m(6) + mrSges(5,1);
t171 = pkin(18) - pkin(22);
t198 = sin(qJ(2));
t201 = cos(qJ(3));
t202 = cos(qJ(2));
t196 = sin(qJ(4));
t200 = cos(qJ(4));
t290 = mrSges(6,1) * t200;
t230 = mrSges(6,2) * t196 - t290;
t271 = t201 * t198;
t265 = pkin(5) * t271;
t341 = -((-t230 + t163) * cos(t121) - t136 * sin(t121)) * (t131 * t202 + t265) - (mrSges(8,1) * cos(t171) - mrSges(8,2) * sin(t171)) * pkin(1) * t202 - (t198 * (mrSges(11,1) * t201 - mrSges(11,2) * t197) + t202 * (pkin(1) * m(11) + mrSges(11,1) * t197 + mrSges(11,2) * t201)) * pkin(4) * cos(t241);
t182 = (pkin(17) - pkin(18));
t134 = 2 * t182;
t210 = 0.2e1 * qJ(3);
t173 = sin(t210);
t211 = 0.2e1 * qJ(2);
t174 = sin(t211);
t176 = cos(t210);
t177 = cos(t211);
t194 = Ifges(7,1) - Ifges(7,2);
t299 = qJD(1) / 0.2e1;
t305 = t177 / 0.2e1;
t306 = t174 / 0.2e1;
t137 = pkin(2) ^ 2 * m(10) - Ifges(9,1) + Ifges(9,2);
t315 = 2 * pkin(19);
t267 = sin(t315);
t268 = cos(t315);
t329 = 2 * Ifges(9,4);
t204 = m(5) + m(6);
t78 = -pkin(5) ^ 2 * t204 + Ifges(11,1) + Ifges(4,1) - Ifges(11,2) - Ifges(4,2);
t5 = (t137 * t268 + t267 * t329 - t78) * qJD(1);
t190 = (Ifges(4,4) + Ifges(11,4));
t7 = (-t137 * t267 + t268 * t329 + (2 * t190)) * qJD(1);
t340 = -(t173 * t5 - t176 * t7) * t305 - (t173 * t7 + t176 * t5) * t306 - ((-0.2e1 * Ifges(7,4) * t174 + t177 * t194) * sin(t134) + (0.2e1 * Ifges(7,4) * t177 + t174 * t194) * cos(t134)) * t299;
t167 = qJD(3) + qJD(2);
t185 = sin(pkin(19));
t187 = cos(pkin(19));
t189 = Ifges(4,5) + Ifges(11,5);
t277 = qJD(1) * pkin(15);
t27 = mrSges(9,2) * t277 + t167 * (pkin(2) * mrSges(10,3) - Ifges(9,5));
t164 = pkin(2) * m(10) + mrSges(9,1);
t28 = Ifges(9,6) * t167 + t164 * t277;
t69 = mrSges(6,1) * t196 + mrSges(6,2) * t200;
t328 = mrSges(5,3) + t69;
t337 = t185 * t28 + t187 * t27 + (pkin(5) * t328 - t189) * t167;
t191 = mrSges(11,2) + mrSges(4,2);
t272 = t197 * t191;
t324 = pkin(5) * t204 + mrSges(11,1) + mrSges(4,1);
t275 = t324 * t201;
t336 = -t272 + t275;
t178 = qJ(3) + qJ(2);
t149 = sin(t178);
t154 = cos(t178);
t335 = -t149 * t324 - t154 * t191;
t170 = qJ(3) + pkin(19);
t147 = qJ(2) + t170;
t123 = sin(t147);
t125 = cos(t147);
t234 = -qJ(2) + t241;
t103 = -qJ(3) + t234;
t233 = qJ(2) + t241;
t104 = qJ(3) + t233;
t107 = 0.2e1 * t147;
t133 = 0.2e1 * t178;
t214 = -((sin(t103) / 0.2e1 - sin(t104) / 0.2e1) * mrSges(11,1) - (cos(t103) / 0.2e1 + cos(t104) / 0.2e1) * mrSges(11,2)) * pkin(4) - t137 * sin(t107) / 0.2e1 + t190 * cos(t133) + Ifges(9,4) * cos(t107) + t78 * sin(t133) / 0.2e1;
t334 = (mrSges(9,2) * t125 + t123 * t164 - t335) * pkin(15) - t214;
t319 = 0.2e1 * t197;
t188 = Ifges(11,6) + Ifges(4,6);
t332 = t167 * t188;
t148 = t211 + t170;
t124 = sin(t148);
t126 = cos(t148);
t139 = sin(t170);
t141 = cos(t170);
t179 = t211 + qJ(3);
t150 = sin(t179);
t155 = cos(t179);
t330 = ((t201 / 0.2e1 - t155 / 0.2e1) * t324 - (t197 / 0.2e1 - t150 / 0.2e1) * t191) * pkin(1) + ((t141 / 0.2e1 - t126 / 0.2e1) * mrSges(10,1) + (t139 / 0.2e1 + t124 / 0.2e1) * mrSges(10,2)) * pkin(2);
t327 = pkin(1) * qJD(2);
t161 = -qJ(4) + t178;
t130 = cos(t161);
t160 = qJ(4) + t178;
t129 = cos(t160);
t287 = mrSges(6,2) * t129;
t128 = sin(t161);
t292 = mrSges(6,1) * t128;
t127 = sin(t160);
t293 = mrSges(6,1) * t127;
t323 = (-mrSges(6,2) * t130 + t287 + t292 + t293) * qJD(1);
t101 = qJ(2) + t121;
t80 = qJ(3) + t101;
t102 = -qJ(2) + t121;
t81 = -qJ(3) + t102;
t217 = (cos(t80) / 0.2e1 - cos(t81) / 0.2e1) * t136 + (sin(t80) / 0.2e1 - sin(t81) / 0.2e1) * t163;
t318 = 0.2e1 * t201;
t317 = pkin(5) / 0.2e1;
t316 = -pkin(9) / 0.2e1;
t66 = -qJ(4) + t81;
t38 = sin(t66);
t314 = t38 / 0.2e1;
t42 = cos(t66);
t313 = -t42 / 0.2e1;
t98 = -qJ(4) + t121;
t85 = -qJ(2) + t98;
t53 = sin(t85);
t312 = -t53 / 0.2e1;
t59 = cos(t85);
t311 = -t59 / 0.2e1;
t310 = mrSges(6,1) * pkin(9);
t309 = mrSges(6,2) * pkin(9);
t308 = -t129 / 0.2e1;
t307 = t130 / 0.2e1;
t304 = pkin(1) * t191;
t298 = qJD(1) / 0.4e1;
t297 = qJD(3) / 0.2e1;
t296 = Ifges(3,4) + Ifges(10,4);
t295 = pkin(5) * qJD(4);
t119 = pkin(9) * qJD(4) - t277;
t294 = mrSges(6,1) * t119;
t288 = mrSges(10,1) * t185;
t284 = mrSges(10,2) * t187;
t209 = 0.2e1 * qJ(4);
t175 = cos(t209);
t283 = Ifges(6,4) * t175;
t281 = pkin(2) * mrSges(10,1) * t187 + pkin(1) * t324;
t280 = mrSges(6,1) * qJD(1);
t279 = mrSges(6,2) * qJD(4);
t278 = pkin(14) * qJD(1);
t172 = sin(t209);
t193 = -Ifges(6,2) + Ifges(6,1);
t273 = t193 * t172;
t270 = t202 * t201;
t266 = pkin(2) * mrSges(10,2) * t185;
t63 = qJ(4) + t80;
t35 = sin(t63);
t263 = t35 / 0.4e1 - t38 / 0.4e1;
t64 = -qJ(4) + t80;
t36 = sin(t64);
t65 = qJ(4) + t81;
t37 = sin(t65);
t262 = t36 / 0.4e1 - t37 / 0.4e1;
t39 = cos(t63);
t261 = t39 / 0.4e1 + t42 / 0.4e1;
t40 = cos(t64);
t41 = cos(t65);
t260 = -t40 / 0.4e1 - t41 / 0.4e1;
t97 = qJ(4) + t121;
t82 = qJ(2) + t97;
t258 = -sin(t82) / 0.4e1 + t53 / 0.4e1;
t84 = -qJ(2) + t97;
t52 = sin(t84);
t83 = qJ(2) + t98;
t257 = sin(t83) / 0.4e1 - t52 / 0.4e1;
t58 = cos(t84);
t255 = cos(t83) / 0.4e1 + t58 / 0.4e1;
t254 = t59 / 0.4e1 + cos(t82) / 0.4e1;
t251 = pkin(5) * t298;
t250 = -t292 / 0.2e1;
t249 = pkin(11) * mrSges(6,1) - Ifges(6,5);
t159 = -pkin(11) * mrSges(6,2) + Ifges(6,6);
t248 = t307 + t308;
t237 = mrSges(6,2) * t251;
t238 = mrSges(6,1) * t251;
t239 = -pkin(5) * t280 / 0.4e1;
t242 = t35 * t238 + t38 * t239 + (t39 + t42) * t237;
t236 = 0.2e1 * t121;
t138 = m(11) + m(4) + m(8) + t204;
t68 = t138 * pkin(1) + mrSges(3,1) + mrSges(10,1);
t235 = pkin(1) ^ 2 * t138 - Ifges(3,1) - Ifges(10,1) + Ifges(3,2) + Ifges(10,2);
t228 = t185 * t27 - t187 * t28;
t195 = mrSges(3,2) + mrSges(10,2);
t227 = t195 * t198 - t202 * t68;
t226 = pkin(5) * t279 * t307 + t37 * t238 + t36 * t239 + t250 * t295 + (t40 + t41) * t237;
t181 = qJ(2) - qJ(4);
t152 = sin(t181);
t157 = cos(t181);
t223 = mrSges(6,1) * t157 / 0.2e1 + mrSges(6,2) * t152 / 0.2e1;
t222 = t287 / 0.2e1 + t293 / 0.2e1;
t220 = t167 * t265 + t202 * (t167 * t302 + t327);
t100 = -qJ(4) + t236;
t61 = 0.2e1 * t97;
t62 = 0.2e1 * t98;
t99 = qJ(4) + t236;
t213 = t69 * t316 - (-t159 + t310) * sin(t99) / 0.4e1 + (t159 + t310) * sin(t100) / 0.4e1 - (t249 + t309) * cos(t100) / 0.4e1 - (-t249 + t309) * cos(t99) / 0.4e1 + (-sin(t62) / 0.8e1 - t172 / 0.4e1 + sin(t61) / 0.8e1) * t193 + (-t175 / 0.2e1 + cos(t61) / 0.4e1 + cos(t62) / 0.4e1) * Ifges(6,4);
t203 = cos(pkin(18));
t199 = sin(pkin(18));
t186 = cos(pkin(20));
t184 = sin(pkin(20));
t180 = qJ(2) + qJ(4);
t169 = qJD(2) - qJD(4);
t168 = qJD(2) + qJD(4);
t165 = qJD(2) + t297;
t162 = qJ(2) + t182;
t158 = cos(t182);
t156 = cos(t180);
t153 = sin(t182);
t151 = sin(t180);
t146 = -qJ(2) + t171;
t145 = qJ(2) + t171;
t144 = -qJD(4) + t167;
t143 = qJD(4) + t167;
t118 = 0.2e1 * t162;
t112 = 0.8e1 * qJD(4) * t159;
t96 = mrSges(7,1) * t278 - qJD(2) * Ifges(7,6);
t95 = mrSges(7,2) * t278 + Ifges(7,5) * qJD(2);
t79 = 0.2e1 * t121;
t75 = cos(t98);
t74 = cos(t97);
t71 = sin(t98);
t70 = sin(t97);
t20 = t304 + (t284 + t288) * pkin(2);
t8 = -t266 + t281;
t6 = pkin(5) * t293 + pkin(1) * (mrSges(6,1) * t156 - mrSges(6,2) * t151);
t4 = t230 * t295 + (-mrSges(9,2) * t185 + t164 * t187 + t324) * t277;
t3 = (t266 + t281) * t201 - (t304 + (-t284 + t288) * pkin(2)) * t197;
t2 = -t197 * t20 + t201 * t8 + t296;
t1 = t20 * t318 + t319 * t8 + t235;
t9 = [((t165 * t124 + t139 * t297) * mrSges(10,2) + (-t126 * t165 + t141 * t297) * mrSges(10,1)) * pkin(2) + (t235 * t306 - t296 * t177 - t194 * sin(t118) / 0.2e1 - Ifges(7,4) * cos(t118) + t227 * pkin(15) + (mrSges(7,1) * cos(t162) - mrSges(7,2) * sin(t162)) * pkin(14)) * qJD(2) + ((t143 * t261 + t144 * t260) * mrSges(6,2) + (t143 * t263 + t144 * t262) * mrSges(6,1)) * pkin(5) + (((t74 / 0.2e1 + t75 / 0.2e1) * mrSges(6,2) + (t70 / 0.2e1 - t71 / 0.2e1) * mrSges(6,1)) * pkin(15) + t213) * qJD(4) + (pkin(5) * t217 - t334) * t167 + ((t191 * t150 - t155 * t324) * t165 + (t275 / 0.2e1 - t272 / 0.2e1) * qJD(3) + (t168 * t258 + t169 * t257) * mrSges(6,2) + (t168 * t254 + t169 * t255) * mrSges(6,1) + ((cos(t101) / 0.2e1 + cos(t102) / 0.2e1) * t163 + (-sin(t101) / 0.2e1 - sin(t102) / 0.2e1) * t136 + (-sin(t146) / 0.2e1 - sin(t145) / 0.2e1) * mrSges(8,2) + (cos(t145) / 0.2e1 + cos(t146) / 0.2e1) * mrSges(8,1) + (cos(t234) / 0.2e1 + cos(t233) / 0.2e1) * m(11) * pkin(4)) * qJD(2)) * pkin(1), (t96 * t158 - t95 * t153 - qJD(2) * (Ifges(3,6) + Ifges(10,6)) - t68 * t277 - (t277 * t324 - t228 + t332) * t197 - (t191 * t277 + t337) * t201) * t202 + (-t95 * t158 - t96 * t153 - qJD(2) * (Ifges(3,5) + Ifges(10,5)) - t201 * t332 + (mrSges(4,3) + mrSges(8,3) + mrSges(11,3) + t328) * t327 + t228 * t318 / 0.2e1 + (t195 - t336) * t277 + t337 * t319 / 0.2e1) * t198 + (t1 * t306 - 0.2e1 * t2 * t305 - t341) * qJD(1) + t340, -t28 * t123 - t27 * t125 + (((t308 - t130 / 0.2e1) * mrSges(6,2) + (-t127 / 0.2e1 + t128 / 0.2e1) * mrSges(6,1)) * pkin(5) - t188 * t149 + (-pkin(5) * mrSges(5,3) + t189) * t154) * t167 + ((mrSges(6,1) * t262 + mrSges(6,2) * t260 + t217) * pkin(5) + t214 + t330) * qJD(1) + t242 + t335 * t277, -(t74 + t75) * t119 * mrSges(6,2) / 0.2e1 + (t69 * pkin(15) + t222 * pkin(5) + ((-t151 / 0.2e1 - t152 / 0.2e1) * mrSges(6,2) + (t156 / 0.2e1 - t157 / 0.2e1) * mrSges(6,1)) * pkin(1) + (t74 - t75) * t249 / 0.2e1) * qJD(4) + (((-t257 + t258) * mrSges(6,2) + (t254 - t255) * mrSges(6,1)) * pkin(1) + t213) * qJD(1) + t226 + t242 + (t112 + 0.8e1 * t294) * t71 / 0.16e2 + (t112 - 0.8e1 * t294) * t70 / 0.16e2; t4 * t271 + (t197 * t4 + (-qJD(4) * t290 + t196 * t279) * pkin(1)) * t202 + (t2 * t177 - t1 * t174 / 0.2e1 + (-(mrSges(7,1) * t202 - mrSges(7,2) * t198) * t158 + (mrSges(7,1) * t198 + mrSges(7,2) * t202) * t153) * pkin(14) + ((-t197 * t198 + t270) * (mrSges(9,2) * t187 + t164 * t185 + t191) - t227) * pkin(15) + t341) * qJD(1) - t340, qJD(3) * (t336 * pkin(1) + (mrSges(10,1) * t141 + mrSges(10,2) * t139) * pkin(2)), t167 * t3, (-t6 / 0.2e1 - t223 * pkin(1) + (mrSges(6,2) * t248 + t250) * pkin(5)) * qJD(1) + (((t313 + t41 / 0.2e1) * mrSges(6,2) + (t314 + t37 / 0.2e1) * mrSges(6,1)) * pkin(5) + ((t52 / 0.2e1 + t312) * mrSges(6,2) + (t311 - t58 / 0.2e1) * mrSges(6,1)) * pkin(1)) * qJD(4); -t222 * t295 + ((-mrSges(6,1) * t263 - mrSges(6,2) * t261 - t217) * pkin(5) - t330 + t334) * qJD(1) + t226, -t3 * qJD(2), 0, (-t323 + ((t41 - t42) * mrSges(6,2) + (t37 + t38) * mrSges(6,1)) * qJD(4)) * t317; t220 * t290 - t196 * (mrSges(6,2) * t220 + t280 * t316) + ((-t159 * t200 + t196 * t249) * sin(t79) + t283) * t299 + ((0.2e1 * t69 * pkin(9) - t273 - 0.2e1 * t283) * cos(t79) + t273) * t298 + (t309 * t200 / 0.2e1 + ((t184 * t199 + t186 * t203) * cos(t166) - (t184 * t203 - t186 * t199) * sin(t166)) * t69 * (-pkin(5) * t270 + t131 * t198 - pkin(15))) * qJD(1), t6 * t299 + (t223 * qJD(1) + ((-t52 / 0.2e1 + t312) * mrSges(6,2) + (t311 + t58 / 0.2e1) * mrSges(6,1)) * qJD(2)) * pkin(1) + ((t128 * t299 + (-t37 / 0.2e1 + t314) * t167) * mrSges(6,1) + ((-t41 / 0.2e1 + t313) * t167 - t248 * qJD(1)) * mrSges(6,2)) * pkin(5), (t323 + ((-t41 - t42) * mrSges(6,2) + (-t37 + t38) * mrSges(6,1)) * t167) * t317, 0;];
Cq = t9;
