% Calculate time derivative of joint inertia matrix for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m1OL_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1OL_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1OL_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:08:45
% EndTime: 2020-05-03 00:09:19
% DurationCPUTime: 6.09s
% Computational Cost: add. (1925->356), mult. (3187->503), div. (0->0), fcn. (600->8), ass. (0->200)
t145 = sin(qJ(4));
t130 = qJD(3) + qJD(4);
t115 = qJD(2) + t130;
t196 = mrSges(6,2) * t115;
t185 = pkin(4) * t196;
t120 = mrSges(6,1) * pkin(6) - Ifges(6,5);
t106 = qJD(5) * t120;
t141 = Ifges(6,1) - Ifges(6,2);
t144 = sin(qJ(5));
t252 = t141 * t144;
t321 = -t115 * t252 - t106;
t311 = t185 + t321;
t149 = cos(qJ(4));
t268 = mrSges(6,1) * qJD(5);
t125 = pkin(4) * t268;
t119 = mrSges(6,2) * pkin(6) - Ifges(6,6);
t131 = qJD(2) + qJD(3);
t337 = qJD(4) + t131;
t215 = t337 * t119;
t243 = qJD(5) * t144;
t36 = -0.4e1 * Ifges(6,4) * t243 - t125 - t215;
t338 = t36 * t149;
t11 = t145 * t311 + t338;
t123 = pkin(3) * t268;
t148 = cos(qJ(5));
t341 = t148 * (-t123 + t11);
t339 = t145 * t36;
t267 = mrSges(6,2) * qJD(5);
t124 = pkin(4) * t267;
t244 = qJD(5) * t141;
t183 = t115 * Ifges(5,6) + t244;
t261 = t115 * t120;
t173 = -t183 + (t124 - t261) * t144;
t283 = mrSges(6,1) * t115;
t186 = pkin(4) * t283;
t245 = qJD(5) * t119;
t73 = t186 + t245;
t86 = t115 * (Ifges(6,4) - Ifges(5,5));
t284 = t73 * t144 + t86;
t313 = t284 * t149;
t9 = -t173 * t145 + t313;
t336 = 0.2e1 * t145;
t212 = t141 * t243;
t198 = -t124 + t212;
t177 = t198 * t148;
t189 = t115 * t148;
t180 = t120 * t189;
t260 = t115 * t144;
t181 = t119 * t260;
t219 = m(6) * pkin(6) + mrSges(6,3);
t88 = t219 * pkin(4) + Ifges(5,4);
t195 = t115 * t88;
t142 = qJD(4) / 0.2e1;
t112 = t142 + t131;
t310 = t148 * mrSges(6,1) - mrSges(6,2) * t144;
t75 = -pkin(4) * m(6) - mrSges(5,1) - t310;
t232 = pkin(3) * t112 * t75;
t265 = qJD(5) * Ifges(6,4);
t307 = t125 * t144 + t265;
t133 = t148 ^ 2;
t111 = t133 * t265;
t96 = 0.2e1 * t111;
t329 = ((t96 + t177 - t307) * t145 - t232 + (t195 + t180 - t181) * t336) * t149;
t276 = t145 * t75;
t10 = t149 * t311 - t339;
t335 = t284 * t145;
t134 = t149 ^ 2;
t253 = t141 * t133;
t217 = t115 * t253;
t225 = Ifges(6,4) * t260;
t228 = (-0.4e1 * t225 - 0.2e1 * t245 - 0.4e1 * t186) * t148 + (-0.2e1 * t106 + 0.4e1 * t185) * t144 + 0.2e1 * t217;
t286 = pkin(3) * t145;
t109 = mrSges(6,1) * t243;
t241 = qJD(5) * t148;
t110 = mrSges(6,2) * t241;
t247 = t109 + t110;
t113 = -mrSges(5,2) + t219;
t299 = -0.2e1 * t113;
t48 = t112 * t299 + t247;
t231 = t48 * t286;
t152 = pkin(6) * mrSges(6,3);
t158 = pkin(4) ^ 2 * m(6);
t162 = pkin(6) ^ 2;
t285 = t162 * m(6);
t71 = -t158 / 0.2e1 + t285 / 0.2e1 + t152 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(6,3) / 0.2e1;
t266 = qJD(4) * t71;
t164 = pkin(3) ^ 2;
t168 = t164 * m(5) - Ifges(4,1) + Ifges(5,1) + Ifges(4,2) - Ifges(5,2) + Ifges(6,2) - Ifges(6,3) + 0.2e1 * t152 - t158;
t63 = (t162 + t164) * m(6) + t168;
t193 = -(-t106 / 0.2e1 + t185) * t144 + (t225 + t245 / 0.2e1 + t186) * t148;
t7 = (-t253 / 0.2e1 - t71) * t115 + t193;
t334 = 0.2e1 * qJD(3) * t63 + 0.8e1 * t7 * t134 + t228 - 0.2e1 * t231 + 0.4e1 * t266;
t293 = mrSges(6,2) * pkin(4);
t178 = (-t252 + t293) * t148;
t282 = mrSges(6,1) * t144;
t330 = -0.2e1 * t133;
t43 = -0.2e1 * (Ifges(6,4) * t330 + pkin(4) * t282 + Ifges(6,4) + t178) * qJD(5);
t220 = qJD(4) * t113 - t247;
t333 = (t276 * qJD(4) + t220 * t149) * pkin(3);
t287 = pkin(3) * t131;
t108 = mrSges(6,2) * t287;
t258 = t115 * t149;
t171 = Ifges(6,4) * t258 + t145 * t244;
t210 = t171 * t330 + t131 * (-Ifges(4,5) + (mrSges(5,3) + t282) * pkin(3));
t322 = t337 * t120;
t304 = (-t124 + t322) * t144 + t183;
t332 = t304 * t145 + (t108 + t10) * t148 + t210;
t331 = 0.2e1 * pkin(3);
t146 = sin(qJ(3));
t59 = t113 * t130 - t247;
t328 = t146 * t59;
t150 = cos(qJ(3));
t275 = t150 * t59;
t320 = -0.16e2 * t111 + 0.8e1 * t307 - 0.8e1 * t177;
t259 = t115 * t145;
t170 = Ifges(6,4) * t259 - t149 * t244;
t227 = pkin(3) * t267;
t179 = Ifges(4,6) * t131 + 0.2e1 * t170 * t133 - t144 * t227;
t317 = t304 * t149 + t179 - t335;
t295 = 2 * qJD(2);
t205 = t131 * t88;
t269 = -0.2e1 * t265 + 0.4e1 * t111;
t8 = t173 * t149 + t335;
t66 = -Ifges(6,5) * t115 + pkin(6) * t283 - t124;
t230 = t115 * t293;
t72 = -Ifges(6,5) * qJD(5) + pkin(6) * t268 - t230;
t190 = t115 * t119;
t74 = t125 + t190;
t303 = (t145 * t72 + t149 * t74 + t123) * t144 + t148 * (t145 * t73 - t149 * t66 + t227) + Ifges(6,3) * t259;
t302 = -2 * pkin(1);
t301 = -4 * pkin(2);
t51 = -t113 * t115 + t247;
t300 = -0.2e1 * t51;
t298 = 0.2e1 * t146;
t147 = sin(qJ(2));
t297 = 0.2e1 * t147;
t296 = 0.2e1 * t148;
t294 = 0.16e2 * t146;
t292 = pkin(3) * t48;
t291 = t220 * pkin(3);
t229 = (0.8e1 * t186 + 0.8e1 * t225 + 0.4e1 * t245) * t148 + (0.4e1 * t106 - 0.8e1 * t230) * t144 - 0.4e1 * t217;
t272 = t71 * t115;
t290 = (t229 - 0.8e1 * t272) * t145 + 0.2e1 * t292;
t289 = pkin(1) * t147;
t288 = pkin(2) * t147;
t280 = Ifges(6,4) * t144;
t279 = pkin(2) * qJD(2);
t278 = qJD(4) * pkin(3);
t270 = 0.2e1 * t198;
t84 = mrSges(6,2) * t148 + t282;
t264 = qJD(5) * t84;
t263 = t106 * t144;
t122 = qJD(2) + qJD(3) / 0.2e1;
t107 = t142 + t122;
t262 = t107 * t113;
t257 = t119 * t144;
t255 = t130 * t146;
t254 = t130 * t150;
t85 = -Ifges(4,4) + t88;
t251 = (Ifges(3,4) + t85) * qJD(2);
t250 = t85 * qJD(2);
t249 = t85 * qJD(3);
t82 = t88 * qJD(4);
t242 = qJD(5) * t146;
t240 = qJD(5) * t150;
t92 = (m(5) + m(6)) * pkin(3) + mrSges(4,1);
t237 = (t131 * t92 - t145 * t51) * t302;
t236 = ((-t110 / 0.2e1 - t109 / 0.2e1 + t262) * t145 + t92 * t122) * t301;
t234 = qJD(2) * t302;
t233 = pkin(1) * t115 * t75;
t226 = -0.16e2 * t146 * t134;
t223 = t149 * t278;
t222 = t75 * t259;
t221 = t75 * t255;
t213 = pkin(6) + t286;
t211 = 0.4e1 * pkin(2) * t107 * t75;
t207 = t145 * t232;
t204 = mrSges(6,1) * t223;
t203 = mrSges(6,2) * t223;
t202 = 0.4e1 * t265 + 0.4e1 * (-0.2e1 * t261 - t198) * t148 + (0.4e1 * t125 + 0.8e1 * t190) * t144 - 0.8e1 * t111 - 0.8e1 * t82;
t117 = -0.2e1 * t125;
t201 = (0.4e1 * t261 + t270) * t148 + (t117 - 0.4e1 * t190) * t144 + 0.4e1 * t82 + t269;
t200 = t111 - t265 / 0.2e1 + (t212 / 0.2e1 - t124 / 0.2e1 + t322) * t148 + (-t125 / 0.2e1 - t215) * t144 + t82;
t174 = (-mrSges(4,2) * t150 - t146 * t92) * qJD(3);
t172 = (t202 - 0.8e1 * t205) * t134 - 0.4e1 * t207 + 0.4e1 * t249 + t201;
t165 = pkin(2) ^ 2;
t153 = m(4) + m(5);
t151 = cos(qJ(2));
t135 = t150 ^ 2;
t89 = pkin(2) * t254;
t60 = pkin(2) * t221;
t47 = t247 - 0.2e1 * t262;
t46 = pkin(2) * t328;
t42 = mrSges(4,2) * t131 - t222;
t38 = mrSges(4,2) * t122 - t107 * t276;
t33 = -t230 - t321;
t5 = (-mrSges(6,1) * t287 - t145 * t66 - t149 * t73) * t148 + (t145 * t74 - t149 * t72 + t108) * t144 - Ifges(6,3) * t258;
t2 = (t145 * t33 + t123 - t338) * t148 + t317;
t1 = t313 + t332;
t3 = [(t290 * t149 + t172 + 0.4e1 * t250) * t135 + ((t233 * t297 + t279 * t299) * t149 + t147 * t237 + 0.2e1 * (mrSges(4,2) - t276) * t279 + (t63 * t295 + 0.4e1 * t329 + t334) * t146) * t150 + (0.4e1 * t205 + t201) * t134 + ((-t75 * t279 + t51 * t289) * t298 + (t228 + 0.4e1 * t272) * t145 + pkin(3) * t300) * t149 + (t42 * t289 + (t113 * t145 + t92) * t279) * t298 + (mrSges(3,1) + (m(6) + t153) * pkin(2)) * t147 * t234 + t222 * t331 + (-t212 - t261) * t296 + 0.2e1 * t115 * t257 - 0.2e1 * t251 - 0.2e1 * t249 - 0.2e1 * t82 - t269 + ((pkin(1) * t300 * t149 + t42 * t302) * t150 + (t47 * t288 + t233) * t298 * t149 + (0.4e1 * t38 * t288 + t237) * t146 + (mrSges(3,2) * t234) + (-0.8e1 * ((t119 * t241 + t263 + (t280 * t296 - t253 - 0.2e1 * t71) * t115 + 0.2e1 * (mrSges(6,1) * t189 - mrSges(6,2) * t260) * pkin(4)) * t134 + t329 - t231 / 0.2e1 + t217 / 0.2e1 + t266 - t193 + t131 * (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + (m(5) / 0.2e1 + m(6) / 0.2e1) * t164 + t71)) * t135 + ((t205 + t200) * t226 + ((t7 * t145 + t292 / 0.4e1) * t294 + t211) * t149 + 0.8e1 * (t131 * t85 + t200 - t207) * t146 + t236) * t150 + 0.4e1 * t329 + (Ifges(3,1) - Ifges(3,2) + t168 + t285 + (t164 - t165) * m(6) - t153 * t165) * t295 + t334) * t147 + (((0.16e2 * t180 - 0.16e2 * t181 + 0.16e2 * t195 - t320) * t134 + (-0.4e1 * t292 + (0.8e1 * t217 + (-0.16e2 * t225 - 0.8e1 * t245) * t148 - 0.8e1 * t263 + 0.16e2 * t272 + 0.16e2 * (t144 * t196 - t148 * t283) * pkin(4)) * t145) * t149 + 0.8e1 * t207 - 0.8e1 * t250 - 0.8e1 * t249 + t202) * t135 + (t7 * t226 + ((t320 * t145 + 0.8e1 * t232) * t146 - 0.2e1 * t47 * pkin(2) + (t144 * t190 - t148 * t261 - t195) * t145 * t294) * t149 + (-0.4e1 * t131 * t63 + t229 + 0.4e1 * t231 - 0.8e1 * t266) * t146 + t38 * t301) * t150 + (t146 * t211 + t290) * t149 + t146 * t236 + 0.4e1 * t251 + t172) * t151) * t151; (((t119 * t243 + t144 * t186 + t86) * t149 + t332) * t150 + t2 * t146 - (qJD(2) * Ifges(3,5))) * t151 + t147 * (t2 * t150 + ((t149 * t33 - t108 + t339) * t148 - t210 - t9) * t146 + (Ifges(3,6) * qJD(2))) + (qJD(2) * (mrSges(4,3) + mrSges(5,3) + t84) * t151 + t147 * t310 * qJD(5)) * pkin(2); t43 + 0.2e1 * t174 * pkin(2) + 0.2e1 * (pkin(2) * t275 + t291 + t60) * t149 + (-t46 + (t89 + t278) * t75) * t336; (t1 * t150 + t146 * (t317 - t341)) * t151 - ((t8 - t179 + t341) * t150 + t1 * t146) * t147; (t60 + 0.2e1 * t291) * t149 + (-t46 + (t89 + 0.2e1 * t278) * t75) * t145 + t43 + (t149 * t275 + t174) * pkin(2); t96 + 0.2e1 * (-(mrSges(6,1) * pkin(4) + t280) * t144 - t178) * qJD(5) + (-t247 * t149 + (t113 * t149 + t276) * qJD(4)) * t331; (0.2e1 * (t146 * t170 - t150 * t171) * t151 + (t146 * t171 + t150 * t170) * t297) * t133 + ((t10 * t150 - t11 * t146) * t151 - t147 * (t10 * t146 + t11 * t150)) * t148 + (-t146 * t8 + t150 * t9) * t151 - (t146 * t9 + t150 * t8) * t147; t43 + t333 + ((t221 + t275) * t149 + (t75 * t254 - t328) * t145) * pkin(2); t117 * t144 + t270 * t148 + t269 + t333; t43; (-pkin(2) * t264 + t5 * t146 - t150 * t303) * t151 + (t146 * t303 + t5 * t150 - t279 * t310) * t147 - pkin(1) * t264; (-t123 * t145 - t106 - t203) * t148 - (-t145 * t227 + t204 - t245) * t144 + (((-mrSges(6,1) * t242 - mrSges(6,2) * t254) * t149 + (-mrSges(6,1) * t240 + mrSges(6,2) * t255) * t145) * t148 - ((mrSges(6,1) * t254 - mrSges(6,2) * t242) * t149 + (-mrSges(6,1) * t255 - mrSges(6,2) * t240) * t145) * t144) * pkin(2); (-t203 - qJD(5) * (mrSges(6,1) * t213 - Ifges(6,5))) * t148 + (-t204 + qJD(5) * (mrSges(6,2) * t213 - Ifges(6,6))) * t144; -qJD(5) * (t120 * t148 - t257); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
