% Calculate matrix of centrifugal and coriolis load on the joints for
% palh1m2DE1
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
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [11x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh1m2DE1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE1_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_coriolismatJ_fixb_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE1_coriolismatJ_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2DE1_coriolismatJ_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2DE1_coriolismatJ_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:57:39
% EndTime: 2020-05-01 20:57:49
% DurationCPUTime: 2.66s
% Computational Cost: add. (2811->291), mult. (3704->393), div. (0->0), fcn. (1368->88), ass. (0->185)
t302 = m(5) + m(6);
t151 = pkin(22) + pkin(21);
t228 = pkin(18) - t151;
t217 = -pkin(20) + t228;
t204 = qJ(4) + t217;
t195 = qJ(2) + t204;
t205 = -qJ(4) + t217;
t196 = -qJ(2) + t205;
t327 = sin(t195) / 0.4e1 - sin(t196) / 0.4e1;
t326 = -cos(t196) / 0.4e1 - cos(t195) / 0.4e1;
t155 = sin(pkin(20));
t156 = cos(pkin(20));
t161 = sin(pkin(18));
t165 = cos(pkin(18));
t163 = cos(qJ(3));
t164 = cos(qJ(2));
t262 = t163 * t164;
t129 = pkin(5) * t262;
t159 = sin(qJ(3));
t139 = pkin(5) * t159 + pkin(1);
t160 = sin(qJ(2));
t58 = t139 * t160 - t129;
t263 = t160 * t163;
t59 = pkin(5) * t263 + t139 * t164;
t18 = t161 * t59 - t165 * t58;
t236 = t161 * t58 + t59 * t165;
t325 = t155 * t18 + t156 * t236;
t87 = t159 * t164 + t263;
t273 = t161 * t87;
t264 = t159 * t160;
t86 = t262 - t264;
t235 = t86 * t165 + t273;
t271 = t165 * t87;
t27 = t161 * t86 - t271;
t324 = t155 * t27 + t156 * t235;
t323 = -t155 * t236 + t156 * t18;
t322 = -t155 * t235 + t156 * t27;
t320 = m(8) + m(11) + m(4) + t302;
t191 = qJ(3) + t195;
t192 = -qJ(3) + t196;
t315 = -sin(t191) / 0.4e1 + sin(t192) / 0.4e1;
t314 = -cos(t191) / 0.4e1 - cos(t192) / 0.4e1;
t152 = qJ(3) + pkin(19);
t237 = qJ(2) + t152;
t131 = sin(t237);
t132 = cos(t237);
t153 = qJ(4) + qJ(2);
t148 = qJ(3) + t153;
t133 = sin(t148);
t260 = qJ(3) + qJ(2);
t149 = -qJ(4) + t260;
t134 = sin(t149);
t136 = cos(t148);
t137 = cos(t149);
t143 = sin(t260);
t146 = cos(t260);
t304 = pkin(5) * m(6);
t250 = t304 / 0.2e1;
t251 = -t304 / 0.2e1;
t278 = rSges(11,3) * m(11);
t287 = rSges(5,3) * m(5);
t288 = rSges(4,3) * m(4);
t296 = m(9) * rSges(9,3);
t305 = pkin(2) * m(10);
t9 = (-pkin(5) * t287 - rSges(11,1) * t278 - rSges(4,1) * t288 + Icges(11,5) + Icges(4,5)) * t146 + (-rSges(9,1) * t296 - rSges(10,3) * t305 + Icges(9,5)) * t132 - (-rSges(9,2) * t296 + Icges(9,6)) * t131 + (t136 + t137) * rSges(6,2) * t251 + (t133 * t251 + t134 * t250) * rSges(6,1) + t143 * (rSges(11,2) * t278 + rSges(4,2) * t288 - Icges(11,6) - Icges(4,6));
t312 = (sin(t152) * rSges(10,2) + cos(t152) * rSges(10,1)) * t305;
t310 = -4 * Icges(6,5);
t309 = 0.2e1 * qJ(4);
t150 = pkin(17) - pkin(18) + qJ(2);
t308 = 0.2e1 * t150;
t307 = pkin(1) / 0.2e1;
t303 = pkin(9) * m(6);
t157 = m(11) * rSges(11,1);
t93 = pkin(5) * t302 + m(4) * rSges(4,1) + t157;
t301 = pkin(1) * t93;
t300 = m(6) * rSges(6,1);
t299 = m(6) * rSges(6,2);
t298 = m(7) * rSges(7,3);
t297 = m(9) * rSges(9,2);
t295 = rSges(3,1) * m(3);
t294 = rSges(6,1) * pkin(9);
t293 = rSges(10,1) * m(10);
t292 = rSges(3,2) * m(3);
t291 = rSges(4,2) * m(4);
t290 = rSges(6,2) * pkin(9);
t289 = rSges(10,2) * m(10);
t284 = cos(t150);
t162 = cos(qJ(4));
t272 = t162 * rSges(6,2);
t158 = sin(qJ(4));
t274 = t158 * rSges(6,1);
t119 = t272 + t274;
t283 = m(6) * t119;
t282 = m(11) * rSges(11,2);
t281 = m(6) * qJD(2);
t166 = -rSges(6,3) - pkin(11);
t280 = rSges(6,1) * t166;
t279 = rSges(6,2) * t166;
t125 = -qJ(2) + t204;
t108 = sin(t125);
t126 = qJ(2) + t205;
t109 = sin(t126);
t110 = cos(t125);
t111 = cos(t126);
t171 = 0.2e1 * qJ(2);
t240 = t171 + t152;
t117 = pkin(2) * sin(t240) * t289;
t135 = sin(t150);
t130 = t282 + t291;
t112 = -qJ(3) + t125;
t113 = qJ(3) + t126;
t247 = -t299 / 0.4e1;
t248 = -t300 / 0.4e1;
t255 = pkin(5) * t300;
t176 = sin(t113) * t255 / 0.4e1 + (sin(t112) * t248 + (cos(t112) + cos(t113)) * t247) * pkin(5);
t206 = qJ(2) + t217;
t197 = qJ(3) + t206;
t207 = -qJ(2) + t217;
t198 = -qJ(3) + t207;
t218 = qJ(2) + t228;
t209 = qJ(3) + t218;
t219 = -qJ(2) + t228;
t210 = -qJ(3) + t219;
t211 = 0.2e1 * t237;
t220 = m(5) * rSges(5,2) + m(6) * t166;
t222 = m(5) * rSges(5,1) + t303;
t223 = 0.2e1 * t260;
t249 = -pkin(4) * m(11) / 0.2e1;
t173 = -t176 + (pkin(2) ^ 2 * m(10) + (rSges(9,1) ^ 2 - rSges(9,2) ^ 2) * m(9) - Icges(9,1) + Icges(9,2)) * sin(t211) / 0.2e1 + (-Icges(11,1) + Icges(11,2) - Icges(4,1) + Icges(4,2) + (rSges(11,1) ^ 2 - rSges(11,2) ^ 2) * m(11) + (rSges(4,1) ^ 2 - rSges(4,2) ^ 2) * m(4) + t302 * pkin(5) ^ 2) * sin(t223) / 0.2e1 + pkin(4) * sin(t210) * t157 / 0.2e1 - (-rSges(9,1) * t297 + Icges(9,4)) * cos(t211) - (-rSges(11,1) * t282 - rSges(4,1) * t291 + Icges(11,4) + Icges(4,4)) * cos(t223) + t315 * t255 + t314 * pkin(5) * t299 + (-sin(t197) / 0.2e1 + sin(t198) / 0.2e1) * pkin(5) * t222 + (cos(t198) / 0.2e1 - cos(t197) / 0.2e1) * pkin(5) * t220 + (t130 * t146 + (m(9) * rSges(9,1) + t305) * t131 + t132 * t297 + t93 * t143) * pkin(15) + (rSges(11,1) * sin(t209) + (cos(t209) + cos(t210)) * rSges(11,2)) * t249;
t194 = pkin(2) * cos(t240) * t293;
t259 = t171 + qJ(3);
t213 = cos(t259) * t301;
t253 = pkin(18) - pkin(22);
t238 = qJ(2) + t253;
t239 = -qJ(2) + t253;
t90 = pkin(1) * t130 * sin(t259);
t1 = -t117 - t90 + t173 - (rSges(3,1) * t292 + rSges(10,1) * t289 - Icges(3,4) - Icges(10,4)) * cos(t171) - (m(7) * rSges(7,1) * rSges(7,2) - Icges(7,4)) * cos(t308) + pkin(1) * t109 * t247 + t213 + t194 + t326 * pkin(1) * t300 + (-rSges(7,1) * t284 + rSges(7,2) * t135) * pkin(14) * m(7) + (Icges(3,1) + Icges(10,1) - Icges(3,2) - Icges(10,2) + (-rSges(10,1) ^ 2 + rSges(10,2) ^ 2) * m(10) + (-rSges(3,1) ^ 2 + rSges(3,2) ^ 2) * m(3) - t320 * pkin(1) ^ 2) * sin(t171) / 0.2e1 + ((-rSges(7,1) ^ 2 + rSges(7,2) ^ 2) * m(7) + Icges(7,1) - Icges(7,2)) * sin(t308) / 0.2e1 + (cos(t219) + cos(t218)) * pkin(1) * t249 + (t111 + t110) * pkin(1) * t248 - (cos(t206) + cos(t207)) * pkin(1) * t222 / 0.2e1 + (sin(t206) + sin(t207)) * t220 * t307 + (-(t289 + t292) * t160 + (pkin(1) * t320 + t293 + t295) * t164) * pkin(15) + (t108 / 0.4e1 + t327) * pkin(1) * t299 + (-(cos(t238) + cos(t239)) * rSges(8,1) / 0.2e1 + (sin(t238) + sin(t239)) * rSges(8,2) / 0.2e1) * pkin(1) * m(8);
t277 = t1 * qJD(1);
t118 = rSges(6,1) * t162 - rSges(6,2) * t158;
t276 = t118 * t59;
t189 = 0.2e1 * t204;
t190 = 0.2e1 * t205;
t225 = 0.2e1 * t217;
t214 = qJ(4) + t225;
t215 = -qJ(4) + t225;
t257 = 0.4e1 * m(6);
t2 = (t274 / 0.2e1 + t272 / 0.2e1) * t303 + ((-cos(t205) / 0.2e1 - cos(t204) / 0.2e1) * rSges(6,2) + (sin(t205) / 0.2e1 - sin(t204) / 0.2e1) * rSges(6,1)) * m(6) * pkin(15) + (rSges(6,1) * t315 + rSges(6,2) * t314) * t304 + ((-t108 / 0.4e1 + t109 / 0.4e1 + t327) * rSges(6,2) + (t110 / 0.4e1 + t111 / 0.4e1 + t326) * rSges(6,1)) * m(6) * pkin(1) + t176 - ((t279 + t294) * m(6) + Icges(6,6)) * sin(t215) / 0.4e1 + ((-t279 + t294) * m(6) - Icges(6,6)) * sin(t214) / 0.4e1 + (cos(t189) / 0.4e1 + cos(t190) / 0.4e1 - cos(t309) / 0.2e1) * (rSges(6,1) * t299 - Icges(6,4)) + (-sin(t309) / 0.4e1 + sin(t189) / 0.8e1 - sin(t190) / 0.8e1) * ((rSges(6,1) ^ 2 - rSges(6,2) ^ 2) * m(6) - Icges(6,1) + Icges(6,2)) - ((-t280 - t290) * t257 + t310) * cos(t214) / 0.16e2 + ((-t280 + t290) * t257 + t310) * cos(t215) / 0.16e2;
t270 = t2 * qJD(1);
t265 = t130 * t159;
t3 = -t117 / 0.2e1 - t90 / 0.2e1 + t194 / 0.2e1 + t173 + t265 * t307 - t163 * t301 / 0.2e1 + t213 / 0.2e1 - t312 / 0.2e1;
t269 = t3 * qJD(1);
t268 = t9 * qJD(3);
t142 = sin(t153);
t154 = qJ(2) - qJ(4);
t144 = sin(t154);
t145 = cos(t153);
t147 = cos(t154);
t180 = (t137 / 0.4e1 - t136 / 0.4e1) * rSges(6,2) + (-t134 / 0.4e1 - t133 / 0.4e1) * rSges(6,1);
t177 = t180 * pkin(5);
t175 = t177 + ((-t144 / 0.4e1 + t142 / 0.4e1) * rSges(6,2) + (-t147 / 0.4e1 - t145 / 0.4e1) * rSges(6,1)) * pkin(1);
t11 = (-t276 / 0.2e1 + t175) * m(6);
t267 = t11 * qJD(1);
t201 = t118 * t87;
t13 = (-t201 / 0.2e1 + t180) * t304;
t266 = t13 * qJD(1);
t246 = m(6) * qJD(4) * t118;
t140 = sin(t151);
t141 = cos(t151);
t241 = (t140 * t322 + t324 * t141) * pkin(5) * t283;
t221 = -t241 / 0.2e1;
t37 = -(m(6) * t280 + Icges(6,5)) * t158 - (t279 * m(6) + Icges(6,6)) * t162;
t84 = pkin(9) * t283;
t216 = t161 * t37 + t165 * t84;
t15 = t312 + (t93 * t163 - t265) * pkin(1);
t178 = -t165 * t264 + t273;
t179 = t161 * t264 + t271;
t174 = (((t155 * t161 + t156 * t165) * t141 - (t155 * t165 - t156 * t161) * t140) * t129 + ((-t155 * t179 + t156 * t178) * t141 - (t155 * t178 + t156 * t179) * t140) * pkin(5)) * t283;
t4 = t221 + t174 / 0.2e1;
t193 = -qJD(2) * t15 + qJD(4) * t4;
t16 = -t161 * t84 + t165 * t37;
t14 = t15 * qJD(3);
t12 = m(6) * t177 + t201 * t250;
t10 = (t276 / 0.2e1 + t175) * m(6);
t5 = -t174 / 0.2e1 + t221;
t6 = [-qJD(2) * t1 - qJD(3) * t3 - qJD(4) * t2, -t277 + t268 + t10 * qJD(4) + ((t144 / 0.2e1 + t142 / 0.2e1) * rSges(6,2) + (t147 / 0.2e1 - t145 / 0.2e1) * rSges(6,1)) * pkin(1) * t281 + (t164 * (rSges(3,3) * t292 + rSges(10,3) * t289 - Icges(3,6) - Icges(10,6)) - (-rSges(7,1) * t298 + Icges(7,5)) * t135 + (rSges(7,2) * t298 - Icges(7,6)) * t284 - (-rSges(3,3) * t295 - rSges(10,3) * t293 + Icges(3,5) + Icges(10,5) + (-rSges(8,3) * m(8) - t278 - t287 - t288) * pkin(1)) * t160 + t9) * qJD(2), t9 * qJD(2) + t12 * qJD(4) + t268 - t269, -t270 + t10 * qJD(2) + t12 * qJD(3) + ((t155 * t16 - t156 * t216) * t141 + (t155 * t216 + t16 * t156) * t140 - (-pkin(15) + t58) * t283) * qJD(4); qJD(4) * t11 + t277, t14, t14 - t193, t267 - t4 * qJD(3) - (t323 * t140 + t325 * t141) * t246; qJD(4) * t13 + t269, t193, 0, t266 + t4 * qJD(2) - (t140 * t324 - t322 * t141) * pkin(5) * t246; -qJD(2) * t11 - qJD(3) * t13 + t270, -t267 - (-t325 * t140 + t323 * t141) * t119 * t281 + t5 * qJD(3), qJD(2) * t5 - qJD(3) * t241 - t266, 0;];
Cq = t6;
