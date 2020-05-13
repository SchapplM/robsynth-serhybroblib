% Calculate time derivative of joint inertia matrix for
% palh3m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [9x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m2TE_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2TE_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_inertiaDJ_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2TE_inertiaDJ_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2TE_inertiaDJ_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2TE_inertiaDJ_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:47:33
% EndTime: 2020-05-07 01:47:34
% DurationCPUTime: 0.88s
% Computational Cost: add. (620->196), mult. (790->242), div. (0->0), fcn. (392->80), ass. (0->108)
t259 = pkin(17) + pkin(18);
t294 = (pkin(15) + t259);
t290 = (pkin(16) + t294);
t257 = qJD(3) + qJD(2);
t279 = m(5) + m(6);
t310 = pkin(8) * mrSges(6,2);
t255 = qJD(4) + qJD(2);
t309 = -t255 / 0.2e1;
t284 = 2 * qJ(2);
t308 = sin(t284);
t256 = qJD(4) - qJD(2);
t244 = qJD(3) - t256;
t307 = pkin(4) * t244;
t272 = sin(pkin(15));
t306 = pkin(8) * t272;
t276 = cos(pkin(15));
t305 = pkin(8) * t276;
t249 = pkin(10) * mrSges(6,2) - Ifges(6,6);
t270 = sin(qJ(3));
t271 = sin(qJ(2));
t303 = t270 * t271;
t275 = cos(qJ(2));
t302 = t270 * t275;
t260 = qJ(3) + qJ(2);
t301 = qJ(4) - qJ(2);
t233 = t279 * pkin(4) + mrSges(4,1) + mrSges(9,1);
t268 = mrSges(4,2) + mrSges(9,2);
t274 = cos(qJ(3));
t194 = (t233 * t270 + t268 * t274) * pkin(1) * qJD(3);
t237 = -pkin(4) * t274 + pkin(1);
t197 = pkin(4) * t303 + t237 * t275;
t300 = qJD(2) * t271;
t299 = qJD(2) * t275;
t297 = pkin(18) + pkin(15);
t296 = 2 * pkin(12);
t243 = qJD(3) + t255;
t295 = pkin(4) * t243 / 0.2e1;
t250 = pkin(10) * mrSges(6,1) - Ifges(6,5);
t293 = 2 * t290;
t292 = -qJ(2) + t294;
t291 = qJ(2) + t294;
t269 = sin(qJ(4));
t273 = cos(qJ(4));
t289 = t273 * mrSges(6,1) - t269 * mrSges(6,2);
t215 = t269 * mrSges(6,1) + t273 * mrSges(6,2);
t206 = t274 * t271 + t302;
t288 = -t274 * t275 + t303;
t225 = -qJ(2) + t290;
t224 = qJ(2) + t290;
t287 = mrSges(6,1) * t306 + t250 * t276;
t220 = -qJ(4) + t225;
t219 = -qJ(4) + t224;
t218 = qJ(4) + t225;
t217 = qJ(4) + t224;
t245 = sin(t260);
t246 = cos(t260);
t251 = qJ(4) + t260;
t252 = qJ(3) - t301;
t286 = (t245 * (Ifges(4,6) + Ifges(9,6)) + (pkin(4) * mrSges(5,3) - Ifges(4,5) - Ifges(9,5)) * t246) * t257 + (cos(t251) * t295 + cos(t252) * t307 / 0.2e1) * mrSges(6,2) + (sin(t251) * t295 - sin(t252) * t307 / 0.2e1) * mrSges(6,1);
t283 = 2 * qJ(4);
t278 = pkin(8) * mrSges(6,1);
t265 = cos(pkin(16));
t264 = sin(pkin(16));
t262 = qJ(2) + qJ(4);
t261 = t284 + qJ(3);
t254 = pkin(8) * m(6) + mrSges(5,1);
t253 = -pkin(15) + pkin(14) - qJ(2);
t248 = -qJ(2) + t297;
t247 = qJ(2) + t297;
t242 = cos(t259);
t241 = sin(t259);
t239 = -pkin(10) * m(6) + mrSges(5,2) - mrSges(6,3);
t238 = 0.2e1 * t260;
t236 = cos(t253);
t235 = sin(t253);
t232 = 0.2e1 * t253;
t231 = -qJ(3) + t292;
t230 = qJ(3) + t291;
t229 = -qJ(4) + t290;
t228 = qJ(4) + t290;
t227 = -qJ(4) + t293;
t226 = qJ(4) + t293;
t222 = -qJ(3) + t225;
t221 = qJ(3) + t224;
t214 = -qJ(3) + t220;
t213 = -qJ(3) + t218;
t212 = qJ(3) + t219;
t211 = qJ(3) + t217;
t210 = 2 * t229;
t209 = 2 * t228;
t208 = t289 * qJD(4);
t204 = mrSges(6,2) * t306 + t249 * t276;
t203 = mrSges(6,1) * t305 - t272 * t250;
t202 = mrSges(6,2) * t305 - t272 * t249;
t196 = -pkin(4) * t302 + t237 * t271;
t193 = t257 * t288;
t192 = t257 * t206;
t191 = t206 * t276 - t272 * t288;
t190 = -t272 * t206 - t276 * t288;
t189 = -t237 * t300 + (t206 * qJD(3) + t270 * t299) * pkin(4);
t188 = t237 * t299 + (t288 * qJD(3) + t270 * t300) * pkin(4);
t187 = -t196 * t272 + t197 * t276;
t186 = t196 * t276 + t197 * t272;
t185 = -t272 * t192 - t193 * t276;
t184 = -t192 * t276 + t272 * t193;
t183 = -t188 * t272 + t189 * t276;
t182 = t188 * t276 + t189 * t272;
t1 = [(((cos(t212) / 0.2e1 + cos(t213) / 0.2e1) * t244 + (-cos(t211) / 0.2e1 - cos(t214) / 0.2e1) * t243) * mrSges(6,2) + ((-sin(t212) / 0.2e1 + sin(t213) / 0.2e1) * t244 + (-sin(t211) / 0.2e1 + sin(t214) / 0.2e1) * t243) * mrSges(6,1)) * pkin(4) + ((t233 * sin(t261) + t268 * cos(t261)) * (0.2e1 * qJD(2) + qJD(3)) + ((cos(t218) / 0.2e1 + cos(t219) / 0.2e1) * t256 + (cos(t217) / 0.2e1 + cos(t220) / 0.2e1) * t255) * mrSges(6,2) + ((sin(t218) / 0.2e1 - sin(t219) / 0.2e1) * t256 + (sin(t217) / 0.2e1 - sin(t220) / 0.2e1) * t255) * mrSges(6,1)) * pkin(1) + ((t233 * t245 + t268 * t246) * t296 + ((-sin(t221) + sin(t222)) * t254 + (cos(t221) - cos(t222)) * t239) * pkin(4) - (pkin(4) ^ 2 * t279 - Ifges(4,1) - Ifges(9,1) + Ifges(4,2) + Ifges(9,2)) * sin(t238) + 0.2e1 * (Ifges(9,4) + Ifges(4,4)) * cos(t238) + ((-cos(t230) - cos(t231)) * mrSges(9,2) + (-sin(t230) + sin(t231)) * mrSges(9,1)) * pkin(3)) * t257 + (-t215 * pkin(8) + (sin(t283) / 0.2e1 - sin(t209) / 0.4e1 + sin(t210) / 0.4e1) * (Ifges(6,2) - Ifges(6,1)) + (-cos(t283) + cos(t209) / 0.2e1 + cos(t210) / 0.2e1) * Ifges(6,4) + ((cos(t228) + cos(t229)) * mrSges(6,2) + (sin(t228) - sin(t229)) * mrSges(6,1)) * pkin(12) + (-t250 - t310) * cos(t226) / 0.2e1 - (-t250 + t310) * cos(t227) / 0.2e1 - (t278 - t249) * sin(t226) / 0.2e1 + (t278 + t249) * sin(t227) / 0.2e1) * qJD(4) + ((Ifges(7,2) - Ifges(7,1)) * sin(t232) - (-Ifges(3,1) + Ifges(3,2)) * t308 + 0.2e1 * Ifges(3,4) * cos(t284) + 0.2e1 * Ifges(7,4) * cos(t232) + (-mrSges(3,1) * t271 - mrSges(3,2) * t275) * t296 + 0.2e1 * (-mrSges(7,1) * t235 + mrSges(7,2) * t236) * pkin(6) + ((sin(t224) - sin(t225)) * t254 + (-t308 * pkin(1) - 0.2e1 * pkin(12) * t271) * (m(4) + m(8) + m(9) + t279) + (cos(t225) - cos(t224)) * t239 + (-cos(t247) + cos(t248)) * mrSges(8,2) + (sin(t247) - sin(t248)) * mrSges(8,1) + (-sin(t292) + sin(t291)) * m(9) * pkin(3)) * pkin(1)) * qJD(2) + t194; (Ifges(3,5) * t275 + Ifges(7,5) * t236 - Ifges(3,6) * t271 + Ifges(7,6) * t235) * qJD(2) + ((-mrSges(4,3) - mrSges(5,3) - mrSges(8,3) - mrSges(9,3)) * t299 + (t256 * cos(-t301) / 0.2e1 + cos(t262) * t309) * mrSges(6,2) + (-t256 * sin(-t301) / 0.2e1 + sin(t262) * t309) * mrSges(6,1)) * pkin(1) + t286; 0.2e1 * t194; t286; t194; 0; -t189 * t289 + ((-(t203 * t265 - t264 * t287) * t269 + (-t202 * t265 + t204 * t264) * t273) * t242 + (-(-t264 * t203 - t287 * t265) * t269 + (t264 * t202 + t204 * t265) * t273) * t241 + (pkin(12) + t197) * t215) * qJD(4); -((t182 * t265 + t264 * t183) * t242 + (-t182 * t264 + t183 * t265) * t241) * t215 - ((t186 * t265 + t187 * t264) * t242 + (-t186 * t264 + t187 * t265) * t241) * t208; (((t184 * t264 + t185 * t265) * t242 + (t184 * t265 - t264 * t185) * t241) * t215 + ((t190 * t264 + t191 * t265) * t242 + (t190 * t265 - t264 * t191) * t241) * t208) * pkin(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
