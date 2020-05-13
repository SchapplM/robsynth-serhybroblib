% Calculate kinetic energy for
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m2DE1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE1_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_energykin_fixb_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE1_energykin_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE1_energykin_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2DE1_energykin_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:55:49
% EndTime: 2020-05-01 20:55:50
% DurationCPUTime: 1.08s
% Computational Cost: add. (872->146), mult. (1668->253), div. (0->0), fcn. (1800->22), ass. (0->71)
t257 = cos(pkin(20));
t262 = sin(pkin(18));
t267 = cos(pkin(18));
t282 = sin(pkin(20));
t232 = t257 * t262 - t267 * t282;
t235 = t267 * t257 + t262 * t282;
t252 = pkin(22) + pkin(21);
t247 = sin(t252);
t248 = cos(t252);
t216 = (t232 * t247 + t235 * t248) * qJD(1);
t260 = sin(qJ(3));
t265 = cos(qJ(3));
t219 = t232 * t265 - t235 * t260;
t222 = t232 * t260 + t235 * t265;
t261 = sin(qJ(2));
t266 = cos(qJ(2));
t210 = t219 * t266 - t222 * t261;
t211 = t219 * t261 + t222 * t266;
t273 = -t210 * t247 - t211 * t248;
t253 = qJD(2) + qJD(3);
t284 = pkin(1) * qJD(2);
t274 = pkin(5) * t253 + t260 * t284;
t276 = t210 * t248 - t211 * t247;
t278 = t265 * t284;
t205 = -t273 * t274 + t276 * t278;
t287 = t205 ^ 2;
t286 = m(11) / 0.2e1;
t215 = (-t232 * t248 + t235 * t247) * qJD(1);
t285 = t215 / 0.2e1;
t244 = pkin(1) * t261 - pkin(15);
t283 = cos(pkin(22));
t281 = sin(pkin(22));
t280 = qJD(2) ^ 2 * pkin(1) ^ 2;
t279 = (t260 ^ 2 + t265 ^ 2) * t280;
t229 = (-t260 * t261 + t265 * t266) * qJD(1);
t227 = -pkin(5) * t229 + t244 * qJD(1);
t277 = pkin(2) * t253;
t234 = -t262 * t283 + t267 * t281;
t270 = qJD(1) ^ 2;
t268 = cos(pkin(17));
t264 = cos(qJ(4));
t263 = sin(pkin(17));
t259 = sin(qJ(4));
t258 = cos(pkin(19));
t256 = sin(pkin(19));
t254 = qJ(3) + qJ(2);
t250 = cos(t254);
t249 = sin(t254);
t238 = t244 ^ 2 * t270;
t237 = t262 * t263 + t267 * t268;
t236 = t262 * t268 - t263 * t267;
t233 = -t256 * t260 + t258 * t265;
t231 = -t262 * t281 - t267 * t283;
t230 = t256 * t265 + t258 * t260;
t228 = (t260 * t266 + t261 * t265) * qJD(1);
t226 = t236 * t261 + t237 * t266;
t225 = -t236 * t266 + t237 * t261;
t224 = t233 * t277;
t223 = t230 * t277;
t221 = t231 * t261 - t234 * t266;
t220 = -t231 * t266 - t234 * t261;
t218 = (-t230 * t261 + t233 * t266) * qJD(1);
t217 = (t230 * t266 + t233 * t261) * qJD(1);
t214 = -pkin(2) * t218 - qJD(1) * pkin(15);
t213 = qJD(4) + t216;
t212 = (-cos(pkin(21)) * t231 - t234 * sin(pkin(21))) * pkin(4) + t244;
t209 = pkin(9) * t216 - pkin(11) * t215 + t227;
t207 = -t273 * t278 - t274 * t276;
t204 = t207 * t264 + t209 * t259;
t203 = -t207 * t259 + t209 * t264;
t1 = m(10) * (t214 ^ 2 + t223 ^ 2 + t224 ^ 2) / 0.2e1 + m(5) * (t207 ^ 2 + t227 ^ 2 + t287) / 0.2e1 + Ifges(4,2) * t229 ^ 2 / 0.2e1 + m(6) * (t203 ^ 2 + t204 ^ 2 + t287) / 0.2e1 + t279 * t286 + m(8) * (t238 + (t220 ^ 2 + t221 ^ 2) * t280) / 0.2e1 + Ifges(9,2) * t218 ^ 2 / 0.2e1 + m(4) * (t238 + t279) / 0.2e1 + (Ifges(4,4) * t229 + Ifges(4,1) * t228 / 0.2e1) * t228 + (Ifges(9,4) * t218 + Ifges(9,1) * t217 / 0.2e1) * t217 + (m(7) * pkin(14) ^ 2 / 0.2e1 + Ifges(2,3) / 0.2e1 + t212 ^ 2 * t286 + (m(3) / 0.2e1 + m(9) / 0.2e1) * pkin(15) ^ 2) * t270 - (-t227 * mrSges(5,1) + t207 * mrSges(5,3) - Ifges(5,2) * t216 / 0.2e1) * t216 + (t203 * mrSges(6,1) - t204 * mrSges(6,2) + Ifges(6,3) * t213 / 0.2e1) * t213 + (t227 * mrSges(5,2) - Ifges(5,4) * t216 + t205 * mrSges(5,3) + Ifges(5,1) * t285 + (Ifges(6,1) * t264 * t285 + t205 * mrSges(6,2) - t203 * mrSges(6,3) + Ifges(6,5) * t213) * t264 + (t205 * mrSges(6,1) - t204 * mrSges(6,3) - Ifges(6,6) * t213 + (-Ifges(6,4) * t264 + Ifges(6,2) * t259 / 0.2e1) * t215) * t259) * t215 + (Ifges(4,5) * t228 + Ifges(9,5) * t217 + Ifges(4,6) * t229 + Ifges(9,6) * t218 + (Ifges(11,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(9,3) / 0.2e1) * t253) * t253 + (t223 * mrSges(10,1) - t224 * mrSges(10,2) + (Ifges(3,3) / 0.2e1 + Ifges(10,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * qJD(2) + ((-t228 * t260 - t229 * t265) * mrSges(4,3) + ((mrSges(11,2) + mrSges(4,2)) * t265 + (mrSges(11,1) + mrSges(4,1)) * t260) * t253) * pkin(1)) * qJD(2) + (t244 * (-mrSges(4,1) * t229 + mrSges(4,2) * t228) - pkin(15) * (-mrSges(9,1) * t218 + mrSges(9,2) * t217) + ((-t212 * mrSges(11,1) + Ifges(11,2) * t250 / 0.2e1) * t250 + (mrSges(8,2) * t244 + Ifges(8,1) * t234 / 0.2e1) * t234 + (pkin(14) * mrSges(7,2) + Ifges(7,1) * t226 / 0.2e1) * t226 + (t212 * mrSges(11,2) + Ifges(11,4) * t250 + Ifges(11,1) * t249 / 0.2e1) * t249 + (-t244 * mrSges(8,1) + Ifges(8,4) * t234 + Ifges(8,2) * t231 / 0.2e1) * t231 + (pkin(14) * mrSges(7,1) - Ifges(7,4) * t226 + Ifges(7,2) * t225 / 0.2e1) * t225) * qJD(1) + (t214 * mrSges(10,2) - t223 * mrSges(10,3) + (-pkin(15) * mrSges(3,2) + (Ifges(10,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t266) * qJD(1)) * t266 + (t214 * mrSges(10,1) - t224 * mrSges(10,3) + (-pkin(15) * mrSges(3,1) + (Ifges(10,2) / 0.2e1 + Ifges(3,2) / 0.2e1) * t261 + (-Ifges(3,4) - Ifges(10,4)) * t266) * qJD(1)) * t261 + t253 * (Ifges(11,5) * t249 + Ifges(11,6) * t250) + (Ifges(7,5) * t226 - Ifges(7,6) * t225 + (Ifges(3,5) + Ifges(10,5)) * t266 + (-Ifges(3,6) - Ifges(10,6)) * t261 + ((t220 * t231 + t221 * t234) * mrSges(8,3) + (-t249 * t260 - t250 * t265) * mrSges(11,3)) * pkin(1)) * qJD(2)) * qJD(1);
T = t1;
