% Calculate kinetic energy for
% palh1m2DE2
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
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m2DE2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_energykin_fixb_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_energykin_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE2_energykin_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2DE2_energykin_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:57:56
% EndTime: 2020-05-02 20:57:58
% DurationCPUTime: 1.19s
% Computational Cost: add. (872->146), mult. (1668->253), div. (0->0), fcn. (1800->22), ass. (0->71)
t256 = cos(pkin(20));
t261 = sin(pkin(18));
t266 = cos(pkin(18));
t280 = sin(pkin(20));
t230 = t256 * t261 - t266 * t280;
t234 = t266 * t256 + t261 * t280;
t249 = pkin(22) + pkin(21);
t244 = sin(t249);
t245 = cos(t249);
t215 = (t230 * t244 + t234 * t245) * qJD(1);
t259 = sin(qJ(3));
t264 = cos(qJ(3));
t218 = t230 * t264 - t234 * t259;
t221 = t230 * t259 + t234 * t264;
t260 = sin(qJ(2));
t265 = cos(qJ(2));
t209 = t218 * t265 - t260 * t221;
t210 = t260 * t218 + t221 * t265;
t272 = -t209 * t244 - t210 * t245;
t250 = qJD(2) + qJD(3);
t281 = pkin(1) * qJD(2);
t273 = pkin(5) * t250 + t259 * t281;
t275 = t209 * t245 - t244 * t210;
t277 = t264 * t281;
t204 = -t272 * t273 + t275 * t277;
t284 = t204 ^ 2;
t283 = m(11) / 0.2e1;
t214 = (-t230 * t245 + t234 * t244) * qJD(1);
t282 = t214 / 0.2e1;
t241 = t260 * pkin(1) - pkin(15);
t279 = qJD(2) ^ 2 * pkin(1) ^ 2;
t278 = (t259 ^ 2 + t264 ^ 2) * t279;
t228 = (-t259 * t260 + t264 * t265) * qJD(1);
t226 = -pkin(5) * t228 + t241 * qJD(1);
t276 = pkin(2) * t250;
t252 = sin(pkin(22));
t255 = cos(pkin(22));
t233 = t252 * t266 - t261 * t255;
t269 = qJD(1) ^ 2;
t267 = cos(pkin(17));
t263 = cos(qJ(4));
t262 = sin(pkin(17));
t258 = sin(qJ(4));
t257 = cos(pkin(19));
t254 = sin(pkin(19));
t251 = qJ(2) + qJ(3);
t247 = cos(t251);
t246 = sin(t251);
t237 = t241 ^ 2 * t269;
t236 = t261 * t262 + t266 * t267;
t235 = t261 * t267 - t262 * t266;
t232 = t252 * t261 + t255 * t266;
t231 = -t254 * t259 + t257 * t264;
t229 = t254 * t264 + t257 * t259;
t227 = (t259 * t265 + t260 * t264) * qJD(1);
t225 = t235 * t260 + t236 * t265;
t224 = -t235 * t265 + t236 * t260;
t223 = t231 * t276;
t222 = t229 * t276;
t220 = -t232 * t260 - t233 * t265;
t219 = t232 * t265 - t233 * t260;
t217 = (-t229 * t260 + t231 * t265) * qJD(1);
t216 = (t229 * t265 + t231 * t260) * qJD(1);
t213 = -pkin(2) * t217 - qJD(1) * pkin(15);
t212 = qJD(4) + t215;
t211 = (cos(pkin(21)) * t232 - t233 * sin(pkin(21))) * pkin(4) + t241;
t208 = pkin(9) * t215 - pkin(11) * t214 + t226;
t206 = -t272 * t277 - t273 * t275;
t203 = t206 * t263 + t208 * t258;
t202 = -t206 * t258 + t208 * t263;
t1 = Ifges(9,2) * t217 ^ 2 / 0.2e1 + t278 * t283 + m(8) * (t237 + (t219 ^ 2 + t220 ^ 2) * t279) / 0.2e1 + m(6) * (t202 ^ 2 + t203 ^ 2 + t284) / 0.2e1 + Ifges(4,2) * t228 ^ 2 / 0.2e1 + m(4) * (t237 + t278) / 0.2e1 + m(10) * (t213 ^ 2 + t222 ^ 2 + t223 ^ 2) / 0.2e1 + m(5) * (t206 ^ 2 + t226 ^ 2 + t284) / 0.2e1 + (Ifges(4,4) * t228 + Ifges(4,1) * t227 / 0.2e1) * t227 + (Ifges(9,4) * t217 + Ifges(9,1) * t216 / 0.2e1) * t216 + (m(7) * pkin(14) ^ 2 / 0.2e1 + Ifges(2,3) / 0.2e1 + t211 ^ 2 * t283 + (m(3) / 0.2e1 + m(9) / 0.2e1) * pkin(15) ^ 2) * t269 - (-t226 * mrSges(5,1) + t206 * mrSges(5,3) - Ifges(5,2) * t215 / 0.2e1) * t215 + (t202 * mrSges(6,1) - t203 * mrSges(6,2) + Ifges(6,3) * t212 / 0.2e1) * t212 + (t204 * mrSges(5,3) + t226 * mrSges(5,2) - Ifges(5,4) * t215 + Ifges(5,1) * t282 + (Ifges(6,1) * t263 * t282 + t204 * mrSges(6,2) - t202 * mrSges(6,3) + Ifges(6,5) * t212) * t263 + (t204 * mrSges(6,1) - t203 * mrSges(6,3) - Ifges(6,6) * t212 + (-Ifges(6,4) * t263 + Ifges(6,2) * t258 / 0.2e1) * t214) * t258) * t214 + (Ifges(4,5) * t227 + Ifges(9,5) * t216 + Ifges(4,6) * t228 + Ifges(9,6) * t217 + (Ifges(11,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(9,3) / 0.2e1) * t250) * t250 + (t222 * mrSges(10,1) - t223 * mrSges(10,2) + (Ifges(3,3) / 0.2e1 + Ifges(10,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * qJD(2) + ((-t259 * t227 - t264 * t228) * mrSges(4,3) + ((mrSges(11,2) + mrSges(4,2)) * t264 + (mrSges(11,1) + mrSges(4,1)) * t259) * t250) * pkin(1)) * qJD(2) + (t241 * (-mrSges(4,1) * t228 + mrSges(4,2) * t227) - pkin(15) * (-mrSges(9,1) * t217 + mrSges(9,2) * t216) + ((-t211 * mrSges(11,1) + Ifges(11,2) * t247 / 0.2e1) * t247 + (t241 * mrSges(8,2) + Ifges(8,1) * t233 / 0.2e1) * t233 + (pkin(14) * mrSges(7,2) + Ifges(7,1) * t225 / 0.2e1) * t225 + (t211 * mrSges(11,2) + Ifges(11,4) * t247 + Ifges(11,1) * t246 / 0.2e1) * t246 + (t241 * mrSges(8,1) - Ifges(8,4) * t233 + Ifges(8,2) * t232 / 0.2e1) * t232 + (pkin(14) * mrSges(7,1) - Ifges(7,4) * t225 + Ifges(7,2) * t224 / 0.2e1) * t224) * qJD(1) + (t213 * mrSges(10,2) - t222 * mrSges(10,3) + (-pkin(15) * mrSges(3,2) + (Ifges(10,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t265) * qJD(1)) * t265 + (t213 * mrSges(10,1) - t223 * mrSges(10,3) + (-pkin(15) * mrSges(3,1) + (Ifges(10,2) / 0.2e1 + Ifges(3,2) / 0.2e1) * t260 + (-Ifges(3,4) - Ifges(10,4)) * t265) * qJD(1)) * t260 + t250 * (Ifges(11,5) * t246 + Ifges(11,6) * t247) + (Ifges(7,5) * t225 - Ifges(7,6) * t224 + (Ifges(3,5) + Ifges(10,5)) * t265 + (-Ifges(3,6) - Ifges(10,6)) * t260 + ((-t219 * t232 + t220 * t233) * mrSges(8,3) + (-t246 * t259 - t247 * t264) * mrSges(11,3)) * pkin(1)) * qJD(2)) * qJD(1);
T = t1;
