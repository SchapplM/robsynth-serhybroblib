% Calculate kinetic energy for
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m2TE_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_energykin_fixb_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2TE_energykin_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2TE_energykin_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2TE_energykin_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:13:26
% EndTime: 2020-05-01 20:13:27
% DurationCPUTime: 0.98s
% Computational Cost: add. (872->146), mult. (1668->253), div. (0->0), fcn. (1800->22), ass. (0->71)
t277 = cos(pkin(20));
t282 = sin(pkin(18));
t287 = cos(pkin(18));
t302 = sin(pkin(20));
t252 = t282 * t277 - t287 * t302;
t255 = t287 * t277 + t282 * t302;
t272 = pkin(22) + pkin(21);
t267 = sin(t272);
t268 = cos(t272);
t236 = (t252 * t267 + t255 * t268) * qJD(1);
t280 = sin(qJ(3));
t285 = cos(qJ(3));
t239 = t252 * t285 - t280 * t255;
t242 = t280 * t252 + t255 * t285;
t281 = sin(qJ(2));
t286 = cos(qJ(2));
t230 = t239 * t286 - t242 * t281;
t231 = t239 * t281 + t242 * t286;
t293 = -t267 * t230 - t231 * t268;
t273 = qJD(2) + qJD(3);
t304 = pkin(1) * qJD(2);
t294 = t273 * pkin(5) + t280 * t304;
t296 = t230 * t268 - t267 * t231;
t298 = t285 * t304;
t225 = -t293 * t294 + t296 * t298;
t307 = t225 ^ 2;
t306 = m(11) / 0.2e1;
t235 = (-t252 * t268 + t255 * t267) * qJD(1);
t305 = t235 / 0.2e1;
t264 = pkin(1) * t281 - pkin(15);
t303 = cos(pkin(22));
t301 = sin(pkin(22));
t300 = qJD(2) ^ 2 * pkin(1) ^ 2;
t299 = (t280 ^ 2 + t285 ^ 2) * t300;
t249 = (-t280 * t281 + t285 * t286) * qJD(1);
t247 = -t249 * pkin(5) + t264 * qJD(1);
t297 = pkin(2) * t273;
t254 = -t282 * t303 + t287 * t301;
t290 = qJD(1) ^ 2;
t288 = cos(pkin(17));
t284 = cos(qJ(4));
t283 = sin(pkin(17));
t279 = sin(qJ(4));
t278 = cos(pkin(19));
t276 = sin(pkin(19));
t274 = qJ(3) + qJ(2);
t270 = cos(t274);
t269 = sin(t274);
t258 = t264 ^ 2 * t290;
t257 = t283 * t282 + t288 * t287;
t256 = t282 * t288 - t287 * t283;
t253 = -t280 * t276 + t285 * t278;
t251 = -t282 * t301 - t303 * t287;
t250 = t285 * t276 + t280 * t278;
t248 = (t280 * t286 + t281 * t285) * qJD(1);
t246 = t281 * t256 + t257 * t286;
t245 = -t256 * t286 + t281 * t257;
t244 = t253 * t297;
t243 = t250 * t297;
t241 = t281 * t251 - t254 * t286;
t240 = -t251 * t286 - t254 * t281;
t238 = (-t250 * t281 + t253 * t286) * qJD(1);
t237 = (t250 * t286 + t253 * t281) * qJD(1);
t234 = -t238 * pkin(2) - qJD(1) * pkin(15);
t233 = qJD(4) + t236;
t232 = (-cos(pkin(21)) * t251 - t254 * sin(pkin(21))) * pkin(4) + t264;
t229 = t236 * pkin(9) - t235 * pkin(11) + t247;
t227 = -t293 * t298 - t296 * t294;
t224 = t284 * t227 + t279 * t229;
t223 = -t279 * t227 + t284 * t229;
t1 = m(4) * (t258 + t299) / 0.2e1 + m(10) * (t234 ^ 2 + t243 ^ 2 + t244 ^ 2) / 0.2e1 + m(5) * (t227 ^ 2 + t247 ^ 2 + t307) / 0.2e1 + Ifges(9,2) * t238 ^ 2 / 0.2e1 + Ifges(4,2) * t249 ^ 2 / 0.2e1 + m(6) * (t223 ^ 2 + t224 ^ 2 + t307) / 0.2e1 + t299 * t306 + m(8) * (t258 + (t240 ^ 2 + t241 ^ 2) * t300) / 0.2e1 + (Ifges(4,4) * t249 + Ifges(4,1) * t248 / 0.2e1) * t248 + (Ifges(9,4) * t238 + Ifges(9,1) * t237 / 0.2e1) * t237 + (t232 ^ 2 * t306 + Ifges(2,3) / 0.2e1 + m(7) * pkin(14) ^ 2 / 0.2e1 + (m(9) / 0.2e1 + m(3) / 0.2e1) * pkin(15) ^ 2) * t290 - (-t247 * mrSges(5,1) + t227 * mrSges(5,3) - Ifges(5,2) * t236 / 0.2e1) * t236 + (t223 * mrSges(6,1) - t224 * mrSges(6,2) + Ifges(6,3) * t233 / 0.2e1) * t233 + (t247 * mrSges(5,2) - Ifges(5,4) * t236 + t225 * mrSges(5,3) + Ifges(5,1) * t305 + (Ifges(6,1) * t284 * t305 + t225 * mrSges(6,2) - t223 * mrSges(6,3) + Ifges(6,5) * t233) * t284 + (t225 * mrSges(6,1) - t224 * mrSges(6,3) - Ifges(6,6) * t233 + (-Ifges(6,4) * t284 + Ifges(6,2) * t279 / 0.2e1) * t235) * t279) * t235 + (Ifges(4,5) * t248 + Ifges(9,5) * t237 + Ifges(4,6) * t249 + Ifges(9,6) * t238 + (Ifges(11,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(9,3) / 0.2e1) * t273) * t273 + (-t244 * mrSges(10,2) + t243 * mrSges(10,1) + (Ifges(3,3) / 0.2e1 + Ifges(7,3) / 0.2e1 + Ifges(10,3) / 0.2e1) * qJD(2) + ((-t280 * t248 - t285 * t249) * mrSges(4,3) + ((mrSges(11,2) + mrSges(4,2)) * t285 + (mrSges(11,1) + mrSges(4,1)) * t280) * t273) * pkin(1)) * qJD(2) + (t264 * (-t249 * mrSges(4,1) + t248 * mrSges(4,2)) - pkin(15) * (-t238 * mrSges(9,1) + t237 * mrSges(9,2)) + ((-t232 * mrSges(11,1) + Ifges(11,2) * t270 / 0.2e1) * t270 + (t264 * mrSges(8,2) + Ifges(8,1) * t254 / 0.2e1) * t254 + (pkin(14) * mrSges(7,2) + Ifges(7,1) * t246 / 0.2e1) * t246 + (t232 * mrSges(11,2) + Ifges(11,4) * t270 + Ifges(11,1) * t269 / 0.2e1) * t269 + (-t264 * mrSges(8,1) + Ifges(8,4) * t254 + Ifges(8,2) * t251 / 0.2e1) * t251 + (pkin(14) * mrSges(7,1) - Ifges(7,4) * t246 + Ifges(7,2) * t245 / 0.2e1) * t245) * qJD(1) + (t234 * mrSges(10,2) - t243 * mrSges(10,3) + (-pkin(15) * mrSges(3,2) + (Ifges(10,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t286) * qJD(1)) * t286 + (t234 * mrSges(10,1) - t244 * mrSges(10,3) + (-pkin(15) * mrSges(3,1) + (Ifges(10,2) / 0.2e1 + Ifges(3,2) / 0.2e1) * t281 + (-Ifges(3,4) - Ifges(10,4)) * t286) * qJD(1)) * t281 + t273 * (Ifges(11,5) * t269 + Ifges(11,6) * t270) + (Ifges(7,5) * t246 - Ifges(7,6) * t245 + (Ifges(3,5) + Ifges(10,5)) * t286 + (-Ifges(3,6) - Ifges(10,6)) * t281 + ((t240 * t251 + t241 * t254) * mrSges(8,3) + (-t269 * t280 - t270 * t285) * mrSges(11,3)) * pkin(1)) * qJD(2)) * qJD(1);
T = t1;
