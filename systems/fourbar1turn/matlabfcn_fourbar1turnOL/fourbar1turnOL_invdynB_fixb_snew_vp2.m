% Calculate vector of inverse dynamics base forces with Newton-Euler for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = fourbar1turnOL_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnOL_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_invdynB_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnOL_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnOL_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:52
% EndTime: 2020-04-12 19:40:53
% DurationCPUTime: 0.53s
% Computational Cost: add. (2159->207), mult. (4521->280), div. (0->0), fcn. (2753->8), ass. (0->80)
t290 = sin(qJ(1));
t294 = cos(qJ(1));
t283 = -t294 * g(1) - t290 * g(2);
t289 = sin(qJ(2));
t293 = cos(qJ(2));
t257 = -t293 * g(3) - t289 * t283;
t272 = (-mrSges(3,1) * t293 + mrSges(3,2) * t289) * qJD(1);
t302 = qJD(1) * qJD(2);
t275 = t289 * qJDD(1) + t293 * t302;
t303 = qJD(1) * t293;
t281 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t303;
t295 = qJD(1) ^ 2;
t251 = (t289 * t293 * t295 + qJDD(2)) * pkin(2) + t257;
t258 = -t289 * g(3) + t293 * t283;
t252 = (-t293 ^ 2 * t295 - qJD(2) ^ 2) * pkin(2) + t258;
t288 = sin(qJ(3));
t292 = cos(qJ(3));
t242 = -t292 * t251 + t288 * t252;
t265 = (t288 * t289 - t292 * t293) * qJD(1);
t300 = t289 * t302;
t277 = t293 * qJDD(1) - t300;
t245 = t265 * qJD(3) - t292 * t275 - t288 * t277;
t266 = (-t288 * t293 - t289 * t292) * qJD(1);
t249 = -t265 * mrSges(4,1) + t266 * mrSges(4,2);
t286 = qJD(2) + qJD(3);
t253 = -t286 * mrSges(4,2) + t265 * mrSges(4,3);
t285 = qJDD(2) + qJDD(3);
t234 = m(4) * t242 + t285 * mrSges(4,1) - t245 * mrSges(4,3) - t266 * t249 + t286 * t253;
t243 = -t288 * t251 - t292 * t252;
t244 = -t266 * qJD(3) + t288 * t275 - t292 * t277;
t254 = t286 * mrSges(4,1) - t266 * mrSges(4,3);
t235 = m(4) * t243 - t285 * mrSges(4,2) + t244 * mrSges(4,3) + t265 * t249 - t286 * t254;
t298 = -t292 * t234 - t288 * t235;
t305 = qJD(1) * t289;
t228 = m(3) * t257 + qJDD(2) * mrSges(3,1) - t275 * mrSges(3,3) + qJD(2) * t281 - t272 * t305 + t298;
t279 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t305;
t229 = m(3) * t258 - qJDD(2) * mrSges(3,2) + t277 * mrSges(3,3) - qJD(2) * t279 + t288 * t234 - t292 * t235 + t272 * t303;
t273 = -t295 * pkin(1) + t283;
t287 = sin(qJ(4));
t291 = cos(qJ(4));
t255 = -t291 * g(3) - t287 * t273;
t271 = (-mrSges(5,1) * t291 + mrSges(5,2) * t287) * qJD(1);
t301 = qJD(1) * qJD(4);
t274 = t287 * qJDD(1) + t291 * t301;
t304 = qJD(1) * t291;
t280 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t304;
t306 = qJD(1) * t287;
t240 = m(5) * t255 + qJDD(4) * mrSges(5,1) - t274 * mrSges(5,3) + qJD(4) * t280 - t271 * t306;
t256 = -t287 * g(3) + t291 * t273;
t276 = t291 * qJDD(1) - t287 * t301;
t278 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t306;
t241 = m(5) * t256 - qJDD(4) * mrSges(5,2) + t276 * mrSges(5,3) - qJD(4) * t278 + t271 * t304;
t225 = m(2) * t283 - t295 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t289 * t228 + t293 * t229 - t287 * t240 + t291 * t241;
t282 = t290 * g(1) - t294 * g(2);
t250 = (-t277 + t300) * pkin(2) - t282;
t296 = m(4) * t250 - t244 * mrSges(4,1) + t245 * mrSges(4,2) - t265 * t253 + t266 * t254;
t270 = -qJDD(1) * pkin(1) - t282;
t297 = -m(5) * t270 + t276 * mrSges(5,1) - t274 * mrSges(5,2) + t280 * t304;
t231 = qJDD(1) * mrSges(2,1) + t277 * mrSges(3,1) - t295 * mrSges(2,2) - t275 * mrSges(3,2) + (m(2) + m(3)) * t282 + (-t287 * t278 - t289 * t279 + t293 * t281) * qJD(1) - t296 + t297;
t308 = t290 * t225 + t294 * t231;
t307 = t291 * t240 + t287 * t241;
t299 = t294 * t225 - t290 * t231;
t264 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t289 + Ifges(3,4) * t293) * qJD(1);
t263 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t287 + Ifges(5,4) * t291) * qJD(1);
t262 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t289 + Ifges(3,2) * t293) * qJD(1);
t261 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t287 + Ifges(5,2) * t291) * qJD(1);
t260 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t289 + Ifges(3,6) * t293) * qJD(1);
t259 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t287 + Ifges(5,6) * t291) * qJD(1);
t248 = Ifges(4,1) * t266 + Ifges(4,4) * t265 + Ifges(4,5) * t286;
t247 = Ifges(4,4) * t266 + Ifges(4,2) * t265 + Ifges(4,6) * t286;
t246 = Ifges(4,5) * t266 + Ifges(4,6) * t265 + Ifges(4,3) * t286;
t237 = mrSges(5,2) * t270 - mrSges(5,3) * t255 + Ifges(5,1) * t274 + Ifges(5,4) * t276 + Ifges(5,5) * qJDD(4) - qJD(4) * t261 + t259 * t304;
t236 = -mrSges(5,1) * t270 + mrSges(5,3) * t256 + Ifges(5,4) * t274 + Ifges(5,2) * t276 + Ifges(5,6) * qJDD(4) + qJD(4) * t263 - t259 * t306;
t233 = mrSges(4,2) * t250 - mrSges(4,3) * t242 + Ifges(4,1) * t245 + Ifges(4,4) * t244 + Ifges(4,5) * t285 + t265 * t246 - t286 * t247;
t232 = -mrSges(4,1) * t250 + mrSges(4,3) * t243 + Ifges(4,4) * t245 + Ifges(4,2) * t244 + Ifges(4,6) * t285 - t266 * t246 + t286 * t248;
t227 = -mrSges(3,2) * t282 - mrSges(3,3) * t257 + Ifges(3,1) * t275 + Ifges(3,4) * t277 + Ifges(3,5) * qJDD(2) - qJD(2) * t262 + t288 * t232 - t292 * t233 + t260 * t303;
t226 = mrSges(3,1) * t282 + mrSges(3,3) * t258 + Ifges(3,4) * t275 + Ifges(3,2) * t277 + Ifges(3,6) * qJDD(2) - pkin(2) * t296 + qJD(2) * t264 - t292 * t232 - t288 * t233 - t260 * t305;
t222 = t295 * Ifges(2,5) - pkin(2) * t298 + mrSges(2,3) * t283 - Ifges(4,3) * t285 - Ifges(5,5) * t274 - Ifges(3,5) * t275 - Ifges(5,6) * t276 - Ifges(3,6) * t277 - mrSges(5,1) * t255 + mrSges(5,2) * t256 - mrSges(3,1) * t257 + mrSges(3,2) * t258 + t265 * t248 - t266 * t247 - Ifges(4,6) * t244 - Ifges(4,5) * t245 - pkin(1) * t307 - mrSges(4,1) * t242 + mrSges(4,2) * t243 + Ifges(2,6) * qJDD(1) - Ifges(3,3) * qJDD(2) - Ifges(5,3) * qJDD(4) + mrSges(2,1) * g(3) + (-t287 * t261 - t289 * t262 + t291 * t263 + t293 * t264) * qJD(1);
t221 = -mrSges(2,2) * g(3) - mrSges(2,3) * t282 + Ifges(2,5) * qJDD(1) - t295 * Ifges(2,6) - t289 * t226 + t293 * t227 - t287 * t236 + t291 * t237;
t1 = [-m(1) * g(1) + t299; -m(1) * g(2) + t308; t293 * t228 + t289 * t229 + (-m(1) - m(2)) * g(3) + t307; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t308 + t294 * t221 - t290 * t222; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t299 + t290 * t221 + t294 * t222; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t282 - mrSges(2,2) * t283 + t289 * t227 + t293 * t226 + t287 * t237 + t291 * t236 + pkin(1) * (-t278 * t306 + t297);];
tauB = t1;
