% Calculate Gravitation load with newton euler on the joints for
% palh2m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m2IC_gravloadJ_floatb_twist_snew_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2IC_gravloadJ_floatb_twist_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2IC_gravloadJ_floatb_twist_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2IC_gravloadJ_floatb_twist_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2IC_gravloadJ_floatb_twist_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2IC_gravloadJ_floatb_twist_snew_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 06:52:55
% EndTime: 2020-05-03 06:52:55
% DurationCPUTime: 0.46s
% Computational Cost: add. (1060->62), mult. (1214->77), div. (0->0), fcn. (1034->12), ass. (0->45)
t1301 = m(6) + m(7);
t1317 = m(5) + t1301;
t1291 = sin(qJ(4));
t1297 = cos(qJ(4));
t1289 = sin(qJ(6));
t1295 = cos(qJ(6));
t1277 = pkin(3) * m(7) + t1295 * mrSges(7,1) - mrSges(7,2) * t1289 + mrSges(6,1);
t1288 = mrSges(7,3) + mrSges(6,2);
t1290 = sin(qJ(5));
t1296 = cos(qJ(5));
t1324 = -t1290 * t1277 - t1288 * t1296 - mrSges(5,2);
t1328 = pkin(5) * t1301 + t1277 * t1296 - t1288 * t1290 + mrSges(5,1);
t1333 = t1324 * t1291 + t1328 * t1297;
t1260 = pkin(2) * t1317 + mrSges(4,1) + t1333;
t1264 = t1328 * t1291 - t1324 * t1297;
t1262 = mrSges(4,2) + t1264;
t1292 = sin(qJ(3));
t1298 = cos(qJ(3));
t1338 = -t1260 * t1292 - t1262 * t1298;
t1294 = sin(qJ(1));
t1300 = cos(qJ(1));
t1282 = -t1294 * g(1) + t1300 * g(2);
t1336 = -t1292 * t1264 + t1333 * t1298;
t1334 = -mrSges(3,2) + t1338;
t1330 = -t1334 + t1338;
t1254 = t1260 * t1298 - t1262 * t1292;
t1310 = m(4) + t1317;
t1302 = pkin(4) * t1310 + mrSges(3,1) + t1254;
t1329 = t1302 - t1254;
t1299 = cos(qJ(2));
t1326 = mrSges(2,1) + (m(3) + t1310) * pkin(1) + t1302 * t1299;
t1304 = t1300 * g(1) + t1294 * g(2);
t1293 = sin(qJ(2));
t1276 = mrSges(7,1) * t1289 + mrSges(7,2) * t1295 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t1274 = -t1293 * g(3) - t1304 * t1299;
t1273 = -t1299 * g(3) + t1304 * t1293;
t1266 = t1292 * t1273 + t1298 * t1274;
t1265 = t1298 * t1273 - t1292 * t1274;
t1258 = t1291 * t1265 + t1297 * t1266;
t1257 = t1297 * t1265 - t1291 * t1266;
t1256 = t1264 * t1298 + t1292 * t1333;
t1250 = t1290 * t1257 + t1296 * t1258;
t1248 = t1295 * t1250 - t1289 * t1282;
t1247 = -t1289 * t1250 - t1295 * t1282;
t1 = [-(-t1276 * t1300 - t1326 * t1294) * g(1) - (-t1294 * t1276 + t1326 * t1300) * g(2) - t1282 * t1293 * t1334; (t1330 * t1293 - t1329 * t1299) * g(3) + t1304 * (t1329 * t1293 + t1330 * t1299); -(-t1256 * t1293 + t1336 * t1299) * g(3) + mrSges(6,2) * t1250 + (-t1289 * t1247 + t1295 * t1248) * mrSges(7,3) - t1277 * (t1296 * t1257 - t1290 * t1258) + t1304 * (t1256 * t1299 + t1336 * t1293); mrSges(7,1) * t1247 - mrSges(7,2) * t1248;];
taug = t1(:);
