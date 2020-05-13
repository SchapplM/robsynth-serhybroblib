% Calculate Gravitation load with newton euler on the joints for
% palh1m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:49
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m2IC_gravloadJ_floatb_twist_snew_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2IC_gravloadJ_floatb_twist_snew_vp2: qJ has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2IC_gravloadJ_floatb_twist_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2IC_gravloadJ_floatb_twist_snew_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2IC_gravloadJ_floatb_twist_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2IC_gravloadJ_floatb_twist_snew_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:48:53
% EndTime: 2020-05-02 23:48:55
% DurationCPUTime: 1.52s
% Computational Cost: add. (1128->151), mult. (1292->165), div. (19->10), fcn. (1079->42), ass. (0->83)
t1394 = sin(qJ(5));
t1454 = cos(qJ(5));
t1363 = -pkin(9) * m(6) - mrSges(6,1) * t1454 + mrSges(6,2) * t1394 - mrSges(5,1);
t1381 = pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3);
t1395 = sin(qJ(4));
t1403 = cos(qJ(4));
t1349 = t1363 * t1403 - t1381 * t1395;
t1470 = m(5) + m(6);
t1343 = t1470 * pkin(5) + mrSges(4,1) - t1349;
t1396 = sin(qJ(3));
t1404 = cos(qJ(3));
t1348 = t1395 * t1363 + t1381 * t1403;
t1469 = -mrSges(4,2) + t1348;
t1336 = t1343 * t1396 - t1469 * t1404;
t1390 = sin(qJ(9));
t1399 = cos(qJ(9));
t1364 = m(10) * pkin(2) - mrSges(10,1) * t1399 + t1390 * mrSges(10,2) + mrSges(9,1);
t1391 = sin(qJ(8));
t1400 = cos(qJ(8));
t1464 = t1390 * mrSges(10,1) + mrSges(10,2) * t1399 - mrSges(9,2);
t1344 = t1364 * t1400 + t1391 * t1464;
t1392 = sin(qJ(7));
t1401 = cos(qJ(7));
t1386 = sin(pkin(19));
t1387 = cos(pkin(19));
t1367 = mrSges(11,1) * t1386 - mrSges(11,2) * t1387;
t1368 = mrSges(11,1) * t1387 + mrSges(11,2) * t1386;
t1388 = sin(qJ(10));
t1389 = cos(qJ(10));
t1455 = pkin(4) * m(11);
t1461 = t1367 * t1389 - t1368 * t1388 - t1386 * t1455 + mrSges(8,2);
t1463 = -t1367 * t1388 - t1368 * t1389 + t1387 * t1455 + mrSges(8,1);
t1467 = -m(8) - m(11) - m(4) - t1470;
t1326 = t1467 * pkin(1) + t1461 * t1392 - t1463 * t1401 - mrSges(3,1) - t1336 - t1344;
t1335 = t1343 * t1404 + t1396 * t1469;
t1460 = t1391 * t1364 - t1400 * t1464;
t1327 = -t1392 * t1463 - t1401 * t1461 - mrSges(3,2) + t1335 - t1460;
t1397 = sin(qJ(2));
t1405 = cos(qJ(2));
t1475 = t1397 * t1326 + t1327 * t1405;
t1393 = sin(qJ(6));
t1402 = cos(qJ(6));
t1474 = -pkin(14) * m(7) + t1402 * mrSges(7,1) - t1393 * mrSges(7,2) + (m(9) + m(10) + m(3) - t1467) * pkin(15) + mrSges(2,1) + t1475;
t1398 = sin(qJ(1));
t1406 = cos(qJ(1));
t1427 = g(1) * t1406 + t1398 * g(2);
t1472 = -t1348 * t1396 + t1349 * t1404;
t1440 = qJ(4) + pkin(18);
t1441 = pkin(19) + qJ(3);
t1428 = t1440 + t1441;
t1449 = qJ(6) - qJ(2);
t1422 = t1428 + t1449;
t1416 = -0.2e1 * qJ(7) - pkin(20) + t1422;
t1423 = t1428 - t1449;
t1418 = pkin(20) + t1423;
t1465 = cos(qJ(10) - t1416) + cos(qJ(10) - t1418);
t1448 = qJ(10) + qJ(7);
t1337 = t1348 * t1404 + t1396 * t1349;
t1447 = (-(t1397 * t1337 - t1472 * t1405) * g(3) - t1427 * (t1337 * t1405 + t1472 * t1397)) / pkin(10);
t1351 = -t1405 * g(3) + t1427 * t1397;
t1352 = t1397 * g(3) + t1427 * t1405;
t1339 = -t1392 * t1351 + t1401 * t1352;
t1341 = t1401 * t1351 + t1392 * t1352;
t1360 = t1386 * t1389 - t1387 * t1388;
t1361 = -t1386 * t1388 - t1387 * t1389;
t1331 = t1361 * t1339 - t1360 * t1341;
t1332 = t1360 * t1339 + t1361 * t1341;
t1324 = mrSges(11,1) * t1331 - mrSges(11,2) * t1332;
t1446 = t1324 / pkin(8);
t1439 = qJ(7) + pkin(20);
t1434 = -qJ(8) + pkin(17) + qJ(3);
t1429 = t1439 - t1449;
t1420 = -qJ(7) + t1423;
t1419 = -qJ(7) + t1422;
t1413 = 0.1e1 / pkin(3);
t1376 = cos(t1429);
t1375 = cos(qJ(9) - t1434);
t1370 = -g(1) * t1398 + g(2) * t1406;
t1350 = t1394 * mrSges(6,1) + mrSges(6,2) * t1454 - mrSges(2,2) + mrSges(11,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3) + mrSges(10,3);
t1340 = t1400 * t1351 + t1391 * t1352;
t1338 = -t1391 * t1351 + t1400 * t1352;
t1333 = t1403 * (t1396 * t1351 - t1404 * t1352) + t1395 * (t1404 * t1351 + t1396 * t1352);
t1 = [-(t1398 * t1350 + t1474 * t1406) * g(2) - (t1350 * t1406 - t1474 * t1398) * g(1); -t1475 * g(3) + (-(pkin(1) * (sin(qJ(10) - t1449) + sin(qJ(10) + t1449)) + (-sin(-qJ(10) + t1429) + sin(qJ(10) + t1429)) * pkin(3)) * pkin(4) * t1447 + ((-(cos(t1420) + cos(t1419)) * pkin(1) + (-cos(t1416) - cos(t1418)) * pkin(3)) * pkin(4) + ((cos(-qJ(10) + t1420) + cos(-qJ(10) + t1419)) * pkin(1) + t1465 * pkin(3)) * pkin(8)) * t1446) / t1465 * t1413 + (pkin(1) * sin(t1439) / pkin(7) * (mrSges(7,1) * (-t1402 * g(3) + t1393 * t1427) - mrSges(7,2) * (-t1393 * g(3) - t1402 * t1427)) + (-pkin(3) * t1376 - pkin(1) * cos(-t1449)) * t1413 * (mrSges(8,1) * t1339 - mrSges(8,2) * t1341 + (t1387 * (t1331 * t1361 + t1332 * t1360) + t1386 * (-t1331 * t1360 + t1332 * t1361)) * t1455 + t1324)) / t1376 - t1427 * (t1326 * t1405 - t1397 * t1327); -(t1335 * t1405 - t1336 * t1397) * g(3) + (t1375 * (-(-t1344 * t1397 - t1405 * t1460) * g(3) + t1427 * (t1344 * t1405 - t1460 * t1397)) + (-pkin(12) * t1375 + pkin(2) * cos(t1434)) / pkin(12) * (mrSges(10,1) * (-t1399 * t1338 + t1390 * t1340) - mrSges(10,2) * (-t1390 * t1338 - t1399 * t1340))) * pkin(6) / pkin(2) / t1390 - pkin(10) * t1447 + t1427 * (t1335 * t1397 + t1405 * t1336) + (-cos(t1441 - t1448) * t1447 + sin(t1440) * t1446) / cos(t1428 - t1448) * pkin(5); mrSges(6,1) * (-t1394 * t1333 + t1454 * t1370) - mrSges(6,2) * (t1454 * t1333 + t1394 * t1370);];
taug = t1(:);
