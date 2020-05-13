% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh3m1TE
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in palh3m1TE_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = palh3m1TE_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1TE_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1TE_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_jacobiaD_rot_sym_varpar: pkin has to be [19x1] (double)');
JaD_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:03
	% EndTime: 2020-04-18 09:52:03
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:53:57
	% EndTime: 2020-04-18 09:53:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:55:55
	% EndTime: 2020-04-18 10:06:46
	% DurationCPUTime: 578.57s
	% Computational Cost: add. (12986363->354), mult. (20105870->650), div. (805951->22), fcn. (12596361->21), ass. (0->334)
	t1354 = sin(qJ(2));
	t1356 = sin(pkin(16));
	t1358 = cos(qJ(2));
	t1360 = cos(pkin(16));
	t1301 = t1354 * t1356 - t1358 * t1360;
	t1286 = pkin(5) * t1301;
	t1280 = (-0.2e1 * t1286 + pkin(1)) * pkin(1);
	t1364 = -pkin(6) - pkin(2);
	t1258 = (pkin(5) - t1364) * (pkin(5) + t1364) + t1280;
	t1363 = -pkin(6) + pkin(2);
	t1259 = (pkin(5) - t1363) * (pkin(5) + t1363) + t1280;
	t1249 = t1259 * t1258;
	t1075 = sqrt(-t1249);
	t1378 = pkin(5) ^ 2;
	t1273 = t1280 + t1378;
	t1264 = pkin(2) ^ 2 - pkin(6) ^ 2 + t1273;
	t1261 = pkin(5) * t1264;
	t1281 = t1354 * t1360 + t1358 * t1356;
	t1253 = t1281 * t1261;
	t1283 = -t1286 + pkin(1);
	t1381 = 0.1e1 / pkin(2);
	t1248 = t1381 * (t1283 * t1075 + t1253);
	t1268 = 0.1e1 / t1273;
	t1423 = t1268 / 0.2e1;
	t1237 = t1248 * t1423;
	t1279 = pkin(5) * t1281;
	t1270 = t1075 * t1279;
	t1247 = t1381 * (t1283 * t1264 - t1270);
	t1243 = t1268 * t1247;
	t1353 = sin(qJ(3));
	t1357 = cos(qJ(3));
	t1437 = -t1357 / 0.2e1;
	t1216 = t1353 * t1237 + t1243 * t1437;
	t1422 = t1353 / 0.2e1;
	t1217 = t1357 * t1237 + t1243 * t1422;
	t1323 = pkin(18) + pkin(19);
	t1312 = sin(t1323);
	t1313 = cos(t1323);
	t1201 = t1313 * t1216 + t1312 * t1217;
	t1198 = pkin(3) * t1201;
	t1193 = (0.2e1 * t1198 + pkin(4)) * pkin(4);
	t1362 = -pkin(8) - pkin(10);
	t1188 = (pkin(3) - t1362) * (pkin(3) + t1362) + t1193;
	t1361 = -pkin(8) + pkin(10);
	t1189 = (pkin(3) - t1361) * (pkin(3) + t1361) + t1193;
	t1177 = t1189 * t1188;
	t1074 = sqrt(-t1177);
	t1275 = pkin(1) * t1279;
	t1244 = -0.2e1 * t1283 * t1275 - t1253;
	t1351 = pkin(5) * t1075;
	t1271 = t1301 * t1351;
	t1366 = pkin(1) * pkin(5);
	t1442 = (t1258 + t1259) * t1366;
	t1037 = 0.2e1 * t1442 * t1281;
	t1370 = 0.1e1 / t1075;
	t1230 = t1370 * t1037;
	t1424 = -t1230 / 0.2e1;
	t1208 = t1381 * (t1279 * t1424 + t1244 + t1271);
	t1204 = t1208 * t1423;
	t1254 = t1301 * t1261;
	t1278 = t1378 * t1281 ^ 2;
	t1375 = -0.2e1 * pkin(1);
	t1445 = t1278 * t1375 - t1254 - t1270;
	t1209 = t1381 * (t1283 * t1230 / 0.2e1 + t1445);
	t1207 = t1268 * t1209;
	t1367 = 0.1e1 / t1273 ^ 2;
	t1411 = t1366 * t1367;
	t1245 = t1247 * t1411;
	t1239 = t1281 * t1245;
	t1246 = t1248 * t1411;
	t1241 = t1281 * t1246;
	t1214 = -t1357 * t1239 + t1353 * t1241;
	t1180 = -t1353 * t1207 / 0.2e1 + t1357 * t1204 - t1214;
	t1215 = t1353 * t1239 + t1357 * t1241;
	t1421 = t1357 / 0.2e1;
	t1181 = t1353 * t1204 + t1207 * t1421 + t1215;
	t1203 = -t1312 * t1216 + t1313 * t1217;
	t1152 = (t1313 * t1180 - t1312 * t1181) * qJD(2) - t1203 * qJD(3);
	t1150 = pkin(3) * t1152;
	t1450 = t1074 * t1150;
	t1446 = pkin(3) * t1203;
	t1449 = t1074 * t1446;
	t1374 = -0.2e1 * pkin(4);
	t1196 = t1201 * t1374;
	t1377 = pkin(8) ^ 2;
	t1379 = pkin(4) ^ 2;
	t1380 = pkin(3) ^ 2;
	t1191 = -pkin(3) * t1196 + pkin(10) ^ 2 - t1377 + t1379 + t1380;
	t1448 = t1191 * t1150;
	t1035 = t1379 + (-t1196 + pkin(3)) * pkin(3);
	t1034 = 0.1e1 / t1035 ^ 2;
	t1372 = 0.1e1 / t1035;
	t1447 = 0.4e1 * t1372 * t1034 * t1379 * t1380;
	t1420 = -t1370 / 0.2e1;
	t1265 = t1279 * t1420;
	t1444 = t1037 * t1265 + t1271;
	t1418 = t1372 / 0.2e1;
	t1146 = pkin(4) * t1152;
	t1144 = pkin(3) * t1146;
	t1402 = t1188 + t1189;
	t1117 = t1402 * t1144;
	t1116 = 0.2e1 * t1117;
	t1371 = 0.1e1 / t1074;
	t1115 = t1371 * t1116;
	t1427 = -pkin(3) / 0.2e1;
	t1114 = t1115 * t1427;
	t1440 = (-t1312 * t1180 - t1313 * t1181) * qJD(2) + t1201 * qJD(3);
	t1149 = pkin(3) * t1440;
	t1376 = 0.1e1 / pkin(10);
	t1194 = t1198 + pkin(4);
	t1430 = -0.2e1 * t1194;
	t1102 = t1376 * (-t1074 * t1149 - t1203 * t1114 + t1144 * t1430 - t1448);
	t1145 = t1376 * t1152;
	t1167 = t1194 * t1191 + t1449;
	t1133 = t1167 * t1145;
	t1365 = pkin(3) * pkin(4);
	t1415 = t1034 * t1365;
	t1443 = t1102 * t1415 + t1133 * t1447;
	t1050 = t1354 * t1353 - t1358 * t1357;
	t1441 = t1358 * t1353 + t1354 * t1357;
	t1428 = -0.2e1 * t1380;
	t1143 = t1146 * t1428;
	t1425 = t1194 / 0.2e1;
	t1103 = t1376 * (t1115 * t1425 - t1143 * t1203 + t1191 * t1149 - t1450);
	t1190 = pkin(3) * t1191;
	t1168 = t1194 * t1074 - t1203 * t1190;
	t1134 = t1168 * t1145;
	t1321 = t1376 * t1418;
	t1342 = t1380 * t1374;
	t1439 = t1103 * t1415 + t1134 * t1447 + (t1342 * t1440 + t1114) * t1321;
	t1382 = pkin(1) ^ 2;
	t1438 = 0.4e1 * t1268 * t1367 * t1382 * t1278;
	t1436 = t1371 / 0.2e1;
	t1166 = t1376 * t1168;
	t1162 = t1166 * t1418;
	t1165 = t1376 * t1167;
	t1163 = t1372 * t1165;
	t1348 = sin(pkin(17));
	t1349 = cos(pkin(17));
	t1025 = -t1349 * t1163 / 0.2e1 + t1348 * t1162;
	t1026 = t1349 * t1162 + t1348 * t1163 / 0.2e1;
	t1359 = cos(qJ(1));
	t1048 = t1050 * t1359;
	t1049 = t1441 * t1359;
	t1011 = t1025 * t1048 + t1026 * t1049;
	t1072 = sin(qJ(4));
	t1073 = cos(qJ(4));
	t1355 = sin(qJ(1));
	t1314 = t1355 * t1073;
	t1288 = -t1011 * t1072 + t1314;
	t1433 = t1288 * qJD(4);
	t1410 = qJD(2) + qJD(3);
	t1227 = t1381 * (t1244 + t1444);
	t1221 = t1268 * t1227;
	t1419 = t1370 / 0.2e1;
	t1276 = t1283 * t1419;
	t1231 = t1381 * (t1037 * t1276 + t1445);
	t1226 = t1268 * t1231;
	t1184 = t1221 * t1422 + t1226 * t1421 + t1215;
	t1431 = -t1221 * t1437 - t1226 * t1422 - t1214;
	t1315 = t1355 * t1072;
	t1005 = t1011 * t1073 + t1315;
	t1000 = 0.1e1 / t1005 ^ 2;
	t1317 = t1359 * t1072;
	t1047 = t1050 * t1355;
	t1038 = t1047 * qJD(1) - t1410 * t1049;
	t1046 = t1441 * t1355;
	t1039 = -t1046 * qJD(1) - t1410 * t1048;
	t1397 = t1102 * t1418 + t1133 * t1415;
	t1409 = t1103 * t1418 + t1134 * t1415;
	t978 = t1409 * t1348 - t1397 * t1349;
	t979 = t1397 * t1348 + t1409 * t1349;
	t954 = -t1025 * t1038 + t1026 * t1039 + t1048 * t978 + t1049 * t979;
	t953 = qJD(1) * t1317 + t954 * t1073 + t1433;
	t999 = 0.1e1 / t1005;
	t1336 = t999 * t1000 * t953;
	t1429 = -0.2e1 * t1336;
	t1373 = 0.4e1 * pkin(4);
	t1164 = -((pkin(3) + pkin(10)) * (pkin(3) - pkin(10)) + t1193) * t1446 * t1373 + 0.4e1 * t1203 * t1377 * t1365;
	t1161 = t1371 * t1164;
	t1426 = -t1161 / 0.2e1;
	t1024 = t1048 * t1026;
	t1012 = t1025 * t1049 - t1024;
	t1017 = -t1025 * t1050 - t1026 * t1441;
	t1302 = t1025 * t1046 - t1026 * t1047;
	t991 = atan2(t1302, t1017);
	t986 = sin(t991);
	t987 = cos(t991);
	t971 = t1017 * t987 + t1302 * t986;
	t969 = 0.1e1 / t971 ^ 2;
	t1335 = t1012 * t969;
	t1015 = 0.1e1 / t1017 ^ 2;
	t1043 = t1410 * t1050;
	t1044 = t1410 * t1441;
	t958 = -t1025 * t1044 + t1026 * t1043 - t1050 * t978 - t1441 * t979;
	t1417 = t1015 * t958;
	t1331 = t1038 * t1026 - t1048 * t979;
	t1416 = t1034 * t1144;
	t1414 = 0.1e1 / t1249 * t1230;
	t1412 = t1275 * t1367;
	t1176 = -t1312 * t1184 + t1313 * t1431;
	t1172 = pkin(4) * t1176;
	t1171 = pkin(3) * t1172;
	t1138 = t1402 * t1171;
	t1137 = 0.2e1 * t1371 * t1138;
	t1175 = -t1313 * t1184 - t1312 * t1431;
	t1173 = pkin(3) * t1175;
	t1174 = pkin(3) * t1176;
	t1110 = t1376 * (-t1172 * t1203 * t1428 - t1074 * t1174 + t1137 * t1425 + t1191 * t1173);
	t1160 = t1166 * t1415;
	t1408 = t1110 * t1418 + t1176 * t1160;
	t1195 = t1203 * t1342;
	t1124 = t1376 * (t1161 * t1425 + t1201 * t1190 + t1203 * t1195 + t1449);
	t1406 = t1124 * t1418 - t1203 * t1160;
	t1159 = t1165 * t1415;
	t1405 = -t1349 * t1159 + t1348 * t1160;
	t1404 = t1348 * t1159 + t1349 * t1160;
	t1403 = 0.2e1 * t1402 * t1365;
	t1401 = -t1365 * t1430 + t1190;
	t1352 = pkin(3) * t1074;
	t1398 = -t1352 - t1195;
	t1136 = t1137 * t1427;
	t1109 = t1376 * (-t1074 * t1173 - t1136 * t1203 + t1171 * t1430 - t1191 * t1174);
	t1396 = t1109 * t1418 + t1176 * t1159;
	t1123 = t1376 * (-t1201 * t1352 + t1401 * t1203 - t1426 * t1446);
	t1395 = t1123 * t1418 - t1203 * t1159;
	t1386 = -t1301 * t1246 + t1248 * t1438 + (t1209 + t1231) * t1412;
	t1387 = -t1301 * t1245 + t1247 * t1438 + (t1208 + t1227) * t1412;
	t1169 = t1387 * t1353 + t1386 * t1357;
	t1170 = t1386 * t1353 - t1387 * t1357;
	t1269 = qJD(2) * t1278;
	t1277 = qJD(2) * t1301;
	t1036 = -0.8e1 * t1382 * t1269 - 0.2e1 * t1442 * t1277;
	t1229 = qJD(2) * t1037;
	t1228 = pkin(5) * t1229;
	t1225 = t1281 * t1228;
	t1032 = t1225 * t1420 + t1283 * t1229 * t1414 / 0.4e1 + t1036 * t1276 + (0.6e1 * t1301 * pkin(1) * t1378 * t1281 - t1253 + t1444) * qJD(2);
	t1033 = t1301 * t1228 * t1419 - t1225 * t1414 / 0.4e1 + t1036 * t1265 + 0.4e1 * pkin(1) * t1269 - (t1283 * t1375 + t1424) * pkin(5) * t1277 + (t1281 * t1351 + t1254) * qJD(2);
	t1266 = t1381 * t1268;
	t1260 = t1266 * t1422;
	t1262 = t1357 * t1266;
	t1182 = -t1033 * t1262 / 0.2e1 + t1032 * t1260 + t1184 * qJD(3);
	t1183 = t1032 * t1262 / 0.2e1 + t1033 * t1260 + t1431 * qJD(3);
	t1018 = -t1313 * t1182 - t1312 * t1183 + (-t1312 * t1169 - t1313 * t1170) * qJD(2);
	t1147 = t1380 * t1152;
	t1135 = t1176 * t1147;
	t1368 = t1371 / t1177;
	t1393 = t1368 * t1117 * t1138 + (t1403 * t1018 - 0.8e1 * t1379 * t1135) * t1436;
	t1139 = t1203 * t1147;
	t1392 = t1368 * t1164 * t1116 / 0.4e1 + (0.8e1 * t1379 * t1139 + t1403 * t1440) * t1436;
	t1125 = pkin(3) * (t1312 * t1182 - t1313 * t1183 + (-t1313 * t1169 + t1312 * t1170) * qJD(2));
	t1389 = (t1398 * t1018 + t1191 * t1125 + t1152 * t1136 + t1175 * t1143 + t1393 * t1194) * t1321 + t1110 * t1416 + t1439 * t1176;
	t1388 = (t1201 * t1143 + t1150 * t1426 + t1392 * t1194 + t1398 * t1440 - t1448) * t1321 + t1124 * t1416 - t1439 * t1203;
	t1385 = (-t1401 * t1018 - t1074 * t1125 + t1175 * t1114 + t1135 * t1373 + t1136 * t1440 + t1393 * t1446) * t1321 + t1109 * t1416 + t1443 * t1176;
	t1384 = (t1201 * t1114 - t1139 * t1373 + t1149 * t1426 + t1392 * t1446 - t1401 * t1440 + t1450) * t1321 + t1123 * t1416 - t1443 * t1203;
	t968 = 0.1e1 / t971;
	t1014 = 0.1e1 / t1017;
	t1007 = t1012 ^ 2;
	t1303 = -t1017 * t986 + t1302 * t987;
	t1324 = t1302 * t1015;
	t1040 = t1049 * qJD(1) - t1410 * t1047;
	t1041 = -t1048 * qJD(1) - t1410 * t1046;
	t957 = t1025 * t1040 + t1026 * t1041 + t1046 * t978 - t1047 * t979;
	t1292 = t1014 * t957 - t958 * t1324;
	t1006 = t1302 ^ 2;
	t990 = t1006 * t1015 + 0.1e1;
	t988 = 0.1e1 / t990;
	t942 = t1292 * t988;
	t933 = t1303 * t942 + t957 * t986 + t958 * t987;
	t1346 = t933 * t968 * t969;
	t955 = t1025 * t1039 + t1049 * t978 + t1331;
	t962 = t1007 * t969 + 0.1e1;
	t1347 = 0.2e1 * (-t1007 * t1346 + t955 * t1335) / t962 ^ 2;
	t1326 = t1000 * t1288;
	t1316 = t1359 * t1073;
	t952 = -qJD(1) * t1316 + t1005 * qJD(4) + t954 * t1072;
	t998 = t1288 ^ 2;
	t985 = t1000 * t998 + 0.1e1;
	t1345 = 0.2e1 * (-t952 * t1326 - t998 * t1336) / t985 ^ 2;
	t1332 = t1014 * t1417;
	t1344 = 0.2e1 * (-t1006 * t1332 + t957 * t1324) / t990 ^ 2;
	t1341 = t1404 * t1018 + t1385 * t1348 + t1389 * t1349 + t978;
	t944 = t1405 * t1018 + t1389 * t1348 - t1385 * t1349;
	t1340 = t944 - t979;
	t1339 = t1384 * t1348 + t1388 * t1349 + t1404 * t1440 + t978;
	t951 = t1388 * t1348 - t1384 * t1349 + t1405 * t1440;
	t1338 = t951 - t979;
	t983 = 0.1e1 / t985;
	t1337 = t1000 * t983;
	t1334 = t1012 * t986;
	t1333 = t1012 * t987;
	t1330 = t1396 * t1348 + t1408 * t1349 + t1025;
	t1329 = t1395 * t1348 + t1406 * t1349 + t1025;
	t981 = t1408 * t1348 - t1396 * t1349;
	t1328 = t1026 - t981;
	t992 = t1406 * t1348 - t1395 * t1349;
	t1327 = t1026 - t992;
	t1325 = t1302 * t1014;
	t1322 = -0.2e1 * t1346;
	t1320 = t1000 * t1345;
	t1319 = t952 * t1337;
	t1318 = t1014 * t1344;
	t1311 = t1012 * t1322;
	t1310 = t1288 * t1429;
	t1309 = 0.2e1 * t1302 * t1332;
	t963 = -t1328 * t1046 - t1330 * t1047;
	t966 = t1328 * t1050 - t1330 * t1441;
	t1291 = -t1014 * t963 + t966 * t1324;
	t972 = -t1327 * t1046 - t1329 * t1047;
	t975 = t1327 * t1050 - t1329 * t1441;
	t1290 = -t1014 * t972 + t975 * t1324;
	t1008 = -t1025 * t1047 - t1026 * t1046;
	t1003 = t1008 * t1073 + t1317;
	t1289 = -t1008 * t1072 + t1316;
	t1287 = -t1072 * t999 - t1073 * t1326;
	t1285 = t986 + (t987 * t1325 - t986) * t988;
	t1284 = t1287 * t983;
	t1250 = t1287 * t1345 + ((qJD(4) * t999 + t1310) * t1073 + (-t1073 * t952 + (-t953 - t1433) * t1072) * t1000) * t983;
	t974 = t1329 * t1048 + t1327 * t1049;
	t965 = t1330 * t1048 + t1328 * t1049;
	t960 = 0.1e1 / t962;
	t956 = t1025 * t1041 - t1026 * t1040 - t1046 * t979 - t1047 * t978;
	t949 = t1285 * t1012;
	t948 = t1290 * t988;
	t946 = t1291 * t988;
	t941 = t1329 * t1043 + t1327 * t1044 - t1338 * t1050 - t1339 * t1441;
	t940 = -t1327 * t1040 + t1329 * t1041 + t1338 * t1046 - t1339 * t1047;
	t938 = -t1303 * t948 + t972 * t986 + t975 * t987;
	t937 = t1330 * t1043 + t1328 * t1044 - t1340 * t1050 - t1341 * t1441;
	t936 = -t1328 * t1040 + t1330 * t1041 + t1340 * t1046 - t1341 * t1047;
	t934 = -t1303 * t946 + t963 * t986 + t966 * t987;
	t931 = t1290 * t1344 + (t975 * t1309 + t1014 * t940 + (-t1302 * t941 - t957 * t975 - t958 * t972) * t1015) * t988;
	t930 = t1291 * t1344 + (t966 * t1309 + t1014 * t936 + (-t1302 * t937 - t957 * t966 - t958 * t963) * t1015) * t988;
	t1 = [t1014 * t955 * t988 + (-t988 * t1417 - t1318) * t1012, t930, t931, 0; (-t1302 * t968 - t949 * t1335) * t1347 + (t957 * t968 + (-t1302 * t933 + t949 * t955) * t969 + (((-t942 * t988 * t1325 + t1344) * t986 + (-t1302 * t1318 + t942 + (t1292 - t942) * t988) * t987) * t1335 + t1322 * t949 + t1285 * t969 * t955) * t1012) * t960, (-t934 * t1335 - t965 * t968) * t1347 + ((-t1330 * t1038 + t1328 * t1039 + t1341 * t1048 - t1340 * t1049) * t968 + t934 * t1311 + (-t965 * t933 + t934 * t955 + (t1302 * t930 - t946 * t957 + t937 + (t1017 * t946 + t963) * t942) * t1333 + (-t1017 * t930 + t946 * t958 + t936 + (t1302 * t946 - t966) * t942) * t1334) * t969) * t960, (-t938 * t1335 - t968 * t974) * t1347 + ((-t1329 * t1038 + t1327 * t1039 + t1339 * t1048 - t1338 * t1049) * t968 + t938 * t1311 + (-t974 * t933 + t938 * t955 + (t1302 * t931 - t948 * t957 + t941 + (t1017 * t948 + t972) * t942) * t1333 + (-t1017 * t931 + t948 * t958 + t940 + (t1302 * t948 - t975) * t942) * t1334) * t969) * t960, 0; (-t1288 * t1320 - t1319) * t1003 - (-t953 * t1337 - t999 * t1345) * t1289 + ((qJD(1) * t1314 + t1003 * qJD(4) + t956 * t1072) * t999 + (-qJD(1) * t1315 + t1289 * qJD(4) + t956 * t1073) * t1326 + t1003 * t1310) * t983, -(-t1038 * t981 + t1330 * t1039 + t1048 * t944 + t1341 * t1049 + t1331) * t1284 + t1250 * (t1048 * t981 + t1330 * t1049 - t1024), -(-t1038 * t992 + t1329 * t1039 + t1048 * t951 + t1339 * t1049 + t1331) * t1284 + t1250 * (t1048 * t992 + t1329 * t1049 - t1024), -t1345 - (0.2e1 * t1319 - (t983 * t1429 - t1320) * t1288) * t1288;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:07
	% EndTime: 2020-04-18 09:52:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:09
	% EndTime: 2020-04-18 09:52:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:53:32
	% EndTime: 2020-04-18 09:53:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
end