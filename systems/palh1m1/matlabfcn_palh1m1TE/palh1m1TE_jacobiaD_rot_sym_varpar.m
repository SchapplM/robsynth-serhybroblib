% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh1m1TE
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
%   Wie in palh1m1TE_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = palh1m1TE_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m1TE_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1TE_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_jacobiaD_rot_sym_varpar: pkin has to be [23x1] (double)');
JaD_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:59
	% EndTime: 2020-04-13 14:17:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:19:09
	% EndTime: 2020-04-13 14:19:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:20:22
	% EndTime: 2020-04-13 14:29:27
	% DurationCPUTime: 562.58s
	% Computational Cost: add. (12986363->358), mult. (20105870->653), div. (805951->22), fcn. (12596361->21), ass. (0->336)
	t1341 = sin(qJ(2));
	t1342 = sin(pkin(19));
	t1344 = cos(qJ(2));
	t1345 = cos(pkin(19));
	t1282 = t1341 * t1345 - t1342 * t1344;
	t1277 = pkin(7) * t1282;
	t1271 = (-0.2e1 * t1277 + pkin(1)) * pkin(1);
	t1362 = pkin(7) ^ 2;
	t1263 = t1271 + t1362;
	t1253 = pkin(3) ^ 2 - pkin(8) ^ 2 + t1263;
	t1349 = -pkin(8) - pkin(3);
	t1247 = (pkin(7) - t1349) * (pkin(7) + t1349) + t1271;
	t1348 = -pkin(8) + pkin(3);
	t1248 = (pkin(7) - t1348) * (pkin(7) + t1348) + t1271;
	t1238 = t1248 * t1247;
	t1064 = sqrt(-t1238);
	t1272 = -t1341 * t1342 - t1344 * t1345;
	t1270 = pkin(7) * t1272;
	t1260 = t1064 * t1270;
	t1274 = -t1277 + pkin(1);
	t1365 = 0.1e1 / pkin(3);
	t1236 = t1365 * (t1253 * t1274 + t1260);
	t1258 = 0.1e1 / t1263;
	t1407 = t1258 / 0.2e1;
	t1226 = t1236 * t1407;
	t1250 = pkin(7) * t1253;
	t1243 = t1272 * t1250;
	t1237 = t1365 * (t1064 * t1274 - t1243);
	t1232 = t1258 * t1237;
	t1340 = sin(qJ(3));
	t1343 = cos(qJ(3));
	t1205 = t1340 * t1226 + t1343 * t1232 / 0.2e1;
	t1421 = -t1340 / 0.2e1;
	t1206 = t1226 * t1343 + t1232 * t1421;
	t1304 = pkin(23) + pkin(22);
	t1298 = sin(t1304);
	t1299 = cos(t1304);
	t1189 = t1205 * t1299 + t1206 * t1298;
	t1186 = pkin(4) * t1189;
	t1182 = (0.2e1 * t1186 + pkin(5)) * pkin(5);
	t1347 = -pkin(9) - pkin(11);
	t1177 = (pkin(4) - t1347) * (pkin(4) + t1347) + t1182;
	t1346 = -pkin(9) + pkin(11);
	t1178 = (pkin(4) - t1346) * (pkin(4) + t1346) + t1182;
	t1166 = t1178 * t1177;
	t1063 = sqrt(-t1166);
	t1266 = pkin(1) * t1270;
	t1414 = -0.2e1 * t1274;
	t1233 = -t1266 * t1414 + t1243;
	t1338 = pkin(7) * t1064;
	t1259 = t1282 * t1338;
	t1351 = pkin(1) * pkin(7);
	t1425 = (t1247 + t1248) * t1351;
	t1028 = -0.2e1 * t1425 * t1272;
	t1355 = 0.1e1 / t1064;
	t1219 = t1355 * t1028;
	t1408 = t1219 / 0.2e1;
	t1197 = t1365 * (t1270 * t1408 + t1233 + t1259);
	t1194 = t1258 * t1197;
	t1242 = t1282 * t1250;
	t1269 = t1362 * t1272 ^ 2;
	t1264 = pkin(1) * t1269;
	t1430 = t1260 - 0.2e1 * t1264 - t1242;
	t1198 = t1365 * (t1274 * t1408 + t1430);
	t1195 = t1198 * t1407;
	t1352 = 0.1e1 / t1263 ^ 2;
	t1398 = t1351 * t1352;
	t1234 = t1236 * t1398;
	t1229 = t1272 * t1234;
	t1235 = t1237 * t1398;
	t1231 = t1272 * t1235;
	t1204 = t1229 * t1340 + t1231 * t1343;
	t1406 = t1340 / 0.2e1;
	t1169 = t1194 * t1406 + t1195 * t1343 - t1204;
	t1203 = t1229 * t1343 - t1231 * t1340;
	t1405 = -t1343 / 0.2e1;
	t1170 = t1194 * t1405 + t1195 * t1340 + t1203;
	t1417 = -t1205 * t1298 + t1206 * t1299;
	t1428 = t1417 * qJD(3) + (t1169 * t1299 - t1170 * t1298) * qJD(2);
	t1433 = pkin(4) * t1428;
	t1435 = t1063 * t1433;
	t1359 = -0.2e1 * pkin(5);
	t1185 = t1189 * t1359;
	t1361 = pkin(9) ^ 2;
	t1363 = pkin(5) ^ 2;
	t1364 = pkin(4) ^ 2;
	t1180 = -pkin(4) * t1185 + pkin(11) ^ 2 - t1361 + t1363 + t1364;
	t1434 = t1180 * t1433;
	t1431 = pkin(4) * t1417;
	t1432 = t1063 * t1431;
	t1292 = t1341 * t1340;
	t1296 = t1344 * t1343;
	t1042 = t1296 - t1292;
	t1062 = cos(qJ(1));
	t1384 = t1042 * t1062;
	t1404 = t1355 / 0.2e1;
	t1254 = t1270 * t1404;
	t1429 = t1028 * t1254 + t1259;
	t1026 = t1363 + (-t1185 + pkin(4)) * pkin(4);
	t1025 = 0.1e1 / t1026 ^ 2;
	t1357 = 0.1e1 / t1026;
	t1427 = 0.4e1 * t1357 * t1025 * t1363 * t1364;
	t1403 = t1357 / 0.2e1;
	t1041 = -t1340 * t1344 - t1341 * t1343;
	t1397 = qJD(2) + qJD(3);
	t1035 = t1397 * t1041;
	t1135 = pkin(5) * t1428;
	t1133 = pkin(4) * t1135;
	t1389 = t1177 + t1178;
	t1106 = t1389 * t1133;
	t1105 = -0.2e1 * t1106;
	t1356 = 0.1e1 / t1063;
	t1104 = t1356 * t1105;
	t1411 = -pkin(4) / 0.2e1;
	t1103 = t1104 * t1411;
	t1424 = (t1169 * t1298 + t1170 * t1299) * qJD(2) + t1189 * qJD(3);
	t1137 = pkin(4) * t1424;
	t1360 = 0.1e1 / pkin(11);
	t1183 = t1186 + pkin(5);
	t1415 = -0.2e1 * t1183;
	t1091 = t1360 * (-t1063 * t1137 - t1103 * t1417 - t1133 * t1415 + t1434);
	t1134 = t1360 * t1428;
	t1156 = t1180 * t1183 + t1432;
	t1122 = t1156 * t1134;
	t1350 = pkin(4) * pkin(5);
	t1400 = t1025 * t1350;
	t1426 = t1091 * t1400 - t1122 * t1427;
	t1412 = -0.2e1 * t1364;
	t1132 = t1135 * t1412;
	t1409 = t1183 / 0.2e1;
	t1092 = t1360 * (t1104 * t1409 + t1132 * t1417 + t1137 * t1180 + t1435);
	t1179 = pkin(4) * t1180;
	t1157 = t1063 * t1183 - t1179 * t1417;
	t1123 = t1157 * t1134;
	t1303 = t1360 * t1403;
	t1326 = t1364 * t1359;
	t1423 = t1092 * t1400 - t1123 * t1427 + (t1326 * t1424 + t1103) * t1303;
	t1420 = t1356 / 0.2e1;
	t1216 = t1365 * (t1233 + t1429);
	t1210 = t1258 * t1216;
	t1267 = t1274 * t1404;
	t1220 = t1365 * (t1028 * t1267 + t1430);
	t1215 = t1258 * t1220;
	t1174 = t1210 * t1405 + t1215 * t1406 + t1203;
	t1416 = -t1210 * t1421 - t1215 * t1405 - t1204;
	t1413 = -0.8e1 * t1363;
	t1358 = 0.4e1 * pkin(5);
	t1153 = -((pkin(4) + pkin(11)) * (pkin(4) - pkin(11)) + t1182) * t1431 * t1358 + 0.4e1 * t1417 * t1361 * t1350;
	t1150 = t1356 * t1153;
	t1410 = -t1150 / 0.2e1;
	t1155 = t1360 * t1157;
	t1151 = t1155 * t1403;
	t1154 = t1360 * t1156;
	t1152 = t1357 * t1154;
	t1335 = sin(pkin(21));
	t1336 = cos(pkin(21));
	t1016 = t1335 * t1152 / 0.2e1 + t1336 * t1151;
	t1015 = t1384 * t1016;
	t1017 = -t1336 * t1152 / 0.2e1 + t1335 * t1151;
	t1040 = t1041 * t1062;
	t1002 = t1017 * t1040 - t1015;
	t1008 = -t1016 * t1041 - t1017 * t1042;
	t1060 = sin(qJ(1));
	t1290 = t1060 * t1292;
	t1037 = -t1060 * t1296 + t1290;
	t1038 = t1041 * t1060;
	t1283 = t1016 * t1037 + t1017 * t1038;
	t982 = atan2(t1283, t1008);
	t977 = sin(t982);
	t978 = cos(t982);
	t962 = t1008 * t978 + t1283 * t977;
	t960 = 0.1e1 / t962 ^ 2;
	t1320 = t1002 * t960;
	t1006 = 0.1e1 / t1008 ^ 2;
	t1385 = -qJD(3) * t1292 + t1296 * t1397;
	t1034 = -qJD(2) * t1292 + t1385;
	t1381 = t1091 * t1403 - t1122 * t1400;
	t1396 = t1092 * t1403 - t1123 * t1400;
	t969 = t1335 * t1381 + t1336 * t1396;
	t970 = t1335 * t1396 - t1336 * t1381;
	t949 = t1016 * t1034 - t1017 * t1035 - t1041 * t969 - t1042 * t970;
	t1402 = t1006 * t949;
	t1029 = t1037 * qJD(1) + t1035 * t1062;
	t1316 = -t1016 * t1029 - t1384 * t969;
	t1401 = t1025 * t1133;
	t1399 = t1266 * t1352;
	t1164 = t1174 * t1298 - t1299 * t1416;
	t1161 = pkin(5) * t1164;
	t1160 = pkin(4) * t1161;
	t1127 = t1389 * t1160;
	t1126 = 0.2e1 * t1356 * t1127;
	t1162 = pkin(4) * t1164;
	t1165 = t1174 * t1299 + t1298 * t1416;
	t1163 = pkin(4) * t1165;
	t1099 = t1360 * (-t1161 * t1412 * t1417 - t1063 * t1162 + t1126 * t1409 + t1163 * t1180);
	t1149 = t1155 * t1400;
	t1395 = t1099 * t1403 + t1149 * t1164;
	t1184 = t1417 * t1326;
	t1113 = t1360 * (t1150 * t1409 + t1179 * t1189 + t1184 * t1417 + t1432);
	t1393 = t1113 * t1403 - t1149 * t1417;
	t1148 = t1154 * t1400;
	t1392 = t1148 * t1335 + t1149 * t1336;
	t1391 = -t1148 * t1336 + t1149 * t1335;
	t1390 = 0.2e1 * t1389 * t1350;
	t1388 = t1350 * t1415 - t1179;
	t1339 = pkin(4) * t1063;
	t1382 = -t1339 - t1184;
	t1125 = t1126 * t1411;
	t1098 = t1360 * (-t1063 * t1163 - t1125 * t1417 + t1160 * t1415 - t1162 * t1180);
	t1380 = t1098 * t1403 + t1148 * t1164;
	t1112 = t1360 * (-t1189 * t1339 - t1388 * t1417 - t1410 * t1431);
	t1379 = t1112 * t1403 - t1148 * t1417;
	t1377 = 0.1e1 / t1238 * t1219 / 0.4e1;
	t1036 = t1258 * t1352;
	t1265 = pkin(1) ^ 2 * t1269;
	t1369 = 0.4e1 * t1036 * t1237 * t1265 - t1282 * t1235 + (-t1198 - t1220) * t1399;
	t1257 = -0.4e1 * t1265;
	t1370 = t1036 * t1236 * t1257 + t1282 * t1234 + (t1197 + t1216) * t1399;
	t1158 = t1340 * t1370 - t1343 * t1369;
	t1159 = t1340 * t1369 + t1343 * t1370;
	t1268 = qJD(2) * t1282;
	t1027 = 0.2e1 * qJD(2) * t1257 - 0.2e1 * t1268 * t1425;
	t1218 = qJD(2) * t1028;
	t1217 = pkin(7) * t1218;
	t1214 = t1272 * t1217;
	t1023 = t1274 * t1218 * t1377 + t1027 * t1267 + t1214 * t1404 + (-0.6e1 * pkin(1) * t1272 * t1282 * t1362 + t1243 + t1429) * qJD(2);
	t1024 = t1282 * t1217 * t1404 + t1214 * t1377 + t1027 * t1254 - (-t1219 / 0.2e1 + pkin(1) * t1414) * pkin(7) * t1268 + (-t1272 * t1338 + t1242 + 0.4e1 * t1264) * qJD(2);
	t1255 = t1365 * t1258;
	t1249 = t1255 * t1405;
	t1251 = t1340 * t1255;
	t1171 = t1024 * t1249 + t1023 * t1251 / 0.2e1 + t1416 * qJD(3);
	t1172 = -t1024 * t1251 / 0.2e1 + t1023 * t1249 + t1174 * qJD(3);
	t1009 = t1299 * t1172 + t1298 * t1171 + (t1158 * t1299 + t1159 * t1298) * qJD(2);
	t1136 = t1364 * t1428;
	t1124 = t1164 * t1136;
	t1353 = t1356 / t1166;
	t1376 = -t1353 * t1106 * t1127 + (t1009 * t1390 - t1124 * t1413) * t1420;
	t1128 = t1417 * t1136;
	t1375 = t1353 * t1153 * t1105 / 0.4e1 + (t1128 * t1413 + t1390 * t1424) * t1420;
	t1114 = pkin(4) * (t1299 * t1171 - t1298 * t1172 + (-t1158 * t1298 + t1159 * t1299) * qJD(2));
	t1372 = (t1009 * t1382 + t1180 * t1114 - t1125 * t1428 - t1165 * t1132 + t1183 * t1376) * t1303 - t1099 * t1401 + t1423 * t1164;
	t1371 = (-t1189 * t1132 + t1183 * t1375 + t1382 * t1424 - t1410 * t1433 + t1434) * t1303 - t1113 * t1401 - t1423 * t1417;
	t1368 = (t1009 * t1388 - t1063 * t1114 + t1165 * t1103 - t1124 * t1358 + t1125 * t1424 + t1376 * t1431) * t1303 - t1098 * t1401 + t1426 * t1164;
	t1367 = (t1189 * t1103 + t1128 * t1358 + t1137 * t1410 + t1375 * t1431 + t1388 * t1424 - t1435) * t1303 - t1112 * t1401 - t1426 * t1417;
	t959 = 0.1e1 / t962;
	t1004 = t1016 * t1040 + t1017 * t1384;
	t1061 = cos(qJ(4));
	t1059 = sin(qJ(4));
	t1309 = t1059 * t1060;
	t996 = t1004 * t1061 + t1309;
	t990 = 0.1e1 / t996;
	t1005 = 0.1e1 / t1008;
	t991 = 0.1e1 / t996 ^ 2;
	t1284 = -t1008 * t977 + t1283 * t978;
	t1310 = t1283 * t1006;
	t1031 = -qJD(1) * t1384 - t1038 * t1397;
	t1032 = -qJD(1) * t1040 - qJD(2) * t1290 + t1060 * t1385;
	t948 = t1016 * t1031 - t1017 * t1032 + t1037 * t969 + t1038 * t970;
	t1280 = t1005 * t948 - t1310 * t949;
	t997 = t1283 ^ 2;
	t981 = t1006 * t997 + 0.1e1;
	t979 = 0.1e1 / t981;
	t933 = t1280 * t979;
	t924 = t1284 * t933 + t948 * t977 + t949 * t978;
	t1333 = t924 * t959 * t960;
	t1030 = -qJD(1) * t1038 - t1384 * t1397;
	t945 = t1017 * t1030 + t1040 * t970 + t1316;
	t998 = t1002 ^ 2;
	t953 = t960 * t998 + 0.1e1;
	t1334 = 0.2e1 * (t1320 * t945 - t1333 * t998) / t953 ^ 2;
	t1307 = t1060 * t1061;
	t995 = t1004 * t1059 - t1307;
	t1328 = t991 * t995;
	t1308 = t1059 * t1062;
	t1321 = qJD(4) * t995;
	t946 = t1016 * t1030 + t1017 * t1029 + t1040 * t969 + t1384 * t970;
	t944 = qJD(1) * t1308 + t946 * t1061 - t1321;
	t1329 = t944 * t990 * t991;
	t1306 = t1061 * t1062;
	t943 = -qJD(1) * t1306 + qJD(4) * t996 + t946 * t1059;
	t989 = t995 ^ 2;
	t976 = t989 * t991 + 0.1e1;
	t1332 = 0.2e1 * (t1328 * t943 - t1329 * t989) / t976 ^ 2;
	t1317 = t1005 * t1402;
	t1331 = 0.2e1 * (t1310 * t948 - t1317 * t997) / t981 ^ 2;
	t1330 = t943 * t991;
	t1325 = t1009 * t1392 + t1335 * t1368 + t1336 * t1372 + t970;
	t935 = t1009 * t1391 + t1335 * t1372 - t1336 * t1368;
	t1324 = -t935 + t969;
	t1323 = t1335 * t1367 + t1336 * t1371 + t1392 * t1424 + t970;
	t942 = t1335 * t1371 - t1336 * t1367 + t1391 * t1424;
	t1322 = -t942 + t969;
	t1319 = t1002 * t977;
	t1318 = t1002 * t978;
	t973 = t1335 * t1395 - t1336 * t1380;
	t1315 = t1016 - t973;
	t984 = t1335 * t1393 - t1336 * t1379;
	t1314 = t1016 - t984;
	t1313 = t1335 * t1380 + t1336 * t1395 + t1017;
	t1312 = t1335 * t1379 + t1336 * t1393 + t1017;
	t1311 = t1283 * t1005;
	t1305 = -0.2e1 * t1333;
	t1302 = t1005 * t1331;
	t1301 = 0.2e1 * t995 * t1329;
	t1300 = t1002 * t1305;
	t1297 = 0.2e1 * t1283 * t1317;
	t1000 = -t1016 * t1038 + t1017 * t1037;
	t994 = t1000 * t1061 + t1308;
	t993 = t1000 * t1059 - t1306;
	t1281 = -t1059 * t990 + t1061 * t1328;
	t954 = t1037 * t1313 - t1038 * t1315;
	t957 = -t1041 * t1313 + t1042 * t1315;
	t1279 = -t1005 * t954 + t1310 * t957;
	t963 = t1037 * t1312 - t1038 * t1314;
	t966 = -t1041 * t1312 + t1042 * t1314;
	t1278 = -t1005 * t963 + t1310 * t966;
	t974 = 0.1e1 / t976;
	t1276 = t1281 * t974;
	t1275 = t977 + (t1311 * t978 - t977) * t979;
	t1239 = t1281 * t1332 + ((-t944 + t1321) * t991 * t1059 + (qJD(4) * t990 + t1301 - t1330) * t1061) * t974;
	t964 = t1040 * t1314 + t1312 * t1384;
	t955 = t1040 * t1315 + t1313 * t1384;
	t951 = 0.1e1 / t953;
	t947 = t1016 * t1032 + t1017 * t1031 + t1037 * t970 - t1038 * t969;
	t940 = t1275 * t1002;
	t939 = t1278 * t979;
	t937 = t1279 * t979;
	t932 = t1034 * t1312 + t1035 * t1314 - t1041 * t1323 + t1042 * t1322;
	t931 = t1031 * t1312 + t1032 * t1314 + t1037 * t1323 - t1038 * t1322;
	t929 = -t1284 * t939 + t963 * t977 + t966 * t978;
	t928 = t1034 * t1313 + t1035 * t1315 - t1041 * t1325 + t1042 * t1324;
	t927 = t1031 * t1313 + t1032 * t1315 + t1037 * t1325 - t1038 * t1324;
	t925 = -t1284 * t937 + t954 * t977 + t957 * t978;
	t922 = t1278 * t1331 + (t966 * t1297 + t1005 * t931 + (-t1283 * t932 - t948 * t966 - t949 * t963) * t1006) * t979;
	t921 = t1279 * t1331 + (t957 * t1297 + t1005 * t927 + (-t1283 * t928 - t948 * t957 - t949 * t954) * t1006) * t979;
	t1 = [t1005 * t945 * t979 + (-t1402 * t979 - t1302) * t1002, t921, t922, 0; (-t1283 * t959 - t1320 * t940) * t1334 + (t948 * t959 + (-t1283 * t924 + t940 * t945) * t960 + (((-t1311 * t933 * t979 + t1331) * t977 + (-t1283 * t1302 + t933 + (t1280 - t933) * t979) * t978) * t1320 + t1305 * t940 + t1275 * t960 * t945) * t1002) * t951, (-t1320 * t925 - t955 * t959) * t1334 + ((t1029 * t1313 + t1030 * t1315 + t1040 * t1324 + t1325 * t1384) * t959 + t925 * t1300 + (-t955 * t924 + t925 * t945 + (t1283 * t921 - t937 * t948 + t928 + (t1008 * t937 + t954) * t933) * t1318 + (-t1008 * t921 + t937 * t949 + t927 + (t1283 * t937 - t957) * t933) * t1319) * t960) * t951, (-t1320 * t929 - t959 * t964) * t1334 + ((t1029 * t1312 + t1030 * t1314 + t1040 * t1322 + t1323 * t1384) * t959 + t929 * t1300 + (-t964 * t924 + t929 * t945 + (t1283 * t922 - t939 * t948 + t932 + (t1008 * t939 + t963) * t933) * t1318 + (-t1008 * t922 + t939 * t949 + t931 + (t1283 * t939 - t966) * t933) * t1319) * t960) * t951, 0; (t1328 * t994 - t990 * t993) * t1332 + ((qJD(1) * t1307 + qJD(4) * t994 + t947 * t1059) * t990 + t994 * t1301 + (-t993 * t944 - (-qJD(1) * t1309 - qJD(4) * t993 + t1061 * t947) * t995 - t994 * t943) * t991) * t974, -(t1029 * t973 + t1030 * t1313 + t1040 * t1325 + t1384 * t935 + t1316) * t1276 + t1239 * (t1040 * t1313 + t1384 * t973 - t1015), -(t1029 * t984 + t1030 * t1312 + t1040 * t1323 + t1384 * t942 + t1316) * t1276 + t1239 * (t1040 * t1312 + t1384 * t984 - t1015), -t1332 + (0.2e1 * t974 * t1330 + (-0.2e1 * t1329 * t974 - t1332 * t991) * t995) * t995;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:18:02
	% EndTime: 2020-04-13 14:18:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:18:03
	% EndTime: 2020-04-13 14:18:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:18:01
	% EndTime: 2020-04-13 14:18:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiaD_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:18:03
	% EndTime: 2020-04-13 14:18:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiaD_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:18:54
	% EndTime: 2020-04-13 14:18:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
end