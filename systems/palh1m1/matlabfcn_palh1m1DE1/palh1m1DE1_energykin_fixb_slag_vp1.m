% Calculate kinetic energy for
% palh1m1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m1DE1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(23,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m1DE1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_energykin_fixb_slag_vp1: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE1_energykin_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1DE1_energykin_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m1DE1_energykin_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-13 14:50:14
% EndTime: 2020-04-13 14:52:14
% DurationCPUTime: 125.46s
% Computational Cost: add. (2862629->705), mult. (4363601->1202), div. (195296->33), fcn. (2746324->46), ass. (0->449)
t1358 = -2 * pkin(1);
t1204 = sin(qJ(1));
t1191 = t1204 ^ 2;
t1210 = cos(qJ(1));
t1192 = t1210 ^ 2;
t1357 = (-pkin(2) - pkin(13));
t1356 = (-pkin(2) + pkin(13));
t1355 = -pkin(3) - pkin(8);
t1354 = -pkin(8) + pkin(3);
t1353 = -pkin(9) - pkin(11);
t1352 = -pkin(9) + pkin(11);
t1222 = pkin(7) ^ 2;
t1230 = pkin(1) ^ 2;
t1203 = sin(qJ(2));
t1205 = sin(pkin(19));
t1209 = cos(qJ(2));
t1211 = cos(pkin(19));
t1172 = t1203 * t1211 - t1205 * t1209;
t1329 = pkin(7) * t1172;
t1295 = t1329 * t1358 + t1230;
t1156 = t1222 + t1295;
t1292 = pkin(3) ^ 2 - pkin(8) ^ 2;
t1141 = t1156 + t1292;
t1165 = pkin(1) - t1329;
t1175 = t1203 * t1205 + t1209 * t1211;
t1169 = t1175 * qJD(2);
t1168 = t1172 * qJD(2);
t1134 = (pkin(7) - t1355) * (pkin(7) + t1355) + t1295;
t1135 = (pkin(7) - t1354) * (pkin(7) + t1354) + t1295;
t1233 = sqrt(-t1135 * t1134);
t1289 = pkin(1) * pkin(7) * t1169;
t1312 = 0.2e1 * (t1134 + t1135) * t1289 / t1233;
t1280 = -t1312 / 0.2e1;
t1239 = t1168 * t1233 + t1175 * t1280;
t1081 = ((t1165 * t1358 - t1141) * t1169 + t1239) * pkin(7);
t1152 = 0.1e1 / t1156;
t1202 = sin(qJ(3));
t1227 = 0.1e1 / pkin(3);
t1104 = pkin(7) * t1141 * t1175 + t1165 * t1233;
t1208 = cos(qJ(3));
t1301 = t1208 * t1104;
t1305 = t1175 * t1233;
t1103 = -pkin(7) * t1305 + t1141 * t1165;
t1304 = t1202 * t1103;
t1244 = t1304 / 0.2e1 + t1301 / 0.2e1;
t1273 = 0.1e1 / t1156 ^ 2 * t1289;
t1302 = t1208 * t1103;
t1303 = t1202 * t1104;
t1342 = -t1208 / 0.2e1;
t1287 = -0.2e1 * t1169 * t1175;
t1306 = t1169 * t1233;
t1082 = t1165 * t1312 / 0.2e1 + t1222 * pkin(1) * t1287 + (-t1141 * t1168 - t1306) * pkin(7);
t1346 = t1082 / 0.2e1;
t1019 = ((-t1302 + t1303) * t1273 + (qJD(3) * t1244 + t1081 * t1342 + t1202 * t1346) * t1152) * t1227;
t1243 = -t1302 / 0.2e1 + t1303 / 0.2e1;
t1347 = -t1081 / 0.2e1;
t1020 = ((-t1301 - t1304) * t1273 + (qJD(3) * t1243 + t1082 * t1342 + t1202 * t1347) * t1152) * t1227;
t1190 = pkin(23) + pkin(22);
t1186 = sin(t1190);
t1187 = cos(t1190);
t991 = t1019 * t1186 + t1020 * t1187;
t1351 = pkin(4) * t991;
t1225 = pkin(4) ^ 2;
t1224 = pkin(5) ^ 2;
t1307 = t1152 * t1227;
t1091 = t1244 * t1307;
t1092 = t1243 * t1307;
t1069 = -t1091 * t1187 + t1092 * t1186;
t1336 = pkin(5) * t1069;
t1297 = -0.2e1 * pkin(4) * t1336 + t1224;
t1065 = t1225 + t1297;
t1294 = pkin(9) ^ 2 - pkin(11) ^ 2;
t1059 = t1065 + t1294;
t1066 = -pkin(4) + t1336;
t1057 = (pkin(4) - t1353) * (pkin(4) + t1353) + t1297;
t1058 = (pkin(4) - t1352) * (pkin(4) + t1352) + t1297;
t1231 = sqrt(-t1058 * t1057);
t1070 = t1091 * t1186 + t1092 * t1187;
t1335 = pkin(5) * t1070;
t976 = t1059 * t1335 - t1066 * t1231;
t1350 = pkin(5) * t976;
t1063 = 0.1e1 / t1065;
t1348 = t1063 / 0.2e1;
t1196 = sin(pkin(20));
t1200 = cos(pkin(20));
t1247 = t1196 * t1208 + t1200 * t1202;
t1331 = pkin(6) * t1247;
t1288 = pkin(1) * t1331;
t1162 = 0.2e1 * t1288;
t1223 = pkin(6) ^ 2;
t1293 = t1223 + t1230;
t1145 = t1162 + t1293;
t1143 = 0.1e1 / t1145;
t1345 = t1143 / 0.2e1;
t1193 = sin(pkin(23));
t1344 = t1193 / 0.2e1;
t1195 = sin(pkin(21));
t1343 = t1195 / 0.2e1;
t1212 = cos(pkin(18));
t1341 = t1212 / 0.2e1;
t1229 = 0.1e1 / pkin(2);
t1339 = t1229 / 0.2e1;
t1171 = t1196 * t1202 - t1200 * t1208;
t1164 = t1171 * qJD(3);
t1338 = pkin(1) * t1164;
t1337 = pkin(1) * t1203;
t1189 = qJD(2) * t1204;
t1176 = qJD(3) * t1204 + t1189;
t1334 = pkin(5) * t1176;
t1299 = -qJD(2) - qJD(3);
t1177 = t1299 * t1210;
t1333 = pkin(5) * t1177;
t1228 = pkin(2) ^ 2;
t1275 = -pkin(13) ^ 2 + t1293;
t1139 = t1162 + t1228 + t1275;
t1161 = -pkin(1) - t1331;
t1296 = t1162 + t1223;
t1132 = ((pkin(1) - t1357) * (pkin(1) + t1357)) + t1296;
t1133 = ((pkin(1) - t1356) * (pkin(1) + t1356)) + t1296;
t1311 = t1133 * t1132;
t1232 = sqrt(-t1311);
t1330 = pkin(6) * t1171;
t1101 = t1139 * t1330 - t1161 * t1232;
t1332 = pkin(6) * t1101;
t1328 = pkin(16) * t1204;
t1327 = Icges(3,4) * t1203;
t1326 = Icges(3,4) * t1209;
t1140 = t1156 - t1292;
t1166 = pkin(1) * t1172 - pkin(7);
t1102 = -pkin(1) * t1305 - t1140 * t1166;
t1105 = pkin(1) * t1140 * t1175 - t1166 * t1233;
t1206 = sin(pkin(18));
t1221 = 0.1e1 / pkin(8);
t1308 = t1152 * t1221;
t1093 = (t1102 * t1341 - t1206 * t1105 / 0.2e1) * t1308;
t1094 = (t1105 * t1341 + t1102 * t1206 / 0.2e1) * t1308;
t1078 = atan2(t1094, t1093);
t1074 = sin(t1078);
t1325 = Icges(7,4) * t1074;
t1075 = cos(t1078);
t1324 = Icges(7,4) * t1075;
t1291 = pkin(5) * t1351;
t1323 = 0.2e1 / t1231 * (t1057 + t1058) * t1291;
t1080 = ((0.2e1 * pkin(7) * t1166 - t1140) * t1169 + t1239) * pkin(1);
t1083 = t1166 * t1280 + t1230 * pkin(7) * t1287 + (-t1140 * t1168 - t1306) * pkin(1);
t1090 = 0.1e1 / t1093 ^ 2;
t1264 = t1212 * t1273;
t1265 = t1206 * t1273;
t1277 = t1152 * t1341;
t1309 = t1152 * t1206;
t972 = ((t1083 * t1277 + t1105 * t1264 + t1080 * t1309 / 0.2e1 + t1102 * t1265) / t1093 - (t1080 * t1277 + t1102 * t1264 - t1083 * t1309 / 0.2e1 - t1105 * t1265) * t1094 * t1090) / (t1090 * t1094 ^ 2 + 0.1e1) * t1221;
t1322 = t1204 * t972;
t1321 = t1210 * t972;
t1320 = t1231 * t991;
t1197 = cos(pkin(23));
t1088 = (-t1197 * t1103 / 0.2e1 + t1104 * t1344) * t1307;
t1087 = 0.1e1 / t1088 ^ 2;
t1089 = (t1197 * t1104 / 0.2e1 + t1103 * t1344) * t1307;
t1266 = t1197 * t1273;
t1267 = t1193 * t1273;
t1278 = t1152 * t1344;
t1310 = t1152 * t1197;
t971 = ((t1081 * t1278 + t1103 * t1267 + t1104 * t1266 + t1310 * t1346) / t1088 - (t1082 * t1278 - t1103 * t1266 + t1104 * t1267 + t1310 * t1347) * t1089 * t1087) / (t1087 * t1089 ^ 2 + 0.1e1) * t1227;
t1319 = -qJD(2) - t971;
t968 = t971 * t1204 + t1189;
t1318 = qJD(2) * t1210;
t1317 = (Icges(7,5) * t1074 + Icges(7,6) * t1075) * qJD(1);
t1199 = cos(pkin(21));
t1316 = t1063 * t1199;
t1217 = 0.1e1 / pkin(11);
t1315 = t1063 * t1217;
t1314 = t1070 * t1231;
t1290 = pkin(6) * t1338;
t1313 = 0.2e1 * (t1132 + t1133) * t1290 / t1232;
t1300 = t1210 * qJD(1);
t1100 = -t1139 * t1161 - t1232 * t1330;
t1099 = 0.1e1 / t1100 ^ 2;
t1144 = 0.1e1 / t1145 ^ 2;
t1163 = t1247 * qJD(3);
t1281 = -t1313 / 0.2e1;
t1012 = 0.2e1 * (((t1161 * t1281 + (t1139 * t1163 - t1164 * t1232) * pkin(6)) * t1345 + (-t1143 * t1171 * t1223 + t1144 * t1332) * t1338) / t1100 - ((-t1164 * t1139 - t1163 * t1232 + t1171 * t1281) * t1345 + (t1100 * t1144 + t1143 * t1161) * t1338) * t1099 * t1332) * pkin(2) / (t1099 * t1101 ^ 2 + 0.1e1) * t1145 * t1229;
t1298 = -qJD(2) - t1012;
t1009 = t1012 * t1204 + t1189;
t1060 = t1065 - t1294;
t1067 = -pkin(4) * t1069 + pkin(5);
t975 = -pkin(4) * t1314 + t1060 * t1067;
t977 = pkin(4) * t1060 * t1070 + t1067 * t1231;
t958 = (t975 * t1343 + t977 * t1199 / 0.2e1) * t1315;
t959 = (-t975 * t1199 / 0.2e1 + t977 * t1343) * t1315;
t1286 = atan2(t958, t959);
t1285 = -t1323 / 0.2e1;
t1064 = 0.1e1 / t1065 ^ 2;
t1274 = t1064 * t1291;
t1268 = t1199 * t1274;
t1269 = t1195 * t1274;
t1283 = t1063 * t1343;
t990 = t1019 * t1187 - t1020 * t1186;
t1241 = t1070 * t1285 - t990 * t1231;
t940 = ((-0.2e1 * pkin(5) * t1067 - t1060) * t991 + t1241) * pkin(4);
t941 = t1067 * t1323 / 0.2e1 - 0.2e1 * t1225 * t991 * t1335 + (t1060 * t990 - t1320) * pkin(4);
t957 = 0.1e1 / t959 ^ 2;
t880 = ((t940 * t1283 + t975 * t1269 + t941 * t1316 / 0.2e1 + t977 * t1268) / t959 - (-t940 * t1316 / 0.2e1 - t975 * t1268 + t941 * t1283 + t977 * t1269) * t958 * t957) / (t957 * t958 ^ 2 + 0.1e1) * t1217;
t877 = t880 * t1204 + t1176;
t1284 = t1209 * t1189;
t1219 = 0.1e1 / pkin(9);
t1282 = t1219 * t1348;
t1279 = t1143 * t1339;
t1276 = 0.1e1 / pkin(13) * t1339;
t1271 = cos(t1286);
t878 = (-t880 + t1299) * t1210;
t1174 = -t1202 * t1203 + t1208 * t1209;
t1157 = t1174 * t1204;
t1263 = -pkin(5) * t1157 - t1328;
t1262 = -rSges(3,1) * t1203 - rSges(3,2) * t1209;
t1261 = rSges(7,1) * t1075 - rSges(7,2) * t1074;
t1260 = -Icges(3,1) * t1203 - t1326;
t1259 = Icges(7,1) * t1075 - t1325;
t1258 = -Icges(3,2) * t1209 - t1327;
t1257 = -Icges(7,2) * t1074 + t1324;
t1256 = -Icges(3,5) * t1203 - Icges(3,6) * t1209;
t1255 = Icges(7,5) * t1075 - Icges(7,6) * t1074;
t1041 = -Icges(7,6) * t1210 + t1204 * t1257;
t1043 = -Icges(7,5) * t1210 + t1204 * t1259;
t1254 = -t1041 * t1074 + t1043 * t1075;
t1042 = Icges(7,6) * t1204 + t1210 * t1257;
t1044 = Icges(7,5) * t1204 + t1210 * t1259;
t1253 = -t1042 * t1074 + t1044 * t1075;
t1076 = atan2(t1089, t1088);
t1071 = sin(t1076);
t1072 = cos(t1076);
t1047 = -t1071 * t1209 - t1072 * t1203;
t1252 = t1071 * t1203 - t1072 * t1209;
t1086 = atan2(t1101 * t1279, t1100 * t1279);
t1084 = sin(t1086);
t1085 = cos(t1086);
t1061 = -t1084 * t1209 - t1085 * t1203;
t1251 = t1084 * t1203 - t1085 * t1209;
t1148 = -Icges(3,6) * t1210 + t1204 * t1258;
t1150 = -Icges(3,5) * t1210 + t1204 * t1260;
t1250 = t1148 * t1209 + t1150 * t1203;
t1149 = Icges(3,6) * t1204 + t1210 * t1258;
t1151 = Icges(3,5) * t1204 + t1210 * t1260;
t1249 = -t1149 * t1209 - t1151 * t1203;
t1179 = -Icges(3,2) * t1203 + t1326;
t1180 = Icges(3,1) * t1209 - t1327;
t1248 = -t1179 * t1209 - t1180 * t1203;
t1173 = t1202 * t1209 + t1203 * t1208;
t1246 = (-t1191 - t1192) * qJD(2) * t1337;
t1245 = -pkin(1) * t1209 * t1318 + qJD(1) * t1204 * t1337;
t1242 = t1173 * t1333 + t1245;
t1050 = Icges(7,2) * t1075 + t1325;
t1051 = Icges(7,1) * t1074 + t1324;
t1240 = (-t1050 * t1074 + t1051 * t1075) * qJD(1);
t974 = -pkin(5) * t1314 - t1059 * t1066;
t1238 = atan2(t976 * t1282, t974 * t1282);
t1185 = pkin(16) * t1300;
t1237 = t1185 + (-t1203 * t1300 - t1284) * pkin(1);
t1159 = t1174 * t1210;
t1236 = t1157 * t1334 - t1159 * t1333 + t1246;
t1235 = sin(t1238);
t1234 = qJD(1) * t1159 * pkin(5) - t1173 * t1334 + t1237;
t1207 = cos(qJ(4));
t1201 = sin(qJ(4));
t1198 = cos(pkin(22));
t1194 = sin(pkin(22));
t1183 = rSges(2,1) * t1210 - rSges(2,2) * t1204;
t1182 = rSges(3,1) * t1209 - rSges(3,2) * t1203;
t1181 = rSges(2,1) * t1204 + rSges(2,2) * t1210;
t1178 = Icges(3,5) * t1209 - Icges(3,6) * t1203;
t1160 = t1173 * t1210;
t1158 = t1173 * t1204;
t1155 = rSges(3,3) * t1204 + t1210 * t1262;
t1154 = -rSges(3,3) * t1210 + t1204 * t1262;
t1147 = Icges(3,3) * t1204 + t1210 * t1256;
t1146 = -Icges(3,3) * t1210 + t1204 * t1256;
t1138 = t1228 - t1275 - 0.2e1 * t1288;
t1137 = 0.1e1 / t1138 ^ 2;
t1131 = rSges(4,1) * t1173 + rSges(4,2) * t1174;
t1129 = Icges(4,1) * t1173 + Icges(4,4) * t1174;
t1128 = Icges(4,4) * t1173 + Icges(4,2) * t1174;
t1127 = Icges(4,5) * t1173 + Icges(4,6) * t1174;
t1126 = qJD(1) * t1155 - t1182 * t1189 + t1185;
t1125 = -t1182 * t1318 + (-t1154 - t1328) * qJD(1);
t1124 = rSges(4,1) * t1159 - rSges(4,2) * t1160 + rSges(4,3) * t1204;
t1123 = rSges(4,1) * t1157 - rSges(4,2) * t1158 - rSges(4,3) * t1210;
t1122 = Icges(4,1) * t1159 - Icges(4,4) * t1160 + Icges(4,5) * t1204;
t1121 = Icges(4,1) * t1157 - Icges(4,4) * t1158 - Icges(4,5) * t1210;
t1120 = Icges(4,4) * t1159 - Icges(4,2) * t1160 + Icges(4,6) * t1204;
t1119 = Icges(4,4) * t1157 - Icges(4,2) * t1158 - Icges(4,6) * t1210;
t1118 = Icges(4,5) * t1159 - Icges(4,6) * t1160 + Icges(4,3) * t1204;
t1117 = Icges(4,5) * t1157 - Icges(4,6) * t1158 - Icges(4,3) * t1210;
t1116 = (t1154 * t1204 + t1155 * t1210) * qJD(2);
t1110 = atan2(t1232 * t1276, t1138 * t1276);
t1108 = cos(t1110);
t1107 = sin(t1110);
t1098 = qJD(1) * t1124 - t1131 * t1176 + t1237;
t1097 = t1131 * t1177 + (-t1123 - t1328) * qJD(1) + t1245;
t1096 = t1123 * t1176 - t1124 * t1177 + t1246;
t1079 = (0.1e1 / t1138 * t1313 / 0.2e1 - 0.2e1 * t1232 * t1137 * t1290) / (-t1137 * t1311 + 0.1e1);
t1056 = t1061 * t1210;
t1055 = t1251 * t1210;
t1054 = t1061 * t1204;
t1053 = t1251 * t1204;
t1052 = rSges(7,1) * t1074 + rSges(7,2) * t1075;
t1046 = rSges(7,3) * t1204 + t1210 * t1261;
t1045 = -rSges(7,3) * t1210 + t1204 * t1261;
t1040 = Icges(7,3) * t1204 + t1210 * t1255;
t1039 = -Icges(7,3) * t1210 + t1204 * t1255;
t1038 = t1047 * t1210;
t1037 = t1252 * t1210;
t1036 = t1047 * t1204;
t1035 = t1252 * t1204;
t1034 = -rSges(9,1) * t1251 + rSges(9,2) * t1061;
t1033 = -Icges(9,1) * t1251 + Icges(9,4) * t1061;
t1032 = -Icges(9,4) * t1251 + Icges(9,2) * t1061;
t1031 = -Icges(9,5) * t1251 + Icges(9,6) * t1061;
t1028 = rSges(9,1) * t1056 + rSges(9,2) * t1055 + rSges(9,3) * t1204;
t1027 = rSges(9,1) * t1054 + rSges(9,2) * t1053 - rSges(9,3) * t1210;
t1026 = Icges(9,1) * t1056 + Icges(9,4) * t1055 + Icges(9,5) * t1204;
t1025 = Icges(9,1) * t1054 + Icges(9,4) * t1053 - Icges(9,5) * t1210;
t1024 = Icges(9,4) * t1056 + Icges(9,2) * t1055 + Icges(9,6) * t1204;
t1023 = Icges(9,4) * t1054 + Icges(9,2) * t1053 - Icges(9,6) * t1210;
t1022 = Icges(9,5) * t1056 + Icges(9,6) * t1055 + Icges(9,3) * t1204;
t1021 = Icges(9,5) * t1054 + Icges(9,6) * t1053 - Icges(9,3) * t1210;
t1018 = -t1061 * t1107 + t1108 * t1251;
t1017 = -t1061 * t1108 - t1107 * t1251;
t1016 = -t1055 * t1107 - t1056 * t1108;
t1015 = -t1055 * t1108 + t1056 * t1107;
t1014 = -t1053 * t1107 - t1054 * t1108;
t1013 = -t1053 * t1108 + t1054 * t1107;
t1010 = t1298 * t1210;
t1008 = -rSges(8,1) * t1252 + rSges(8,2) * t1047;
t1007 = -Icges(8,1) * t1252 + Icges(8,4) * t1047;
t1006 = -Icges(8,4) * t1252 + Icges(8,2) * t1047;
t1005 = -Icges(8,5) * t1252 + Icges(8,6) * t1047;
t1004 = (-t1047 * t1194 - t1198 * t1252) * pkin(4);
t1003 = rSges(8,1) * t1038 + rSges(8,2) * t1037 + rSges(8,3) * t1204;
t1002 = rSges(8,1) * t1036 + rSges(8,2) * t1035 - rSges(8,3) * t1210;
t1001 = Icges(8,1) * t1038 + Icges(8,4) * t1037 + Icges(8,5) * t1204;
t1000 = Icges(8,1) * t1036 + Icges(8,4) * t1035 - Icges(8,5) * t1210;
t999 = Icges(8,4) * t1038 + Icges(8,2) * t1037 + Icges(8,6) * t1204;
t998 = Icges(8,4) * t1036 + Icges(8,2) * t1035 - Icges(8,6) * t1210;
t997 = Icges(8,5) * t1038 + Icges(8,6) * t1037 + Icges(8,3) * t1204;
t996 = Icges(8,5) * t1036 + Icges(8,6) * t1035 - Icges(8,3) * t1210;
t995 = (-t1037 * t1194 + t1038 * t1198) * pkin(4);
t994 = (-t1035 * t1194 + t1036 * t1198) * pkin(4);
t993 = (-t1079 + t1298) * t1210;
t992 = t1079 * t1204 + t1009;
t989 = rSges(10,1) * t1018 + rSges(10,2) * t1017;
t988 = Icges(10,1) * t1018 + Icges(10,4) * t1017;
t987 = Icges(10,4) * t1018 + Icges(10,2) * t1017;
t986 = Icges(10,5) * t1018 + Icges(10,6) * t1017;
t985 = rSges(10,1) * t1016 + rSges(10,2) * t1015 + rSges(10,3) * t1204;
t984 = rSges(10,1) * t1014 + rSges(10,2) * t1013 - rSges(10,3) * t1210;
t983 = Icges(10,1) * t1016 + Icges(10,4) * t1015 + Icges(10,5) * t1204;
t982 = Icges(10,1) * t1014 + Icges(10,4) * t1013 - Icges(10,5) * t1210;
t981 = Icges(10,4) * t1016 + Icges(10,2) * t1015 + Icges(10,6) * t1204;
t980 = Icges(10,4) * t1014 + Icges(10,2) * t1013 - Icges(10,6) * t1210;
t979 = Icges(10,5) * t1016 + Icges(10,6) * t1015 + Icges(10,3) * t1204;
t978 = Icges(10,5) * t1014 + Icges(10,6) * t1013 - Icges(10,3) * t1210;
t973 = 0.1e1 / t974 ^ 2;
t969 = t1319 * t1210;
t967 = qJD(1) * t1028 - t1009 * t1034 + t1185;
t966 = t1010 * t1034 + (-t1027 - t1328) * qJD(1);
t964 = t1009 * t1027 - t1010 * t1028;
t963 = -t1052 * t1322 + (-pkin(15) * t1210 + t1046) * qJD(1);
t962 = -t1052 * t1321 + (pkin(15) * t1204 - t1045) * qJD(1);
t960 = cos(t1238);
t955 = qJD(1) * t1003 - t1008 * t968 + t1237;
t954 = t1008 * t969 + (-t1002 - t1328) * qJD(1) + t1245;
t953 = (t1045 * t1204 + t1046 * t1210) * t972;
t952 = qJD(1) * t985 - t989 * t992 + t1185 + (qJD(1) * t1056 + t1009 * t1251) * pkin(2);
t951 = -pkin(2) * t1010 * t1251 + t989 * t993 + (-pkin(2) * t1054 - t1328 - t984) * qJD(1);
t950 = t1002 * t968 - t1003 * t969 + t1246;
t948 = sin(t1286);
t947 = -t1194 * t1235 - t1198 * t960;
t946 = t1194 * t960 - t1198 * t1235;
t942 = t984 * t992 - t985 * t993 + (t1009 * t1054 - t1010 * t1056) * pkin(2);
t939 = t1173 * t1271 + t1174 * t948;
t938 = t1173 * t948 - t1174 * t1271;
t937 = qJD(4) * t938 + qJD(1);
t936 = t1159 * t1271 - t1160 * t948;
t935 = t1159 * t948 + t1160 * t1271;
t934 = t1157 * t1271 - t1158 * t948;
t933 = t1157 * t948 + t1158 * t1271;
t932 = t1201 * t1204 + t1207 * t936;
t931 = -t1201 * t936 + t1204 * t1207;
t930 = -t1201 * t1210 + t1207 * t934;
t929 = -t1201 * t934 - t1207 * t1210;
t928 = t1047 * t946 - t1252 * t947;
t927 = t1047 * t947 + t1252 * t946;
t926 = t1037 * t946 + t1038 * t947;
t925 = t1037 * t947 - t1038 * t946;
t924 = t1035 * t946 + t1036 * t947;
t923 = t1035 * t947 - t1036 * t946;
t922 = pkin(10) * t939 + pkin(12) * t938;
t921 = rSges(5,1) * t939 - rSges(5,2) * t938;
t920 = Icges(5,1) * t939 - Icges(5,4) * t938;
t919 = Icges(5,4) * t939 - Icges(5,2) * t938;
t918 = Icges(5,5) * t939 - Icges(5,6) * t938;
t917 = pkin(10) * t936 + pkin(12) * t935;
t916 = pkin(10) * t934 + pkin(12) * t933;
t915 = rSges(5,1) * t936 - rSges(5,2) * t935 + rSges(5,3) * t1204;
t914 = rSges(5,1) * t934 - rSges(5,2) * t933 - rSges(5,3) * t1210;
t913 = Icges(5,1) * t936 - Icges(5,4) * t935 + Icges(5,5) * t1204;
t912 = Icges(5,1) * t934 - Icges(5,4) * t933 - Icges(5,5) * t1210;
t911 = Icges(5,4) * t936 - Icges(5,2) * t935 + Icges(5,6) * t1204;
t910 = Icges(5,4) * t934 - Icges(5,2) * t933 - Icges(5,6) * t1210;
t909 = Icges(5,5) * t936 - Icges(5,6) * t935 + Icges(5,3) * t1204;
t908 = Icges(5,5) * t934 - Icges(5,6) * t933 - Icges(5,3) * t1210;
t907 = rSges(11,1) * t928 + rSges(11,2) * t927;
t906 = Icges(11,1) * t928 + Icges(11,4) * t927;
t905 = Icges(11,4) * t928 + Icges(11,2) * t927;
t904 = Icges(11,5) * t928 + Icges(11,6) * t927;
t903 = rSges(11,1) * t926 + rSges(11,2) * t925 + rSges(11,3) * t1204;
t902 = rSges(11,1) * t924 + rSges(11,2) * t923 - rSges(11,3) * t1210;
t901 = Icges(11,1) * t926 + Icges(11,4) * t925 + Icges(11,5) * t1204;
t900 = Icges(11,1) * t924 + Icges(11,4) * t923 - Icges(11,5) * t1210;
t899 = Icges(11,4) * t926 + Icges(11,2) * t925 + Icges(11,6) * t1204;
t898 = Icges(11,4) * t924 + Icges(11,2) * t923 - Icges(11,6) * t1210;
t897 = Icges(11,5) * t926 + Icges(11,6) * t925 + Icges(11,3) * t1204;
t896 = Icges(11,5) * t924 + Icges(11,6) * t923 - Icges(11,3) * t1210;
t895 = 0.2e1 * (((t1066 * t1285 + (t1059 * t990 - t1320) * pkin(5)) * t1348 + (-t1063 * t1070 * t1224 + t1064 * t1350) * t1351) / t974 - ((-t991 * t1059 + t1241) * t1348 + (t1063 * t1066 + t1064 * t974) * t1351) * t973 * t1350) * pkin(9) * t1065 * t1219 / (t973 * t976 ^ 2 + 0.1e1);
t894 = (-t895 + t1319) * t1210;
t893 = t1204 * t895 + t968;
t892 = rSges(6,3) * t938 + (rSges(6,1) * t1207 - rSges(6,2) * t1201) * t939;
t891 = Icges(6,5) * t938 + (Icges(6,1) * t1207 - Icges(6,4) * t1201) * t939;
t890 = Icges(6,6) * t938 + (Icges(6,4) * t1207 - Icges(6,2) * t1201) * t939;
t889 = Icges(6,3) * t938 + (Icges(6,5) * t1207 - Icges(6,6) * t1201) * t939;
t888 = rSges(6,1) * t932 + rSges(6,2) * t931 + rSges(6,3) * t935;
t887 = rSges(6,1) * t930 + rSges(6,2) * t929 + rSges(6,3) * t933;
t886 = Icges(6,1) * t932 + Icges(6,4) * t931 + Icges(6,5) * t935;
t885 = Icges(6,1) * t930 + Icges(6,4) * t929 + Icges(6,5) * t933;
t884 = Icges(6,4) * t932 + Icges(6,2) * t931 + Icges(6,6) * t935;
t883 = Icges(6,4) * t930 + Icges(6,2) * t929 + Icges(6,6) * t933;
t882 = Icges(6,5) * t932 + Icges(6,6) * t931 + Icges(6,3) * t935;
t881 = Icges(6,5) * t930 + Icges(6,6) * t929 + Icges(6,3) * t933;
t876 = qJD(4) * t933 + t878;
t875 = qJD(4) * t935 + t877;
t874 = -pkin(1) * t1284 - t1004 * t968 - t893 * t907 + t1185 + (-t1210 * t1337 + t903 + t995) * qJD(1);
t873 = t1004 * t969 + t894 * t907 + (-t902 - t994 - t1328) * qJD(1) + t1245;
t872 = qJD(1) * t915 - t877 * t921 + t1234;
t871 = t878 * t921 + (t1263 - t914) * qJD(1) + t1242;
t870 = t893 * t902 - t894 * t903 + t968 * t994 - t969 * t995 + t1246;
t869 = t877 * t914 - t878 * t915 + t1236;
t868 = qJD(1) * t917 - t875 * t892 - t877 * t922 + t888 * t937 + t1234;
t867 = t876 * t892 + t878 * t922 - t887 * t937 + (t1263 - t916) * qJD(1) + t1242;
t866 = t875 * t887 - t876 * t888 + t877 * t916 - t878 * t917 + t1236;
t1 = (m(2) * (t1181 ^ 2 + t1183 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + t1177 * ((-t1118 * t1210 - t1120 * t1158 + t1122 * t1157) * t1176 + (-t1117 * t1210 - t1119 * t1158 + t1121 * t1157) * t1177 + (-t1127 * t1210 - t1128 * t1158 + t1129 * t1157) * qJD(1)) / 0.2e1 + t1176 * ((t1118 * t1204 - t1120 * t1160 + t1122 * t1159) * t1176 + (t1117 * t1204 - t1119 * t1160 + t1121 * t1159) * t1177 + (t1127 * t1204 - t1128 * t1160 + t1129 * t1159) * qJD(1)) / 0.2e1 + m(11) * (t870 ^ 2 + t873 ^ 2 + t874 ^ 2) / 0.2e1 + (((-t1149 * t1203 + t1151 * t1209) * t1204 - (-t1148 * t1203 + t1150 * t1209) * t1210) * qJD(2) + ((t1042 * t1075 + t1044 * t1074) * t1204 - (t1041 * t1075 + t1043 * t1074) * t1210) * t972 + (t1120 * t1174 + t1122 * t1173) * t1176 + (t1119 * t1174 + t1121 * t1173) * t1177 + (t1024 * t1061 - t1026 * t1251) * t1009 + (t1023 * t1061 - t1025 * t1251) * t1010 + (-t1001 * t1252 + t1047 * t999) * t968 + (-t1000 * t1252 + t1047 * t998) * t969 + (t1017 * t981 + t1018 * t983) * t992 + (t1017 * t980 + t1018 * t982) * t993 + (-t911 * t938 + t913 * t939) * t877 + (-t910 * t938 + t912 * t939) * t878 + (t899 * t927 + t901 * t928) * t893 + (t898 * t927 + t900 * t928) * t894 + (-t1179 * t1203 + t1180 * t1209 + t1050 * t1075 + t1051 * t1074 + t1174 * t1128 + t1173 * t1129 + t1061 * t1032 - t1251 * t1033 + t1047 * t1006 - t1252 * t1007 + t1017 * t987 + t1018 * t988 - t919 * t938 + t920 * t939 + t905 * t927 + t906 * t928) * qJD(1)) * qJD(1) / 0.2e1 - ((-t1210 * t1178 + t1204 * t1248) * qJD(1) + (t1192 * t1146 + (t1249 * t1204 + (-t1147 + t1250) * t1210) * t1204) * qJD(2)) * t1318 / 0.2e1 + ((t1204 * t1178 + t1210 * t1248) * qJD(1) + (t1191 * t1147 + (t1250 * t1210 + (-t1146 + t1249) * t1204) * t1210) * qJD(2)) * t1189 / 0.2e1 + t877 * ((t1204 * t909 - t911 * t935 + t913 * t936) * t877 + (t1204 * t908 - t910 * t935 + t912 * t936) * t878 + (t1204 * t918 - t919 * t935 + t920 * t936) * qJD(1)) / 0.2e1 + t893 * ((t1204 * t897 + t899 * t925 + t901 * t926) * t893 + (t1204 * t896 + t898 * t925 + t900 * t926) * t894 + (t1204 * t904 + t905 * t925 + t906 * t926) * qJD(1)) / 0.2e1 + m(8) * (t950 ^ 2 + t954 ^ 2 + t955 ^ 2) / 0.2e1 + m(3) * (t1116 ^ 2 + t1125 ^ 2 + t1126 ^ 2) / 0.2e1 + t992 * ((t1015 * t981 + t1016 * t983 + t1204 * t979) * t992 + (t1015 * t980 + t1016 * t982 + t1204 * t978) * t993 + (t1015 * t987 + t1016 * t988 + t1204 * t986) * qJD(1)) / 0.2e1 + t968 * ((t1001 * t1038 + t1037 * t999 + t1204 * t997) * t968 + (t1000 * t1038 + t1037 * t998 + t1204 * t996) * t969 + (t1005 * t1204 + t1006 * t1037 + t1007 * t1038) * qJD(1)) / 0.2e1 + t1009 * ((t1204 * t1022 + t1055 * t1024 + t1056 * t1026) * t1009 + (t1021 * t1204 + t1023 * t1055 + t1025 * t1056) * t1010 + (t1031 * t1204 + t1032 * t1055 + t1033 * t1056) * qJD(1)) / 0.2e1 + t937 * ((t882 * t875 + t881 * t876 + t889 * t937) * t938 + ((-t1201 * t884 + t1207 * t886) * t875 + (-t1201 * t883 + t1207 * t885) * t876 + (-t1201 * t890 + t1207 * t891) * t937) * t939) / 0.2e1 + t1010 * ((-t1022 * t1210 + t1024 * t1053 + t1026 * t1054) * t1009 + (-t1210 * t1021 + t1053 * t1023 + t1054 * t1025) * t1010 + (-t1031 * t1210 + t1032 * t1053 + t1033 * t1054) * qJD(1)) / 0.2e1 + t894 * ((-t1210 * t897 + t899 * t923 + t901 * t924) * t893 + (-t1210 * t896 + t898 * t923 + t900 * t924) * t894 + (-t1210 * t904 + t905 * t923 + t906 * t924) * qJD(1)) / 0.2e1 + t993 * ((t1013 * t981 + t1014 * t983 - t1210 * t979) * t992 + (t1013 * t980 + t1014 * t982 - t1210 * t978) * t993 + (t1013 * t987 + t1014 * t988 - t1210 * t986) * qJD(1)) / 0.2e1 + t875 * ((t882 * t935 + t884 * t931 + t886 * t932) * t875 + (t881 * t935 + t883 * t931 + t885 * t932) * t876 + (t889 * t935 + t890 * t931 + t891 * t932) * t937) / 0.2e1 + t876 * ((t882 * t933 + t884 * t929 + t886 * t930) * t875 + (t881 * t933 + t883 * t929 + t885 * t930) * t876 + (t889 * t933 + t890 * t929 + t891 * t930) * t937) / 0.2e1 - ((t1039 * t1321 - t1317) * t1210 + (t1240 + (t1253 * t1204 + (-t1040 - t1254) * t1210) * t972) * t1204) * t1321 / 0.2e1 + m(5) * (t869 ^ 2 + t871 ^ 2 + t872 ^ 2) / 0.2e1 + t878 * ((-t1210 * t909 - t911 * t933 + t913 * t934) * t877 + (-t1210 * t908 - t910 * t933 + t912 * t934) * t878 + (-t1210 * t918 - t919 * t933 + t920 * t934) * qJD(1)) / 0.2e1 + t969 * ((t1001 * t1036 + t1035 * t999 - t1210 * t997) * t968 + (t1000 * t1036 + t1035 * t998 - t1210 * t996) * t969 + (-t1005 * t1210 + t1006 * t1035 + t1007 * t1036) * qJD(1)) / 0.2e1 + m(7) * (t953 ^ 2 + t962 ^ 2 + t963 ^ 2) / 0.2e1 + m(4) * (t1096 ^ 2 + t1097 ^ 2 + t1098 ^ 2) / 0.2e1 + ((t1040 * t1322 + t1317) * t1204 + (t1240 + (-t1254 * t1210 + (-t1039 + t1253) * t1204) * t972) * t1210) * t1322 / 0.2e1 + m(6) * (t866 ^ 2 + t867 ^ 2 + t868 ^ 2) / 0.2e1 + m(9) * (t964 ^ 2 + t966 ^ 2 + t967 ^ 2) / 0.2e1 + m(10) * (t942 ^ 2 + t951 ^ 2 + t952 ^ 2) / 0.2e1;
T = t1;
