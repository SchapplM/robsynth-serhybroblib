% Calculate kinetic energy for
% picker2Dm1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
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
% Datum: 2020-05-10 08:43
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = picker2Dm1TE_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(9,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1TE_energykin_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'picker2Dm1TE_energykin_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1TE_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1TE_energykin_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1TE_energykin_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm1TE_energykin_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 00:13:02
% EndTime: 2020-05-10 00:13:57
% DurationCPUTime: 44.95s
% Computational Cost: add. (712376->731), mult. (1986703->1272), div. (24860->44), fcn. (541762->14), ass. (0->544)
t1066 = (pkin(3) ^ 2);
t1371 = -2 * t1066;
t1075 = (pkin(7) ^ 2);
t1050 = -2 * t1075;
t1018 = sin(qJ(2));
t1019 = sin(qJ(1));
t1265 = t1018 * t1019;
t1343 = pkin(1) * pkin(3);
t1162 = t1265 * t1343;
t1136 = qJD(2) * t1162;
t1020 = cos(qJ(2));
t1021 = cos(qJ(1));
t1261 = t1020 * t1021;
t1175 = qJD(1) * t1261;
t1137 = t1175 * t1343;
t1370 = -0.6e1 * t1137 + 0.6e1 * t1136;
t1369 = -t1137 + t1136;
t1061 = (pkin(4) ^ 2);
t1072 = (pkin(1) ^ 2);
t1070 = t1072 ^ 2;
t1263 = t1019 * t1020;
t1207 = pkin(3) * t1263;
t1354 = -4 * pkin(7);
t1234 = t1075 - t1061;
t960 = t1072 + t1234;
t1211 = pkin(3) * t960 * t1354;
t1274 = 4 * pkin(1);
t971 = pkin(3) * t1018;
t953 = t971 + pkin(7);
t1296 = t1021 * t953;
t1238 = -t1066 + t1075;
t982 = t1018 ^ 2;
t1293 = t1066 * t982;
t1348 = 0.2e1 * t1018;
t1251 = pkin(7) * t1348;
t955 = pkin(3) * t1251;
t910 = t1238 + t955 + 0.2e1 * t1293;
t986 = t1021 ^ 2;
t1306 = t910 * t986;
t1366 = 0.2e1 * pkin(3);
t919 = t955 + t960;
t1160 = pkin(1) * t1207;
t937 = -0.2e1 * t1160;
t885 = t937 + t919;
t1254 = t1066 * t1075;
t1362 = 4 * t1072;
t987 = t1066 * t1362;
t958 = t987 - 4 * t1254;
t815 = t958 * t982 + t1018 * t1211 - t1070 - (t1075 - (t1366 + pkin(4)) * pkin(4)) * (t1075 + (t1366 - pkin(4)) * pkin(4)) + (t1050 + (2 * t1061) - (4 * t1066) - 0.4e1 * t1306) * t1072 + (t1207 * t919 - t1296 * t885) * t1274;
t1076 = sqrt(t815);
t1023 = 0.15e2 * t1070;
t1030 = 18 * t1075;
t1031 = -2 * t1061;
t1033 = -6 * t1061;
t1035 = 2 * t1066;
t1074 = t1075 ^ 2;
t1044 = 3 * t1074;
t1083 = pkin(3) * t1066;
t1063 = t1083 ^ 2;
t1084 = t1066 ^ 2;
t1290 = t1072 * t986;
t1212 = 0.12e2 * t1290;
t1213 = 0.12e2 * t1293;
t1048 = 3 * t1075;
t1246 = 15 * t1072 + t1048;
t1273 = 6 * pkin(1);
t981 = t1018 * t982;
t1283 = t1083 * t981;
t1328 = pkin(7) * t1021;
t1163 = -0.2e1 * t1207;
t1079 = pkin(1) * t1072;
t985 = t1021 * t986;
t1285 = t1079 * t985;
t1224 = pkin(7) * t1285;
t1250 = 0.4e1 * t1328;
t1258 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t1248 = t1072 / 0.3e1 + t1075;
t890 = -0.4e1 / 0.9e1 * t1160 + 0.4e1 / 0.9e1 * t1066 - t1061 / 0.9e1 + t1248;
t1006 = 0.2e1 / 0.3e1 * t1066;
t1122 = t1075 - t1160;
t997 = -t1061 / 0.6e1;
t904 = t1006 + t1122 + t997;
t1005 = 0.4e1 / 0.3e1 * t1066;
t977 = t1072 + t1075;
t999 = -t1061 / 0.3e1;
t1180 = t999 + t977;
t939 = t1005 + t1180;
t969 = -t1072 / 0.3e1 + t1075;
t825 = 0.4e1 * t1224 + t939 * t1258 + 0.6e1 * t890 * t1290 + (t1163 * t969 + t1250 * t904) * pkin(1);
t973 = pkin(1) * t1021;
t1231 = 0.4e1 * t973;
t1040 = 6 * t1072;
t1147 = -0.4e1 * t1160;
t1001 = -0.2e1 / 0.3e1 * t1061;
t1049 = 2 * t1075;
t1170 = t1001 + t1006 + t1049;
t1169 = t1001 + t977;
t1276 = t1084 + (t1006 + t1169) * t977;
t840 = t939 * t1147 + t1276 + (t1040 + t1170) * t1066;
t898 = t937 + t939;
t957 = t1238 * t1362;
t827 = pkin(7) * t1231 * t898 + t957 * t986 + t840;
t1259 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t1009 = -t1066 / 0.3e1;
t967 = t1009 + t1075;
t912 = t967 * t937;
t859 = t1259 * t939 + t912;
t988 = 0.10e2 / 0.3e1 * t1072;
t860 = (t988 + t1170) * t1066 + t1276;
t1218 = 0.4e1 * t1290;
t1232 = 0.2e1 * t973;
t956 = pkin(7) * t1232;
t979 = -3 * t1072 + t1075;
t911 = t956 + t1218 + t979;
t1185 = 0.8e1 * t1224;
t976 = -3 * t1066 + t1075;
t927 = t976 * t1185;
t1167 = t1066 + t977;
t950 = -t1061 + t1167;
t954 = t973 + pkin(7);
t970 = t977 ^ 2;
t797 = t825 * t1213 + t927 + t859 * t1212 + t1063 + ((-t1061 + t1066 + t1246) * t1084) + (t1023 + ((t1030 + t1033 + 6 * t1066) * t1072) + t1044 + ((t1031 + t1035) * t1075)) * t1066 + (t970 * t950) + (0.8e1 * t1283 * t911 + 0.6e1 * t827 * t971) * t954 + (-t1207 * t860 + t1328 * t840) * t1273;
t1200 = t954 * t1283;
t1159 = -0.8e1 * t1200;
t1222 = -0.4e1 * t1293;
t1299 = t1018 * t954;
t972 = pkin(3) * t1020;
t1205 = t1072 * t972;
t962 = 3 * t1066 + t977;
t1297 = t1019 * t962;
t1041 = 3 * t1072;
t1242 = t1041 + t1075;
t961 = t1066 + t1242;
t926 = t961 * t972;
t1322 = t1019 * pkin(1);
t930 = t972 - t1322;
t839 = -0.2e1 * t986 * t1205 + t926 + (0.2e1 * t1328 * t930 - t1297) * pkin(1);
t1230 = 0.2e1 * t972;
t1203 = pkin(7) * t1230;
t1111 = -t1019 * t1079 + t1205;
t920 = 0.2e1 * t1111;
t1042 = 2 * t1072;
t975 = t1042 + t1066;
t978 = -t1072 + t1075;
t842 = t978 * t972 + t920 * t986 + (t1019 * t975 + t1021 * t1203) * pkin(1);
t879 = -pkin(1) * t1297 + t926;
t1176 = t1019 * t1259;
t1344 = 4 * t1079;
t1345 = 0.4e1 * t1070;
t1361 = 8 * t1072;
t952 = pkin(3) * t1345 + t1083 * t1361;
t882 = t1020 * t952 + t1176 * t1344;
t1039 = 0.5e1 * t1070;
t1235 = t1074 + t1084;
t1025 = 10 * t1072;
t1245 = t1025 + t1049;
t1252 = t1075 * t1072;
t907 = (t1066 * t1245) + t1039 + t1235 + (6 * t1252);
t1034 = 5 * t1084;
t1046 = 6 * t1075;
t916 = t1034 + (t1025 + t1046) * t1066 + t970;
t806 = t842 * t1222 + t882 * t986 + (-0.4e1 * t839 * t1299 + (-t907 + t1185) * t1020) * pkin(3) + (-0.4e1 * t879 * t1328 + (t916 + t1159) * t1019) * pkin(1);
t922 = t971 + t954;
t793 = t1076 * t806 + t797 * t922;
t1368 = t793 ^ 2;
t1300 = 2 * pkin(7);
t1367 = 2 * pkin(1);
t1233 = pkin(7) * t973;
t1202 = 0.6e1 * t1233;
t1002 = -0.3e1 / 0.2e1 * t1061;
t1253 = t1075 * t1061;
t1060 = t1061 ^ 2;
t1281 = -t1084 / 0.2e1 + t1060 / 0.2e1;
t1146 = t1044 - (3 * t1253) + t1281;
t1237 = t1070 + t1074;
t1277 = t1063 + t977 * ((t1002 + t1049) * t1072 - 0.3e1 / 0.2e1 * t1253 + t1237 + t1281);
t1241 = t1049 - t1061;
t1109 = (t1241 * t1072) + t1060 / 0.6e1 + t1237 - t1253;
t1105 = 0.5e1 / 0.6e1 * t1084 + t1109;
t868 = (t988 + t1241) * t1066 + t1105;
t824 = -0.6e1 * t868 * t1160 + t1277 + (t1023 + ((t1030 - 9 * t1061) * t1072) + t1146) * t1066 + (t1002 + t1246) * t1084;
t1164 = -0.4e1 * t1207;
t1000 = -t1061 / 0.2e1;
t944 = t1000 + t1167;
t1139 = t944 * t1164;
t1126 = pkin(1) * t1139;
t843 = t1126 + ((t1040 + t1241) * t1066) + t1105;
t861 = t1259 * t944 + t912;
t810 = t1202 * t843 + t1212 * t861 + t824 + t927;
t1365 = 0.8e1 * t810;
t1007 = t1066 / 0.3e1;
t1181 = t1007 + t1061 / 0.3e1 + t1049;
t1364 = 0.10e2 / 0.3e1 * t1084 + (-t1072 + t1181) * t1371 + (t1009 + t1180) * t1050;
t949 = t1035 + t1041 + t1234;
t880 = t949 + t955 + t1147;
t1206 = pkin(3) * t1261;
t893 = t1019 * t953 + t1206;
t936 = -pkin(1) + t1207;
t802 = -t880 * t1296 + t893 * t1076 + (t1251 * t936 + t1263 * t949) * pkin(3) + (-0.2e1 * t1306 + (0.2e1 * t982 - 0.4e1) * t1066 - t960) * pkin(1);
t1334 = 0.2e1 / t802;
t1363 = 0.1e1 / t1334;
t1360 = t976 / 0.2e1;
t998 = -t1061 / 0.4e1;
t1359 = t1066 / 0.2e1 + t998;
t1358 = qJD(1) - qJD(2);
t1357 = t1048 - t1061 - t1066;
t1138 = qJD(1) * t1162;
t1267 = qJD(2) * t1020;
t1191 = t954 * t1267;
t1153 = pkin(3) * t1191;
t1356 = t1138 - t1153;
t1178 = qJD(2) * t1265;
t1106 = t1175 - t1178;
t1355 = -8 * pkin(7);
t1353 = -4 * pkin(1);
t1352 = -2 * pkin(2);
t1053 = 2 * pkin(2);
t1022 = cos(pkin(9));
t918 = t1035 + t919;
t844 = t1232 * t910 + t918 * t953;
t848 = t1021 * t918 + (0.4e1 * t986 - 0.2e1) * t953 * pkin(1);
t891 = -t936 + t1296;
t805 = t1019 * t844 + t1076 * t891 + t848 * t972;
t1057 = pkin(5) ^ 2;
t1015 = sin(pkin(8));
t1016 = cos(pkin(8));
t909 = -t1015 * t1019 - t1016 * t1021;
t1331 = pkin(5) * t909;
t896 = -0.2e1 * pkin(1) * t1331;
t1278 = t1057 + t896;
t873 = -(t1367 + pkin(5)) * pkin(5) + t1278;
t874 = pkin(5) * (t1367 - pkin(5)) + t1278;
t1078 = sqrt(-t873 * t874);
t887 = t896 + 0.2e1 * t1057;
t894 = -pkin(1) + t1331;
t908 = t1015 * t1021 - t1016 * t1019;
t821 = pkin(5) * t887 * t908 - t1078 * t894;
t1312 = t805 * t821;
t1286 = t1078 * t908;
t819 = -pkin(5) * t1286 - t887 * t894;
t1315 = t802 * t819;
t1114 = t1315 / 0.4e1 + t1312 / 0.4e1;
t1067 = 0.1e1 / pkin(3);
t881 = t1072 + t1278;
t877 = 0.1e1 / t881;
t1294 = 0.1e1 / pkin(5) * t877;
t1216 = pkin(1) * t1296;
t857 = t1167 + t937 + t955 + 0.2e1 * t1216;
t855 = 0.1e1 / t857;
t1148 = t1067 * t855 * t1294;
t1100 = t1114 * t1148;
t1333 = sin(pkin(9));
t1311 = t819 * t805;
t1314 = t802 * t821;
t1113 = t1314 / 0.4e1 - t1311 / 0.4e1;
t787 = t1113 * t1148;
t781 = t1022 * t1100 - t1333 * t787;
t1351 = 0.1e1 / t781;
t1350 = 0.1e1 / t819;
t1349 = 0.1e1 / t781 ^ 2;
t1347 = -0.8e1 * t1021;
t1346 = -0.2e1 * t1021;
t1301 = 0.2e1 * t1369;
t947 = qJD(2) * t1203;
t1201 = t947 + t1301;
t1271 = qJD(1) * t1019;
t856 = 0.1e1 / t857 ^ 2;
t1309 = ((qJD(2) * t1206 - t1271 * t953) * t1367 + t1201) * t856;
t1194 = -t1309 / 0.2e1;
t1323 = pkin(3) * t1067;
t1336 = t855 / 0.2e1;
t1269 = qJD(1) * t1021;
t1179 = t1019 * t1269;
t1143 = t910 * t1179;
t1266 = qJD(2) * t1066;
t1141 = t1018 * t1020 * t1266;
t1128 = 0.4e1 * t1141;
t902 = t947 + t1128;
t1307 = t902 * t986;
t1313 = ((0.8e1 * t1143 - 0.4e1 * t1307) * t1072 + (t1348 * t958 + t1211) * t1267 + (-0.4e1 * t1201 * t1296 + (0.8e1 * pkin(7) * t1020 ^ 2 * t1266 + 0.4e1 * qJD(1) * t885 * t953) * t1019 + 0.4e1 * (t919 * t1175 + (-t1261 * t885 - t1265 * t919) * qJD(2)) * pkin(3)) * pkin(1)) / t1076;
t1196 = t1313 / 0.2e1;
t1260 = t1021 * t1076;
t1262 = t1019 * t1076;
t1270 = qJD(1) * t1020;
t788 = t891 * t1196 + (t844 * t1021 - t1262 * t953) * qJD(1) + (t1021 * t902 - t1271 * t910) * t1019 * t1367 + ((-t1260 + (-t918 - 0.8e1 * t1216) * t1019) * t1270 + ((-t848 + t1262) * t1018 + (t1260 + (t1300 * t953 + t918) * t1019 + (t1367 * t986 - pkin(1) + t1328) * t1020 * t1366) * t1020) * qJD(2)) * pkin(3);
t1102 = 0.4e1 * t1106;
t1099 = t1102 * t973;
t1264 = t1018 * t1021;
t1119 = -t1263 + t1264;
t789 = t893 * t1196 + (-t947 * t1021 + (t1019 * t880 + t1260) * qJD(1)) * t953 + (t1128 + 0.4e1 * t1143 - 0.2e1 * t1307) * pkin(1) + ((t1021 * t949 - t1262) * t1270 + (-t1076 * t1119 - t1261 * t880 - t1265 * t949) * qJD(2) + (t1106 * t971 + t1267 * t936) * t1300 + t953 * t1099) * pkin(3);
t1096 = t802 ^ 2;
t801 = 0.1e1 / t1096;
t804 = t805 ^ 2;
t754 = qJD(1) + ((t1194 * t805 + t1336 * t788) * t1334 - 0.2e1 * (t1194 * t802 + t1336 * t789) * t805 * t801) / (t801 * t804 + 0.1e1) * t857 * t1323;
t1325 = pkin(1) * qJD(1);
t1217 = t855 * t1325;
t1155 = t1067 * t1217;
t800 = t1155 * t1363;
t1342 = -pkin(2) * t754 / 0.2e1 - t800 / 0.2e1;
t791 = 0.1e1 / t793;
t1341 = t791 / 0.4e1;
t1340 = t802 / 0.4e1;
t1339 = -t805 / 0.2e1;
t1338 = t805 / 0.4e1;
t1337 = t819 / 0.4e1;
t1047 = 4 * t1075;
t1183 = t1006 + 0.2e1 / 0.3e1 * t1061 + t1047;
t1184 = t1005 + 0.4e1 / 0.3e1 * t1061 + t1050;
t989 = -0.20e2 / 0.3e1 * t1072;
t1335 = t1084 / 0.2e1 - (t989 + t1183) * t1066 / 0.2e1 + 0.3e1 / 0.2e1 * t1070 - t1184 * t1072 / 0.2e1 - t1074 / 0.2e1;
t878 = 0.1e1 / t881 ^ 2;
t900 = t908 * qJD(1);
t1308 = t878 * t900;
t1161 = pkin(1) * t855 * t1308;
t1193 = t900 * t1367;
t1310 = (t873 + t874) * pkin(5) * t1193 / t1078;
t1195 = t1310 / 0.2e1;
t901 = t909 * qJD(1);
t1110 = -t901 * t1078 + t1195 * t908;
t798 = (-(t1367 * t894 - t887) * t900 + t1110) * pkin(5);
t1229 = 0.2e1 * t900 * t908;
t1287 = t1078 * t900;
t799 = t894 * t1195 + t1057 * pkin(1) * t1229 + (t887 * t901 + t1287) * pkin(5);
t750 = (-(t1315 / 0.2e1 + t1312 / 0.2e1) * t1161 + (-t1114 * t1309 + (t789 * t1337 + t798 * t1340 + t788 * t821 / 0.4e1 + t799 * t1338) * t855) * t1294) * t1067;
t751 = (-(-t1314 / 0.2e1 + t1311 / 0.2e1) * t1161 + (t1113 * t1309 + (-t789 * t821 / 0.4e1 - t802 * t799 / 0.4e1 + t798 * t1338 + t788 * t1337) * t855) * t1294) * t1067;
t735 = t1022 * t750 + t1333 * t751;
t1332 = pkin(2) * t735;
t1330 = pkin(6) * t781;
t783 = -t1022 * t787 - t1100 * t1333;
t1329 = pkin(6) * t783;
t1013 = t1072 / 0.2e1;
t1327 = -t1076 / 0.4e1;
t1326 = t1076 / 0.4e1;
t1324 = pkin(1) * t986;
t1055 = pkin(6) ^ 2;
t778 = t1330 * t1053;
t1279 = t1055 + t778;
t771 = -(t1053 + pkin(6)) * pkin(6) + t1279;
t772 = pkin(6) * (t1053 - pkin(6)) + t1279;
t1316 = t771 * t772;
t1077 = sqrt(-t1316);
t1321 = 0.2e1 * (-t771 - t772) * pkin(6) * t1332 / t1077;
t1320 = t735 * t783;
t1125 = t1083 * t982 * t1191;
t1115 = -0.24e2 * t1125;
t1116 = t969 * t1136;
t1117 = t969 * t1137;
t1284 = t1079 * t986;
t1223 = pkin(7) * t1284;
t1134 = t1223 * t1271;
t1120 = -0.24e2 * t1134;
t1140 = t1072 * t1179;
t1127 = -0.24e2 * t1140;
t1268 = qJD(2) * t1018;
t1208 = pkin(3) * t1268;
t1129 = t1208 * t1285;
t1192 = qJD(1) * t1283;
t1154 = pkin(1) * t1192;
t1131 = t1019 * t1154;
t1142 = t1259 * t1269;
t1210 = pkin(1) * t1271;
t1165 = pkin(7) * t1210;
t1149 = -0.6e1 * t1165;
t1151 = -0.2e1 * t1179;
t1209 = pkin(1) * t1269;
t1152 = t962 * t1209;
t1166 = t1290 * t1343;
t1177 = qJD(2) * t1264;
t1215 = pkin(3) * t1299;
t1302 = 0.4e1 * t1369 * t939;
t1103 = t1106 * pkin(3);
t1098 = pkin(1) * t967 * t1103 * t1290;
t1303 = t976 * t1120 - 0.24e2 * t1098;
t865 = t1072 * t1103 * t1250;
t917 = pkin(3) * t1267 - t1210;
t948 = -0.2e1 * t1165;
t983 = t1019 ^ 2;
t1319 = ((t954 * t1154 * t1347 + t983 * t1192 * t1361 + t1115 * t1322 + (0.2e1 * (-t1072 * t1208 - t1079 * t1269) * t986 + t920 * t1151 - t978 * t1208 + (t975 * t1269 + (-qJD(1) * t1263 - t1177) * pkin(7) * t1366) * pkin(1)) * t1222 - 0.8e1 * t842 * t1141 - 0.4e1 * (-t1152 + (0.4e1 * t1020 * t1140 + (-t961 + 0.2e1 * t1290) * t1268) * pkin(3) + ((-t1208 - t1209) * t973 - t930 * t1210) * t1300) * t1215 + t1129 * t1355 + t1120 * t972 + (t1142 * t1344 - t1268 * t952) * t986 + t882 * t1151 - 0.4e1 * (-t1208 * t961 - t1152) * t1233 + 0.4e1 * t879 * t1165 + t907 * t1208 + t916 * t1209 + 0.4e1 * t1356 * t839) * t1076 + t806 * t1196 + (0.8e1 * (t948 - 0.8e1 * t1140) * t1200 + (-0.12e2 * t1134 + 0.6e1 * (-0.4e1 / 0.9e1 * t1175 + 0.4e1 / 0.9e1 * t1178) * t1166 - 0.12e2 * t890 * t1140 - t865 - 0.4e1 * t904 * t1165 - 0.2e1 * t1117 + 0.2e1 * t1116) * t1213 + 0.24e2 * t825 * t1141 + 0.6e1 * (t957 * t1151 + (t1021 * t1301 - t1271 * t898) * pkin(7) * t1274 + t1302) * t1215 + t859 * t1127 + t1302 * t1202 + t840 * t1149 + t1303 + (-0.8e1 * t1131 + 0.24e2 * t1125) * t911 + t1370 * t860 - 0.6e1 * t1356 * t827) * t922 + t797 * t917) / t1368;
t1088 = pkin(7) * t1075;
t1240 = t1060 - t1084;
t1112 = 6 * t1074 + t1240 - 6 * t1253;
t1171 = 0.32e2 * t1224;
t1305 = t970 * (-t1066 + t960);
t1172 = 0.16e2 * t1224;
t984 = t986 ^ 2;
t1292 = t1070 * t984;
t1221 = 0.8e1 * t1292;
t1249 = 0.6e1 * t1328;
t1280 = 0.4e1 / 0.7e1 * t1075 - t1061 / 0.7e1;
t933 = t1075 + t1066 / 0.4e1 + t1072 / 0.4e1 - t1061 / 0.8e1;
t836 = -0.32e2 / 0.21e2 * t933 * t1160 + 0.5e1 / 0.42e2 * t1084 + (0.16e2 / 0.21e2 * t1072 + t1280) * t1066 + t1070 / 0.7e1 + t1280 * t1072 + t1074 - 0.3e1 / 0.7e1 * t1253 + t1060 / 0.42e2;
t1011 = 0.4e1 / 0.3e1 * t1072;
t935 = t1248 + t1359;
t838 = -0.8e1 / 0.3e1 * t935 * t1160 + 0.5e1 / 0.18e2 * t1084 + (t1011 + t999) * t1066 + t1074 - t1070 / 0.3e1 + t1060 / 0.18e2 + (t1005 + 0.2e1 / 0.3e1 * t1072 + t1001) * t1075;
t899 = -t1084 / 0.6e1 + t1109;
t1247 = t1013 + t1075;
t903 = -0.2e1 / 0.3e1 * t1160 + t998 + t1247;
t964 = (t1047 + t1061) * t1072;
t1010 = -0.2e1 / 0.3e1 * t1066;
t968 = t1010 + t1075;
t809 = t968 * t1221 + t903 * t1172 + 0.14e2 * t836 * t1290 + (t978 * t1084) + (-0.10e2 / 0.3e1 * t1070 + (2 * t1074) - t1253 + t964) * t1066 + t899 * t1258 + (t1139 * t969 + t1249 * t838) * pkin(1);
t1045 = 8 * t1075;
t886 = t1079 * t1164 + t987 + t1345 + ((t1031 + t1045) * t1072);
t897 = -t1072 + t1122 + t1359;
t826 = t1185 + t886 * t986 + t944 * t979 + (t1163 * t1258 + t1250 * t897) * pkin(1);
t1236 = t1074 - t1070;
t828 = t967 * t1126 - t1063 + (-t988 - t1234) * t1084 + (t964 + t1084 / 0.6e1 - t1060 / 0.6e1 + t1236) * t1066 + t899 * t1075;
t1032 = -5 * t1061;
t1038 = 0.7e1 * t1070;
t832 = (t1002 + t1048 + (7 * t1072)) * t1084 + (t1038 + ((t1032 + 10 * t1075) * t1072) + t1146) * t1066 + t1277;
t959 = -12 * pkin(7) * t1079 + t1088 * t1274;
t966 = -0.8e1 * t1070 + (12 * t1252);
t850 = t1021 * t959 + t966 * t986 + t1172 + t1221 + t1237 - (6 * t1252);
t863 = t1259 * t937 + t944 * t976;
t928 = 0.16e2 * (t1235 - 6 * t1254) * t1070;
t974 = -30 * t1061 + 60 * t1075;
t980 = t982 ^ 2;
t785 = t928 * t984 + t863 * t1171 + 0.24e2 * t828 * t1290 + (t1031 + t1047 + 28 * t1072) * t1063 + (t950 * t1305) + (0.32e2 * t1283 * t826 + t1365 * t971) * t954 + ((t1033 * t1074) + (t1040 * t1112) + (t1049 * t1240) + t974 * t1070 + (28 * t1079 ^ 2) + (4 * t1088 ^ 2) + 0.24e2 * t809 * t982) * t1066 + 0.8e1 * (-t1207 * t832 + t1328 * t824) * pkin(1) + ((t1072 * t974) + 0.16e2 * t850 * t980 + 0.70e2 * t1070 + t1084 + t1112) * t1084;
t1282 = t1084 * t980;
t905 = 0.4e1 / 0.3e1 * t1290 + t956 + t969;
t1156 = -0.24e2 * t905 * t1282;
t1214 = -0.12e2 * t1293;
t1220 = -0.6e1 * t1290;
t1291 = t1070 * t985;
t1190 = t1019 * t1291;
t931 = 0.7e1 / 0.6e1 * t1066 + t997 + t1247;
t1182 = t1007 + t1075 + t997;
t934 = t1011 + t1182;
t858 = -t1322 * t931 + t934 * t972;
t864 = (t1066 * t978) - 0.5e1 / 0.3e1 * t1070 + t1181 * t1072 + t1075 * (t999 + t967);
t941 = t1072 + t1182;
t963 = t1035 + t978;
t870 = t941 * t972 - t963 * t1322 / 0.2e1;
t811 = t1190 * t1354 + t858 * t1218 + (-0.8e1 / 0.3e1 * t1292 + t864) * t972 + (t1019 * t1335 + t1250 * t870) * pkin(1);
t1204 = t1079 * t972;
t1219 = -0.4e1 * t1290;
t1244 = t1031 + t1371;
t1168 = t1046 + t1244;
t869 = t1034 + ((t1025 + t1168) * t1066) + (t1010 + t1169) * t977;
t1298 = t1019 * t869;
t854 = t1084 + (t1001 + t1010 + t1245) * t1066 + t1039 + (t1168 * t1072) + t1075 * (t1001 + t968);
t847 = t854 * t972;
t938 = 0.8e1 / 0.3e1 * t1066 + t1180;
t940 = t1006 + t999 + t1242;
t862 = -t1322 * t938 + t940 * t972;
t945 = 0.5e1 / 0.6e1 * t1066 + t1013 + t997;
t876 = pkin(1) * t1176 + t1230 * t945;
t812 = t985 * t1204 * t1355 + t876 * t1219 + t847 + (t1250 * t862 - t1298) * pkin(1);
t829 = -pkin(1) * t1298 + t847;
t867 = -(3 * t1084) + (t989 + t1184) * t1066 + t1183 * t1072 + t1236;
t830 = t1322 * t1364 + t867 * t972;
t1243 = t1032 - 5 * t1066;
t946 = t1357 * t1025;
t831 = t1063 + ((21 * t1072 + t1357) * t1084) + ((t1075 * t1244) + t1044 + 0.35e2 * t1070 + t946) * t1066 + (t1038 + ((t1045 + t1243) * t1072) + (t1075 * (-t1066 + t1234))) * t977;
t921 = 0.4e1 * t1111;
t929 = t972 + 0.2e1 * t1322;
t943 = t1000 + t961;
t833 = t979 * t972 + t921 * t986 + (t1019 * t943 + t1328 * t929) * t1367;
t841 = 0.7e1 * t1063 + ((35 * t1072 + 15 * t1075 + t1243) * t1084) + (0.21e2 * t1070 + t946 + (9 * t1074) + ((t1033 - 6 * t1066) * t1075)) * t1066 + t1305;
t932 = t1075 + 0.5e1 / 0.2e1 * t1066 + 0.3e1 / 0.2e1 * t1072 + t1000;
t872 = t1322 * t1360 + t932 * t972;
t794 = t872 * t1172 + t833 * t1159 + t811 * t1214 + t830 * t1220 + (-0.6e1 * t812 * t1299 + (0.24e2 * t1292 * t967 - t831) * t1020) * pkin(3) + (-0.6e1 * t829 * t1328 + (t841 + t1156) * t1019) * pkin(1);
t770 = t1076 * t794 + t785 * t922;
t1318 = 0.1e1 / t770 ^ 2 * t1368;
t1317 = t770 * t791;
t1158 = t805 * t1217;
t1132 = t1158 / 0.4e1;
t1056 = 0.1e1 / pkin(6);
t1082 = pkin(2) ^ 2;
t775 = t1082 + t1279;
t1295 = t1056 / t775;
t776 = t778 + 0.2e1 * t1055;
t777 = -pkin(2) - t1330;
t738 = t1077 * t1329 - t776 * t777;
t739 = -t1077 * t777 - t1329 * t776;
t726 = (t1067 * t1132 * t739 + t1342 * t738) * t1295;
t1062 = 0.1e1 / pkin(4);
t1068 = 0.1e1 / pkin(3) ^ 2;
t1133 = -t1158 / 0.4e1;
t1186 = t1067 * (pkin(3) * t754 + t800) / 0.2e1;
t1255 = t1062 * t1076;
t741 = t1062 * t1186 * t1317 + t1068 * t1133 * t1255;
t1304 = t1370 * t868;
t1289 = 0.1e1 / pkin(1) * t877;
t1288 = t1077 * t735;
t1257 = t1056 / pkin(2);
t1256 = t1062 * t1068;
t1227 = 0.1e1 / t775 ^ 2 * t1332;
t1226 = pkin(5) * t1308;
t1199 = -t1319 / 0.4e1;
t1198 = t791 * t1340;
t1197 = t791 * t1338;
t1189 = t1019 * t1282;
t1188 = t855 * t1256;
t1187 = t1295 / 0.2e1;
t1174 = t1077 * t1257;
t1150 = pkin(3) * t1177;
t1145 = qJD(1) * t1190;
t1144 = t1084 * t981 * t1267;
t1130 = t1208 * t1292;
t736 = -t1022 * t751 + t1333 * t750;
t737 = 0.1e1 / t738 ^ 2;
t716 = t754 + 0.2e1 * (((-t777 * t1321 / 0.2e1 + t1055 * t1320 * t1352 + (t736 * t776 + t1288) * pkin(6)) * t1187 - t739 * t1227) / t738 - (-t738 * t1227 + (t783 * t1321 / 0.2e1 - t1077 * t736 + (t1352 * t777 + t776) * t735) * pkin(6) * t1187) * t739 * t737) * pkin(6) / (t737 * t739 ^ 2 + 0.1e1) * t775;
t1101 = t1106 * t1367;
t1104 = -t1208 * t854 - t1209 * t869;
t1121 = -0.48e2 * t1134;
t747 = (t1156 * t1209 - 0.24e2 * pkin(1) * (-0.8e1 / 0.3e1 * t1140 + t948) * t1189 - 0.96e2 * t905 * t1144 * t1322 + ((-t979 + t1219) * t1208 + 0.2e1 * (pkin(1) * t943 - t1019 * t921 - 0.2e1 * t1284) * t1269 + (-0.2e1 * t1150 + (-0.2e1 * t1019 * t929 + 0.4e1 * t1324) * qJD(1)) * pkin(1) * pkin(7)) * t1159 + (0.8e1 / 0.3e1 * t1130 + (-t1208 * t934 - t1209 * t931) * t1218 - t864 * t1208 + t1209 * t1335 + (0.32e2 / 0.3e1 * t1291 * t972 + t1072 * t858 * t1347) * t1271 + (t941 * t1150 * t1353 + ((0.12e2 * t983 * t986 - 0.4e1 * t984) * t1070 + (-0.4e1 * t1019 * t870 - 0.2e1 * t1324 * t963) * pkin(1)) * qJD(1)) * pkin(7)) * t1214 - 0.24e2 * t811 * t1141 - 0.6e1 * ((-0.4e1 * (pkin(1) * t1142 - 0.2e1 * t1208 * t945) * t986 + 0.8e1 * t876 * t1179) * t1072 + (0.8e1 * t1129 + (-t1208 * t940 - t1209 * t938) * t1231 + (0.24e2 * t1204 * t986 + t1353 * t862) * t1271) * pkin(7) + t1104) * t1215 + (-t932 * t1208 + t1209 * t1360) * t1172 + t872 * t1121 + (-t1208 * t867 + t1209 * t1364) * t1220 + 0.12e2 * t830 * t1140 - 0.6e1 * t1104 * t1233 + 0.6e1 * t829 * t1165 + t831 * t1208 + t841 * t1209 + (-0.96e2 * t1145 * t972 - 0.24e2 * t1130) * t967 + (0.8e1 * t1131 + t1115) * t833 + 0.6e1 * t1356 * t812) * t1076 + t794 * t1196 + (0.16e2 * (t1346 * t966 - 0.48e2 * t1223 - 0.32e2 * t1291 - t959) * qJD(1) * t1189 + 0.64e2 * t850 * t1144 + 0.32e2 * (-t865 + (t886 * t1346 + (t1353 * t897 - 0.24e2 * t1284) * pkin(7)) * t1271 + (-t1101 * t1258 - t1102 * t1284) * pkin(3)) * t1200 + 0.24e2 * (-0.32e2 * t968 * t1145 + (-0.2e1 / 0.3e1 * t1175 + 0.2e1 / 0.3e1 * t1178) * t1172 * t1343 + t903 * t1121 - 0.28e2 * t836 * t1140 + (-0.8e1 / 0.3e1 * t1175 + 0.8e1 / 0.3e1 * t1178) * t935 * pkin(3) * t1072 * t1249 + t838 * t1149 + 0.4e1 * (-t1117 + t1116) * t944 - 0.64e2 / 0.3e1 * t1106 * t933 * t1166) * t1293 + 0.48e2 * t809 * t1141 - 0.8e1 * t810 * t1138 + 0.8e1 * (t861 * t1127 + (-pkin(3) * t1099 * t944 - t1271 * t843) * pkin(7) * t1273 + t1303 + t1304) * t1215 + t1153 * t1365 - 0.4e1 * t928 * t985 * t1271 - pkin(3) * t1101 * t1171 * t1259 - 0.96e2 * t863 * t1134 - 0.96e2 * t944 * t1098 - 0.48e2 * t828 * t1140 + 0.8e1 * t1304 * t1233 - 0.8e1 * t824 * t1165 + 0.8e1 * t1369 * t832 + (-0.32e2 * t1131 + 0.96e2 * t1125) * t826) * t922 + t785 * t917;
t724 = t754 + (0.1e1 / t770 * t793 * t1196 - 0.2e1 * (t747 * t791 / 0.2e1 - t770 * t1319 / 0.2e1) * pkin(4) * t1255 * t1318 * t1323) / (t1318 * t815 + 0.1e1);
t1118 = t1261 + t1265;
t1108 = t1197 * t770 + t1326 * t802;
t1107 = t1198 * t770 + t1327 * t805;
t1097 = t819 ^ 2;
t895 = -pkin(1) * t909 + pkin(5);
t888 = t896 + t1042;
t853 = t1358 * t1118;
t852 = t1358 * t1119;
t822 = pkin(1) * t888 * t908 + t1078 * t895;
t820 = -pkin(1) * t1286 + t888 * t895;
t818 = t821 ^ 2;
t817 = 0.1e1 / t820 ^ 2;
t816 = 0.1e1 / t1097;
t768 = (-((-t895 * t1310 / 0.2e1 + t1072 * pkin(5) * t1229 + (t888 * t901 + t1287) * pkin(1)) * t1289 / 0.2e1 - t822 * t1226) / t820 - (t820 * t1226 - (-(-0.2e1 * pkin(5) * t895 - t888) * t900 + t1110) * pkin(1) * t1289 / 0.2e1) * t822 * t817) / (t817 * t822 ^ 2 + 0.1e1) * t881 * t1367;
t767 = qJD(1) + (t799 * t1350 * t1294 + (-t798 * t816 * t1294 - (-t819 * t816 + t1350) * t878 * t1193) * t821) * t881 / (t816 * t818 + 0.1e1) * pkin(5);
t762 = t1108 * t1188;
t761 = t1107 * t1188;
t756 = (t783 * t1339 - t781 * t802 / 0.2e1) * t1155;
t755 = (t781 * t1339 + t1363 * t783) * t1155;
t745 = t1118 * t762 - t1119 * t761;
t744 = -t1118 * t761 - t1119 * t762;
t743 = 0.1e1 / t745 ^ 2;
t740 = (t1068 * t1132 * t1317 + t1076 * t1186) * t1062;
t729 = (-t1108 * t1309 + (t789 * t1326 + t802 * t1313 / 0.8e1 + t747 * t1197 + (t1199 * t805 + t1341 * t788) * t770) * t855) * t1256;
t728 = (-t1107 * t1309 + (t747 * t1198 + t788 * t1327 - t805 * t1313 / 0.8e1 + (t1199 * t802 + t1341 * t789) * t770) * t855) * t1256;
t725 = (t1067 * t1133 * t738 + t1342 * t739) * t1295;
t723 = pkin(4) * t724 + t741;
t722 = (t1320 * t1349 + t1351 * t736) / (t1349 * t783 ^ 2 + 0.1e1) + t754;
t721 = -t723 * t744 + t740 * t745;
t720 = t723 * t745 + t740 * t744;
t717 = (-(-t1118 * t728 - t1119 * t729 - t761 * t852 + t762 * t853) / t745 - (-t1118 * t729 + t1119 * t728 - t761 * t853 - t762 * t852) * t744 * t743) / (t743 * t744 ^ 2 + 0.1e1) + t724;
t715 = pkin(6) * t716 + t726;
t714 = (-t1351 * t1321 / 0.4e1 + t1349 * t1288 / 0.2e1) / (0.1e1 - 0.1e1 / pkin(6) ^ 2 / t1082 * t1349 * t1316 / 0.4e1) * t1257 + t716;
t713 = t781 * t725 - t715 * t1174 / 0.2e1;
t712 = t725 * t1174 / 0.2e1 + t781 * t715;
t1 = m(4) * (t725 ^ 2 + t726 ^ 2) / 0.2e1 + qJD(2) ^ 2 * Ifges(8,3) / 0.2e1 + t768 ^ 2 * Ifges(6,3) / 0.2e1 + m(5) * (t740 ^ 2 + t741 ^ 2) / 0.2e1 + m(10) * (t712 ^ 2 + t713 ^ 2) / 0.2e1 + m(7) * (t755 ^ 2 + t756 ^ 2) / 0.2e1 + m(11) * (t720 ^ 2 + t721 ^ 2) / 0.2e1 + (mrSges(3,2) * t1155 * t1339 + mrSges(3,1) * t800 + Ifges(3,3) * t754 / 0.2e1) * t754 + (t741 * mrSges(5,1) - t740 * mrSges(5,2) + Ifges(5,3) * t724 / 0.2e1) * t724 + (t756 * mrSges(7,1) - t755 * mrSges(7,2) + Ifges(7,3) * t722 / 0.2e1) * t722 + (t720 * mrSges(11,1) - t721 * mrSges(11,2) + Ifges(11,3) * t717 / 0.2e1) * t717 + (t726 * mrSges(4,1) - t725 * mrSges(4,2) + Ifges(4,3) * t716 / 0.2e1) * t716 + (t712 * mrSges(10,1) - t713 * mrSges(10,2) + Ifges(10,3) * t714 / 0.2e1) * t714 + ((-t819 * mrSges(9,1) / 0.2e1 + t821 * mrSges(9,2) / 0.2e1) * t1294 * t1325 + Ifges(9,3) * t767 / 0.2e1) * t767 + (Ifges(2,3) / 0.2e1 + (m(3) * (t804 / 0.4e1 + t1096 / 0.4e1) * t856 * t1068 + m(9) * (t818 / 0.4e1 + t1097 / 0.4e1) * t878 / pkin(5) ^ 2) * t1013) * qJD(1) ^ 2;
T = t1;
