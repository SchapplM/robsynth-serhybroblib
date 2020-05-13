% Calculate vector of inverse dynamics joint torques with newton euler and ic for
% palh2m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
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
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh2m2IC_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2IC_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2IC_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'palh2m2IC_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2IC_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2IC_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2IC_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2IC_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2IC_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 06:51:38
% EndTime: 2020-05-03 06:52:55
% DurationCPUTime: 75.87s
% Computational Cost: add. (15144->1343), mult. (24588->1630), div. (0->0), fcn. (12748->12), ass. (0->744)
t585 = sin(qJ(6));
t591 = cos(qJ(6));
t451 = mrSges(7,1) * t591 - mrSges(7,2) * t585;
t410 = pkin(3) * m(7) + mrSges(6,1) + t451;
t592 = cos(qJ(5));
t368 = t410 * t592;
t597 = m(6) + m(7);
t583 = mrSges(7,3) + mrSges(6,2);
t586 = sin(qJ(5));
t980 = t583 * t586;
t431 = -pkin(5) * t597 - mrSges(5,1) + t980;
t255 = t431 - t368;
t367 = t586 * t410;
t495 = t583 * t592;
t284 = t495 + t367;
t276 = mrSges(5,2) + t284;
t587 = sin(qJ(4));
t593 = cos(qJ(4));
t1224 = -t255 * t587 + t276 * t593;
t1227 = t255 * t593;
t1230 = -t276 * t587 - t1227;
t594 = cos(qJ(3));
t1241 = t1230 * t594;
t588 = sin(qJ(3));
t1243 = t1224 * t588 - t1241;
t346 = t367 + mrSges(5,2);
t1189 = t346 * pkin(4);
t1202 = t1189 * t588;
t1203 = pkin(2) * t431;
t382 = pkin(2) * t410;
t1076 = pkin(4) * t588;
t882 = t583 * t1076;
t293 = -t382 + t882;
t266 = t293 * t592;
t106 = -t266 - t1203 - t1202;
t390 = pkin(5) * t410;
t1103 = 0.2e1 * t390;
t541 = pkin(2) * t583;
t885 = t410 * t1076;
t288 = t541 + t885;
t263 = t288 * t587;
t574 = t591 ^ 2;
t582 = Ifges(7,1) - Ifges(7,2);
t981 = t582 * t574;
t1041 = pkin(3) * mrSges(7,1);
t1060 = Ifges(7,4) * t585;
t990 = (t1041 + t1060) * t591;
t1129 = t981 - 0.2e1 * t990;
t915 = pkin(3) * t585;
t501 = mrSges(7,2) * t915;
t486 = -0.2e1 * t501;
t581 = pkin(3) ^ 2 * m(7);
t724 = Ifges(7,1) + Ifges(6,3) + t486 + t581 - t1129;
t481 = pkin(5) * t980;
t660 = pkin(5) ^ 2;
t556 = t597 * t660;
t928 = -0.2e1 * t481 + t556;
t695 = Ifges(5,3) + t724 + t928;
t1184 = pkin(2) * t346;
t1201 = t431 * t1076;
t931 = t1201 - t1184;
t1242 = pkin(4) * t1241 + t106 * t593 + t931 * t587 + t695 + (-t263 + t1103) * t592;
t570 = m(5) + t597;
t1228 = m(4) + t570;
t1239 = pkin(4) * t1228 + mrSges(3,1);
t1150 = (pkin(3) * mrSges(7,2) - t582 * t585) * t591;
t502 = mrSges(7,1) * t915;
t1196 = -0.2e1 * Ifges(7,4) * t574 + Ifges(7,4) + t1150 + t502;
t271 = Ifges(6,5) + t1196;
t1058 = Ifges(7,6) * t591;
t1059 = Ifges(7,5) * t585;
t452 = t1058 + t1059;
t443 = -Ifges(6,6) + t452;
t405 = t443 * t592;
t154 = t271 * t586 - t405;
t151 = Ifges(5,6) + t154;
t1030 = t151 * t587;
t1080 = mrSges(6,3) * pkin(5);
t404 = t586 * t443;
t153 = t271 * t592 + t404;
t1064 = mrSges(7,1) * t585;
t518 = pkin(5) * t1064;
t1062 = mrSges(7,2) * t591;
t519 = pkin(5) * t1062;
t1145 = Ifges(5,5) + t518 + t519 - t1080 + t153;
t66 = t1145 * t593 - t1030;
t778 = t346 + t495;
t1115 = t587 * t778 + t1227;
t675 = t1103 * t592 + t695;
t1240 = -pkin(2) * t1115 + t675;
t105 = mrSges(4,2) + t1224;
t479 = pkin(2) * t570 + mrSges(4,1);
t99 = t479 - t1115;
t96 = t99 * t594;
t46 = -t105 * t588 + t96;
t1238 = t46 + t1239;
t1232 = t1228 * pkin(4) ^ 2;
t1237 = Ifges(3,2) - Ifges(3,1) + t1232;
t1008 = t288 * t592;
t1128 = -t931 + t1008;
t1236 = t587 * t106 + t1128 * t593;
t1231 = pkin(2) * t1224;
t920 = pkin(2) * t587;
t1207 = t431 * t920;
t1167 = -0.2e1 * t1207;
t865 = Ifges(6,2) + Ifges(7,3) + t581;
t1117 = -Ifges(6,1) - Ifges(7,2) + t865 - t1129;
t1216 = 0.2e1 * t501 - t1117;
t221 = 0.4e1 * t1216;
t193 = t221 * t586;
t549 = pkin(5) * t583;
t158 = t193 - 0.4e1 * t549;
t810 = t410 * t920;
t333 = 0.2e1 * t810;
t453 = Ifges(7,5) * t591 - Ifges(7,6) * t585;
t1113 = -mrSges(7,3) * pkin(3) + Ifges(6,4) + t453;
t575 = t592 ^ 2;
t1002 = t1113 * t575;
t1165 = 0.2e1 * t1002;
t1109 = t1165 - t1113;
t732 = 0.4e1 * t1109;
t879 = pkin(5) * t367;
t1134 = -0.4e1 * t879 + (4 * Ifges(5,4)) + t732;
t576 = t593 ^ 2;
t736 = 0.8e1 * t1109;
t706 = 0.8e1 * t879 - (8 * Ifges(5,4)) - t736;
t220 = 0.8e1 * t1216;
t185 = t220 * t586;
t946 = -t185 + 0.8e1 * t549;
t697 = (t946 * t592 + t706) * t576 + t1134;
t1226 = t1167 + (t333 + t158) * t592 + t697;
t589 = sin(qJ(2));
t1075 = pkin(4) * t589;
t1204 = pkin(1) * t431;
t311 = t589 * t1189;
t1218 = -t311 - t1204;
t1225 = -t1075 * t495 + t1204 + t1218;
t503 = pkin(2) * t1062;
t504 = pkin(2) * t1064;
t584 = -mrSges(5,3) - mrSges(6,3);
t540 = t584 * pkin(2);
t54 = Ifges(4,5) + t503 + t504 + t540 + t66;
t68 = t1145 * t587 + t151 * t593;
t65 = Ifges(4,6) + t68;
t31 = t54 * t588 + t65 * t594;
t1188 = m(3) + t1228;
t101 = t105 * t594;
t444 = t479 * t588;
t429 = t444 + mrSges(3,2);
t1195 = t1115 * t588 - t101 - t429;
t595 = cos(qJ(2));
t1223 = t1188 * pkin(1) + t1195 * t589 + t1238 * t595 + mrSges(2,1);
t753 = -0.2e1 * t1196;
t1222 = -0.2e1 * t1203;
t1221 = 0.2e1 * t1203;
t1011 = t1196 * t575;
t363 = t1113 * t586;
t330 = -0.16e2 * t363;
t225 = -t330 + 0.8e1 * t390;
t201 = t225 * t587;
t314 = 0.2e1 * t1184;
t508 = 0.2e1 * t541;
t1220 = (t201 + t508) * t592 + t314;
t223 = 0.2e1 * t1216;
t1208 = t223 * t586;
t1219 = t1208 - t549;
t690 = -(4 * Ifges(4,4)) + t697;
t313 = 0.4e1 * t1184;
t671 = -(4 * Ifges(5,1)) + (4 * Ifges(5,2)) - 0.8e1 * t481 + 0.4e1 * t556 + t221;
t126 = t671 * t587;
t601 = Ifges(7,2) / 0.2e1;
t553 = t601 - Ifges(7,1) / 0.2e1;
t1107 = t501 + Ifges(6,1) / 0.2e1 - Ifges(6,2) / 0.2e1 - t581 / 0.2e1 - t553 * t574 - t990;
t598 = Ifges(7,3) / 0.2e1;
t218 = -t1107 + t598 - Ifges(7,2) / 0.2e1;
t987 = t575 * t587;
t852 = t218 * t987;
t159 = 0.16e2 * t852;
t958 = t126 + t159;
t824 = t313 + t958;
t942 = t201 + 0.4e1 * t541;
t1217 = (0.4e1 * t810 + t158) * t592 + (t592 * t942 + t824) * t593 + t690 - 0.4e1 * t1207;
t209 = -t1203 + t1202;
t1213 = -0.2e1 * t1184;
t509 = -0.2e1 * t541;
t621 = pkin(2) ^ 2;
t988 = t570 * t621;
t1158 = Ifges(4,3) + t988;
t683 = t695 + t1158;
t1063 = mrSges(4,2) * t588;
t903 = pkin(4) * t1063;
t95 = pkin(4) * t96;
t1215 = ((t509 - t885) * t587 + t1103) * t592 + t683 + t95 - t903 + (t1201 + t1213) * t587;
t919 = pkin(2) * t588;
t807 = t431 * t919;
t1214 = 0.4e1 * t807;
t1212 = pkin(4) * t431;
t1210 = -t46 + t1238;
t47 = t588 * t99 + t101;
t1209 = -t1195 - t47;
t197 = t223 * t575;
t1159 = -Ifges(4,2) - t988;
t1206 = Ifges(4,1) + t1159;
t419 = t452 * t586;
t454 = t1062 + t1064;
t440 = pkin(5) * t454;
t270 = 0.2e1 * t419 + t440;
t247 = t270 * t592;
t1205 = t247 - t1196;
t961 = t592 * t593;
t977 = t586 * t587;
t1131 = t961 - t977;
t997 = t1131 * t594;
t1162 = Ifges(5,2) - Ifges(5,1);
t1198 = -t556 - t1162;
t1199 = -0.2e1 * t1198 + t223;
t543 = pkin(1) * t589;
t1141 = t543 - 0.2e1 * t919;
t1001 = t410 * t1141;
t179 = t218 * t586;
t150 = t179 + t549 / 0.2e1;
t972 = t587 * t588;
t853 = t150 * t972;
t1197 = 0.8e1 * t853 + t1001;
t1088 = -0.2e1 * t589;
t890 = 0.2e1 * t977;
t1192 = -0.2e1 * pkin(1);
t1191 = -2 * qJD(2);
t966 = t589 * t150;
t1190 = -0.16e2 * t966;
t799 = t346 * t920;
t812 = t588 * t543;
t450 = -pkin(2) + t812;
t995 = t450 * t346;
t342 = -0.4e1 * t363;
t227 = -t342 + t1103;
t205 = t227 * t587;
t994 = t450 * t583;
t1186 = t205 - t994;
t752 = 0.2e1 * t1196;
t1185 = t752 * t575 + t1205;
t1050 = t65 * t588;
t30 = t54 * t594 - t1050;
t322 = pkin(1) * t346;
t1111 = t486 + t1117;
t670 = t586 * t1111;
t828 = t591 * t1041;
t416 = Ifges(7,3) - t501 + t828;
t1000 = t416 * t586;
t278 = pkin(5) * t451 + t453 * t586;
t173 = t416 * t592 + t278;
t421 = t453 * t592;
t1180 = (-t421 + t1000) * t593 + t173 * t587;
t599 = -Ifges(7,3) / 0.2e1;
t816 = t981 + t553;
t907 = -0.2e1 * t1060;
t250 = t599 + (t907 - t1041) * t591 + t501 + t816;
t119 = t250 * t592 - t278;
t144 = t250 * t586 + t421;
t1179 = t119 * t587 + t144 * t593;
t590 = sin(qJ(1));
t596 = cos(qJ(1));
t1121 = g(1) * t596 + g(2) * t590;
t971 = t587 * t592;
t975 = t586 * t593;
t1130 = -t975 - t971;
t998 = t1130 * t594;
t234 = t1131 * t588 - t998;
t728 = -t1130 * t588 - t997;
t1176 = t589 * t234 + t595 * t728;
t341 = -0.4e1 * t1002;
t755 = 0.2e1 * t1113;
t735 = t341 + t755;
t1135 = 0.2e1 * t879 - (2 * Ifges(5,4)) + t735;
t943 = -t1208 + 0.2e1 * t549;
t1175 = (t943 * t592 + t1135) * t576;
t1074 = pkin(4) * t594;
t1174 = t1074 * t1224 + t1236;
t550 = pkin(4) * t583;
t1173 = (t410 * t543 - t550) * t588 - t382;
t813 = t431 * t543;
t1172 = (-t813 - t1189) * t588 + t1203;
t1005 = t293 * t587;
t1143 = -t390 + t263;
t1028 = t1143 * t586;
t1123 = -t293 * t586 + t1008;
t283 = t368 - t980;
t876 = (t283 * t587 + t284 * t593) * t1074;
t1171 = t1123 * t593 - (-t549 + t1005) * t592 - t1028 + t876;
t139 = t158 * t592;
t1170 = (t139 + t1134) * t576 + t1135;
t714 = -t882 - t382;
t781 = -t714 * t592 + t209;
t1168 = 0.4e1 * t781;
t577 = t594 ^ 2;
t982 = t577 * t589;
t1166 = -0.8e1 * t982;
t1160 = t1175 - Ifges(4,4);
t339 = -0.8e1 * t363;
t226 = t339 - 0.4e1 * t390;
t1154 = t226 * t587;
t391 = pkin(4) * t410;
t1101 = -0.2e1 * t391;
t805 = t583 * t919;
t300 = t1101 + 0.4e1 * t805;
t267 = t300 * t587;
t1153 = t267 * t592;
t939 = t1216 * t586 - t549;
t1152 = t592 * t939;
t976 = t586 * t592;
t847 = t1196 * t976;
t763 = t587 * t847;
t836 = t452 * t987;
t1072 = pkin(5) * t586;
t877 = t454 * t1072;
t302 = t452 + t877;
t973 = t587 * t302;
t1151 = (0.8e1 * t763 - 0.8e1 * t836 + 0.4e1 * t973) * t593;
t236 = -Ifges(5,4) + t1113 + t879;
t899 = -0.2e1 * t1002;
t688 = -t236 - t899;
t803 = t586 * t920;
t366 = t454 * t803;
t269 = -0.2e1 * t440 - 0.4e1 * t419;
t754 = -0.4e1 * t1196;
t720 = t269 * t592 + t754 * t575 + t752;
t922 = pkin(2) * t454;
t974 = t587 * t1196;
t1146 = t366 + t720 * t576 + ((0.4e1 * t586 * t974 - t922) * t592 - 0.4e1 * t836 + 0.2e1 * t973) * t593 + t1185;
t343 = -0.2e1 * t363;
t229 = t343 - t390;
t301 = -0.2e1 * t805 + t391;
t1085 = 0.2e1 * t592;
t992 = t452 * t575;
t1040 = (-0.2e1 * t992 + (t1085 * t1196 + t440) * t586 + t452) * t576;
t1127 = t1040 + t992;
t1139 = -0.2e1 * t847 + 0.2e1 * t1127 - t302;
t917 = pkin(2) * t592;
t921 = pkin(2) * t586;
t1138 = -pkin(2) * t1130 - (t1076 * t592 + t921) * t593 - t587 * t917;
t991 = t454 * t592;
t1137 = -pkin(3) * t991 - t419 - t440;
t375 = t431 * t1075;
t909 = 0.2e1 * t1212;
t1136 = -t589 * t909 + t375;
t1126 = (0.2e1 * Ifges(6,5) + t752) * t592 + 0.2e1 * t404;
t1124 = (qJD(4) + qJD(3));
t1120 = 0.2e1 * t1206;
t1116 = (-t1141 * t583 - t391) * t592 - t346 * t1141 + t1212;
t560 = m(6) / 0.2e1 + m(7) / 0.2e1;
t1110 = t481 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1 - t560 * t660;
t1105 = -0.8e1 * t236;
t1104 = 0.2e1 * t382;
t446 = t543 - t919;
t1098 = -0.2e1 * t446;
t1097 = 0.2e1 * t446;
t1096 = -0.2e1 * t450;
t1095 = 0.2e1 * t550;
t1094 = -0.2e1 * t586;
t1093 = 0.2e1 * t586;
t1092 = -0.2e1 * t587;
t1091 = 0.2e1 * t587;
t1090 = 0.8e1 * t587;
t1089 = 0.4e1 * t588;
t1087 = 0.2e1 * t589;
t1086 = -0.2e1 * t592;
t1083 = -2 * qJD(6);
t1081 = mrSges(4,2) * pkin(4);
t1077 = pkin(4) * t454;
t1073 = pkin(5) * t284;
t1067 = t479 * pkin(4);
t1029 = t154 * t587;
t82 = t153 * t587 + t154 * t593;
t1054 = t588 * t82;
t81 = t153 * t593 - t1029;
t36 = t588 * t81 + t594 * t82;
t1057 = t36 * t589;
t98 = t234 * t595 - t589 * t728;
t83 = (t591 * t907 + t598 + t816) * t98 * t1083;
t1066 = (((t1126 * t593 - 0.2e1 * t1029) * t594 - 0.2e1 * t1054) * t595 - 0.2e1 * t1057) * qJD(5) + t83;
t968 = t588 * t589;
t384 = (t594 * t595 - t968) * qJD(1);
t963 = t589 * t594;
t385 = (t588 * t595 + t963) * qJD(1);
t199 = t384 * t593 - t385 * t587;
t200 = t384 * t587 + t385 * t593;
t619 = qJD(1) ^ 2;
t422 = -pkin(1) * t619 - t1121;
t962 = t589 * t595;
t238 = -g(3) * t595 - t422 * t589 + (t619 * t962 + qJDD(2)) * pkin(4);
t578 = t595 ^ 2;
t618 = qJD(2) ^ 2;
t252 = -g(3) * t589 + t422 * t595 + (-t578 * t619 - t618) * pkin(4);
t886 = qJDD(2) + qJDD(3);
t69 = t594 * t238 - t588 * t252 + (t384 * t385 + t886) * pkin(2);
t770 = qJDD(4) + t886;
t569 = qJD(2) + qJD(3);
t85 = t238 * t588 + t252 * t594 + (-t384 ^ 2 - (t569 ^ 2)) * pkin(2);
t18 = -t587 * t85 + t593 * t69 + (t199 * t200 + t770) * pkin(5);
t538 = qJD(2) + t1124;
t29 = t587 * t69 + t593 * t85 + (-t199 ^ 2 - (t538 ^ 2)) * pkin(5);
t1065 = t586 * t18 + t592 * t29;
t491 = qJD(5) + t538;
t913 = qJD(1) * t491;
t43 = -qJDD(1) * t1176 - t98 * t913;
t1021 = t229 * t592;
t127 = t218 + t1110;
t217 = t1216 * t575;
t58 = t127 + t217 + t1021;
t1056 = t58 * t576;
t1055 = t588 * t68;
t1053 = t589 * t31;
t34 = t588 * t66 + t594 * t68;
t1052 = t589 * t34;
t228 = -t342 + t390;
t1022 = t228 * t592;
t1051 = (-t197 + t1216 - t481 + t1022) * t576;
t1049 = t1185 * t576;
t986 = t575 * t588;
t851 = t218 * t986;
t161 = 0.16e2 * t851;
t230 = t390 + 0.2e1 * t363;
t1020 = t230 * t592;
t219 = t601 + t599 + t1107;
t59 = t1111 * t575 + t1020 - t1110 + t219;
t984 = t576 * t588;
t861 = t59 * t984;
t1047 = t161 - 0.16e2 * t861;
t165 = -0.8e1 * t851;
t1046 = t165 + 0.8e1 * t861;
t985 = t575 * t589;
t850 = t218 * t985;
t167 = -0.8e1 * t850;
t983 = t576 * t589;
t1044 = 0.8e1 * t59 * t983 + t167;
t676 = t1216 + t928 + t1162;
t1043 = t217 + (t227 * t592 - t197 + t676) * t576;
t829 = t1113 * t985;
t759 = t588 * t829;
t305 = 0.16e2 * t759;
t171 = t549 + t670;
t152 = t171 * t592;
t80 = t152 - t688;
t841 = t576 * t968;
t771 = t80 * t841;
t1042 = t305 + 0.16e2 * t771;
t1036 = t1224 * t594;
t1033 = t127 * t589;
t1031 = t150 * t588;
t1027 = t171 * t587;
t210 = -t1184 - t1201;
t1026 = t210 * t589;
t811 = t346 * t919;
t215 = 0.4e1 * t811 + t909;
t1025 = t215 * t587;
t1024 = t218 * t575;
t1023 = t219 * t575;
t231 = Ifges(4,4) + t236;
t1019 = t231 * t588;
t1018 = t236 * t587;
t241 = t755 + t879;
t1017 = t241 * t587;
t1009 = t283 * t593;
t289 = -t541 + t885;
t1007 = t289 * t589;
t806 = t583 * t920;
t297 = t390 - t806;
t1004 = t297 * t592;
t1003 = t346 * t587;
t996 = t446 * t587;
t993 = t450 * t587;
t979 = t583 * t587;
t978 = t586 * t1196;
t970 = t588 * t127;
t969 = t588 * t230;
t967 = t588 * t593;
t160 = -0.32e2 * t852;
t960 = ((8 * Ifges(5,1)) - (8 * Ifges(5,2)) - t220 + 0.16e2 * t481 - 0.8e1 * t556) * t587 + t160;
t672 = 0.4e1 * t481 - t1199;
t125 = t672 * t587;
t163 = -0.8e1 * t852;
t959 = t125 + t163;
t128 = t587 * t1190;
t925 = pkin(1) * t588;
t796 = t410 * t925;
t331 = -0.2e1 * t796;
t957 = t128 + t331;
t130 = 0.32e2 * t853;
t809 = t410 * t919;
t956 = t130 - 0.8e1 * t809;
t134 = t676 * t587;
t954 = t134 + 0.4e1 * t852;
t222 = 0.4e1 * t1111;
t371 = -0.2e1 * t382;
t532 = 0.4e1 * t549;
t953 = (t222 * t586 + t532) * t587 + t371;
t472 = -0.4e1 * t481;
t952 = (t472 + t221) * t587 + t159;
t835 = t587 * t968;
t155 = 0.16e2 * t230 * t835;
t542 = pkin(1) * t583;
t511 = -0.2e1 * t542;
t951 = t155 + t511;
t834 = t1113 * t986;
t760 = t587 * t834;
t303 = -0.32e2 * t760;
t848 = t236 * t972;
t950 = 0.16e2 * t848 + t303;
t304 = 0.16e2 * t760;
t947 = -0.8e1 * t848 + t304;
t945 = -t185 + t532;
t534 = -0.2e1 * t549;
t944 = t193 + t534;
t941 = t1154 + t509;
t940 = t205 + t541;
t938 = -t1087 * t714 + t331;
t937 = t300 * t589 + t511;
t837 = t1113 * t987;
t324 = 0.4e1 * t837;
t935 = t324 - 0.2e1 * t1018;
t537 = -0.2e1 * t550;
t934 = -0.4e1 * t809 + t537;
t356 = -0.2e1 * t885;
t933 = t356 + t508;
t792 = -t920 / 0.2e1;
t932 = -t390 / 0.2e1 + t583 * t792;
t383 = pkin(1) * t410;
t918 = pkin(2) * t589;
t916 = pkin(2) * t593;
t914 = pkin(3) * t586;
t912 = (-t1063 + t1239) * pkin(1);
t455 = t479 * pkin(1);
t911 = qJD(1) * qJD(2);
t910 = qJD(1) * qJD(6);
t102 = (0.2e1 * Ifges(5,5) + 0.2e1 * t518 + 0.2e1 * t519 - 0.2e1 * t1080 + t1126) * t593;
t906 = (((t102 - 0.2e1 * t1030) * t594 - 0.2e1 * t1055) * t595 - 0.2e1 * t1052) * qJD(4) + t1066;
t905 = 0.4e1 * t587;
t904 = 0.2e1 * t593;
t901 = -0.4e1 * t1023;
t900 = 0.4e1 * t236 * t589;
t898 = -0.2e1 * t996;
t897 = 0.2e1 * t996;
t896 = 0.2e1 * t993;
t893 = 0.4e1 * t983;
t892 = -0.4e1 * t982;
t888 = -0.2e1 * t961;
t559 = m(6) / 0.4e1 + m(7) / 0.4e1;
t884 = t410 * t1075;
t880 = pkin(4) * t961;
t878 = pkin(5) * t368;
t873 = pkin(4) * t444;
t871 = -0.2e1 * t925;
t870 = 0.2e1 * t925;
t866 = t429 * t1192;
t107 = (m(5) / 0.2e1 + t560) * t621 - Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + t127;
t141 = (t1208 + t534) * t587;
t862 = (t107 - t799 + (-t806 + t229) * t592 + ((t141 + t382) * t592 - t1203 + t935) * t593 + t1043) * t982;
t864 = -0.8e1 * t862 + t1044;
t863 = qJD(5) * t1073;
t859 = t80 * t984;
t666 = t672 - t1120;
t858 = t666 * t588 + t1046;
t857 = (0.4e1 * Ifges(4,1) - 0.4e1 * Ifges(4,2) + t671 - 0.4e1 * t988) * t588 + t1047;
t856 = mrSges(4,2) * t1192 + t1042;
t855 = mrSges(4,2) * t543;
t149 = t179 + t549 / 0.4e1;
t854 = t149 * t971;
t849 = t219 * t987;
t846 = t410 * t977;
t845 = t1141 * t1003;
t844 = t1141 * t979;
t843 = t410 * t993;
t842 = t450 * t979;
t840 = t586 * t972;
t839 = t270 * t971;
t838 = t575 * t974;
t833 = t589 * t967;
t764 = t588 * t850;
t136 = 0.32e2 * t587 * t764;
t317 = -0.2e1 * t322;
t761 = t127 * t835;
t827 = -0.16e2 * t761 + t136 + t317;
t826 = t1213 + t959;
t825 = t125 + 0.8e1 * t849 - t1184;
t823 = 0.2e1 * t1201 + t958;
t822 = t1184 + t954;
t821 = 0.8e1 * t807 + t950;
t758 = t587 * t829;
t306 = 0.16e2 * t758;
t358 = t431 * t871;
t820 = -0.8e1 * t589 * t1018 + t306 - t358;
t355 = 0.2e1 * t885;
t819 = t355 + t941;
t818 = 0.4e1 * t1018 - 0.8e1 * t837 + t1221;
t815 = t587 * t383;
t804 = t410 * t921;
t802 = t588 * t918;
t800 = t592 * t916;
t795 = qJD(4) * t1231;
t794 = t589 * t911;
t793 = t922 / 0.2e1;
t789 = t592 * t18 - t29 * t586;
t91 = t199 * t586 + t200 * t592;
t79 = -t491 * t585 + t591 * t91;
t784 = t583 * t897;
t782 = t236 + t152;
t90 = t592 * t199 - t200 * t586;
t282 = t346 * t870;
t775 = 0.4e1 * t809;
t460 = t583 * t870;
t773 = t375 - t322;
t772 = t1105 + 0.16e2 * t1002;
t769 = 0.2e1 * t1081 + t858;
t765 = mrSges(4,2) * t812;
t762 = t592 * t835;
t364 = -0.2e1 * t1201;
t756 = t364 + t826;
t751 = t410 * t803;
t750 = t583 * t803;
t749 = (t1104 - t882) * t592 - t1202;
t745 = -t542 - 0.2e1 * t884;
t744 = -t542 - t884;
t318 = 0.2e1 * t1189;
t212 = t318 - 0.4e1 * t807;
t296 = t1095 + t775;
t739 = -pkin(1) + 0.2e1 * t802;
t738 = t587 * t914 - pkin(2);
t734 = -0.2e1 * pkin(5) * t991 + t753;
t733 = 0.16e2 * t1109;
t78 = -t491 * t591 - t585 * t91;
t729 = -t1137 * t587 + (-t452 * t592 + t454 * t914) * t593;
t725 = 0.2e1 * t778;
t723 = t744 * t588;
t719 = 0.4e1 * t879 - t736;
t716 = -0.2e1 * t879 + t732;
t457 = t481 / 0.2e1;
t156 = t457 + t218;
t713 = -0.2e1 * t1022 + 0.4e1 * t156 + 0.4e1 * t1051;
t712 = t1020 - t127 + 0.2e1 * t1024;
t711 = qJD(6) * t491;
t710 = -0.16e2 * t854 - 0.4e1 * t1017;
t708 = t212 + t947;
t707 = (t592 * t945 + t719) * t576 + t716;
t701 = -mrSges(4,3) + t454 + t584;
t178 = -0.16e2 * t670;
t700 = ((t178 - 0.16e2 * t549) * t592 - 0.16e2 * t879 + (16 * Ifges(5,4)) + t733) * t576 + t706;
t694 = t1130 * t1077;
t649 = 2 * Ifges(4,4);
t691 = t649 + t1170;
t689 = t284 * t916 + (t549 + t810) * t592 + t297 * t586;
t687 = t712 * t905;
t327 = -0.32e2 * t363;
t198 = (t327 - 0.16e2 * t390) * t587;
t685 = (((t198 - 0.8e1 * t541) * t592 - 0.8e1 * t1184 + t960) * t593 + (-0.8e1 * t810 + t946) * t592 + 0.8e1 * t1207 + (8 * Ifges(4,4)) + t700) * t577 + t690;
t260 = 0.2e1 * t843;
t407 = 0.2e1 * t994;
t682 = -t431 * t896 + ((t1154 + t407) * t592 + 0.2e1 * t995 + t959) * t593 + (t260 + t943) * t592 + t1170;
t677 = t688 + t1160;
t14 = t1217 * t577;
t674 = t14 + 0.2e1 * t873 - (2 * Ifges(3,4)) + t691;
t673 = t472 + t1199;
t667 = t676 + t1206;
t194 = t221 * t575;
t664 = -t197 + (t226 * t592 + t194 + t672) * t576 + t667;
t663 = t431 * t993 + (-t843 + t939) * t592 + t688;
t662 = (-t1060 / 0.2e1 - t1041 / 0.2e1) * t591 + (-Ifges(7,2) / 0.4e1 + Ifges(7,1) / 0.4e1) * t574 + t559 * t660 - Ifges(5,1) / 0.4e1 + Ifges(6,1) / 0.4e1 + Ifges(5,2) / 0.4e1 - Ifges(6,2) / 0.4e1 + Ifges(7,2) / 0.4e1 - Ifges(7,3) / 0.4e1 + t501 / 0.2e1 - t581 / 0.4e1;
t617 = qJD(3) ^ 2;
t616 = qJD(4) ^ 2;
t615 = qJD(5) ^ 2;
t614 = qJD(6) ^ 2;
t573 = -0.2e1 * t1081;
t536 = -0.8e1 * t549;
t510 = -0.4e1 * t541;
t497 = -0.2e1 * t903;
t480 = qJDD(5) + t770;
t470 = 0.2e1 * t481;
t464 = -0.4e1 * t805;
t462 = 0.8e1 * t805;
t441 = pkin(4) * t451;
t434 = qJDD(1) * t595 - t794;
t427 = t455 * t1088;
t424 = t587 * t594 + t967;
t423 = -t593 * t594 + t972;
t396 = 0.2e1 * t1204;
t393 = -0.4e1 * t1203;
t392 = mrSges(2,2) - mrSges(3,3) + t701;
t378 = -0.2e1 * t390;
t376 = 0.4e1 * t390;
t372 = -0.2e1 * t383;
t370 = 0.4e1 * t382;
t362 = t382 / 0.2e1;
t359 = -0.8e1 * t807;
t338 = t840 * t1077;
t335 = 0.8e1 * t809;
t334 = -0.4e1 * t810;
t319 = -0.2e1 * t1189;
t316 = -0.4e1 * t1184;
t299 = t462 - 0.4e1 * t391;
t298 = t462 + t1101;
t291 = t370 + 0.2e1 * t882;
t286 = -0.2e1 * t1202;
t281 = -0.4e1 * t811;
t280 = 0.8e1 * t811;
t275 = t410 * t1098;
t256 = t714 * t1091;
t214 = t280 + 0.4e1 * t1212;
t213 = t280 + t909;
t208 = t393 + 0.2e1 * t1202;
t202 = (-t330 + t376) * t587;
t192 = 0.8e1 * t1019;
t191 = -0.4e1 * t589 * t230;
t189 = -0.4e1 * t969;
t187 = 0.8e1 * t969;
t177 = (t460 + 0.2e1 * t1007) * t587;
t175 = t215 * t589;
t169 = t209 * t1087;
t146 = (t296 * t589 + t372) * t587;
t131 = -0.16e2 * t853;
t129 = t588 * t1190;
t122 = 0.4e1 * t970;
t89 = qJD(6) + t90;
t84 = -t550 + t1197;
t67 = t119 * t593 - t144 * t587;
t56 = t1230 * t588 + t1036;
t42 = t98 * qJDD(1) - t1176 * t913;
t40 = qJDD(6) + t43;
t35 = t594 * t81 - t1054;
t33 = t594 * t66 - t1055;
t27 = Ifges(3,6) + t31;
t25 = qJD(6) * t78 + t591 * t42 - t585 * t480;
t24 = -qJD(6) * t79 - t585 * t42 - t591 * t480;
t23 = t701 * pkin(4) + Ifges(3,5) + t30;
t22 = Ifges(7,1) * t79 + Ifges(7,4) * t78 + Ifges(7,5) * t89;
t21 = Ifges(7,4) * t79 + Ifges(7,2) * t78 + Ifges(7,6) * t89;
t20 = Ifges(7,5) * t79 + Ifges(7,6) * t78 + Ifges(7,3) * t89;
t9 = t33 * t589 + t34 * t595;
t7 = -qJDD(1) * pkin(1) - t590 * g(1) + t596 * g(2) + (-(-t423 * t595 - t424 * t589) * qJDD(1) + (-(t423 * t589 - t424 * t595) * qJD(1) + t200) * t538) * pkin(5) + (-t434 + t794) * pkin(4) + (t491 * t91 - t43) * pkin(3) + (t588 * (qJDD(1) * t589 + t595 * t911) - t594 * t434 + (qJD(3) + t569) * t385) * pkin(2);
t6 = t30 * t589 + t31 * t595;
t5 = t23 * t589 + t27 * t595;
t4 = (-t491 ^ 2 - t90 ^ 2) * pkin(3) + t1065;
t3 = (t90 * t91 + t480) * pkin(3) + t789;
t2 = t4 * t591 - t585 * t7;
t1 = -t4 * t585 - t591 * t7;
t8 = [(t35 * t589 + t36 * t595) * qJDD(5) + (t33 * t595 - t1052) * t616 + (((t1137 * t593 - t452 * t971 + t738 * t454) * t594 - t1077 + t729 * t588) * t595 + t729 * t963 - t1137 * t833 + t452 * t762 - (t738 * t968 + pkin(1)) * t454) * t614 + (t35 * t595 - t1057) * t615 - (-t1223 * t590 - t392 * t596) * g(1) + ((((((t178 + t536) * t592 - 0.8e1 * t879 + t733) * t576 + (t160 + ((t327 - 0.8e1 * t390) * t587 + t510) * t592 + (0.8e1 * t481 - t220) * t587 - 0.4e1 * t804) * t593 + (t334 + t945) * t592 + 0.4e1 * t750 + t719) * t577 + ((t1094 * t301 + t592 * t934 + t303) * t593 + t161 + t1153 + t550 * t890) * t594 + ((t202 + t933) * t592 - t714 * t1093 + t952) * t593 + (-t256 + t944) * t592 + t289 * t890 + t707) * t578 + (((t1094 * t383 + t592 * t937 + t136) * t593 + t305 + t146 * t592 + t542 * t890) * t594 + (t592 * t938 + t306) * t593 + t167 + t177 * t592) * t595 + (((t202 + t508) * t592 + 0.2e1 * t804 + t952) * t593 + (t333 + t944) * t592 - 0.2e1 * t750 + t707) * t577 + ((t1097 * t980 + t275 * t592 + t304) * t593 + t165 + t592 * t784 + t846 * t1097) * t594 + (t592 * t944 + t716) * t576 + (t163 + ((t339 + t378) * t587 + t407) * t592 + (t470 - t223) * t587 + 0.2e1 * t450 * t367) * t593 + (t260 + 0.4e1 * t179) * t592 + t842 * t1094 + t735 + (0.2e1 * (t542 * t975 + t586 * t815) * t595 + ((-0.8e1 * t1051 + (0.32e2 * t854 + 0.8e1 * t1017) * t593 + 0.4e1 * t1022 + 0.4e1 * t751 - 0.8e1 * t156) * t578 + t593 * t710 + t713) * t594) * t588 + (-0.8e1 * (t1051 + (t324 + (t1219 * t587 + t362) * t592 - t1017 - t583 * t921 / 0.2e1) * t593 + t217 + (t343 + t932) * t592 - t751 / 0.2e1 + t156) * t577 + 0.2e1 * (t301 * t977 + t550 * t975) * t594 + (t1093 * t289 + t710) * t593 + t714 * t890 + t713 + (0.8e1 * (-t1219 * t592 + t241 + t341) * t576 + (-0.16e2 * t156 * t587 + 0.8e1 * t228 * t971 + 0.4e1 * t804) * t593 - 0.16e2 * t149 * t592 - 0.4e1 * t241) * t594 * t588) * t962) * qJD(5) + (((((0.8e1 * t1150 + 0.8e1 * t502 + (-0.16e2 * t574 + 0.8e1) * Ifges(7,4)) * t575 + (0.4e1 * t440 + 0.8e1 * t419) * t592 + t754) * t576 - t1151 - 0.2e1 * t366 + t720) * t577 + t1146) * t578 + t1146 * t577 + t1049 + (-0.2e1 * t763 + 0.2e1 * t836 - t973) * t593 - t1011 - t452 * t976 + ((0.2e1 * t577 * t800 + (-pkin(4) * t977 + t880) * t594) * t578 + (-t739 * t997 + t1130 * (-t918 + t925)) * t595 + pkin(3) - t1131 * t450 + t446 * t998) * t454 + (((-0.4e1 * t1040 + (-t454 * t921 - t587 * t753 - 0.4e1 * t838 - 0.2e1 * t839) * t904 - 0.4e1 * t992 + (t454 * t920 - 0.2e1 * t978) * t1086 + 0.2e1 * t302) * t594 + t694) * t578 + ((0.2e1 * t838 + t839 - t974) * t904 + t1139) * t594) * t588 + (-0.4e1 * ((t586 * t793 + (0.2e1 * t1011 + t1205) * t587) * t593 + (t587 * t793 - t978) * t592 - t877 / 0.2e1 - t1059 / 0.2e1 - t1058 / 0.2e1 + t1127) * t577 + t694 * t594 + (-t269 * t971 + 0.4e1 * t838 - 0.2e1 * t974) * t593 + t338 + t1139 + ((0.4e1 * t1011 - 0.4e1 * t1049 + 0.2e1 * t247 + t753 + t1151) * t594 - t454 * t880) * t588) * t962) * t1083) * qJD(1) + (((pkin(2) * t451 + t173 * t593 - t416 * t977 + t453 * t971) * t594 + t441 - t1180 * t588) * t595 - t1180 * t963 - t173 * t833 - t453 * t762 + t835 * t1000 + (pkin(1) - t802) * t451) * qJDD(6) + t9 * qJDD(4) - (t1223 * t596 - t392 * t590) * g(2) + ((((t1111 + t470 + (t673 + (t376 + 0.8e1 * t363) * t592 + t222 * t575) * t576 - t799 + (t382 * t592 - t1203 + (0.8e1 * t1002 - 0.4e1 * t236 + t139) * t587) * t593 + (t378 - t806 + t342) * t592 + t197 + t1198) * t892 + (((t155 + t937) * t592 + t175 + t827) * t593 + (t146 + t129) * t592 + (t212 * t589 + t396) * t587 + t968 * t1105 + t1042) * t594 + ((t128 + t938) * t592 + t169 + t820) * t593 + (t177 + t191) * t592 + (t282 + 0.2e1 * t1026) * t587 + 0.4e1 * t1033 + t1044) * t595 + t682 + ((((t198 + t510) * t592 + t316 + t960) * t593 + (t334 + t946) * t592 + 0.4e1 * t1207 + t700) * t577 + (((t130 + t934) * t592 + t1214 + t319 + t950) * t593 + (t267 + t187) * t592 + t1025 - 0.8e1 * t970 + t1047) * t594 + ((t201 + t933) * t592 + t314 + t823) * t593 + (-t256 + t158) * t592 + t209 * t1091 + t697) * t578 + (((t131 + t275) * t592 - t431 * t1098 + t947) * t593 + (t189 + t784) * t592 + t346 * t897 + t122 + t1046) * t594 + ((t958 + t1220) * t593 + t1226) * t577) * qJD(1) + t1066) * qJD(4) + t6 * qJDD(3) + (t486 + t865 + t589 * t866 + t346 * t896 + ((t1096 * t410 + t141) * t592 - t431 * t1096 + t935) * t593 + (0.4e1 * t859 + (-t446 * t725 + t588 * t687) * t593 + 0.4e1 * t834 + (t410 * t898 - 0.4e1 * t1031) * t592 - t431 * t898 - 0.2e1 * t855 - 0.2e1 * t1019) * t594 + (t497 + (-0.8e1 * t859 + ((0.2e1 * t391 + t464) * t592 + t281 - 0.2e1 * t1212 - 0.8e1 * t712 * t972) * t593 - 0.8e1 * t834 + (-t296 * t587 + 0.8e1 * t1031) * t592 - t212 * t587 + 0.4e1 * t1019 + 0.2e1 * t1067) * t594 + (t194 + t666 - 0.4e1 * t799 + (t370 * t592 + t393 + ((t536 - 0.8e1 * t670) * t592 + t772) * t587) * t593 + (-0.4e1 * t806 + t226) * t592 + (-t220 * t575 + t225 * t592 + t671) * t576) * t577 + t210 * t1092 + (t1092 * t289 + t227) * t592 + ((-0.2e1 * t882 + t953) * t592 + t286 + t818) * t593 + t664 + t1237) * t578 + (-t363 + t842) * t1085 + (((t940 * t592 + t822) * t593 + (t810 + t939) * t592 - t1207 + t677) * t892 + (-0.8e1 * t59 * t841 + (t383 * t1085 - 0.2e1 * t1204 + (-t296 * t592 - t212 + (0.16e2 * t150 * t592 - t772) * t972) * t589) * t593 + 0.8e1 * t764 + (0.4e1 * t230 * t968 + t587 * t937) * t592 + (t175 + t317) * t587 + ((t673 + t1120) * t588 + t573) * t589 + 0.2e1 * t455) * t594 + t80 * t893 + (-t725 * t925 + (t1086 * t289 - 0.2e1 * t210 + t687) * t589) * t593 + 0.4e1 * t829 + (t587 * t938 - 0.4e1 * t966) * t592 + (t169 - t358) * t587 + (Ifges(3,4) - t231 - t873) * t1087 + 0.2e1 * t912) * t595 + Ifges(3,1) + Ifges(5,1) + Ifges(2,3) + (0.2e1 * t799 + (t592 * t953 + t818) * t593 + (0.2e1 * t806 + t227) * t592 + t664) * t577 + 0.2e1 * t828 + t1043 - t1159 + t1188 * pkin(1) ^ 2) * qJDD(1) + (t30 * t595 - t1053) * t617 + qJD(5) * t83 + (t23 * t595 - t27 * t589) * t618 + t5 * qJDD(2) + ((((0.2e1 * t540 + 0.2e1 * Ifges(4,5) + t102 + 0.2e1 * t503 + 0.2e1 * t504 + (0.2e1 * t405 + (-0.2e1 * Ifges(6,5) + t753) * t586 - 0.2e1 * Ifges(5,6)) * t587) * t594 - 0.2e1 * t1050) * t595 - 0.2e1 * t1053) * qJD(3) + (t674 + t912 * t1088 + ((((t299 * t589 + t951) * t592 + t214 * t589 + t827) * t593 + (((t335 + 0.4e1 * t550) * t589 + t372) * t587 + t129) * t592 + ((t359 + 0.4e1 * t1189) * t589 + t396) * t587 + (-t192 - 0.4e1 * t1067) * t589 + t856) * t594 + (t589 * t1168 + t592 * t957 + t820) * t593 + ((t460 + 0.4e1 * t1007) * t587 + t191) * t592 + (t282 + 0.4e1 * t1026) * t587 + (0.2e1 * Ifges(3,1) - 0.2e1 * Ifges(3,2) + t666 + 0.4e1 * t903 - 0.2e1 * t1232) * t589 + t866 + t864) * t595 + ((t708 + 0.2e1 * t813) * t593 + t189 * t592 + t427 + t84 * t888 - t1116 * t1091 + t769) * t594 + (t1091 * t1173 + t943) * t592 + t1172 * t1091 + (t685 + t139 - 0.4e1 * t873 + t587 * t1168 + (((-0.4e1 * t550 + t956) * t592 - 0.4e1 * t1189 + t821) * t593 + (t299 * t587 + t187) * t592 + t214 * t587 - 0.4e1 * t1081 + t857) * t594 + ((-0.4e1 * t885 + t942) * t592 + 0.4e1 * t1201 + t824) * t593 + (4 * Ifges(3,4))) * t578 + ((t460 * t589 + t819) * t592 + t589 * t282 + t756) * t593) * qJD(1) + t906) * qJD(2) + ((0.2e1 * t765 + t649 + t682 + t14 + ((((t298 * t589 + t951) * t592 + t213 * t589 + t827) * t593 + (((t335 + t1095) * t589 + t372) * t587 + t129) * t592 + ((t359 + t318) * t589 + t396) * t587 + (-t192 - 0.2e1 * t1067) * t589 + t856) * t594 + ((t291 * t589 + t957) * t592 + t208 * t589 + t820) * t593 + (((t355 + t510) * t589 + t460) * t587 + t191) * t592 + ((t364 + t316) * t589 + t282) * t587 + (t666 + 0.2e1 * t903) * t589 + t479 * t871 + t864) * t595 + (((t131 - 0.2e1 * t1001) * t592 + 0.2e1 * t1141 * t431 + t947) * t593 + (t189 + 0.2e1 * t844) * t592 + 0.2e1 * t845 + t427 + t107 * t1089 + t1046) * t594 + (t685 - 0.2e1 * t873 + (((t537 + t956) * t592 + t319 + t821) * t593 + (t298 * t587 + t187) * t592 + t213 * t587 + t573 + t857) * t594 + ((t356 + t942) * t592 + t313 + t823) * t593 + (t291 * t587 + t158) * t592 + t208 * t587) * t578) * qJD(1) + t906) * qJD(3); 0.2e1 * qJD(3) * t795 + ((-0.2e1 * t876 - t1123 * t904 + (t534 + 0.2e1 * t1005) * t592 + 0.2e1 * t1028) * qJD(5) - 0.2e1 * t1174 * t569) * qJD(4) + (-t1174 + t1231) * t616 + (-t1171 + t689) * t615 - 0.2e1 * (-t1124 * t689 + t1171 * t569) * qJD(5) + (-t1240 + t1242) * qJDD(4) + (-0.2e1 * t1004 - 0.2e1 * (-t1003 - t1227) * pkin(2) - t683 + (t1222 + t749) * t593 + t1215) * qJDD(3) + (t675 + t497 + t1232 - t1128 * t1091 + t293 * t888 + Ifges(3,3) + 0.2e1 * t95 + (-t749 + t286) * t593 + t1158 - t1215) * qJDD(2) + (t5 - t6) * qJDD(1) + (t1209 * t589 - t1210 * t595) * g(3) + ((t840 + t998) * pkin(4) + t1138) * t614 * t451 + t441 * t910 * t1087 + (qJD(3) * t1191 - t617 - t618) * pkin(4) * t47 + (t338 + (pkin(4) * t998 + t1138) * t454) * qJDD(6) + (-t1004 - (-t846 + t1009) * pkin(2) + (-t284 * t587 + t1009) * t1074 + (-t288 * t586 - t266) * t593 - t1143 * t592 + t293 * t977) * qJDD(5) + (-qJD(5) * t689 - t795) * t1191 + (-t663 + t765 + t677 + t1239 * t543 + t1152 + (-t855 - t1067) * t588 + Ifges(3,4) + (t995 - t954 + (-t346 * t543 + t1212) * t588 + t822 + ((-t543 * t583 - t391) * t588 + t940 - t1186) * t592) * t593 + (t781 * t1092 + t674 - t691 - t873 + (t1201 - t826 + t756 + (-t885 - t941 + t819) * t592) * t593) * t578 + (t844 * t592 - t1081 + t845 + (-t1025 - t1081 + t769 - t858 - t1153) * t578 + (-t1189 + (t84 - t1197) * t592 + (t1214 - t1189 + t708 - t947 + (-t550 - t775 + t296) * t592) * t578) * t593) * t594 + (-t1217 * t578 + t1167 + 0.2e1 * t1207) * t577 + ((-t662 + t346 * t792 + (m(5) / 0.4e1 + t559) * t621 + (-t363 + t932) * t592 + t59 * t576 + t1023 + t457 - Ifges(4,1) / 0.4e1 + Ifges(4,2) / 0.4e1) * t1166 - t1075 * t1063 - t107 * t1088 + 0.4e1 * t862 + (t667 + t1237) * t589 + pkin(1) * mrSges(3,2) + ((t362 * t592 - t1203 / 0.2e1) * t1166 + t1225 * t588) * t593 + (t589 * t1067 + ((t744 - t745) * t592 + t1136) * t593) * t594) * t595 - t1160 + (-(((t464 + t391) * t592 + t281 - t1212) * t594 + (t371 - t882) * t592 + t1221 - t1202) * t578 + t1116 * t594 - t1173 * t592 + ((t1165 - t236 + t1152) * t593 * t1166 - t723 * t592 + t1225 * t594 + (t592 * t745 - t1136) * t588) * t595 - t1172) * t587) * t619 + t1121 * (t1209 * t595 + t1210 * t589) + (-t734 + t753 + 0.2e1 * (t800 - pkin(4) * t997 + (t1076 * t586 - t917) * t593 + (pkin(4) * t972 - pkin(5)) * t592) * t454) * t711; (t663 + (((-t410 * t917 / 0.4e1 + t1203 / 0.4e1 + (t782 + t899) * t587) * t593 - t481 / 0.2e1 + t662 + t799 / 0.4e1 + (t806 / 0.4e1 + t363 + t390 / 0.2e1) * t592 + t1056 + t1024) * t1166 + (-0.8e1 * t771 + (t542 * t592 + 0.8e1 * t761 + (pkin(4) * t368 + (0.16e2 * t849 + (t1090 * t229 + t509) * t592 + t1213) * t588) * t589 - t773) * t593 - 0.8e1 * t759 + (t815 + ((-0.2e1 * t809 - t550) * t587 + t171 * t1089) * t589) * t592 + (t431 * t739 - t311) * t587 + t588 * t900) * t594 + t58 * t893 + (-0.8e1 * t758 + (t796 + (t714 + 0.4e1 * t1027) * t589) * t592 + t587 * t900 + t1218 * t588 + t431 * t918) * t593 + t589 * t901 + ((t583 * t918 + t723) * t587 + t229 * t1088) * t592 + (t346 * t918 + t588 * t773) * t587 - 0.2e1 * t1033) * t595 + (-0.8e1 * t760 * t593 - t1115 * t446 + (t593 * t782 * t905 - 0.2e1 * t1021 + 0.4e1 * t1056 - 0.2e1 * t127 + t901) * t588) * t594 + (((t1154 - t541) * t592 + t825) * t593 + (-t810 + t943) * t592 + t1207 + t1170) * t577 + t1175 + (((-0.16e2 * t849 + t126 + t1220) * t593 + t1226) * t577 + (-0.8e1 * t58 * t984 + (t1189 + t550 * t592 + t304 + ((t1104 - 0.8e1 * t1027) * t592 - t236 * t1090 + t1222) * t588) * t593 + 0.8e1 * t219 * t986 + (t229 * t1089 + t301 * t587) * t592 + (t1213 * t588 - t1212) * t587 + t122) * t594 + ((t1154 + t289) * t592 - t1201 + t825) * t593 + (t587 * t714 + t943) * t592 - t209 * t587 + t1170) * t578 + (t1186 * t592 + t134 - 0.4e1 * t849 - t995) * t593) * t619 + t1240 * qJDD(3) + 0.2e1 * (qJD(3) * t1231 - t863) * qJD(2) + t617 * t1231 + (pkin(4) * t1036 + t1236) * t618 - mrSges(6,1) * t789 + mrSges(6,2) * t1065 + (-t481 + t724 + t878) * qJDD(5) - 0.2e1 * t1124 * t863 + t734 * t711 + (-t1072 * t451 - t453) * t614 - t615 * t1073 - t302 * qJDD(6) + (t695 + 0.2e1 * t878) * qJDD(4) - pkin(3) * (m(7) * t3 + t25 * mrSges(7,2) - t24 * mrSges(7,1) + t79 * (mrSges(7,1) * t89 - mrSges(7,3) * t79) - t78 * (-mrSges(7,2) * t89 + mrSges(7,3) * t78)) - Ifges(6,5) * t42 - Ifges(6,6) * t43 + t9 * qJDD(1) + t591 * (-mrSges(7,1) * t3 + mrSges(7,3) * t2 + Ifges(7,4) * t25 + Ifges(7,2) * t24 + Ifges(7,6) * t40 - t20 * t79 + t22 * t89) + ((-0.2e1 * t1179 * t594 - 0.2e1 * t588 * t67) * t595 + (-t1179 * t588 + t67 * t594) * t1088) * t910 + t585 * (mrSges(7,2) * t3 - mrSges(7,3) * t1 + Ifges(7,1) * t25 + Ifges(7,4) * t24 + Ifges(7,5) * t40 + t20 * t78 - t21 * t89) + t1242 * qJDD(2) - Ifges(6,3) * t480 + t90 * (Ifges(6,1) * t91 + Ifges(6,4) * t90 + Ifges(6,5) * t491) - t91 * (Ifges(6,4) * t91 + Ifges(6,2) * t90 + Ifges(6,6) * t491) + t1121 * (-t1243 * t589 + t56 * t595) - (-t1243 * t595 - t56 * t589) * g(3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + Ifges(7,5) * t25 + Ifges(7,6) * t24 + Ifges(7,3) * t40 + t21 * t79 - t22 * t78;];
tau = t8(:);
