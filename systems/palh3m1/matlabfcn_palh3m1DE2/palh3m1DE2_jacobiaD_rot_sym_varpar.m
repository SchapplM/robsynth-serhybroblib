% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh3m1DE2
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
%   Wie in palh3m1DE2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = palh3m1DE2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1DE2_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1DE2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_jacobiaD_rot_sym_varpar: pkin has to be [19x1] (double)');
JaD_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:42
	% EndTime: 2020-04-20 16:20:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:42
	% EndTime: 2020-04-20 16:20:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:42
	% EndTime: 2020-04-20 16:20:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:42
	% EndTime: 2020-04-20 16:20:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:23:42
	% EndTime: 2020-04-20 16:23:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:26:03
	% EndTime: 2020-04-20 16:40:58
	% DurationCPUTime: 733.44s
	% Computational Cost: add. (20094077->329), mult. (30878438->688), div. (1290831->26), fcn. (19372407->24), ass. (0->314)
	t1000 = cos(pkin(17));
	t1009 = pkin(4) ^ 2;
	t1002 = sin(qJ(3));
	t1012 = pkin(1) ^ 2;
	t1190 = sin(qJ(2));
	t1191 = sin(pkin(16));
	t1192 = cos(qJ(2));
	t1193 = cos(pkin(16));
	t989 = t1190 * t1191 - t1192 * t1193;
	t1187 = pkin(5) * t989;
	t1128 = -0.2e1 * pkin(1) * t1187 + t1012;
	t1205 = -pkin(6) - pkin(2);
	t978 = (pkin(5) - t1205) * (pkin(5) + t1205) + t1128;
	t1204 = -pkin(6) + pkin(2);
	t979 = (pkin(5) - t1204) * (pkin(5) + t1204) + t1128;
	t1153 = t979 * t978;
	t1014 = sqrt(-t1153);
	t990 = t1190 * t1193 + t1192 * t1191;
	t1186 = pkin(5) * t990;
	t1008 = pkin(5) ^ 2;
	t984 = t1008 + t1128;
	t980 = pkin(2) ^ 2 - pkin(6) ^ 2 + t984;
	t985 = pkin(1) - t1187;
	t972 = t1014 * t985 + t980 * t1186;
	t1133 = t972 * t1002;
	t1005 = cos(qJ(3));
	t1138 = t1014 * t990;
	t971 = -pkin(5) * t1138 + t985 * t980;
	t1134 = t971 * t1005;
	t1054 = -t1133 / 0.2e1 + t1134 / 0.2e1;
	t1011 = 0.1e1 / pkin(2);
	t981 = 0.1e1 / t984;
	t1139 = t1011 * t981;
	t962 = t1054 * t1139;
	t1132 = t972 * t1005;
	t1135 = t971 * t1002;
	t1053 = t1132 / 0.2e1 + t1135 / 0.2e1;
	t963 = t1053 * t1139;
	t996 = pkin(18) + pkin(19);
	t994 = sin(t996);
	t995 = cos(t996);
	t953 = -t962 * t995 + t963 * t994;
	t1188 = pkin(3) * t953;
	t1208 = -2 * pkin(4);
	t1129 = -t1188 * t1208 + t1009;
	t1203 = -pkin(8) - pkin(10);
	t943 = (pkin(3) - t1203) * (pkin(3) + t1203) + t1129;
	t1202 = -pkin(8) + pkin(10);
	t944 = (pkin(3) - t1202) * (pkin(3) + t1202) + t1129;
	t1154 = t944 * t943;
	t1013 = sqrt(-t1154);
	t1065 = -t962 * t994 - t963 * t995;
	t1189 = pkin(3) * t1065;
	t1010 = pkin(3) ^ 2;
	t949 = t1010 + t1129;
	t945 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t949;
	t950 = pkin(4) + t1188;
	t925 = -t1013 * t1189 + t945 * t950;
	t1137 = t925 * t1000;
	t1007 = 0.1e1 / pkin(10);
	t946 = 0.1e1 / t949;
	t1141 = t1007 * t946;
	t926 = t1013 * t950 + t945 * t1189;
	t999 = sin(pkin(17));
	t1155 = t926 * t999;
	t919 = (-t1137 / 0.2e1 + t1155 / 0.2e1) * t1141;
	t1136 = t926 * t1000;
	t1156 = t925 * t999;
	t920 = (t1136 / 0.2e1 + t1156 / 0.2e1) * t1141;
	t906 = qJ(2) + qJ(3) + atan2(t920, t919);
	t904 = sin(t906);
	t900 = t904 ^ 2;
	t905 = cos(t906);
	t902 = 0.1e1 / t905 ^ 2;
	t1160 = t900 * t902;
	t1063 = t1133 - t1134;
	t1206 = pkin(1) * pkin(5);
	t982 = 0.1e1 / t984 ^ 2;
	t1120 = t982 * t1206;
	t1180 = t1005 / 0.2e1;
	t1182 = -t1002 / 0.2e1;
	t1185 = t985 * pkin(1);
	t1089 = t980 + 0.2e1 * t1185;
	t988 = t989 * qJD(2);
	t1130 = t988 * t1014;
	t1061 = 0.2e1 * (t978 + t979) * t1206;
	t987 = t990 * qJD(2);
	t973 = t987 * t1061;
	t976 = 0.1e1 / t1014;
	t1214 = t973 * t976;
	t964 = -t1186 * t1214 / 0.2e1;
	t958 = t964 + (-t1089 * t987 + t1130) * pkin(5);
	t1213 = t1008 * t990;
	t1102 = t987 * t1213;
	t1080 = pkin(1) * t1102;
	t1196 = t976 / 0.2e1;
	t1097 = t985 * t1196;
	t1131 = t987 * t1014;
	t959 = t973 * t1097 - 0.2e1 * t1080 + (-t988 * t980 - t1131) * pkin(5);
	t934 = ((-t1053 * qJD(3) + t958 * t1180 + t959 * t1182) * t981 - t1063 * t987 * t1120) * t1011;
	t1062 = t1132 + t1135;
	t1050 = t987 * t1062;
	t1181 = t1002 / 0.2e1;
	t935 = (t1050 * t1120 + (t1054 * qJD(3) + t959 * t1180 + t958 * t1181) * t981) * t1011;
	t930 = t934 * t995 - t935 * t994;
	t1218 = t1013 * t930;
	t1001 = sin(qJ(4));
	t1004 = cos(qJ(4));
	t1006 = cos(qJ(1));
	t1121 = t1004 * t1006;
	t1003 = sin(qJ(1));
	t1124 = t1001 * t1003;
	t898 = -t905 * t1121 + t1124;
	t892 = 0.1e1 / t898 ^ 2;
	t1122 = t1003 * t1004;
	t1123 = t1001 * t1006;
	t896 = t905 * t1123 + t1122;
	t1162 = t892 * t896;
	t891 = 0.1e1 / t898;
	t1057 = t1001 * t891 + t1004 * t1162;
	t1040 = t1057 * t904;
	t890 = t896 ^ 2;
	t876 = t890 * t892 + 0.1e1;
	t874 = 0.1e1 / t876;
	t1036 = t874 * t1040;
	t1217 = -0.4e1 * t981 * t982;
	t947 = 0.1e1 / t949 ^ 2;
	t1119 = pkin(3) * pkin(4) * t947;
	t1037 = (-t1137 + t1155) * t1119;
	t1184 = -t1000 / 0.2e1;
	t1195 = t999 / 0.2e1;
	t1088 = t950 * t1208 - t945;
	t1197 = -t1065 / 0.2e1;
	t937 = 0.1e1 / t1013;
	t1099 = t937 * t1197;
	t1067 = -t934 * t994 - t935 * t995;
	t1212 = t1067 * t1013;
	t1127 = 2 * pkin(4);
	t1060 = pkin(3) * (t943 + t944) * t1127;
	t921 = t930 * t1060;
	t886 = (t1088 * t930 + t921 * t1099 - t1212) * pkin(3);
	t1082 = t1010 * t1065 * t1208;
	t1198 = t937 / 0.2e1;
	t1100 = t950 * t1198;
	t887 = t921 * t1100 + t930 * t1082 + (t1067 * t945 - t1218) * pkin(3);
	t869 = ((t886 * t1184 + t887 * t1195) * t946 + t930 * t1037) * t1007;
	t916 = 0.1e1 / t919 ^ 2;
	t1216 = t869 * t916;
	t1215 = t904 * t1160;
	t1078 = qJD(4) * t905 + qJD(1);
	t1143 = t1003 * t904;
	t881 = atan2(t1143, t905);
	t877 = sin(t881);
	t1146 = t1003 * t877;
	t1125 = qJD(1) * t1006;
	t1093 = t904 * t1125;
	t1144 = t1003 * t900;
	t1158 = t916 * t920;
	t1038 = (t1136 + t1156) * t1119;
	t1183 = t1000 / 0.2e1;
	t870 = ((t887 * t1183 + t886 * t1195) * t946 + t930 * t1038) * t1007;
	t918 = t920 ^ 2;
	t909 = t916 * t918 + 0.1e1;
	t907 = 0.1e1 / t909;
	t915 = 0.1e1 / t919;
	t855 = qJD(2) + qJD(3) + (-t869 * t1158 + t870 * t915) * t907;
	t1147 = t1003 * t855;
	t997 = t1003 ^ 2;
	t882 = t997 * t1160 + 0.1e1;
	t879 = 0.1e1 / t882;
	t901 = 0.1e1 / t905;
	t848 = ((t905 * t1147 + t1093) * t901 + t855 * t902 * t1144) * t879;
	t1211 = -t855 * t1146 + t848 * t877;
	t878 = cos(t881);
	t1210 = t878 * t1125 - t848 * t1146 + t855 * t877;
	t974 = t990 * t1061;
	t1074 = t974 * t1196 + t980;
	t960 = (t989 * t1014 + (-t1074 - 0.2e1 * t1185) * t990) * pkin(5);
	t1209 = (t972 * qJD(3) - t958) * t990 - t960 * t987 + t971 * t988;
	t1109 = t877 * t1143;
	t867 = t878 * t905 + t1109;
	t864 = 0.1e1 / t867;
	t1207 = -0.2e1 * t874;
	t865 = 0.1e1 / t867 ^ 2;
	t1201 = -t921 / 0.2e1;
	t1179 = pkin(1) * t1008;
	t961 = t974 * t1097 - 0.2e1 * t990 ^ 2 * t1179 + (-t980 * t989 - t1138) * pkin(5);
	t1055 = t961 * t1180 + t960 * t1181;
	t1085 = t990 * t1120;
	t941 = (t1055 * t981 + t1062 * t1085) * t1011;
	t1056 = t960 * t1180 + t961 * t1182;
	t942 = (-t1056 * t981 + t1063 * t1085) * t1011;
	t933 = -t941 * t994 - t942 * t995;
	t922 = t933 * t1060;
	t1200 = -t922 / 0.2e1;
	t931 = t1065 * t1060;
	t1199 = -t931 / 0.2e1;
	t1194 = t879 - 0.1e1;
	t1069 = t1125 * t1144;
	t1163 = t878 * t904;
	t843 = (t1003 * t848 - t855) * t1163 + (t1093 + (-t848 + t1147) * t905) * t877;
	t1118 = -0.2e1 * t843 * t864 * t865;
	t998 = t1006 ^ 2;
	t1159 = t900 * t998;
	t1170 = t855 * t905;
	t862 = t865 * t1159 + 0.1e1;
	t1178 = (t1118 * t1159 + 0.2e1 * (t1170 * t904 * t998 - t1069) * t865) / t862 ^ 2;
	t1077 = qJD(1) * t905 + qJD(4);
	t1142 = t1006 * t904;
	t1035 = t1077 * t1003 + t855 * t1142;
	t854 = t1035 * t1004 + t1078 * t1123;
	t893 = t891 * t892;
	t1173 = t854 * t893;
	t853 = t1035 * t1001 - t1078 * t1121;
	t1177 = 0.2e1 * (-t853 * t1162 - t890 * t1173) / t876 ^ 2;
	t1039 = 0.2e1 * t855 * (t1215 + t904) * t901;
	t1175 = (t1039 * t997 + 0.2e1 * t902 * t1069) / t882 ^ 2;
	t1174 = t853 * t892;
	t1171 = t855 * t904;
	t1165 = t915 * t1216;
	t1169 = 0.2e1 * (t870 * t1158 - t918 * t1165) / t909 ^ 2;
	t1168 = t864 * t905;
	t1167 = t865 * t904;
	t1166 = t865 * t905;
	t1164 = t877 * t905;
	t1161 = t900 * t901;
	t1157 = t921 * t937 / t1154;
	t1151 = qJD(1) * t904;
	t1149 = t1000 * t946;
	t1145 = t1003 * t878;
	t1140 = t1010 * t930;
	t1126 = qJD(1) * t1003;
	t1117 = -0.2e1 * t893 * t896;
	t1115 = t904 * t1178;
	t1114 = t904 * t1175;
	t1113 = 0.1e1 / t1153 * t974 * t1214;
	t1112 = t1194 * t904;
	t1108 = t1003 * t879 * t901;
	t1106 = t1006 * t1178;
	t1105 = t1006 * t855 * t864;
	t1104 = t865 * t1142;
	t1103 = t1006 * t1166;
	t1101 = t1157 / 0.4e1;
	t1098 = t946 * t1195;
	t1096 = 0.1e1 + t1160;
	t1095 = t864 * t1126;
	t1092 = t1009 * t1140;
	t1091 = -t1149 / 0.2e1;
	t1090 = t1149 / 0.2e1;
	t1087 = t1010 * t1127;
	t1084 = t907 * t1119;
	t1083 = 0.4e1 * pkin(4) * t1140;
	t1079 = -0.8e1 * t1092;
	t1076 = t900 * t1108;
	t1075 = -t1065 * t1157 / 0.4e1;
	t1071 = t1012 * t1102;
	t1070 = t946 * t947 * t1092;
	t1068 = t1096 * t879;
	t1064 = t1000 * t1070;
	t1059 = t1003 * t1068;
	t1058 = 0.4e1 * t999 * t1070;
	t1052 = t933 * t1058;
	t1051 = t1065 * t1058;
	t1049 = t1068 * t1125;
	t1048 = -0.4e1 * t925 * t1064;
	t1047 = 0.4e1 * t926 * t1064;
	t1046 = -t915 * t1169 - t907 * t1216;
	t1045 = t1006 * t1118 - t865 * t1126;
	t932 = -t941 * t995 + t942 * t994;
	t888 = (-t932 * t1013 + t1088 * t933 + t922 * t1099) * pkin(3);
	t1034 = t961 * t987 - t972 * t988 + (qJD(3) * t971 + t959) * t990;
	t970 = -t1061 * t988 - 0.8e1 * t1071;
	t939 = t964 + (t1113 / 0.4e1 + t970 * t1196) * t985 + (0.2e1 * t987 * t989 + 0.4e1 * t988 * t990) * t1179 + (-t1074 * t987 + t1130) * pkin(5);
	t940 = 0.4e1 * t1080 + (t1131 - t990 * t1113 / 0.4e1 + t1089 * t988 + (t989 * t973 / 0.2e1 + t988 * t974 / 0.2e1 - t990 * t970 / 0.2e1) * t976) * pkin(5);
	t923 = (-t1063 * t1071 * t1217 + (-t940 * t1005 / 0.2e1 + t939 * t1181 + t1055 * qJD(3)) * t981 + (t1034 * t1002 + t1209 * t1005) * t1120) * t1011;
	t924 = (-t1012 * t1050 * t1213 * t1217 + (t1056 * qJD(3) + t939 * t1180 + t940 * t1181) * t981 + (-t1209 * t1002 + t1034 * t1005) * t1120) * t1011;
	t911 = -t923 * t995 - t924 * t994;
	t1044 = t886 * t933 + t888 * t930 + t911 * t925;
	t913 = (-t953 * t1013 + t1065 * t1088 + t931 * t1099) * pkin(3);
	t1043 = t1065 * t886 + t1067 * t925 + t913 * t930;
	t889 = t922 * t1100 + t933 * t1082 + (-t1013 * t933 + t932 * t945) * pkin(3);
	t1042 = t887 * t933 + t889 * t930 + t911 * t926;
	t914 = t931 * t1100 + t1065 * t1082 + (-t1013 * t1065 + t945 * t953) * pkin(3);
	t1041 = t1065 * t887 + t1067 * t926 + t914 * t930;
	t1033 = t1126 * t1036;
	t1032 = t1158 * t1169 + (0.2e1 * t920 * t1165 - t870 * t916) * t907;
	t1031 = t879 * t1039 - t1096 * t1175;
	t1030 = -t1040 * t1177 + (t1057 * t1170 + ((-qJD(4) * t896 - t854) * t892 * t1001 + (qJD(4) * t891 + t854 * t1117 - t1174) * t1004) * t904) * t874;
	t912 = t1060 * t1067 + t1065 * t1079;
	t910 = t923 * t994 - t924 * t995;
	t895 = t905 * t1122 + t1123;
	t894 = t905 * t1124 - t1121;
	t885 = ((t914 * t1183 + t913 * t1195) * t946 + t1065 * t1038) * t1007;
	t884 = ((t913 * t1184 + t914 * t1195) * t946 + t1065 * t1037) * t1007;
	t883 = t1060 * t911 + t933 * t1079;
	t873 = ((t889 * t1183 + t888 * t1195) * t946 + t933 * t1038) * t1007;
	t872 = ((t888 * t1184 + t889 * t1195) * t946 + t933 * t1037) * t1007;
	t871 = t1065 * t1083 + (t1218 + t931 * t1075 + t1088 * t1067 + (t1067 * t1199 + t912 * t1197 + t953 * t1201) * t937) * pkin(3);
	t868 = (t931 * t1101 + t912 * t1198) * t950 + (-0.2e1 * t1065 * t1067 - t930 * t953) * t1087 + (-t930 * t945 - t1212 + (t1065 * t1201 + t930 * t1199) * t937) * pkin(3);
	t863 = t933 * t1083 + (-t910 * t1013 + t922 * t1075 + t1088 * t911 + (t1067 * t1200 + t883 * t1197 + t932 * t1201) * t937) * pkin(3);
	t860 = 0.1e1 / t862;
	t859 = 0.1e1 + (-t884 * t1158 + t885 * t915) * t907;
	t858 = (t922 * t1101 + t883 * t1198) * t950 + (-t1065 * t911 - t1067 * t933 - t930 * t932) * t1087 + (t910 * t945 - t911 * t1013 + (t930 * t1200 + t933 * t1201) * t937) * pkin(3);
	t856 = 0.1e1 + (-t872 * t1158 + t873 * t915) * t907;
	t852 = (t878 * t1076 - t877 * t1112) * t1006;
	t851 = t859 * t1059;
	t850 = t856 * t1059;
	t846 = -t851 * t1164 - t859 * t1163 + (t851 * t1163 + t859 * t1164) * t1003;
	t845 = t1046 * t885 + t1032 * t884 + (((t1047 * t1065 + t925 * t1051 + t868 * t1090 + t871 * t1098) * t915 - (t1048 * t1065 + t926 * t1051 + t871 * t1091 + t868 * t1098) * t1158) * t907 + ((-t1041 * t1158 + t1043 * t915) * t999 + (t1041 * t915 + t1043 * t1158) * t1000) * t1084) * t1007;
	t844 = -t850 * t1164 - t856 * t1163 + (t850 * t1163 + t856 * t1164) * t1003;
	t842 = t1046 * t873 + t1032 * t872 + (((t1047 * t933 + t1052 * t925 + t1090 * t858 + t1098 * t863) * t915 - (t1048 * t933 + t1052 * t926 + t1091 * t863 + t1098 * t858) * t1158) * t907 + ((-t1042 * t1158 + t1044 * t915) * t999 + (t1042 * t915 + t1044 * t1158) * t1000) * t1084) * t1007;
	t840 = t859 * t1049 + (t1031 * t859 + t1068 * t845) * t1003;
	t839 = t856 * t1049 + (t1031 * t856 + t1068 * t842) * t1003;
	t1 = [-t1108 * t1151 + (t1068 * t855 - t1114 * t901) * t1006, t839, t840, 0; (-t864 * t1115 + (t855 * t1168 + (-qJD(1) * t852 - t843) * t1167) * t860) * t1003 + (-t865 * t1115 * t852 + (((-t848 * t1076 - t1194 * t1170 + t1114) * t877 + (-t848 * t1112 + (-t1161 * t1175 + (0.2e1 * t904 + t1215) * t879 * t855) * t1003) * t878) * t1104 + (t1118 * t904 + t1166 * t855) * t852 + (t864 + ((-t997 + t998) * t879 * t878 * t1161 + t1194 * t1109) * t865) * t1151) * t860) * t1006, (-t1167 * t844 + t1168 * t856) * t1106 + ((t856 * t1095 + (-t842 * t864 + (t843 * t856 + t844 * t855) * t865) * t1006) * t905 + (t856 * t1105 + (t1145 * t839 + t1210 * t850 + t1211 * t856 - t842 * t878) * t1104 + t1045 * t844 + ((t1003 * t842 + t1125 * t856 - t839) * t877 + (-t848 * t850 - t855 * t856 + (t848 * t856 + t850 * t855) * t1003) * t878) * t1103) * t904) * t860, (-t1167 * t846 + t1168 * t859) * t1106 + ((t859 * t1095 + (-t845 * t864 + (t843 * t859 + t846 * t855) * t865) * t1006) * t905 + (t859 * t1105 + (t1145 * t840 + t1210 * t851 + t1211 * t859 - t845 * t878) * t1104 + t1045 * t846 + ((t1003 * t845 + t1125 * t859 - t840) * t877 + (-t848 * t851 - t855 * t859 + (t848 * t859 + t851 * t855) * t1003) * t878) * t1103) * t904) * t860, 0; (-t1162 * t895 - t891 * t894) * t1177 + (-t895 * t1174 + (t1117 * t895 - t892 * t894) * t854 + ((-t1001 * t1171 + t1078 * t1004) * t891 + (-t1078 * t1001 - t1004 * t1171) * t1162) * t1003 + t1057 * t1006 * t1077) * t874, -t856 * t1033 + (t1030 * t856 + t1036 * t842) * t1006, -t859 * t1033 + (t1030 * t859 + t1036 * t845) * t1006, -t1177 + (t1174 * t1207 + (t1173 * t1207 - t1177 * t892) * t896) * t896;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:50
	% EndTime: 2020-04-20 16:20:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:50
	% EndTime: 2020-04-20 16:20:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:22:25
	% EndTime: 2020-04-20 16:22:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
end