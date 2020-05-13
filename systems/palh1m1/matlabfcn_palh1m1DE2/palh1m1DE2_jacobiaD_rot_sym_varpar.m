% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh1m1DE2
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
%   Wie in palh1m1DE2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = palh1m1DE2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m1DE2_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1DE2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_jacobiaD_rot_sym_varpar: pkin has to be [23x1] (double)');
JaD_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:35
	% EndTime: 2020-04-15 18:49:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:35
	% EndTime: 2020-04-15 18:49:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:35
	% EndTime: 2020-04-15 18:49:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:35
	% EndTime: 2020-04-15 18:49:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:51:50
	% EndTime: 2020-04-15 18:51:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:53:45
	% EndTime: 2020-04-15 19:06:37
	% DurationCPUTime: 692.72s
	% Computational Cost: add. (20094671->332), mult. (30878438->695), div. (1290831->26), fcn. (19372407->24), ass. (0->316)
	t1008 = 0.1e1 / pkin(11);
	t1011 = pkin(4) ^ 2;
	t1010 = pkin(5) ^ 2;
	t1005 = cos(qJ(3));
	t1013 = pkin(1) ^ 2;
	t1002 = sin(qJ(2));
	t1007 = cos(pkin(19));
	t1188 = sin(pkin(19));
	t1189 = cos(qJ(2));
	t1046 = -t1002 * t1007 + t1189 * t1188;
	t1185 = pkin(7) * t1046;
	t1207 = -2 * pkin(1);
	t1128 = -t1185 * t1207 + t1013;
	t1204 = -pkin(8) - pkin(3);
	t977 = (pkin(7) - t1204) * (pkin(7) + t1204) + t1128;
	t1203 = -pkin(8) + pkin(3);
	t978 = (pkin(7) - t1203) * (pkin(7) + t1203) + t1128;
	t1149 = t978 * t977;
	t1015 = sqrt(-t1149);
	t989 = t1002 * t1188 + t1189 * t1007;
	t1184 = pkin(7) * t989;
	t1009 = pkin(7) ^ 2;
	t983 = t1009 + t1128;
	t979 = pkin(3) ^ 2 - pkin(8) ^ 2 + t983;
	t984 = pkin(1) + t1185;
	t971 = t1015 * t984 + t979 * t1184;
	t1137 = t1005 * t971;
	t1001 = sin(qJ(3));
	t1132 = t1015 * t989;
	t970 = -pkin(7) * t1132 + t979 * t984;
	t1145 = t1001 * t970;
	t1057 = t1145 / 0.2e1 + t1137 / 0.2e1;
	t1012 = 0.1e1 / pkin(3);
	t980 = 0.1e1 / t983;
	t1133 = t1012 * t980;
	t962 = t1057 * t1133;
	t1138 = t1005 * t970;
	t1144 = t1001 * t971;
	t1055 = t1138 / 0.2e1 - t1144 / 0.2e1;
	t963 = t1055 * t1133;
	t995 = pkin(23) + pkin(22);
	t993 = sin(t995);
	t994 = cos(t995);
	t951 = t962 * t994 + t963 * t993;
	t1187 = pkin(4) * t951;
	t1206 = -2 * pkin(5);
	t1129 = -t1187 * t1206 + t1010;
	t948 = t1011 + t1129;
	t945 = 0.1e1 / t948;
	t1135 = t1008 * t945;
	t1202 = -pkin(9) - pkin(11);
	t942 = (pkin(4) - t1202) * (pkin(4) + t1202) + t1129;
	t1201 = -pkin(9) + pkin(11);
	t943 = (pkin(4) - t1201) * (pkin(4) + t1201) + t1129;
	t1151 = t943 * t942;
	t1014 = sqrt(-t1151);
	t1080 = t962 * t993 - t994 * t963;
	t1186 = pkin(4) * t1080;
	t944 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t948;
	t949 = pkin(5) + t1187;
	t924 = t1014 * t949 + t944 * t1186;
	t999 = cos(pkin(21));
	t1152 = t924 * t999;
	t923 = -t1014 * t1186 + t944 * t949;
	t998 = sin(pkin(21));
	t1155 = t923 * t998;
	t917 = (t1155 / 0.2e1 + t1152 / 0.2e1) * t1135;
	t1153 = t924 * t998;
	t1154 = t923 * t999;
	t918 = (-t1154 / 0.2e1 + t1153 / 0.2e1) * t1135;
	t904 = qJ(2) + qJ(3) + atan2(t917, t918);
	t902 = sin(t904);
	t898 = t902 ^ 2;
	t903 = cos(t904);
	t900 = 0.1e1 / t903 ^ 2;
	t1159 = t898 * t900;
	t1066 = -t1138 + t1144;
	t1205 = pkin(1) * pkin(7);
	t981 = 0.1e1 / t983 ^ 2;
	t1119 = t981 * t1205;
	t1037 = t1066 * t1119;
	t1182 = -t1005 / 0.2e1;
	t1183 = t1001 / 0.2e1;
	t1120 = t984 * t1207;
	t1087 = -t979 + t1120;
	t986 = t1046 * qJD(2);
	t1131 = t986 * t1015;
	t975 = 0.1e1 / t1015;
	t1195 = -t975 / 0.2e1;
	t1064 = 0.2e1 * (t977 + t978) * t1205;
	t987 = t989 * qJD(2);
	t972 = t987 * t1064;
	t964 = t972 * t1184 * t1195;
	t957 = t964 + (t1087 * t987 - t1131) * pkin(7);
	t1213 = t1009 * t987;
	t1101 = t989 * t1213;
	t1078 = pkin(1) * t1101;
	t1194 = t975 / 0.2e1;
	t1094 = t984 * t1194;
	t1130 = t987 * t1015;
	t958 = t972 * t1094 - 0.2e1 * t1078 + (t986 * t979 - t1130) * pkin(7);
	t933 = (t987 * t1037 + (t1057 * qJD(3) + t957 * t1182 + t958 * t1183) * t980) * t1012;
	t1065 = t1137 + t1145;
	t1181 = t1005 / 0.2e1;
	t934 = (t1065 * t987 * t1119 + (t1055 * qJD(3) + t958 * t1181 + t957 * t1183) * t980) * t1012;
	t928 = t933 * t993 - t934 * t994;
	t1216 = t1014 * t928;
	t1000 = sin(qJ(4));
	t1004 = cos(qJ(4));
	t1006 = cos(qJ(1));
	t1121 = t1004 * t1006;
	t1003 = sin(qJ(1));
	t1124 = t1000 * t1003;
	t895 = t903 * t1121 + t1124;
	t890 = 0.1e1 / t895 ^ 2;
	t1122 = t1003 * t1004;
	t1123 = t1000 * t1006;
	t894 = t1123 * t903 - t1122;
	t1162 = t890 * t894;
	t889 = 0.1e1 / t895;
	t1059 = -t1000 * t889 + t1004 * t1162;
	t1041 = t1059 * t902;
	t888 = t894 ^ 2;
	t874 = t888 * t890 + 0.1e1;
	t872 = 0.1e1 / t874;
	t1036 = t872 * t1041;
	t946 = 0.1e1 / t948 ^ 2;
	t1118 = pkin(4) * pkin(5) * t946;
	t1040 = (t1153 - t1154) * t1118;
	t1192 = -t999 / 0.2e1;
	t1193 = t998 / 0.2e1;
	t1086 = t949 * t1206 - t944;
	t1196 = -t1080 / 0.2e1;
	t936 = 0.1e1 / t1014;
	t1098 = t936 * t1196;
	t1081 = t994 * t933 + t934 * t993;
	t1212 = t1081 * t1014;
	t1127 = 2 * pkin(5);
	t1063 = pkin(4) * (t942 + t943) * t1127;
	t919 = t928 * t1063;
	t884 = (t1086 * t928 + t1098 * t919 - t1212) * pkin(4);
	t1082 = t1011 * t1080 * t1206;
	t1197 = t936 / 0.2e1;
	t1099 = t949 * t1197;
	t885 = t919 * t1099 + t928 * t1082 + (t1081 * t944 - t1216) * pkin(4);
	t868 = ((t884 * t1192 + t885 * t1193) * t945 + t928 * t1040) * t1008;
	t915 = 0.1e1 / t918 ^ 2;
	t1215 = t868 * t915;
	t1214 = t902 * t1159;
	t1211 = -qJD(4) * t903 + qJD(1);
	t1139 = t1003 * t902;
	t877 = atan2(-t1139, -t903);
	t875 = sin(t877);
	t1142 = t1003 * t875;
	t1125 = qJD(1) * t1006;
	t1090 = t902 * t1125;
	t1140 = t1003 * t898;
	t1157 = t915 * t917;
	t1039 = (t1152 + t1155) * t1118;
	t1191 = t999 / 0.2e1;
	t867 = ((t885 * t1191 + t884 * t1193) * t945 + t928 * t1039) * t1008;
	t913 = t917 ^ 2;
	t907 = t913 * t915 + 0.1e1;
	t905 = 0.1e1 / t907;
	t914 = 0.1e1 / t918;
	t853 = qJD(2) + qJD(3) + (-t1157 * t868 + t867 * t914) * t905;
	t1143 = t1003 * t853;
	t996 = t1003 ^ 2;
	t880 = t1159 * t996 + 0.1e1;
	t878 = 0.1e1 / t880;
	t899 = 0.1e1 / t903;
	t846 = (-(-t1143 * t903 - t1090) * t899 + t853 * t900 * t1140) * t878;
	t1210 = t853 * t1142 - t846 * t875;
	t876 = cos(t877);
	t1209 = -t876 * t1125 + t846 * t1142 - t853 * t875;
	t1180 = pkin(1) * t1009;
	t973 = t989 * t1064;
	t960 = t973 * t1094 - 0.2e1 * t989 ^ 2 * t1180 + (t1046 * t979 - t1132) * pkin(7);
	t1208 = t960 * t987 + t971 * t986 + (qJD(3) * t970 + t958) * t989;
	t1109 = t875 * t1139;
	t865 = -t876 * t903 - t1109;
	t862 = 0.1e1 / t865;
	t863 = 0.1e1 / t865 ^ 2;
	t1200 = -t919 / 0.2e1;
	t1050 = t989 * t1065;
	t1073 = t973 * t1195 - t979;
	t959 = (-t1046 * t1015 + (t1073 + t1120) * t989) * pkin(7);
	t1058 = t960 * t1181 + t959 * t1183;
	t940 = (-t1050 * t1119 - t1058 * t980) * t1012;
	t1056 = t959 * t1182 + t960 * t1183;
	t941 = (t989 * t1037 + t1056 * t980) * t1012;
	t930 = t940 * t994 + t941 * t993;
	t920 = t930 * t1063;
	t1199 = -t920 / 0.2e1;
	t929 = t1080 * t1063;
	t1198 = -t929 / 0.2e1;
	t1190 = t878 - 0.1e1;
	t1070 = t1125 * t1140;
	t997 = t1006 ^ 2;
	t1158 = t898 * t997;
	t1170 = t853 * t903;
	t1163 = t876 * t902;
	t841 = (-t1003 * t846 + t853) * t1163 + (-t1090 + (t846 - t1143) * t903) * t875;
	t1178 = t841 * t862 * t863;
	t860 = t1158 * t863 + 0.1e1;
	t1179 = 0.2e1 * (-t1158 * t1178 + (t902 * t997 * t1170 - t1070) * t863) / t860 ^ 2;
	t1076 = -qJD(1) * t903 + qJD(4);
	t1136 = t1006 * t902;
	t1104 = t853 * t1136;
	t852 = t1211 * t1123 + (t1003 * t1076 - t1104) * t1004;
	t891 = t889 * t890;
	t1173 = t852 * t891;
	t1054 = t903 * t1124 + t1121;
	t851 = qJD(1) * t1054 - t895 * qJD(4) + t1000 * t1104;
	t1177 = 0.2e1 * (-t1162 * t851 - t888 * t1173) / t874 ^ 2;
	t1038 = 0.2e1 * t853 * (t1214 + t902) * t899;
	t1175 = (t1038 * t996 + 0.2e1 * t1070 * t900) / t880 ^ 2;
	t1174 = t851 * t890;
	t1171 = t853 * t902;
	t1165 = t914 * t1215;
	t1169 = 0.2e1 * (t1157 * t867 - t1165 * t913) / t907 ^ 2;
	t1168 = t862 * t903;
	t1167 = t863 * t902;
	t1166 = t863 * t903;
	t1164 = t875 * t903;
	t893 = -t1122 * t903 + t1123;
	t1161 = t893 * t894;
	t1160 = t898 * t899;
	t1156 = t919 * t936 / t1151;
	t1150 = t945 * t999;
	t1148 = qJD(1) * t902;
	t1141 = t1003 * t876;
	t1134 = t1011 * t928;
	t1126 = qJD(1) * t1003;
	t1117 = 0.2e1 * t1178;
	t1116 = -0.2e1 * t1173;
	t1115 = t902 * t1179;
	t1114 = t902 * t1175;
	t1113 = t975 / t1149 * t973 * t972;
	t1112 = t1190 * t902;
	t1108 = t1003 * t878 * t899;
	t1106 = t1006 * t1179;
	t1105 = t1006 * t853 * t862;
	t1103 = t863 * t1136;
	t1102 = t1006 * t1166;
	t1100 = t1156 / 0.4e1;
	t1097 = t945 * t1193;
	t1096 = -t1150 / 0.2e1;
	t1095 = t1150 / 0.2e1;
	t1093 = 0.1e1 + t1159;
	t1092 = t862 * t1126;
	t1088 = t1010 * t1134;
	t1085 = t1011 * t1127;
	t1084 = t905 * t1118;
	t1083 = 0.4e1 * pkin(5) * t1134;
	t1077 = -0.8e1 * t1088;
	t1075 = t898 * t1108;
	t1074 = -t1080 * t1156 / 0.4e1;
	t1072 = t1013 * t1101;
	t1071 = t945 * t946 * t1088;
	t1069 = t1093 * t878;
	t1068 = 0.4e1 * t1071;
	t1062 = t1003 * t1069;
	t1061 = t924 * t1068;
	t1060 = t998 * t1068;
	t1053 = -0.4e1 * t1071 * t1154;
	t1052 = t1080 * t1060;
	t1051 = t930 * t1061;
	t1049 = t1069 * t1125;
	t1048 = -t914 * t1169 - t905 * t1215;
	t1047 = t1006 * t1117 + t863 * t1126;
	t931 = -t940 * t993 + t941 * t994;
	t886 = (-t931 * t1014 + t1086 * t930 + t1098 * t920) * pkin(4);
	t1035 = -t959 * t987 - t970 * t986 + (qJD(3) * t971 - t957) * t989;
	t969 = t986 * t1064 - 0.8e1 * t1072;
	t938 = t964 + (t1113 / 0.4e1 + t969 * t1194) * t984 + (-0.2e1 * t1046 * t987 - 0.4e1 * t986 * t989) * t1180 + (t1073 * t987 - t1131) * pkin(7);
	t939 = 0.4e1 * t1078 + (t1130 - t989 * t1113 / 0.4e1 + t1087 * t986 + (-t1046 * t972 / 0.2e1 - t986 * t973 / 0.2e1 - t989 * t969 / 0.2e1) * t975) * pkin(7);
	t982 = t980 * t981;
	t921 = (0.4e1 * t1066 * t982 * t1072 + (t1058 * qJD(3) + t939 * t1182 + t938 * t1183) * t980 + (t1208 * t1001 + t1035 * t1005) * t1119) * t1012;
	t922 = (-0.4e1 * t982 * t1013 * t1050 * t1213 + (-t1001 * t939 / 0.2e1 + t938 * t1182 + t1056 * qJD(3)) * t980 + (t1035 * t1001 - t1208 * t1005) * t1119) * t1012;
	t909 = t921 * t993 + t922 * t994;
	t1045 = t884 * t930 + t886 * t928 + t909 * t923;
	t911 = (-t951 * t1014 + t1080 * t1086 + t1098 * t929) * pkin(4);
	t1044 = t1080 * t884 + t1081 * t923 + t911 * t928;
	t887 = t920 * t1099 + t930 * t1082 + (-t1014 * t930 + t931 * t944) * pkin(4);
	t1043 = t885 * t930 + t887 * t928 + t909 * t924;
	t912 = t929 * t1099 + t1080 * t1082 + (-t1014 * t1080 + t944 * t951) * pkin(4);
	t1042 = t1080 * t885 + t1081 * t924 + t912 * t928;
	t1034 = t1126 * t1036;
	t1033 = t1157 * t1169 + (0.2e1 * t917 * t1165 - t867 * t915) * t905;
	t1032 = t878 * t1038 - t1093 * t1175;
	t1031 = -t1041 * t1177 + (t1059 * t1170 + ((-qJD(4) * t894 + t852) * t890 * t1000 + (-qJD(4) * t889 + t894 * t1116 - t1174) * t1004) * t902) * t872;
	t910 = t1063 * t1081 + t1077 * t1080;
	t908 = t921 * t994 - t922 * t993;
	t883 = ((t911 * t1192 + t912 * t1193) * t945 + t1080 * t1040) * t1008;
	t882 = ((t912 * t1191 + t911 * t1193) * t945 + t1080 * t1039) * t1008;
	t881 = t1063 * t909 + t1077 * t930;
	t871 = ((t886 * t1192 + t887 * t1193) * t945 + t930 * t1040) * t1008;
	t870 = ((t887 * t1191 + t886 * t1193) * t945 + t930 * t1039) * t1008;
	t869 = t1080 * t1083 + (t1216 + t929 * t1074 + t1086 * t1081 + (t1081 * t1198 + t910 * t1196 + t951 * t1200) * t936) * pkin(4);
	t866 = (t929 * t1100 + t910 * t1197) * t949 + (-0.2e1 * t1080 * t1081 - t928 * t951) * t1085 + (-t1212 - t928 * t944 + (t1080 * t1200 + t928 * t1198) * t936) * pkin(4);
	t861 = t930 * t1083 + (-t908 * t1014 + t920 * t1074 + t1086 * t909 + (t1081 * t1199 + t881 * t1196 + t931 * t1200) * t936) * pkin(4);
	t858 = 0.1e1 / t860;
	t857 = (t920 * t1100 + t881 * t1197) * t949 + (-t1080 * t909 - t1081 * t930 - t928 * t931) * t1085 + (t908 * t944 - t909 * t1014 + (t928 * t1199 + t930 * t1200) * t936) * pkin(4);
	t856 = 0.1e1 + (-t1157 * t883 + t882 * t914) * t905;
	t854 = 0.1e1 + (-t1157 * t871 + t870 * t914) * t905;
	t850 = (-t1075 * t876 + t1112 * t875) * t1006;
	t849 = t856 * t1062;
	t848 = t854 * t1062;
	t844 = t849 * t1164 + t856 * t1163 + (-t1163 * t849 - t1164 * t856) * t1003;
	t843 = t1048 * t882 + t1033 * t883 + (((t1061 * t1080 * t999 + t1052 * t923 + t1095 * t866 + t1097 * t869) * t914 - (t1052 * t924 + t1053 * t1080 + t1096 * t869 + t1097 * t866) * t1157) * t905 + ((t1042 * t914 + t1044 * t1157) * t999 + (-t1042 * t1157 + t1044 * t914) * t998) * t1084) * t1008;
	t842 = t848 * t1164 + t854 * t1163 + (-t1163 * t848 - t1164 * t854) * t1003;
	t840 = t1048 * t870 + t1033 * t871 + (((t1060 * t923 * t930 + t1051 * t999 + t1095 * t857 + t1097 * t861) * t914 - (t1051 * t998 + t1053 * t930 + t1096 * t861 + t1097 * t857) * t1157) * t905 + ((t1043 * t914 + t1045 * t1157) * t999 + (-t1043 * t1157 + t1045 * t914) * t998) * t1084) * t1008;
	t838 = t856 * t1049 + (t1032 * t856 + t1069 * t843) * t1003;
	t837 = t854 * t1049 + (t1032 * t854 + t1069 * t840) * t1003;
	t1 = [-t1108 * t1148 + (t1069 * t853 - t1114 * t899) * t1006, t837, t838, 0; (t862 * t1115 + (-t853 * t1168 + (qJD(1) * t850 + t841) * t1167) * t858) * t1003 + (t863 * t1115 * t850 + (-((t846 * t1075 + t1190 * t1170 - t1114) * t875 + (t846 * t1112 + (t1160 * t1175 + (-0.2e1 * t902 - t1214) * t878 * t853) * t1003) * t876) * t1103 + (t1117 * t902 - t1166 * t853) * t850 + (-t862 + (-(t996 - t997) * t878 * t876 * t1160 + t1190 * t1109) * t863) * t1148) * t858) * t1006, (t1167 * t842 - t1168 * t854) * t1106 + ((-t854 * t1092 + (t840 * t862 + (-t841 * t854 - t842 * t853) * t863) * t1006) * t903 + (-t854 * t1105 - (-t1141 * t837 + t1209 * t848 + t1210 * t854 + t840 * t876) * t1103 + t1047 * t842 - ((-t1003 * t840 - t1125 * t854 + t837) * t875 + (t846 * t848 + t853 * t854 + (-t846 * t854 - t848 * t853) * t1003) * t876) * t1102) * t902) * t858, (t1167 * t844 - t1168 * t856) * t1106 + ((-t856 * t1092 + (t843 * t862 + (-t841 * t856 - t844 * t853) * t863) * t1006) * t903 + (-t856 * t1105 - (-t1141 * t838 + t1209 * t849 + t1210 * t856 + t843 * t876) * t1103 + t1047 * t844 - ((-t1003 * t843 - t1125 * t856 + t838) * t875 + (t846 * t849 + t853 * t856 + (-t846 * t856 - t849 * t853) * t1003) * t876) * t1102) * t902) * t858, 0; (t1054 * t889 + t1161 * t890) * t1177 + (t893 * t1174 + (t1054 * t890 + 0.2e1 * t1161 * t891) * t852 + ((t1000 * t1171 + t1211 * t1004) * t889 - (-t1211 * t1000 + t1004 * t1171) * t1162) * t1003 - t1059 * t1006 * t1076) * t872, -t854 * t1034 + (t1031 * t854 + t1036 * t840) * t1006, -t856 * t1034 + (t1031 * t856 + t1036 * t843) * t1006, -t1177 + (-0.2e1 * t872 * t1174 + (t1116 * t872 - t890 * t1177) * t894) * t894;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:43
	% EndTime: 2020-04-15 18:49:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:42
	% EndTime: 2020-04-15 18:49:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:40
	% EndTime: 2020-04-15 18:49:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiaD_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:43
	% EndTime: 2020-04-15 18:49:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiaD_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:50:55
	% EndTime: 2020-04-15 18:50:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
end