% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = palh1m2DE2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m2DE2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_jacobiRD_rot_sym_varpar: pkin has to be [22x1] (double)');
JRD_rot=NaN(9,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0; -t31, 0, 0, 0; 0, 0, 0, 0; t31, 0, 0, 0; -t30, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->8), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t32 = sin(qJ(1));
	t39 = qJD(1) * t32;
	t34 = cos(qJ(1));
	t38 = qJD(1) * t34;
	t31 = sin(qJ(2));
	t37 = qJD(2) * t31;
	t33 = cos(qJ(2));
	t36 = qJD(2) * t33;
	t35 = qJD(2) * t34;
	t30 = -t32 * t37 + t33 * t38;
	t29 = t31 * t38 + t32 * t36;
	t28 = t31 * t35 + t33 * t39;
	t27 = t31 * t39 - t33 * t35;
	t1 = [t29, t28, 0, 0; t27, -t30, 0, 0; 0, -t36, 0, 0; t30, -t27, 0, 0; t28, t29, 0, 0; 0, t37, 0, 0; -t39, 0, 0, 0; t38, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:29
	% EndTime: 2020-05-02 21:08:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (61->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t68 = qJ(2) + qJ(3);
	t65 = sin(t68);
	t67 = qJD(2) + qJD(3);
	t76 = t67 * t65;
	t66 = cos(t68);
	t75 = t67 * t66;
	t69 = sin(qJ(1));
	t74 = t67 * t69;
	t70 = cos(qJ(1));
	t73 = t67 * t70;
	t72 = qJD(1) * t69;
	t71 = qJD(1) * t70;
	t64 = t65 * t74 - t66 * t71;
	t63 = t65 * t71 + t66 * t74;
	t62 = t65 * t73 + t66 * t72;
	t61 = t65 * t72 - t66 * t73;
	t1 = [t64, t61, t61, 0; -t62, -t63, -t63, 0; 0, -t76, -t76, 0; t63, t62, t62, 0; t61, t64, t64, 0; 0, -t75, -t75, 0; -t72, 0, 0, 0; t71, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:34
	% EndTime: 2020-05-02 21:08:36
	% DurationCPUTime: 2.11s
	% Computational Cost: add. (18937->50), mult. (35154->89), div. (516->4), fcn. (50250->15), ass. (0->65)
	t326 = pkin(22) + pkin(21);
	t324 = sin(t326);
	t325 = cos(t326);
	t327 = cos(pkin(20));
	t331 = sin(pkin(18));
	t374 = sin(pkin(20));
	t375 = cos(pkin(18));
	t321 = t327 * t331 - t375 * t374;
	t322 = t375 * t327 + t331 * t374;
	t328 = sin(qJ(3));
	t332 = cos(qJ(3));
	t313 = t321 * t332 - t322 * t328;
	t329 = sin(qJ(2));
	t333 = cos(qJ(2));
	t378 = t328 * t321 + t322 * t332;
	t344 = t313 * t333 - t378 * t329;
	t394 = t329 * t313 + t333 * t378;
	t398 = t344 * t324 + t325 * t394;
	t288 = 0.1e1 / t398;
	t289 = 0.1e1 / t398 ^ 2;
	t376 = t378 * qJD(3);
	t388 = t313 * qJD(3);
	t397 = qJD(2) * t394 + t329 * t388 + t376 * t333;
	t404 = t344 * qJD(2) - t329 * t376 + t388 * t333;
	t408 = -t397 * t324 + t404 * t325;
	t382 = t408 * t289;
	t366 = t288 * t382;
	t399 = t394 * t324 - t344 * t325;
	t403 = -0.2e1 * t399;
	t416 = -t399 * t366 * t403 - t408 * t288;
	t282 = t324 * t404 + t397 * t325;
	t287 = t399 ^ 2;
	t286 = t287 * t289 + 0.1e1;
	t284 = 0.1e1 / t286;
	t396 = t289 * t399;
	t275 = qJD(2) + qJD(3) + (-t282 * t288 + t396 * t408) * t284;
	t341 = t288 * t398 + t399 * t396;
	t407 = -t341 * t284 + 0.1e1;
	t415 = t275 * t407;
	t281 = qJ(2) + qJ(3) + atan2(t399, -t398);
	t279 = sin(t281);
	t330 = sin(qJ(1));
	t350 = qJD(1) * t330;
	t280 = cos(t281);
	t334 = cos(qJ(1));
	t367 = t280 * t334;
	t337 = -t275 * t367 + t279 * t350;
	t414 = t337 * t407;
	t349 = qJD(1) * t334;
	t368 = t280 * t330;
	t338 = t275 * t368 + t279 * t349;
	t413 = t338 * t407;
	t369 = t279 * t334;
	t339 = t275 * t369 + t280 * t350;
	t412 = t339 * t407;
	t370 = t279 * t330;
	t340 = t275 * t370 - t280 * t349;
	t411 = t340 * t407;
	t410 = t280 * t415;
	t409 = t279 * t415;
	t392 = t282 * t396;
	t405 = 0.2e1 * t341 * (-t287 * t366 + t392) / t286 ^ 2;
	t274 = (t398 * t382 - 0.2e1 * t392 + t416) * t284 + t405;
	t273 = t405 + ((t282 * t403 + t398 * t408) * t289 + t416) * t284;
	t1 = [t340, -t273 * t369 + t414, -t274 * t369 + t414, 0; -t339, -t273 * t370 - t413, -t274 * t370 - t413, 0; 0, t273 * t280 - t409, t274 * t280 - t409, 0; t338, -t273 * t367 + t412, -t274 * t367 + t412, 0; t337, -t273 * t368 + t411, -t274 * t368 + t411, 0; 0, -t273 * t279 - t410, -t274 * t279 - t410, 0; -t350, 0, 0, 0; t349, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:50
	% EndTime: 2020-05-02 21:08:54
	% DurationCPUTime: 3.71s
	% Computational Cost: add. (32622->77), mult. (60268->146), div. (882->4), fcn. (86488->17), ass. (0->90)
	t857 = pkin(22) + pkin(21);
	t855 = sin(t857);
	t856 = cos(t857);
	t858 = cos(pkin(20));
	t863 = sin(pkin(18));
	t931 = sin(pkin(20));
	t932 = cos(pkin(18));
	t852 = t863 * t858 - t932 * t931;
	t853 = t932 * t858 + t863 * t931;
	t860 = sin(qJ(3));
	t865 = cos(qJ(3));
	t844 = t852 * t865 - t860 * t853;
	t861 = sin(qJ(2));
	t866 = cos(qJ(2));
	t935 = t860 * t852 + t853 * t865;
	t883 = t844 * t866 - t861 * t935;
	t951 = t861 * t844 + t866 * t935;
	t955 = t883 * t855 + t856 * t951;
	t819 = 0.1e1 / t955;
	t820 = 0.1e1 / t955 ^ 2;
	t933 = t935 * qJD(3);
	t945 = t844 * qJD(3);
	t954 = qJD(2) * t951 + t861 * t945 + t933 * t866;
	t961 = t883 * qJD(2) - t861 * t933 + t945 * t866;
	t965 = -t954 * t855 + t961 * t856;
	t939 = t965 * t820;
	t921 = t819 * t939;
	t956 = t951 * t855 - t883 * t856;
	t960 = -0.2e1 * t956;
	t975 = -t956 * t921 * t960 - t965 * t819;
	t812 = qJ(2) + qJ(3) + atan2(t956, -t955);
	t810 = sin(t812);
	t859 = sin(qJ(4));
	t864 = cos(qJ(4));
	t813 = t855 * t961 + t954 * t856;
	t818 = t956 ^ 2;
	t817 = t818 * t820 + 0.1e1;
	t815 = 0.1e1 / t817;
	t953 = t820 * t956;
	t806 = qJD(2) + qJD(3) + (-t813 * t819 + t953 * t965) * t815;
	t811 = cos(t812);
	t867 = cos(qJ(1));
	t922 = t811 * t867;
	t887 = t806 * t922;
	t898 = qJD(4) * t867;
	t862 = sin(qJ(1));
	t903 = qJD(1) * t862;
	t879 = t819 * t955 + t956 * t953;
	t964 = -t879 * t815 + 0.1e1;
	t974 = (-t864 * t887 + (t859 * t898 + t864 * t903) * t810) * t964;
	t924 = t811 * t862;
	t888 = t806 * t924;
	t900 = qJD(4) * t859;
	t902 = qJD(1) * t867;
	t973 = (-t864 * t888 + (t862 * t900 - t864 * t902) * t810) * t964;
	t972 = (t859 * t887 + (-t859 * t903 + t864 * t898) * t810) * t964;
	t899 = qJD(4) * t864;
	t971 = (t859 * t888 + (t859 * t902 + t862 * t899) * t810) * t964;
	t930 = t806 * t810;
	t895 = t864 * t930;
	t970 = (-t811 * t900 - t895) * t964;
	t969 = (-t811 * t899 + t859 * t930) * t964;
	t926 = t810 * t867;
	t894 = t806 * t926;
	t968 = (-t811 * t903 - t894) * t964;
	t927 = t810 * t862;
	t896 = t806 * t927;
	t967 = (t811 * t902 - t896) * t964;
	t966 = t806 * t811 * t964;
	t949 = t813 * t953;
	t962 = 0.2e1 * t879 * (-t818 * t921 + t949) / t817 ^ 2;
	t925 = t811 * t859;
	t923 = t811 * t864;
	t906 = t862 * t864;
	t905 = t864 * t867;
	t893 = t859 * t927;
	t892 = t859 * t926;
	t891 = t810 * t906;
	t890 = t810 * t905;
	t885 = qJD(4) * t811 - qJD(1);
	t884 = qJD(1) * t811 - qJD(4);
	t881 = t885 * t859;
	t874 = t884 * t862 + t894;
	t805 = -t884 * t905 + (t881 + t895) * t862;
	t804 = t885 * t906 + (t884 * t867 - t896) * t859;
	t803 = t874 * t864 + t867 * t881;
	t802 = t874 * t859 - t885 * t905;
	t801 = (t955 * t939 - 0.2e1 * t949 + t975) * t815 + t962;
	t800 = t962 + ((t813 * t960 + t955 * t965) * t820 + t975) * t815;
	t1 = [t805, -t800 * t890 + t974, -t801 * t890 + t974, t802; -t803, -t800 * t891 + t973, -t801 * t891 + t973, -t804; 0, t800 * t923 + t970, t801 * t923 + t970, -t806 * t925 - t810 * t899; t804, t800 * t892 + t972, t801 * t892 + t972, t803; t802, t800 * t893 + t971, t801 * t893 + t971, t805; 0, -t800 * t925 + t969, -t801 * t925 + t969, -t806 * t923 + t810 * t900; -t810 * t902 - t888, t800 * t922 + t968, t801 * t922 + t968, 0; -t810 * t903 + t887, t800 * t924 + t967, t801 * t924 + t967, 0; 0, t800 * t810 + t966, t801 * t810 + t966, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:29
	% EndTime: 2020-05-02 21:08:29
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (63->12), mult. (154->20), div. (0->0), fcn. (198->8), ass. (0->21)
	t76 = sin(pkin(18));
	t77 = sin(pkin(17));
	t80 = cos(pkin(18));
	t81 = cos(pkin(17));
	t72 = t76 * t81 - t80 * t77;
	t73 = t77 * t76 + t81 * t80;
	t74 = sin(qJ(2));
	t78 = cos(qJ(2));
	t68 = t72 * t78 - t74 * t73;
	t66 = t68 * qJD(2);
	t75 = sin(qJ(1));
	t85 = qJD(1) * t75;
	t79 = cos(qJ(1));
	t84 = qJD(1) * t79;
	t67 = -t74 * t72 - t73 * t78;
	t83 = -t75 * t66 + t67 * t84;
	t82 = -t79 * t66 - t67 * t85;
	t65 = t67 * qJD(2);
	t64 = -t65 * t75 - t68 * t84;
	t63 = -t79 * t65 + t68 * t85;
	t1 = [t64, t82, 0, 0; -t63, t83, 0, 0; 0, t65, 0, 0; -t83, t63, 0, 0; t82, t64, 0, 0; 0, -t66, 0, 0; -t85, 0, 0, 0; t84, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiRD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:29
	% EndTime: 2020-05-02 21:08:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->8), mult. (18->12), div. (0->0), fcn. (10->6), ass. (0->8)
	t45 = pkin(18) - pkin(22);
	t43 = qJ(1) + t45;
	t47 = qJD(1) * sin(t43);
	t44 = -qJ(1) + t45;
	t46 = qJD(1) * cos(t44);
	t40 = qJD(1) * cos(t43) / 0.2e1;
	t39 = -qJD(1) * sin(t44) / 0.2e1;
	t1 = [t46 / 0.2e1 + t40, 0, 0, 0; t39 + t47 / 0.2e1, 0, 0, 0; 0, 0, 0, 0; -t47 / 0.2e1 + t39, 0, 0, 0; t40 - t46 / 0.2e1, 0, 0, 0; 0, 0, 0, 0; -qJD(1) * sin(qJ(1)), 0, 0, 0; qJD(1) * cos(qJ(1)), 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiRD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:29
	% EndTime: 2020-05-02 21:08:30
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (706->26), mult. (1602->52), div. (254->4), fcn. (2066->9), ass. (0->41)
	t110 = cos(pkin(19));
	t112 = cos(qJ(3));
	t126 = sin(pkin(19));
	t132 = sin(qJ(3));
	t106 = t132 * t110 + t112 * t126;
	t103 = 0.1e1 / t106 ^ 2;
	t107 = t112 * t110 - t132 * t126;
	t105 = t107 ^ 2;
	t100 = t105 * t103 + 0.1e1;
	t133 = 0.1e1 / t106;
	t101 = t106 * qJD(3);
	t102 = t107 * qJD(3);
	t124 = t103 * t107;
	t98 = 0.1e1 / t100;
	t92 = qJD(2) + (t101 * t133 + t102 * t124) * t98;
	t97 = qJ(2) + atan2(-t107, t106);
	t95 = sin(t97);
	t131 = t92 * t95;
	t96 = cos(t97);
	t130 = t92 * t96;
	t125 = t102 * t133 * t103;
	t129 = (-t101 * t124 - t105 * t125) / t100 ^ 2;
	t111 = sin(qJ(1));
	t128 = t111 * t92;
	t113 = cos(qJ(1));
	t127 = t113 * t92;
	t123 = qJD(1) * t111;
	t122 = qJD(1) * t113;
	t121 = t95 * t123;
	t120 = t96 * t123;
	t119 = t95 * t122;
	t118 = t96 * t122;
	t87 = -0.2e1 * t129 + 0.2e1 * (-t101 * t103 * t98 + (-t103 * t129 - t98 * t125) * t107) * t107;
	t93 = t100 * t98;
	t117 = t93 * t130 + t87 * t95;
	t116 = t93 * t131 - t87 * t96;
	t91 = -t95 * t128 + t118;
	t90 = t96 * t128 + t119;
	t89 = t95 * t127 + t120;
	t88 = -t96 * t127 + t121;
	t1 = [t90, t89, t116 * t113 + t93 * t120, 0; t88, -t91, t116 * t111 - t93 * t118, 0; 0, -t130, -t117, 0; t91, -t88, t117 * t113 - t93 * t121, 0; t89, t90, t117 * t111 + t93 * t119, 0; 0, t131, t116, 0; -t123, 0, 0, 0; t122, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiRD_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:30
	% EndTime: 2020-05-02 21:08:31
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (1454->32), mult. (3234->57), div. (542->4), fcn. (4204->10), ass. (0->40)
	t198 = cos(pkin(19));
	t200 = cos(qJ(3));
	t222 = sin(pkin(19));
	t223 = sin(qJ(3));
	t192 = t198 * t223 + t200 * t222;
	t225 = 0.1e1 / t192 ^ 2;
	t233 = t192 * t225;
	t232 = qJD(3) * t233;
	t194 = t198 * t200 - t222 * t223;
	t191 = t194 ^ 2;
	t211 = t191 * t225;
	t229 = 0.1e1 + t211;
	t231 = 0.1e1 / t229;
	t207 = t194 * t232;
	t188 = 0.1e1 / t192;
	t224 = t188 * t225;
	t177 = qJ(2) + atan2(-t194, t192) + atan2(-t194, -t192);
	t175 = sin(t177);
	t221 = qJD(2) * t175;
	t176 = cos(t177);
	t220 = qJD(2) * t176;
	t199 = sin(qJ(1));
	t219 = t175 * t199;
	t201 = cos(qJ(1));
	t218 = t175 * t201;
	t217 = t176 * t199;
	t216 = t176 * t201;
	t185 = t194 * qJD(3);
	t215 = t185 * t224;
	t212 = t224 * t191;
	t209 = qJD(1) * t199;
	t208 = qJD(1) * t201;
	t204 = t188 * t192 + t211;
	t170 = -qJD(2) * t219 + t176 * t208;
	t168 = -qJD(2) * t218 - t176 * t209;
	t169 = -qJD(2) * t217 - t175 * t208;
	t167 = -qJD(2) * t216 + t175 * t209;
	t172 = (-t204 + t229) * t231;
	t166 = 0.2e1 * t204 / t229 ^ 2 * (-t185 * t212 - t207) + (0.2e1 * (-t194 * t215 - t232) * t194 + 0.4e1 * t207 + (-t188 + 0.2e1 * t212 + t233) * t185 + 0.2e1 * t191 * t215) * t231;
	t1 = [t169, t168, t166 * t216 + t168 * t172, 0; -t167, t170, t166 * t217 + t170 * t172, 0; 0, t220, t166 * t175 + t172 * t220, 0; -t170, t167, -t166 * t218 + t167 * t172, 0; t168, t169, -t166 * t219 + t169 * t172, 0; 0, -t221, t166 * t176 - t172 * t221, 0; -t209, 0, 0, 0; t208, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiRD_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:33
	% EndTime: 2020-05-02 21:08:35
	% DurationCPUTime: 1.42s
	% Computational Cost: add. (12153->61), mult. (21768->116), div. (822->8), fcn. (29886->16), ass. (0->70)
	t370 = pkin(22) + pkin(21);
	t368 = sin(t370);
	t369 = cos(t370);
	t415 = sin(qJ(3));
	t416 = sin(pkin(18));
	t417 = cos(qJ(3));
	t418 = cos(pkin(18));
	t362 = -t415 * t416 - t417 * t418;
	t363 = -t415 * t418 + t417 * t416;
	t372 = sin(qJ(2));
	t374 = cos(qJ(2));
	t390 = t362 * t374 - t372 * t363;
	t420 = t372 * t362 + t374 * t363;
	t437 = t368 * t420 - t390 * t369;
	t321 = 0.1e1 / t437 ^ 2;
	t430 = t368 * t390 + t369 * t420;
	t440 = t430 ^ 2 * t321;
	t358 = t362 * qJD(3);
	t359 = t363 * qJD(3);
	t380 = t390 * qJD(2) + t358 * t374 - t372 * t359;
	t381 = qJD(2) * t420 + t358 * t372 + t359 * t374;
	t314 = t368 * t380 + t381 * t369;
	t320 = 0.1e1 / t437;
	t438 = t314 * t320;
	t439 = t438 * t440;
	t317 = 0.1e1 + t440;
	t315 = 0.1e1 / t317;
	t406 = t321 * t430;
	t435 = -t320 * t437 - t406 * t430;
	t307 = t435 * t315;
	t436 = -t368 * t381 + t369 * t380;
	t304 = -0.2e1 * t435 * (t406 * t436 - t439) / t317 ^ 2 + (-t438 + 0.2e1 * t439 + (t314 * t437 - 0.2e1 * t430 * t436) * t321) * t315;
	t371 = sin(pkin(22));
	t414 = cos(pkin(22));
	t360 = t416 * t371 + t414 * t418;
	t361 = t418 * t371 - t416 * t414;
	t345 = t372 * t360 + t361 * t374;
	t341 = 0.1e1 / t345 ^ 2;
	t421 = t360 * t374 - t361 * t372;
	t433 = t341 * t421 ^ 2;
	t340 = 0.1e1 / t345;
	t432 = t340 * t433;
	t339 = t345 * qJD(2);
	t405 = t341 * t421;
	t426 = t339 * t405;
	t336 = 0.1e1 + t433;
	t334 = 0.1e1 / t336;
	t338 = t421 * qJD(2);
	t305 = -qJD(2) + (t338 * t405 + t339 * t340) * t334 + (-t314 * t406 + t320 * t436) * t315;
	t311 = -qJ(2) - atan2(t421, t345) + pkin(21) - atan2(-t430, t437);
	t309 = sin(t311);
	t413 = t305 * t309;
	t310 = cos(t311);
	t412 = t305 * t310;
	t373 = sin(qJ(1));
	t411 = t309 * t373;
	t375 = cos(qJ(1));
	t410 = t309 * t375;
	t409 = t310 * t373;
	t408 = t310 * t375;
	t397 = qJD(1) * t373;
	t396 = qJD(1) * t375;
	t386 = t340 * t345 + t433;
	t385 = t305 * t411 - t310 * t396;
	t384 = t305 * t410 + t310 * t397;
	t383 = t305 * t409 + t309 * t396;
	t382 = -t305 * t408 + t309 * t397;
	t306 = t386 * t334 + t307 - 0.1e1;
	t303 = (-t426 + (-t341 * t345 - 0.2e1 * t432) * t338 + (t340 * t421 - t345 * t405) * qJD(2)) * t334 - 0.2e1 * t386 / t336 ^ 2 * (-t338 * t432 - t426) + t304;
	t1 = [t383, -t303 * t408 + t384 * t306, -t304 * t408 + t384 * t307, 0; t382, -t303 * t409 + t385 * t306, -t304 * t409 + t385 * t307, 0; 0, t303 * t309 + t306 * t412, t304 * t309 + t307 * t412, 0; t385, -t303 * t410 + t382 * t306, -t304 * t410 + t382 * t307, 0; -t384, -t303 * t411 - t383 * t306, -t304 * t411 - t383 * t307, 0; 0, -t303 * t310 + t306 * t413, -t304 * t310 + t307 * t413, 0; -t397, 0, 0, 0; t396, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
end