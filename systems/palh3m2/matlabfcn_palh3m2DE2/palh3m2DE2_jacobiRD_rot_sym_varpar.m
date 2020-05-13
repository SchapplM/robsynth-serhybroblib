% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh3m2DE2
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
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = palh3m2DE2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE2_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2DE2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_jacobiRD_rot_sym_varpar: pkin has to be [18x1] (double)');
JRD_rot=NaN(9,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:22
	% EndTime: 2020-05-07 02:13:22
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:23
	% EndTime: 2020-05-07 02:13:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0; -t31, 0, 0, 0; 0, 0, 0, 0; t31, 0, 0, 0; -t30, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:23
	% EndTime: 2020-05-07 02:13:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t32 = sin(qJ(1));
	t39 = qJD(1) * t32;
	t34 = cos(qJ(1));
	t38 = qJD(1) * t34;
	t31 = sin(qJ(2));
	t37 = qJD(2) * t31;
	t33 = cos(qJ(2));
	t36 = qJD(2) * t33;
	t35 = qJD(2) * t34;
	t30 = t32 * t37 - t33 * t38;
	t29 = t31 * t38 + t32 * t36;
	t28 = t31 * t35 + t33 * t39;
	t27 = t31 * t39 - t33 * t35;
	t1 = [t30, t27, 0, 0; -t28, -t29, 0, 0; 0, -t37, 0, 0; t29, t28, 0, 0; t27, t30, 0, 0; 0, -t36, 0, 0; -t39, 0, 0, 0; t38, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:24
	% EndTime: 2020-05-07 02:13:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (57->10), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t77 = qJD(2) + qJD(3);
	t79 = sin(qJ(1));
	t84 = t77 * t79;
	t80 = cos(qJ(1));
	t83 = t77 * t80;
	t82 = qJD(1) * t79;
	t81 = qJD(1) * t80;
	t78 = qJ(2) + qJ(3);
	t76 = cos(t78);
	t75 = sin(t78);
	t74 = t77 * t76;
	t73 = t77 * t75;
	t72 = -t75 * t84 + t76 * t81;
	t71 = t75 * t81 + t76 * t84;
	t70 = t75 * t83 + t76 * t82;
	t69 = -t75 * t82 + t76 * t83;
	t1 = [t72, t69, t69, 0; t70, t71, t71, 0; 0, t73, t73, 0; -t71, -t70, -t70, 0; t69, t72, t72, 0; 0, t74, t74, 0; -t82, 0, 0, 0; t81, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:27
	% EndTime: 2020-05-07 02:13:29
	% DurationCPUTime: 1.93s
	% Computational Cost: add. (18937->52), mult. (35154->90), div. (516->4), fcn. (50250->15), ass. (0->66)
	t312 = pkin(17) + pkin(18);
	t310 = sin(t312);
	t311 = cos(t312);
	t313 = sin(pkin(16));
	t314 = cos(pkin(16));
	t318 = sin(pkin(15));
	t322 = cos(pkin(15));
	t308 = t313 * t322 + t314 * t318;
	t309 = -t313 * t318 + t314 * t322;
	t315 = sin(qJ(3));
	t319 = cos(qJ(3));
	t302 = t315 * t308 - t309 * t319;
	t316 = sin(qJ(2));
	t320 = cos(qJ(2));
	t331 = t308 * t319 + t315 * t309;
	t372 = t302 * t316 - t331 * t320;
	t373 = t302 * t320 + t316 * t331;
	t375 = t373 * t310 + t311 * t372;
	t280 = t375 ^ 2;
	t380 = -t372 * t310 + t311 * t373;
	t282 = 0.1e1 / t380 ^ 2;
	t279 = t280 * t282 + 0.1e1;
	t277 = 0.1e1 / t279;
	t281 = 0.1e1 / t380;
	t353 = t282 * t375;
	t300 = t331 * qJD(3);
	t371 = t302 * qJD(3);
	t325 = qJD(2) * t373 + t316 * t300 + t320 * t371;
	t379 = -qJD(2) * t372 + t300 * t320 - t316 * t371;
	t376 = t310 * t379 + t311 * t325;
	t384 = -t325 * t310 + t379 * t311;
	t268 = qJD(2) + qJD(3) + (-t281 * t376 + t353 * t384) * t277;
	t330 = t281 * t380 + t375 * t353;
	t387 = -t330 * t277 + 0.1e1;
	t394 = t268 * t387;
	t274 = qJ(2) + qJ(3) + atan2(t375, -t380);
	t272 = sin(t274);
	t317 = sin(qJ(1));
	t339 = qJD(1) * t317;
	t273 = cos(t274);
	t321 = cos(qJ(1));
	t357 = t273 * t321;
	t326 = t268 * t357 - t272 * t339;
	t393 = t326 * t387;
	t338 = qJD(1) * t321;
	t358 = t273 * t317;
	t327 = t268 * t358 + t272 * t338;
	t392 = t327 * t387;
	t359 = t272 * t321;
	t328 = t268 * t359 + t273 * t339;
	t391 = t328 * t387;
	t360 = t272 * t317;
	t329 = -t268 * t360 + t273 * t338;
	t390 = t329 * t387;
	t389 = t273 * t394;
	t388 = t272 * t394;
	t370 = t384 * t282;
	t356 = t281 * t370;
	t381 = t376 * t353;
	t385 = 0.2e1 * t330 * (-t280 * t356 + t381) / t279 ^ 2;
	t383 = -0.2e1 * t375;
	t382 = t384 * t281;
	t335 = t356 * t383;
	t267 = ((-t282 * t376 - t335) * t375 - t382 + t380 * t370 - t381) * t277 + t385;
	t266 = t385 + (-t382 - t375 * t335 + (t376 * t383 + t380 * t384) * t282) * t277;
	t1 = [t329, t266 * t359 + t393, t267 * t359 + t393, 0; t328, t266 * t360 + t392, t267 * t360 + t392, 0; 0, -t266 * t273 + t388, -t267 * t273 + t388, 0; -t327, t266 * t357 - t391, t267 * t357 - t391, 0; t326, t266 * t358 + t390, t267 * t358 + t390, 0; 0, t266 * t272 + t389, t267 * t272 + t389, 0; -t339, 0, 0, 0; t338, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:48
	% EndTime: 2020-05-07 02:13:52
	% DurationCPUTime: 3.41s
	% Computational Cost: add. (32622->78), mult. (60268->144), div. (882->4), fcn. (86488->17), ass. (0->92)
	t901 = pkin(17) + pkin(18);
	t899 = sin(t901);
	t900 = cos(t901);
	t906 = sin(qJ(2));
	t911 = cos(qJ(2));
	t902 = sin(pkin(16));
	t903 = cos(pkin(16));
	t908 = sin(pkin(15));
	t913 = cos(pkin(15));
	t897 = t902 * t913 + t903 * t908;
	t898 = -t902 * t908 + t903 * t913;
	t905 = sin(qJ(3));
	t910 = cos(qJ(3));
	t928 = t897 * t910 + t905 * t898;
	t982 = -t897 * t905 + t898 * t910;
	t983 = t928 * t906 - t982 * t911;
	t992 = t982 * t906 + t911 * t928;
	t990 = t899 * t983 - t900 * t992;
	t869 = t990 ^ 2;
	t995 = t992 * t899 + t983 * t900;
	t871 = 0.1e1 / t995 ^ 2;
	t868 = t869 * t871 + 0.1e1;
	t870 = 0.1e1 / t995;
	t966 = t871 * t990;
	t925 = t870 * t995 + t990 * t966;
	t889 = t928 * qJD(3);
	t991 = t982 * qJD(3);
	t920 = t983 * qJD(2) + t889 * t906 - t911 * t991;
	t997 = t992 * qJD(2) + t889 * t911 + t906 * t991;
	t1002 = -t920 * t899 + t997 * t900;
	t986 = t1002 * t871;
	t969 = t870 * t986;
	t994 = t899 * t997 + t900 * t920;
	t998 = t994 * t966;
	t1013 = 0.2e1 * t925 * (-t869 * t969 + t998) / t868 ^ 2;
	t866 = 0.1e1 / t868;
	t1003 = -t925 * t866 + 0.1e1;
	t863 = qJ(2) + qJ(3) + atan2(t990, -t995);
	t861 = sin(t863);
	t904 = sin(qJ(4));
	t909 = cos(qJ(4));
	t857 = qJD(2) + qJD(3) + (t1002 * t966 - t870 * t994) * t866;
	t862 = cos(t863);
	t912 = cos(qJ(1));
	t970 = t862 * t912;
	t940 = t857 * t970;
	t947 = qJD(4) * t912;
	t907 = sin(qJ(1));
	t951 = qJD(1) * t907;
	t1012 = t1003 * (t909 * t940 + (-t904 * t947 - t909 * t951) * t861);
	t971 = t862 * t909;
	t941 = t857 * t971;
	t949 = qJD(4) * t904;
	t950 = qJD(1) * t912;
	t1011 = t1003 * (t907 * t941 + (-t907 * t949 + t909 * t950) * t861);
	t1010 = t1003 * (-t904 * t940 + (t904 * t951 - t909 * t947) * t861);
	t972 = t862 * t907;
	t942 = t857 * t972;
	t948 = qJD(4) * t909;
	t1009 = t1003 * (-t904 * t942 + (-t904 * t950 - t907 * t948) * t861);
	t978 = t857 * t861;
	t944 = t909 * t978;
	t1008 = t1003 * (t862 * t949 + t944);
	t1007 = t1003 * (t862 * t948 - t904 * t978);
	t974 = t861 * t912;
	t943 = t857 * t974;
	t1006 = t1003 * (t862 * t951 + t943);
	t975 = t861 * t907;
	t945 = t857 * t975;
	t1005 = t1003 * (-t862 * t950 + t945);
	t1004 = t1003 * t857 * t862;
	t1000 = -t1002 * t870 + t995 * t986;
	t999 = t994 * t871;
	t931 = qJD(1) * t862 + qJD(4);
	t979 = t931 * t907 + t943;
	t973 = t862 * t904;
	t953 = t907 * t909;
	t939 = t904 * t975;
	t938 = t904 * t974;
	t937 = t861 * t953;
	t936 = t909 * t974;
	t934 = 0.2e1 * t990 * t969;
	t932 = qJD(4) * t862 + qJD(1);
	t927 = t932 * t912;
	t926 = t931 * t912;
	t856 = t909 * t926 + (-t932 * t904 - t944) * t907;
	t855 = t932 * t953 + (t926 - t945) * t904;
	t854 = t904 * t927 + t979 * t909;
	t853 = -t979 * t904 + t909 * t927;
	t852 = ((t934 - t999) * t990 - t998 + t1000) * t866 + t1013;
	t851 = t1013 + ((t934 - 0.2e1 * t999) * t990 + t1000) * t866;
	t1 = [t856, t851 * t936 + t1012, t852 * t936 + t1012, t853; t854, t851 * t937 + t1011, t852 * t937 + t1011, t855; 0, -t851 * t971 + t1008, -t852 * t971 + t1008, t857 * t973 + t861 * t948; -t855, -t851 * t938 + t1010, -t852 * t938 + t1010, -t854; t853, -t851 * t939 + t1009, -t852 * t939 + t1009, t856; 0, t851 * t973 + t1007, t852 * t973 + t1007, -t861 * t949 + t941; t861 * t950 + t942, -t851 * t970 + t1006, -t852 * t970 + t1006, 0; t861 * t951 - t940, -t851 * t972 + t1005, -t852 * t972 + t1005, 0; 0, -t851 * t861 - t1004, -t852 * t861 - t1004, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:52
	% EndTime: 2020-05-07 02:13:52
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (63->11), mult. (154->20), div. (0->0), fcn. (198->8), ass. (0->21)
	t78 = sin(pkin(15));
	t79 = sin(pkin(14));
	t82 = cos(pkin(15));
	t83 = cos(pkin(14));
	t74 = t78 * t83 - t82 * t79;
	t75 = t78 * t79 + t82 * t83;
	t76 = sin(qJ(2));
	t80 = cos(qJ(2));
	t86 = t76 * t74 - t75 * t80;
	t68 = t86 * qJD(2);
	t77 = sin(qJ(1));
	t88 = qJD(1) * t77;
	t81 = cos(qJ(1));
	t87 = qJD(1) * t81;
	t69 = -t80 * t74 - t76 * t75;
	t85 = t77 * t68 + t69 * t87;
	t84 = t81 * t68 - t69 * t88;
	t67 = t69 * qJD(2);
	t66 = -t67 * t77 + t86 * t87;
	t65 = -t81 * t67 - t86 * t88;
	t1 = [t66, t84, 0, 0; -t65, t85, 0, 0; 0, t67, 0, 0; -t85, t65, 0, 0; t84, t66, 0, 0; 0, t68, 0, 0; -t88, 0, 0, 0; t87, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiRD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:53
	% EndTime: 2020-05-07 02:13:53
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (1553->27), mult. (3288->56), div. (252->4), fcn. (4590->11), ass. (0->41)
	t120 = sin(pkin(18));
	t121 = cos(pkin(18));
	t124 = sin(pkin(15));
	t127 = cos(pkin(15));
	t118 = t127 * t120 + t124 * t121;
	t119 = -t124 * t120 + t127 * t121;
	t122 = sin(qJ(2));
	t125 = cos(qJ(2));
	t115 = t122 * t118 - t125 * t119;
	t112 = 0.1e1 / t115 ^ 2;
	t116 = t125 * t118 + t122 * t119;
	t138 = t116 ^ 2 * t112;
	t111 = 0.1e1 / t115;
	t108 = 0.1e1 + t138;
	t106 = 0.1e1 / t108;
	t109 = t115 * qJD(2);
	t110 = t116 * qJD(2);
	t140 = t112 * t116;
	t100 = qJD(2) + (-t109 * t111 - t110 * t140) * t106;
	t134 = t111 * t115 + t138;
	t101 = -t106 * t134 + 0.1e1;
	t145 = t100 * t101;
	t105 = qJ(2) + atan2(t116, t115);
	t103 = sin(t105);
	t123 = sin(qJ(1));
	t144 = t103 * t123;
	t126 = cos(qJ(1));
	t143 = t103 * t126;
	t104 = cos(t105);
	t142 = t104 * t123;
	t141 = t104 * t126;
	t139 = t111 * t138;
	t137 = qJD(1) * t123;
	t136 = qJD(1) * t126;
	t135 = t109 * t140;
	t133 = t100 * t144 - t104 * t136;
	t132 = t100 * t143 + t104 * t137;
	t131 = t100 * t142 + t103 * t136;
	t130 = -t100 * t141 + t103 * t137;
	t99 = 0.2e1 * t134 / t108 ^ 2 * (-t110 * t139 - t135) + (0.2e1 * t135 + (t112 * t115 - t111 + 0.2e1 * t139) * t110) * t106;
	t1 = [t133, t101 * t130 - t143 * t99, 0, 0; -t132, -t101 * t131 - t144 * t99, 0, 0; 0, -t103 * t145 + t99 * t104, 0, 0; t131, t101 * t132 - t141 * t99, 0, 0; t130, t101 * t133 - t142 * t99, 0, 0; 0, -t99 * t103 - t104 * t145, 0, 0; -t137, 0, 0, 0; t136, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiRD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:56
	% EndTime: 2020-05-07 02:13:57
	% DurationCPUTime: 1.29s
	% Computational Cost: add. (11985->59), mult. (21426->114), div. (804->8), fcn. (29366->16), ass. (0->70)
	t344 = pkin(17) + pkin(18);
	t342 = sin(t344);
	t343 = cos(t344);
	t349 = sin(pkin(15));
	t350 = cos(qJ(3));
	t353 = cos(pkin(15));
	t392 = sin(qJ(3));
	t338 = t350 * t349 + t392 * t353;
	t347 = sin(qJ(2));
	t351 = cos(qJ(2));
	t396 = t392 * t349 - t350 * t353;
	t397 = t347 * t338 + t351 * t396;
	t398 = t338 * t351 - t347 * t396;
	t403 = t342 * t398 + t343 * t397;
	t301 = 0.1e1 / t403 ^ 2;
	t404 = -t342 * t397 + t343 * t398;
	t409 = t404 ^ 2 * t301;
	t334 = t396 * qJD(3);
	t335 = t338 * qJD(3);
	t312 = qJD(2) * t398 - t334 * t347 + t335 * t351;
	t394 = qJD(2) * t397 + t334 * t351 + t335 * t347;
	t294 = t312 * t343 - t342 * t394;
	t300 = 0.1e1 / t403;
	t402 = t294 * t300;
	t408 = t402 * t409;
	t297 = 0.1e1 + t409;
	t295 = 0.1e1 / t297;
	t384 = t301 * t404;
	t407 = -t300 * t403 - t384 * t404;
	t287 = t407 * t295;
	t293 = t312 * t342 + t394 * t343;
	t284 = -0.2e1 * t407 * (-t293 * t384 - t408) / t297 ^ 2 + (-t402 + 0.2e1 * t408 + (0.2e1 * t293 * t404 + t294 * t403) * t301) * t295;
	t345 = sin(pkin(18));
	t346 = cos(pkin(18));
	t336 = t353 * t345 + t349 * t346;
	t337 = -t349 * t345 + t353 * t346;
	t323 = t347 * t336 - t351 * t337;
	t320 = 0.1e1 / t323 ^ 2;
	t324 = t351 * t336 + t347 * t337;
	t381 = t324 ^ 2 * t320;
	t319 = 0.1e1 / t323;
	t315 = 0.1e1 + t381;
	t313 = 0.1e1 / t315;
	t317 = t323 * qJD(2);
	t318 = t324 * qJD(2);
	t383 = t320 * t324;
	t285 = -qJD(2) + (t317 * t319 + t318 * t383) * t313 + (-t293 * t300 - t294 * t384) * t295;
	t291 = -qJ(2) - atan2(t324, t323) + pkin(17) - atan2(-t404, t403);
	t289 = sin(t291);
	t391 = t285 * t289;
	t290 = cos(t291);
	t390 = t285 * t290;
	t348 = sin(qJ(1));
	t389 = t289 * t348;
	t352 = cos(qJ(1));
	t388 = t289 * t352;
	t387 = t290 * t348;
	t386 = t290 * t352;
	t382 = t319 * t381;
	t374 = qJD(1) * t348;
	t373 = qJD(1) * t352;
	t372 = t317 * t383;
	t367 = t319 * t323 + t381;
	t362 = -t285 * t389 + t290 * t373;
	t361 = t285 * t388 + t290 * t374;
	t360 = t285 * t387 + t289 * t373;
	t359 = -t285 * t386 + t289 * t374;
	t286 = t367 * t313 + t287 - 0.1e1;
	t283 = (-0.2e1 * t372 + (-t320 * t323 + t319 - 0.2e1 * t382) * t318) * t313 - 0.2e1 * t367 / t315 ^ 2 * (-t318 * t382 - t372) + t284;
	t1 = [t362, t283 * t388 - t359 * t286, t284 * t388 - t359 * t287, 0; t361, t283 * t389 + t360 * t286, t284 * t389 + t360 * t287, 0; 0, t283 * t290 - t286 * t391, t284 * t290 - t287 * t391, 0; t360, -t283 * t386 + t361 * t286, -t284 * t386 + t361 * t287, 0; t359, -t283 * t387 - t362 * t286, -t284 * t387 - t362 * t287, 0; 0, t283 * t289 + t286 * t390, t284 * t289 + t287 * t390, 0; -t374, 0, 0, 0; t373, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
end