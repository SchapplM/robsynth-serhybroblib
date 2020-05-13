% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh1m1OL
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% 
% Output:
% JaD_transl [3x13]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = palh1m1OL_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),uint8(0),zeros(3,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_jacobiaD_transl_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1OL_jacobiaD_transl_sym_varpar: qJD has to be [13x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh1m1OL_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1OL_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_jacobiaD_transl_sym_varpar: pkin has to be [20x1] (double)');
JaD_transl=NaN(3,13);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:18
	% EndTime: 2020-04-15 19:46:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:18
	% EndTime: 2020-04-15 19:46:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:18
	% EndTime: 2020-04-15 19:46:18
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (17->14), mult. (60->29), div. (0->0), fcn. (38->4), ass. (0->12)
	t17 = sin(qJ(1));
	t26 = qJD(1) * t17;
	t19 = cos(qJ(1));
	t25 = qJD(1) * t19;
	t24 = qJD(2) * t17;
	t23 = qJD(2) * t19;
	t16 = sin(qJ(2));
	t18 = cos(qJ(2));
	t22 = r_i_i_C(1) * t18 - r_i_i_C(2) * t16;
	t21 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18 - pkin(15);
	t20 = t22 * qJD(2);
	t1 = [t22 * t24 + (-r_i_i_C(3) * t17 + t21 * t19) * qJD(1), (-t16 * t26 + t18 * t23) * r_i_i_C(2) + (t16 * t23 + t18 * t26) * r_i_i_C(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t19 * t20 + (r_i_i_C(3) * t19 + t21 * t17) * qJD(1), (t16 * t25 + t18 * t24) * r_i_i_C(2) + (t16 * t24 - t18 * t25) * r_i_i_C(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:18
	% EndTime: 2020-04-15 19:46:19
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (77->25), mult. (110->37), div. (0->0), fcn. (71->6), ass. (0->26)
	t37 = qJD(2) + qJD(3);
	t38 = qJ(2) + qJ(3);
	t36 = cos(t38);
	t56 = r_i_i_C(2) * t36;
	t35 = sin(t38);
	t58 = r_i_i_C(1) * t35;
	t47 = t56 + t58;
	t45 = t47 * t37;
	t41 = cos(qJ(2));
	t59 = pkin(1) * t41;
	t60 = qJD(2) * t59 + t45;
	t57 = r_i_i_C(2) * t35;
	t55 = t36 * t37;
	t40 = sin(qJ(1));
	t54 = qJD(1) * t40;
	t42 = cos(qJ(1));
	t53 = qJD(1) * t42;
	t39 = sin(qJ(2));
	t52 = qJD(2) * t39;
	t51 = r_i_i_C(1) * t55;
	t50 = t37 * t57;
	t49 = qJD(1) * t56;
	t46 = t39 * pkin(1) - r_i_i_C(1) * t36 - pkin(15) + t57;
	t44 = t40 * t49 + t54 * t58 + (t50 - t51) * t42;
	t30 = t40 * t50;
	t1 = [t60 * t40 + (-r_i_i_C(3) * t40 + t46 * t42) * qJD(1), (t41 * t54 + t42 * t52) * pkin(1) + t44, t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t60 * t42 + (r_i_i_C(3) * t42 + t46 * t40) * qJD(1), t30 + (pkin(1) * t52 - t51) * t40 + (-t47 - t59) * t53, -t42 * t49 + t30 + (-t35 * t53 - t40 * t55) * r_i_i_C(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t60, -t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:18
	% EndTime: 2020-04-15 19:46:19
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (196->35), mult. (166->47), div. (0->0), fcn. (107->8), ass. (0->36)
	t54 = cos(qJ(2));
	t50 = qJD(2) + qJD(3);
	t46 = qJD(4) + t50;
	t51 = qJ(2) + qJ(3);
	t49 = qJ(4) + t51;
	t45 = cos(t49);
	t74 = r_i_i_C(2) * t45;
	t44 = sin(t49);
	t76 = r_i_i_C(1) * t44;
	t60 = t74 + t76;
	t58 = t60 * t46;
	t47 = sin(t51);
	t77 = pkin(5) * t47;
	t56 = -t50 * t77 - t58;
	t71 = pkin(1) * qJD(2);
	t79 = -t54 * t71 + t56;
	t73 = t45 * t46;
	t66 = r_i_i_C(1) * t73;
	t48 = cos(t51);
	t72 = t48 * t50;
	t78 = -pkin(5) * t72 - t66;
	t75 = r_i_i_C(2) * t44;
	t53 = sin(qJ(1));
	t70 = qJD(1) * t53;
	t55 = cos(qJ(1));
	t69 = qJD(1) * t55;
	t65 = t46 * t75;
	t63 = qJD(1) * t74;
	t64 = t53 * t63 + t55 * t65 + t70 * t76;
	t52 = sin(qJ(2));
	t61 = t52 * t71 + t78;
	t59 = t52 * pkin(1) - pkin(5) * t48 - r_i_i_C(1) * t45 - pkin(15) + t75;
	t57 = -t55 * t66 + t64;
	t43 = -t54 * pkin(1) - t77;
	t38 = t53 * t65;
	t1 = [-t79 * t53 + (-r_i_i_C(3) * t53 + t59 * t55) * qJD(1), -t43 * t70 + t61 * t55 + t64, (t47 * t70 - t55 * t72) * pkin(5) + t57, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0; t79 * t55 + (r_i_i_C(3) * t55 + t59 * t53) * qJD(1), t38 + t61 * t53 + (t43 - t60) * t69, t38 + t78 * t53 + (-t60 - t77) * t69, -t55 * t63 + t38 + (-t44 * t69 - t53 * t73) * r_i_i_C(1), 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t79, t56, -t58, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:20
	% EndTime: 2020-04-15 19:46:20
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (570->68), mult. (520->101), div. (0->0), fcn. (384->10), ass. (0->63)
	t271 = qJ(2) + qJ(3);
	t269 = qJ(4) + t271;
	t264 = sin(t269);
	t275 = cos(qJ(5));
	t330 = r_i_i_C(1) * t275 + pkin(9);
	t291 = t330 * t264;
	t272 = sin(qJ(5));
	t312 = qJD(5) * t275;
	t265 = cos(t269);
	t270 = qJD(2) + qJD(3);
	t266 = qJD(4) + t270;
	t320 = t265 * t266;
	t332 = t264 * t312 + t272 * t320;
	t327 = pkin(11) + r_i_i_C(3);
	t305 = t327 * t265;
	t276 = cos(qJ(2));
	t321 = pkin(1) * qJD(2);
	t307 = t276 * t321;
	t267 = sin(t271);
	t325 = pkin(5) * t270;
	t311 = t267 * t325;
	t324 = pkin(9) * t264;
	t331 = (t305 - t324) * t266 - t307 - t311;
	t313 = qJD(5) * t272;
	t298 = t264 * t313;
	t328 = r_i_i_C(1) * t298 + t332 * r_i_i_C(2);
	t268 = cos(t271);
	t290 = t330 * t265;
	t306 = t327 * t264;
	t280 = (-t290 - t306) * t266 - t268 * t325;
	t326 = pkin(5) * t267;
	t322 = r_i_i_C(2) * t272;
	t274 = sin(qJ(1));
	t319 = t266 * t274;
	t318 = t266 * t275;
	t277 = cos(qJ(1));
	t317 = t266 * t277;
	t316 = t275 * t277;
	t315 = qJD(1) * t274;
	t314 = qJD(1) * t277;
	t309 = t264 * t322;
	t308 = qJD(1) * t322;
	t304 = t327 * t274;
	t303 = t264 * t318;
	t293 = qJD(5) * t265 - qJD(1);
	t292 = qJD(1) * t265 - qJD(5);
	t289 = t330 * t277;
	t288 = t328 * t277 + t315 * t291;
	t287 = t293 * t272;
	t286 = t277 * t264 * t308 + t328 * t274 + t314 * t305;
	t285 = -t305 - t309;
	t273 = sin(qJ(2));
	t284 = qJD(1) * (t273 * pkin(1) - pkin(5) * t268 - pkin(9) * t265 - pkin(15) - t306);
	t283 = t264 * t317 + t292 * t274;
	t281 = t273 * t321 + t280;
	t279 = -t265 * r_i_i_C(2) * t312 + (-t265 * t313 - t303) * r_i_i_C(1) + t327 * t320 + (-t324 + t309) * t266;
	t278 = t279 - t311;
	t263 = -t276 * pkin(1) - t326;
	t245 = -t292 * t316 + (t287 + t303) * t274;
	t244 = t293 * t275 * t274 + (-t264 * t319 + t292 * t277) * t272;
	t243 = t283 * t275 + t277 * t287;
	t242 = t283 * t272 - t293 * t316;
	t1 = [t245 * r_i_i_C(1) + t244 * r_i_i_C(2) - t331 * t274 + t277 * t284, (-t263 + t285) * t315 + t281 * t277 + t288, (t285 + t326) * t315 + t280 * t277 + t288, (-t274 * t308 - t327 * t317) * t264 + (-qJD(1) * t304 - t266 * t289) * t265 + t288, t242 * r_i_i_C(1) + t243 * r_i_i_C(2), 0, 0, 0, 0, 0, 0, 0, 0; -t243 * r_i_i_C(1) + t242 * r_i_i_C(2) + t274 * t284 + t331 * t277, (t263 - t291) * t314 + t281 * t274 + t286, (-t291 - t326) * t314 + t280 * t274 + t286, -t290 * t319 + (-qJD(1) * t289 - t266 * t304) * t264 + t286, -t244 * r_i_i_C(1) + t245 * r_i_i_C(2), 0, 0, 0, 0, 0, 0, 0, 0; 0, t278 - t307, t278, t279, (-t265 * t318 + t298) * r_i_i_C(2) - t332 * r_i_i_C(1), 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:18
	% EndTime: 2020-04-15 19:46:18
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (17->14), mult. (60->29), div. (0->0), fcn. (38->4), ass. (0->12)
	t17 = sin(qJ(1));
	t26 = qJD(1) * t17;
	t19 = cos(qJ(1));
	t25 = qJD(1) * t19;
	t24 = qJD(6) * t17;
	t23 = qJD(6) * t19;
	t16 = sin(qJ(6));
	t18 = cos(qJ(6));
	t22 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18;
	t21 = -r_i_i_C(1) * t18 + r_i_i_C(2) * t16 + pkin(14);
	t20 = t22 * qJD(6);
	t1 = [t22 * t24 + (-r_i_i_C(3) * t17 + t21 * t19) * qJD(1), 0, 0, 0, 0, (t16 * t23 + t18 * t26) * r_i_i_C(2) + (t16 * t26 - t18 * t23) * r_i_i_C(1), 0, 0, 0, 0, 0, 0, 0; -t19 * t20 + (r_i_i_C(3) * t19 + t21 * t17) * qJD(1), 0, 0, 0, 0, (t16 * t24 - t18 * t25) * r_i_i_C(2) + (-t16 * t25 - t18 * t24) * r_i_i_C(1), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:18
	% EndTime: 2020-04-15 19:46:19
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (77->20), mult. (110->34), div. (0->0), fcn. (71->6), ass. (0->22)
	t42 = qJ(2) + qJ(7);
	t39 = sin(t42);
	t40 = cos(t42);
	t66 = -r_i_i_C(1) * t40 + r_i_i_C(2) * t39;
	t65 = t66 * qJD(1);
	t41 = qJD(2) + qJD(7);
	t59 = t40 * t41;
	t60 = t39 * t41;
	t64 = r_i_i_C(1) * t60 + r_i_i_C(2) * t59;
	t45 = cos(qJ(2));
	t52 = qJD(2) * t45 * pkin(1);
	t63 = -t41 * t66 + t52;
	t58 = qJD(1) * t45;
	t43 = sin(qJ(2));
	t57 = qJD(2) * t43;
	t51 = -r_i_i_C(1) * t59 + r_i_i_C(2) * t60;
	t49 = t43 * pkin(1) + r_i_i_C(1) * t39 + r_i_i_C(2) * t40 - pkin(15);
	t44 = sin(qJ(1));
	t46 = cos(qJ(1));
	t48 = t64 * t44 + t65 * t46;
	t47 = -t65 * t44 + t64 * t46;
	t1 = [t63 * t44 + (-r_i_i_C(3) * t44 + t49 * t46) * qJD(1), (t44 * t58 + t46 * t57) * pkin(1) + t47, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0; -t63 * t46 + (r_i_i_C(3) * t46 + t49 * t44) * qJD(1), (t44 * t57 - t46 * t58) * pkin(1) + t48, 0, 0, 0, 0, t48, 0, 0, 0, 0, 0, 0; 0, t51 - t52, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:18
	% EndTime: 2020-04-15 19:46:18
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (68->17), mult. (88->29), div. (0->0), fcn. (56->4), ass. (0->16)
	t29 = qJD(2) + qJD(8);
	t31 = sin(qJ(1));
	t38 = t29 * t31;
	t32 = cos(qJ(1));
	t37 = t29 * t32;
	t36 = qJD(1) * t31;
	t35 = qJD(1) * t32;
	t30 = qJ(2) + qJ(8);
	t27 = sin(t30);
	t28 = cos(t30);
	t34 = r_i_i_C(1) * t28 - r_i_i_C(2) * t27;
	t33 = r_i_i_C(1) * t27 + r_i_i_C(2) * t28 - pkin(15);
	t26 = t34 * t29;
	t25 = (t27 * t35 + t28 * t38) * r_i_i_C(2) + (t27 * t38 - t28 * t35) * r_i_i_C(1);
	t24 = (-t27 * t36 + t28 * t37) * r_i_i_C(2) + (t27 * t37 + t28 * t36) * r_i_i_C(1);
	t1 = [t34 * t38 + (-r_i_i_C(3) * t31 + t33 * t32) * qJD(1), t24, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0; -t32 * t26 + (r_i_i_C(3) * t32 + t33 * t31) * qJD(1), t25, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0; 0, -t26, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiaD_transl_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:18
	% EndTime: 2020-04-15 19:46:18
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (187->26), mult. (152->37), div. (0->0), fcn. (98->6), ass. (0->30)
	t50 = qJ(2) + qJ(8);
	t47 = cos(t50);
	t49 = qJD(2) + qJD(8);
	t45 = qJD(9) + t49;
	t48 = qJ(9) + t50;
	t43 = sin(t48);
	t64 = r_i_i_C(2) * t43;
	t44 = cos(t48);
	t65 = r_i_i_C(1) * t44;
	t57 = (-t64 + t65) * t45;
	t67 = pkin(2) * t49;
	t70 = -t47 * t67 + t57;
	t68 = pkin(2) * t47;
	t66 = r_i_i_C(1) * t43;
	t63 = r_i_i_C(2) * t44;
	t51 = sin(qJ(1));
	t62 = t45 * t51;
	t52 = cos(qJ(1));
	t61 = t45 * t52;
	t60 = qJD(1) * t51;
	t59 = qJD(1) * t52;
	t55 = -t63 - t66;
	t46 = sin(t50);
	t54 = pkin(2) * t46 - pkin(15) + t55;
	t53 = t55 * t45 + t46 * t67;
	t41 = t59 * t65;
	t40 = t60 * t64;
	t37 = t41 + (-t64 - t68) * t59 + t53 * t51;
	t36 = t40 + (-t65 + t68) * t60 + t53 * t52;
	t1 = [-t70 * t51 + (-r_i_i_C(3) * t51 + t54 * t52) * qJD(1), t36, 0, 0, 0, 0, 0, t36, -t61 * t63 + t40 + (-t43 * t61 - t44 * t60) * r_i_i_C(1), 0, 0, 0, 0; t70 * t52 + (r_i_i_C(3) * t52 + t54 * t51) * qJD(1), t37, 0, 0, 0, 0, 0, t37, -t62 * t66 + t41 + (-t43 * t59 - t44 * t62) * r_i_i_C(2), 0, 0, 0, 0; 0, t70, 0, 0, 0, 0, 0, t70, t57, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiaD_transl_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 19:46:18
	% EndTime: 2020-04-15 19:46:19
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (248->35), mult. (166->47), div. (0->0), fcn. (107->8), ass. (0->36)
	t61 = cos(qJ(2));
	t58 = -qJD(7) - qJD(2);
	t56 = -qJD(10) + t58;
	t57 = pkin(19) - qJ(7) - qJ(2);
	t55 = -qJ(10) + t57;
	t51 = sin(t55);
	t81 = r_i_i_C(2) * t51;
	t52 = cos(t55);
	t82 = r_i_i_C(1) * t52;
	t67 = t81 + t82;
	t65 = t67 * t56;
	t54 = cos(t57);
	t83 = pkin(4) * t54;
	t63 = t58 * t83 - t65;
	t77 = pkin(1) * qJD(2);
	t85 = -t61 * t77 + t63;
	t79 = t51 * t56;
	t73 = r_i_i_C(1) * t79;
	t53 = sin(t57);
	t78 = t53 * t58;
	t84 = pkin(4) * t78 - t73;
	t80 = r_i_i_C(2) * t52;
	t60 = sin(qJ(1));
	t76 = qJD(1) * t60;
	t62 = cos(qJ(1));
	t75 = qJD(1) * t62;
	t72 = t56 * t80;
	t70 = qJD(1) * t81;
	t71 = t60 * t72 + t62 * t70 + t75 * t82;
	t59 = sin(qJ(2));
	t68 = t59 * t77 + t84;
	t66 = t59 * pkin(1) - pkin(4) * t53 + r_i_i_C(1) * t51 - pkin(15) - t80;
	t64 = -t60 * t73 + t71;
	t49 = -t61 * pkin(1) - t83;
	t45 = t62 * t72;
	t1 = [-t85 * t60 + (-r_i_i_C(3) * t60 + t66 * t62) * qJD(1), t45 + t68 * t62 + (-t49 - t67) * t76, 0, 0, 0, 0, t45 + t84 * t62 + (-t67 + t83) * t76, 0, 0, -t60 * t70 + t45 + (-t52 * t76 - t62 * t79) * r_i_i_C(1), 0, 0, 0; t85 * t62 + (r_i_i_C(3) * t62 + t66 * t60) * qJD(1), t49 * t75 + t68 * t60 + t71, 0, 0, 0, 0, (-t54 * t75 + t60 * t78) * pkin(4) + t64, 0, 0, t64, 0, 0, 0; 0, t85, 0, 0, 0, 0, t63, 0, 0, -t65, 0, 0, 0;];
	JaD_transl = t1;
end