% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh3m1OL
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% JaD_transl [3x10]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:16
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = palh3m1OL_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),uint8(0),zeros(3,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1OL_jacobiaD_transl_sym_varpar: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m1OL_jacobiaD_transl_sym_varpar: qJD has to be [10x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh3m1OL_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1OL_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1OL_jacobiaD_transl_sym_varpar: pkin has to be [16x1] (double)');
JaD_transl=NaN(3,10);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:35
	% EndTime: 2020-04-20 17:16:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:35
	% EndTime: 2020-04-20 17:16:36
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:36
	% EndTime: 2020-04-20 17:16:36
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (17->14), mult. (60->29), div. (0->0), fcn. (38->4), ass. (0->12)
	t17 = sin(qJ(1));
	t26 = qJD(1) * t17;
	t19 = cos(qJ(1));
	t25 = qJD(1) * t19;
	t24 = qJD(2) * t17;
	t23 = qJD(2) * t19;
	t16 = sin(qJ(2));
	t18 = cos(qJ(2));
	t22 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18;
	t21 = -r_i_i_C(1) * t18 + r_i_i_C(2) * t16 - pkin(12);
	t20 = t22 * qJD(2);
	t1 = [t22 * t24 + (-r_i_i_C(3) * t17 + t21 * t19) * qJD(1), (t16 * t23 + t18 * t26) * r_i_i_C(2) + (t16 * t26 - t18 * t23) * r_i_i_C(1), 0, 0, 0, 0, 0, 0, 0, 0; -t19 * t20 + (r_i_i_C(3) * t19 + t21 * t17) * qJD(1), (t16 * t24 - t18 * t25) * r_i_i_C(2) + (-t16 * t25 - t18 * t24) * r_i_i_C(1), 0, 0, 0, 0, 0, 0, 0, 0; 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:36
	% EndTime: 2020-04-20 17:16:36
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (77->25), mult. (110->40), div. (0->0), fcn. (71->6), ass. (0->27)
	t46 = qJD(2) + qJD(3);
	t47 = qJ(2) + qJ(3);
	t45 = cos(t47);
	t65 = r_i_i_C(2) * t45;
	t44 = sin(t47);
	t66 = r_i_i_C(1) * t44;
	t54 = t65 + t66;
	t48 = sin(qJ(2));
	t67 = pkin(1) * t48;
	t55 = qJD(2) * t67;
	t68 = t54 * t46 - t55;
	t64 = t44 * t46;
	t63 = t45 * t46;
	t62 = r_i_i_C(1) * t64 + r_i_i_C(2) * t63;
	t49 = sin(qJ(1));
	t61 = qJD(1) * t49;
	t51 = cos(qJ(1));
	t60 = qJD(1) * t51;
	t50 = cos(qJ(2));
	t59 = qJD(2) * t50;
	t58 = r_i_i_C(1) * t63;
	t57 = r_i_i_C(2) * t64;
	t56 = qJD(1) * t66;
	t53 = -t50 * pkin(1) + r_i_i_C(1) * t45 - r_i_i_C(2) * t44 - pkin(12);
	t52 = t51 * t56 + t60 * t65 + (-t57 + t58) * t49;
	t38 = t51 * t58;
	t1 = [-t68 * t49 + (-r_i_i_C(3) * t49 + t53 * t51) * qJD(1), t38 + (-pkin(1) * t59 - t57) * t51 + (-t54 + t67) * t61, -t49 * t56 + t38 + (-t45 * t61 - t51 * t64) * r_i_i_C(2), 0, 0, 0, 0, 0, 0, 0; t68 * t51 + (r_i_i_C(3) * t51 + t53 * t49) * qJD(1), (-t48 * t60 - t49 * t59) * pkin(1) + t52, t52, 0, 0, 0, 0, 0, 0, 0; 0, -t55 + t62, t62, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:36
	% EndTime: 2020-04-20 17:16:36
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (196->37), mult. (166->50), div. (0->0), fcn. (107->8), ass. (0->37)
	t60 = qJ(2) + qJ(3);
	t58 = qJ(4) + t60;
	t53 = sin(t58);
	t59 = qJD(2) + qJD(3);
	t55 = qJD(4) + t59;
	t80 = t53 * t55;
	t71 = r_i_i_C(2) * t80;
	t57 = cos(t60);
	t78 = t57 * t59;
	t85 = pkin(4) * t78 - t71;
	t56 = sin(t60);
	t83 = pkin(4) * t56;
	t52 = t59 * t83;
	t61 = sin(qJ(2));
	t76 = pkin(1) * qJD(2);
	t42 = -t61 * t76 + t52;
	t54 = cos(t58);
	t81 = r_i_i_C(2) * t54;
	t82 = r_i_i_C(1) * t53;
	t67 = t81 + t82;
	t84 = t67 * t55 + t42;
	t79 = t54 * t55;
	t77 = r_i_i_C(1) * t80 + r_i_i_C(2) * t79;
	t62 = sin(qJ(1));
	t75 = qJD(1) * t62;
	t64 = cos(qJ(1));
	t74 = qJD(1) * t64;
	t72 = r_i_i_C(1) * t79;
	t69 = qJD(1) * t82;
	t70 = t62 * t72 + t64 * t69 + t74 * t81;
	t63 = cos(qJ(2));
	t68 = -t63 * t76 + t85;
	t66 = -t63 * pkin(1) + pkin(4) * t57 + r_i_i_C(1) * t54 - r_i_i_C(2) * t53 - pkin(12);
	t65 = -t62 * t71 + t70;
	t51 = -t61 * pkin(1) + t83;
	t45 = t64 * t72;
	t1 = [-t84 * t62 + (-r_i_i_C(3) * t62 + t66 * t64) * qJD(1), t45 + t68 * t64 + (-t51 - t67) * t75, t45 + t85 * t64 + (-t67 - t83) * t75, -t62 * t69 + t45 + (-t54 * t75 - t64 * t80) * r_i_i_C(2), 0, 0, 0, 0, 0, 0; t84 * t64 + (r_i_i_C(3) * t64 + t66 * t62) * qJD(1), t51 * t74 + t68 * t62 + t70, (t56 * t74 + t62 * t78) * pkin(4) + t65, t65, 0, 0, 0, 0, 0, 0; 0, t42 + t77, t52 + t77, t77, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:37
	% EndTime: 2020-04-20 17:16:38
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (570->73), mult. (520->104), div. (0->0), fcn. (384->10), ass. (0->66)
	t313 = pkin(10) + r_i_i_C(3);
	t263 = qJD(2) + qJD(3);
	t264 = qJ(2) + qJ(3);
	t260 = sin(t264);
	t312 = pkin(4) * t260;
	t256 = t263 * t312;
	t259 = qJD(4) + t263;
	t262 = qJ(4) + t264;
	t258 = cos(t262);
	t296 = t313 * t258;
	t266 = sin(qJ(2));
	t309 = pkin(1) * qJD(2);
	t297 = t266 * t309;
	t257 = sin(t262);
	t311 = pkin(8) * t257;
	t317 = (t296 - t311) * t259 - t256 + t297;
	t267 = sin(qJ(1));
	t282 = qJD(1) * t258 + qJD(5);
	t270 = cos(qJ(1));
	t306 = t259 * t270;
	t293 = t257 * t306;
	t314 = t282 * t267 + t293;
	t310 = pkin(8) * t258;
	t308 = t259 * t267;
	t268 = cos(qJ(5));
	t307 = t259 * t268;
	t261 = cos(t264);
	t305 = t261 * t263;
	t304 = qJD(1) * t267;
	t303 = qJD(1) * t270;
	t265 = sin(qJ(5));
	t302 = qJD(5) * t265;
	t301 = qJD(5) * t268;
	t300 = pkin(4) * t305;
	t299 = r_i_i_C(2) * t257 * t265;
	t298 = r_i_i_C(2) * t301;
	t295 = t257 * t308;
	t294 = t257 * t307;
	t292 = t258 * t259 * t265;
	t291 = t258 * t308;
	t290 = t258 * t307;
	t289 = t257 * t303;
	t287 = t257 * t302;
	t286 = r_i_i_C(1) * t290;
	t285 = r_i_i_C(2) * t292;
	t284 = r_i_i_C(1) * t287;
	t283 = qJD(5) * t258 + qJD(1);
	t281 = t283 * t270;
	t280 = t282 * t270;
	t279 = t268 * r_i_i_C(1) * t289 + t267 * t286 + t313 * t295 + (t289 + t291) * pkin(8);
	t278 = -t296 - t299;
	t277 = t257 * t301 + t292;
	t269 = cos(qJ(2));
	t276 = qJD(1) * (-t269 * pkin(1) + pkin(4) * t261 + t313 * t257 - pkin(12) + t310);
	t275 = t258 * t298 + (t278 + t311) * t259 + (t258 * t302 + t294) * r_i_i_C(1);
	t274 = -t277 * r_i_i_C(2) - t284;
	t273 = t256 + t275;
	t272 = t270 * t286 + t306 * t310 + t299 * t304 + (-t270 * t298 - pkin(8) * t304 + (-t268 * t304 - t270 * t302) * r_i_i_C(1)) * t257 + t313 * (t258 * t304 + t293);
	t271 = -t270 * t285 + t272;
	t255 = -t266 * pkin(1) + t312;
	t238 = -t269 * t309 + t300;
	t234 = t268 * t280 + (-t283 * t265 - t294) * t267;
	t233 = t283 * t268 * t267 + (t280 - t295) * t265;
	t232 = t265 * t281 + t314 * t268;
	t231 = -t314 * t265 + t268 * t281;
	t1 = [t234 * r_i_i_C(1) - t233 * r_i_i_C(2) + t317 * t267 + t270 * t276, -t255 * t304 + (t238 - t285) * t270 + t272, (-t260 * t304 + t270 * t305) * pkin(4) + t271, t271, t231 * r_i_i_C(1) - t232 * r_i_i_C(2), 0, 0, 0, 0, 0; t232 * r_i_i_C(1) + t231 * r_i_i_C(2) + t267 * t276 - t317 * t270, (t238 + t274) * t267 + (t255 + t278) * t303 + t279, (t274 + t300) * t267 + (t278 + t312) * t303 + t279, -t267 * t284 - t296 * t303 + (-t265 * t291 + (-t265 * t303 - t267 * t301) * t257) * r_i_i_C(2) + t279, t233 * r_i_i_C(1) + t234 * r_i_i_C(2), 0, 0, 0, 0, 0; 0, t273 - t297, t273, t275, (-t287 + t290) * r_i_i_C(2) + t277 * r_i_i_C(1), 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:36
	% EndTime: 2020-04-20 17:16:36
	% DurationCPUTime: 0.07s
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
	t21 = -r_i_i_C(1) * t18 + r_i_i_C(2) * t16 + pkin(6);
	t20 = t22 * qJD(6);
	t1 = [t22 * t24 + (-r_i_i_C(3) * t17 + t21 * t19) * qJD(1), 0, 0, 0, 0, (t16 * t23 + t18 * t26) * r_i_i_C(2) + (t16 * t26 - t18 * t23) * r_i_i_C(1), 0, 0, 0, 0; -t19 * t20 + (r_i_i_C(3) * t19 + t21 * t17) * qJD(1), 0, 0, 0, 0, (t16 * t24 - t18 * t25) * r_i_i_C(2) + (-t16 * t25 - t18 * t24) * r_i_i_C(1), 0, 0, 0, 0; 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:36
	% EndTime: 2020-04-20 17:16:36
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (77->25), mult. (110->37), div. (0->0), fcn. (71->6), ass. (0->26)
	t36 = qJD(2) + qJD(7);
	t37 = qJ(2) + qJ(7);
	t35 = cos(t37);
	t55 = r_i_i_C(2) * t35;
	t34 = sin(t37);
	t57 = r_i_i_C(1) * t34;
	t46 = t55 + t57;
	t44 = t46 * t36;
	t38 = sin(qJ(2));
	t58 = pkin(1) * t38;
	t59 = qJD(2) * t58 + t44;
	t56 = r_i_i_C(2) * t34;
	t54 = t35 * t36;
	t39 = sin(qJ(1));
	t53 = qJD(1) * t39;
	t41 = cos(qJ(1));
	t52 = qJD(1) * t41;
	t40 = cos(qJ(2));
	t51 = qJD(2) * t40;
	t50 = r_i_i_C(1) * t54;
	t49 = t36 * t56;
	t48 = qJD(1) * t55;
	t45 = -t40 * pkin(1) - r_i_i_C(1) * t35 - pkin(12) + t56;
	t43 = t39 * t48 + t53 * t57 + (t49 - t50) * t41;
	t29 = t39 * t49;
	t1 = [t59 * t39 + (-r_i_i_C(3) * t39 + t45 * t41) * qJD(1), (t38 * t53 - t41 * t51) * pkin(1) + t43, 0, 0, 0, 0, t43, 0, 0, 0; -t59 * t41 + (r_i_i_C(3) * t41 + t45 * t39) * qJD(1), t29 + (-pkin(1) * t51 - t50) * t39 + (-t46 - t58) * t52, 0, 0, 0, 0, -t41 * t48 + t29 + (-t34 * t52 - t39 * t54) * r_i_i_C(1), 0, 0, 0; 0, -t59, 0, 0, 0, 0, -t44, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:36
	% EndTime: 2020-04-20 17:16:36
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (248->36), mult. (166->46), div. (0->0), fcn. (107->8), ass. (0->34)
	t51 = sin(qJ(2));
	t49 = pkin(15) - qJ(7) - qJ(2);
	t45 = sin(t49);
	t50 = -qJD(7) - qJD(2);
	t48 = -qJD(8) + t50;
	t47 = -qJ(8) + t49;
	t44 = cos(t47);
	t71 = r_i_i_C(2) * t44;
	t43 = sin(t47);
	t74 = r_i_i_C(1) * t43;
	t62 = (-t71 + t74) * t48;
	t75 = pkin(3) * t50;
	t57 = -t45 * t75 + t62;
	t68 = pkin(1) * qJD(2);
	t77 = -t51 * t68 + t57;
	t46 = cos(t49);
	t72 = r_i_i_C(2) * t43;
	t73 = r_i_i_C(1) * t44;
	t61 = -t72 - t73;
	t55 = t46 * t75 + t61 * t48;
	t76 = pkin(3) * t45;
	t52 = sin(qJ(1));
	t70 = t48 * t52;
	t54 = cos(qJ(1));
	t69 = t48 * t54;
	t67 = qJD(1) * t52;
	t66 = qJD(1) * t54;
	t53 = cos(qJ(2));
	t59 = -t53 * pkin(1) - pkin(3) * t46 - pkin(12) - t61;
	t56 = -t53 * t68 + t55;
	t42 = -t51 * pkin(1) + t76;
	t41 = t66 * t71;
	t40 = t67 * t74;
	t1 = [-t77 * t52 + (-r_i_i_C(3) * t52 + t59 * t54) * qJD(1), t40 + (-t42 - t71) * t67 + t56 * t54, 0, 0, 0, 0, t40 + (-t71 - t76) * t67 + t55 * t54, -t69 * t73 + t40 + (-t43 * t69 - t44 * t67) * r_i_i_C(2), 0, 0; t77 * t54 + (r_i_i_C(3) * t54 + t59 * t52) * qJD(1), t41 + (t42 - t74) * t66 + t56 * t52, 0, 0, 0, 0, t41 + (-t74 + t76) * t66 + t55 * t52, -t70 * t72 + t41 + (-t43 * t66 - t44 * t70) * r_i_i_C(1), 0, 0; 0, t77, 0, 0, 0, 0, t57, t62, 0, 0;];
	JaD_transl = t1;
end