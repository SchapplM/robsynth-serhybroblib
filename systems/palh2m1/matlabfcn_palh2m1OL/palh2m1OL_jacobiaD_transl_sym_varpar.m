% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh2m1OL
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = palh2m1OL_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh2m1OL_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m1OL_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_jacobiaD_transl_sym_varpar: pkin has to be [6x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:40
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (17->13), mult. (60->29), div. (0->0), fcn. (38->4), ass. (0->12)
	t17 = sin(qJ(1));
	t26 = qJD(1) * t17;
	t19 = cos(qJ(1));
	t25 = qJD(1) * t19;
	t24 = qJD(2) * t17;
	t23 = qJD(2) * t19;
	t16 = sin(qJ(2));
	t18 = cos(qJ(2));
	t22 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18;
	t21 = -r_i_i_C(1) * t18 + r_i_i_C(2) * t16 - pkin(1);
	t20 = t22 * qJD(2);
	t1 = [t17 * t20 + (r_i_i_C(3) * t17 + t19 * t21) * qJD(1), (t16 * t23 + t18 * t26) * r_i_i_C(2) + (t16 * t26 - t18 * t23) * r_i_i_C(1), 0, 0, 0; -t22 * t23 + (-r_i_i_C(3) * t19 + t17 * t21) * qJD(1), (t16 * t24 - t18 * t25) * r_i_i_C(2) + (-t16 * t25 - t18 * t24) * r_i_i_C(1), 0, 0, 0; 0, t20, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (77->25), mult. (110->40), div. (0->0), fcn. (71->6), ass. (0->27)
	t41 = qJD(2) + qJD(3);
	t42 = qJ(2) + qJ(3);
	t40 = cos(t42);
	t60 = r_i_i_C(2) * t40;
	t39 = sin(t42);
	t61 = r_i_i_C(1) * t39;
	t49 = t60 + t61;
	t43 = sin(qJ(2));
	t62 = pkin(2) * t43;
	t50 = qJD(2) * t62;
	t63 = t41 * t49 + t50;
	t59 = t39 * t41;
	t58 = t40 * t41;
	t57 = r_i_i_C(1) * t59 + r_i_i_C(2) * t58;
	t44 = sin(qJ(1));
	t56 = qJD(1) * t44;
	t46 = cos(qJ(1));
	t55 = qJD(1) * t46;
	t45 = cos(qJ(2));
	t54 = qJD(2) * t45;
	t53 = r_i_i_C(1) * t58;
	t52 = r_i_i_C(2) * t59;
	t51 = qJD(1) * t60;
	t48 = -t45 * pkin(2) - r_i_i_C(1) * t40 + r_i_i_C(2) * t39 - pkin(1);
	t47 = t44 * t51 + t56 * t61 + (t52 - t53) * t46;
	t32 = t44 * t52;
	t1 = [t63 * t44 + (r_i_i_C(3) * t44 + t46 * t48) * qJD(1), (t43 * t56 - t46 * t54) * pkin(2) + t47, t47, 0, 0; -t63 * t46 + (-r_i_i_C(3) * t46 + t44 * t48) * qJD(1), t32 + (-pkin(2) * t54 - t53) * t44 + (-t49 - t62) * t55, -t46 * t51 + t32 + (-t39 * t55 - t44 * t58) * r_i_i_C(1), 0, 0; 0, t50 + t57, t57, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (196->37), mult. (166->50), div. (0->0), fcn. (107->8), ass. (0->37)
	t57 = qJ(2) + qJ(3);
	t55 = qJ(4) + t57;
	t51 = cos(t55);
	t56 = qJD(2) + qJD(3);
	t52 = qJD(4) + t56;
	t76 = t51 * t52;
	t69 = r_i_i_C(1) * t76;
	t54 = cos(t57);
	t75 = t54 * t56;
	t82 = -pkin(3) * t75 - t69;
	t53 = sin(t57);
	t80 = pkin(3) * t53;
	t49 = t56 * t80;
	t58 = sin(qJ(2));
	t73 = pkin(2) * qJD(2);
	t39 = -t58 * t73 - t49;
	t78 = r_i_i_C(2) * t51;
	t50 = sin(t55);
	t79 = r_i_i_C(1) * t50;
	t64 = t78 + t79;
	t81 = t64 * t52 - t39;
	t77 = t50 * t52;
	t74 = r_i_i_C(1) * t77 + r_i_i_C(2) * t76;
	t59 = sin(qJ(1));
	t72 = qJD(1) * t59;
	t61 = cos(qJ(1));
	t71 = qJD(1) * t61;
	t68 = r_i_i_C(2) * t77;
	t66 = qJD(1) * t78;
	t67 = t59 * t66 + t61 * t68 + t72 * t79;
	t60 = cos(qJ(2));
	t65 = -t60 * t73 + t82;
	t63 = -pkin(2) * t60 - pkin(3) * t54 - r_i_i_C(1) * t51 + r_i_i_C(2) * t50 - pkin(1);
	t62 = -t61 * t69 + t67;
	t48 = -pkin(2) * t58 - t80;
	t41 = t59 * t68;
	t1 = [t81 * t59 + (r_i_i_C(3) * t59 + t63 * t61) * qJD(1), -t48 * t72 + t65 * t61 + t67, (t53 * t72 - t61 * t75) * pkin(3) + t62, t62, 0; -t81 * t61 + (-r_i_i_C(3) * t61 + t63 * t59) * qJD(1), t41 + t65 * t59 + (t48 - t64) * t71, t41 + t82 * t59 + (-t64 - t80) * t71, -t61 * t66 + t41 + (-t50 * t71 - t59 * t76) * r_i_i_C(1), 0; 0, -t39 + t74, t49 + t74, t74, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:42
	% EndTime: 2020-05-03 00:53:43
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (570->67), mult. (520->100), div. (0->0), fcn. (384->10), ass. (0->61)
	t269 = qJ(2) + qJ(3);
	t267 = qJ(4) + t269;
	t262 = sin(t267);
	t273 = cos(qJ(5));
	t325 = r_i_i_C(1) * t273 + pkin(4);
	t289 = t325 * t262;
	t263 = cos(t267);
	t268 = qJD(2) + qJD(3);
	t264 = qJD(4) + t268;
	t270 = sin(qJ(5));
	t308 = qJD(5) * t273;
	t327 = t263 * t264 * t270 + t262 * t308;
	t322 = pkin(6) + r_i_i_C(3);
	t303 = t322 * t263;
	t265 = sin(t269);
	t320 = pkin(3) * t268;
	t261 = t265 * t320;
	t271 = sin(qJ(2));
	t316 = pkin(2) * qJD(2);
	t305 = t271 * t316;
	t319 = pkin(4) * t262;
	t326 = (t303 - t319) * t264 - t261 - t305;
	t309 = qJD(5) * t270;
	t296 = t262 * t309;
	t323 = r_i_i_C(1) * t296 + t327 * r_i_i_C(2);
	t266 = cos(t269);
	t288 = t325 * t263;
	t304 = t322 * t262;
	t276 = (-t288 - t304) * t264 - t266 * t320;
	t321 = pkin(3) * t265;
	t317 = r_i_i_C(2) * t270;
	t272 = sin(qJ(1));
	t315 = t264 * t272;
	t314 = t264 * t273;
	t275 = cos(qJ(1));
	t313 = t264 * t275;
	t312 = t273 * t275;
	t311 = qJD(1) * t272;
	t310 = qJD(1) * t275;
	t306 = qJD(1) * t317;
	t302 = t322 * t272;
	t301 = t262 * t314;
	t291 = qJD(5) * t263 + qJD(1);
	t290 = qJD(1) * t263 + qJD(5);
	t287 = t325 * t275;
	t286 = t323 * t275 + t311 * t289;
	t285 = t291 * t270;
	t284 = t275 * t262 * t306 + t323 * t272 + t310 * t303;
	t283 = -t262 * t317 - t303;
	t274 = cos(qJ(2));
	t282 = qJD(1) * (-t274 * pkin(2) - pkin(3) * t266 - pkin(4) * t263 - pkin(1) - t304);
	t281 = t262 * t313 + t290 * t272;
	t279 = -t274 * t316 + t276;
	t278 = t263 * r_i_i_C(2) * t308 + (t283 + t319) * t264 + (t263 * t309 + t301) * r_i_i_C(1);
	t277 = t261 + t278;
	t260 = -t271 * pkin(2) - t321;
	t241 = -t290 * t312 + (t285 + t301) * t272;
	t240 = t291 * t273 * t272 + (-t262 * t315 + t290 * t275) * t270;
	t239 = t281 * t273 + t275 * t285;
	t238 = t281 * t270 - t291 * t312;
	t1 = [t241 * r_i_i_C(1) + t240 * r_i_i_C(2) - t326 * t272 + t275 * t282, (-t260 + t283) * t311 + t279 * t275 + t286, (t283 + t321) * t311 + t276 * t275 + t286, (-t272 * t306 - t322 * t313) * t262 + (-qJD(1) * t302 - t264 * t287) * t263 + t286, t238 * r_i_i_C(1) + t239 * r_i_i_C(2); -t239 * r_i_i_C(1) + t238 * r_i_i_C(2) + t272 * t282 + t326 * t275, (t260 - t289) * t310 + t279 * t272 + t284, (-t289 - t321) * t310 + t276 * t272 + t284, -t288 * t315 + (-qJD(1) * t287 - t264 * t302) * t262 + t284, -t240 * r_i_i_C(1) + t241 * r_i_i_C(2); 0, t277 + t305, t277, t278, (t263 * t314 - t296) * r_i_i_C(2) + t327 * r_i_i_C(1);];
	JaD_transl = t1;
end