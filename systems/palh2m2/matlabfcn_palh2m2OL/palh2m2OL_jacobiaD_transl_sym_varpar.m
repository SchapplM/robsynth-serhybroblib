% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh2m2OL
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = palh2m2OL_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh2m2OL_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m2OL_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_jacobiaD_transl_sym_varpar: pkin has to be [5x1] (double)');
JaD_transl=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
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
	t22 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18;
	t21 = -r_i_i_C(1) * t18 + r_i_i_C(2) * t16 - pkin(1);
	t20 = t22 * qJD(2);
	t1 = [t22 * t24 + (-r_i_i_C(3) * t17 + t19 * t21) * qJD(1), (t16 * t23 + t18 * t26) * r_i_i_C(2) + (t16 * t26 - t18 * t23) * r_i_i_C(1), 0, 0, 0, 0; -t19 * t20 + (r_i_i_C(3) * t19 + t17 * t21) * qJD(1), (t16 * t24 - t18 * t25) * r_i_i_C(2) + (-t16 * t25 - t18 * t24) * r_i_i_C(1), 0, 0, 0, 0; 0, -t20, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (77->25), mult. (110->37), div. (0->0), fcn. (71->6), ass. (0->26)
	t36 = qJD(2) + qJD(3);
	t37 = qJ(2) + qJ(3);
	t35 = cos(t37);
	t55 = r_i_i_C(2) * t35;
	t34 = sin(t37);
	t57 = r_i_i_C(1) * t34;
	t46 = t55 + t57;
	t44 = t46 * t36;
	t38 = sin(qJ(2));
	t58 = pkin(4) * t38;
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
	t45 = -t40 * pkin(4) - r_i_i_C(1) * t35 - pkin(1) + t56;
	t43 = t39 * t48 + t53 * t57 + (t49 - t50) * t41;
	t29 = t39 * t49;
	t1 = [t59 * t39 + (-r_i_i_C(3) * t39 + t45 * t41) * qJD(1), (t38 * t53 - t41 * t51) * pkin(4) + t43, t43, 0, 0, 0; -t59 * t41 + (r_i_i_C(3) * t41 + t45 * t39) * qJD(1), t29 + (-pkin(4) * t51 - t50) * t39 + (-t46 - t58) * t52, -t41 * t48 + t29 + (-t34 * t52 - t39 * t54) * r_i_i_C(1), 0, 0, 0; 0, -t59, -t44, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:58
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (196->35), mult. (166->48), div. (0->0), fcn. (107->8), ass. (0->36)
	t51 = sin(qJ(2));
	t50 = qJ(2) + qJ(3);
	t46 = sin(t50);
	t49 = qJD(2) + qJD(3);
	t45 = qJD(4) + t49;
	t48 = qJ(4) + t50;
	t44 = cos(t48);
	t74 = r_i_i_C(2) * t44;
	t43 = sin(t48);
	t76 = r_i_i_C(1) * t43;
	t59 = t74 + t76;
	t57 = t59 * t45;
	t70 = pkin(2) * t49;
	t55 = -t46 * t70 - t57;
	t72 = pkin(4) * qJD(2);
	t78 = -t51 * t72 + t55;
	t47 = cos(t50);
	t73 = t44 * t45;
	t67 = r_i_i_C(1) * t73;
	t77 = -t47 * t70 - t67;
	t75 = r_i_i_C(2) * t43;
	t71 = pkin(2) * t46;
	t52 = sin(qJ(1));
	t69 = qJD(1) * t52;
	t54 = cos(qJ(1));
	t68 = qJD(1) * t54;
	t66 = t45 * t75;
	t64 = qJD(1) * t74;
	t65 = t52 * t64 + t54 * t66 + t69 * t76;
	t53 = cos(qJ(2));
	t60 = -t53 * t72 + t77;
	t58 = -pkin(2) * t47 - pkin(4) * t53 - r_i_i_C(1) * t44 - pkin(1) + t75;
	t56 = -t54 * t67 + t65;
	t40 = -pkin(4) * t51 - t71;
	t37 = t52 * t66;
	t1 = [-t78 * t52 + (-r_i_i_C(3) * t52 + t58 * t54) * qJD(1), -t40 * t69 + t60 * t54 + t65, (-t47 * t49 * t54 + t46 * t69) * pkin(2) + t56, t56, 0, 0; t78 * t54 + (r_i_i_C(3) * t54 + t58 * t52) * qJD(1), t37 + t60 * t52 + (t40 - t59) * t68, t37 + t77 * t52 + (-t59 - t71) * t68, -t54 * t64 + t37 + (-t43 * t68 - t52 * t73) * r_i_i_C(1), 0, 0; 0, t78, t55, -t57, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:58
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (390->46), mult. (230->56), div. (0->0), fcn. (148->10), ass. (0->44)
	t63 = sin(qJ(2));
	t62 = qJ(2) + qJ(3);
	t58 = sin(t62);
	t61 = qJD(2) + qJD(3);
	t57 = qJD(4) + t61;
	t53 = qJD(5) + t57;
	t60 = qJ(4) + t62;
	t56 = qJ(5) + t60;
	t52 = cos(t56);
	t89 = r_i_i_C(2) * t52;
	t51 = sin(t56);
	t91 = r_i_i_C(1) * t51;
	t72 = t89 + t91;
	t70 = t72 * t53;
	t54 = sin(t60);
	t92 = pkin(5) * t54;
	t68 = -t57 * t92 - t70;
	t85 = pkin(2) * t61;
	t67 = -t58 * t85 + t68;
	t86 = pkin(4) * qJD(2);
	t94 = -t63 * t86 + t67;
	t88 = t52 * t53;
	t80 = r_i_i_C(1) * t88;
	t55 = cos(t60);
	t87 = t55 * t57;
	t93 = -pkin(5) * t87 - t80;
	t59 = cos(t62);
	t73 = -t59 * t85 + t93;
	t90 = r_i_i_C(2) * t51;
	t64 = sin(qJ(1));
	t84 = qJD(1) * t64;
	t66 = cos(qJ(1));
	t83 = qJD(1) * t66;
	t79 = t53 * t90;
	t77 = qJD(1) * t89;
	t78 = t64 * t77 + t66 * t79 + t84 * t91;
	t65 = cos(qJ(2));
	t74 = -t65 * t86 + t73;
	t48 = -pkin(2) * t58 - t92;
	t71 = -pkin(2) * t59 - t65 * pkin(4) - pkin(5) * t55 - r_i_i_C(1) * t52 - pkin(1) + t90;
	t69 = -t66 * t80 + t78;
	t46 = t64 * t79;
	t45 = -t63 * pkin(4) + t48;
	t1 = [-t94 * t64 + (-r_i_i_C(3) * t64 + t71 * t66) * qJD(1), -t45 * t84 + t74 * t66 + t78, -t48 * t84 + t73 * t66 + t78, (t54 * t84 - t66 * t87) * pkin(5) + t69, t69, 0; t94 * t66 + (r_i_i_C(3) * t66 + t71 * t64) * qJD(1), t46 + t74 * t64 + (t45 - t72) * t83, t46 + t73 * t64 + (t48 - t72) * t83, t46 + t93 * t64 + (-t72 - t92) * t83, -t66 * t77 + t46 + (-t51 * t83 - t64 * t88) * r_i_i_C(1), 0; 0, t94, t67, t68, -t70, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:59
	% EndTime: 2020-05-03 06:34:59
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (871->76), mult. (590->100), div. (0->0), fcn. (434->12), ass. (0->64)
	t265 = qJ(2) + qJ(3);
	t263 = qJ(4) + t265;
	t259 = qJ(5) + t263;
	t255 = cos(t259);
	t264 = qJD(2) + qJD(3);
	t260 = qJD(4) + t264;
	t256 = qJD(5) + t260;
	t266 = sin(qJ(6));
	t269 = cos(qJ(6));
	t254 = sin(t259);
	t305 = qJD(6) * t254;
	t321 = t255 * t256 * t266 + t269 * t305;
	t313 = r_i_i_C(3) * t255;
	t319 = r_i_i_C(1) * t269 + pkin(3);
	t320 = t254 * t319 + t313;
	t258 = cos(t263);
	t286 = t319 * t256;
	t282 = t255 * t286;
	t315 = pkin(5) * t260;
	t275 = -t258 * t315 - t282;
	t295 = t266 * t305;
	t318 = r_i_i_C(1) * t295 + t321 * r_i_i_C(2);
	t262 = cos(t265);
	t308 = pkin(2) * t264;
	t279 = -t262 * t308 + t275;
	t261 = sin(t265);
	t297 = t261 * t308;
	t267 = sin(qJ(2));
	t312 = pkin(4) * qJD(2);
	t301 = t267 * t312;
	t257 = sin(t263);
	t304 = t257 * t315;
	t317 = (pkin(3) * t254 + t313) * t256 + t297 + t301 + t304;
	t316 = pkin(5) * t257;
	t311 = t254 * t256;
	t310 = t256 * t269;
	t271 = cos(qJ(1));
	t309 = t269 * t271;
	t268 = sin(qJ(1));
	t307 = qJD(1) * t268;
	t306 = qJD(1) * t271;
	t302 = r_i_i_C(2) * t254 * t266;
	t300 = t268 * t311;
	t299 = t271 * t311;
	t290 = qJD(1) * t302;
	t288 = qJD(6) * t255 + qJD(1);
	t287 = qJD(1) * t255 + qJD(6);
	t285 = r_i_i_C(3) * t300 + t318 * t268 + t271 * t290;
	t251 = -pkin(2) * t261 - t316;
	t283 = t288 * t266;
	t281 = r_i_i_C(3) * t299 + t318 * t271 + t320 * t307;
	t270 = cos(qJ(2));
	t280 = -t270 * t312 + t279;
	t278 = qJD(1) * (-pkin(2) * t262 - pkin(3) * t255 - t270 * pkin(4) - pkin(5) * t258 + r_i_i_C(3) * t254 - pkin(1));
	t276 = t287 * t268 + t299;
	t274 = t256 * t302 - t254 * t286 + (-r_i_i_C(3) * t256 + (-r_i_i_C(1) * t266 - r_i_i_C(2) * t269) * qJD(6)) * t255;
	t273 = t274 - t304;
	t272 = t273 - t297;
	t241 = -t267 * pkin(4) + t251;
	t234 = -t287 * t309 + (t254 * t310 + t283) * t268;
	t233 = t288 * t269 * t268 + (t287 * t271 - t300) * t266;
	t232 = t276 * t269 + t271 * t283;
	t231 = t276 * t266 - t288 * t309;
	t1 = [t234 * r_i_i_C(1) + t233 * r_i_i_C(2) + t268 * t317 + t271 * t278, (-t241 - t302) * t307 + t280 * t271 + t281, (-t251 - t302) * t307 + t279 * t271 + t281, (-t302 + t316) * t307 + t275 * t271 + t281, -t268 * t290 - t271 * t282 + t281, t231 * r_i_i_C(1) + t232 * r_i_i_C(2); -t232 * r_i_i_C(1) + t231 * r_i_i_C(2) + t268 * t278 - t271 * t317, t280 * t268 + (t241 - t320) * t306 + t285, t279 * t268 + (t251 - t320) * t306 + t285, t275 * t268 + (-t320 - t316) * t306 + t285, -t268 * t282 - t306 * t320 + t285, -t233 * r_i_i_C(1) + t234 * r_i_i_C(2); 0, t272 - t301, t272, t273, t274, (-t255 * t310 + t295) * r_i_i_C(2) - t321 * r_i_i_C(1);];
	JaD_transl = t1;
end